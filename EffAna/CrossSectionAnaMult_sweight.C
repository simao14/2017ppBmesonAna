#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>
#include "../henri2022/parameter.h" 
//#include "tnp_weight_lowptPbPb.h"

//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, 
		std::vector<std::vector<double> > numbers , std::string caption, int useerr, std::vector<std::vector<double> > error ={{0}}){
	
	std::ofstream file_check;
	std::ofstream file;

	//Begin Document
                                                                                    
	file.open(filename + ".tex");
	file_check.open(filename + "_check.tex");


	file_check << "\\documentclass{article}" << std::endl;
	//file << "\\usepackage[utf8]{inputenc}" << std::endl;     
	file_check << "\\usepackage{rotating}" << std::endl;                                                                                   
	// file_check << "\\usepackage{cancel}" << std::endl;
	file_check << "\\usepackage{geometry}" << std::endl;
	file_check << "\\usepackage{booktabs}" << std::endl;
	file_check << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;

	file_check << "\\begin{document}" << std::endl;
	// Create table                                                                                                                                                                                                                                                
	std::string col="c";
	for(int i=1; i<n_col; i++)
		col+="|c";

		file_check << "\\begin{sidewaystable}"<< std::endl;
		file_check << "\\begin{tabular}{"+col+"}" << std::endl;
		file_check << "\\toprule" << std::endl;
		file << "\\begin{tabular}{"+col+"}" << std::endl;
		file << "\\toprule" << std::endl;
	
	for(int c=0; c<n_col-1; c++){
		file << col_name[c] << " & ";
		file_check << col_name[c] << " & ";
	}
	
	file << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;
	file_check << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;
	if (useerr==0){
		for(int i=1; i<n_lin; i++)
		{	
			file << labels[i-1] << " & ";
			file_check << labels[i-1] << " & ";
			for(int c=1; c<n_col-1; c++){
				file << std::setprecision(3)<<  numbers[i-1][c-1]<< " \\% & ";
				file_check << std::setprecision(3)<<  numbers[i-1][c-1]<< " \\% & ";
										}
				file << std::setprecision(3)<<  numbers[i-1][n_col-2]<< " \\% \\\\" << std::endl;
				file_check << std::setprecision(3)<<  numbers[i-1][n_col-2]<< " \\% \\\\" << std::endl; 
		}
	}
	if (useerr==1) {
		for(int i=1; i<n_lin; i++)
		{	
			file << labels[i-1] << " & ";
			file_check << labels[i-1] << " & ";
			for(int c=1; c<n_col-1; c++){
				file << std::setprecision(3)<<  numbers[i-1][c-1] << "$\\pm$" <<std::setprecision(1)<< error[i-1][c-1]<< " \\% & ";
				file_check << std::setprecision(3)<<  numbers[i-1][c-1] << "$\\pm$" <<std::setprecision(1)<< error[i-1][c-1]<< " \\% & ";
										}
				file << std::setprecision(3)<<  numbers[i-1][n_col-2] << "$\\pm$" <<std::setprecision(1)<< error[i-1][n_col-2]<< " \\% \\\\" << std::endl;
				file_check << std::setprecision(3)<<  numbers[i-1][n_col-2] << "$\\pm$" << std::setprecision(1)<<error[i-1][n_col-2]<< " \\% \\\\" << std::endl; 
		}
	}
	if (useerr==2) {
		for(int i=1; i<n_lin; i++)
		{	
			file << labels[i-1] << " & ";
			file_check << labels[i-1] << " & ";
			for(int c=1; c<n_col-1; c++){
				file <<  numbers[i-1][c-1] << " & ";
				file_check <<  numbers[i-1][c-1] << " & ";
										}
				file <<  numbers[i-1][n_col-2] <<  " \\\\" << std::endl;
				file_check <<  numbers[i-1][n_col-2] <<  " \\\\" << std::endl; 
		}
	}

	file << "\\bottomrule" << std::endl;
	file_check << "\\bottomrule" << std::endl;

	//End Table                                                                                                                                    
	file << "\\end{tabular}" << std::endl;
	file_check << "\\end{tabular}" << std::endl;
	file_check << "\\caption{"+caption+"}" << std::endl;
	file_check << "\\end{sidewaystable}"<< std::endl;

	//file_check << "\\end{table}" << std::endl;
	//End document                                                                                                                                 

	file_check << "\\end{document}" << std::endl;

	file.close();
	file_check.close();
	
	system(("pdflatex " + filename+ "_check.tex").c_str());
	system(("open " + filename + "_check.pdf").c_str());
}


void CrossSectionAnaMult_sweight(int DoTnP,int whichvar,int meson_n, int BsBP=0, int usemc=0){

	TString var_n;
	TString var_N;
	TString var_n2;
	TString bsbpbins = "";
	if(BsBP==1 && meson_n==0 ) {bsbpbins = "_BsBPBINS";}

	double BRchain;
	int NCand;
	if (meson_n == 0 && BsBP == 0){
		var_n="BP";
		var_n2="BP";
		var_N="B^{+}";
		NCand = 10;
		BRchain = 6.02061e-5;
	} 
	else if (meson_n == 0 && BsBP == 1 ) {
		var_n="BP";
		var_n2="BPBs";
		var_N="B^{+}";
		NCand = 10;
		BRchain = 6.02061e-5;
	}
	else {
		var_n="Bs";
		var_n2="Bs";
		var_N="B^{0}_{s}";
		NCand = 15;
		BRchain = 3.1189e-5;
	}
	
	gSystem->mkdir( var_n2.Data() , true);
	gSystem->mkdir(Form("%s/EffFinal",var_n2.Data()) ,true);
	gSystem->mkdir(Form("%s/FinalFiles",var_n2.Data()) ,true);
	gSystem->mkdir(Form("%s/Trash",var_n2.Data()) ,true);

	int NBins;
	//const int NBins = 6;

	int TnP = 1;

	const double epsilon = 1e-7;

	

	//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);


	TString FileName;
	TFile * fin;
	TFile * fin_0;
	TFile * fin_1;
	TFile * fin_2;
	TFile * fin_3;
	TFile * fin_4;
	TFile * fin_5;
	TFile * fin_6;
	TFile * fin_7;
	TTree * EffInfoTree;
	TTree * EffInfoTree_0;
	TTree * EffInfoTree_1;
	TTree * EffInfoTree_2;
	TTree * EffInfoTree_3;
	TTree * EffInfoTree_4;
	TTree * EffInfoTree_5;
	TTree * EffInfoTree_6;
	TTree * EffInfoTree_7;
	TTree * root;
	
	FileName =Form("/data3/smcosta/data/%sData_nom_BDT.root",var_n.Data());

	if (usemc==0){
		if (meson_n==0 && BsBP == 0){
			if (whichvar==0){
				fin_0 = new TFile("/data3/smcosta/presel/forSim/pt/5.0_7.0_ntKp_w_sWeights_Bpt.root");
				fin_1 = new TFile("/data3/smcosta/presel/forSim/pt/7.0_10.0_ntKp_w_sWeights_Bpt.root");
				fin_2 = new TFile("/data3/smcosta/presel/forSim/pt/10.0_15.0_ntKp_w_sWeights_Bpt.root");
				fin_3 = new TFile("/data3/smcosta/presel/forSim/pt/15.0_20.0_ntKp_w_sWeights_Bpt.root");
				fin_4 = new TFile("/data3/smcosta/presel/forSim/pt/20.0_30.0_ntKp_w_sWeights_Bpt.root");
				fin_5 = new TFile("/data3/smcosta/presel/forSim/pt/30.0_50.0_ntKp_w_sWeights_Bpt.root");
				fin_6 = new TFile("/data3/smcosta/presel/forSim/pt/50.0_60.0_ntKp_w_sWeights_Bpt.root");
			} 
			else if (whichvar == 1){
				fin_0 = new TFile("/data3/smcosta/presel/forSim/By/0.0_0.5_ntKp_w_sWeights_By.root");
				fin_1 = new TFile("/data3/smcosta/presel/forSim/By/0.5_1.0_ntKp_w_sWeights_By.root");
				fin_2 = new TFile("/data3/smcosta/presel/forSim/By/1.0_1.5_ntKp_w_sWeights_By.root");
				fin_3 = new TFile("/data3/smcosta/presel/forSim/By/1.5_2.0_ntKp_w_sWeights_By.root");
				fin_4 = new TFile("/data3/smcosta/presel/forSim/By/2.0_2.4_ntKp_w_sWeights_By.root");
			}
			else if (whichvar == 2){
				fin_0 = new TFile("/data3/smcosta/presel/forSim/nmul/0.0_20.0_ntKp_w_sWeights_nmult.root");
				fin_1 = new TFile("/data3/smcosta/presel/forSim/nmul/20.0_30.0_ntKp_w_sWeights_nmult.root");
				fin_2 = new TFile("/data3/smcosta/presel/forSim/nmul/30.0_40.0_ntKp_w_sWeights_nmult.root");
				fin_3 = new TFile("/data3/smcosta/presel/forSim/nmul/40.0_50.0_ntKp_w_sWeights_nmult.root");
				fin_4 = new TFile("/data3/smcosta/presel/forSim/nmul/50.0_60.0_ntKp_w_sWeights_nmult.root");
				fin_5 = new TFile("/data3/smcosta/presel/forSim/nmul/50.0_60.0_ntKp_w_sWeights_nmult.root");
				fin_6 = new TFile("/data3/smcosta/presel/forSim/nmul/60.0_70.0_ntKp_w_sWeights_nmult.root");
				fin_7 = new TFile("/data3/smcosta/presel/forSim/nmul/70.0_100.0_ntKp_w_sWeights_nmult.root");
			}

		} else if (meson_n ==0 && BsBP == 1){

			if (whichvar==0){
				fin_0 = new TFile("/data3/smcosta/presel/forSim/pt/7.0_10.0_ntKp_w_sWeights_Bpt.root");
				fin_1 = new TFile("/data3/smcosta/presel/forSim/pt/10.0_15.0_ntKp_w_sWeights_Bpt.root");
				fin_2 = new TFile("/data3/smcosta/presel/forSim/pt/15.0_20.0_ntKp_w_sWeights_Bpt.root");
				fin_3 = new TFile("/data3/smcosta/presel/forSim/pt/20.0_50.0_ntKp_w_sWeights_Bpt.root");
			} 
			else if (whichvar == 1){
				fin_0 = new TFile("/data3/smcosta/presel/forSim/By/0.0_0.5_ntKp_w_sWeights_By_BsBP.root");
				fin_1 = new TFile("/data3/smcosta/presel/forSim/By/0.5_1.0_ntKp_w_sWeights_By_BsBP.root");
				fin_2 = new TFile("/data3/smcosta/presel/forSim/By/1.0_1.5_ntKp_w_sWeights_By_BsBP.root");
				fin_3 = new TFile("/data3/smcosta/presel/forSim/By/1.5_2.0_ntKp_w_sWeights_By_BsBP.root");
				fin_4 = new TFile("/data3/smcosta/presel/forSim/By/2.0_2.4_ntKp_w_sWeights_By_BsBP.root");
			}
			else if (whichvar == 2){
				fin_0 = new TFile("/data3/smcosta/presel/forSim/nmul/0.0_20.0_ntKp_w_sWeights_nMult_BsBP.root");
				fin_1 = new TFile("/data3/smcosta/presel/forSim/nmul/20.0_30.0_ntKp_w_sWeights_nMult_BsBP.root");
				fin_2 = new TFile("/data3/smcosta/presel/forSim/nmul/30.0_40.0_ntKp_w_sWeights_nMult_BsBP.root");
				fin_3 = new TFile("/data3/smcosta/presel/forSim/nmul/40.0_50.0_ntKp_w_sWeights_nMult_BsBP.root");
				fin_4 = new TFile("/data3/smcosta/presel/forSim/nmul/50.0_60.0_ntKp_w_sWeights_nMult_BsBP.root");
				fin_5 = new TFile("/data3/smcosta/presel/forSim/nmul/50.0_60.0_ntKp_w_sWeights_nMult_BsBP.root");
				fin_6 = new TFile("/data3/smcosta/presel/forSim/nmul/60.0_70.0_ntKp_w_sWeights_nMult_BsBP.root");
				fin_7 = new TFile("/data3/smcosta/presel/forSim/nmul/70.0_100.0_ntKp_w_sWeights_nMult_BsBP.root");
			}
		} else {
			if (whichvar==0){
				fin_0 = new TFile("/data3/smcosta/presel/forSim/pt/7.0_10.0_ntphi_w_sWeights_Bpt.root");
				fin_1 = new TFile("/data3/smcosta/presel/forSim/pt/10.0_15.0_ntphi_w_sWeights_Bpt.root");
				fin_2 = new TFile("/data3/smcosta/presel/forSim/pt/15.0_20.0_ntphi_w_sWeights_Bpt.root");
				fin_3 = new TFile("/data3/smcosta/presel/forSim/pt/20.0_50.0_ntphi_w_sWeights_Bpt.root");
			} 
			else if (whichvar == 1){
				fin_0 = new TFile("/data3/smcosta/presel/forSim/By/0.0_0.5_ntphi_w_sWeights_By.root");
				fin_1 = new TFile("/data3/smcosta/presel/forSim/By/0.5_1.0_ntphi_w_sWeights_By.root");
				fin_2 = new TFile("/data3/smcosta/presel/forSim/By/1.0_1.5_ntphi_w_sWeights_By.root");
				fin_3 = new TFile("/data3/smcosta/presel/forSim/By/1.5_2.0_ntphi_w_sWeights_By.root");
				fin_4 = new TFile("/data3/smcosta/presel/forSim/By/2.0_2.4_ntphi_w_sWeights_By.root");
			}
			else if (whichvar == 2){
				fin_0 = new TFile("/data3/smcosta/presel/forSim/nmul/0.0_20.0_ntphi_w_sWeights_nmult.root");
				fin_1 = new TFile("/data3/smcosta/presel/forSim/nmul/20.0_30.0_ntphi_w_sWeights_nmult.root");
				fin_2 = new TFile("/data3/smcosta/presel/forSim/nmul/30.0_40.0_ntphi_w_sWeights_nmult.root");
				fin_3 = new TFile("/data3/smcosta/presel/forSim/nmul/40.0_50.0_ntphi_w_sWeights_nmult.root");
				fin_4 = new TFile("/data3/smcosta/presel/forSim/nmul/50.0_60.0_ntphi_w_sWeights_nmult.root");
				fin_5 = new TFile("/data3/smcosta/presel/forSim/nmul/50.0_60.0_ntphi_w_sWeights_nmult.root");
				fin_6 = new TFile("/data3/smcosta/presel/forSim/nmul/60.0_70.0_ntphi_w_sWeights_nmult.root");
				fin_7 = new TFile("/data3/smcosta/presel/forSim/nmul/70.0_100.0_ntphi_w_sWeights_nmult.root");
			}
			
		} 

		


	}
	else {
		FileName = Form("/data3/tasheng/presel/%sMC_nom.root",var_n.Data());	
	}

	fin = new TFile(FileName.Data());

	

	if (meson_n==0){

		EffInfoTree = (TTree * ) fin->Get("ntKp");

		if (whichvar == 0 && BsBP == 0){
			EffInfoTree_0 = (TTree * ) fin_0->Get("ntKp");
			EffInfoTree_1 = (TTree * ) fin_1->Get("ntKp");
			EffInfoTree_2 = (TTree * ) fin_2->Get("ntKp");
			EffInfoTree_3 = (TTree * ) fin_3->Get("ntKp");
			EffInfoTree_4 = (TTree * ) fin_4->Get("ntKp");
			EffInfoTree_5 = (TTree * ) fin_5->Get("ntKp");
			EffInfoTree_6 = (TTree * ) fin_6->Get("ntKp");
		} 
		else if (whichvar == 0 && BsBP == 1){
			EffInfoTree_0 = (TTree * ) fin_0->Get("ntKp");
			EffInfoTree_1 = (TTree * ) fin_1->Get("ntKp");
			EffInfoTree_2 = (TTree * ) fin_2->Get("ntKp");
			EffInfoTree_3 = (TTree * ) fin_3->Get("ntKp");
			
		}
		else if (whichvar == 1){
			EffInfoTree_0 = (TTree * ) fin_0->Get("ntKp");
			EffInfoTree_1 = (TTree * ) fin_1->Get("ntKp");
			EffInfoTree_2 = (TTree * ) fin_2->Get("ntKp");
			EffInfoTree_3 = (TTree * ) fin_3->Get("ntKp");
			EffInfoTree_4 = (TTree * ) fin_4->Get("ntKp");
		}
		else if (whichvar == 2){
			EffInfoTree_0 = (TTree * ) fin_0->Get("ntKp");
			EffInfoTree_1 = (TTree * ) fin_1->Get("ntKp");
			EffInfoTree_2 = (TTree * ) fin_2->Get("ntKp");
			EffInfoTree_3 = (TTree * ) fin_3->Get("ntKp");
			EffInfoTree_4 = (TTree * ) fin_4->Get("ntKp");
			EffInfoTree_5 = (TTree * ) fin_5->Get("ntKp");
			EffInfoTree_6 = (TTree * ) fin_6->Get("ntKp");
			EffInfoTree_7 = (TTree * ) fin_7->Get("ntKp");
		}
	}
	else {
		EffInfoTree = (TTree * ) fin->Get("ntphi");
		if (whichvar == 0){
			EffInfoTree_0 = (TTree * ) fin_0->Get("ntphi");
			EffInfoTree_1 = (TTree * ) fin_1->Get("ntphi");
			EffInfoTree_2 = (TTree * ) fin_2->Get("ntphi");
			EffInfoTree_3 = (TTree * ) fin_3->Get("ntphi");
		}
		else if (whichvar == 1){
			EffInfoTree_0 = (TTree * ) fin_0->Get("ntphi");
			EffInfoTree_1 = (TTree * ) fin_1->Get("ntphi");
			EffInfoTree_2 = (TTree * ) fin_2->Get("ntphi");
			EffInfoTree_3 = (TTree * ) fin_3->Get("ntphi");
			EffInfoTree_4 = (TTree * ) fin_4->Get("ntphi");
		}
		else if (whichvar == 2){
			EffInfoTree_0 = (TTree * ) fin_0->Get("ntphi");
			EffInfoTree_1 = (TTree * ) fin_1->Get("ntphi");
			EffInfoTree_2 = (TTree * ) fin_2->Get("ntphi");
			EffInfoTree_3 = (TTree * ) fin_3->Get("ntphi");
			EffInfoTree_4 = (TTree * ) fin_4->Get("ntphi");
			EffInfoTree_5 = (TTree * ) fin_5->Get("ntphi");
			EffInfoTree_6 = (TTree * ) fin_6->Get("ntphi");
			EffInfoTree_7 = (TTree * ) fin_7->Get("ntphi");
		}
	}

	fin->cd();

	int NEvents = EffInfoTree->GetEntries();
	
	
	Double_t BmassNew[NCand];
	Double_t BptNew[NCand];
	Double_t ByNew[NCand];
	Float_t BEff[NCand];
	Float_t BEffErr[NCand];

	Double_t sweight_0[NCand];
	Double_t sweight_1[NCand];
	Double_t sweight_2[NCand];
	Double_t sweight_3[NCand];
	Double_t sweight_4[NCand];
	Double_t sweight_5[NCand];
	Double_t sweight_6[NCand];
	Double_t sweight_7[NCand];
	

	Long64_t nMult;
	Long64_t trackSelection;

	Float_t BEffInv[NCand];
	Float_t BEffInvErr[NCand];
	Float_t BEffInv1D[NCand];
	Float_t BEffInvErr1D[NCand];

	Float_t BEffInvFit[NCand];
	Float_t BEffInvErrFit[NCand];

	Float_t BEffInvBDTWeighted[NCand];
	Float_t BEffInvErrBDTWeighted[NCand];



	Float_t BEffInvUp[NCand];
	Float_t BEffInvErrUp[NCand];
	Float_t BEffInvDown[NCand];
	Float_t BEffInvErrDown[NCand];

	
	EffInfoTree->SetBranchAddress("Bmass",BmassNew);
	EffInfoTree->SetBranchAddress("By",ByNew);
	EffInfoTree->SetBranchAddress("Bpt",BptNew);
	EffInfoTree->SetBranchAddress("nMult",&nMult); 
	EffInfoTree->SetBranchAddress("track", &trackSelection);
	EffInfoTree_0->SetBranchAddress("sWeight", &sweight_0);
	EffInfoTree_1->SetBranchAddress("sWeight", &sweight_1);
	EffInfoTree_2->SetBranchAddress("sWeight", &sweight_2);
	EffInfoTree_3->SetBranchAddress("sWeight", &sweight_3);
		if (meson_n==0 && whichvar == 0 && BsBP == 0){
			EffInfoTree_4->SetBranchAddress("sWeight", &sweight_4);
			EffInfoTree_5->SetBranchAddress("sWeight", &sweight_5);
			EffInfoTree_6->SetBranchAddress("sWeight", &sweight_6);
		}
		else if (whichvar == 1) {
			EffInfoTree_4->SetBranchAddress("sWeight", &sweight_4);
		}
		else if (whichvar == 2) {
			EffInfoTree_4->SetBranchAddress("sWeight", &sweight_4);
			EffInfoTree_5->SetBranchAddress("sWeight", &sweight_5);
			EffInfoTree_6->SetBranchAddress("sWeight", &sweight_6);
			EffInfoTree_7->SetBranchAddress("sWeight", &sweight_7);
		}

	Float_t var;
	Float_t sw;
	TString var_M;
	TString var_file;

	if (whichvar==1){
		var_M="y";
		var_file="y";
		NBins= nyBins_both;} 
	if (whichvar==2){
		var_M="Mult";
		var_file="Mult";
		NBins=nmBins_both;}
	if (whichvar==0 && meson_n==0 && BsBP==0){
		var_M="p_{T}";
		var_file="pt";
		NBins=nptBinsBP;}
	if ((whichvar==0 && meson_n==1)|| (whichvar==0 && meson_n==0 && BsBP==1) ){
		var_M="p_{T}";
		var_file="pt";
		NBins=nptBins;}

	double ptBins[NBins + 1];
	double Ntrk[NBins];
	double Counts[NBins];
	int ck[NBins];
	double CountsTight[NBins];
	double CountsLoose[NBins];
	double SumCounts[NBins];
	double SumCountsErr[NBins];
	double SumCountsTight[NBins];
	double SumCountsLoose[NBins];
	double NewEff[NBins];
	double NewEffErr[NBins];
	double NewEffTight[NBins];
	double NewEffLoose[NBins];
	double NewEffReal[NBins];
	double NewEffRealErr[NBins];
	double SumCountsUp[NBins];
	double SumCountsErrUp[NBins];
	double SumCountsDown[NBins];
	double SumCountsErrDown[NBins];
	double SumCountsEff[NBins];
	double SumCountsEffErr[NBins];
	
	
	
	
	double lumi = 302.3;
	if (whichvar==1){
	for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ybinsvec[i];             //taken from parameter.h
		}
	} 
	if (whichvar==2){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  nmbinsvec[i];             //taken from parameter.h
		}
	}
	if (whichvar==0 && meson_n==0 && BsBP==0){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvecBP[i];             //taken from parameter.h
		}
	}
	if ((whichvar==0 && meson_n==1)|| (whichvar==0 && meson_n==0 && BsBP==1)){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvec[i];             //taken from parameter.h
		}
	}

	for(int i = 0; i < NBins; i++){
		Counts[i] = 0;
		Ntrk[i] = 0;
		ck[i]=0;
		CountsTight[i] = 0;
		CountsLoose[i] = 0;
		SumCounts[i] = 0;
		SumCountsErr[i] = 0;
		SumCountsEff[i] = 0;
		SumCountsEffErr[i] = 0;
		SumCountsUp[i] = 0;
		SumCountsErrUp[i] = 0;
		SumCountsDown[i] = 0;
		SumCountsErrDown[i] = 0;
		SumCountsTight[i] = 0;
		SumCountsLoose[i] = 0;
	}


	int EtaBin;
	int PtBin;

	double trgtnp1;
	double trktnp1;
	double muidtnp1;

	double trgtnp1systup;
	double trgtnp1systdown;
	double trgtnp1statup;
	double trgtnp1statdown;

	double trktnp1systup;
	double trktnp1systdown;
	double trktnp1statup;
	double trktnp1statdown;

	double muidtnp1systup;
	double muidtnp1systdown;
	double muidtnp1statup;
	double muidtnp1statdown;
	double tnptotal1;
	double tnptotal1up;
	double tnptotal1down;
	double tnptotal1systup;
	double tnptotal1systdown;
	double tnptotal1statup;
	double tnptotal1statdown;
	double trgtnp2;
	double trktnp2;
	double muidtnp2;
	double trgtnp2systup;
	double trgtnp2systdown;
	double trgtnp2statup;
	double trgtnp2statdown;
	double trktnp2systup;
	double trktnp2systdown;
	double trktnp2statup;
	double trktnp2statdown;
	double muidtnp2systup;
	double muidtnp2systdown;
	double muidtnp2statup;
	double muidtnp2statdown;
	double tnptotal2;
	double tnptotal2up;
	double tnptotal2down;
	double tnptotal2systup;
	double tnptotal2systdown;
	double tnptotal2statup;
	double tnptotal2statdown;

	double tnpabssystup;
	double tnpabssystdown;


	//pt Binned Correction//
	TFile * fin1DEff; 

	if(DoTnP == 0) fin1DEff = new TFile(Form("%s/NewEff2DMaps/EffFineNoTnP.root",var_n.Data()));
	if(DoTnP == 1) fin1DEff = new TFile(Form("%s/NewEff2DMaps/EffFineBDT.root",var_n.Data()));	


	fin1DEff->cd();


	//Add 2D eff calculations
	
	TH2D * invEff2D;

	invEff2D= (TH2D *) fin1DEff->Get("invEff2DY");
	//invEff2D= (TH2D *) fin1DEff->Get("invEff2DBptSyst");
	TH2D * invEffTrkTight = (TH2D *) fin1DEff->Get("invEffTrkTight");
	TH2D * invEffTrkLoose = (TH2D *) fin1DEff->Get("invEffTrkLoose");
	TH2D * DrawinvEff2D= (TH2D *) fin1DEff->Get("invEff2D");
	TH2D * DrawinvEff2DY= (TH2D *) fin1DEff->Get("invEff2DY");

	
	/*
	if (whichvar==3){invEff2D= (TH2D *) fin1DEff->Get("invEff2D");}
	else {invEff2D= (TH2D *) fin1DEff->Get("invEff2DY");}
	*/

	Float_t EffInvTrkTight;
    Float_t EffInvTrkLoose;

	int XBin;
	int YBin;
	double ymax=1.5;
	double pthigh ;
	double ptlow;
	double BMass;

	if (meson_n==0){
		BMass=5.27932;
		ptlow=5;
		pthigh = 60;}
	else {
		BMass=5.3663;
		ptlow=7;
		pthigh = 50;}
	
	if(BsBP==1){ // ENSURE SAME PT REGION for both mesons
		ptlow=7;
		pthigh=50;}

	for( int i = 0; i < NEvents; i++){

		EffInfoTree->GetEntry(i);
		int j=0;
		if (whichvar==1){var=abs(ByNew[j]);}
		if (whichvar==2){var=nMult;}
		if (whichvar==0){var=BptNew[j];}
		
		for(int k = 0; k < NBins; k++){

			if (k==0){EffInfoTree_0->GetEntry(ck[k]); sw=sweight_0[j];}
			if (k==1){EffInfoTree_1->GetEntry(ck[k]); sw=sweight_1[j];}
			if (k==2){EffInfoTree_2->GetEntry(ck[k]); sw=sweight_2[j];}
			if (k==3){EffInfoTree_3->GetEntry(ck[k]); sw=sweight_3[j];}
			if (k==4){EffInfoTree_4->GetEntry(ck[k]); sw=sweight_4[j];}
			if (k==5){EffInfoTree_5->GetEntry(ck[k]); sw=sweight_5[j];}
			if (k==6){EffInfoTree_6->GetEntry(ck[k]); sw=sweight_6[j];}
			if (k==7){EffInfoTree_7->GetEntry(ck[k]); sw=sweight_7[j];}
			//sw = 1 ;
			if(BmassNew[j]>5 && BmassNew[j]<6 && var > ptBins[k] && var < ptBins[k+1] && TMath::Abs(ByNew[j]) < 2.4  && ((BptNew[j] > ptlow && BptNew[j] < 10 && abs(ByNew[j]) > ymax )||(BptNew[j] > 10 && BptNew[j]<pthigh))){
				
				ck[k]=ck[k]+1;
				
				XBin = invEff2D->GetXaxis()->FindBin( BptNew[j]);
				YBin = invEff2D->GetYaxis()->FindBin( TMath::Abs(ByNew[j]));
			
				BEffInv[j] = invEff2D->GetBinContent(XBin,YBin);
				BEffInvErr[j] = invEff2D->GetBinError(XBin,YBin);

				if (EffInvTrkLoose > 0) {
					SumCountsLoose[k] += EffInvTrkLoose*sw;
					CountsLoose[k]+=sw;
				}
				
				if (trackSelection > 1 && EffInvTrkTight > 0) {
					SumCountsTight[k] += EffInvTrkTight +sw;
					CountsTight[k]+=sw;
				} 
				
				if(trackSelection>0 && BEffInv[j] > 0){
					SumCounts[k] = SumCounts[k] + BEffInv[j]*sw;
					SumCountsErr[k] = SumCountsErr[k] + BEffInvErr[j] * BEffInvErr[j]*sw*sw;
					Counts[k] = Counts[k] + sw;
					Ntrk[k] += nMult;
			}
			}
		}
	}


	TH1D * hInvEff = new TH1D("hInvEff","",NBins,ptBins);
	hInvEff->GetXaxis()->SetTitle(Form("%s %s [GeV/c]",var_N.Data(),var_file.Data()));
	hInvEff->GetYaxis()->SetTitle("<1/ #alpha #times #epsilon >");
	hInvEff->GetYaxis()->SetTitleOffset(1.4);
	hInvEff->GetXaxis()->CenterTitle();
	hInvEff->GetYaxis()->CenterTitle();
	hInvEff->SetMarkerColor(1);
	hInvEff->SetLineColor(1);
	hInvEff->SetMarkerStyle(20);

	TH1D * hEffInv = new TH1D("hEffInv","",NBins,ptBins);
	hEffInv->GetXaxis()->SetTitle(Form("%s %s [GeV/c]",var_N.Data(),var_M.Data()));
	hEffInv->GetYaxis()->SetTitle("<1/ #alpha #times #epsilon >^{-1}");
	hEffInv->GetYaxis()->SetTitleOffset(1.4);
	hEffInv->GetXaxis()->CenterTitle();
	hEffInv->GetYaxis()->CenterTitle();
	hEffInv->SetMarkerColor(1);
	hEffInv->SetLineColor(1);
	hEffInv->SetMarkerStyle(20);

	TH1D *hInvEffTight = (TH1D*) hInvEff->Clone("hInvEffTight");
	hInvEffTight->GetYaxis()->SetTitle("<1/(Eff * Acc)> - Tight track selection");

 	TH1D *hInvEffLoose = (TH1D*) hInvEff->Clone("hInvEffLoose");
	hInvEffLoose->GetYaxis()->SetTitle("<1/(Eff * Acc)> - Loose track selection");
	
	TH1D * hNtrk = new TH1D("hNtrk","",NBins,ptBins);
	hNtrk->GetXaxis()->SetTitle(Form("%s %s",var_N.Data(),var_M.Data()));
	hNtrk->GetYaxis()->SetTitle("N_{trk} count >");
	hNtrk->GetYaxis()->SetTitleOffset(1.4);
	hNtrk->GetXaxis()->CenterTitle();
	hNtrk->GetYaxis()->CenterTitle();
	hNtrk->SetMarkerColor(1);
	hNtrk->SetLineColor(1);
	hNtrk->SetMarkerStyle(20);

	TH1D * hcounts = new TH1D("hcounts","",NBins,ptBins);
	hcounts->GetXaxis()->SetTitle(Form("%s %s",var_N.Data(),var_M.Data()));
	hcounts->GetYaxis()->SetTitle("event count >");
	hcounts->GetYaxis()->SetTitleOffset(1.4);
	hcounts->GetXaxis()->CenterTitle();
	hcounts->GetYaxis()->CenterTitle();
	hcounts->SetMarkerColor(1);
	hcounts->SetLineColor(1);
	hcounts->SetMarkerStyle(20);


	

	for(int i = 0; i < NBins; i++){

		
		NewEff[i] = SumCounts[i]/Counts[i];
		NewEffErr[i] = TMath::Sqrt(SumCountsErr[i])/Counts[i];

		NewEffTight[i] = SumCountsTight[i]/CountsTight[i];
		NewEffLoose[i] = SumCountsLoose[i]/CountsLoose[i];

		hInvEff->SetBinContent(i+1,NewEff[i]);
		hInvEff->SetBinError(i+1,NewEffErr[i]);

		hEffInv->SetBinContent(i+1,1/NewEff[i]);
		hEffInv->SetBinError(i+1,NewEffErr[i]/(NewEff[i] * NewEff[i]));

		hInvEffTight->SetBinContent(i+1, NewEffTight[i]);
		hInvEffTight->SetBinError(i+1, epsilon);
    	// TODO: currently loose is the same as nominal
		hInvEffLoose->SetBinContent(i+1, NewEff[i]);
		hInvEffLoose->SetBinError(i+1, epsilon);
		
		hNtrk->SetBinContent(i+1,Ntrk[i]);
		hcounts->SetBinContent(i+1,Counts[i]);

		//	cout << "Real eff = " << SumCountsReal[i]/Counts[i] << endl;
		//cout << "Counts = " << Counts[i] << endl;

		cout << "--------------------------------------------------------------------------------------------------------------" << endl;

		cout << "   NewEff = " << NewEff[i] << "     NewEffErr = " << NewEffErr[i] << "  Fractional = " << NewEffErr[i]/NewEff[i] << endl;
		
		cout << "--------------------------------------------------------------------------------------------------------------" << endl;

		//	cout << "   NewEff = " << NewEffUp[i] << "     NewEffErr = " << NewEffErrUp[i] << "  Fractional = " << NewEffErrUp[i]/NewEffUp[i] << endl;

		//NewEffErr[i] = 0; //Remove Error on Efficiency Correction//
	}


	TFile * RawYield = new TFile(Form("../henri2022/ROOTfiles/yields_%s_binned_%s%s.root",var_n.Data(),var_file.Data(), bsbpbins.Data()));
	RawYield->cd();
	TH1D * hPt = (TH1D *) RawYield->Get("hPt");

	TFile * RawYieldTight;
	TH1D * hPtTight;
	RawYieldTight = new TFile(Form("../henri2022/ROOTfiles/yields_%s_binned_%s_trk%s.root",var_n.Data(),var_file.Data(), bsbpbins.Data()));
	hPtTight = (TH1D *) RawYieldTight->Get("hPt");

	double RawCount;
	double RawCountErr;
	double RawCountTight;
	double CorrYield= 0;
	double CorrYieldErr= 0;

	double CorrYieldDiff[NBins];
	double CorrYieldDiffErr[NBins];

	for(int i = 0; i < NBins;i++){

		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);

		cout << "RawCount = " << RawCount << "  RawCountErr = " << RawCountErr << " NewEff[i] =   " << NewEff[i] << "  NewEffErr[i] =  " << NewEffErr[i] << endl; 

		cout << "CORR YIELD PT:  " <<  RawCount *   NewEff[i]/(BRchain*2 * lumi) << endl;
		CorrYield = RawCount * (ptBins[i+1] - ptBins[i]) *  NewEff[i]  + CorrYield;
		CorrYieldErr = ((RawCountErr * (ptBins[i+1] - ptBins[i]) *  NewEff[i]) *(RawCountErr * (ptBins[i+1] - ptBins[i]) *  NewEff[i]) + (RawCount * (ptBins[i+1] - ptBins[i]) *  NewEffErr[i]) * (RawCount * (ptBins[i+1] - ptBins[i]) *  NewEffErr[i]))  + CorrYieldErr;

		cout << "PrintYield = " << RawCount* (ptBins[i+1] - ptBins[i]) *  NewEff[i] << endl;

	}

	
	TH1D * CorrDiffHis = new TH1D("hPtSigma","",NBins,ptBins);
	CorrDiffHis->GetXaxis()->SetTitle(Form("%s [GeV/c]",var_file.Data()));
	CorrDiffHis->GetYaxis()->SetTitle(Form("d #sigma/d %s (pb GeV^{-1} c)",var_M.Data()));

	CorrDiffHis->GetYaxis()->SetTitleOffset(1.3);
	CorrDiffHis->GetXaxis()->SetLabelSize(0.02);
	CorrDiffHis->GetYaxis()->SetLabelSize(0.02);
	CorrDiffHis->GetXaxis()->CenterTitle();
	CorrDiffHis->GetYaxis()->CenterTitle();

	TH1D * CorrDiffHisBin = new TH1D("CorrDiffHisBin","",NBins,ptBins);
	CorrDiffHisBin->GetXaxis()->SetTitle("nMult");
	CorrDiffHisBin->GetYaxis()->SetTitle("#sigma (pb)");

	CorrDiffHisBin->GetYaxis()->SetTitleOffset(1.3);
	CorrDiffHisBin->GetXaxis()->CenterTitle();
	CorrDiffHisBin->GetYaxis()->CenterTitle();

	CorrDiffHisBin->SetMarkerColor(kRed);
	CorrDiffHisBin->SetLineColor(kRed);	
	CorrDiffHisBin->SetMarkerSize(1);
	CorrDiffHisBin->SetMarkerStyle(20);

	TH1D * CorrDiffHisBin_tight = new TH1D("CorrDiffHisBin_tight","",NBins,ptBins);
	CorrDiffHisBin_tight->GetXaxis()->SetTitle("nMult");
	CorrDiffHisBin_tight->GetYaxis()->SetTitle("#sigma (pb)");

	CorrDiffHisBin_tight->GetYaxis()->SetTitleOffset(1.3);
	CorrDiffHisBin_tight->GetXaxis()->CenterTitle();
	CorrDiffHisBin_tight->GetYaxis()->CenterTitle();

	CorrDiffHisBin_tight->SetMarkerColor(kRed);
	CorrDiffHisBin_tight->SetLineColor(kRed);	
	CorrDiffHisBin_tight->SetMarkerSize(1);
	CorrDiffHisBin_tight->SetMarkerStyle(20);

	TH1D * InvEff1DHisTight = new TH1D("InvEff1DHisTight","",NBins,ptBins);
	InvEff1DHisTight->GetXaxis()->SetTitle("nMult");
	InvEff1DHisTight->GetYaxis()->SetTitle("#sigma (pb)");

	InvEff1DHisTight->GetYaxis()->SetTitleOffset(1.3);
	InvEff1DHisTight->GetXaxis()->CenterTitle();
	InvEff1DHisTight->GetYaxis()->CenterTitle();

	InvEff1DHisTight->SetMarkerColor(kRed);
	InvEff1DHisTight->SetLineColor(kRed);	
	InvEff1DHisTight->SetMarkerSize(1);
	InvEff1DHisTight->SetMarkerStyle(20);

	TH1D * InvEff1DHisLoose = new TH1D("InvEff1DHisLoose","",NBins,ptBins);
	InvEff1DHisLoose->GetXaxis()->SetTitle("nMult");
	InvEff1DHisLoose->GetYaxis()->SetTitle("#sigma (pb)");

	InvEff1DHisLoose->GetYaxis()->SetTitleOffset(1.3);
	InvEff1DHisLoose->GetXaxis()->CenterTitle();
	InvEff1DHisLoose->GetYaxis()->CenterTitle();

	InvEff1DHisLoose->SetMarkerColor(kRed);
	InvEff1DHisLoose->SetLineColor(kRed);	
	InvEff1DHisLoose->SetMarkerSize(1);
	InvEff1DHisLoose->SetMarkerStyle(20);

	TH1D * CorrDiffHisLow = new TH1D("hPtSigmaLow","",NBins,ptBins);
	CorrDiffHisLow->GetXaxis()->SetTitle(Form("%s [GeV/c]",var_file.Data()));
	CorrDiffHisLow->GetYaxis()->SetTitle(Form("d #sigma/d %s (pb GeV^{-1} c)",var_M.Data()));

	CorrDiffHisLow->GetYaxis()->SetTitleOffset(1.3);
	CorrDiffHisLow->GetXaxis()->SetLabelSize(0.02);
	CorrDiffHisLow->GetYaxis()->SetLabelSize(0.02);
	CorrDiffHisLow->GetXaxis()->CenterTitle();
	CorrDiffHisLow->GetYaxis()->CenterTitle();

	TH1D * CorrDiffHisTight = (TH1D*) CorrDiffHis->Clone("hPtSigma_tight");

	for(int i = 0; i < NBins;i++){
		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);
		CorrYieldDiff[i] = (RawCount *  NewEff[i])/(BRchain*2* lumi);
		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr *  NewEff[i]) *(RawCountErr  *  NewEff[i]) + (RawCount *  NewEffErr[i]) * (RawCount  *  NewEffErr[i]))/(BRchain*2* lumi);
		CorrDiffHis->SetBinContent(i+1,CorrYieldDiff[i]);
		CorrDiffHis->SetBinError(i+1,CorrYieldDiffErr[i]);

		RawCountTight = hPtTight->GetBinContent(i+1);
		CorrDiffHisTight->SetBinContent(i+1, (RawCountTight *  NewEffTight[i]) / (BRchain*2* lumi) );
		CorrDiffHisTight->SetBinError(i+1, epsilon);
	}

	//CorrDiffHis->SetTitle("(Preliminary) %s #rightarrow J/#psi K^{+} p_{T} Differential Cross Section in pp");

	int startLow;
	int endLow;

	if (meson_n==0 && whichvar==0 && BsBP==0){startLow=0; endLow=5;}
	if ((meson_n==1 && whichvar==0)|| (meson_n==0 && whichvar==0 && BsBP==1)){startLow=0; endLow=3;}
	if (whichvar==1){startLow=1; endLow=1;}

	for (int i=startLow; i<NBins-endLow;i++){
		CorrDiffHisLow->SetBinContent(i+1,CorrDiffHis->GetBinContent(i+1));
		CorrDiffHisLow->SetBinError(i+1,CorrDiffHis->GetBinError(i+1));
	}


	CorrDiffHis->SetMarkerColor(kBlack);
	CorrDiffHis->SetMarkerSize(1);
	CorrDiffHis->SetMarkerStyle(20);

	CorrDiffHisLow->SetMarkerColor(kBlue);
	CorrDiffHisLow->SetMarkerSize(1);
	CorrDiffHisLow->SetMarkerStyle(20);

	
	



	TFile * foutCorr;
	if(DoTnP == 0 && usemc==0)	foutCorr = new TFile(Form("%s/FinalFiles/%sPPCorrYield%sNoTnP%s.root",var_n2.Data(),var_n.Data(),var_file.Data(),bsbpbins.Data()),"RECREATE");
	if(DoTnP == 1 && usemc==0)	foutCorr = new  TFile(Form("%s/FinalFiles/%sPPCorrYield%s%s.root",var_n2.Data(),var_n.Data(),var_file.Data(),bsbpbins.Data()),"RECREATE");
	if(DoTnP == 0 && usemc==1)	foutCorr = new TFile(Form("%s/FinalFiles/%sPPCorrYield%sNoTnPMC%s.root",var_n2.Data(),var_n.Data(),var_file.Data(),bsbpbins.Data()),"RECREATE");
	if(DoTnP == 1 && usemc==1)	foutCorr = new  TFile(Form("%s/FinalFiles/%sPPCorrYield%sMC%s.root",var_n2.Data(),var_n.Data(),var_file.Data(),bsbpbins.Data()),"RECREATE");

	TH1D * Eff1DHisvar;
	TH1D * InvEff1DHisvarTight;
	TH1D * InvEff1DHisvarLoose;
	if (whichvar==1){
		Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisY");
		//Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisYFid");
		//Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisYFid10");

		InvEff1DHisvarTight=(TH1D *) fin1DEff->Get("InvEff1DHisTightY");
		InvEff1DHisvarLoose=(TH1D *) fin1DEff->Get("InvEff1DHisLooseY");
	    }
	if (whichvar==2){
		Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisMult");
		//Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisMultFid");
		//Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisMultFid10");

		InvEff1DHisvarTight=(TH1D *) fin1DEff->Get("InvEff1DHisTightMult");
		InvEff1DHisvarLoose=(TH1D *) fin1DEff->Get("InvEff1DHisLooseMult");
		}
	if (whichvar==0){
		Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHis");
		//Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisFid");
		//Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisFid10");

		InvEff1DHisvarTight=(TH1D *) fin1DEff->Get("InvEff1DHisTight");
		InvEff1DHisvarLoose=(TH1D *) fin1DEff->Get("InvEff1DHisLoose");
	}

	float Eff1D[NBins];
	float Eff1DErr[NBins];
	float Eff1DTight[NBins];
	float Eff1DTightErr[NBins];
	double XTemp;
	double YTemp;

	for(int i = 0; i < NBins;i++){
		
		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);
		Eff1D[i] = Eff1DHisvar->GetBinContent(i+1);
		Eff1DErr[i] = Eff1DHisvar->GetBinError(i+1);

		Eff1DTight[i] = InvEff1DHisvarTight->GetBinContent(i+1);
		Eff1DTightErr[i] = InvEff1DHisvarTight->GetBinError(i+1);

		CorrYieldDiff[i] = (RawCount /  Eff1D[i])/(BRchain*2* lumi);
		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr /  Eff1D[i]) *(RawCountErr  /  Eff1D[i]) + (RawCount /Eff1D[i] *  Eff1DErr[i]) * (RawCount /Eff1D[i] *  Eff1DErr[i]))/(BRchain*2* lumi);

		CorrDiffHisBin->SetBinContent(i+1,CorrYieldDiff[i]);
		CorrDiffHisBin->SetBinError(i+1,CorrYieldDiffErr[i]);

		CorrYieldDiff[i] = (RawCount *  Eff1DTight[i])/(BRchain*2* lumi);
		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr *  Eff1DTight[i]) *(RawCountErr  *  Eff1DTight[i]) + (RawCount *  Eff1DTightErr[i]) * (RawCount  *  Eff1DTightErr[i]))/(BRchain*2* lumi);

		CorrDiffHisBin_tight->SetBinContent(i+1,CorrYieldDiff[i]);
		CorrDiffHisBin_tight->SetBinError(i+1,CorrYieldDiffErr[i]);

		InvEff1DHisLoose->SetBinContent(i+1,InvEff1DHisvarLoose->GetBinContent(i+1));
		InvEff1DHisLoose->SetBinError(i+1,InvEff1DHisvarLoose->GetBinError(i+1));

		InvEff1DHisTight->SetBinContent(i+1,InvEff1DHisvarTight->GetBinContent(i+1));
		InvEff1DHisTight->SetBinError(i+1,InvEff1DHisvarTight->GetBinError(i+1));
	}

	double dmax = -10000;
	for(int i = 0; i < NBins;i++){
		if (dmax < hInvEff->GetBinContent(i+1)){dmax=hInvEff->GetBinContent(i+1);}
	}

	hInvEff->SetMaximum(dmax*1.5);
	TCanvas *c = new TCanvas("c","c",700,700);
	c->cd();

	hInvEff->Draw("ep");
	
	if (usemc==0) {c->SaveAs(Form("%s/EffFinal/ReAnaEff_%dBins_%s.pdf",var_n2.Data(),NBins,var_file.Data()));}
	else {c->SaveAs(Form("%s/EffFinal/ReAnaEff_%dBins_%s_MC.pdf",var_n2.Data(),NBins,var_file.Data()));}
	
	Eff1DHisvar->Draw("ep");

	if (usemc==0) {c->SaveAs(Form("%s/EffFinal/ReAnaEff1D_%dBins_%s.pdf",var_n2.Data(),NBins,var_file.Data()));}
	else {c->SaveAs(Form("%s/EffFinal/ReAnaEff1D_%dBins_%s_MC.pdf",var_n2.Data(),NBins,var_file.Data()));}
	
	TH2D * invAcc2D=(TH2D *) fin1DEff->Get("invAcc2D");
	invAcc2D->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	invAcc2D->GetYaxis()->SetTitle("rapidity");

	invAcc2D->Draw("pcolz");
	c->SaveAs(Form("%s/EffFinal/acc_2Dmap%s.pdf",var_n2.Data(),bsbpbins.Data()));


	DrawinvEff2D->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	DrawinvEff2D->GetYaxis()->SetTitle("|y|");
	DrawinvEff2D->GetZaxis()->SetLabelSize(0.02);

	DrawinvEff2D->Draw("pcolz");
	c->SaveAs(Form("%s/EffFinal/totaleff_2Dmap_%s%s.pdf",var_n2.Data(),var_n2.Data(), bsbpbins.Data()));
	


	DrawinvEff2DY->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	DrawinvEff2DY->GetYaxis()->SetTitle("|y|");
	DrawinvEff2DY->GetZaxis()->SetLabelSize(0.02);
	c->SetLogz();
	DrawinvEff2DY->Draw("pcolz");

    c->SaveAs(Form("%s/EffFinal/totaleff_Fid_2Dmap_%s%s.pdf",var_n2.Data(),var_n2.Data(),bsbpbins.Data()));
	
	
	CorrDiffHis->SetMinimum(0);
	CorrDiffHisLow->SetMinimum(0);
	CorrDiffHis->Draw("ep");
	CorrDiffHisLow->Draw("epsame");

	TLegend* leg3 = new TLegend(0.55,0.65,0.85,0.83,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetFillStyle(0);
	leg3->AddEntry(CorrDiffHis,"Cross Section","PL");
	if(whichvar == 0) {leg3->AddEntry(CorrDiffHisLow,"Cross Section (|y| > 1.5)","PL");}
	if(whichvar == 1) {leg3->AddEntry(CorrDiffHisLow,"Cross Section (p_{T} > 10 GeV/c)","PL");}
	leg3->Draw("same");

	if (usemc==0){c->SaveAs(Form("%s/EffFinal/xsection_%s_%s.pdf",var_n2.Data(),var_file.Data(),var_n2.Data()));}
	else {c->SaveAs(        Form("%s/EffFinal/xsection_%s_%s_MC.pdf",var_n2.Data(),var_file.Data(),var_n2.Data()));}
	

	double cmax = -10000;
	double cmin =  10000;
	for(int i = 0; i < NBins;i++){
		if (cmax < hEffInv->GetBinContent(i+1)){cmax=hEffInv->GetBinContent(i+1);}
		if (cmax < Eff1DHisvar->GetBinContent(i+1)){cmax=Eff1DHisvar->GetBinContent(i+1);}
		if (cmin > hEffInv->GetBinContent(i+1)){cmin=hEffInv->GetBinContent(i+1);}
		if (cmin > Eff1DHisvar->GetBinContent(i+1)){cmin=Eff1DHisvar->GetBinContent(i+1);}
	}

	hEffInv->GetYaxis()->SetRangeUser( cmin*0.9, cmax*1.1);
	hEffInv->GetYaxis()->SetTitle("#alpha #times #epsilon");
	hEffInv->SetMarkerColor(kRed+1);
	hEffInv->SetLineColor(kRed+1);
	hEffInv->Draw("ep");
	if (whichvar!=2){
		Eff1DHisvar->Draw("epsame");
	}
	

	TLegend *leg = new TLegend(0.15,0.82,0.4,0.87,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.020);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(1);

	leg->AddEntry(hEffInv,"2D < 1/ #alpha #times #epsilon >^{-1}","pl");
	if (whichvar!=2){
		leg->AddEntry(Eff1DHisvar,"1D","pl");
	}
	leg->Draw("SAME");

		if (usemc==0){c->SaveAs(Form("%s/EffFinal/ReAnaEffcomp_%dBins_%s.pdf",var_n2.Data(),NBins,var_file.Data()));}
		else {c->SaveAs(Form("%s/EffFinal/ReAnaEffcomp_%dBins_%s_MC.pdf",var_n2.Data(),NBins,var_file.Data()));}
	



	double tusk;
	double d4c;

	std::vector<double> cv2D;
	std::vector<double> cvinv;

	
	std::vector<double> eff_val1d;
	std::vector<double> eff_val2d;
	std::vector<double> eff_err1d;
	std::vector<double> eff_err2d;
	
	std::vector<double> eff_inv_val1d;
	std::vector<double> eff_inv_val2dinv;
	std::vector<double> eff_inv_err1d;
	std::vector<double> eff_inv_err2dinv;
	
	
	std::vector<std::vector<double>> comparison_values;
	std::vector<std::vector<double>> comparison_inv_values;
	std::vector<std::vector<double>> eff_val;
	std::vector<std::vector<double>> eff_err;
	std::vector<std::vector<double>> eff_inv_val;
	std::vector<std::vector<double>> eff_inv_err;
	std::vector<std::string> labels = {"2D"};
	std::vector<std::string> labelseff = {"1D","2D"};
	std::vector<std::string> col_name;
	std::vector<std::string> col_name_eff;
	
	
	string name;
	TString whichvarname;
	col_name.push_back("Relative difference 1D vs:");
	col_name_eff.push_back("Eff values");
	if(whichvar==0){name="$<p_T<$"; whichvarname="pt";} else if(whichvar==1){name="$<y<$";whichvarname="y";} else if(whichvar==2){name="$<nTrks<$";whichvarname="nMult";}
	for(int i=0;i<NBins;i++){
		std::ostringstream clabel;
		clabel<<ptBins[i]<<name<<ptBins[i+1];
		std::string label1 = clabel.str();
		col_name.push_back(label1);
		col_name_eff.push_back(label1);

		cv2D.push_back(abs(1/NewEff[i]-Eff1D[i])/Eff1D[i]*100);
		cvinv.push_back(abs(NewEff[i]-1/Eff1D[i])*Eff1D[i]*100);

		eff_val1d.push_back(Eff1D[i]);
		eff_val2d.push_back(1/NewEff[i]);
		
		eff_err1d.push_back(Eff1DErr[i]);
		eff_err2d.push_back(NewEffErr[i]/(NewEff[i]*NewEff[i]));
		
		eff_inv_val1d.push_back(1/Eff1D[i]);
		eff_inv_val2dinv.push_back(NewEff[i]);

		eff_inv_err1d.push_back(Eff1DErr[i]/(Eff1D[i]*Eff1D[i]));
		eff_inv_err2dinv.push_back(NewEffErr[i]);
	}
	
	comparison_values.push_back(cv2D);

	comparison_inv_values.push_back(cvinv);

	eff_val.push_back(eff_val1d);
	eff_val.push_back(eff_val2d);
	
	eff_err.push_back(eff_err1d);
	eff_err.push_back(eff_err2d);
	
	eff_inv_val.push_back(eff_inv_val1d);
	eff_inv_val.push_back(eff_inv_val2dinv);
	
	eff_inv_err.push_back(eff_inv_err1d);
	eff_inv_err.push_back(eff_inv_err2dinv);
	

	latex_table(Form("1D2Dcomparisons_%s",whichvarname.Data()), NBins+1,  2,  col_name , labels , comparison_values, "1D vs 2D efficiency comparisons",0);
	latex_table(Form("1D2Dcomparisons_inv_%s",whichvarname.Data()), NBins+1,  2,  col_name , labels , comparison_inv_values, "1D vs 2D efficiency comparisons",0);
	latex_table(Form("Effyields_%s",whichvarname.Data()), NBins+1,  3,  col_name_eff , labelseff , eff_val , "Eff Values per bin",1, eff_err);
	latex_table(Form("Effyields_inv_%s",whichvarname.Data()), NBins+1,  3,  col_name_eff , labelseff , eff_inv_val , "Eff Values per bin",1, eff_inv_err);
	
	std::vector<std::string> filetype ={"_check.aux", "_check.log",".tex","_check.tex"};
	for (int j=0;j<(int)(filetype.size());j++){
				rename(("1D2Dcomparisons_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n2.Data())+"/Trash/1D2Dcomparisons_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				rename(("Effyields_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n2.Data())+"/Trash/Effyields_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				rename(("1D2Dcomparisons_inv_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n2.Data())+"/Trash/1D2Dcomparisons_inv_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				rename(("Effyields_inv_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n2.Data())+"/Trash/Effyields_inv_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				
			}
	rename(("1D2Dcomparisons_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n2.Data())+"/EffFinal/1D2Dcomparisons_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("Effyields_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n2.Data())+"/EffFinal/Effyields_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("1D2Dcomparisons_inv_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n2.Data())+"/EffFinal/1D2Dcomparisons_inv_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("Effyields_inv_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n2.Data())+"/EffFinal/Effyields_inv_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	



	foutCorr->cd();

	CorrDiffHis->Write();
	hNtrk->Write();
	hcounts->Write();
	CorrDiffHisBin->Write();
	CorrDiffHisBin_tight->Write();
	CorrDiffHisTight->Write();
	hPtTight->SetName("hPtTight");
	hPtTight->Write();
	hInvEffTight->Write();
	hInvEffLoose->Write();
	hInvEff->Write();
	Eff1DHisvar->Write();
	InvEff1DHisTight->Write();
	InvEff1DHisLoose->Write();
	hPt->Write();

	foutCorr->Close();

}