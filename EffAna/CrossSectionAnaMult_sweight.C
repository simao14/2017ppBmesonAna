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
	if (meson_n == 0){
		var_n="BP";
		var_n2="Bp";
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
	
	gSystem->mkdir( var_n.Data() , true);
	gSystem->mkdir(Form("%s/EffFinal",var_n.Data()) ,true);
	gSystem->mkdir(Form("%s/FinalFiles",var_n.Data()) ,true);
	gSystem->mkdir(Form("%s/Trash",var_n.Data()) ,true);

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
	TTree * EffInfoTree;
	TTree * EffInfoTree_0;
	TTree * EffInfoTree_1;
	TTree * EffInfoTree_2;
	TTree * EffInfoTree_3;
	TTree * EffInfoTree_4;
	TTree * EffInfoTree_5;
	TTree * EffInfoTree_6;
	TTree * root;

	if (usemc==0){
		if (meson_n==0){
			FileName ="/data3/smcosta/presel/ntKp_w_sWeights.root";
			fin_0 = new TFile("/data3/smcosta/presel/5_7_ntKp_w_sWeights.root");
			fin_1 = new TFile("/data3/smcosta/presel/7_10_ntKp_w_sWeights.root");
			fin_2 = new TFile("/data3/smcosta/presel/10_15_ntKp_w_sWeights.root");
			fin_3 = new TFile("/data3/smcosta/presel/15_20_ntKp_w_sWeights.root");
			fin_4 = new TFile("/data3/smcosta/presel/20_30_ntKp_w_sWeights.root");
			fin_5 = new TFile("/data3/smcosta/presel/30_50_ntKp_w_sWeights.root");
			fin_6 = new TFile("/data3/smcosta/presel/50_60_ntKp_w_sWeights.root");
		} else {
			FileName = "/data3/smcosta/presel/ntphi_w_sWeights.root";
			fin_0 = new TFile("/data3/smcosta/presel/7_10_ntphi_w_sWeights.root");
			fin_1 = new TFile("/data3/smcosta/presel/10_15_ntphi_w_sWeights.root");
			fin_2 = new TFile("/data3/smcosta/presel/15_20_ntphi_w_sWeights.root");
			fin_3 = new TFile("/data3/smcosta/presel/20_50_ntphi_w_sWeights.root");
		}
	}
	else {
		if (meson_n==0){
			FileName = Form("/data3/tasheng/presel/%sMC_nom.root",var_n.Data());
		} else {
			FileName = Form("/data3/tasheng/presel/%sMC_nom.root",var_n.Data());
		}
		
	}
	fin = new TFile(FileName.Data());

	

	if (meson_n==0){
		EffInfoTree = (TTree * ) fin->Get("ntKp");
		EffInfoTree_0 = (TTree * ) fin_0->Get("ntKp");
		EffInfoTree_1 = (TTree * ) fin_1->Get("ntKp");
		EffInfoTree_2 = (TTree * ) fin_2->Get("ntKp");
		EffInfoTree_3 = (TTree * ) fin_3->Get("ntKp");
		EffInfoTree_4 = (TTree * ) fin_4->Get("ntKp");
		EffInfoTree_5 = (TTree * ) fin_5->Get("ntKp");
		EffInfoTree_6 = (TTree * ) fin_6->Get("ntKp");
		}
	else {
		EffInfoTree = (TTree * ) fin->Get("ntphi");
		EffInfoTree_0 = (TTree * ) fin_0->Get("ntphi");
		EffInfoTree_1 = (TTree * ) fin_1->Get("ntphi");
		EffInfoTree_2 = (TTree * ) fin_2->Get("ntphi");
		EffInfoTree_3 = (TTree * ) fin_3->Get("ntphi");
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
	

	Int_t nMult;
	Double_t trackSelection;

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
	//EffInfoTree->SetBranchAddress("nMult",&nMult); 
	EffInfoTree->SetBranchAddress("track", &trackSelection);
	EffInfoTree_0->SetBranchAddress("sWeight", &sweight_0);
	EffInfoTree_1->SetBranchAddress("sWeight", &sweight_1);
	EffInfoTree_2->SetBranchAddress("sWeight", &sweight_2);
	EffInfoTree_3->SetBranchAddress("sWeight", &sweight_3);
	if (meson_n==0){
		EffInfoTree_4->SetBranchAddress("sWeight", &sweight_4);
		EffInfoTree_5->SetBranchAddress("sWeight", &sweight_5);
		EffInfoTree_6->SetBranchAddress("sWeight", &sweight_6);
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

	int Counts[NBins];
	int ck[NBins];
	double Counts_weighted[NBins];
	double Counts_weighted_mass[NBins];
	int CountsTight[NBins];
	int CountsLoose[NBins];
	double SumCounts[NBins];
	double SumCountsErr[NBins];
	double SumCounts_weighted[NBins];
	double SumCountsErr_weighted[NBins];
	double SumCounts_weighted_mass[NBins];
	double SumCountsErr_weighted_mass[NBins];
	double SumCountsTight[NBins];
	double SumCountsLoose[NBins];
	double NewEff[NBins];
	double NewEffErr[NBins];
	double NewEff_weighted[NBins];
	double NewEffErr_weighted[NBins];
	double NewEff_weighted_mass[NBins];
	double NewEffErr_weighted_mass[NBins];
	double NewEffTight[NBins];
	double NewEffLoose[NBins];
	double NewEffReal[NBins];
	double NewEffRealErr[NBins];
	double NewEffReal_weighted[NBins];
	double NewEffRealErr_weighted[NBins];
	double SumCountsUp[NBins];
	double SumCountsErrUp[NBins];
	double SumCountsDown[NBins];
	double SumCountsErrDown[NBins];
	double SumCountsEff[NBins];
	double SumCountsEffErr[NBins];
	double SumCountsEff_weighted[NBins];
	double SumCountsEffErr_weighted[NBins];
	
	
	
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
		ck[i]=0;
		Counts_weighted[i] = 0;
		Counts_weighted_mass[i] = 0;
		CountsTight[i] = 0;
		CountsLoose[i] = 0;
		SumCounts[i] = 0;
		SumCountsErr[i] = 0;
		SumCounts_weighted[i] = 0;
		SumCountsErr_weighted[i] = 0;
		SumCounts_weighted_mass[i] = 0;
		SumCountsErr_weighted_mass[i] = 0;
		SumCountsEff[i] = 0;
		SumCountsEffErr[i] = 0;
		SumCountsEff_weighted[i] = 0;
		SumCountsEffErr_weighted[i] = 0;
		SumCountsUp[i] = 0;
		SumCountsErrUp[i] = 0;
		SumCountsDown[i] = 0;
		SumCountsErrDown[i] = 0;
		SumCountsTight[i] = 0;
		SumCountsLoose[i] = 0;
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
	TH1D * invEff1D;

	invEff2D= (TH2D *) fin1DEff->Get("invEff2DY");
	//invEff2D= (TH2D *) fin1DEff->Get("invEff2DBptSyst");
	TH2D * invEffTrkTight = (TH2D *) fin1DEff->Get("invEffTrkTight");
	TH2D * invEffTrkLoose = (TH2D *) fin1DEff->Get("invEffTrkLoose");
	TH2D * DrawinvEff2D= (TH2D *) fin1DEff->Get("invEff2D");
	TH2D * DrawinvEff2DY= (TH2D *) fin1DEff->Get("invEff2DY");

	if (whichvar==1){ invEff1D = (TH1D *) fin1DEff->Get("invEff1DY");} 
	else if (whichvar==0){ invEff1D = (TH1D *) fin1DEff->Get("invEff1DFGpt");}
	/*
	if (whichvar==3){invEff2D= (TH2D *) fin1DEff->Get("invEff2D");}
	else {invEff2D= (TH2D *) fin1DEff->Get("invEff2DY");}
	*/

	Float_t EffInvTrkTight;
    Float_t EffInvTrkLoose;

	int XBin;
	int YBin;
	int DBin;
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
		if (whichvar==1){var=ByNew[j];}
		if (whichvar==2){var=nMult;}
		if (whichvar==0){var=BptNew[j];}
		
		for(int k = 0; k < NBins; k++){

			if (k==0){EffInfoTree_0->GetEntry(ck[k]);   sw=sweight_0[j];}
			if (k==1){EffInfoTree_1->GetEntry(ck[k]);  sw=sweight_1[j];}
			if (k==2){EffInfoTree_2->GetEntry(ck[k]); sw=sweight_2[j];}
			if (k==3){EffInfoTree_3->GetEntry(ck[k]); sw=sweight_3[j];}
			if (k==4){EffInfoTree_4->GetEntry(ck[k]); sw=sweight_4[j];}
			if (k==5){EffInfoTree_5->GetEntry(ck[k]); sw=sweight_5[j];}
			if (k==6){EffInfoTree_6->GetEntry(ck[k]); sw=sweight_6[j];}
			
			if(var > ptBins[k] && var < ptBins[k+1]){
				
				ck[k]=ck[k]+1;
				
				XBin = invEff2D->GetXaxis()->FindBin( BptNew[j]);
				YBin = invEff2D->GetYaxis()->FindBin( TMath::Abs(ByNew[j]));
			
				BEffInv[j] = invEff2D->GetBinContent(XBin,YBin);
				BEffInvErr[j] = invEff2D->GetBinError(XBin,YBin);

				DBin = invEff1D->FindBin( TMath::Abs(var));

				BEffInv1D[j] = invEff1D->GetBinContent(DBin);
				BEffInvErr1D[j] = invEff1D->GetBinError(DBin);

				
				SumCounts_weighted_mass[k] = SumCounts_weighted_mass[k] + BEffInv[j]*sw;
				SumCountsErr_weighted_mass[k] = SumCountsErr_weighted_mass[k] + BEffInvErr[j] * BEffInvErr[j]*sw*sw;
				Counts_weighted_mass[k] = Counts_weighted_mass[k] + sw;

				if( (TMath::Abs(BmassNew[j] - BMass) < 0.08) && TMath::Abs(ByNew[j]) < 2.4  && ((BptNew[j] > ptlow && BptNew[j] < 10 && abs(ByNew[j]) > ymax )||(BptNew[j] > 10 && BptNew[j]<pthigh)))
				{
					
					
					if(trackSelection>0 && BEffInv[j] > 0){
						SumCounts[k] = SumCounts[k] + BEffInv[j];
						SumCountsErr[k] = SumCountsErr[k] + BEffInvErr[j] * BEffInvErr[j];
						
						SumCountsEff[k] = SumCountsEff[k] + BEffInv1D[j];
						SumCountsEffErr[k] = SumCountsEffErr[k] + BEffInvErr1D[j] * BEffInvErr1D[j];

						SumCounts_weighted[k] = SumCounts_weighted[k] + BEffInv[j]*sw;
						SumCountsErr_weighted[k] = SumCountsErr_weighted[k] + BEffInvErr[j] * BEffInvErr[j]*sw*sw;
						
						SumCountsEff_weighted[k] = SumCountsEff_weighted[k] + BEffInv1D[j]*sw;
						SumCountsEffErr_weighted[k] = SumCountsEffErr_weighted[k] + BEffInvErr1D[j] * BEffInvErr1D[j]*sw*sw;

						Counts[k] = Counts[k] + 1;
						Counts_weighted[k] = Counts_weighted[k] + sw;
				}
			}
			}
		}
	}


	TH1D * hInvEff = new TH1D("hInvEff","",NBins,ptBins);
	hInvEff->GetXaxis()->SetTitle(Form("%s %s [GeV/c]",var_N.Data(),var_file.Data()));
	hInvEff->GetYaxis()->SetTitle("<1/ #alpha #times #epsilon >^{-1}");
	hInvEff->GetYaxis()->SetTitleOffset(1.4);
	hInvEff->GetXaxis()->CenterTitle();
	hInvEff->GetYaxis()->CenterTitle();
	hInvEff->SetMarkerColor(1);
	hInvEff->SetLineColor(1);
	hInvEff->SetMarkerStyle(20);
	//hInvEff->SetMinimum(0);

	TH1D * hInvEff_weighted = new TH1D("hInvEff_weighted","",NBins,ptBins);
	hInvEff_weighted->GetXaxis()->SetTitle(Form("%s %s [GeV/c]",var_N.Data(),var_file.Data()));
	hInvEff_weighted->GetYaxis()->SetTitle("<1/ #alpha #times #epsilon >^{-1}");
	hInvEff_weighted->GetYaxis()->SetTitleOffset(1.4);
	hInvEff_weighted->GetXaxis()->CenterTitle();
	hInvEff_weighted->GetYaxis()->CenterTitle();
	hInvEff_weighted->SetMarkerColor(1);
	hInvEff_weighted->SetLineColor(1);
	hInvEff_weighted->SetMarkerStyle(20);
	//hInvEff->SetMinimum(0);

	TH1D * hInvEff_weighted_mass = new TH1D("hInvEff_weighted_mass","",NBins,ptBins);
	hInvEff_weighted_mass->GetXaxis()->SetTitle(Form("%s %s [GeV/c]",var_N.Data(),var_file.Data()));
	hInvEff_weighted_mass->GetYaxis()->SetTitle("<1/ #alpha #times #epsilon >^{-1}");
	hInvEff_weighted_mass->GetYaxis()->SetTitleOffset(1.4);
	hInvEff_weighted_mass->GetXaxis()->CenterTitle();
	hInvEff_weighted_mass->GetYaxis()->CenterTitle();
	hInvEff_weighted_mass->SetMarkerColor(1);
	hInvEff_weighted_mass->SetLineColor(1);
	hInvEff_weighted_mass->SetMarkerStyle(20);
	//hInvEff->SetMinimum(0);



	TH1D * hEffInv1D = new TH1D("hEffInv1D","",NBins,ptBins);
	hEffInv1D->GetXaxis()->SetTitle(Form("%s %s [GeV/c]",var_N.Data(),var_file.Data()));
	hEffInv1D->GetYaxis()->SetTitle("<1/ #alpha #times #epsilon >^{-1}");
	hEffInv1D->GetYaxis()->SetTitleOffset(1.4);
	hEffInv1D->GetXaxis()->CenterTitle();
	hEffInv1D->GetYaxis()->CenterTitle();
	hEffInv1D->SetMarkerColor(1);
	hEffInv1D->SetLineColor(1);
	hEffInv1D->SetMarkerStyle(20);

	TH1D * hEffInv1D_weighted = new TH1D("hEffInv1D_weighted","",NBins,ptBins);
	hEffInv1D_weighted->GetXaxis()->SetTitle(Form("%s %s [GeV/c]",var_N.Data(),var_file.Data()));
	hEffInv1D_weighted->GetYaxis()->SetTitle("<1/ #alpha #times #epsilon >^{-1}");
	hEffInv1D_weighted->GetYaxis()->SetTitleOffset(1.4);
	hEffInv1D_weighted->GetXaxis()->CenterTitle();
	hEffInv1D_weighted->GetYaxis()->CenterTitle();
	hEffInv1D_weighted->SetMarkerColor(1);
	hEffInv1D_weighted->SetLineColor(1);
	hEffInv1D_weighted->SetMarkerStyle(20);


	

	

	for(int i = 0; i < NBins; i++){

		NewEff[i] = SumCounts[i]/Counts[i];
		NewEffErr[i] = TMath::Sqrt(SumCountsErr[i])/Counts[i];

		NewEffReal[i] = SumCountsEff[i]/Counts[i];
		NewEffRealErr[i] = TMath::Sqrt(SumCountsEffErr[i])/Counts[i];

		NewEff_weighted[i] = SumCounts_weighted[i]/Counts_weighted[i];
		NewEffErr_weighted[i] = TMath::Sqrt(SumCountsErr_weighted[i])/Counts_weighted[i];

		NewEff_weighted_mass[i] = SumCounts_weighted_mass[i]/Counts_weighted_mass[i];
		NewEffErr_weighted_mass[i] = TMath::Sqrt(SumCountsErr_weighted_mass[i])/Counts_weighted_mass[i];

		NewEffReal_weighted[i] = SumCountsEff_weighted[i]/Counts_weighted[i];
		NewEffRealErr_weighted[i] = TMath::Sqrt(SumCountsEffErr_weighted[i])/Counts_weighted[i];

		hInvEff->SetBinContent(i+1,1/NewEff[i]);
		hInvEff->SetBinError(i+1,NewEffErr[i]/(NewEff[i] * NewEff[i]));
		
		hEffInv1D->SetBinContent(i+1,1/NewEffReal[i]);
		hEffInv1D->SetBinError(i+1,NewEffRealErr[i]/(NewEffReal[i]*NewEffReal[i]));

		hInvEff_weighted->SetBinContent(i+1,1/NewEff_weighted[i]);
		hInvEff_weighted->SetBinError(i+1,NewEffErr_weighted[i]/(NewEff_weighted[i] * NewEff_weighted[i]));

		hInvEff_weighted_mass->SetBinContent(i+1,1/NewEff_weighted_mass[i]);
		hInvEff_weighted_mass->SetBinError(i+1,NewEffErr_weighted_mass[i]/(NewEff_weighted_mass[i] * NewEff_weighted_mass[i]));

		hEffInv1D_weighted->SetBinContent(i+1,1/NewEffReal_weighted[i]);
		hEffInv1D_weighted->SetBinError(i+1,NewEffRealErr_weighted[i]/(NewEffReal_weighted[i]*NewEffReal_weighted[i]));

		

		//	cout << "Real eff = " << SumCountsReal[i]/Counts[i] << endl;
		//cout << "Counts = " << Counts[i] << endl;

		cout << "--------------------------------------------------------------------------------------------------------------" << endl;

		cout << "   NewEff = " << NewEff[i] << "     NewEffErr = " << NewEffErr[i] << "  Fractional = " << NewEffErr[i]/NewEff[i] << endl;
		


		cout << "--------------------------------------------------------------------------------------------------------------" << endl;

		//	cout << "   NewEff = " << NewEffUp[i] << "     NewEffErr = " << NewEffErrUp[i] << "  Fractional = " << NewEffErrUp[i]/NewEffUp[i] << endl;

		//NewEffErr[i] = 0; //Remove Error on Efficiency Correction//
	}


	TFile * foutCorr;
	if(DoTnP == 0 && usemc==0)	foutCorr = new TFile(Form("%s/FinalFiles/weights%sPPCorrYield%sNoTnP%s.root",var_n.Data(),var_n.Data(),var_file.Data(),bsbpbins.Data()),"RECREATE");
	if(DoTnP == 1 && usemc==0)	foutCorr = new  TFile(Form("%s/FinalFiles/weights%sPPCorrYield%s%s.root",var_n.Data(),var_n.Data(),var_file.Data(),bsbpbins.Data()),"RECREATE");
	if(DoTnP == 0 && usemc==1)	foutCorr = new TFile(Form("%s/FinalFiles/weights%sPPCorrYield%sNoTnPMC%s.root",var_n.Data(),var_n.Data(),var_file.Data(),bsbpbins.Data()),"RECREATE");
	if(DoTnP == 1 && usemc==1)	foutCorr = new  TFile(Form("%s/FinalFiles/weights%sPPCorrYield%sMC%s.root",var_n.Data(),var_n.Data(),var_file.Data(),bsbpbins.Data()),"RECREATE");

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
	double XTemp;
	double YTemp;

	for(int i = 0; i < NBins;i++){
		
		Eff1D[i] = Eff1DHisvar->GetBinContent(i+1);
		Eff1DErr[i] = Eff1DHisvar->GetBinError(i+1);
	}

	hInvEff->SetMaximum(NewEff[0]*1.5);
	TCanvas *c = new TCanvas("c","c",700,700);
	c->cd();

	hInvEff->Draw("ep");
	
	if (usemc==0) {c->SaveAs(Form("%s/EffFinal/ReAnaEff_%dBins_%s.pdf",var_n.Data(),NBins,var_file.Data()));}
	else {c->SaveAs(Form("%s/EffFinal/ReAnaEff_%dBins_%s_MC.pdf",var_n.Data(),NBins,var_file.Data()));}
	
	Eff1DHisvar->Draw("ep");
if(BsBP==0){
	if (usemc==0) {c->SaveAs(Form("%s/EffFinal/ReAnaEff1D_%dBins_%s.pdf",var_n.Data(),NBins,var_file.Data()));}
	else {c->SaveAs(Form("%s/EffFinal/ReAnaEff1D_%dBins_%s_MC.pdf",var_n.Data(),NBins,var_file.Data()));}
}
	TH2D * invAcc2D=(TH2D *) fin1DEff->Get("invAcc2D");
	invAcc2D->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	invAcc2D->GetYaxis()->SetTitle("rapidity");

	invAcc2D->Draw("pcolz");
	if (usemc==0){c->SaveAs(Form("%s/EffFinal/acc_2Dmap%s.pdf",var_n.Data(),bsbpbins.Data()));}
	else {c->SaveAs(        Form("%s/EffFinal/acc_2Dmap%s_MC.pdf",var_n.Data(),bsbpbins.Data()));}


	DrawinvEff2D->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	DrawinvEff2D->GetYaxis()->SetTitle("|y|");
	DrawinvEff2D->GetZaxis()->SetLabelSize(0.02);

	DrawinvEff2D->Draw("pcolz");
	if (usemc==0){c->SaveAs(Form("%s/EffFinal/totaleff_2Dmap_%s%s.pdf",var_n.Data(),var_n.Data(), bsbpbins.Data()));}
	else {c->SaveAs(        Form("%s/EffFinal/totaleff_2Dmap_%s%s_MC.pdf",var_n.Data(),var_n.Data(), bsbpbins.Data()));}


	DrawinvEff2DY->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	DrawinvEff2DY->GetYaxis()->SetTitle("|y|");
	DrawinvEff2DY->GetZaxis()->SetLabelSize(0.02);
	c->SetLogz();
	DrawinvEff2DY->Draw("pcolz");

	if (usemc==0){c->SaveAs(Form("%s/EffFinal/totaleff_Fid_2Dmap_%s%s.pdf",var_n.Data(),var_n.Data(),bsbpbins.Data()));}
	else {c->SaveAs(        Form("%s/EffFinal/totaleff_Fid_2Dmap_%s%s_MC.pdf",var_n.Data(),var_n.Data(), bsbpbins.Data()));}

	

	

	double cmax = -10000;
	double cmin =  10000;
	for(int i = 0; i < NBins;i++){
		if (cmax < hEffInv1D->GetBinContent(i+1)){cmax=hEffInv1D->GetBinContent(i+1);}
		if (cmax < hInvEff->GetBinContent(i+1)){cmax=hInvEff->GetBinContent(i+1);}
		if (cmax < Eff1DHisvar->GetBinContent(i+1)){cmax=Eff1DHisvar->GetBinContent(i+1);}
		if (cmin > hEffInv1D->GetBinContent(i+1)){cmin=hEffInv1D->GetBinContent(i+1);}
		if (cmin > hInvEff->GetBinContent(i+1)){cmin=hInvEff->GetBinContent(i+1);}
		if (cmin > Eff1DHisvar->GetBinContent(i+1)){cmin=Eff1DHisvar->GetBinContent(i+1);}
	}
	Eff1DHisvar->GetYaxis()->SetRangeUser( cmin*0.9, cmax*1.1);

	hInvEff->SetMarkerColor(kRed+1);
	hInvEff->SetLineColor(kRed+1);
	hEffInv1D->SetMarkerColor(kBlue+1);
	hEffInv1D->SetLineColor(kBlue+1);
	Eff1DHisvar->Draw("ep");
	hInvEff->Draw("epsame");
	hEffInv1D->Draw("epsame");

	TLegend *leg = new TLegend(0.15,0.82,0.4,0.87,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.020);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(1);

	leg->AddEntry(hInvEff,"2D < 1/ #alpha #times #epsilon >^{-1}","pl");
	leg->AddEntry(hEffInv1D,"1D < 1/ #alpha #times #epsilon >^{-1}","pl");
	leg->AddEntry(Eff1DHisvar,"1D","pl");

	leg->Draw("SAME");

	if(BsBP==0){
		if (usemc==0){c->SaveAs(Form("%s/EffFinal/ReAnaEffcomp_%dBins_%s.pdf",var_n.Data(),NBins,var_file.Data()));}
		else {c->SaveAs(Form("%s/EffFinal/ReAnaEffcomp_%dBins_%s_MC.pdf",var_n.Data(),NBins,var_file.Data()));}
	}


	double cmax_weighted = -10000;
	double cmin_weighted =  10000;
	for(int i = 0; i < NBins;i++){
		if (cmax_weighted < hEffInv1D_weighted->GetBinContent(i+1)){cmax_weighted=hEffInv1D_weighted->GetBinContent(i+1);}
		if (cmax_weighted < hInvEff_weighted->GetBinContent(i+1)){cmax_weighted=hInvEff_weighted->GetBinContent(i+1);}
		if (cmax_weighted < Eff1DHisvar->GetBinContent(i+1)){cmax_weighted=Eff1DHisvar->GetBinContent(i+1);}
		if (cmin_weighted > hEffInv1D_weighted->GetBinContent(i+1)){cmin_weighted=hEffInv1D_weighted->GetBinContent(i+1);}
		if (cmin_weighted > hInvEff_weighted->GetBinContent(i+1)){cmin_weighted=hInvEff_weighted->GetBinContent(i+1);}
		if (cmin_weighted > Eff1DHisvar->GetBinContent(i+1)){cmin_weighted=Eff1DHisvar->GetBinContent(i+1);}
	}
	Eff1DHisvar->GetYaxis()->SetRangeUser( cmin_weighted*0.9, cmax_weighted*1.1);

	hInvEff_weighted->SetMarkerColor(kRed+1);
	hInvEff_weighted->SetLineColor(kRed+1);
	hEffInv1D_weighted->SetMarkerColor(kBlue+1);
	hEffInv1D_weighted->SetLineColor(kBlue+1);
	Eff1DHisvar->Draw("ep");
	hInvEff_weighted->Draw("epsame");
	hEffInv1D_weighted->Draw("epsame");

	TLegend *leg_weighted = new TLegend(0.15,0.82,0.4,0.87,NULL,"brNDC");
	leg_weighted->SetBorderSize(0);
	leg_weighted->SetTextSize(0.020);
	leg_weighted->SetTextFont(42);
	leg_weighted->SetFillStyle(0);
	leg_weighted->SetLineWidth(1);

	leg_weighted->AddEntry(hInvEff_weighted,"2D < 1/ #alpha #times #epsilon >^{-1}","pl");
	leg_weighted->AddEntry(hEffInv1D_weighted,"1D < 1/ #alpha #times #epsilon >^{-1}","pl");
	leg_weighted->AddEntry(Eff1DHisvar,"1D","pl");

	leg_weighted->Draw("SAME");

	if(BsBP==0){
	if (usemc==0){c->SaveAs(Form("%s/EffFinal/ReAnaEffcompweighted_%dBins_%s.pdf",var_n.Data(),NBins,var_file.Data()));}
	else {c->SaveAs(Form("%s/EffFinal/ReAnaEffcompweighted_%dBins_%s_MC.pdf",var_n.Data(),NBins,var_file.Data()));}
	}

	double cmax_both = -10000;
	double cmin_both =  10000;
	for(int i = 0; i < NBins;i++){
		if (cmax_both < hEffInv1D_weighted->GetBinContent(i+1)){cmax_both=hEffInv1D_weighted->GetBinContent(i+1);}
		if (cmax_both < hInvEff_weighted->GetBinContent(i+1)){cmax_both=hInvEff_weighted->GetBinContent(i+1);}
		if (cmax_both < hEffInv1D->GetBinContent(i+1)){cmax_both=hEffInv1D->GetBinContent(i+1);}
		if (cmax_both < hInvEff->GetBinContent(i+1)){cmax_both=hInvEff->GetBinContent(i+1);}

		if (cmin_both > hEffInv1D_weighted->GetBinContent(i+1)){cmin_both=hEffInv1D_weighted->GetBinContent(i+1);}
		if (cmin_both > hInvEff_weighted->GetBinContent(i+1)){cmin_both=hInvEff_weighted->GetBinContent(i+1);}
		if (cmin_both > hEffInv1D->GetBinContent(i+1)){cmin_both=hEffInv1D->GetBinContent(i+1);}
		if (cmin_both > hInvEff->GetBinContent(i+1)){cmin_both=hInvEff->GetBinContent(i+1);}
		
	}
	Eff1DHisvar->GetYaxis()->SetRangeUser( cmin_both*0.9, cmax_both*1.1);

	//hInvEff_weighted->SetMarkerColor(kRed+1);
	//hInvEff_weighted->SetLineColor(kRed+1);
	hInvEff_weighted_mass->SetMarkerColor(kRed+1);
	hInvEff_weighted_mass->SetLineColor(kRed+1);
	hInvEff->SetMarkerColor(kBlue+1);
	hInvEff->SetLineColor(kBlue+1);

	Eff1DHisvar->Draw("ep");
	//hInvEff_weighted->Draw("epsame");
	hInvEff_weighted_mass->Draw("epsame");
	hInvEff->Draw("epsame");

	TLegend *leg_both = new TLegend(0.15,0.82,0.4,0.87,NULL,"brNDC");
	leg_both->SetBorderSize(0);
	leg_both->SetTextSize(0.020);
	leg_both->SetTextFont(42);
	leg_both->SetFillStyle(0);
	leg_both->SetLineWidth(1);

	leg_both->AddEntry(hInvEff,"2D < 1/ #alpha #times #epsilon >^{-1}","pl");
	//leg_both->AddEntry(hInvEff_weighted,"weighted 2D ","pl");
	leg_both->AddEntry(hInvEff_weighted_mass,"2D weighted ","pl");
	leg_both->AddEntry(Eff1DHisvar,"1D ","pl");
	
	leg_both->Draw("SAME");

	if(BsBP==0){
	if (usemc==0){c->SaveAs(Form("%s/EffFinal/ReAnaEffcompboth_%dBins_%s.pdf",var_n.Data(),NBins,var_file.Data()));}
	else {c->SaveAs(Form("%s/EffFinal/ReAnaEffcompboth_%dBins_%s_MC.pdf",var_n.Data(),NBins,var_file.Data()));}
	}

	double tusk;
	double d4c;
	std::vector<double> cv1D;
	std::vector<double> cv2D;
	std::vector<double> cvweight1D;
	std::vector<double> cvweight2D;
	std::vector<double> cvboth2D;
	std::vector<double> cvboth2Dmass;
	std::vector<double> cvboth1D;
	
	std::vector<double> eff_val1d;
	std::vector<double> eff_val1dFG;
	std::vector<double> eff_val2d;
	std::vector<double> eff_valweight1d;
	std::vector<double> eff_valweight2d;
	std::vector<double> eff_valweight2dmass;
	std::vector<double> eff_err1d;
	std::vector<double> eff_err1dFG;
	std::vector<double> eff_err2d;
	std::vector<double> eff_errweight1d;
	std::vector<double> eff_errweight2d;
	std::vector<double> eff_errweight2dmass;
	std::vector<std::vector<double>> comparison_values;
	std::vector<std::vector<double>> comparison_values_weighted;
	std::vector<std::vector<double>> comparison_values_both;
	std::vector<std::vector<double>> eff_val;
	std::vector<std::vector<double>> eff_err;
	std::vector<std::vector<double>> eff_val_weight;
	std::vector<std::vector<double>> eff_err_weight;
	std::vector<std::vector<double>> eff_val_both;
	std::vector<std::vector<double>> eff_err_both;
	std::vector<std::string> labels = {"2D", "1D F-G"};
	std::vector<std::string> labelsboth = {"1D", "2D weighted"};
	std::vector<std::string> labelseff = {"1D","2D", "1D F-G"};
	std::vector<std::string> labelseffboth = {"2D", "2D weighted", "1D"};
	std::vector<std::string> col_name;
	std::vector<std::string> col_name_both;
	std::vector<std::string> col_name_eff;
	
	string name;
	TString whichvarname;
	col_name.push_back("Relative difference 1D vs:");
	col_name_both.push_back("Relative difference 2D vs:");
	col_name_eff.push_back("Eff values");
	if(whichvar==0){name="$<p_T<$"; whichvarname="pt";} else if(whichvar==1){name="$<y<$";whichvarname="y";} else if(whichvar==2){name="$<nTrks<$";whichvarname="nMult";}
	for(int i=0;i<NBins;i++){
		std::ostringstream clabel;
		clabel<<ptBins[i]<<name<<ptBins[i+1];
		std::string label1 = clabel.str();
		col_name.push_back(label1);
		col_name_eff.push_back(label1);
		col_name_both.push_back(label1);

		cv1D.push_back(abs(1/NewEffReal[i]-Eff1D[i])/Eff1D[i]*100);
		cv2D.push_back(abs(1/NewEff[i]-Eff1D[i])/Eff1D[i]*100);

		cvweight1D.push_back(abs(1/NewEffReal_weighted[i]-Eff1D[i])/Eff1D[i]*100);
		cvweight2D.push_back(abs(1/NewEff_weighted[i]-Eff1D[i])/Eff1D[i]*100);


		cvboth1D.push_back(abs(Eff1D[i]-1/NewEff[i])*NewEff[i]*100);
		cvboth2D.push_back(abs(1/NewEff_weighted[i]-1/NewEff[i])*NewEff[i]*100);
		cvboth2Dmass.push_back(abs(1/NewEff_weighted_mass[i]-1/NewEff[i])*NewEff[i]*100);

		eff_val1d.push_back(Eff1D[i]);
		eff_val2d.push_back(1/NewEff[i]);
		eff_val1dFG.push_back(1/NewEffReal[i]);

		eff_err1d.push_back(Eff1DErr[i]);
		eff_err2d.push_back(NewEffErr[i]/(NewEff[i]*NewEff[i]));
		eff_err1dFG.push_back(NewEffRealErr[i]/(NewEffReal[i]*NewEffReal[i]));

		
		eff_valweight2d.push_back(1/NewEff_weighted[i]);
		eff_valweight2dmass.push_back(1/NewEff_weighted_mass[i]);
		eff_valweight1d.push_back(1/NewEffReal_weighted[i]);

		
		eff_errweight2d.push_back(NewEffErr_weighted[i]/(NewEff_weighted[i]*NewEff_weighted[i]));
		eff_errweight2dmass.push_back(NewEffErr_weighted_mass[i]/(NewEff_weighted_mass[i]*NewEff_weighted_mass[i]));
		eff_errweight1d.push_back(NewEffRealErr[i]/(NewEffReal_weighted[i]*NewEffReal_weighted[i]));
		
	}
	
	comparison_values.push_back(cv2D);
	comparison_values.push_back(cv1D);
	comparison_values_weighted.push_back(cvweight2D);
	comparison_values_weighted.push_back(cvweight1D);

	comparison_values_both.push_back(cvboth1D);
	comparison_values_both.push_back(cvboth2Dmass);

	eff_val.push_back(eff_val1d);
	eff_val.push_back(eff_val2d);
	eff_val.push_back(eff_val1dFG);
	eff_err.push_back(eff_err1d);
	eff_err.push_back(eff_err2d);
	eff_err.push_back(eff_err1dFG);

	eff_val_weight.push_back(eff_val1d);
	eff_val_weight.push_back(eff_valweight2d);
	eff_val_weight.push_back(eff_valweight1d);
	eff_err_weight.push_back(eff_err1d);
	eff_err_weight.push_back(eff_errweight2d);
	eff_err_weight.push_back(eff_errweight1d);

	eff_val_both.push_back(eff_val2d);
	eff_val_both.push_back(eff_valweight2dmass);
	eff_val_both.push_back(eff_val1d);
	
	eff_err_both.push_back(eff_err2d);
	eff_err_both.push_back(eff_errweight2dmass);
	eff_err_both.push_back(eff_err1d);
	
	
	

	latex_table(Form("1D2Dcomparisons_%s",whichvarname.Data()), NBins+1,  3,  col_name , labels , comparison_values, "1D vs 2D efficiency comparisons",0);

	

	latex_table(Form("Effyields_%s",whichvarname.Data()), NBins+1,  4,  col_name_eff , labelseff , eff_val , "Eff Values per bin",1, eff_err);

	
	latex_table(Form("1D2Dcomparisons_weighted_%s",whichvarname.Data()), NBins+1,  3,  col_name , labels , comparison_values_weighted, "1D vs 2D efficiency comparisons",0);
	
	latex_table(Form("Effyields_weighted_%s",whichvarname.Data()), NBins+1,  4,  col_name_eff , labelseff , eff_val_weight , "Eff Values per bin",1, eff_err_weight);

	
	latex_table(Form("1D2Dcomparisons_both_%s",whichvarname.Data()), NBins+1,  3,  col_name_both , labelsboth , comparison_values_both, "2D method comparisons",0);
	
	latex_table(Form("Effyields_both_%s",whichvarname.Data()), NBins+1,  4,  col_name_eff , labelseffboth , eff_val_both , "Eff Values per bin",1, eff_err_both);

	
	

	std::vector<std::string> filetype ={"_check.aux", "_check.log",".tex","_check.tex"};
	for (int j=0;j<(int)(filetype.size());j++){
				rename(("1D2Dcomparisons_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n.Data())+"/Trash/1D2Dcomparisons_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				rename(("1D2Dcomparisons_weighted_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n.Data())+"/Trash/1D2Dcomparisons_weighted_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				rename(("1D2Dcomparisons_both_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n.Data())+"/Trash/1D2Dcomparisons_both_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				rename(("Effyields_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n.Data())+"/Trash/Effyields_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				rename(("Effyields_weighted_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n.Data())+"/Trash/Effyields_weighted_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				rename(("Effyields_both_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),(std::string (var_n.Data())+"/Trash/Effyields_both_"+std::string (whichvarname.Data())+filetype[j]).c_str());
				
			}
	rename(("1D2Dcomparisons_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n.Data())+"/EffFinal/1D2Dcomparisons_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("1D2Dcomparisons_weighted_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n.Data())+"/EffFinal/1D2Dcomparisons_weighted_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("1D2Dcomparisons_both_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n.Data())+"/EffFinal/1D2Dcomparisons_both_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("Effyields_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n.Data())+"/EffFinal/Effyields_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("Effyields_weighted_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n.Data())+"/EffFinal/Effyields_weighted_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("Effyields_both_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),(std::string (var_n.Data())+"/EffFinal/Effyields_both_"+std::string (whichvarname.Data())+"_check.pdf").c_str());



	foutCorr->cd();
	hInvEff->Write();
	Eff1DHisvar->Write();
	foutCorr->Close();

}