#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TSystem.h"
#include <fstream>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "scale.h"
#include <iostream>
#include <fstream>
#include "../../../henri2022/parameter.h" 
#include "CMS_lumi.C" 

using namespace std;
using std::cout;
using std::endl;

// MOST OF THE VARIABLES SAY BP BUT ARE WORKING FOR Bs AS WELL 
// THIS CODE CAN BE EASILY EXTENDED TO ACCOUNT FOR A FUTURE MESON... 
// JUST FOLLOW THE LOGIC
void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, 
		std::vector<std::vector<double> > numbers , std::string caption, int ori){
	
	std::ofstream file_check;
	std::ofstream file;
                                                                                    
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

	if (ori==0){
		for(int i=1; i<n_lin; i++)
		{	
			file << labels[i-1] << " & ";
			file_check << labels[i-1] << " & ";
			for(int c=1; c<n_col-1; c++){
				file << std::setprecision(3)<<  numbers[c-1][i-1]<< " & ";
				file_check << std::setprecision(3)<<  numbers[c-1][i-1]<< " & ";
										}
				file << std::setprecision(3)<<  numbers[n_col-2][i-1]<< "  \\\\" << std::endl;
				file_check << std::setprecision(3)<<  numbers[n_col-2][i-1]<< "  \\\\" << std::endl; 
		}
	} else {
		for(int i=1; i<n_lin; i++)
		{	
			file << labels[i-1] << " & ";
			file_check << labels[i-1] << " & ";
			for(int c=1; c<n_col-1; c++){
				file << std::setprecision(3)<< numbers[i-1][c-1]<< " \\% & ";
				file_check << std::setprecision(3)<< numbers[i-1][c-1]<< " \\% & ";
										}
				file << std::setprecision(3)<<  numbers[i-1][n_col-2]<< " \\% \\\\" << std::endl;
				file_check << std::setprecision(3)<<  numbers[i-1][n_col-2]<< " \\% \\\\" << std::endl; 
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



void Bmeson_Comparisons(int meson_n, int whichvar){
                          
    constexpr bool fidFONLL = true;
	TString B_m ;
	TString t_tree ;
	TString b_m;
	int NBins = 7;
	int NBinsLow = 0  ;
  	int NBinsHigh = 0 ;
	int NBins2015 = 0 ;
	int lowstart  ;
	int lowend    ;
	TString var_n;
	TString var_N;
	TString var_l;
	
	if(meson_n == 0){t_tree = "ntKp"; B_m = "BP"; b_m = "bp";} 
	if(meson_n == 1){t_tree = "ntphi"; 	B_m = "Bs"; b_m = "bs";}

	if(meson_n == 0 && whichvar==0){
		NBins = nptBinsBP;
		var_n="pt";
		var_N="PT";
		var_l="p_{T} [GeV/c]";
		lowstart = 0;
		lowend = 1;
		NBinsLow = 2 ;
		NBinsHigh = 5;
		NBins2015 = 5;
	} 
	if(meson_n == 1 && whichvar==0){
		var_n="pt";
		var_N="PT";
		var_l="p_{T} [GeV/c]";
		lowstart = 0;
		lowend = 0;
		NBins = nptBins;
		NBinsLow = 1 ;
		NBinsHigh = 3;
		NBins2015 = 3;
	}
	if(whichvar==1){
		NBins = nyBins_both;
		var_n="y";
		var_N="Y";
		var_l="|y|";
		lowstart = 0;   //0-0.5-1-1.5
		lowend = 2;
		NBinsLow = 3;
		NBinsHigh = 2;
	}
	if(whichvar==2){
		NBins = nmBins_both;
		var_n="Mult";
		var_N="Mult";
		var_l="Mult";
		lowstart = 100;
		lowend = -1;
		NBinsLow = 0 ;         
		NBinsHigh = 0;
	}

	gSystem->mkdir("Plots/", true);
	gSystem->mkdir(Form("Plots/%s",B_m.Data()), true);
	
	TString InfileB = Form("../../../EffAna/%s/FinalFiles/%sPPCorrYield%s.root",B_m.Data(),B_m.Data(),var_N.Data());
	TFile * FileB= new TFile(InfileB.Data());
	
	double ptBins[NBins+1];
	for(int i = 0; i < NBins + 1; i++){
		if (whichvar==0){
			if (meson_n==0){ ptBins[i] =  ptbinsvecBP[i];} 
			else if (meson_n == 1){ptBins[i] =  ptbinsvec[i];}
		}
		else if (whichvar==1){ ptBins[i] =  ybinsvec[i]; }          
		else if (whichvar==2){ ptBins[i] =  nmbinsvec[i];}
	}
	//center of the bin and its left and right margins
	//THIS NEED TO BE UPDATED, THIS IS WRONG! WE WANT THE MEAN OF THE BIN CENTER!
	float BPXsecPPX[NBins];
	float BXSecPPXErrUp[NBins] ;
	float BXSecPPXErrDown[NBins] ;

	for( int c=0; c < NBins; c++){
		BPXsecPPX[c]=(ptBins[c]+ptBins[c+1])/2;

		BXSecPPXErrUp[c]=(abs(ptBins[c+1]-ptBins[c]))/2;
		BXSecPPXErrDown[c]=(abs(ptBins[c+1]-ptBins[c]))/2;
	}
	//THIS NEED TO BE UPDATED, THIS IS WRONG! WE WANT THE MEAN OF THE BIN CENTER!


  // cross section with 2D Map eff correction
	float BPXsecPPY2D[NBins];
	float BPXSecPPY2DErrUp[NBins];
	float BPXSecPPY2DErrDown[NBins];
	float BPXSecPPY2DErrUpRatio[NBins];
	float BPXSecPPY2DErrDownRatio[NBins];
  
  // cross section with pT < 10 scaled to full y
	float BPXsecPPY2DScaled[NBins];
	float BPXSecPPY2DErrUpScaled[NBins];
	float BPXSecPPY2DErrDownScaled[NBins];

	TH1D * BCross2D = (TH1D *) FileB->Get("hPtSigma");
	for(int i = 0; i < NBins; i++){
		BPXsecPPY2D[i] = BCross2D->GetBinContent(i+1);
		BPXSecPPY2DErrUp[i] = BCross2D->GetBinError(i+1);
		BPXSecPPY2DErrDown[i] = BCross2D->GetBinError(i+1);
		BPXSecPPY2DErrUpRatio[i] = BPXSecPPY2DErrUp[i] / BPXsecPPY2D[i];
		BPXSecPPY2DErrDownRatio[i] = BPXSecPPY2DErrDown[i] / BPXsecPPY2D[i];
		
		BPXsecPPY2DScaled[i] = BCross2D->GetBinContent(i+1);
		BPXSecPPY2DErrUpScaled[i] = BCross2D->GetBinError(i+1);
		BPXSecPPY2DErrDownScaled[i] = BCross2D->GetBinError(i+1);
	}
	
	// cross section with 2D Map eff correction
	float BPXsecPPY1D[NBins];
	float BPXSecPPY1DErrUp[NBins];
	float BPXSecPPY1DErrDown[NBins];
	float BPXSecPPY1DErrUpRatio[NBins];
	float BPXSecPPY1DErrDownRatio[NBins];

	TH1D * BCross1D = (TH1D *) FileB->Get("CorrDiffHisBin");
	for(int i = 0; i < NBins; i++){
		BPXsecPPY1D[i] = BCross1D->GetBinContent(i+1);
		BPXSecPPY1DErrUp[i] = BCross1D->GetBinError(i+1);
		BPXSecPPY1DErrDown[i] = BCross1D->GetBinError(i+1);
		BPXSecPPY1DErrUpRatio[i] = BPXSecPPY1DErrUp[i] / BPXsecPPY1D[i];
		BPXSecPPY1DErrDownRatio[i] = BPXSecPPY1DErrDown[i] / BPXsecPPY1D[i];
		
	}

	//Syst Add Up PP//
  	TString errorFile = Form("../../../2DMapSyst/OutFiles/%sError2D_%s.root", B_m.Data(),var_n.Data());
  	TFile fError(errorFile);

	TH1D * TnPSyst = (TH1D *) fError.Get("TnPSyst");
	TH1D * BptSyst = (TH1D *) fError.Get("BptSyst");
	TH1D * MCDataSyst = (TH1D *) fError.Get("MCDataSyst");
  	if (!MCDataSyst) MCDataSyst = (TH1D *) fError.Get("BDTSyst");

	TString errorFile1D = Form("../../../1DMapSyst/OutFiles/%sError1D_%s.root", B_m.Data(),var_n.Data());
  	TFile fError1D(errorFile1D);

	TH1D * TnPSyst1D = (TH1D *) fError1D.Get("TnPSyst");
	TH1D * BptSyst1D = (TH1D *) fError1D.Get("BptSyst");
	TH1D * MCDataSyst1D = (TH1D *) fError1D.Get("BDTSyst");

	TString pdfErrorFile = Form("../../../syst_error/%s_pdf_%s.root",b_m.Data(),var_n.Data());
	TFile fPdfError(pdfErrorFile);
	TGraph* pdfSyst = (TGraph *) fPdfError.Get(Form("%s_error",b_m.Data()));
	
	TString trackSelErrorFile = Form("../../../syst_error/syst_track_sel_%s.root",var_n.Data());
	TFile fTrackSelError(trackSelErrorFile);
	TGraph* trackSelSyst = (TGraph *) fTrackSelError.Get(Form("%s_track_sel_error", b_m.Data()));

	TString trackSelErrorFile1D = Form("../../../syst_error/syst_track_sel_%s_1D.root",var_n.Data());
	TFile fTrackSelError1D(trackSelErrorFile1D);
	TGraph* trackSelSyst1D = (TGraph *) fTrackSelError1D.Get(Form("%s_track_sel_error", b_m.Data()));

	float BPXSecPPY2DSystUp[NBins];
	float BPXSecPPY2DSystDown[NBins];
	float BPXSecPPY1DSystUp[NBins];
	float BPXSecPPY1DSystDown[NBins];
	float BPXSecPPYSystUpScaled[NBins];
	float BPXSecPPYSystDownScaled[NBins];

  // percent error

  	double B_nu;
  	if (meson_n == 0){ B_nu = 2.4 ;}
  	else { B_nu = 4.8;}
	double BPTrackingSyst[NBins];
	for( int c=0; c < NBins; c++){ BPTrackingSyst[c]= B_nu ;}
	
	float BPMCDataSyst[NBins];
	float BPPDFSyst[NBins];
	float BPTrackSelSyst[NBins];
	float BPPtShapeSyst[NBins];
	float BPTnPSystDown[NBins];
	float BPTnPSystUp[NBins];
	float BP1DMCDataSyst[NBins];
	float BP1DPDFSyst[NBins];
	float BP1DTrackSelSyst[NBins];
	float BP1DPtShapeSyst[NBins];
	float BP1DTnPSystDown[NBins];
	float BP1DTnPSystUp[NBins];

  // Get systematics from input files
    for (auto ibin = 0; ibin < NBins; ++ibin){
		BPMCDataSyst[ibin] = MCDataSyst->GetBinContent(ibin + 1);
		BPPtShapeSyst[ibin] = BptSyst->GetBinContent(ibin + 1);
		BPTnPSystDown[ibin] = TnPSyst->GetBinContent(ibin + 1);

		BP1DMCDataSyst[ibin] = MCDataSyst1D->GetBinContent(ibin + 1);
		BP1DPtShapeSyst[ibin] = BptSyst1D->GetBinContent(ibin + 1);
		BP1DTnPSystDown[ibin] = TnPSyst1D->GetBinContent(ibin + 1);
		
		// TnP systematics are symmetric in the binned pT case
		BPTnPSystUp[ibin] = BPTnPSystDown[ibin];
		BPPDFSyst[ibin] = pdfSyst->GetY()[ibin];
		BPTrackSelSyst[ibin] = trackSelSyst->GetY()[ibin];
		BP1DTrackSelSyst[ibin] = trackSelSyst1D->GetY()[ibin];

		BP1DTnPSystUp[ibin] = BP1DTnPSystDown[ibin];
		BP1DPDFSyst[ibin] = pdfSyst->GetY()[ibin];
		//BP1DTrackSelSyst[ibin] = trackSelSyst->GetY()[ibin];
  	}

  // RMS of all the errors

	float BP2DTotalSystDownRatio[NBins];
	float BP2DTotalSystUpRatio[NBins];
	float BP1DTotalSystDownRatio[NBins];
	float BP1DTotalSystUpRatio[NBins];

	for(int i = 0; i < NBins; i++){
		BP2DTotalSystDownRatio[i] = TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BPMCDataSyst[i], 2) +
                                          TMath::Power(BPPDFSyst[i], 2) + TMath::Power(BPTrackSelSyst[i], 2) +
                                          TMath::Power(BPPtShapeSyst[i], 2) + TMath::Power(BPTnPSystDown[i], 2)) / 100;

    	BP2DTotalSystUpRatio[i] = TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BPMCDataSyst[i], 2) +
                                        TMath::Power(BPPDFSyst[i], 2) + TMath::Power(BPTrackSelSyst[i], 2) +
                                        TMath::Power(BPPtShapeSyst[i], 2) + TMath::Power(BPTnPSystUp[i], 2)) / 100;
		
		BP1DTotalSystDownRatio[i] = TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BP1DMCDataSyst[i], 2) +
                                          TMath::Power(BP1DPDFSyst[i], 2) + TMath::Power(BP1DTrackSelSyst[i], 2) +
                                          TMath::Power(BP1DPtShapeSyst[i], 2) + TMath::Power(BP1DTnPSystDown[i], 2)) / 100;

    	BP1DTotalSystUpRatio[i] = TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BP1DMCDataSyst[i], 2) +
                                        TMath::Power(BP1DPDFSyst[i], 2) + TMath::Power(BP1DTrackSelSyst[i], 2) +
                                        TMath::Power(BP1DPtShapeSyst[i], 2) + TMath::Power(BP1DTnPSystUp[i], 2)) / 100;
	}

  // global uncertainty from branching ratio and luminosity
  // Fixed, copied from the paper draft

	double numb;
	if(meson_n == 0){numb = 0.035;}
	else {numb = 0.077;}
  	vector<float> globUncert(NBins, numb);
	for(int i = 0; i < NBins; i++){
		BPXSecPPY2DSystUp[i] = BPXsecPPY2D[i] * BP2DTotalSystUpRatio[i];
		BPXSecPPY2DSystDown[i] = BPXsecPPY2D[i] * BP2DTotalSystDownRatio[i];
		BPXSecPPY1DSystUp[i] = BPXsecPPY1D[i] * BP1DTotalSystUpRatio[i];
		BPXSecPPY1DSystDown[i] = BPXsecPPY1D[i] * BP1DTotalSystDownRatio[i];
    	

		}

// CREATE THE CANVAS and the pads
	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd(); 
	if (whichvar==0){c->SetLogy();}
	c->SetLeftMargin(0.15);

	//Setup histograms for different purposs
	TH2D * HisEmpty;
	if(meson_n == 0 && whichvar==0) {HisEmpty = new TH2D("HisEmpty","",100,5,60,100,300.0,2000000);} 
	if(meson_n == 1 && whichvar==0) {HisEmpty = new TH2D("HisEmpty","",100,7,50,100,300.0,2000000);}

	if(meson_n == 0 && whichvar==1) {HisEmpty = new TH2D("HisEmpty","",100,0,2.4,100,1250000.0,5000000);}
	if(meson_n == 1 && whichvar==1) {HisEmpty = new TH2D("HisEmpty","",100,0,2.4,100,220000.0,600000);}

	if(meson_n == 0 && whichvar==2) {HisEmpty = new TH2D("HisEmpty","",100,0,100,100,0,4200000);}   // need to adjust range for when we have nmult results
	if(meson_n == 1 && whichvar==2) {HisEmpty = new TH2D("HisEmpty","",100,0,100,100,0,600000);}

	HisEmpty->GetXaxis()->SetTitle(var_l.Data());
	if (whichvar==0) {HisEmpty->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb c/GeV]");}
	else {HisEmpty->GetYaxis()->SetTitle(Form("d#sigma/d%s [pb c/GeV]",var_n.Data()));}
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	//HisEmpty->GetYaxis()->SetTitleOffset(1.8);
	//HisEmpty->GetXaxis()->SetTitleOffset(1.3);	

	TH2D * HisEmpty2;
	if (meson_n == 0){
		HisEmpty2 = new TH2D("HisEmpty2","",100,5,60,100,300.0,3000000);
		HisEmpty2->GetXaxis()->SetTitle(var_l.Data());
	} else {	
		HisEmpty2 = new TH2D("HisEmpty2","",100,7,50,100,300.0,3000000);
		HisEmpty2->GetXaxis()->SetTitle(var_l.Data());
	}
	if (whichvar==0) {HisEmpty2->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb c/GeV]");}
	else {HisEmpty2->GetYaxis()->SetTitle(Form("d#sigma/d%s [pb c/GeV]",var_n.Data()));}
		HisEmpty2->GetXaxis()->CenterTitle();
		HisEmpty2->GetYaxis()->CenterTitle();

// CREATE THE CANVAS and the pads

  // separate plots for different fiducial regions
  	vector<float> BPXsecPPXLow ;
  	vector<float> BPXsecPPXHigh ;
	vector<float> BPXsecPPXErrLow ;
	vector<float> BPXsecPPXErrHigh ;

  	vector<float> BPXsecPPYLow ;
  	vector<float> BPXsecPPYHigh ;
	vector<float> BPXsecPPYErrDownLow ;
	vector<float> BPXsecPPYErrDownHigh ;
	vector<float> BPXsecPPYErrUpLow ;
	vector<float> BPXsecPPYErrUpHigh ;
	vector<float> BPYSystDown_low ;
	vector<float> BPYSystDown_high ;
	vector<float> BPYSystUp_low ;
	vector<float> BPYSystUp_high ;

	vector<float> BP1DXsecPPYLow ;
  	vector<float> BP1DXsecPPYHigh ;
	vector<float> BP1DXsecPPYErrDownLow ;
	vector<float> BP1DXsecPPYErrDownHigh ;
	vector<float> BP1DXsecPPYErrUpLow ;
	vector<float> BP1DXsecPPYErrUpHigh ;
	vector<float> BP1DYSystDown_low ;
	vector<float> BP1DYSystDown_high ;
	vector<float> BP1DYSystUp_low ;
	vector<float> BP1DYSystUp_high ;

	for (int i=0 ; i<NBins; i++){
		if (i >= lowstart && i <= lowend){
			BPXsecPPXLow.push_back(BPXsecPPX[i]);
			BPXsecPPXErrLow.push_back(BXSecPPXErrDown[i]);
			BPXsecPPYLow.push_back(BPXsecPPY2D[i]);
			BPXsecPPYErrDownLow.push_back(BPXSecPPY2DErrDown[i]);
			BPXsecPPYErrUpLow.push_back(BPXSecPPY2DErrUp[i]);
			BPYSystDown_low.push_back(BPXSecPPY2DSystDown[i]);
			BPYSystUp_low.push_back(BPXSecPPY2DSystUp[i]);

			BP1DXsecPPYLow.push_back(BPXsecPPY1D[i]);
			BP1DXsecPPYErrDownLow.push_back(BPXSecPPY1DErrDown[i]);
			BP1DXsecPPYErrUpLow.push_back(BPXSecPPY1DErrUp[i]);
			BP1DYSystDown_low.push_back(BPXSecPPY1DSystDown[i]);
			BP1DYSystUp_low.push_back(BPXSecPPY1DSystUp[i]);
		} else {
			BPXsecPPXHigh.push_back(BPXsecPPX[i]);
			BPXsecPPXErrHigh.push_back(BXSecPPXErrDown[i]);
			BPXsecPPYHigh.push_back(BPXsecPPY2D[i]);
			BPXsecPPYErrDownHigh.push_back(BPXSecPPY2DErrDown[i]);
			BPXsecPPYErrUpHigh.push_back(BPXSecPPY2DErrUp[i]);
			BPYSystDown_high.push_back(BPXSecPPY2DSystDown[i]);
			BPYSystUp_high.push_back(BPXSecPPY2DSystUp[i]);

			BP1DXsecPPYHigh.push_back(BPXsecPPY1D[i]);
			BP1DXsecPPYErrDownHigh.push_back(BPXSecPPY1DErrDown[i]);
			BP1DXsecPPYErrUpHigh.push_back(BPXSecPPY1DErrUp[i]);
			BP1DYSystDown_high.push_back(BPXSecPPY1DSystDown[i]);
			BP1DYSystUp_high.push_back(BPXSecPPY1DSystUp[i]);
		}
	}

	float zero[NBinsLow];
	for (int i=0 ; i<NBinsLow; i++){zero[i]=0;}

	TGraphAsymmErrors *BPRAAGraph_low_just_marker   = new TGraphAsymmErrors(NBinsLow, BPXsecPPXLow.data(), BPXsecPPYLow.data() ,zero, zero, zero, zero);
	TGraphAsymmErrors *BP1DRAAGraph_low_just_marker = new TGraphAsymmErrors(NBinsLow, BPXsecPPXLow.data(), BP1DXsecPPYLow.data() ,zero, zero, zero, zero);

	TGraphAsymmErrors *BPRAAGraph_low = new TGraphAsymmErrors(NBinsLow , BPXsecPPXLow.data() , BPXsecPPYLow.data() , BPXsecPPXErrLow.data() , BPXsecPPXErrLow.data() , BPXsecPPYErrDownLow.data() , BPXsecPPYErrUpLow.data());
	TGraphAsymmErrors *BPRAAGraph     = new TGraphAsymmErrors(NBinsHigh, BPXsecPPXHigh.data(), BPXsecPPYHigh.data(), BPXsecPPXErrHigh.data(), BPXsecPPXErrHigh.data(), BPXsecPPYErrDownHigh.data(), BPXsecPPYErrUpHigh.data());     
	TGraphAsymmErrors *BPRAAGraphSyst_low  = new TGraphAsymmErrors(NBinsLow , BPXsecPPXLow.data() , BPXsecPPYLow.data() , BPXsecPPXErrLow.data() , BPXsecPPXErrLow.data() , BPYSystDown_low.data() , BPYSystUp_low.data());                 											
	TGraphAsymmErrors *BPRAAGraphSyst      = new TGraphAsymmErrors(NBinsHigh, BPXsecPPXHigh.data(), BPXsecPPYHigh.data(), BPXsecPPXErrHigh.data(), BPXsecPPXErrHigh.data(), BPYSystDown_high.data(), BPYSystUp_high.data());                 											
	
	TGraphAsymmErrors *BP1DRAAGraph_low = new TGraphAsymmErrors(NBinsLow , BPXsecPPXLow.data() , BP1DXsecPPYLow.data() , BPXsecPPXErrLow.data() , BPXsecPPXErrLow.data() , BP1DXsecPPYErrDownLow.data() , BP1DXsecPPYErrUpLow.data());
	TGraphAsymmErrors *BP1DRAAGraph     = new TGraphAsymmErrors(NBinsHigh, BPXsecPPXHigh.data(), BP1DXsecPPYHigh.data(), BPXsecPPXErrHigh.data(), BPXsecPPXErrHigh.data(), BP1DXsecPPYErrDownHigh.data(), BP1DXsecPPYErrUpHigh.data());     
	TGraphAsymmErrors *BP1DRAAGraphSyst_low  = new TGraphAsymmErrors(NBinsLow , BPXsecPPXLow.data() , BP1DXsecPPYLow.data() , BPXsecPPXErrLow.data() , BPXsecPPXErrLow.data() , BP1DYSystDown_low.data() , BP1DYSystUp_low.data());                 											
	TGraphAsymmErrors *BP1DRAAGraphSyst      = new TGraphAsymmErrors(NBinsHigh, BPXsecPPXHigh.data(), BP1DXsecPPYHigh.data(), BPXsecPPXErrHigh.data(), BPXsecPPXErrHigh.data(), BP1DYSystDown_high.data(), BP1DYSystUp_high.data());
  // separate plots for different fiducial regions
 
  	cout << endl << "-------------------------------------------------------  "<< Form("%s meson Xsection", B_m.Data()) <<"  -------------------------------------------------------" << endl;

	for(int i=0;i<NBins;i++){		
		cout << "BIN " <<              Form("[%.1f,%.1f]  ",ptBins[i],ptBins[i+1]) << Form("%.0f #pm (STATup) %.0f #pm (SYSTup) %.0f #pm (STATdown) %.0f #pm (SYSTdown) %.0f ",BPXsecPPY2D[i], BPXSecPPY2DErrUp[i], BPXSecPPY2DSystUp[i], BPXSecPPY2DErrDown[i], BPXSecPPY2DSystDown[i]) << endl;
		cout << "(normalized) BIN " << Form("[%.1f,%.1f]  ",ptBins[i],ptBins[i+1]) << Form("%.0f #pm (STATup) %.1f #pm (SYSTup) %.1f #pm (STATdown) %.1f #pm (SYSTdown) %.1f ",BPXsecPPY2D[i], 100*BPXSecPPY2DErrUp[i]/BPXsecPPY2D[i], 100*BPXSecPPY2DSystUp[i]/BPXsecPPY2D[i], 100*BPXSecPPY2DErrDown[i]/BPXsecPPY2D[i], 100*BPXSecPPY2DSystDown[i]/BPXsecPPY2D[i]) << endl;
	}
 
 	cout<< endl << "-------------------------------------------------------  "<< Form("%s meson Xsection", B_m.Data()) <<"  -------------------------------------------------------" << endl;




// CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection 
// CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection 
	
	BPRAAGraph->SetMarkerStyle(20);
	BPRAAGraph->SetMarkerSize(1);
	BPRAAGraph_low_just_marker ->SetMarkerStyle(20);
	BPRAAGraph_low_just_marker ->SetMarkerSize(0.9);
	BPRAAGraph_low_just_marker ->SetMarkerColor(kWhite);
	BPRAAGraph_low ->SetMarkerStyle(24);
	BPRAAGraph_low ->SetMarkerSize(1);

	if (meson_n==0){
		BPRAAGraph_low ->SetMarkerColor(kGreen+2);
		BPRAAGraph_low ->SetLineColor(kGreen+2);
		BPRAAGraphSyst_low ->SetFillColorAlpha(kGreen-7,0.5);
		BPRAAGraphSyst_low ->SetLineColor(kGreen-7);
		BPRAAGraph ->SetMarkerColor(kGreen+2);
		BPRAAGraph ->SetLineColor(kGreen+2);
		BPRAAGraphSyst ->SetFillColorAlpha(kGreen-7,0.5);
		BPRAAGraphSyst ->SetLineColor(kGreen-7);
	} else {
		BPRAAGraph_low ->SetMarkerColor(kBlue+2);
		BPRAAGraph_low ->SetLineColor(kBlue+2);
		BPRAAGraphSyst_low ->SetFillColorAlpha(kBlue-3,0.5);
		BPRAAGraphSyst_low ->SetLineColor(kBlue-3);
		BPRAAGraph->SetLineColor(kBlue+2);
		BPRAAGraph->SetMarkerColor(kBlue+2);
		BPRAAGraphSyst->SetFillColorAlpha(kBlue-3,0.5);
		BPRAAGraphSyst->SetLineColor(kBlue-3);
	}

	BP1DRAAGraph->SetMarkerStyle(20);
	BP1DRAAGraph->SetMarkerSize(1);
	BP1DRAAGraph_low_just_marker ->SetMarkerStyle(20);
	BP1DRAAGraph_low_just_marker ->SetMarkerSize(0.9);
	BP1DRAAGraph_low_just_marker ->SetMarkerColor(kWhite);
	BP1DRAAGraph_low ->SetMarkerStyle(24);
	BP1DRAAGraph_low ->SetMarkerSize(1);
	BP1DRAAGraph_low ->SetMarkerColor(kPink+2);
	BP1DRAAGraph_low ->SetLineColor(kPink+2);
	BP1DRAAGraphSyst_low ->SetFillColorAlpha(kPink-3,0.5);
	BP1DRAAGraphSyst_low ->SetLineColor(kPink-3);
	BP1DRAAGraph->SetLineColor(kPink+2);
	BP1DRAAGraph->SetMarkerColor(kPink+2);
	BP1DRAAGraphSyst->SetFillColorAlpha(kPink-3,0.5);
	BP1DRAAGraphSyst->SetLineColor(kPink-3);
	
	HisEmpty->Draw();

	BPRAAGraphSyst_low->Draw("5same");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph_low->Draw("epSAME");
	BPRAAGraph_low_just_marker->Draw("epSAME");

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.025); 
	lat->SetTextFont(42);

	lat->DrawLatex(0.15,0.91 , "CMS work in progress");
	if (meson_n == 0) {lat->DrawLatex(0.6,0.7 ,Form("2017 pp global Unc. #pm %.1f%%",3.5));} 
	else {	lat->DrawLatex(0.6,0.7,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;}

    TLegend* leged;
	if (meson_n==0) {leged = new TLegend(0.65,0.67,0.9,0.75,NULL,"brNDC");}
	else {leged = new TLegend(0.65,0.6,0.9,0.68,NULL,"brNDC");}
	leged->SetBorderSize(0);
	leged->SetFillStyle(0);
	if (whichvar==0){
		leged->AddEntry((TObject*)0, "y region:", "");
		leged->AddEntry(BPRAAGraph,"|y|<2.4","P");
		leged->AddEntry(BPRAAGraph_low,"|y|>1.5","P");
	} 
	if (whichvar==1){
		leged->AddEntry((TObject*)0, "p_{T} region:", "");
		if (meson_n==0){
			leged->AddEntry(BPRAAGraph,"5<p_{T}<60 GeV/c","P");
		} else {
			leged->AddEntry(BPRAAGraph,"7<p_{T}<50 GeV/c","P");
		}
		leged->AddEntry(BPRAAGraph_low,"p_{T}>10 GeV/c","P");
	} 
	leged->SetTextSize(0.022);
	leged->Draw("same");

	c->SaveAs(Form("Plots/%s/CrossONLYLog_%s.pdf", B_m.Data(), var_n.Data()));


// CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection 
// CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection 

	
// 	COMPARISON OF 1D vs 2D methods
	if (meson_n == 0) {lat->DrawLatex(0.65,0.6 ,Form("2017 pp global Unc. #pm %.1f%%",3.5));} 
	else {	lat->DrawLatex(0.65,0.52,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;}

	BP1DRAAGraph->SetMarkerSize(0.8);
	BP1DRAAGraph_low ->SetMarkerSize(0.8);
	BP1DRAAGraph_low_just_marker ->SetMarkerSize(0.7);
	BPRAAGraph->SetMarkerSize(0.8);
	BPRAAGraph_low ->SetMarkerSize(0.8);
	BPRAAGraph_low_just_marker ->SetMarkerSize(0.7);

	HisEmpty->Draw();

	BP1DRAAGraphSyst_low->Draw("5same");
	BP1DRAAGraphSyst->Draw("5same");
	BP1DRAAGraph->Draw("epSAME");
	BP1DRAAGraph_low->Draw("epSAME");
	BP1DRAAGraph_low_just_marker->Draw("epSAME");

	BPRAAGraphSyst_low->Draw("5same");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph_low->Draw("epSAME");
	BPRAAGraph_low_just_marker->Draw("epSAME");

	if (meson_n==0) {leged = new TLegend(0.65,0.67,0.9,0.75,NULL,"brNDC");}
	else {leged = new TLegend(0.65,0.6,0.9,0.68,NULL,"brNDC");}
	leged->SetBorderSize(0);
	leged->SetFillStyle(0);
	if (whichvar==0){
		leged->AddEntry(BPRAAGraph,"2D","P");
		leged->AddEntry(BPRAAGraph_low,"2D (|y|>1.5)","P");
		leged->AddEntry(BP1DRAAGraph,"1D","P");
		leged->AddEntry(BP1DRAAGraph_low,"1D (|y|>1.5)","P");
	} 
	if (whichvar==1){
		leged->AddEntry(BPRAAGraph,"2D","P");
		leged->AddEntry(BP1DRAAGraph,"1D","P");
		leged->AddEntry(BPRAAGraph_low,"2D (p_{T}>10 GeV/c)","P");
		leged->AddEntry(BP1DRAAGraph_low,"1D (p_{T}>10 GeV/c)","P");
		} 
	leged->SetTextSize(0.020);
	leged->Draw("same");

	c->SaveAs(Form("Plots/%s/Cross1D2Dcomp_%s.pdf", B_m.Data(), var_n.Data()));

// 	COMPARISON OF 1D vs 2D methods


	string name;
	TString whichvarname;
	if(whichvar==0){name="$<p_T<$"; whichvarname="pt";} 
	else if(whichvar==1){name="$<y<$"; whichvarname="y";} 
	else if(whichvar==2){name="$<nTrks<$"; whichvarname="nMult";}
	
	std::vector<std::string> col_name;
	std::vector<std::string> col_name_diff;
	std::vector<std::string> labels = {"","Xsection 1D","Stat error","Syst error", "Xsection 2D","Stat error","Syst error"};
	std::vector<std::string> labels_diff = {"Relative difference 1D vs 2D"};
	std::vector<std::vector<double>> xsec_values;
	std::vector<std::vector<double>> xsec_diff;
	std::vector<double> val1D;
	std::vector<double> stat1D;
	std::vector<double> syst1D;
	std::vector<double> val2D;
	std::vector<double> stat2D;
	std::vector<double> syst2D;
	std::vector<double> xdiff;
	col_name_diff.push_back("");
	
	for(int i=0;i<NBins;i++){
		std::ostringstream clabel;
		clabel<<ptBins[i]<<name<<ptBins[i+1];
		std::string label1 = clabel.str();
		col_name.push_back(label1);
		col_name_diff.push_back(label1);

		val1D.push_back(BPXsecPPY1D[i]);
		stat1D.push_back(BPXSecPPY1DErrDown[i]/BPXsecPPY1D[i]);
		syst1D.push_back(BPXSecPPY1DSystDown[i]/BPXsecPPY1D[i]);
		val2D.push_back(BPXsecPPY2D[i]);
		stat2D.push_back(BPXSecPPY2DErrUp[i]/BPXsecPPY2D[i]);
		syst2D.push_back(BPXSecPPY2DSystDown[i]/BPXsecPPY2D[i]);

		xdiff.push_back(abs(BPXsecPPY2D[i]-BPXsecPPY1D[i])/BPXsecPPY2D[i]*100);
	}

	xsec_values.push_back(val1D);
	xsec_values.push_back(stat1D);
	xsec_values.push_back(syst1D);
	xsec_values.push_back(val2D);
	xsec_values.push_back(stat2D);
	xsec_values.push_back(syst2D);
	xsec_diff.push_back(xdiff);

	gSystem->mkdir("Trash",true);

	latex_table(Form("1D2DXseccomparisons_%s",whichvarname.Data()), 7,  NBins+1,  labels , col_name , xsec_values, "1D vs 2D Cross section comparisons",0);
	latex_table(Form("1D2DXsecdiff_%s",whichvarname.Data()), NBins+1,  2,  col_name_diff , labels_diff , xsec_diff, "1D vs 2D Cross section differences",1);
	std::vector<std::string> filetype ={"_check.aux", "_check.log",".tex","_check.tex"};
	for (int j=0;j<(int)(filetype.size());j++){
		rename(("1D2DXseccomparisons_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),("Trash/1D2DXseccomparisons_"+std::string (whichvarname.Data())+filetype[j]).c_str());
		rename(("1D2DXsecdiff_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),("Trash/1D2DXsecdiff_"+std::string (whichvarname.Data())+filetype[j]).c_str());
	}
	rename(("1D2DXseccomparisons_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),("Plots/"+std::string (B_m.Data())+"/1D2DXseccomparisons_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("1D2DXsecdiff_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),("Plots/"+std::string (B_m.Data())+"/1D2DXsecdiff_"+std::string (whichvarname.Data())+"_check.pdf").c_str());



//  XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb 
//2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 
		float BXsecPPX2015[NBins2015] ;
		float BXSecPPXErrDown2015[NBins2015] ;
		float BXSecPPXErrUp2015[NBins2015] ;
		float BXsecPPY2015[NBins2015] ;
		float BXSecPPYErrDown2015[NBins2015] ;
		float BXSecPPYErrUp2015[NBins2015] ;
		float BXSecPPYSystDown2015[NBins2015] ;
		float BXSecPPYSystUp2015[NBins2015] ;

		float sys_up_pbpb[4] ={0,0,0,0};
		float sys_down_pbpb[4] ={0,0,0,0};
		float sat_up_pbpb[4] ={0,0,0,0};
		float sat_down_pbpb[4] ={0,0,0,0};

if (whichvar==0){


		TGraphAsymmErrors *BPPbPbCrossGraph;
  		TGraphAsymmErrors *BPPbPbCrossGraphSyst;

		if (meson_n == 0){
			for(int i=0; i<4; i++) {
			sys_up_pbpb[i]  = BPXSecPbPbYSystUpRatio[i]*BPXsecPbPbY[i];
			sys_down_pbpb[i]= BPXSecPbPbYSystDownRatio[i]*BPXsecPbPbY[i];
			sat_up_pbpb[i]  = BPXSecPbPbYErrUpRatio[i]*BPXsecPbPbY[i];
			sat_down_pbpb[i]= BPXSecPbPbYErrDownRatio[i]*BPXsecPbPbY[i];
						}
			BPPbPbCrossGraph = new TGraphAsymmErrors(4, abscissae, BPXsecPbPbY,abscissae_x_y, abscissae_x_y,sat_down_pbpb,sat_up_pbpb);
			BPPbPbCrossGraphSyst  = new TGraphAsymmErrors(4, abscissae, BPXsecPbPbY, abscissae_x_y, abscissae_x_y,sys_down_pbpb ,sys_up_pbpb);

		} else {
			for(int i=0; i<4; i++) {
			sys_up_pbpb[i]  = BsXSecPbPbYSystUpPercent[i]*BsXsecPbPbY[i];
			sys_down_pbpb[i]= BsXSecPbPbYSystDownPercent[i]*BsXsecPbPbY[i];
			sat_up_pbpb[i]  = BsXSecPbPbYErrUpPercent[i]*BsXsecPbPbY[i];
			sat_down_pbpb[i]= BsXSecPbPbYErrDownPercent[i]*BsXsecPbPbY[i];
			}
			BPPbPbCrossGraph = new TGraphAsymmErrors(4, abscissae, BsXsecPbPbY,abscissae_x_y, abscissae_x_y,sat_down_pbpb,sat_up_pbpb);
			BPPbPbCrossGraphSyst  = new TGraphAsymmErrors(4, abscissae, BsXsecPbPbY, abscissae_x_y, abscissae_x_y,sys_down_pbpb ,sys_up_pbpb);
		}

			BPPbPbCrossGraph->SetLineColor(kOrange+1);
			BPPbPbCrossGraph->SetMarkerColor(kOrange+1);
			BPPbPbCrossGraph->SetMarkerStyle(21);
			BPPbPbCrossGraph->SetMarkerSize(1);
			BPPbPbCrossGraphSyst->SetFillColorAlpha(kOrange+1,0.5);
			BPPbPbCrossGraphSyst->SetLineColor(kOrange+1);

			HisEmpty->Draw();
			BPPbPbCrossGraph->Draw("ep");	
			BPPbPbCrossGraphSyst->Draw("5same");	
			BPRAAGraphSyst_low->Draw("5same");
			BPRAAGraphSyst->Draw("5same");
			BPRAAGraph->Draw("epSAME");
			BPRAAGraph_low->Draw("epSAME");
			BPRAAGraph_low_just_marker->Draw("epSAME");

			lat->DrawLatex(0.15,0.91 , "CMS work in progress");
				if (meson_n == 0) {
					lat->DrawLatex(0.6,0.65 ,Form("2018 PbPb global Unc. #pm %.1f%%",3.9));
					lat->DrawLatex(0.6,0.7 ,Form("2017 pp global Unc. #pm %.1f%%",3.5)) ;
				} else {
					lat->DrawLatex(0.6,0.65 ,Form("2018 PbPb global Unc. #pm %.1f%%",7.9));
					lat->DrawLatex(0.6,0.7,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;
				}

			TLegend* leg1 = new TLegend(0.65,0.74,0.9,0.85,NULL,"brNDC");
			leg1->SetBorderSize(0);
			leg1->SetFillStyle(0);
			if(meson_n == 0) { leg1->AddEntry((TObject*)0, "B^{+}", "");}
			else {leg1->AddEntry((TObject*)0, "B^{0}_{s}", "");}
			leg1->AddEntry(BPRAAGraph,"2017 pp ","P");
			leg1->AddEntry(BPRAAGraph_low,"2017 pp (|y|>1.5)","P");
			leg1->AddEntry(BPPbPbCrossGraph,"2018 PbPb  ","P");
			leg1->SetTextSize(0.025);
			leg1->Draw("same");

			if (whichvar==0){c->SaveAs(Form("Plots/%s/PbPbPPCrossLog.pdf", B_m.Data()));}

//  XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb 
//2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 


		if(meson_n == 0) { 
			for( int c=0; c <NBins2015; c++){ 
				BXsecPPX2015[c]= vect_BPXsecPPX2015[c] ;
				BXSecPPXErrDown2015[c]= vect_BPXSecPPXErrDown2015[c];
				BXSecPPXErrUp2015[c]= vect_BPXSecPPXErrUp2015[c];
				BXsecPPY2015[c]= vect_BPXsecPPY2015[c];
				BXSecPPYErrDown2015[c]= vect_BPXSecPPYErrDown2015[c];
				BXSecPPYErrUp2015[c]= vect_BPXSecPPYErrUp2015[c];
				BXSecPPYSystDown2015[c]= vect_BPXSecPPYSystDown2015[c];
				BXSecPPYSystUp2015[c] = vect_BPXSecPPYSystUp2015[c];
				}
		} else {
			for( int c=0; c <NBins2015; c++){ 
				BXsecPPX2015[c]= vect_BsXsecPPX2015[c] ;
				BXSecPPXErrDown2015[c]= vect_BsXSecPPXErrDown2015[c];
				BXSecPPXErrUp2015[c]= vect_BsXSecPPXErrUp2015[c];
				BXsecPPY2015[c]= vect_BsXsecPPY2015[c];
				BXSecPPYErrDown2015[c]= vect_BsXSecPPYErrDown2015[c];
				BXSecPPYErrUp2015[c]= vect_BsXSecPPYErrUp2015[c];
				BXSecPPYSystDown2015[c]= vect_BsXSecPPYSystDown2015[c];
				BXSecPPYSystUp2015[c] = vect_BsXSecPPYSystUp2015[c];
				}			
			}

	TGraphAsymmErrors *BPPPCrossGraph2015 = new TGraphAsymmErrors(NBins2015, BXsecPPX2015, BXsecPPY2015,BXSecPPXErrDown2015, BXSecPPXErrUp2015,BXSecPPYErrDown2015,BXSecPPYErrUp2015);
	TGraphAsymmErrors *BPPPCrossGraph2015Syst = new TGraphAsymmErrors(NBins2015, BXsecPPX2015, BXsecPPY2015,BXSecPPXErrDown2015, BXSecPPXErrUp2015,BXSecPPYSystDown2015,BXSecPPYSystUp2015);
	BPPPCrossGraph2015->SetLineColor(kOrange+1);
	BPPPCrossGraph2015->SetMarkerColor(kOrange+1);
	BPPPCrossGraph2015->SetMarkerStyle(21);
	BPPPCrossGraph2015->SetMarkerSize(1);
	BPPPCrossGraph2015Syst->SetLineColor(kOrange+1);
	BPPPCrossGraph2015Syst->SetFillColorAlpha(kOrange+1, 0.5);

	HisEmpty2->Draw();
	BPPPCrossGraph2015Syst->Draw("5same");
	BPPPCrossGraph2015->Draw("epSAME");
	BPRAAGraphSyst_low->Draw("5same");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph_low->Draw("epSAME");
	BPRAAGraph_low_just_marker->Draw("epSAME");

	lat->DrawLatex(0.15,0.91 , "CMS work in progress");
				if (meson_n == 0) {
					lat->DrawLatex(0.62,0.65 ,Form("2015 pp global Unc. #pm %.1f%%",3.8));
					lat->DrawLatex(0.62,0.7 ,Form("2017 pp global Unc. #pm %.1f%%",3.5)) ;
				} else {
					lat->DrawLatex(0.62,0.65 ,Form("2015 pp global Unc. #pm %.1f%%",7.9));
					lat->DrawLatex(0.62,0.7,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;
				}
	
			TLegend* lege = new TLegend(0.65,0.74,0.9,0.85,NULL,"brNDC");
				lege->SetBorderSize(0);
				lege->SetFillStyle(0);
			if(meson_n == 0) { lege->AddEntry((TObject*)0, "B^{+}", "");}
			else {lege->AddEntry((TObject*)0, "B^{0}_{s}", "");}
			lege->AddEntry(BPRAAGraph,"2017 pp ","P");
			lege->AddEntry(BPRAAGraph_low,"2017 pp (|y|>1.5)","P");
			lege->AddEntry(BPPPCrossGraph2015,"2015 pp  ","P");
			lege->SetTextSize(0.025);
			lege->Draw("same");

			c->SaveAs(Form("Plots/%s/pp_2015_CrossLog.pdf", B_m.Data()));
}
//2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 
//  XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb 




// vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL 
if (whichvar==0){

//	gStyle->SetPadTickX(1);
//	gStyle->SetPadTickY(1);
//	gStyle->SetTitleColor(1, "XYZ");
//	gStyle->SetTitleFont(42, "XYZ");
//	gStyle->SetTitleSize(0.15, "XYZ");
 	gStyle->SetTickLength(0.04, "XYZ");
//  gStyle->SetNdivisions(305, "XYZ");
//	gStyle->SetHatchesLineWidth(5);
// 	gStyle->SetHatchesSpacing(0.05);

	TCanvas * cr = new TCanvas("cr","cr",600,600);
	cr->cd();    
	cr->SetLeftMargin(0.15);

	//for comparisons we need this pad
	TPad * MyPadr;
	MyPadr = new TPad("MyPadr","",0,0.25,1,0.95);
	MyPadr->SetBottomMargin(0.);
	MyPadr->SetLogy();
	MyPadr->Draw();

	//for ratios we need 2nd pad
	TPad * MyPadr2;
	MyPadr2 = new TPad("MyPadr2","",0,0.0,1,0.25);
	MyPadr2->SetBottomMargin(0.3);
	MyPadr2->SetTopMargin(0);
	
	MyPadr2->Draw();

	TH2D * HisEmpty3;
	if (meson_n == 0){
		HisEmpty3 = new TH2D("HisEmpty3","",100,5,60,100,0.5,1.5);
		HisEmpty3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	} else {
		HisEmpty3 = new TH2D("HisEmpty3","",100,7,50,100,0.5,1.5);
		HisEmpty3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	}
		HisEmpty3->GetYaxis()->SetTitle("Data/FONLL");
		HisEmpty3->GetXaxis()->CenterTitle();
		HisEmpty3->GetYaxis()->CenterTitle();
		HisEmpty3->GetYaxis()->SetNdivisions(305);
		HisEmpty3->GetYaxis()->SetTitleSize(0.1);
		//HisEmpty3->GetXaxis()->SetLabelSize(0.1);
		HisEmpty3->GetYaxis()->SetLabelSize(0.1);
		HisEmpty3->GetXaxis()->SetTitleSize(0.1);
		HisEmpty3->GetYaxis()->SetTitleOffset(0.4);
		HisEmpty3->GetXaxis()->SetTitleOffset(1.0);
		HisEmpty3->GetXaxis()->SetLabelSize(0.1);

		HisEmpty2->GetXaxis()->SetTitleSize(0.035);
	/*		 


	pull_plotMC->GetXaxis()->SetTickLength(0.16);

	frameMC->GetXaxis()->SetTitleOffset(1.2);
	frameMC->GetXaxis()->SetTitleSize(0.035);
	frameMC->GetXaxis()->SetTitleFont(42);
	frameMC->GetXaxis()->CenterTitle();
	frameMC->GetXaxis()->SetLabelSize(0.035); */

    TFile * finFONLL ;
	if(meson_n == 0){ finFONLL = new TFile("FONLLs/fonllOutput_pp_Bplus_5p03TeV_y2p4.root");}
	else{ finFONLL = new TFile("FONLLs/BsFONLL.root");}
	finFONLL->cd();
	TGraphAsymmErrors *BPFONLL = (TGraphAsymmErrors*) finFONLL->Get("gaeSigmaBplus");
	BPFONLL->SetLineWidth(2);
	BPFONLL->SetLineColor(kRed+2);
	BPFONLL->SetFillStyle(0);
	BPFONLL->SetMarkerStyle(20);
	BPFONLL->SetMarkerSize(1);
	BPFONLL->SetMarkerColor(kRed+2);

	TFile * finFONLL2 ;
	if(meson_n == 0){ finFONLL2 = new TFile("FONLLs/fonllOutput_pp_Bplus_5p03TeV_yFid.root");}
	else{ finFONLL2 = new TFile("FONLLs/BsFONLLFid.root");}
    finFONLL2->cd();
	TGraphAsymmErrors *BFONLL2 = (TGraphAsymmErrors*) finFONLL2->Get("gaeSigmaBplus");

	TGraphAsymmErrors *BFONLLLow= new TGraphAsymmErrors();
	BFONLLLow->SetLineWidth(2);
	BFONLLLow->SetLineColor(kRed-7);
	BFONLLLow->SetFillStyle(0);
	BFONLLLow->SetMarkerStyle(20);
	BFONLLLow->SetMarkerSize(1);
	BFONLLLow->SetMarkerColor(kRed-7);


double XTempChange;
double YTempChange;
double YErrLowTemp;
double YErrHighTemp;

	for(int i = 0; i < NBinsLow; i ++){                      // STILL NEED TO CHANGE FOR OTHER VARIABLES
		BFONLL2->GetPoint(i,XTempChange,YTempChange);
		YErrLowTemp = BFONLL2->GetErrorYlow(i);
		YErrHighTemp = BFONLL2->GetErrorYhigh(i);
		BPFONLL->SetPoint(i,XTempChange,YTempChange);
		BPFONLL->SetPointEYhigh(i,YErrHighTemp);
		BPFONLL->SetPointEYlow(i,YErrLowTemp);
		BFONLLLow->AddPoint(XTempChange,YTempChange);
		BFONLLLow->SetPointEYhigh(i,YErrHighTemp);
		BFONLLLow->SetPointEYlow(i,YErrLowTemp);
		BFONLLLow->SetPointEXhigh(i,BPFONLL->GetErrorXhigh(i));
		BFONLLLow->SetPointEXlow(i,BPFONLL->GetErrorXlow(i));
}
	//BPPPCrossGraph2015Syst->Draw("5same");	
	//BPPPCrossGraph2015->Draw("epSAME");
	MyPadr->cd();
	HisEmpty2->Draw();
	BPRAAGraphSyst_low->Draw("5same");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph_low->Draw("epSAME");
	BPRAAGraph_low_just_marker->Draw("epSAME");
	BPFONLL->Draw("5");
	BFONLLLow->Draw("5");


	lat->SetTextSize(0.035); 
	lat->DrawLatex(0.1,0.91 , "CMS work in progress");
    if (meson_n == 0) {lat->DrawLatex(0.57,0.62 ,Form("2017 pp global Unc. #pm %.1f%%",3.5)) ;}
	else {lat->DrawLatex(0.6,0.62,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;}

			TLegend* leg3 = new TLegend(0.6,0.68,0.9,0.85,NULL,"brNDC");
				leg3->SetBorderSize(0);
				leg3->SetFillStyle(0);
			if(meson_n == 0) { leg3->AddEntry((TObject*)0, "B^{+}", "");}
			else {leg3->AddEntry((TObject*)0, "B^{0}_{s}", "");}
			leg3->AddEntry(BPRAAGraph,"2017 pp ","P");
			leg3->AddEntry(BPRAAGraph_low,"2017 pp (|y|>1.5)","P");
			//leg3->AddEntry(BPPPCrossGraph2015,"2018 PbPb","P");
			leg3->AddEntry(BPFONLL,"FONLL","f");
			leg3->AddEntry(BFONLLLow,"FONLL (|y| > 1.5)","f");
			leg3->Draw("same");
			MyPadr->Update();


	//Ratio
	float Ratio1Y[NBins2015];
	float Ratio1YErr[NBins2015];
	float Ratio2Y[NBins2015];
	float Ratio2YErr[NBins2015];

	for(int i = 1; i < NBins2015; i++){
		Ratio2Y[i] = BPXsecPPY2D[i+1]/BXsecPPY2015[i];
		Ratio2YErr[i] = Ratio2Y[i] *TMath::Sqrt(BPXSecPPY2DErrDown[i+1]/BPXsecPPY2D[i+1] * BPXSecPPY2DErrDown[i+1]/BPXsecPPY2D[i+1] + BXSecPPYErrDown2015[i]/BXsecPPY2015[i] * BXSecPPYErrDown2015[i]/BXsecPPY2015[i] );
			}

	//These vectors are just for BP
  std::vector<float> RatioDataYLow(1);
  std::vector<float> RatioDataYLowErr(1);
  	//These vectors are just for BP

if (meson_n == 0){
	Ratio2YErr[0] = 0.00001;
	Ratio2Y[0] = -0.1;
	Ratio1YErr[0] = 0.00001;
	Ratio1Y[0] = -0.1;
	// low pt bin
  RatioDataYLow[0] = BPXsecPPYLow[1] / BXsecPPY2015[0];
  RatioDataYLowErr[0] = RatioDataYLow[0] * TMath::Sqrt(pow(BPXSecPPY2DErrDown[1]/BPXsecPPY2D[1], 2) + pow(BXSecPPYErrDown2015[0]/BXsecPPY2015[0], 2));
}else{
	Ratio1Y[0] = -1;
	Ratio1YErr[0] = 0.001;
	Ratio2Y[0] = -1;
	Ratio2YErr[0] = 0.001;
}

	TGraphAsymmErrors *Ratio2 = new TGraphAsymmErrors(NBins2015, BXsecPPX2015, Ratio2Y,BXSecPPXErrDown2015, BXSecPPXErrUp2015,Ratio2YErr,Ratio2YErr);
		Ratio2->SetLineColor(kOrange+1);
		Ratio2->SetMarkerStyle(21);
		Ratio2->SetMarkerSize(1);
		Ratio2->SetMarkerColor(kOrange+1);


  	TGraphAsymmErrors *RatioDataLow;
  	if (meson_n ==0){
		RatioDataLow = new TGraphAsymmErrors(1, BPXsecPPXLow.data(), RatioDataYLow.data(), BPXsecPPXErrLow.data(), BPXsecPPXErrLow.data(), RatioDataYLowErr.data(), RatioDataYLowErr.data());
		RatioDataLow->SetLineColor(kOrange+1);
		RatioDataLow->SetMarkerStyle(25);
		RatioDataLow->SetMarkerSize(1);
		RatioDataLow->SetMarkerColor(kOrange+1);
		}

	MyPadr2->cd();
	TLine * Unity2 = new TLine(0,0,0,0);
	if (meson_n == 0){Unity2 = new TLine(5,1,60,1);}
	else {Unity2 = new TLine(7,1,50,1);}
	Unity2->SetLineWidth(2);
	Unity2->SetLineStyle(2);
	Unity2->SetLineColor(1);
	Unity2->Draw("SAME");
	HisEmpty3->Draw();
	MyPadr2->Update();

	//FONLL
	double XTempFONLL;
	double YTempFONLL;
	float Ratio4Y[NBins];
	float Ratio4YErr[NBins];
	float FONLLY[NBins];
	float FONLLYErr[NBins];
	for(int i = 0; i < NBins; i++){
		BPFONLL->GetPoint(i,XTempFONLL,YTempFONLL);
		FONLLY[i] = YTempFONLL;
		FONLLYErr[i] = BPFONLL->GetErrorYhigh (i);
		Ratio4Y[i] = BPXsecPPY2DScaled[i]/FONLLY[i];
		Ratio4YErr[i] = Ratio4Y[i] *TMath::Sqrt(BPXSecPPY2DErrDown[i]/BPXsecPPY2D[i] * BPXSecPPY2DErrDown[i]/BPXsecPPY2D[i] + FONLLYErr[i]/FONLLY[i] * FONLLYErr[i]/FONLLY[i] );
		}

 // Get ratio plots

float binlow[NBinsLow];
float glbSystUp;
float glbSystDown;
float bl_low[NBinsLow];
float bl_low_yStatL[NBinsLow];
float bl_low_yStatH[NBinsLow];
float bl_low_xErrL[NBinsLow];
float bl_low_xErrH[NBinsLow];
float bl_low_ySystL[NBinsLow];
float bl_low_ySystH[NBinsLow];
float binhigh[NBins-NBinsLow];
float bl_high[NBins-NBinsLow];
float bl_high_yStatL[NBins-NBinsLow];
float bl_high_yStatH[NBins-NBinsLow];
float bl_high_xErrL[NBins-NBinsLow];
float bl_high_xErrH[NBins-NBinsLow];
float bl_high_ySystL[NBins-NBinsLow];
float bl_high_ySystH[NBins-NBinsLow];
int NBinsLow2015;
if (meson_n==0){NBinsLow2015=1;} else {NBinsLow2015=0;}
float binlow_2015[NBinsLow2015];
float bl_low_2015[NBinsLow2015];
float bl_low_2015_yStatL[NBinsLow2015];
float bl_low_2015_yStatH[NBinsLow2015];
float bl_low_2015_xErrL[NBinsLow2015];
float bl_low_2015_xErrH[NBinsLow2015];
float bl_low_2015_ySystL[NBinsLow2015];
float bl_low_2015_ySystH[NBinsLow2015];
float binhigh_2015[NBins2015-NBinsLow2015];
float bl_high_2015[NBins2015-NBinsLow2015];
float bl_high_2015_yStatL[NBins2015-NBinsLow2015];
float bl_high_2015_yStatH[NBins2015-NBinsLow2015];
float bl_high_2015_xErrL[NBins2015-NBinsLow2015];
float bl_high_2015_xErrH[NBins2015-NBinsLow2015];
float bl_high_2015_ySystL[NBins2015-NBinsLow2015];
float bl_high_2015_ySystH[NBins2015-NBinsLow2015];


for (int i=0;i<NBins;++i){
	
	if( i<NBinsLow){
		binlow[i]=BPXsecPPX[i];
		glbSystUp=globUncert[i]*100;
		glbSystDown=globUncert[i]*100;
		bl_low[i]=BPXsecPPY2DScaled[i];
		bl_low_yStatL[i]=BPXSecPPY2DErrDownScaled[i];
		bl_low_yStatH[i]=BPXSecPPY2DErrUpScaled[i];
		bl_low_xErrL[i]=BPXsecPPX[i]-ptBins[i];
		bl_low_xErrH[i]=ptBins[i + 1]-BPXsecPPX[i];
		bl_low_ySystL[i]=BP2DTotalSystDownRatio[i]*BPXsecPPY2DScaled[i];
		bl_low_ySystH[i]=BP2DTotalSystUpRatio[i]*BPXsecPPY2DScaled[i];
	} 
	else {
		binhigh[i-NBinsLow]=BPXsecPPX[i];
		bl_high[i-NBinsLow]=BPXsecPPY2DScaled[i];
		bl_high_yStatL[i-NBinsLow]=BPXSecPPY2DErrDownScaled[i];
		bl_high_yStatH[i-NBinsLow]=BPXSecPPY2DErrUpScaled[i];
		bl_high_xErrL[i-NBinsLow]=BPXsecPPX[i]-ptBins[i];
		bl_high_xErrH[i-NBinsLow]=ptBins[i + 1]-BPXsecPPX[i];
		bl_high_ySystL[i-NBinsLow]=BP2DTotalSystDownRatio[i]*BPXsecPPY2DScaled[i];
		bl_high_ySystH[i-NBinsLow]=BP2DTotalSystUpRatio[i]*BPXsecPPY2DScaled[i];
	}
}
double ptBins2015[NBins2015+1];
if (meson_n==0){
	for(auto i=0;i<NBins2015+1;++i){
		ptBins2015[i]=vect_ptBins2015bp[i];
	}
}
else {for(auto i=0;i<NBins2015+1;++i){
		ptBins2015[i]=vect_ptBins2015bs[i];
	}}
for (int i=0;i<NBins2015;++i){
	
	if( i<NBinsLow2015){
		binlow_2015[i]=BXsecPPX2015[i];
		bl_low_2015[i]=BXsecPPY2015[i];
		bl_low_2015_yStatL[i]=BXSecPPYErrDown2015[i];
		bl_low_2015_yStatH[i]=BXSecPPYErrUp2015[i];
		bl_low_2015_xErrL[i]=BXsecPPX2015[i]-ptBins2015[i];
		bl_low_2015_xErrH[i]=ptBins2015[i + 1]-BXsecPPX2015[i];
		bl_low_2015_ySystL[i]=BXSecPPYSystDown2015[i];
		bl_low_2015_ySystH[i]=BXSecPPYSystDown2015[i];	
	} 
	else {
		binhigh_2015[i-NBinsLow2015]=BXsecPPX2015[i];
		bl_high_2015[i-NBinsLow2015]=BXsecPPY2015[i];
		bl_high_2015_yStatL[i-NBinsLow2015]=BXSecPPYErrDown2015[i];
		bl_high_2015_yStatH[i-NBinsLow2015]=BXSecPPYErrUp2015[i];
		bl_high_2015_xErrL[i-NBinsLow2015]=BXsecPPX2015[i]-ptBins2015[i];
		bl_high_2015_xErrH[i-NBinsLow2015]=ptBins2015[i + 1]-BXsecPPX2015[i];
		bl_high_2015_ySystL[i-NBinsLow2015]=BXSecPPYSystDown2015[i];
		bl_high_2015_ySystH[i-NBinsLow2015]=BXSecPPYSystUp2015[i];
	}
}
  vector<float> BXsec;
  vector<float> BXsecStat;
  vector<float> BXsecSyst;
  vector<float> BXsec2015;
  vector<float> BXsecStat2015;
  vector<float> BXsecSyst2015;
  vector<float> FONLL;
  vector<float> FONLLUp;
  vector<float> FONLLDown;

 
  for (auto i = 0; i < NBinsLow; ++i) {
    BXsec.push_back(bl_low[i]);
    BXsecStat.push_back(bl_low_yStatL[i]);
    BXsecSyst.push_back(bl_low_ySystL[i]);
  }
  for (auto i = 0; i < NBinsHigh; ++i) {
    BXsec.push_back(bl_high[i]);
    BXsecStat.push_back(bl_high_yStatL[i]);
    BXsecSyst.push_back(bl_high_ySystL[i]);
  }
  for (auto i = 0; i < NBinsLow2015; ++i) {
    BXsec2015.push_back(bl_low_2015[i]);
    BXsecStat2015.push_back(bl_low_2015_yStatL[i]);
    BXsecSyst2015.push_back(bl_low_2015_ySystL[i]);
  }
  for (auto i = 0; i < NBins2015-NBinsLow2015; ++i) {
    BXsec2015.push_back(bl_high_2015[i]);
    BXsecStat2015.push_back(bl_high_2015_yStatL[i]);
    BXsecSyst2015.push_back(bl_high_2015_ySystL[i]);
  }

  float RatioBs[NBins];
  float RatioBsStat[NBins];
  float RatioBsSyst[NBins];
  float RatioBsFonErrHigh[NBins];
  float RatioBsFonErrLow[NBins];
  float RatioBs2015[NBins2015];
  float RatioBsStat2015[NBins2015];
  float RatioBsSyst2015[NBins2015];

std::vector<float> Unity(NBins, 1);

for (auto i = 0; i < NBins; ++i) {
  
    BPFONLL->GetPoint(i, XTempFONLL, YTempFONLL);
    FONLL.push_back(YTempFONLL);
    FONLLUp.push_back(BPFONLL->GetErrorYhigh(i));
    FONLLDown.push_back(BPFONLL->GetErrorYlow(i));

    RatioBs[i] = BXsec[i] / FONLL[i];
    RatioBsStat[i] = BXsecStat[i] / FONLL[i];
    RatioBsSyst[i] = BXsecSyst[i] / FONLL[i];
    RatioBsFonErrHigh[i] = FONLLUp[i] / FONLL[i];
    RatioBsFonErrLow[i] = FONLLDown[i] / FONLL[i];
  }

int start2015;
if(meson_n==0){start2015=0;} else {start2015=1;}

for (auto i=start2015;i<NBins2015;i++){
  	
	RatioBs2015[i] = BXsec2015[i] / BXsec[i+1];
    RatioBsStat2015[i] = BXsecStat2015[i] / BXsec[i+1];
    RatioBsSyst2015[i] = BXsecSyst2015[i] / BXsec[i+1];
}
  if(meson_n==0){for (int i=0; i<NBinsLow; ++i){BPFONLL->RemovePoint(0);}}

  //DATA
  TGraphAsymmErrors gRatioBs_low(NBinsLow, binlow, RatioBs, bl_low_xErrL, bl_low_xErrH, RatioBsStat, RatioBsStat);

  TGraphAsymmErrors *BPRAAGraph_low_just_m = new TGraphAsymmErrors(NBinsLow, binlow, RatioBs ,zero, zero, zero, zero);

  TGraphAsymmErrors gRatioBs_high(NBinsHigh,binhigh,RatioBs + NBinsLow,bl_high_xErrL, bl_high_xErrH,RatioBsStat + NBinsLow, RatioBsStat + NBinsLow);
  TGraphAsymmErrors gRatioBs_syst_low(NBinsLow,binlow,RatioBs,bl_low_xErrL, bl_low_xErrH,RatioBsSyst, RatioBsSyst);
  TGraphAsymmErrors gRatioBs_syst_high(NBinsHigh,binhigh,RatioBs + NBinsLow,bl_high_xErrL, bl_high_xErrH,RatioBsSyst + NBinsLow, RatioBsSyst + NBinsLow);
  
  //FONLL
  TGraphAsymmErrors gRatioBs_Fon_low(NBinsLow,binlow,Unity.data(),bl_low_xErrL, bl_low_xErrH,RatioBsFonErrLow, RatioBsFonErrHigh);
  TGraphAsymmErrors gRatioBs_Fon_high(NBinsHigh,binhigh,Unity.data(),bl_high_xErrL, bl_high_xErrH,RatioBsFonErrLow + NBinsLow,RatioBsFonErrHigh + NBinsLow);

  //2015 currently not being ploted
  TGraphAsymmErrors gRatioBs2015_high(NBins2015-NBinsLow2015,binhigh_2015,RatioBs2015 + NBinsLow2015,bl_high_2015_xErrL, bl_high_2015_xErrH,RatioBsStat2015 + NBinsLow2015, RatioBsStat2015 + NBinsLow2015);
  TGraphAsymmErrors gRatioBs2015_syst_high(NBins2015-NBinsLow2015,binhigh_2015,RatioBs2015 + NBinsLow2015,bl_high_2015_xErrL, bl_high_2015_xErrH,RatioBsSyst2015 + NBinsLow2015, RatioBsSyst2015 + NBinsLow2015);
  gRatioBs2015_high.SetMarkerColor(kOrange+2);
  gRatioBs2015_high.SetLineColor(kOrange+2);
  gRatioBs2015_high.SetMarkerStyle(20);
  gRatioBs2015_syst_high.SetFillColorAlpha(kOrange+1, 0.5);

int color_mark =  kBlue + 2;
int color_syst = kBlue -3;
if(meson_n==0){	
color_mark = kGreen +2;
color_syst = kGreen -7;
}

BPRAAGraph_low_just_m ->SetMarkerStyle(20);
BPRAAGraph_low_just_m ->SetMarkerSize(0.9);
BPRAAGraph_low_just_m ->SetMarkerColor(kWhite);

  gRatioBs_syst_low.SetFillColorAlpha(color_syst, 0.5);
  gRatioBs_syst_low.SetLineColor(color_syst);
  gRatioBs_low.SetMarkerStyle(24);
  gRatioBs_low.SetMarkerColor(color_mark);
  gRatioBs_low.SetLineColor(color_mark);
  gRatioBs_syst_high.SetFillColorAlpha(color_syst, 0.5);
  gRatioBs_syst_high.SetLineColor(color_syst);
  gRatioBs_high.SetMarkerStyle(20);
  gRatioBs_high.SetMarkerColor(color_mark);
  gRatioBs_high.SetLineColor(color_mark);


  gRatioBs_Fon_low.SetLineColor(kRed-7); 
  gRatioBs_Fon_low.SetFillStyle(0);
  gRatioBs_Fon_low.SetLineWidth(2);
  gRatioBs_Fon_high.SetLineColor(kRed+2);
  gRatioBs_Fon_high.SetFillStyle(0);
  gRatioBs_Fon_high.SetLineWidth(2);


  gRatioBs_syst_low.Draw("5");
  gRatioBs_syst_high.Draw("5");
  BPRAAGraph_low_just_m->Draw("p");
  gRatioBs_low.Draw("ep");
  gRatioBs_high.Draw("ep");
  gRatioBs_Fon_high.Draw("5");
  gRatioBs_Fon_low.Draw("5");
  Unity2->Draw("SAME");

	MyPadr2->Update();
	//CMS_lumi(MyPadr,19011,0);
	MyPadr->Update();
  
  	cr->SetLogy();   
	if (whichvar==0){cr->SaveAs(Form("Plots/%s/CrossCompLog.pdf", B_m.Data()));}
	//FONLL
}




















  // summary of errors (in ratio, not percent)
  gSystem->mkdir("../../../MakeFinalPlots/NominalPlots/CrossSection/dataSource/" ,true );
  string outFile;
  if(meson_n == 0){ outFile = "../../../MakeFinalPlots/NominalPlots/CrossSection/dataSource/corryield_pt_Bp_New.txt";}
  else {outFile = "../../../MakeFinalPlots/NominalPlots/CrossSection/dataSource/corryield_pt_Bs_New.txt";}
  ofstream out;
  out.open(outFile);
  unsigned columnWidth = 14;
  out << std::left << std::setw(columnWidth) <<
    "ptmin" << std::setw(columnWidth) << "ptmax" << std::setw(columnWidth) <<
    "central val" << std::setw(columnWidth) <<
    "statUp" << std::setw(columnWidth) << "statDown" << std::setw(columnWidth) <<
    "systUp" << std::setw(columnWidth) << "systDown" << std::setw(columnWidth) <<
    "glbUp" << std::setw(columnWidth) << "glbDown" << std::setw(columnWidth) <<
    "abscissae" << endl;
  for (auto i = 0; i < NBins; ++i ) {
    out << std::setw(columnWidth) <<
      ptBins[i] << std::setw(columnWidth) << ptBins[i + 1] << std::setw(columnWidth) <<
      setprecision(0) << std::fixed << BPXsecPPY2D[i] << std::setw(columnWidth) <<
      setprecision(2) << std::defaultfloat <<
      BPXSecPPY2DErrUpRatio[i] << std::setw(columnWidth) <<
      BPXSecPPY2DErrDownRatio[i] << std::setw(columnWidth) <<
      BP2DTotalSystDownRatio[i] << std::setw(columnWidth) <<
      BP2DTotalSystDownRatio[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      setprecision(3) << BPXsecPPX[i] << "\n";
  }
  out.close();

}


						