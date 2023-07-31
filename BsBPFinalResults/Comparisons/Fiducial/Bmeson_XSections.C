#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TSystem.h"
#include <fstream>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include "../../../henri2022/parameter.h" 
#include "CMS_lumi.C" 

using namespace std;
using std::cout;
using std::endl;

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













































void Bmeson_XSections(TString meson_n, TString whichvar, int BsBPBins = 0){
                          
	int NBins = 7;
	int NBinsLow = 0  ;
  	int NBinsHigh = 0 ;
	int NBins2015 = 0 ;
	int lowend    ;
	TString var_l;
	TString Path_to_bin_Center = "../../../henri2022/ROOTfiles/" ;
	TString bsbpbins = "";

	if(whichvar == "pt"){
	
		if (meson_n=="BP" && BsBPBins==0){
			NBins = nptBinsBP;
			lowend = 1;
			NBinsLow = 2 ;
			NBinsHigh = 5;
			NBins2015 = 5;
			Path_to_bin_Center = Path_to_bin_Center + "yields_Bp_binned_pt.root" ;
		} 
		else if (BsBPBins==1 || meson_n == "Bs"){
			NBins = nptBins;
			lowend = 0;
			NBinsLow = 1 ;
			NBinsHigh = 3;
			NBins2015 = 3;
			
			if (meson_n == "BP" ){
				bsbpbins="_BsBPBINS";
				Path_to_bin_Center = Path_to_bin_Center + "yields_Bp_MatchingBINS_pt_BsBPBINS.root" ;
			} else {
				Path_to_bin_Center = Path_to_bin_Center + "yields_Bs_binned_pt.root" ;
			}
		}
		var_l="p_{T} [GeV/c]";
	
	} if(whichvar=="y"){
		NBins = nyBins_both;
		var_l="|y|";
		lowend = 2;
		NBinsLow = 3;
		NBinsHigh = 2;
		if(meson_n == "BP"){Path_to_bin_Center = Path_to_bin_Center + "yields_Bp_binned_y.root" ;} 
		else if (meson_n == "Bs") {Path_to_bin_Center =  Path_to_bin_Center + "yields_Bs_binned_y.root";} 
	
	} if(whichvar=="Mult"){
		NBins = nmBins_both;
		var_l="Mult";
		lowend = 100;
		NBinsLow = nmBins_both;         
		if(meson_n == "BP") {Path_to_bin_Center = Path_to_bin_Center + "yields_Bp_binned_Mult.root" ;} 
		else if (meson_n == "Bs") {Path_to_bin_Center = Path_to_bin_Center + "yields_Bs_binned_Mult.root" ;} 
	}

	gSystem->mkdir("Plots/", true);
	gSystem->mkdir(Form("Plots/%s",meson_n.Data()), true);
	
	TString InfileB = Form("../../../EffAna/%s/FinalFiles/%sPPCorrYield%s.root",meson_n.Data(),meson_n.Data(),whichvar.Data());
	TFile * FileB= new TFile(InfileB.Data());
	
	// BINS
	double ptBins[NBins+1];
	for(int i = 0; i < NBins + 1; i++){
		if (whichvar=="pt"){
			if (meson_n=="BP" && BsBPBins==0){ ptBins[i] =  ptbinsvecBP[i];} 
			else if (meson_n=="Bs" || BsBPBins=1){ptBins[i] =  ptbinsvec[i];}
		}
		else if (whichvar=="y"){ ptBins[i] =  ybinsvec[i]; }          
		else if (whichvar=="Mult"){ ptBins[i] =  nmbinsvec[i];}
	}

	//center of the bin and its left and right margins
	TFile* Dif_Plot = TFile::Open(Path_to_bin_Center, "READ");
	TMultiGraph* multiGraph = dynamic_cast<TMultiGraph*>(Dif_Plot->Get("TG"));
    TGraphAsymmErrors* Agraph = dynamic_cast<TGraphAsymmErrors*>(multiGraph->GetListOfGraphs()->At(0));
	float XsecPP_X[NBins];
	float XsecPP_X_BinRight[NBins] ;
	float XsecPP_X_BinLeft[NBins] ;
	for( int c=0; c < NBins; c++){
		XsecPP_X[c]= Agraph->GetX()[c];
		XsecPP_X_BinLeft[c] = XsecPP_X[c] - ptBins[c];
		XsecPP_X_BinRight[c]= ptBins[c+1] - XsecPP_X[c];
	}
    Dif_Plot->Close();
	//center of the bin and its left and right margins
	// BINS

  // cross section with 2D Map eff correction
	float XsecPP_Y[NBins];
	float XsecPP_Y_StatUp[NBins];
	float XsecPP_Y_StatDown[NBins];
	float BPXSecPPY2DErrUpRatio[NBins];
	float BPXSecPPY2DErrDownRatio[NBins];
  
	float BPXsecPPY2DScaled[NBins];
	float BPXSecPPY2DErrUpScaled[NBins];
	float BPXSecPPY2DErrDownScaled[NBins];

	TH1D * BCross2D = (TH1D *) FileB->Get("hPtSigma");
	for(int i = 0; i < NBins; i++){
		XsecPP_Y[i] = BCross2D->GetBinContent(i+1);
		XsecPP_Y_StatUp[i] = BCross2D->GetBinError(i+1);
		XsecPP_Y_StatDown[i] = BCross2D->GetBinError(i+1);
		BPXSecPPY2DErrUpRatio[i] = XsecPP_Y_StatUp[i] / XsecPP_Y[i];
		BPXSecPPY2DErrDownRatio[i] = XsecPP_Y_StatDown[i] / XsecPP_Y[i];
		
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
  	TString errorFile = Form("../../../2DMapSyst/OutFiles/%sError2D_%s%s.root", meson_n.Data(),whichvar.Data(),bsbpbins);
  	TFile fError(errorFile);

	TH1D * TnPSyst = (TH1D *) fError.Get("TnPSyst");
	TH1D * BptSyst = (TH1D *) fError.Get("BptSyst");
	TH1D * MCDataSyst = (TH1D *) fError.Get("MCDataSyst");
  	if (!MCDataSyst) MCDataSyst = (TH1D *) fError.Get("BDTSyst");

	TString errorFile1D = Form("../../../1DMapSyst/OutFiles/%sError1D_%s%s.root", meson_n.Data(),whichvar.Data(),bsbpbins);
  	TFile fError1D(errorFile1D);

	TH1D * TnPSyst1D = (TH1D *) fError1D.Get("TnPSyst");
	TH1D * BptSyst1D = (TH1D *) fError1D.Get("BptSyst");
	TH1D * MCDataSyst1D = (TH1D *) fError1D.Get("BDTSyst");

	TString pdfErrorFile = Form("../../../syst_error/%s_pdf_%s%s.root",meson_n.Data(),whichvar.Data(),bsbpbins);
	TFile fPdfError(pdfErrorFile);
	TGraph* pdfSyst = (TGraph *) fPdfError.Get(Form("%s_error",meson_n.Data()));
	
	TString trackSelErrorFile = Form("../../../syst_error/syst_track_sel_%s%s.root",whichvar.Data(),bsbpbins);
	TFile fTrackSelError(trackSelErrorFile);
	TGraph* trackSelSyst = (TGraph *) fTrackSelError.Get(Form("%s_track_sel_error%s", meson_n.Data()));

	TString trackSelErrorFile1D = Form("../../../syst_error/syst_track_sel_%s%s_1D.root",whichvar.Data(), bsbpbins);
	TFile fTrackSelError1D(trackSelErrorFile1D);
	TGraph* trackSelSyst1D = (TGraph *) fTrackSelError1D.Get(Form("%s_track_sel_error", meson_n.Data()));

	float XsecPP_Y_SystUp[NBins];
	float XsecPP_Y_SystDown[NBins];
	float BPXSecPPY1DSystUp[NBins];
	float BPXSecPPY1DSystDown[NBins];
	float BPXSecPPYSystUpScaled[NBins];
	float BPXSecPPYSystDownScaled[NBins];

  // percent error

  	double B_nu;
  	if (meson_n=="BP"){ B_nu = 2.4 ;}
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
	if(meson_n=="BP"){numb = 0.035;}
	else {numb = 0.077;}
  	vector<float> globUncert(NBins, numb);
	for(int i = 0; i < NBins; i++){
		XsecPP_Y_SystUp[i] = XsecPP_Y[i] * BP2DTotalSystUpRatio[i];
		XsecPP_Y_SystDown[i] = XsecPP_Y[i] * BP2DTotalSystDownRatio[i];
		BPXSecPPY1DSystUp[i] = BPXsecPPY1D[i] * BP1DTotalSystUpRatio[i];
		BPXSecPPY1DSystDown[i] = BPXsecPPY1D[i] * BP1DTotalSystDownRatio[i];
		}

// CREATE THE CANVAS and the pads
	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd(); 
	if (whichvar=="pt"){c->SetLogy();}
	c->SetLeftMargin(0.15);

	//Setup histograms for different purposs
	TH2D * HisEmpty;
	if(whichvar=="pt") {HisEmpty = new TH2D("HisEmpty","",100,ptBins[0],ptBins[NBins],100,300.0,3000000);} 

	if(meson_n=="BP" && whichvar=="y") {HisEmpty = new TH2D("HisEmpty","",100,ptBins[0],ptBins[NBins],100,1250000.0,10000000);}
	if(meson_n=="Bs" && whichvar=="y") {HisEmpty = new TH2D("HisEmpty","",100,ptBins[0],ptBins[NBins],100,220000.0,1250000);}

	if(meson_n=="BP" && whichvar=="Mult") {HisEmpty = new TH2D("HisEmpty","",100,ptBins[0],ptBins[NBins],100,0,4200000);}   // need to adjust range for when we have nmult results
	if(meson_n=="Bs" && whichvar=="Mult") {HisEmpty = new TH2D("HisEmpty","",100,ptBins[0],ptBins[NBins],100,0,600000);}

	HisEmpty->GetXaxis()->SetTitle(var_l.Data());
	if (whichvar=="pt") {HisEmpty->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb c/GeV]");}
	else {HisEmpty->GetYaxis()->SetTitle(Form("d#sigma/d%s [pb c/GeV]",whichvar.Data()));}
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
// CREATE THE CANVAS and the pads

  // separate plots for different fiducial regions
  	vector<float> XsecPP_X_Low ;
	vector<float> XsecPP_X_High ;
	vector<float> XsecPP_X_BinL_Low ;
	vector<float> XsecPP_X_BinR_Low ;
	vector<float> XsecPP_X_BinL_High ;
	vector<float> XsecPP_X_BinR_High ;

  	vector<float> XsecPP_Y_Low ;
  	vector<float> XsecPP_Y_High ;
	vector<float> XsecPP_Y_StatDown_Low ;
	vector<float> XsecPP_Y_StatDown_High ;
	vector<float> XsecPP_Y_StatUp_Low ;
	vector<float> XsecPP_Y_StatUp_High ;
	vector<float> XsecPP_Y_SystDown_Low ;
	vector<float> XsecPP_Y_SystDown_High ;
	vector<float> XsecPP_Y_SystUp_Low ;
	vector<float> XsecPP_Y_SystUp_High ;

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
		if ( i <= lowend){
			XsecPP_X_Low.push_back(XsecPP_X[i]);
			XsecPP_X_BinL_Low.push_back(XsecPP_X_BinLeft[i]);
			XsecPP_X_BinR_Low.push_back(XsecPP_X_BinRight[i]);

			XsecPP_Y_Low.push_back(XsecPP_Y[i]);
			XsecPP_Y_StatDown_Low.push_back(XsecPP_Y_StatDown[i]);
			XsecPP_Y_StatUp_Low.push_back(XsecPP_Y_StatUp[i]);
			XsecPP_Y_SystDown_Low.push_back(XsecPP_Y_SystDown[i]);
			XsecPP_Y_SystUp_Low.push_back(XsecPP_Y_SystUp[i]);

			BP1DXsecPPYLow.push_back(BPXsecPPY1D[i]);
			BP1DXsecPPYErrDownLow.push_back(BPXSecPPY1DErrDown[i]);
			BP1DXsecPPYErrUpLow.push_back(BPXSecPPY1DErrUp[i]);
			BP1DYSystDown_low.push_back(BPXSecPPY1DSystDown[i]);
			BP1DYSystUp_low.push_back(BPXSecPPY1DSystUp[i]);

		} else {
			XsecPP_X_High.push_back(XsecPP_X[i]);
			XsecPP_X_BinL_High.push_back(XsecPP_X_BinLeft[i]);
			XsecPP_X_BinR_High.push_back(XsecPP_X_BinRight[i]);

			XsecPP_Y_High.push_back(XsecPP_Y[i]);
			XsecPP_Y_StatDown_High.push_back(XsecPP_Y_StatDown[i]);
			XsecPP_Y_StatUp_High.push_back(XsecPP_Y_StatUp[i]);
			XsecPP_Y_SystDown_High.push_back(XsecPP_Y_SystDown[i]);
			XsecPP_Y_SystUp_High.push_back(XsecPP_Y_SystUp[i]);

			BP1DXsecPPYHigh.push_back(BPXsecPPY1D[i]);
			BP1DXsecPPYErrDownHigh.push_back(BPXSecPPY1DErrDown[i]);
			BP1DXsecPPYErrUpHigh.push_back(BPXSecPPY1DErrUp[i]);
			BP1DYSystDown_high.push_back(BPXSecPPY1DSystDown[i]);
			BP1DYSystUp_high.push_back(BPXSecPPY1DSystUp[i]);
		}
	}

	float zero[NBinsLow];
	for (int i=0 ; i<NBinsLow; i++){zero[i]=0;}

	TGraphAsymmErrors *BPRAAGraph_low_just_marker   = new TGraphAsymmErrors(NBinsLow, XsecPP_X_Low.data(), XsecPP_Y_Low.data() ,zero, zero, zero, zero);
	TGraphAsymmErrors *BP1DRAAGraph_low_just_marker = new TGraphAsymmErrors(NBinsLow, XsecPP_X_Low.data(), BP1DXsecPPYLow.data() ,zero, zero, zero, zero);

	TGraphAsymmErrors *BPRAAGraph_low = new TGraphAsymmErrors(NBinsLow , XsecPP_X_Low.data() , XsecPP_Y_Low.data() , XsecPP_X_BinL_Low.data() , XsecPP_X_BinR_Low.data() , XsecPP_Y_StatDown_Low.data() , XsecPP_Y_StatUp_Low.data());
	TGraphAsymmErrors *BPRAAGraph     = new TGraphAsymmErrors(NBinsHigh, XsecPP_X_High.data(), XsecPP_Y_High.data(), XsecPP_X_BinL_High.data(), XsecPP_X_BinR_High.data(), XsecPP_Y_StatDown_High.data(), XsecPP_Y_StatUp_High.data());     
	TGraphAsymmErrors *BPRAAGraphSyst_low  = new TGraphAsymmErrors(NBinsLow , XsecPP_X_Low.data() , XsecPP_Y_Low.data() , XsecPP_X_BinL_Low.data() , XsecPP_X_BinR_Low.data() , XsecPP_Y_SystDown_Low.data() , XsecPP_Y_SystUp_Low.data());                 											
	TGraphAsymmErrors *BPRAAGraphSyst      = new TGraphAsymmErrors(NBinsHigh, XsecPP_X_High.data(), XsecPP_Y_High.data(), XsecPP_X_BinL_High.data(), XsecPP_X_BinR_High.data(), XsecPP_Y_SystDown_High.data(), XsecPP_Y_SystUp_High.data());                 											
	
	TGraphAsymmErrors *BP1DRAAGraph_low = new TGraphAsymmErrors(NBinsLow , XsecPP_X_Low.data() , BP1DXsecPPYLow.data() , XsecPP_X_BinL_Low.data() , XsecPP_X_BinR_Low.data() , BP1DXsecPPYErrDownLow.data() , BP1DXsecPPYErrUpLow.data());
	TGraphAsymmErrors *BP1DRAAGraph     = new TGraphAsymmErrors(NBinsHigh, XsecPP_X_High.data(), BP1DXsecPPYHigh.data(), XsecPP_X_BinL_High.data(), XsecPP_X_BinR_High.data(), BP1DXsecPPYErrDownHigh.data(), BP1DXsecPPYErrUpHigh.data());     
	TGraphAsymmErrors *BP1DRAAGraphSyst_low  = new TGraphAsymmErrors(NBinsLow , XsecPP_X_Low.data() , BP1DXsecPPYLow.data() , XsecPP_X_BinL_Low.data() , XsecPP_X_BinR_Low.data() , BP1DYSystDown_low.data() , BP1DYSystUp_low.data());                 											
	TGraphAsymmErrors *BP1DRAAGraphSyst      = new TGraphAsymmErrors(NBinsHigh, XsecPP_X_High.data(), BP1DXsecPPYHigh.data(), XsecPP_X_BinL_High.data(), XsecPP_X_BinR_High.data(), BP1DYSystDown_high.data(), BP1DYSystUp_high.data());
  // separate plots for different fiducial regions
 
  	cout << endl << "-------------------------------------------------------  "<< Form("%s meson Xsection", meson_n.Data()) <<"  -------------------------------------------------------" << endl;

	for(int i=0;i<NBins;i++){		
		cout << "BIN " <<              Form("[%.1f,%.1f]  ",ptBins[i],ptBins[i+1]) << Form("%.0f #pm (STATup) %.0f #pm (SYSTup) %.0f #pm (STATdown) %.0f #pm (SYSTdown) %.0f ",XsecPP_Y[i], XsecPP_Y_StatUp[i], XsecPP_Y_SystUp[i], XsecPP_Y_StatDown[i], XsecPP_Y_SystDown[i]) << endl;
		cout << "(normalized) BIN " << Form("[%.1f,%.1f]  ",ptBins[i],ptBins[i+1]) << Form("%.0f #pm (STATup) %.1f #pm (SYSTup) %.1f #pm (STATdown) %.1f #pm (SYSTdown) %.1f ",XsecPP_Y[i], 100*XsecPP_Y_StatUp[i]/XsecPP_Y[i], 100*XsecPP_Y_SystUp[i]/XsecPP_Y[i], 100*XsecPP_Y_StatDown[i]/XsecPP_Y[i], 100*XsecPP_Y_SystDown[i]/XsecPP_Y[i]) << endl;
	}
 
 	cout<< endl << "-------------------------------------------------------  "<< Form("%s meson Xsection", meson_n.Data()) <<"  -------------------------------------------------------" << endl;

// CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection 
// CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection 
	
	BPRAAGraph->SetMarkerStyle(20);
	BPRAAGraph->SetMarkerSize(1);
	BPRAAGraph_low_just_marker ->SetMarkerStyle(20);
	BPRAAGraph_low_just_marker ->SetMarkerSize(0.9);
	BPRAAGraph_low_just_marker ->SetMarkerColor(kWhite);
	BPRAAGraph_low ->SetMarkerStyle(24);
	BPRAAGraph_low ->SetMarkerSize(1);

	if (meson_n=="BP"){
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
	if (meson_n=="BP") {lat->DrawLatex(0.6,0.7 ,Form("2017 pp global Unc. #pm %.1f%%",3.5));} 
	else {	lat->DrawLatex(0.6,0.7,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;}

    TLegend* leged;
	if (meson_n=="BP") {leged = new TLegend(0.65,0.6,0.9,0.7,NULL,"brNDC");}
	else {leged = new TLegend(0.65,0.6,0.9,0.7,NULL,"brNDC");}
	leged->SetBorderSize(0);
	leged->SetFillStyle(0);
	if (whichvar=="pt"){
		leged->AddEntry((TObject*)0, "y region:", "");
		leged->AddEntry(BPRAAGraph,"|y|<2.4","P");
		leged->AddEntry(BPRAAGraph_low,"|y|>1.5","P");
	} 
	if (whichvar=="y"){
		leged->AddEntry((TObject*)0, "p_{T} region:", "");
		if (meson_n=="BP"){leged->AddEntry(BPRAAGraph,"5<p_{T}<60 GeV/c","P");}
		else {leged->AddEntry(BPRAAGraph,"7<p_{T}<50 GeV/c","P");}
		leged->AddEntry(BPRAAGraph_low,"p_{T}>10 GeV/c","P");
	} 

	leged->SetTextSize(0.022);
	leged->Draw("same");

	c->SaveAs(Form("Plots/%s_Xsection_%s%s.pdf", meson_n.Data(), whichvar.Data()),bsbpbins);


// CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection 
// CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection 

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// 	COMPARISON OF 1D vs 2D methods
	if(BsBPBins==0){
	if (meson_n=="BP") {lat->DrawLatex(0.65,0.6 ,Form("2017 pp global Unc. #pm %.1f%%",3.5));} 
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

	if (meson_n=="BP") {leged = new TLegend(0.65,0.67,0.9,0.75,NULL,"brNDC");}
	else {leged = new TLegend(0.65,0.6,0.9,0.68,NULL,"brNDC");}
	leged->SetBorderSize(0);
	leged->SetFillStyle(0);
	if (whichvar=="pt"){
		leged->AddEntry(BPRAAGraph,"2D","P");
		leged->AddEntry(BPRAAGraph_low,"2D (|y|>1.5)","P");
		leged->AddEntry(BP1DRAAGraph,"1D","P");
		leged->AddEntry(BP1DRAAGraph_low,"1D (|y|>1.5)","P");
	} 
	if (whichvar=="y"){
		leged->AddEntry(BPRAAGraph,"2D","P");
		leged->AddEntry(BP1DRAAGraph,"1D","P");
		leged->AddEntry(BPRAAGraph_low,"2D (p_{T}>10 GeV/c)","P");
		leged->AddEntry(BP1DRAAGraph_low,"1D (p_{T}>10 GeV/c)","P");
		} 
	leged->SetTextSize(0.020);
	leged->Draw("same");

	c->SaveAs(Form("Plots/%s/Bmeson_1D2Dcomp_%s.pdf", meson_n.Data(), whichvar.Data()));
// 	COMPARISON OF 1D vs 2D methods

	string name;
	TString whichvarname;
	if(whichvar=="pt"){name="$<p_T<$"; whichvarname="pt";} 
	else if(whichvar=="y"){name="$<y<$"; whichvarname="y";} 
	else if(whichvar=="Mult"){name="$<nTrks<$"; whichvarname="nMult";}
	
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
		val2D.push_back(XsecPP_Y[i]);
		stat2D.push_back(XsecPP_Y_StatUp[i]/XsecPP_Y[i]);
		syst2D.push_back(XsecPP_Y_SystDown[i]/XsecPP_Y[i]);

		xdiff.push_back(abs(XsecPP_Y[i]-BPXsecPPY1D[i])/XsecPP_Y[i]*100);
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
	rename(("1D2DXseccomparisons_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),("Plots/"+std::string (meson_n.Data())+"/1D2DXseccomparisons_"+std::string (whichvarname.Data())+"_check.pdf").c_str());
	rename(("1D2DXsecdiff_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),("Plots/"+std::string (meson_n.Data())+"/1D2DXsecdiff_"+std::string (whichvarname.Data())+"_check.pdf").c_str());

	}
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

		float BXsecPPX2015[NBins2015] ;
		float BXSecPPXErrDown2015[NBins2015] ;
		float BXSecPPXErrUp2015[NBins2015] ;
		float BXsecPPY2015[NBins2015] ;
		float BXSecPPYErrDown2015[NBins2015] ;
		float BXSecPPYErrUp2015[NBins2015] ;
		float BXSecPPYSystDown2015[NBins2015] ;
		float BXSecPPYSystUp2015[NBins2015] ;


if (whichvar=="pt" && BsBPBins == 0){
//  XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb 

	TGraphAsymmErrors *BPPbPbCrossGraph;
  	TGraphAsymmErrors *BPPbPbCrossGraphSyst;

	// get the values for the histograms (check parameter.h)
	for(int i=0; i<4; i++) {

		XsecPbPb_XL_BP[i] = XsecPbPb_X_BP[i] - BINS_PbPb[i];
		XsecPbPb_XR_BP[i]= BINS_PbPb[i+1] - XsecPbPb_X_BP[i];
		XsecPbPb_XL_Bs[i] = XsecPbPb_X_Bs[i] - BINS_PbPb[i];
		XsecPbPb_XR_Bs[i]= BINS_PbPb[i+1] - XsecPbPb_X_Bs[i];

		XSecPbPb_BP_Y_StatUpRatio[i]  = XSecPbPb_BP_Y_StatUpRatio[i]*XsecPbPb_Y_BP[i];
		XSecPbPb_BP_Y_StatDownRatio[i]= XSecPbPb_BP_Y_StatDownRatio[i]*XsecPbPb_Y_BP[i];
		XSecPbPb_BP_Y_SystUpRatio[i]  = XSecPbPb_BP_Y_SystUpRatio[i]*XsecPbPb_Y_BP[i];
		XSecPbPb_BP_Y_SystDownRatio[i]= XSecPbPb_BP_Y_SystDownRatio[i]*XsecPbPb_Y_BP[i];
		XSecPbPb_Bs_Y_StatUpRatio[i]  = XSecPbPb_Bs_Y_StatUpRatio[i]*XsecPbPb_Y_Bs[i];
		XSecPbPb_Bs_Y_StatDownRatio[i]= XSecPbPb_Bs_Y_StatDownRatio[i]*XsecPbPb_Y_Bs[i];
		XSecPbPb_Bs_Y_SystUpRatio[i]  = XSecPbPb_Bs_Y_SystUpRatio[i]*XsecPbPb_Y_Bs[i];
		XSecPbPb_Bs_Y_SystDownRatio[i]= XSecPbPb_Bs_Y_SystDownRatio[i]*XsecPbPb_Y_Bs[i];
	}
	// get the value for the histograms (check parameter.h)

		if (meson_n=="BP"){
			BPPbPbCrossGraph      = new TGraphAsymmErrors(4, XsecPbPb_X_BP, XsecPbPb_Y_BP, XsecPbPb_XL_BP, XsecPbPb_XR_BP, XSecPbPb_BP_Y_StatDownRatio, XSecPbPb_BP_Y_StatUpRatio);
			BPPbPbCrossGraphSyst  = new TGraphAsymmErrors(4, XsecPbPb_X_BP, XsecPbPb_Y_BP, XsecPbPb_XL_BP, XsecPbPb_XR_BP, XSecPbPb_BP_Y_SystDownRatio, XSecPbPb_BP_Y_SystUpRatio);
		} else {
			BPPbPbCrossGraph      = new TGraphAsymmErrors(4, XsecPbPb_X_Bs, XsecPbPb_Y_Bs, XsecPbPb_XL_Bs, XsecPbPb_XR_Bs, XSecPbPb_Bs_Y_StatDownRatio, XSecPbPb_Bs_Y_StatUpRatio);
			BPPbPbCrossGraphSyst  = new TGraphAsymmErrors(4, XsecPbPb_X_Bs, XsecPbPb_Y_Bs, XsecPbPb_XL_Bs, XsecPbPb_XR_Bs, XSecPbPb_Bs_Y_SystDownRatio ,XSecPbPb_Bs_Y_SystUpRatio);
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
				if (meson_n=="BP") {
					lat->DrawLatex(0.6,0.65 ,Form("2018 PbPb global Unc. #pm %.1f%%",3.9));
					lat->DrawLatex(0.6,0.7 ,Form("2017 pp global Unc. #pm %.1f%%",3.5)) ;
				} else {
					lat->DrawLatex(0.6,0.65 ,Form("2018 PbPb global Unc. #pm %.1f%%",7.9));
					lat->DrawLatex(0.6,0.7,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;
				}

			TLegend* leg1 = new TLegend(0.65,0.74,0.9,0.85,NULL,"brNDC");
			leg1->SetBorderSize(0);
			leg1->SetFillStyle(0);
			if(meson_n=="BP") { leg1->AddEntry((TObject*)0, "B^{+}", "");}
			else {leg1->AddEntry((TObject*)0, "B^{0}_{s}", "");}
			leg1->AddEntry(BPRAAGraph,"2017 pp ","P");
			leg1->AddEntry(BPRAAGraph_low,"2017 pp (|y|>1.5)","P");
			leg1->AddEntry(BPPbPbCrossGraph,"2018 PbPb  ","P");
			leg1->SetTextSize(0.025);
			leg1->Draw("same");

			c->SaveAs(Form("Plots/%s_Xsection_%s_vsPbPb.pdf", meson_n.Data(),whichvar.Data()));

	//  XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb 
	//2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 

		if(meson_n=="BP") { 
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

	HisEmpty->Draw();
	BPPPCrossGraph2015Syst->Draw("5same");
	BPPPCrossGraph2015->Draw("epSAME");
	BPRAAGraphSyst_low->Draw("5same");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph_low->Draw("epSAME");
	BPRAAGraph_low_just_marker->Draw("epSAME");

	lat->DrawLatex(0.15,0.91 , "CMS work in progress");
				if (meson_n=="BP") {
					lat->DrawLatex(0.62,0.65 ,Form("2015 pp global Unc. #pm %.1f%%",3.8));
					lat->DrawLatex(0.62,0.7 ,Form("2017 pp global Unc. #pm %.1f%%",3.5)) ;
				} else {
					lat->DrawLatex(0.62,0.65 ,Form("2015 pp global Unc. #pm %.1f%%",7.9));
					lat->DrawLatex(0.62,0.7,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;
				}
	
			TLegend* lege = new TLegend(0.65,0.74,0.9,0.85,NULL,"brNDC");
				lege->SetBorderSize(0);
				lege->SetFillStyle(0);
			if(meson_n=="BP") { lege->AddEntry((TObject*)0, "B^{+}", "");}
			else {lege->AddEntry((TObject*)0, "B^{0}_{s}", "");}
			lege->AddEntry(BPRAAGraph,"2017 pp ","P");
			lege->AddEntry(BPRAAGraph_low,"2017 pp (|y|>1.5)","P");
			lege->AddEntry(BPPPCrossGraph2015,"2015 pp  ","P");
			lege->SetTextSize(0.025);
			lege->Draw("same");

			c->SaveAs(Form("Plots/%s_Xsection_%s_vs2015.pdf", meson_n.Data(),whichvar.Data()));

	//2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 
	//////////////////////////////////////////////////////////////////////////////////

	// vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL 
	// vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL 

 	gStyle->SetTickLength(0.04, "XYZ");
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
	if (meson_n=="BP"){
		HisEmpty3 = new TH2D("HisEmpty3","",100,ptBins[0],ptBins[NBins],100,0.5,1.5);
		HisEmpty3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	} else {
		HisEmpty3 = new TH2D("HisEmpty3","",100,ptBins[0],ptBins[NBins],100,0.5,1.5);
		HisEmpty3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	}
		HisEmpty3->GetYaxis()->SetTitle("Data/FONLL");
		HisEmpty3->GetXaxis()->CenterTitle();
		HisEmpty3->GetYaxis()->CenterTitle();
		HisEmpty3->GetYaxis()->SetNdivisions(305);
		HisEmpty3->GetYaxis()->SetTitleSize(0.1);
		HisEmpty3->GetYaxis()->SetLabelSize(0.1);
		HisEmpty3->GetXaxis()->SetTitleSize(0.1);
		HisEmpty3->GetYaxis()->SetTitleOffset(0.4);
		HisEmpty3->GetXaxis()->SetTitleOffset(1.0);
		HisEmpty3->GetXaxis()->SetLabelSize(0.1);

    TFile * finFONLL ;
	if(meson_n=="BP"){ finFONLL = new TFile("FONLLs/fonllOutput_pp_Bplus_5p03TeV_y2p4.root");}
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
	if(meson_n=="BP"){ finFONLL2 = new TFile("FONLLs/fonllOutput_pp_Bplus_5p03TeV_yFid.root");}
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

	for(int i = 0; i < NBinsLow; i ++){                   
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

	MyPadr->cd();
	HisEmpty->Draw();
	BPRAAGraphSyst_low->Draw("5same");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph_low->Draw("epSAME");
	BPRAAGraph_low_just_marker->Draw("epSAME");
	BPFONLL->Draw("5");
	BFONLLLow->Draw("5");

	lat->SetTextSize(0.035); 
	lat->DrawLatex(0.1,0.91 , "CMS work in progress");
    if (meson_n=="BP") {lat->DrawLatex(0.57,0.62 ,Form("2017 pp global Unc. #pm %.1f%%",3.5)) ;}
	else {lat->DrawLatex(0.6,0.62,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;}

			TLegend* leg3 = new TLegend(0.6,0.68,0.9,0.85,NULL,"brNDC");
				leg3->SetBorderSize(0);
				leg3->SetFillStyle(0);
			if(meson_n=="BP") { leg3->AddEntry((TObject*)0, "B^{+}", "");}
			else {leg3->AddEntry((TObject*)0, "B^{0}_{s}", "");}
			leg3->AddEntry(BPRAAGraph,"2017 pp ","P");
			leg3->AddEntry(BPRAAGraph_low,"2017 pp (|y|>1.5)","P");
			leg3->AddEntry(BPFONLL,"FONLL","f");
			leg3->AddEntry(BFONLLLow,"FONLL (|y| > 1.5)","f");
			leg3->Draw("same");
			MyPadr->Update();

	MyPadr2->cd();
	TLine * Unity2 = new TLine(0,0,0,0);
	if (meson_n=="BP"){Unity2 = new TLine(5,1,60,1);}
	else {Unity2 = new TLine(7,1,50,1);}
	Unity2->SetLineWidth(2);
	Unity2->SetLineStyle(2);
	Unity2->SetLineColor(1);
	Unity2->Draw("SAME");
	HisEmpty3->Draw();
	MyPadr2->Update();


 // Get ratio plots
float binlow[NBinsLow];
float bl_low[NBinsLow];
float bl_low_yStatL[NBinsLow];
float bl_low_xErrL[NBinsLow];
float bl_low_xErrH[NBinsLow];
float bl_low_ySystL[NBinsLow];
float binhigh[NBins-NBinsLow];
float bl_high[NBins-NBinsLow];
float bl_high_yStatL[NBins-NBinsLow];
float bl_high_xErrL[NBins-NBinsLow];
float bl_high_xErrH[NBins-NBinsLow];
float bl_high_ySystL[NBins-NBinsLow];

for (int i=0;i<NBins;++i){
	
	if( i<NBinsLow){
		binlow[i]=XsecPP_X[i];
		bl_low[i]=BPXsecPPY2DScaled[i];
		bl_low_yStatL[i]=BPXSecPPY2DErrDownScaled[i];
		bl_low_xErrL[i]=XsecPP_X[i]-ptBins[i];
		bl_low_xErrH[i]=ptBins[i + 1]-XsecPP_X[i];
		bl_low_ySystL[i]=BP2DTotalSystDownRatio[i]*BPXsecPPY2DScaled[i];
	} 
	else {
		binhigh[i-NBinsLow]=XsecPP_X[i];
		bl_high[i-NBinsLow]=BPXsecPPY2DScaled[i];
		bl_high_yStatL[i-NBinsLow]=BPXSecPPY2DErrDownScaled[i];
		bl_high_xErrL[i-NBinsLow]=XsecPP_X[i]-ptBins[i];
		bl_high_xErrH[i-NBinsLow]=ptBins[i + 1]-XsecPP_X[i];
		bl_high_ySystL[i-NBinsLow]=BP2DTotalSystDownRatio[i]*BPXsecPPY2DScaled[i];
	}
}

  vector<float> BXsec;
  vector<float> BXsecStat;
  vector<float> BXsecSyst;
 
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

  	std::vector<float> Unity(NBins, 1);

	//FONLL
  	vector<float> FONLL;
  	vector<float> FONLLUp;
  	vector<float> FONLLDown;
	double XTempFONLL;
	double YTempFONLL;
 	float FONLL_Norm_ratio[NBins];
  	float RatioBsStat[NBins];
  	float RatioBsSyst[NBins];
  	float RatioBsFonErrHigh[NBins];
  	float RatioBsFonErrLow[NBins];

	//BINS AND RESPECTIVE CONTENT NORMALIZED BY FONLL
	for (auto i = 0; i < NBins; ++i) {
		BPFONLL->GetPoint(i, XTempFONLL, YTempFONLL);
		FONLL.push_back(YTempFONLL);
		FONLLUp.push_back(BPFONLL->GetErrorYhigh(i));
		FONLLDown.push_back(BPFONLL->GetErrorYlow(i));
		FONLL_Norm_ratio[i] = BXsec[i] / FONLL[i];
		RatioBsStat[i] = BXsecStat[i] / FONLL[i];
		RatioBsSyst[i] = BXsecSyst[i] / FONLL[i];
		RatioBsFonErrHigh[i] = FONLLUp[i] / FONLL[i];
		RatioBsFonErrLow[i] = FONLLDown[i] / FONLL[i];
	}

	//DATA
	TGraphAsymmErrors gRatioBs_low(NBinsLow, binlow, FONLL_Norm_ratio, bl_low_xErrL, bl_low_xErrH, RatioBsStat, RatioBsStat);
	TGraphAsymmErrors gRatioBs_syst_low(NBinsLow,binlow,FONLL_Norm_ratio,bl_low_xErrL, bl_low_xErrH,RatioBsSyst, RatioBsSyst);
	TGraphAsymmErrors *BPRAAGraph_low_just_m = new TGraphAsymmErrors(NBinsLow, binlow, FONLL_Norm_ratio ,zero, zero, zero, zero);
	TGraphAsymmErrors gRatioBs_high(NBinsHigh,binhigh,FONLL_Norm_ratio + NBinsLow,bl_high_xErrL, bl_high_xErrH,RatioBsStat + NBinsLow, RatioBsStat + NBinsLow);
	TGraphAsymmErrors gRatioBs_syst_high(NBinsHigh,binhigh,FONLL_Norm_ratio + NBinsLow,bl_high_xErrL, bl_high_xErrH,RatioBsSyst + NBinsLow, RatioBsSyst + NBinsLow);
	
	//FONLL
	TGraphAsymmErrors gRatioBs_Fon_low(NBinsLow,binlow,Unity.data(),bl_low_xErrL, bl_low_xErrH,RatioBsFonErrLow, RatioBsFonErrHigh);
	TGraphAsymmErrors gRatioBs_Fon_high(NBinsHigh,binhigh,Unity.data(),bl_high_xErrL, bl_high_xErrH,RatioBsFonErrLow + NBinsLow,RatioBsFonErrHigh + NBinsLow);

	int color_mark =  kBlue + 2;
	int color_syst = kBlue -3;
	if(meson_n=="BP"){	
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
	cr->SaveAs(Form("Plots/%s_Xsection_%s_vsFONL.pdf", meson_n.Data(),whichvar.Data()));
	//FONLL
}
// vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL 
// vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL 


if(meson_n=="Bs" || BsBPBins == 1){}
// Save Histogram for the Ratio of Bmesons
	gSystem->mkdir("./ROOTFiles",true); 
	TFile *ratio_f= new TFile(Form("./ROOTFiles/%s_Xsection_%s.root", meson_n.Data(), whichvar.Data()),"recreate");
	ratio_f->cd();

	TMultiGraph* mg = new TMultiGraph();
	TGraphAsymmErrors* gr_staterr = new TGraphAsymmErrors(NBins, XsecPP_X, XsecPP_Y, XsecPP_X_BinLeft, XsecPP_X_BinRight, XsecPP_Y_StatDown, XsecPP_Y_StatUp);
	gr_staterr->SetName("Y_stat");
	gr_staterr->SetLineColor(1); 
	mg->Add(gr_staterr, "stat");

	TGraphAsymmErrors* gr_systerr = new TGraphAsymmErrors(NBins, XsecPP_X, XsecPP_Y, nullptr, nullptr, XsecPP_Y_SystDown, XsecPP_Y_SystUp);
	gr_systerr->SetName("Y_syst");
	gr_systerr->SetLineColor(2);
	mg->Add(gr_systerr,"syst");

	 if(whichvar=="y"){
		mg->GetXaxis()->SetTitle("Rapidity (y)");
		mg->GetYaxis()->SetTitle(Form("d#sigma/d%s [pb c/GeV]",whichvar.Data()));
		mg->GetXaxis()->SetLimits(0,2.4);
	 }
	 else if(whichvar=="pt"){
		mg->GetXaxis()->SetTitle("p_{T}");
		mg->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb c/GeV]");
		if (meson_n=="BP"){ mg->GetXaxis()->SetLimits(0 ,80); }
		if (meson_n=="Bs"){ mg->GetXaxis()->SetLimits(0 ,60); }
	 }
	 else if(whichvar=="Mult"){
		 mg->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 mg->GetYaxis()->SetTitle(Form("d#sigma/d%s [pb c/GeV]",whichvar.Data()));
		 mg->GetXaxis()->SetLimits(0, 100);
	 }

	 mg->Write(Form("%s_Xsection_%s", meson_n.Data(), whichvar.Data()));
	 
	 ratio_f->Close();
// Save Histogram for the Ratio of Bmesons






if(BsBPBins==0){
  // summary of errors (in ratio, not percent)
  gSystem->mkdir("../../../MakeFinalPlots/NominalPlots/CrossSection/dataSource/" ,true );
  string outFile;
  if(meson_n=="BP"){ outFile = "../../../MakeFinalPlots/NominalPlots/CrossSection/dataSource/corryield_pt_Bp_New.txt";}
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
      setprecision(0) << std::fixed << XsecPP_Y[i] << std::setw(columnWidth) <<
      setprecision(2) << std::defaultfloat <<
      BPXSecPPY2DErrUpRatio[i] << std::setw(columnWidth) <<
      BPXSecPPY2DErrDownRatio[i] << std::setw(columnWidth) <<
      BP2DTotalSystDownRatio[i] << std::setw(columnWidth) <<
      BP2DTotalSystDownRatio[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      setprecision(3) << XsecPP_X[i] << "\n";
  }
  out.close();
}


 }


						