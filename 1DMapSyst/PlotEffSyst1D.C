#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "../henri2022/parameter.h" 

using namespace std;
using std::cout;
using std::endl;

void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, 
		std::vector<std::vector<double> > numbers , std::string caption){
	
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

void PlotEffSyst1D(int Opt, int whichvar){

	TString var_n;
	TString var_l;
	TString var_N;
	int NBins = 7;

	if(whichvar==0 && Opt==0){
		var_n="pt";
		var_l="p_{T}";
		var_N="PT";
		NBins = nptBinsBP;
	}
	if(whichvar==0 && Opt==1){
		var_n="pt";
		var_l="p_{T}";
		var_N="PT";
		NBins = nptBins;
	}

	if(whichvar==1){
		var_n="y";
		var_l="Rapidity";
		var_N="Y";
		NBins = nyBins_both;
	}
	if(whichvar==2){
		var_n="Mult";
		var_l="Multiplicity";
		var_N="Mult";
		NBins = nmBins_both;
	}

	double ptBins[NBins+1];
	if (whichvar==1){
	for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ybinsvec[i];             
		}
	} 
	if (whichvar==2){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  nmbinsvec[i];             
		}
	}
	if (whichvar==0 && Opt==0){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvecBP[i];             
		}
	}
	if (whichvar==0 && Opt!=0){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvec[i];             
		}
	}

	TCanvas * c = new TCanvas("c","c",700,700);
	c->cd();

	gStyle->SetOptStat(0);

	TString infile;
  	TString outFile;
	TString BmesonName;
	if(Opt == 0) {
	BmesonName =  "BP" ;
    infile =  Form("OutFiles/BPSyst1D_%s.root",var_n.Data());
    outFile =  Form("OutFiles/BPError1D_%s.root",var_n.Data());
  } else if(Opt == 1) {
	BmesonName =  "Bs";
    infile =  Form("OutFiles/BsSyst1D_%s.root",var_n.Data());
    outFile =  Form("OutFiles/BsError1D_%s.root",var_n.Data());
  }

	gSystem->mkdir("EffSystPlots", true);
	gSystem->mkdir(Form("EffSystPlots/%s",BmesonName.Data()), true);

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TH1D * Eff1DHis = (TH1D * ) fin->Get("Eff1DHisInv");
	TH1D * Eff1DHisTnPUp = (TH1D * ) fin->Get("Eff1DTnPUpSystHis");
	TH1D * Eff1DHisTnPDown = (TH1D * ) fin->Get("Eff1DTnPDownSystHis");
	TH1D * Eff1DHisBpt = (TH1D * ) fin->Get("Eff1DBptHis");
	TH1D * Eff1DHisBDT = (TH1D * ) fin->Get("Eff1DBDTHis");

  	TFile fout(outFile, "recreate");

	//Draw Eff2DHis
	TCanvas * EFFplot  = new TCanvas("cSyst","cSyst",600,600);
	EFFplot->cd();
	Eff1DHis->SetMarkerStyle(20);
	Eff1DHis->SetMarkerSize(1);
	Eff1DHis->SetMarkerColor(kBlack);
	Eff1DHis->SetLineColor(kBlack);
	Eff1DHis->Draw("ep");
	EFFplot->SaveAs(Form("EffSystPlots/%s/%s_EFFplot.pdf",BmesonName.Data(),var_n.Data()));


	//Draw Systematic Uncertainties
	TCanvas * cSyst  = new TCanvas("cSyst","cSyst",600,600);
	cSyst->cd();
	Eff1DHis->SetMarkerStyle(20);
	Eff1DHis->SetMarkerSize(1);
	Eff1DHis->SetMarkerColor(kBlack);
	Eff1DHis->SetLineColor(kBlack);
	Eff1DHisTnPUp->SetMarkerStyle(20);
	Eff1DHisTnPUp->SetMarkerSize(1);
	Eff1DHisTnPUp->SetMarkerColor(kRed);
	Eff1DHisTnPUp->SetLineColor(kRed);
	Eff1DHisTnPDown->SetMarkerStyle(20);
	Eff1DHisTnPDown->SetMarkerSize(1);
	Eff1DHisTnPDown->SetMarkerColor(kBlue);
	Eff1DHisTnPDown->SetLineColor(kBlue);
	Eff1DHisTnPUp->Draw("ep");
	Eff1DHis->Draw("epSAME");
	Eff1DHisTnPDown->Draw("epSAME");

	TLegend* leg = new TLegend(0.62,0.77,0.89,0.89,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.025);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(Eff1DHis,"Nominal","PL");
	leg->AddEntry(Eff1DHisTnPUp,"TnP Variation Up","PL");
	leg->AddEntry(Eff1DHisTnPDown,"TnP Variation Down","PL");
	leg->Draw("same");

	cSyst->SaveAs(Form("EffSystPlots/%s/%s_TnPSystComp.pdf",BmesonName.Data(),var_n.Data()));


	Eff1DHisBDT->SetMarkerStyle(20);
	Eff1DHisBDT->SetMarkerSize(1);
	Eff1DHisBDT->SetMarkerColor(kRed);
	Eff1DHisBDT->SetLineColor(kRed);
	Eff1DHis->Draw("ep");
	Eff1DHisBDT->Draw("epSAME");

	TLegend* leg2 = new TLegend(0.7,0.8,0.89,0.89,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.025);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->AddEntry(Eff1DHis,"Nominal","PL");
	leg2->AddEntry(Eff1DHisBDT,"BDT Weighted","PL");
	leg2->Draw("same");

	cSyst->SaveAs(Form("EffSystPlots/%s/%s_MCDataSystComp.pdf",BmesonName.Data(),var_n.Data()));

	Eff1DHisBpt->SetMarkerStyle(20);
	Eff1DHisBpt->SetMarkerSize(1);
	Eff1DHisBpt->SetMarkerColor(kRed);
	Eff1DHisBpt->SetLineColor(kRed);

	if (Opt==0){
		Eff1DHisBpt->Draw("ep");
		Eff1DHis->Draw("epSAME");
	} else if (Opt==1) {
		Eff1DHis->Draw("ep");
		Eff1DHisBpt->Draw("epSame");
	}

	TLegend* leg3 = new TLegend(0.7,0.8,0.89,0.89,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.025);
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->AddEntry(Eff1DHis,"Nominal","PL");
	leg3->AddEntry(Eff1DHisBpt,"Bpt Weighted","PL");
	leg3->Draw("same");

	cSyst->SaveAs(Form("EffSystPlots/%s/%s_BptSystComp.pdf",BmesonName.Data(),var_n.Data()));

	//Done drawing 


	TH1D * TnPSyst = (TH1D *) Eff1DHisTnPUp->Clone("TnPSyst");
	TnPSyst->GetYaxis()->SetTitle("TnP Systematic Uncertainties (%)");
	TnPSyst->GetYaxis()->SetTitleOffset(1.3);
	TnPSyst->SetLineColor(kBlack);
	TnPSyst->SetMarkerColor(kBlack);	
	TnPSyst->Reset();

	TH1D * BptSyst = (TH1D *) Eff1DHisBpt->Clone("BptSyst");
	BptSyst->GetYaxis()->SetTitle(Form("B-meson %s Shape Systematic Uncertainties (%s)",var_l.Data(),"%"));
	BptSyst->GetYaxis()->SetTitleOffset(1.3);
	BptSyst->SetLineColor(kBlack);
	BptSyst->SetMarkerColor(kBlack);	
	BptSyst->Reset();

	TH1D * BDTSyst = (TH1D *) Eff1DHisBDT->Clone("BDTSyst");
	BDTSyst->GetYaxis()->SetTitle("MC/Data Discrepancy Systematic Uncertainties (%)");
	BDTSyst->GetYaxis()->SetTitleOffset(1.3);
	BDTSyst->SetLineColor(kBlack);
	BDTSyst->SetMarkerColor(kBlack);	
	BDTSyst->Reset();

	c->cd();


	float SystValue;
	float SystValueError;

	for(int i = 0; i < Eff1DHisTnPUp->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisTnPUp->GetBinContent(i+1) -  Eff1DHis->GetBinContent(i+1) )/Eff1DHis->GetBinContent(i+1) * 100;
		cout << "TnP Syst: " << SystValue << endl;
		TnPSyst->SetBinContent(i+1,SystValue);
		TnPSyst->SetBinError(i+1,0.001);
	}

	for(int i = 0; i < Eff1DHisTnPUp->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBpt->GetBinContent(i+1) -  Eff1DHis->GetBinContent(i+1) )/Eff1DHis->GetBinContent(i+1) * 100;
		BptSyst->SetBinContent(i+1,SystValue);
		BptSyst->SetBinError(i+1,0.001);
		cout << "Bpt Syst: " << SystValue << endl;
	}

	for(int i = 0; i < Eff1DHisTnPUp->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBDT->GetBinContent(i+1) -  Eff1DHis->GetBinContent(i+1) )/Eff1DHis->GetBinContent(i+1) * 100;
		BDTSyst->SetBinContent(i+1,SystValue);
		BDTSyst->SetBinError(i+1,0.001);
		cout << "BDT Syst: " << SystValue << endl;
    	cout << SystValue << "\n";

	}

	TnPSyst->Draw("ep");
	c->SaveAs(Form("EffSystPlots/%s/%s_TnPSystRatio.pdf",BmesonName.Data(),var_n.Data()));
	BptSyst->Draw("ep");
	c->SaveAs(Form("EffSystPlots/%s/%s_BptSysRatio.pdf",BmesonName.Data(),var_n.Data()));
	BDTSyst->Draw("ep");
	c->SaveAs(Form("EffSystPlots/%s/%s_MCDataSystRatio.pdf",BmesonName.Data(),var_n.Data()));

  	TnPSyst->Write();
  	BptSyst->Write();
  	BDTSyst->Write();

	//eff with syst (comparison with 2D)
	TString InfileB = Form("../EffAna/%s/FinalFiles/%sPPCorrYield%s.root",BmesonName.Data(),BmesonName.Data(),var_N.Data());
	TFile * FileB= new TFile(InfileB.Data());

	TString errorFile = Form("../2DMapSyst/OutFiles/%sError2D_%s.root", BmesonName.Data(),var_n.Data());
  	TFile fError(errorFile);

	TString trackSelErrorFile = Form("../syst_error/syst_track_sel_%s.root",var_n.Data());
	TFile *fTrackSelError(trackSelErrorFile);

	TH1D * Eff2DHis = (TH1D *) FileB->Get("hInvEff");
	TH1D * TnPSyst2D = (TH1D *) fError.Get("TnPSyst");
	TH1D * BptSyst2D = (TH1D *) fError.Get("BptSyst");
	TH1D * MCDataSyst2D = (TH1D *) fError.Get("BDTSyst");

	TGraphAsymmErrors* trackSelSyst =  dynamic_cast<TGraphAsymmErrors*>( fTrackSelError->Get( Form("%s_track_sel_error", BmesonName.Data()) ) );
	TGraphAsymmErrors* trackSelSystMC= dynamic_cast<TGraphAsymmErrors*>( fTrackSelError->Get( Form("%s_track_sel_error", BmesonName.Data()) ) );

	float BP1DEffX[NBins];
	float B1DEffXErrUp[NBins] ;
	float BP1DEffY[NBins];
	float BP1DEffYErrUp[NBins];

	float BP2DEffY[NBins];
	float BP2DEffYErrUp[NBins];

	for( int c=0; c <NBins; c++){ 
		BP1DEffX[c]=(ptBins[c]+ptBins[c+1])/2;
		B1DEffXErrUp[c]=(abs(ptBins[c+1]-ptBins[c]))/2;
		BP1DEffY[c] = Eff1DHis->GetBinContent(c+1);
		BP1DEffYErrUp[c] = Eff1DHis->GetBinError(c+1);
		BP2DEffY[c] = Eff2DHis->GetBinContent(c+1);
		BP2DEffYErrUp[c] = Eff2DHis->GetBinError(c+1);
		cout<< "1D eff:"<<BP1DEffY[c]<<endl;
		cout<< "2D eff:"<<BP2DEffY[c]<<endl;
	}

	float BP1DEffYSystUp[NBins];
	float BP1DEffYSystDown[NBins];

	float BP2DEffYSystUp[NBins];
	float BP2DEffYSystDown[NBins];

	double B_nu;
  	if (Opt == 0){ B_nu = 2.4 ;}
  	else { B_nu = 4.8;}
	double BPTrackingSyst[NBins];
	for( int c=0; c < NBins; c++){BPTrackingSyst[c]= B_nu ;}
	float BPPDFSyst[NBins];
	
	float BPMCDataSyst1D[NBins];
	float BPPtShapeSyst1D[NBins];
	float BPTnPSystDown1D[NBins];
	float BPTnPSystUp1D[NBins];

	float BPMCDataSyst2D[NBins];
	float BPPtShapeSyst2D[NBins];
	float BPTnPSystDown2D[NBins];
	float BPTnPSystUp2D[NBins];
	float BPTrackSelSyst[NBins];
	float BPTrackSelSystMC[NBins];

  // Get systematics from input files
    for (auto ibin = 0; ibin < NBins; ++ibin){
		BPMCDataSyst1D[ibin] = BDTSyst->GetBinContent(ibin + 1);
		BPPtShapeSyst1D[ibin] = BptSyst->GetBinContent(ibin + 1);
		BPTnPSystDown1D[ibin] = TnPSyst->GetBinContent(ibin + 1);
		BPTnPSystUp1D[ibin] = BPTnPSystDown1D[ibin];

	    BPMCDataSyst2D[ibin] = MCDataSyst2D->GetBinContent(ibin + 1);
		BPPtShapeSyst2D[ibin] = BptSyst2D->GetBinContent(ibin + 1);
		BPTnPSystDown2D[ibin] = TnPSyst2D->GetBinContent(ibin + 1);
		BPTnPSystUp2D[ibin] = BPTnPSystDown2D[ibin];
		BPTrackSelSyst[ibin] = trackSelSyst->GetY()[ibin];
		BPTrackSelSystMC[ibin] = trackSelSystMC->GetY()[ibin];
  	}

  // RMS of all the errors
	float BPTotalSystDownRatio1D[NBins];
	float BPTotalSystUpRatio1D[NBins];

	float BPTotalSystDownRatio2D[NBins];
	float BPTotalSystUpRatio2D[NBins];

	for(int i = 0; i < NBins; i++){

		BPTotalSystDownRatio1D[i] = TMath::Sqrt(TMath::Power(BPMCDataSyst1D[i], 2) + TMath::Power(BPPtShapeSyst1D[i], 2) + TMath::Power(BPTnPSystDown1D[i], 2)) / 100;
    	BPTotalSystUpRatio1D[i] = TMath::Sqrt(TMath::Power(BPMCDataSyst1D[i], 2) + TMath::Power(BPPtShapeSyst1D[i], 2) + TMath::Power(BPTnPSystUp1D[i], 2)) / 100;

		//BPTotalSystDownRatio1D[i] = TMath::Sqrt(TMath::Power(TMath::Power(BPTrackingSyst[i], 2) + BPMCDataSyst1D[i], 2) + TMath::Power(BPTrackSelSystMC[i], 2) + TMath::Power(BPPtShapeSyst1D[i], 2) + TMath::Power(BPTnPSystDown1D[i], 2)) / 100;
    	//BPTotalSystUpRatio1D[i] = TMath::Sqrt(TMath::Power(TMath::Power(BPTrackingSyst[i], 2) + BPMCDataSyst1D[i], 2) + TMath::Power(BPTrackSelSystMC[i], 2) + TMath::Power(BPPtShapeSyst1D[i], 2) + TMath::Power(BPTnPSystUp1D[i], 2)) / 100;

		//BPTotalSystDownRatio2D[i] = TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BPMCDataSyst2D[i], 2) + TMath::Power(BPTrackSelSyst[i], 2) + TMath::Power(BPPtShapeSyst2D[i], 2) + TMath::Power(BPTnPSystDown2D[i], 2)) / 100;
    	//BPTotalSystUpRatio2D[i] = TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BPMCDataSyst2D[i], 2) + TMath::Power(BPTrackSelSyst[i], 2) + TMath::Power(BPPtShapeSyst2D[i], 2) + TMath::Power(BPTnPSystUp2D[i], 2)) / 100;

		BPTotalSystDownRatio2D[i] = TMath::Sqrt(TMath::Power(BPMCDataSyst2D[i], 2) + TMath::Power(BPPtShapeSyst2D[i], 2) + TMath::Power(BPTnPSystDown2D[i], 2)) / 100;
    	BPTotalSystUpRatio2D[i] = TMath::Sqrt(TMath::Power(BPMCDataSyst2D[i], 2) + TMath::Power(BPPtShapeSyst2D[i], 2) + TMath::Power(BPTnPSystUp2D[i], 2)) / 100;	
	}

	for(int i = 0; i < NBins; i++){
		BP1DEffYSystUp[i] = BP1DEffY[i] * BPTotalSystUpRatio1D[i];
		BP1DEffYSystDown[i] = BP1DEffY[i] * BPTotalSystDownRatio1D[i];

		BP2DEffYSystUp[i] = BP2DEffY[i] * BPTotalSystUpRatio2D[i];
		BP2DEffYSystDown[i] = BP2DEffY[i] * BPTotalSystDownRatio2D[i];
	}

	TGraphAsymmErrors *BP1DEffGraph = new TGraphAsymmErrors(NBins, BP1DEffX ,BP1DEffY , B1DEffXErrUp, B1DEffXErrUp, BP1DEffYErrUp, BP1DEffYErrUp); 
	TGraphAsymmErrors *BP1DEffGraphSyst = new TGraphAsymmErrors(NBins, BP1DEffX , BP1DEffY, B1DEffXErrUp, B1DEffXErrUp, BP1DEffYSystDown, BP1DEffYSystUp);

	TGraphAsymmErrors *BP2DEffGraph = new TGraphAsymmErrors(NBins, BP1DEffX ,BP2DEffY , B1DEffXErrUp, B1DEffXErrUp, BP2DEffYErrUp, BP2DEffYErrUp); 
	TGraphAsymmErrors *BP2DEffGraphSyst = new TGraphAsymmErrors(NBins, BP1DEffX , BP2DEffY, B1DEffXErrUp, B1DEffXErrUp, BP2DEffYSystDown, BP2DEffYSystUp);

	gStyle->SetOptStat(0);
	cSyst->cd(); 
	cSyst->SetLeftMargin(0.15);

	//Setup histograms for different purposs
	TH2D * HisEmpty;
	if(Opt == 0 && whichvar==0) {HisEmpty = new TH2D("HisEmpty","",100,5,60,100,0,60);} 
	if(Opt == 1 && whichvar==0) {HisEmpty = new TH2D("HisEmpty","",100,7,50,100,0,60);}

	if(Opt == 0 && whichvar==1) {HisEmpty = new TH2D("HisEmpty","",100,-2.4,2.4,100,0,60);}
	if(Opt == 1 && whichvar==1) {HisEmpty = new TH2D("HisEmpty","",100,-2.4,2.4,100,0,60);}

	if(Opt == 0 && whichvar==2) {HisEmpty = new TH2D("HisEmpty","",100,0,100,100,0,60);}   // need to adjust range for when we have nmult results
	if(Opt == 1 && whichvar==2) {HisEmpty = new TH2D("HisEmpty","",100,0,100,100,0,60);}

	HisEmpty->GetXaxis()->SetTitle(var_l.Data());
	HisEmpty->GetYaxis()->SetTitle("1/(#alpha #times #epsilon)");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();

	BP1DEffGraph->SetMarkerStyle(20);
	BP1DEffGraph->SetMarkerSize(1);
	BP1DEffGraph ->SetMarkerColor(kGreen+2);
	BP1DEffGraph ->SetLineColor(kGreen+2);
	BP1DEffGraphSyst ->SetFillColorAlpha(kGreen-7,0.5);
	BP1DEffGraphSyst ->SetLineColor(kGreen-7);

	BP2DEffGraph->SetMarkerStyle(20);
	BP2DEffGraph->SetMarkerSize(1);
	BP2DEffGraph->SetLineColor(kBlue+2);
	BP2DEffGraph->SetMarkerColor(kBlue+2);
	BP2DEffGraphSyst->SetFillColorAlpha(kBlue-3,0.5);
	BP2DEffGraphSyst->SetLineColor(kBlue-3);
	
	HisEmpty->Draw();
	BP2DEffGraphSyst->Draw("5same");
	BP2DEffGraph->Draw("epSAME");
	BP1DEffGraphSyst->Draw("5same");
	BP1DEffGraph->Draw("epSAME");


	
	TLegend* leged;
	leged = new TLegend(0.65,0.6,0.9,0.68,NULL,"brNDC");
	leged->SetBorderSize(0);
	leged->SetFillStyle(0);
	leged->AddEntry(BP1DEffGraph,"1D Eff.","P");
	leged->AddEntry(BP2DEffGraph,"2D Eff.","P");
	leged->SetTextSize(0.022);
	leged->Draw("same");
	

	cSyst->SaveAs(Form("EffSystPlots/%s/Effcomp_%s.pdf", BmesonName.Data(), var_n.Data()));
	string name;
	TString whichvarname;
	if(whichvar==0){name="$<p_T<$"; whichvarname="pt";} else if(whichvar==1){name="$<y<$";whichvarname="y";} else if(whichvar==2){name="$<nTrks<$";whichvarname="nMult";}
	std::vector<std::string> col_name;
	std::vector<std::string> labels = {"","1D value","Stat error","Syst error", "2D value","Stat error","Syst error"};
	std::vector<std::vector<double>> eff_values;
	std::vector<double> val1D;
	std::vector<double> stat1D;
	std::vector<double> syst1D;
	std::vector<double> val2D;
	std::vector<double> stat2D;
	std::vector<double> syst2D;
	for(int i=0;i<NBins;i++){

		std::ostringstream clabel;
		clabel<<ptBins[i]<<name<<ptBins[i+1];
		std::string label1 = clabel.str();
		col_name.push_back(label1);

		val1D.push_back(BP1DEffY[i]);
		stat1D.push_back(BP1DEffYErrUp[i]/BP1DEffY[i]);
		syst1D.push_back(BP1DEffYSystDown[i]/BP1DEffY[i]);

		val2D.push_back(BP2DEffY[i]);
		stat2D.push_back(BP2DEffYErrUp[i]/BP2DEffY[i]);
		syst2D.push_back(BP2DEffYSystDown[i]/BP2DEffY[i]);

	}

	eff_values.push_back(val1D);
	eff_values.push_back(stat1D);
	eff_values.push_back(syst1D);
	eff_values.push_back(val2D);
	eff_values.push_back(stat2D);
	eff_values.push_back(syst2D);

	gSystem->mkdir("Trash",true);

	latex_table(Form("1D2Dcomparisons_%s",whichvarname.Data()), 7,  NBins+1,  labels , col_name , eff_values, "1D vs 2D efficiency comparisons");

	std::vector<std::string> filetype ={"_check.aux", "_check.log",".tex","_check.tex"};
	for (int j=0;j<(int)(filetype.size());j++){
				rename(("1D2Dcomparisons_"+ std::string (whichvarname.Data()) +filetype[j]).c_str(),("Trash/1D2Dcomparisons_"+std::string (whichvarname.Data())+filetype[j]).c_str());
			}
	rename(("1D2Dcomparisons_"+ std::string (whichvarname.Data()) +"_check.pdf").c_str(),("EffSystPlots/"+std::string (BmesonName.Data())+"/1D2Dcomparisons_"+std::string (whichvarname.Data())+"_check.pdf").c_str());


}


