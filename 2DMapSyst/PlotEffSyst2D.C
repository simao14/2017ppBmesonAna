#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;
using std::cout;
using std::endl;

void PlotEffSyst2D(int Opt){

	TCanvas * c = new TCanvas("c","c",700,700);
	c->cd();

	gStyle->SetOptStat(0);

	TString infile;
  TString outFile;

	if(Opt == 0) {
    infile =  "OutFiles/BPSyst2D.root";
    outFile =  "OutFiles/BPError2D.root";
  } else if(Opt == 1) {
    infile =  "OutFiles/BsSyst2D.root";
    outFile =  "OutFiles/BsError2D.root";
  }

	TString BmesonName;
	if(Opt == 0) BmesonName =  "BP";
	if(Opt == 1) BmesonName =  "Bs";
	gSystem->mkdir(Form("SystPlots/%s/Pt",BmesonName.Data()), true);

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TH1D * Eff1DHis = (TH1D * ) fin->Get("Eff2DHis");
	TH1D * Eff1DHisTnPUp = (TH1D * ) fin->Get("Eff2DTnPUpSystHis");
	TH1D * Eff1DHisTnPDown = (TH1D * ) fin->Get("Eff2DTnPDownSystHis");
	TH1D * Eff1DHisBpt = (TH1D * ) fin->Get("Eff2DBptHis");
	TH1D * Eff1DHisBDT = (TH1D * ) fin->Get("Eff2DBDTHis");

  TFile fout(outFile, "recreate");

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

	TLegend* leg = new TLegend(0.80,0.75,0.89,0.89,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.025);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(Eff1DHis,"Nominal","PL");
	leg->AddEntry(Eff1DHisTnPUp,"T&P Variation Up","PL");
	leg->AddEntry(Eff1DHisTnPDown,"T&P Variation Down","PL");
	leg->Draw("same");

	cSyst->SaveAs(Form("SystPlots/%s/Pt/TnPSystComp.png",BmesonName.Data()));


	Eff1DHisBDT->SetMarkerStyle(20);
	Eff1DHisBDT->SetMarkerSize(1);
	Eff1DHisBDT->SetMarkerColor(kRed);
	Eff1DHisBDT->SetLineColor(kRed);
	Eff1DHis->Draw("ep");
	Eff1DHisBDT->Draw("epSAME");

	TLegend* leg2 = new TLegend(0.80,0.75,0.89,0.89,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.025);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->AddEntry(Eff1DHis,"Nominal","PL");
	leg2->AddEntry(Eff1DHisBDT,"BDT Weighted","PL");
	leg2->Draw("same");

	cSyst->SaveAs(Form("SystPlots/%s/Pt/MCDataSystComp.png",BmesonName.Data()));

	Eff1DHisBpt->SetMarkerStyle(20);
	Eff1DHisBpt->SetMarkerSize(1);
	Eff1DHisBpt->SetMarkerColor(kRed);
	Eff1DHisBpt->SetLineColor(kRed);
	Eff1DHis->Draw("ep");
	Eff1DHisBpt->Draw("epSAME");

	TLegend* leg3 = new TLegend(0.80,0.75,0.89,0.89,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.025);
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->AddEntry(Eff1DHis,"Nominal","PL");
	leg3->AddEntry(Eff1DHisBpt,"Bpt Weighted","PL");
	leg3->Draw("same");

	cSyst->SaveAs(Form("SystPlots/%s/Pt/BptSystComp.png",BmesonName.Data()));

	//Done drawing 


	TH1D * TnPSyst = (TH1D *) Eff1DHisTnPUp->Clone("TnPSyst");
	TnPSyst->GetYaxis()->SetTitle("TnP Systematic Uncertainties (%)");
	TnPSyst->GetYaxis()->SetTitleOffset(1.3);
	TnPSyst->SetLineColor(kBlack);
	TnPSyst->SetMarkerColor(kBlack);	
	TnPSyst->Reset();

	TH1D * BptSyst = (TH1D *) Eff1DHisBpt->Clone("BptSyst");
	BptSyst->GetYaxis()->SetTitle("B-meson p_{T} Shape Systematic Uncertainties (%)");
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
	c->SaveAs(Form("SystPlots/%s/Pt/TnPSystRatio.png",BmesonName.Data()));
	BptSyst->Draw("ep");
	c->SaveAs(Form("SystPlots/%s/Pt/BptSysRatio.png",BmesonName.Data()));
	BDTSyst->Draw("ep");
	c->SaveAs(Form("SystPlots/%s/Pt/MCDataSystRatio.png",BmesonName.Data()));

  TnPSyst->Write();
  BptSyst->Write();
  BDTSyst->Write();


/*
	TH1D * Eff1DHisMult = (TH1D * ) fin->Get("Eff1DHisMult");
	TH1D * Eff1DHisTnPUpMult = (TH1D * ) fin->Get("Eff1DHisTnPUpMult");
	TH1D * Eff1DHisBptMult = (TH1D * ) fin->Get("Eff1DHisBptMult");
	TH1D * Eff1DHisBDTMult = (TH1D * ) fin->Get("Eff1DHisBDTMult");
	TH1D * Eff1DHisTnPDownMult = (TH1D * ) fin->Get("Eff1DHisTnPDownMult");



	//Draw Systematic Uncertainties Comparison

		cSyst->cd();
		Eff1DHisMult->SetMarkerStyle(20);
		Eff1DHisMult->SetMarkerSize(1);
		Eff1DHisMult->SetMarkerColor(kBlack);
		Eff1DHisMult->SetLineColor(kBlack);


		Eff1DHisTnPUpMult->SetMarkerStyle(20);
		Eff1DHisTnPUpMult->SetMarkerSize(1);
		Eff1DHisTnPUpMult->SetMarkerColor(kRed);
		Eff1DHisTnPUpMult->SetLineColor(kRed);


		Eff1DHisTnPDownMult->SetMarkerStyle(20);
		Eff1DHisTnPDownMult->SetMarkerSize(1);
		Eff1DHisTnPDownMult->SetMarkerColor(kBlue);
		Eff1DHisTnPDownMult->SetLineColor(kBlue);



		Eff1DHisTnPUpMult->Draw("ep");
		Eff1DHisMult->Draw("epSAME");
		Eff1DHisTnPDownMult->Draw("epSAME");

		TLegend* legMult = new TLegend(0.50,0.35,0.80,0.60,NULL,"brNDC");
		legMult->SetBorderSize(0);
		legMult->SetTextSize(0.040);
		legMult->SetTextFont(42);
		legMult->SetFillStyle(0);
		legMult->SetLineWidth(3);
		legMult->AddEntry(Eff1DHisMult,"Nominal","PL");
		legMult->AddEntry(Eff1DHisTnPUpMult,"T&P Variation Up","PL");
		legMult->AddEntry(Eff1DHisTnPDownMult,"T&P Variation Down","PL");
		legMult->Draw("same");

		cSyst->SaveAs(Form("SystPlots/%s/Mult/TnPSystComp.png",BmesonName.Data()));



		Eff1DHisBDTMult->SetMarkerStyle(20);
		Eff1DHisBDTMult->SetMarkerSize(1);
		Eff1DHisBDTMult->SetMarkerColor(kRed);
		Eff1DHisBDTMult->SetLineColor(kRed);

		Eff1DHisMult->Draw("ep");
		Eff1DHisBDTMult->Draw("epSAME");

		TLegend* leg2Mult = new TLegend(0.50,0.35,0.80,0.60,NULL,"brNDC");
		leg2Mult->SetBorderSize(0);
		leg2Mult->SetTextSize(0.040);
		leg2Mult->SetTextFont(42);
		leg2Mult->SetFillStyle(0);
		leg2Mult->SetLineWidth(3);
		leg2Mult->AddEntry(Eff1DHis,"Nominal","PL");
		leg2Mult->AddEntry(Eff1DHisBDT,"BDT Weighted","PL");
		leg2Mult->Draw("same");



		cSyst->SaveAs(Form("SystPlots/%s/Mult/MCDataSystComp.png",BmesonName.Data()));






		Eff1DHisBptMult->SetMarkerStyle(20);
		Eff1DHisBptMult->SetMarkerSize(1);
		Eff1DHisBptMult->SetMarkerColor(kRed);
		Eff1DHisBptMult->SetLineColor(kRed);

		Eff1DHisMult->Draw("ep");
		Eff1DHisBptMult->Draw("epSAME");

		TLegend* leg3Mult = new TLegend(0.50,0.35,0.80,0.60,NULL,"brNDC");
		leg3Mult->SetBorderSize(0);
		leg3Mult->SetTextSize(0.040);
		leg3Mult->SetTextFont(42);
		leg3Mult->SetFillStyle(0);
		leg3Mult->SetLineWidth(3);
		leg3Mult->AddEntry(Eff1DHis,"Nominal","PL");
		leg3Mult->AddEntry(Eff1DHisBpt,"Bpt Weighted","PL");
		leg3Mult->Draw("same");

		
		cSyst->SaveAs(Form("SystPlots/%s/Mult/BptSystComp.png",BmesonName.Data()));



	//Done Drawing



	TH1D * TnPSystMult = (TH1D *) Eff1DHisTnPUp->Clone("TnPSyst");
	TnPSystMult->GetYaxis()->SetTitle("TnP Systematic Uncertainties (%)");
	TnPSystMult->GetYaxis()->SetTitleOffset(1.3);

	TnPSystMult->SetLineColor(kBlack);
	TnPSystMult->SetMarkerColor(kBlack);	
	TnPSystMult->Reset();


	TH1D * BptSystMult = (TH1D *) Eff1DHisBpt->Clone("BptSystMult");
	BptSystMult->GetYaxis()->SetTitle("B-meson p_{T} Shape Systematic Uncertainties (%)");
	BptSystMult->GetYaxis()->SetTitleOffset(1.3);

	BptSystMult->SetLineColor(kBlack);
	BptSystMult->SetMarkerColor(kBlack);	
	BptSystMult->Reset();

	TH1D * BDTSystMult = (TH1D *) Eff1DHisBDT->Clone("BDTSystMult");
	BDTSystMult->GetYaxis()->SetTitle("BDT Systematic Uncertainties (%)");
	BDTSystMult->GetYaxis()->SetTitleOffset(1.3);

	BDTSystMult->SetLineColor(kBlack);
	BDTSystMult->SetMarkerColor(kBlack);	
	BDTSystMult->Reset();


	c->cd();
	
	for(int i = 0; i < Eff1DHisTnPUpMult->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisTnPUpMult->GetBinContent(i+1) -  Eff1DHisMult->GetBinContent(i+1) )/Eff1DHisMult->GetBinContent(i+1) * 100;




		cout << "TnP Syst: " << SystValue << endl;

		TnPSystMult->SetBinContent(i+1,SystValue);
		TnPSystMult->SetBinError(i+1,0.001);



	}




	for(int i = 0; i < Eff1DHisTnPUpMult->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBptMult->GetBinContent(i+1) -  Eff1DHisMult->GetBinContent(i+1) )/Eff1DHisMult->GetBinContent(i+1) * 100;

		BptSystMult->SetBinContent(i+1,SystValue);
		BptSystMult->SetBinError(i+1,0.001);

		cout << "Bpt Syst: " << SystValue << endl;



	}




	for(int i = 0; i < Eff1DHisTnPUpMult->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBDTMult->GetBinContent(i+1) -  Eff1DHisMult->GetBinContent(i+1) )/Eff1DHisMult->GetBinContent(i+1) * 100;

		BDTSystMult->SetBinContent(i+1,SystValue);
		BDTSystMult->SetBinError(i+1,0.001);

		cout << "BDT Syst: " << SystValue << endl;



	}








	TnPSystMult->Draw("ep");

	c->SaveAs(Form("SystPlots/%s/Mult/TnPSystRatio.png",BmesonName.Data()));



	BptSystMult->Draw("ep");
	c->SaveAs(Form("SystPlots/%s/Mult/BptSystRatio.png",BmesonName.Data()));




	BDTSystMult->Draw("ep");
	c->SaveAs(Form("SystPlots/%s/Mult/MCDataSystRatio.png",BmesonName.Data()));

*/

}


