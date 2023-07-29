#include "TROOT.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
using namespace std;

using std::cout;
using std::endl;

void PlotBsBPRatio(int Opt){

	gSystem->mkdir("Plots/2DMap", true);

	TLatex* tex = new TLatex();
	tex->SetTextFont(42);
	tex->SetTextSize(0.03);
	tex->SetLineWidth(2);

	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetLeftMargin(0.16);

	TString InfileBP = "../../BP/EffAna/FinalFiles/BPPPCorrYieldMult.root";

	TString FileName;
	if(Opt == 0) FileName = "CorrDiffHisBin";
	if(Opt == 1) FileName = "hPtSigma";

	TString Method;
	if(Opt == 0) Method = "Binned";
	if(Opt == 1) Method = "2DMap";

	TFile * FileBP = new TFile(InfileBP.Data());
	TH1D * BPCross = (TH1D *) FileBP->Get(FileName.Data());
	BPCross->SetMarkerStyle(20);
	BPCross->SetMarkerSize(1);
	BPCross->SetMarkerColor(1);
	BPCross->SetLineColor(1);

	//B+ PP//
	const int NBins = 10;
	float BPXsecPPY[NBins];
	float BPXsecPPX[NBins] = {7.5,20,27.5,32.5,37.5,45,57.5,72.5,90,115};
	float BPXSecPPYErrUp[NBins];
	float BPXSecPPYErrDown[NBins];
	float BPXSecPPYErrUpPercent[NBins];
	float BPXSecPPYErrDownPercent[NBins];
	float BPXSecPPXErrUp[NBins] = {7.5,5,2.5,2.5,2.5,5,7.5,7.5,10,15};
	float BPXSecPPXErrDown[NBins] = {7.5,5,2.5,2.5,2.5,5,7.5,7.5,10,15};


	for(int i = 0; i < NBins; i++){
		BPXsecPPY[i] = BPCross->GetBinContent(i+1);
		BPXSecPPYErrUp[i] = BPCross->GetBinError(i+1);
		BPXSecPPYErrDown[i] = BPCross->GetBinError(i+1);
		BPXSecPPYErrUpPercent[i] = BPXSecPPYErrUp[i]/BPXsecPPY[i];
		BPXSecPPYErrDownPercent[i] = BPXSecPPYErrDown[i]/BPXsecPPY[i];

		cout << "BPXsecPPY = " << BPXsecPPY[i] << endl;
	}

	TH2D * HisEmpty = new TH2D("HisEmpty","",100,0,130,100,20.0,600000);
	HisEmpty->GetXaxis()->SetTitle("Event Multiplicity");
	HisEmpty->GetYaxis()->SetTitle("#Delta #sigma (pb)");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	HisEmpty->GetYaxis()->SetTitleOffset(1.4);
	HisEmpty->Draw();

	TGraphAsymmErrors *BPPPCrossGraph = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY,BPXSecPPXErrDown, BPXSecPPXErrUp,BPXSecPPYErrDown,BPXSecPPYErrUp);
	BPPPCrossGraph->SetLineColor(kBlue+2);
	BPPPCrossGraph->SetMarkerStyle(21);
	BPPPCrossGraph->SetMarkerSize(1);
	BPPPCrossGraph->SetMarkerColor(kBlue+2);

	BPPPCrossGraph->Draw("ep");

	c->SaveAs(Form("Plots/%s/BPPPCrossMult.png",Method.Data()));
	c->SetLogy();

	tex->DrawLatex(0.9,1.8,"p_{T} < 100 GeV/c");
	tex->DrawLatex(1.2,0.79,"|y| < 2.4 for 10 < p_{T} < 100 GeV/c");
	c->SaveAs(Form("Plots/%s/BPPPCrossMultLog.png",Method.Data()));

	//Bs PP

	TString InfileBs = "../../Bs/EffAna/FinalFiles/BsPPCorrYieldMult.root";
	TFile * FileBs = new TFile(InfileBs.Data());
	TH1D * BsCross = (TH1D *) FileBs->Get(FileName.Data());
	BsCross->SetMarkerStyle(20);
	BsCross->SetMarkerSize(1);
	BsCross->SetMarkerColor(1);
	BsCross->SetLineColor(1);

	float BsXsecPPY[NBins];
	float BsXsecPPX[NBins] = {7.5,20,27.5,32.5,37.5,45,57.5,72.5,90,115};
	float BsXSecPPYErrUp[NBins];
	float BsXSecPPYErrDown[NBins];
	float BsXSecPPYErrUpPercent[NBins];
	float BsXSecPPYErrDownPercent[NBins];
	float BsXSecPPXErrUp[NBins] = {7.5,5,2.5,2.5,2.5,5,7.5,7.5,10,15};
	float BsXSecPPXErrDown[NBins] = {7.5,5,2.5,2.5,2.5,5,7.5,7.5,10,15};


	for(int i = 0; i < NBins; i++){
		BsXsecPPY[i] = BsCross->GetBinContent(i+1);
		BsXSecPPYErrUp[i] = BsCross->GetBinError(i+1);
		BsXSecPPYErrDown[i] = BsCross->GetBinError(i+1);
		BsXSecPPYErrUpPercent[i] = BsXSecPPYErrUp[i]/BPXsecPPY[i];
		BsXSecPPYErrDownPercent[i] = BsXSecPPYErrDown[i]/BPXsecPPY[i];
		cout << "BsXsecPPY = " << BsXsecPPY[i] << endl;
	}

	TCanvas * c2 = new TCanvas("c","c",600,600);
	c2->cd();
	TH2D * HisEmpty2 = new TH2D("HisEmpty2","",100,0,130,100,20.0,50000);
	HisEmpty2->GetXaxis()->SetTitle("Event Multiplicity");
	HisEmpty2->GetYaxis()->SetTitle("#Delta #sigma (pb)");
	HisEmpty2->GetXaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->SetTitleOffset(1.4);
	HisEmpty2->Draw();

	TGraphAsymmErrors *BsPPCrossGraph = new TGraphAsymmErrors(NBins, BsXsecPPX, BsXsecPPY,BsXSecPPXErrDown, BsXSecPPXErrUp,BsXSecPPYErrDown,BsXSecPPYErrUp);
	BsPPCrossGraph->SetLineColor(kGreen+2);
	BsPPCrossGraph->SetMarkerStyle(20);
	BsPPCrossGraph->SetMarkerSize(1);
	BsPPCrossGraph->SetMarkerColor(kGreen+2);
	BsPPCrossGraph->Draw("ep");

	c2->SaveAs(Form("Plots/%s/BsPPCrossMult.png",Method.Data()));
	c2->SetLogy();
	c2->SaveAs(Form("Plots/%s/BsPPCrossMultLog.png",Method.Data()));

	TCanvas * c3 = new TCanvas("c","c",600,600);
	c3->cd();
	
	TH2D * HisEmpty3 = new TH2D("HisEmpty3","",100,0,130,100,20.0,600000);
	HisEmpty3->GetXaxis()->SetTitle("Event Multiplicity");
	HisEmpty3->GetYaxis()->SetTitle("#Delta #sigma (pb)");
	HisEmpty3->GetXaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->SetTitleOffset(1.4);
	HisEmpty3->Draw();

	BsPPCrossGraph->Draw("ep");
	BPPPCrossGraph->Draw("epSAME");

	TLegend* leg = new TLegend(0.65,0.65,0.95,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(BsPPCrossGraph,"B^{0}_{s}","PL");
	leg->AddEntry(BPPPCrossGraph,"B^{+}","PL");
	leg->Draw("same");

	c3->SaveAs(Form("Plots/%s/BsBPMult.png",Method.Data()));
	c3->SetLogy();
	c3->SaveAs(Form("Plots/%s/BsBPMultLog.png",Method.Data()));


	//Bs B+ Ratio PP//

	TH1D * BsBPRatio = (TH1D *) BsCross->Clone("BsBPRatio");
	BsBPRatio->Sumw2();
	BPCross->Sumw2();

	BsBPRatio->Divide(BPCross);

	float BsBPPPY[NBins];
	float BsBPPPX[NBins] = {7.5,20,27.5,32.5,37.5,45,57.5,72.5,90,115};

	float BsBPPPYErrUp[NBins];
	float BsBPPPYErrDown[NBins];


	float BsBPPPYErrUpPercent[NBins];
	float BsBPPPYErrDownPercent[NBins];


	float BsBPPPXErrUp[NBins] = {7.5,5,2.5,2.5,2.5,5,7.5,7.5,10,15};
	float BsBPPPXErrDown[NBins] = {7.5,5,2.5,2.5,2.5,5,7.5,7.5,10,15};


	for(int i = 0; i < NBins; i++){

		BsBPPPY[i] = BsBPRatio->GetBinContent(i+1);
		BsBPPPYErrUp[i] = BsBPRatio->GetBinError(i+1);
		BsBPPPYErrDown[i] = BsBPRatio->GetBinError(i+1);
		BsBPPPYErrUpPercent[i] = BsBPPPYErrUp[i]/BsBPPPY[i];
		BsBPPPYErrDownPercent[i] = BsBPPPYErrDown[i]/BsBPPPY[i];

		cout << "BsBPPPY = " << BsBPPPY[i] << endl;

	}


	TGraphAsymmErrors *BsBPRatioGraph = new TGraphAsymmErrors(NBins, BsBPPPX, BsBPPPY,BsBPPPXErrDown, BsBPPPXErrUp,BsBPPPYErrDown,BsBPPPYErrUp);




	
	TCanvas * c4 = new TCanvas("c","c",600,600);
	c4->cd();

	TH2D * HisEmpty4 = new TH2D("HisEmpty4","",100,0,130,100,0.0,0.4);
	HisEmpty4->GetXaxis()->SetTitle("Event Multiplicity");
	HisEmpty4->GetYaxis()->SetTitle("B^{0}_{s}/B^{+}");
	HisEmpty4->GetXaxis()->CenterTitle();
	HisEmpty4->GetYaxis()->CenterTitle();
	HisEmpty4->GetYaxis()->SetTitleOffset(1.4);
	HisEmpty4->Draw();


	BsBPRatioGraph->SetLineColor(kRed+2);
	BsBPRatioGraph->SetMarkerStyle(41);
	BsBPRatioGraph->SetMarkerSize(1);
	BsBPRatioGraph->SetMarkerColor(kRed+2);

	BsBPRatioGraph->Draw("ep");


	c4->SaveAs(Form("Plots/%s/BsBPRatio.png",Method.Data()));
	



}

