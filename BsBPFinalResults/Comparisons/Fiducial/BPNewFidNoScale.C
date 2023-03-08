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
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"

//#include "his.h"
using namespace std;
using std::cout;
using std::endl;

// scaled lower data pt bins to full rapidity
constexpr bool scaleMC = false;
// swap the lower fonll bins with full rapidity
constexpr bool fidFONLL = true;

void BPNewFidNoScale(){

	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);

	c->cd();
	c->SetLeftMargin(0.16);

	TString InfileBP = "../../../BP/EffAna/FinalFiles/BPPPCorrYieldPT.root";
	TFile * FileBP = new TFile(InfileBP.Data());
	TH1D * BPCross = (TH1D *) FileBP->Get("CorrDiffHisBin");
	BPCross->SetMarkerStyle(20);
	BPCross->SetMarkerSize(1);
	BPCross->SetMarkerColor(1);
	BPCross->SetLineColor(1);

	//B+ PbPb//
	
	const int NBins = 7;
	float BPXsecPPY[NBins];
	float BPXsecPPX[NBins] = {6,8.5,12.5,17.5,25,40,55};
	float BPXSecPPYErrUp[NBins];
	float BPXSecPPYErrDown[NBins];
	float BPXSecPPYErrUpPercent[NBins];
	float BPXSecPPYErrDownPercent[NBins];


	for(int i = 0; i < NBins; i++){

		BPXsecPPY[i] = BPCross->GetBinContent(i+1);
		BPXSecPPYErrUp[i] = BPCross->GetBinError(i+1);
		BPXSecPPYErrDown[i] = BPCross->GetBinError(i+1);
		BPXSecPPYErrUpPercent[i] = BPXSecPPYErrUp[i]/BPXsecPPY[i];
		BPXSecPPYErrDownPercent[i] = BPXSecPPYErrDown[i]/BPXsecPPY[i];
		cout << "BPXsecPPY[i] = " << BPXsecPPY[i]  << "   Stat[i]  = " << BPXSecPPYErrUpPercent[i] << endl;

	}


	float BPXsecPPY2D[NBins];
	float BPXsecPPY2DScaled[NBins];
	float BPXSecPPY2DErrUp[NBins];
	float BPXSecPPY2DErrDown[NBins];
	TH1D * BPCross2D = (TH1D *) FileBP->Get("hPtSigma");
	BPCross2D->SetMarkerStyle(20);
	BPCross2D->SetMarkerSize(1);
	BPCross2D->SetMarkerColor(1);
	BPCross2D->SetLineColor(1);






	for(int i = 0; i < NBins; i++){

		BPXsecPPY2D[i] = BPCross2D->GetBinContent(i+1);
		BPXSecPPY2DErrUp[i] = BPCross2D->GetBinError(i+1);
		BPXSecPPY2DErrDown[i] = BPCross2D->GetBinError(i+1);
		BPXsecPPY2DScaled[i] = BPCross2D->GetBinContent(i+1);
	}

  std::vector<double> factor = {2.532, 3.927};
  if (scaleMC) {
    for (auto i = 0; i < factor.size(); ++i) {
      BPXsecPPY2D[i] = BPXsecPPY2D[i] * factor[i];
    }
  }




	float BPXSecPPXErrUp[NBins] = {1,1.5,2.5,2.5,5,10,5};
	float BPXSecPPXErrDown[NBins] = {1,1.5,2.5,2.5,5,10,5};
	float BPXsecPbPbY[NBins] = {4.82132e+06/11.1,311668,270167,64384.4,208537/11.1,28700.6/11.1,7000.73/11.1};
	float BPXsecPbPbX[NBins] = {6,8.73,12.4,17.2,25,40,55};
	float BPXSecPbPbXErrUp[NBins] = {1,1.27,2.6,2.8,5,10,5};
	float BPXSecPbPbXErrDown[NBins] = {1,1.23,2.4,2.2,5,10,5};
	float BPXSecPbPbYErrUpPercent[NBins] = {0.278198,0.159,0.041,0.0654,0.0690334,0.104543,0.24575};
	float BPXSecPbPbYErrDownPercent[NBins] = {0.278198,0.145,0.0795,0.065,0.0690334,0.104543,0.24575};
	float BPXSecPbPbYErrUp[NBins];
	float BPXSecPbPbYErrDown[NBins];

	for(int i = 0; i < NBins; i++){
		BPXSecPbPbYErrUp[i] = BPXSecPbPbYErrUpPercent[i] * BPXsecPbPbY[i];
		BPXSecPbPbYErrDown[i] = BPXSecPbPbYErrDownPercent[i] * BPXsecPbPbY[i];
	}




	//Syst Add Up PP//



	float BPXSecPPYSystUp[NBins];
	float BPXSecPPYSystDown[NBins];
	float BPTrackingSyst[NBins] = {0.05,0.05,0.05,0.05,0.05,0.05,0.05};
	float BPMDDataSyst[NBins] = {1.92761, 1.43481, 0.301207, 0.131291, 0, 0, 0};
	// float BPPDFSyst[NBins] = {0.0212,0.0122,0.0091,0.0411,0.0120,0.0150,0.0089};
	float BPPDFSyst[NBins] = {0};
	// float BPPtShapeSyst[NBins] = {0.0096,0.00814,0.000408,0.000617,0.000833,0.00138,0.000640};
	float BPPtShapeSyst[NBins] = {0.00538802, 0.00302489, 0.00115669, 0.00031239, 0.00069791, 0.00274824, 0.00864261};
	// float BPTnPSystDown[NBins] = {0.00454,0.00476,0.00382,0.00354,0.00429,0.00581,0.00613};
	float BPTnPSystDown[NBins] = { 0.46582, 0.44548, 0.37401, 0.36262, 0.4261, 0.62824, 0.60984};
	// float BPTnPSystUp[NBins] = {0.00454,0.00476,0.00382,0.00354,0.00429,0.00581,0.00613};
	float BPTnPSystUp[NBins] = { 0.46582, 0.44548, 0.37401, 0.36262, 0.4261, 0.62824, 0.60984};

	float BPTotalSystDown[NBins];
	float BPTotalSystUp[NBins];

	for(int i = 0; i < NBins; i++){

	
		BPTotalSystDown[i] = TMath::Sqrt(BPTrackingSyst[i] * BPTrackingSyst[i] + BPMDDataSyst[i] * BPMDDataSyst[i] + BPPDFSyst[i] * BPPDFSyst[i] + BPPtShapeSyst[i] * BPPtShapeSyst[i] + BPTnPSystDown[i] * BPTnPSystDown[i]) / 100;
		BPTotalSystUp[i] = TMath::Sqrt(BPTrackingSyst[i] * BPTrackingSyst[i] + BPMDDataSyst[i] * BPMDDataSyst[i] + BPPDFSyst[i] * BPPDFSyst[i] + BPPtShapeSyst[i] * BPPtShapeSyst[i] + BPTnPSystUp[i] * BPTnPSystUp[i]) / 100;


	}
  // global uncertainty probably means cross section and BR
  // no, it means uncertanty from branching ratio and luminosity
  std::vector<float> globUncert(NBins);
	for(int i = 0; i < NBins; i++){

		BPXSecPPYSystUp[i] = BPXsecPPY[i] * ( BPTotalSystUp[i]);
		BPXSecPPYSystDown[i] = BPXsecPPY[i] * (BPTotalSystDown[i] );
    globUncert[i] = sqrt(std::pow(BPXSecPPYErrUpPercent[i], 2) + std::pow(BPTotalSystDown[i], 2));

		cout << "i = " << i << "     BP syst[i] = " << BPTotalSystDown[i] <<
		"    glob: " << globUncert[i] << "\n";

	}


	//PbPb

	float BPXSecPbPbYSystUpPercent[NBins] = {0.3577,0.1404,0.1714,0.0775,0.0858,0.0715,0.1253};
	float BPXSecPbPbYSystDownPercent[NBins] = {0.3210,0.1359,0.1705,0.0761,0.0843,0.0699,0.1220};
	float BPXSecPbPbYSystUp[NBins];
	float BPXSecPbPbYSystDown[NBins];

	for(int i = 0; i < NBins; i++){

		BPXSecPbPbYSystDown[i] = (BPXSecPbPbYSystDownPercent[i]) * BPXsecPbPbY[i];
		BPXSecPbPbYSystUp[i] = (BPXSecPbPbYSystUpPercent[i]) * BPXsecPbPbY[i];

	}





	//Setup the Syst



	TH2D * HisEmpty = new TH2D("HisEmpty","",100,5,60,100,100.0,2000000);
	HisEmpty->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty->GetYaxis()->SetTitle("d#sigma/dp_{T} (pb c/GeV)");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	HisEmpty->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty->GetXaxis()->SetTitleOffset(1.3);		
	HisEmpty->Draw();

	TGraphAsymmErrors *BPPPCrossGraph = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY,BPXSecPPXErrDown, BPXSecPPXErrUp,BPXSecPPYErrDown,BPXSecPPYErrUp);
	TGraphAsymmErrors *BPPPCrossGraphSyst  = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY, BPXSecPPXErrDown, BPXSecPPXErrUp, BPXSecPPYSystDown,BPXSecPPYSystUp);
	TGraphAsymmErrors *BPPPCrossGraph2D = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY2D, BPXSecPPXErrDown, BPXSecPPXErrUp,BPXSecPPY2DErrDown,BPXSecPPY2DErrUp);
	// TGraphAsymmErrors *BPPPCrossGraph2DScaled = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY2DScaled, BPXSecPPXErrDown, BPXSecPPXErrUp,BPXSecPPY2DErrDown,BPXSecPPY2DErrUp);
	TGraphAsymmErrors *BPPbPbCrossGraph = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY,BPXSecPbPbXErrDown, BPXSecPbPbXErrUp,BPXSecPbPbYErrDown,BPXSecPbPbYErrUp);
  	TGraphAsymmErrors *BPPbPbCrossGraphSyst    = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY, BPXSecPbPbXErrDown, BPXSecPbPbXErrUp, BPXSecPbPbYSystDown,BPXSecPbPbYSystUp);
 
	BPPbPbCrossGraph->SetLineColor(kGreen+2);
//	BPPbPbCrossGraph->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraph->SetMarkerStyle(20);
	BPPbPbCrossGraph->SetMarkerSize(1);
	BPPbPbCrossGraph->SetMarkerColor(kGreen+2);
	BPPPCrossGraph->SetLineColor(kBlue+2);
//	BPPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraph->SetMarkerStyle(21);
	BPPPCrossGraph->SetMarkerSize(1);
	BPPPCrossGraph->SetMarkerColor(kBlue+2);
	BPPPCrossGraphSyst->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraphSyst->SetLineColor(kBlue-9);
	BPPbPbCrossGraphSyst->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraphSyst->SetLineColor(kGreen-9);
	BPPPCrossGraph->Draw("ep");	
	BPPPCrossGraphSyst->Draw("5same");	


	// c->SaveAs("Plots/BP/BPCrossONLY.png");
	// c->SaveAs("Plots/BP/BPCrossONLY.pdf");
	// c->SetLogy();
	// c->SaveAs("Plots/BP/BPCrossONLYLog.png");
	// c->SaveAs("Plots/BP/BPCrossONLYLog.pdf");


	TCanvas * c2New = new TCanvas("c2New","c2New",600,600);
	c2New->cd();
	HisEmpty->Draw();

	BPPbPbCrossGraph->Draw("epsame");	
	BPPPCrossGraph->Draw("epsame");	

	TLegend* leg = new TLegend(0.30,0.60,0.60,0.80,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(BPPbPbCrossGraph,"2018 PbPb 5.02 TeV","PL");
	leg->AddEntry(BPPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg->Draw("same");
	BPPPCrossGraphSyst->Draw("5same");	
	BPPbPbCrossGraphSyst->Draw("5same");	
	
	// c2New->SaveAs("Plots/BP/BPPbPbPPCross.png");
	// c2New->SetLogy();
	// c2New->SaveAs("Plots/BP/BPPbPbPPCrossLog.png");
	//2015 Reference//


	TCanvas * cRatio = new TCanvas("cRatio","cRatio",800,1200);
    TPad * MyPad1;
	MyPad1 = new TPad("MyPad1","",0,0.5,1,1.0);
	MyPad1->Draw();
   
	TPad * MyPad2;
	MyPad2 = new TPad("MyPad2","",0,0.25,1,0.50);
	// MyPad2->Draw();

	TPad * MyPad3;

	// MyPad3 = new TPad("MyPad3","",0,0.00,1,0.25);
	MyPad3 = new TPad("MyPad2","",0,0.25,1,0.50);
	MyPad3->Draw();


	MyPad1->cd();



	TH2D * HisEmpty2 = new TH2D("HisEmpty2","",100,5,60,100,100.0,30000000);
	HisEmpty2->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty2->GetYaxis()->SetTitle("d#sigma/d p_{T} (pb c/GeV)");
	HisEmpty2->GetXaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->CenterTitle();
	HisEmpty2->SetTitle("B^{+} Cross Section With Fiducial Region");
	HisEmpty2->GetYaxis()->SetTitleOffset(1.2);
	HisEmpty2->Draw();

	// BPPPCrossGraph->Draw("ep");
	// BPPPCrossGraphSyst->Draw("5same");

	const int NBins2015 = 5;
	float BPXsecPPX2015[NBins2015] = {8.5,12.5,17.5,25,40};
	float BPXSecPPXErrDown2015[NBins2015] = {1.5,2.5,2.5,5,10};
	float BPXSecPPXErrUp2015[NBins2015] = {1.5,2.5,2.5,5,10};
	float BPXsecPPY2015[NBins2015] = {2610000,744000,197000,46500,5300};
	float BPXSecPPYErrDown2015[NBins2015] = {170000,29000,9000,2400,500};
	float BPXSecPPYErrUp2015[NBins2015] = {170000,29000,9000,2400,500};
	float BPXSecPPYSystDown2015[NBins2015] = {230000,59000,15000,3500,400};
	float BPXSecPPYSystUp2015[NBins2015] = {230000,59000,15000,3500,400};

	TGraphAsymmErrors *BPPPCrossGraph2015 = new TGraphAsymmErrors(NBins2015, BPXsecPPX2015, BPXsecPPY2015,BPXSecPPXErrDown2015, BPXSecPPXErrUp2015,BPXSecPPYErrDown2015,BPXSecPPYErrUp2015);
	TGraphAsymmErrors *BPPPCrossGraph2015Syst = new TGraphAsymmErrors(NBins2015, BPXsecPPX2015, BPXsecPPY2015,BPXSecPPXErrDown2015, BPXSecPPXErrUp2015,BPXSecPPYSystDown2015,BPXSecPPYSystUp2015);
	BPPPCrossGraph2015Syst->SetFillColorAlpha(kGreen-9+2,0.5);
	BPPPCrossGraph2015Syst->SetLineColor(kGreen-9+2);
	BPPPCrossGraph2015->SetLineColor(kGreen+2);
	BPPPCrossGraph2015->SetMarkerStyle(33);
	BPPPCrossGraph2015->SetMarkerSize(1);
	BPPPCrossGraph2015->SetMarkerColor(kGreen+2);
	// BPPPCrossGraph2015->Draw("epSAME");
	// BPPPCrossGraph2015Syst->Draw("5same");	

	BPPPCrossGraph2D->SetLineColor(kOrange+1);
	BPPPCrossGraph2D->SetMarkerStyle(34);
	BPPPCrossGraph2D->SetMarkerSize(1);
	BPPPCrossGraph2D->SetMarkerColor(kOrange+1);

	TFile * finFONLL = new TFile("FONLLs/fonllOutput_pp_Bplus_5p03TeV_y2p4.root");
	finFONLL->cd();
	TGraphAsymmErrors *BPFONLL = (TGraphAsymmErrors*) finFONLL->Get("gaeSigmaBplus");
	BPFONLL->SetLineColor(kRed+2);
	BPFONLL->SetMarkerStyle(20);
	BPFONLL->SetMarkerSize(1);
	BPFONLL->SetMarkerColor(kRed+2);
	BPFONLL->Draw("epSAME");

	TFile * finFONLL2 = new TFile("FONLLs/fonllOutput_pp_Bplus_5p03TeV_yFid.root");
	finFONLL2->cd();
	TGraphAsymmErrors *BPFONLL2 = (TGraphAsymmErrors*) finFONLL2->Get("gaeSigmaBplus");
	BPFONLL2->SetLineColor(kRed+2);
	BPFONLL2->SetMarkerStyle(20);
	BPFONLL2->SetMarkerSize(1);
	BPFONLL2->SetMarkerColor(kRed+2);
	// BPFONLL2->Draw("epSAME");

	double XTempChange;
	double YTempChange;
	double YErrLowTemp;
	double YErrHighTemp;

  if (fidFONLL) {
    for(int i = 0; i < 2; i ++){
      BPFONLL2->GetPoint(i,XTempChange,YTempChange);
      YErrLowTemp = BPFONLL2->GetErrorYlow(i);
      YErrHighTemp = BPFONLL2->GetErrorYhigh(i);
      BPFONLL->SetPoint(i,XTempChange,YTempChange);
      BPFONLL->SetPointEYhigh(i,YErrHighTemp);
      BPFONLL->SetPointEYlow(i,YErrLowTemp);
    }
  }

	BPPPCrossGraph2D->Draw("epSAME");
	TLegend* leg3 = new TLegend(0.37,0.50,0.70,0.80,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.040);
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->SetLineWidth(3);
	// leg3->AddEntry(BPPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg3->AddEntry(BPPPCrossGraph2D,"2017 pp 5.02 TeV (2D Map)","PL");	
	// leg3->AddEntry(BPPPCrossGraph2015,"2015 pp 5.02 TeV","PL");
	leg3->AddEntry(BPFONLL,"FONLL Calculations","PL");
	leg3->Draw("same");

	MyPad1->Update();

	//Ratio

	float Ratio1Y[NBins2015];
	float Ratio1YErr[NBins2015];
	float Ratio2Y[NBins2015];
	float Ratio2YErr[NBins2015];

	for(int i = 1; i < NBins2015; i++){

		Ratio1Y[i] = BPXsecPPY[i+1]/BPXsecPPY2015[i];
		Ratio1YErr[i] = Ratio1Y[i] * TMath::Sqrt(BPXSecPPYErrDown[i+1]/BPXsecPPY[i+1] * BPXSecPPYErrDown[i+1]/BPXsecPPY[i+1] + BPXSecPPYErrDown2015[i]/BPXsecPPY2015[i] * BPXSecPPYErrDown2015[i]/BPXsecPPY2015[i]);

	//	cout << "BPXSecPPYErrDown2015[i]/BPXsecPPY2015[i] = " << BPXSecPPYErrDown[i+1]/BPXsecPPY[i+1]  << endl;
	//	cout << "sqrt(BPXSecPPYErrDown[i+1]/BPXsecPPY[i+1] * BPXSecPPYErrDown[i+1]/BPXsecPPY[i+1] + BPXSecPPYErrDown2015[i]/BPXsecPPY2015[i] * BPXSecPPYErrDown2015[i]/BPXsecPPY2015[i]) = " << sqrt(BPXSecPPYErrDown[i+1]/BPXsecPPY[i+1] * BPXSecPPYErrDown[i+1]/BPXsecPPY[i+1] + BPXSecPPYErrDown2015[i]/BPXsecPPY2015[i] * BPXSecPPYErrDown2015[i]/BPXsecPPY2015[i]) << endl;

		Ratio2Y[i] = BPXsecPPY2D[i+1]/BPXsecPPY2015[i];
		Ratio2YErr[i] = Ratio2Y[i] *TMath::Sqrt(BPXSecPPY2DErrDown[i+1]/BPXsecPPY2D[i+1] * BPXSecPPY2DErrDown[i+1]/BPXsecPPY2D[i+1] + BPXSecPPYErrDown2015[i]/BPXsecPPY2015[i] * BPXSecPPYErrDown2015[i]/BPXsecPPY2015[i] );
	}

	Ratio2YErr[0] = 0.00001;
	Ratio2Y[0] = -0.1;
	Ratio1YErr[0] = 0.00001;
	Ratio1Y[0] = -0.1;

	cRatio->cd();



	MyPad2->cd();



	TH2D * HisEmpty3 = new TH2D("HisEmpty3","",100,5,60,100,0,2);
	HisEmpty3->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty3->GetYaxis()->SetTitle("2017/2015 Data");
	HisEmpty3->GetXaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->SetTitleOffset(1.2);
	HisEmpty3->Draw();


	TGraphAsymmErrors *Ratio1 = new TGraphAsymmErrors(NBins2015, BPXsecPPX2015, Ratio1Y ,BPXSecPPXErrDown2015, BPXSecPPXErrUp2015,Ratio1YErr,Ratio1YErr);
	
	TGraphAsymmErrors *Ratio2 = new TGraphAsymmErrors(NBins2015, BPXsecPPX2015, Ratio2Y,BPXSecPPXErrDown2015, BPXSecPPXErrUp2015,Ratio2YErr,Ratio2YErr);

	// Ratio1->SetLineColor(kBlue+2);
	// Ratio1->SetMarkerStyle(21);
	// Ratio1->SetMarkerSize(1);
	// Ratio1->SetMarkerColor(kBlue+2);
	// Ratio1->Draw("epSAME");




	Ratio2->SetLineColor(kOrange+1);
	Ratio2->SetMarkerStyle(34);
	Ratio2->SetMarkerSize(1);
	Ratio2->SetMarkerColor(kOrange+1);
	// Ratio2->Draw("epSAME");


	TLine * Unity2 = new TLine(5,1,60,1);
	Unity2->SetLineWidth(2);
	Unity2->SetLineStyle(2);
	Unity2->SetLineColor(1);
	// Unity2->Draw("SAME");




	// MyPad2->Update();


	//FONLL




	float Ratio3Y[NBins];
	float Ratio3YErr[NBins];
	
	float Ratio4Y[NBins];
	float Ratio4YErr[NBins];

	float FONLLY[NBins];
	float FONLLYErr[NBins];


	double XTempFONLL;
	double YTempFONLL;

	for(int i = 0; i < NBins; i++){



		
		BPFONLL->GetPoint(i,XTempFONLL,YTempFONLL);
		FONLLY[i] = YTempFONLL;
		FONLLYErr[i] = BPFONLL->GetErrorYhigh (i);


		cout << "FONLLY[i]  = " << FONLLY[i] << endl;

		Ratio3Y[i] = BPXsecPPY[i]/FONLLY[i];
		Ratio3YErr[i] = Ratio3Y[i] * sqrt(BPXSecPPYErrDown[i]/BPXsecPPY[i] * BPXSecPPYErrDown[i]/BPXsecPPY[i] + FONLLYErr[i]/FONLLY[i] * FONLLYErr[i]/FONLLY[i]);

		Ratio4Y[i] = BPXsecPPY2D[i]/FONLLY[i];
		Ratio4YErr[i] = Ratio4Y[i] *TMath::Sqrt(BPXSecPPY2DErrDown[i]/BPXsecPPY2D[i] * BPXSecPPY2DErrDown[i]/BPXsecPPY2D[i] + FONLLYErr[i]/FONLLY[i] * FONLLYErr[i]/FONLLY[i] );
	

		//cout << "Ratio3YErr[i]  = " << Ratio3YErr[i]  << "   Ratio4YErr[i] = " << Ratio4YErr[i] << endl;

	}

	cRatio->cd();



	MyPad3->cd();



	TH2D * HisEmpty4 = new TH2D("HisEmpty4","",100,5,60,100,0,2);
	HisEmpty4->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty4->GetYaxis()->SetTitle("2017 Data/FONLL");
	HisEmpty4->GetXaxis()->CenterTitle();
	HisEmpty4->GetYaxis()->CenterTitle();
	HisEmpty4->GetYaxis()->SetTitleOffset(1.2);
	HisEmpty4->Draw();


	TGraphAsymmErrors *Ratio3 = new TGraphAsymmErrors(NBins, BPXsecPPX, Ratio3Y ,BPXSecPPXErrDown, BPXSecPPXErrUp,Ratio3YErr,Ratio3YErr);

	TGraphAsymmErrors *Ratio4 = new TGraphAsymmErrors(NBins, BPXsecPPX, Ratio4Y,BPXSecPPXErrDown, BPXSecPPXErrUp,Ratio4YErr,Ratio4YErr);

	// Ratio3->SetLineColor(kBlue+2);
	// Ratio3->SetMarkerStyle(21);
	// Ratio3->SetMarkerSize(1);
	// Ratio3->SetMarkerColor(kBlue+2);
	// Ratio3->Draw("epSAME");




	Ratio4->SetLineColor(kOrange+1);
	Ratio4->SetMarkerStyle(34);
	Ratio4->SetMarkerSize(1);
	Ratio4->SetMarkerColor(kOrange+1);
	Ratio4->Draw("epSAME");



	Unity2->Draw("SAME");




	MyPad2->Update();








	// cRatio->SaveAs("Plots/BP/BPCrossComp.png");
//	cRatio->SetLogy();

	//MyPad1->SetLogy();
	
	MyPad1->SetLogy();
	MyPad1->Update();

	cRatio->SaveAs("Plots/BP/BPCrossCompLogNoScale.png");

	cRatio->SaveAs("Plots/BP/BPCrossCompLogNoScale.pdf");




	//FONLL







  std::vector<int> ptbins = {5, 7, 10, 15, 20, 30, 50, 60};
  std::vector<float> abscissae = {6.0, 8.75, 12.5, 17.5, 25, 40, 55};

gSystem->mkdir("../../../MakeFinalPlots/NominalPlots/CrossSection/dataSource/" ,true );
  string outFile = "../../../MakeFinalPlots/NominalPlots/CrossSection/dataSource/corryield_pt_Bp_New.txt";
  ofstream out;
  out.open(outFile);
  out << "ptmin\t ptmax\t central val\t stat\t systUp\t systDown\t glbUp\t glbDown" << endl;
  for (auto i = 0; i < NBins; ++i ) {
    out << ptbins[i] << "\t" << ptbins[i + 1] << "\t" << BPXsecPPY[i] << "\t" <<
      BPXSecPPYErrUpPercent[i] << "\t" << BPXSecPPYErrDownPercent[i] << "\t" <<
      BPTotalSystDown[i] << "\t" << BPTotalSystDown[i] << "\t" << globUncert[i] <<
      "\t" << globUncert[i] << "\t" << abscissae[i] << "\n";
  }
  out.close();

}
