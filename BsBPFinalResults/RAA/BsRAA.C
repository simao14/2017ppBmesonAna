#include "TROOT.h"
#include "TStyle.h"
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
#include <iomanip>
#include <fstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TLine.h"
//#include "tnp_weight_lowptPbPb.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void BsRAA(){




	gStyle->SetOptStat(0);

	TString InfileBs = "../../Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root";
	//TString InfileBs = "BsPPCorrYieldPT.root";

	TFile * FileBs = new TFile(InfileBs.Data());

	TH1D * BsCross = (TH1D *) FileBs->Get("hPtSigma");
	BsCross->SetMarkerStyle(20);
	BsCross->SetMarkerSize(1);
	BsCross->SetMarkerColor(1);
	BsCross->SetLineColor(1);
	BsCross->Draw("ep");

	//B+ PbPb//
	
	const int NBins = 4;
	float BsXsecPPY[NBins];
	float BsXsecPPX[NBins] = {8.5,12.5,17.5,35};
	float BsXSecPPYErrUp[NBins];
	float BsXSecPPYErrDown[NBins];
	float BsXSecPPYErrUpPercent[NBins];
	float BsXSecPPYErrDownPercent[NBins];

	for(int i = 0; i < NBins; i++){
		BsXsecPPY[i] = BsCross->GetBinContent(i+1);
		BsXSecPPYErrUp[i] = BsCross->GetBinError(i+1);
		BsXSecPPYErrDown[i] = BsCross->GetBinError(i+1);
		BsXSecPPYErrUpPercent[i] = BsXSecPPYErrUp[i]/BsXsecPPY[i];
		BsXSecPPYErrDownPercent[i] = BsXSecPPYErrDown[i]/BsXsecPPY[i];
	}

	float BsXSecPPXErrUp[NBins] = {1.5,2.5,2.5,15};
	float BsXSecPPXErrDown[NBins] = {1.5,2.5,2.5,15};
	float BsXsecPbPbY[NBins] = {160432,75523.7,25354.5,2272.18};
	float BsXsecPbPbX[NBins] = {8.75,12.6,17.4,27.3};
	float BsXSecPbPbXErrUp[NBins] = {1.25,2.4,2.4,22.7};
	float BsXSecPbPbXErrDown[NBins] = {1.25,2.6,2.6,7.3};
	float BsXSecPbPbYErrUpPercent[NBins] = {0.513,0.224,0.216,0.216};
	float BsXSecPbPbYErrDownPercent[NBins] = {0.483,0.256,0.207,0.163};
	float BsXSecPbPbYErrUp[NBins];
	float BsXSecPbPbYErrDown[NBins];

	for(int i = 0; i < NBins; i++){
		BsXSecPbPbYErrUp[i] = BsXSecPbPbYErrUpPercent[i] * BsXsecPbPbY[i];
		BsXSecPbPbYErrDown[i] = BsXSecPbPbYErrDownPercent[i] * BsXsecPbPbY[i];
	}

	//Syst//

	float BsXSecPPYSystUp[NBins];
	float BsXSecPPYSystDown[NBins];

  TString errorFile = "../../2DMapSyst/OutFiles/BsError2D.root";
  TFile fError(errorFile);

	TH1D * TnPSyst = (TH1D *) fError.Get("TnPSyst");
	TH1D * BptSyst = (TH1D *) fError.Get("BptSyst");
	TH1D * MCDataSyst = (TH1D *) fError.Get("MCDataSyst");
  	if (!MCDataSyst) MCDataSyst = (TH1D *) fError.Get("BDTSyst");

  TString pdfErrorFile = "../../bs_pdf.root";
  TFile fPdfError(pdfErrorFile);
  TGraph* pdfSyst = (TGraph *) fPdfError.Get("bs_error");

  TString trackSelErrorFile = "../../syst_track_sel.root";
  TFile fTrackSelError(trackSelErrorFile);
  TGraph* trackSelSyst = (TGraph *) fTrackSelError.Get("bs_track_sel_error");

  // percent error
	float BsTrackingSyst[NBins] = {[0 ... NBins - 1] = 10};
	float BsMCDataSyst[NBins];
	float BsPtShapeSyst[NBins];
	float BsPDFSyst[NBins];
	float BsTrackSelSyst[NBins];
	float BsTnPSystDown[NBins];
	float BsTnPSystUp[NBins];

  for (auto ibin = 0; ibin < NBins; ++ibin) {
    BsMCDataSyst[ibin] = MCDataSyst->GetBinContent(ibin + 1);
    BsPtShapeSyst[ibin] = BptSyst->GetBinContent(ibin + 1);
    BsTnPSystDown[ibin] = TnPSyst->GetBinContent(ibin + 1);
    // TnP systematics are symmetric in the binned pT case
    BsTnPSystUp[ibin] = BsTnPSystDown[ibin];
    BsPDFSyst[ibin] = pdfSyst->GetY()[ibin];
    BsTrackSelSyst[ibin] = trackSelSyst->GetY()[ibin];
  }

  // RMS of all the errors
	float BsTotalSystDownRatio[NBins];
	float BsTotalSystUpRatio[NBins];

	for(int i = 0; i < NBins; i++){
		BsTotalSystDownRatio[i] = TMath::Sqrt(TMath::Power(BsTrackingSyst[i], 2) + TMath::Power(BsMCDataSyst[i], 2) +
                                          TMath::Power(BsPDFSyst[i], 2) + TMath::Power(BsTrackSelSyst[i], 2) +
                                          TMath::Power(BsPtShapeSyst[i], 2) + TMath::Power(BsTnPSystDown[i], 2)) / 100;
		BsTotalSystUpRatio[i] = BsTotalSystDownRatio[i];
	}

	for(int i = 0; i < NBins; i++){
		BsXSecPPYSystUp[i] = BsXsecPPY[i] * BsTotalSystUpRatio[i];
		BsXSecPPYSystDown[i] = BsXsecPPY[i] * BsTotalSystDownRatio[i];
	}

	//PbPb
	float BsXSecPbPbYSystUpPercent[NBins] = {0.4564,0.1482,0.1218,0.1647};
	float BsXSecPbPbYSystDownPercent[NBins] = {0.4564,0.1454,0.1210,0.1640};
	float BsXSecPbPbYSystUp[NBins];
	float BsXSecPbPbYSystDown[NBins];

	for(int i = 0; i < NBins; i++){
		BsXSecPbPbYSystDown[i] = (BsXSecPbPbYSystDownPercent[i]) * BsXsecPbPbY[i];
		BsXSecPbPbYSystUp[i] = (BsXSecPbPbYSystUpPercent[i]) * BsXsecPbPbY[i];
	}

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetLeftMargin(0.16);

	TH2D * HisEmpty = new TH2D("HisEmpty","",100,7,50,100,200.0,350000);
	HisEmpty->GetXaxis()->SetTitle("B^{0}_{s} p_{T} [GeV/c]");
	HisEmpty->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	HisEmpty->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty->Draw();

	TGraphAsymmErrors *BsPPCrossGraph = new TGraphAsymmErrors(NBins, BsXsecPPX, BsXsecPPY,BsXSecPPXErrDown, BsXSecPPXErrUp,BsXSecPPYErrDown,BsXSecPPYErrUp);
	TGraphAsymmErrors *BPPbPbCrossGraph = new TGraphAsymmErrors(NBins, BsXsecPbPbX, BsXsecPbPbY,BsXSecPbPbXErrDown, BsXSecPbPbXErrUp,BsXSecPbPbYErrDown,BsXSecPbPbYErrUp);

	BPPbPbCrossGraph->SetLineColor(kGreen+2);
//	BPPbPbCrossGraph->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraph->SetMarkerStyle(20);
	BPPbPbCrossGraph->SetMarkerSize(1);
	BPPbPbCrossGraph->SetMarkerColor(kGreen+2);
	BsPPCrossGraph->SetLineColor(kBlue-9);
//	BsPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BsPPCrossGraph->SetMarkerStyle(21);
	BsPPCrossGraph->SetMarkerSize(1);
	BsPPCrossGraph->SetMarkerColor(kBlue-9);
	TGraphAsymmErrors *BPPPCrossGraphSyst  = new TGraphAsymmErrors(NBins, BsXsecPPX, BsXsecPPY, BsXSecPPXErrDown, BsXSecPPXErrUp, BsXSecPPYSystDown,BsXSecPPYSystUp);
  	TGraphAsymmErrors *BPPbPbCrossGraphSyst    = new TGraphAsymmErrors(NBins, BsXsecPbPbX, BsXsecPbPbY, BsXSecPbPbXErrDown, BsXSecPbPbXErrUp, BsXSecPbPbYSystDown,BsXSecPbPbYSystUp);
	BPPPCrossGraphSyst->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraphSyst->SetLineColor(kBlue-9);
	BPPbPbCrossGraphSyst->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraphSyst->SetLineColor(kGreen-9);

	BPPbPbCrossGraph->Draw("epsame");	
	BsPPCrossGraph->Draw("epsame");	
	BPPPCrossGraphSyst->Draw("5same");	
	BPPbPbCrossGraphSyst->Draw("5same");	

	TLegend* leg = new TLegend(0.75,0.80,0.9,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->AddEntry(BPPbPbCrossGraph,"2018 PbPb","PL");
	leg->AddEntry(BsPPCrossGraph,"2017 pp","PL");
	leg->Draw("same");


	//c->SaveAs("RAAPlots/BsPbPbPPCross.png");
	c->SetLogy();
	c->SaveAs("RAAPlots/BsPbPbPPCrossLog.pdf");


	//2015 References
	TH2D * HisEmpty2 = new TH2D("HisEmpty2","",100,5,60,100,100.0,3000000);
	HisEmpty2->GetXaxis()->SetTitle("B^{0}_{s} p_{T} [GeV/c]");
	HisEmpty2->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty2->GetXaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty2->Draw();
	BsPPCrossGraph->Draw("ep");
	const int NBins2015PP = 3;
	float BsXsecPPX2015[NBins2015PP] = {11,17.5,35.0};
	float BsXSecPPXErrDown2015[NBins2015PP] = {4,2.5,15};
	float BsXSecPPXErrUp2015[NBins2015PP] = {4,2.5,15};
	float BsXsecPPY2015[NBins2015PP] = {316000,34100,3830};
	float BsXSecPPYErrDown2015[NBins2015PP] = {37000,6300,670};
	float BsXSecPPYErrUp2015[NBins2015PP] = {37000,6300,670};


	TGraphAsymmErrors *BsPPCrossGraph2015 = new TGraphAsymmErrors(NBins2015PP, BsXsecPPX2015, BsXsecPPY2015,BsXSecPPXErrDown2015, BsXSecPPXErrUp2015,BsXSecPPYErrDown2015,BsXSecPPYErrUp2015);
	BsPPCrossGraph2015->SetLineColor(kOrange+2);
	BsPPCrossGraph2015->SetMarkerStyle(25);
	BsPPCrossGraph2015->SetMarkerSize(1);
	BsPPCrossGraph2015->SetMarkerColor(kOrange+2);
	BsPPCrossGraph2015->Draw("epSAME");

/*
	TFile * finFONLL = new TFile("../BsFONLL.root");
	finFONLL->cd();
	TGraphAsymmErrors *BPFONLL = (TGraphAsymmErrors*) finFONLL->Get("gaeSigmaBplus");
	BPFONLL->SetLineColor(kBlue-9);
	BPFONLL->SetMarkerStyle(20);
	BPFONLL->SetMarkerSize(1);
	BPFONLL->SetMarkerColor(kBlue-9);
	BPFONLL->Draw("epSAME");

	

	TLegend* leg3 = new TLegend(0.50,0.55,0.80,0.80,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.040);
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->SetLineWidth(3);
	leg3->AddEntry(BsPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg3->AddEntry(BsPPCrossGraph2015,"2015 pp 5.02 TeV","PL");
	leg3->AddEntry(BPFONLL,"FONLL Calculations","PL");
	leg3->Draw("same");


*/
	//c->SaveAs("RAAPlots/Bs/BsCrossComp.png");
	c->SetLogy();

	c->SaveAs("RAAPlots/BsCrossCompLog.pdf");

/*


	//Introduce Non-Fiducial
	//
	//
	TH2D * HisEmpty6 = new TH2D("HisEmpty6","",100,5,60,100,100.0,8000000);
	HisEmpty6->GetXaxis()->SetTitle("B^{+} p_{T} [GeV/c]");
	HisEmpty6->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty6->GetXaxis()->CenterTitle();
	HisEmpty6->GetYaxis()->CenterTitle();
	HisEmpty6->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty6->Draw();

	

	TString InfileBP2 = "FinalFiles/BsPPCorrYieldPT2.root";


	TFile * FileBP2 = new TFile(InfileBP2.Data());
	TH1D * BsCross2 = (TH1D *) FileBP2->Get("CorrDiffHisBin");
	BsCross2->SetMarkerStyle(20);
	BsCross2->SetMarkerSize(1);
	BsCross2->SetMarkerColor(1);
	BsCross2->SetLineColor(1);


	


	float BsXsecPPY2[NBins];


	float BsXSecPPYErrUpPercent2[NBins];
	float BsXSecPPYErrDownPercent2[NBins];
	float BsXSecPPYErrUp2[NBins];
	float BsXSecPPYErrDown2[NBins];

	for(int i = 0; i < NBins; i++){

		BsXsecPPY2[i] = BsCross2->GetBinContent(i+1);
		BsXSecPPYErrUp2[i] = BsCross2->GetBinError(i+1);
		BsXSecPPYErrDown2[i] = BsCross2->GetBinError(i+1);
		BsXSecPPYErrUpPercent2[i] = BsXSecPPYErrUp2[i]/BsXsecPPY2[i];
		BsXSecPPYErrDownPercent2[i] = BsXSecPPYErrDown2[i]/BsXsecPPY2[i];
		
	}

	TGraphAsymmErrors *BsPPCrossGraph2 = new TGraphAsymmErrors(NBins, BsXsecPPX, BsXsecPPY2,BsXSecPPXErrDown, BsXSecPPXErrUp,BsXSecPPYErrDown2,BsXSecPPYErrUp2);


	BsPPCrossGraph2->SetLineColor(kGreen+2);
//	BsPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BsPPCrossGraph2->SetMarkerStyle(20);
	BsPPCrossGraph2->SetMarkerSize(1);
	BsPPCrossGraph2->SetMarkerColor(kGreen+2);

	BsPPCrossGraph2->Draw("ep");
	BsPPCrossGraph->Draw("epSAME");	

	TLegend* leg5 = new TLegend(0.45,0.55,0.75,0.80,NULL,"brNDC");
	leg5->SetBorderSize(0);
	leg5->SetTextSize(0.040);
	leg5->SetTextFont(42);
	leg5->SetFillStyle(0);
	leg5->SetLineWidth(3);
	leg5->AddEntry(BsPPCrossGraph,"Fiducial Region Used","PL");
	leg5->AddEntry(BsPPCrossGraph2,"No Fiducial Region","PL");
	leg5->Draw("same");





	c->SaveAs("RAAPlots/Bs/BsCrossLog.png");

*/
//DONE


	TCanvas * c2 = new TCanvas("c","c",600,600);

	c2->cd();

	c2->SetLeftMargin(0.16);


	TLine * Unity = new TLine(7,1,50,1);
	Unity->SetLineWidth(2);
	Unity->SetLineStyle(2);
	Unity->SetLineColor(1);



	TH2D * HisEmptyRAA = new TH2D("HisEmptyRAA","",100,7,50,100,0,2);
	HisEmptyRAA->GetXaxis()->SetTitle("B^{0}_{s} p_{T} [GeV/c]");
	HisEmptyRAA->GetYaxis()->SetTitle("RAA = #frac{1}{T_{AA}} #frac{dN_{PbPb}/dp_{T}}{d #sigma_{pp}/d p_{T}}");
	HisEmptyRAA->GetXaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->SetTitleOffset(1.8);
	HisEmptyRAA->GetXaxis()->SetTitleOffset(1.3);
	
	HisEmptyRAA->Draw();

	float BsRAAY[NBins];
	float BsRAAX[NBins] = {8.73,12.4,17.2,27.3};
	float BsRAAXErrUp[NBins] = {1.27,2.6,2.8,22.7};
	float BsRAAXErrDown[NBins] = {1.73,2.4,2.2,7.3};
	float BsRAAYErrUp[NBins];
	float BsRAAYErrDown[NBins];


	//2015 Data Sets
	
	const int NBins2015 = 2;

	float BsRAAY2015[NBins2015] = {1.51,0.87};
	float BsRAAX2015[NBins2015] = {11.0,32.5};


	float BsRAAXErrUp2015[NBins2015] = {4.0,17.5};
	float BsRAAXErrDown2015[NBins2015] = {4.0,17.5};

	float BsRAAYErrUp2015[NBins2015] = {0.61,0.30};
	float BsRAAYErrDown2015[NBins2015] = {0.61,0.30};


	float BsRAAYSystUp2015[NBins2015] = {0.50,0.17};
	float BsRAAYSystDown2015[NBins2015] = {0.50,0.17};


	float BsRAAYSystUp[NBins] ;
	float BsRAAYSystDown[NBins];

  float BsRAAYErrUpRatio[NBins];
  float BsRAAYErrDownRatio[NBins];
  float BsRAAYSystUpRatio[NBins];
  float BsRAAYSystDownRatio[NBins];

	for(int i = 0; i < NBins; i++){

		BsRAAY[i] = BsXsecPbPbY[i] / BsXsecPPY[i];
	
		cout << "i = " << i <<  "    BsRAAY[i] = " << BsRAAY[i] << endl;

		BsRAAYErrUp[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYErrUpPercent[i] * BsXSecPbPbYErrUpPercent[i] + BsXSecPPYErrUpPercent[i] * BsXSecPPYErrUpPercent[i]);
		BsRAAYErrDown[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYErrDownPercent[i] * BsXSecPbPbYErrDownPercent[i] + BsXSecPPYErrDownPercent[i] * BsXSecPPYErrDownPercent[i]);
		BsRAAYSystDown[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYSystDownPercent[i] * BsXSecPbPbYSystDownPercent[i] + BsTotalSystDownRatio[i] * BsTotalSystDownRatio[i]);
		BsRAAYSystUp[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYSystUpPercent[i] * BsXSecPbPbYSystUpPercent[i] + BsTotalSystUpRatio[i] * BsTotalSystUpRatio[i]);

    BsRAAYErrUpRatio[i] = BsRAAYErrUp[i]/BsRAAY[i];
    BsRAAYErrDownRatio[i] = BsRAAYErrDown[i]/BsRAAY[i];
    BsRAAYSystUpRatio[i] = BsRAAYSystUp[i]/BsRAAY[i];
    BsRAAYSystDownRatio[i] = BsRAAYSystDown[i]/BsRAAY[i];

		cout << "i = " << i <<  "    BsRAAY[i] = " << BsRAAY[i] <<
      "   BsRAAYErrUp[i]  = "  << BsRAAYErrDownRatio[i] <<
      "  BsRAAYErrDown[i] =   " << BsRAAYErrDownRatio[i] <<
      "   BsRAAYSystUp[i]  = "  << BsRAAYSystUpRatio[i]   <<
      "  BsRAAYSystDown[i] =   " << BsRAAYSystDownRatio[i] << endl;
			
	}

	TGraphAsymmErrors *BsRAAGraph = new TGraphAsymmErrors(NBins, BsRAAX, BsRAAY,BsRAAXErrDown, BsRAAXErrUp,BsRAAYErrDown,BsRAAYErrUp);

	BsRAAGraph->SetName("BsRAAGraph");
	TGraphAsymmErrors *BsRAAGraphSyst = new TGraphAsymmErrors(NBins, BsRAAX, BsRAAY,BsRAAXErrDown, BsRAAXErrUp,BsRAAYSystDown,BsRAAYSystUp);
	TGraphAsymmErrors *BsRAAGraph2015 = new TGraphAsymmErrors(NBins2015, BsRAAX2015, BsRAAY2015,BsRAAXErrDown2015, BsRAAXErrUp2015,BsRAAYErrDown2015,BsRAAYErrUp2015);
	TGraphAsymmErrors *BsRAAGraphSyst2015 = new TGraphAsymmErrors(NBins2015, BsRAAX2015, BsRAAY2015,BsRAAXErrDown2015, BsRAAXErrUp2015,BsRAAYSystDown2015,BsRAAYSystUp2015);

	BsRAAGraph->SetLineColor(kBlue-9);
	BsRAAGraph->SetMarkerStyle(20);
	BsRAAGraph->SetMarkerSize(1);
	BsRAAGraph->SetMarkerColor(kBlue-9);
	BsRAAGraphSyst->SetFillColorAlpha(kBlue-9,0.5);
	BsRAAGraphSyst->SetLineColor(kBlue-9);
	BsRAAGraph2015->SetLineColor(kOrange+1);
	BsRAAGraph2015->SetMarkerColor(kOrange+1);
	BsRAAGraph2015->SetMarkerStyle(21);
	BsRAAGraph2015->SetMarkerSize(1);
	BsRAAGraphSyst2015->SetFillColorAlpha(kOrange+1,0.5);
	BsRAAGraphSyst2015->SetLineColor(kOrange+1);

	HisEmptyRAA->Draw();
	BsRAAGraph->Draw("ep");
	BsRAAGraphSyst->Draw("5same");
	Unity->Draw("SAME");
	c2->SaveAs("RAAPlots/BsRAA.pdf");

	HisEmptyRAA->Draw();

	//BsRAAGraph2015->Draw("epSAME");
	BsRAAGraphSyst2015->Draw("5same");
	BsRAAGraph2015->Draw("epSAME");
	BsRAAGraphSyst->Draw("5same");
	BsRAAGraph->Draw("ep");

	Unity->Draw("SAME");



	TLegend* leg2 = new TLegend(0.65,0.80,0.9,0.85,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	leg2->AddEntry(BsRAAGraph,"2018 PbPb + 2017 pp","PL");
	leg2->AddEntry(BsRAAGraph2015,"2015 PbPb + 2015 pp","PL");
	leg2->Draw("same");

	c2->SaveAs("RAAPlots/BsRAAComparison.pdf");


	TFile * fout = new TFile("OutFile/BsRAA.root","RECREATE");
	fout->cd();
	BsRAAGraph->Write();
	fout->Close();


  std::vector<float> globUncert(NBins, 0.077);
  // summary of errors (in ratio, not percent)
  std::vector<int> ptbins = {7, 10, 15, 20, 50};
  std::vector<float> abscissae = {8.5, 12.5, 17.5, 35.0};

		gSystem->mkdir("../../MakeFinalPlots/NominalPlots/RAA/dataSource/" ,true );

  string outFile = "../../MakeFinalPlots/NominalPlots/RAA/dataSource/RAA_pt_Bs_New.txt";
  ofstream out;
  out.open(outFile);
  unsigned columnWidth = 14;
  out << std::left << std::setw(columnWidth) <<
    "ptmin" << std::setw(columnWidth) << "ptmax" << std::setw(columnWidth) <<
    "central_val" << std::setw(columnWidth) <<
    "statUp" << std::setw(columnWidth) << "statDown" << std::setw(columnWidth) <<
    "systUp" << std::setw(columnWidth) << "systDown" << std::setw(columnWidth) <<
    "glbUp" << std::setw(columnWidth) << "glbDown" << std::setw(columnWidth) <<
    "abscissae" << endl;
  for (auto i = 0; i < NBins; ++i ) {
    out << std::setw(columnWidth) <<
      ptbins[i] << std::setw(columnWidth) << ptbins[i + 1] << std::setw(columnWidth) <<
      setprecision(3) << BsRAAY[i] << std::setw(columnWidth) <<
      setprecision(2) << std::defaultfloat <<
      BsRAAYErrUpRatio[i] << std::setw(columnWidth) <<
      BsRAAYErrDownRatio[i] << std::setw(columnWidth) <<
      BsRAAYSystUpRatio[i] << std::setw(columnWidth) <<
      BsRAAYSystDownRatio[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      setprecision(3) << abscissae[i] << "\n";
  }
  out.close();
}
