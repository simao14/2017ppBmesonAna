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
#include "../../../henri2022/parameter.h" 

using namespace std;
using std::cout;
using std::endl;

// MOST OF THE VARIABLES SAY BP BUT ARE WORKING FOR Bs AS WELL 
// THIS CODE CAN BE EASILY EXTENDED TO ACCOUNT FOR A FUTURE MESON... 
// JUST FOLLOW THE LOGIC

void Bmeson_Comparisons(int meson_n){


	// scaled lower data pt bins to full rapidity
   	 	//constexpr bool scaleMC = false;                              
    // swap the lower fonll bins with full rapidity
    constexpr bool fidFONLL = true;

	TString B_m ;
	TString t_tree ;
	TString b_m;
	int NBins = 7;
    vector<double> scaledPt;
	int NBinsLow ;
  	int NBinsHigh ;
	int NBins2015 ;
	double hist3max;
	
	if(meson_n == 0){
		NBins = nptBinsBP;
		B_m = "BP";
		b_m = "bp";
		t_tree = "ntKp";
		NBinsLow = 2 ;
		NBinsHigh = 5;
		NBins2015 = 5;
		hist3max = 1.2;
		scaledPt = {5, 7, 10};
	} else {
		NBins = nptBins;
		B_m = "Bs";
		b_m = "bs";
		t_tree = "ntphi";
		NBinsLow = 1 ;
		NBinsHigh = 3;
		NBins2015 = 3;
		hist3max = 1.5;
		scaledPt = {7, 10};
	}

	

	gSystem->mkdir("Plots/", true);
	TString InfileB = Form("../../../%s/EffAna/FinalFiles/%sPPCorrYieldPT.root",B_m.Data(),B_m.Data());
	TFile * FileB= new TFile(InfileB.Data());
	TH1D * BCross = (TH1D *) FileB->Get("CorrDiffHisBin");
	BCross->SetMarkerStyle(20);
	BCross->SetMarkerSize(1);
	BCross->SetMarkerColor(1);
	BCross->SetLineColor(1);
	
	double ptBins[NBins+1];


	if(meson_n==0) { for( int c=0; c <NBins+1; c++){ ptBins[c]=ptbinsvecBP[c];}}
	else{            for( int c=0; c <NBins+1; c++){ ptBins[c]=ptbinsvec[c];}}
	//center of the bin
	float BPXsecPPX[NBins];
	for( int c=0; c <NBins; c++){ BPXsecPPX[c]=(ptBins[c]+ptBins[c+1])/2;}

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

  vector<double> factor = scaleFactor(Form("/data3/tasheng/presel/%sMC_nom.root", B_m.Data()), t_tree.Data(), scaledPt);
  int factor_start;
  if (meson_n==0){factor_start=1;}     
  else{factor_start=0;}
  /*if(meson_n==1){
  for (auto i = factor_start; i < factor.size(); ++i) {
    cout << "applying scaling factor: " << factor[i] << "\n";
    BPXsecPPY2DScaled[i] *= factor[i];
	BPXSecPPY2DErrUpScaled[i] *= factor[i];
    BPXSecPPY2DErrDownScaled[i] *= factor[i];
  		}
  }*/
  
float BXSecPPXErrUp[NBins] ;
float BXSecPPXErrDown[NBins] ;



if (meson_n == 0){
	vector<float> vect_BXSecPPXErrUp{1,1.5,2.5,2.5,5,10,5};
	vector<float> vect_BXSecPPXErrDown{1,1.5,2.5,2.5,5,10,5};
	
	for( int c=0; c <NBins; c++){ 
		BXSecPPXErrUp[c]=vect_BXSecPPXErrUp[c];
		BXSecPPXErrDown[c]=vect_BXSecPPXErrDown[c];
		}
} else {
	vector<float> vect_BXSecPPXErrUp {1.5,2.5,2.5,15};
	vector<float> vect_BXSecPPXErrDown {1.5,2.5,2.5,15};
	for( int c=0; c <NBins; c++){ 
		BXSecPPXErrUp[c]=vect_BXSecPPXErrUp[c];
		BXSecPPXErrDown[c]=vect_BXSecPPXErrDown[c];
		}
}




// FOR BP ONLY FOR BP ONLY (for now)
// Values of PbPb yields


	float BPXsecPbPbY[NBins] ;
	float BPXsecPbPbX[NBins] ;
	float BPXSecPbPbXErrUp[NBins] ;
	float BPXSecPbPbXErrDown[NBins] ;
	float BPXSecPbPbYErrUpRatio[NBins] ;
	float BPXSecPbPbYErrDownRatio[NBins] ;
	float BPXSecPbPbYSystUpRatio[NBins];
	float BPXSecPbPbYSystDownRatio[NBins] ;
	float BPXSecPbPbYErrUp[NBins] ;
	float BPXSecPbPbYErrDown[NBins] ;
	float BPXSecPbPbYSystUp[NBins] ;
	float BPXSecPbPbYSystDown[NBins] ;
if(meson_n == 0){
	for( int c=0; c <NBins; c++){ 
	 BPXsecPbPbY[c]= vect_BPXsecPbPbY[c] ;
	 BPXsecPbPbX[c]= vect_BPXsecPbPbX[c];
	 BPXSecPbPbXErrUp[c]= vect_BPXSecPbPbXErrUp[c];
	 BPXSecPbPbXErrDown[c]= vect_BPXSecPbPbXErrDown[c];
	 BPXSecPbPbYErrUpRatio[c]= vect_BPXSecPbPbYErrUpRatio[c];
	 BPXSecPbPbYErrDownRatio[c]= vect_BPXSecPbPbYErrDownRatio[c];
	 BPXSecPbPbYSystUpRatio[c]= vect_BPXSecPbPbYSystUpRatio[c];
	 BPXSecPbPbYSystDownRatio[c] = vect_BPXSecPbPbYSystDownRatio[c];
	 }
	for(int i = 0; i < NBins; i++){
		BPXSecPbPbYErrUp[i] = BPXSecPbPbYErrUpRatio[i] * BPXsecPbPbY[i];
		BPXSecPbPbYErrDown[i] = BPXSecPbPbYErrDownRatio[i] * BPXsecPbPbY[i];
		BPXSecPbPbYSystUp[i] = (BPXSecPbPbYSystUpRatio[i]) * BPXsecPbPbY[i];
		BPXSecPbPbYSystDown[i] = (BPXSecPbPbYSystDownRatio[i]) * BPXsecPbPbY[i];
		}
}



// Values of PbPb yields
// FOR BP ONLY FOR BP ONLY (for now)

	//Syst Add Up PP//
  	TString errorFile = Form("../../../2DMapSyst/OutFiles/%sError2D.root", B_m.Data());
  	TFile fError(errorFile);

	TH1D * TnPSyst = (TH1D *) fError.Get("TnPSyst");
	TH1D * BptSyst = (TH1D *) fError.Get("BptSyst");
	TH1D * MCDataSyst = (TH1D *) fError.Get("MCDataSyst");
  	if (!MCDataSyst) MCDataSyst = (TH1D *) fError.Get("BDTSyst");

	TString pdfErrorFile = Form("../../../%s_pdf.root",b_m.Data());
	TFile fPdfError(pdfErrorFile);
	TGraph* pdfSyst = (TGraph *) fPdfError.Get(Form("%s_error",b_m.Data()));
	TString trackSelErrorFile = "../../../syst_track_sel.root";
	TFile fTrackSelError(trackSelErrorFile);
	TGraph* trackSelSyst = (TGraph *) fTrackSelError.Get(Form("%s_track_sel_error", b_m.Data()));

	float BPXSecPPYSystUp[NBins];
	float BPXSecPPYSystDown[NBins];
	float BPXSecPPYSystUpScaled[NBins];
	float BPXSecPPYSystDownScaled[NBins];

  // percent error
  	int B_nu;
  	if (meson_n == 0){ B_nu = 5 ;}
  	else { B_nu =10 ;}
	float BPTrackingSyst[NBins];
	for( int c=0; c < NBins; c++){BPTrackingSyst[c]= B_nu ;}
	float BPMCDataSyst[NBins];
	float BPPDFSyst[NBins];
	float BPTrackSelSyst[NBins];
	float BPPtShapeSyst[NBins];
	float BPTnPSystDown[NBins];
	float BPTnPSystUp[NBins];

  // Get systematics from input files
    for (auto ibin = 0; ibin < NBins; ++ibin){
    BPMCDataSyst[ibin] = MCDataSyst->GetBinContent(ibin + 1);
    BPPtShapeSyst[ibin] = BptSyst->GetBinContent(ibin + 1);
    BPTnPSystDown[ibin] = TnPSyst->GetBinContent(ibin + 1);
    // TnP systematics are symmetric in the binned pT case
    BPTnPSystUp[ibin] = BPTnPSystDown[ibin];
    BPPDFSyst[ibin] = pdfSyst->GetY()[ibin];
    BPTrackSelSyst[ibin] = trackSelSyst->GetY()[ibin];
  	}

  // RMS of all the errors
	float BPTotalSystDownRatio[NBins];
	float BPTotalSystUpRatio[NBins];

	for(int i = 0; i < NBins; i++){
		BPTotalSystDownRatio[i] = TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BPMCDataSyst[i], 2) +
                                          TMath::Power(BPPDFSyst[i], 2) + TMath::Power(BPTrackSelSyst[i], 2) +
                                          TMath::Power(BPPtShapeSyst[i], 2) + TMath::Power(BPTnPSystDown[i], 2)) / 100;
    	BPTotalSystUpRatio[i] = TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BPMCDataSyst[i], 2) +
                                        TMath::Power(BPPDFSyst[i], 2) + TMath::Power(BPTrackSelSyst[i], 2) +
                                        TMath::Power(BPPtShapeSyst[i], 2) + TMath::Power(BPTnPSystUp[i], 2)) / 100;
	}

  // global uncertainty from branching ratio and luminosity
  // Fixed, copied from the paper draft
	double numb;
	if(meson_n == 0){numb = 0.035;}
	else {numb = 0.077;}
  vector<float> globUncert(NBins, numb);
	for(int i = 0; i < NBins; i++){
		BPXSecPPYSystUp[i] = BPXsecPPY2D[i] * BPTotalSystUpRatio[i];
		BPXSecPPYSystDown[i] = BPXsecPPY2D[i] * BPTotalSystDownRatio[i];
    	BPXSecPPYSystDownScaled[i] = BPXsecPPY2DScaled[i] * BPTotalSystDownRatio[i];
    	BPXSecPPYSystUpScaled[i] = BPXsecPPY2DScaled[i] * BPTotalSystUpRatio[i];
		}

	// CREATE THE CANVAS
	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",700,700);
	c->cd();    
	c->SetLeftMargin(0.15);
 
	//Setup the Syst
	TH2D * HisEmpty;
	if(meson_n == 0) { 
		HisEmpty = new TH2D("HisEmpty","",100,5,60,100,200.0,2000000);
		HisEmpty->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	} else {
		HisEmpty = new TH2D("HisEmpty","",100,7,50,100,200.0,2000000);
		HisEmpty->GetXaxis()->SetTitle("p_{T} [GeV/c]");
		}
	HisEmpty->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb c/GeV]");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	HisEmpty->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty->GetXaxis()->SetTitleOffset(1.3);	
	HisEmpty->Draw();

  // separate plots for different fiducial regions
  vector<double> BPXsecPPXLow = {BPXsecPPX, BPXsecPPX + NBinsLow};
  vector<double> BPXsecPPXHigh = {BPXsecPPX + NBinsLow, BPXsecPPX + NBins};
  vector<double> BPXsecPPXErrDownLow = {BXSecPPXErrDown, BXSecPPXErrDown + NBinsLow};
  vector<double> BPXsecPPXErrDownHigh = {BXSecPPXErrDown + NBinsLow, BXSecPPXErrDown + NBins};
  vector<double> BPXsecPPXErrUpLow = {BXSecPPXErrUp, BXSecPPXErrUp + NBinsLow};
  vector<double> BPXsecPPXErrUpHigh = {BXSecPPXErrUp + NBinsLow, BXSecPPXErrUp + NBins};
  vector<double> BPXsecPPYLow = {BPXsecPPY2DScaled, BPXsecPPY2DScaled + NBinsLow};
  vector<double> BPXsecPPYHigh = {BPXsecPPY2D + NBinsLow, BPXsecPPY2D + NBins};
  vector<double> BPXsecPPYErrDownLow = {BPXSecPPY2DErrDownScaled, BPXSecPPY2DErrDownScaled + NBinsLow};
  vector<double> BPXsecPPYErrDownHigh = {BPXSecPPY2DErrDown + NBinsLow, BPXSecPPY2DErrDown + NBins};
  vector<double> BPXsecPPYErrUpLow = {BPXSecPPY2DErrUpScaled, BPXSecPPY2DErrUpScaled + NBinsLow};
  vector<double> BPXsecPPYErrUpHigh = {BPXSecPPY2DErrUp + NBinsLow, BPXSecPPY2DErrUp + NBins};

  // separate plots for different fiducial regions
	TGraphAsymmErrors *BPPPCrossGraph2D = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY2D, BXSecPPXErrDown, BXSecPPXErrUp, BPXSecPPY2DErrDown, BPXSecPPY2DErrUp);
	TGraphAsymmErrors *BPPPCrossGraph2DLow = new TGraphAsymmErrors(NBinsLow, BPXsecPPXLow.data(), BPXsecPPYLow.data(), BPXsecPPXErrDownLow.data(), BPXsecPPXErrUpLow.data(), BPXsecPPYErrDownLow.data(), BPXsecPPYErrUpLow.data());
	TGraphAsymmErrors *BPPPCrossGraph2DHigh = new TGraphAsymmErrors(NBinsHigh, BPXsecPPXHigh.data(), BPXsecPPYHigh.data(), BPXsecPPXErrDownHigh.data(), BPXsecPPXErrUpHigh.data(), BPXsecPPYErrDownHigh.data(), BPXsecPPYErrUpHigh.data());                                                   
	TGraphAsymmErrors *BPPPCrossGraph2DSyst  = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY2D, BXSecPPXErrDown, BXSecPPXErrUp, BPXSecPPYSystDown,BPXSecPPYSystUp);                 											
	TGraphAsymmErrors *BPPPCrossGraph2DScaledSyst  = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY2DScaled, BXSecPPXErrDown, BXSecPPXErrUp, BPXSecPPYSystDownScaled, BPXSecPPYSystUpScaled);
  // separate plots for different fiducial regions
 
 /*
  double Xchange;
  double Ychange;
  if (meson_n==0){            
	BPPPCrossGraph2DLow->GetPoint(0,Xchange,Ychange);
	BPPPCrossGraph2DHigh->AddPoint(Xchange,Ychange);
	BPPPCrossGraph2DLow->RemovePoint(0);
    }
	*/

	// CrossSection 
	BPPPCrossGraph2D->SetMarkerStyle(20);
	BPPPCrossGraph2D->SetMarkerSize(1);
	if (meson_n==0){
		BPPPCrossGraph2D->SetLineColor(kGreen-9);
		BPPPCrossGraph2D->SetMarkerColor(kGreen-9);
		BPPPCrossGraph2DSyst->SetFillColorAlpha(kGreen-9,0.5);
		BPPPCrossGraph2DSyst->SetLineColor(kGreen-9);
	} else {
		BPPPCrossGraph2D->SetLineColor(kBlue-9);
		BPPPCrossGraph2D->SetMarkerColor(kBlue-9);
		BPPPCrossGraph2DSyst->SetFillColorAlpha(kBlue-9,0.5);
		BPPPCrossGraph2DSyst->SetLineColor(kBlue-9);

	}
	
	TLegend* leged = new TLegend(0.8,0.70,0.95,0.95,NULL,"brNDC");
	leged->SetBorderSize(0);
	leged->SetTextSize(0.025);     
	leged->SetTextFont(42);
	leged->SetFillStyle(0);
	if(meson_n == 0) {leged->AddEntry(BPPPCrossGraph2D,"B^{+}","P");}
	else {leged->AddEntry(BPPPCrossGraph2D,"B^{0}_{s}","P");}
	leged->Draw();
	BPPPCrossGraph2D->Draw("ep");	
	BPPPCrossGraph2DSyst->Draw("5same");	
	//c->SaveAs(Form("Plots/%s/%sCrossONLY.png",B_m.Data(), B_m.Data()));
	c->SetLogy();
	c->SaveAs(Form("Plots/%sCrossONLYLog.pdf", B_m.Data()));
	// CrossSection (log scale) 



// FOR BP ONLY (for now), PbPb plots for comparison
 	if (meson_n == 0){
	TGraphAsymmErrors *BPPbPbCrossGraph;
  	TGraphAsymmErrors *BPPbPbCrossGraphSyst;
	BPPbPbCrossGraph = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY,BPXSecPbPbXErrDown, BPXSecPbPbXErrUp,BPXSecPbPbYErrDown,BPXSecPbPbYErrUp);
  	BPPbPbCrossGraphSyst  = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY, BPXSecPbPbXErrDown, BPXSecPbPbXErrUp, BPXSecPbPbYSystDown,BPXSecPbPbYSystUp);
	BPPbPbCrossGraph->SetLineColor(kOrange+1);
	BPPbPbCrossGraph->SetMarkerStyle(21);
	BPPbPbCrossGraph->SetMarkerSize(1);
	BPPbPbCrossGraph->SetMarkerColor(kOrange+1);
	BPPbPbCrossGraphSyst->SetFillColorAlpha(kOrange+1,0.5);
	BPPbPbCrossGraphSyst->SetLineColor(kOrange+1);

	TCanvas * c2New = new TCanvas("c2New","c2New",700,700);
	c2New->cd();
	c2New->SetLeftMargin(0.15);
	HisEmpty->Draw();

	TLegend* leg = new TLegend(0.7,0.8,0.95,0.9,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.025);     
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(BPPbPbCrossGraph,"2018 PbPb","P");
	leg->AddEntry(BPPPCrossGraph2D,"2017 pp","P");
	leg->Draw();
	BPPbPbCrossGraph->Draw("ep");	
	BPPPCrossGraph2D->Draw("ep");	
	BPPbPbCrossGraphSyst->Draw("5same");	
	BPPPCrossGraph2DSyst->Draw("5same");

	//c2New->SaveAs("Plots/BP/BPPbPbPPCross.png");
	c2New->SetLogy();
	c2New->SaveAs(Form("Plots/%sPbPbPPCrossLog.pdf", B_m.Data()));
	}
 // FOR BP ONLY (for now), PbPb plots for comparison

	//2015 Reference (	Big plots  )
	TCanvas * cRatio = new TCanvas("cRatio","cRatio",700,800);
	gStyle->SetPadTopMargin(0.08); 
	gStyle->SetPadBottomMargin(0.12); 
	gStyle->SetPadLeftMargin(0.16); 
 	gStyle->SetPadRightMargin(0.04);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetTitleColor(1, "XYZ");
	gStyle->SetTitleFont(42, "XYZ");
	//gStyle->SetTitleSize(0.06, "XYZ");
 	gStyle->SetTickLength(0.03, "XYZ");
  	gStyle->SetNdivisions(510, "XYZ");
	gStyle->SetHatchesLineWidth(5);
 	gStyle->SetHatchesSpacing(0.05);

	TPad * MyPad1;
	MyPad1 = new TPad("MyPad1","",0,0.45,1,1.0);
	MyPad1->SetBottomMargin(0);
	MyPad1->Draw();

	TPad * MyPad2;
	MyPad2 = new TPad("MyPad2","",0,0.25,1,0.45);
	MyPad2->SetBottomMargin(0.3);
	MyPad2->SetTopMargin(0);

	MyPad2->Draw();


	/*TPad * MyPad3;
	MyPad3 = new TPad("MyPad3","",0,0.00,1,0.25);
	MyPad3->SetBottomMargin(0.2);
	MyPad3->SetTopMargin(0);
	MyPad3->Draw();*/


	MyPad1->cd();
	TH2D * HisEmpty2;
	if (meson_n == 0){
		HisEmpty2 = new TH2D("HisEmpty2","",100,5,60,100,500.0,2300000);
		HisEmpty2->GetXaxis()->SetTitle("p_{T} [GeV/c]");}
	else {	
		HisEmpty2 = new TH2D("HisEmpty2","",100,7,50,100,500.0,2300000);
		HisEmpty2->GetXaxis()->SetTitle("p_{T} [GeV/c]");}
	HisEmpty2->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb c/GeV]");
	HisEmpty2->GetYaxis()->SetTitleSize(40);
	HisEmpty2->GetXaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->CenterTitle();
	//HisEmpty2->SetTitle("B^{+} Cross Section With Fiducial Region");
	HisEmpty2->GetYaxis()->SetTitleOffset(1.0);
	HisEmpty2->GetXaxis()->SetTitleOffset(1.2);
	HisEmpty2->GetYaxis()->SetTitleSize(0.06);
	HisEmpty2->Draw();

		float BXsecPPX2015[NBins2015] ;
		float BXSecPPXErrDown2015[NBins2015] ;
		float BXSecPPXErrUp2015[NBins2015] ;
		float BXsecPPY2015[NBins2015] ;
		float BXSecPPYErrDown2015[NBins2015] ;
		float BXSecPPYErrUp2015[NBins2015] ;
		float BXSecPPYSystDown2015[NBins2015] ;
		float BXSecPPYSystUp2015[NBins2015] ;
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

		/*if(meson_n==0){  
			int i=0;
			BXsecPPY2015[i] *= (1/factor[1+i]);
			BXSecPPYErrDown2015[i] *= (1/factor[1+i]);
			BXSecPPYErrUp2015[i] *= (1/factor[1+i]);
			}*/

	TGraphAsymmErrors *BPPPCrossGraph2015 = new TGraphAsymmErrors(NBins2015, BXsecPPX2015, BXsecPPY2015,BXSecPPXErrDown2015, BXSecPPXErrUp2015,BXSecPPYErrDown2015,BXSecPPYErrUp2015);
	TGraphAsymmErrors *BPPPCrossGraph2015Syst = new TGraphAsymmErrors(NBins2015, BXsecPPX2015, BXsecPPY2015,BXSecPPXErrDown2015, BXSecPPXErrUp2015,BXSecPPYSystDown2015,BXSecPPYSystUp2015);
	BPPPCrossGraph2015->SetLineColor(kOrange+1);
	BPPPCrossGraph2015->SetMarkerColor(kOrange+1);
	BPPPCrossGraph2015->SetMarkerStyle(20);
	BPPPCrossGraph2015->SetMarkerSize(1);

	BPPPCrossGraph2DLow->SetMarkerStyle(25);
	BPPPCrossGraph2DLow->SetMarkerSize(1);
	BPPPCrossGraph2DHigh->SetMarkerStyle(21);
	BPPPCrossGraph2DHigh->SetMarkerSize(1);
	
	BPPPCrossGraph2015Syst->SetLineColor(kOrange+1);
	BPPPCrossGraph2015Syst->SetFillColorAlpha(kOrange+1, 0.5);

	if (meson_n==0){
		BPPPCrossGraph2DLow->SetLineColor(kGreen-9);
		BPPPCrossGraph2DLow->SetMarkerColor(kGreen-9);
		BPPPCrossGraph2DHigh->SetLineColor(kGreen-9);
		BPPPCrossGraph2DHigh->SetMarkerColor(kGreen-9);
		BPPPCrossGraph2DScaledSyst->SetFillColorAlpha(kGreen-9, 0.5);
		BPPPCrossGraph2DScaledSyst->SetLineColor(kGreen-9);
	} else{
		BPPPCrossGraph2DLow->SetLineColor(kBlue-9);
		BPPPCrossGraph2DLow->SetMarkerColor(kBlue-9);
		BPPPCrossGraph2DHigh->SetLineColor(kBlue-9);
		BPPPCrossGraph2DHigh->SetMarkerColor(kBlue-9);
		BPPPCrossGraph2DScaledSyst->SetFillColorAlpha(kBlue-9, 0.5);
		BPPPCrossGraph2DScaledSyst->SetLineColor(kBlue-9);
	}


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
//if(meson_n==0){          //HENRIQUE STUFF
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
//}

	/*double XTempChange2;
	double YTempChange2;
	double YErrLowTemp2;
	double YErrHighTemp2;
	TGraphAsymmErrors *BPPPCrossGraph2015Low= new TGraphAsymmErrors();
	TGraphAsymmErrors *BPPPCrossGraph2015SystLow= new TGraphAsymmErrors();
	BPPPCrossGraph2015SystLow->SetLineColor(kViolet+1);
	BPPPCrossGraph2015SystLow->SetFillColorAlpha(kViolet+1, 0.5);
	BPPPCrossGraph2015Low->SetLineColor(kViolet+1);
	BPPPCrossGraph2015Low->SetMarkerColor(kViolet+1);
	BPPPCrossGraph2015Low->SetMarkerStyle(20);
	BPPPCrossGraph2015Low->SetMarkerSize(1);
	if(meson_n==0){
		int i=0;
		BPPPCrossGraph2015->GetPoint(i,XTempChange2,YTempChange2);
		BPPPCrossGraph2015Low->AddPoint(XTempChange2,YTempChange2);
		BPPPCrossGraph2015Low->SetPointEYhigh(i,BPPPCrossGraph2015->GetErrorYhigh(i));
		BPPPCrossGraph2015Low->SetPointEYlow(i,BPPPCrossGraph2015->GetErrorYlow(i));
		BPPPCrossGraph2015Low->SetPointEXhigh(i,BPPPCrossGraph2015->GetErrorXhigh(i));
		BPPPCrossGraph2015Low->SetPointEXlow(i,BPPPCrossGraph2015->GetErrorXlow(i));
		BPPPCrossGraph2015->RemovePoint(0);

		BPPPCrossGraph2015Syst->GetPoint(i,XTempChange2,YTempChange2);
		BPPPCrossGraph2015SystLow->AddPoint(XTempChange2,YTempChange2);
		BPPPCrossGraph2015SystLow->SetPointEYhigh(i,BPPPCrossGraph2015Syst->GetErrorYhigh(i));
		BPPPCrossGraph2015SystLow->SetPointEYlow(i,BPPPCrossGraph2015Syst->GetErrorYlow(i));
		BPPPCrossGraph2015SystLow->SetPointEXhigh(i,BPPPCrossGraph2015Syst->GetErrorXhigh(i));
		BPPPCrossGraph2015SystLow->SetPointEXlow(i,BPPPCrossGraph2015Syst->GetErrorXlow(i));
		BPPPCrossGraph2015Syst->RemovePoint(0);

	}
	
	BPPPCrossGraph2015Low->Draw("epSAME");
	BPPPCrossGraph2015SystLow->Draw("epSAME");*/
	BPPPCrossGraph2015Syst->Draw("5same");	
	BPPPCrossGraph2015->Draw("epSAME");
	BPPPCrossGraph2DScaledSyst->Draw("5same");
	BPPPCrossGraph2DLow->Draw("epSAME");
	BPPPCrossGraph2DHigh->Draw("epSAME");
	BPFONLL->Draw("5");
	BFONLLLow->Draw("5");

	TLegend* leg3 = new TLegend(0.75,0.64,0.9,0.85,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.025);     
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->AddEntry(BPPPCrossGraph2DLow,"2017 pp (|y| < 1.5)","PL");	
	leg3->AddEntry(BPPPCrossGraph2DHigh,"2017 pp","PL");	
	leg3->AddEntry(BPPPCrossGraph2015,"2015 ","PL");
	//leg3->AddEntry(BPPPCrossGraph2015Low,"2015 pp 5.02 TeV(scaled)","PL");
	leg3->AddEntry(BPFONLL,"FONLL","f");
	leg3->AddEntry(BFONLLLow,"FONLL (|y| < 1.5)","f");
	leg3->Draw("same");
	MyPad1->Update();

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
  std::vector<double> RatioDataYLow(1);
  std::vector<double> RatioDataYLowErr(1);
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
	cRatio->cd();
	MyPad2->cd();
	
	TLegend* leg4 = new TLegend(0.8,0.3,0.95,0.66,NULL,"brNDC");
	leg4->SetBorderSize(0);
	leg4->SetTextSize(0.025);     
	leg4->SetTextFont(42);
	leg4->SetFillStyle(0);

TH2D * HisEmpty3;
if (meson_n == 0){
	HisEmpty3 = new TH2D("HisEmpty3","",100,5,60,100,0,2);
	HisEmpty3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
} else {
	HisEmpty3 = new TH2D("HisEmpty3","",100,7,50,100,0,2);
	HisEmpty3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
}
	HisEmpty3->GetYaxis()->SetTitle("Data/FONLL");
	HisEmpty3->GetXaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->SetTitleOffset(0.5);
	HisEmpty3->GetYaxis()->SetLabelSize(0.1);
	HisEmpty3->GetXaxis()->SetLabelSize(0.1);
	HisEmpty3->GetXaxis()->SetTitleSize(0.12);
	HisEmpty3->GetYaxis()->SetTitleSize(0.12);
	HisEmpty3->GetYaxis()->SetRangeUser(0.4, 1.5);
	HisEmpty3->Draw();

	TGraphAsymmErrors *Ratio2 = new TGraphAsymmErrors(NBins2015, BXsecPPX2015, Ratio2Y,BXSecPPXErrDown2015, BXSecPPXErrUp2015,Ratio2YErr,Ratio2YErr);
		Ratio2->SetLineColor(kOrange+1);
		Ratio2->SetMarkerStyle(21);
		Ratio2->SetMarkerSize(1);
		Ratio2->SetMarkerColor(kOrange+1);
		//leg4->AddEntry(Ratio2,"2017/2015","PL");	
		//Ratio2->Draw("epSAME");

  	TGraphAsymmErrors *RatioDataLow;
  	if (meson_n ==0){
		RatioDataLow = new TGraphAsymmErrors(1, BPXsecPPXLow.data() + 1, RatioDataYLow.data(),BPXsecPPXErrDownLow.data() + 1, BPXsecPPXErrUpLow.data() + 1, RatioDataYLowErr.data(), RatioDataYLowErr.data());
		RatioDataLow->SetLineColor(kOrange+1);
		RatioDataLow->SetMarkerStyle(25);
		RatioDataLow->SetMarkerSize(1);
		RatioDataLow->SetMarkerColor(kOrange+1);
		//RatioDataLow->Draw("epSAME");
		//leg4->AddEntry(RatioDataLow,"2017/2015 (scaled to |y| < 2.4)","PL");	
		}

	TLine * Unity2 = new TLine(0,0,0,0);
	if (meson_n == 0){Unity2 = new TLine(5,1,60,1);}
	else {Unity2 = new TLine(7,1,50,1);}
	Unity2->SetLineWidth(2);
	Unity2->SetLineStyle(2);
	Unity2->SetLineColor(1);
	Unity2->Draw("SAME");

	MyPad2->Update();

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

	//cRatio->cd();
	//MyPad3->cd();

TH2D * HisEmpty4;
if (meson_n == 0){
	HisEmpty4 = new TH2D("HisEmpty4","",100,5,60,100,0,2);
	HisEmpty4->GetXaxis()->SetTitle("B^{+} p_{T} [GeV/c]");
} else {
	HisEmpty4 = new TH2D("HisEmpty4","",100,7,50,100,0,2);
	HisEmpty4->GetXaxis()->SetTitle("B^{0}_{s} p_{T} [GeV/c]");
}
	HisEmpty4->GetYaxis()->SetTitle("2017 Data/FONLL");
	HisEmpty4->GetXaxis()->CenterTitle();
	HisEmpty4->GetYaxis()->CenterTitle();
	HisEmpty4->GetYaxis()->SetTitleOffset(1.0);
	HisEmpty4->GetYaxis()->SetTitleSize(40);
	HisEmpty4->GetYaxis()->SetRangeUser(0.4, 1.5);
	//HisEmpty4->Draw();

 // Get ratio plots

double binlow[NBinsLow];
double glbSystUp;
double glbSystDown;
double bl_low[NBinsLow];
double bl_low_yStatL[NBinsLow];
double bl_low_yStatH[NBinsLow];
double bl_low_xErrL[NBinsLow];
double bl_low_xErrH[NBinsLow];
double bl_low_ySystL[NBinsLow];
double bl_low_ySystH[NBinsLow];
double binhigh[NBins-NBinsLow];
double bl_high[NBins-NBinsLow];
double bl_high_yStatL[NBins-NBinsLow];
double bl_high_yStatH[NBins-NBinsLow];
double bl_high_xErrL[NBins-NBinsLow];
double bl_high_xErrH[NBins-NBinsLow];
double bl_high_ySystL[NBins-NBinsLow];
double bl_high_ySystH[NBins-NBinsLow];
int NBinsLow2015;
if (meson_n==0){NBinsLow2015=1;} else {NBinsLow2015=0;}
double binlow_2015[NBinsLow2015];
double bl_low_2015[NBinsLow2015];
double bl_low_2015_yStatL[NBinsLow2015];
double bl_low_2015_yStatH[NBinsLow2015];
double bl_low_2015_xErrL[NBinsLow2015];
double bl_low_2015_xErrH[NBinsLow2015];
double bl_low_2015_ySystL[NBinsLow2015];
double bl_low_2015_ySystH[NBinsLow2015];
double binhigh_2015[NBins2015-NBinsLow2015];
double bl_high_2015[NBins2015-NBinsLow2015];
double bl_high_2015_yStatL[NBins2015-NBinsLow2015];
double bl_high_2015_yStatH[NBins2015-NBinsLow2015];
double bl_high_2015_xErrL[NBins2015-NBinsLow2015];
double bl_high_2015_xErrH[NBins2015-NBinsLow2015];
double bl_high_2015_ySystL[NBins2015-NBinsLow2015];
double bl_high_2015_ySystH[NBins2015-NBinsLow2015];


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
		bl_low_ySystL[i]=BPTotalSystDownRatio[i]*BPXsecPPY2DScaled[i];
		bl_low_ySystH[i]=BPTotalSystUpRatio[i]*BPXsecPPY2DScaled[i];
	} 
	else {
		binhigh[i-NBinsLow]=BPXsecPPX[i];
		bl_high[i-NBinsLow]=BPXsecPPY2DScaled[i];
		bl_high_yStatL[i-NBinsLow]=BPXSecPPY2DErrDownScaled[i];
		bl_high_yStatH[i-NBinsLow]=BPXSecPPY2DErrUpScaled[i];
		bl_high_xErrL[i-NBinsLow]=BPXsecPPX[i]-ptBins[i];
		bl_high_xErrH[i-NBinsLow]=ptBins[i + 1]-BPXsecPPX[i];
		bl_high_ySystL[i-NBinsLow]=BPTotalSystDownRatio[i]*BPXsecPPY2DScaled[i];
		bl_high_ySystH[i-NBinsLow]=BPTotalSystUpRatio[i]*BPXsecPPY2DScaled[i];
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
  vector<double> BXsec;
  vector<double> BXsecStat;
  vector<double> BXsecSyst;
  vector<double> BXsec2015;
  vector<double> BXsecStat2015;
  vector<double> BXsecSyst2015;
  vector<double> FONLL;
  vector<double> FONLLUp;
  vector<double> FONLLDown;

 
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

  double RatioBs[NBins];
  double RatioBsStat[NBins];
  double RatioBsSyst[NBins];
  double RatioBsFonErrHigh[NBins];
  double RatioBsFonErrLow[NBins];
  double RatioBs2015[NBins2015];
  double RatioBsStat2015[NBins2015];
  double RatioBsSyst2015[NBins2015];

std::vector<double> Unity(NBins, 1);

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
  // low pT graph
  TGraphAsymmErrors *RatioFonLow = new TGraphAsymmErrors(NBinsLow, BPXsecPPX, Ratio4Y,
                                                         BXSecPPXErrDown, BXSecPPXErrUp,
                                                         Ratio4YErr,Ratio4YErr);
  // high pT graph
  TGraphAsymmErrors *RatioFonHigh = new TGraphAsymmErrors(NBinsHigh, BPXsecPPX + NBinsLow,
                                                          Ratio4Y + NBinsLow,
                                                          BXSecPPXErrDown + NBinsLow,
                                                          BXSecPPXErrUp + NBinsLow,
                                                          Ratio4YErr + NBinsLow,
                                                          Ratio4YErr + NBinsLow);

	

  TGraphAsymmErrors gRatioBs_low(NBinsLow, binlow,RatioBs,bl_low_xErrL, bl_low_xErrH,RatioBsStat, RatioBsStat);
  TGraphAsymmErrors gRatioBs_high(NBinsHigh,binhigh,RatioBs + NBinsLow,bl_high_xErrL, bl_high_xErrH,RatioBsStat + NBinsLow, RatioBsStat + NBinsLow);
  TGraphAsymmErrors gRatioBs_syst_low(NBinsLow,binlow,RatioBs,bl_low_xErrL, bl_low_xErrH,RatioBsSyst, RatioBsSyst);
  TGraphAsymmErrors gRatioBs_syst_high(NBinsHigh,binhigh,RatioBs + NBinsLow,bl_high_xErrL, bl_high_xErrH,RatioBsSyst + NBinsLow, RatioBsSyst + NBinsLow);
  TGraphAsymmErrors gRatioBs_Fon_low(NBinsLow,binlow,Unity.data(),bl_low_xErrL, bl_low_xErrH,RatioBsFonErrLow, RatioBsFonErrHigh);
  TGraphAsymmErrors gRatioBs_Fon_high(NBinsHigh,binhigh,Unity.data(),bl_high_xErrL, bl_high_xErrH,RatioBsFonErrLow + NBinsLow,RatioBsFonErrHigh + NBinsLow);

if(meson_n==0){								  
  TGraphAsymmErrors gRatioBs2015_low(NBinsLow2015,binlow_2015,RatioBs2015,bl_low_2015_xErrL, bl_low_2015_xErrH,RatioBsStat2015, RatioBsStat2015);
  TGraphAsymmErrors gRatioBs2015_syst_low(NBinsLow2015,binlow_2015,RatioBs2015,bl_low_2015_xErrL, bl_low_2015_xErrH,RatioBsSyst2015, RatioBsSyst2015);
  gRatioBs2015_low.SetMarkerStyle(20);
  gRatioBs2015_low.SetMarkerColor(kViolet+1);
  gRatioBs2015_syst_low.SetFillColorAlpha(kViolet+1, 0.5);
  //gRatioBs2015_syst_low.Draw("5");
  //gRatioBs2015_low.Draw("ep");

}

  TGraphAsymmErrors gRatioBs2015_high(NBins2015-NBinsLow2015,binhigh_2015,RatioBs2015 + NBinsLow2015,bl_high_2015_xErrL, bl_high_2015_xErrH,RatioBsStat2015 + NBinsLow2015, RatioBsStat2015 + NBinsLow2015);
  TGraphAsymmErrors gRatioBs2015_syst_high(NBins2015-NBinsLow2015,binhigh_2015,RatioBs2015 + NBinsLow2015,bl_high_2015_xErrL, bl_high_2015_xErrH,RatioBsSyst2015 + NBinsLow2015, RatioBsSyst2015 + NBinsLow2015);
	
  int hcolor = kBlue - 9;
  double halpha = 0.5;
  if (meson_n == 0) { hcolor = kGreen - 9;}

  gRatioBs_low.SetMarkerStyle(25);
  gRatioBs_low.SetMarkerColor(hcolor);
  gRatioBs_low.SetLineColor(hcolor);
  gRatioBs_high.SetMarkerStyle(21);
  gRatioBs_high.SetMarkerColor(hcolor);
  gRatioBs_high.SetLineColor(hcolor);
  gRatioBs2015_high.SetMarkerColor(kOrange+1);
  gRatioBs2015_high.SetLineColor(kOrange+1);
  gRatioBs2015_high.SetMarkerStyle(20);

  TLatex *lat = new TLatex();
  lat->SetNDC();
    lat->SetTextSize(0.1); 
    if (meson_n == 0) {lat->DrawLatex(0.2,0.85 ,Form("B^{+} global Uncertainty: #pm %.1f%%",3.5)) ;}
	 else {lat->DrawLatex(0.2,0.85,Form("B_{s}^{0} Global Uncertainty: #pm %.1f%%",7.7)) ;}

  gRatioBs2015_syst_high.SetFillColorAlpha(kOrange+1, 0.5);
  gRatioBs_syst_low.SetFillColorAlpha(hcolor, halpha);
  gRatioBs_syst_high.SetFillColorAlpha(hcolor, halpha);
  gRatioBs_Fon_low.SetLineColor(kRed-7); 
  gRatioBs_Fon_low.SetFillStyle(0);
  gRatioBs_Fon_high.SetLineColor(kRed+2);
  gRatioBs_Fon_high.SetFillStyle(0);
  gRatioBs_Fon_low.SetLineWidth(2);
  gRatioBs_Fon_high.SetLineWidth(2);
  gRatioBs2015_high.Draw("ep");
  gRatioBs2015_syst_high.Draw("5");
  gRatioBs_low.Draw("ep");
  gRatioBs_high.Draw("ep");
  gRatioBs_syst_low.Draw("5");
  gRatioBs_syst_high.Draw("5");
  gRatioBs_Fon_high.Draw("5");
  gRatioBs_Fon_low.Draw("5");

	RatioFonHigh->SetLineColor(kRed+2);
	RatioFonHigh->SetMarkerStyle(21);
	RatioFonHigh->SetMarkerSize(1);
	RatioFonHigh->SetMarkerColor(kRed+2);
	//RatioFonHigh->Draw("epSAME");
	RatioFonLow->SetLineColor(kRed+2);
	RatioFonLow->SetMarkerStyle(25);
	RatioFonLow->SetMarkerSize(1);
	RatioFonLow->SetMarkerColor(kRed+2);
	//RatioFonLow->Draw("epSAME");
	Unity2->Draw("SAME");
	leg4->AddEntry(RatioFonLow,"2017/FONLL (scaled to |y| < 2.4)","PL");
	leg4->AddEntry(RatioFonHigh,"2017/FONLL","PL");
	//leg4->Draw("same");
	MyPad2->Update();

	//cRatio->SaveAs(Form("Plots/%s/%sCrossComp.png", B_m.Data(), B_m.Data()));
	MyPad1->SetLogy();
	MyPad1->Update();
	cRatio->SaveAs(Form("Plots/%sCrossCompLog.pdf", B_m.Data()));
	//FONLL


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
      BPTotalSystDownRatio[i] << std::setw(columnWidth) <<
      BPTotalSystDownRatio[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      setprecision(3) << BPXsecPPX[i] << "\n";
  }
  out.close();

}


						