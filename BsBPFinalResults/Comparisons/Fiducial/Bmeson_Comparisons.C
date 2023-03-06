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

//MOST OF THE VARIABLES SAY BP BUT ARE WORKING FOR Bs AS WELL 
// THIS CODE CAN BE EASILY EXTENDED TO ACCOUNT FOR A FUTURE MESON... 
//JUST FOLLOW THE LOGIC

void Bmeson_Comparisons(int meson_n ){

	TString B_m ;
	int NBins = 7;
    vector<double> scaledPt;
	TString t_tree ;
	TString b_m;

	int NBinsLow ;
  	int NBinsHigh ;
	int NBins2015 ;
	if(meson_n == 0){
		NBins = nptBinsBP;
		B_m = "BP";
		b_m = "bp";
		scaledPt = {5, 7, 10};
		t_tree = "ntKp";
		NBinsLow = 2 ;
		NBinsHigh = 5;
		NBins2015 = 5;
	} else {
		NBins = nptBins;
		B_m = "Bs";
		b_m = "bs";
		scaledPt = {7, 10};
		t_tree = "ntphi";
		NBinsLow = 1 ;
		NBinsHigh = 3;
		NBins2015 = 3;
	}

	gSystem->mkdir(Form("Plots/%s", B_m.Data()), true);
	TString InfileB = Form("../../../%s/EffAna/FinalFiles/%sPPCorrYieldPT.root",B_m.Data(),B_m.Data());
	TFile * FileB= new TFile(InfileB.Data());
	double ptBins[NBins+1];
	float BPXSecPPYErrUp[NBins];
	float BPXSecPPYErrDown[NBins];

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
  for (auto i = 0; i < factor.size(); ++i) {
    cout << "applying scaling factor: " << factor[i] << "\n";
    BPXsecPPY2DScaled[i] *= factor[i];
	BPXSecPPY2DErrUpScaled[i] *= factor[i];
    BPXSecPPY2DErrDownScaled[i] *= factor[i];
  		}

float BXSecPPXErrUp[NBins] ;
float BXSecPPXErrDown[NBins] ;

if (meson_n == 0){
	vector<float> vect_BXSecPPXErrUp{1,1.5,2.5,2.5,5,10,5};
	vector<float> vect_BXSecPPXErrDown{1,1.5,2.5,2.5,5,10,5};
	for( int c=0; c <NBins; c++){ 
		BXSecPPXErrUp[c]=vect_BXSecPPXErrUp[c];
		BXSecPPXErrDown[c]=vect_BXSecPPXErrDown[c];}
} else {
	vector<float> vect_BXSecPPXErrUp {1.5,2.5,2.5,15};
	vector<float> vect_BXSecPPXErrDown {1.5,2.5,2.5,15};
	for( int c=0; c <NBins; c++){ 
		BXSecPPXErrUp[c]=vect_BXSecPPXErrUp[c];
		BXSecPPXErrDown[c]=vect_BXSecPPXErrDown[c];}
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
		HisEmpty = new TH2D("HisEmpty","",100,5,60,100,100.0,2000000);
		HisEmpty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	} else {
		HisEmpty = new TH2D("HisEmpty","",100,7,50,100,100.0,2000000);
		HisEmpty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		}
	HisEmpty->GetYaxis()->SetTitle("d#sigma/dp_{T} (pb c/GeV)");
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

	// CrossSection 
	BPPPCrossGraph2D->SetLineColor(kBlue+2);
	BPPPCrossGraph2D->SetMarkerStyle(21);
	BPPPCrossGraph2D->SetMarkerSize(1);
	BPPPCrossGraph2D->SetMarkerColor(kBlue+2);
	BPPPCrossGraph2DSyst->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraph2DSyst->SetLineColor(kBlue-9);

	TLegend* leged = new TLegend(0.75,0.70,0.95,0.95,NULL,"brNDC");
	leged->SetBorderSize(0);
	leged->SetTextSize(0.025);     
	leged->SetTextFont(42);
	leged->SetFillStyle(0);
	if(meson_n == 0) {leged->AddEntry(BPPPCrossGraph2D,"B^{#pm}","PL");}
	else {leged->AddEntry(BPPPCrossGraph2D,"B^{0}_{s}","PL");}
	leged->Draw();
	BPPPCrossGraph2D->Draw("ep");	
	BPPPCrossGraph2DSyst->Draw("5same");	
	//c->SaveAs(Form("Plots/%s/%sCrossONLY.png",B_m.Data(), B_m.Data()));
	c->SetLogy();
	c->SaveAs(Form("Plots/%s/%sCrossONLYLog.pdf", B_m.Data(), B_m.Data()));
	// CrossSection (log scale) 

// FOR BP ONLY (for now), PbPb plots for comparison
 	if (meson_n == 0){
	TGraphAsymmErrors *BPPbPbCrossGraph;
  	TGraphAsymmErrors *BPPbPbCrossGraphSyst;
	BPPbPbCrossGraph = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY,BPXSecPbPbXErrDown, BPXSecPbPbXErrUp,BPXSecPbPbYErrDown,BPXSecPbPbYErrUp);
  	BPPbPbCrossGraphSyst  = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY, BPXSecPbPbXErrDown, BPXSecPbPbXErrUp, BPXSecPbPbYSystDown,BPXSecPbPbYSystUp);
	BPPbPbCrossGraph->SetLineColor(kGreen+2);
	BPPbPbCrossGraph->SetMarkerStyle(21);
	BPPbPbCrossGraph->SetMarkerSize(1);
	BPPbPbCrossGraph->SetMarkerColor(kGreen+2);
	BPPbPbCrossGraphSyst->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraphSyst->SetLineColor(kGreen-9);
	
	TCanvas * c2New = new TCanvas("c2New","c2New",700,700);
	c2New->cd();
	c2New->SetLeftMargin(0.15);
	HisEmpty->Draw();

	TLegend* leg = new TLegend(0.55,0.8,0.9,0.9,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.025);     
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(BPPbPbCrossGraph,"2018 PbPb 5.02 TeV","PL");
	leg->AddEntry(BPPPCrossGraph2D,"2017 pp 5.02 TeV","PL");
	leg->Draw();
	BPPbPbCrossGraph->Draw("ep");	
	BPPPCrossGraph2D->Draw("ep");	
	BPPbPbCrossGraphSyst->Draw("5same");	
	BPPPCrossGraph2DSyst->Draw("5same");

	//c2New->SaveAs("Plots/BP/BPPbPbPPCross.png");
	c2New->SetLogy();
	c2New->SaveAs(Form("Plots/%s/%sPbPbPPCrossLog.pdf", B_m.Data(), B_m.Data()));
	}
 // FOR BP ONLY (for now), PbPb plots for comparison

	//2015 Reference (	Big plots  )
	TCanvas * cRatio = new TCanvas("cRatio","cRatio",700,800);
	TPad * MyPad1;
	MyPad1 = new TPad("MyPad1","",0,0.4,1,1.0);
	MyPad1->Draw();
	TPad * MyPad2;
	MyPad2 = new TPad("MyPad2","",0,0.2,1,0.4);
	MyPad2->Draw();
	TPad * MyPad3;
	MyPad3 = new TPad("MyPad3","",0,0.00,1,0.2);
	MyPad3->Draw();

	MyPad1->cd();
	TH2D * HisEmpty2;
	if (meson_n == 0){
		HisEmpty2 = new TH2D("HisEmpty2","",100,5,60,100,100.0,30000000);
		HisEmpty2->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");}
	else {	
		HisEmpty2 = new TH2D("HisEmpty2","",100,7,50,100,100.0,30000000);
		HisEmpty2->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");}
	HisEmpty2->GetYaxis()->SetTitle("d#sigma/d p_{T} (pb c/GeV)");
	HisEmpty2->GetXaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->CenterTitle();
	//HisEmpty2->SetTitle("B^{+} Cross Section With Fiducial Region");
	HisEmpty2->GetYaxis()->SetTitleOffset(1.2);
	HisEmpty2->GetXaxis()->SetTitleOffset(1.2);
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
			vector<float> vect_BXsecPPX2015{8.5,12.5,17.5,25,40};
			vector<float> vect_BXSecPPXErrDown2015{1.5,2.5,2.5,5,10};
			vector<float> vect_BXSecPPXErrUp2015{1.5,2.5,2.5,5,10};
			vector<float> vect_BXsecPPY2015{2610000,744000,197000,46500,5300};
			vector<float> vect_BXSecPPYErrDown2015{170000,29000,9000,2400,500};
			vector<float> vect_BXSecPPYErrUp2015{170000,29000,9000,2400,500};
			vector<float> vect_BXSecPPYSystDown2015{230000,59000,15000,3500,400};
			vector<float> vect_BXSecPPYSystUp2015{230000,59000,15000,3500,400};
			for( int c=0; c <NBins; c++){ 
				BXsecPPX2015[c]= vect_BXsecPPX2015[c] ;
				BXSecPPXErrDown2015[c]= vect_BXSecPPXErrDown2015[c];
				BXSecPPXErrUp2015[c]= vect_BXSecPPXErrUp2015[c];
				BXsecPPY2015[c]= vect_BXsecPPY2015[c];
				BXSecPPYErrDown2015[c]= vect_BXSecPPYErrDown2015[c];
				BXSecPPYErrUp2015[c]= vect_BXSecPPYErrUp2015[c];
				BXSecPPYSystDown2015[c]= vect_BXSecPPYSystDown2015[c];
				BXSecPPYSystUp2015[c] = vect_BXSecPPYSystUp2015[c];
				}
		} else {
			vector<float> vect_BXsecPPX2015{11,17.5,35.0};
			vector<float> vect_BXSecPPXErrDown2015{4,2.5,15};
			vector<float> vect_BXSecPPXErrUp2015{4,2.5,15};
			vector<float> vect_BXsecPPY2015{316000,34100,3830};
			vector<float> vect_BXSecPPYErrDown2015{37000,6300,670};
			vector<float> vect_BXSecPPYErrUp2015{37000,6300,670};
			vector<float> vect_BXSecPPYSystDown2015{62000,3200,360};
			vector<float> vect_BXSecPPYSystUp2015{62000,3200,360};
			for( int c=0; c <NBins; c++){ 
				BXsecPPX2015[c]= vect_BXsecPPX2015[c] ;
				BXSecPPXErrDown2015[c]= vect_BXSecPPXErrDown2015[c];
				BXSecPPXErrUp2015[c]= vect_BXSecPPXErrUp2015[c];
				BXsecPPY2015[c]= vect_BXsecPPY2015[c];
				BXSecPPYErrDown2015[c]= vect_BXSecPPYErrDown2015[c];
				BXSecPPYErrUp2015[c]= vect_BXSecPPYErrUp2015[c];
				BXSecPPYSystDown2015[c]= vect_BXSecPPYSystDown2015[c];
				BXSecPPYSystUp2015[c] = vect_BXSecPPYSystUp2015[c];
				}			
			}

	TGraphAsymmErrors *BPPPCrossGraph2015 = new TGraphAsymmErrors(NBins2015, BXsecPPX2015, BXsecPPY2015,BXSecPPXErrDown2015, BXSecPPXErrUp2015,BXSecPPYErrDown2015,BXSecPPYErrUp2015);
	TGraphAsymmErrors *BPPPCrossGraph2015Syst = new TGraphAsymmErrors(NBins2015, BXsecPPX2015, BXsecPPY2015,BXSecPPXErrDown2015, BXSecPPXErrUp2015,BXSecPPYSystDown2015,BXSecPPYSystUp2015);

	BPPPCrossGraph2015Syst->SetFillColorAlpha(kGreen-9+2,0.5);
	BPPPCrossGraph2015Syst->SetLineColor(kGreen-9+2);
	BPPPCrossGraph2015->SetLineColor(kGreen+2);
	BPPPCrossGraph2015->SetMarkerStyle(33);
	BPPPCrossGraph2015->SetMarkerSize(1);
	BPPPCrossGraph2015->SetMarkerColor(kGreen+2);
	BPPPCrossGraph2DLow->SetLineColor(kOrange+1);
	BPPPCrossGraph2DLow->SetMarkerStyle(25);
	BPPPCrossGraph2DLow->SetMarkerSize(1);
	BPPPCrossGraph2DLow->SetMarkerColor(kOrange+1);
	BPPPCrossGraph2DHigh->SetLineColor(kOrange+1);
	BPPPCrossGraph2DHigh->SetMarkerStyle(34);
	BPPPCrossGraph2DHigh->SetMarkerSize(1);
	BPPPCrossGraph2DHigh->SetMarkerColor(kOrange+1);
	BPPPCrossGraph2DScaledSyst->SetFillColorAlpha(kOrange+1, 0.3);
	BPPPCrossGraph2DScaledSyst->SetLineColor(kOrange+1);

    TFile * finFONLL ;
	if(meson_n == 0){ finFONLL = new TFile("FONLLs/fonllOutput_pp_Bplus_5p03TeV_y2p4.root");}
	else{ finFONLL = new TFile("FONLLs/BsFONLL.root");}
	finFONLL->cd();
	TGraphAsymmErrors *BPFONLL = (TGraphAsymmErrors*) finFONLL->Get("gaeSigmaBplus");
	BPFONLL->SetLineColor(kRed+2);
	BPFONLL->SetMarkerStyle(20);
	BPFONLL->SetMarkerSize(1);
	BPFONLL->SetMarkerColor(kRed+2);
	BPFONLL->Draw("epSAME");
	
	BPPPCrossGraph2015->Draw("epSAME");
	BPPPCrossGraph2DScaledSyst->Draw("5SAME");
	BPPPCrossGraph2015Syst->Draw("5same");	
	BPPPCrossGraph2DLow->Draw("epSAME");
	BPPPCrossGraph2DHigh->Draw("epSAME");

	TLegend* leg3 = new TLegend(0.55,0.64,0.8,0.85,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.025);     
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->AddEntry(BPPPCrossGraph2DLow,"2017 pp 5.02 TeV (scaled to |y| < 2.4)","PL");	
	leg3->AddEntry(BPPPCrossGraph2DHigh,"2017 pp 5.02 TeV","PL");	
	leg3->AddEntry(BPPPCrossGraph2015,"2015 pp 5.02 TeV","PL");
	leg3->AddEntry(BPFONLL,"FONLL Calculations","PL");
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

TH2D * HisEmpty3;
if (meson_n == 0){
	HisEmpty3 = new TH2D("HisEmpty3","",100,5,60,100,0,2);
	HisEmpty3->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
} else {
	HisEmpty3 = new TH2D("HisEmpty3","",100,7,50,100,0,2);
	HisEmpty3->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
}
	HisEmpty3->GetYaxis()->SetTitle("2017/2015 Data");
	HisEmpty3->GetXaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->SetTitleOffset(1.2);
	HisEmpty3->Draw();

	TGraphAsymmErrors *Ratio2 = new TGraphAsymmErrors(NBins2015, BXsecPPX2015, Ratio2Y,BXSecPPXErrDown2015, BXSecPPXErrUp2015,Ratio2YErr,Ratio2YErr);
		Ratio2->SetLineColor(kOrange+1);
		Ratio2->SetMarkerStyle(34);
		Ratio2->SetMarkerSize(1);
		Ratio2->SetMarkerColor(kOrange+1);
		Ratio2->Draw("epSAME");

  	TGraphAsymmErrors *RatioDataLow;
  	if (meson_n ==0){
		RatioDataLow = new TGraphAsymmErrors(1, BPXsecPPXLow.data() + 1, RatioDataYLow.data(),BPXsecPPXErrDownLow.data() + 1, BPXsecPPXErrUpLow.data() + 1, RatioDataYLowErr.data(), RatioDataYLowErr.data());
		RatioDataLow->SetLineColor(kOrange+1);
		RatioDataLow->SetMarkerStyle(25);
		RatioDataLow->SetMarkerSize(1);
		RatioDataLow->SetMarkerColor(kOrange+1);
		RatioDataLow->Draw("epSAME");
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
		Ratio4Y[i] = BPXsecPPY2DScaled[i]/FONLLY[i];
		Ratio4YErr[i] = Ratio4Y[i] *TMath::Sqrt(BPXSecPPY2DErrDown[i]/BPXsecPPY2D[i] * BPXSecPPY2DErrDown[i]/BPXsecPPY2D[i] + FONLLYErr[i]/FONLLY[i] * FONLLYErr[i]/FONLLY[i] );
		}

	cRatio->cd();
	MyPad3->cd();

TH2D * HisEmpty4;
if (meson_n == 0){
	HisEmpty4 = new TH2D("HisEmpty4","",100,5,60,100,0,2);
	HisEmpty4->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
} else {
	HisEmpty4 = new TH2D("HisEmpty4","",100,7,50,100,0,2);
	HisEmpty4->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
}
	HisEmpty4->GetYaxis()->SetTitle("2017 Data/FONLL");
	HisEmpty4->GetXaxis()->CenterTitle();
	HisEmpty4->GetYaxis()->CenterTitle();
	HisEmpty4->GetYaxis()->SetTitleOffset(1.2);
	HisEmpty4->Draw();

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


	RatioFonHigh->SetLineColor(kOrange+1);
	RatioFonHigh->SetMarkerStyle(34);
	RatioFonHigh->SetMarkerSize(1);
	RatioFonHigh->SetMarkerColor(kOrange+1);
	RatioFonHigh->Draw("epSAME");
	RatioFonLow->SetLineColor(kOrange+1);
	RatioFonLow->SetMarkerStyle(25);
	RatioFonLow->SetMarkerSize(1);
	RatioFonLow->SetMarkerColor(kOrange+1);
	RatioFonLow->Draw("epSAME");
	Unity2->Draw("SAME");
	MyPad2->Update();

	//cRatio->SaveAs(Form("Plots/%s/%sCrossComp.png", B_m.Data(), B_m.Data()));
	MyPad1->SetLogy();
	MyPad1->Update();
	cRatio->SaveAs(Form("Plots/%s/%sCrossCompLog.pdf", B_m.Data(), B_m.Data()));
	//FONLL


  // summary of errors (in ratio, not percent)
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

