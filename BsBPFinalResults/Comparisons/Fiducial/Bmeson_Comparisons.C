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
#include "CMS_lumi.C" 

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

// CREATE THE CANVAS and the pads
	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd(); 
	c->SetLogy();   
	c->SetLeftMargin(0.15);

	//Setup histograms for different purposs
	TH2D * HisEmpty;
	if(meson_n == 0) { 
		HisEmpty = new TH2D("HisEmpty","",100,5,60,100,300.0,2000000);
		HisEmpty->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	} else {
		HisEmpty = new TH2D("HisEmpty","",100,7,50,100,300.0,2000000);
		HisEmpty->GetXaxis()->SetTitle("p_{T} [GeV/c]");
		}
	HisEmpty->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb c/GeV]");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	//HisEmpty->GetYaxis()->SetTitleOffset(1.8);
	//HisEmpty->GetXaxis()->SetTitleOffset(1.3);	

	TH2D * HisEmpty2;
	if (meson_n == 0){
		HisEmpty2 = new TH2D("HisEmpty2","",100,5,60,100,300.0,3000000);
		HisEmpty2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	} else {	
		HisEmpty2 = new TH2D("HisEmpty2","",100,7,50,100,300.0,3000000);
		HisEmpty2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	}
		HisEmpty2->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb c/GeV]");
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
		vector<float>  BPYSystDown_low ;
		vector<float>  BPYSystDown_high ;
		vector<float>  BPYSystUp_low ;
		vector<float>  BPYSystUp_high ;

	if(meson_n == 0) { 
			BPXsecPPXLow =  {BPXsecPPX[0],BPXsecPPX[1]};
			BPXsecPPXHigh = {BPXsecPPX[2], BPXsecPPX[3], BPXsecPPX[4], BPXsecPPX[5], BPXsecPPX[6]};
				BPXsecPPXErrLow = {BXSecPPXErrDown[0],BXSecPPXErrDown[1]};
				BPXsecPPXErrHigh = {BXSecPPXErrDown[2], BXSecPPXErrDown[3], BXSecPPXErrDown[4], BXSecPPXErrDown[5], BXSecPPXErrDown[6]};

			BPXsecPPYLow = {BPXsecPPY2D[0],BPXsecPPY2D[1]};
			BPXsecPPYHigh = {BPXsecPPY2D[2],BPXsecPPY2D[3],BPXsecPPY2D[4],BPXsecPPY2D[5],BPXsecPPY2D[6]};
				BPXsecPPYErrDownLow = {BPXSecPPY2DErrDown[0],BPXSecPPY2DErrDown[1]};
				BPXsecPPYErrDownHigh = {BPXSecPPY2DErrDown[2], BPXSecPPY2DErrDown[3],BPXSecPPY2DErrDown[4], BPXSecPPY2DErrDown[5], BPXSecPPY2DErrDown[6]};
				BPXsecPPYErrUpLow= {BPXSecPPY2DErrUp[0],BPXSecPPY2DErrUp[1]};
				BPXsecPPYErrUpHigh = {BPXSecPPY2DErrUp[2], BPXSecPPY2DErrUp[3],BPXSecPPY2DErrUp[4], BPXSecPPY2DErrUp[5], BPXSecPPY2DErrUp[6]};
				BPYSystDown_low = {BPXSecPPYSystDown[0],BPXSecPPYSystDown[1]};
				BPYSystDown_high = {BPXSecPPYSystDown[2],BPXSecPPYSystDown[3],BPXSecPPYSystDown[4],BPXSecPPYSystDown[5],BPXSecPPYSystDown[6]};
				BPYSystUp_low = {BPXSecPPYSystUp[0],BPXSecPPYSystUp[1]};
				BPYSystUp_high = {BPXSecPPYSystUp[2],BPXSecPPYSystUp[3],BPXSecPPYSystUp[4],BPXSecPPYSystUp[5],BPXSecPPYSystUp[6]};
	} else {
			BPXsecPPXLow =  {BPXsecPPX[0]};
			BPXsecPPXHigh = {BPXsecPPX[1], BPXsecPPX[2], BPXsecPPX[3]};
				BPXsecPPXErrLow = {BXSecPPXErrDown[0]};
				BPXsecPPXErrHigh = {BXSecPPXErrDown[1], BXSecPPXErrDown[2], BXSecPPXErrDown[3]};

			BPXsecPPYLow = {BPXsecPPY2D[0]};
			BPXsecPPYHigh = {BPXsecPPY2D[1],BPXsecPPY2D[2],BPXsecPPY2D[3]};
				BPXsecPPYErrDownLow = {BPXSecPPY2DErrDown[0]};
				BPXsecPPYErrDownHigh = {BPXSecPPY2DErrDown[1], BPXSecPPY2DErrDown[2], BPXSecPPY2DErrDown[3] };
				BPXsecPPYErrUpLow= {BPXSecPPY2DErrUp[0]};
				BPXsecPPYErrUpHigh = {BPXSecPPY2DErrUp[1], BPXSecPPY2DErrUp[2], BPXSecPPY2DErrUp[3]};
				BPYSystDown_low = {BPXSecPPYSystDown[0]};
				BPYSystDown_high = {BPXSecPPYSystDown[1],BPXSecPPYSystDown[2],BPXSecPPYSystDown[3]};
				BPYSystUp_low = {BPXSecPPYSystUp[0]};
				BPYSystUp_high = {BPXSecPPYSystUp[1],BPXSecPPYSystUp[2],BPXSecPPYSystUp[3]};
		}



	float zero1[1] = {0};
	float zero2[2] = {0,0};


	TGraphAsymmErrors *BPRAAGraph_low_just_marker ;
	if(meson_n != 0) { BPRAAGraph_low_just_marker = new TGraphAsymmErrors(NBinsLow, BPXsecPPXLow.data(), BPXsecPPYLow.data() ,zero1, zero1, zero1, zero1);} 
	else {BPRAAGraph_low_just_marker = new TGraphAsymmErrors(NBinsLow, BPXsecPPXLow.data(), BPXsecPPYLow.data() ,zero2, zero2, zero2, zero2);}
	
	TGraphAsymmErrors *BPRAAGraph_low = new TGraphAsymmErrors(NBinsLow , BPXsecPPXLow.data() , BPXsecPPYLow.data() , BPXsecPPXErrLow.data() , BPXsecPPXErrLow.data() , BPXsecPPYErrDownLow.data() , BPXsecPPYErrUpLow.data());
	TGraphAsymmErrors *BPRAAGraph     = new TGraphAsymmErrors(NBinsHigh, BPXsecPPXHigh.data(), BPXsecPPYHigh.data(), BPXsecPPXErrHigh.data(), BPXsecPPXErrHigh.data(), BPXsecPPYErrDownHigh.data(), BPXsecPPYErrUpHigh.data());     
	TGraphAsymmErrors *BPRAAGraphSyst_low  = new TGraphAsymmErrors(NBinsLow , BPXsecPPXLow.data() , BPXsecPPYLow.data() , BPXsecPPXErrLow.data() , BPXsecPPXErrLow.data() , BPYSystDown_low.data() , BPYSystUp_low.data());                 											
	TGraphAsymmErrors *BPRAAGraphSyst      = new TGraphAsymmErrors(NBinsHigh, BPXsecPPXHigh.data(), BPXsecPPYHigh.data(), BPXsecPPXErrHigh.data(), BPXsecPPXErrHigh.data(), BPYSystDown_high.data(), BPYSystUp_high.data());                 											
  // separate plots for different fiducial regions
 
  	cout << endl << "-------------------------------------------------------  "<< Form("%s meson Xsection", B_m.Data()) <<"  -------------------------------------------------------" << endl;

	for(int i=0;i<NBins;i++){		
		cout << "BIN " << Form("[%.1f,%.1f]  ",ptBins[i],ptBins[i+1]) << Form("%f #pm (STATup) %.1f #pm (SYSTup) %.1f #pm %.1f (STATdown) %.1f #pm (SYSTdown) %.1f ",BPXsecPPY2D[i],BPXSecPPY2DErrUp[i],BPXSecPPYSystUp[i],BPXSecPPY2DErrDown[i],BPXSecPPYSystDown[i]) << endl;
		cout << "(normalized) BIN " << Form("[%.1f,%.1f]  ",ptBins[i],ptBins[i+1]) << Form("%.1f #pm (STATup) %.1f #pm (SYSTup) %.1f #pm %.1f (STATdown) %.1f #pm (SYSTdown) %.1f ",BPXsecPPY2D[i],BPXSecPPY2DErrUp[i]/BPXsecPPY2D[i],BPXSecPPYSystUp[i]/BPXsecPPY2D[i],BPXSecPPY2DErrDown[i]/BPXsecPPY2D[i],BPXSecPPYSystDown[i]/BPXsecPPY2D[i]) << endl;
	}
 
 	cout<< endl << "-------------------------------------------------------  "<< Form("%s meson Xsection", B_m.Data()) <<"  -------------------------------------------------------" << endl;



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
	
	HisEmpty->Draw();
	BPRAAGraphSyst_low->Draw("5same");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph_low->Draw("epSAME");
	BPRAAGraph_low_just_marker->Draw("epSAME");
	c->SetLogy();

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.02); 
	lat->SetTextFont(42);
	if (meson_n == 0) {lat->DrawLatex(0.6,0.7 ,Form("2017 pp global Unc. #pm %.1f%%",3.5));} 
	else {	lat->DrawLatex(0.6,0.7,Form("2017 pp Global Unc. #pm %.1f%%",7.7)) ;}
	
    TLegend* leged = new TLegend(0.65,0.77,0.9,0.85,NULL,"brNDC");
	leged->SetBorderSize(0);
	leged->SetFillStyle(0);
	leged->AddEntry((TObject*)0, "y region:", "");
	leged->AddEntry(BPRAAGraph,"|y|<2.4","P");
	leged->AddEntry(BPRAAGraph_low,"|y|>1.5","P");
	leged->SetTextSize(0.025);
	leged->Draw("same");

	c->SaveAs(Form("Plots/%sCrossONLYLog.pdf", B_m.Data()));
// CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection CrossSection 


//  XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb XSEC vs PbPb 

		float sys_up_pbpb[4] ={0,0,0,0};
		float sys_down_pbpb[4] ={0,0,0,0};
		float sat_up_pbpb[4] ={0,0,0,0};
		float sat_down_pbpb[4] ={0,0,0,0};
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

			c->SaveAs(Form("Plots/%sPbPbPPCrossLog.pdf", B_m.Data()));
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

			c->SaveAs(Form("Plots/%spp_2015_CrossLog.pdf", B_m.Data()));
//2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 2015 Reference 























// vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL vs FONL 

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

  //DATA
  TGraphAsymmErrors gRatioBs_low(NBinsLow, binlow, RatioBs, bl_low_xErrL, bl_low_xErrH, RatioBsStat, RatioBsStat);
  TGraphAsymmErrors *BPRAAGraph_low_just_m ;
  if(meson_n != 0) { BPRAAGraph_low_just_m= new TGraphAsymmErrors(NBinsLow, binlow, RatioBs ,zero1, zero1, zero1, zero1);} 
  else {BPRAAGraph_low_just_m             = new TGraphAsymmErrors(NBinsLow, binlow, RatioBs ,zero2, zero2, zero2, zero2);}
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
	CMS_lumi(MyPadr,19011,0);
	MyPadr->Update();
  
  	cr->SetLogy();   
	cr->SaveAs(Form("Plots/%sCrossCompLog.pdf", B_m.Data()));
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


						