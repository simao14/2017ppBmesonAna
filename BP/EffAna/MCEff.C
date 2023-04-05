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
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "../../parameter.h"


using namespace std;

using std::cout;
using std::endl;

bool reweightPtOnY = true;


void  MCEff(int DoTnP, int Rescale){
	gSystem->mkdir( "Syst" , true);
	gSystem->mkdir( "NewEff2DMaps" , true);
	gSystem->mkdir( "1DEffPlots" , true);
	gSystem->mkdir( "TnPHis" , true);
	gSystem->mkdir( "MuonInfoPlots" , true);
	gSystem->mkdir( "Eff2DMapTnP" , true);
	gSystem->mkdir( "Plot1DEfficiency/Pt" , true);
	gSystem->mkdir( "Plot1DEfficiency/Mult" , true);

	gStyle->SetOptStat(0);

	int ptmin = 10;
	int ptmax = 50;

	TString infile;
	infile = "/data3/tasheng/presel/output/BP_MC_BDTs_nom_tnp.root";

	TFile * fin = new TFile(infile.Data());

	fin->cd();

	TTree * ntKp = (TTree * ) fin->Get("Bfinder/ntKp");
	//	TTree * BDT = (TTree * ) fin->Get("BDT");
	TTree * ntHi = (TTree * ) fin->Get("hiEvtAnalyzer/HiTree");
	TTree * ntSkim = (TTree * ) fin->Get("skimanalysis/HltTree");
	TTree * ntHlt = (TTree *) fin->Get("hltanalysis/HltTree");
	//	TTree * TnPInfo = (TTree * ) fin->Get("TnPInfo");
	//	TTree * CentWeightTree =	(TTree * ) fin->Get("CentWeightTree");
	TTree * ntGen = (TTree * ) fin->Get("Bfinder/ntGen");
	
//	TString BDT1Name = "BDT_pt_3_5";
	TString BDT2Name = "BDT_pt_5_7";
	TString BDT3Name = "BDT_pt_7_10";
	TString BDT4Name = "BDT_pt_10_15";
	TString BDT5Name = "BDT_pt_15_20";
	TString BDT6Name = "BDT_pt_20_50";
//	TString BDT7Name = "BDT_pt_2_3";
//	TString BDT8Name = "BDT_pt_1_2";

	if(Rescale == 1){
//	 BDT1Name = "BDT_pt_New_3_5";
	 BDT2Name = "BDT_pt_New_5_7";
	 BDT3Name = "BDT_pt_New_7_10";
	 BDT4Name = "BDT_pt_New_10_15";
	 BDT5Name = "BDT_pt_New_15_20";
	 BDT6Name = "BDT_pt_New_20_50";
//	 BDT7Name = "BDT_pt_New_2_3";
//	 BDT8Name = "BDT_pt_New_1_2";

	}

	
//	TTree * BDT1 = (TTree *) fin->Get(BDT1Name.Data());
	TTree * BDT2 = (TTree *) fin->Get(BDT2Name.Data());
	TTree * BDT3 = (TTree *) fin->Get(BDT3Name.Data());
	TTree * BDT4 = (TTree *) fin->Get(BDT4Name.Data());
	TTree * BDT5 = (TTree *) fin->Get(BDT5Name.Data());
	TTree * BDT6 = (TTree *) fin->Get(BDT6Name.Data());
//	TTree * BDT7 = (TTree *) fin->Get("BDT_pt_2_3");
//	TTree * BDT7 = (TTree *) fin->Get(BDT7Name.Data());
//	TTree * BDT8 = (TTree *) fin->Get(BDT8Name.Data());
	


	TTree * root = (TTree * ) fin->Get("Bfinder/root");
	TTree * TnPInfo = (TTree *) fin->Get("TnPInfo");

	Int_t nMult;

	root->SetBranchAddress("EvtInfo.nMult",&nMult);

	const int NCand = 13000;


	int run;
	int lumi;
	int evt;
	int hiBin;
	Float_t PVz;
	Int_t pclusterCompatibilityFilter;
	Int_t pprimaryVertexFilter;
	Int_t phfCoincFilter2Th4;
	Int_t   Bsize;
	Float_t Btrk1Pt[NCand];
	Float_t Btrk2Pt[NCand];
	Float_t Btrk1PtErr[NCand];
	Float_t Btrk2PtErr[NCand];
	Float_t Bchi2cl[NCand];
	Float_t BsvpvDistance[NCand];
	Float_t BsvpvDisErr[NCand];
	Float_t Bpt[NCand];
	Float_t Btrk1Eta[NCand];
	Float_t Btrk2Eta[NCand];
	Float_t By[NCand];
	Bool_t Bmu1isTriggered[NCand];
	Bool_t Bmu2isTriggered[NCand];
	Float_t Bmass[NCand];
	Float_t Bmumumass[NCand];
	Float_t Bmu1eta[NCand];
	Float_t Bmu1pt[NCand];
	Float_t Bmu2eta[NCand];
	Float_t Bmu2pt[NCand];

	//	Float_t Bmu1phi[NCand];
	//	Float_t Bmu2phi[NCand];

	Bool_t Bmu1TMOneStationTight[NCand];
	Int_t Bmu1InPixelLayer[NCand];
	Int_t Bmu1InStripLayer[NCand];
	Bool_t Bmu2TMOneStationTight[NCand];	
	Int_t Bmu2InPixelLayer[NCand];
	Int_t Bmu2InStripLayer[NCand];
	Bool_t Bmu1isGlobalMuon[NCand];
	Bool_t Bmu2isGlobalMuon[NCand];
	Bool_t Bmu1isTrackerMuon[NCand];
	Bool_t Bmu2isTrackerMuon[NCand];
	Float_t Bmu1dxyPV[NCand];
	Float_t Bmu2dxyPV[NCand];
	Float_t Bmu1dzPV[NCand];
	Float_t Bmu2dzPV[NCand];
	Bool_t Btrk1highPurity[NCand];
	Bool_t Btrk2highPurity[NCand];
	Float_t Btktkmass[NCand];
	Float_t Btrk1PixelHit[NCand];
	Float_t Btrk2PixelHit[NCand];
	Float_t Btrk1StripHit[NCand];
	Float_t Btrk2StripHit[NCand];
	Float_t Btrk1Chi2ndf[NCand];
	Float_t Btrk2Chi2ndf[NCand];
	Float_t Btrk1nStripLayer[NCand];
	Float_t Btrk2nStripLayer[NCand];
	Float_t Btrk1nPixelLayer[NCand];
	Float_t Btrk2nPixelLayer[NCand];
	Float_t Bgen[NCand];
	Float_t Bgenpt[NCand];
	Float_t Bgeny[NCand];


	//	Float_t pthatweight;

	Float_t pthat;
	Float_t weight;
	Float_t Bdtheta[NCand];
//	Double_t BDT_pt_3_5[NCand];
	Double_t BDT_pt_5_7[NCand];
	Double_t BDT_pt_7_10[NCand];
	Double_t BDT_pt_10_15[NCand];
	Double_t BDT_pt_15_20[NCand];
	Double_t BDT_pt_20_50[NCand];
//	Double_t BDT_pt_2_3[NCand];
//	Double_t BDT_pt_2_3[NCand];
//	Double_t BDT_pt_1_2[NCand];
	
	/*

	   BDT->SetBranchAddress("BDT_5_7",BDT_pt_5_7);
	   BDT->SetBranchAddress("BDT_7_10",BDT_pt_7_10);
	   BDT->SetBranchAddress("BDT_10_15",BDT_pt_10_15);
	   BDT->SetBranchAddress("BDT_15_20",BDT_pt_15_20);
	   BDT->SetBranchAddress("BDT_22_30",BDT_pt_22_30);
	   BDT->SetBranchAddress("BDT_30_40",BDT_pt_30_40);
	   BDT->SetBranchAddress("BDT_40_50",BDT_pt_40_50);
	   BDT->SetBranchAddress("BDT_50_60",BDT_pt_50_60);

	   BDT->SetBranchAddress("run",&run);
	   BDT->SetBranchAddress("evt",&evt);
	   BDT->SetBranchAddress("lumi",&lumi);

*/

	ntHi->SetBranchAddress("hiBin",&hiBin);
	ntHi->SetBranchAddress("pthat",&pthat);
	ntHi->SetBranchAddress("weight",&weight);

	int HBHENoiseFilterResult;
	int pPAprimaryVertexFilter;
	int pBeamScrapingFilter;

	ntSkim->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
	ntSkim->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
	ntSkim->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
	ntKp->SetBranchAddress("Bsize",&Bsize);
	ntKp->SetBranchAddress("PVz",&PVz);
	ntKp->SetBranchAddress("Btrk1Pt",Btrk1Pt);
	//	ntKp->SetBranchAddress("Btrk2Pt",Btrk2Pt);
	ntKp->SetBranchAddress("Btrk1PtErr",Btrk1PtErr);
	//	ntKp->SetBranchAddress("Btrk2PtErr",Btrk2PtErr);
	ntKp->SetBranchAddress("Bchi2cl",Bchi2cl);
	ntKp->SetBranchAddress("BsvpvDistance",BsvpvDistance);
	ntKp->SetBranchAddress("BsvpvDisErr",BsvpvDisErr);
	ntKp->SetBranchAddress("Bpt",Bpt);
	ntKp->SetBranchAddress("By",By);
	ntKp->SetBranchAddress("Btrk1Eta",Btrk1Eta);
	//ntKp->SetBranchAddress("Btrk2Eta",Btrk2Eta);
	ntKp->SetBranchAddress("Bmass",Bmass);
	ntKp->SetBranchAddress("Bdtheta",Bdtheta);
	ntKp->SetBranchAddress("Bmu1isTriggered",Bmu1isTriggered);
	ntKp->SetBranchAddress("Bmu2isTriggered",Bmu2isTriggered);
	ntKp->SetBranchAddress("Bmumumass",Bmumumass);
	ntKp->SetBranchAddress("Bmu1eta",Bmu1eta);
	ntKp->SetBranchAddress("Bmu2eta",Bmu2eta);
	ntKp->SetBranchAddress("Bmu1pt",Bmu1pt);
	ntKp->SetBranchAddress("Bmu2pt",Bmu2pt);

	//	ntKp->SetBranchAddress("Bmu1phi",Bmu1phi);
	//	ntKp->SetBranchAddress("Bmu2phi",Bmu2phi);

	ntKp->SetBranchAddress("Bmu1TMOneStationTight",Bmu1TMOneStationTight);
	ntKp->SetBranchAddress("Bmu1InPixelLayer",Bmu1InPixelLayer);
	ntKp->SetBranchAddress("Bmu1InStripLayer",Bmu1InStripLayer);
	ntKp->SetBranchAddress("Bmu2TMOneStationTight",Bmu2TMOneStationTight);
	ntKp->SetBranchAddress("Bmu2InPixelLayer",Bmu2InPixelLayer);
	ntKp->SetBranchAddress("Bmu2InStripLayer",Bmu2InStripLayer);
	ntKp->SetBranchAddress("Bmu1isGlobalMuon",Bmu1isGlobalMuon);
	ntKp->SetBranchAddress("Bmu2isGlobalMuon",Bmu2isGlobalMuon);
	ntKp->SetBranchAddress("Bmu1isTrackerMuon",Bmu1isTrackerMuon);
	ntKp->SetBranchAddress("Bmu2isTrackerMuon",Bmu2isTrackerMuon);


	//	Int_t HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1;


	//	ntHlt->SetBranchAddress("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",&HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1);



	ntKp->SetBranchAddress("Bmu1dxyPV",Bmu1dxyPV);
	ntKp->SetBranchAddress("Bmu2dxyPV",Bmu2dxyPV);
	ntKp->SetBranchAddress("Bmu1dzPV",Bmu1dzPV);
	ntKp->SetBranchAddress("Bmu2dzPV",Bmu2dzPV);
	ntKp->SetBranchAddress("Btrk1highPurity",Btrk1highPurity);
	ntKp->SetBranchAddress("Btrk2highPurity",Btrk2highPurity);
	ntKp->SetBranchAddress("Btktkmass",Btktkmass);
	ntKp->SetBranchAddress("Btrk1PixelHit",Btrk1PixelHit);
	ntKp->SetBranchAddress("Btrk2PixelHit",Btrk2PixelHit);
	ntKp->SetBranchAddress("Btrk1StripHit",Btrk1StripHit);
	ntKp->SetBranchAddress("Btrk2StripHit",Btrk2StripHit);
	ntKp->SetBranchAddress("Btrk1Chi2ndf",Btrk1Chi2ndf);
	ntKp->SetBranchAddress("Btrk2Chi2ndf",Btrk2Chi2ndf);
	ntKp->SetBranchAddress("Bgen",Bgen);
	ntKp->SetBranchAddress("Bgenpt",Bgenpt);
	ntKp->SetBranchAddress("Bgeny",Bgeny);
	ntKp->SetBranchAddress("Btrk1nStripLayer",Btrk1nStripLayer);	
	ntKp->SetBranchAddress("Btrk2nStripLayer",Btrk2nStripLayer);
	ntKp->SetBranchAddress("Btrk1nPixelLayer",Btrk1nPixelLayer);
	ntKp->SetBranchAddress("Btrk2nPixelLayer",Btrk2nPixelLayer);

	if(Rescale == 0){


//	BDT1->SetBranchAddress("BDT_pt_3_5",BDT_pt_3_5);
	BDT2->SetBranchAddress("BDT_pt_5_7",BDT_pt_5_7);
	BDT3->SetBranchAddress("BDT_pt_7_10",BDT_pt_7_10);
	BDT4->SetBranchAddress("BDT_pt_10_15",BDT_pt_10_15);
	BDT5->SetBranchAddress("BDT_pt_15_20",BDT_pt_15_20);
	BDT6->SetBranchAddress("BDT_pt_20_50",BDT_pt_20_50);
//	BDT7->SetBranchAddress("BDT_pt_2_3",BDT_pt_2_3);

//	BDT7->SetBranchAddress("BDT_pt_2_3",BDT_pt_2_3);
//	BDT8->SetBranchAddress("BDT_pt_1_2",BDT_pt_1_2);

	}
	if(Rescale == 1){

//	BDT1->SetBranchAddress("BDT_pt_New_3_5",BDT_pt_3_5);
	BDT2->SetBranchAddress("BDT_pt_New_5_7",BDT_pt_5_7);
	BDT3->SetBranchAddress("BDT_pt_New_7_10",BDT_pt_7_10);
	BDT4->SetBranchAddress("BDT_pt_New_10_15",BDT_pt_10_15);
	BDT5->SetBranchAddress("BDT_pt_New_15_20",BDT_pt_15_20);
	BDT6->SetBranchAddress("BDT_pt_New_20_50",BDT_pt_20_50);
//	BDT7->SetBranchAddress("BDT_pt_New_2_3",BDT_pt_2_3);
//	BDT8->SetBranchAddress("BDT_pt_New_1_2",BDT_pt_1_2);

	}
	Bool_t Bmu1SoftMuID[NCand];
	Bool_t Bmu2SoftMuID[NCand];
	Bool_t Bmu1isAcc[NCand];
	Bool_t Bmu2isAcc[NCand];

	ntKp->SetBranchAddress("Bmu1SoftMuID",Bmu1SoftMuID);
	ntKp->SetBranchAddress("Bmu2SoftMuID",Bmu2SoftMuID);
	ntKp->SetBranchAddress("Bmu1isAcc",Bmu1isAcc);
	ntKp->SetBranchAddress("Bmu2isAcc",Bmu2isAcc);

	Int_t Gsize;
	Float_t Gy[NCand];
	Float_t Gpt[NCand];
	Int_t GisSignal[NCand];
	Int_t GcollisionId[NCand];
	Int_t GpdgId[NCand];
	Float_t Gmu1pt[NCand];
	Float_t Gmu1eta[NCand];
	Float_t Gmu1phi[NCand];
	Float_t Gmu2pt[NCand];
	Float_t Gmu2eta[NCand];
	Float_t Gmu2phi[NCand];
	Float_t Gtk1pt[NCand];
	Float_t Gtk1eta[NCand];
	Float_t Gtk1phi[NCand];

	ntGen->SetBranchAddress("Gsize",&Gsize);
	ntGen->SetBranchAddress("Gy",Gy);
	ntGen->SetBranchAddress("Gpt",Gpt);
	ntGen->SetBranchAddress("GisSignal",GisSignal);
	ntGen->SetBranchAddress("GcollisionId",GcollisionId);
	ntGen->SetBranchAddress("GpdgId",GpdgId);
	ntGen->SetBranchAddress("Gmu1pt",Gmu1pt);
	ntGen->SetBranchAddress("Gmu1eta",Gmu1eta);
	ntGen->SetBranchAddress("Gmu1phi",Gmu1phi);
	ntGen->SetBranchAddress("Gmu2pt",Gmu2pt);
	ntGen->SetBranchAddress("Gmu2eta",Gmu2eta);
	ntGen->SetBranchAddress("Gmu2phi",Gmu2phi);
	ntGen->SetBranchAddress("Gtk1pt",Gtk1pt);
	ntGen->SetBranchAddress("Gtk1eta",Gtk1eta);
	ntGen->SetBranchAddress("Gtk1phi",Gtk1phi);

	Int_t HLT_HIL1DoubleMu0_v1;
	ntHlt->SetBranchAddress("HLT_HIL1DoubleMu0_v1",&HLT_HIL1DoubleMu0_v1);

	double CentWeight;

	//CentWeightTree->SetBranchAddress("CentWeight",&CentWeight);


	double muid1[NCand];
	double muid2[NCand];
	double trk1[NCand];
	double trk2[NCand];
	double trg1[NCand];
	double trg2[NCand];

	/*
	   TnPInfo->SetBranchAddress("muid1",muid1);
	   TnPInfo->SetBranchAddress("muid2",muid2);
	   TnPInfo->SetBranchAddress("trk1",trk1);
	   TnPInfo->SetBranchAddress("trk2",trk2);
	   TnPInfo->SetBranchAddress("trg1",trg1);
	   TnPInfo->SetBranchAddress("trg2",trg2);
	   */




	//Syst Purpose//
	double muid1statup[NCand];
	double trk1statup[NCand];
	double trg1statup[NCand];

	double muid1statdown[NCand];
	double trk1statdown[NCand];
	double trg1statdown[NCand];


	double muid1systup[NCand];
	double trk1systup[NCand];
	double trg1systup[NCand];

	double muid1systdown[NCand];
	double trk1systdown[NCand];
	double trg1systdown[NCand];


	double muid2statup[NCand];
	double trk2statup[NCand];
	double trg2statup[NCand];

	double muid2statdown[NCand];
	double trk2statdown[NCand];
	double trg2statdown[NCand];


	double muid2systup[NCand];
	double trk2systup[NCand];
	double trg2systup[NCand];

	double muid2systdown[NCand];
	double trk2systdown[NCand];
	double trg2systdown[NCand];



	double muid1syst;
	double muid1stat;
	double muid2syst;
	double muid2stat;

	double trk1syst;
	double trk1stat;
	double trk2syst;
	double trk2stat;

	double trg1syst;
	double trg1stat;
	double trg2syst;
	double trg2stat;


	double tnptotal1syst;
	double tnptotal1stat;


	double tnptotal2syst;
	double tnptotal2stat;


	double tnptotal1err;
	double tnptotal2err;

	double tnptotalerr;


	/*



	   TnPInfo->SetBranchAddress("muid1statup",muid1statup);
	   TnPInfo->SetBranchAddress("trk1statup",trk1statup);
	   TnPInfo->SetBranchAddress("trg1statup",trg1statup);
	   TnPInfo->SetBranchAddress("muid1statdown",muid1statdown);
	   TnPInfo->SetBranchAddress("trk1statdown",trk1statdown);
	   TnPInfo->SetBranchAddress("trg1statdown",trg1statdown);
	   TnPInfo->SetBranchAddress("muid1systup",muid1systup);
	   TnPInfo->SetBranchAddress("trk1systup",trk1systup);
	   TnPInfo->SetBranchAddress("trg1systup",trg1systup);
	   TnPInfo->SetBranchAddress("muid1systdown",muid1systdown);
	   TnPInfo->SetBranchAddress("trk1systdown",trk1systdown);
	   TnPInfo->SetBranchAddress("trg1systdown",trg1systdown);



	   TnPInfo->SetBranchAddress("muid2statup",muid2statup);
	   TnPInfo->SetBranchAddress("trk2statup",trk2statup);
	   TnPInfo->SetBranchAddress("trg2statup",trg2statup);
	   TnPInfo->SetBranchAddress("muid2statdown",muid2statdown);
	   TnPInfo->SetBranchAddress("trk2statdown",trk2statdown);
	   TnPInfo->SetBranchAddress("trg2statdown",trg2statdown);
	   TnPInfo->SetBranchAddress("muid2systup",muid2systup);
	   TnPInfo->SetBranchAddress("trk2systup",trk2systup);
	   TnPInfo->SetBranchAddress("trg2systup",trg2systup);
	   TnPInfo->SetBranchAddress("muid2systdown",muid2systdown);
	   TnPInfo->SetBranchAddress("trk2systdown",trk2systdown);
	   TnPInfo->SetBranchAddress("trg2systdown",trg2systdown);
	   */


	float Factor = 4.0;


	int NEvents = ntKp->GetEntries();

	const int yBinN = 5;
  std::vector<double> yBins ({0.0, 0.5, 1.0, 1.5, 2.0, 2.4});
	double yBinning[yBinN+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.4};

  // create a vector of pT binning with specified regional widths
  auto createBins = [] (std::vector<double> edges,
                        std::vector<double> binWidth) {
    std::vector<double> bins = {edges[0]};
    while (bins.back() < edges.back()) {
      auto iRegion = std::upper_bound(edges.begin(), edges.end(), bins.back())
        - edges.begin() - 1;
      bins.push_back(bins.back() + binWidth.at(iRegion));
    }
    return bins;
  };

  auto bptBinVec = createBins({0, 10, 40, 50, 60}, {1/8., 1/4., 1/2., 1});
  auto BptBinning = bptBinVec.data();
  const int BptBin = bptBinVec.size() - 1;

	double PVzWeight;


	double EventWeight;
	double TnPWeight;
	double muidWeight;
	double trkWeight;
	double TotalWeight;
	double muidtrkWeight;


	double TotalWeightSystUp;
	double TotalWeightSystDown;

	double TotalWeightMuidUp;
	double TotalWeightMuidDown;
	double TotalWeightTrkUp;
	double TotalWeightTrkDown;
	double TotalWeightTrgUp;
	double TotalWeightTrgDown;


	double muid1total;
	double muid2total;

	double trk1total;
	double trk2total;

	double trg1total;
	double trg2total;


	double muidtotalerr;
	double trktotalerr;
	double trgtotalerr;


	TH2D * NoWeightHis = new TH2D("NoWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * EvtWeightHis = new TH2D("EvtWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * muidWeightHis = new TH2D("muidWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * trkWeightHis = new TH2D("trkWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * muidtrkWeightHis = new TH2D("muidtrkWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHis = new TH2D("TnPWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisSystUp = new TH2D("TnPWeightHisSystUp","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisSystDown = new TH2D("TnPWeightHisSystDown","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * BDTWeightHisSyst = new TH2D("BDTWeightHisSyst","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * BptWeightHisSyst = new TH2D("BptWeightHisSyst","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TrkLooseHis = (TH2D*) TnPWeightHis->Clone("TrkLooseHis");
	TH2D * TrkTightHis = (TH2D*) TnPWeightHis->Clone("TrkTightHis");
	TH2D * TnPWeightHisMuidUp = new TH2D("TnPWeightHisMuidUp","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisMuidDown = new TH2D("TnPWeightHisMuidDown","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisTrkUp = new TH2D("TnPWeightHisTrkUp","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisTrkDown = new TH2D("TnPWeightHisTrkDown","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisTrgUp = new TH2D("TnPWeightHisTrgUp","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisTrgDown = new TH2D("TnPWeightHisTrgDown","",BptBin,BptBinning,yBinN,yBinning);

	//Gen//
	TH2D * NoWeightGenHis = new TH2D("NoWeightGenHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * EvtWeightGenHis = new TH2D("EvtWeightGenHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * BptWeightGenHis = new TH2D("BptWeightGenHis","",BptBin,BptBinning,yBinN,yBinning);	
	TH2D * EvtWeightGenAccHis = new TH2D("EvtWeightGenAccHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * NoWeightGenAccHis = new TH2D("NoWeightGenAccHis","",BptBin,BptBinning,yBinN,yBinning);


	TH1D * Bmu1ptHis = new TH1D("Bmu1ptHis","",200,0,50);
	Bmu1ptHis->GetXaxis()->SetTitle("Bmu1pt (GeV/c)");	
	Bmu1ptHis->GetYaxis()->SetTitle("Counts");
	Bmu1ptHis->GetXaxis()->CenterTitle();	
	Bmu1ptHis->GetYaxis()->CenterTitle();
	Bmu1ptHis->GetXaxis()->SetTitleOffset(1.2);	
	Bmu1ptHis->GetYaxis()->SetTitleOffset(1.5);



	const int NPtBins = 10;
	double PtBin[NPtBins + 1] = {0,1,2,3,5,7,10,15,20,50,100};
	

//	const int NPtBins1D = 11;
	//double  PtBin1D[NPtBins1D + 1] = {0,1,2,3,5,7,10,15,20,30,50,100};
	/*

	const int NPtBins1D = 10;
	double  PtBin1D[NPtBins1D + 1] = {1,2,3,5,7,10,15,20,30,50,100};
	*/
	
//	const int NPtBins1D = 4;
//	double  PtBin1D[NPtBins1D + 1] = {7,10,15,20,50};

	
	const int NPtBins1D = 7;
	double  PtBin1D[NPtBins1D + 1] = {5,7,10,15,20,30,50,60};

	const int NMultiBin = 7;
	double  MultiBin1D[NMultiBin + 1] = {0,20,30,40,50,60,70,100};

	const int YBin = 8;
	double  YBin1D[YBin + 1] = {-2.4,-1.5,-1.0,-0.5,0.0 ,0.5, 1.0, 1.5, 2.4};

//	TH1D * Eff1DRECOHis = new TH1D("Eff1DRECOHis","",NPtBins,PtBin);
	TH1D * Eff1DRECOHis = new TH1D("Eff1DRECOHis","",NPtBins1D,PtBin1D);

	Eff1DRECOHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DRECOHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHis->GetXaxis()->CenterTitle();	
	Eff1DRECOHis->GetYaxis()->CenterTitle();
	Eff1DRECOHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHis->GetYaxis()->SetTitleOffset(1.5);


	//TnP Varied

	TH1D * Eff1DRECOHisTnPUp = new TH1D("Eff1DRECOHisTnPUp","",NPtBins1D,PtBin1D);

	Eff1DRECOHisTnPUp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DRECOHisTnPUp->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisTnPUp->GetXaxis()->CenterTitle();	
	Eff1DRECOHisTnPUp->GetYaxis()->CenterTitle();
	Eff1DRECOHisTnPUp->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisTnPUp->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff1DRECOHisTnPDown = new TH1D("Eff1DRECOHisTnPDown","",NPtBins1D,PtBin1D);

	Eff1DRECOHisTnPDown->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DRECOHisTnPDown->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisTnPDown->GetXaxis()->CenterTitle();	
	Eff1DRECOHisTnPDown->GetYaxis()->CenterTitle();
	Eff1DRECOHisTnPDown->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisTnPDown->GetYaxis()->SetTitleOffset(1.5);




	//BDT Weighted
	TH1D * Eff1DRECOHisBDT = new TH1D("Eff1DRECOHisBDT","",NPtBins1D,PtBin1D);

	Eff1DRECOHisBDT->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DRECOHisBDT->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisBDT->GetXaxis()->CenterTitle();	
	Eff1DRECOHisBDT->GetYaxis()->CenterTitle();
	Eff1DRECOHisBDT->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisBDT->GetYaxis()->SetTitleOffset(1.5);
	


	TFile * finBDTWeight = new TFile("BDTWeights/BPw.root");

	TH1D * weights_BDT_pt_5_7 = (TH1D * ) finBDTWeight->Get("weights_BDT_pt_5_7");
	TH1D * weights_BDT_pt_7_10 = (TH1D * ) finBDTWeight->Get("weights_BDT_pt_7_10");
	TH1D * weights_BDT_pt_10_15 = (TH1D * ) finBDTWeight->Get("weights_BDT_pt_10_15");
	TH1D * weights_BDT_pt_15_20 = (TH1D * ) finBDTWeight->Get("weights_BDT_pt_15_20");
	TH1D * weights_BDT_pt_20_50 = (TH1D * ) finBDTWeight->Get("weights_BDT_pt_20_50");

	
	TH1D * Eff1DRECOHisBpt = new TH1D("Eff1DRECOHisBpt","",NPtBins1D,PtBin1D);

	Eff1DRECOHisBpt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DRECOHisBpt->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisBpt->GetXaxis()->CenterTitle();	
	Eff1DRECOHisBpt->GetYaxis()->CenterTitle();
	Eff1DRECOHisBpt->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisBpt->GetYaxis()->SetTitleOffset(1.5);




	TH1D * Eff1DRECOHisTnPUpY = new TH1D("Eff1DRECOHisTnPUpY","",YBin,YBin1D);

	Eff1DRECOHisTnPUpY->GetXaxis()->SetTitle("rapidity");
	Eff1DRECOHisTnPUpY->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisTnPUpY->GetXaxis()->CenterTitle();	
	Eff1DRECOHisTnPUpY->GetYaxis()->CenterTitle();
	Eff1DRECOHisTnPUpY->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisTnPUpY->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff1DRECOHisTnPDownY = new TH1D("Eff1DRECOHisTnPDownY","",YBin,YBin1D);

	Eff1DRECOHisTnPDownY->GetXaxis()->SetTitle("rapidity");
	Eff1DRECOHisTnPDownY->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisTnPDownY->GetXaxis()->CenterTitle();	
	Eff1DRECOHisTnPDownY->GetYaxis()->CenterTitle();
	Eff1DRECOHisTnPDownY->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisTnPDownY->GetYaxis()->SetTitleOffset(1.5);	

	//Mult Stuffs

	TH1D * Eff1DRECOHisTnPUpMult = new TH1D("Eff1DRECOHisTnPUp","",NMultiBin,MultiBin1D);

	Eff1DRECOHisTnPUpMult->GetXaxis()->SetTitle("Multiplicity");
	Eff1DRECOHisTnPUpMult->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisTnPUpMult->GetXaxis()->CenterTitle();	
	Eff1DRECOHisTnPUpMult->GetYaxis()->CenterTitle();
	Eff1DRECOHisTnPUpMult->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisTnPUpMult->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff1DRECOHisTnPDownMult = new TH1D("Eff1DRECOHisTnPDownMult","",NMultiBin,MultiBin1D);

	Eff1DRECOHisTnPDownMult->GetXaxis()->SetTitle("Multiplicity");
	Eff1DRECOHisTnPDownMult->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisTnPDownMult->GetXaxis()->CenterTitle();	
	Eff1DRECOHisTnPDownMult->GetYaxis()->CenterTitle();
	Eff1DRECOHisTnPDownMult->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisTnPDownMult->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff1DRECOHisBDTMult = new TH1D("Eff1DRECOHisBDTMult","",NMultiBin,MultiBin1D);

	Eff1DRECOHisBDTMult->GetXaxis()->SetTitle("Multiplicity");
	Eff1DRECOHisBDTMult->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisBDTMult->GetXaxis()->CenterTitle();	
	Eff1DRECOHisBDTMult->GetYaxis()->CenterTitle();
	Eff1DRECOHisBDTMult->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisBDTMult->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff1DRECOHisBptMult = new TH1D("Eff1DRECOHisBptMult","",NMultiBin,MultiBin1D);

	Eff1DRECOHisBptMult->GetXaxis()->SetTitle("Multiplicity");
	Eff1DRECOHisBptMult->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHisBptMult->GetXaxis()->CenterTitle();	
	Eff1DRECOHisBptMult->GetYaxis()->CenterTitle();
	Eff1DRECOHisBptMult->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHisBptMult->GetYaxis()->SetTitleOffset(1.5);



	float BptWeight;

//	TF1 * BptWFunc = new TF1("BptWFunc","8.646057/(x*x) +0.425834*TMath::Log(x) -0.148249",0,100);

//	TF1 * BptWFunc = new TF1("BptWFunc","20801.474609/(x*x*x*x*x) -92.810738/(x*x) + 1.475384 ",0,100);

	// TF1 * BptWFunc = new TF1("BptWFunc"," 968.466885/x**(4.494808) + 0.800573 + 0.003677 * x",0,100);
	// TF1 * BptWFunc = new TF1("BptWFunc","1.000000/(x*x) +0.435893*TMath::Log(x) - 0.116910",0,100);
	// TF1 * BptWFunc = new TF1("BptWFunc","1.070585/x**(8.245110) + 0.833796 + 0.016723 * x",0,100);
	TF1 * BptWFunc = new TF1("BptWFunc","10.577117/x**(1.906323) + 0.654119 + 0.012688 * x",0,100);

  TFile fBptWeight("../../NewBptStudies/ResultFile/BptWeight_BP.root");
  std::map<int, TF1*> BptWtF;
  if (reweightPtOnY) {
    for (auto iy = 0; iy < yBinN; ++iy) {
      BptWtF[iy] = (TF1*) fBptWeight.Get(TString::Format("BptWeight_y%d", iy));
    	}
  }

	int BDTWeightBin;
	float BDTWeight;

	TH1D * Eff1DRECOMultHis = new TH1D("Eff1DRECOMultHis","",NMultiBin,MultiBin1D);

	Eff1DRECOMultHis->GetXaxis()->SetTitle("Multiplicity");
	Eff1DRECOMultHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOMultHis->GetXaxis()->CenterTitle();	
	Eff1DRECOMultHis->GetYaxis()->CenterTitle();
	Eff1DRECOMultHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOMultHis->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff1DRECOYHis = new TH1D("Eff1DRECOYHis","",YBin,YBin1D);

	Eff1DRECOYHis->GetXaxis()->SetTitle("Rapidity");
	Eff1DRECOYHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOYHis->GetXaxis()->CenterTitle();	
	Eff1DRECOYHis->GetYaxis()->CenterTitle();
	Eff1DRECOYHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOYHis->GetYaxis()->SetTitleOffset(1.5);

//	TH1D * Eff1DGENHis = new TH1D("Eff1DGENHis","",NPtBins,PtBin);
	TH1D * Eff1DGENHis = new TH1D("Eff1DGENHis","",NPtBins1D,PtBin1D);
	
	Eff1DGENHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DGENHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENHis->GetXaxis()->CenterTitle();	
	Eff1DGENHis->GetYaxis()->CenterTitle();
	Eff1DGENHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENHis->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff1DGENHisGpt = new TH1D("Eff1DGENHisGpt","",NPtBins1D,PtBin1D);
	
	Eff1DGENHisGpt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DGENHisGpt->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENHisGpt->GetXaxis()->CenterTitle();	
	Eff1DGENHisGpt->GetYaxis()->CenterTitle();
	Eff1DGENHisGpt->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENHisGpt->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff1DGENAccHis = new TH1D("Eff1DGENAccHis","",NPtBins1D,PtBin1D);
	
	Eff1DGENAccHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DGENAccHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENAccHis->GetXaxis()->CenterTitle();	
	Eff1DGENAccHis->GetYaxis()->CenterTitle();
	Eff1DGENAccHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENAccHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff1DGENAccYHis = new TH1D("Eff1DGENAccYHis","",YBin,YBin1D);
	
	Eff1DGENAccYHis->GetXaxis()->SetTitle("rapidity");
	Eff1DGENAccYHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENAccYHis->GetXaxis()->CenterTitle();	
	Eff1DGENAccYHis->GetYaxis()->CenterTitle();
	Eff1DGENAccYHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENAccYHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff1DGENMultHis = new TH1D("Eff1DGENMultHis","",NMultiBin,MultiBin1D);

	Eff1DGENMultHis->GetXaxis()->SetTitle("Multiplicity");
	Eff1DGENMultHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENMultHis->GetXaxis()->CenterTitle();	
	Eff1DGENMultHis->GetYaxis()->CenterTitle();
	Eff1DGENMultHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENMultHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff1DGENYHis = new TH1D("Eff1DGENYHis","",YBin,YBin1D);

	Eff1DGENYHis->GetXaxis()->SetTitle("Rapidity");
	Eff1DGENYHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENYHis->GetXaxis()->CenterTitle();	
	Eff1DGENYHis->GetYaxis()->CenterTitle();
	Eff1DGENYHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENYHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff1DGENMultHisGpt = new TH1D("Eff1DGENMultHisGpt","",NMultiBin,MultiBin1D);

	Eff1DGENMultHisGpt->GetXaxis()->SetTitle("Multiplicity");
	Eff1DGENMultHisGpt->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENMultHisGpt->GetXaxis()->CenterTitle();	
	Eff1DGENMultHisGpt->GetYaxis()->CenterTitle();
	Eff1DGENMultHisGpt->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENMultHisGpt->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff1DGENAccMultHis = new TH1D("Eff1DGENAccMultHis","",NMultiBin,MultiBin1D);

	Eff1DGENAccMultHis->GetXaxis()->SetTitle("Multiplicity");
	Eff1DGENAccMultHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENAccMultHis->GetXaxis()->CenterTitle();	
	Eff1DGENAccMultHis->GetYaxis()->CenterTitle();
	Eff1DGENAccMultHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENAccMultHis->GetYaxis()->SetTitleOffset(1.5);



	TH1D * Bmu2ptHis = new TH1D("Bmu2ptHis","",200,0,50);
	Bmu2ptHis->GetXaxis()->SetTitle("Bmu2pt (GeV/c)");	
	Bmu2ptHis->GetYaxis()->SetTitle("Counts");
	Bmu2ptHis->GetXaxis()->CenterTitle();	
	Bmu2ptHis->GetYaxis()->CenterTitle();
	Bmu2ptHis->GetXaxis()->SetTitleOffset(1.2);	
	Bmu2ptHis->GetYaxis()->SetTitleOffset(1.5);

	const int NPtbin = 200;
	double PtMin = 0;
	double PtMax = 100;
	double PtStep = (PtMax - PtMin)/NPtbin;

	const int NEtabin = 4;
	double Etabinning[NEtabin + 1] = {0,1.2,1.8,2.1,2.4};



	TH2D * Bmu1pteta = new TH2D("Bmu1pteta","",NPtbin,PtMin,PtMax,NEtabin,Etabinning);
	Bmu1pteta->GetXaxis()->SetTitle("Bmu1pt");	
	Bmu1pteta->GetYaxis()->SetTitle("Bmu1eta");
	Bmu1pteta->GetXaxis()->CenterTitle();	
	Bmu1pteta->GetYaxis()->CenterTitle();
	Bmu1pteta->GetXaxis()->SetTitleOffset(1.2);	
	Bmu1pteta->GetYaxis()->SetTitleOffset(1.5);



	TH2D * Bmu2pteta = new TH2D("Bmu2pteta","",NPtbin,PtMin,PtMax,NEtabin,Etabinning);
	Bmu2pteta->GetXaxis()->SetTitle("Bmu2pt");	
	Bmu2pteta->GetYaxis()->SetTitle("Bmu2eta");
	Bmu2pteta->GetXaxis()->CenterTitle();	
	Bmu2pteta->GetYaxis()->CenterTitle();
	Bmu2pteta->GetXaxis()->SetTitleOffset(1.2);	
	Bmu2pteta->GetYaxis()->SetTitleOffset(1.5);



	TH1D * Bmu1etaHis = new TH1D("Bmu1etaHis","",100,0,3);
	Bmu1etaHis->GetXaxis()->SetTitle("Bmu1eta");	
	Bmu1etaHis->GetYaxis()->SetTitle("Counts");
	Bmu1etaHis->GetXaxis()->CenterTitle();	
	Bmu1etaHis->GetYaxis()->CenterTitle();
	Bmu1etaHis->GetXaxis()->SetTitleOffset(1.2);	
	Bmu1etaHis->GetYaxis()->SetTitleOffset(1.5);




	TH1D * Bmu2etaHis = new TH1D("Bmu2etaHis","",100,0,3);
	Bmu2etaHis->GetXaxis()->SetTitle("Bmu2eta");	
	Bmu2etaHis->GetYaxis()->SetTitle("Counts");
	Bmu2etaHis->GetXaxis()->CenterTitle();	
	Bmu2etaHis->GetYaxis()->CenterTitle();
	Bmu2etaHis->GetXaxis()->SetTitleOffset(1.2);	
	Bmu2etaHis->GetYaxis()->SetTitleOffset(1.5);




	TH1D * Bmu1TrgSF = new TH1D("Bmu1TrgSF","",100,0.7,1.3);
	Bmu1TrgSF->GetXaxis()->SetTitle("Bmu1TrgSF");	
	Bmu1TrgSF->GetYaxis()->SetTitle("Counts");
	Bmu1TrgSF->GetXaxis()->CenterTitle();	
	Bmu1TrgSF->GetYaxis()->CenterTitle();
	Bmu1TrgSF->GetXaxis()->SetTitleOffset(1.2);	
	Bmu1TrgSF->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Bmu2TrgSF = new TH1D("Bmu2TrgSF","",100,0.7,1.3);
	Bmu2TrgSF->GetXaxis()->SetTitle("Bmu2TrgSF");	
	Bmu2TrgSF->GetYaxis()->SetTitle("Counts");
	Bmu2TrgSF->GetXaxis()->CenterTitle();	
	Bmu2TrgSF->GetYaxis()->CenterTitle();
	Bmu2TrgSF->GetXaxis()->SetTitleOffset(1.2);	
	Bmu2TrgSF->GetYaxis()->SetTitleOffset(1.5);

	
	Float_t TnPNominal[NCand];
	Float_t TnPMu1Nominal[NCand];
	Float_t TnPMu2Nominal[NCand];
	Float_t TnPError[NCand];

	TnPInfo->SetBranchAddress("TnPNominal",TnPNominal);
	TnPInfo->SetBranchAddress("TnPMu1Nominal",TnPMu1Nominal);
	TnPInfo->SetBranchAddress("TnPMu2Nominal",TnPMu2Nominal);
	TnPInfo->SetBranchAddress("TnPError",TnPError);

	
	TH1D * TnPSFHis = new TH1D("TnPSFHis","",100,0.8,1.2);
	TnPSFHis->GetXaxis()->SetTitle("SF^{#mu_{1}} #times SF^{#mu_{2}}");
	TnPSFHis->GetYaxis()->SetTitle("Counts");
	TnPSFHis->GetXaxis()->CenterTitle();	
	TnPSFHis->GetYaxis()->CenterTitle();
	TnPSFHis->GetXaxis()->SetTitleOffset(1.2);	
	TnPSFHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * TnPSFHisMu1 = new TH1D("TnPSFHisMu1","",100,0.8,1.2);
	TnPSFHisMu1->GetXaxis()->SetTitle("SF^{#mu_{1}}");
	TnPSFHisMu1->GetYaxis()->SetTitle("Counts");
	TnPSFHisMu1->GetXaxis()->CenterTitle();	
	TnPSFHisMu1->GetYaxis()->CenterTitle();
	TnPSFHisMu1->GetXaxis()->SetTitleOffset(1.2);	
	TnPSFHisMu1->GetYaxis()->SetTitleOffset(1.5);

	TH1D * TnPSFHisMu2 = new TH1D("TnPSFHisMu2","",100,0.8,1.2);
	TnPSFHisMu2->GetXaxis()->SetTitle("SF^{#mu_{2}}");
	TnPSFHisMu2->GetYaxis()->SetTitle("Counts");
	TnPSFHisMu2->GetXaxis()->CenterTitle();	
	TnPSFHisMu2->GetYaxis()->CenterTitle();
	TnPSFHisMu2->GetXaxis()->SetTitleOffset(1.2);	
	TnPSFHisMu2->GetYaxis()->SetTitleOffset(1.5);



	//	TFile * fout = new TFile("WeightCheck.root","RECREATE");

	TString outfileName;

	if(Rescale == 0){
		if(DoTnP == 0) outfileName = "NewEff2DMaps/EffFineNoTnP.root";
		if(DoTnP == 1) outfileName = "NewEff2DMaps/EffFineBDT.root";
	}
	
	if(Rescale == 1){

		if(DoTnP == 1) outfileName = "NewEff2DMaps/EffFineBDTNew.root";
	}

	TFile * fout = new TFile(outfileName.Data(),"RECREATE");

	fout->cd();


	//NEvents = 10;

	for(int i = 0; i < NEvents; i++){


		if(i%10000==0) cout << "Now Working on = " << i  <<  endl;

		ntKp->GetEntry(i);
		ntSkim->GetEntry(i);
		ntHi->GetEntry(i);
		ntHlt->GetEntry(i);
		//	BDT->GetEntry(i);
		//	CentWeightTree->GetEntry(i);
		//	TnPInfo->GetEntry(i);
		ntGen->GetEntry(i);

		//BDT//
//		BDT1->GetEntry(i);
		BDT2->GetEntry(i);
		BDT3->GetEntry(i);
		BDT4->GetEntry(i);
		BDT5->GetEntry(i);
		BDT6->GetEntry(i);
//		BDT7->GetEntry(i);
//		BDT8->GetEntry(i);
	
		root->GetEntry(i);

	
		TnPInfo->GetEntry(i);


		for(int j = 0; j < Bsize; j++){
      bool passBDT = ((Bpt[j] > 5 && Bpt[j] < 7 && BDT_pt_5_7[j] > 0.08)
                      || (Bpt[j] > 7 && Bpt[j] < 10 && BDT_pt_7_10[j] > 0.07)
                      || (Bpt[j] > 10 && Bpt[j] < 15 && BDT_pt_10_15[j] > 0)
                      || (Bpt[j] > 15 && Bpt[j] < 20 && BDT_pt_15_20[j] > 0.02)
                      || (Bpt[j] > 20 && Bpt[j] < 50 && BDT_pt_20_50[j] > 0.04)
                      || (Bpt[j] > 50 && Bpt[j] < 60));

      bool preselection = (pPAprimaryVertexFilter == 1
                           && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1)
        && (Bmu1isTriggered[j] == 1 && Bmu2isTriggered[j] == 1)
        && (Btrk1Pt[j] > 0.5 && Bchi2cl[j] > 0.05 && BsvpvDistance[j]/BsvpvDisErr[j] > 2.0
            && Bpt[j] > 2 && abs(Btrk1Eta[j]-0.0) < 2.4  &&
            (TMath::Abs(By[j])<2.4 && TMath::Abs(Bmumumass[j]-3.096916)<0.15
             &&((abs(Bmu1eta[j])<1.2&&Bmu1pt[j]>3.5)
             ||(abs(Bmu1eta[j])>1.2&&abs(Bmu1eta[j])<2.1 &&Bmu1pt[j]>(5.47-1.89*abs(Bmu1eta[j])))
             ||(abs(Bmu1eta[j])>2.1&&abs(Bmu1eta[j])<2.4&&Bmu1pt[j]>1.5))
             &&((abs(Bmu2eta[j])<1.2&&Bmu2pt[j]>3.5)
             ||(abs(Bmu2eta[j])>1.2&&abs(Bmu2eta[j])<2.1&&Bmu2pt[j]>(5.47-1.89*abs(Bmu2eta[j])))
             ||(abs(Bmu2eta[j])>2.1&&abs(Bmu2eta[j])<2.4&&Bmu2pt[j]>1.5))
             &&Bmu1TMOneStationTight[j]&&Bmu2TMOneStationTight[j]
             &&Bmu1InPixelLayer[j]>0&&(Bmu1InPixelLayer[j]+Bmu1InStripLayer[j])>5
             &&Bmu2InPixelLayer[j]>0&& (Bmu2InPixelLayer[j]+Bmu2InStripLayer[j])>5
             &&Bmu1dxyPV[j]<0.3&&Bmu2dxyPV[j]<0.3
             &&Bmu1dzPV[j]<20&&Bmu2dzPV[j]<20
             &&Bmu1isTrackerMuon[j]&&Bmu2isTrackerMuon[j]
             &&Bmu1isGlobalMuon[j]&&Bmu2isGlobalMuon[j]
             &&Btrk1highPurity[j]
             &&abs(Btrk1Eta[j])<2.4&&Btrk1Pt[j]>0.5)
            && (Btrk1PixelHit[j] + Btrk1StripHit[j] > 10) && (abs(PVz)<15));

      // tracking selection variation
      auto passTracking = [=] (Tracking cut) {
        return (Btrk1PtErr[j]/Btrk1Pt[j] < ptErr[cut] &&
                Btrk1Chi2ndf[j]/(Btrk1nStripLayer[j]+Btrk1nPixelLayer[j])
                < chi2Nlayer[cut]);
      };

      if ((Bgen[j] == 23333) && preselection
          && passTracking(Tracking::loose) && passBDT){

				
				//EventWeight = pthat * weight;
				EventWeight = weight;
				//Turn off TnP weight for now
				
				if(DoTnP == 0){

				muid1[j] = 1;
				muid2[j] = 1;
				trk1[j] = 1;
				trk2[j] = 1;
				trg1[j] = 1;
				trg2[j] = 1;
			
				TnPWeight = muid1[j] * trk1[j] * trg1[j] * muid2[j] * trk2[j] * trg2[j];

				}

				if(DoTnP == 1) TnPWeight = TnPNominal[j];

				muidWeight = EventWeight * muid1[j] * muid2[j];
				trkWeight = EventWeight * trk1[j] * trk2[j];

				muidtrkWeight = EventWeight * muid1[j] * muid2[j] * trk1[j] * trk2[j];
				//TotalWeight = EventWeight * TnPWeight;
				TotalWeight=1;


        if (passTracking(Tracking::tight)) {
          TrkTightHis->Fill(Bpt[j],abs(By[j]),TotalWeight);
        }
        if (passTracking(Tracking::standard)) {
          TnPWeightHis->Fill(Bpt[j],abs(By[j]),TotalWeight);
          TnPSFHis->Fill(TnPWeight);
          TnPSFHisMu1->Fill(TnPMu1Nominal[j]);
          TnPSFHisMu2->Fill(TnPMu2Nominal[j]);
          NoWeightHis->Fill(Bpt[j],abs(By[j]),1);
          EvtWeightHis->Fill(Bpt[j],abs(By[j]),EventWeight);
          muidWeightHis->Fill(Bpt[j],abs(By[j]),muidWeight);
          trkWeightHis->Fill(Bpt[j],abs(By[j]),trkWeight);
          muidtrkWeightHis->Fill(Bpt[j],abs(By[j]),muidtrkWeight);
          Eff1DRECOHis->Fill(Bpt[j],TotalWeight);
          Eff1DRECOMultHis->Fill(nMult,TotalWeight);
		  Eff1DRECOYHis->Fill(By[j],TotalWeight);
        }
        TrkLooseHis->Fill(Bpt[j],abs(By[j]),TotalWeight);

        // The following systematics are only related to standard track sel
        if (!passTracking(Tracking::standard)) {
          continue;
        }


				//Now Everything is About the Error//
				if(muid1systup[j] >= muid1systdown[j]) muid1syst = muid1systup[j];
				if(muid1systdown[j] > muid1systup[j]) muid1syst = muid1systdown[j];
				if(muid2systup[j] >= muid2systdown[j]) muid2syst = muid2systup[j];
				if(muid2systdown[j] > muid2systup[j]) muid2syst = muid2systdown[j];

				if(muid1statup[j] >= muid1statdown[j]) muid1stat = muid1statup[j];
				if(muid1statdown[j] > muid1statup[j]) muid1stat = muid1statdown[j];
				if(muid2statup[j] >= muid2statdown[j]) muid2stat = muid2statup[j];
				if(muid2statdown[j] > muid2statup[j]) muid2stat = muid2statdown[j];

				//Trk
				if(trk1systup[j] >= trk1systdown[j]) trk1syst = trk1systup[j];
				if(trk1systdown[j] > trk1systup[j]) trk1syst = trk1systdown[j];
				if(trk2systup[j] >= trk2systdown[j]) trk2syst = trk2systup[j];
				if(trk2systdown[j] > trk2systup[j]) trk2syst = trk2systdown[j];

				if(trk1statup[j] >= trk1statdown[j]) trk1stat = trk1statup[j];
				if(trk1statdown[j] > trk1statup[j]) trk1stat = trk1statdown[j];
				if(trk2statup[j] >= trk2statdown[j]) trk2stat = trk2statup[j];
				if(trk2statdown[j] > trk2statup[j]) trk2stat = trk2statdown[j];



				//Trg
				if(trg1systup[j] >= trg1systdown[j]) trg1syst = trg1systup[j];
				if(trg1systdown[j] > trg1systup[j]) trg1syst = trg1systdown[j];
				if(trg2systup[j] >= trg2systdown[j]) trg2syst = trg2systup[j];
				if(trg2systdown[j] > trg2systup[j]) trg2syst = trg2systdown[j];

				if(trg1statup[j] >= trg1statdown[j]) trg1stat = trg1statup[j];
				if(trg1statdown[j] > trg1statup[j]) trg1stat = trg1statdown[j];
				if(trg2statup[j] >= trg2statdown[j]) trg2stat = trg2statup[j];
				if(trg2statdown[j] > trg2statup[j]) trg2stat = trg2statdown[j];


				tnptotal1syst = sqrt(muid1syst/muid1[j] * muid1syst/muid1[j] + trk1syst/trk1[j] * trk1syst/trk1[j] + trg1syst/trg1[j] * trg1syst/trg1[j]);
				tnptotal1stat = sqrt(muid1stat/muid1[j] * muid1stat/muid1[j] + trk1stat/trk1[j] * trk1stat/trk1[j] + trg1stat/trg1[j] * trg1stat/trg1[j]);


				tnptotal2syst = sqrt(muid2syst/muid2[j] * muid2syst/muid2[j] + trk2syst/trk2[j] * trk2syst/trk2[j] + trg2syst/trg2[j] * trg2syst/trg2[j]);
				tnptotal2stat = sqrt(muid2stat/muid2[j] * muid2stat/muid2[j] + trk2stat/trk2[j] * trk2stat/trk2[j] + trg2stat/trg2[j] * trg2stat/trg2[j]);


				tnptotal1err = sqrt(tnptotal1stat * tnptotal1stat + tnptotal1syst * tnptotal1syst);
				tnptotal2err = sqrt(tnptotal2stat * tnptotal2stat + tnptotal2syst * tnptotal2syst);	

				tnptotalerr = sqrt(tnptotal1err * tnptotal1err + tnptotal2err * tnptotal2err);

				//	cout << "tnptotalerr = " << tnptotalerr << endl;
			 
				if(DoTnP == 1) tnptotalerr = TnPError[j];

				TotalWeightSystUp = TotalWeight * ( 1 + tnptotalerr);
				TotalWeightSystDown = TotalWeight * ( 1 - tnptotalerr);

				TnPWeightHisSystUp->Fill(Bpt[j],abs(By[j]),TotalWeightSystUp);
				TnPWeightHisSystDown->Fill(Bpt[j],abs(By[j]),TotalWeightSystDown);

				Eff1DRECOHisTnPUp->Fill(Bpt[j],TotalWeightSystUp);
				Eff1DRECOHisTnPDown->Fill(Bpt[j],TotalWeightSystDown);

				Eff1DRECOHisTnPUpMult->Fill(nMult,TotalWeightSystUp);
				Eff1DRECOHisTnPDownMult->Fill(nMult,TotalWeightSystDown);
				
				Eff1DRECOHisTnPUpY->Fill(By[j],TotalWeightSystUp);
				Eff1DRECOHisTnPDownY->Fill(By[j],TotalWeightSystDown);

				muid1total =  sqrt(muid1syst/muid1[j] * muid1syst/muid1[j] + muid1stat/muid1[j] * muid1stat/muid1[j]);
				muid2total =  sqrt(muid2syst/muid2[j] * muid2syst/muid2[j] + muid2stat/muid2[j] * muid2stat/muid2[j]);

				trk1total =  sqrt(trk1syst/trk1[j] * trk1syst/trk1[j] + trk1stat/trk1[j] * trk1stat/trk1[j]);
				trk2total =	 sqrt(trk2syst/trk2[j] * trk2syst/trk2[j] + trk2stat/trk2[j] * trk2stat/trk2[j]);

				trg1total =  sqrt(trg1syst/trg1[j] * trg1syst/trg1[j] + trg1stat/trg1[j] * trg1stat/trg1[j]);
				trg2total =	sqrt(trg2syst/trg2[j] * trg2syst/trg2[j] + trg2stat/trg2[j] * trg2stat/trg2[j]);


				muidtotalerr = sqrt(muid1total * muid1total + muid2total * muid2total);
				trktotalerr = sqrt(trk1total * trk1total + trk2total * trk2total);
				trgtotalerr = sqrt(trg1total * trg1total + trg2total * trg2total);

				//		cout << "muid1total = " << muid1total << "    muid2total = " <<muid2total << "  muidtotalerr = " << muidtotalerr  << endl;


				TotalWeightMuidUp = TotalWeight * ( 1 + muidtotalerr);
				TotalWeightMuidDown = TotalWeight * ( 1 - muidtotalerr);

				TotalWeightTrkUp = TotalWeight * ( 1 + trktotalerr);
				TotalWeightTrkDown = TotalWeight * ( 1 - trktotalerr);

				TotalWeightTrgUp = TotalWeight * ( 1 + trgtotalerr);
				TotalWeightTrgDown = TotalWeight * ( 1 - trgtotalerr);

				//For the Moment//

				TotalWeightMuidUp = 1;
				TotalWeightMuidDown = 1;


				TotalWeightTrkUp = 1;
				TotalWeightTrkDown = 1;
	
				TotalWeightTrgUp = 1;
				TotalWeightTrgDown = 1;

				TnPWeightHisMuidUp->Fill(Bpt[j],abs(By[j]),TotalWeightMuidUp);
				TnPWeightHisMuidDown->Fill(Bpt[j],abs(By[j]),TotalWeightMuidDown);


				//cout << "TotalWeightMuidUp = " << TotalWeightMuidUp  << "   TotalWeightMuidDown = " << TotalWeightMuidDown << endl;


				TnPWeightHisTrkUp->Fill(Bpt[j],abs(By[j]),TotalWeightTrkUp);
				TnPWeightHisTrkDown->Fill(Bpt[j],abs(By[j]),TotalWeightTrkDown);



				TnPWeightHisTrgUp->Fill(Bpt[j],abs(By[j]),TotalWeightTrgUp);
				TnPWeightHisTrgDown->Fill(Bpt[j],abs(By[j]),TotalWeightTrgDown);


				if(Bpt[j] < ptmax && Bpt[j] > ptmin){
					Bmu1ptHis->Fill(Bmu1pt[j],TotalWeight);
					Bmu2ptHis->Fill(Bmu2pt[j],TotalWeight);
					Bmu1etaHis->Fill(abs(Bmu1eta[j]),TotalWeight);
					Bmu2etaHis->Fill(abs(Bmu2eta[j]),TotalWeight);
					Bmu1TrgSF->Fill(trg1[j],TotalWeight);
					Bmu2TrgSF->Fill(trg2[j],TotalWeight);

					Bmu1pteta->Fill(Bmu1pt[j],abs(Bmu1eta[j]),TotalWeight);
					Bmu2pteta->Fill(Bmu2pt[j],abs(Bmu2eta[j]),TotalWeight);
				}
	
				//BDTWeight
	
				BDTWeight = 1;

				if(Bpt[j] < 7 && Bpt[j] > 5){
					BDTWeightBin = weights_BDT_pt_5_7->GetXaxis()->FindBin(BDT_pt_5_7[j]);
					BDTWeight = weights_BDT_pt_5_7->GetBinContent(BDTWeightBin);
				}	

				if(Bpt[j] < 10 && Bpt[j] > 7){
					BDTWeightBin = weights_BDT_pt_7_10->GetXaxis()->FindBin(BDT_pt_7_10[j]);
					BDTWeight = weights_BDT_pt_7_10->GetBinContent(BDTWeightBin);
				}	
				if(Bpt[j] < 15 && Bpt[j] > 10){
					BDTWeightBin = weights_BDT_pt_10_15->GetXaxis()->FindBin(BDT_pt_10_15[j]);
					BDTWeight = weights_BDT_pt_10_15->GetBinContent(BDTWeightBin);
					
				}
				
				if(Bpt[j] < 20 && Bpt[j] > 15){
					BDTWeightBin = weights_BDT_pt_15_20->GetXaxis()->FindBin(BDT_pt_15_20[j]);
					BDTWeight = weights_BDT_pt_15_20->GetBinContent(BDTWeightBin);
				
				}
				
				if(Bpt[j] < 50 && Bpt[j] > 20){
					BDTWeightBin = weights_BDT_pt_20_50->GetXaxis()->FindBin(BDT_pt_20_50[j]);
					BDTWeight = weights_BDT_pt_20_50->GetBinContent(BDTWeightBin);

				}

				if(Bpt[j] < 60 && Bpt[j] > 50) BDTWeight = 1;

			
				Eff1DRECOHisBDT->Fill(Bpt[j],TotalWeight * BDTWeight);
				BDTWeightHisSyst->Fill(Bpt[j],abs(By[j]),TotalWeight * BDTWeight);

			
        auto iY = std::upper_bound(yBins.begin(), yBins.end(), abs(Bgeny[j]))
          - yBins.begin() - 1;
				BptWeight = BptWtF[iY]->Eval(Bgenpt[j]);

				Eff1DRECOHisBpt->Fill(Bpt[j],TotalWeight * BptWeight);
				BptWeightHisSyst->Fill(Bpt[j],abs(By[j]),TotalWeight * BptWeight);
			


				Eff1DRECOHisBDTMult->Fill(nMult,TotalWeight * BDTWeight);				
				Eff1DRECOHisBptMult->Fill(nMult,TotalWeight * BptWeight);

			}





		}



  }




		cout << "Now Loop Gen" << endl;
		

		
		for(int i = 0; i < NEvents; i++){


			ntGen->GetEntry(i);
			ntHi->GetEntry(i);
			//CentWeightTree->GetEntry(i);
			ntKp->GetEntry(i);
			root->GetEntry(i);

			// PVzWeight = 1;
		//	PVzWeight = (0.0132 * TMath::Exp((PVz -  0.7538) * (PVz -  0.7538)/(2* 6.024* 6.024)))/(  0.0137 * TMath::Exp((PVz -  0.7538) * (PVz -  0.7538)) );
				
			CentWeight = 1;
			
			 

			//PVzWeight = (0.0132 * TMath::Exp((PVz -  0.7538) * (PVz -  0.7538)/(2* 6.024* 6.024)))/(  0.0137 * TMath::Exp((PVz -  0.7538) * (PVz -  0.7538)) );

			//	PVzWeight = (0.163562 * TMath::Exp(- 0.021039 * (PVz - 0.426587)*(PVz - 0.426587)))/(0.159629 * TMath::Exp(- 0.020014 * (PVz - 0.589381)*(PVz - 0.589381)));

			//	PVzWeight = (TMath::Gaus(PVz,0.432315,4.874300)/(sqrt(2*3.14159)*4.874300))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989));
      PVzWeight = (0.013245 * TMath::Exp(-(PVz-0.753876)*(PVz-0.753876)/(2 * 6.023671 * 6.023671)))/(0.013790 * TMath::Exp(-(PVz-0.608178)*(PVz-0.608178)/(2 * 5.785230 * 5.785230)));
			//	EventWeight = PVzWeight * CentWeight * pthat * weight;
			
			//EventWeight = PVzWeight * CentWeight * weight;
			EventWeight = 1;	


			for(int j = 0; j < Gsize; j++){

		
			//	EventWeight = PVzWeight * BptWeight * weight;
				




	//			if((TMath::Abs(Gy[j])<2.4 && TMath::Abs(GpdgId[j])==521 && GisSignal[j]==1 && GcollisionId[j]==0 ) ){
				if( (TMath::Abs(Gy[j])<2.4 && TMath::Abs(GpdgId[j])==521 && GisSignal[j]==1 && GcollisionId[j]==0 )  && ((Gpt[j]>5 && Gpt[j]<10 && TMath::Abs(Gy[j])>1.5) || (Gpt[j]>10))  ){

          auto iY = std::upper_bound(yBins.begin(), yBins.end(), abs(Gy[j]))
            - yBins.begin() - 1;
          BptWeight = BptWtF[iY]->Eval(Gpt[j]);

					NoWeightGenHis->Fill(Gpt[j],abs(Gy[j]),1);
					EvtWeightGenHis->Fill(Gpt[j],abs(Gy[j]),EventWeight);
					BptWeightGenHis->Fill(Gpt[j],abs(Gy[j]),EventWeight * BptWeight);					
					Eff1DGENHis->Fill(Gpt[j],EventWeight);
					Eff1DGENMultHis->Fill(nMult,EventWeight);
					Eff1DGENYHis->Fill(By[j],EventWeight);
					Eff1DGENHisGpt->Fill(Gpt[j],EventWeight * BptWeight);
					Eff1DGENMultHisGpt->Fill(nMult,EventWeight * BptWeight);
				
				}

				if((TMath::Abs(Gy[j])<2.4 && TMath::Abs(GpdgId[j])==521 && GisSignal[j]==1 && GcollisionId[j]==0 && TMath::Abs(Gtk1eta[j])<2.4 && ((TMath::Abs(Gmu1eta[j])<1.2 && Gmu1pt[j]>3.5) || (TMath::Abs(Gmu1eta[j])>1.2 && TMath::Abs(Gmu1eta[j])<2.1 && Gmu1pt[j]>5.47-1.89*TMath::Abs(Gmu1eta[j])) || (TMath::Abs(Gmu1eta[j])>2.1 && TMath::Abs(Gmu2eta[j])<2.4 && Gmu1pt[j]>1.5)) && ((TMath::Abs(Gmu2eta[j])<1.2 && Gmu2pt[j]>3.5) || (TMath::Abs(Gmu2eta[j])>1.2 && TMath::Abs(Gmu2eta[j])<2.1 && Gmu2pt[j]>5.47-1.89*TMath::Abs(Gmu2eta[j])) || (TMath::Abs(Gmu2eta[j])>2.1 && TMath::Abs(Gmu2eta[j])<2.4 && Gmu2pt[j]>1.5)) && ((Gpt[j]>5 && Gpt[j]<10 && TMath::Abs(Gy[j])>1.5) || (Gpt[j]>10))) ){

					NoWeightGenAccHis->Fill(Gpt[j],abs(Gy[j]),EventWeight);
					EvtWeightGenAccHis->Fill(Gpt[j],abs(Gy[j]),EventWeight);
					Eff1DGENAccHis->Fill(Gpt[j],EventWeight);
					Eff1DGENAccYHis->Fill(Gy[j],EventWeight);
					Eff1DGENAccMultHis->Fill(nMult,EventWeight);

				}
			}
		}

		cout << "START MAKING HIS BRO" << endl;

		//TOTAL EFF
		TH1D * Eff1DHis = (TH1D * ) Eff1DRECOHis->Clone("Eff1DHis");
		Eff1DHis->Sumw2();
		Eff1DGENHis->Sumw2();
		Eff1DHis->Divide(Eff1DGENHis);

		//SEL EFF
		TH1D * Sel1DHis = (TH1D * ) Eff1DRECOHis->Clone("Sel1DHis");
		Sel1DHis->Sumw2();
		Eff1DGENAccHis->Sumw2();
		Sel1DHis->Divide(Eff1DGENAccHis);
		
		//ACEPTANCE

		TH1D * Acc1DHis = (TH1D * ) Eff1DGENAccHis->Clone("Acc1DHis");
		Acc1DHis->Sumw2();
		Eff1DGENHis->Sumw2();
		Acc1DHis->Divide(Eff1DGENHis);
	

		//Save 1D Eff Plots/

		Eff1DHis->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
		Eff1DHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
		Eff1DHis->GetYaxis()->SetTitleOffset(1.4);
		Eff1DHis->GetXaxis()->CenterTitle();
		Eff1DHis->GetYaxis()->CenterTitle();
	

		Eff1DHis->SetMarkerStyle(20);
		Eff1DHis->SetMarkerSize(1);
		Eff1DHis->SetMarkerColor(kBlack);
		Eff1DHis->SetLineColor(kBlack);


		Sel1DHis->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
		Sel1DHis->GetYaxis()->SetTitle("#epsilon");
		Sel1DHis->GetYaxis()->SetTitleOffset(1.4);
		Sel1DHis->GetXaxis()->CenterTitle();
		Sel1DHis->GetYaxis()->CenterTitle();


		Sel1DHis->SetMarkerStyle(20);
		Sel1DHis->SetMarkerSize(1);
		Sel1DHis->SetMarkerColor(kBlack);
		Sel1DHis->SetLineColor(kBlack);

	


		Acc1DHis->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
		Acc1DHis->GetYaxis()->SetTitle("#alpha");
		Acc1DHis->GetYaxis()->SetTitleOffset(1.4);
		Acc1DHis->GetXaxis()->CenterTitle();
		Acc1DHis->GetYaxis()->CenterTitle();


		Acc1DHis->SetMarkerStyle(20);
		Acc1DHis->SetMarkerSize(1);
		Acc1DHis->SetMarkerColor(kBlack);
		Acc1DHis->SetLineColor(kBlack);




		//Now Syst
		



		TH1D * Eff1DHisTnPUp = (TH1D * ) Eff1DRECOHisTnPUp->Clone("Eff1DHisTnPUp");
		Eff1DHisTnPUp->Sumw2();
		Eff1DGENHis->Sumw2();
		Eff1DHisTnPUp->Divide(Eff1DGENHis);


		TH1D * Eff1DHisTnPDown = (TH1D * ) Eff1DRECOHisTnPDown->Clone("Eff1DHisTnPDown");
		Eff1DHisTnPDown->Sumw2();
		Eff1DGENHis->Sumw2();
		Eff1DHisTnPDown->Divide(Eff1DGENHis);



		TH1D * Eff1DHisBDT = (TH1D * ) Eff1DRECOHisBDT->Clone("Eff1DHisBDT");
		Eff1DHisBDT->Sumw2();
		Eff1DGENHis->Sumw2();
		Eff1DHisBDT->Divide(Eff1DGENHis);



		TH1D * Eff1DHisBpt = (TH1D * ) Eff1DRECOHisBpt->Clone("Eff1DHisBpt");
		Eff1DHisBpt->Sumw2();
		Eff1DGENHisGpt->Sumw2();
		Eff1DHisBpt->Divide(Eff1DGENHisGpt);


		//Draw Syst//

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

		TLegend* leg = new TLegend(0.50,0.35,0.80,0.60,NULL,"brNDC");
		leg->SetBorderSize(0);
		leg->SetTextSize(0.040);
		leg->SetTextFont(42);
		leg->SetFillStyle(0);
		leg->SetLineWidth(3);
		leg->AddEntry(Eff1DHis,"Nominal","PL");
		leg->AddEntry(Eff1DHisTnPUp,"T&P Variation Up","PL");
		leg->AddEntry(Eff1DHisTnPDown,"T&P Variation Down","PL");
		leg->Draw("same");

		cSyst->SaveAs("Syst/TnPSyst.png");



		Eff1DHisBDT->SetMarkerStyle(20);
		Eff1DHisBDT->SetMarkerSize(1);
		Eff1DHisBDT->SetMarkerColor(kRed);
		Eff1DHisBDT->SetLineColor(kRed);

		Eff1DHis->Draw("ep");
		Eff1DHisBDT->Draw("epSAME");

		TLegend* leg2 = new TLegend(0.50,0.35,0.80,0.60,NULL,"brNDC");
		leg2->SetBorderSize(0);
		leg2->SetTextSize(0.040);
		leg2->SetTextFont(42);
		leg2->SetFillStyle(0);
		leg2->SetLineWidth(3);
		leg2->AddEntry(Eff1DHis,"Nominal","PL");
		leg2->AddEntry(Eff1DHisBDT,"BDT Weighted","PL");
		leg2->Draw("same");


		cSyst->SaveAs("Syst/BDTWeighted.png");






		Eff1DHisBpt->SetMarkerStyle(20);
		Eff1DHisBpt->SetMarkerSize(1);
		Eff1DHisBpt->SetMarkerColor(kRed);
		Eff1DHisBpt->SetLineColor(kRed);

		Eff1DHis->Draw("ep");
		Eff1DHisBpt->Draw("epSAME");

		TLegend* leg3 = new TLegend(0.50,0.35,0.80,0.60,NULL,"brNDC");
		leg3->SetBorderSize(0);
		leg3->SetTextSize(0.040);
		leg3->SetTextFont(42);
		leg3->SetFillStyle(0);
		leg3->SetLineWidth(3);
		leg3->AddEntry(Eff1DHis,"Nominal","PL");
		leg3->AddEntry(Eff1DHisBpt,"Bpt Weighted","PL");
		leg3->Draw("same");

		cSyst->SaveAs("Syst/BptWeighted.png");














		for(int i = 0; i < NPtBins1D; i++){

			float TnPSystValue =  (Eff1DHisTnPUp->GetBinContent(i+1) - Eff1DHis->GetBinContent(i+1))/Eff1DHis->GetBinContent(i+1);


			cout << "i = " << i << "   TnPSystValue = " << TnPSystValue << endl;
		}

		for(int i = 0; i < NPtBins1D; i++){

			float BDTSystValue =  abs(Eff1DHisBDT->GetBinContent(i+1) - Eff1DHis->GetBinContent(i+1))/Eff1DHis->GetBinContent(i+1);

			cout << "i = " << i << "   BDTSystValue = " << BDTSystValue << endl;
			
		}




		TH1D * Eff1DHisMult = (TH1D * ) Eff1DRECOMultHis->Clone("Eff1DHisMult");
		Eff1DHisMult->Sumw2();
		Eff1DGENMultHis->Sumw2();
		Eff1DHisMult->Divide(Eff1DGENMultHis);

		TH1D * Eff1DHisY = (TH1D * ) Eff1DRECOYHis->Clone("Eff1DHisY");
		Eff1DHisY->Sumw2();
		Eff1DGENYHis->Sumw2();
		Eff1DHisY->Divide(Eff1DGENYHis);

		TH1D * Sel1DHisMult = (TH1D * ) Eff1DRECOMultHis->Clone("Sel1DHisMult");
		Sel1DHisMult->Sumw2();
		Eff1DGENAccMultHis->Sumw2();
		Sel1DHisMult->Divide(Eff1DGENAccMultHis);
	
		TH1D * Acc1DHisMult = (TH1D * ) Eff1DGENAccMultHis->Clone("Acc1DHisMult");
		Acc1DHisMult->Sumw2();
		Eff1DGENMultHis->Sumw2();
		Acc1DHisMult->Divide(Eff1DGENMultHis);
	








		//Save 1D Eff Plots/

		Eff1DHisMult->GetXaxis()->SetTitle("B^{+} nMult");
		Eff1DHisMult->GetYaxis()->SetTitle("#alpha #times #epsilon");
		Eff1DHisMult->GetYaxis()->SetTitleOffset(1.4);
		Eff1DHisMult->GetXaxis()->CenterTitle();
		Eff1DHisMult->GetYaxis()->CenterTitle();
	

		Eff1DHisMult->SetMarkerStyle(20);
		Eff1DHisMult->SetMarkerSize(1);
		Eff1DHisMult->SetMarkerColor(kBlack);
		Eff1DHisMult->SetLineColor(kBlack);

		Eff1DHisY->GetXaxis()->SetTitle("B^{+} y (GeV/c)");
		Eff1DHisY->GetYaxis()->SetTitle("#alpha #times #epsilon");
		Eff1DHisY->GetYaxis()->SetTitleOffset(1.4);
		Eff1DHisY->GetXaxis()->CenterTitle();
		Eff1DHisY->GetYaxis()->CenterTitle();
	

		Eff1DHisY->SetMarkerStyle(20);
		Eff1DHisY->SetMarkerSize(1);
		Eff1DHisY->SetMarkerColor(kBlack);
		Eff1DHisY->SetLineColor(kBlack);


		Sel1DHisMult->GetXaxis()->SetTitle("B^{+} nMult ");
		Sel1DHisMult->GetYaxis()->SetTitle("#epsilon");
		Sel1DHisMult->GetYaxis()->SetTitleOffset(1.4);
		Sel1DHisMult->GetXaxis()->CenterTitle();
		Sel1DHisMult->GetYaxis()->CenterTitle();


		Sel1DHisMult->SetMarkerStyle(20);
		Sel1DHisMult->SetMarkerSize(1);
		Sel1DHisMult->SetMarkerColor(kBlack);
		Sel1DHisMult->SetLineColor(kBlack);

	


		Acc1DHisMult->GetXaxis()->SetTitle("B^{+} nMult ");
		Acc1DHisMult->GetYaxis()->SetTitle("#alpha");
		Acc1DHisMult->GetYaxis()->SetTitleOffset(1.4);
		Acc1DHisMult->GetXaxis()->CenterTitle();
		Acc1DHisMult->GetYaxis()->CenterTitle();


		Acc1DHisMult->SetMarkerStyle(20);
		Acc1DHisMult->SetMarkerSize(1);
		Acc1DHisMult->SetMarkerColor(kBlack);
		Acc1DHisMult->SetLineColor(kBlack);







		//Now Syst Mult
		



		TH1D * Eff1DHisTnPUpMult = (TH1D * ) Eff1DRECOHisTnPUpMult->Clone("Eff1DHisTnPUpMult");
		Eff1DHisTnPUpMult->Sumw2();
		Eff1DGENMultHis->Sumw2();
		Eff1DHisTnPUpMult->Divide(Eff1DGENMultHis);


		TH1D * Eff1DHisTnPDownMult = (TH1D * ) Eff1DRECOHisTnPDownMult->Clone("Eff1DHisTnPDownMult");
		Eff1DHisTnPDownMult->Sumw2();
		Eff1DGENMultHis->Sumw2();
		Eff1DHisTnPDownMult->Divide(Eff1DGENMultHis);

		TH1D * Eff1DHisTnPUpY = (TH1D * ) Eff1DRECOHisTnPUpY->Clone("Eff1DHisTnPUpY");
		Eff1DHisTnPUpY->Sumw2();
		Eff1DGENYHis->Sumw2();
		Eff1DHisTnPUpY->Divide(Eff1DGENYHis);

		TH1D * Eff1DHisTnPDownY = (TH1D * ) Eff1DRECOHisTnPDownY->Clone("Eff1DHisTnPDownY");
		Eff1DHisTnPDownY->Sumw2();
		Eff1DGENYHis->Sumw2();
		Eff1DHisTnPDownY->Divide(Eff1DGENYHis);

		TH1D * Eff1DHisBDTMult = (TH1D * ) Eff1DRECOHisBDTMult->Clone("Eff1DHisBDTMult");
		Eff1DHisBDTMult->Sumw2();
		Eff1DGENMultHis->Sumw2();
		Eff1DHisBDTMult->Divide(Eff1DGENMultHis);

		TH1D * Eff1DHisBptMult = (TH1D * ) Eff1DRECOHisBptMult->Clone("Eff1DHisBptMult");
		Eff1DHisBptMult->Sumw2();
		Eff1DGENMultHisGpt->Sumw2();
		Eff1DHisBptMult->Divide(Eff1DGENMultHisGpt);


		//Draw Syst//


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

		Eff1DHisTnPUpY->SetMarkerStyle(20);
		Eff1DHisTnPUpY->SetMarkerSize(1);
		Eff1DHisTnPUpY->SetMarkerColor(kRed);
		Eff1DHisTnPUpY->SetLineColor(kRed);


		Eff1DHisTnPDownY->SetMarkerStyle(20);
		Eff1DHisTnPDownY->SetMarkerSize(1);
		Eff1DHisTnPDownY->SetMarkerColor(kBlue);
		Eff1DHisTnPDownY->SetLineColor(kBlue);

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

		cSyst->SaveAs("Syst/TnPSystMult.png");



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


		cSyst->SaveAs("Syst/BDTWeightedMult.png");






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

		cSyst->SaveAs("Syst/BptWeightedMult.png");

		TCanvas * cSyst_y  = new TCanvas("cSyst_y","cSyst_y",600,600);
		cSyst_y->cd();

		Eff1DHisTnPUpY->SetMarkerStyle(20);
		Eff1DHisTnPUpY->SetMarkerSize(1);
		Eff1DHisTnPUpY->SetMarkerColor(kRed);
		Eff1DHisTnPUpY->SetLineColor(kRed);

		Eff1DHisTnPDownY->SetMarkerStyle(20);
		Eff1DHisTnPDownY->SetMarkerSize(1);
		Eff1DHisTnPDownY->SetMarkerColor(kBlue);
		Eff1DHisTnPDownY->SetLineColor(kBlue);

		Eff1DHisTnPUpMult->Draw("ep");
		Eff1DHisMult->Draw("epSAME");
		Eff1DHisTnPDownMult->Draw("epSAME");

		TLegend* legY = new TLegend(0.50,0.35,0.80,0.60,NULL,"brNDC");
		legY->SetBorderSize(0);
		legY->SetTextSize(0.040);
		legY->SetTextFont(42);
		legY->SetFillStyle(0);
		legY->SetLineWidth(3);
		legY->AddEntry(Eff1DHisY,"Nominal","PL");
		legY->AddEntry(Eff1DHisTnPUpY,"T&P Variation Up","PL");
		legY->AddEntry(Eff1DHisTnPDownY,"T&P Variation Down","PL");
		legY->Draw("same");

		cSyst_y->SaveAs("Syst/TnPSystY.png");

		//2D shits

		TH2D * invAcc2D = (TH2D * ) EvtWeightGenHis->Clone("invAcc2D");
		invAcc2D->Sumw2();
		invAcc2D->Divide(EvtWeightGenAccHis);

		TH2D * invEffonly2D = (TH2D * ) EvtWeightGenAccHis->Clone("invEffonly2D");
		invEffonly2D->Sumw2();
		invEffonly2D->Divide(TnPWeightHis);

		TH2D * invEff2D = (TH2D * ) EvtWeightGenHis->Clone("invEff2D");
		invEff2D->Sumw2();
		invEff2D->Divide(TnPWeightHis);

    // Track selection variation
		TH2D * invEffTrkLoose = (TH2D * ) EvtWeightGenHis->Clone("invEffTrkLoose");
		invEffTrkLoose->Sumw2();
		invEffTrkLoose->Divide(TrkLooseHis);

		TH2D * invEffTrkTight = (TH2D * ) EvtWeightGenHis->Clone("invEffTrkTight");
		invEffTrkTight->Sumw2();
		invEffTrkTight->Divide(TrkTightHis);

		//Systematics




		TH2D * invEff2DTnPSystUp = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTnPSystUp");
		invEff2DTnPSystUp->Sumw2();
		invEff2DTnPSystUp->Divide(TnPWeightHisSystDown);

		TH2D * invEff2DTnPSystDown = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTnPSystDown");
		invEff2DTnPSystDown->Sumw2();
		invEff2DTnPSystDown->Divide(TnPWeightHisSystUp);


		TH2D * invEff2DMuidUp = (TH2D * ) EvtWeightGenHis->Clone("invEff2DMuidUp");
		invEff2DMuidUp->Sumw2();
		invEff2DMuidUp->Divide(TnPWeightHisMuidDown);

		TH2D * invEff2DMuidDown = (TH2D * ) EvtWeightGenHis->Clone("invEff2DMuidDown");
		invEff2DMuidDown->Sumw2();
		invEff2DMuidDown->Divide(TnPWeightHisMuidUp);


		TH2D * invEff2DTrkUp = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTrkUp");
		invEff2DTrkUp->Sumw2();
		invEff2DTrkUp->Divide(TnPWeightHisTrkDown);

		TH2D * invEff2DTrkDown = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTrkDown");
		invEff2DTrkDown->Sumw2();
		invEff2DTrkDown->Divide(TnPWeightHisTrkUp);



		TH2D * invEff2DTrgUp = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTrgUp");
		invEff2DTrgUp->Sumw2();
		invEff2DTrgUp->Divide(TnPWeightHisTrgDown);

		TH2D * invEff2DTrgDown = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTrgDown");
		invEff2DTrgDown->Sumw2();
		invEff2DTrgDown->Divide(TnPWeightHisTrgUp);


		TH2D * invEff2DBDTSyst = (TH2D * ) EvtWeightGenHis->Clone("invEff2DBDTSyst");
		invEff2DBDTSyst->Sumw2();
		BDTWeightHisSyst->Sumw2();
		invEff2DBDTSyst->Divide(BDTWeightHisSyst);


		TH2D * invEff2DBptSyst = (TH2D * ) BptWeightGenHis->Clone("invEff2DBptSyst");
		invEff2DBptSyst->Sumw2();
		BptWeightHisSyst->Sumw2();
		invEff2DBptSyst->Divide(BptWeightHisSyst);



		NoWeightHis->Write();
		EvtWeightHis->Write();
		muidWeightHis->Write();
		trkWeightHis->Write();
		muidtrkWeightHis->Write();
		TnPWeightHis->Write();

		NoWeightGenHis->Write();
		EvtWeightGenHis->Write();

		NoWeightGenAccHis->Write();
		EvtWeightGenAccHis->Write();

		invAcc2D->Write();
		invEffonly2D->Write();
		invEff2D->Write();
		invEffTrkTight->Write();
		invEffTrkLoose->Write();

    	// debugging purpose
    	BDTWeightHisSyst->Write();
    	BptWeightHisSyst->Write();

		invEff2DTnPSystUp->Write();
		invEff2DTnPSystDown->Write();

		invEff2DMuidUp->Write();
		invEff2DMuidDown->Write();

		invEff2DTrkUp->Write();
		invEff2DTrkDown->Write();

		invEff2DTrgUp->Write();
		invEff2DTrgDown->Write();
	

		Eff1DRECOHis->Write();
		Eff1DGENHis->Write();
		Eff1DHis->Write();


		Eff1DRECOMultHis->Write();
		Eff1DGENMultHis->Write();
		Eff1DHisMult->Write();
		

		Eff1DRECOYHis->Write();
		Eff1DGENYHis->Write();
		Eff1DHisY->Write();

		TFile * fout2 = new TFile(Form("BPMuonInfoPlots_%d_%d.root",ptmin,ptmax),"RECREATE");
		fout2->cd();

		TCanvas *c = new TCanvas("c","c",600,600);
		c->cd();

		Eff1DHis->Draw("ep");
		c->SaveAs("1DEffPlots/Eff1DHis.png");

		for(int i = 0; i < NPtBins; i ++){

			cout << "pT =  " << PtBin[i] << " - " << PtBin[i+1] << "  Eff:  " << Eff1DHis->GetBinContent(i+1)   << "   Eff Error: " << Eff1DHis->GetBinError(i+1)<< endl;


		}



		TnPSFHis->Draw();
		c->SaveAs("TnPHis/TnPSFHis.png");

		TnPSFHisMu1->Draw();
		c->SaveAs("TnPHis/TnPSFHisMu1.png");

		TnPSFHisMu2->Draw();
		c->SaveAs("TnPHis/TnPSFHisMu2.png");



		Bmu1ptHis->Draw();
		c->SaveAs("MuonInfoPlots/Bmu1ptHis.png");

		Bmu2ptHis->Draw();
		c->SaveAs("MuonInfoPlots/Bmu2ptHis.png");


		Bmu1etaHis->Draw();
		c->SaveAs("MuonInfoPlots/Bmu1etaHis.png");

		Bmu2etaHis->Draw();
		c->SaveAs("MuonInfoPlots/Bmu2etaHis.png");

		Bmu1TrgSF->Draw();
		c->SaveAs("MuonInfoPlots/Bmu1TrgSF.png");

		Bmu2TrgSF->Draw();
		c->SaveAs("MuonInfoPlots/Bmu2TrgSF.png");

		Bmu1ptHis->Write();
		Bmu2ptHis->Write();
		Bmu1etaHis->Write();
		Bmu2etaHis->Write();
		Bmu1TrgSF->Write();
		Bmu2TrgSF->Write();

		Bmu1pteta->Write();
		Bmu2pteta->Write();
		fout2->Close();


		c->SetLogz();

		invEff2DTnPSystUp->GetYaxis()->SetTitle("B |y|");
		invEff2DTnPSystUp->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
		invEff2DTnPSystUp->GetXaxis()->CenterTitle();
		invEff2DTnPSystUp->GetYaxis()->CenterTitle();
		invEff2DTnPSystUp->GetYaxis()->SetTitleOffset(1.2);
		invEff2DTnPSystUp->SetTitle("");

		invEff2DTnPSystUp->Draw("COLZ");
		c->SaveAs("Eff2DMapTnP/Eff2D_Up.png");




		invEff2DTnPSystDown->GetYaxis()->SetTitle("B |y|");
		invEff2DTnPSystDown->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
		invEff2DTnPSystDown->GetXaxis()->CenterTitle();
		invEff2DTnPSystDown->GetYaxis()->CenterTitle();
		invEff2DTnPSystDown->GetYaxis()->SetTitleOffset(1.2);
		invEff2DTnPSystDown->SetTitle("");

		invEff2DTnPSystDown->Draw("COLZ");
		c->SaveAs("Eff2DMapTnP/Eff2D_Down.png");


		TCanvas * c1DSave = new TCanvas("c1DSave","c1DSave",600,600);
		c1DSave->cd();
		

		Acc1DHis->Draw("ep");
		c1DSave->SaveAs("Plot1DEfficiency/Pt/Acc1DHis.png");
	
		Sel1DHis->Draw("ep");
		c1DSave->SaveAs("Plot1DEfficiency/Pt/Sel1DHis.png");

		Eff1DHis->Draw("ep");
		c1DSave->SaveAs("Plot1DEfficiency/Pt/Eff1DHis.png");

		Acc1DHisMult->Draw("ep");
		c1DSave->SaveAs("Plot1DEfficiency/Mult/Acc1DHis.png");
	
		Sel1DHisMult->Draw("ep");
		c1DSave->SaveAs("Plot1DEfficiency/Mult/Sel1DHis.png");

		Eff1DHisMult->Draw("ep");
		c1DSave->SaveAs("Plot1DEfficiency/Mult/Eff1DHis.png");



		TFile * foutSyst = new TFile("NewEff2DMaps/BPSyst.root","RECREATE");
		foutSyst->cd();
		Eff1DHis->Write();
		Eff1DHisTnPUp->Write();
		Eff1DHisTnPDown->Write();
		Eff1DHisBpt->Write();
		Eff1DHisBDT->Write();

		Eff1DHisY->Write();
		Eff1DHisMult->Write();
		Eff1DHisTnPUpMult->Write();
		Eff1DHisTnPDownMult->Write();
		Eff1DHisBptMult->Write();
		Eff1DHisBDTMult->Write();
		foutSyst->Close();



		TFile * foutSyst2D = new TFile("NewEff2DMaps/BPSyst2D.root","RECREATE");

		invEff2D->Write();
		invEff2DTnPSystUp->Write();
		invEff2DTnPSystDown->Write();
		invEff2DBDTSyst->Write();
		invEff2DBptSyst->Write();

		foutSyst2D->Close();

	

		fout->Close();
		fin->Close();





	}
