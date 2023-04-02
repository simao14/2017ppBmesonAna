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
#include "../../henri2022/parameter.h" 
//#include "tnp_weight_lowptPbPb.h"

//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void CrossSectionAnaMult(int DoTnP,int whichvar, int usemc=0){

	gSystem->mkdir("EffFinal" ,true);
	gSystem->mkdir("FinalFiles" ,true);

	int NBins;
	//const int NBins = 6;

	int TnP = 1;


	double BRchain = 6.02061e-5;

	//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);


	TString FileName;
	TFile * fin;
	TTree * EffInfoTree;
	TTree * root;
	int NCand;


	if (usemc==0){
	FileName = "/data3/tasheng/presel/BPData_nom.root";
	fin = new TFile(FileName.Data());
	EffInfoTree = (TTree * ) fin->Get("ntKp");
	NCand=10;
	}
	else {
	FileName = "/data3/tasheng/presel/output/BP_MC_BDTs_nom_tnp.root";
	fin = new TFile(FileName.Data());
	EffInfoTree = (TTree * ) fin->Get("Bfinder/ntKp");
	root = (TTree * ) fin->Get("Bfinder/root");
	NCand=13000;
	}
	
	fin->cd();

	int NEvents = EffInfoTree->GetEntries();
	if (usemc==1){NEvents=NEvents/100;}
	
	Int_t BsizeNew;
	Int_t runNew;
	Int_t lumiNew;
	Int_t evtNew;
	Float_t BmassNew[NCand];
	Float_t BptNew[NCand];
	Float_t ByNew[NCand];
	Float_t BEff[NCand];
	Float_t BEffErr[NCand];

	Int_t nMult;
	Int_t trackSelection;

	Float_t BEffInv[NCand];
	Float_t BEffInvErr[NCand];
	Float_t BEffInv1D[NCand];
	Float_t BEffInvErr1D[NCand];

	Float_t BEffInvFit[NCand];
	Float_t BEffInvErrFit[NCand];

	Float_t BEffInvBDTWeighted[NCand];
	Float_t BEffInvErrBDTWeighted[NCand];



	Float_t BEffInvUp[NCand];
	Float_t BEffInvErrUp[NCand];
	Float_t BEffInvDown[NCand];
	Float_t BEffInvErrDown[NCand];

	EffInfoTree->SetBranchAddress("Bsize",&BsizeNew);
	EffInfoTree->SetBranchAddress("Bmass",BmassNew);
	EffInfoTree->SetBranchAddress("By",ByNew);
	EffInfoTree->SetBranchAddress("Bpt",BptNew);
	if (usemc==0){
			EffInfoTree->SetBranchAddress("nMult",&nMult); 
			EffInfoTree->SetBranchAddress("track", &trackSelection);
			}
	else {root->SetBranchAddress("EvtInfo.nMult",&nMult);}
	

	Float_t var;
	TString var_m;
	
	if (whichvar==0){
		var_m="y"; 
		NBins=8;} 
	if (whichvar==1){
		var_m="Mult";
		NBins=7;}
	if (whichvar==2){
		var_m="pt";
		NBins=7;}

	double ptBins[NBins + 1];

	int Counts[NBins];
	double SumCounts[NBins];
	double SumCountsErr[NBins];
	double NewEff[NBins];
	double NewEffErr[NBins];

	double NewEffReal[NBins];
	double NewEffRealErr[NBins];



	double SumCountsUp[NBins];
	double SumCountsErrUp[NBins];

	double SumCountsDown[NBins];
	double SumCountsErrDown[NBins];


	double SumCountsEff[NBins];
	double SumCountsEffErr[NBins];



	double SumCountsSyst[NBins];
	double SumCountsSystErr[NBins];
	double NewEffSyst[NBins];
	double NewEffSystErr[NBins];


	//	double CorrectionFactor[NBins];

	double NewEffUp[NBins];
	double NewEffErrUp[NBins];
	double NewEffDown[NBins];
	double NewEffErrDown[NBins];


	double lumi = 302.3;
	if (whichvar==0){
	for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ybinsvec[i];             //taken from parameter.h
		}
	} 
	if (whichvar==1){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  nmbinsvec[i];             //taken from parameter.h
		}
	}
	if (whichvar==2){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvecBP[i];             //taken from parameter.h
		}
	}
	for(int i = 0; i < NBins; i++){
		Counts[i] = 0;
		SumCounts[i] = 0;
		SumCountsErr[i] = 0;
		SumCountsEff[i] = 0;
		SumCountsEffErr[i] = 0;
		SumCountsSyst[i] = 0;
		SumCountsSystErr[i] = 0;
		SumCountsUp[i] = 0;
		SumCountsErrUp[i] = 0;
		SumCountsDown[i] = 0;
		SumCountsErrDown[i] = 0;

	}





	int EtaBin;
	int PtBin;

	double trgtnp1;
	double trktnp1;
	double muidtnp1;

	double trgtnp1systup;
	double trgtnp1systdown;
	double trgtnp1statup;
	double trgtnp1statdown;

	double trktnp1systup;
	double trktnp1systdown;
	double trktnp1statup;
	double trktnp1statdown;

	double muidtnp1systup;
	double muidtnp1systdown;
	double muidtnp1statup;
	double muidtnp1statdown;

	double tnptotal1;
	double tnptotal1up;
	double tnptotal1down;

	double tnptotal1systup;
	double tnptotal1systdown;
	double tnptotal1statup;
	double tnptotal1statdown;

	double trgtnp2;
	double trktnp2;
	double muidtnp2;

	double trgtnp2systup;
	double trgtnp2systdown;
	double trgtnp2statup;
	double trgtnp2statdown;


	double trktnp2systup;
	double trktnp2systdown;
	double trktnp2statup;
	double trktnp2statdown;

	double muidtnp2systup;
	double muidtnp2systdown;
	double muidtnp2statup;
	double muidtnp2statdown;

	double tnptotal2;
	double tnptotal2up;
	double tnptotal2down;

	double tnptotal2systup;
	double tnptotal2systdown;
	double tnptotal2statup;
	double tnptotal2statdown;

	double tnpabssystup;
	double tnpabssystdown;


	//pt Binned Correction//
	TFile * fin1DEff; 

	if(DoTnP == 0) fin1DEff = new TFile("NewEff2DMaps/EffFineNoTnP.root");
	if(DoTnP == 1) fin1DEff = new TFile("NewEff2DMaps/EffFineBDT.root");	


	fin1DEff->cd();


	//Add 2D eff calculations
	
	TH2D * invEff2D = (TH2D *) fin1DEff->Get("invEff2D");
	TH2D * invEffTrkTight = (TH2D *) fin1DEff->Get("invEffTrkTight");
	TH2D * invEffTrkLoose = (TH2D *) fin1DEff->Get("invEffTrkLoose");

	Float_t EffInvTrkTight;
  	Float_t EffInvTrkLoose;
	int XBin;
	int YBin;
	if (usemc==0){
		for( int i = 0; i < NEvents; i++){

			EffInfoTree->GetEntry(i);
			if (whichvar==0){var=ByNew[0]; trackSelection=1;}
			if (whichvar==1){var=nMult; trackSelection=1;}
			if (whichvar==2){var=BptNew[0];}
			int j=0;
				for(int k = 0; k < NBins; k++){

					if(var > ptBins[k] && var < ptBins[k+1] && TMath::Abs(BmassNew[j] - 5.27932) < 0.08 &&  TMath::Abs(ByNew[j]) < 2.4  && ((BptNew[j] > 5 && BptNew[j] < 10 && abs(ByNew[j]) > 1.5 )||(BptNew[j] > 10)))
					{

						
						XBin = invEff2D->GetXaxis()->FindBin( BptNew[j]);
						YBin = invEff2D->GetYaxis()->FindBin( TMath::Abs(ByNew[j]));
						BEffInv[j] = invEff2D->GetBinContent(XBin,YBin);

						BEffInvErr[j] = invEff2D->GetBinError(XBin,YBin);
						BEff[j] = 1.0/invEff2D->GetBinContent(XBin,YBin);

						BEffErr[j] = BEffInvErr[j]/(BEffInv[j] * BEffInv[j]);

						if(trackSelection>0 && BEffInv[j] > 0){
							SumCounts[k] = SumCounts[k] + BEffInv[j];
							SumCountsErr[k] = SumCountsErr[k] + BEffInvErr[j] * BEffInvErr[j];
							SumCountsEff[k] = SumCountsEff[k] + BEff[j];
							SumCountsEffErr[k] = SumCountsEffErr[k] + BEffErr[j] * BEffErr[j];
							SumCountsSyst[k] = 	SumCountsSyst[k]  + BEffInvBDTWeighted[j];
							SumCountsSystErr[k] = 	SumCountsSystErr[k]  + BEffInvErrBDTWeighted[j] * BEffInvErrBDTWeighted[j];

							SumCountsUp[k] = SumCountsUp[k] + BEffInvUp[j];
							SumCountsErrUp[k] = SumCountsErrUp[k] + BEffInvErrUp[j] * BEffInvErrUp[j];

							SumCountsDown[k] = SumCountsDown[k] + BEffInvDown[j];
							SumCountsErrDown[k] = SumCountsErrUp[k] + BEffInvErrDown[j] * BEffInvErrDown[j];

							Counts[k] = Counts[k] + 1;
						}
					}
				}
			}
	}
	else {
		for( int i = 0; i < NEvents; i++){
			EffInfoTree->GetEntry(i);
			root->GetEntry(i);
			for (int j = 0 ; j< BsizeNew ;j++){
				if (whichvar==0){var=ByNew[j];}
				if (whichvar==1){var=nMult;}
				if (whichvar==2){var=BptNew[j];}
				for(int k = 0; k < NBins; k++){
					if(var > ptBins[k] && var < ptBins[k+1] && TMath::Abs(BmassNew[j] - 5.27932) < 0.08 &&  TMath::Abs(ByNew[j]) < 2.4  && ((BptNew[j] > 5 && BptNew[j] < 10 && abs(ByNew[j]) > 1.5 )||(BptNew[j] > 10)))
					{
						XBin = invEff2D->GetXaxis()->FindBin( BptNew[j]);
						YBin = invEff2D->GetYaxis()->FindBin( TMath::Abs(ByNew[j]));
						BEffInv[j] = invEff2D->GetBinContent(XBin,YBin);

						BEffInvErr[j] = invEff2D->GetBinError(XBin,YBin);
						BEff[j] = 1.0/invEff2D->GetBinContent(XBin,YBin);

						BEffErr[j] = BEffInvErr[j]/(BEffInv[j] * BEffInv[j]);

						if(BEffInv[j] > 0){
							SumCounts[k] = SumCounts[k] + BEffInv[j];
							SumCountsErr[k] = SumCountsErr[k] + BEffInvErr[j] * BEffInvErr[j];
							SumCountsEff[k] = SumCountsEff[k] + BEff[j];
							SumCountsEffErr[k] = SumCountsEffErr[k] + BEffErr[j] * BEffErr[j];
							SumCountsSyst[k] = 	SumCountsSyst[k]  + BEffInvBDTWeighted[j];
							SumCountsSystErr[k] = 	SumCountsSystErr[k]  + BEffInvErrBDTWeighted[j] * BEffInvErrBDTWeighted[j];

							SumCountsUp[k] = SumCountsUp[k] + BEffInvUp[j];
							SumCountsErrUp[k] = SumCountsErrUp[k] + BEffInvErrUp[j] * BEffInvErrUp[j];

							SumCountsDown[k] = SumCountsDown[k] + BEffInvDown[j];
							SumCountsErrDown[k] = SumCountsErrUp[k] + BEffInvErrDown[j] * BEffInvErrDown[j];

							Counts[k] = Counts[k] + 1;
						}
					}
				}
			}
		  }
		}


	TH1D * hInvEff = new TH1D("hInvEff","",NBins,ptBins);


	hInvEff->GetXaxis()->SetTitle(Form("B^{+} %s (GeV/c)",var_m.Data()));
	hInvEff->GetYaxis()->SetTitle("<1/(Eff * Acc)>");
	hInvEff->GetYaxis()->SetTitleOffset(1.4);
	hInvEff->GetXaxis()->CenterTitle();
	hInvEff->GetYaxis()->CenterTitle();
	hInvEff->SetMarkerColor(1);
	hInvEff->SetLineColor(1);
	hInvEff->SetMarkerStyle(20);

	hInvEff->SetMinimum(0);


	TH1D * hInvEffSyst = new TH1D("hInvEffSyst","",NBins,ptBins);

	hInvEffSyst->GetXaxis()->SetTitle(Form("B^{+} %s (GeV/c)",var_m.Data()));
	hInvEffSyst->GetYaxis()->SetTitle("<1/(Eff * Acc)> - BDT Data-MC Weighted");
	hInvEffSyst->GetYaxis()->SetTitleOffset(1.4);
	hInvEffSyst->GetXaxis()->CenterTitle();
	hInvEffSyst->GetYaxis()->CenterTitle();
	hInvEffSyst->SetMarkerColor(1);
	hInvEffSyst->SetLineColor(2);
	hInvEffSyst->SetMarkerStyle(20);

	hInvEffSyst->SetMinimum(0);

	TH1D * hEff = new TH1D("hEff","",NBins,ptBins);


	hEff->GetXaxis()->SetTitle(Form("B^{+} %s (GeV/c)",var_m.Data()));
	hEff->GetYaxis()->SetTitle("<(Eff * Acc)>");
	hEff->GetYaxis()->SetTitleOffset(1.4);
	hEff->GetXaxis()->CenterTitle();
	hEff->GetYaxis()->CenterTitle();
	hEff->SetMarkerColor(1);
	hEff->SetLineColor(1);
	hEff->SetMarkerStyle(20);

	hEff->SetMinimum(0);

	TH1D * hEffInv = new TH1D("hEffInv","",NBins,ptBins);


	hEffInv->GetXaxis()->SetTitle(Form("B^{+} %s (GeV/c)",var_m.Data()));
	hEffInv->GetYaxis()->SetTitle("<(Eff * Acc)> 2D Inv");
	hEffInv->GetYaxis()->SetTitleOffset(1.4);
	hEffInv->GetXaxis()->CenterTitle();
	hEffInv->GetYaxis()->CenterTitle();
	hEffInv->SetMarkerColor(1);
	hEffInv->SetLineColor(1);
	hEffInv->SetMarkerStyle(20);

	hEffInv->SetMinimum(0);


	TH1D * hInvEffUp = new TH1D("hInvEffUp","",NBins,ptBins);


	hInvEffUp->GetXaxis()->SetTitle(Form("B^{+} %s (GeV/c)",var_m.Data()));
	hInvEffUp->GetYaxis()->SetTitle("<1./(Eff * Acc)>");
	hInvEffUp->GetYaxis()->SetTitleOffset(1.4);
	hInvEffUp->GetXaxis()->CenterTitle();
	hInvEffUp->GetYaxis()->CenterTitle();
	hInvEffUp->SetMarkerColor(1);
	hInvEffUp->SetLineColor(1);
	hInvEffUp->SetMarkerStyle(20);
	hInvEffUp->SetMinimum(0);



	TH1D * hInvEffDown = new TH1D("hInvEffDown","",NBins,ptBins);


	hInvEffDown->GetXaxis()->SetTitle(Form("B^{+} %s (GeV/c)",var_m.Data()));
	hInvEffDown->GetYaxis()->SetTitle("<1/(Eff * Acc)>");
	hInvEffDown->GetYaxis()->SetTitleOffset(1.4);
	hInvEffDown->GetXaxis()->CenterTitle();
	hInvEffDown->GetYaxis()->CenterTitle();
	hInvEffDown->SetMarkerColor(1);
	hInvEffDown->SetLineColor(1);
	hInvEffDown->SetMarkerStyle(20);
	hInvEffDown->SetMinimum(0);



	for(int i = 0; i < NBins; i++){


		NewEff[i] = SumCounts[i]/Counts[i];
		NewEffErr[i] = TMath::Sqrt(SumCountsErr[i])/Counts[i];


		NewEffUp[i] = SumCountsUp[i]/Counts[i];
		NewEffErrUp[i] = TMath::Sqrt(SumCountsErrUp[i])/Counts[i];



		NewEffDown[i] = SumCountsDown[i]/Counts[i];
		NewEffErrDown[i] = TMath::Sqrt(SumCountsErrDown[i])/Counts[i];


		NewEffReal[i] = SumCountsEff[i]/Counts[i];
		NewEffRealErr[i] = TMath::Sqrt(SumCountsEffErr[i])/Counts[i];


		hInvEff->SetBinContent(i+1,NewEff[i]);
		hInvEff->SetBinError(i+1,NewEffErr[i]);

		hEff->SetBinContent(i+1,NewEffReal[i]);
		hEff->SetBinError(i+1,NewEffRealErr[i]);

		hEffInv->SetBinContent(i+1,1/NewEff[i]);
		hEffInv->SetBinError(i+1,NewEffErr[i]/(NewEff[i] * NewEff[i]));

		NewEffSyst[i] = SumCountsSyst[i]/Counts[i];
		NewEffSystErr[i] = TMath::Sqrt(SumCountsSystErr[i])/Counts[i];


		hInvEffSyst->SetBinContent(i+1,	NewEffSyst[i]);
		hInvEffSyst->SetBinError(i+1, NewEffSystErr[i]);



		hInvEffUp->SetBinContent(i+1,NewEffUp[i]);
		hInvEffUp->SetBinError(i+1,NewEffErrUp[i]);


		hInvEffDown->SetBinContent(i+1,NewEffDown[i]);
		hInvEffDown->SetBinError(i+1,NewEffErrDown[i]);

		//	cout << "Real eff = " << SumCountsReal[i]/Counts[i] << endl;
		//cout << "Counts = " << Counts[i] << endl;
		cout << "Count =  " <<  Counts[i] << "   NewEff = " << NewEff[i] << "     NewEffErr = " << NewEffErr[i] << endl;
		cout << "Count =  " <<  Counts[i] << "   NewEffSyst = " << NewEffSyst[i] << "     NewEffSystErr = " << NewEffSystErr[i] << endl;



		cout << "-----------------------------------------------------------------------------------------------" << endl;

		cout << "   NewEff = " << NewEff[i] << "     NewEffErr = " << NewEffErr[i] << "  Fractional = " << NewEffErr[i]/NewEff[i] << endl;
		//	cout << "   NewEff = " << NewEffUp[i] << "     NewEffErr = " << NewEffErrUp[i] << "  Fractional = " << NewEffErrUp[i]/NewEffUp[i] << endl;



		//NewEffErr[i] = 0; //Remove Error on Efficiency Correction//
	}



	TFile * RawYield = new TFile(Form("../../henri2022/ROOTfiles/yields_Bp_binned_%s.root",var_m.Data()));
	RawYield->cd();
	TH1D * hPt = (TH1D *) RawYield->Get("hPt");


	double RawCount;
	double RawCountErr;

	double CorrYield= 0;
	double CorrYieldErr= 0;

	double CorrYieldDiff[NBins];
	double CorrYieldDiffErr[NBins];

	for(int i = 0; i < NBins;i++){

		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);

		cout << "RawCount = " << RawCount << "  RawCountErr = " << RawCountErr << " NewEff[i] =   " << NewEff[i] << "  NewEffErr[i] =  " << NewEffErr[i] << endl; 

		cout << "CORR YIELD PT:  " <<  RawCount *   NewEff[i]/(BRchain*2 * lumi) << endl;
		CorrYield = RawCount * (ptBins[i+1] - ptBins[i]) *  NewEff[i]  + CorrYield;
		CorrYieldErr = ((RawCountErr * (ptBins[i+1] - ptBins[i]) *  NewEff[i]) *(RawCountErr * (ptBins[i+1] - ptBins[i]) *  NewEff[i]) + (RawCount * (ptBins[i+1] - ptBins[i]) *  NewEffErr[i]) * (RawCount * (ptBins[i+1] - ptBins[i]) *  NewEffErr[i]))  + CorrYieldErr;

		cout << "PrintYield = " << RawCount* (ptBins[i+1] - ptBins[i]) *  NewEff[i] << endl;

	}

	TH1D * CorrDiffHis = new TH1D("hPtSigma","",NBins,ptBins);
	CorrDiffHis->GetXaxis()->SetTitle(Form("%s (GeV/c)",var_m.Data()));
	CorrDiffHis->GetYaxis()->SetTitle("d #sigma/d p_{T} (pb GeV^{-1} c)");

	CorrDiffHis->GetYaxis()->SetTitleOffset(1.3);
	CorrDiffHis->GetXaxis()->CenterTitle();
	CorrDiffHis->GetYaxis()->CenterTitle();




	for(int i = 0; i < NBins;i++){
		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);
		CorrYieldDiff[i] = (RawCount *  NewEff[i])/(BRchain*2* lumi);
		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr *  NewEff[i]) *(RawCountErr  *  NewEff[i]) + (RawCount *  NewEffErr[i]) * (RawCount  *  NewEffErr[i]))/(BRchain*2* lumi);
		CorrDiffHis->SetBinContent(i+1,CorrYieldDiff[i]);
		CorrDiffHis->SetBinError(i+1,CorrYieldDiffErr[i]);

	}

	
	//CorrDiffHis->SetTitle("(Preliminary) B^{+} #rightarrow J/#psi K^{+} p_{T} Differential Cross Section in pp");

	CorrDiffHis->SetMarkerColor(kBlack);
	CorrDiffHis->SetMarkerSize(1);
	CorrDiffHis->SetMarkerStyle(20);

	TH1D * CorrDiffHisReal = new TH1D("hPtSigmaReal","",NBins,ptBins);
	CorrDiffHisReal->GetXaxis()->SetTitle(Form("%s (GeV/c)",var_m.Data()));
	CorrDiffHisReal->GetYaxis()->SetTitle("#sigma (pb)");

	CorrDiffHisReal->GetYaxis()->SetTitleOffset(1.3);
	CorrDiffHisReal->GetXaxis()->CenterTitle();
	CorrDiffHisReal->GetYaxis()->CenterTitle();




	for(int i = 0; i < NBins;i++){
		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);
		CorrYieldDiff[i] = (RawCount /  NewEffReal[i])/(BRchain*2* lumi);
		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr /  NewEffReal[i]) *(RawCountErr  /  NewEffReal[i]) + (RawCount /NewEffReal[i] *  NewEffRealErr[i]) * (RawCount/NewEffReal[i]  *  NewEffRealErr[i]))/(BRchain*2* lumi);
		CorrDiffHisReal->SetBinContent(i+1,CorrYieldDiff[i]);
		CorrDiffHisReal->SetBinError(i+1,CorrYieldDiffErr[i]);

	}

	
	//CorrDiffHisReal->SetTitle("(Preliminary) B^{+} #rightarrow J/#psi K^{+} p_{T} Differential Cross Section in pp");

	CorrDiffHisReal->SetMarkerColor(kBlack);
	CorrDiffHisReal->SetMarkerSize(1);
	CorrDiffHisReal->SetMarkerStyle(20);

	TFile * foutCorr;
	if(DoTnP == 0 && usemc==0)	foutCorr = new TFile(Form("FinalFiles/BPPPCorrYield%sNoTnP.root",var_m.Data()),"RECREATE");
	if(DoTnP == 1 && usemc==0)	foutCorr = new  TFile(Form("FinalFiles/BPPPCorrYield%s.root",var_m.Data()),"RECREATE");
	if(DoTnP == 0 && usemc==1)	foutCorr = new TFile(Form("FinalFiles/BPPPCorrYield%sNoTnPMC.root",var_m.Data()),"RECREATE");
	if(DoTnP == 1 && usemc==1)	foutCorr = new  TFile(Form("FinalFiles/BPPPCorrYield%sMC.root",var_m.Data()),"RECREATE");

	TH1D * Eff1DHisvar;
	if (whichvar==0){
		Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisY");
	    }
	if (whichvar==1){
		Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHisMult");
		}
	if (whichvar==2){
		Eff1DHisvar=(TH1D *) fin1DEff->Get("Eff1DHis");
	}



	TH1D * CorrDiffHisBin = new TH1D("CorrDiffHisBin","",NBins,ptBins);
	CorrDiffHisBin->GetXaxis()->SetTitle("nMult");
	CorrDiffHisBin->GetYaxis()->SetTitle("#sigma (pb)");

	CorrDiffHisBin->GetYaxis()->SetTitleOffset(1.3);
	CorrDiffHisBin->GetXaxis()->CenterTitle();
	CorrDiffHisBin->GetYaxis()->CenterTitle();

	CorrDiffHisBin->SetMarkerColor(kRed);
	CorrDiffHisBin->SetLineColor(kRed);	
	CorrDiffHisBin->SetMarkerSize(1);
	CorrDiffHisBin->SetMarkerStyle(20);

	float Eff1D[NBins];
	float Eff1DErr[NBins];

	for(int i = 0; i < NBins;i++){
		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);
		Eff1D[i] = Eff1DHisvar->GetBinContent(i+1);
		Eff1DErr[i] = Eff1DHisvar->GetBinError(i+1);

		//		CorrYieldDiff[i] = (RawCount *  Eff1D[i])/(BRchain*2* lumi);
		//		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr *  Eff1D[i]) *(RawCountErr  *  Eff1D[i]) + (RawCount *  Eff1DErr[i]) * (RawCount  *  Eff1DErr[i]))/(BRchain*2* lumi);
		CorrYieldDiff[i] = (RawCount /  Eff1D[i])/(BRchain*2* lumi);
		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr /  Eff1D[i]) *(RawCountErr  /  Eff1D[i]) + (RawCount /Eff1D[i] *  Eff1DErr[i]) * (RawCount /Eff1D[i] *  Eff1DErr[i]))/(BRchain*2* lumi);

		CorrDiffHisBin->SetBinContent(i+1,CorrYieldDiff[i]);
		CorrDiffHisBin->SetBinError(i+1,CorrYieldDiffErr[i]);

	}
	
	hInvEff->SetMaximum(NewEff[0]*1.5);
	TCanvas *c = new TCanvas("c","c",600,600);
	c->cd();

	hEffInv->Draw("ep");
	
  	//c->BuildLegend(0.6, 0.6, 0.9, 0.8);
  	
	gSystem->mkdir("EffFinal",true);
	//c->SaveAs(Form("EffFinal/ReAnaEff_%dBins.png",NBins));
	if (usemc==0){c->SaveAs(Form("EffFinal/ReAnaEff_%dBins_%s.pdf",NBins,var_m.Data()));}
	else {c->SaveAs(Form("EffFinal/ReAnaEff_%dBins_%s_MC.pdf",NBins,var_m.Data()));}
	hEff->Draw("ep");

	//c->SaveAs(Form("EffFinal/ReAnaEffReal_%dBins.png",NBins));
	if (usemc==0){c->SaveAs(Form("EffFinal/ReAnaEffReal_%dBins_%s.pdf",NBins,var_m.Data()));}
	else {c->SaveAs(Form("EffFinal/ReAnaEffReal_%dBins_%s_MC.pdf",NBins,var_m.Data()));}
	Eff1DHisvar->Draw("ep");

	if (usemc==0){c->SaveAs(Form("EffFinal/ReAnaEff1D_%dBins_%s.pdf",NBins,var_m.Data()));}
	else {c->SaveAs(Form("EffFinal/ReAnaEff1D_%dBins_%s_MC.pdf",NBins,var_m.Data()));}

	TH2D * invAcc2D=(TH2D *) fin1DEff->Get("invAcc2D");
	invAcc2D->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	invAcc2D->GetYaxis()->SetTitle("rapidity");

	invAcc2D->Draw("pcol");

	if (usemc==0){c->SaveAs("EffFinal/acc_2Dmap.pdf");}
	else {c->SaveAs("EffFinal/acc_2Dmap.pdf");}

	invEff2D->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	invEff2D->GetYaxis()->SetTitle("rapidity");

	invEff2D->Draw("pcol");

	if (usemc==0){c->SaveAs("EffFinal/totaleff_2Dmap.pdf");}
	else {c->SaveAs("EffFinal/totaleff_2Dmap.pdf");}

	hEffInv->SetMarkerColor(kRed+1);
	hEffInv->SetLineColor(kRed+1);
	hEff->SetMarkerColor(kBlue+1);
	hEff->SetLineColor(kBlue+1);
	hEffInv->Draw("ep");
	hEff->Draw("ep");
	Eff1DHisvar->Draw("ep");

	TLegend *leg = new TLegend(0.30,0.65,0.75,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);

	leg->AddEntry(hEffInv,"2D Inv.","pl");
	leg->AddEntry(hEff,"2D","pl");
	leg->AddEntry(Eff1DHisvar,"1D","pl");

	leg->Draw("SAME");

	//if (usemc==0){c->SaveAs(Form("EffFinal/ReAnaEffcomp_%dBins_%s.pdf",NBins,var_m.Data()));}
	//else {c->SaveAs(Form("EffFinal/ReAnaEffcomp_%dBins_%s_MC.pdf",NBins,var_m.Data()));}




	foutCorr->cd();
	CorrDiffHis->Write();
	CorrDiffHisBin->Write();
	CorrDiffHisReal->Write();
	hInvEff->Write();
	hEff->Write();
	hEffInv->Write();
	Eff1DHisvar->Write();
	hPt->Write();
	foutCorr->Close();

}
