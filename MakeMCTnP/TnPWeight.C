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
#include "tnp_weight.h"

using namespace std;

using std::cout;
using std::endl;


void TnPWeight(int Opt){


	const int NCand = 8000;

	TString infile;
	TString TreeName;
	TString outfile;

	if(Opt == 0){
	//	infile = "../UnskimmedSamples/OfficialMC/BPMC.root";
		// infile = "../../bmva/TMVA/BP/sample/BPMC_5_60.root";
		// infile = "../../dat/BP_MC_all.root";
		infile = "/data3/smcosta/data/BPMC_nom_BDT.root";
		TreeName = "ntKp";
		outfile = "BPTnPInfo.root";

	}
	if(Opt == 1){

		// infile = "../UnskimmedSamples/OfficialMC/BsMC.root";
		// infile = "../../bmva/TMVA/Bs/sample/BsMC_7_50.root";
		// infile = "../../dat/Bs_MC_all.root";
		infile = "/data3/smcosta/data/BsMC_nom_BDT.root";
		TreeName = "ntphi";
		//	outfile = "BsTnPInfo.root";
		outfile = "BsTnPInfo.root";
	}

	Double_t Bmu1pt;
	Double_t Bmu2pt;
	Double_t Bmu1eta;
	Double_t Bmu2eta;
	Long64_t Bsize;
	Int_t BsizeNew;

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	cout << "Pass 1" << endl;
	cout << "TreeName = " << TreeName.Data() << endl;

	TTree * t = (TTree * ) fin->Get(TreeName.Data());


	t->SetBranchAddress("Bmu1pt",&Bmu1pt);
	t->SetBranchAddress("Bmu2pt",&Bmu2pt);
	t->SetBranchAddress("Bmu1eta",&Bmu1eta);
	t->SetBranchAddress("Bmu2eta",&Bmu2eta);
	t->SetBranchAddress("Bsize",&Bsize);


	Int_t BsizeTnP;
	Float_t TnPMu1Nominal;
	Float_t TnPMu1StatError;
	Float_t TnPMu1SystError;
	Float_t TnPMu1Error;
	Float_t TnPMu2Nominal;
	Float_t TnPMu2StatError;
	Float_t TnPMu2SystError;
	Float_t TnPMu2Error;
	Float_t TnPNominal;
	Float_t TnPStatError;
	Float_t TnPSystError;
	Float_t TnPError;


	TFile * fout = new TFile(outfile.Data(),"RECREATE");
	fout->cd();

	TTree * TnPInfo = new TTree("TnPInfo","TnPInfo");
	//TnPInfo->Branch("BsizeTnP",BsizeTnP,"BsizeTnP/F");

	TnPInfo->Branch("BsizeNew",&BsizeNew,"BsizeNew/I");
	TnPInfo->Branch("TnPMu1Nominal",&TnPMu1Nominal,"TnPMu1Nominal/F");
	TnPInfo->Branch("TnPMu1StatError",&TnPMu1StatError,"TnPMu1StatError/F");
	TnPInfo->Branch("TnPMu1SystError",&TnPMu1SystError,"TnPMu1SystError/F");
	TnPInfo->Branch("TnPMu1Error",&TnPMu1Error,"TnPMu1Error/F");
	TnPInfo->Branch("TnPMu2Nominal",&TnPMu2Nominal,"TnPMu2Nominal/F");
	TnPInfo->Branch("TnPMu2StatError",&TnPMu2StatError,"TnPMu2StatError/F");
	TnPInfo->Branch("TnPMu2SystError",&TnPMu2SystError,"TnPMu2SystError/F");
	TnPInfo->Branch("TnPMu2Error",&TnPMu2Error,"TnPMu2Error/F");
	TnPInfo->Branch("TnPNominal",&TnPNominal,"TnPNominal/F");
	TnPInfo->Branch("TnPStatError",&TnPStatError,"TnPStatError/F");
	TnPInfo->Branch("TnPSystError",&TnPSystError,"TnPSystError/F");
	TnPInfo->Branch("TnPError",&TnPError,"TnPError/F");

	int NEvent = t->GetEntries();

	//NEvent = 10;

	for(int i = 0; i < NEvent; i++){

		if(i%100000 == 0) cout << "Event = " << i << endl;

		t->GetEntry(i);
		//MuonInfoTree->GetEntry(i);
		//		BsizeTnP = Bsize;
//		cout  << "i = " << i << "  Bmu1pt =" << Bmu1pt << "  Bmu1eta = " << Bmu1eta << endl;
		BsizeNew = Bsize;

		auto TnPSF1 = tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(Bmu1pt,Bmu1eta);
		auto TnPSF2 = tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(Bmu2pt,Bmu2eta);
		
		//cout << "Bmu1pt[j] = " << Bmu1pt[j] << " Bmu1eta[j] =  " << Bmu1eta[j] << endl;
		//cout << "Bmu2pt[j] = " << Bmu2pt[j] << " Bmu2eta[j] =  " << Bmu2eta[j] << endl;

		TnPMu1Nominal =  std::get<0>(TnPSF1);
		TnPMu1StatError = std::get<1>(TnPSF1);
		TnPMu1SystError = std::get<2>(TnPSF1);
		TnPMu1Error = std::get<3>(TnPSF1);


		TnPMu2Nominal =  std::get<0>(TnPSF2);
		TnPMu2StatError = std::get<1>(TnPSF2);
		TnPMu2SystError = std::get<2>(TnPSF2);
		TnPMu2Error = std::get<3>(TnPSF2);

//		cout << "TnPNominal = " << TnPNominal << "   TnPMu1Nominal = " << TnPMu1Nominal <<  "   TnPMu2Nominal = " << TnPMu2Nominal  << endl;

		TnPNominal =  TnPMu1Nominal * TnPMu2Nominal;
		TnPStatError = TnPNominal * TMath::Sqrt(TnPMu1StatError/TnPMu1Nominal * TnPMu1StatError/TnPMu1Nominal + TnPMu2StatError/TnPMu2Nominal * TnPMu2StatError/TnPMu2Nominal);
		TnPSystError = TnPNominal * TMath::Sqrt(TnPMu1SystError/TnPMu1Nominal * TnPMu1SystError/TnPMu1Nominal + TnPMu2SystError/TnPMu2Nominal * TnPMu2SystError/TnPMu2Nominal);
		TnPError = TnPNominal * TMath::Sqrt(TnPMu1Error/TnPMu1Nominal * TnPMu1Error/TnPMu1Nominal + TnPMu2Error/TnPMu2Nominal * TnPMu2Error/TnPMu2Nominal);



		TnPInfo->Fill();	


	}

	fout->Write();
	fout->Close();

	fin->Close();
}


