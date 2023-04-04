#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "../henri2022/parameter.h" 

using namespace std;


// 0 for B+ 
// 1 for B_s
void CalEffSystB(int meson_n ){

	TString B_m ;
	int NBins = 7;
	TString t_tree ;
	int NCand = 10 ;
	double b_m_mass ;

	if(meson_n == 0){
		NBins = nptBinsBP;
		B_m = "BP";
		t_tree = "ntKp";
		NCand = 10;
		b_m_mass = 5.27932 ;
	} else {
		NBins = nptBins;
		B_m = "Bs";
		t_tree = "ntphi";
		NCand = 15;
		b_m_mass = 5.3663 ;
	}

	//MOMENTUM BINS
	double ptBins[NBins + 1];
	if(meson_n==0) { for( int c=0; c <NBins+1; c++){ ptBins[c]=ptbinsvecBP[c];}}
	else{            for( int c=0; c <NBins+1; c++){ ptBins[c]=ptbinsvec[c];  }}

//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	TString FileName = Form("/data3/tasheng/presel/%sData_nom.root",B_m.Data());
	TFile * fin = new TFile(FileName.Data());
	fin->cd();
	TTree * EffInfoTree = (TTree * ) fin->Get(t_tree.Data());
	int NEvents = EffInfoTree->GetEntries();

	Int_t BsizeNew;
	Float_t BmassNew[NCand];
	Float_t ByNew[NCand];
	Float_t BptNew[NCand];
	EffInfoTree->SetBranchAddress("Bsize", &BsizeNew);
	EffInfoTree->SetBranchAddress("Bmass", BmassNew);
	EffInfoTree->SetBranchAddress("By", ByNew);
	EffInfoTree->SetBranchAddress("Bpt", BptNew);
	Float_t BEffInv[NCand];
	Float_t BEffInvErr[NCand];
	Float_t BEffInvTnPUp[NCand];
	Float_t BEffInvTnPDown[NCand];
	Float_t BEffInvBDT[NCand];
	Float_t BEffInvBpt[NCand];
	
	int Counts[NBins];
	double SumCounts[NBins];
	double SumCountsErr[NBins];
	double NewEff[NBins];
	double NewEffErr[NBins];

	//Syst Collection
	double SumCountsTnPUpSyst[NBins];
	double SumCountsTnPDownSyst[NBins];
	double SumCountsBDTSyst[NBins];
	double SumCountsBptSyst[NBins];

	// Initialize in zero
	for(int i = 0; i < NBins; i++){
		Counts[i] = 0;
		SumCounts[i] = 0;
		SumCountsErr[i] = 0;

		SumCountsTnPUpSyst[i] = 0;
		SumCountsTnPDownSyst[i] = 0;
		SumCountsBDTSyst[i] = 0;
		SumCountsBptSyst[i] = 0;
	}

	TFile * finSyst2D = new TFile(Form("../%s/EffAna/NewEff2DMaps/%sSyst2D.root",B_m.Data(), B_m.Data()));
	TH2D * invEff2D = (TH2D *) finSyst2D->Get("invEff2D");
	TH2D * invEff2DTnPSystUp = (TH2D *) finSyst2D->Get("invEff2DTnPSystUp");
	TH2D * invEff2DTnPSystDown = (TH2D *) finSyst2D->Get("invEff2DTnPSystDown");
	TH2D * invEff2DBDTSyst = (TH2D *) finSyst2D->Get("invEff2DBDTSyst");
	TH2D * invEff2DBptSyst = (TH2D *) finSyst2D->Get("invEff2DBptSyst");

	int XBin;
	int YBin;
	for( int i = 0; i < NEvents; i++){
		EffInfoTree->GetEntry(i);

		for(int j = 0; j < BsizeNew; j++){
			for(int k = 0; k < NBins; k++){
				
				if(BptNew[j] > ptBins[k] && BptNew[j] < ptBins[k+1] && TMath::Abs(BmassNew[j] - b_m_mass) < 0.08 && TMath::Abs(ByNew[j]) < 2.4 && ((BptNew[j] > 5 && BptNew[j] < 10 && abs(ByNew[j]) > 1.5)||(BptNew[j] > 10))){
					XBin = invEff2D->GetXaxis()->FindBin( BptNew[j]);
					YBin = invEff2D->GetYaxis()->FindBin( TMath::Abs(ByNew[j]));
					
					BEffInv[j] = invEff2D->GetBinContent(XBin,YBin);
					BEffInvErr[j] = invEff2D->GetBinError(XBin,YBin);
					BEffInvTnPUp[j] = invEff2DTnPSystUp->GetBinContent(XBin,YBin);
					BEffInvTnPDown[j] = invEff2DTnPSystDown->GetBinContent(XBin,YBin);
					BEffInvBDT[j] = invEff2DBDTSyst->GetBinContent(XBin,YBin);
					BEffInvBpt[j] = invEff2DBptSyst->GetBinContent(XBin,YBin);

					if(BEffInv[j] > 0){
						SumCounts[k] = SumCounts[k] + BEffInv[j];
						SumCountsErr[k] = SumCountsErr[k] + BEffInvErr[j] * BEffInvErr[j];
						SumCountsTnPUpSyst[k] = BEffInvTnPUp[j] + SumCountsTnPUpSyst[k];
						SumCountsTnPDownSyst[k] = BEffInvTnPDown[j] + SumCountsTnPDownSyst[k];
						SumCountsBDTSyst[k] = BEffInvBDT[j] + SumCountsBDTSyst[k];
						SumCountsBptSyst[k] = BEffInvBpt[j] + SumCountsBptSyst[k];
						Counts[k] = Counts[k] + 1;
					}
				}
			}
		}
	}

	double EffTnPUp[NBins];
	double EffTnPDown[NBins];
	double EffBDT[NBins];
	double EffBpt[NBins];
	for(int i = 0; i < NBins; i++){
		NewEff[i] = SumCounts[i]/Counts[i];
		NewEffErr[i] = TMath::Sqrt(SumCountsErr[i])/Counts[i];
		EffTnPUp[i] = SumCountsTnPUpSyst[i]/Counts[i];
		EffTnPDown[i] = SumCountsTnPDownSyst[i]/Counts[i];
		EffBDT[i] = SumCountsBDTSyst[i]/Counts[i];
		EffBpt[i] = SumCountsBptSyst[i]/Counts[i];

	cout << "----------------------------------------------------------------------------------------" << endl;
	cout << "Count =  " <<  Counts[i] << "   NewEff = " << NewEff[i] << "     NewEffErr = " << NewEffErr[i] << "  Fractional = " << NewEffErr[i]/NewEff[i] << endl;

	}

	TH1D * Eff2DHis = new TH1D("Eff2DHis","",NBins,ptBins);
	Eff2DHis->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	Eff2DHis->GetYaxis()->SetTitle("<1/(#alpha#epsilon)>");
	Eff2DHis->GetXaxis()->CenterTitle();	
	Eff2DHis->GetYaxis()->CenterTitle();
	Eff2DHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff2DTnPUpSystHis = new TH1D("Eff2DTnPUpSystHis","",NBins,ptBins);
	Eff2DTnPUpSystHis->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	Eff2DTnPUpSystHis->GetYaxis()->SetTitle("<1/(#alpha#epsilon)>");
	Eff2DTnPUpSystHis->GetXaxis()->CenterTitle();	
	Eff2DTnPUpSystHis->GetYaxis()->CenterTitle();
	Eff2DTnPUpSystHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DTnPUpSystHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff2DTnPDownSystHis = new TH1D("Eff2DTnPDownSystHis","",NBins,ptBins);
	Eff2DTnPDownSystHis->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	Eff2DTnPDownSystHis->GetYaxis()->SetTitle("<1/(#alpha#epsilon)>");
	Eff2DTnPDownSystHis->GetXaxis()->CenterTitle();	
	Eff2DTnPDownSystHis->GetYaxis()->CenterTitle();
	Eff2DTnPDownSystHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DTnPDownSystHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff2DBDTHis = new TH1D("Eff2DBDTHis","",NBins,ptBins);
	Eff2DBDTHis->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	Eff2DBDTHis->GetYaxis()->SetTitle("<1/(#alpha#epsilon)>");
	Eff2DBDTHis->GetXaxis()->CenterTitle();	
	Eff2DBDTHis->GetYaxis()->CenterTitle();
	Eff2DBDTHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DBDTHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff2DBptHis = new TH1D("Eff2DBptHis","",NBins,ptBins);
	Eff2DBptHis->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	Eff2DBptHis->GetYaxis()->SetTitle("<1/(#alpha#epsilon)>");
	Eff2DBptHis->GetXaxis()->CenterTitle();	
	Eff2DBptHis->GetYaxis()->CenterTitle();
	Eff2DBptHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DBptHis->GetYaxis()->SetTitleOffset(1.5);

	for(int i = 0; i < NBins; i++){
		Eff2DHis->SetBinContent(i+1, NewEff[i]);
		Eff2DHis->SetBinError(i+1, NewEffErr[i]);
		Eff2DTnPUpSystHis->SetBinContent(i+1, EffTnPUp[i]);		
		Eff2DTnPUpSystHis->SetBinError(i+1, NewEffErr[i]);
		Eff2DTnPDownSystHis->SetBinContent(i+1, EffTnPDown[i]);		
		Eff2DTnPDownSystHis->SetBinError(i+1, NewEffErr[i]);
		Eff2DBDTHis->SetBinContent(i+1, EffBDT[i]);		
		Eff2DBDTHis->SetBinError(i+1, NewEffErr[i]);
		Eff2DBptHis->SetBinContent(i+1, EffBpt[i]);		
		Eff2DBptHis->SetBinError(i+1, NewEffErr[i]);
	}
	
	gSystem->mkdir("OutFiles", true);
	TFile * fout = new TFile(Form("OutFiles/%sSyst2D.root",B_m.Data()),"RECREATE");
	fout->cd();
	Eff2DHis->Write();
	Eff2DTnPUpSystHis->Write();
	Eff2DTnPDownSystHis->Write();
	Eff2DBDTHis->Write();
	Eff2DBptHis->Write();
	fout->Close();
}
