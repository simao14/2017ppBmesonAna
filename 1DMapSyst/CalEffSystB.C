#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "../henri2022/parameter.h" 

using namespace std;


// 0 for B+ 
// 1 for B_s
void CalEffSystB(int meson_n, int whichvar){

	TString B_m ;
	int NBins;
	double b_m_mass ;
	TString var_n;
	TString var_l;

	if(meson_n == 0){
		B_m = "BP";
	} 
	if(meson_n == 1){
		B_m = "Bs";
	}


	if(meson_n == 0 && whichvar==0){
		NBins = nptBinsBP;
		var_n="pt";
		var_l="p_{T} [GeV/c]";
		
	}
	if(meson_n == 1 && whichvar==0){
		NBins = nptBins;
		var_n="pt";
		var_l="p_{T} [GeV/c]";
	}
	if(whichvar==1){
		NBins = nyBins_both;
		var_n="y";
		var_l="Rapidity";
		
	}
	if(whichvar==2){
		NBins = nmBins_both;
		var_n="Mult";
		var_l="Multiplicity";
		
	}

	// BINS
	double ptBins[NBins + 1];
	if (whichvar==1){
	for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ybinsvec[i];             
		}
	} 
	if (whichvar==2){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  nmbinsvec[i];             
		}
	}
	if (whichvar==0 && meson_n==0){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvecBP[i];             
		}
	}
	if (whichvar==0 && meson_n!=0){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvec[i];             
		}
	}

//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	double NewEff[NBins];
	double NewEffErr[NBins];

	TFile * finSyst1D = new TFile(Form("../EffAna/%s/NewEff2DMaps/%sSyst.root",B_m.Data(), B_m.Data()));
	TH1D * Eff1D;
	TH1D * Eff1DTnPSystUp;
	TH1D * Eff1DTnPSystDown;
	TH1D * Eff1DBDTSyst;
	TH1D * Eff1DBptSyst;

	if (whichvar==0) {
		Eff1D = (TH1D *) finSyst1D->Get("Eff1DHis");
		Eff1DTnPSystUp = (TH1D *) finSyst1D->Get("Eff1DHisTnPUp");
		Eff1DTnPSystDown = (TH1D *) finSyst1D->Get("Eff1DHisTnPDown");
		Eff1DBDTSyst = (TH1D *) finSyst1D->Get("Eff1DHisBDT");
		Eff1DBptSyst = (TH1D *) finSyst1D->Get("Eff1DHisBpt");
	} else if (whichvar==1) {
		Eff1D = (TH1D *) finSyst1D->Get("Eff1DHisY");
		Eff1DTnPSystUp = (TH1D *) finSyst1D->Get("Eff1DHisTnPUpY");
		Eff1DTnPSystDown = (TH1D *) finSyst1D->Get("Eff1DHisTnPDownY");
		Eff1DBDTSyst = (TH1D *) finSyst1D->Get("Eff1DHisBDTY");
		Eff1DBptSyst = (TH1D *) finSyst1D->Get("Eff1DHisBptY");
	}	else if (whichvar==2) {
		Eff1D = (TH1D *) finSyst1D->Get("Eff1DHisMult");
		Eff1DTnPSystUp = (TH1D *) finSyst1D->Get("Eff1DHisTnPUpMult");
		Eff1DTnPSystDown = (TH1D *) finSyst1D->Get("Eff1DHisTnPDownMult");
		Eff1DBDTSyst = (TH1D *) finSyst1D->Get("Eff1DHisBDTMult");
		Eff1DBptSyst = (TH1D *) finSyst1D->Get("Eff1DHisBptMult");
	}

	double EffTnPUp[NBins];
	double EffTnPDown[NBins];
	double EffBDT[NBins];
	double EffBpt[NBins];
	for(int i = 0; i < NBins; i++){
		NewEff[i] = Eff1D->GetBinContent(i+1);
		NewEffErr[i] = Eff1D->GetBinError(i+1);
		EffTnPUp[i] = Eff1DTnPSystUp->GetBinContent(i+1);
		EffTnPDown[i] = Eff1DTnPSystDown->GetBinContent(i+1);
		EffBDT[i] = Eff1DBDTSyst->GetBinContent(i+1);
		EffBpt[i] = Eff1DBptSyst->GetBinContent(i+1);
	}

	TH1D * Eff1DHisInv = new TH1D("Eff1DHisInv","",NBins,ptBins);
	Eff1DHisInv->GetXaxis()->SetTitle(var_l.Data());
	Eff1DHisInv->GetYaxis()->SetTitle("1/(#alpha #times #epsilon)");
	Eff1DHisInv->GetXaxis()->CenterTitle();	
	Eff1DHisInv->GetYaxis()->CenterTitle();
	Eff1DHisInv->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DHisInv->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff1DTnPUpSystHis = new TH1D("Eff1DTnPUpSystHis","",NBins,ptBins);
	Eff1DTnPUpSystHis->GetXaxis()->SetTitle(var_l.Data());
	Eff1DTnPUpSystHis->GetYaxis()->SetTitle("1/(#alpha #times #epsilon)");
	Eff1DTnPUpSystHis->GetXaxis()->CenterTitle();	
	Eff1DTnPUpSystHis->GetYaxis()->CenterTitle();
	Eff1DTnPUpSystHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DTnPUpSystHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff1DTnPDownSystHis = new TH1D("Eff1DTnPDownSystHis","",NBins,ptBins);
	Eff1DTnPDownSystHis->GetXaxis()->SetTitle(var_l.Data());
	Eff1DTnPDownSystHis->GetYaxis()->SetTitle("1/(#alpha #times #epsilon)");
	Eff1DTnPDownSystHis->GetXaxis()->CenterTitle();	
	Eff1DTnPDownSystHis->GetYaxis()->CenterTitle();
	Eff1DTnPDownSystHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DTnPDownSystHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff1DBDTHis = new TH1D("Eff1DBDTHis","",NBins,ptBins);
	Eff1DBDTHis->GetXaxis()->SetTitle(var_l.Data());
	Eff1DBDTHis->GetYaxis()->SetTitle("1/(#alpha #times #epsilon)");
	Eff1DBDTHis->GetXaxis()->CenterTitle();	
	Eff1DBDTHis->GetYaxis()->CenterTitle();
	Eff1DBDTHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DBDTHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff1DBptHis = new TH1D("Eff1DBptHis","",NBins,ptBins);
	Eff1DBptHis->GetXaxis()->SetTitle(var_l.Data());
	Eff1DBptHis->GetYaxis()->SetTitle("1/(#alpha #times #epsilon)");
	Eff1DBptHis->GetXaxis()->CenterTitle();	
	Eff1DBptHis->GetYaxis()->CenterTitle();
	Eff1DBptHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DBptHis->GetYaxis()->SetTitleOffset(1.5);

	for(int i = 0; i < NBins; i++){
		Eff1DHisInv->SetBinContent(i+1, 1/NewEff[i]);
		Eff1DHisInv->SetBinError(i+1, NewEffErr[i]/(NewEff[i] * NewEff[i]));
		Eff1DTnPUpSystHis->SetBinContent(i+1, 1/EffTnPUp[i]);		
		Eff1DTnPUpSystHis->SetBinError(i+1, NewEffErr[i]/(EffTnPUp[i]*EffTnPUp[i]));
		Eff1DTnPDownSystHis->SetBinContent(i+1, 1/EffTnPDown[i]);		
		Eff1DTnPDownSystHis->SetBinError(i+1, NewEffErr[i]/(EffTnPDown[i]*EffTnPDown[i]));
		Eff1DBDTHis->SetBinContent(i+1, 1/EffBDT[i]);		
		Eff1DBDTHis->SetBinError(i+1, NewEffErr[i]/(EffBDT[i]*EffBDT[i]));
		Eff1DBptHis->SetBinContent(i+1, 1/EffBpt[i]);		
		Eff1DBptHis->SetBinError(i+1, NewEffErr[i]/(EffBpt[i]*EffBpt[i]));
	}
	
	gSystem->mkdir("OutFiles", true);
	TFile * fout = new TFile(Form("OutFiles/%sSyst1D_%s.root",B_m.Data(),var_n.Data()),"RECREATE");
	fout->cd();
	Eff1DHisInv->Write();
	Eff1DTnPUpSystHis->Write();
	Eff1DTnPDownSystHis->Write();
	Eff1DBDTHis->Write();
	Eff1DBptHis->Write();
	fout->Close();
}
