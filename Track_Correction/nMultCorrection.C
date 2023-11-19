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
#include "trackingEfficiency2017pp.h"
using namespace std;

using std::cout;
using std::endl;

void nMultCorrection(int meson_n){

TString var;
if (meson_n==0){var="BP";}
else {var="Bs";}

TFile* of = new TFile(Form("/data3/tasheng/presel/%sData_nom.root",var.Data()));

TTree* ot;

if (meson_n==0){ot = (TTree*) of->Get("ntKp");}
else {ot = (TTree*) of->Get("ntphi");}

int NCand=15;

Int_t nMult;
float Btrk1Pt[NCand];
float Btrk1Eta[NCand];

ot->SetBranchAddress("Btrk1Pt",&Btrk1Pt);
ot->SetBranchAddress("nMult",&nMult); 
ot->SetBranchAddress("Btrk1Eta", &Btrk1Eta);

ot->SetBranchStatus("*",1);

TFile* nf = new TFile(Form("/data3/smcosta/data/%sData_nom_Multcorr.root",var.Data()),"recreate");
TTree* nt = ot->CloneTree(0);


TrkEff2017pp trkEff =  TrkEff2017pp(false, "");

for (auto i=0; i<ot->GetEntries(); i++){

    ot->GetEntry(i);
    nMult = nMult * trkEff.getCorrection(Btrk1Pt[0],Btrk1Eta[0]);
    nt->Fill();
}
nt->AutoSave();

}