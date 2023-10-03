#include <TFile.h>        // To work with ROOT files
#include <TTree.h>        // To work with TTree objects
#include <TH2D.h>         // To create and work with 2D histograms
#include <TCanvas.h>      // To create and work with canvas
#include <TROOT.h>        // To work with ROOT environment settings
#include <iostream>       // For standard input/output operations
   
   
   
#define BS_MASS 5.36682  //"Bmass > 5.35182 && Bmass < 5.38182"
#define BP_MASS 5.27915


void DataView() {
    TFile* file = TFile::Open("~/Downloads/BPMC_nom_sm.root"); 
    TTree* tree = (TTree*)file->Get("ntKp");


    TCanvas *canvas3D = new TCanvas("canvas3D", "3D Histogram", 800, 600);
    canvas3D->SetTheta(30);
    canvas3D->SetPhi(30);
    gStyle->SetOptStat(0);

    tree->Draw("abs(By):Bpt >>+ hist2D(50,5,60,50,1.4,2.1)", "Bmass > 5.26915 && Bmass < 5.28915", "colz");

    // Retrieve the 2D histogram from the current pad
    TH2D *hist2D = (TH2D*)gPad->GetPrimitive("hist2D");
    TLine *line = new TLine(hist2D->GetXaxis()->GetXmin(), 1.5, hist2D->GetXaxis()->GetXmax(), 1.5);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    TLine *line1 = new TLine(hist2D->GetXaxis()->GetXmin(), 2, hist2D->GetXaxis()->GetXmax(), 2);
    line1->SetLineColor(kRed);
    line1->SetLineWidth(2);
    line1->Draw();
    line->Draw();

    // Keep the canvas on the screen
    //canvas3D->Draw();
}


