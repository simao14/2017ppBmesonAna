#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TSystem.h"
#include <fstream>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include "../../../henri2022/parameter.h" 
#include "CMS_lumi.C" 

using namespace std;
using std::cout;
using std::endl;

void divideTGraphsInFiles(TString inputFile1, TString inputFile2, TString whichvar);



void Bmeson_Ratio(){

    //First Retrive the BP Xsection with matching bins

    //First Retrive the BP Xsection with matching bins

// fs/fu 
// pT Fragmentation Fraction

TString BP_FILE_pt = "./ROOTFiles/BP_Xsection_pt.root";
TString Bs_FILE_pt = "./ROOTFiles/Bs_Xsection_pt.root";
divideTGraphsInFiles(Bs_FILE_pt, BP_FILE_pt, "p_{T}") ;

// pT Fragmentation Fraction
// y Fragmentation Fraction

TString BP_FILE_y = "./ROOTFiles/BP_Xsection_y.root";
TString Bs_FILE_y = "./ROOTFiles/Bs_Xsection_y.root";
divideTGraphsInFiles(Bs_FILE_y, BP_FILE_y, "|y|") ;

// y Fragmentation Fraction
// fs/fu 


}












void divideTGraphsInFiles(TString inputFile1, TString inputFile2, TString whichvar) {
    //inputFile1 SHOULD BE Bmeson NUMERATOR
    //inputFile2 SHOULD BE Bmeson DENOMINATOR      
    TString Name1= "";
    TString Name2= "";

    if (whichvar == "p_{T}"){
        Name2= "BP_Xsection_pt";
        Name1= "Bs_Xsection_pt";
    } else {
        Name2= "BP_Xsection_y";
        Name1= "Bs_Xsection_y";
    }

    TFile* file1 = TFile::Open(inputFile1, "READ");
	TMultiGraph *mg1= (TMultiGraph*) file1->Get(Name1);
    TGraphAsymmErrors* Y_stat_B1 = dynamic_cast<TGraphAsymmErrors*>(mg1->GetListOfGraphs()->FindObject("Y_stat"));
    TGraphAsymmErrors* Y_syst_B1 = dynamic_cast<TGraphAsymmErrors*>(mg1->GetListOfGraphs()->FindObject("Y_syst"));

    TFile* file2 = TFile::Open(inputFile2, "READ");
	TMultiGraph *mg2= (TMultiGraph*) file2->Get(Name2);
    TGraphAsymmErrors* Y_stat_B2 = dynamic_cast<TGraphAsymmErrors*>(mg2->GetListOfGraphs()->FindObject("Y_stat"));
    TGraphAsymmErrors* Y_syst_B2 = dynamic_cast<TGraphAsymmErrors*>(mg2->GetListOfGraphs()->FindObject("Y_syst"));

    if (!Y_stat_B1 || !Y_syst_B1 || !Y_stat_B2 || !Y_syst_B2) {
        std::cout << "Error: One or more TGraphAsymmErrors not found in the multi-graphs." << std::endl;
        return;
    }

    // Retrieve values
    // FOR NOW THERE IS NO BP with matching number of bins so i=0 i=1 i=2 (first 3 bins) are being considered only
    double X_POS[3] ;
    double Frag_f[3];
    double Frag_f_Stat_U[3];
    double Frag_f_Syst_U[3];

    Double_t B_X1, B_Y1, B_X2, B_Y2, Y1_stat, Y1_syst, Y2_stat, Y2_syst ;

    int dummy = 0;
    if (whichvar == "|y|")
    {
        dummy = 1;
    }

    for (int i = 0; i < 3; ++i) {

        Y_stat_B1->GetPoint(i, B_X1, B_Y1);           // X and Yield values of NUMERATOR Bmeson
        Y1_stat = Y_stat_B1->GetErrorYhigh(i);   // Yield Statistical Unc. of NUMERATOR Bmeson
        Y1_syst = Y_syst_B1->GetErrorYhigh(i);   // Yield Systematic Unc. of NUMERATOR Bmeson
       
        Y_stat_B2->GetPoint(i+1-dummy, B_X2, B_Y2);           // X and Yield values of DENOMINATOR Bmeson
        Y2_stat = Y_stat_B2->GetErrorYhigh(i+1-dummy);   // Yield Statistical Unc. of DENOMINATOR Bmeson
        Y2_syst = Y_syst_B2->GetErrorYhigh(i+1-dummy);   // Yield Systematic Unc. of DENOMINATOR Bmeson

        X_POS[i] = (B_X1+B_X2)/2;
        Frag_f[i] = B_Y1 / B_Y2;
        Frag_f_Stat_U[i] = fabs(Frag_f[i]) * sqrt(pow(Y1_stat / B_Y1, 2) + pow(Y2_stat / B_Y2, 2));  //Propagate STAT. UNC.
        Frag_f_Syst_U[i] = fabs(Frag_f[i]) * sqrt(pow(Y1_syst / B_Y1, 2) + pow(Y2_syst / B_Y2, 2));  //Propagate SYST. UNC.

        cout << i << "-th BIN FRAGratio: " <<  Frag_f[i] << " +/- " << Frag_f_Stat_U[i] << " +/- " << Frag_f_Syst_U[i] << endl;
        cout << "Bin Xposition " << B_X1 <<"   " <<B_X2 << "  --->  "<< X_POS[i] << endl;
        //graph1->GetErrorXlow(i)
    }

	//center of the bin and its left and right margins
	double ptBins[3];
	double XsecPP_X_BinLeft[3] ;
    double XsecPP_X_BinRight[3] ;

	for(int i = 0; i < 3 + 1; i++){
		if (whichvar=="p_{T}"){ ptBins[i] = ptbinsvec[i];}
		else if (whichvar=="|y|"){ ptBins[i] =  ybinsvec[i];}          
		else if (whichvar=="BMult"){ ptBins[i] =  nmbinsvec[i];}
	}
	for( int c=0; c < 3; c++){
		XsecPP_X_BinLeft[c] = X_POS[c] - ptBins[c];
		XsecPP_X_BinRight[c]= ptBins[c+1] - X_POS[c];
	}
	//center of the bin and its left and right margins






    // Clean up
    file1->Close();
    file2->Close();

    // Plot the result
	TH2D * HisEmpty = new TH2D("HisEmpty","",100,0,2.4,100,0,0.8); 
    TString units = "";
    if( whichvar == "p_{T}"){ 
        units = "[GeV/c]" ;
        HisEmpty = new TH2D("HisEmpty","",100,7,50,100,0,0.8); }
	HisEmpty->GetXaxis()->SetTitle(Form("%s %s",whichvar.Data(), units.Data()));
	HisEmpty->GetYaxis()->SetTitle("B^{0}_{s}/B^{+}");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
    HisEmpty->SetStats(0);

    TCanvas* canvas = new TCanvas("canvas", "Result of Division", 600, 600);
    canvas->cd();

    // Create a new TGraphAsymmErrors for the result of the division
    TGraphAsymmErrors* FragRatio_stat = new TGraphAsymmErrors(3, X_POS, Frag_f, XsecPP_X_BinLeft, XsecPP_X_BinRight, Frag_f_Stat_U, Frag_f_Stat_U);     
	FragRatio_stat->SetLineColor(1); 
    TGraphAsymmErrors* FragRatio_syst = new TGraphAsymmErrors(3, X_POS, Frag_f, nullptr, nullptr, Frag_f_Syst_U, Frag_f_Syst_U);     
	FragRatio_syst->SetLineColor(2); 

	HisEmpty->Draw();
	FragRatio_syst->Draw("5same");
	FragRatio_stat->Draw("epSAME");

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.025); 
	lat->SetTextFont(42);
	lat->DrawLatex(0.1,0.91 , "CMS work in progress");
    //lat->DrawLatex(0.57,0.62 ,Form("2017 pp global Unc. #pm %.1f%%",0000)) ;

	TLegend* leg3 = new TLegend(0.6,0.68,0.9,0.85,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetFillStyle(0);
	leg3->AddEntry((TObject*)0, "B^{0}_{s}/B^{+}", "");
	leg3->AddEntry(FragRatio_stat,"2017 pp ","P");
	//leg3->AddEntry(,"2017 pp (|y|>1.5)","P");
	leg3->SetTextSize(0.022);
	leg3->Draw("same");

    gSystem->mkdir("./Plots/Fragmentation",true); 
	if( whichvar == "p_{T}"){ canvas->SaveAs(Form("./Plots/Fragmentation/fs_fu_pT.pdf"));}
    else { canvas->SaveAs(Form("./Plots/Fragmentation/fs_fu_%s.pdf", whichvar.Data() ));}



}
