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

void divideTGraphsInFiles(TString inputFile1, TString inputFile2, TString whichvar, int DRatios = 0 );



void Bmeson_Ratio(){
//First Retrive the BP Xsection with matching bins
//First Retrive the BP Xsection with matching bins

// fs/fu 
// pT Fragmentation Fraction

TString BP_FILE_pt = "./ROOTFiles/BP_Xsection_pt_BsBPBINS.root";
TString Bs_FILE_pt = "./ROOTFiles/Bs_Xsection_pt.root";
divideTGraphsInFiles(Bs_FILE_pt, BP_FILE_pt, "p_{T}", 1) ;

// pT Fragmentation Fraction
// y Fragmentation Fraction

TString BP_FILE_y = "./ROOTFiles/BP_Xsection_y_BsBPBINS.root";
TString Bs_FILE_y = "./ROOTFiles/Bs_Xsection_y.root";
divideTGraphsInFiles(Bs_FILE_y, BP_FILE_y, "|y|") ;

// y Fragmentation Fraction
// fs/fu 

}




void divideTGraphsInFiles(TString inputFile1, TString inputFile2, TString whichvar, int DRatios = 0) {
    //inputFile1 SHOULD BE Bmeson NUMERATOR
    //inputFile2 SHOULD BE Bmeson DENOMINATOR      
    TString Name1= "";
    TString Name2= "";
    int NumberBin =0 ;

    if (whichvar == "p_{T}"){
        Name2= "BP_Xsection_pt";
        Name1= "Bs_Xsection_pt";
        NumberBin = nptBins;
    } else {
        Name2= "BP_Xsection_y";
        Name1= "Bs_Xsection_y";
        NumberBin = nyBins_both;
    }

    TFile* file1 = TFile::Open(inputFile1, "READ");
	TMultiGraph *mg1= (TMultiGraph*) file1->Get(Name1);
    TGraphAsymmErrors* Y_stat_B1 = dynamic_cast<TGraphAsymmErrors*>(mg1->GetListOfGraphs()->FindObject("Y_stat"));
    TGraphAsymmErrors* Y_syst_B1 = dynamic_cast<TGraphAsymmErrors*>(mg1->GetListOfGraphs()->FindObject("Y_syst"));

    TFile* file2 = TFile::Open(inputFile2, "READ");
	TMultiGraph *mg2= (TMultiGraph*) file2->Get(Name2);
    TGraphAsymmErrors* Y_stat_B2 = dynamic_cast<TGraphAsymmErrors*>(mg2->GetListOfGraphs()->FindObject("Y_stat"));
    TGraphAsymmErrors* Y_syst_B2 = dynamic_cast<TGraphAsymmErrors*>(mg2->GetListOfGraphs()->FindObject("Y_syst"));

    double X_POS[NumberBin+1] ;
    double Frag_f[NumberBin+1];
    double Frag_f_Stat_U[NumberBin+1];
    double Frag_f_Syst_U[NumberBin+1];

    Double_t B_X1, B_Y1, B_X2, B_Y2, Y1_stat, Y1_syst, Y2_stat, Y2_syst ;
    for (int i = 0; i < NumberBin+1; ++i) {

        Y_stat_B1->GetPoint(i, B_X1, B_Y1);      // X and Yield values of NUMERATOR Bmeson
        Y1_stat = Y_stat_B1->GetErrorYhigh(i);   // Yield Statistical Unc. of NUMERATOR Bmeson
        Y1_syst = Y_syst_B1->GetErrorYhigh(i);   // Yield Systematic Unc. of NUMERATOR Bmeson

        Y_stat_B2->GetPoint(i, B_X2, B_Y2);    // X and Yield values of DENOMINATOR Bmeson
        Y2_stat = Y_stat_B2->GetErrorYhigh(i);   // Yield Statistical Unc. of DENOMINATOR Bmeson
        Y2_syst = Y_syst_B2->GetErrorYhigh(i);   // Yield Systematic Unc. of DENOMINATOR Bmeson

        X_POS[i] = (B_X1+B_X2)/2;     //mean of both X1 x2
        Frag_f[i] = B_Y1 / B_Y2;
        Frag_f_Stat_U[i] = fabs(Frag_f[i]) * sqrt(pow(Y1_stat / B_Y1, 2) + pow(Y2_stat / B_Y2, 2));  //Propagate STAT. UNC.
        Frag_f_Syst_U[i] = fabs(Frag_f[i]) * sqrt(pow(Y1_syst / B_Y1, 2) + pow(Y2_syst / B_Y2, 2));  //Propagate SYST. UNC.

        cout << i << "-th BIN FRAGratio: " <<  Frag_f[i] << " +/- " << Frag_f_Stat_U[i] << " +/- " << Frag_f_Syst_U[i] << endl;
        cout << "Bin Xposition " << B_X1 <<"   " <<B_X2 << "  --->  "<< X_POS[i] << endl;
    }

	//center of the bin and its left and right margins
	double ptBins[NumberBin+1];
	double XsecPP_X_BinLeft[NumberBin+1] ;
    double XsecPP_X_BinRight[NumberBin+1] ;

	for(int i = 0; i < NumberBin + 1; i++){
		if (whichvar=="p_{T}"){ ptBins[i] = ptbinsvec[i];}
		else if (whichvar=="|y|"){ ptBins[i] =  ybinsvec[i];}          
		else if (whichvar=="BMult"){ ptBins[i] =  nmbinsvec[i];}
	}
	for( int c=0; c < NumberBin+1; c++){
		XsecPP_X_BinLeft[c] = X_POS[c] - ptBins[c];
		XsecPP_X_BinRight[c]= ptBins[c+1] - X_POS[c];
	}
	//center of the bin and its left and right margins

    // Clean up
    file1->Close();
    file2->Close();

    // Plot the result
	TH2D * HisEmpty = new TH2D("HisEmpty","",100,0,2.4,100,0,0.5); 
    TString units = "";
    if( whichvar == "p_{T}"){ 
        units = "[GeV/c]" ;
        HisEmpty = new TH2D("HisEmpty","",100,7,50,100,0,0.5); }
	HisEmpty->GetXaxis()->SetTitle(Form("%s %s",whichvar.Data(), units.Data()));
	HisEmpty->GetYaxis()->SetTitle("B^{0}_{s}/B^{+}");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
    HisEmpty->SetStats(0);

    TCanvas* canvas = new TCanvas("canvas", "Result of Division", 700, 700);
    canvas->cd();

    // Create a new TGraphAsymmErrors for the result of the division
    TGraphAsymmErrors* FragRatio_stat = new TGraphAsymmErrors(NumberBin, X_POS, Frag_f, XsecPP_X_BinLeft, XsecPP_X_BinRight, Frag_f_Stat_U, Frag_f_Stat_U);     
	FragRatio_stat->SetLineColor(1); 
    TGraphAsymmErrors* FragRatio_syst = new TGraphAsymmErrors(NumberBin, X_POS, Frag_f, nullptr, nullptr, Frag_f_Syst_U, Frag_f_Syst_U);     
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





    // Fragmentation Fraction Double Ratios
    if(DRatios==1 && NumberBin==4 ){

        double X_POS_DR[NumberBin+1] ;
        double Frag_f_DR[NumberBin+1];
        double Frag_f_Stat_U_DR[NumberBin+1];
        double Frag_f_Syst_U_DR[NumberBin+1];
	    double XsecPP_X_BinLeft_DR[NumberBin+1] ;
        double XsecPP_X_BinRight_DR[NumberBin+1] ;

    //compute the Fr fraction double ratio
    for (int i = 0; i < NumberBin+1; ++i) {

        X_POS_DR = (FRfr2018_X[i]+X_POS[i])/2;     //mean of both X positions ?
        Frag_f_DR[i] = FRfr2018_Y[i] / Frag_f[i];
        Frag_f_Stat_U_DR[i] = fabs(Frag_f[i]) * sqrt(pow(FRfr2018_Y_StatUpRatio[i] / FRfr2018_Y[i], 2) + pow(Frag_f_Stat_U[i] / Frag_f[i], 2));  //Propagate STAT. UNC.
        Frag_f_Syst_U_DR[i] = fabs(Frag_f[i]) * sqrt(pow(FRfr2018_Y_SystUpRatio[i] / FRfr2018_Y[i], 2) + pow(Frag_f_Syst_U[i] / Frag_f[i], 2));  //Propagate SYST. UNC.

        cout << i << "-th BIN DOUBLE FRAGratio: " <<  Frag_f_DR[i] << " +/- " << Frag_f_Stat_U_DR[i] << " +/- " << Frag_f_Syst_U_DR[i] << endl;
    }
    //compute the Fr fraction double ratio
    //left and right limits of the bin
	for( int c=0; c < NumberBin+1; c++){
		XsecPP_X_BinLeft_DR[c] = X_POS_DR[c] - ptBins[c];
		XsecPP_X_BinRight_DR[c]= ptBins[c+1] - X_POS_DR[c];
	}
    //left and right limits of the bin

    // Plot the result
	TH2D * HisEmpty3 = new TH2D("HisEmpty","",100,7,50,100,0.5,1.5); 
	HisEmpty3->GetXaxis()->SetTitle("p_{T}");
	HisEmpty3->GetYaxis()->SetTitle("fsfu Double Ratio ");
	HisEmpty3->GetXaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->CenterTitle();
    HisEmpty3->SetStats(0);

    TCanvas* canvas3 = new TCanvas("canvas3", "Result of Division", 700, 700);
    canvas3->cd();

    // Create a new TGraphAsymmErrors for the result of the division
    TGraphAsymmErrors* FragRatio_stat_DR = new TGraphAsymmErrors(NumberBin, X_POS_DR, Frag_f_DR, XsecPP_X_BinLeft_DR, XsecPP_X_BinRight_DR, Frag_f_Stat_U_DR, Frag_f_Stat_U_DR);     
	FragRatio_stat_DR->SetLineColor(1); 
    TGraphAsymmErrors* FragRatio_syst_DR = new TGraphAsymmErrors(NumberBin, X_POS_DR, Frag_f_DR, nullptr, nullptr, Frag_f_Syst_U_DR, Frag_f_Syst_U_DR);     
	FragRatio_syst->SetLineColor(2); 

	HisEmpty3->Draw();
	FragRatio_syst_DR->Draw("5same");
	FragRatio_stat_DR->Draw("epSAME");

	TLatex *lat2 = new TLatex();
	lat2->SetNDC();
	lat2->SetTextSize(0.025); 
	lat2->SetTextFont(42);
	lat2->DrawLatex(0.1,0.91 , "CMS work in progress");
    //lat->DrawLatex(0.57,0.62 ,Form("2017 pp global Unc. #pm %.1f%%",0000)) ;

	TLegend* leg31 = new TLegend(0.6,0.68,0.9,0.85,NULL,"brNDC");
	leg31->SetBorderSize(0);
	leg31->SetFillStyle(0);
	leg31->AddEntry((TObject*)0, "B^{0}_{s}/B^{+}", "");
	leg31->AddEntry(FragRatio_stat_DR,"2017 pp ","P");
	//leg3->AddEntry(,"2017 pp (|y|>1.5)","P");
	leg31->SetTextSize(0.022);
	leg31->Draw("same");

	canvas3->SaveAs(Form("./Plots/Fragmentation/DOUBLER_fs_fu_pT.pdf"));


    }

}
