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
        cout << "Bin Xposition " << B_X1 <<"   " << B_X2 << "  --->  "<< X_POS[i] << endl;
    }

	//center of the bin and its left and right margins
	double ptBins[NumberBin+1];
	double X_BinLeft[NumberBin+1] ;
    double X_BinRight[NumberBin+1] ;

    int binlow = 0;

	for(int i = 0; i < NumberBin + 1; i++){
		if (whichvar=="p_{T}"){ ptBins[i] = ptbinsvec[i];}
		else if (whichvar=="|y|"){ ptBins[i] =  ybinsvec[i];}          
		else if (whichvar=="BMult"){ ptBins[i] =  nmbinsvec[i];}
	}
	for( int c=0; c < NumberBin+1; c++){
		X_BinLeft[c] = X_POS[c] - ptBins[c];
		X_BinRight[c]= ptBins[c+1] - X_POS[c];
        if ((whichvar=="p_{T}") && X_POS[c]< 10){binlow += 1;}
        else if(whichvar == "|y|" && abs(X_POS[c]) < 1.5){binlow += 1;}
	}
    int binhigh = NumberBin- binlow;
    cout << "Number of lowBINS" << binlow << endl;
	cout << "Number of lowBINS" << binhigh << endl;	
    
    //center of the bin and its left and right margins

    //Separate according to the Fid Region
    double X_POS_ycut[binlow+1] ;
    double X_BinLeft_ycut[binlow+1] ;
    double X_BinRight_ycut[binlow+1] ;
    double X_POS_yall[binhigh+1] ;
    double X_BinLeft_yall[binhigh+1] ;
    double X_BinRight_yall[binhigh+1] ;

    double Frag_f_ycut[binlow+1] ;
    double Frag_f_Stat_U_ycut[binlow+1] ;
    double Frag_f_Syst_U_ycut [binlow+1] ;
    double Frag_f_yall[binhigh+1] ;
    double Frag_f_Stat_U_yall[binhigh+1] ;
    double Frag_f_Syst_U_yall[binhigh+1] ;

	for( int c=0; c < NumberBin+1; c++){
        if(c<binlow){
            X_POS_ycut[c] = X_POS[c] ;
            X_BinLeft_ycut[c] = X_BinLeft[c] ;
            X_BinRight_ycut[c] = X_BinRight[c];
            Frag_f_ycut[c] = Frag_f[c];
            Frag_f_Stat_U_ycut[c] = Frag_f_Stat_U[c];
            Frag_f_Syst_U_ycut[c] = Frag_f_Syst_U[c];

        } else {
            X_POS_yall[c-binlow] = X_POS[c];
            X_BinLeft_yall [c-binlow] = X_BinLeft[c];
            X_BinRight_yall[c-binlow] = X_BinRight[c];
            Frag_f_yall[c-binlow] = Frag_f[c];
            Frag_f_Stat_U_yall[c-binlow] = Frag_f_Stat_U[c];
            Frag_f_Syst_U_yall[c-binlow] = Frag_f_Syst_U[c];
        }
	}
    //Separate according to the Fid Region

    // Clean up
    file1->Close();
    file2->Close();

    //CANVAS
	gStyle->SetOptStat(0);
	TCanvas * canvas = new TCanvas("c","c",600,600);
	canvas->cd(); 
	canvas->SetLeftMargin(0.15);

    // Plot the result
	TH2D * HisEmpty = new TH2D("HisEmpty","",100,0,2.4,100,0,0.5); 
    TString units = "";
    if( whichvar == "p_{T}"){ 
        units = "[GeV/c]" ;
        HisEmpty = new TH2D("HisEmpty","",100,7,50,100,0,0.5); 
    }
	HisEmpty->GetXaxis()->SetTitle(Form("%s %s",whichvar.Data(), units.Data()));
	HisEmpty->GetYaxis()->SetTitle("f_{s}/f_{u}");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
    HisEmpty->SetStats(0);


    // Create a new TGraphAsymmErrors for the result of the division
	TGraphAsymmErrors *FragRatio_Graph_low_marker= new TGraphAsymmErrors(binlow, X_POS_ycut, Frag_f_ycut ,nullptr, nullptr, nullptr, nullptr);
	TGraphAsymmErrors *FragRatio_Graph_Stat_low  = new TGraphAsymmErrors(binlow, X_POS_ycut, Frag_f_ycut , X_BinLeft_ycut , X_BinRight_ycut , Frag_f_Stat_U_ycut , Frag_f_Stat_U_ycut);
	TGraphAsymmErrors *FragRatio_Graph_Syst_low  = new TGraphAsymmErrors(binlow, X_POS_ycut, Frag_f_ycut , X_BinLeft_ycut , X_BinRight_ycut , Frag_f_Syst_U_ycut , Frag_f_Syst_U_ycut);                 											
    TGraphAsymmErrors *FragRatio_Graph_Stat      = new TGraphAsymmErrors(binhigh, X_POS_yall, Frag_f_yall, X_BinLeft_yall, X_BinRight_yall, Frag_f_Stat_U_yall, Frag_f_Stat_U_yall);     
	TGraphAsymmErrors *FragRatio_Graph_Syst      = new TGraphAsymmErrors(binhigh, X_POS_yall, Frag_f_yall, X_BinLeft_yall, X_BinRight_yall, Frag_f_Syst_U_yall, Frag_f_Syst_U_yall);                 											

	FragRatio_Graph_Stat->SetMarkerStyle(20);
	FragRatio_Graph_Stat->SetMarkerSize(1);
	FragRatio_Graph_low_marker ->SetMarkerStyle(20);
	FragRatio_Graph_low_marker ->SetMarkerSize(0.9);
	FragRatio_Graph_low_marker ->SetMarkerColor(kWhite);
	FragRatio_Graph_Stat_low ->SetMarkerStyle(24);
	FragRatio_Graph_Stat_low ->SetMarkerSize(1);

	FragRatio_Graph_Stat_low ->SetMarkerColor(kCyan+2);
	FragRatio_Graph_Stat_low ->SetLineColor(kCyan+2);
	FragRatio_Graph_Syst_low ->SetFillColorAlpha(kCyan-7,0.5);
	FragRatio_Graph_Syst_low ->SetLineColor(kCyan-7);
	FragRatio_Graph_Stat ->SetMarkerColor(kCyan+2);
	FragRatio_Graph_Stat ->SetLineColor(kCyan+2);
	FragRatio_Graph_Syst ->SetFillColorAlpha(kCyan-7,0.5);
	FragRatio_Graph_Syst ->SetLineColor(kCyan-7);

	HisEmpty->Draw();
	FragRatio_Graph_Syst_low->Draw("5same");
	FragRatio_Graph_Syst->Draw("5same");
	FragRatio_Graph_Stat->Draw("epSAME");
	FragRatio_Graph_Stat_low->Draw("epSAME");
	FragRatio_Graph_low_marker->Draw("epSAME");

    //LABELS
	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.025); 
	lat->SetTextFont(42);
	lat->DrawLatex(0.15,0.91 , "CMS work in progress");
	TLegend* leged = new TLegend(0.65,0.74,0.9,0.85,NULL,"brNDC");
	leged->SetBorderSize(0);
	leged->SetFillStyle(0);
	if (whichvar=="p_{t}"){
		leged->AddEntry((TObject*)0, "y region:", "");
		leged->AddEntry(FragRatio_Graph_Stat,"|y|<2.4","P");
		leged->AddEntry(FragRatio_Graph_Stat_low,"|y|>1.5","P");
	} 
	if (whichvar=="|y|"){
		leged->AddEntry((TObject*)0, "p_{T} region:", "");
        leged->AddEntry(FragRatio_Graph_Stat,"7<p_{T}<50 GeV/c","P");
		leged->AddEntry(FragRatio_Graph_Stat_low,"10 < p_{T} < 50 GeV/c","P");
	} 
    leged->SetTextSize(0.022);
	leged->Draw("same");
    //LABELS

    gSystem->mkdir("./Plots/Fragmentation",true); 
	if( whichvar == "p_{T}"){ canvas->SaveAs(Form("./Plots/Fragmentation/fs_fu_pT.pdf"));}
    else { canvas->SaveAs(Form("./Plots/Fragmentation/fs_fu_%s.pdf", whichvar.Data() ));}





    // Fragmentation Fraction Double Ratios
    if(DRatios==1){

        double X_POS_DR[NumberBin+1] ;
        double Frag_f_DR[NumberBin+1];
        double Frag_f_Stat_U_DR[NumberBin+1];
        double Frag_f_Syst_U_DR[NumberBin+1];
	    double X_BinLeft_DR[NumberBin+1] ;
        double X_BinRight_DR[NumberBin+1] ;

    //compute the Fr fraction double ratio
    for (int i = 0; i < NumberBin+1; ++i) {

        X_POS_DR[i] = (FRfr2018_X[i]+X_POS[i])/2;     //mean of both X positions ?
        Frag_f_DR[i] = FRfr2018_Y[i] / Frag_f[i];
        Frag_f_Stat_U_DR[i] = fabs(Frag_f_DR[i]) * sqrt( pow(FRfr2018_Y_StatUpRatio[i] / FRfr2018_Y[i], 2) + pow(Frag_f_Stat_U[i] / Frag_f[i], 2));  //Propagate STAT. UNC.
        Frag_f_Syst_U_DR[i] = fabs(Frag_f_DR[i]) * sqrt( pow(FRfr2018_Y_SystUpRatio[i] / FRfr2018_Y[i], 2) + pow(Frag_f_Syst_U[i] / Frag_f[i], 2));  //Propagate SYST. UNC.
        cout << i << "-th BIN DOUBLE FRAGratio: " <<  Frag_f_DR[i] << " +/- " << Frag_f_Stat_U_DR[i] << " +/- " << Frag_f_Syst_U_DR[i] << endl;
		
        X_BinLeft_DR[i] = X_POS_DR[i] - ptBins[i];
		X_BinRight_DR[i]= ptBins[i+1] - X_POS_DR[i];
	}

    //Separate according to the Fid Region
    double X_POS_DR_ycut[binlow] ;
    double X_BinLeft_DR_ycut[binlow] ;
    double X_BinRight_DR_ycut[binlow] ;
    double Frag_f_DR_ycut[binlow] ;
    double Frag_f_Stat_U_DR_ycut[binlow] ;
    double Frag_f_Syst_U_DR_ycut [binlow] ;
    double X_POS_DR_yall[binhigh] ;
    double X_BinLeft_DR_yall[binhigh] ;
    double X_BinRight_DR_yall[binhigh] ;
    double Frag_f_DR_yall[binhigh] ;
    double Frag_f_Stat_U_DR_yall[binhigh] ;
    double Frag_f_Syst_U_DR_yall[binhigh] ;

	for( int c=0; c < NumberBin+1; c++){
        if(c<binlow){
            X_POS_DR_ycut[c] = X_POS_DR[c] ;
            X_BinLeft_DR_ycut[c] = X_BinLeft_DR[c] ;
            X_BinRight_DR_ycut[c] = X_BinRight_DR[c];
            Frag_f_DR_ycut[c] = Frag_f_DR[c];
            Frag_f_Stat_U_DR_ycut[c] = Frag_f_Stat_U_DR[c];
            Frag_f_Syst_U_DR_ycut[c] = Frag_f_Syst_U_DR[c];

        } else {
            X_POS_DR_yall[c-binlow] = X_POS_DR[c];
            X_BinLeft_DR_yall [c-binlow] = X_BinLeft_DR[c];
            X_BinRight_DR_yall[c-binlow] = X_BinRight_DR[c];
            Frag_f_DR_yall[c-binlow] = Frag_f_DR[c];
            Frag_f_Stat_U_DR_yall[c-binlow] = Frag_f_Stat_U_DR[c];
            Frag_f_Syst_U_DR_yall[c-binlow] = Frag_f_Syst_U_DR[c];
        }
	}
    //Separate according to the Fid Region

    // Plot the result
	TH2D * HisEmpty3 = new TH2D("HisEmpty","",100,7,50,100,0.7,3.5); 
	HisEmpty3->GetXaxis()->SetTitle("p_{T}");
	HisEmpty3->GetYaxis()->SetTitle("fsfu Double Ratio ");
	HisEmpty3->GetXaxis()->CenterTitle();
	HisEmpty3->GetYaxis()->CenterTitle();
    HisEmpty3->SetStats(0);

    TCanvas* canvas3 = new TCanvas("canvas3", "Result of Division", 700, 700);
    canvas3->cd();

    // Create a new TGraphAsymmErrors for the result of the division
    TGraphAsymmErrors* FragRatio_stat_DR = new TGraphAsymmErrors(NumberBin, X_POS_DR, Frag_f_DR, X_BinLeft_DR, X_BinRight_DR, Frag_f_Stat_U_DR, Frag_f_Stat_U_DR);     
	FragRatio_stat_DR->SetLineColor(1); 
    TGraphAsymmErrors* FragRatio_syst_DR = new TGraphAsymmErrors(NumberBin, X_POS_DR, Frag_f_DR, nullptr, nullptr, Frag_f_Syst_U_DR, Frag_f_Syst_U_DR);     
	FragRatio_syst_DR->SetLineColor(2); 

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
