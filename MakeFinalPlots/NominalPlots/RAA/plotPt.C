/*
   Macro to plot xsec and ratio vs pt for Bs and Bp)

Input: txt files in inputDir, with 10 columns:
ptmin, ptmax, central val, statUp, statDown, systUp, systDown, glbUp, glbDown, abscissae

Output: xsec vs pt, ratio vs pt.

*/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <Riostream.h>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"

#include "TPaveStats.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "CMS_lumi.C"
#include "tdrstyle.C"
#include "DrawLHCb.C"

#include "auxiliaryPt.h"
#include "auxiliaryRef.h"

//#include "outsideSource/lhcb.C"
#endif
using namespace std;

void plotPt(bool bSavePlots       = 1,
		bool bDoDebug         = 0, //  figure out if things are read properly
		bool whichPlot        = 1, //0 is x-sec, 1 is for ratio
		bool drawRef          = 1, //draw Ref (for ratio only
		bool drawlhcb         = 1,
		const char* inputDir  = "dataSource", //inptu txt files
		const char* outputDir = "figs")// where the output figures will be
		{
		gSystem->mkdir(Form("./%s/png",outputDir), kTRUE);
		gSystem->mkdir(Form("./%s/pdf",outputDir), kTRUE);

		//   gROOT->ProcessLine(".x lhcb.C");
		//  lhcb();
		//set the style
		setTDRStyle();

		//samples:
		const unsigned int nMes      = 2;
		const char* inputFileType[2] = {"RAA_pt", "RAA_pt"};
		const char* mesonName[nMes]  = {"Bs", "Bp"};

		Int_t endMes = nMes;

    std::vector<int> nlines = {4, 7};
		for (Int_t ib=0; ib<endMes; ib++){
			ifstream in;
			string inputFileName = Form("%s/%s_%s_New.txt",inputDir,inputFileType[whichPlot],mesonName[ib]);
			//if(whichPlot==1) inputFileName = Form("%s/%s.txt",inputDir,inputFileType[whichPlot]);
	//		if(whichPlot==1)  inputFileName = Form("%s/%s_%s.txt",inputDir,inputFileType[whichPlot],mesonName[ib]);		
			cout << "########## Input file name: " << inputFileName << endl;

			in.open(inputFileName.c_str());
			if (!in.is_open()) {
				cout << "input file " << inputFileName << " cannot be open" << endl;
				continue;
			}
			double x[20];

			//get first line, the header, and discard it
			string tmpstrg;
			getline(in,tmpstrg);//ignore first line/ the header

      cout << "start reading" << endl;
		int nEntry=0;
		for (auto iline = 0; iline < nlines[ib]; ++iline){
        	in >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5] >> x[6] >> x[7] >> x[8] >> x[9];
			glbSystDown = x[8]*100;
			glbSystUp   = x[7]*100;

        cout << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ", " << x[4] << ", " << x[5] << ", " << x[6] << ", " << x[7] << ", " << x[8] << ", " << x[9] << endl;

					if(ib==0){//bs
						if(nEntry==0){
						binLow[nEntry]    = x[9];							
						glbSystUpBs  = x[8]*100;
						glbSystDownBs   = x[7]*100;
						bs_low[nEntry] = x[2]; //central value
						//stat uncert
						bs_low_yStatL[nEntry] = x[4]*x[2];
						bs_low_yStatH[nEntry] =  x[3]*x[2];
						//bin width
						bs_low_xErrL[nEntry] = x[9]-x[0];
						bs_low_xErrH[nEntry] = x[1]-x[9];
						//systm. uncert
						bs_low_ySystL[nEntry] = x[6]*x[2];
						bs_low_ySystH[nEntry] = x[5]*x[2];
						cout << "nEntry = " << nEntry << "   x[1] =  " << x[1] << "   x[9] = " << x[9] << "  bs_low_xErrH[0] =  " << bs_low_xErrH[0] << endl;

						if(bDoDebug){
							cout << "Central  \t statUncertL \t systUncerL " << endl;
							cout<< bs_low[nEntry] << "\t" << bs_low_yStatL[nEntry] << "\t" <<   bs_low_ySystL[nEntry] << endl;
							cout<<"==========================================="<<endl;
						}
						}
					}
					else{//bp
						//central value	
						if(nEntry<2){
						binLowNew[nEntry]    = x[9];	
						glbSystUpBp = x[8]*100;
						glbSystDownBp   = x[7]*100;

						bpl_low[nEntry] = x[2];
						//stat uncert
						bpl_low_yStatL[nEntry] = x[4]*x[2];
						bpl_low_yStatH[nEntry] =  x[3]*x[2];
						//bin width
						bpl_low_xErrL[nEntry] = x[9]-x[0];
						bpl_low_xErrH[nEntry] = x[1]-x[9];
						//systm. uncert
						bpl_low_ySystL[nEntry] = x[6]*x[2];
						bpl_low_ySystH[nEntry] = x[5]*x[2];

						if(bDoDebug){
							cout << "Central  \t statUncertL  \t systUncertL " << endl;
							cout<< bpl_low[nEntry] << "\t" << bpl_low_yStatL[nEntry] << "\t" <<   bpl_low_ySystL[nEntry] << endl;
							cout<<"==========================================="<<endl;
						}
					}
					
				}
				// high-pt bins, |y|<2.4
					if(ib==0){//bs
						if(nEntry>0){
						//central value
						binHigh[nEntry-1] = x[9];										
						bs_high[nEntry-1] = x[2];
						//stat uncert
						bs_high_yStatL[nEntry-1] = x[4]*x[2];
						bs_high_yStatH[nEntry-1] = x[3]*x[2];
						//bin width
						bs_high_xErrL[nEntry-1] = x[9]-x[0];
						bs_high_xErrH[nEntry-1] = x[1]-x[9];
						//systm. uncert
						bs_high_ySystL[nEntry-1] = x[6]*x[2];
						bs_high_ySystH[nEntry-1] = x[5]*x[2];
						//	cout << "Bs Bro: " << " bs_high_ySystL = " << bs_high_ySystL[nEntry-1] << "  bs_high_ySystH = " << bs_high_ySystH[nEntry-1] << endl;

						if(bDoDebug){
							if(nEntry==1) cout << "Central  \t statUncertL  \t systUncertL "<< endl;
							cout<< bs_high[nEntry-2] << "\t" << bs_high_yStatL[nEntry-2] << "\t" << bs_high_ySystL[nEntry-2] << endl;
						}
						}
					}else{//bp
						if(nEntry>1){
						//central value
						binHighNew[nEntry-2] = x[9];				
						bpl_high[nEntry-2] = x[2];
						//stat uncert
						bpl_high_yStatL[nEntry-2] = x[4]*x[2];
						bpl_high_yStatH[nEntry-2] = x[3]*x[2];
						//bin width
						bpl_high_xErrL[nEntry-2] = x[9]-x[0];
						bpl_high_xErrH[nEntry-2] = x[1]-x[9];
						//systm. uncert
						bpl_high_ySystL[nEntry-2] = x[6]*x[2];
						bpl_high_ySystH[nEntry-2] = x[5]*x[2];

            cout << "ysyst:" << bpl_high_ySystL[nEntry-2] << ", ratio: " <<
              bpl_high_ySystL[nEntry-2] / bpl_high[nEntry-2] << "\n";


						}
					
				}//high-pt bins
				nEntry++;
			}//reading input file line by line

			cout<<"@@@@@@@@@@@ Finished meson "<<mesonName[ib]<<endl;
			if(bDoDebug){
				if(ib==0) cout<<"Element_low = "<< binLow[0] <<"\t bs_low[0]= " << bs_low[0] << "\t statUncertL = "<<bs_low_yStatL[0]<<endl;

				if(ib==1)
				{
					for(int i=0; i<3; i++){
						cout<<"Element_high " << i << "\t binHigh[i] = "<< binHigh[i] <<"\t bs_high[i]= " << bs_high[i] << "\t statUncertL = "<<bs_high_yStatL[i]<<"\t systUncertL = "<<bs_high_ySystL[i]<< endl;
					}
				}
			}//bDebug
			in.close();//close input file

		}//for each meson,ib
//		bs_low_xErrH[0] = 1.0;
//		bs_low_xErrL[1] = 1.73;
		
//		bs_high_xErrH[4] = 5.0;


		//----------------------------------------------------------------
		// gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
		// Bs
		TGraphAsymmErrors *pgBs_low = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_yStatL,bs_low_yStatH);
		TGraphAsymmErrors *pgBs_high= new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
				bs_high_xErrL, bs_high_xErrH,
				bs_high_yStatL,bs_high_yStatH);
		TGraphAsymmErrors *pgBs_lowWhite = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_yStatL,bs_low_yStatH);   //AddWhite//
		

		// Bplus
		TGraphAsymmErrors *pgBpl_low = new TGraphAsymmErrors(nBinsLowNew, binLowNew, bpl_low,
				bpl_low_xErrL, bpl_low_xErrH,
				bpl_low_yStatL,bpl_low_yStatH);
		
		TGraphAsymmErrors *pgBpl_lowWhite = new TGraphAsymmErrors(nBinsLowNew, binLowNew, bpl_low,
				bpl_low_xErrL, bpl_low_xErrH,
				bpl_low_yStatL,bpl_low_yStatH);
		
		TGraphAsymmErrors *pgBpl_high= new TGraphAsymmErrors(nBinsHighNew,binHighNew,bpl_high,
				bpl_high_xErrL, bpl_high_xErrH,
				bpl_high_yStatL,bpl_high_yStatH);

		for(int i = 0; i < 5; i++){

			cout << "binLowHigh = " << binHighNew[i] << "   bs_high_xErrL = " << bs_high_xErrL[i] << "   bs_high_xErrH = " << bs_high_xErrH[i] << endl;
		}



		//Systmeatic uncertainty
		// Bs
		TGraphAsymmErrors *pgBs_syst_low    = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_ySystL,bs_low_ySystH);
		TGraphAsymmErrors *pgBs_syst_high   = new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
				bs_high_xErrL, bs_high_xErrH,
				bs_high_ySystL,bs_high_ySystH);

		// Bp
		TGraphAsymmErrors *pgBpl_syst_low    = new TGraphAsymmErrors(nBinsLowNew, binLowNew,  bpl_low,
				bpl_low_xErrL, bpl_low_xErrH,
				bpl_low_ySystL, bpl_low_ySystH);
		TGraphAsymmErrors *pgBpl_syst_high   = new TGraphAsymmErrors(nBinsHighNew,binHighNew, bpl_high,
				bpl_high_xErrL,bpl_high_xErrH,
				bpl_high_ySystL,bpl_high_ySystH);


		// pgBs_syst_high->Draw("APL");
		//==========================================
		//------------------------------------------
		TGraphAsymmErrors *pgRatio_low    = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_yStatL,bs_low_yStatH);

		TGraphAsymmErrors *pgRatio_lowWhite    = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_yStatL,bs_low_yStatH);

		TGraphAsymmErrors *pgRatio_high   = new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
				bs_high_xErrL, bs_high_xErrH,
				bs_high_yStatL,bs_high_yStatH);

		TGraphAsymmErrors *pgRatio_syst_low = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_ySystL,bs_low_ySystH);

		TGraphAsymmErrors *pgRatio_syst_lowWhite = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_ySystL,bs_low_ySystH);
		
		TGraphAsymmErrors *pgRatio_syst_high= new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
				bs_high_xErrL, bs_high_xErrH,
				bs_high_ySystL,bs_high_ySystH);

		//================== reference
		TGraphAsymmErrors * FragBand = new TGraphAsymmErrors(BandBin,BandX,BandY,BandXErr,BandXErr,BandYErr,BandYErr);
		FragBand->SetName("BandErr");
		FragBand->SetMarkerStyle(20);
		FragBand->SetMarkerSize(0.8);
		FragBand->SetFillColorAlpha(kGreen,0.5);
		FragBand->SetFillStyle(3004);
		FragBand->SetLineWidth(2);
		FragBand->SetLineColor(kGreen);

		TGraphAsymmErrors *LHCb7TeVRef = new TGraphAsymmErrors(LHCb7TeVNPoints, LHCb7TeVX, LHCb7TeVY,
				LHCb7TeVXErrDown, LHCb7TeVXErrUp,
				LHCb7TeVYErr,LHCb7TeVYErr);

		LHCb7TeVRef->SetName("LHCb7TeVRef");

		//==================
		// // **************** marker setup
		// Bs
		// marker style
		pgBs_low->SetMarkerStyle(markerLow[0]);
		pgBs_lowWhite->SetMarkerStyle(markerHigh[0]);
		pgBs_high->SetMarkerStyle(markerHigh[0]);

		pgBpl_lowWhite->SetMarkerStyle(markerHigh[1]);
		pgBpl_low->SetMarkerStyle(markerLow[1]);
		pgBpl_high->SetMarkerStyle(markerHigh[1]);
		pgRatio_low->SetMarkerStyle(markerRatio[0]);
		pgRatio_lowWhite->SetMarkerStyle(markerRatio[1]);
		pgRatio_high->SetMarkerStyle(markerRatio[1]);

		// marker size
		pgBs_low->SetMarkerSize(markerSizeLow[0]);
		pgBs_lowWhite->SetMarkerSize(markerSizeLow[0]);
		pgBs_high->SetMarkerSize(markerSizeHigh[0]);
		pgBpl_lowWhite->SetMarkerSize(markerSizeLow[1]*0.86);
		pgBpl_low->SetMarkerSize(markerSizeLow[1]);
		pgBpl_high->SetMarkerSize(markerSizeHigh[1]);
		pgRatio_low->SetMarkerSize(markerSizeRatio[0]);
		pgRatio_lowWhite->SetMarkerSize(markerSizeRatio[0]);
		pgRatio_high->SetMarkerSize(markerSizeRatio[1]);

		// marker color
		pgBs_low->SetMarkerColor(colorLow[0]);
		pgBs_high->SetMarkerColor(colorHigh[0]);
	
		pgBs_lowWhite->SetMarkerColor(kWhite);
		pgBpl_lowWhite->SetMarkerColor(kWhite);
	
		pgBs_low->SetMarkerColor(kBlue+2);
		pgBs_high->SetMarkerColor(kBlue+2);

		//cout << "Color = " << TColor::GetColor("#C32148") << endl; 

//		pgBs_low->SetMarkerColor(924);
		pgBpl_low->SetMarkerColor(colorLow[1]);
		pgBpl_high->SetMarkerColor(colorHigh[1]);
		pgBpl_low->SetMarkerColor(kGreen+3);	
		pgBpl_high->SetMarkerColor(kGreen+3);

		pgRatio_lowWhite->SetMarkerColor(kWhite);
		pgRatio_low->SetMarkerColor(kRed + 2);
		pgRatio_high->SetMarkerColor(kRed + 2);

		// line color
		pgBs_low->SetLineColor(colorLow[0]);
		pgBs_high->SetLineColor(colorHigh[0]);
		pgBs_low->SetLineColor(kBlue+2);	
		pgBs_high->SetLineColor(kBlue+2);

		pgBpl_low->SetLineColor(colorLow[1]);
		pgBpl_high->SetLineColor(colorHigh[1]);
		pgBpl_low->SetLineColor(kGreen+3);
		pgBpl_high->SetLineColor(kGreen+3);


		pgRatio_lowWhite->SetLineColor(kRed + 2);
		pgRatio_low->SetLineColor(kRed + 2);
		pgRatio_high->SetLineColor(kRed + 2);


		// systematic boxes
		pgBs_syst_low->SetFillColorAlpha(kBlue-9,0.5);
		pgBs_syst_high->SetFillColorAlpha(kBlue-9,0.5);
		pgBpl_syst_low->SetFillColorAlpha(kGreen-9,0.5);
		pgBpl_syst_high->SetFillColorAlpha(kGreen-9,0.5);
		pgRatio_syst_low->SetFillColorAlpha(colorRatio[0],0.2);
		pgRatio_lowWhite->SetFillColor(kWhite);
		pgRatio_syst_high->SetFillColorAlpha(colorRatio[1],0.2);

		//Reference 
		LHCb7TeVRef->SetMarkerColor(kBlue);
		LHCb7TeVRef->SetLineColor(kBlue);
		LHCb7TeVRef->SetMarkerStyle(20);
		LHCb7TeVRef->SetMarkerSize(markerSizeLow[1]);


		//================
		//=============================================
		// read TAMU
		std::ifstream TAMUBsBP("outsideSource/TAMUPT.dat");
		const int NBinsTAMU = 150;
		double TAMUBsBPPt[NBinsTAMU];
		double TAMUBsBPPtErr[NBinsTAMU];
		double TAMUBsBPRatio[NBinsTAMU];
		double TAMUBsBPRatioErr[NBinsTAMU];
		for(int i = 0; i < NBinsTAMU; i++){
			TAMUBsBP >> TAMUBsBPPt[i] >> TAMUBsBPRatio[i];
		}
		for(int i = 0; i < NBinsTAMU; i++){
			TAMUBsBPPtErr[i] = 0.1;
			TAMUBsBPRatioErr[i] = 0.005;
			//     cout  << "TAMUBsBPPt[i] = "  <<  TAMUBsBPPt[i] << "     TAMUBsBPRatio[i] = " <<  TAMUBsBPRatio[i] << endl;
		}
		TGraphErrors* TAMUTheory = new TGraphErrors(NBinsTAMU,TAMUBsBPPt,TAMUBsBPRatio,TAMUBsBPPtErr,TAMUBsBPRatioErr);
		TAMUTheory->SetName("TAMUTheory");
		/*
		   TAMUTheory->SetMarkerStyle(20);
		   TAMUTheory->SetMarkerSize(0.0);
		   TAMUTheory->SetFillColor(kOrange+3);
		   TAMUTheory->SetFillStyle(3002);
		   TAMUTheory->SetLineColor(kOrange);
		   TAMUTheory->SetLineWidth(3);
		   */

		TAMUTheory->SetLineWidth(3);
//		TAMUTheory->SetMarkerColor(kOrange);
//		TAMUTheory->SetLineColor(kOrange);

		TAMUTheory->SetMarkerColor(kOrange + 1);   //Stronger TAMU Color
		TAMUTheory->SetLineColor(kOrange + 1);
		//=============================================

		// read CAO
		std::ifstream CAOBsBP("outsideSource/CAOPT.dat");
		const int NBinsCAO = 40;
		double CAOBsBPPt[NBinsCAO];
		double CAOBsBPPtErr[NBinsCAO];
		double CAOBsBPRatio[NBinsCAO];
		double CAOBsBPRatioErr[NBinsCAO];
		for(int i = 0; i < NBinsCAO; i++){
			CAOBsBP >> CAOBsBPPt[i] >> CAOBsBPRatio[i];
		}
		for(int i = 0; i < NBinsCAO; i++){
			CAOBsBPPtErr[i] = 0.1;
			CAOBsBPRatioErr[i] = 0.005;
			//     cout  << "CAOBsBPPt[i] = "  <<  CAOBsBPPt[i] << "     CAOBsBPRatio[i] = " <<  CAOBsBPRatio[i] << endl;
		}
		TGraphErrors* CAOTheory = new TGraphErrors(NBinsCAO,CAOBsBPPt,CAOBsBPRatio,CAOBsBPPtErr,CAOBsBPRatioErr);
		CAOTheory->SetName("CAOTheory");
		CAOTheory->SetMarkerStyle(20);
		CAOTheory->SetLineWidth(3);
		CAOTheory->SetMarkerColor(kGreen+1);
		CAOTheory->SetLineColor(kGreen+1);

		//-------------------------------------------
		TF1 *f4 = new TF1("f4","1",5,60);
		f4->SetLineWidth(1);
		f4->SetLineColor(1);
		f4->SetLineStyle(1);
		f4->GetYaxis()->SetTitle(yAxName[whichPlot]);
		f4->GetXaxis()->SetTitle(xAxName[0]);
		f4->GetXaxis()->CenterTitle(kTRUE);
		f4->GetYaxis()->CenterTitle();
		
		if(whichPlot==1){
			f4->GetYaxis()->SetTitleSize(0.06*0.83);
			f4->GetXaxis()->SetTitleSize(0.06*0.83);
			f4->GetYaxis()->SetTitleOffset(1.40);
			f4->GetXaxis()->SetTitleOffset(1.20);
						}

		if(whichPlot==0){
			f4->GetYaxis()->SetTitleSize(0.06*0.80);
			f4->GetYaxis()->SetTitleOffset(1.5);
			f4->GetXaxis()->SetTitleOffset(1.05);
			f4->GetXaxis()->SetTitleSize(f4->GetXaxis()->GetTitleSize() * 0.77);
			f4->GetXaxis()->SetTitleOffset(1.18);
			f4->GetXaxis()->SetTitleSize(0.06*0.83);  //Unify Textsize
			f4->GetYaxis()->SetTitleSize(0.06*0.83);
			f4->GetYaxis()->SetTitleOffset(1.40);
			f4->GetXaxis()->SetTitleOffset(1.20);
			cout << "Offset = " << f4->GetYaxis()->GetTitleOffset() << endl;
						}

		f4->GetYaxis()->SetRangeUser(1e2,2e7);
		//if(whichPlot==1) f4->GetYaxis()->SetRangeUser(0.0,1.8);
		if(whichPlot==1) f4->GetYaxis()->SetRangeUser(0.0,4.0);
		if(whichPlot==0) f4->GetYaxis()->SetRangeUser(100.0,2000000); //1.8 -> 0.9
		
		f4->GetYaxis()->SetRangeUser(0.0,2); 
		f4->SetLineColor(0);
		//f4->GetXaxis()->SetNdivisions(-6);

		//---------------- general stuff
		TLatex *lat = new TLatex();
		lat->SetNDC();

		// // ##################################################### x-sec canvas
		TCanvas *pc1 = new TCanvas("pc1","pc1");
		f4->SetLineStyle(2);
		f4->Draw();// axis
		//pc1->SetBottomMargin(0.10);
	
	//	cout << "pc1->GetBottomMargin() = " << pc1->GetBottomMargin() << endl;
	
		pc1->SetBottomMargin(0.135); //Enlarge Margin

		CMS_lumi(pc1,19011,0);

		//Added to fix//
		pgBs_syst_high->SetLineColor(kBlue-9);
		pgBs_syst_low->SetLineColor(kBlue-9);
		pgBpl_syst_high->SetLineColor(kGreen-9);
		pgBpl_syst_low->SetLineColor(kGreen-9);



		if(whichPlot==0)// x-section
		{
			//gPad->SetLogy();
			pgBs_syst_high->SetLineWidth(1);
			pgBs_syst_high->Draw("5same");
			pgBs_syst_low->Draw("5same");	
		//	pgBs_low->SetMarkerStyle(24);
			
			
			pgBs_lowWhite->Draw("P");
			pgBs_low->Draw("P");
			pgBs_high->Draw("P");
			pgBpl_syst_low->Draw("5same");
			pgBpl_syst_high->Draw("5same");
			pgBs_lowWhite->Draw("P");
			pgBpl_lowWhite->Draw("P");
			pgBpl_low->Draw("P");
			pgBpl_high->Draw("P");
			pgBpl_lowWhite->Draw("P");
			TLine * Unity = new TLine(5,1,60,1);
			Unity->SetLineWidth(2);
			Unity->SetLineStyle(2);
			Unity->SetLineColor(1);
			Unity->Draw("SAME");

		}
		else{//ratio plot
			if(drawRef){
				//		FragBand->Draw("5same");
				TAMUTheory->Draw("l");
				CAOTheory->Draw("l");
				LHCb7TeVRef->Draw("P");
			}

			pgRatio_syst_low->SetLineColor(colorRatio[0]);
			pgRatio_syst_high->SetLineColor(colorRatio[1]);
			pgRatio_syst_low->Draw("5same");
			pgRatio_low->Draw("P");
			pgRatio_lowWhite->Draw("P");
			pgRatio_syst_high->Draw("5same");
			pgRatio_high->Draw("P");
		}

		//supplemental info on plot:
		if(whichPlot==0){

			lat->SetTextFont(42);
			lat->SetTextSize(ltxSetTextSize2 * 1.3);
			lat->SetTextSize(ltxSetTextSize4); //Enlarge Labels
			
			//lat->DrawLatex(xsec_ltxText1_xStart,xsec_ltxText1_yStart,"Cent. 0-90%");
			//lat->DrawLatex(xsec_ltxText1_xStart,xsec_ltxText1_yStart,"Centrality 0-90%"); //Expand Cent.
			//lat->DrawLatex(xsec_ltxText1_xStart + 0.25,xsec_ltxText1_yStart,"Centrality 0-90%"); // Move Cent
			//lat->DrawLatex(xsec_ltxText1_xStart + 0.34,xsec_ltxText1_yStart - 0.22,"Centrality 0-90%"); // Move Cent
			cout << "ltxSetTextSize4 = " << ltxSetTextSize4 << endl;

	//		lat->DrawLatex(xsec_ltxText1_xStart,xsec_ltxText1_yStart-0.65,Form("B_{s}^{0} global uncert.: #pm %.1f %%",glbSystUpBs));
	//		lat->DrawLatex(xsec_ltxText1_xStart,xsec_ltxText1_yStart-0.70,Form("B^{+} global uncert.: #pm %.1f %%",glbSystUpBp));

		//	lat->DrawLatex(xsec_ltxText1_xStart + 0.02,xsec_ltxText1_yStart-0.65+0.037,Form("B_{s}^{0} global uncertainty: #pm %.1f%%",glbSystUpBs));
		//	lat->DrawLatex(xsec_ltxText1_xStart + 0.02,xsec_ltxText1_yStart-0.70 + 0.037,Form("B^{+} global uncertainty: #pm %.1f%%",glbSystUpBp));
			lat->DrawLatex(xsec_ltxText1_xStart + 0.30,xsec_ltxText1_yStart - 0.30+ 0.037,Form("global uncertainty: #pm %.1f%%",glbSystUpBp));

		//	cout << "X1 = " << xsec_ltxText1_xStart + 0.02 << "    Y1 = " << xsec_ltxText1_yStart-0.65+0.037 << endl;

		//	cout << "X2 = " << xsec_ltxText1_xStart + 0.02 << "    Y2 = " << xsec_ltxText1_yStart-0.65+0.037 << endl;

			// lat->DrawLatex(xsec_ltxText1_xStart,xsec_ltxText1_yStart-0.7,Form("- %.2f %%",glbSystDown));

			// legend
			//			TLegend *legXSec = new TLegend(legXsec_xLowStart,legXsec_y,legXsec_xLowEnd,legXsec_y+0.15,"B_{s}^{0}                    B^{+}","brNDC");
			//			legXSec->SetBorderSize(0);


			//		TLegend *legXSec = new TLegend(legXsec_xLowStart+0.01,legXsec_y-0.05,legXsec_xLowEnd+0.18,legXsec_y+0.05,"B_{s}^{0}   B^{+}","brNDC");
			//		legXSec->SetBorderSize(0);

			//	TLegend *legXSec = new TLegend(legXsec_xLowStart,legXsec_y,legXsec_xLowEnd,legXsec_y+0.15,"B_{s}^{0}                    B^{+}",



			//		TLegend *legXSec = new TLegend(legXsec_xLowStart+0.01,legXsec_y-0.05,legXsec_xLowEnd+0.18,legXsec_y+0.05,NULL,"brNDC");
			//		legXSec->SetBorderSize(0);



			double ShiftX = 0.05;
			double ShiftY = 0.13;

			lat->SetTextSize(0.05);
			lat->SetTextSize(0.05 * ltxSetTextSize4/ltxSetTextSize2);  //Enlarge Labels
			lat->SetTextSize(0.048 * 1.15);  //Enlarge Labels
		//	lat->SetTextSize(0.025);  //Enlarge Labels
		
			//lat->DrawLatex(legXsec_xLowStart- ShiftX,legXsec_y,"#bf{B_{s}^{0}              B^{+}}");

		//	lat->DrawLatex(legXsec_xLowStart+ShiftX,legXsec_y,"#bf{B_{s}^{0}}");
		//	lat->DrawLatex(legXsec_xLowStart+ ShiftX + 0.080,legXsec_y,"#bf{B^{+}}");


			lat->DrawLatex(legXsec_xLowStart+ShiftX+0.02 + 0.009,legXsec_y + 0.07,"#bf{B_{s}^{0}}"); //Enlarge Label
			lat->DrawLatex(legXsec_xLowStart+ ShiftX + 0.125 + 0.003,legXsec_y+ 0.07,"#bf{B^{+}}");  //Enlarge Label


		//	lat->DrawLatex(legXsec_xLowStart+ShiftX+0.02 - 0.01,legXsec_y + 0.08,"#bf{Rebinned}"); //Enlarge Label + Shift up
		//	lat->DrawLatex(legXsec_xLowStart+ ShiftX + 0.125 + 0.00,legXsec_y  + 0.08,"#bf{Submitted}");  //Enlarge Label + Shift up
			

			cout << "Bs B+ Y Location = " << legXsec_y + 0.08 << endl;


			lat->SetTextSize(ltxSetTextSize2 * 1.3);
			lat->SetTextSize(ltxSetTextSize4); //Enlarge Labels

			
		//	lat->DrawLatex(legXsec_xLowStart-0.15,legXsec_y+0.062-ShiftY,"1.5 < |y| < 2.4");
		//	lat->DrawLatex(legXsec_xLowStart-0.10,legXsec_y+0.017-ShiftY,"|y| < 2.4 ");




//			lat->DrawLatex(legXsec_xLowStart-0.18,legXsec_y+0.062-ShiftY,"1.5 < |y| < 2.4"); //Enlarge Label
//			lat->DrawLatex(legXsec_xLowStart-0.13,legXsec_y+0.017-ShiftY,"|y| < 2.4 "); //Enlarge Label


			lat->DrawLatex(legXsec_xLowStart-0.18-0.02,legXsec_y+0.062-ShiftY + 0.08,"1.5 < |y| < 2.4"); //Enlarge Label + Shift up
			lat->DrawLatex(legXsec_xLowStart-0.13-0.02,legXsec_y+0.017-ShiftY + 0.08,"|y| < 2.4 "); //Enlarge Label + Shift up




//			TLegend *legXSec = new TLegend(legXsec_xLowStart + ShiftX + 0.04,legXsec_y-ShiftY,legXsec_xLowEnd+ShiftX-0.10+0.06,legXsec_y+0.15-ShiftY,"                    ","brNDC");
			TLegend *legXSec = new TLegend(legXsec_xLowStart + ShiftX + 0.04,legXsec_y-ShiftY + 0.08,legXsec_xLowEnd+ShiftX-0.10+0.06,legXsec_y+0.15-ShiftY + 0.08,"                    ","brNDC");
		
			legXSec->SetBorderSize(0);

			cout << "legXsec_xLowStart + ShiftX + 0.04 = " << legXsec_xLowStart + ShiftX + 0.04 << endl;
			cout << "legXsec_y+0.15-ShiftY + 0.08 = " << legXsec_y+0.15-ShiftY + 0.08 << endl;
			cout << "legXsec_y-ShiftY + 0.08 = " << legXsec_y-ShiftY + 0.08 << endl;
			legXSec->SetTextSize(ltxSetTextSize2);
			legXSec->SetLineColor(1);
			legXSec->SetLineStyle(1);
			legXSec->SetLineWidth(1);
			legXSec->SetFillColor(19);
			legXSec->SetFillStyle(0);
			legXSec->SetTextFont(42);
			legXSec->SetNColumns(2);
			legXSec->SetColumnSeparation(0.0);
			//legXSec->AddEntry(pgBs_low,"1.5 < |y| < 2.4","p");
			legXSec->AddEntry(pgBs_low," ","p");
			legXSec->SetTextFont(42);
			legXSec->AddEntry(pgBpl_low," ","p");
			//legXSec->AddEntry(pgBs_high,"|y| < 2.4","p");
			legXSec->AddEntry(pgBs_high," ","p");
			legXSec->SetTextFont(42);
			legXSec->AddEntry(pgBpl_high," ","p");
				


			legXSec->Draw();

		}else{
			//			lat->SetTextSize(ltxSetTextSize1 * 0.64);

			lat->SetTextSize(ltxSetTextSize2*1.7);
			lat->SetTextSize(ltxSetTextSize2*1.7/1.3 *1.17 * ltxSetTextSize4/ltxSetTextSize2); //Enlarge Labels
			lat->SetTextSize(ltxSetTextSize2*1.7/1.3 *1.17 * ltxSetTextSize4/ltxSetTextSize2 /1.05); //Enlarge Labels - Reduce by 8% to match
			

			lat->SetTextFont(42);
			//lat->DrawLatex(ratio_ltxText1_xStart,ratio_ltxText1_yStart,"#bf{#frac{B_{s}^{0}}{B^{+}}}");

			//cout << "ltxSetTextSize1 = " << ltxSetTextSize1 << endl;

			lat->SetTextFont(42);
			lat->SetTextSize(ltxSetTextSize2 * 1.3);
			lat->SetTextSize(ltxSetTextSize4*1.08/1.08); //Enlarge Labels
	
			//lat->DrawLatex(ratio_ltxText2_xStart-0.04,ratio_ltxText2_yStart,"Cent. 0-90%");
	//		lat->DrawLatex(ratio_ltxText2_xStart-0.04,ratio_ltxText2_yStart,"Cent. 0-90%");
	//		lat->DrawLatex(ratio_ltxText2_xStart-0.06,ratio_ltxText2_yStart,"Centrality 0-90%"); //Expand Cent.

			lat->SetTextFont(42);
			lat->SetTextSize(ltxSetTextSize2 * 1.3);
			lat->SetTextSize(ltxSetTextSize4*1.08); //Enlarge Labels
			lat->SetTextSize(ltxSetTextSize4*1.08/1.08); //Enlarge Labels - Reduce 8%
	
		//	lat->DrawLatex(legRatio_xLowStart-0.14,legRatio_y-0.15 - 0.30,Form("global uncert.: #pm %.1f %%",glbSystDown));
		//	lat->DrawLatex(legRatio_xLowStart-0.20,legRatio_y-0.15 - 0.30,Form("global uncertainty: #pm %.1f%%",glbSystDown));

			cout << "legRatio_xLowStart = " << legRatio_xLowStart << "   legRatio_y = " << legRatio_y << endl; 

			// legend
			TLegend *legRatio = new TLegend(legRatio_xLowStart-0.04,legRatio_y+0.08,legRatio_xLowEnd-0.04,legRatio_y+0.18,NULL,"brNDC");
			legRatio->SetBorderSize(0);
			legRatio->SetTextFont(42);
			legRatio->SetTextSize(ltxSetTextSize2 * 1.3);
			legRatio->SetTextSize(ltxSetTextSize4*1.08); //Enlarge Labels
			legRatio->SetTextSize(ltxSetTextSize4*1.08/1.08); //Enlarge Labels - Reduce by 6% to match

			
			legRatio->SetLineColor(1);
			legRatio->SetLineStyle(1);
			legRatio->SetLineWidth(1);
			legRatio->SetFillColor(19);
			legRatio->SetFillStyle(0);
			TLegendEntry *entry1 = legRatio->AddEntry("pgRatio_low","1.5 < |y| < 2.4","p");
			entry1->SetTextFont(42);
			entry1->SetMarkerStyle(markerRatio[0]);
			entry1->SetMarkerSize(1.7);
			entry1->SetMarkerColor(kRed+3);
			
			entry1->SetFillStyle(1001);

			TLegendEntry *entry2 = legRatio->AddEntry("pgRatio_high","|y| < 2.4","p");
			entry2->SetTextFont(42);
			entry2->SetMarkerStyle(markerRatio[1]);
			entry2->SetMarkerColor(kRed+3);
		
			entry2->SetMarkerSize(1.7);
			entry2->SetFillStyle(1001);

			legRatio->Draw();
			//---------------
			if(drawRef)
			{
				TLegend *legRatioRef = new TLegend(legRatioRef_xLowStart-0.06 + 0.09,legRatioRef_y-0.09,legRatioRef_xLowEnd-0.03+ 0.09,legRatio_y+0.11,NULL,"brNDC");
				legRatioRef->SetBorderSize(0);
				legRatioRef->SetTextFont(42);

				legRatioRef->SetTextSize(ltxSetTextSize2*1.3);
				legRatioRef->SetTextSize(ltxSetTextSize4*1.08); //Enlarge Labels
				legRatioRef->SetTextSize(ltxSetTextSize4*1.08/1.08); //Enlarge Labels - Reduced by 8% to Match
		
				legRatioRef->SetLineColor(1);
				legRatioRef->SetLineStyle(1);
				legRatioRef->SetLineWidth(1);
				legRatioRef->SetFillColor(19);
				legRatioRef->SetFillStyle(0);
				//		TLegendEntry *entry1Ref = legRatioRef->AddEntry("FragBand","f_{s}/f_{u} reference: PDG","P");

				/*
				   TLegendEntry *entry1Ref = legRatioRef->AddEntry("FragBand","f_{s}/f_{u} LCHb 13TeV","P");
				   entry1Ref->SetTextFont(42);
				   entry1Ref->SetFillStyle(1001);
				   entry1Ref->SetMarkerStyle(25);
				   entry1Ref->SetMarkerSize(1.4);
				   entry1Ref->SetMarkerColor(kGreen);
				   entry1Ref->SetLineWidth(5);
				   */

				//		TLegendEntry *entry4Ref = legRatioRef->AddEntry(pgRatio_high,"PbPb: CMS 5.02 TeV","p");
				//		entry4Ref->SetTextFont(42);
				//		entry4Ref->SetLineColor(colorRatio[1]);
				//		entry4Ref->SetLineWidth(3);

				TLegendEntry *entry2Ref = legRatioRef->AddEntry("TAMUTheory","PbPb: TAMU","l");
				entry2Ref->SetTextFont(42);
				entry2Ref->SetLineColor(kOrange+1);
				entry2Ref->SetLineWidth(3);


			//	TLegendEntry *entry5Ref = legRatioRef->AddEntry("CAOTheory","PbPb: Cao, Sun, Ko (Cent. 0-80%)","l");
			//	TLegendEntry *entry5Ref = legRatioRef->AddEntry("CAOTheory","PbPb: Langevin (Cent. 0-80%)","l");
				TLegendEntry *entry5Ref = legRatioRef->AddEntry("CAOTheory","PbPb: Langevin (Centrality 0-80%)","l"); //Expand Cent.

				entry5Ref->SetTextFont(42);
				entry5Ref->SetLineColor(kGreen+1);
				entry5Ref->SetLineWidth(3);



				if(drawlhcb){


					TLegendEntry *entry6Ref = legRatioRef->AddEntry("LHCb7TeVRef","pp: LHCb 7 TeV","Pl");
					entry6Ref->SetTextFont(42);
					entry6Ref->SetLineColor(kBlue);
					entry6Ref->SetLineWidth(3);


				}


				legRatioRef->Draw();

			}


		}


		// gPad->RedrawAxis();
		pc1->Update();

		if(bSavePlots)
		{
			if (whichPlot==0)
			{
				pc1->SaveAs(Form("%s/pdf/BRAA_vsPt.pdf",outputDir));
				pc1->SaveAs(Form("%s/png/BRAA_vsPt.png",outputDir));
			}else{
				pc1->SaveAs(Form("%s/pdf/ratio_vsPt_ref%d_%d.pdf",outputDir,drawRef,drawlhcb));
				pc1->SaveAs(Form("%s/png/ratio_vsPt_ref%d_%d.png",outputDir,drawRef,drawlhcb));
			}
		}

		}
