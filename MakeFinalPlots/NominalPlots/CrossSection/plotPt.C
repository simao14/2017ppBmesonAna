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
#include "TPad.h"

#include "CMS_lumi.C"
#include "tdrstyle.C"
#include "DrawLHCb.C"

#include "auxiliaryPt.h"
#include "auxiliaryRef.h"
#include "../../../parameter.h"

//#include "outsideSource/lhcb.C"
#endif
using namespace std;

void makePlot(int mes,
              TString fonFull,
              TString fonFid,
              TString outputDir,
              TGraphAsymmErrors* gLow,
              TGraphAsymmErrors* gSystLow,
              TGraphAsymmErrors* gLowWhite,
              TGraphAsymmErrors* gHigh,
              TGraphAsymmErrors* gSystHigh,
              unsigned int nBins,
              unsigned int nBinsLow,
              double* xLow,
              double* xHigh,
              double* xErrLeftLow,
              double* xErrRightLow,
              double* xErrLeftHigh,
              double* xErrRightHigh,
              double* yLow,
              double* yStatDownLow,
              double* yStatUpLow,
              double* ySystDownLow,
              double* ySystUpLow,
              double* yHigh,
              double* yStatDownHigh,
              double* yStatUpHigh,
              double* ySystDownHigh,
              double* ySystUpHigh,
              double glbSyst);

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
		const char* inputFileType[2] = {"corryield_pt", "RAA_pt"};
		const char* mesonName[nMes]  = {"Bs", "Bp"};

		Int_t endMes = nMes;

    std::vector<int> nlines = {4, 7};
		for (Int_t ib=0; ib<endMes; ib++){
			ifstream in;
			string inputFileName = Form("%s/%s_%s_New.txt",inputDir,inputFileType[whichPlot],mesonName[ib]);
			//if(whichPlot==1) inputFileName = Form("%s/%s.txt",inputDir,inputFileType[whichPlot]);
	  	//	string inputFileName = Form("%s/%s_%s_New.txt",inputDir,inputFileType[whichPlot],mesonName[ib]);		
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
			// while(in >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5] >> x[6] >> x[7] >> x[8] >> x[9])
			// {
			for (auto iline = 0; iline < nlines[ib]; ++iline)
			{
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

		//----------------------------------------------------------------
		// gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
		// Bs
		TGraphAsymmErrors *pgBs_low = new TGraphAsymmErrors(nBinsLow, binLow, bs_low, bs_low_xErrL, bs_low_xErrH, bs_low_yStatL,bs_low_yStatH);
		TGraphAsymmErrors *pgBs_high= new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high, bs_high_xErrL, bs_high_xErrH, bs_high_yStatL,bs_high_yStatH);
		TGraphAsymmErrors *pgBs_lowWhite = new TGraphAsymmErrors(nBinsLow, binLow, bs_low, bs_low_xErrL, bs_low_xErrH, bs_low_yStatL,bs_low_yStatH);   //AddWhite//
		// Bplus
		TGraphAsymmErrors *pgBpl_low = new TGraphAsymmErrors(nBinsLowNew, binLowNew, bpl_low, bpl_low_xErrL, bpl_low_xErrH, bpl_low_yStatL,bpl_low_yStatH);
		TGraphAsymmErrors *pgBpl_lowWhite = new TGraphAsymmErrors(nBinsLowNew, binLowNew, bpl_low, bpl_low_xErrL, bpl_low_xErrH, bpl_low_yStatL,bpl_low_yStatH);
		TGraphAsymmErrors *pgBpl_high= new TGraphAsymmErrors(nBinsHighNew,binHighNew,bpl_high, bpl_high_xErrL, bpl_high_xErrH, bpl_high_yStatL,bpl_high_yStatH);

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
				bpl_high_ySystL, bpl_high_ySystH);

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
		pgBs_syst_high->SetLineColor(kBlue-9);
		pgBs_syst_low->SetLineColor(kBlue-9);
		pgBpl_syst_high->SetLineColor(kGreen-9);
		pgBpl_syst_low->SetLineColor(kGreen-9);


    TString bFonDir = "../../../BsBPFinalResults/Comparisons/Fiducial/FONLLs/";

    makePlot(0,
             bFonDir + "BsFONLL.root",
             bFonDir + "BsFONLLFid.root",
             outputDir,
             pgBs_low, pgBs_syst_low, pgBs_lowWhite,
             pgBs_high, pgBs_syst_high,
             nBinsLow + nBinsHigh,
             nBinsLow,
             binLow, binHigh,
             bs_low_xErrL, bs_low_xErrH,
             bs_high_xErrL, bs_high_xErrH,
             bs_low, bs_low_yStatL, bs_low_yStatH, bs_low_ySystL, bs_low_ySystH,
             bs_high, bs_high_yStatL, bs_high_yStatH, bs_high_ySystL, bs_high_ySystH,
             glbSystUpBs);

    makePlot(1, bFonDir + "fonllOutput_pp_Bplus_5p03TeV_y2p4.root",
             bFonDir + "fonllOutput_pp_Bplus_5p03TeV_yFid.root",
             outputDir,
             pgBpl_low, pgBpl_syst_low, pgBpl_lowWhite,
             pgBpl_high, pgBpl_syst_high,
             nBinsLowNew + nBinsHighNew,
             nBinsLowNew,
             binLowNew, binHighNew,
             bpl_low_xErrL, bpl_low_xErrH,
             bpl_high_xErrL, bpl_high_xErrH,
             bpl_low, bpl_low_yStatL, bpl_low_yStatH, bpl_low_ySystL, bpl_low_ySystH,
             bpl_high, bpl_high_yStatL, bpl_high_yStatH, bpl_high_ySystL, bpl_high_ySystH,
             glbSystUpBp);
    }

// compute ratio and plot them
void makePlot(int mes,
              TString fonFull,
              TString fonFid,
              TString outputDir,
              TGraphAsymmErrors* gLow,
              TGraphAsymmErrors* gSystLow,
              TGraphAsymmErrors* gLowWhite,
              TGraphAsymmErrors* gHigh,
              TGraphAsymmErrors* gSystHigh,
              unsigned int nBins,
              unsigned int nBinsLow,
              double* xLow,
              double* xHigh,
              double* xErrLeftLow,
              double* xErrRightLow,
              double* xErrLeftHigh,
              double* xErrRightHigh,
              double* yLow,
              double* yStatDownLow,
              double* yStatUpLow,
              double* ySystDownLow,
              double* ySystUpLow,
              double* yHigh,
              double* yStatDownHigh,
              double* yStatUpHigh,
              double* ySystDownHigh,
              double* ySystUpHigh,
              double glbSyst) {
  cout << "nbins" << nBins << "\n";

  unsigned int nBinsHigh = nBins - nBinsLow;
  int whichPlot = 0;

  TFile * finFONLLBs = new TFile(fonFull);
  TFile * finFONLLBs2 = new TFile(fonFid);

  TGraphAsymmErrors *BsFONLL = (TGraphAsymmErrors*) finFONLLBs->Get("gaeSigmaBplus");
  BsFONLL->SetLineColor(kRed+2);
  BsFONLL->SetLineWidth(2);
  BsFONLL->SetFillStyle(0);

  finFONLLBs2->cd();
  TGraphAsymmErrors *BsFONLL2 = (TGraphAsymmErrors*) finFONLLBs2->Get("gaeSigmaBplus");
  BsFONLL2->SetLineColor(kRed+2);
  BsFONLL->SetFillColorAlpha(kRed+2, 0.5);

  TGraphAsymmErrors *BsFONLLscaled =  new TGraphAsymmErrors();

  BsFONLLscaled->SetLineColor(kRed-7);
  BsFONLLscaled->SetLineStyle(7);
  BsFONLLscaled->SetLineWidth(2);
  BsFONLLscaled->SetFillStyle(0);
  BsFONLLscaled->SetFillColorAlpha(kRed-7, 0.5);

  double XTempChange;
  double YTempChange;
  double YErrLowTemp;
  double YErrHighTemp;

  for(int i = 0; i < nBinsLow; i ++){
    BsFONLL2->GetPoint(i,XTempChange,YTempChange);
    YErrLowTemp = BsFONLL2->GetErrorYlow(i);
    YErrHighTemp = BsFONLL2->GetErrorYhigh(i);
    BsFONLLscaled->AddPoint(XTempChange,YTempChange);
    BsFONLLscaled->SetPointEYhigh(i,YErrHighTemp);
    BsFONLLscaled->SetPointEYlow(i,YErrLowTemp);
    BsFONLLscaled->SetPointEXhigh(i,BsFONLL->GetErrorXhigh(i));
    BsFONLLscaled->SetPointEXlow(i,BsFONLL->GetErrorXlow(i));
  }

  // Get ratio plots
  vector<double> BsXsec;
  vector<double> BsXsecStat;
  vector<double> BsXsecSyst;
  vector<double> FONLL;
  vector<double> FONLLUp;
  vector<double> FONLLDown;

  double XTempFONLL;
  double YTempFONLL;
  for (auto i = 0; i < nBinsLow; ++i) {
    BsXsec.push_back(yLow[i]);
    BsXsecStat.push_back(yStatDownLow[i]);
    BsXsecSyst.push_back(ySystDownLow[i]);
  }
  for (auto i = 0; i < nBinsHigh; ++i) {
    BsXsec.push_back(yHigh[i]);
    BsXsecStat.push_back(yStatDownHigh[i]);
    BsXsecSyst.push_back(ySystDownHigh[i]);
  }

  double RatioBs[nBins];
  double RatioBsStat[nBins];
  double RatioBsSyst[nBins];
  double RatioBsFonErrHigh[nBins];
  double RatioBsFonErrLow[nBins];


  std::vector<double> Unity(nBins, 1);

for (auto i = 0; i < nBins; ++i) {
  if(i>=nBinsLow){
    BsFONLL->GetPoint(i, XTempFONLL, YTempFONLL);
    FONLL.push_back(YTempFONLL);
    FONLLUp.push_back(BsFONLL->GetErrorYhigh(i));
    FONLLDown.push_back(BsFONLL->GetErrorYlow(i));
  } else {
    BsFONLLscaled->GetPoint(i, XTempFONLL, YTempFONLL);
    FONLL.push_back(YTempFONLL);
    FONLLUp.push_back(BsFONLLscaled->GetErrorYhigh(i));
    FONLLDown.push_back(BsFONLLscaled->GetErrorYlow(i));
  }
    cout << "bs xsex:" << BsXsec[i] << "\n";
    cout << "FONLL:" << FONLL[i] << "\n";
    RatioBs[i] = BsXsec[i] / FONLL[i];
    RatioBsStat[i] = BsXsecStat[i] / FONLL[i];
    RatioBsSyst[i] = BsXsecSyst[i] / FONLL[i];
    cout << "ratio:" << RatioBs[i] << "\n";
    RatioBsFonErrHigh[i] = FONLLUp[i] / FONLL[i];
    RatioBsFonErrLow[i] = FONLLDown[i] / FONLL[i];
  }

for (int i=0; i<nBinsLow; ++i){BsFONLL->RemovePoint(0);} 

  TGraphAsymmErrors gRatioBs_low(nBinsLow,
                                 xLow,
                                 RatioBs,
                                 xErrLeftLow, xErrRightLow,
                                 RatioBsStat, RatioBsStat);

  TGraphAsymmErrors gRatioBs_high(nBinsHigh,
                                  xHigh,
                                  RatioBs + nBinsLow,
                                  xErrLeftHigh, xErrRightHigh,
                                  RatioBsStat + nBinsLow, RatioBsStat + nBinsLow);

  TGraphAsymmErrors gRatioBs_syst_low(nBinsLow,
                                      xLow,
                                      RatioBs,
                                      xErrLeftLow, xErrRightLow,
                                      RatioBsSyst, RatioBsSyst);

  TGraphAsymmErrors gRatioBs_syst_high(nBinsHigh,
                                       xHigh,
                                       RatioBs + nBinsLow,
                                       xErrLeftHigh, xErrRightHigh,
                                       RatioBsSyst + nBinsLow, RatioBsSyst + nBinsLow);

  TGraphAsymmErrors gRatioBs_Fon_low(nBinsLow,
                                     xLow,
                                     Unity.data(),
                                     xErrLeftLow, xErrRightLow,
                                     RatioBsFonErrLow, RatioBsFonErrHigh);

  TGraphAsymmErrors gRatioBs_Fon_high(nBinsHigh,
                                      xHigh,
                                      Unity.data(),
                                      xErrLeftHigh, xErrRightHigh,
                                      RatioBsFonErrLow + nBinsLow,
                                      RatioBsFonErrHigh + nBinsLow);

  gRatioBs_low.SetMarkerStyle(markerLow[mes]);
  gRatioBs_high.SetMarkerStyle(markerHigh[mes]);

  int hcolor = kBlue - 9;
  double halpha = 0.5;
  if (mes == 1) {
    hcolor = kGreen - 9;
  }
  gRatioBs_syst_low.SetFillColorAlpha(hcolor, halpha);
  gRatioBs_syst_high.SetFillColorAlpha(hcolor, halpha);
  gRatioBs_Fon_low.SetLineColor(kRed-7);
  gRatioBs_Fon_low.SetLineStyle(7);
  gRatioBs_Fon_high.SetLineColor(kRed+2);
  gRatioBs_Fon_low.SetLineWidth(2);
  gRatioBs_Fon_high.SetLineWidth(2);

  // gRatioBs_Fon_low.SetFillColorAlpha(kRed+2, 0.5);
  // gRatioBs_Fon_high.SetFillColorAlpha(kRed+2, 0.5);



  //---------------- general stuff
  TLatex *lat = new TLatex();
  lat->SetNDC();

  // // ##################################################### x-sec canvas
  TCanvas *pc1 = new TCanvas("pc1","pc1");
  //pc1->SetBottomMargin(0.10);
  
  //	cout << "pc1->GetBottomMargin() = " << pc1->GetBottomMargin() << endl;
  
  pc1->SetBottomMargin(0.135); //Enlarge Margin

  //Added to fix//


  // TCanvas * cRatio = new TCanvas("cRatio","cRatio",800, 1200);
  TPad * dataPad = new TPad("MyPad1","",0,0.3,1, 0.96);
  dataPad->SetBottomMargin(0);
  dataPad->SetLogy();
  dataPad->Draw();

  TPad * ratioPad = new TPad("MyPad2","",0,0.0,1,0.3);
  ratioPad->SetTopMargin(0);
  ratioPad->SetBottomMargin(0.3);
  ratioPad->Draw();

  //supplemental info on plot:
  lat->SetTextFont(42);
  lat->SetTextSize(ltxSetTextSize2 * 1.3);
  lat->SetTextSize(ltxSetTextSize4); //Enlarge Labels

  cout << "ltxSetTextSize4 = " << ltxSetTextSize4 << endl;

  dataPad->cd();
  double xmin = 5;
  double xmax = 60;
  if (mes == 0) {
    xmin = 7;
    xmax = 50;
  }
  //-------------------------------------------
  TF1 *f4 = new TF1("f4","1", xmin, xmax);
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
  if(whichPlot==1) f4->GetYaxis()->SetRangeUser(0.0,0.90);
  if(whichPlot==0) f4->GetYaxis()->SetRangeUser(100.0,2000000); //1.8 -> 0.9

  //f4->GetXaxis()->SetNdivisions(-6);
  f4->Draw();// axis

  TH2D * HisEmpty = new TH2D("HisEmpty","",100, xmin, xmax, 100, 100.0, 2000000);
  // HisEmpty->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
  // HisEmpty->GetYaxis()->SetTitle("d#sigma/dp_{T} (pb c/GeV)");
  // HisEmpty->GetXaxis()->CenterTitle();
  // HisEmpty->GetYaxis()->CenterTitle();
  // HisEmpty->GetYaxis()->SetTitleOffset(1.8);
  // HisEmpty->GetXaxis()->SetTitleOffset(1.3);		
  // HisEmpty->Draw();

  // TH1D test0("test0", "", 10, 0, 10);
  // test0.Fill(2);
  // test0.Draw();
  BsFONLL->Draw("5");
  BsFONLLscaled->Draw("5");
  gSystHigh->Draw("5");
  gHigh->Draw("P");
  gSystLow->Draw("5");
  gLow->Draw("P");
  gLowWhite->Draw("P");
  if (mes == 0) {lat->DrawLatex(xsec_ltxText1_xStart + 0.02,xsec_ltxText1_yStart-0.65+0.037,Form("B_{s}^{0} global uncertainty: #pm %.1f%%",glbSystUpBs));
  } else { lat->DrawLatex(xsec_ltxText1_xStart + 0.02,xsec_ltxText1_yStart-0.70 + 0.037,Form("B^{+} global uncertainty: #pm %.1f%%",glbSystUpBp)); }

  double ShiftX = 0.05;
  double ShiftY = 0.13;
  lat->SetTextSize(0.05);
  lat->SetTextSize(0.05 * ltxSetTextSize4/ltxSetTextSize2);  //Enlarge Labels
      
  lat->SetTextSize(0.048 * 1.15);  //Enlarge Labels
  if (mes == 0) {
    lat->DrawLatex(legXsec_xLowStart+ShiftX+0.02 + 0.009,legXsec_y + 0.07,"#bf{B_{s}^{0}}"); //Enlarge Label
  } else {
    // lat->DrawLatex(legXsec_xLowStart+ ShiftX + 0.125 + 0.003,legXsec_y+ 0.07,"#bf{B^{+}}");  //Enlarge Label
    lat->DrawLatex(legXsec_xLowStart+ ShiftX + 0.02 + 0.009,legXsec_y+ 0.07,"#bf{B^{+}}");  //Enlarge Label
  }

  cout << "Bs B+ Y Location = " << legXsec_y + 0.08 << endl;

  lat->SetTextSize(ltxSetTextSize2 * 1.3);
  lat->SetTextSize(ltxSetTextSize4); //Enlarge Labels
  lat->DrawLatex(legXsec_xLowStart-0.18-0.02,legXsec_y+0.062-ShiftY + 0.08,"1.5 < |y| < 2.4"); //Enlarge Label + Shift up
  lat->DrawLatex(legXsec_xLowStart-0.13-0.02,legXsec_y+0.017-ShiftY + 0.08,"|y| < 2.4 "); //Enlarge Label + Shift up

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
  legXSec->AddEntry(gLow,"Data","p");
  // legXSec->SetTextFont(42);
  // legXSec->AddEntry(pgBpl_low," ","p");
  legXSec->AddEntry(BsFONLL,"FONLL","f");
  // legXSec->AddEntry(pgBpl_low," ","");
  //legXSec->AddEntry(pgBs_high,"|y| < 2.4","p");
  legXSec->AddEntry(gHigh," ","p");
  legXSec->SetTextFont(42);
  //legXSec->AddEntry(gHigh," ","");
  legXSec->AddEntry(BsFONLLscaled,"FONLL(|y|<1.5)","f");
  legXSec->Draw("same");

  ratioPad->cd();
  // TH2D * HisEmpty4 = new TH2D("HisEmpty4","",100, xmin, xmax, 100, 0.5, 1.5);
  // HisEmpty4->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
  // HisEmpty4->GetYaxis()->SetTitle("Data/FONLL");
  // HisEmpty4->GetXaxis()->CenterTitle();
  // HisEmpty4->GetYaxis()->CenterTitle();
  // HisEmpty4->GetYaxis()->SetTitleOffset(0);
  // HisEmpty4->GetYaxis()->SetTitleSize(0.1);
  // HisEmpty4->GetYaxis()->SetLabelSize(0.1);
  // HisEmpty4->GetXaxis()->SetTitleSize(0.1);
  // HisEmpty4->GetXaxis()->SetLabelSize(0.1);

  double rRange = 0.55;
  TF1 fRatio("fRatio","1", xmin, xmax);
  fRatio.GetYaxis()->SetRangeUser(1 - rRange, 1 + rRange);
  fRatio.GetYaxis()->SetTitle("Data / FONLL");
  fRatio.GetXaxis()->SetTitle(xAxName[0]);
  fRatio.GetXaxis()->CenterTitle();
  fRatio.GetYaxis()->CenterTitle();
  fRatio.GetXaxis()->SetTitleSize(0.1);  //Unify Textsize
  fRatio.GetXaxis()->SetTitleOffset(1.2);
  fRatio.GetYaxis()->SetTitleSize(0.1);
  fRatio.GetYaxis()->SetTitleOffset(0.68);
  fRatio.GetYaxis()->SetLabelSize(0.10);
  fRatio.GetXaxis()->SetLabelSize(0.10);

  // HisEmpty4->Draw();
  fRatio.Draw();

  gRatioBs_Fon_low.Draw("5");
  gRatioBs_Fon_high.Draw("5");
  gRatioBs_syst_low.Draw("5");
  gRatioBs_syst_high.Draw("5");
  gRatioBs_low.Draw("ep");
  gRatioBs_high.Draw("ep");
  

  auto gRatio_white = (TGraphAsymmErrors*) gRatioBs_low.Clone();
  gRatio_white->SetMarkerSize(markerSizeRatio[0]);
  gRatio_white->SetMarkerSize(markerSizeLow[1]*0.86);
  gRatio_white->SetMarkerColor(kWhite);
  // gRatio_white->Draw("P");

  TLine * UnityLine = new TLine(xmin, 1, xmax, 1);
  UnityLine->SetLineWidth(2);
  UnityLine->SetLineStyle(2);
  UnityLine->SetLineColor(1);
  UnityLine->Draw("SAME");

  dataPad->Update();
  ratioPad->Update();

  CMS_lumi(pc1,19011,0);

  pc1->Update();


  if (mes == 0) {
    pc1->SaveAs(Form("%s/pdf/xsec_vsPt_Bs.pdf", outputDir.Data()));
    pc1->SaveAs(Form("%s/png/xsec_vsPt_Bs.png", outputDir.Data()));
  } else {
    pc1->SaveAs(Form("%s/pdf/xsec_vsPt_BP.pdf", outputDir.Data()));
    pc1->SaveAs(Form("%s/png/xsec_vsPt_BP.png", outputDir.Data()));
  }

}
