#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "../henri2022/parameter.h" 

using namespace std;
using std::cout;
using std::endl;

void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, 
		std::vector<std::vector<double> > numbers , std::string caption){
	
	std::ofstream file_check;
	std::ofstream file;

	//Begin Document
                                                                                    
	file.open(filename + ".tex");
	file_check.open(filename + "_check.tex");


	file_check << "\\documentclass{article}" << std::endl;
	//file << "\\usepackage[utf8]{inputenc}" << std::endl;     
	file_check << "\\usepackage{rotating}" << std::endl;                                                                                   
	// file_check << "\\usepackage{cancel}" << std::endl;
	file_check << "\\usepackage{geometry}" << std::endl;
	file_check << "\\usepackage{booktabs}" << std::endl;
	file_check << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;

	file_check << "\\begin{document}" << std::endl;
	// Create table                                                                                                                                                                                                                                                
	std::string col="c";
	for(int i=1; i<n_col; i++)
		col+="|c";

		file_check << "\\begin{sidewaystable}"<< std::endl;
		file_check << "\\begin{tabular}{"+col+"}" << std::endl;
		file_check << "\\toprule" << std::endl;
		file << "\\begin{tabular}{"+col+"}" << std::endl;
		file << "\\toprule" << std::endl;
	
	for(int c=0; c<n_col-1; c++){
		file << col_name[c] << " & ";
		file_check << col_name[c] << " & ";
	}
	
	file << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;
	file_check << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;

	for(int i=1; i<n_lin; i++)
	{	
		file << labels[i-1] << " & ";
		file_check << labels[i-1] << " & ";
		for(int c=1; c<n_col-1; c++){
			file << std::setprecision(3)<<  numbers[c-1][i-1]<< " & ";
			file_check << std::setprecision(3)<<  numbers[c-1][i-1]<< " & ";
									}
			file << std::setprecision(3)<<  numbers[n_col-2][i-1]<< "  \\\\" << std::endl;
			file_check << std::setprecision(3)<<  numbers[n_col-2][i-1]<< "  \\\\" << std::endl; 
	}
	
	file << "\\bottomrule" << std::endl;
	file_check << "\\bottomrule" << std::endl;

	//End Table                                                                                                                                    
	file << "\\end{tabular}" << std::endl;
	file_check << "\\end{tabular}" << std::endl;
	file_check << "\\caption{"+caption+"}" << std::endl;
	file_check << "\\end{sidewaystable}"<< std::endl;

	//file_check << "\\end{table}" << std::endl;
	//End document                                                                                                                                 

	file_check << "\\end{document}" << std::endl;

	file.close();
	file_check.close();
	
	system(("pdflatex " + filename+ "_check.tex").c_str());
	system(("open " + filename + "_check.pdf").c_str());
}






















void PlotEffSyst2D(TString meson_n, TString whichvar, int BsBP=0, int usemc=0){

	TString var_l;
	int NBins = 7;
	TString bsbpbins = "";
	if (BsBP==1){bsbpbins = "_BsBPBINS";}

	if (whichvar =="pt" ){
		if(meson_n == "BP" && BsBP==0){
			NBins = nptBinsBP;
			var_l="p_{T} [GeV/c]";
		
		} else if(meson_n == "Bs" || BsBP==1){
			NBins = nptBins;
			var_l="p_{T} [GeV/c]";
		}

	} else if(whichvar =="y"){
		NBins = nyBins_both;
		var_l="Rapidity";
	} else if(whichvar =="Mult"){
		NBins = nmBins_both;
		var_l="Multiplicity";
	}

	double ptBins[NBins+1];
	if (whichvar =="y"){
	for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ybinsvec[i];             
		}
	} 
	if (whichvar =="Mult"){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  nmbinsvec[i];             
		}
	}
	if (whichvar =="pt" && meson_n == "BP" && BsBP==0){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvecBP[i];             
		}
	}
	if (whichvar =="pt" && (meson_n=="Bs"||BsBP==1)){
		for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvec[i];             
		}
	}

	TCanvas * c = new TCanvas("c","c",700,700);
	c->cd();

	gStyle->SetOptStat(0);

	TString infile;
  	TString outFile;

	if(meson_n == "BP" && usemc==0) {
    infile =  Form("OutFiles/BPSyst2D_%s%s.root",whichvar.Data(),bsbpbins.Data());
    outFile =  Form("OutFiles/BPError2D_%s%s.root",whichvar.Data(),bsbpbins.Data());
  } else if(meson_n == "Bs" && usemc==0) {
    infile =  Form("OutFiles/BsSyst2D_%s%s.root",whichvar.Data(),bsbpbins.Data());
    outFile =  Form("OutFiles/BsError2D_%s%s.root",whichvar.Data(),bsbpbins.Data());
  } else if(meson_n == "BP" && usemc==1) {
    infile =  Form("OutFiles/BPSyst2D_%s%s_MC.root",whichvar.Data(),bsbpbins.Data());
    outFile =  Form("OutFiles/BPError2D_%s%s_MC.root",whichvar.Data(),bsbpbins.Data());
  } else if(meson_n == "Bs" && usemc==1) {
    infile =  Form("OutFiles/BsSyst2D_%s%s_MC.root",whichvar.Data(),bsbpbins.Data());
    outFile =  Form("OutFiles/BsError2D_%s%s_MC.root",whichvar.Data(),bsbpbins.Data());
  }

	gSystem->mkdir("EffSystPlots", true);
	gSystem->mkdir(Form("EffSystPlots/%s",meson_n.Data()), true);

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TH1D * Eff1DHis = (TH1D * ) fin->Get("Eff2DHis");
	TH1D * Eff1DHisTnPUp = (TH1D * ) fin->Get("Eff2DTnPUpSystHis");
	TH1D * Eff1DHisTnPDown = (TH1D * ) fin->Get("Eff2DTnPDownSystHis");
	TH1D * Eff1DHisBpt = (TH1D * ) fin->Get("Eff2DBptHis");
	TH1D * Eff1DHisBDT = (TH1D * ) fin->Get("Eff2DBDTHis");

  TFile fout(outFile, "recreate");

	//Draw Eff2DHis
	TCanvas * EFFplot  = new TCanvas("cSyst","cSyst",600,600);
	EFFplot->cd();
	Eff1DHis->SetMarkerStyle(20);
	Eff1DHis->SetMarkerSize(1);
	Eff1DHis->SetMarkerColor(kBlack);
	Eff1DHis->SetLineColor(kBlack);
	Eff1DHis->Draw("ep");

	if (usemc==0){EFFplot->SaveAs(Form("EffSystPlots/%s/%s_EFFplot%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
	else {EFFplot->SaveAs(Form("EffSystPlots/%s/%s_EFFplot_MC%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}

	//Draw Systematic Uncertainties
	TCanvas * cSyst  = new TCanvas("cSyst","cSyst",600,600);
	cSyst->cd();
	Eff1DHis->SetMarkerStyle(20);
	Eff1DHis->SetMarkerSize(1);
	Eff1DHis->SetMarkerColor(kBlack);
	Eff1DHis->SetLineColor(kBlack);
	Eff1DHisTnPUp->SetMarkerStyle(20);
	Eff1DHisTnPUp->SetMarkerSize(1);
	Eff1DHisTnPUp->SetMarkerColor(kRed);
	Eff1DHisTnPUp->SetLineColor(kRed);
	Eff1DHisTnPDown->SetMarkerStyle(20);
	Eff1DHisTnPDown->SetMarkerSize(1);
	Eff1DHisTnPDown->SetMarkerColor(kBlue);
	Eff1DHisTnPDown->SetLineColor(kBlue);
	Eff1DHisTnPUp->Draw("ep");
	Eff1DHis->Draw("epSAME");
	Eff1DHisTnPDown->Draw("epSAME");

	TLegend* leg = new TLegend(0.62,0.77,0.89,0.89,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.025);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(Eff1DHis,"Nominal","PL");
	leg->AddEntry(Eff1DHisTnPUp,"TnP Variation Up","PL");
	leg->AddEntry(Eff1DHisTnPDown,"TnP Variation Down","PL");
	leg->Draw("same");
	if (usemc==0){cSyst->SaveAs(Form("EffSystPlots/%s/%s_TnPSystComp%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
	else {cSyst->SaveAs(Form("EffSystPlots/%s/%s_TnPSystComp_MC%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
	Eff1DHisBDT->SetMarkerStyle(20);
	Eff1DHisBDT->SetMarkerSize(1);
	Eff1DHisBDT->SetMarkerColor(kRed);
	Eff1DHisBDT->SetLineColor(kRed);
	Eff1DHis->Draw("ep");
	Eff1DHisBDT->Draw("epSAME");

	TLegend* leg2 = new TLegend(0.7,0.8,0.89,0.89,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.025);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->AddEntry(Eff1DHis,"Nominal","PL");
	leg2->AddEntry(Eff1DHisBDT,"BDT Weighted","PL");
	leg2->Draw("same");
	if (usemc==0) {cSyst->SaveAs(Form("EffSystPlots/%s/%s_MCDataSystComp%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
	else {cSyst->SaveAs(Form("EffSystPlots/%s/%s_MCDataSystComp_MC%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
	Eff1DHisBpt->SetMarkerStyle(20);
	Eff1DHisBpt->SetMarkerSize(1);
	Eff1DHisBpt->SetMarkerColor(kRed);
	Eff1DHisBpt->SetLineColor(kRed);
	Eff1DHis->Draw("ep");
	Eff1DHisBpt->Draw("epSAME");

	TLegend* leg3 = new TLegend(0.7,0.8,0.89,0.89,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.025);
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->AddEntry(Eff1DHis,"Nominal","PL");
	leg3->AddEntry(Eff1DHisBpt,"Bpt Weighted","PL");
	leg3->Draw("same");
	if (usemc==0){cSyst->SaveAs(Form("EffSystPlots/%s/%s_BptSystComp%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
	else {cSyst->SaveAs(Form("EffSystPlots/%s/%s_BptSystComp_MC%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
	//Done drawing  (only draw for BsBP=0)


	TH1D * TnPSyst = (TH1D *) Eff1DHisTnPUp->Clone("TnPSyst");
	TnPSyst->GetYaxis()->SetTitle("TnP Systematic Uncertainties (%)");
	TnPSyst->GetYaxis()->SetTitleOffset(1.3);
	TnPSyst->SetLineColor(kBlack);
	TnPSyst->SetMarkerColor(kBlack);	
	TnPSyst->Reset();

	TH1D * BptSyst = (TH1D *) Eff1DHisBpt->Clone("BptSyst");
	BptSyst->GetYaxis()->SetTitle(Form("B-meson %s Shape Systematic Uncertainties (%s)",var_l.Data(),"%"));
	BptSyst->GetYaxis()->SetTitleOffset(1.3);
	BptSyst->SetLineColor(kBlack);
	BptSyst->SetMarkerColor(kBlack);	
	BptSyst->Reset();

	TH1D * BDTSyst = (TH1D *) Eff1DHisBDT->Clone("BDTSyst");
	BDTSyst->GetYaxis()->SetTitle("MC/Data Discrepancy Systematic Uncertainties (%)");
	BDTSyst->GetYaxis()->SetTitleOffset(1.3);
	BDTSyst->SetLineColor(kBlack);
	BDTSyst->SetMarkerColor(kBlack);	
	BDTSyst->Reset();

	c->cd();

	float SystValue;
	float SystValueError;

	for(int i = 0; i < Eff1DHisTnPUp->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisTnPUp->GetBinContent(i+1) -  Eff1DHis->GetBinContent(i+1) )/Eff1DHis->GetBinContent(i+1) * 100;
		cout << "TnP Syst: " << SystValue << endl;
		TnPSyst->SetBinContent(i+1,SystValue);
		TnPSyst->SetBinError(i+1,0.001);
	}

	for(int i = 0; i < Eff1DHisTnPUp->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBpt->GetBinContent(i+1) -  Eff1DHis->GetBinContent(i+1) )/Eff1DHis->GetBinContent(i+1) * 100;
		BptSyst->SetBinContent(i+1,SystValue);
		BptSyst->SetBinError(i+1,0.001);
		cout << "Bpt Syst: " << SystValue << endl;
	}

	for(int i = 0; i < Eff1DHisTnPUp->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBDT->GetBinContent(i+1) -  Eff1DHis->GetBinContent(i+1) )/Eff1DHis->GetBinContent(i+1) * 100;
		BDTSyst->SetBinContent(i+1,SystValue);
		BDTSyst->SetBinError(i+1,0.001);
		cout << "BDT Syst: " << SystValue << endl;
    	cout << SystValue << "\n";

	}

	TnPSyst->Draw("ep");
		if (usemc==0){c->SaveAs(Form("EffSystPlots/%s/%s_TnPSystRatio%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
		else {c->SaveAs(Form("EffSystPlots/%s/%s_TnPSystRatio_MC%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
		BptSyst->Draw("ep");
		if (usemc==0){c->SaveAs(Form("EffSystPlots/%s/%s_BptSysRatio%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
		else {c->SaveAs(Form("EffSystPlots/%s/%s_BptSysRatio_MC%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
		BDTSyst->Draw("ep");
		if (usemc==0){c->SaveAs(Form("EffSystPlots/%s/%s_MCDataSystRatio%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}
		else {c->SaveAs(Form("EffSystPlots/%s/%s_MCDataSystRatio_MC%s.pdf",meson_n.Data(),whichvar.Data(),bsbpbins.Data()));}

  TnPSyst->Write();
  BptSyst->Write();
  BDTSyst->Write();

	float BP1DEffX[NBins];
	float B1DEffXErrUp[NBins] ;
	float BP1DEffY[NBins];
	float BP1DEffYErrUp[NBins];

	for( int c=0; c <NBins; c++){ 
		BP1DEffX[c]=(ptBins[c]+ptBins[c+1])/2;
		B1DEffXErrUp[c]=(abs(ptBins[c+1]-ptBins[c]))/2;
		BP1DEffY[c] = Eff1DHis->GetBinContent(c+1);
		BP1DEffYErrUp[c] = Eff1DHis->GetBinError(c+1);
	}

	float BP1DEffYSystUp[NBins];
	float BP1DEffYSystDown[NBins];

	float BP2DEffYSystUp[NBins];
	float BP2DEffYSystDown[NBins];

	double B_nu;
  	if (meson_n == "BP"){ B_nu = 2.4 ;}
  	else { B_nu = 4.8;}
	double BPTrackingSyst[NBins];
	for( int c=0; c < NBins; c++){BPTrackingSyst[c]= B_nu ;}
	float BPPDFSyst[NBins];
	
	float BPMCDataSyst1D[NBins];
	float BPPtShapeSyst1D[NBins];
	float BPTnPSystDown1D[NBins];
	float BPTnPSystUp1D[NBins];

  // Get systematics from input files
    for (auto ibin = 0; ibin < NBins; ++ibin){
		BPMCDataSyst1D[ibin] = BDTSyst->GetBinContent(ibin + 1);
		BPPtShapeSyst1D[ibin] = BptSyst->GetBinContent(ibin + 1);
		BPTnPSystDown1D[ibin] = TnPSyst->GetBinContent(ibin + 1);
		BPTnPSystUp1D[ibin] = BPTnPSystDown1D[ibin];
  	}

  // RMS of all the errors
	float BPTotalSystDownRatio1D[NBins];
	float BPTotalSystUpRatio1D[NBins];


	for(int i = 0; i < NBins; i++){

		BPTotalSystDownRatio1D[i] = TMath::Sqrt(TMath::Power(BPMCDataSyst1D[i], 2) + TMath::Power(BPPtShapeSyst1D[i], 2) + TMath::Power(BPTnPSystDown1D[i], 2)) / 100;
    	BPTotalSystUpRatio1D[i] = TMath::Sqrt(TMath::Power(BPMCDataSyst1D[i], 2) + TMath::Power(BPPtShapeSyst1D[i], 2) + TMath::Power(BPTnPSystUp1D[i], 2)) / 100;

		//BPTotalSystDownRatio1D[i] = TMath::Sqrt(TMath::Power(TMath::Power(BPTrackingSyst[i], 2) + BPMCDataSyst1D[i], 2) + TMath::Power(BPTrackSelSystMC[i], 2) + TMath::Power(BPPtShapeSyst1D[i], 2) + TMath::Power(BPTnPSystDown1D[i], 2)) / 100;
    	//BPTotalSystUpRatio1D[i] = TMath::Sqrt(TMath::Power(TMath::Power(BPTrackingSyst[i], 2) + BPMCDataSyst1D[i], 2) + TMath::Power(BPTrackSelSystMC[i], 2) + TMath::Power(BPPtShapeSyst1D[i], 2) + TMath::Power(BPTnPSystUp1D[i], 2)) / 100;
	
	}

	for(int i = 0; i < NBins; i++){
		BP1DEffYSystUp[i] = BP1DEffY[i] * BPTotalSystUpRatio1D[i];
		BP1DEffYSystDown[i] = BP1DEffY[i] * BPTotalSystDownRatio1D[i];

	}

	TGraphAsymmErrors *BP1DEffGraph = new TGraphAsymmErrors(NBins, BP1DEffX ,BP1DEffY , B1DEffXErrUp, B1DEffXErrUp, BP1DEffYErrUp, BP1DEffYErrUp); 
	TGraphAsymmErrors *BP1DEffGraphSyst = new TGraphAsymmErrors(NBins, BP1DEffX , BP1DEffY, B1DEffXErrUp, B1DEffXErrUp, BP1DEffYSystDown, BP1DEffYSystUp);

	
	gStyle->SetOptStat(0);
	cSyst->cd(); 
	cSyst->SetLeftMargin(0.15);

	//Setup histograms for different purposs
	TH2D * HisEmpty;
	if(meson_n == "BP" && whichvar == "pt") {HisEmpty = new TH2D("HisEmpty","",100,5,60,100,0,60);} 
	if(meson_n == "Bs" && whichvar == "pt") {HisEmpty = new TH2D("HisEmpty","",100,7,50,100,0,60);}
	if(whichvar == "y") {HisEmpty = new TH2D("HisEmpty","",100,0,2.4,100,0,50);}
	if(meson_n == "BP" && whichvar == "Mult") {HisEmpty = new TH2D("HisEmpty","",100,0,100,100,0,60);}   // need to adjust range for when we have nmult results
	if(meson_n == "Bs" && whichvar == "Mult") {HisEmpty = new TH2D("HisEmpty","",100,0,100,100,0,60);}

	HisEmpty->GetXaxis()->SetTitle(var_l.Data());
	HisEmpty->GetYaxis()->SetTitle("1/(#alpha #times #epsilon)");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();

	BP1DEffGraph->SetMarkerStyle(20);
	BP1DEffGraph->SetMarkerSize(1);
	BP1DEffGraph ->SetMarkerColor(kGreen+2);
	BP1DEffGraph ->SetLineColor(kGreen+2);
	BP1DEffGraphSyst ->SetFillColorAlpha(kGreen-7,0.5);
	BP1DEffGraphSyst ->SetLineColor(kGreen-7);

	
	
	HisEmpty->Draw();
	BP1DEffGraphSyst->Draw("5same");
	BP1DEffGraph->Draw("epSAME");
	cSyst->SaveAs(Form("EffSystPlots/%s/Effcomp_%s%s.pdf", meson_n.Data(), whichvar.Data(),bsbpbins.Data()));

	
	


	string name;
	TString whichvarname;
	if(whichvar == "pt"){name="$<p_T<$"; whichvarname="pt";} 
	else if(whichvar == "y"){name="$<y<$";whichvarname="y";} 
	else if(whichvar == "Mult"){name="$<nTrks<$";whichvarname="nMult";}
	std::vector<std::string> col_name;
	std::vector<std::string> labels = {"","Eff value","Stat error","Syst error"};
	std::vector<std::vector<double>> eff_values;
	std::vector<double> val1D;
	std::vector<double> stat1D;
	std::vector<double> syst1D;
	for(int i=0;i<NBins;i++){

		std::ostringstream clabel;
		clabel<<ptBins[i]<<name<<ptBins[i+1];
		std::string label1 = clabel.str();
		col_name.push_back(label1);

		val1D.push_back(BP1DEffY[i]);
		stat1D.push_back(BP1DEffYErrUp[i]/BP1DEffY[i]);
		syst1D.push_back(BP1DEffYSystDown[i]/BP1DEffY[i]);

	}

	eff_values.push_back(val1D);
	eff_values.push_back(stat1D);
	eff_values.push_back(syst1D);


	gSystem->mkdir("Trash",true);

	
	latex_table(Form("2DEfftable_%s%s",whichvarname.Data(),bsbpbins.Data()), 4,  NBins+1,  labels , col_name , eff_values, "Efficiency values and errors");

	std::vector<std::string> filetype ={"_check.aux", "_check.log",".tex","_check.tex"};
	for (int j=0;j<(int)(filetype.size());j++){
				rename(("2DEfftable_"+ std::string (whichvarname.Data()) + std::string (bsbpbins.Data()) +filetype[j]).c_str(),("Trash/2DEfftable_"+std::string (whichvarname.Data())  + std::string (bsbpbins.Data()) +filetype[j]).c_str());
			}
	rename(("2DEfftable_"+ std::string (whichvarname.Data()) + std::string (bsbpbins.Data()) +"_check.pdf").c_str(),("EffSystPlots/"+std::string (meson_n.Data())+"/2DEfftable_"+std::string (whichvarname.Data())  + std::string (bsbpbins.Data()) +"_check.pdf").c_str());
	


}


