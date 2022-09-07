#include <cmath> 

void yield_ratio(int syst,TString varExp){

	


	TFile *y1= new TFile(Form("./results/ntphi_%s_ratio.root",varExp.Data()),"read");
	TFile *y2= new TFile(Form("./results/ntKp_%s_ratio.root",varExp.Data()),"read");

	TMultiGraph *TMG1= (TMultiGraph*) y1->Get("TG");
	TMultiGraph *TMG2= (TMultiGraph*) y2->Get("TG");

	auto tl1 = TMG1->GetListOfGraphs();
	auto tl2 = TMG2->GetListOfGraphs();


	auto TG1 = (TGraphAsymmErrors *) tl1->At(0);
  	auto TG2 = (TGraphAsymmErrors *) tl2->At(0);
	
	const int _nBins = (TG1->GetN());
	TCanvas* c=new TCanvas();
	TLegend* leg_ratio=new TLegend(0.7,0.7,0.9,0.9);
	TMultiGraph* m_ratio=new TMultiGraph();
	double x[_nBins];
	double y[_nBins];
	double ey[_nBins];
	double ex_l[_nBins];
	double ex_h[_nBins];

	
	for (int i=0;i<_nBins;i++){
		y[i]=TG1->GetY()[i]/TG2->GetY()[i];
		x[i]=(TG1->GetX()[i]+TG2->GetX()[i])/2;
		ey[i]=abs(y[i])*sqrt(pow(TG1->GetErrorY(i)/TG1->GetY()[i],2)+pow(TG2->GetErrorY(i)/TG2->GetY()[i],2));
		ex_l[i]=sqrt(pow(TG1->GetErrorXlow(i),2)+pow(TG2->GetErrorXlow(i),2))/2;
		ex_h[i]=sqrt(pow(TG1->GetErrorXhigh(i),2)+pow(TG2->GetErrorXhigh(i),2))/2;
	}

	if (syst==1){
		auto TG1_syst = (TGraphAsymmErrors *) tl1->At(1);
  		auto TG2_syst = (TGraphAsymmErrors *) tl2->At(1);
  		double x_syst[_nBins];
		double y_syst[_nBins];
		double ey_syst[_nBins];
		double ex_l_syst[_nBins];
		double ex_h_syst[_nBins];
	
		for (int i=0;i<_nBins;i++){
			y_syst[i]=TG1_syst->GetY()[i]/TG2_syst->GetY()[i];
			x_syst[i]=(TG1_syst->GetX()[i]+TG2_syst->GetX()[i])/2;
			ey_syst[i]=abs(y_syst[i])*sqrt(pow(TG1_syst->GetErrorY(i)/TG1_syst->GetY()[i],2)+pow(TG2_syst->GetErrorY(i)/TG2_syst->GetY()[i],2));
			ex_l_syst[i]=sqrt(pow(TG1_syst->GetErrorXlow(i),2)+pow(TG2_syst->GetErrorXlow(i),2))/2;
			ex_h_syst[i]=sqrt(pow(TG1_syst->GetErrorXhigh(i),2)+pow(TG2_syst->GetErrorXhigh(i),2))/2;
	}		
		TGraphAsymmErrors *g_ratio_syst= new TGraphAsymmErrors (_nBins,x_syst,y_syst,ex_l_syst,ex_h_syst,ey_syst,ey_syst);
		g_ratio_syst->SetLineColor(2);
		m_ratio->Add(g_ratio_syst);
		leg_ratio->AddEntry(g_ratio_syst,"Systematic Uncertainty", "e");
	}

	
	TGraphAsymmErrors *g_ratio= new TGraphAsymmErrors (_nBins,x,y,ex_l,ex_h,ey,ey);
	g_ratio->SetMarkerColor(1);
	g_ratio->SetMarkerStyle(21);
	
	m_ratio->Add(g_ratio);
	m_ratio->GetYaxis()->SetTitle("B_{s}/B^{+} yield ratio");
	
	double yield_max=0;
	for(int i = 0; i < _nBins; i++){
		if(y[i] > yield_max){
		yield_max = y[i];
		}
	}
	m_ratio->GetYaxis()->SetRangeUser(0, yield_max * 1.5);

	if(varExp == "By"){
		 m_ratio->GetXaxis()->SetTitle("Rapidity (y)");
		 //m_ratio->GetYaxis()->SetTitle("B_{s}/B^{+} yield ratio");
		 m_ratio->GetXaxis()->SetLimits(-2.4 ,2.4);
	 }
	 if(varExp == "Bpt"){
		 m_ratio->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 //m_ratio->GetYaxis()->SetTitle("dN_{S}/dp_{T}");
		//if (tree == "ntKp"){ m_ratio->GetXaxis()->SetLimits(0 ,80); }
		//if (tree == "ntphi"){ m_ratio->GetXaxis()->SetLimits(0 ,60); }
	 }
	 if(varExp == "nMult"){
		 m_ratio->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 //m_ratio->GetYaxis()->SetTitle("dN_{S}/dMult");
		 m_ratio->GetXaxis()->SetLimits(0, 110);
	 }

	m_ratio->Draw("ap");
	leg_ratio->AddEntry(g_ratio, "Statistical Uncertainty", "e");
	leg_ratio->SetBorderSize(0);
	leg_ratio->SetFillStyle(0);
	leg_ratio->SetTextSize(0);
	leg_ratio->Draw();
	c->SaveAs(Form("./results/%s_ratioplot.png",varExp.Data())); 
	return;
}