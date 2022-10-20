#include <cmath> 

void yield_cs(int syst,TString varExp,TString tree){

	
	double brafrac_ntKp = 0.00102;
	double brafrac_ntKp_err = 0.000019;

	double brafrac_ntphi = 0.00079;
	double brafrac_ntphi_err = 0.00007;
	
	double B;
	double B_err;
	if(tree == "ntphi"){
		B = brafrac_ntphi;
		B_err = brafrac_ntphi_err;
	}
	if(tree == "ntKp"){
		B = brafrac_ntKp;
		B_err = brafrac_ntKp_err;
	}



	double lumi = 302.3;
	
	TString* var;
	TString* par;
	if(varExp == "By"){var = new TString("Y");}
	if(varExp == "nMult"){var = new TString("Mult");}
	if(varExp == "Bpt"){var = new TString("PT");}

	if(tree == "ntphi"){par = new TString("Bs");}
	if(tree == "ntKp"){par = new TString("BP");}

	TFile *diff_f = new TFile(Form("./results/%s_%s_ratio.root",tree.Data(),varExp.Data()),"read");
	TFile *eff_f = new TFile(Form("./FinalFiles/%sPPCorrYield%s.root",par->Data(),var->Data()),"read");

	TMultiGraph *TMG_diff = (TMultiGraph*) diff_f->Get("TG");
	TH1D *TH_eff = (TH1D*) eff_f->Get("hInvEff");

	auto tl_diff = TMG_diff->GetListOfGraphs();


	auto TG_diff = (TGraphAsymmErrors *) tl_diff->At(0);
	
	const int _nBins = (TG_diff->GetN());
	TCanvas* c=new TCanvas();
	TLegend* leg_diff=new TLegend(0.7,0.7,0.9,0.9);
	TMultiGraph* m_diff=new TMultiGraph();
	double x[_nBins];
	double y[_nBins];
	double ey[_nBins];
	double ex_l[_nBins];
	double ex_h[_nBins];
	
	for (int i=0;i<_nBins;i++){
		y[i]=TG_diff->GetY()[i] * (TH_eff->GetBinContent(i+1) / (B*lumi));
		x[i]=TG_diff->GetX()[i];
		ey[i]=abs(y[i]) * (TG_diff->GetErrorY(i)/TG_diff->GetY()[i]);
		ex_l[i]=TG_diff->GetErrorXlow(i);
		ex_h[i]=TG_diff->GetErrorXhigh(i);
	}

	if (syst==1){
		auto TG_syst = (TGraphAsymmErrors *) tl_diff->At(1);
  		double x_syst[_nBins];
		double y_syst[_nBins];
		double ex_syst_l[_nBins];
		double ex_syst_h[_nBins];
		double ey_syst[_nBins];
	
		for (int i=0;i<_nBins;i++){
			y_syst[i]=TG_syst->GetY()[i] * (TH_eff->GetBinContent(i+1) / (B*lumi));
			x_syst[i]=TG_syst->GetX()[i];
			ey_syst[i]=abs(y_syst[i])*sqrt(pow(TG_syst->GetErrorY(i)/TG_syst->GetY()[i],2)+pow(TH_eff->GetBinError(i+1)/TH_eff->GetBinContent(i+1),2) + pow(B_err/B,2) + pow(0.019,2));
		ex_syst_l[i]=TG_syst->GetErrorXlow(i);
		ex_syst_h[i]=TG_syst->GetErrorXhigh(i);

		}
		
		TGraphAsymmErrors *g_diff_syst= new TGraphAsymmErrors (_nBins,x_syst,y_syst,ex_syst_l,ex_syst_h,ey_syst,ey_syst);
		g_diff_syst->SetLineColor(2);
		m_diff->Add(g_diff_syst);
		leg_diff->AddEntry(g_diff_syst,"Systematic Uncertainty", "e");
	}

	
	TGraphAsymmErrors *g_diff= new TGraphAsymmErrors (_nBins,x,y,ex_l,ex_h,ey,ey);
	g_diff->SetMarkerColor(1);
	//g_diff->SetMarkerStyle(21);
	
	m_diff->Add(g_diff);
	m_diff->GetYaxis()->SetTitle("Differential Cross Section");
	
	double yield_max=0;
	for(int i = 0; i < _nBins-1; i++){
		if(y[i] > yield_max){
		yield_max = y[i];
		}
	}
	m_diff->GetYaxis()->SetRangeUser(0.01, yield_max * 1.5);
	m_diff->GetYaxis()->SetTitle(Form("d#sigma/d%s", varExp.Data()));
	if(varExp == "By"){
		 m_diff->GetXaxis()->SetTitle("Rapidity (y)");
		 m_diff->GetYaxis()->SetTitle("d#sigma/dy (pb)");
		 m_diff->GetXaxis()->SetLimits(-2.4 ,2.4);
	 }
	 if(varExp == "Bpt"){
		 m_diff->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 c->SetLogy();
		 //m_diff->GetYaxis()->SetTitle("dN_{S}/dp_{T}");
		//if (tree == "ntKp"){ m_diff->GetXaxis()->SetLimits(0 ,80); }
		//if (tree == "ntphi"){ m_diff->GetXaxis()->SetLimits(0 ,60); }
	 }
	 if(varExp == "nMult"){
		 m_diff->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 m_diff->GetYaxis()->SetTitle("d#sigma/dMult (pb)");
		 m_diff->GetXaxis()->SetLimits(0, 110);
	 }

	m_diff->Draw("ap");
	leg_diff->AddEntry(g_diff, "Statistical Uncertainty", "e");
	leg_diff->SetBorderSize(0);
	leg_diff->SetFillStyle(0);
	leg_diff->SetTextSize(0);
	leg_diff->Draw();
	c->SaveAs(Form("./results/%s_%s_diffcrosssectionplot.png",tree.Data(),varExp.Data())); 
	return;
}
