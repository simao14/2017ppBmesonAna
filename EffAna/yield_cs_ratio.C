#include <cmath> 

void yield_cs_ratio(TString varExp){

	gSystem->mkdir("ratios",true);
	
	TString* var;
	if(varExp == "By"){var = new TString("Y");}
	if(varExp == "nMult"){var = new TString("Mult");}
	if(varExp == "Bpt"){var = new TString("PT");}

	TFile *cross_bs = new TFile(Form("Bs/FinalFiles/BsPPCorrYield%s.root",var->Data()),"read");
	TFile *cross_bp = new TFile(Form("BP/FinalFiles/BPPPCorrYield%s.root",var->Data()),"read");

	TH1D *TH_cross_bs = (TH1D*) cross_bs->Get("hPtSigma");

	TH1D *TH_cross_bp = (TH1D*) cross_bp->Get("hPtSigma");

	TCanvas* c=new TCanvas();
	c->cd();
	TLegend* leg_ratio=new TLegend(0.7,0.7,0.9,0.9);
	TH1D * cs_ratio = (TH1D*) TH_cross_bs->Clone("cs_ratio");
	cs_ratio->Divide(TH_cross_bp);

	cs_ratio->SetMarkerStyle(20);
	cs_ratio->SetMarkerSize(1);
	cs_ratio->SetMarkerColor(kBlack);
	cs_ratio->SetLineColor(kBlack);

	if(varExp == "By"){
		 cs_ratio->GetXaxis()->SetTitle("Rapidity (y)");
		 cs_ratio->GetYaxis()->SetTitle("B_{s}/B^{+} d #sigma/d y ratio");
		 cs_ratio->GetXaxis()->SetLimits(-2.4 ,2.4);
	 }
	 if(varExp == "Bpt"){
		 cs_ratio->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 cs_ratio->GetYaxis()->SetTitle("B_{s}/B^{+} d #sigma/d p_{T} ratio");
	 }
	 if(varExp == "nMult"){
		 cs_ratio->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 cs_ratio->GetYaxis()->SetTitle("B_{s}/B^{+} d #sigma/d Mult ratio");
		 cs_ratio->GetXaxis()->SetLimits(0, 110);
	 }

	cs_ratio->Draw("ep");
	//leg_ratio->AddEntry(g_ratio, "Statistical Uncertainty", "e");
	leg_ratio->SetBorderSize(0);
	leg_ratio->SetFillStyle(0);
	leg_ratio->SetTextSize(0);
	//leg_ratio->Draw();
	c->SaveAs(Form("ratios/%s_cs_ratioplot.pdf",varExp.Data())); 
	return;
}