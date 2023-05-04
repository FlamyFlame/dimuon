void hist_helper(TH1* h, std::string xtitle = ""){
	h->SetStats(0);
	h->Scale(1.,"width");
	if (xtitle != "") h->GetXaxis()->SetTitle(xtitle.c_str());
	h->GetYaxis()->SetLabelFont(43);
	h->GetYaxis()->SetLabelSize(32);
	h->GetYaxis()->SetLabelOffset(0.01);
	h->GetYaxis()->SetTitleFont(43);
	h->GetYaxis()->SetTitleSize(32);
	h->GetYaxis()->SetTitleOffset(1.8);
	h->GetXaxis()->SetLabelFont(43);
	h->GetXaxis()->SetLabelSize(32);
	h->GetXaxis()->SetLabelOffset(0.02);
	h->GetXaxis()->SetTitleFont(43);  
	h->GetXaxis()->SetTitleSize(32);
	h->GetXaxis()->SetTitleOffset(1);
	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.7);
}


void plot_crossx_bb_cc(){
  	TCanvas* c = new TCanvas("c","c",1200,800);
  	gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.135);

    TLegend* l = new TLegend(0.2,0.7,0.5,0.85);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    l->SetTextSize(l->GetTextSize()*2.4);
    // l->SetMarkerSize(l->GetMarkerSize()*3);
    // l->SetMarkerStyle(20);
    l->SetMargin(0.2);
    l->SetTextColor(1);

  	TFile* f[2];
  	f[0] = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_bb.root");
  	f[1] = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_cc.root");
  	
  	TH1D* h[2];
    h[0] = (TH1D*) f[0]->Get("h_crossx");
    h[1] = (TH1D*) f[1]->Get("h_crossx");

  	hist_helper(h[0],"#sigma (pb^{-1})");
  	hist_helper(h[1],"#sigma (pb^{-1})");

  	gPad->SetLogx();
  	gPad->SetLogy();

    // l->AddEntry(h[0],"bb","lp");
    // l->AddEntry(h[1],"cc","lp");
    l->AddEntry(h[0],"bb","lp");
    l->AddEntry(h[1],"cc","lp");

    h[0]->SetMarkerColor(kBlack);
    h[0]->SetLineColor(kBlack);
    h[1]->SetMarkerColor(kRed);
    h[1]->SetLineColor(kRed);

    // h[0]->GetYaxis()->SetRangeUser(0,1220);

  	h[1]->Draw("E");
  	h[0]->Draw("E,same");
    l->Draw();

  	c->SaveAs("plots/mc_truth/crossx_bb_cc_compr.png");
  	c->Close();
  	delete c;


}