void hist_helper(TH1* h, float norm = 1){
	h->SetStats(0);
	h->Scale(norm,"width");
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


void dphi_same_b_all_op_compr(){
  	TCanvas* c = new TCanvas("c","c",1200,800);
  	gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.135);

    TLegend* l = new TLegend(0.2,0.7,0.5,0.87);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    l->SetTextSize(l->GetTextSize()*3);
    l->SetMargin(0.02);
    l->SetTextColor(1);

  	TFile* f = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_bb.root");
  	TH1D* h[2];

    h[0] = (TH1D*) f->Get("h_dphi_bb_op_near");
    h[1] = (TH1D*) f->Get("h_dphi_bb_op_near_from_same_b");

	h[0]->Rebin(2);
	h[1]->Rebin(2);

  	hist_helper(h[0]);
  	hist_helper(h[1]);

  	h[0]->SetTitle("bb sample, near side");

    l->AddEntry(h[0],"opposite sign, all","lp");
    l->AddEntry(h[1],"opposite sign, muons from same b only","lp");

    h[0]->SetMarkerColor(kBlack);
    h[0]->SetLineColor(kBlack);
    h[1]->SetMarkerColor(kRed);
    h[1]->SetLineColor(kRed);

    h[0]->GetYaxis()->SetRangeUser(0,1220);

  	h[0]->Draw("E");
  	h[1]->Draw("E,same");
    l->Draw();

  	c->SaveAs("plots/mc_truth/dphi_same_b_all_op_compr.png");
  	c->Close();
  	delete c;


}