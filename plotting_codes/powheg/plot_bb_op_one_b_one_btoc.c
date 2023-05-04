void hist_helper(TH1* h, bool scale_to_unity, float norm = 1.){
	h->SetStats(0);
	// if (scale_to_unity) h->Scale(1./h->Integral());
	// else h->Scale(norm,"width");
	h->GetYaxis()->SetLabelFont(43);
	h->GetYaxis()->SetLabelSize(32);
	h->GetYaxis()->SetLabelOffset(0.01);
	h->GetYaxis()->SetTitleFont(43);
	h->GetYaxis()->SetTitleSize(32);
	h->GetYaxis()->SetTitleOffset(1.5);
	h->GetXaxis()->SetLabelFont(43);
	h->GetXaxis()->SetLabelSize(32);
	h->GetXaxis()->SetLabelOffset(0.02);
	h->GetXaxis()->SetTitleFont(43);  
	h->GetXaxis()->SetTitleSize(32);
	h->GetXaxis()->SetTitleOffset(1);
	h->SetMarkerStyle(20);
	h->SetMarkerSize(1.2);
}


void plot_bb_op_one_b_one_btoc(){
  	TCanvas* c = new TCanvas("c","c",2200,1000);
  	c->Divide(2,1);

    // TLegend* l = new TLegend(0.2,0.7,0.5,0.87);
    // l->SetBorderSize(0);
    // l->SetFillStyle(0);
    // l->SetTextFont(42);
    // l->SetTextSize(l->GetTextSize()*3);
    // l->SetMargin(0.2);
    // l->SetTextColor(1);

  	TFile* f = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_bb.root");
  	TH1D* h[2];

    h[0] = (TH1D*) f->Get("h_bb_op_one_b_one_btoc_near");
    h[1] = (TH1D*) f->Get("h_bb_op_one_b_one_btoc_away");

  	// hist_helper(h[0],true);
  	// hist_helper(h[1],true);
  	hist_helper(h[0],false);
  	hist_helper(h[1],false);

  	// h[0]->SetTitle("bb sample, same sign, both from b, from same ancestors, near side");
  	h[0]->SetTitle("near side");
  	h[1]->SetTitle("away side");

    // l->AddEntry(h[0],"near side","lp");
    // l->AddEntry(h[1],"away side","lp");

    // h[0]->SetMarkerColor(kRed);
    // h[0]->SetLineColor(kRed);
    // h[1]->SetMarkerColor(kBlue);
    // h[1]->SetLineColor(kBlue);

    c->cd(1);
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.16);
  	h[0]->Draw("E");

  	c->cd(2);
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.16);
  	h[1]->Draw("E");

  	c->SaveAs("plots/mc_truth/bb_op_one_b_one_btoc.png");
  	c->Close();
  	delete c;

}