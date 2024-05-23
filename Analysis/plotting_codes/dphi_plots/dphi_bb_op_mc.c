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


void dphi_bb_op_mc(){
  	TCanvas* c = new TCanvas("c","c",1200,800);
  	gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.135);

    TLegend* l = new TLegend(0.2,0.7,0.5,0.87);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    l->SetTextSize(l->GetTextSize()*3);
    l->SetMargin(0.2);
    l->SetTextColor(1);

  	TFile* f[3];
  	TH2D* h2d;
  	TH1D* h[3];

  	f[0] = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_bb.root");
  	f[1] = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_bb.root");
  	f[2] = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/histograms_mc_truth_bb.root");

    h[0] = (TH1D*) f[0]->Get("h_dphi_bb_op_near_both_from_b");
    h[1] = (TH1D*) f[1]->Get("h_dphi_bb_op_near_one_b_one_btoc");
  	h2d = (TH2D*) f[2]->Get("h_Deta_Dphi_dr3_sign2_gapcut1");
    h[2] = (TH1D*) h2d->ProjectionX("h_dphi");

	h[0]->Rebin(2);
	h[1]->Rebin(2);
	h[2]->Rebin(2);

  	hist_helper(h[0]);
  	hist_helper(h[1]);
  	hist_helper(h[2]);


    l->AddEntry(h[0],"both direct b","lp");
    l->AddEntry(h[1],"both direct b + one b, one b-to-c","lp");
    l->AddEntry(h[2],"all pairs","lp");

    h[0]->SetMarkerColor(kRed);
    h[0]->SetLineColor(kRed);
    h[1]->SetMarkerColor(kBlue);
    h[1]->SetLineColor(kBlue);
    h[2]->SetMarkerColor(kBlack);
    h[2]->SetLineColor(kBlack);

  	h[0]->SetTitle("bb sample, opposite sign, near side");
    h[0]->GetYaxis()->SetRangeUser(0,1220);

  	h[0]->Draw();
  	h[1]->Add(h[0]);
  	h[1]->Draw("same");
  	h[2]->Draw("same");
    l->Draw();

  	c->SaveAs("plots/mc_truth/dphi_bb_op_mc.png");
  	c->Close();
  	delete c;


}