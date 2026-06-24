// Clean D_OS vs D_SS low-mass overlay for the interview talk.
// Source: pp24 template-fit histos (no_res_cut), cross-section weighted.
#include "AtlasStyle.C"
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

void make_lowmass_overlay(){
  SetAtlasStyle();

  const char* fpath = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                      "histograms_real_pairs_pp_2024_2mu4_nominal_template_fit.root";
  TFile* f = TFile::Open(fpath);
  TH1D* hOS = (TH1D*) f->Get("h1d_crossx_minv_0_4_op_dsigma");
  TH1D* hSS = (TH1D*) f->Get("h1d_crossx_minv_0_4_ss_dsigma");

  hOS->SetLineColor(kBlack);   hOS->SetMarkerColor(kBlack);   hOS->SetMarkerStyle(20); hOS->SetMarkerSize(0.7);
  hSS->SetLineColor(kRed+1);   hSS->SetMarkerColor(kRed+1);   hSS->SetMarkerStyle(24); hSS->SetMarkerSize(0.7);
  hOS->SetLineWidth(2);        hSS->SetLineWidth(2);
  hOS->SetStats(0);            hSS->SetStats(0);
  hOS->SetTitle("");

  hOS->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  hOS->GetYaxis()->SetTitle("d#sigma/dm_{#mu#mu}  [arb.]");

  // y range for log scale: floor above zero, headroom for the J/psi spike
  double ymax = hOS->GetMaximum() * 3.0;
  hOS->GetYaxis()->SetRangeUser(0.5, ymax);

  TCanvas* c = new TCanvas("c_lowmass","c_lowmass",700,560);
  c->SetLogy();
  hOS->Draw("E1");
  hSS->Draw("E1 same");

  // signal window markers [1.08, 2.9]
  c->Update();
  TLine* l1 = new TLine(1.08, 0.5, 1.08, ymax); l1->SetLineStyle(2); l1->SetLineColor(kGray+2); l1->Draw();
  TLine* l2 = new TLine(2.90, 0.5, 2.90, ymax); l2->SetLineStyle(2); l2->SetLineColor(kGray+2); l2->Draw();

  TLegend* leg = new TLegend(0.55,0.74,0.92,0.90);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.038);
  leg->AddEntry(hOS, "Opposite-sign (data)", "lp");
  leg->AddEntry(hSS, "Same-sign (data)", "lp");
  leg->Draw();

  TLatex* atlas = new TLatex(); atlas->SetNDC(); atlas->SetTextFont(72); atlas->SetTextSize(0.050);
  atlas->DrawLatex(0.20, 0.88, "ATLAS");
  TLatex* lab = new TLatex(); lab->SetNDC(); lab->SetTextFont(42); lab->SetTextSize(0.042);
  lab->DrawLatex(0.345, 0.88, "Internal");
  TLatex* sub = new TLatex(); sub->SetNDC(); sub->SetTextFont(42); sub->SetTextSize(0.034);
  sub->DrawLatex(0.20, 0.82, "pp 2024, single-b selection (no mass cut)");
  TLatex* win = new TLatex(); win->SetNDC(); win->SetTextFont(42); win->SetTextSize(0.030); win->SetTextColor(kGray+3);
  win->DrawLatex(0.43, 0.20, "signal window");

  const char* out = "/usatlas/u/yuhanguo/usatlasdata/ew_onia_slides/talk_figs/FIG-LOWMASS-OVERLAY.png";
  c->SaveAs(out);
  printf("wrote %s\n", out);
}
