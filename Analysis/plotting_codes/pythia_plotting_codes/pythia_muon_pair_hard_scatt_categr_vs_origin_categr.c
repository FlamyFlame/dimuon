#include <TROOT.h>
// #include <TChain.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include "vector"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

#include<algorithm>
// #include "string"
// #include <vector>
// #include "time.h"
// #include "struct_hist.h"


const int nSigns = 2;
const int nflavors = 3;
const int nDphi = 2;

std::string flavors[nflavors] = {"bb","bc","cc"};
std::string py_path = "/usatlas/u/yuhanguo/usatlasdata/pythia/";
// std::string fname = "muon_pairs_pythia.root";
std::string fname = "muon_pairs_pythia_0429.root";
std::vector<std::vector<std::vector<std::string>>> hist_names = {};

std::string labels[nflavors][nSigns][nDphi] = {
  {{"both from b, same sign, #Delta #phi < #pi/2","both from b, same sign, #Delta #phi #geq #pi/2"}, {"both from b, opposite sign, #Delta #phi < #pi/2","both from b, opposite sign, #Delta #phi #geq #pi/2"}}, 
  {{"one from b one from c, same sign, #Delta #phi < #pi/2","one from b one from c, same sign, #Delta #phi #geq #pi/2"}, {"one from b one from c, opposite sign, #Delta #phi < #pi/2","one from b one from c, opposite sign, #Delta #phi #geq #pi/2"}}, 
  {{"both from c, same sign, #Delta #phi < #pi/2","both from c, same sign, #Delta #phi #geq #pi/2"}, {"both from c, opposite sign, #Delta #phi < #pi/2","both from c, opposite sign, #Delta #phi #geq #pi/2"}}
};

void hist_helper(TH1* h, char * title){
  h->SetStats(0);
  h->SetTitle(title);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(38);
  h->GetYaxis()->SetLabelOffset(0.01);
  // h->GetYaxis()->SetTitleFont(43);
  // h->GetYaxis()->SetTitleSize(32);
  // h->GetYaxis()->SetTitleOffset(2);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(38);
  h->GetXaxis()->SetLabelOffset(0.01);
  // h->GetXaxis()->SetTitleFont(43);  
  // h->GetXaxis()->SetTitleSize(32);
  // h->GetXaxis()->SetTitleOffset(1);
}

void pythia_muon_pair_hard_scatt_categr_vs_origin_categr(){

  TH2D* h;
  TFile* f = TFile::Open((py_path + fname).c_str());

  hist_names.clear();
  hist_names.push_back({{"h_both_from_b_hard_scatt_categr_vs_origin_categr_sign1_near","h_both_from_b_hard_scatt_categr_vs_origin_categr_sign1_away"},{"h_both_from_b_hard_scatt_categr_vs_origin_categr_sign2_near","h_both_from_b_hard_scatt_categr_vs_origin_categr_sign2_away"}});
  hist_names.push_back({{"h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr_sign1_near","h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr_sign1_away"},{"h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr_sign2_near","h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr_sign2_away"}});
  hist_names.push_back({{"h_both_from_c_hard_scatt_categr_vs_origin_categr_sign1_near","h_both_from_c_hard_scatt_categr_vs_origin_categr_sign1_away"},{"h_both_from_c_hard_scatt_categr_vs_origin_categr_sign2_near","h_both_from_c_hard_scatt_categr_vs_origin_categr_sign2_away"}});

  for (int imc = 0; imc < nflavors; imc++){

    TCanvas* c = new TCanvas(Form("c-%s",flavors[imc].c_str()),Form("c-%s",flavors[imc].c_str()),2900,2000);
    c->Divide(2,2);

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){
      for (int lphi = 0; lphi < nDphi; lphi++){

        c->cd(ksign * nSigns + lphi + 1);
        h = (TH2D*) f->Get(hist_names[imc][ksign][lphi].c_str());
        h->Scale(1000);
        hist_helper(h, Form("%s (%.2f nb)",(labels[imc][ksign][lphi]).c_str(), h->Integral()));
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.135);

        // TLegend* l = new TLegend(0.65,0.74,0.87,0.87);
        // l->SetBorderSize(0);
        // l->SetFillStyle(0);
        // l->SetTextFont(38);
        // l->SetTextSize(l->GetTextSize()*3);
        // l->SetMargin(0.02);
        // l->SetTextColor(1);
        // l->AddEntry(h, labels[imc][ksign][lphi].c_str(),"lp");
        h->Draw("colz");
        // l->Draw(); 
      }
    }

    // if (mode == 1) c->SaveAs(Form("plots/mc_truth/prt_grouping/pythia_muon_pair_dp_ancestor_groups_%s_unweighted.png",flavors[imc].c_str()));
    // else           c->SaveAs(Form("plots/mc_truth/prt_grouping/pythia_muon_pair_dp_ancestor_groups_%s.png",flavors[imc].c_str()));
    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/pythia/plots/ancestor_tracing/pythia_muon_pair_hard_scatt_categr_vs_origin_categr_%s.png",flavors[imc].c_str()));
    c->Close();
    delete c;
  }
}





