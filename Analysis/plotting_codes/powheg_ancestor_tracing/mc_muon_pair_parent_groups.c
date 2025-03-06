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
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/plotting_codes/powheg_ancestor_tracing/powheg_settings.h"
// #include "string"
// #include <vector>
// #include "time.h"
// #include "struct_hist.h"

std::vector<std::vector<std::string>> hist_names = {{"h_parent_groups_sign1_near","h_parent_groups_sign1_away"},{"h_parent_groups_sign2_near","h_parent_groups_sign2_away"}};

std::string labels[nMCmodes][nSigns][nDphi] = {
  {{"bb, same sign, #Delta #phi < #pi/2","bb, same sign, #Delta #phi #geq #pi/2"}, {"bb, opposite sign, #Delta #phi < #pi/2","bb, opposite sign, #Delta #phi #geq #pi/2"}}, 
  {{"cc, same sign, #Delta #phi < #pi/2","cc, same sign, #Delta #phi #geq #pi/2"}, {"cc, opposite sign, #Delta #phi < #pi/2","cc, opposite sign, #Delta #phi #geq #pi/2"}}
};

// void hist_helper(TH1* h, std::string title){
//   h->SetStats(0);
//   h->SetTitle(title.c_str());
//   h->GetYaxis()->SetLabelFont(43);
//   h->GetYaxis()->SetLabelSize(38);
//   h->GetYaxis()->SetLabelOffset(0.01);
//   // h->GetYaxis()->SetTitleFont(43);
//   // h->GetYaxis()->SetTitleSize(32);
//   // h->GetYaxis()->SetTitleOffset(2);
//   h->GetXaxis()->SetLabelFont(43);
//   h->GetXaxis()->SetLabelSize(38);
//   h->GetXaxis()->SetLabelOffset(0.01);
//   // h->GetXaxis()->SetTitleFont(43);  
//   // h->GetXaxis()->SetTitleSize(32);
//   // h->GetXaxis()->SetTitleOffset(1);
// }

void mc_muon_pair_parent_groups(){

  TH2D* h;
  TFile* f[nMCmodes][nBatches];

  for (int imc = 0; imc < nMCmodes; imc++){

    TCanvas* c = new TCanvas(Form("c-%s",mcmodes[imc].c_str()),Form("c-%s",mcmodes[imc].c_str()),2900,2000);
    c->Divide(2,2);

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){
      for (int lphi = 0; lphi < nDphi; lphi++){

        c->cd(ksign * nSigns + lphi + 1);

        f[imc][0] = TFile::Open((mcdir + sub_dirs[imc] + fnames[imc][0]).c_str());
        h = (TH2D*) f[imc][0]->Get(hist_names[ksign][lphi].c_str());

        for (int jbatch = 1; jbatch < nBatches; jbatch++){
          f[imc][jbatch] = TFile::Open((mcdir + sub_dirs[imc] + fnames[imc][jbatch]).c_str());
          TH2D* h_cur_batch = (TH2D*) f[imc][jbatch]->Get(hist_names[ksign][lphi].c_str());
          h->Add(h_cur_batch);
        }

        h->Scale(1./4999000);
        hist_helper(h,labels[imc][ksign][lphi]);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.135);

        // TLegend* l = new TLegend(0.65,0.74,0.87,0.87);
        // l->SetBorderSize(0);
        // l->SetFillStyle(0);
        // l->SetTextFont(38);
        // l->SetTextSize(l->GetTextSize()*3);
        // l->SetMargin(0.2);
        // l->SetTextColor(1);
        // l->AddEntry(h, labels[imc][ksign][lphi].c_str(),"lp");
        h->Draw("colz");
        // l->Draw(); 
      }
    }

    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/powheg/prt_grouping/muon_pair_parent_groups%s.png",mcmodes[imc].c_str()));
    c->Close();
    delete c;
  }
}




