//  ──────────────────────────────────────────────────────────────────────
//  TrigEffPlotter.C          (compile once with .L TrigEffPlotter.C+)
//  ──────────────────────────────────────────────────────────────────────
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <iostream>
#include <vector>
#include <map>

class TrigEffPlotter {
public:
   /// ctor:  file name  |  1-D vars               |  optional log settings
   ///                   |  2-D vars (optional)
   TrigEffPlotter(const std::string&  fileName,
                  const std::vector<std::string>& vars1D,
                  const std::map<std::string,std::pair<bool,bool>>& logopt_ = {},
                  const std::vector<std::string>& vars2D_ = {})
      : fFile(nullptr), variables(vars1D), var2Ds(vars2D_), logopt(logopt_)
   { initialize(fileName); }

   /// loop over all variables (1-D and 2-D)
   void Run() {
      if(!fFile || fFile->IsZombie()) return;
      gStyle->SetOptStat(0);
      for(const auto& v : variables) plotOne(v);
      for(const auto& v : var2Ds )  plot2DOne(v);
   }

private:
   std::string data_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_2024/";
   TFile*                          fFile;
   std::vector<std::string>        variables;   // 1-D
   std::vector<std::string>        var2Ds;      // 2-D  (NEW)
   std::map<std::string,std::pair<bool,bool>> logopt;

   void initialize(const std::string& fn) {
      fFile = TFile::Open(fn.c_str(),"READ");
      if(!fFile || fFile->IsZombie())
         std::cerr<<"[TrigEffPlotter] ERROR: cannot open "<<fn<<std::endl;
   }

   // ---------- helpers ----------------------------------------------------
   template<typename T>
   T* getObj(const std::string& n) {
      T* h = dynamic_cast<T*>(fFile->Get(n.c_str()));
      if(!h) std::cerr<<"[TrigEffPlotter] WARNING: "<<n<<" not found\n";
      return h;
   }
   void setStyle(TH1* h,int col,int m) {
      if(!h) return;
      h->SetMarkerStyle(m); h->SetMarkerColor(col);
      h->SetLineColor(col); h->SetLineWidth(2);
   }
   void scaleErr(TH1* h, double frac = 0.1) {
      if(!h) return;
      for(int b = 1; b <= h->GetNbinsX(); ++b)
         h->SetBinError(b, frac * h->GetBinContent(b));   // 10 % error bars
   }
   void applyLog(const std::pair<bool,bool>& xy){
      gPad->SetLogx(xy.first); gPad->SetLogy(xy.second);
   }
   std::pair<bool,bool> xyFor(const std::string& v) const {
      auto it = logopt.find(v);
      return (it==logopt.end()) ? std::make_pair(false,false) : it->second;
   }

   // ---------- 1-D plotting (unchanged from previous reply) ---------------
   void plotOne(const std::string& var){
      const auto xy = xyFor(var);
      // names --------------------------------------------------------------
      std::string h_ss_1 = Form("h_%s_mu4_mu4noL1_given_mu4_sign1",var.c_str());
      std::string h_ss_2 = Form("h_%s_2mu4_given_mu4_sign1",      var.c_str());
      std::string h_os_1 = Form("h_%s_mu4_mu4noL1_given_mu4_sign2",var.c_str());
      std::string h_os_2 = Form("h_%s_2mu4_given_mu4_sign2",      var.c_str());
      TH1* h1_ss = getObj<TH1>(h_ss_1);
      TH1* h2_ss = getObj<TH1>(h_ss_2);
      TH1* h1_os = getObj<TH1>(h_os_1);
      TH1* h2_os = getObj<TH1>(h_os_2);
      setStyle(h1_ss,kRed+1,20);  setStyle(h2_ss,kAzure+2,24);
      setStyle(h1_os,kRed+1,20);  setStyle(h2_os,kAzure+2,24);

      scaleErr(h1_ss);  scaleErr(h2_ss); // scale error down by 1/10
      scaleErr(h1_os);  scaleErr(h2_os); // scale error down by 1/10

      // ── 1st canvas: 2 pads ---------------------------------------------
      TCanvas* c1 = new TCanvas(Form("c1_%s",var.c_str()),var.c_str(),1200,500);
      c1->Divide(2,1);

      c1->cd(1); applyLog(xy);
      if(h1_ss){  h1_ss->SetTitle("same sign"); 
                  h1_ss->GetYaxis()->SetRangeUser(0.,1.);
                  h1_ss->Draw("E1"); }
      if(h2_ss)  h2_ss->Draw("E1 SAME");
      { TLegend* L=new TLegend(0.52,0.17,0.84,0.31 );
        L->AddEntry(h1_ss,"mu4_mu4noL1","lep");
        L->AddEntry(h2_ss,"2mu4","lep"); 
        L->AddEntry("","errorbar/10",""); 
        L->Draw(); }

      c1->cd(2); applyLog(xy);
      if(h1_os){  h1_os->SetTitle("opposite sign"); 
                  h1_os->GetYaxis()->SetRangeUser(0.,1.);
                  h1_os->Draw("E1"); }
      if(h2_os)  h2_os->Draw("E1 SAME");
      { TLegend* L=new TLegend(0.52,0.17,0.84,0.31 );
        L->AddEntry(h1_os,"mu4_mu4noL1","lep");
        L->AddEntry(h2_os,"2mu4","lep"); 
        L->AddEntry("","errorbar/10",""); 
        L->Draw(); }

      c1->SaveAs(Form("%strig_effcy_plots/%s_trig_effcy.png", data_dir.c_str(), var.c_str()));

      // ── 2nd canvas: w_sig_sel ------------------------------------------
      std::string h_sel_1 = Form("h_%s_mu4_mu4noL1_given_mu4_w_sig_sel_sign2",var.c_str());
      std::string h_sel_2 = Form("h_%s_2mu4_given_mu4_w_sig_sel_sign2",      var.c_str());
      TH1* h1_sel = getObj<TH1>(h_sel_1);
      TH1* h2_sel = getObj<TH1>(h_sel_2);
      setStyle(h1_sel,kRed+1,20); setStyle(h2_sel,kAzure+2,24);
      
      scaleErr(h1_sel); scaleErr(h2_sel);  // scale error down by 1/10

      TCanvas* c2 = new TCanvas(Form("c2_%s",var.c_str()),var.c_str(),700,600);
      applyLog(xy);
      if(h1_sel){
         h1_sel->GetYaxis()->SetRangeUser(0.,1.);
         h1_sel->Draw("E1");
      }
      if(h2_sel) h2_sel->Draw("E1 SAME");
      { TLegend* L=new TLegend(0.52,0.17,0.84,0.31 );
        L->AddEntry(h1_sel,"mu4_mu4noL1","lep");
        L->AddEntry(h2_sel,"2mu4","lep"); 
        L->AddEntry("","errorbar/10",""); 
        L->Draw(); }

      c2->SaveAs(Form("%strig_effcy_plots/%s_trig_effcy_w_sig_sel.png", data_dir.c_str(), var.c_str()));
   }

   // ---------- NEW: 2-D plotting -----------------------------------------
   void plot2DOne(const std::string& var){
      const auto xy = xyFor(var);
      // first 2×2 canvas ---------------------------------------------------
      struct HN{ std::string tag; std::string title; };
      const HN names[4] = {
         {Form("h_%s_mu4_mu4noL1_given_mu4_sign1",var.c_str()), "mu4_mu4noL1, same sign"},
         {Form("h_%s_2mu4_given_mu4_sign1",      var.c_str()), "2mu4, same sign"},
         {Form("h_%s_mu4_mu4noL1_given_mu4_sign2",var.c_str()), "mu4_mu4noL1, opposite sign"},
         {Form("h_%s_2mu4_given_mu4_sign2",      var.c_str()), "2mu4, opposite sign"}
      };

      TCanvas* c = new TCanvas(Form("c2d_%s",var.c_str()),var.c_str(),1100,900);
      c->Divide(2,2);

      for(int i=0;i<4;++i){
         c->cd(i+1);
         gPad->SetRightMargin(0.15);   // leave room for color bar
         applyLog(xy);                 // z-axis never log
         TH2* h = getObj<TH2>(names[i].tag);
         if(h){
            h->SetTitle(names[i].title.c_str());
            h->Draw("COLZ");
         }
      }
      c->SaveAs(Form("%strig_effcy_plots/%s_trig_effcy.png", data_dir.c_str(), var.c_str()));

      // second 1×2 canvas (signal-selection, sign2 only) -------------------
      std::string h1n = Form("h_%s_mu4_mu4noL1_given_mu4_w_sig_sel_sign2",var.c_str());
      std::string h2n = Form("h_%s_2mu4_given_mu4_w_sig_sel_sign2",      var.c_str());
      TH2* h1 = getObj<TH2>(h1n);
      TH2* h2 = getObj<TH2>(h2n);

      TCanvas* cB = new TCanvas(Form("c2dB_%s",var.c_str()),var.c_str(),900,450);
      cB->Divide(2,1);

      cB->cd(1); gPad->SetRightMargin(0.15); applyLog(xy);
      if(h1){ h1->SetTitle("mu4_mu4noL1, signal selection"); h1->Draw("COLZ"); }

      cB->cd(2); gPad->SetRightMargin(0.15); applyLog(xy);
      if(h2){ h2->SetTitle("2mu4, signal selection");        h2->Draw("COLZ"); }

      cB->SaveAs(Form("%strig_effcy_plots/%s_trig_effcy_w_sig_sel.png", data_dir.c_str(), var.c_str()));
   }
};



void trig_effcy_plot(){
   std::vector<std::string> var1Ds = {"Deta", "Deta_zoomin", "Dphi", "Dphi_zoomin", "DR", "DR_zoomin", "pt2nd", "minv_zoomin", "pair_pt_log"};
   std::vector<std::string> var2Ds = {"pt2nd_vs_q_eta_2nd", "pair_eta_vs_pair_pT", "Deta_Dphi", "eta1_eta2", "eta_avg_Deta", "eta_avg_Dphi", "minv_pair_pt_log"};
   std::map<std::string,std::pair<bool,bool>> logaxes = {
      {"pt2nd",{true,false}},
      {"pair_pt_log",{true,false}},
      {"pt2nd_vs_q_eta_2nd",{false,true}},
      {"pair_eta_vs_pair_pT",{true,false}},
      {"minv_pair_pt_log",{true,true}},
   };

   TrigEffPlotter plotter("/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_2024/histograms_real_pairs_pp_2024_single_mu4.root",var1Ds,logaxes,var2Ds);
   plotter.Run();
}

