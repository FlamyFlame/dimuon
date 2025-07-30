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
#include "TGraphAsymmErrors.h"


class TrigEffPlotter {
public:
   /// ctor:  file name  |  1-D vars               |  optional log settings
   ///                   |  2-D vars (optional)
   // TrigEffPlotter(const std::string&  fileName,
   TrigEffPlotter(int  runYear,
                  const std::vector<std::string>& vars1D,
                  bool plot_mu4_mu4noL1_excl_ = false,
                  bool plot_sepr_ = false,
                  bool plot_good_accept_ = false,
                  bool good_accept_means_bad_ = false,
                  const std::map<std::string,std::pair<bool,bool>>& logopt_ = {},
                  const std::vector<std::string>& vars2D_ = {},
                  const std::vector<std::string>& vars2D_for_profile = {})
      : fFile(nullptr), plot_mu4_mu4noL1_excl(plot_mu4_mu4noL1_excl_), plot_sepr(plot_sepr_), plot_good_accept(plot_good_accept_), good_accept_means_bad(good_accept_means_bad_), var1Ds(vars1D), var2Ds(vars2D_), var2Ds_for_profile(vars2D_for_profile), logopt(logopt_)
   {
      isRun2pp = (runYear % 2000 == 17);
      switch (runYear % 2000){
      case 17:
         data_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_run2/";
         infile_name = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_run2/histograms_real_pairs_pp_2017_single_mu4.root";
         break;
      case 24:
         data_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_2024/";
         infile_name = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_2024/histograms_real_pairs_pp_2024_single_mu4.root";
         break;
      default:
         std::cout << "first argument to constructor, runYear, must be 17 or 24 for pp data!" << std::endl;
         throw std::exception();
      }
      initialize();
   }

   /// loop over all var1Ds (1-D and 2-D)
   void Run() {
      if(!fFile || fFile->IsZombie()) return;
      gStyle->SetOptStat(0);
      for(const auto& v : var1Ds) plotOne(v);
      for(const auto& v : var2Ds )  plot2DOne(v);
      for(const auto& v : var2Ds_for_profile )  plotProfileOne(v);
   }

private:
   bool isRun2pp;
   std::string data_dir;
   std::string infile_name;

   bool plot_mu4_mu4noL1_excl;
   bool plot_sepr;
   bool plot_good_accept;
   bool good_accept_means_bad;

   std::string sepr_suffix;

   TFile*                          fFile;
   std::vector<std::string>        var1Ds;   // 1-D
   std::vector<std::string>        var2Ds;      // 2-D  (NEW)
   std::vector<std::string>        var2Ds_for_profile;      // 2-D  (NEW)
   std::map<std::string,std::pair<bool,bool>> logopt;
   std::vector<std::string> vars_for_good_accept = {"Deta", "Deta_zoomin", "Dphi", "Dphi_zoomin", "DR", "DR_zoomin", "minv_zoomin", "pair_pt_log"};
   
   void initialize() {
      sepr_suffix = plot_sepr? "_sepr" : "";
      if (plot_good_accept){
         if (good_accept_means_bad) sepr_suffix = "_bad_accept";
         else                       sepr_suffix = "_good_accept";
      }

      if (plot_sepr && plot_good_accept){
         std::cerr << "ERROR: plot_sepr and plot_good_accept CANNOT both be true! Return without plotting!" << std::endl;
         throw std::exception();
      }
      fFile = TFile::Open(infile_name.c_str(),"READ");
      if(!fFile || fFile->IsZombie())
         std::cerr<<"[TrigEffPlotter] ERROR: cannot open "<<infile_name<<std::endl;
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
   void setStyle(TGraph* h,int col,int m) {
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

   void setTGraphBinning(TGraphAsymmErrors* g, TH1* h_pass, TH1* h_total){
      int n = h_pass->GetNbinsX();
      int ip = 0;  // index of the valid point in the TGraph

      for (int i = 1; i <= n; ++i) {
          double denom = h_total->GetBinContent(i);
          if (denom <= 0) continue;  // this bin was skipped by BayesDivide

          double xlow  = h_pass->GetXaxis()->GetBinLowEdge(i);
          double xhigh = h_pass->GetXaxis()->GetBinUpEdge(i);
          double xcenter = 0.5 * (xlow + xhigh);
          double halfwidth = 0.5 * (xhigh - xlow);

          g->SetPoint(ip, xcenter, g->GetY()[ip]);  // keep existing Y
          g->SetPointError(ip, halfwidth, halfwidth,
                               g->GetEYlow()[ip], g->GetEYhigh()[ip]);
          ++ip;
      }
   }

   void adjustLogXRange(TGraphAsymmErrors* g, TH1* h_ref){
      double xmin = h_ref->GetXaxis()->GetBinLowEdge(1);
      double xmax = h_ref->GetXaxis()->GetBinLowEdge(h_ref->GetNbinsX() + 1);
      g->GetXaxis()->SetRangeUser(xmin, xmax);
   }

   // ---------- 1-D plotting (unchanged from previous reply) ---------------
   void plotOne(const std::string& var){
      const auto xy = xyFor(var);

      bool var_available_good_accept = (std::find(vars_for_good_accept.begin(), vars_for_good_accept.end(), var) != vars_for_good_accept.end()); // variable is suitable for good-acceptance comparison
      if (plot_good_accept && !var_available_good_accept) return;

      // ------------------------------------------------------------------------------------------------
      // -------------------------- FULL RANGE (NO SIGNAL SELECTION) ------------------------------------
      // ------------------------------------------------------------------------------------------------

      // ─────────────────────── hist names ───────────────────────
      std::string h_name_ss_ref        = Form("h_%s_mu4_sign1",            var.c_str());
      std::string h_name_ss_1          = Form("h_%s_mu4_mu4noL1_sign1",    var.c_str());
      std::string h_name_ss_2          = Form("h_%s_2mu4_sign1",           var.c_str());

      std::string h_name_op_ref        = Form("h_%s_mu4_sign2",            var.c_str());
      std::string h_name_op_1          = Form("h_%s_mu4_mu4noL1_sign2",    var.c_str());
      std::string h_name_op_2          = Form("h_%s_2mu4_sign2",           var.c_str());

      // “sepr” version (only needed when plot_sepr == true)
      std::string h_name_ss_ref_sepr   = Form("h_%s_mu4_sepr_sign1",       var.c_str());
      std::string h_name_ss_1_sepr     = Form("h_%s_mu4_mu4noL1_sepr_sign1",var.c_str());
      std::string h_name_ss_2_sepr     = Form("h_%s_2mu4_sepr_sign1",      var.c_str());

      std::string h_name_op_ref_sepr   = Form("h_%s_mu4_sepr_sign2",       var.c_str());
      std::string h_name_op_1_sepr     = Form("h_%s_mu4_mu4noL1_sepr_sign2",var.c_str());
      std::string h_name_op_2_sepr     = Form("h_%s_2mu4_sepr_sign2",      var.c_str());

      // “sepr” version (only needed when plot_sepr == true)
      std::string h_name_ss_ref_good_accept   = Form("h_%s_mu4_good_accept_sign1",       var.c_str());
      std::string h_name_ss_1_good_accept     = Form("h_%s_mu4_mu4noL1_good_accept_sign1",var.c_str());
      std::string h_name_ss_2_good_accept     = Form("h_%s_2mu4_good_accept_sign1",      var.c_str());

      std::string h_name_op_ref_good_accept   = Form("h_%s_mu4_good_accept_sign2",       var.c_str());
      std::string h_name_op_1_good_accept     = Form("h_%s_mu4_mu4noL1_good_accept_sign2",var.c_str());
      std::string h_name_op_2_good_accept     = Form("h_%s_2mu4_good_accept_sign2",      var.c_str());

      // ─────────────────────── get histograms ───────────────────────
      TH1 *h0_ss = getObj<TH1>(h_name_ss_ref), *h1_ss = getObj<TH1>(h_name_ss_1),
          *h2_ss = getObj<TH1>(h_name_ss_2),
          *h0_op = getObj<TH1>(h_name_op_ref), *h1_op = getObj<TH1>(h_name_op_1),
          *h2_op = getObj<TH1>(h_name_op_2);

      TH1 *h0_ss_sepr=nullptr,*h1_ss_sepr=nullptr,*h2_ss_sepr=nullptr,
          *h0_op_sepr=nullptr,*h1_op_sepr=nullptr,*h2_op_sepr=nullptr;

      TH1 *h0_ss_good_accept=nullptr,*h1_ss_good_accept=nullptr,*h2_ss_good_accept=nullptr,
          *h0_op_good_accept=nullptr,*h1_op_good_accept=nullptr,*h2_op_good_accept=nullptr;

      if (plot_sepr) {
         h0_ss_sepr = getObj<TH1>(h_name_ss_ref_sepr);
         h1_ss_sepr = getObj<TH1>(h_name_ss_1_sepr);
         h2_ss_sepr = getObj<TH1>(h_name_ss_2_sepr);

         h0_op_sepr = getObj<TH1>(h_name_op_ref_sepr);
         h1_op_sepr = getObj<TH1>(h_name_op_1_sepr);
         h2_op_sepr = getObj<TH1>(h_name_op_2_sepr);
      }

      if (plot_good_accept) {
         h0_ss_good_accept = getObj<TH1>(h_name_ss_ref_good_accept);
         h1_ss_good_accept = getObj<TH1>(h_name_ss_1_good_accept);
         h2_ss_good_accept = getObj<TH1>(h_name_ss_2_good_accept);

         h0_op_good_accept = getObj<TH1>(h_name_op_ref_good_accept);
         h1_op_good_accept = getObj<TH1>(h_name_op_1_good_accept);
         h2_op_good_accept = getObj<TH1>(h_name_op_2_good_accept);
      }

      // (basic null checks for the mandatory histograms omitted for brevity)

      // ─────────────────────── efficiencies ───────────────────────
      auto mkGR = [this](TH1* num, TH1* den) {
         auto g = new TGraphAsymmErrors();
         g->BayesDivide(num, den);
         // setTGraphBinning(g,num,den);
         return g;
      };

      TGraphAsymmErrors *gr1_ss      = mkGR(h1_ss,h0_ss),
                        *gr2_ss      = mkGR(h2_ss,h0_ss),
                        *gr1_op      = mkGR(h1_op,h0_op),
                        *gr2_op      = mkGR(h2_op,h0_op),
                        *gr1_ss_sepr = nullptr, *gr2_ss_sepr = nullptr,
                        *gr1_op_sepr = nullptr, *gr2_op_sepr = nullptr,
                        *gr1_ss_good_accept = nullptr, *gr2_ss_good_accept = nullptr,
                        *gr1_op_good_accept = nullptr, *gr2_op_good_accept = nullptr;

      if (plot_sepr) {
         if (!isRun2pp){
            gr1_ss_sepr = mkGR(h1_ss_sepr,h0_ss_sepr);
            gr1_op_sepr = mkGR(h1_op_sepr,h0_op_sepr);
         }  

         gr2_ss_sepr    = mkGR(h2_ss_sepr,h0_ss_sepr);
         gr2_op_sepr    = mkGR(h2_op_sepr,h0_op_sepr);
      }

      if (plot_good_accept) {
         if (!isRun2pp){
            gr1_ss_good_accept = mkGR(h1_ss_good_accept,h0_ss_good_accept);
            gr1_op_good_accept = mkGR(h1_op_good_accept,h0_op_good_accept);
         }  

         gr2_ss_good_accept    = mkGR(h2_ss_good_accept,h0_ss_good_accept);
         gr2_op_good_accept    = mkGR(h2_op_good_accept,h0_op_good_accept);
      }

      // ─────────────────────── styles ───────────────────────
      auto setStyle = [](TGraphAsymmErrors* g, Color_t col, Style_t m, Width_t w=2, Style_t ls=1){
         if(!g) return;
         g->SetMarkerColor(col); g->SetLineColor(col);
         g->SetMarkerStyle(m);   g->SetLineWidth(w); g->SetLineStyle(ls);
      };

      setStyle(gr2_ss,      kAzure+2, 24);                           // solid
      setStyle(gr1_ss,      kRed+1,   20);
      setStyle(gr2_op,      kAzure+2, 24);
      setStyle(gr1_op,      kRed+1,   20);

      setStyle(gr2_ss_sepr, kAzure+1, 25, 2, 2);                     // dashed, open marker
      setStyle(gr1_ss_sepr, kRed+2,   24, 2, 2);
      setStyle(gr2_op_sepr, kAzure+1, 25, 2, 2);
      setStyle(gr1_op_sepr, kRed+2,   24, 2, 2);

      setStyle(gr2_ss_good_accept, kAzure+1, 25, 2, 2);                     // dashed, open marker
      setStyle(gr1_ss_good_accept, kRed+2,   24, 2, 2);
      setStyle(gr2_op_good_accept, kAzure+1, 25, 2, 2);
      setStyle(gr1_op_good_accept, kRed+2,   24, 2, 2);

      // ─────────────────────── canvas ───────────────────────
      TCanvas* c1 = new TCanvas(Form("c1_%s",var.c_str()),var.c_str(),1200,500);
      c1->Divide(2,1);

      // helper lambdas -----------------------------------------------------------
      auto drawPad = [&](TPad* pad,
                         TGraphAsymmErrors* g_ref,             // always drawn
                         TGraphAsymmErrors* g_alt,             // drawn if !isRun2pp
                         TGraphAsymmErrors* g_ref_sepr,        // drawn if plot_sepr
                         TGraphAsymmErrors* g_alt_sepr,        // drawn if (plot_sepr && !isRun2pp)
                         TGraphAsymmErrors* g_ref_good_accept,        // drawn if plot_good_accept
                         TGraphAsymmErrors* g_alt_good_accept,        // drawn if (plot_good_accept && !isRun2pp)
                         TH1* h_ref,                           // histogram to help set pad properties
                         const char* title){
         pad->cd(); applyLog(xy);

         // 1) reference plain graph sets the frame
         g_ref->Draw("APL");
         g_ref->SetTitle(title);
         g_ref->GetXaxis()->SetTitle(h_ref->GetXaxis()->GetTitle());
         g_ref->GetYaxis()->SetTitle("#epsilon");
         g_ref->GetYaxis()->SetRangeUser(0.,1.);
         if (xy.first) adjustLogXRange(g_ref, h_ref);

         // 2) overlay others
         if (!isRun2pp) g_alt->Draw("PL SAME");
         if (plot_sepr) {
            g_ref_sepr->Draw("PL SAME");
            if (!isRun2pp) g_alt_sepr->Draw("PL SAME");
         }
         if (plot_good_accept) {
            g_ref_good_accept->Draw("PL SAME");
            if (!isRun2pp) g_alt_good_accept->Draw("PL SAME");
         }

         // 3) legend
         int nLines = 1                               // g_ref
                    + (!isRun2pp ? 1 : 0)             // g_alt
                    + (plot_sepr ? 1 : 0)             // g_ref_sepr
                    + (plot_sepr && !isRun2pp ? 1:0) // g_alt_sepr
                    + (plot_good_accept ? 1 : 0)             // g_ref_good_accept
                    + (plot_good_accept && !isRun2pp ? 1:0); // g_alt_good_accept
         double leg_ymax = 0.17 + 0.065 * nLines;
         
         double right_move, up_move = 0.;
         if (var == "Deta" || var =="DR"){
            right_move = -0.05;
         }
         if (var =="DR"){
            right_move = -0.12;
         }
         if (var == "Deta_zoomin" || var == "Dphi" || var == "Dphi_zoomin" || var =="minv_zoomin"){
            right_move = 0.04;
            // up_move = -0.02;
         }
         TLegend* L = new TLegend(0.52 + right_move, 0.14 + up_move, 0.84 + right_move, leg_ymax + up_move);
         L->SetBorderSize(0);
         if (!isRun2pp) {
            L->AddEntry(g_alt,        "mu4_mu4noL1",     "lep");
            if (plot_sepr) L->AddEntry(g_alt_sepr,"mu4_mu4noL1 (sepr)","lep");
            if (plot_good_accept){
               if (good_accept_means_bad) L->AddEntry(g_alt_good_accept,"mu4_mu4noL1 (bad accept)","lep");
               else                       L->AddEntry(g_alt_good_accept,"mu4_mu4noL1 (good accept)","lep");
            }
         }
         L->AddEntry(g_ref,          "2mu4",               "lep");
         if (plot_sepr) L->AddEntry(g_ref_sepr,"2mu4 (sepr)",          "lep");
         if (plot_good_accept){
            if (good_accept_means_bad) L->AddEntry(g_ref_good_accept,"2mu4 (bad accept)",          "lep");
            else                       L->AddEntry(g_ref_good_accept,"2mu4 (good accept)",          "lep");
         }
         L->Draw();
      };

      // ───────────── same-sign pad ─────────────
      drawPad((TPad*)c1->cd(1),
              gr2_ss, gr1_ss, gr2_ss_sepr, gr1_ss_sepr, gr2_ss_good_accept, gr1_ss_good_accept, h1_ss,
              "same sign");

      // ───────────── opposite-sign pad ─────────────
      drawPad((TPad*)c1->cd(2),
              gr2_op, gr1_op, gr2_op_sepr, gr1_op_sepr, gr2_op_good_accept, gr1_op_good_accept, h1_op,
              "opposite sign");

      c1->SaveAs(Form("%strig_effcy_plots%s/%s_trig_effcy%s.png", data_dir.c_str(), sepr_suffix.c_str(), var.c_str(), sepr_suffix.c_str()));

      // ------------------------------------------------------------------------------------------------
      // --------------------------- WITH SINGLE-B SIGNAL SELECTION -------------------------------------
      // ------------------------------------------------------------------------------------------------

      if (plot_sepr || plot_good_accept) return; // return without with-single-b-signal-selection efficiency plots

      // -------------------------- get histograms ------------------------------------

      std::string h_name_sel_0 = Form("h_%s_mu4_w_single_b_sig_sel",        var.c_str());
      std::string h_name_sel_1 = Form("h_%s_mu4_mu4noL1_w_single_b_sig_sel",var.c_str());
      std::string h_name_sel_2 = Form("h_%s_2mu4_w_single_b_sig_sel",       var.c_str());
      TH1* h0_sel = getObj<TH1>(h_name_sel_0);
      TH1* h1_sel = getObj<TH1>(h_name_sel_1);
      TH1* h2_sel = getObj<TH1>(h_name_sel_2);

      if (!h0_sel){ cout << "Function h_" + var + "_mu4_w_single_b_sig_sel NOT FOUND!"; return;}
      if (!h1_sel){ cout << "Function h_" + var + "_mu4_mu4noL1_w_single_b_sig_sel NOT FOUND!"; return;}
      if (!h2_sel){ cout << "Function h_" + var + "_2mu4_w_single_b_sig_sel NOT FOUND!"; return;}

      // -------------------------- calculate efficiencies ------------------------------------

      TGraphAsymmErrors *gr1_sel      = mkGR(h1_sel,h0_sel),
                        *gr2_sel      = mkGR(h2_sel,h0_sel);
      
      // -------------------------- draw canvas w/ 1 pad + save ------------------------------------
      TCanvas* c2 = new TCanvas(Form("c2_%s",var.c_str()),var.c_str(),700,500);
      applyLog(xy);
      setStyle(gr1_sel,kRed+1,20); setStyle(gr2_sel,kAzure+2,24);
      gr2_sel->Draw("APL");
      gr2_sel->GetXaxis()->SetTitle(h1_sel->GetXaxis()->GetTitle()); 
      gr2_sel->GetYaxis()->SetTitle("#epsilon"); 
      gr2_sel->GetYaxis()->SetRangeUser(0.,1.);
      if (xy.first){ // if log x
         adjustLogXRange(gr2_sel, h2_sel);
      }
      if (!isRun2pp) gr1_sel->Draw("PL SAME");

      double right_move, up_move = 0.;
      double leg_ymax = isRun2pp? 0.24 : 0.31;
      if (var =="DR_zoomin"){
         right_move = -0.1;
      }

      { TLegend* L=new TLegend(0.52 + right_move,0.17,0.84 + right_move,leg_ymax);
        if (!isRun2pp) L->AddEntry(gr1_sel,"mu4_mu4noL1","lep");
        L->AddEntry(gr2_sel,"2mu4","lep"); 
        L->Draw(); }

      c2->SaveAs(Form("%strig_effcy_plots/%s_trig_effcy_w_sig_sel.png", data_dir.c_str(), var.c_str()));
   }

   // ---------- 2-D plotting -----------------------------------------
   void plot2DOne(const std::string& var)
   {
      if (plot_good_accept) return;

      const auto xy = xyFor(var);

      // ----------------------------------------------------------------------
      // 1)   PASSED  /  TOTAL   lists, conditioned on Run-2 pp
      // ----------------------------------------------------------------------
      struct HN { std::string hname; std::string title; };

      std::vector<HN> h_passed, h_total;

      if (isRun2pp)
      {
         // ─── only the 2μ4 ratios ───────────────────────────────────────────
         h_passed = {
            {Form("h_%s_2mu4%s_sign1", var.c_str(), sepr_suffix.c_str()), "2mu4, same sign"},
            {Form("h_%s_2mu4%s_sign2", var.c_str(), sepr_suffix.c_str()), "2mu4, opposite sign"}
         };
         h_total  = {
            {Form("h_%s_mu4%s_sign1",  var.c_str(), sepr_suffix.c_str()), "mu4, same sign"},
            {Form("h_%s_mu4%s_sign2",  var.c_str(), sepr_suffix.c_str()), "mu4, opposite sign"}
         };
      }
      else
      {
         // ─── full heavy-ion set (2 × 2) ────────────────────────────────────
         h_passed = {
            {Form("h_%s_mu4_mu4noL1%s_sign1", var.c_str(), sepr_suffix.c_str()), "mu4_mu4noL1, same sign"},
            {Form("h_%s_mu4_mu4noL1%s_sign2", var.c_str(), sepr_suffix.c_str()), "mu4_mu4noL1, opposite sign"},
            {Form("h_%s_2mu4%s_sign1",        var.c_str(), sepr_suffix.c_str()), "2mu4, same sign"},
            {Form("h_%s_2mu4%s_sign2",        var.c_str(), sepr_suffix.c_str()), "2mu4, opposite sign"}
         };
         h_total  = {
            {Form("h_%s_mu4%s_sign1",         var.c_str(), sepr_suffix.c_str()), "mu4, same sign"},
            {Form("h_%s_mu4%s_sign2",         var.c_str(), sepr_suffix.c_str()), "mu4, opposite sign"}
         };
      }

      // ----------------------------------------------------------------------
      // 2)   FIRST CANVAS (global efficiencies)
      // ----------------------------------------------------------------------
      const int nPads = static_cast<int>(h_passed.size());      // 2 or 4
      const int nCols = 2;
      const int nRows = (isRun2pp ? 1 : 2);

      TCanvas* c = new TCanvas(Form("c2d_%s",var.c_str()),var.c_str(),
                               1100, (isRun2pp ? 500 : 900));
      c->Divide(nCols, nRows);

      for (int i = 0; i < nPads; ++i)
      {
         c->cd(i + 1);
         gPad->SetRightMargin(0.15);
         applyLog(xy);

         TH2* h_pass = getObj<TH2>(h_passed[i].hname);
         TH2* h_tot  = getObj<TH2>(h_total[i % 2].hname);  // 0→same, 1→opp

         if (!h_pass) { std::cout << "Histogram " << h_passed[i].hname
                                  << " NOT FOUND!\n"; return; }
         if (!h_tot)  { std::cout << "Histogram " << h_total[i / 2].hname
                                  << " NOT FOUND!\n"; return; }

         h_pass->Divide(h_tot);
         h_pass->SetTitle(h_passed[i].title.c_str());
         h_pass->Draw("COLZ");
      }

      c->SaveAs(Form("%strig_effcy_plots%s/%s_trig_effcy%s.png", data_dir.c_str(), sepr_suffix.c_str(), var.c_str(), sepr_suffix.c_str()));


      // ----------------------------------------------------------------------
      // 3)   SECOND CANVAS (signal-selection, sign 2 only)
      // ----------------------------------------------------------------------

      if (plot_sepr) return;

      std::string h2d_name_mu4  = Form("h_%s_mu4_w_single_b_sig_sel",   var.c_str());
      std::string h2d_name_2mu4 = Form("h_%s_2mu4_w_single_b_sig_sel", var.c_str());
      std::string h2d_name_mu4noL1 = Form("h_%s_mu4_mu4noL1_w_single_b_sig_sel",
                                          var.c_str());

      TH2* h_mu4   = getObj<TH2>(h2d_name_mu4);
      TH2* h_2mu4  = getObj<TH2>(h2d_name_2mu4);
      TH2* h_mu4noL1 = (isRun2pp ? nullptr :
                        getObj<TH2>(h2d_name_mu4noL1));

      // mandatory histograms
      if (!h_mu4 ) { std::cout << "Histogram " << h2d_name_mu4  << " NOT FOUND!\n"; return; }
      if (!h_2mu4) { std::cout << "Histogram " << h2d_name_2mu4 << " NOT FOUND!\n"; return; }
      if (!isRun2pp && !h_mu4noL1) {
          std::cout << "Histogram " << h2d_name_mu4noL1 << " NOT FOUND!\n"; return;
      }

      if (isRun2pp)
      {
         // ─── 1 × 1 canvas ──────────────────────────────────────────────────
         TCanvas* cB = new TCanvas(Form("c2dB_%s",var.c_str()),var.c_str(),500,450);
         cB->cd(); gPad->SetRightMargin(0.15); applyLog(xy);

         h_2mu4->Divide(h_mu4);
         h_2mu4->SetTitle("2mu4, signal selection");
         h_2mu4->Draw("COLZ");

         cB->SaveAs(Form("%strig_effcy_plots/%s_trig_effcy_w_sig_sel.png",
                         data_dir.c_str(), var.c_str()));
      }
      else
      {
         // ─── 1 × 2 canvas (heavy-ion / Run-3) ──────────────────────────────
         TCanvas* cB = new TCanvas(Form("c2dB_%s",var.c_str()),var.c_str(),900,450);
         cB->Divide(2,1);

         // pad 1 : mu4_mu4noL1 / mu4
         cB->cd(1); gPad->SetRightMargin(0.15); applyLog(xy);
         h_mu4noL1->Divide(h_mu4);
         h_mu4noL1->SetTitle("mu4_mu4noL1, signal selection");
         h_mu4noL1->Draw("COLZ");

         // pad 2 : 2mu4 / mu4
         cB->cd(2); gPad->SetRightMargin(0.15); applyLog(xy);
         h_2mu4->Divide(h_mu4);
         h_2mu4->SetTitle("2mu4, signal selection");
         h_2mu4->Draw("COLZ");

         cB->SaveAs(Form("%strig_effcy_plots/%s_trig_effcy_w_sig_sel.png",
                         data_dir.c_str(), var.c_str()));
      }
   }
   // ─────────────────────────────────────────────────────────────────────────────
   //   Plot the 2-D “profile” histograms (no division)
   //   * canvas layout mirrors the efficiency version
   //   * saves   …/trig_effcy_plots…/<var>_profile[ _sepr].png
   // ─────────────────────────────────────────────────────────────────────────────
   void plotProfileOne(const std::string& var)
   {
       if (plot_good_accept) return;

       const auto xy = xyFor(var);

       // ------------------------------------------------------------------ 1)
       // Histograms to draw in the first canvas (global view)
       // ------------------------------------------------------------------
       struct HN { std::string hname; std::string title; };

       std::vector<HN> h_list;
       h_list = {
           {Form("h_%s_mu4%s_sign1", var.c_str(), sepr_suffix.c_str()),
            "mu4, same sign"},
           {Form("h_%s_mu4%s_sign2", var.c_str(), sepr_suffix.c_str()),
            "mu4, opposite sign"}
       };

       // first canvas geometry
       const int nPads = static_cast<int>(h_list.size());  // always 2 entries
       const int nCols = 2;
       const int nRows = (isRun2pp ? 1 : 2);               // 2×1  or  2×2

       TCanvas* c = new TCanvas(Form("cProf_%s", var.c_str()), var.c_str(),
                                1100, (isRun2pp ? 500 : 900));
       c->Divide(nCols, nRows);

       // book‐keeping: if heavy-ion call wants 4 pads, duplicate the two hists
       if (!isRun2pp) h_list.insert(h_list.end(), h_list.begin(), h_list.end());

       for (int i = 0; i < nCols * nRows; ++i) {
           c->cd(i + 1);
           gPad->SetRightMargin(0.15);
           gPad->SetLogz();
           applyLog(xy);

           if (i < static_cast<int>(h_list.size())) {
               TH2* h = getObj<TH2>(h_list[i].hname);
               if (!h) { std::cout << "Histogram " << h_list[i].hname
                                   << " NOT FOUND!\n"; return; }
               h->SetTitle(h_list[i].title.c_str());
               h->Draw("COLZ");

               TH1* h_pfx = h->ProfileX(Form("%s_pfx",h->GetName()), 1, -1, "s");
               h_pfx->SetLineColor(kRed);
               h_pfx->SetMarkerColor(kRed);
               h_pfx->SetLineWidth(2);
               h_pfx->Draw("E,same");   // E1 = draw points with vertical error bars

           }
       }

       c->SaveAs(Form("%strig_effcy_plots%s/%s_profile%s.png",
                      data_dir.c_str(), sepr_suffix.c_str(),
                      var.c_str(),      sepr_suffix.c_str()));

       // ------------------------------------------------------------------ 2)
       // Second canvas: signal-selection profile (only one histogram exists)
       // ------------------------------------------------------------------
       if (plot_sepr) return;               // same behaviour as efficiency version

       std::string h2d_name_mu4 = Form("h_%s_mu4_w_single_b_sig_sel", var.c_str());
       TH2* h_mu4 = getObj<TH2>(h2d_name_mu4);
       if (!h_mu4) {
           std::cout << "Histogram " << h2d_name_mu4 << " NOT FOUND!\n";
           return;
       }

       // canvas geometry   (Run-2 pp → 1×1 ;  otherwise 1×2, second pad left blank)
       if (isRun2pp) {
           TCanvas* cB = new TCanvas(Form("cProfB_%s", var.c_str()),
                                     var.c_str(), 500, 450);
           cB->cd(); gPad->SetRightMargin(0.15); applyLog(xy);
           gPad->SetLogz();
           h_mu4->SetTitle("mu4, signal selection");
           h_mu4->Draw("COLZ");

           TH1* h_pfx = h_mu4->ProfileX(Form("%s_pfx",h_mu4->GetName()), 1, -1, "s");
           h_pfx->SetLineColor(kRed);
           h_pfx->SetMarkerColor(kRed);
           h_pfx->SetLineWidth(2);
           h_pfx->Draw("E,same");   // E1 = draw points with vertical error bars

           cB->SaveAs(Form("%strig_effcy_plots/%s_profile_w_sig_sel.png",
                           data_dir.c_str(), var.c_str()));
       } else {
           TCanvas* cB = new TCanvas(Form("cProfB_%s", var.c_str()),
                                     var.c_str(), 900, 450);
           cB->Divide(2, 1);

           // pad 1 : mu4  (signal selection)
           cB->cd(1); gPad->SetRightMargin(0.15); applyLog(xy);
           gPad->SetLogz();
           h_mu4->SetTitle("mu4, signal selection");
           h_mu4->Draw("COLZ");

           TH1* h_pfx = h_mu4->ProfileX(Form("%s_pfx",h_mu4->GetName()), 1, -1, "s");
           h_pfx->SetLineColor(kRed);
           h_pfx->SetMarkerColor(kRed);
           h_pfx->SetLineWidth(2);
           h_pfx->Draw("E,same");   // E1 = draw points with vertical error bars

           // pad 2 intentionally left blank for symmetry with efficiency canvas
           cB->cd(2);
           gPad->SetRightMargin(0.15);
           cB->Update();

           cB->SaveAs(Form("%strig_effcy_plots/%s_profile_w_sig_sel.png",
                           data_dir.c_str(), var.c_str()));
       }
   }

};



void trig_effcy_plot(){
   std::vector<std::string> var1Ds = {"Deta", "Deta_zoomin", "Dphi", "Dphi_zoomin", "DR", "DR_zoomin", "pt2nd", "minv_zoomin", "pair_pt_log"};
   std::vector<std::string> var2Ds = {"pt2nd_vs_q_eta_2nd", "pair_eta_vs_pair_pT", "Deta_Dphi", "eta1_eta2", "eta_avg_Deta", "eta_avg_Dphi", "minv_pair_pt_log",
                                      // "Deta_vs_pT_1st", "Deta_zoomin_vs_pT_1st", "Dphi_vs_pT_1st", "Dphi_zoomin_vs_pT_1st", "DR_vs_pT_1st", "DR_zoomin_vs_pT_1st", "minv_zoomin_vs_pT_1st", "pair_pt_log_vs_pT_1st", "pt2nd_vs_pT_1st", 
                                      // "Deta_vs_pair_pT", "Deta_zoomin_vs_pair_pT", "Dphi_vs_pair_pT", "Dphi_zoomin_vs_pair_pT", "DR_vs_pair_pT", "DR_zoomin_vs_pair_pT", "minv_zoomin_vs_pair_pT", "pair_pt_log_vs_pair_pT", "pt2nd_vs_pair_pT",
                                     };
   std::vector<std::string> var2Ds_for_profile = {
      "Deta_vs_pT_1st", "Deta_zoomin_vs_pT_1st", "Dphi_vs_pT_1st", "Dphi_zoomin_vs_pT_1st", "DR_vs_pT_1st", "DR_zoomin_vs_pT_1st", "minv_zoomin_vs_pT_1st", "pair_pt_log_vs_pT_1st", "pt2nd_vs_pT_1st", 
      "Deta_vs_pair_pT", "Deta_zoomin_vs_pair_pT", "Dphi_vs_pair_pT", "Dphi_zoomin_vs_pair_pT", "DR_vs_pair_pT", "DR_zoomin_vs_pair_pT", "minv_zoomin_vs_pair_pT", "pair_pt_log_vs_pair_pT", "pt2nd_vs_pair_pT",
   };


   std::map<std::string,std::pair<bool,bool>> logaxes = {
      {"pt2nd",{true,false}},
      {"pair_pt_log",{true,false}},
      {"pt2nd_vs_q_eta_2nd",{false,true}},
      {"pair_eta_vs_pair_pT",{true,false}},
      {"minv_pair_pt_log",{true,true}},
      {"Deta_vs_pT_1st",{true,false}},
      {"Deta_zoomin_vs_pT_1st",{true,false}},
      {"Dphi_vs_pT_1st",{true,false}},
      {"Dphi_zoomin_vs_pT_1st",{true,false}},
      {"DR_vs_pT_1st",{true,false}},
      {"DR_zoomin_vs_pT_1st",{true,false}},
      {"minv_zoomin_vs_pT_1st",{true,false}},
      {"pair_pt_log_vs_pT_1st",{true,true}},
      {"pt2nd_vs_pT_1st",{true,true}},
      {"Deta_vs_pair_pT",{true,false}},
      {"Deta_zoomin_vs_pair_pT",{true,false}},
      {"Dphi_vs_pair_pT",{true,false}},
      {"Dphi_zoomin_vs_pair_pT",{true,false}},
      {"DR_vs_pair_pT",{true,false}},
      {"DR_zoomin_vs_pair_pT",{true,false}},
      {"minv_zoomin_vs_pair_pT",{true,false}},
      {"pair_pt_log_vs_pair_pT",{true,true}},
      {"pt2nd_vs_pair_p",{true,true}},
   };

   // plot_mu4_mu4noL1_excl_, plot_sepr_, plot_resn_

   // TrigEffPlotter plotter(17,var1Ds,false,false,logaxes,var2Ds,var2Ds_for_profile);
   // TrigEffPlotter plotter(17,var1Ds,false,true,logaxes,var2Ds,var2Ds_for_profile);
   // TrigEffPlotter plotter(24,var1Ds,false,false,false,logaxes,var2Ds,var2Ds_for_profile);
   TrigEffPlotter plotter(24,var1Ds,false,false,true,true,logaxes,var2Ds,var2Ds_for_profile);
   // TrigEffPlotter plotter(24,var1Ds,false,true,logaxes,var2Ds,var2Ds_for_profile);
   plotter.Run();
}

