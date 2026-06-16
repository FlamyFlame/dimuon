// =============================================================================
// build_run2_reco_eff_placeholder.C
//
// PLACEHOLDER reco-efficiency builder. Digitizes the Run 2 *single-muon*
// reconstruction efficiencies and writes them to a ROOT file for use as a
// temporary reco-eff correction in the Run 3 crossx pipeline, until the proper
// Run 3 *pair* efficiency epsilon_reco(pair pT, pair eta, dR) from Pythia
// fullsim + HIJING overlay becomes available (see
// docs/tracking/reco_eff_placeholder_run2.md, task_05).
//
// Sources (Medium muons):
//   PbPb : ATL-COM-PHYS-2021-1094 (Run 2 dimuon note), Appendix F.2 -- single-muon
//          eps_reco(pT, q*eta) per centrality interval, HIJING overlay, 5.02 TeV.
//   pp   : HION-2019-58 / arXiv:2109.00411 (Run 2 HF-muon R_AA note), Fig. 31 --
//          data-driven Medium-muon eps_reco^pp(pT) in barrel (0.1<|eta|<1.05)
//          and endcap (1.3<|eta|<2.1). This is the source the dimuon note itself
//          cites for its pp efficiencies.
//
// Values are eyeball-digitized anchor points; precision is intentionally modest
// (this is a flagged placeholder). The reproduction plots are cross-checked
// against the original PDF panels (/review-plot).
//
// Run:  root -l -b -q build_run2_reco_eff_placeholder.C
// =============================================================================

#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>
#include <string>
#include <vector>
#include <utility>
#include <iostream>

#include "../../Utilities/proj_range_to_suffix.cxx"  // pairToSuffix definition

// q*eta slices -- MUST match CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap_run2
static const std::vector<std::pair<float,float>> kQEtaSlices = {
    {-2.4f,-2.0f},{-2.0f,-1.5f},{-1.5f,-1.0f},{-1.0f,-0.5f},{-0.5f,0.5f},
    {0.5f,1.0f},{1.0f,1.5f},{1.5f,2.0f},{2.0f,2.4f}
};

// F.2 centrality intervals (note appendix), 7 of them
static const std::vector<std::pair<int,int>> kCtrIntervals = {
    {0,10},{10,20},{20,30},{30,40},{40,50},{50,60},{60,80}
};

// pT anchor points [GeV]
static const std::vector<double> kPt = {4.0,4.5,5.0,6.0,7.0,8.0,10.0,13.0,16.0,19.0};
static const int kNpt = (int)kPt.size();

// Base single-muon eps_reco for the 0-10% interval, per q*eta slice [9][10].
// Digitized from F.2 (Medium muons). Forward/backward slices ~flat; central
// slices show the low-pT turn-on.
static const double kEffBase[9][10] = {
    {0.78,0.79,0.80,0.81,0.81,0.82,0.82,0.82,0.83,0.83}, // [-2.4,-2.0]
    {0.76,0.77,0.78,0.78,0.79,0.79,0.79,0.79,0.80,0.80}, // [-2.0,-1.5]
    {0.62,0.70,0.78,0.85,0.87,0.88,0.89,0.90,0.90,0.91}, // [-1.5,-1.0]
    {0.80,0.85,0.88,0.90,0.91,0.92,0.93,0.93,0.94,0.95}, // [-1.0,-0.5]
    {0.55,0.62,0.68,0.78,0.82,0.84,0.85,0.86,0.86,0.87}, // [-0.5, 0.5]
    {0.45,0.58,0.68,0.82,0.88,0.90,0.92,0.93,0.93,0.94}, // [ 0.5, 1.0]
    {0.45,0.58,0.70,0.84,0.87,0.88,0.89,0.89,0.90,0.90}, // [ 1.0, 1.5]
    {0.76,0.77,0.78,0.79,0.80,0.80,0.81,0.81,0.82,0.82}, // [ 1.5, 2.0]
    {0.78,0.79,0.80,0.81,0.81,0.82,0.82,0.82,0.83,0.84}  // [ 2.0, 2.4]
};

// Mild central->peripheral additive increase (note: "small but systematic,
// increasing from central to peripheral"). Applied per centrality interval.
static const double kCtrBump[7] = {0.00,0.01,0.02,0.03,0.035,0.04,0.05};

// pp Medium-muon eps_reco (HF R_AA Fig. 31), barrel and endcap.
static const double kEffPpBarrel[10] = {0.77,0.85,0.92,0.95,0.96,0.96,0.97,0.98,0.98,0.98};
static const double kEffPpEndcap[10] = {0.60,0.72,0.85,0.92,0.93,0.94,0.95,0.96,0.96,0.97};

static std::string ctrKey(int lo, int hi){ return "_ctr" + std::to_string(lo) + "_" + std::to_string(hi); }

static TGraph* makeGraph(const std::string& name, const double* eff, double bump){
    TGraph* g = new TGraph(kNpt);
    g->SetName(name.c_str());
    for (int i=0;i<kNpt;++i){
        double v = eff[i] + bump;
        if (v > 1.0) v = 1.0;
        if (v < 0.01) v = 0.01;
        g->SetPoint(i, kPt[i], v);
    }
    return g;
}

void build_run2_reco_eff_placeholder(){
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);

    const std::string out_root = "../../EfficiencyCorrs/EffFiles/run2_reco_eff_placeholder.root";
    const std::string plot_dir = "../../plots/reco_effcy/placeholder/";
    gSystem->mkdir(plot_dir.c_str(), kTRUE);

    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");

    // --- PbPb: TGraph per (centrality interval, q*eta slice) ---
    for (size_t ic=0; ic<kCtrIntervals.size(); ++ic){
        const int lo = kCtrIntervals[ic].first, hi = kCtrIntervals[ic].second;
        for (size_t is=0; is<kQEtaSlices.size(); ++is){
            std::string nm = "gr_reco_eff_medium_pbpb" + ctrKey(lo,hi)
                           + "_q_eta_" + pairToSuffix(kQEtaSlices[is]);
            TGraph* g = makeGraph(nm, kEffBase[is], kCtrBump[ic]);
            g->Write();
        }
    }
    // --- pp: barrel + endcap ---
    makeGraph("gr_reco_eff_medium_pp_barrel", kEffPpBarrel, 0.0)->Write();
    makeGraph("gr_reco_eff_medium_pp_endcap", kEffPpEndcap, 0.0)->Write();

    fout->Close();
    std::cout << "[build] wrote " << out_root << std::endl;

    // --- Reproduction plots: one 3x3 canvas per centrality (F.2 layout) ---
    TFile* fin = TFile::Open(out_root.c_str(), "READ");
    for (size_t ic=0; ic<kCtrIntervals.size(); ++ic){
        const int lo = kCtrIntervals[ic].first, hi = kCtrIntervals[ic].second;
        TCanvas c("c","c",1200,1000);
        c.Divide(3,3);
        for (size_t is=0; is<kQEtaSlices.size(); ++is){
            c.cd(is+1);
            gPad->SetGridx(); gPad->SetGridy();
            std::string nm = "gr_reco_eff_medium_pbpb" + ctrKey(lo,hi)
                           + "_q_eta_" + pairToSuffix(kQEtaSlices[is]);
            TGraph* g = (TGraph*)fin->Get(nm.c_str());
            g->SetMarkerStyle(24); g->SetMarkerSize(0.9);
            g->SetTitle("");
            g->GetXaxis()->SetLimits(4,20);
            g->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
            g->GetYaxis()->SetTitle("Efficiency");
            g->GetYaxis()->SetRangeUser(0,1.6);
            g->Draw("APL");
            TLatex t; t.SetNDC(); t.SetTextSize(0.05);
            t.DrawLatex(0.18,0.85, Form("%.1f<q*#eta<%.1f, Medium-#mu",
                        kQEtaSlices[is].first, kQEtaSlices[is].second));
            t.DrawLatex(0.18,0.79, Form("%d-%d%%  (PLACEHOLDER)", lo, hi));
        }
        std::string out = plot_dir + Form("reco_eff_placeholder_pbpb_ctr%d_%d.png", lo, hi);
        c.SaveAs(out.c_str());
    }
    // --- pp reproduction canvas ---
    {
        TCanvas c("cpp","cpp",1000,500);
        c.Divide(2,1);
        const char* names[2] = {"gr_reco_eff_medium_pp_barrel","gr_reco_eff_medium_pp_endcap"};
        const char* labs[2]  = {"Medium muon, 0.10<#eta<1.05","Medium muon, 1.30<#eta<2.10"};
        for (int k=0;k<2;++k){
            c.cd(k+1); gPad->SetGridx(); gPad->SetGridy();
            TGraph* g = (TGraph*)fin->Get(names[k]);
            g->SetMarkerStyle(20); g->SetMarkerColor(kBlack); g->SetLineColor(kBlack);
            g->SetTitle(""); g->GetXaxis()->SetLimits(4,20);
            g->GetXaxis()->SetTitle("p_{T} [GeV]"); g->GetYaxis()->SetTitle("Reconstruction efficiency");
            g->GetYaxis()->SetRangeUser(0,1.2); g->Draw("APL");
            TLatex t; t.SetNDC(); t.SetTextSize(0.045);
            t.DrawLatex(0.30,0.40, labs[k]);
            t.DrawLatex(0.30,0.33, "pp placeholder (HF R_{AA} Fig.31)");
        }
        c.SaveAs((plot_dir + "reco_eff_placeholder_pp.png").c_str());
    }
    fin->Close();
    std::cout << "[build] wrote reproduction plots to " << plot_dir << std::endl;
}
