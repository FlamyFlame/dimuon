// =============================================================================
// build_run2_reco_eff_placeholder.C
//
// PLACEHOLDER reco-efficiency builder. Writes the Run 2 *single-muon*
// reconstruction efficiencies to a ROOT file for use as a temporary reco-eff
// correction in the Run 3 crossx pipeline, until the proper Run 3 *pair*
// efficiency epsilon_reco(pair pT, pair eta, dR) from Pythia fullsim + HIJING
// overlay becomes available (see docs/tracking/reco_eff_placeholder_run2.md,
// task_05).
//
// WORKING POINTS: this builder writes BOTH the Medium and the Tight key sets
// into a single output file, so the consumer (RDFBasedHistFillingData::
// EvaluateSingleMuonRecoEffPlaceholder) can pick the WP-matched keys at runtime
// (nominal WP = Tight; see docs/muon_wp_registry.md §4). Only the muon quality
// bit differs between the two source files; every other convention is identical.
//
// Sources (per WP):
//   PbPb : EXACT Run 2 fits supplied by colleague (the ones used in the Run 2
//          internal note),
//          EfficiencyCorrs/EffFiles/MuonRecoEffcyRun2MC_{medium,tight}.root.
//          Logistic TF1 fits  tf1_eff_fit_cent{C}_eta{E}  of the single-muon
//          eps_reco(pT) per (centrality, q*eta) bin. Both files share the SAME
//          171-fit structure (verified).
//   pp   : HION-2019-58 / arXiv:2109.00411 (Run 2 HF-muon R_AA note), Fig. 31 --
//          data-driven MEDIUM-muon eps_reco^pp(pT) in barrel (0.1<|eta|<1.05)
//          and endcap (1.3<|eta|<2.1).
//          *** PHYSICS GAP (Tight pp): there is NO Tight pp reco-eff source ***
//          The dimuon note App. F.1 (Tight) / F.2 (Medium) are BOTH PbPb-only;
//          the only pp curve anywhere is Fig.31 MEDIUM. We therefore write the
//          Tight pp keys as an INTERIM reuse of the Medium Fig.31 values (NO
//          invented numbers), clearly labelled. The orchestrator/user decides
//          the eventual Tight pp source (documented fallback: most-peripheral
//          60-80% PbPb Tight). See reco_eff_placeholder_run2.md:77-92, memory
//          project_pp_reco_eff_placeholder.
//
// The PbPb fits are written out as TF1 objects (clones of the colleague's fits)
// named tf1_reco_eff_{wp}_pbpb_ctr{lo}_{hi}_q_eta_{suffix}, and the downstream
// RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder evaluates each
// fit AT the exact muon pT (no resampling). pp stays as digitized TGraphs
// named gr_reco_eff_{wp}_pp_{barrel,endcap}.
//
// Run:  root -l -b -q build_run2_reco_eff_placeholder.C
// =============================================================================

#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>
#include <string>
#include <vector>
#include <utility>
#include <iostream>

#include "../../Utilities/proj_range_to_suffix.cxx"  // pairToSuffix definition

// q*eta slices -- MUST match CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap_run2
// AND the colleague's GetEtaBin ordering (eta bin index E == slice index here).
static const std::vector<std::pair<float,float>> kQEtaSlices = {
    {-2.4f,-2.0f},{-2.0f,-1.5f},{-1.5f,-1.0f},{-1.0f,-0.5f},{-0.5f,0.5f},
    {0.5f,1.0f},{1.0f,1.5f},{1.5f,2.0f},{2.0f,2.4f}
};

// F.2 centrality intervals (7) and the colleague file's matching cent index
// (from muon_reco_effcy_run2.txt: original bins 0-10, then coarse bins 11..17 =
// {0,100},{0,10},{10,20},{20,40},{40,60},{0,20},{80,100}).
//   0-10 -> 12, 10-20 -> 13, 20-30 -> 4, 30-40 -> 5, 40-50 -> 6, 50-60 -> 7, 60-80 -> 8
static const std::vector<std::pair<int,int>> kCtrIntervals = {
    {0,10},{10,20},{20,30},{30,40},{40,50},{50,60},{60,80}
};
static const int kCtrColleagueIdx[7] = { 12, 13, 4, 5, 6, 7, 8 };

// pp Medium-muon eps_reco (HF R_AA Fig. 31), barrel and endcap.
// NOTE: MEDIUM only -- no Tight pp source exists; Tight pp reuses these (INTERIM).
static const std::vector<double> kPtPp = {4.0,4.5,5.0,6.0,7.0,8.0,10.0,13.0,16.0,19.0};
static const double kEffPpBarrel[10] = {0.77,0.85,0.92,0.95,0.96,0.96,0.97,0.98,0.98,0.98};
static const double kEffPpEndcap[10] = {0.60,0.72,0.85,0.92,0.93,0.94,0.95,0.96,0.96,0.97};

static std::string ctrKey(int lo, int hi){ return "_ctr" + std::to_string(lo) + "_" + std::to_string(hi); }

// pp graph from eyeball anchor array.
static TGraph* makePpGraph(const std::string& name, const double* eff){
    const int n = (int)kPtPp.size();
    TGraph* g = new TGraph(n);
    g->SetName(name.c_str());
    for (int i=0;i<n;++i){
        double v = eff[i];
        if (v > 1.0) v = 1.0;
        if (v < 0.01) v = 0.01;
        g->SetPoint(i, kPtPp[i], v);
    }
    return g;
}

// Build one WP's key set into the already-open output file. Returns nWritten.
static int buildOneWP(const std::string& wp, TFile* fout, int& nMissingOut){
    const std::string in_root = "../../EfficiencyCorrs/EffFiles/MuonRecoEffcyRun2MC_" + wp + ".root";
    TFile* fin_fits = TFile::Open(in_root.c_str(), "READ");
    if (!fin_fits || fin_fits->IsZombie()){
        std::cerr << "[build] FATAL: cannot open " << in_root << std::endl; nMissingOut = -1; return 0;
    }

    // --- PbPb: clone the colleague's TF1 fit per (centrality interval, q*eta slice) ---
    int nWritten = 0, nMissing = 0;
    for (size_t ic=0; ic<kCtrIntervals.size(); ++ic){
        const int lo = kCtrIntervals[ic].first, hi = kCtrIntervals[ic].second;
        const int cidx = kCtrColleagueIdx[ic];
        for (size_t is=0; is<kQEtaSlices.size(); ++is){
            const int eidx = (int)is;  // slice index == colleague eta bin index
            std::string fitName = Form("tf1_eff_fit_cent%d_eta%d", cidx, eidx);
            TF1* fit = (TF1*)fin_fits->Get(fitName.c_str());
            std::string nm = "tf1_reco_eff_" + wp + "_pbpb" + ctrKey(lo,hi)
                           + "_q_eta_" + pairToSuffix(kQEtaSlices[is]);
            if (!fit){
                std::cerr << "[build] MISSING fit " << fitName << " for " << nm << std::endl;
                ++nMissing; continue;
            }
            TF1* out = (TF1*)fit->Clone(nm.c_str());
            fout->cd();
            out->Write();
            ++nWritten;
        }
    }
    // --- pp: barrel + endcap. Medium = Fig.31; Tight = SAME values (INTERIM, no pp Tight source) ---
    fout->cd();
    makePpGraph("gr_reco_eff_" + wp + "_pp_barrel", kEffPpBarrel)->Write();
    makePpGraph("gr_reco_eff_" + wp + "_pp_endcap", kEffPpEndcap)->Write();

    fin_fits->Close();
    nMissingOut = nMissing;
    return nWritten;
}

// Reproduction plots for one WP (7 PbPb 3x3 canvases + 1 pp), read back from fout path.
static void makePlotsOneWP(const std::string& wp, const std::string& out_root, const std::string& plot_dir){
    const bool pp_interim = (wp == "tight");  // pp Tight is an interim reuse of Medium Fig.31
    TFile* fin = TFile::Open(out_root.c_str(), "READ");
    for (size_t ic=0; ic<kCtrIntervals.size(); ++ic){
        const int lo = kCtrIntervals[ic].first, hi = kCtrIntervals[ic].second;
        TCanvas c("c","c",1200,1000);
        c.Divide(3,3);
        for (size_t is=0; is<kQEtaSlices.size(); ++is){
            c.cd(is+1);
            gPad->SetGridx(); gPad->SetGridy();
            std::string nm = "tf1_reco_eff_" + wp + "_pbpb" + ctrKey(lo,hi)
                           + "_q_eta_" + pairToSuffix(kQEtaSlices[is]);
            TF1* g = (TF1*)fin->Get(nm.c_str());
            if (!g) continue;
            g->SetLineColor(kAzure+1); g->SetLineWidth(2);
            g->SetTitle("");
            g->SetRange(4,19);
            g->GetXaxis()->SetLimits(4,20);
            g->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
            g->GetYaxis()->SetTitle("Efficiency");
            g->GetYaxis()->SetRangeUser(0,1.6);
            g->Draw("L");
            TLatex t; t.SetNDC(); t.SetTextSize(0.05);
            t.DrawLatex(0.18,0.85, Form("%.1f<q*#eta<%.1f, %s-#mu",
                        kQEtaSlices[is].first, kQEtaSlices[is].second, wp.c_str()));
            t.DrawLatex(0.18,0.79, Form("%d-%d%%  (Run 2 fit, PLACEHOLDER)", lo, hi));
        }
        std::string out = plot_dir + Form("reco_eff_placeholder_%s_pbpb_ctr%d_%d.png", wp.c_str(), lo, hi);
        c.SaveAs(out.c_str());
    }
    // --- pp reproduction canvas ---
    {
        TCanvas c("cpp","cpp",1000,500);
        c.Divide(2,1);
        const std::string bn = "gr_reco_eff_" + wp + "_pp_barrel";
        const std::string en = "gr_reco_eff_" + wp + "_pp_endcap";
        const char* names[2] = {bn.c_str(), en.c_str()};
        const char* labs[2]  = {"0.10<#eta<1.05","1.30<#eta<2.10"};
        for (int k=0;k<2;++k){
            c.cd(k+1); gPad->SetGridx(); gPad->SetGridy();
            TGraph* g = (TGraph*)fin->Get(names[k]);
            if (!g) continue;
            g->SetMarkerStyle(20); g->SetMarkerColor(kBlack); g->SetLineColor(kBlack);
            g->SetTitle(""); g->GetXaxis()->SetLimits(4,20);
            g->GetXaxis()->SetTitle("p_{T} [GeV]"); g->GetYaxis()->SetTitle("Reconstruction efficiency");
            g->GetYaxis()->SetRangeUser(0,1.2); g->Draw("APL");
            TLatex t; t.SetNDC(); t.SetTextSize(0.045);
            t.DrawLatex(0.30,0.40, Form("%s muon, %s", wp.c_str(), labs[k]));
            t.DrawLatex(0.30,0.33, "pp placeholder (HF R_{AA} Fig.31)");
            if (pp_interim){
                t.SetTextColor(kRed+1);
                t.DrawLatex(0.30,0.26, "INTERIM: Medium Fig.31 (no Tight pp source)");
            }
        }
        c.SaveAs((plot_dir + "reco_eff_placeholder_" + wp + "_pp.png").c_str());
    }
    fin->Close();
}

void build_run2_reco_eff_placeholder(){
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);

    const std::string out_root = "../../EfficiencyCorrs/EffFiles/run2_reco_eff_placeholder.root";
    const std::string plot_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/reco_effcy_placeholder/";
    gSystem->mkdir(plot_dir.c_str(), kTRUE);

    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()){ std::cerr << "[build] FATAL: cannot create " << out_root << std::endl; return; }

    int totWritten = 0;
    for (const std::string& wp : {std::string("medium"), std::string("tight")}){
        int nMissing = 0;
        int nWritten = buildOneWP(wp, fout, nMissing);
        if (nMissing < 0){ std::cerr << "[build] ABORT: missing input for WP=" << wp << std::endl; fout->Close(); return; }
        totWritten += nWritten;
        std::cout << "[build] WP=" << wp << ": wrote " << nWritten
                  << " PbPb fits (" << nMissing << " missing) + 2 pp graphs" << std::endl;
    }
    fout->Close();
    std::cout << "[build] wrote " << out_root << " (" << totWritten
              << " PbPb fits total across both WPs, + 4 pp graphs)" << std::endl;

    // --- Reproduction plots (per WP) ---
    for (const std::string& wp : {std::string("medium"), std::string("tight")})
        makePlotsOneWP(wp, out_root, plot_dir);
    std::cout << "[build] wrote reproduction plots to " << plot_dir << std::endl;
    std::cout << "[build] NOTE: pp Tight keys are an INTERIM reuse of Medium Fig.31 "
                 "(no Tight pp reco-eff source exists) -- see file header + "
                 "docs/tracking/reco_eff_placeholder_run2.md." << std::endl;
}
