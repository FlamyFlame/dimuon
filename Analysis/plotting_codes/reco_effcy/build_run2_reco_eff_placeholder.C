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
// Sources (Medium muons):
//   PbPb : EXACT Run 2 fits supplied by colleague (the ones used in the Run 2
//          internal note), EfficiencyCorrs/EffFiles/MuonRecoEffcyRun2MC_medium.root.
//          Logistic TF1 fits  tf1_eff_fit_cent{C}_eta{E}  of the single-muon
//          eps_reco(pT) per (centrality, q*eta) bin. These REPLACE the earlier
//          eyeball-digitized F.2 anchor arrays (Step 11, 2026-06-18).
//   pp   : HION-2019-58 / arXiv:2109.00411 (Run 2 HF-muon R_AA note), Fig. 31 --
//          data-driven Medium-muon eps_reco^pp(pT) in barrel (0.1<|eta|<1.05)
//          and endcap (1.3<|eta|<2.1). UNCHANGED (colleague file is PbPb-only).
//
// The PbPb fits are written out as TF1 objects (clones of the colleague's fits)
// named tf1_reco_eff_medium_pbpb_ctr{lo}_{hi}_q_eta_{suffix}, and the downstream
// RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder evaluates each fit
// AT the exact muon pT (no resampling). pp stays as digitized TGraphs.
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

// pp Medium-muon eps_reco (HF R_AA Fig. 31), barrel and endcap (UNCHANGED).
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

void build_run2_reco_eff_placeholder(){
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);

    const std::string in_root  = "../../EfficiencyCorrs/EffFiles/MuonRecoEffcyRun2MC_medium.root";
    const std::string out_root = "../../EfficiencyCorrs/EffFiles/run2_reco_eff_placeholder.root";
    const std::string plot_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/reco_effcy_placeholder/";
    gSystem->mkdir(plot_dir.c_str(), kTRUE);

    TFile* fin_fits = TFile::Open(in_root.c_str(), "READ");
    if (!fin_fits || fin_fits->IsZombie()){
        std::cerr << "[build] FATAL: cannot open " << in_root << std::endl; return;
    }

    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");

    // --- PbPb: clone the colleague's TF1 fit per (centrality interval, q*eta slice) ---
    // The downstream lookup evaluates these fits at the exact muon pT.
    int nWritten = 0, nMissing = 0;
    for (size_t ic=0; ic<kCtrIntervals.size(); ++ic){
        const int lo = kCtrIntervals[ic].first, hi = kCtrIntervals[ic].second;
        const int cidx = kCtrColleagueIdx[ic];
        for (size_t is=0; is<kQEtaSlices.size(); ++is){
            const int eidx = (int)is;  // slice index == colleague eta bin index
            std::string fitName = Form("tf1_eff_fit_cent%d_eta%d", cidx, eidx);
            TF1* fit = (TF1*)fin_fits->Get(fitName.c_str());
            std::string nm = "tf1_reco_eff_medium_pbpb" + ctrKey(lo,hi)
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
    // --- pp: barrel + endcap (UNCHANGED; digitized points -> TGraph) ---
    fout->cd();
    makePpGraph("gr_reco_eff_medium_pp_barrel", kEffPpBarrel)->Write();
    makePpGraph("gr_reco_eff_medium_pp_endcap", kEffPpEndcap)->Write();

    fout->Close();
    fin_fits->Close();
    std::cout << "[build] wrote " << out_root << " (" << nWritten
              << " PbPb fits, " << nMissing << " missing, + 2 pp graphs)" << std::endl;

    // --- Reproduction plots: one 3x3 canvas per centrality (F.2 layout) ---
    TFile* fin = TFile::Open(out_root.c_str(), "READ");
    for (size_t ic=0; ic<kCtrIntervals.size(); ++ic){
        const int lo = kCtrIntervals[ic].first, hi = kCtrIntervals[ic].second;
        TCanvas c("c","c",1200,1000);
        c.Divide(3,3);
        for (size_t is=0; is<kQEtaSlices.size(); ++is){
            c.cd(is+1);
            gPad->SetGridx(); gPad->SetGridy();
            std::string nm = "tf1_reco_eff_medium_pbpb" + ctrKey(lo,hi)
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
            t.DrawLatex(0.18,0.85, Form("%.1f<q*#eta<%.1f, Medium-#mu",
                        kQEtaSlices[is].first, kQEtaSlices[is].second));
            t.DrawLatex(0.18,0.79, Form("%d-%d%%  (Run 2 fit, PLACEHOLDER)", lo, hi));
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
