// Abnormal-run investigation (docs/tracking/data_per_run_signal_yield_qa.md R2d): muon-pair distributions of ONE run
// vs the REST of its dataset, unit-area shapes + ratio, from the NOMINAL muon-pair NTuple (no re-derivation; Pb+Pb
// event selection already applied at NTuple stage; bad runs/LBs filtered again here as a no-op guard).
// Base selection = the crossx signal selection used in plot_signal_yield_per_lumi_per_run (OS = sign2, Tight,
// ParamsSet signal cuts, Pb+Pb centrality [0,80)); each signal-cut variable is shown N-1 (its own cut dropped).
// Binnings are the canonical ones (ParamsSet / CommonEffcyConfig), never retyped.
// Also writes <outcsv>: per-run counts for ALL runs (signal OS, signal SS, OS Tight before signal cuts, mean centrality).
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "RDFBasedHistFilling/CommonEffcyConfig.h"
#include "NTupleProcessingCode/PbPbEventSelConfig.h"

namespace pdr {
inline std::vector<double> EdgesFromRanges(const QEtaBinning& r) {
    std::vector<double> e; for (auto& p : r) { if (e.empty() || std::abs(e.back() - p.first) > 1e-6) e.push_back(p.first); e.push_back(p.second); }
    return e;
}
struct Var { std::string name, expr, cut_drop, xtitle; std::vector<double> edges; bool logx, pbpb_only, both_muons; };
}

void pair_dist_run_vs_rest(const std::string& tag, const std::string& label, const std::string& file, bool is_pbpb,
                           int year, std::vector<int> runs, const std::string& plotbase, const std::string& outcsv) {
    using namespace pdr;
    ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0);
    ParamsSet ps; CommonEffcyConfig ce;

    // signal cuts, individually (N-1)
    std::map<std::string, std::string> cut = {
        {"minv",    ParamsSet::SignalMinvCutExpr("minv")},
        {"pair_pt", ParamsSet::SignalPairPtCutExpr("pair_pt")},
        {"gap",     ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && " + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta")},
        {"pair_eta",ParamsSet::PairFiducialEtaCutExpr("pair_eta")},
    };
    auto sig_except = [&](const std::string& drop) {
        std::string s; for (auto& kv : cut) if (kv.first != drop) s += (s.empty() ? "" : " && ") + kv.second;
        return s;
    };
    const std::string ctr = "avg_centrality >= " + std::to_string(ParamsSet::ctrbins.front())
                          + " && avg_centrality < " + std::to_string(ParamsSet::ctrbins.back());

    std::vector<double> mwide; for (int i = 0; i <= ps.minv_nbins[0]; ++i) mwide.push_back(ps.minv_max[0] * i / ps.minv_nbins[0]);
    std::vector<double> ctrb(ParamsSet::ctrbins.begin(), ParamsSet::ctrbins.end());
    std::vector<double> phib; for (int i = 0; i <= 16; ++i) phib.push_back(-M_PI + 2 * M_PI * i / 16);  // QA axis only
    std::vector<double> drb(ce.dr_bins_edges_for_reco_effcy.begin(), ce.dr_bins_edges_for_reco_effcy.end());
    const std::vector<Var> vars = {
        {"pair_pt",   "pair_pt",  "pair_pt",  "p_{T}^{#mu#mu} [GeV]",         ps.pair_pt_coarse_bins, true,  false, false},
        {"pair_eta",  "pair_eta", "pair_eta", "#eta^{#mu#mu}",                EdgesFromRanges(ce.pair_eta_proj_ranges_coarse_incl_gap), false, false, false},
        {"minv",      "minv",     "minv",     "m_{#mu#mu} [GeV]",             ParamsSet::SignalMinvRecoEffcyBins(), false, false, false},
        {"minv_wide", "minv",     "minv",     "m_{#mu#mu} [GeV]",             mwide,                  false, false, false},
        {"mu_pt",     "mu_pt",    "",         "single-muon p_{T} [GeV]",      ps.single_mu_pt_coarse_bins, true, false, true},
        {"mu_qeta",   "mu_qeta",  "gap",      "single-muon q#upoint#eta",     EdgesFromRanges(ce.q_eta_proj_ranges_coarse_incl_gap), false, false, true},
        {"mu_phi",    "mu_phi",   "",         "single-muon #phi",             phib,                   false, false, true},
        {"dr",        "dr",       "",         "#DeltaR_{#mu#mu}",             drb,                    false, false, false},
        {"centrality","avg_centrality","",    "centrality [%]",               ctrb,                   false, true,  false},
    };

    const auto bad_runs = PbPbBadRuns(year); const auto bad_lbs = PbPbBadLBRanges(year);
    auto good = [bad_runs, bad_lbs, is_pbpb](UInt_t r, UInt_t lb) { return !is_pbpb || !PbPbIsBadRunLB(bad_runs, bad_lbs, (int)r, lb); };

    auto base = [&](const char* tree) {
        ROOT::RDF::RNode n = ROOT::RDataFrame(tree, file).Filter(good, {"run_number", "lb"}).Filter("pair_pass_tight");
        if (is_pbpb) n = n.Filter(ctr);
        return n.Define("mu_pt",   "ROOT::RVecF{m1.pt, m2.pt}")
                .Define("mu_qeta", "ROOT::RVecF{float(m1.charge * m1.eta), float(m2.charge * m2.eta)}")
                .Define("mu_phi",  "ROOT::RVecF{m1.phi, m2.phi}");
    };
    ROOT::RDF::RNode os = base("muon_pair_tree_sign2"), ss = base("muon_pair_tree_sign1");

    // ---- per-run table, ALL runs
    auto os_sig = os.Filter(sig_except(""));
    auto ss_sig = ss.Filter(sig_except(""));
    auto r_os_all = os.Take<UInt_t>("run_number");
    auto r_os_sig = os_sig.Take<UInt_t>("run_number");
    auto r_ss_sig = ss_sig.Take<UInt_t>("run_number");
    auto c_os_sig = is_pbpb ? os_sig.Take<Int_t>("avg_centrality") : os_sig.Define("zero", "0").Take<Int_t>("zero");

    // ---- histograms: per target run and rest
    struct HH { ROOT::RDF::RResultPtr<TH1D> run, rest; };
    std::map<int, std::map<std::string, HH>> hh;
    for (int r : runs) for (auto& v : vars) {
        if (v.pbpb_only && !is_pbpb) continue;
        auto src = os.Filter(sig_except(v.cut_drop));
        const std::string isr = "run_number == " + std::to_string(r);
        ROOT::RDF::TH1DModel m(Form("h_%s_%d", v.name.c_str(), r), "", (int)v.edges.size() - 1, v.edges.data());
        ROOT::RDF::TH1DModel m2(Form("h_%s_%d_rest", v.name.c_str(), r), "", (int)v.edges.size() - 1, v.edges.data());
        hh[r][v.name] = {src.Filter(isr).Histo1D(m, v.expr), src.Filter("!(" + isr + ")").Histo1D(m2, v.expr)};
    }

    // per-run CSV
    std::map<unsigned, long> nos_all, nos, nss; std::map<unsigned, double> csum;
    for (auto r : *r_os_all) nos_all[r]++;
    for (size_t i = 0; i < r_os_sig->size(); ++i) { nos[(*r_os_sig)[i]]++; csum[(*r_os_sig)[i]] += (*c_os_sig)[i]; }
    for (auto r : *r_ss_sig) nss[r]++;
    std::ofstream o(outcsv);
    o << "run,n_os_tight_all,n_os_signal,n_ss_signal,mean_centrality_signal\n";
    for (auto& kv : nos_all) o << kv.first << "," << kv.second << "," << nos[kv.first] << "," << nss[kv.first] << ","
                               << (nos[kv.first] ? csum[kv.first] / nos[kv.first] : 0.) << "\n";

    // ---- plots
    for (int r : runs) {
        const std::string dir = plotbase + "/" + tag + "_run" + std::to_string(r) + "/pairs";
        gSystem->mkdir(dir.c_str(), true);
        for (auto& v : vars) {
            if (v.pbpb_only && !is_pbpb) continue;
            TH1D* a = (TH1D*)hh[r][v.name].run->Clone(); TH1D* b = (TH1D*)hh[r][v.name].rest->Clone();
            const double na = a->Integral(), nb = b->Integral();
            if (na <= 0 || nb <= 0) { std::cout << "empty " << v.name << " run " << r << "\n"; continue; }
            double chi2 = 0; int ndf = 0, igood = 0; const double p = a->Chi2TestX(b, chi2, ndf, igood, "UU");
            std::cout << Form("SHAPE %s run %d %-10s chi2/ndf = %5.1f/%2d  p = %.3g  (N_run = %.0f)\n",
                              tag.c_str(), r, v.name.c_str(), chi2, ndf, p, na);
            a->Scale(1. / na, "width"); b->Scale(1. / nb, "width");

            TCanvas c(Form("c_%s_%d", v.name.c_str(), r), "", 700, 700);
            TPad up("up", "", 0, 0.32, 1, 1); up.SetBottomMargin(0.02); up.SetLeftMargin(0.14); up.Draw();
            TPad dn("dn", "", 0, 0, 1, 0.32); dn.SetTopMargin(0.02); dn.SetBottomMargin(0.32); dn.SetLeftMargin(0.14); dn.Draw();
            up.cd(); if (v.logx) up.SetLogx();
            b->SetLineColor(kBlack); b->SetFillColor(kGray); b->SetFillStyle(1001);
            a->SetLineColor(kRed + 1); a->SetMarkerColor(kRed + 1); a->SetMarkerStyle(20);
            const double ymax = std::max(b->GetBinContent(b->GetMaximumBin()), a->GetBinContent(a->GetMaximumBin()) + a->GetBinError(a->GetMaximumBin()));
            b->SetMaximum(1.45 * ymax); b->SetMinimum(0);
            b->GetYaxis()->SetTitle("normalised to unit area / bin width"); b->GetYaxis()->SetTitleSize(0.05);
            b->GetXaxis()->SetLabelSize(0);
            b->Draw("HIST"); a->Draw("E1 SAME");
            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.042);
            tl.DrawLatex(0.17, 0.85, (label + ", OS Tight, signal selection" + (v.cut_drop.empty() ? "" : " (N-1)")).c_str());
            tl.DrawLatex(0.17, 0.79, Form("#chi^{2}/ndf = %.1f/%d, p = %.2g", chi2, ndf, p));
            TLegend lg(0.60, 0.60, 0.88, 0.75); lg.SetBorderSize(0); lg.SetFillStyle(0);
            lg.AddEntry(a, Form("run %d (%.0f pairs)", r, na), "lep");
            lg.AddEntry(b, Form("other runs (%.0f pairs)", nb), "f"); lg.Draw();

            dn.cd(); if (v.logx) dn.SetLogx();
            TH1D* rt = (TH1D*)a->Clone(); rt->Divide(b);
            rt->SetMinimum(0); rt->SetMaximum(2.5);
            rt->GetYaxis()->SetTitle("run / others"); rt->GetYaxis()->SetNdivisions(505);
            rt->GetYaxis()->SetTitleSize(0.10); rt->GetYaxis()->SetTitleOffset(0.6); rt->GetYaxis()->SetLabelSize(0.09);
            rt->GetXaxis()->SetTitle(v.xtitle.c_str()); rt->GetXaxis()->SetTitleSize(0.11); rt->GetXaxis()->SetLabelSize(0.09);
            rt->Draw("E1");
            TLine one(v.edges.front(), 1, v.edges.back(), 1); one.SetLineStyle(2); one.Draw();
            c.SaveAs((dir + "/" + v.name + "_run" + std::to_string(r) + ".png").c_str());
        }
    }
}
