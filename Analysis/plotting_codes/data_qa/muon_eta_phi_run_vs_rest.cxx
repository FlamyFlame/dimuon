// Abnormal-run investigation (docs/tracking/data_per_run_signal_yield_qa.md): single-muon eta vs phi of ONE run vs
// the rest of its dataset, from the NOMINAL muon-pair NTuple (no single-muon tree exists; muons are taken from the
// pairs and each muon is counted ONCE per event, key = run, LB, BCID, entry, muon index).
// Selections:  "all_tight"  = muons of Tight pairs, OS + SS, no signal cuts (statistics for detector holes)
//              "signal"     = muons of OS Tight signal-region pairs (as in plot_signal_yield_per_lumi_per_run)
// Pb+Pb: centrality [0,80) and bad runs/LBs filtered (no-op on a regenerated NTuple).
// Panels: run, other runs (unit area), ratio run/others.  The eta-phi grid is a QA display axis only.
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "NTupleProcessingCode/PbPbEventSelConfig.h"

void muon_eta_phi_run_vs_rest(const std::string& tag, const std::string& label, const std::string& file, bool is_pbpb,
                              int year, int run, const std::string& plotbase) {
    ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0);
    const std::string sig =
        ParamsSet::SignalMinvCutExpr("minv") + " && " + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
        + ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && " + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta")
        + " && " + ParamsSet::PairFiducialEtaCutExpr("pair_eta");
    const std::string ctr = "avg_centrality >= " + std::to_string(ParamsSet::ctrbins.front())
                          + " && avg_centrality < " + std::to_string(ParamsSet::ctrbins.back());
    const auto bad_runs = PbPbBadRuns(year); const auto bad_lbs = PbPbBadLBRanges(year);
    auto good = [bad_runs, bad_lbs, is_pbpb](UInt_t r, UInt_t lb) { return !is_pbpb || !PbPbIsBadRunLB(bad_runs, bad_lbs, (int)r, lb); };

    const int neta = 20, nphi = 16;   // 0.25 x pi/8
    struct Sel { std::string name, title; std::vector<std::string> trees; bool signal; };
    const std::vector<Sel> sels = {{"all_tight", "muons of Tight pairs (OS+SS)", {"muon_pair_tree_sign1", "muon_pair_tree_sign2"}, false},
                                   {"signal",    "muons of OS Tight signal-region pairs", {"muon_pair_tree_sign2"}, true}};
    const std::string dir = plotbase + "/" + tag + "_run" + std::to_string(run) + "/muons";
    gSystem->mkdir(dir.c_str(), true);

    for (auto& s : sels) {
        TH2D hr("hr", ";#eta^{#mu};#phi^{#mu}", neta, -2.5, 2.5, nphi, -M_PI, M_PI);
        TH2D ho("ho", ";#eta^{#mu};#phi^{#mu}", neta, -2.5, 2.5, nphi, -M_PI, M_PI);
        std::set<std::tuple<unsigned, unsigned, unsigned, int, int>> seen;
        for (auto& tr : s.trees) {
            ROOT::RDF::RNode n = ROOT::RDataFrame(tr, file).Filter(good, {"run_number", "lb"}).Filter("pair_pass_tight");
            if (is_pbpb) n = n.Filter(ctr);
            if (s.signal) n = n.Filter(sig);
            auto R = n.Take<UInt_t>("run_number"); auto L = n.Take<UInt_t>("lb"); auto B = n.Take<UInt_t>("bcid");
            auto E = n.Take<Int_t>("m1.ev_num");
            auto i1 = n.Take<Int_t>("m1.ind"), i2 = n.Take<Int_t>("m2.ind");
            auto e1 = n.Take<Float_t>("m1.eta"), e2 = n.Take<Float_t>("m2.eta");
            auto p1 = n.Take<Float_t>("m1.phi"), p2 = n.Take<Float_t>("m2.phi");
            for (size_t k = 0; k < R->size(); ++k) {
                const int idx[2] = {(*i1)[k], (*i2)[k]}; const float eta[2] = {(*e1)[k], (*e2)[k]}, phi[2] = {(*p1)[k], (*p2)[k]};
                for (int m = 0; m < 2; ++m) {
                    if (!seen.insert({(*R)[k], (*L)[k], (*B)[k], (*E)[k], idx[m]}).second) continue;
                    ((int)(*R)[k] == run ? hr : ho).Fill(eta[m], phi[m]);
                }
            }
        }
        const double nr = hr.Integral(), no = ho.Integral();
        double chi2 = 0; int ndf = 0, igood = 0; const double pval = hr.Chi2TestX(&ho, chi2, ndf, igood, "UU");
        std::cout << Form("%s run %d [%s]: %.0f muons in run, %.0f in other runs; chi2/ndf = %.1f/%d, p = %.3g\n",
                          tag.c_str(), run, s.name.c_str(), nr, no, chi2, ndf, pval);
        TH2D* ratio = (TH2D*)hr.Clone("ratio");
        ratio->Divide(&ho); ratio->Scale(no / nr);
        for (int i = 1; i <= neta; ++i) for (int j = 1; j <= nphi; ++j) if (ho.GetBinContent(i, j) <= 0) ratio->SetBinContent(i, j, 0);
        TH2D* hrn = (TH2D*)hr.Clone("hrn"); hrn->Scale(1. / nr);
        TH2D* hon = (TH2D*)ho.Clone("hon"); hon->Scale(1. / no);
        const double zmax = std::max(hrn->GetMaximum(), hon->GetMaximum());
        hrn->SetMaximum(zmax); hon->SetMaximum(zmax); hrn->SetMinimum(0); hon->SetMinimum(0);
        hrn->GetZaxis()->SetTitle("fraction of muons"); hon->GetZaxis()->SetTitle("fraction of muons");
        ratio->GetZaxis()->SetTitle("run / others (unit area)");
        ratio->SetMinimum(0); ratio->SetMaximum(2);

        ratio->SetMarkerSize(0.6); gStyle->SetPaintTextFormat(".1f");
        TCanvas c(Form("c_%s", s.name.c_str()), "", 1800, 560); c.Divide(3, 1, 0.002, 0.002);
        TH2D* hs[3] = {hrn, hon, ratio};
        const std::string ttl[3] = {Form("run %d (%.0f muons)", run, nr), Form("other runs (%.0f muons)", no), Form("run / other runs, #chi^{2}/ndf = %.0f/%d, p = %.2g", chi2, ndf, pval)};
        for (int p = 0; p < 3; ++p) {
            c.cd(p + 1); gPad->SetRightMargin(0.18); gPad->SetLeftMargin(0.12); gPad->SetTopMargin(0.13); gPad->SetBottomMargin(0.12);
            hs[p]->GetZaxis()->SetTitleOffset(1.4); hs[p]->SetContour(99); hs[p]->Draw(p == 2 ? "COLZ TEXT" : "COLZ");
            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.045);
            tl.DrawLatex(0.12, 0.94, (label + ", " + s.title).c_str());
            tl.DrawLatex(0.12, 0.89, ttl[p].c_str());
        }
        c.SaveAs((dir + "/muon_eta_phi_" + s.name + "_run" + std::to_string(run) + ".png").c_str());
    }
}
