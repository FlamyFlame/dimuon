// plot_sig_accept_cutflow_above_60GeV.cxx
//
// Diagnostic: cumulative and per-cut signal acceptance for op-sign, from_same_b pairs.
// Overlays two pT selections on the same canvas:
//   - all pair pT  (black markers)
//   - pair pT > 60 GeV  (red markers)
//
// Two plots per MC:
//   prefix_sig_accept_above_60GeV_accum.png
//       #bins = #cuts + 1; bin 0 = measured range (A_sig denominator); bin i = crossx after
//       the first i mass cuts / crossx in the measured range (last bin = integrated A_sig)
//   prefix_sig_accept_above_60GeV.png
//       #bins = #cuts; bin i = crossx_i / crossx_{i-1}
//
// Run (ACLiC recommended for struct member access in JIT filters):
//   root -l -b -q 'plotting_codes/single_b_analysis/plot_sig_accept_cutflow_above_60GeV.cxx+'

#include "../../MuonObjectsParamsAndHelpers/MuonPairPythia.h"
#include "../../MuonObjectsParamsAndHelpers/MuonPairPowheg.h"
#include "../../Utilities/MC_helpers.h"
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../Utilities/SignalAcceptanceDefinition.h"
#include "TString.h"   // Form(), used to build the cut labels from the ParamsSet values

#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>

// Cutflow of the truth signal acceptance A_sig (2026-10-05, user: the measured-range cuts are
// NOT acceptance). The starting point ("measured range") is the A_sig DENOMINATOR,
// SignalAcceptanceDefinition::DenomExpr() = from_same_b && truth pair pT > 9 && |truth eta^pair| < 2.2
// -- these cuts DEFINE the measurement and are never shown as an acceptance loss. The cuts below are
// exactly what A_sig corrects: the truth mass window. So the last cumulative bin IS the integrated
// A_sig. NO detector-gap cut (it is in eps_reco). Thresholds READ FROM ParamsSet, labels built from
// the values (docs/tracking/signal_acceptance_mass_window.md).
static const std::vector<std::pair<std::string, std::string>>& kCuts() {
    static const std::vector<std::pair<std::string, std::string>> cuts = {
        {Form("truth_minv > %g", ParamsSet::signal_minv_min), Form("m_{#mu#mu}>%g", ParamsSet::signal_minv_min)},
        {Form("truth_minv < %g", ParamsSet::signal_minv_max), Form("m_{#mu#mu}<%g", ParamsSet::signal_minv_max)},
    };
    return cuts;
}

// Compute both cutflows (all pT and pT>60) in a single RDF event-loop pass.
// Returns {sums_all, sums_60}, each of length ncuts+1.
static std::pair<std::vector<double>, std::vector<double>>
RunCutflowPair(ROOT::RDF::RNode df, const std::string& weight_col)
{
    auto df_all = df.Filter(SignalAcceptanceDefinition::DenomExpr());
    auto df_60  = df.Filter(SignalAcceptanceDefinition::DenomExpr() + " && truth_pair_pt > 60.0");

    std::vector<ROOT::RDF::RResultPtr<double>> ptrs_all, ptrs_60;
    auto dfcur_all = df_all;
    auto dfcur_60  = df_60;

    // Collect all result pointers before any dereference so RDF evaluates in one pass.
    ptrs_all.push_back(dfcur_all.Sum<double>(weight_col));
    ptrs_60.push_back(dfcur_60.Sum<double>(weight_col));
    for (auto& [expr, lbl] : kCuts()) {
        dfcur_all = dfcur_all.Filter(expr);
        dfcur_60  = dfcur_60.Filter(expr);
        ptrs_all.push_back(dfcur_all.Sum<double>(weight_col));
        ptrs_60.push_back(dfcur_60.Sum<double>(weight_col));
    }

    std::vector<double> sums_all, sums_60;
    for (auto& p : ptrs_all) sums_all.push_back(*p);
    for (auto& p : ptrs_60)  sums_60.push_back(*p);
    // Integrated A_sig = after the mass window / measured range (bin 0).
    const size_t k = 0;
    if (sums_all.at(k) > 0 && sums_60.at(k) > 0)
        std::cout << "[INFO] integrated A_sig (mass window | measured range): all pT = "
                  << sums_all.back() / sums_all.at(k) << ", pT>60 = " << sums_60.back() / sums_60.at(k) << "\n";
    return {sums_all, sums_60};
}

// Style helpers
static void StyleHist(TH1D& h, Color_t col, Style_t marker) {
    h.SetLineWidth(2);
    h.SetLineColor(col);
    h.SetMarkerStyle(marker);
    h.SetMarkerColor(col);
    h.SetMarkerSize(1.3);
}

static void MakePlots(const std::vector<double>& sums_all,
                      const std::vector<double>& sums_60,
                      const std::string& out_dir,
                      const std::string& prefix,
                      const std::string& mc_label)
{
    gStyle->SetOptStat(0);
    const int ncuts = static_cast<int>(kCuts().size());

    const double norm_all = sums_all.at(0);
    const double norm_60  = sums_60.at(0);
    if (norm_all <= 0 || norm_60 <= 0) {
        std::cerr << "[WARN] zero norm, skipping " << prefix << "\n";
        return;
    }

    // x-axis labels
    std::vector<std::string> lbls;
    // first bin = the measured range itself (A_sig denominator), labelled by its physical cuts
    lbls.push_back(Form("p_{T}^{pair}>%g, |#eta^{pair}|<%g", ParamsSet::signal_pair_pt_min, ParamsSet::pair_eta_fiducial_max));
    for (auto& [expr, lbl] : kCuts()) lbls.push_back(lbl);

    // ========================== Accum plot ==========================
    TH1D h_acc_all(("h_acc_all_" + prefix).c_str(), "", ncuts + 1, -0.5, ncuts + 0.5);
    TH1D h_acc_60 (("h_acc_60_"  + prefix).c_str(), "", ncuts + 1, -0.5, ncuts + 0.5);
    for (int i = 0; i <= ncuts; ++i) {
        h_acc_all.SetBinContent(i + 1, sums_all.at(i) / norm_all);
        h_acc_60 .SetBinContent(i + 1, sums_60 .at(i) / norm_60);
        h_acc_all.GetXaxis()->SetBinLabel(i + 1, lbls.at(i).c_str());
        h_acc_60 .GetXaxis()->SetBinLabel(i + 1, lbls.at(i).c_str());
    }
    {
        TCanvas c(("c_acc_" + prefix).c_str(), "", 1000, 580);
        c.SetBottomMargin(0.28); c.SetLeftMargin(0.14); c.SetRightMargin(0.05);

        StyleHist(h_acc_all, kBlack,   20);
        StyleHist(h_acc_60,  kRed + 1, 24);

        h_acc_all.GetYaxis()->SetTitle("cumul. fraction (normalized to measured range)");
        h_acc_all.GetYaxis()->SetRangeUser(0., 1.4);
        h_acc_all.GetXaxis()->SetLabelSize(0.048);
        h_acc_all.GetXaxis()->LabelsOption("v");
        h_acc_all.Draw("HIST P");
        h_acc_60 .Draw("HIST P SAME");

        TLegend leg(0.55, 0.74, 0.93, 0.90);
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.040);
        leg.AddEntry(&h_acc_all, "all p_{T}^{pair}",       "lep");
        leg.AddEntry(&h_acc_60,  "p_{T}^{pair} > 60 GeV", "lep");
        leg.Draw();

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.040);
        lat.DrawLatex(0.16, 0.93, mc_label.c_str());
        lat.DrawLatex(0.16, 0.88, "Truth single-b");

        const std::string out = out_dir + "/" + prefix + "_sig_accept_above_60GeV_accum.png";
        c.SaveAs(out.c_str());
        std::cout << "[INFO] Saved: " << out << "\n";
    }

    // ========================== Per-cut plot ==========================
    TH1D h_per_all(("h_per_all_" + prefix).c_str(), "", ncuts, 0.5, ncuts + 0.5);
    TH1D h_per_60 (("h_per_60_"  + prefix).c_str(), "", ncuts, 0.5, ncuts + 0.5);
    for (int i = 1; i <= ncuts; ++i) {
        const double r_all = (sums_all.at(i-1) > 0) ? sums_all.at(i) / sums_all.at(i-1) : 0.;
        const double r_60  = (sums_60 .at(i-1) > 0) ? sums_60 .at(i) / sums_60 .at(i-1) : 0.;
        h_per_all.SetBinContent(i, r_all);
        h_per_60 .SetBinContent(i, r_60);
        h_per_all.GetXaxis()->SetBinLabel(i, lbls.at(i).c_str());
        h_per_60 .GetXaxis()->SetBinLabel(i, lbls.at(i).c_str());
    }
    {
        TCanvas c(("c_per_" + prefix).c_str(), "", 1000, 580);
        c.SetBottomMargin(0.28); c.SetLeftMargin(0.14); c.SetRightMargin(0.05);

        StyleHist(h_per_all, kBlack,   20);
        StyleHist(h_per_60,  kRed + 1, 24);

        h_per_all.GetYaxis()->SetTitle("per-cut fraction (given all prior cuts passed)");
        h_per_all.GetYaxis()->SetRangeUser(0., 1.4);
        h_per_all.GetXaxis()->SetLabelSize(0.048);
        h_per_all.GetXaxis()->LabelsOption("v");
        h_per_all.Draw("HIST P");
        h_per_60 .Draw("HIST P SAME");

        TLegend leg(0.55, 0.74, 0.93, 0.90);
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.040);
        leg.AddEntry(&h_per_all, "all p_{T}^{pair}",       "lep");
        leg.AddEntry(&h_per_60,  "p_{T}^{pair} > 60 GeV", "lep");
        leg.Draw();

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.040);
        lat.DrawLatex(0.16, 0.93, mc_label.c_str());
        lat.DrawLatex(0.16, 0.88, "Truth single-b");

        const std::string out = out_dir + "/" + prefix + "_sig_accept_above_60GeV.png";
        c.SaveAs(out.c_str());
        std::cout << "[INFO] Saved: " << out << "\n";
    }
}

void plot_sig_accept_cutflow_above_60GeV()
{
    const std::string pythia_file =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/"
        "muon_pairs_pythia_5p36TeV_no_data_resonance_cuts.root";

    // cc truth file does not yet exist as a merged single file; use bb only.
    const std::vector<std::string> powheg_files = {
        // NO-resonance-cut NTuple: A_sig needs ALL masses in the denominator; the nominal POWHEG
        // truth NTuple is resonance-vetoed (docs/tracking/signal_acceptance_mass_window.md R2).
        "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_evgen_truth_full_sample/muon_pairs_powheg_bb_truth_no_data_resonance_cuts.root",
    };

    const std::string pythia_out =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pythia";
    const std::string powheg_out =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/powheg";
    gSystem->mkdir(pythia_out.c_str(), true);
    gSystem->mkdir(powheg_out.c_str(), true);

    // ---- Pythia ----
    std::cout << "\n[INFO] Pythia cutflow...\n";
    {
        auto df = ROOT::RDataFrame("muon_pair_tree_sign2",
                                   std::vector<std::string>{pythia_file});
        auto [sums_all, sums_60] = RunCutflowPair(df, "weight");
        std::cout << "[INFO] Pythia all-pT sums:"; for (double s : sums_all) std::cout << " " << s; std::cout << "\n";
        std::cout << "[INFO] Pythia  >60 sums:";  for (double s : sums_60)  std::cout << " " << s; std::cout << "\n";
        MakePlots(sums_all, sums_60, pythia_out, "pythia",
                  "Pythia 8.3, pp #sqrt{s_{NN}} = 5.36 TeV");
    }

    // ---- Powheg ----
    std::cout << "\n[INFO] Powheg cutflow...\n";
    {
        const double norm = SumMetaNentriesBeforeFilter(powheg_files);
        std::cout << "[INFO] Powheg nentries_before_cuts_sum = " << norm << "\n";
        auto df = ROOT::RDataFrame("muon_pair_tree_sign2", powheg_files)
                      .Define("weight_norm", Form("weight / %.15g", norm));
        auto [sums_all, sums_60] = RunCutflowPair(df, "weight_norm");
        std::cout << "[INFO] Powheg all-pT sums:"; for (double s : sums_all) std::cout << " " << s; std::cout << "\n";
        std::cout << "[INFO] Powheg  >60 sums:";  for (double s : sums_60)  std::cout << " " << s; std::cout << "\n";
        MakePlots(sums_all, sums_60, powheg_out, "powheg",
                  "POWHEG+Pythia8, pp #sqrt{s_{NN}} = 5.36 TeV");
    }
}
