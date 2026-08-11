// =================================================================================================
// write_mc_pair_statistics_tables.cxx
//
// HOW MUCH SAMPLE does each (pair pT, pair eta) cell of the MC DeltaR correction actually have?
// This macro answers that in CSV form, for both pair charge combinations, at both working points.
//
// IT DOES NOT TOUCH THE MUON-PAIR TREES. It reads the six TH2Ds that FillMCTrigEffHists.cxx books
// ON THE STEP-3 DENOMINATOR NODE (see its "STATISTICS BOOKKEEPING" block, filled right after the
// `step3 selection` Filter). That is the whole point: those histograms are one entry per selected
// pair under EXACTLY the selection the DeltaR correction is measured with, so this table cannot
// drift from it. A macro that re-opened the trees and re-derived the selection would be free to.
//
//   h_mc_paircount_vs_pt_eta_{ss,os}  raw pair counts (UNWEIGHTED) -> the statistical sample size
//   h_mc_pairsumw_vs_pt_eta_{ss,os}   sum of the per-pair MC weight = the cell cross section, in nb
//   h_mc_pairsumw2_vs_pt_eta_{ss,os}  sum of weight^2 -> stat. error on the cross section = sqrt()
//
// Sign convention (Analysis README / repo-wide): muon_pair_tree_sign1 = SAME sign  -> "ss"
//                                                muon_pair_tree_sign2 = OPPOSITE sign -> "os"
//
// SIX CSVs per working point, in three sets of two (counts | cross section):
//   pair_pt_counts.csv                   pair_pt_crossx.csv                 (eta-INTEGRATED, 2 cols)
//   same_sign_pt_vs_eta_counts.csv       same_sign_pt_vs_eta_crossx.csv     (pT = cols, eta = rows)
//   opposite_sign_pt_vs_eta_counts.csv   opposite_sign_pt_vs_eta_crossx.csv (idem)
//
// Every table is written by ONE lambda that is called twice -- once with std::cout, once with the
// output file stream -- so the log and the file cannot diverge (the plot_mc_trig_eff.cxx
// write_table_A convention).
//
// BINNINGS: never retyped. Row and column labels are formatted from the histogram axes themselves.
// PATHS: never hardcoded. Input dir / file label / plot root all come from GetDrCorrSample().
//
// Compile/run (ACLiC, from this directory):
//   root -l -b -q 'write_mc_pair_statistics_tables.cxx+("pp_full", true)'    // Tight
//   root -l -b -q 'write_mc_pair_statistics_tables.cxx+("pp_full", false)'   // Medium
// =================================================================================================

#include <TAxis.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

// Sample identity (input dir, file label, plot root) + the WP file suffix + the pair-pT binning
// token, shared with the fill and the plot stages so this table can only ever describe the file
// those stages produced.
#include "dr_correction_sample_cfg.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"
#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"

namespace {

std::string Fmt(const char* f, double a)           { char b[256]; snprintf(b, sizeof b, f, a);       return b; }
std::string Fmt(const char* f, double a, double b_){ char b[256]; snprintf(b, sizeof b, f, a, b_);   return b; }

// The canvas headline is a ROOT-latex string; a CSV comment is plain text. Only the two
// beam-energy tokens need translating -- deliberately targeted, so an unexpected latex construct
// stays visible rather than being silently mangled.
std::string PlainText(std::string s)
{
    const std::vector<std::pair<std::string, std::string>> subs = {
        {"#sqrt{s_{NN}}", "sqrt(s_NN)"}, {"#sqrt{s}", "sqrt(s)"}};
    for (const auto& kv : subs)
        for (size_t p = s.find(kv.first); p != std::string::npos; p = s.find(kv.first, p))
            s.replace(p, kv.first.size(), kv.second);
    return s;
}

// Missing key = HARD FAILURE. A silently-zero table would look exactly like a sample with no
// statistics, which is the one thing this table exists to measure.
TH2D* Get2D(TFile* f, const std::string& name)
{
    TH2D* h = dynamic_cast<TH2D*>(f->Get(name.c_str()));
    if (!h)
        throw std::runtime_error(
            "write_mc_pair_statistics_tables: histogram '" + name + "' is MISSING from " +
            f->GetName() + ".\n  The Step-3 statistics histograms are booked by "
            "FillMCTrigEffHists.cxx (STATISTICS BOOKKEEPING block); a file filled before that "
            "block was added does not contain them. Re-run the filler for this sample/WP.");
    h->SetDirectory(nullptr);
    return h;
}

// All six must live on the SAME (pair pT, pair eta) grid, or a cell of one table would not be the
// same cell of the next.
void CheckSameAxes(const TH2D* ref, const TH2D* h)
{
    const bool same = h->GetNbinsX() == ref->GetNbinsX() && h->GetNbinsY() == ref->GetNbinsY();
    if (!same)
        throw std::runtime_error(std::string("write_mc_pair_statistics_tables: binning mismatch "
                                 "between '") + ref->GetName() + "' and '" + h->GetName() + "'");
    for (int i = 1; i <= ref->GetNbinsX() + 1; ++i)
        if (std::fabs(h->GetXaxis()->GetBinLowEdge(i) - ref->GetXaxis()->GetBinLowEdge(i)) > 1e-6)
            throw std::runtime_error(std::string("write_mc_pair_statistics_tables: pair-pT edge ")
                                     + std::to_string(i) + " differs between '" + ref->GetName()
                                     + "' and '" + h->GetName() + "'");
    for (int i = 1; i <= ref->GetNbinsY() + 1; ++i)
        if (std::fabs(h->GetYaxis()->GetBinLowEdge(i) - ref->GetYaxis()->GetBinLowEdge(i)) > 1e-6)
            throw std::runtime_error(std::string("write_mc_pair_statistics_tables: pair-eta edge ")
                                     + std::to_string(i) + " differs between '" + ref->GetName()
                                     + "' and '" + h->GetName() + "'");
}

// The repo convention: one writer, called with std::cout and with the file, so the two can never
// say different things.
void Emit(const std::string& path, const std::function<void(std::ostream&)>& write)
{
    std::cout << "\n----- " << path << " -----\n";
    write(std::cout);
    std::ofstream ofs(path);
    if (!ofs)
        throw std::runtime_error("write_mc_pair_statistics_tables: cannot open '" + path +
                                 "' for writing");
    write(ofs);
}

}  // namespace

void write_mc_pair_statistics_tables(const std::string& sample = "pp_full",
                                     bool use_tight_wp = true)
{
    gROOT->SetBatch(kTRUE);

    const DrCorrSample cfg    = GetDrCorrSample(sample, use_tight_wp);
    const std::string  wp_suf = DrCorrWpSuffix(use_tight_wp);        // "" | "_medium_wp"
    const std::string  wp_txt = use_tight_wp ? "Tight" : "Medium";

    // Input: the Step-3 fill output. Same construction as plot_mc_trig_eff.cxx, so the two stages
    // read the same file for the same (sample, WP, pair-pT binning).
    const std::string in_path = cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf
                              + MCTrigEffPairPt::FileSuffix() + "_step3.root";
    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie())
        throw std::runtime_error("write_mc_pair_statistics_tables: cannot open " + in_path);

    struct SignSet { std::string key, label, file_stem; TH2D *cnt, *sw, *sw2; };
    std::vector<SignSet> signs = {
        {"ss", "same sign",     "same_sign",     nullptr, nullptr, nullptr},
        {"os", "opposite sign", "opposite_sign", nullptr, nullptr, nullptr}
    };
    for (auto& s : signs) {
        s.cnt = Get2D(fin, "h_mc_paircount_vs_pt_eta_" + s.key);
        s.sw  = Get2D(fin, "h_mc_pairsumw_vs_pt_eta_"  + s.key);
        s.sw2 = Get2D(fin, "h_mc_pairsumw2_vs_pt_eta_" + s.key);
    }
    fin->Close();

    const TH2D* ref = signs[0].cnt;
    for (const auto& s : signs) { CheckSameAxes(ref, s.cnt); CheckSameAxes(ref, s.sw);
                                  CheckSameAxes(ref, s.sw2); }

    const int npt  = ref->GetNbinsX();
    const int neta = ref->GetNbinsY();
    const TAxis* ax = ref->GetXaxis();
    const TAxis* ay = ref->GetYaxis();

    // Consistency between the binning this file was FILLED with and the one this process would
    // select: they are both driven by MCTRIGEFF_PAIRPT_4BIN, so a mismatch means the environment
    // changed between fill and table and the file name is lying.
    if (npt != MCTrigEffPairPt::NBins())
        throw std::runtime_error("write_mc_pair_statistics_tables: " + in_path + " has "
            + std::to_string(npt) + " pair-pT bins but this process is configured for "
            + std::to_string(MCTrigEffPairPt::NBins()) + " (" + MCTrigEffPairPt::Describe() + ")");

    // Nothing to tabulate = hard failure, never a table of zeros.
    double grand = 0.;
    for (const auto& s : signs) grand += s.cnt->Integral();
    if (!(grand > 0.))
        throw std::runtime_error("write_mc_pair_statistics_tables: the statistics histograms in "
                                 + in_path + " are EMPTY inside the binned range -- refusing to "
                                 "write a table of zeros");

    // -------- output directory -------------------------------------------------------------
    // cfg.out_base is ".../<beam>_trigger_efficiency/mc_based<variant>/". These tables are a
    // sibling product set, not part of the mc_based plot tree, so they get their own TOP-LEVEL
    // directory under the same plot root -- and, per the VARIANT LAYOUT block in
    // dr_correction_sample_cfg.h, the working point lives in that top-level name.
    // The pair-pT token is taken from the ACTUAL number of bins in the file, so the directory
    // name can never claim a binning the numbers inside it do not have.
    std::string base = cfg.out_base;
    if (!base.empty() && base.back() == '/') base.pop_back();
    const size_t slash = base.find_last_of('/');
    if (slash == std::string::npos || base.compare(slash + 1, 8, "mc_based") != 0)
        throw std::runtime_error("write_mc_pair_statistics_tables: unexpected out_base '"
            + cfg.out_base + "' (expected its last component to start with 'mc_based')");
    const std::string out_dir = base.substr(0, slash + 1) + "mc_statistics_pt"
                              + std::to_string(npt) + "bins"
                              + (use_tight_wp ? "" : "_medium") + "/";
    gSystem->mkdir(out_dir.c_str(), kTRUE);

    // -------- labels, formatted from the axes (never retyped) ------------------------------
    auto pt_lab  = [&](int ix) { return Fmt("%.1f-%.1f", ax->GetBinLowEdge(ix), ax->GetBinUpEdge(ix)); };
    auto eta_lab = [&](int iy) { return Fmt("%.1f to %.1f", ay->GetBinLowEdge(iy), ay->GetBinUpEdge(iy)); };

    // -------- provenance header, identical on all six files --------------------------------
    // The Step-3 selection wording follows mc_trigger_efficiency.md §3.0/§3.3 and the
    // `step3 selection` Filter in RDFBasedHistFilling/FillMCTrigEffHists.cxx. The gap-cut windows
    // are printed from ParamsSet so the header cannot state windows the fill did not apply.
    auto header = [&](std::ostream& os, const std::string& quantity) {
        os << "# MC muon-pair statistics of the Step-3 DeltaR-correction sample "
              "(mc_trigger_efficiency.md 3.3).\n";
        os << "# sample=" << cfg.key << "  label=" << cfg.mc_label
           << "  identity=" << PlainText(cfg.sample_text) << "\n";
        os << "# muon working point = " << wp_txt << " (required of BOTH legs)\n";
        os << "# source = " << in_path << "\n";
        os << "#   histograms h_mc_{paircount,pairsumw,pairsumw2}_vs_pt_eta_{ss,os}, booked by\n";
        os << "#   RDFBasedHistFilling/FillMCTrigEffHists.cxx on the Step-3 DENOMINATOR node, i.e.\n";
        os << "#   one entry per selected pair under exactly the selection quoted below.\n";
        os << "# sign convention: muon_pair_tree_sign1 = SAME sign (ss), "
              "muon_pair_tree_sign2 = OPPOSITE sign (os).\n";
        os << "#\n";
        if (quantity == "counts") {
            os << "# QUANTITY: counts = raw number of reconstructed muon pairs, UNWEIGHTED and\n";
            os << "#   DIMENSIONLESS. This is the statistical sample size of the cell -- it is NOT\n";
            os << "#   a yield: the MC pT-hat slices are mixed with very different weights.\n";
        } else {
            os << "# QUANTITY: crossx = sum over the pairs of the cell of the per-pair MC weight\n";
            os << "#   w = sigma_slice * eps_filt * r_isospin / N_slice, i.e. the CROSS SECTION of\n";
            os << "#   the cell. UNITS: **nb** (nanobarn). The AMI cross sections this weight is\n";
            os << "#   built from are in nb, NOT pb (Analysis/docs/ami_weights.md); to compare with\n";
            os << "#   pp data, whose luminosity is in pb^-1, MULTIPLY THESE NUMBERS BY 1000.\n";
            os << "#   STATISTICAL ERROR = sqrt(sum of w^2) over the same cells, reported in its\n";
            os << "#   own column (see the column header line below). Same units, nb.\n";
        }
        os << "#\n";
        os << "# SELECTION (Step-3 selection; mc_trigger_efficiency.md 3.0/3.3, code:\n";
        os << "#   RDFBasedHistFilling/FillMCTrigEffHists.cxx, Filter \"step3 selection\").\n";
        os << "#   Input tree: the mc_trig muon-pair trees (no data resonance cuts applied).\n";
        os << "#   Required of BOTH muons of the pair:\n";
        os << "#     - the generic muon cuts of the data analysis at the " << wp_txt
           << " working point\n";
        os << "#       (combined quality + WP bit + IDCuts + MuonCuts, dp/p < 0.12 ONE-SIDED,\n";
        os << "#        |d0| < 2 mm, |z0 sin(theta)| < 2 mm) -- the GENERIC selection, not the\n";
        os << "#        signal selection;\n";
        os << "#     - pT > 4 GeV and |eta| < 2.4;\n";
        os << "#     - truth fiducial: truth pT > 4 GeV and |truth eta| < 2.4;\n";
        os << "#     - forward low-pT veto: (pT > 7 GeV) OR (q*eta > -2);\n";
        os << "#     - q*eta fiducial gap cut: " << ParamsSet::FiducialGapCutExpr("q*eta") << ";\n";
        if (cfg.key == "overlay")
            os << "#     - event centrality in [0,5) %.\n";
        os << "#   NO TRIGGER REQUIREMENT (these are the Step-3 denominator pairs).\n";
        os << "#   NO SIGNAL-SELECTION CUT: no pair-pT, no pair q*eta, no DeltaR, no m_mumu and no\n";
        os << "#   resonance veto. Those define the measurement, not the muon.\n";
        os << "#\n";
        os << "# BINNED RANGE -- pairs OUTSIDE it are NOT COUNTED anywhere in this table:\n";
        os << "#   pair pT must be in [" << Fmt("%.4g", ax->GetXmin()) << ", "
           << Fmt("%.4g", ax->GetXmax()) << "] GeV and pair eta in ["
           << Fmt("%.4g", ay->GetXmin()) << ", " << Fmt("%.4g", ay->GetXmax()) << "].\n";
        os << "#   This is why these totals are SMALLER than the totals of the DeltaR histograms\n";
        os << "#   filled from the same node, which have no pair-pT / pair-eta range restriction.\n";
        for (const auto& s : signs) {
            const double in_r  = s.cnt->Integral();
            const double all_r = s.cnt->Integral(0, npt + 1, 0, neta + 1);
            os << "#   " << s.label << ": " << Fmt("%.0f", in_r) << " of " << Fmt("%.0f", all_r)
               << " selected pairs are inside the binned range ("
               << Fmt("%.2f", all_r > 0. ? 100. * in_r / all_r : 0.) << " %).\n";
        }
        os << "#\n";
        os << "# Bin edges, exactly as used (the labels below are these edges rounded to 0.1):\n";
        os << "#   pair pT [GeV] :";
        for (int ix = 1; ix <= npt + 1; ++ix) os << " " << Fmt("%.4f", ax->GetBinLowEdge(ix));
        os << "\n#   pair eta      :";
        for (int iy = 1; iy <= neta + 1; ++iy) os << " " << Fmt("%.4f", ay->GetBinLowEdge(iy));
        os << "\n#\n";
    };

    // ============================ SET 1: pair pT, pair eta INTEGRATED =======================
    // One row per pair-pT bin, one data column per sign, summed over ALL pair-eta bins.
    Emit(out_dir + "pair_pt_counts.csv", [&](std::ostream& os) {
        header(os, "counts");
        os << "# Rows = pair-pT bins; columns = pair charge combination. Pair eta INTEGRATED over\n";
        os << "# the full binned range above (all " << neta << " bins summed).\n";
        os << "pair_pT_GeV,same_sign_pairs,opposite_sign_pairs\n";
        for (int ix = 1; ix <= npt; ++ix) {
            os << pt_lab(ix);
            for (const auto& s : signs)
                os << "," << Fmt("%.0f", s.cnt->Integral(ix, ix, 1, neta));
            os << "\n";
        }
        os << "# total," << Fmt("%.0f", signs[0].cnt->Integral())
           << "," << Fmt("%.0f", signs[1].cnt->Integral()) << "\n";
    });

    Emit(out_dir + "pair_pt_crossx.csv", [&](std::ostream& os) {
        header(os, "crossx");
        os << "# Rows = pair-pT bins. Pair eta INTEGRATED over the full binned range above.\n";
        os << "# Each sign contributes TWO columns: the cross section and its statistical error,\n";
        os << "# both in nb.\n";
        os << "pair_pT_GeV,same_sign_nb,same_sign_stat_err_nb,"
              "opposite_sign_nb,opposite_sign_stat_err_nb\n";
        for (int ix = 1; ix <= npt; ++ix) {
            os << pt_lab(ix);
            for (const auto& s : signs)
                os << "," << Fmt("%.6g", s.sw->Integral(ix, ix, 1, neta))
                   << "," << Fmt("%.4g", std::sqrt(s.sw2->Integral(ix, ix, 1, neta)));
            os << "\n";
        }
        os << "# total";
        for (const auto& s : signs)
            os << "," << Fmt("%.6g", s.sw->Integral())
               << "," << Fmt("%.4g", std::sqrt(s.sw2->Integral()));
        os << "\n";
    });

    // ==================== SETS 2 and 3: pair pT (columns) x pair eta (rows) =================
    for (const auto& s : signs) {
        Emit(out_dir + s.file_stem + "_pt_vs_eta_counts.csv", [&](std::ostream& os) {
            header(os, "counts");
            os << "# " << s.label << " pairs only. Rows = pair-eta bins, columns = pair-pT bins\n";
            os << "# [GeV]. Cell = raw pair count in that (pair pT, pair eta) cell.\n";
            os << "pair_eta\\pair_pT_GeV";
            for (int ix = 1; ix <= npt; ++ix) os << "," << pt_lab(ix);
            os << "\n";
            for (int iy = 1; iy <= neta; ++iy) {
                os << eta_lab(iy);
                for (int ix = 1; ix <= npt; ++ix)
                    os << "," << Fmt("%.0f", s.cnt->GetBinContent(ix, iy));
                os << "\n";
            }
            os << "# total (all cells) = " << Fmt("%.0f", s.cnt->Integral()) << "\n";
        });

        Emit(out_dir + s.file_stem + "_pt_vs_eta_crossx.csv", [&](std::ostream& os) {
            header(os, "crossx");
            os << "# " << s.label << " pairs only. Rows = pair-eta bins, columns = pair-pT bins\n";
            os << "# [GeV]. EACH pair-pT bin occupies TWO adjacent columns: '<range>' is the cross\n";
            os << "# section of the cell in nb, '<range>_err' is its statistical error in nb.\n";
            os << "pair_eta\\pair_pT_GeV";
            for (int ix = 1; ix <= npt; ++ix) os << "," << pt_lab(ix) << "," << pt_lab(ix) << "_err";
            os << "\n";
            for (int iy = 1; iy <= neta; ++iy) {
                os << eta_lab(iy);
                for (int ix = 1; ix <= npt; ++ix)
                    os << "," << Fmt("%.6g", s.sw->GetBinContent(ix, iy))
                       << "," << Fmt("%.4g", std::sqrt(s.sw2->GetBinContent(ix, iy)));
                os << "\n";
            }
            os << "# total (all cells) = " << Fmt("%.6g", s.sw->Integral())
               << " +- " << Fmt("%.4g", std::sqrt(s.sw2->Integral())) << " nb\n";
        });
    }

    std::cout << "\n[write_mc_pair_statistics_tables] " << cfg.key << ", " << wp_txt
              << " WP: 6 CSVs written to " << out_dir << "\n";

    for (auto& s : signs) { delete s.cnt; delete s.sw; delete s.sw2; }
}
