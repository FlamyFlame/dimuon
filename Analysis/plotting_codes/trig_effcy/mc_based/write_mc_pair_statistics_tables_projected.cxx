// =================================================================================================
// write_mc_pair_statistics_tables_projected.cxx
//
// PROJECTED-STATISTICS companion to write_mc_pair_statistics_tables.cxx (2026-09-03 user request;
// docs/tracking/pythia_pp24_pthat_stats_projection.md). Answers: on the pair-pT x pair-eta
// statistics grid of the Step-3 DeltaR-correction sample, what would each cell's RAW PAIR COUNT
// look like if the two highest pT-hat slices (kn4 = 70-125 GeV, kn5 = 125-300 GeV) had
// PtHatKn45Projected::kNTarget events each instead of the kNCurrent they were actually produced
// with -- a decision aid for an MC production request, NOT a new measurement.
//
// SOURCE (never re-derived): the SAME Step-3 statistics histograms
// write_mc_pair_statistics_tables.cxx reads (h_mc_paircount_vs_pt_eta_{ss,os}, booked by
// FillMCTrigEffHists.cxx on the Step-3 denominator node), PLUS the new per-slice counterparts
// h_mc_paircount_vs_pt_eta_{ss,os}_kn{4,5} that the SAME fill stage now ALSO books (opt-in
// `book_kn_stats` flag) from the per-pT-hat-slice trees the mc_trig pair file already carries.
//
// PROJECTION (tracking doc Physics Procedure §4), per (pair-pT bin, pair-eta GROUP, sign):
//   count_proj = count_all_slices - count_kn4 - count_kn5 + sf*(count_kn4 + count_kn5)
// i.e. remove the two slices' CURRENT contribution and put back their PROJECTED one; kn0-3
// untouched. sf = PtHatKn45Projected::kSf (SAME constant the projected crossx plot uses).
//
// BINNING: pair-pT keeps the filled 8-bin axis (ParamsSet::pair_pt_coarse_bins) untouched. Pair-eta
// is MERGED into the 3 sign-independent |eta^pair| bins 0-1, 1-2, 2-2.4 via MakeDrEtaGroups
// (dr_correction_cell_groups.h) -- the SAME utility just added for the dR-correction fit cells
// (mode "nocorr_etamerge"); never re-invented here, per .claude/CLAUDE.md 'Binnings'.
//
// FOUR CSVs per (sample, WP): for each sign (same-sign, opposite-sign),
//   {sign}_pt_vs_abs_eta_counts_projected.csv  the projected count (formula above)
//   {sign}_pt_vs_abs_eta_counts_current.csv    the SAME cells, SAME binning, un-projected (all
//                                               six pT-hat slices as actually produced) -- the
//                                               direct baseline for a side-by-side comparison.
//
// Compile/run (ACLiC, from this directory):
//   root -l -b -q 'write_mc_pair_statistics_tables_projected.cxx+("pp_full", true)'   // Tight
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

#include "dr_correction_sample_cfg.h"
#include "dr_correction_cell_groups.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"
#include "../../../Utilities/PtHatKn45ProjectedStats.h"
#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"

namespace {

std::string Fmt(const char* f, double a)            { char b[256]; snprintf(b, sizeof b, f, a);    return b; }
std::string Fmt(const char* f, double a, double b_) { char b[256]; snprintf(b, sizeof b, f, a, b_); return b; }

TH2D* Get2D(TFile* f, const std::string& name)
{
    TH2D* h = dynamic_cast<TH2D*>(f->Get(name.c_str()));
    if (!h)
        throw std::runtime_error(
            "write_mc_pair_statistics_tables_projected: histogram '" + name + "' is MISSING from "
            + f->GetName() + ".\n  The per-pTHat-slice statistics histograms are booked by "
              "FillMCTrigEffHists.cxx's book_kn_stats option; rerun Step 3 with book_kn_stats=true "
              "for this sample/WP (docs/tracking/pythia_pp24_pthat_stats_projection.md).");
    h->SetDirectory(nullptr);
    return h;
}

void CheckSameAxes(const TH2D* ref, const TH2D* h)
{
    const bool same = h->GetNbinsX() == ref->GetNbinsX() && h->GetNbinsY() == ref->GetNbinsY();
    if (!same)
        throw std::runtime_error(std::string("write_mc_pair_statistics_tables_projected: binning "
                                 "mismatch between '") + ref->GetName() + "' and '" + h->GetName() + "'");
    for (int i = 1; i <= ref->GetNbinsX() + 1; ++i)
        if (std::fabs(h->GetXaxis()->GetBinLowEdge(i) - ref->GetXaxis()->GetBinLowEdge(i)) > 1e-6)
            throw std::runtime_error("write_mc_pair_statistics_tables_projected: pair-pT edge "
                                     + std::to_string(i) + " differs between '" + ref->GetName()
                                     + "' and '" + h->GetName() + "'");
    for (int i = 1; i <= ref->GetNbinsY() + 1; ++i)
        if (std::fabs(h->GetYaxis()->GetBinLowEdge(i) - ref->GetYaxis()->GetBinLowEdge(i)) > 1e-6)
            throw std::runtime_error("write_mc_pair_statistics_tables_projected: pair-eta edge "
                                     + std::to_string(i) + " differs between '" + ref->GetName()
                                     + "' and '" + h->GetName() + "'");
}

// Sum a TH2D's bin content over pair-pT bin `ix` (1-based, identity axis) and the pair-eta GROUP
// `g` (1-based) of `G`, i.e. over every (possibly folded) source sub-range G.ranges[g-1] holds.
double GroupSum(const TH2D* h, int ix, const DrAxisGroups& G, int g)
{
    double s = 0.;
    for (const auto& r : G.ranges[g - 1])
        s += h->Integral(ix, ix, r.first, r.second);
    return s;
}

void Emit(const std::string& path, const std::function<void(std::ostream&)>& write)
{
    std::cout << "\n----- " << path << " -----\n";
    write(std::cout);
    std::ofstream ofs(path);
    if (!ofs)
        throw std::runtime_error("write_mc_pair_statistics_tables_projected: cannot open '" + path
                                 + "' for writing");
    write(ofs);
}

}  // namespace

void write_mc_pair_statistics_tables_projected(const std::string& sample = "pp_full",
                                               bool use_tight_wp = true)
{
    gROOT->SetBatch(kTRUE);

    const DrCorrSample cfg    = GetDrCorrSample(sample, use_tight_wp);
    const std::string  wp_suf = DrCorrWpSuffix(use_tight_wp);
    const std::string  wp_txt = use_tight_wp ? "Tight" : "Medium";
    const double        sf    = PtHatKn45Projected::kSf;

    const std::string in_path = cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf
                              + MCTrigEffPairPt::FileSuffix() + "_step3.root";
    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie())
        throw std::runtime_error("write_mc_pair_statistics_tables_projected: cannot open " + in_path);

    struct SignSet { std::string key, label, file_stem; TH2D *all, *kn4, *kn5; };
    std::vector<SignSet> signs = {
        {"ss", "same sign",     "same_sign",     nullptr, nullptr, nullptr},
        {"os", "opposite sign", "opposite_sign", nullptr, nullptr, nullptr}
    };
    for (auto& s : signs) {
        s.all = Get2D(fin, "h_mc_paircount_vs_pt_eta_" + s.key);
        s.kn4 = Get2D(fin, "h_mc_paircount_vs_pt_eta_" + s.key + "_kn4");
        s.kn5 = Get2D(fin, "h_mc_paircount_vs_pt_eta_" + s.key + "_kn5");
    }
    fin->Close();

    const TH2D* ref = signs[0].all;
    for (const auto& s : signs) { CheckSameAxes(ref, s.all); CheckSameAxes(ref, s.kn4); CheckSameAxes(ref, s.kn5); }

    const int npt = ref->GetNbinsX();
    if (npt != MCTrigEffPairPt::NBins())
        throw std::runtime_error("write_mc_pair_statistics_tables_projected: " + in_path + " has "
            + std::to_string(npt) + " pair-pT bins but this process is configured for "
            + std::to_string(MCTrigEffPairPt::NBins()) + " (" + MCTrigEffPairPt::Describe() + ")");

    // Pair-eta MERGE into the 3 sign-independent |eta^pair| bins -- reused verbatim, never
    // re-invented (dr_correction_cell_groups.h; .claude/CLAUDE.md 'Binnings').
    const DrAxisGroups Geta = MakeDrEtaGroups(ref->GetYaxis(), /*merge_eta=*/true);
    std::cout << DrGroupsDescribe(Geta, "pair eta", " GeV") << "\n";

    const TAxis* ax = ref->GetXaxis();
    auto pt_lab = [&](int ix) { return Fmt("%.1f-%.1f", ax->GetBinLowEdge(ix), ax->GetBinUpEdge(ix)); };
    auto eta_lab = [&](int g) { return Fmt("%.1f to %.1f |eta|", Geta.edges[g - 1], Geta.edges[g]); };

    const std::string out_dir = cfg.mc_dir;   // sits beside the source Step-3 file, like its input

    auto header = [&](std::ostream& os, const std::string& sign_label) {
        os << "# PROJECTED same-sign/opposite-sign muon-pair statistics of the Step-3 "
              "DeltaR-correction sample (mc_trigger_efficiency.md 3.3), on the SAME selection as\n";
        os << "# write_mc_pair_statistics_tables.cxx, with the pair-eta axis MERGED into the 3\n";
        os << "# sign-independent |pair eta| bins and the two highest pT-hat slices' contribution\n";
        os << "# PROJECTED to " << Fmt("%.0f", PtHatKn45Projected::kNTarget) << " events each "
              "(currently " << Fmt("%.0f", PtHatKn45Projected::kNCurrent) << " each, sf = "
           << Fmt("%.5f", sf) << "). Central-value cross sections are UNCHANGED by this "
              "projection -- see docs/tracking/pythia_pp24_pthat_stats_projection.md sect 1. This "
              "is a DECISION AID for an MC production request, not a measurement.\n";
        os << "#\n";
        os << "# " << sign_label << " pairs only. QUANTITY: PROJECTED raw pair count (UNWEIGHTED,\n";
        os << "#   DIMENSIONLESS) = count_all_slices - count_kn4 - count_kn5 + sf*(count_kn4 + "
              "count_kn5).\n";
        os << "# sample=" << cfg.key << "  label=" << cfg.mc_label << "  WP=" << wp_txt
           << " (required of BOTH legs)\n";
        os << "# source = " << in_path << "\n";
        os << "# sign convention: muon_pair_tree_sign1 = SAME sign, muon_pair_tree_sign2 = "
              "OPPOSITE sign.\n";
        os << "# pair-pT bin edges [GeV]:";
        for (int ix = 1; ix <= npt + 1; ++ix) os << " " << Fmt("%.4f", ax->GetBinLowEdge(ix));
        os << "\n";
        os << "# pair-eta GROUP edges (|eta^pair|, folded neg+pos source bins):";
        for (double e : Geta.edges) os << " " << Fmt("%.4f", e);
        os << "\n#\n";
    };

    // CURRENT (un-projected) companion, same selection/binning/eta-merge, all six pT-hat slices
    // as actually produced -- for a DIRECT comparison against the *_projected.csv cell by cell.
    // Not a projection: no sf anywhere here, this is simply GroupSum(s.all, ...) on its own.
    auto header_current = [&](std::ostream& os, const std::string& sign_label) {
        os << "# CURRENT (un-projected, all six pT-hat slices as actually produced) same-sign/\n";
        os << "# opposite-sign muon-pair statistics of the Step-3 DeltaR-correction sample\n";
        os << "# (mc_trigger_efficiency.md 3.3), on the SAME selection as\n";
        os << "# write_mc_pair_statistics_tables.cxx, with the pair-eta axis MERGED into the 3\n";
        os << "# sign-independent |pair eta| bins -- the direct baseline for\n";
        os << "# *_pt_vs_abs_eta_counts_projected.csv (same cells, same binning, sf NOT applied).\n";
        os << "#\n";
        os << "# " << sign_label << " pairs only. QUANTITY: raw pair count (UNWEIGHTED, "
              "DIMENSIONLESS).\n";
        os << "# sample=" << cfg.key << "  label=" << cfg.mc_label << "  WP=" << wp_txt
           << " (required of BOTH legs)\n";
        os << "# source = " << in_path << "\n";
        os << "# sign convention: muon_pair_tree_sign1 = SAME sign, muon_pair_tree_sign2 = "
              "OPPOSITE sign.\n";
        os << "# pair-pT bin edges [GeV]:";
        for (int ix = 1; ix <= npt + 1; ++ix) os << " " << Fmt("%.4f", ax->GetBinLowEdge(ix));
        os << "\n";
        os << "# pair-eta GROUP edges (|eta^pair|, folded neg+pos source bins):";
        for (double e : Geta.edges) os << " " << Fmt("%.4f", e);
        os << "\n#\n";
    };

    for (const auto& s : signs) {
        Emit(out_dir + s.file_stem + "_pt_vs_abs_eta_counts_current.csv", [&](std::ostream& os) {
            header_current(os, s.label);
            os << "# Rows = |pair eta| GROUP; columns = pair-pT bins [GeV]. Cell = CURRENT raw\n";
            os << "# pair count, all six pT-hat slices as actually produced (no projection).\n";
            os << "pair_abs_eta\\pair_pT_GeV";
            for (int ix = 1; ix <= npt; ++ix) os << "," << pt_lab(ix);
            os << "\n";
            double grand_cur = 0.;
            for (int g = 1; g <= Geta.n; ++g) {
                os << eta_lab(g);
                for (int ix = 1; ix <= npt; ++ix) {
                    const double all = GroupSum(s.all, ix, Geta, g);
                    grand_cur += all;
                    os << "," << Fmt("%.0f", all);
                }
                os << "\n";
            }
            os << "# total: current=" << Fmt("%.0f", grand_cur) << "\n";
        });
    }

    for (const auto& s : signs) {
        Emit(out_dir + s.file_stem + "_pt_vs_abs_eta_counts_projected.csv", [&](std::ostream& os) {
            header(os, s.label);
            os << "# Rows = |pair eta| GROUP; columns = pair-pT bins [GeV]. Cell = PROJECTED raw\n";
            os << "# pair count (see the formula above). The DIRECT, cell-by-cell CURRENT\n";
            os << "# (un-projected) companion, same selection/binning/eta-merge, is the sibling\n";
            os << "# file *_pt_vs_abs_eta_counts_current.csv.\n";
            os << "pair_abs_eta\\pair_pT_GeV";
            for (int ix = 1; ix <= npt; ++ix) os << "," << pt_lab(ix);
            os << "\n";
            double grand_cur = 0., grand_proj = 0.;
            for (int g = 1; g <= Geta.n; ++g) {
                os << eta_lab(g);
                for (int ix = 1; ix <= npt; ++ix) {
                    const double all = GroupSum(s.all, ix, Geta, g);
                    const double k4  = GroupSum(s.kn4, ix, Geta, g);
                    const double k5  = GroupSum(s.kn5, ix, Geta, g);
                    const double proj = all - k4 - k5 + sf * (k4 + k5);
                    grand_cur  += all;
                    grand_proj += proj;
                    os << "," << Fmt("%.1f", proj);
                }
                os << "\n";
            }
            os << "# total: current=" << Fmt("%.0f", grand_cur) << "  projected="
               << Fmt("%.1f", grand_proj) << "  (x" << Fmt("%.3f", grand_proj / grand_cur) << ")\n";
        });
    }

    std::cout << "\n[write_mc_pair_statistics_tables_projected] " << cfg.key << ", " << wp_txt
              << " WP: 4 CSVs written to " << out_dir << "\n";

    for (auto& s : signs) { delete s.all; delete s.kn4; delete s.kn5; }
}
