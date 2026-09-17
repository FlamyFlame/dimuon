// =================================================================================================
// pp24 SECONDARY-VERTEX STATISTICS
//
// What fraction of good pp24 muon pairs comes from a vertex other than the primary one?
//
// This is a STATISTICS RECORD, not a diagnostic and not a systematic: the all-vertex procedure is
// the CORRECT one for pp. pp24 has ~4 inelastic collisions per bunch crossing and the pp
// luminosity is defined over ALL of them, so a cross-section that counted only primary-vertex
// pairs would put an under-counted numerator over an all-collision L_int.
// Physics + method: docs/tracking/pp24_all_vertex_pairs.md (§1, §3d, §3e).
//
// PROVENANCE. Reads the NTuple-processing OUTPUT pair trees -- never the raw NTUPs -- so every cut
// (quality/WP, pT > 4.5, |eta| < 2.4, one-sided Delta-p/p, the same-vertex impact-parameter cut, the
// 2mu4 trigger match, the opposite-sign resonance veto) is exactly the nominal one, applied once,
// by the ntuple stage. See .claude/CLAUDE.md "NTuple-Processing Provenance".
//
// TWO FRACTIONS, deliberately reported side by side. They are different questions and conflating
// them would misstate the size of the change:
//   f_sec  = N(matched_vtx_ind != 0) / N   -- pairs whose OWN best-matching vertex is secondary.
//   f_new  = N(!pass_primary_vtx)   / N   -- pairs that fail the primary vertex outright, i.e.
//                                            the ones the all-vertex rule ADDS to the analysis.
// f_sec >= f_new always: a pair can pass w.r.t. the primary and still match a secondary vertex
// more closely.
//
// BINNING. Pair pT in the canonical coarse axis ParamsSet::pair_pt_coarse_bins
// (N_COARSE_PAIR_PT_BINS = 8 logarithmic bins, 9-150 GeV), READ from ParamsSet and never retyped
// (.claude/CLAUDE.md §Binnings). Pairs below 9 GeV and above 150 GeV are reported on their own
// rows rather than folded into the first/last bin, so the "all pairs" total is complete.
//
// POPULATIONS. Two, because the answer differs and both matter:
//   ALL PAIRS      -- everything in the output tree, the "good muon pairs" of the request.
//   SIGNAL REGION  -- the single-b OS crossx selection this feeds
//                     (Tight WP, 1.08 < m_uu < 2.9 GeV, pair pT > ParamsSet::signal_pair_pt_min
//                     (9 GeV since 2026-09-08), both muons outside every
//                      ParamsSet::single_mu_fiducial_gap_cuts window, |eta^pair| < 2.2).
//
// Usage:  root -l -b -q 'pp24_secondary_vertex_stats.cxx+()'
//         root -l -b -q 'pp24_secondary_vertex_stats.cxx+(false)'   // Medium WP
// Output: <plots>/single_b_analysis/pp24/pp24_secondary_vertex_fraction.csv  (+ a console table)
// =================================================================================================

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TSystem.h"
#include "ROOT/RDataFrame.hxx"
#include "TEfficiency.h"

#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../MuonObjectsParamsAndHelpers/DatasetTriggerMap.h"

namespace {

// Binomial error on a fraction k/n. The two samples are NESTED (k is a subset of n), so this is
// the right error, not the uncorrelated-ratio one.
//
// At the endpoints the normal-approximation form collapses to exactly 0, which in a published
// table reads as an exact measurement of zero rather than as "no events yet" -- 0/16 is not
// 0.000 +- 0.000. There the 68 % Clopper-Pearson interval is used instead and its one-sided
// width reported, so an empty high-pT bin shows an honest uncertainty.
double FracErr(double k, double n) {
    if (n <= 0) return 0.0;
    if (k <= 0)  return TEfficiency::ClopperPearson((int)n, 0,      0.6827, true);
    if (k >= n)  return 1.0 - TEfficiency::ClopperPearson((int)n, (int)n, 0.6827, false);
    const double p = k / n;
    return std::sqrt(p * (1.0 - p) / n);
}

struct Row {
    std::string label;
    double n{0}, n_sec{0}, n_new{0};
};

// At k = 0 / k = N the error is a ONE-SIDED Clopper-Pearson width, so it must not be printed as
// "+-": "0.000 +- 10.869 %" claims a two-sided interval that runs below zero.
const char* ErrSign(double k, double n) {
    if (n <= 0)          return "    ";  // 0/0 has no interval at all, one-sided or otherwise
    if (k <= 0)          return " +  ";
    if (k >= n)          return " -  ";
    return " +- ";
}

void PrintRow(std::ostream& os, const Row& r) {
    auto pct = [](double a, double b){ return b > 0 ? 100.0 * a / b : 0.0; };
    os << std::left << std::setw(22) << r.label << std::right
       << std::setw(10) << (long long)r.n
       << std::setw(9)  << (long long)r.n_sec
       << "  " << std::fixed << std::setprecision(3) << std::setw(6) << pct(r.n_sec, r.n)
       << ErrSign(r.n_sec, r.n) << std::setw(6) << 100.0 * FracErr(r.n_sec, r.n) << " %"
       << std::setw(9)  << (long long)r.n_new
       << "  " << std::setw(6) << pct(r.n_new, r.n)
       << ErrSign(r.n_new, r.n) << std::setw(6) << 100.0 * FracErr(r.n_new, r.n) << " %" << std::endl;
}

}  // namespace

// `input_override` exists only so the macro can be exercised on a `_test` NTuple output before
// the full 12-batch rerun lands; leave it empty for the nominal file.
void pp24_secondary_vertex_stats(bool use_tight_wp = true,
                                 const std::string& input_override = "") {
    const int run_year = 24;
    const std::string data_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024";
    const std::string out_dir  =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pp24";
    gSystem->mkdir(out_dir.c_str(), true);

    const ParamsSet pms{};
    const auto& pt_bins = pms.pair_pt_coarse_bins;              // CANONICAL, never retyped
    const int    n_pt   = ParamsSet::N_COARSE_PAIR_PT_BINS;
    if ((int)pt_bins.size() != n_pt + 1)
        throw std::runtime_error("pair_pt_coarse_bins size does not match N_COARSE_PAIR_PT_BINS");

    const std::string trig = DatasetTriggerMap::GetTrigger(run_year, "pp");
    const std::string path = input_override.empty()
        ? (data_dir + "/muon_pairs_pp_20" + std::to_string(run_year) + "_" + trig
           + "_mindR_0_02.root")
        : input_override;
    if (gSystem->AccessPathName(path.c_str()))
        throw std::runtime_error("Missing NTuple-processing output: " + path
                                 + " -- run the pp24 NTuple stage first.");

    // sign2 = OPPOSITE sign (_op); sign1 = same sign. The signal region is opposite-sign, and the
    // OS tree is the one the ntuple stage resonance-vetoes.
    const char* tree_os = "muon_pair_tree_sign2";
    const char* tree_ss = "muon_pair_tree_sign1";

    // WORKING POINT -- read this before quoting any number from here.
    // The pair TREE is filled at MEDIUM: DimuonDataAlgCoreT.h:363 has `requireTight = false` and
    // the nominal pp24 run script does not override it, so PassCuts_DataCore requires quality
    // bit 8 (Medium), not bit 16 (Tight). Tight is therefore an ADDITIONAL filter applied here,
    // via the `pair_pass_tight` column the ntuple stage writes.
    // Consequence, and the reason this is a named variable rather than an inline ternary: the WP
    // must be applied to EVERY population in the table, not only to the signal region. Before
    // 2026-09-08 it was applied only inside `signal_cuts`, so the "all pairs" rows were Medium
    // while the signal-region rows were Tight and the header called the whole table "Tight WP" --
    // two working points in one table under one label.
    const std::string wp_cut  = use_tight_wp ? "pair_pass_tight" : "";  // "" => tree WP (Medium)
    const std::string wp_name = use_tight_wp ? "Tight" : "Medium";

    // Compose two cut strings, tolerating an empty one (the Medium case).
    auto and_cut = [](const std::string& a, const std::string& b) {
        if (a.empty()) return b;
        if (b.empty()) return a;
        return "(" + a + ") && (" + b + ")";
    };

    // Signal region, built from ParamsSet so no gap-window or pair-eta number is ever retyped.
    // Mirrors `signal_cuts` in RDFBasedHistFilling/RDFBasedHistFillingPP.cxx, in the pair tree's
    // dotted column names.
    const std::string signal_cuts = and_cut(wp_cut,
        ParamsSet::SignalMinvCutExpr("minv") + " && " + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
        + ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && "
        + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta") + " && "
        + ParamsSet::PairFiducialEtaCutExpr("pair_eta"));

    auto count = [](ROOT::RDF::RNode df, const std::string& cut) {
        Row r;
        auto sel = cut.empty() ? df : ROOT::RDF::RNode(df.Filter(cut));
        r.n     = *sel.Count();
        r.n_sec = *sel.Filter("matched_vtx_ind != 0").Count();
        r.n_new = *sel.Filter("!pass_primary_vtx").Count();
        return r;
    };

    // ---- ALL PAIRS (OS + SS), the "good muon pairs" of the request ------------------------
    ROOT::RDataFrame df_os(tree_os, path);
    ROOT::RDataFrame df_ss(tree_ss, path);

    if (!df_os.HasColumn("matched_vtx_ind") || !df_os.HasColumn("pass_primary_vtx"))
        throw std::runtime_error(
            path + " has no matched_vtx_ind / pass_primary_vtx branch: it predates the "
            "all-vertex impact-parameter selection. Rerun the pp24 NTuple stage "
            "(docs/tracking/pp24_all_vertex_pairs.md).");

    std::vector<Row> rows;

    Row all_os = count(df_os, wp_cut);  all_os.label = "all OS pairs";
    Row all_ss = count(df_ss, wp_cut);  all_ss.label = "all SS pairs";
    Row all_both{"all pairs (OS+SS)",
                 all_os.n + all_ss.n, all_os.n_sec + all_ss.n_sec, all_os.n_new + all_ss.n_new};
    rows.push_back(all_both);
    rows.push_back(all_os);
    rows.push_back(all_ss);

    // ---- per canonical coarse pair-pT bin, OS+SS ------------------------------------------
    // Built with ostringstream, NOT snprintf into a fixed buffer: `signal_cuts` is several
    // hundred characters (it carries every gap window), and a truncated cut string does not
    // fail loudly -- RDF just refuses to JIT a half-expression.
    // The edges are serialised into the cut strings at setprecision(10), i.e. ~1e-10 GeV from
    // the exact double. That is NOT a retyped binning and cannot mis-tile: the SAME stream
    // produces bin i's upper edge and bin i+1's lower edge, so the half-open [lo, hi) bins meet
    // at a bit-identical threshold, and the residual is far below the float granularity of
    // pair_pt. `check_complete` below enforces the tiling regardless.
    auto pt_cut = [](double lo, double hi) {
        std::ostringstream o;
        o << std::setprecision(10) << "pair_pt >= " << lo << " && pair_pt < " << hi;
        return o.str();
    };
    auto bin_label = [](double lo, double hi) {
        std::ostringstream o;
        o << std::fixed << std::setprecision(2) << lo << "-" << hi;
        return o.str();
    };
    auto sum_os_ss = [&](const std::string& cut, const std::string& label) {
        const std::string c = and_cut(wp_cut, cut);   // the WP applies to every row, not just
        Row a = count(df_os, c), b = count(df_ss, c); // to the signal region
        return Row{label, a.n + b.n, a.n_sec + b.n_sec, a.n_new + b.n_new};
    };

    std::vector<Row> pt_rows;
    {
        std::ostringstream lo_cut;
        lo_cut << std::setprecision(10) << "pair_pt < " << pt_bins.front();
        pt_rows.push_back(sum_os_ss(lo_cut.str(), Form("pair pT < %g (below)", pt_bins.front())));

        for (int i = 0; i < n_pt; ++i)
            pt_rows.push_back(sum_os_ss(pt_cut(pt_bins[i], pt_bins[i + 1]),
                                        bin_label(pt_bins[i], pt_bins[i + 1])));

        std::ostringstream hi_cut;
        hi_cut << std::setprecision(10) << "pair_pt >= " << pt_bins.back();
        pt_rows.push_back(sum_os_ss(hi_cut.str(), "pair pT >= 150 (above)"));
    }

    // ---- signal region (OS only, by definition) -------------------------------------------
    Row sig = count(df_os, signal_cuts);
    sig.label = "signal region (OS)";

    std::vector<Row> sig_pt_rows;
    for (int i = 0; i < n_pt; ++i) {
        Row r = count(df_os, "(" + signal_cuts + ") && " + pt_cut(pt_bins[i], pt_bins[i + 1]));
        r.label = bin_label(pt_bins[i], pt_bins[i + 1]);
        sig_pt_rows.push_back(r);
    }
    {
        // The signal region already requires pair_pt > signal_pair_pt_min, so it has no below-axis row -- but it
        // DOES have pairs above 150 GeV, and without this row the eight bins would silently sum
        // to less than the integrated signal total, contradicting the completeness the CSV
        // header claims. (Measured 2026-09-08: 3 pairs. Small, but a table must not assert
        // something untrue about itself.)
        std::ostringstream hi_cut;
        hi_cut << std::setprecision(10) << "(" << signal_cuts << ") && pair_pt >= "
               << pt_bins.back();
        Row r = count(df_os, hi_cut.str());
        r.label = "pair pT >= 150 (above)";
        sig_pt_rows.push_back(r);
    }

    // Completeness assertion, not a comment: the per-bin rows must exhaust the integrated
    // totals, or one of the cut strings has drifted from the axis.
    auto check_complete = [](const std::vector<Row>& rows, const Row& total, const char* what) {
        double n = 0;
        for (const auto& r : rows) n += r.n;
        if (n != total.n)
            throw std::runtime_error(std::string("pp24_secondary_vertex_stats: the ") + what
                + " pair-pT rows sum to " + std::to_string((long long)n) + " but the integrated "
                "total is " + std::to_string((long long)total.n) + ". The binning rows no longer "
                "tile the sample -- check the below-/above-axis cuts against pair_pt_coarse_bins.");
    };
    check_complete(pt_rows,     all_both, "all-pairs");
    check_complete(sig_pt_rows, sig,      "signal-region");

    // ---- console table ---------------------------------------------------------------------
    auto header = [](const std::string& title) {
        std::cout << "\n" << title << "\n"
                  << std::left << std::setw(22) << "population" << std::right
                  << std::setw(10) << "N"
                  << std::setw(9)  << "N_sec" << std::setw(19) << "f_sec"
                  << std::setw(9)  << "N_new" << std::setw(19) << "f_new" << std::endl;
    };

    std::cout << "\npp24 secondary-vertex statistics, " << wp_name << " WP"
              << (use_tight_wp ? "  (the pair tree is filled at MEDIUM; Tight is applied here "
                                 "via pair_pass_tight, to EVERY row)" : "  (the tree WP)")
              << "\n"
              << "  f_sec = pairs whose best-matching vertex is NOT the primary\n"
              << "  f_new = pairs that fail the primary vertex outright (ADDED by the change)\n"
              << "  input: " << path << std::endl;

    header("--- integrated ---");
    for (const auto& r : rows) PrintRow(std::cout, r);
    PrintRow(std::cout, sig);

    header("--- all pairs (OS+SS), canonical coarse pair-pT bins [GeV] ---");
    for (const auto& r : pt_rows) PrintRow(std::cout, r);

    header("--- signal region (OS), canonical coarse pair-pT bins [GeV] ---");
    for (const auto& r : sig_pt_rows) PrintRow(std::cout, r);

    // ---- CSV --------------------------------------------------------------------------------
    const std::string csv = out_dir + "/pp24_secondary_vertex_fraction.csv";
    std::ofstream out(csv);
    if (!out) throw std::runtime_error("Cannot write " + csv);

    if (!input_override.empty())
        out << "# WARNING: produced from a NON-NOMINAL input file (" << path << ").\n";
    out << "# pp 2024 data, " << wp_name << " WP -- applied to EVERY population below, not "
           "only to the signal region. NOTE: the pair tree is filled at MEDIUM "
           "(DimuonDataAlgCoreT requireTight = false), so Tight is an additional "
           "`pair_pass_tight` filter applied at this stage. Fraction of good muon pairs "
           "coming from a "
           "SECONDARY (pile-up) vertex. pp24 has ~4 collisions per bunch crossing and the pp "
           "luminosity counts all of them, so pairs are selected against ALL reconstructed "
           "vertices: both muons must pass |d0| < " << pms.d0cut << " mm and "
           "|z0 sin(theta)| < " << pms.z0cut << " mm w.r.t. the SAME track-bearing vertex "
           "(vtx_ntrk >= 2). Read from the NTuple-processing output pair trees, so every other "
           "cut is the nominal one. See docs/tracking/pp24_all_vertex_pairs.md.\n"
        << "# f_sec = N(matched_vtx_ind != 0)/N -- the pair's own best-matching vertex is "
           "secondary. f_new = N(!pass_primary_vtx)/N -- the pair fails the primary vertex "
           "outright, i.e. it is ADDED by the all-vertex rule. f_sec >= f_new by construction.\n"
        << "# pair-pT rows use ParamsSet::pair_pt_coarse_bins (N_COARSE_PAIR_PT_BINS = "
        << n_pt << ", logarithmic, 9-150 GeV); the below-/above-axis rows keep the totals "
           "complete. Signal region = OS + " << wp_name << " + 1.08 < m_uu < 2.9 GeV + "
           "pair pT > 9 GeV + both muons outside every ParamsSet::single_mu_fiducial_gap_cuts "
           "window + |eta^pair| < " << ParamsSet::pair_eta_fiducial_max << ".\n"
        << "# The pair-pT edges in the `selection` column are ROUNDED to 2 dp for legibility; "
           "the cuts themselves use the exact ParamsSet edges. Errors are nested-binomial, "
           "except at k = 0 or k = N, where the 68 % Clopper-Pearson width is given and is "
           "ONE-SIDED -- upward at k = 0, downward at k = N (a bin with no secondary-vertex "
           "pairs is not a measurement of exactly zero, and its interval does not extend below "
           "it). These columns carry no sign, so read them with this line.\n"
        << "population,selection,N,N_sec,f_sec,f_sec_err,N_new,f_new,f_new_err\n";

    auto write_row = [&](const std::string& pop, const Row& r) {
        auto frac = [](double a, double b){ return b > 0 ? a / b : 0.0; };
        out << pop << ",\"" << r.label << "\"," << (long long)r.n << "," << (long long)r.n_sec
            << "," << std::fixed << std::setprecision(6) << frac(r.n_sec, r.n)
            << "," << FracErr(r.n_sec, r.n)
            << "," << (long long)r.n_new
            << "," << frac(r.n_new, r.n) << "," << FracErr(r.n_new, r.n) << "\n";
    };

    for (const auto& r : rows)        write_row("integrated", r);
    write_row("integrated", sig);
    for (const auto& r : pt_rows)     write_row("all_pairs_vs_pair_pt", r);
    for (const auto& r : sig_pt_rows) write_row("signal_region_vs_pair_pt", r);
    out.close();

    // ---- Two FOCUSED CSVs -------------------------------------------------------------------
    // The master CSV above holds everything; these two carry the headline tables on their own,
    // for dropping straight into the note or a talk. They are emitted from the SAME Row objects
    // in the SAME pass, so they cannot disagree with the master or with each other -- which is
    // the whole reason they are produced here rather than derived by hand afterwards.
    const std::string common_hdr =
        "# pp 2024 data, " + wp_name + " WP, applied to EVERY population in this file (the pair "
        "tree is filled at MEDIUM, so Tight is an additional pair_pass_tight filter). "
        "f_sec = the pair's own best-matching vertex is not "
        "the primary. f_new = the pair fails the primary vertex outright, i.e. it is ADDED by "
        "the all-vertex rule; f_sec >= f_new by construction. Both muons must pass |d0| and "
        "|z0 sin(theta)| w.r.t. the SAME track-bearing vertex (vtx_ntrk >= 2). Read from the "
        "NTuple-processing output pair trees. Errors are nested-binomial, except at k = 0 or "
        "k = N where the ONE-SIDED 68 % Clopper-Pearson width is given (upward at k = 0, "
        "downward at k = N); these columns carry no sign. "
        "See docs/tracking/pp24_all_vertex_pairs.md.\n";

    auto write_focused = [&](const std::string& path,
                             const std::string& extra_hdr,
                             const std::string& first_col,
                             const std::vector<std::pair<std::string, Row>>& body) {
        std::ofstream o(path);
        if (!o) throw std::runtime_error("Cannot write " + path);
        if (!input_override.empty())
            o << "# WARNING: produced from a NON-NOMINAL input file (" << path << ").\n";
        o << common_hdr << extra_hdr;
        o << first_col << ",N,N_sec,f_sec,f_sec_err,N_new,f_new,f_new_err\n";
        auto frac = [](double a, double b){ return b > 0 ? a / b : 0.0; };
        for (const auto& kv : body) {
            const Row& r = kv.second;
            o << "\"" << kv.first << "\"," << (long long)r.n << "," << (long long)r.n_sec
              << "," << std::fixed << std::setprecision(6) << frac(r.n_sec, r.n)
              << "," << FracErr(r.n_sec, r.n)
              << "," << (long long)r.n_new
              << "," << frac(r.n_new, r.n) << "," << FracErr(r.n_new, r.n) << "\n";
        }
        o.close();
        std::cout << "CSV written to " << path << std::endl;
    };

    // (1) INTEGRATED: the whole pair sample against the single-b signal region.
    write_focused(out_dir + "/pp24_secondary_vertex_integrated.csv",
        "# INTEGRATED over all pair pT. 'all pairs (OS+SS)' is every pair in the output trees; "
        "'signal region (OS)' additionally requires " + wp_name + " + 1.08 < m_uu < 2.9 GeV + "
        "pair pT > 9 GeV + both muons outside every ParamsSet::single_mu_fiducial_gap_cuts "
        "window + |eta^pair| < " + std::to_string(ParamsSet::pair_eta_fiducial_max) + ".\n",
        "population",
        { {"all pairs (OS+SS)", all_both},
          {"all OS pairs",      all_os},
          {"all SS pairs",      all_ss},
          {"signal region (OS)", sig} });

    // (2) SIGNAL REGION vs pair pT -- the table that carries the physics: f_new falls by a
    // factor ~9 across the axis, so this is a SHAPE change to dsigma/dpair-pT, not a
    // normalisation. Every bin is emitted, including the sparse high-pT ones and the
    // above-axis row, so the rows still exhaust the integrated total (check_complete above).
    {
        std::vector<std::pair<std::string, Row>> body;
        for (const auto& r : sig_pt_rows) body.emplace_back(r.label, r);
        write_focused(out_dir + "/pp24_secondary_vertex_vs_pair_pt_signal.csv",
            "# SIGNAL REGION (OS) only, split by the canonical coarse pair-pT axis "
            "ParamsSet::pair_pt_coarse_bins (N_COARSE_PAIR_PT_BINS = "
            + std::to_string(n_pt) + ", logarithmic, 9-150 GeV), read from ParamsSet and never "
            "retyped. Bin edges in the first column are ROUNDED to 2 dp for legibility; the cuts "
            "use the exact edges. The last row is the above-axis overflow, kept so the rows sum "
            "to the integrated signal-region total.\n",
            "pair_pt_range_GeV", body);
    }

    // ---- Two BRIEF CSVs, for presenting -----------------------------------------------------
    // Same numbers again, stripped to what the console prints: no counts-and-provenance columns,
    // percentages rather than fractions, and ONE header line instead of five. Emitted from the
    // same Row objects in the same pass as everything else, so the brief and detailed files
    // cannot disagree.
    //
    // The one header line is not decoration: f_sec and f_new are two similar-looking percentages
    // and a bare table of them is ambiguous, which is exactly the confusion this file set has
    // already had to correct once (the WP label). One line is the minimum that keeps it honest.
    auto write_brief = [&](const std::string& path,
                           const std::string& first_col,
                           const std::string& scope,
                           const std::vector<std::pair<std::string, Row>>& body) {
        std::ofstream o(path);
        if (!o) throw std::runtime_error("Cannot write " + path);
        o << "# pp24 " << wp_name << " WP, " << scope
          << ". f_sec = pair's best-matching vertex is not the primary; "
             "f_new = pair fails the primary outright (ADDED by the all-vertex rule). "
             "Percentages. Errors binomial; one-sided 68% CL where f = 0.\n";
        o << first_col << ",N,f_sec_percent,f_sec_err,f_new_percent,f_new_err\n";
        auto pct = [](double a, double b){ return b > 0 ? 100.0 * a / b : 0.0; };
        for (const auto& kv : body) {
            const Row& r = kv.second;
            o << "\"" << kv.first << "\"," << (long long)r.n
              << "," << std::fixed << std::setprecision(3) << pct(r.n_sec, r.n)
              << "," << 100.0 * FracErr(r.n_sec, r.n)
              << "," << pct(r.n_new, r.n)
              << "," << 100.0 * FracErr(r.n_new, r.n) << "\n";
        }
        o.close();
        std::cout << "Brief CSV written to " << path << std::endl;
    };

    write_brief(out_dir + "/pp24_secondary_vertex_integrated_brief.csv",
                "population", "integrated over all pair pT",
                { {"all pairs (OS+SS)", all_both},
                  {"all OS pairs",      all_os},
                  {"all SS pairs",      all_ss},
                  {"signal region (OS)", sig} });

    {
        std::vector<std::pair<std::string, Row>> body;
        for (const auto& r : sig_pt_rows) body.emplace_back(r.label, r);
        write_brief(out_dir + "/pp24_secondary_vertex_vs_pair_pt_signal_brief.csv",
                    "pair_pt_GeV", "single-b signal region (OS), canonical coarse pair-pT bins",
                    body);
    }

    std::cout << "\nMaster CSV written to " << csv << std::endl;
}
