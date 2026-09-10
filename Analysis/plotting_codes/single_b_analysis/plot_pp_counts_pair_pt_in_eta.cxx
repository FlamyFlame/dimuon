// =================================================================================================
// plot_pp_counts_pair_pt_in_eta.cxx
//
// RAW pp24 PAIR COUNTS in the single-b signal region: pair pT on the x axis, the 9 canonical
// pair-eta ranges as subplots, plus the same numbers as a CSV table (pair pT in COLUMNS, pair eta
// in ROWS).
//
// WHY a counts plot at all. Every point of `pp24_crossx_pair_pt_in_eta_subplots.png` is
// (1/L) * SUM 1/(eps_trig * eps_reco), so a cell can sit at a perfectly healthy dsigma/dpT while
// resting on a handful of pairs. The statistical reach of the measurement -- where the error bars
// are Poisson-dominated, how finely the result can be binned, which cells must be merged -- is a
// property of the RAW COUNT, and this is it.
//
//   input : dimuon_data/pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root
//           h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts  (nominal 9 -> 150 GeV axis)
//           h2d_counts_pt_120_pair_eta_binned_w_signal_cuts   (opt-in  9 -> 120 GeV axis)
//           = UNWEIGHTED OS pair count, filled from the very same RDF node as the crossx
//             histogram it mirrors (Tight WP + `signal_cuts`: 1.08 < m < 2.9 GeV,
//             pair pT > ParamsSet::signal_pair_pt_min, both muons outside every
//             ParamsSet::single_mu_fiducial_gap_cuts window), on the SAME axes
//             (RDFBasedHistFillingPP.cxx).
//
// BINNING: identical to `pp24_crossx_pair_pt_in_eta_subplots.png` by construction -- the SAME
// histogram (selected through SingleBCrossxPlotterBase::PtAxisHist, so the counts figure can never
// sit on a different pair-pT axis than the cross-section figure beside it) and the same nine
// panels (CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap), projected by the SAME
// SingleBCrossxPlotterBase::DrawPairPtByEta code path. Nothing is retyped here.
//
// NOT width-scaled (`differential = false`): a count is a count. That is the whole distinction
// from the crossx figure beside it.
//
// Output: plots/single_b_analysis/pp24/          (nominal, pair pT 9 -> 150 GeV)
//         plots/single_b_analysis/pp24_pt_120/   (opt-in alternative, 9 -> 120 GeV)
//           pp_counts_pair_pt_in_eta_subplots.png
//           pp_counts_pair_pt_in_eta.csv
// Usage:  root -l -b -q 'plot_pp_counts_pair_pt_in_eta.cxx+()'
// =================================================================================================

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

#include "SingleBCrossxPlotterBase.cxx"
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"   // fiducial-cut expressions, for provenance

// Repo rule: every plot set and its code expose a Medium/Tight working-point config var and
// default to TIGHT (docs/muon_wp_registry.md). Here it switches the DATA input file, which is the
// only way a WP variant is meaningful: the pp24 crossx producer writes a Medium run to a DISTINCT
// output (`RDFBasedHistFillingData.cxx`, `if (!isTight) out_file_suffix += "_medium_wp"`), so
// picking that file switches the SELECTED SPECTRUM -- and hence the count -- coherently.
enum class CountsWP { Tight, Medium };

inline CountsWP ParseCountsWP(const std::string& s)
{
    if (s == "tight"  || s == "Tight"  || s == "TIGHT")  return CountsWP::Tight;
    if (s == "medium" || s == "Medium" || s == "MEDIUM") return CountsWP::Medium;
    throw std::runtime_error("plot_pp_counts_pair_pt_in_eta: unknown working point '" + s
                             + "'; expected \"tight\" (nominal) or \"medium\".");
}
inline const char* CountsWPName(CountsWP wp){ return wp == CountsWP::Tight ? "tight" : "medium"; }
// The Tight nominal keeps the UN-suffixed filename every other crossx consumer reads.
inline std::string CountsWPFileSuffix(CountsWP wp){ return wp == CountsWP::Tight ? "" : "_medium_wp"; }

class PPSignalCountsPlotter : public SingleBCrossxPlotterBase {
    int run_year;
    CountsWP wp;

public:
    PPSignalCountsPlotter(const std::string& in_path, int yr, CountsWP wp_in)
        : SingleBCrossxPlotterBase(in_path,
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pp"
              + std::to_string(yr % 100)),
          run_year(yr % 100), wp(wp_in) {}

    // The CSV. Uses the SAME panel projection as DrawPairPtByEta -- FindBin(lo+1e-6) ..
    // FindBin(hi-1e-6) on the fine pair-eta axis -- and then PROVES it was exact: with
    // ParamsSet::N_PAIR_ETA_CROSSX_BINS = 48 (width 0.1) every coarse boundary is a bin edge, so
    // the nine panels must partition the axis and their sum must equal the 2D integral. If the
    // axis is ever moved off that property the panels silently overlap (that bug cost +20.04 %
    // in Aug 2026), so it is asserted here rather than trusted.
    void WriteCsv(const std::string& h2_name, const std::string& csv_name) {
        CheckHistogramExists(h2_name, "TH2D");
        TH2D* h2 = dynamic_cast<TH2D*>(GetHistObject(h2_name));
        if (!h2) throw std::runtime_error("Failed to retrieve TH2D: " + h2_name);

        const TAxis* ax = h2->GetXaxis();
        const int npt = ax->GetNbins();

        std::vector<std::vector<double>> rows(q_eta_bins.size(), std::vector<double>(npt, 0.));
        std::vector<double> col_tot(npt, 0.);
        double grand = 0.;

        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            const auto yb = PairEtaPanels::Bins(h2->GetYaxis(), q_eta_bins[ieta],
                                               "PPSignalCountsPlotter::WriteCsv");
            const int y1 = yb.first, y2 = yb.second;
            for (int ix = 1; ix <= npt; ++ix) {
                double n = 0.;
                for (int iy = y1; iy <= y2; ++iy) n += h2->GetBinContent(ix, iy);
                rows[ieta][ix - 1] = n;
                col_tot[ix - 1] += n;
                grand += n;
            }
        }

        // Exactness check. Since 2026-09-07 the nine panels tile the SURVIVING region
        // |eta^pair| < ParamsSet::pair_eta_fiducial_max = 2.2, NOT the whole +-2.4 fine axis --
        // so this compares the panel sum against the 2D integral only because the producer
        // applied the pair-level gap cut, leaving |eta^pair| > 2.2 empty. Panel-vs-axis
        // ALIGNMENT is now checked separately and unconditionally by PairEtaPanels::Bins above;
        // a residual here therefore means the input FILE is stale (filled before the pair-level
        // cut existed), not that a boundary is off a bin edge.
        const double total_2d = h2->Integral();
        if (std::fabs(grand - total_2d) > 1e-6 * std::max(1.0, total_2d)) {
            throw std::runtime_error(
                "PPSignalCountsPlotter::WriteCsv: the 9 pair-eta panels do not account for the "
                "whole fine pair-eta axis (panel sum " + std::to_string(grand) + " vs 2D integral "
                + std::to_string(total_2d) + "). The panels cover |eta^pair| < "
                + std::to_string(ParamsSet::pair_eta_fiducial_max) + "; a nonzero residual outside "
                "that means this input file was filled BEFORE the pair-level gap cut was added -- "
                "rerun the pp24 crossx hist filling.");
        }

        const std::string path = output_dir + "/" + csv_name;
        std::ofstream out(path);
        if (!out) throw std::runtime_error("Cannot write " + path);

        // Provenance line. The CSV travels separately from the PNG, and a 9x15 grid of bare
        // integers says nothing about what was counted. `#`-prefixed, which every CSV reader
        // worth using can skip (pandas: comment="#").
        out << "# pp 20" << run_year << " data, " << CountsWPName(wp)
            << " WP: RAW opposite-sign muon-pair COUNTS in the single-b signal region "
               "(1.08 < m_uu < 2.9 GeV, pair pT > " << ParamsSet::signal_pair_pt_min
            << " GeV, both muons outside every "
               "ParamsSet::single_mu_fiducial_gap_cuts window, and |eta^pair| < "
            << ParamsSet::pair_eta_fiducial_max
            << "). Unweighted, NOT efficiency "
               "corrected, NOT background subtracted. Columns = pair pT [GeV] ("
            << (use_pt_bins_120 ? "ParamsSet::pT_bins_120" : "ParamsSet::pT_bins_150")
            << ", " << npt << " log bins " << std::fixed << std::setprecision(0)
            << ax->GetBinLowEdge(1) << " - " << ax->GetBinUpEdge(npt)
            << " GeV); rows = the 9 "
               "CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap panels.\n";
        out << "pair_eta_range";
        for (int ix = 1; ix <= npt; ++ix)
            out << ",\"" << std::fixed << std::setprecision(2) << ax->GetBinLowEdge(ix)
                << "-" << ax->GetBinUpEdge(ix) << "\"";
        out << ",all_pt\n";

        out << std::setprecision(1);
        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            double row_tot = 0.;
            out << "\"" << std::setprecision(1) << q_eta_bins[ieta].first
                << " to " << q_eta_bins[ieta].second << "\"";
            for (int ix = 0; ix < npt; ++ix) {
                out << "," << std::setprecision(0) << rows[ieta][ix];
                row_tot += rows[ieta][ix];
            }
            out << "," << std::setprecision(0) << row_tot << "\n";
        }

        out << "\"all_eta\"";
        for (int ix = 0; ix < npt; ++ix) out << "," << std::setprecision(0) << col_tot[ix];
        out << "," << std::setprecision(0) << grand << "\n";
        out.close();

        std::cout << "[INFO] Saved: " << path << "  (total signal-region OS pairs = "
                  << std::setprecision(0) << grand << ")" << std::endl;
    }

    void Run() override {
        // Same axis discipline as the cross-section figure this one mirrors: ONE pair-pT axis per
        // invocation, into its own directory (pp24/ nominal, pp24_pt_120/ alternative).
        output_dir += PtAxisDirSuffix();
        if (!Init()) return;

        const std::string trig_label_pp = DatasetTriggerMap::GetTriggerLabel(run_year, "pp");
        const std::string h2 = PtAxisHist("h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts");
        // Built from run_year and the WP, never hard-coded: a "PP 2024, tight WP" literal would
        // silently mislabel any other year or working point the macro is asked for.
        const std::string info_line = Form("PP 20%d, %s WP", run_year, CountsWPName(wp));

        DrawPairPtByEta(h2, info_line, trig_label_pp,
                        "pp_counts_pair_pt_in_eta_subplots.png",
                        "N_{pairs}", /*differential=*/false);
        // Same asymmetry as the draw methods: absent in the partial "pt_120" alternative -> skip;
        // absent in the nominal view -> WriteCsv throws, because that would be a producer bug.
        if (!SkipInAltView(h2, "pp_counts_pair_pt_in_eta.csv"))
            WriteCsv(h2, "pp_counts_pair_pt_in_eta.csv");
    }
};

// `wp` = "tight" (nominal) or "medium" (the WP systematic variant).
// `also_pt_120 = true` additionally refreshes the opt-in 9 -> 120 GeV view in pp24_pt_120/, in the
// same invocation as the nominal 9 -> 150 GeV one, so the two cannot drift apart.
void plot_pp_counts_pair_pt_in_eta(int run_year = 24, const std::string& input_file = "",
                                   const char* wp = "tight", bool also_pt_120 = false)
{
    const CountsWP muon_wp = ParseCountsWP(wp);
    std::cout << "plot_pp_counts_pair_pt_in_eta: muon WP = " << CountsWPName(muon_wp) << std::endl;

    std::string in_path = input_file;
    const int yr = run_year % 100;
    if (in_path.empty()) {
        const std::string& trig = DatasetTriggerMap::GetTrigger(yr, "pp");
        const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_20"
                                 + std::to_string(yr) + "/";
        const std::string fname_base = base + "histograms_real_pairs_pp_20"
                                      + std::to_string(yr) + "_" + trig;
        const std::string wp_sfx = CountsWPFileSuffix(muon_wp);
        std::vector<std::string> candidates = {
            fname_base + "_nominal" + wp_sfx + ".root",
            fname_base + "_coarse_q_eta_bin" + wp_sfx + ".root",
            fname_base + wp_sfx + ".root",
        };
        for (const auto& c : candidates) {
            if (!gSystem->AccessPathName(c.c_str())) { in_path = c; break; }
        }
        if (in_path.empty()) {
            std::string msg =
                "plot_pp_counts_pair_pt_in_eta: input file not found.\n"
                "  Expected (trigger fixed by DatasetTriggerMap for pp 20" +
                std::to_string(yr) + " -> " + trig + "):\n"
                "  " + candidates.front();
            if (!wp_sfx.empty())
                msg += "\n  The Medium working point is the WP SYSTEMATIC variant. Produce it"
                       "\n  first by rerunning the pp24 data RDF crossx stage with"
                       "\n  RDFBasedHistFillingData::isTight = false (it writes the _medium_wp"
                       "\n  output; the Tight nominal is untouched).";
            throw std::runtime_error(msg);
        }
    }

    PPSignalCountsPlotter pl(in_path, yr, muon_wp);   // DEFAULT: pT_bins_150 -> pp24/
    pl.Run();

    if (also_pt_120) {
        PPSignalCountsPlotter pl120(in_path, yr, muon_wp);
        pl120.use_pt_bins_120 = true;                 // OPT-IN: pT_bins_120 -> pp24_pt_120/
        pl120.Run();
    }
}
