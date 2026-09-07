// =================================================================================================
// plot_pp24_muon_pt45_diagnostic.cxx
//
// TEMPORARY DIAGNOSTIC — muon reconstructed pT > 4.5 GeV cut study, pp24.
// docs/tracking/muon_pt45_cut_diagnostic.md. Compares RAW (unweighted, uncorrected) pp24
// single-b signal-region OS pair counts under the current muon reconstructed-pT>4 GeV cut
// (h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts, PRE-EXISTING, untouched) vs the candidate
// pT>4.5 GeV cut (h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_mupt45, the temporary
// additive histogram in RDFBasedHistFillingPP.cxx). Same input file, same axes, same signal
// region as plot_pp_counts_pair_pt_in_eta.cxx — this is its two-series diagnostic sibling.
//
// If the 4.5 GeV cut is adopted, the permanent change belongs in NTuple processing
// (NTupleProcessingCode/DimuonDataAlgCoreT.c:596 for data), every data/MC result must be rerun,
// and this macro + its RDFBasedHistFillingPP.cxx block must be DELETED.
//
// Output: plots/single_b_analysis/pp24/muon_pt45_diagnostic/
//           pair_eta_dependence_mupt_45_vs_40.png   (pair-pT integrated)
//           pair_pt_dependence_mupt_45_vs_40.png    (pair-eta integrated)
//           pair_pt_in_eta_subplots_mupt_45_vs_40.png
//           muon_pt45_diagnostic_summary.txt        (total counts + % decrease)
// Usage:  root -l -b -q 'plot_pp24_muon_pt45_diagnostic.cxx+()'
// =================================================================================================

#include <fstream>
#include <stdexcept>
#include <string>

#include "SingleBCrossxPlotterBase.cxx"

class PPMuonPt45DiagPlotter : public SingleBCrossxPlotterBase {
    int run_year;

public:
    PPMuonPt45DiagPlotter(const std::string& in_path, int yr)
        : SingleBCrossxPlotterBase(in_path,
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pp"
              + std::to_string(yr % 100) + "/muon_pt45_diagnostic"),
          run_year(yr % 100) {}

    void Run() override {
        if (!Init()) return;

        const std::string trig_label_pp = DatasetTriggerMap::GetTriggerLabel(run_year, "pp");
        const std::string h2_base = "h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts";
        const std::string h2_a = h2_base;                    // muon pT > 4 GeV (current)
        const std::string h2_b = h2_base + "_diag_mupt45";    // muon pT > 4.5 GeV (candidate)
        const std::string label_a = "muon p_{T} > 4 GeV (current)";
        const std::string label_b = "muon p_{T} > 4.5 GeV (candidate)";
        const std::string info_line1 = Form("PP 20%d, tight WP", run_year);
        const std::string info_line2 = "single-b signal region";

        DrawPairEtaIntegratedTwoSeries(h2_a, h2_b, label_a, label_b, info_line1, info_line2,
            "pair_eta_dependence_mupt_45_vs_40.png", "N_{pairs}");
        DrawPairPtIntegratedTwoSeries(h2_a, h2_b, label_a, label_b, info_line1, info_line2,
            "pair_pt_dependence_mupt_45_vs_40.png", "N_{pairs}");
        // Shorter labels for the 9-panel grid: each subplot's legend box is narrow, and
        // "muon p_{T} > 4.5 GeV (candidate)" clips past the frame there (review 2026-09-06
        // self-check) — the marker color already carries current-vs-candidate via the other two
        // (single-panel) PNGs' full labels.
        DrawPairPtByEtaTwoSeries(h2_a, h2_b, "p_{T} > 4 GeV", "p_{T} > 4.5 GeV", info_line1, info_line2,
            "pair_pt_in_eta_subplots_mupt_45_vs_40.png", "N_{pairs}");

        TH2D* ha = dynamic_cast<TH2D*>(GetHistObject(h2_a));
        TH2D* hb = dynamic_cast<TH2D*>(GetHistObject(h2_b));
        const double n_a = ha->Integral(0, -1, 0, -1);
        const double n_b = hb->Integral(0, -1, 0, -1);
        const double pct = (n_a > 0) ? 100.0 * (n_a - n_b) / n_a : 0.0;

        const std::string path = output_dir + "/muon_pt45_diagnostic_summary.txt";
        std::ofstream out(path);
        if (!out) throw std::runtime_error("Cannot write " + path);
        out << "pp 20" << run_year << " data, tight WP, single-b signal region\n";
        out << "RAW opposite-sign pair counts (unweighted, no trigger/reco correction):\n";
        out << "  muon pT > 4.0 GeV (current):   " << n_a << "\n";
        out << "  muon pT > 4.5 GeV (candidate): " << n_b << "\n";
        out << "  percentage decrease: " << pct << " %\n";
        out.close();
        std::cout << "[INFO] Saved: " << path << std::endl;
        std::cout << "[INFO] pp20" << run_year << " muon pT 4.0->4.5 GeV: " << n_a << " -> " << n_b
                  << "  (" << pct << " % decrease)" << std::endl;
    }
};

void plot_pp24_muon_pt45_diagnostic(int run_year = 24, const std::string& input_file = "")
{
    std::string in_path = input_file;
    const int yr = run_year % 100;
    if (in_path.empty()) {
        const std::string& trig = DatasetTriggerMap::GetTrigger(yr, "pp");
        const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_20"
                                 + std::to_string(yr) + "/";
        in_path = base + "histograms_real_pairs_pp_20" + std::to_string(yr) + "_" + trig
                 + "_nominal.root";
        if (gSystem->AccessPathName(in_path.c_str())) {
            throw std::runtime_error(
                "plot_pp24_muon_pt45_diagnostic: input file not found: " + in_path +
                "\n  Run RDFBasedHistFilling/run_crossx_hist_filling_pp24.sh first.");
        }
    }

    PPMuonPt45DiagPlotter pl(in_path, yr);
    pl.Run();
}
