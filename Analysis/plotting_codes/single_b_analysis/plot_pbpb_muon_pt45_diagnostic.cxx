// =================================================================================================
// plot_pbpb_muon_pt45_diagnostic.cxx
//
// TEMPORARY DIAGNOSTIC — muon reconstructed pT > 4.5 GeV cut study, Pb+Pb 2023+2024+2025 combined.
// docs/tracking/muon_pt45_cut_diagnostic.md. Compares RAW (unweighted, uncorrected) Pb+Pb
// single-b signal-region OS pair counts under the current muon reconstructed-pT>4 GeV cut vs the
// candidate pT>4.5 GeV cut, on the SAME pair-pT x pair-eta axes as pp24
// (ParamsSet::pT_bins_120 x ParamsSet::N_PAIR_ETA_CROSSX_BINS).
//
// Input: the standalone RDFBasedHistFilling/fill_pbpb_muon_pt45_diag_counts.cxx output --
//   histograms_pbpb_23_24_25_muon_pt45_diag_counts.root
//   h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_baseline  (muon pT>4 GeV, current)
//   h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_mupt45    (muon pT>4.5 GeV, candidate)
// already combined across the 3 years (plain TChain union -- correct for an unweighted count, see
// that file's header comment) and already centrality-inclusive (no per-centrality split was
// needed for this diagnostic). NOT read from RDFBasedHistFillingPbPb's own crossx output: as of
// 2026-09 that class's trigger-weighted event loop throws for PbPb (pre-existing, unrelated bug,
// docs/tracking/pp_trig_eff_highpt_jump.md, ACTIVE, blocked on a user decision) which poisons
// EVERY histogram sharing that event loop -- see fill_pbpb_muon_pt45_diag_counts.cxx for the full
// story. TODO (not this diagnostic's scope): that bug should eventually be fixed so PbPb crossx
// hist filling (and its own raw-counts histogram) works again.
//
// If the 4.5 GeV cut is adopted, the permanent change belongs in NTuple processing
// (NTupleProcessingCode/DimuonDataAlgCoreT.c:596 for data), every data/MC result must be rerun,
// and this macro + fill_pbpb_muon_pt45_diag_counts.cxx must be DELETED.
//
// Output: plots/single_b_analysis/pbpb_23_24_25_combined/muon_pt45_diagnostic/
//           pair_eta_dependence_mupt_45_vs_40.png   (pair-pT integrated)
//           pair_pt_dependence_mupt_45_vs_40.png    (pair-eta integrated)
//           pair_pt_in_eta_subplots_mupt_45_vs_40.png
//           muon_pt45_diagnostic_summary.txt        (total counts + % decrease)
// Usage:  root -l -b -q 'plot_pbpb_muon_pt45_diagnostic.cxx+()'
// =================================================================================================

#include <fstream>
#include <stdexcept>
#include <string>

#include "SingleBCrossxPlotterBase.cxx"

class PbPbMuonPt45DiagPlotter : public SingleBCrossxPlotterBase {
public:
    explicit PbPbMuonPt45DiagPlotter(const std::string& in_path)
        : SingleBCrossxPlotterBase(in_path,
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/"
              "pbpb_23_24_25_combined/muon_pt45_diagnostic") {}

    void Run() override {
        if (!Init()) return;

        const std::string h2_a = "h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_baseline";
        const std::string h2_b = "h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_mupt45";
        const std::string label_a = "muon p_{T} > 4 GeV (current)";
        const std::string label_b = "muon p_{T} > 4.5 GeV (candidate)";
        const std::string info_line1 = "Pb+Pb 2023+2024+2025 combined, tight WP";
        const std::string info_line2 = "single-b signal region";

        DrawPairEtaIntegratedTwoSeries(h2_a, h2_b, label_a, label_b, info_line1, info_line2,
            "pair_eta_dependence_mupt_45_vs_40.png", "N_{pairs}");
        DrawPairPtIntegratedTwoSeries(h2_a, h2_b, label_a, label_b, info_line1, info_line2,
            "pair_pt_dependence_mupt_45_vs_40.png", "N_{pairs}");
        // Shorter labels for the 9-panel grid — see plot_pp24_muon_pt45_diagnostic.cxx for why.
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
        out << "Pb+Pb 2023+2024+2025 combined data, tight WP, single-b signal region\n";
        out << "RAW opposite-sign pair counts (unweighted, no trigger/reco correction):\n";
        out << "  muon pT > 4.0 GeV (current):   " << n_a << "\n";
        out << "  muon pT > 4.5 GeV (candidate): " << n_b << "\n";
        out << "  percentage decrease: " << pct << " %\n";
        out.close();
        std::cout << "[INFO] Saved: " << path << std::endl;
        std::cout << "[INFO] PbPb 23+24+25 muon pT 4.0->4.5 GeV: " << n_a << " -> " << n_b
                  << "  (" << pct << " % decrease)" << std::endl;
    }
};

void plot_pbpb_muon_pt45_diagnostic(const std::string& input_file = "")
{
    std::string in_path = input_file;
    if (in_path.empty()) {
        in_path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/"
                  "histograms_pbpb_23_24_25_muon_pt45_diag_counts.root";
        if (gSystem->AccessPathName(in_path.c_str())) {
            throw std::runtime_error(
                "plot_pbpb_muon_pt45_diagnostic: input file not found: " + in_path +
                "\n  Run RDFBasedHistFilling/fill_pbpb_muon_pt45_diag_counts.cxx first.");
        }
    }

    PbPbMuonPt45DiagPlotter pl(in_path);
    pl.Run();
}
