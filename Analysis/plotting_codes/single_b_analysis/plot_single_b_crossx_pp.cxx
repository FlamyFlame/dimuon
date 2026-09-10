#include "SingleBCrossxPlotterBase.cxx"

class SingleBCrossxPlotterPP : public SingleBCrossxPlotterBase {
    int run_year;
public:
    SingleBCrossxPlotterPP(const std::string& in_path, int yr)
        : SingleBCrossxPlotterBase(in_path,
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pp"
              + std::to_string(yr % 100)),
          run_year(yr % 100) {}

    void Run() override {
        // One invocation draws exactly ONE pair-pT axis, into its own directory: the nominal
        // 9 -> 150 GeV view (default) into pp24/, the opt-in 9 -> 120 GeV view into pp24_pt_120/.
        // Every histogram below is named ONCE, in its canonical unsuffixed form, and mapped onto
        // the selected family by PtAxisHist() -- so the whole figure set is guaranteed to sit on
        // ONE axis. Mixing them (a 9 -> 120 minv panel beside a 9 -> 150 spectrum) is exactly the
        // silent failure .claude/CLAUDE.md §Binnings exists to prevent.
        output_dir += PtAxisDirSuffix();
        if (!Init()) return;

        const std::string trig_label_pp = DatasetTriggerMap::GetTriggerLabel(run_year, "pp");

        Save2DColz(PtAxisHist("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts"),
                   "pp24_crossx_pair_pt_pair_eta.png",
                   "d^{2}#sigma/dp_{T}d#eta [pb GeV^{-1}]");
        Save2DColz(PtAxisHist("h2d_crossx_pair_pt_minv_w_signal_cuts"),
                   "pp24_crossx_pair_pt_minv.png",
                   "d^{2}#sigma/dp_{T}dm_{#mu#mu} [pb GeV^{-1} GeV^{-1}]");
        Save2DColz(PtAxisHist("h2d_crossx_pair_pt_dr_w_signal_cuts"),
                   "pp24_crossx_pair_pt_dr.png",
                   "d^{2}#sigma/dp_{T}d#DeltaR [pb GeV^{-1}]");
        DrawPairPtByEtaWithDrLines(
            PtAxisHist("h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts"),
            "PP 2024, tight WP", trig_label_pp,
            "pp24_crossx_pair_pt_in_eta_subplots_dr_lines.png",
            "d#sigma/dp_{T} [pb GeV^{-1}]");
        DrawPairPtByEta(
            PtAxisHist("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts"),
            "PP 2024, tight WP", trig_label_pp,
            "pp24_crossx_pair_pt_in_eta_subplots.png",
            "d#sigma/dp_{T} [pb GeV^{-1}]");
    }
};

// `also_pt_120 = true` additionally refreshes the OPT-IN 9 -> 120 GeV alternative view in
// pp24_pt_120/, in the SAME invocation as the nominal 9 -> 150 GeV one. Refreshing them together
// is the point: when only one of the two directories was regenerated they drifted apart for
// months (2026-04-17 to 2026-08-04, back when the 150 view was the opt-in one).
void plot_single_b_crossx_pp(
    int run_year = 24,
    const std::string& input_file = "",
    bool also_pt_120 = false)
{
    std::string in_path = input_file;
    const int yr = run_year % 100;
    if (in_path.empty()) {
        const std::string& trig = DatasetTriggerMap::GetTrigger(yr, "pp");
        const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_20"
                                 + std::to_string(yr) + "/";
        const std::string fname_base = base + "histograms_real_pairs_pp_20"
                                      + std::to_string(yr) + "_" + trig;
        std::vector<std::string> candidates = {
            fname_base + "_nominal.root",
            fname_base + "_coarse_q_eta_bin.root",
            fname_base + ".root",
        };
        for (const auto& c : candidates) {
            if (!gSystem->AccessPathName(c.c_str())) { in_path = c; break; }
        }
        if (in_path.empty()) {
            throw std::runtime_error(
                "plot_single_b_crossx_pp: input file not found.\n"
                "  Expected (trigger fixed by DatasetTriggerMap for pp 20" +
                std::to_string(yr) + " -> " + trig + "):\n"
                "  " + candidates.front());
        }
    }

    SingleBCrossxPlotterPP pl(in_path, yr);   // DEFAULT: pT_bins_150 (9 -> 150 GeV) -> pp24/
    pl.Run();

    if (also_pt_120) {
        SingleBCrossxPlotterPP pl120(in_path, yr);
        pl120.use_pt_bins_120 = true;         // OPT-IN: pT_bins_120 (9 -> 120 GeV) -> pp24_pt_120/
        pl120.Run();
    }
}
