#include "SingleBCrossxPlotterBase.cxx"

class SingleBCrossxPlotterPP : public SingleBCrossxPlotterBase {
    int run_year;
public:
    bool use_pt_bins_150 = false;

    SingleBCrossxPlotterPP(const std::string& in_path, int yr)
        : SingleBCrossxPlotterBase(in_path,
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pp"
              + std::to_string(yr % 100)),
          run_year(yr % 100) {}

    void Run() override {
        if (!Init()) return;

        const std::string trig_label_pp = DatasetTriggerMap::GetTriggerLabel(run_year, "pp");

        // --- standard pT_bins_120 plots ---
        Save2DColz("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts", "pp24_crossx_pair_pt_pair_eta.png",
                   "d^{2}#sigma/dp_{T}d#eta [pb GeV^{-1}]");
        Save2DColz("h2d_crossx_pair_pt_minv_w_signal_cuts", "pp24_crossx_pair_pt_minv.png",
                   "d^{2}#sigma/dp_{T}dm_{#mu#mu} [pb GeV^{-1} GeV^{-1}]");
        Save2DColz("h2d_crossx_pair_pt_dr_w_signal_cuts", "pp24_crossx_pair_pt_dr.png",
                   "d^{2}#sigma/dp_{T}d#DeltaR [pb GeV^{-1}]");
        DrawPairPtByEtaWithDrLines(
            "h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts",
            "PP 2024, medium WP", trig_label_pp,
            "pp24_crossx_pair_pt_in_eta_subplots_dr_lines.png",
            "d#sigma/dp_{T} [pb GeV^{-1}]");
        DrawPairPtByEta(
            "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts",
            "PP 2024, medium WP", trig_label_pp,
            "pp24_crossx_pair_pt_in_eta_subplots.png",
            "d#sigma/dp_{T} [pb GeV^{-1}]");

        if (use_pt_bins_150) {
            const std::string base_out = output_dir;
            output_dir = base_out + "_pt_150";
            gSystem->mkdir(output_dir.c_str(), true);
            Save2DColz("h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts", "pp24_crossx_pair_pt_pair_eta.png",
                       "d^{2}#sigma/dp_{T}d#eta [pb GeV^{-1}]");
            DrawPairPtByEtaWithDrLines(
                "h3d_crossx_dr_vs_pair_eta_vs_pt_150_w_signal_cuts",
                "PP 2024, medium WP", trig_label_pp,
                "pp24_crossx_pair_pt_in_eta_subplots_dr_lines.png",
                "d#sigma/dp_{T} [pb GeV^{-1}]");
            DrawPairPtByEta(
                "h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts",
                "PP 2024, medium WP", trig_label_pp,
                "pp24_crossx_pair_pt_in_eta_subplots.png",
                "d#sigma/dp_{T} [pb GeV^{-1}]");
            output_dir = base_out;
        }
    }
};

void plot_single_b_crossx_pp(
    int run_year = 24,
    const std::string& input_file = "",
    bool use_pt_bins_150 = false)
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

    SingleBCrossxPlotterPP pl(in_path, yr);
    pl.use_pt_bins_150 = use_pt_bins_150;
    pl.Run();
}
