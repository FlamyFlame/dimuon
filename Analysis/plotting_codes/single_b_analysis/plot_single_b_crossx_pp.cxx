#include "SingleBCrossxPlotterBase.cxx"

class SingleBCrossxPlotterPP : public SingleBCrossxPlotterBase {
public:
    explicit SingleBCrossxPlotterPP(const std::string& in_path)
        : SingleBCrossxPlotterBase(in_path, "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis") {}

    void Run() override {
        if (!Init()) return;

        Save2DColz("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts", "pp24_crossx_pair_pt_pair_eta_colz.png");
        Save2DColz("h2d_crossx_pair_pt_minv_w_signal_cuts", "pp24_crossx_pair_pt_minv_colz.png");
        Save2DColz("h2d_crossx_pair_pt_dr_w_signal_cuts", "pp24_crossx_pair_pt_dr_colz.png");

        DrawPairPtByEtaWithDrLines(
            "h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts",
            "PP 2024, medium WP, mu4_mu4noL1",
            "pp24_crossx_pair_pt_in_eta_subplots_dr_lines.png");
    }
};

void plot_single_b_crossx_pp(
    const std::string& input_file = "")
{
    std::string in_path = input_file;
    if (in_path.empty()) {
        const std::vector<std::string> candidates = {
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_mu4_mu4noL1_coarse_q_eta_bin.root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_mu4_mu4noL1_mu4_mu4noL1_coarse_q_eta_bin.root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_mu4_mu4noL1_mindR_noCut.root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_mu4_mu4noL1.root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_single_mu4.root"
        };
        for (const auto& p : candidates) {
            if (!gSystem->AccessPathName(p.c_str())) {
                in_path = p;
                break;
            }
        }
        if (in_path.empty()) {
            in_path = candidates.front();
        }
    }

    SingleBCrossxPlotterPP pl(in_path);
    pl.Run();
}
