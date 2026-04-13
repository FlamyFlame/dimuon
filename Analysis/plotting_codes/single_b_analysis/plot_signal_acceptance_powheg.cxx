#include "SignalAcceptancePlotter.cxx"

class SignalAcceptancePlotterPowheg : public SignalAcceptancePlotter {
public:
    SignalAcceptancePlotterPowheg(const std::string& in_path)
        : SignalAcceptancePlotter(in_path,
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/powheg") {}

    void Run() override {
        if (use_pt_bins_150) output_dir += "_pt_150";
        if (!Init()) return;

        Draw2DColz("powheg_sig_accept_pair_pt_vs_pair_eta.png");

        DrawAcceptancePtByEta(
            "powheg_sig_accept_pair_pt_in_eta_subplots.png",
            "POWHEG pp #sqrt{s_{NN}} = 5.36 TeV",
            "Truth single-b #rightarrow #mu^{+}#mu^{-}");
    }
};

void plot_signal_acceptance_powheg(const std::string& input_file = "",
                                    bool use_pt_bins_150 = false)
{
    std::string in_path = input_file;
    if (in_path.empty()) {
        in_path = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/histograms_powheg_truth.root";
    }

    SignalAcceptancePlotterPowheg pl(in_path);
    pl.use_pt_bins_150 = use_pt_bins_150;
    pl.Run();
}
