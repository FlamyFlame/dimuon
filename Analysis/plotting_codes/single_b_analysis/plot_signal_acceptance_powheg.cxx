#include "SignalAcceptancePlotter.cxx"

class SignalAcceptancePlotterPowheg : public SignalAcceptancePlotter {
public:
    SignalAcceptancePlotterPowheg(const std::string& in_path)
        : SignalAcceptancePlotter(in_path,
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/powheg") {}

    void Run() override {
        output_dir += PtAxisDirSuffix();
        if (!Init()) return;

        Draw2DColz("powheg_sig_accept_pair_pt_vs_pair_eta.png");

        DrawAcceptancePtByEta(
            "powheg_sig_accept_pair_pt_in_eta_subplots.png",
            "POWHEG pp #sqrt{s_{NN}} = 5.36 TeV",
            "Truth single-b #rightarrow #mu^{+}#mu^{-}");
    }
};

// `also_pt_120 = true` additionally refreshes the opt-in 9 -> 120 GeV view in powheg_pt_120/.
void plot_signal_acceptance_powheg(const std::string& input_file = "",
                                    bool also_pt_120 = false)
{
    std::string in_path = input_file;
    if (in_path.empty()) {
        in_path = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/histograms_powheg_truth.root";
    }

    SignalAcceptancePlotterPowheg pl(in_path);   // DEFAULT: pT_bins_150 -> powheg/
    pl.Run();

    if (also_pt_120) {
        SignalAcceptancePlotterPowheg pl120(in_path);
        pl120.use_pt_bins_120 = true;            // OPT-IN: pT_bins_120 -> powheg_pt_120/
        pl120.Run();
    }
}
