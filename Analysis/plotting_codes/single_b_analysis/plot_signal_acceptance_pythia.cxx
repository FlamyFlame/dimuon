#include "SignalAcceptancePlotter.cxx"

class SignalAcceptancePlotterPythia : public SignalAcceptancePlotter {
public:
    SignalAcceptancePlotterPythia(const std::string& in_path)
        : SignalAcceptancePlotter(in_path,
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pythia") {}

    void Run() override {
        if (use_pt_bins_150) output_dir += "_pt_150";
        if (!Init()) return;

        Draw2DColz("pythia_sig_accept_pair_pt_vs_pair_eta.png");

        DrawAcceptancePtByEta(
            "pythia_sig_accept_pair_pt_in_eta_subplots.png",
            "Pythia 8.3 pp #sqrt{s_{NN}} = 5.36 TeV",
            "Truth single-b #rightarrow #mu^{+}#mu^{-}");
    }
};

void plot_signal_acceptance_pythia(
    double ecom = 5.36,
    bool with_data_resonance_cuts = false,
    const std::string& input_file = "",
    bool use_pt_bins_150 = false)
{
    std::string in_path = input_file;
    if (in_path.empty()) {
        std::string ecom_tag;
        std::string ecom_subdir;
        if (std::abs(ecom - 5.36) < 0.01) {
            ecom_tag   = "5p36TeV";
            ecom_subdir = "pythia_5p36TeV";
        } else if (std::abs(ecom - 5.02) < 0.01) {
            ecom_tag   = "5p02TeV";
            ecom_subdir = "pythia_5TeV";
        } else {
            throw std::runtime_error("plot_signal_acceptance_pythia: ecom must be 5.02 or 5.36 TeV.");
        }
        const std::string cut_suffix = with_data_resonance_cuts
            ? "_with_data_resonance_cuts"
            : "_no_data_resonance_cuts";
        in_path = "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/"
                  + ecom_subdir + "/histograms_pythia_" + ecom_tag + cut_suffix + ".root";
    }

    SignalAcceptancePlotterPythia pl(in_path);
    pl.use_pt_bins_150 = use_pt_bins_150;
    pl.Run();
}
