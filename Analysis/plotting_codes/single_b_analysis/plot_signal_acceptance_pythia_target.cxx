#include "SignalAcceptancePlotter.cxx"

// TARGET signal acceptance (docs/tracking/signal_acceptance_mass_window.md PD4/D5-D7) from the private
// no-muon-pT-filter Pythia sample: denominator T_target (single-b OS, truth pair pT > 9,
// |eta^pair| < 2.2, both |eta_mu| <= 2.5, no muon pT cut). One plot set per factor:
//   total A_sig, the single-muon pT cut, the single-muon eta cut, both muon cuts, the mass window.
// Output: its own directory (never the interim pythia/ plots).
class SignalAcceptancePlotterPythiaTarget : public SignalAcceptancePlotter {
public:
    std::string file_tag;   // "" for the total, else the factor tag

    SignalAcceptancePlotterPythiaTarget(const std::string& in_path, const std::string& tag, const std::string& ytitle)
        : SignalAcceptancePlotter(in_path,
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pythia_target")
    {
        num_tag         = tag;
        file_tag        = tag.empty() ? "total" : tag;
        y_title         = ytitle;
        expected_marker = SignalAcceptanceDefinition::TargetValue();
    }

    void Run() override {
        output_dir += PtAxisDirSuffix();
        if (!Init()) return;
        const std::string mc1 = "Pythia 8 A14 pp #sqrt{s} = 5.36 TeV";
        const std::string mc2 = "Truth single-b #rightarrow #mu^{+}#mu^{-}";
        const std::string base = "pythia_target_sig_accept_" + file_tag;
        Draw2DColz(base + "_pair_pt_vs_pair_eta.png");
        DrawAcceptancePtByEta(base + "_pair_pt_in_eta_subplots.png", mc1, mc2);
        if (!use_pt_bins_120) DrawAcceptanceEtaByPt(base + "_pair_eta_in_pt_subplots.png", mc1, mc2);
    }
};

void plot_signal_acceptance_pythia_target(
    const std::string& input_file =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_private_nomupt/histograms_pythia_nomupt_muon_fid_target_no_data_resonance_cuts.root")
{
    // Labels from the ParamsSet constants the numerators use -- never retyped.
    const std::string pt_cut   = Form("p_{T}^{#mu} > %g GeV", ParamsSet::SINGLE_MU_PT_MIN);
    const std::string eta_cut  = Form("|#eta^{#mu}| < %g", ParamsSet::SINGLE_MU_ETA_MAX);
    const std::string mass_cut = Form("%g < m_{#mu#mu} < %g GeV", ParamsSet::signal_minv_min, ParamsSet::signal_minv_max);
    const std::vector<std::pair<std::string, std::string>> factors = {
        {"",      "#alpha_{signal}"},
        {"mupt",  "#alpha(" + pt_cut + ")"},
        {"mueta", "#alpha(" + eta_cut + ")"},
        {"mu",    "#alpha(" + pt_cut + ", " + eta_cut + ")"},
        {"mass",  "#alpha(" + mass_cut + ")"}};
    for (const auto& [tag, ytitle] : factors) {
        SignalAcceptancePlotterPythiaTarget pl(input_file, tag, ytitle);
        pl.Run();
    }
}
