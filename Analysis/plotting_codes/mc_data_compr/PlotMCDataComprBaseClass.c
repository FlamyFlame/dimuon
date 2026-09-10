#include "PlotMCDataComprBaseClass.h"

#include <stdexcept>

// -------------------------------------------------------------------------------------------
// THE observable table. `data_var` / `mc_var` / `powheg_var` are the ONLY place a producer
// histogram name is spelled out; the key builders in the header assemble the rest. Nothing here
// retypes a BINNING: the axes come from the histograms, which the producers book by name from
// `hist_binning_map` / `ParamsSet` on ALL THREE sides, and AssertSameAxis1D() enforces that they
// really did.
//
// `powheg_var` is the Pythia `mc_var` minus its `truth_` prefix -- the POWHEG truth producer
// names its columns without it. An EMPTY `powheg_var` means the sample has no counterpart for
// that observable (pair eta and pair pT are SIGNAL-family observables; POWHEG appears only in
// the generic family).
//
// AXIS TITLES: `#DeltaR`, not `#Delta R`. The space is not ignored by TLatex -- it renders as
// "Δ R", which is wrong for a single symbol. Same for #Delta#phi / #Delta#eta.
// -------------------------------------------------------------------------------------------
const std::vector<PlotMCDataComprBaseClass::McDataObservable>&
PlotMCDataComprBaseClass::Observables(){
    static const std::vector<McDataObservable> v = {
        // kin           x title              y symbol       y unit         data_var       mc_var                     powheg_var           logx
        {"DR",           "#DeltaR",           "#DeltaR",     "pb",          "DR",          "truth_dr_ppbin",          "dr_ppbin",          false},
        {"DR_zoomin",    "#DeltaR",           "#DeltaR",     "pb",          "DR_zoomin",   "truth_dr_zoomin_ppbin",   "dr_zoomin_ppbin",   false},
        {"Dphi",         "#Delta#phi",        "#Delta#phi",  "pb",          "Dphi",        "truth_dphi_ppbin",        "dphi_ppbin",        false},
        {"Dphi_zoomin",  "#Delta#phi",        "#Delta#phi",  "pb",          "Dphi_zoomin", "truth_dphi_zoomin_ppbin", "dphi_zoomin_ppbin", false},
        {"Deta_zoomin",  "#Delta#eta",        "#Delta#eta",  "pb",          "Deta_zoomin", "truth_deta_zoomin_ppbin", "deta_zoomin_ppbin", false},
        {"minv_zoomin",  "m_{#mu#mu} [GeV]",  "m_{#mu#mu}",  "pb GeV^{-1}", "minv_zoomin", "truth_minv_zoomin_ppbin", "minv_zoomin_ppbin", false},
        {"pair_eta",     "#eta^{pair}",       "#eta^{pair}", "pb",          "pair_eta",    "truth_pair_eta_crossx",   "",                  false},
        // pair pT: BOTH sides are on the NOMINAL 9 -> 150 GeV axis. The data variable
        // `pair_pt_150` (RDFBasedHistFillingPP signal_region_1d_vars) and the MC variable
        // `truth_pair_pt_log_150` (var1D_{pythia,powheg}_fullsim.json) are both bound to
        // "binning": "pT_bins_150"; neither has a second, alternative-axis member, so there is
        // nothing to select between here. The "_150" is a leftover naming token, not an
        // opt-in marker -- the opt-in 9 -> 120 GeV view is the separate "pt_120" family and it
        // does not reach these 1D spectra.
        {"pair_pt",      "p_{T}^{pair} [GeV]","p_{T}^{pair}","pb GeV^{-1}", "pair_pt_150", "truth_pair_pt_log_150",   "",                  true}
    };
    return v;
}

const PlotMCDataComprBaseClass::McDataObservable&
PlotMCDataComprBaseClass::Observable(const std::string& kin){
    for (const auto& o : Observables()) if (o.kin == kin) return o;
    throw std::runtime_error("PlotMCDataComprBaseClass::Observable: unknown observable '" + kin
                             + "'. Add it to Observables() in PlotMCDataComprBaseClass.c.");
}

void PlotMCDataComprBaseClass::initialize(){
    std::cout << "Base class initialization function (muon WP = "
              << McDataComprConfig::WPName(muon_wp) << ")" << std::endl;

    pythia_with_data_resonance_cuts_suffix = pythia_with_data_resonance_cuts? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";

    // Input paths and file names come from McDataComprConfig -- ONE place, shared with
    // plot_mc_data_pair_pt_in_eta.cxx, which does not derive from this class.
    dt_paths = {
        McDataComprConfig::PowhegDir(),
        // The pp24-condition FULLSIM FULL sample (pp beam only, isospin weight 1, AMI-weighted),
        // NOT the private sample and NOT the Pythia truth full sample -- see the header comment.
        McDataComprConfig::PythiaDir(),
        McDataComprConfig::DataDir()
    };
    fnames = {
        McDataComprConfig::PowhegFileName(),
        McDataComprConfig::PythiaFileName(pythia_with_data_resonance_cuts),
        McDataComprConfig::DataFileName(muon_wp)
    };

    // Normalizations, in one place.
    //   Pythia: nb -> pb only. The old `TRUTH_COMBINE_RENORM = 1e3` hand factor belonged to the
    //           private truth-combined sample and is gone with it.
    //   POWHEG: exactly 1 -- its weight is already a pb cross-section. See PowhegNormFactor()
    //           for why it must NEVER be given Pythia's x1000.
    //   Data:   1/L_int, except in the SIGNAL family, where the histograms are already
    //           d(sigma) in pb (the crossx weight carries 1/L) and a second division would be
    //           a silent double count.
    norm_factor[DataType::pythia] = (float)NbToPb();
    norm_factor[DataType::powheg] = (float)PowhegNormFactor();
    norm_factor[DataType::pp_2024_2mu4] =
        data_already_dsigma ? 1.f : (float)PPBaseClass::GetCrossxFactor(24, "2mu4");

    for (int idt = 0; idt < s_nDtTypes; idt++){
        const std::string path = dt_paths[idt] + fnames[idt];
        std::cout << "opening file: " << path << std::endl;
        McDataComprConfig::AssertInputExists(path, is_data[idt] ? "data" : "MC");
        f[idt] = TFile::Open(path.c_str());
        if (!f[idt] || f[idt]->IsZombie()){
            throw std::runtime_error("Cannot open file: " + path);
        }
    }

    set_legend_position();
}

void PlotMCDataComprBaseClass::set_legend_position(){
    // With a LOG y axis there IS no corner that is empty for every shape: dR fills the top, dphi
    // is U-shaped and fills the bottom corners. So on log y the sentinel is deliberately LEFT IN
    // PLACE and the caller builds a default-constructed TLegend, which ROOT auto-places in a gap
    // it finds itself (the same fix used for the gap-cut figures, muon_gap_cuts_acceptance.md
    // F13). A position the caller set explicitly still wins. The auto-placed box is then pushed
    // clear of the top frame line by McDataComprRatio::FixLegendTopOverlap().
    if (logy) return;
    if (legend_position_same_sign[0] < -999.){
        legend_position_same_sign = {0.5,0.6,0.95,0.89};
    }
    if (legend_position_opp_sign[0] < -999.){
        legend_position_opp_sign = {0.5,0.6,0.93,0.89};
    }
}

std::string PlotMCDataComprBaseClass::OutputPath(const std::string& filename) const {
    std::string dir = PlotBaseDir();
    if (!output_subdir.empty()){
        dir += output_subdir;
        if (dir.back() != '/') dir += '/';
    }
    gSystem->mkdir(dir.c_str(), kTRUE);
    return dir + filename;
}

void PlotMCDataComprBaseClass::AssertSameAxis1D(const TH1* a, const TH1* b,
                                                const std::string& ctx){
    if (!a || !b) return;
    const TAxis* xa = a->GetXaxis();
    const TAxis* xb = b->GetXaxis();
    if (xa->GetNbins() != xb->GetNbins())
        throw std::runtime_error(ctx + ": axis bin count differs ("
                                 + std::to_string(xa->GetNbins()) + " vs "
                                 + std::to_string(xb->GetNbins()) + ")");
    for (int i = 1; i <= xa->GetNbins() + 1; ++i){
        const double ea = xa->GetBinLowEdge(i);
        const double eb = xb->GetBinLowEdge(i);
        // Relative tolerance: one side is built as (nbins, min, max) and the other from a
        // generated edge vector, and those agree only to double rounding.
        if (std::fabs(ea - eb) > 1e-9 * std::max(1.0, std::fabs(ea)))
            throw std::runtime_error(ctx + ": axis edge " + std::to_string(i) + " differs ("
                                     + std::to_string(ea) + " vs " + std::to_string(eb) + ")");
    }
}

TH1D* PlotMCDataComprBaseClass::GetHist1D(int idt, const std::string& base_name) const {
    if (base_name.empty())
        throw std::runtime_error("GetHist1D: empty histogram name requested for " + fnames[idt]
                                 + " (this observable has no counterpart in that sample)");
    TH1D* h = (TH1D*)f[idt]->Get(base_name.c_str());
    if (!h) {
        throw std::runtime_error("GetHist1D: histogram '" + base_name + "' not found in " + fnames[idt]);
    }
    TH1D* c = (TH1D*)h->Clone(Form("h_%s_%d", base_name.c_str(), idt));
    c->SetDirectory(nullptr);
    return c;
}
