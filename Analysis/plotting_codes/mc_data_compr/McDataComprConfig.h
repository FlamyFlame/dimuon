#pragma once

#include <stdexcept>
#include <string>

#include <TSystem.h>

// =================================================================================================
// McDataComprConfig -- the ONE place where the mc_data_compr plot set names its input files and
// its muon working point.
//
// It lives at namespace scope, not as members of PlotMCDataComprBaseClass, for the same reason
// McDataComprColors.h does: `plot_mc_data_pair_pt_in_eta.cxx` draws the same SIGNAL family without
// deriving from that class, and a protected member is unreachable from there (which is how the
// colour constants once broke Stage 8 of pipeline_pp_crossx.sh).
//
// WORKING POINT (repo rule: every plot set exposes a Medium/Tight config var, default TIGHT).
// What it actually switches, and what it CANNOT switch, is spelled out at MuonWP below -- read it
// before quoting a "Medium" version of any of these plots.
// =================================================================================================
namespace McDataComprConfig {

// -------------------------------------------------------------------------------------------
// The muon working point. NOMINAL = Tight (2026-07-07 decision; docs/muon_wp_registry.md).
//
// SWITCHES: the DATA input file only. The pp24 crossx producer routes a Medium run to a
//   DISTINCT output (`RDFBasedHistFillingData.cxx:137`, `if (!isTight) out_file_suffix +=
//   "_medium_wp"`), because Medium changes BOTH the selected spectrum (quality bit 8 instead of
//   16) AND the efficiencies it is corrected with (the Medium turn-on fits and the Medium
//   reco-eff placeholder). Picking the Medium file therefore switches the whole data side
//   coherently -- selection and correction together -- which is the only way a WP variant is
//   meaningful.
//
// CANNOT SWITCH: the MC side, and this is physics, not a missing feature. Every MC histogram in
//   this plot set is a TRUTH-quantity histogram (`*_pass_signal_truth`, `*_gapcut_truth`, in
//   truth kinematics, with the gap cut on truth q*eta). A reconstruction working point is a
//   property of a RECONSTRUCTED muon; it does not exist for a truth pair. The Pythia fullsim and
//   POWHEG files accordingly carry no WP token in their names. So a Medium run of this plot set
//   compares the Medium DATA against the SAME MC curve -- which is exactly what the WP systematic
//   asks for.
// -------------------------------------------------------------------------------------------
enum class MuonWP { Tight, Medium };

inline MuonWP ParseWP(const std::string& s){
    if (s == "tight" || s == "Tight" || s == "TIGHT") return MuonWP::Tight;
    if (s == "medium" || s == "Medium" || s == "MEDIUM") return MuonWP::Medium;
    throw std::runtime_error("McDataComprConfig::ParseWP: unknown working point '" + s
                             + "'; expected \"tight\" (nominal) or \"medium\".");
}

inline const char* WPName(MuonWP wp){ return wp == MuonWP::Tight ? "Tight" : "Medium"; }

// The Tight nominal keeps the UN-suffixed filename that every other crossx consumer reads.
inline std::string WPFileSuffix(MuonWP wp){ return wp == MuonWP::Tight ? "" : "_medium_wp"; }

// --- input directories ---------------------------------------------------------------------
inline const char* DataDir(){
    return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/";
}
// The pp24-condition FULLSIM FULL sample: pp beam only, isospin weight 1, AMI-weighted.
inline const char* PythiaDir(){
    return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/";
}
inline const char* PowhegDir(){
    return "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
}

// --- input files -----------------------------------------------------------------------------
inline std::string DataFileName(MuonWP wp){
    return "histograms_real_pairs_pp_2024_2mu4_nominal" + WPFileSuffix(wp) + ".root";
}
inline std::string DataFile(MuonWP wp){ return std::string(DataDir()) + DataFileName(wp); }

inline std::string PythiaFileName(bool with_data_resonance_cuts){
    return std::string("histograms_pythia_fullsim_pp24")
         + (with_data_resonance_cuts ? "_with_data_resonance_cuts" : "_no_data_resonance_cuts")
         + "_full.root";
}
inline std::string PythiaFile(bool with_data_resonance_cuts){
    return std::string(PythiaDir()) + PythiaFileName(with_data_resonance_cuts);
}

// POWHEG truth, ONE combined file (2026-08-25 rework). The plot set used to open the two 2025
// legacy files `{bb,cc}_evgen_truth_full_sample/histograms_mc_truth_{bb,cc}_combined.root`, which
// the current producer (`RDFBasedHistFillingPowhegTruth.cxx`) no longer writes at all.
inline std::string PowhegFileName(){ return "histograms_powheg_truth.root"; }
inline std::string PowhegFile(){ return std::string(PowhegDir()) + PowhegFileName(); }

// A missing Medium data file is the one failure mode a user will actually hit, so say what to do
// about it instead of letting TFile::Open return a null pointer.
inline void AssertInputExists(const std::string& path, const std::string& what){
    if (gSystem->AccessPathName(path.c_str())){
        std::string msg = what + " input does not exist: " + path;
        if (path.find("_medium_wp") != std::string::npos)
            msg += "\n  The Medium working point is the WP SYSTEMATIC variant. Produce it first by"
                   "\n  rerunning the pp24 data RDF crossx stage with RDFBasedHistFillingData::"
                   "isTight = false"
                   "\n  (it writes the _medium_wp output; the Tight nominal is untouched).";
        throw std::runtime_error(msg);
    }
}

}  // namespace McDataComprConfig
