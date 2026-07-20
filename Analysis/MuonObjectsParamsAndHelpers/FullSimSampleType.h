#pragma once
#include <string>
#include <stdexcept>

// `noovl` = the r17663 NO-OVERLAY diagnostic sample: Pythia pp collisions reconstructed
// with the SAME PbPb23-conditions pass as the HIJING overlay (Athena 24.0.58,
// OFLCOND-MC23-SDR-RUN3-05, ConditionsRunNumber=460000, L1 MC_HI_run3_v1) but with
// `Digitization.PileUp=False` and no HIJING input. Single pT-hat slice (pTH8_14, DSID
// 802781). It exists ONLY to separate "conditions/L1 configuration" from "HIJING
// occupancy" as the cause of the forward-endcap L1 anomaly
// (mc_trigger_efficiency.md R8; round-4 outcome tree). It is NOT an overlay: no HIJING,
// no centrality, no overlay Extras.
enum class FullSimSampleType { pp, hijing, zmumu, data, noovl };

inline std::string FullSimSampleSuffix(FullSimSampleType t) {
    switch (t) {
    case FullSimSampleType::pp:     return "";
    case FullSimSampleType::hijing: return "_hijing";
    case FullSimSampleType::zmumu:  return "_zmumu";
    case FullSimSampleType::data:   return "_data";
    case FullSimSampleType::noovl:  return "_noovl";
    }
    throw std::runtime_error("FullSimSampleSuffix: unknown type");
}

// Input directory. `is_test_sample` selects the small TEST sample vs the FULL production.
// It is the SAME switch that selects the isospin treatment (FullSimSampleUsesFourBeams below)
// -- deliberately one flag, so the input path and the isospin weight can never disagree.
inline std::string FullSimSampleInputDir(FullSimSampleType t, bool is_test_sample) {
    switch (t) {
    case FullSimSampleType::pp:
        return is_test_sample
            ? "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/"
            : "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/";
    case FullSimSampleType::hijing:
        return is_test_sample
            ? "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/"
            : "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_full_sample/";
    case FullSimSampleType::zmumu:
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_zmumu_overlay_test_sample/";
    case FullSimSampleType::data:
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_data_overlay_test_sample/";
    case FullSimSampleType::noovl:
        // Only one production exists (10 k events, single slice) -- the same directory
        // either way, so the isTestSample switch cannot point it anywhere wrong.
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_no_overlay_test_sample/";
    }
    throw std::runtime_error("FullSimSampleInputDir: unknown type");
}

// "Overlay" = a minimum-bias/underlying event was overlaid on the Pythia signal, so the
// event carries a heavy-ion environment (centrality, FCal, HIJING truth). `noovl` is
// reconstructed with PbPb conditions but has NO overlaid event -> NOT an overlay.
inline bool FullSimSampleIsOverlay(FullSimSampleType t) {
    return t != FullSimSampleType::pp && t != FullSimSampleType::noovl;
}

// ISOSPIN CONTENT OF THE SIMULATED COLLISION SYSTEM -- the single rule, in one place.
//
//   pp CONDITIONS fullsim simulates pp COLLISIONS
//       -> ONE beam (pp), isospin weight 1. There is nothing to isospin-average.
//   PbPb CONDITIONS HIJING overlay simulates Pb+Pb, whose nucleons are a p/n mix
//       -> FOUR beams {pp,pn,np,nn}, combined with the Pb ratio 4:6:6:9 (Z=82, N=126).
//
// The two TEST samples are exceptions in OPPOSITE directions -- they are what happens to
// exist on disk, not what the physics wants:
//   pp24 TEST sample        : produced with 4 isospin beams BY MISTAKE  -> 4 beams
//   HIJING overlay TEST smp : only the pp beam was produced             -> 1 beam
// so the exception is exactly "is this a test sample?", and the rule collapses to a XOR:
inline bool FullSimSampleUsesFourBeams(FullSimSampleType t, bool is_test_sample) {
    // r17663: pp COLLISIONS, and only the pp-beam DSID (802781) was ever produced ->
    // one beam, isospin weight 1, regardless of the test/full switch. (It is a
    // single-slice diagnostic anyway: one global weight cancels in every ratio.)
    if (t == FullSimSampleType::noovl) return false;
    return FullSimSampleIsOverlay(t) != is_test_sample;
}

// NTUP file tag.  Frozen: it is baked into the skimmed NTUP file names on disk and
// into the grid output-dataset names (SkimCode/run_pythia_fullsim_HIJING_overlay/
// grid_sub*.sh).  "PP24" here is a legacy misnomer for the HIJING overlay -- the
// overlay is Pb+Pb (see FullSimSampleLabel) -- but renaming it would orphan the
// existing NTUPs and grid datasets.
inline std::string FullSimSampleFileTag(FullSimSampleType t) {
    switch (t) {
    case FullSimSampleType::pp:     return "FullSimPP24";
    case FullSimSampleType::hijing: return "FullSimHIJINGOverlayPP24";
    case FullSimSampleType::zmumu:  return "FullSimZmumuOverlayPP24";
    case FullSimSampleType::data:   return "FullSimDataOverlayPP24";
    // r17663: skimmed with the overlay run mode, so it inherits the overlay basename,
    // distinguished by the _r17663 suffix (grid_sub_r17663_nooverlay.sh).
    case FullSimSampleType::noovl:  return "FullSimHIJINGOverlayPP24_r17663";
    }
    throw std::runtime_error("FullSimSampleFileTag: unknown type");
}

// Output-file / plot-directory label.  The HIJING overlay simulates Pb+Pb collisions,
// never pp: the test sample (r17618 / r17662) is reconstructed with Pb+Pb 2023
// conditions (ConditionsRunNumber=460000), and the full sample now in production will
// use Pb+Pb 2024 conditions -> it will be labelled "hijing_overlay_pbpb24".
inline std::string FullSimSampleLabel(FullSimSampleType t) {
    switch (t) {
    case FullSimSampleType::pp:     return "pp24";
    case FullSimSampleType::hijing: return "hijing_overlay_pbpb23";
    case FullSimSampleType::zmumu:  return "zmumu_overlay_pp24";
    case FullSimSampleType::data:   return "data_overlay_pp24";
    case FullSimSampleType::noovl:  return "r17663_no_overlay";
    }
    throw std::runtime_error("FullSimSampleLabel: unknown type");
}

inline std::string FullSimSamplePlotDir(FullSimSampleType t) {
    switch (t) {
    case FullSimSampleType::pp:     return "pp24";
    case FullSimSampleType::hijing: return "hijing_overlay_pbpb23";
    case FullSimSampleType::zmumu:  return "zmumu_overlay_pp24";
    case FullSimSampleType::data:   return "data_overlay_pp24";
    case FullSimSampleType::noovl:  return "r17663_no_overlay";
    }
    throw std::runtime_error("FullSimSamplePlotDir: unknown type");
}
