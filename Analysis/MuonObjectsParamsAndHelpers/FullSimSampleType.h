#pragma once
#include <string>
#include <stdexcept>

enum class FullSimSampleType { pp, hijing, zmumu, data };

inline std::string FullSimSampleSuffix(FullSimSampleType t) {
    switch (t) {
    case FullSimSampleType::pp:     return "";
    case FullSimSampleType::hijing: return "_hijing";
    case FullSimSampleType::zmumu:  return "_zmumu";
    case FullSimSampleType::data:   return "_data";
    }
    throw std::runtime_error("FullSimSampleSuffix: unknown type");
}

inline std::string FullSimSampleInputDir(FullSimSampleType t) {
    switch (t) {
    case FullSimSampleType::pp:
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
    case FullSimSampleType::hijing:
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/";
    case FullSimSampleType::zmumu:
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_zmumu_overlay_test_sample/";
    case FullSimSampleType::data:
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_data_overlay_test_sample/";
    }
    throw std::runtime_error("FullSimSampleInputDir: unknown type");
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
    }
    throw std::runtime_error("FullSimSampleFileTag: unknown type");
}

inline bool FullSimSampleIsOverlay(FullSimSampleType t) {
    return t != FullSimSampleType::pp;
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
    }
    throw std::runtime_error("FullSimSampleLabel: unknown type");
}

inline std::string FullSimSamplePlotDir(FullSimSampleType t) {
    switch (t) {
    case FullSimSampleType::pp:     return "pp24";
    case FullSimSampleType::hijing: return "hijing_overlay_pbpb23";
    case FullSimSampleType::zmumu:  return "zmumu_overlay_pp24";
    case FullSimSampleType::data:   return "data_overlay_pp24";
    }
    throw std::runtime_error("FullSimSamplePlotDir: unknown type");
}
