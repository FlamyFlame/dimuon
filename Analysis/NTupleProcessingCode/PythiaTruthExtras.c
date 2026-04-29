#include "PythiaTruthExtras.h"
#include "PythiaAlgCoreT.h"
#include "TChain.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <cmath>

// ===========================================================================
// Init helpers
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::InitParamsExtra() {
    self().perform_truth = true;
    ResonanceNameMap();
}

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::ResonanceNameMap() {
    resonance_ids = {113, 221, 223, 331, 333, 443, 553};
    resonance_id_to_name_and_crossx_map[113] = {"rho0",    0.};
    resonance_id_to_name_and_crossx_map[221] = {"eta",     0.};
    resonance_id_to_name_and_crossx_map[223] = {"omega",   0.};
    resonance_id_to_name_and_crossx_map[331] = {"eta'",    0.};
    resonance_id_to_name_and_crossx_map[333] = {"phi",     0.};
    resonance_id_to_name_and_crossx_map[443] = {"J/psi",   0.};
    resonance_id_to_name_and_crossx_map[553] = {"Upsilon", 0.};
}

// ---------------------------------------------------------------------------
// InitInputExtra — set truth particle branches on all relevant chains
// ---------------------------------------------------------------------------

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::InitInputExtra() {

    // Ensure vector<vector<int>> CollectionProxy exists before SetBranchAddress.
    // Without this, ROOT returns -3 ("no compiled CollectionProxy") in sessions
    // where PythiaTruthExtras.h was not JIT-compiled by Cling.
    if (!isPrivate)
        gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");

    // Helper lambda: bind truth particle branches onto one chain
    auto bind_truth = [&](TChain* ch) {
        if (!ch) return;
        ch->SetBranchStatus("truth_id",       1); ch->SetBranchAddress("truth_id",      &truth_id);
        ch->SetBranchStatus("truth_barcode",  1); ch->SetBranchAddress("truth_barcode", &truth_barcode);
        ch->SetBranchStatus("truth_status",   1); ch->SetBranchAddress("truth_status",  &truth_status);

        if (isPrivate) {
            ch->SetBranchStatus("truth_m",   1); ch->SetBranchAddress("truth_m",   &truth_m);
            ch->SetBranchStatus("truth_pt",  1); ch->SetBranchAddress("truth_pt",  &truth_pt);
            ch->SetBranchStatus("truth_eta", 1); ch->SetBranchAddress("truth_eta", &truth_eta);
            ch->SetBranchStatus("truth_phi", 1); ch->SetBranchAddress("truth_phi", &truth_phi);
            ch->SetBranchStatus("truth_mother1", 1); ch->SetBranchAddress("truth_mother1", &truth_mother1);
            ch->SetBranchStatus("truth_mother2", 1); ch->SetBranchAddress("truth_mother2", &truth_mother2);
        } else {
            ch->SetBranchStatus("truth_m",   1); ch->SetBranchAddress("truth_m",   &truth_m_f);
            ch->SetBranchStatus("truth_pt",  1); ch->SetBranchAddress("truth_pt",  &truth_pt_f);
            ch->SetBranchStatus("truth_eta", 1); ch->SetBranchAddress("truth_eta", &truth_eta_f);
            ch->SetBranchStatus("truth_phi", 1); ch->SetBranchAddress("truth_phi", &truth_phi_f);
            ch->SetBranchStatus("truth_parents", 1); ch->SetBranchAddress("truth_parents", &truth_parents);
        }
    };

    if (isPrivate) {
        bind_truth(self().evChain);
    } else {
        for (int ikin = 0; ikin < self().nKinRanges; ikin++)
            for (int ibeam = 0; ibeam < self().nBeamTypes; ibeam++)
                bind_truth(self().evChains_kn_beam[ikin][ibeam]);
    }
}

// ---------------------------------------------------------------------------
// InitTempVariablesExtra
// ---------------------------------------------------------------------------

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::InitTempVariablesExtra() {
    m1_history          = &m1_history_storage;
    m2_history          = &m2_history_storage;
    m1_history_particle = &m1_history_particle_storage;
    m2_history_particle = &m2_history_particle_storage;
    m1_first_hf_hadron_prt_pt_eta_phi_m = &m1_first_hf_hadron_prt_pt_eta_phi_m_storage;
    m2_first_hf_hadron_prt_pt_eta_phi_m = &m2_first_hf_hadron_prt_pt_eta_phi_m_storage;
    CrossxClear();
}

// ---------------------------------------------------------------------------
// InitOutputHistsExtra
// ---------------------------------------------------------------------------

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::InitOutputHistsExtra() {
    h_numMuonPairs = new TH1D("h_numMuonPairs", "h_numMuonPairs", 6, 0, 6);

    for (int isign = 0; isign < ParamsSet::nSigns; isign++) {
        for (int jdphi = 0; jdphi < 2; jdphi++) {
            h_parent_groups[isign][jdphi] = new TH2D(
                Form("h_parent_groups_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_parent_groups_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                nParentGroups, 0, nParentGroups, nParentGroups, 0, nParentGroups);

            h_hard_scatt_categr_vs_origin_categr[isign][jdphi] = new TH2D(
                Form("h_hard_scatt_categr_vs_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_hard_scatt_categr_vs_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                nHFMuonPairOriginCategories, 0, nHFMuonPairOriginCategories,
                nHardScattCategories, 0, nHardScattCategories);

            h_both_from_b_hard_scatt_categr_vs_origin_categr[isign][jdphi] = new TH2D(
                Form("h_both_from_b_hard_scatt_categr_vs_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_both_from_b_hard_scatt_categr_vs_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                nHFMuonPairOriginCategories, 0, nHFMuonPairOriginCategories,
                nHardScattCategories, 0, nHardScattCategories);

            h_both_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi] = new TH2D(
                Form("h_both_from_c_hard_scatt_categr_vs_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_both_from_c_hard_scatt_categr_vs_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                nHFMuonPairOriginCategories, 0, nHFMuonPairOriginCategories,
                nHardScattCategories, 0, nHardScattCategories);

            h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi] = new TH2D(
                Form("h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                nHFMuonPairOriginCategories, 0, nHFMuonPairOriginCategories,
                nHardScattCategories, 0, nHardScattCategories);

            h_muon_pair_origin_categr[isign][jdphi] = new TH1D(
                Form("h_muon_pair_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_muon_pair_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                nHFMuonPairOriginCategories, 0, nHFMuonPairOriginCategories);

            h_both_from_b_muon_pair_origin_categr[isign][jdphi] = new TH1D(
                Form("h_both_from_b_muon_pair_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_both_from_b_muon_pair_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                nHFMuonPairOriginCategories, 0, nHFMuonPairOriginCategories);

            h_both_from_c_muon_pair_origin_categr[isign][jdphi] = new TH1D(
                Form("h_both_from_c_muon_pair_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_both_from_c_muon_pair_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                nHFMuonPairOriginCategories, 0, nHFMuonPairOriginCategories);

            h_one_from_b_one_from_c_muon_pair_origin_categr[isign][jdphi] = new TH1D(
                Form("h_one_from_b_one_from_c_muon_pair_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_one_from_b_one_from_c_muon_pair_origin_categr_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                nHFMuonPairOriginCategories, 0, nHFMuonPairOriginCategories);

            h_Qsplit_gs_ISR_one_hard_scatt[isign][jdphi] = new TH1D(
                Form("h_Qsplit_gs_ISR_one_hard_scatt_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_Qsplit_gs_ISR_one_hard_scatt_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                80, 0, 80);

            h_Qsplit_gs_FSR[isign][jdphi] = new TH1D(
                Form("h_Qsplit_gs_FSR_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_Qsplit_gs_FSR_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                80, 0, 80);

            h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[isign][jdphi] = new TH1D(
                Form("h_Qsplit_to_mHard_gs_ISR_one_hard_scatt_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_Qsplit_to_mHard_gs_ISR_one_hard_scatt_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                50, 0, 1.);

            h_Qsplit_to_mHard_gs_FSR[isign][jdphi] = new TH1D(
                Form("h_Qsplit_to_mHard_gs_FSR_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                Form("h_Qsplit_to_mHard_gs_FSR_sign%d%s", isign+1, jdphi ? "_away" : "_near"),
                50, 0, 1.);

            // Enable Sumw2
            h_parent_groups[isign][jdphi]->Sumw2();
            h_hard_scatt_categr_vs_origin_categr[isign][jdphi]->Sumw2();
            h_both_from_b_hard_scatt_categr_vs_origin_categr[isign][jdphi]->Sumw2();
            h_both_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi]->Sumw2();
            h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi]->Sumw2();
            h_muon_pair_origin_categr[isign][jdphi]->Sumw2();
            h_both_from_b_muon_pair_origin_categr[isign][jdphi]->Sumw2();
            h_both_from_c_muon_pair_origin_categr[isign][jdphi]->Sumw2();
            h_one_from_b_one_from_c_muon_pair_origin_categr[isign][jdphi]->Sumw2();
            h_Qsplit_gs_ISR_one_hard_scatt[isign][jdphi]->Sumw2();
            h_Qsplit_gs_FSR[isign][jdphi]->Sumw2();
            h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[isign][jdphi]->Sumw2();
            h_Qsplit_to_mHard_gs_FSR[isign][jdphi]->Sumw2();
        }
    }
}

// ---------------------------------------------------------------------------
// InitOutputExtra — open log files
// ---------------------------------------------------------------------------

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::InitOutputExtra() {
    std::string diag_dir = self().py_dir + "ancestor_tracing_output_diagnostic_files/";
    const std::string& bsuf = self().batch_suffix;
    const std::string sign_labels[ParamsSet::nSigns] = {"same_sign", "op_sign"};
    const std::string dphis[2] = {"_near", "_away"};

    if (print_bad_warnings) {
        m_very_bad_warning_file = new std::ofstream(diag_dir + "very_bad_warnings_pythia" + bsuf + ".txt");
        *m_very_bad_warning_file << "Test writing into warning files [beginning]" << std::endl;
        m_hard_scattering_warning_file = new std::ofstream(diag_dir + "hard_scattering_warnings_pythia" + bsuf + ".txt");
    }
    m_crossx_summary_file = new std::ofstream(diag_dir + "crossx_summary_pythia" + bsuf + ".txt");

    if (print_low_minv_resonances)
        m_very_low_minv_resonance_file = new std::ofstream(diag_dir + "very_low_minv_resonances_pythia" + bsuf + ".txt");

    if (print_unspecified_parent)
        m_unspecified_parent_file = new std::ofstream(diag_dir + "unspecified_parents_pythia" + bsuf + ".txt");

    if (print_FE)
        m_FE_file = new std::ofstream(diag_dir + "FE_pythia" + bsuf + ".txt");

    if (print_HF_pair_origin_others_history)
        m_HF_pair_origin_others_category_file = new std::ofstream(diag_dir + "HF_pair_origin_others_category" + bsuf + ".txt");

    if (print_other_flavor_history)
        m_other_flavor_category_file = new std::ofstream(diag_dir + "other_flavor_category" + bsuf + ".txt");

    if (print_prt_history) {
        for (int isign = 0; isign < ParamsSet::nSigns; isign++)
            for (int jdphi = 0; jdphi < 2; jdphi++) {
                m_b_parent_file[isign][jdphi] = new std::ofstream(
                    Form("%sb_parents_%s%s%s.txt", diag_dir.c_str(), sign_labels[isign].c_str(), dphis[jdphi].c_str(), bsuf.c_str()));
                m_c_parent_file[isign][jdphi] = new std::ofstream(
                    Form("%sc_parents_%s%s%s.txt", diag_dir.c_str(), sign_labels[isign].c_str(), dphis[jdphi].c_str(), bsuf.c_str()));
            }
    }

    if (debug_print_history_nevents > 0)
        m_debug_history_file = new std::ofstream(diag_dir + "debug_history" + bsuf + ".txt");
    if (debug_print_bhadron_id_nevents > 0)
        m_debug_bhadron_id_file = new std::ofstream(diag_dir + "debug_bhadron_id" + bsuf + ".txt");
}

// ---------------------------------------------------------------------------
// HistAdjustExtra — set bin labels
// ---------------------------------------------------------------------------

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::HistAdjustExtra() {
    for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++) {
        for (int lphi = 0; lphi < 2; lphi++) {
            for (int iprt = 0; iprt < nParentGroups; iprt++) {
                h_parent_groups[ksign][lphi]->GetXaxis()->SetBinLabel(iprt+1, parentGroupLabels[iprt].c_str());
                h_parent_groups[ksign][lphi]->GetYaxis()->SetBinLabel(iprt+1, parentGroupLabels[iprt].c_str());
            }
            for (int iorigin = 0; iorigin < nHFMuonPairOriginCategories; iorigin++) {
                const char* lbl = HFMuonPairOriginCategoryLabels[iorigin].c_str();
                h_muon_pair_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1, lbl);
                h_both_from_b_muon_pair_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1, lbl);
                h_both_from_c_muon_pair_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1, lbl);
                h_one_from_b_one_from_c_muon_pair_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1, lbl);
                h_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1, lbl);
                h_both_from_b_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1, lbl);
                h_both_from_c_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1, lbl);
                h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1, lbl);
            }
            for (int ihard = 0; ihard < nHardScattCategories; ihard++) {
                const char* lbl = HardScattCategoryLabels[ihard].c_str();
                h_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetYaxis()->SetBinLabel(ihard+1, lbl);
                h_both_from_b_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetYaxis()->SetBinLabel(ihard+1, lbl);
                h_both_from_c_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetYaxis()->SetBinLabel(ihard+1, lbl);
                h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetYaxis()->SetBinLabel(ihard+1, lbl);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// FinalizeExtra
// ---------------------------------------------------------------------------

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::FinalizeExtra() {
    WriteCrossxSummary();

    // Always print truncated-history stats to stdout for quick inspection
    std::cout << "[TruthExtras] #muons from B-hadron: " << m_n_bhadron_muons << std::endl;
    if (m_n_bhadron_muons > 0) {
        double frac = 100. * m_n_truncated_bhadron_muons / m_n_bhadron_muons;
        std::cout << "[TruthExtras] #muons with truncated B-hadron history: "
                  << m_n_truncated_bhadron_muons << " (" << frac << "%)" << std::endl;
        if (m_n_truncated_bhadron_muons > 0)
            std::cout << "[TruthExtras]   eldest B-hadron barcode range: ["
                      << m_truncated_bc_min << ", " << m_truncated_bc_max << "]" << std::endl;
    }

    if (m_crossx_summary_file) { m_crossx_summary_file->close(); delete m_crossx_summary_file; m_crossx_summary_file = nullptr; }

    if (print_bad_warnings) {
        if (m_very_bad_warning_file) {
            *m_very_bad_warning_file << "Test writing into warning files [end]" << std::endl;
            m_very_bad_warning_file->close(); delete m_very_bad_warning_file; m_very_bad_warning_file = nullptr;
        }
        if (m_hard_scattering_warning_file) {
            m_hard_scattering_warning_file->close(); delete m_hard_scattering_warning_file; m_hard_scattering_warning_file = nullptr;
        }
    }
    if (print_low_minv_resonances && m_very_low_minv_resonance_file) {
        m_very_low_minv_resonance_file->close(); delete m_very_low_minv_resonance_file; m_very_low_minv_resonance_file = nullptr;
    }
    if (print_unspecified_parent && m_unspecified_parent_file) {
        m_unspecified_parent_file->close(); delete m_unspecified_parent_file; m_unspecified_parent_file = nullptr;
    }
    if (print_FE && m_FE_file) {
        m_FE_file->close(); delete m_FE_file; m_FE_file = nullptr;
    }
    if (print_HF_pair_origin_others_history && m_HF_pair_origin_others_category_file) {
        m_HF_pair_origin_others_category_file->close(); delete m_HF_pair_origin_others_category_file; m_HF_pair_origin_others_category_file = nullptr;
    }
    if (print_other_flavor_history && m_other_flavor_category_file) {
        m_other_flavor_category_file->close(); delete m_other_flavor_category_file; m_other_flavor_category_file = nullptr;
    }
    if (print_prt_history) {
        for (int isign = 0; isign < ParamsSet::nSigns; isign++)
            for (int jdphi = 0; jdphi < 2; jdphi++) {
                if (m_b_parent_file[isign][jdphi]) { m_b_parent_file[isign][jdphi]->close(); delete m_b_parent_file[isign][jdphi]; m_b_parent_file[isign][jdphi] = nullptr; }
                if (m_c_parent_file[isign][jdphi]) { m_c_parent_file[isign][jdphi]->close(); delete m_c_parent_file[isign][jdphi]; m_c_parent_file[isign][jdphi] = nullptr; }
            }
    }
}

// ===========================================================================
// Kinematic accessors
// ===========================================================================

template <class PairT, class Derived>
int PythiaTruthExtras<PairT, Derived>::GetParticleIndex(int barcode) const {
    if (barcode < 0)
        throw std::runtime_error("GetParticleIndex: negative barcode " + std::to_string(barcode));

    if (isPrivate) {
        if (!truth_status)
            throw std::runtime_error("GetParticleIndex: truth_status is null for private sample.");
        if (barcode >= static_cast<int>(truth_status->size()))
            throw std::runtime_error("GetParticleIndex: private barcode out of range: " + std::to_string(barcode));
        return barcode;
    }

    if (!truth_barcode)
        throw std::runtime_error("GetParticleIndex: truth_barcode is null for non-private sample.");

    int bc_back = truth_barcode->empty() ? -1 : truth_barcode->back();
    if (barcode_lookup_source != truth_barcode
        || barcode_lookup_size != truth_barcode->size()
        || barcode_lookup_back != bc_back) {
        barcode_to_index_cache.clear();
        for (size_t i = 0; i < truth_barcode->size(); ++i)
            barcode_to_index_cache.emplace(truth_barcode->at(i), static_cast<int>(i));
        barcode_lookup_source = truth_barcode;
        barcode_lookup_size = truth_barcode->size();
        barcode_lookup_back = bc_back;
    }

    auto it = barcode_to_index_cache.find(barcode);
    if (it == barcode_to_index_cache.end())
        throw std::runtime_error("GetParticleIndex: barcode not found in truth_barcode: " + std::to_string(barcode));
    return it->second;
}

template <class PairT, class Derived>
double PythiaTruthExtras<PairT, Derived>::TruthPtAt(size_t i) const {
    const int idx = GetParticleIndex(static_cast<int>(i));
    if (isPrivate && truth_pt)    return truth_pt->at(idx);
    if (!isPrivate && truth_pt_f) return static_cast<double>(truth_pt_f->at(idx));
    return 0.;
}

template <class PairT, class Derived>
double PythiaTruthExtras<PairT, Derived>::TruthEtaAt(size_t i) const {
    const int idx = GetParticleIndex(static_cast<int>(i));
    if (isPrivate && truth_eta)    return truth_eta->at(idx);
    if (!isPrivate && truth_eta_f) return static_cast<double>(truth_eta_f->at(idx));
    return 0.;
}

template <class PairT, class Derived>
double PythiaTruthExtras<PairT, Derived>::TruthPhiAt(size_t i) const {
    const int idx = GetParticleIndex(static_cast<int>(i));
    if (isPrivate && truth_phi)    return truth_phi->at(idx);
    if (!isPrivate && truth_phi_f) return static_cast<double>(truth_phi_f->at(idx));
    return 0.;
}

template <class PairT, class Derived>
double PythiaTruthExtras<PairT, Derived>::TruthMAt(size_t i) const {
    const int idx = GetParticleIndex(static_cast<int>(i));
    if (isPrivate && truth_m)    return truth_m->at(idx);
    if (!isPrivate && truth_m_f) return static_cast<double>(truth_m_f->at(idx));
    return 0.;
}

// ===========================================================================
// Cross-section tracking helpers
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::CrossxClear() {
    pair_counter = 0;
    total_crossx = 0.;
    from_resonance_total_crossx = 0.;
    from_same_b_total_crossx = 0.;
    either_from_tau_total_crossx = 0.;
    both_incoming_total_crossx = 0.;
    both_incoming_FE_QQ_from_same_b_total_crossx = 0.;
    FE_total_crossx = 0.;
    from_same_gluon_spitting_total_crossx = 0.;
    hard_QQ_scatt_total_crossx = 0.;
    FE_from_same_b_total_crossx = 0.;
    FE_from_same_GS_total_crossx = 0.;
    FE_from_diff_ancestors_total_crossx = 0.;
    FE_from_same_ancestors_not_same_b_or_gs_total_crossx = 0.;
    skipped_event_crossx = 0.;
}

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::PerPairCrossxUpdate() {
    if (!truth_id) return;
    PairT* mpair = self().mpairRef().get();
    if (!mpair) return;

    total_crossx += mpair->weight;
    pair_counter++;

    const int idx3 = GetParticleIndex(3);
    const int idx4 = GetParticleIndex(4);
    bool parton1_hq = (abs(truth_id->at(idx3)) == 4 || abs(truth_id->at(idx3)) == 5);
    bool parton2_hq = (abs(truth_id->at(idx4)) == 4 || abs(truth_id->at(idx4)) == 5);
    bool parton1_glq = (abs(truth_id->at(idx3)) <= 3 || abs(truth_id->at(idx3)) == 21);
    bool parton2_glq = (abs(truth_id->at(idx4)) <= 3 || abs(truth_id->at(idx4)) == 21);
    bool flavor_excitation = (parton1_hq && parton2_glq) || (parton2_hq && parton1_glq);
    bool hard_QQ_scatt = parton1_hq && parton2_hq;
    bool both_incoming = m1_ancestor_is_incoming && m2_ancestor_is_incoming;

    if (m1_from_tau || m2_from_tau) either_from_tau_total_crossx += mpair->weight;
    if (hard_QQ_scatt) hard_QQ_scatt_total_crossx += mpair->weight;

    if (flavor_excitation) {
        FE_total_crossx += mpair->weight;
        if (mpair->from_same_b)
            FE_from_same_b_total_crossx += mpair->weight;
        else if (from_same_gluon_photon_splitting_or_both_HQ_incoming)
            FE_from_same_GS_total_crossx += mpair->weight;
        else if (!mpair->from_same_ancestors)
            FE_from_diff_ancestors_total_crossx += mpair->weight;
        else {
            FE_from_same_ancestors_not_same_b_or_gs_total_crossx += mpair->weight;
            if (print_FE && m_FE_file) {
                *m_FE_file << "Event #: " << mpair->m1.ev_num
                           << ". FE, same ancestors, not same b or GS." << std::endl;
                if (truth_id->size() > 6)
                    *m_FE_file << truth_id->at(GetParticleIndex(3)) << " " << truth_id->at(GetParticleIndex(4)) << " "
                               << truth_id->at(GetParticleIndex(5)) << " " << truth_id->at(GetParticleIndex(6)) << " "
                               << mpair->from_same_ancestors << std::endl;
                PrintHistory(m_FE_file, false, mpair->from_same_ancestors);
            }
        }
    }

    if (mpair->from_same_b) from_same_b_total_crossx += mpair->weight;
    if (from_same_gluon_photon_splitting_or_both_HQ_incoming) from_same_gluon_spitting_total_crossx += mpair->weight;

    if (both_incoming) {
        both_incoming_total_crossx += mpair->weight;
        if ((flavor_excitation && mpair->from_same_b) || hard_QQ_scatt)
            both_incoming_FE_QQ_from_same_b_total_crossx += mpair->weight;
    }
}

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::WriteCrossxSummary() {
    if (!m_crossx_summary_file) return;
    *m_crossx_summary_file << "#muon pairs: " << pair_counter << std::endl;
    *m_crossx_summary_file << "#muons from B-hadron: " << m_n_bhadron_muons << std::endl;
    if (m_n_bhadron_muons > 0) {
        double frac = 100. * m_n_truncated_bhadron_muons / m_n_bhadron_muons;
        *m_crossx_summary_file
            << "#muons from B-hadron with truncated history (no b-quark found): "
            << m_n_truncated_bhadron_muons << " (" << frac << "%)" << std::endl;
        if (m_n_truncated_bhadron_muons > 0)
            *m_crossx_summary_file
                << "  eldest B-hadron barcode range for truncated muons: ["
                << m_truncated_bc_min << ", " << m_truncated_bc_max << "]" << std::endl;
    }
    *m_crossx_summary_file << "Total crossx: " << total_crossx << std::endl;
    *m_crossx_summary_file << "Total crossx of skipped HF events: " << skipped_event_crossx << std::endl;
    *m_crossx_summary_file << "Total crossx [from resonance]: " << from_resonance_total_crossx << std::endl;
    *m_crossx_summary_file << "Total crossx [from same b]: " << from_same_b_total_crossx << std::endl;
    *m_crossx_summary_file << "Total crossx [BOTH INCOMING]: " << both_incoming_total_crossx << std::endl;
    *m_crossx_summary_file << "Total crossx [both incoming] - FE & same b: " << both_incoming_FE_QQ_from_same_b_total_crossx << std::endl;
    *m_crossx_summary_file << "Total crossx [same GS]: " << from_same_gluon_spitting_total_crossx << std::endl;
    *m_crossx_summary_file << "Total crossx [QQ']: " << hard_QQ_scatt_total_crossx << std::endl;
    *m_crossx_summary_file << "Total crossx [FE]: " << FE_total_crossx << std::endl;
    *m_crossx_summary_file << "  [FE - same b]: " << FE_from_same_b_total_crossx << std::endl;
    *m_crossx_summary_file << "  [FE - same GS]: " << FE_from_same_GS_total_crossx << std::endl;
    *m_crossx_summary_file << "  [FE - diff ancestors]: " << FE_from_diff_ancestors_total_crossx << std::endl;
    *m_crossx_summary_file << "  [FE - others]: " << FE_from_same_ancestors_not_same_b_or_gs_total_crossx << std::endl;
    *m_crossx_summary_file << "Total crossx [either from tau]: " << either_from_tau_total_crossx << std::endl;
    *m_crossx_summary_file << "For pairs from resonance decay:" << std::endl;
    for (const auto& res : resonance_id_to_name_and_crossx_map)
        *m_crossx_summary_file << "  crossx [" << res.second.first << " decay]: " << res.second.second << std::endl;
}

// ===========================================================================
// MuonPairTagsReinit
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::MuonPairTagsReinit() {
    m1_c_tag = false;  m2_c_tag = false;
    m1_osc = false;    m2_osc = false;
    m1_from_tau = false; m2_from_tau = false;
    skip_event_origin_analysis = false;
    from_same_gluon_photon_splitting_or_both_HQ_incoming = false;
    m1_ancestor_is_incoming = false;
    m2_ancestor_is_incoming = false;
    m1_from_hard_scatt_before_gs = false; m1_from_hard_scatt_after_gs = false;
    m2_from_hard_scatt_before_gs = false; m2_from_hard_scatt_after_gs = false;
    m1_hard_scatt_in_bar1 = -10; m2_hard_scatt_in_bar1 = -10;
    m1_eldest_bhadron_barcode = -10; m2_eldest_bhadron_barcode = -10;
    m1_hard_scatt_Q = -10.; m2_hard_scatt_Q = -10.;
    m1_resonance_barcode = -10; m2_resonance_barcode = -10;
    m1_youngest_non_chadron_parent_barcode = -10;
    m2_youngest_non_chadron_parent_barcode = -10;

    PairT* mpair = self().mpairRef().get();
    if (!mpair) return;

    mpair->m1_parent_group    = -10;
    mpair->m2_parent_group    = -10;
    mpair->Qsplit             = -10.;
    mpair->mHard_relevant     = -10.;
    mpair->m1_hard_scatt_category = -10;
    mpair->m2_hard_scatt_category = -10;
    mpair->muon_pair_origin_category  = -10;
    mpair->muon_pair_flavor_category  = other_flavor;
    mpair->from_same_b        = false;
    mpair->from_same_ancestors = false;
    mpair->m1_from_pdf        = false;
    mpair->m2_from_pdf        = false;

    if (m1_history) { for (auto& v : *m1_history) v.clear(); m1_history->clear(); }
    if (m2_history) { for (auto& v : *m2_history) v.clear(); m2_history->clear(); }
    if (m1_history_particle) { for (auto& v : *m1_history_particle) v.clear(); m1_history_particle->clear(); }
    if (m2_history_particle) { for (auto& v : *m2_history_particle) v.clear(); m2_history_particle->clear(); }

    if (m1_first_hf_hadron_prt_pt_eta_phi_m) m1_first_hf_hadron_prt_pt_eta_phi_m->clear();
    if (m2_first_hf_hadron_prt_pt_eta_phi_m) m2_first_hf_hadron_prt_pt_eta_phi_m->clear();

    m1_parton_ancestor_ids.clear();  m2_parton_ancestor_ids.clear();
    m1_parton_ancestor_bars.clear(); m2_parton_ancestor_bars.clear();
    m1_multi_hf_quark_ids.clear();   m2_multi_hf_quark_ids.clear();
}

// ===========================================================================
// PerformTruthPairAnalysis
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::PerformTruthPairAnalysis() {
    MuonPairAncestorTracing();
}

// ===========================================================================
// PrintHistory
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::PrintHistory(std::ostream* f, bool print_single,
                                                      bool muon1_sameancestor, bool print_category) {
    if (!f) return;
    PairT* mpair = self().mpairRef().get();
    if (!mpair) return;

    *f << "Event #: " << mpair->m1.ev_num << std::endl;

    if (print_single) {
        if (print_category) {
            int hsc = muon1_sameancestor ? mpair->m1_hard_scatt_category : mpair->m2_hard_scatt_category;
            std::string hname = (hsc > 0 && hsc < nHardScattCategories) ? HardScattCategoryLabels[hsc] : "no hard scattering";
            *f << "hard scattering category: " << hname << std::endl;
            std::string pname = (mpair->muon_pair_origin_category > 0 && mpair->muon_pair_origin_category < nHFMuonPairOriginCategories)
                                ? HFMuonPairOriginCategoryLabels[mpair->muon_pair_origin_category] : "undetermined";
            *f << "pair origin: " << pname << std::endl;
        }
        auto* hp = muon1_sameancestor ? m1_history_particle : m2_history_particle;
        if (!hp) return;
        *f << "id (barcode, status)" << std::endl;
        for (auto& step : *hp) {
            for (auto& p : step) *f << p.id << " (" << p.barcode << ", " << p.status << ") ";
            *f << "<--- ";
        }
        *f << std::endl << std::endl;
    } else {
        if (print_category) {
            *f << (mpair->from_same_ancestors ? "same ancestors" : "different ancestors") << std::endl;
            std::string pname = (mpair->muon_pair_origin_category > 0 && mpair->muon_pair_origin_category < nHFMuonPairOriginCategories)
                                ? HFMuonPairOriginCategoryLabels[mpair->muon_pair_origin_category] : "undetermined";
            *f << "pair origin: " << pname << std::endl;
        }
        *f << "id (barcode, status)" << std::endl;
        if (m1_history_particle) {
            for (auto& step : *m1_history_particle) {
                for (auto& p : step) *f << p.id << " (" << p.barcode << ", " << p.status << ") ";
                *f << "<--- ";
            }
            *f << std::endl << std::endl;
        }
        if (m2_history_particle) {
            for (auto& step : *m2_history_particle) {
                for (auto& p : step) *f << p.id << " (" << p.barcode << ", " << p.status << ") ";
                *f << "<--- ";
            }
            *f << std::endl << std::endl;
        }
    }
}

// ===========================================================================
// GetPtEtaPhiMFromBarcode
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::GetPtEtaPhiMFromBarcode(int barcode, std::vector<float>* out) {
    out->clear();
    out->push_back(static_cast<float>(TruthPtAt(barcode)));
    out->push_back(static_cast<float>(TruthEtaAt(barcode)));
    out->push_back(static_cast<float>(TruthPhiAt(barcode)));
    out->push_back(static_cast<float>(TruthMAt(barcode)));
}

// ===========================================================================
// ParentGrouping
// ===========================================================================

template <class PairT, class Derived>
int PythiaTruthExtras<PairT, Derived>::ParentGrouping(std::vector<int>& parent_ids, bool c_tag, bool prev_is_lepton) {

    if (parent_ids.size() == 2 && abs(parent_ids[0]) <= 5
        && parent_ids[1] == -parent_ids[0] && prev_is_lepton)
        return prt_drell_yan;

    int pid = abs(parent_ids[0]) % 10000;

    auto it = std::find(resonance_ids.begin(), resonance_ids.end(), pid);
    if (it != resonance_ids.end()) return resonance_decay;

    if ((pid >= 500 && pid < 600) || (pid >= 5000 && pid < 6000))
        return c_tag ? b_to_c : direct_b;

    if (c_tag) return direct_c;

    if ((pid >= 100 && pid < 400) || (pid >= 1000 && pid < 4000))
        return s_light;

    if (pid == 22) return single_photon;

    std::cout << "Not in any parent group. IDs: ";
    for (auto id : parent_ids) std::cout << id << " ";
    std::cout << std::endl;
    return -1;
}

// ===========================================================================
// UpdateCurParents
// ===========================================================================

template <class PairT, class Derived>
std::pair<int,int> PythiaTruthExtras<PairT, Derived>::UpdateCurParents(
    bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids,
    bool before_gs, bool& prev_out_hard_scatt, int hf_quark_index)
{
    PairT* mpair = self().mpairRef().get();

    bool cur_is_quark_gluon = (abs(cur_prt_ids[0]) <= 5 || cur_prt_ids[0] == 21);
    if (cur_prt_bars.size() > 1 && !cur_is_quark_gluon) {
        if (print_unspecified_parent && m_unspecified_parent_file) {
            *m_unspecified_parent_file << "Event#: " << (mpair ? mpair->m1.ev_num : -1) << std::endl;
            *m_unspecified_parent_file << "More than one parent at hadronic level:" << std::endl;
            for (int id : cur_prt_ids) *m_unspecified_parent_file << id << " ";
            *m_unspecified_parent_file << std::endl << std::endl;
        }
    }

    int prev_first_prt_id  = cur_prt_ids[0];
    int prev_first_prt_bar = cur_prt_bars[0];
    int index = (hf_quark_index >= 0) ? hf_quark_index : 0;

    int mother1 = 0, mother2 = 0;
    const int cur_idx = GetParticleIndex(cur_prt_bars[index]);

    cur_prt_bars.clear();
    cur_prt_ids.clear();

    if (isPrivate) {
        mother1 = truth_mother1->at(cur_idx);
        mother2 = truth_mother2->at(cur_idx);

        int status = abs(truth_status->at(cur_idx));

        if (mother1 == 0) {
            std::cout << "Error:: Mother1 cannot be 0" << std::endl;
            PrintHistory(&std::cout, true, isMuon1);
            std::cout << "bar: " << prev_first_prt_bar << ", mother1 " << mother1 << ", mother2 " << mother2 << std::endl;
            return {0, 0};
        }

        if (mother2 == mother1 || mother2 == 0) {
            cur_prt_bars.push_back(mother1);
            cur_prt_ids.push_back(truth_id->at(GetParticleIndex(mother1)) % 10000);
            if ((mother1 == 1 || mother1 == 2) && before_gs) {
                if (isMuon1) m1_ancestor_is_incoming = true;
                else         m2_ancestor_is_incoming = true;
            }
        } else if (((status >= 81 && status <= 86) || (status >= 101 && status <= 106)) && mother1 < mother2) {
            for (int cm = mother1; cm <= mother2; cm++) {
                cur_prt_bars.push_back(cm);
                cur_prt_ids.push_back(truth_id->at(GetParticleIndex(cm)) % 10000);
            }
        } else {
            cur_prt_bars.push_back(mother1);
            cur_prt_ids.push_back(truth_id->at(GetParticleIndex(mother1)) % 10000);
            cur_prt_bars.push_back(mother2);
            cur_prt_ids.push_back(truth_id->at(GetParticleIndex(mother2)) % 10000);
        }
    } else {
        if (!truth_parents || cur_idx >= (int)truth_parents->size()) {
            std::cout << "Error:: parent list unavailable" << std::endl;
            PrintHistory(&std::cout, true, isMuon1);
            std::cout << "bar: " << prev_first_prt_bar << std::endl;
            return {0, 0};
        }

        const auto& parents = truth_parents->at(cur_idx);
        if (parents.empty()) {
            std::cout << "Error:: parent list empty" << std::endl;
            PrintHistory(&std::cout, true, isMuon1);
            std::cout << "bar: " << prev_first_prt_bar << std::endl;
            return {0, 0};
        }

        for (int parent_barcode : parents) {
            int parent_idx = -1;
            try {
                parent_idx = GetParticleIndex(parent_barcode);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "UpdateCurParents(non-private): failed to find parent barcode "
                    + std::to_string(parent_barcode)
                    + " for child barcode " + std::to_string(prev_first_prt_bar)
                    + ". " + e.what()
                );
            }

            cur_prt_bars.push_back(parent_barcode);
            cur_prt_ids.push_back(truth_id->at(parent_idx) % 10000);
            if ((parent_barcode == 1 || parent_barcode == 2) && before_gs) {
                if (isMuon1) m1_ancestor_is_incoming = true;
                else         m2_ancestor_is_incoming = true;
            }
        }
    }

    if ((cur_prt_bars.empty() || cur_prt_ids.empty()) && (print_bad_warnings && m_very_bad_warning_file)) {
        *m_very_bad_warning_file << "Error:: parent list empty." << std::endl;
        if (mpair) *m_very_bad_warning_file << "Event: " << mpair->m1.ev_num << std::endl;
    }

    // Build profile for history
    std::vector<Particle> cur_prt_profiles;
    for (int bar : cur_prt_bars) {
        float pt  = static_cast<float>(TruthPtAt(bar));
        float eta = static_cast<float>(TruthEtaAt(bar));
        float phi = static_cast<float>(TruthPhiAt(bar));
        int id_p  = truth_id->at(GetParticleIndex(bar));
        int st_p  = truth_status->at(GetParticleIndex(bar));
        cur_prt_profiles.push_back({pt, eta, phi, bar, id_p, st_p});
    }

    if (isMuon1) {
        m1_history->push_back(cur_prt_ids);
        m1_history_particle->push_back(cur_prt_profiles);
    } else {
        m2_history->push_back(cur_prt_ids);
        m2_history_particle->push_back(cur_prt_profiles);
    }

    // Handle "prev_out_hard_scatt" case — we just traced back from an outgoing hard-scatter particle
    if (prev_out_hard_scatt) {
        const int cur0_idx = cur_prt_bars.empty() ? -1 : GetParticleIndex(cur_prt_bars[0]);
        if ((cur_prt_bars.size() != 2 || (abs(truth_status->at(cur0_idx)) != 21 && abs(truth_status->at(cur0_idx)) != 31))
            && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)) {
            *m_very_bad_warning_file << "Current parents should be incoming particles of a hard scattering." << std::endl;
            PrintHistory(m_very_bad_warning_file, true, isMuon1);
        }
        int catgr = (cur_prt_bars.size() >= 2) ? HardAnalysisCategr(cur_prt_bars[0], cur_prt_bars[1]) : -1;

        bool bars_ok = false;
        if (cur_prt_bars.size() >= 2) {
            try {
                (void)GetParticleIndex(cur_prt_bars[0] + 2);
                (void)GetParticleIndex(cur_prt_bars[1] + 2);
                bars_ok = true;
            } catch (...) {
                bars_ok = false;
            }
        }
        if (isMuon1) {
            if (bars_ok) {
                vm1out1.SetPtEtaPhiM(TruthPtAt(cur_prt_bars[0]+2), TruthEtaAt(cur_prt_bars[0]+2),
                                      TruthPhiAt(cur_prt_bars[0]+2), TruthMAt(cur_prt_bars[0]+2));
                vm1out2.SetPtEtaPhiM(TruthPtAt(cur_prt_bars[1]+2), TruthEtaAt(cur_prt_bars[1]+2),
                                      TruthPhiAt(cur_prt_bars[1]+2), TruthMAt(cur_prt_bars[1]+2));
            }
        } else {
            if (bars_ok) {
                vm2out1.SetPtEtaPhiM(TruthPtAt(cur_prt_bars[0]+2), TruthEtaAt(cur_prt_bars[0]+2),
                                      TruthPhiAt(cur_prt_bars[0]+2), TruthMAt(cur_prt_bars[0]+2));
                vm2out2.SetPtEtaPhiM(TruthPtAt(cur_prt_bars[1]+2), TruthEtaAt(cur_prt_bars[1]+2),
                                      TruthPhiAt(cur_prt_bars[1]+2), TruthMAt(cur_prt_bars[1]+2));
            }
        }

        if (isMuon1) {
            if (before_gs) m1_from_hard_scatt_before_gs = true;
            else           m1_from_hard_scatt_after_gs  = true;
            m1_hard_scatt_in_bar1 = cur_prt_bars[0];
            if (mpair) mpair->m1_hard_scatt_category = catgr;
            TLorentzVector vhard = vm1out1 + vm1out2;
            m1_hard_scatt_Q = vhard.M();
        } else {
            if (before_gs) m2_from_hard_scatt_before_gs = true;
            else           m2_from_hard_scatt_after_gs  = true;
            m2_hard_scatt_in_bar1 = cur_prt_bars[0];
            if (mpair) mpair->m2_hard_scatt_category = catgr;
            TLorentzVector vhard = vm2out1 + vm2out2;
            m2_hard_scatt_Q = vhard.M();
        }
        prev_out_hard_scatt = false;
    }

    // Check if current parent is outgoing from hard scatter (status 23 or 33)
    if (!cur_prt_bars.empty() &&
        (abs(truth_status->at(GetParticleIndex(cur_prt_bars[0]))) == 23 || abs(truth_status->at(GetParticleIndex(cur_prt_bars[0]))) == 33)) {
        if ((isMuon1 && m1_hard_scatt_in_bar1 > 0) || (!isMuon1 && m2_hard_scatt_in_bar1 > 0)) {
            if (print_bad_warnings && m_very_bad_warning_file)
                *m_very_bad_warning_file << "More than one hard scattering in muon history chain." << std::endl;
        }
        prev_out_hard_scatt = true;
    }

    // Record oscillation
    if (!cur_prt_ids.empty() && cur_prt_ids[0] == -prev_first_prt_id
        && abs(prev_first_prt_id) != 4 && abs(prev_first_prt_id) != 5) {
        if (isMuon1) m1_osc = true;
        else         m2_osc = true;
    }

    // Check for resonance
    if (!cur_prt_ids.empty()) {
        auto it_res = std::find(resonance_ids.begin(), resonance_ids.end(), abs(cur_prt_ids[0]) % 10000);
        if (it_res != resonance_ids.end()) {
            if (isMuon1) {
                m1_resonance_barcode = cur_prt_bars[0];
                if (mpair) mpair->m1_parent_group = static_cast<int>(single_muon_parent_group::resonance_decay);
            } else {
                m2_resonance_barcode = cur_prt_bars[0];
                if (mpair) mpair->m2_parent_group = static_cast<int>(single_muon_parent_group::resonance_decay);
            }
            return {0, 0};
        }
    }

    // c-tag
    if (!cur_prt_ids.empty() && ((abs(cur_prt_ids[0]) >= 400 && abs(cur_prt_ids[0]) < 500) ||
                                  (abs(cur_prt_ids[0]) >= 4000 && abs(cur_prt_ids[0]) < 5000))) {
        if (isMuon1) m1_c_tag = true;
        else         m2_c_tag = true;
    }

    // muon to hadron stage: record last HF hadron
    if (mpair && (abs(prev_first_prt_id) == 13 || abs(prev_first_prt_id) == 15) && !cur_prt_bars.empty()) {
        int cpid = abs(cur_prt_ids[0]);
        if ((cpid >= 400 && cpid < 600) || (cpid >= 4000 && cpid < 6000)) {
            if (isMuon1) GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m1_last_hf_hadron_prt_pt_eta_phi_m);
            else         GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m2_last_hf_hadron_prt_pt_eta_phi_m);
        }
    }

    // last b-flavored hadron
    if (!cur_prt_bars.empty()) {
        bool prev_is_b = ((abs(prev_first_prt_id) >= 500 && abs(prev_first_prt_id) < 600) ||
                          (abs(prev_first_prt_id) >= 5000 && abs(prev_first_prt_id) < 6000));
        bool cur_is_b  = ((abs(cur_prt_ids[0]) >= 500 && abs(cur_prt_ids[0]) < 600) ||
                          (abs(cur_prt_ids[0]) >= 5000 && abs(cur_prt_ids[0]) < 6000));
        if (!prev_is_b && cur_is_b && mpair) {
            if (isMuon1) GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m1_last_b_hadron_prt_pt_eta_phi_m);
            else         GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m2_last_b_hadron_prt_pt_eta_phi_m);
        }
    }

    // from hadron to quark/gluon stage: record first HF hadron
    if (!cur_prt_bars.empty()) {
        bool prev_is_hf = ((abs(prev_first_prt_id) >= 400 && abs(prev_first_prt_id) < 600) ||
                           (abs(prev_first_prt_id) >= 4000 && abs(prev_first_prt_id) < 6000));
        if (prev_is_hf && ((abs(cur_prt_ids[0]) < 4000 && abs(cur_prt_ids[0]) % 100 <= 5) || cur_prt_ids[0] == 21)) {
            if (isMuon1) GetPtEtaPhiMFromBarcode(prev_first_prt_bar, m1_first_hf_hadron_prt_pt_eta_phi_m);
            else         GetPtEtaPhiMFromBarcode(prev_first_prt_bar, m2_first_hf_hadron_prt_pt_eta_phi_m);
            return {prev_first_prt_bar, prev_first_prt_id};
        }
    }

    return {0, 0};
}

// ===========================================================================
// FindHeavyQuarks
// ===========================================================================

template <class PairT, class Derived>
int PythiaTruthExtras<PairT, Derived>::FindHeavyQuarks(
    std::vector<int>& cur_prt_ids, std::vector<int>& cur_prt_bars,
    int quark_type, bool isMuon1, int prev_hq_bar, int hadron_child_id)
{
    if ((quark_type != 4 && quark_type != 5) && print_bad_warnings && m_very_bad_warning_file) {
        *m_very_bad_warning_file << "Error: quark_type must be 4 or 5." << std::endl;
    }
    if (cur_prt_ids.empty() && print_bad_warnings && m_very_bad_warning_file) {
        *m_very_bad_warning_file << "Error: parent list empty in FindHeavyQuarks." << std::endl;
        return -1;
    }

    int status = truth_status->at(GetParticleIndex(cur_prt_bars[0]));

    if (cur_prt_ids.size() == 1 && (cur_prt_ids[0] == 2212 || cur_prt_ids[0] == 2112))
        return -2;

    int index_candidate = -1;
    auto it_q    = std::find(cur_prt_ids.begin(), cur_prt_ids.end(),  quark_type);
    auto it_qbar = std::find(cur_prt_ids.begin(), cur_prt_ids.end(), -quark_type);

    if (it_q != cur_prt_ids.end()) {
        index_candidate = it_q - cur_prt_ids.begin();
        if (it_qbar != cur_prt_ids.end()) {
            int sign_correction = (abs(hadron_child_id) >= 500 && abs(hadron_child_id) < 600) ? -1 : +1;
            if (hadron_child_id * sign_correction > 0)
                index_candidate = it_q - cur_prt_ids.begin();
            else if (hadron_child_id * sign_correction < 0)
                index_candidate = it_qbar - cur_prt_ids.begin();
            else if (abs(status) == 21 || abs(status) == 31) {
                index_candidate = (prev_hq_bar - 2 == cur_prt_bars[0]) ? 0 : 1;
                if (truth_id->at(GetParticleIndex(prev_hq_bar)) != cur_prt_ids[index_candidate]) {
                    if (self().debug_mode){                        
                        std::cout << "HQ mismatch in FindHeavyQuarks." << std::endl;
                        PrintHistory(&std::cout, true, isMuon1);
                    }
                }
            } else if (print_bad_warnings && m_very_bad_warning_file) {
                *m_very_bad_warning_file << "Both Q and Qbar in parent vec, not incoming to HS, not immediate parent of hadron." << std::endl;
                PrintHistory(m_very_bad_warning_file, true, isMuon1);
            }
        }

        // Two quarks of same type
        if (std::find(it_q+1, cur_prt_ids.end(), quark_type) != cur_prt_ids.end()) {
            if (cur_prt_ids.size() != 2 || (abs(status) != 21 && abs(status) != 31) || it_qbar != cur_prt_ids.end()) {
                skip_event_origin_analysis = true;
                return index_candidate;
            }
            index_candidate = (prev_hq_bar - 2 == cur_prt_bars[0]) ? 0 : 1;
        }
    } else if (it_qbar != cur_prt_ids.end()) {
        index_candidate = it_qbar - cur_prt_ids.begin();
        auto it_qbar2 = std::find(it_qbar+1, cur_prt_ids.end(), -quark_type);
        if (it_qbar2 != cur_prt_ids.end()) {
            if (cur_prt_ids.size() != 2 || (abs(status) != 21 && abs(status) != 31)) {
                std::cout << "More than one same-sign Qbar in parent vec, not from HS." << std::endl;
                PrintHistory(&std::cout, true, isMuon1);
                return index_candidate;
            }
            index_candidate = (prev_hq_bar - 2 == cur_prt_bars[0]) ? 0 : 1;
        }
    }
    return index_candidate;
}

// ===========================================================================
// HardAnalysisCategr
// ===========================================================================

template <class PairT, class Derived>
int PythiaTruthExtras<PairT, Derived>::HardAnalysisCategr(int in_bar1, int in_bar2) {
    PairT* mpair = self().mpairRef().get();

    int s1 = abs(truth_status->at(GetParticleIndex(in_bar1)));
    int s2 = abs(truth_status->at(GetParticleIndex(in_bar2)));
    int id1 = abs(truth_id->at(GetParticleIndex(in_bar1)));
    int id2 = abs(truth_id->at(GetParticleIndex(in_bar2)));

    if ((in_bar2 != in_bar1 + 1 || (s1 != 21 && s1 != 31) || (s2 != 21 && s2 != 31))
        && print_bad_warnings && m_hard_scattering_warning_file)
        *m_hard_scattering_warning_file << "Invalid hard-scatter barcodes!" << std::endl;

    bool in1_hq = (id1 == 4 || id1 == 5);
    bool in2_hq = (id2 == 4 || id2 == 5);
    bool in1_glq = (id1 <= 3 || id1 == 21);
    bool in2_glq = (id2 <= 3 || id2 == 21);

    int out_bar1 = in_bar1 + 2;
    int out_bar2 = in_bar2 + 2;
    int out_id1 = abs(truth_id->at(GetParticleIndex(out_bar1)));
    int out_id2 = abs(truth_id->at(GetParticleIndex(out_bar2)));

    bool fe = (in1_hq && in2_glq) || (in2_hq && in1_glq);
    if (fe) return flavor_excit;

    if (in1_hq && in2_hq && truth_id->at(GetParticleIndex(in_bar1)) == truth_id->at(GetParticleIndex(out_bar1))
        && truth_id->at(GetParticleIndex(in_bar2)) == truth_id->at(GetParticleIndex(out_bar2)))
        return hq_scatt;

    bool final_QQbar = ((out_id1 == 4 || out_id1 == 5)
                        && truth_id->at(GetParticleIndex(out_bar2)) == -truth_id->at(GetParticleIndex(out_bar1)));
    if (final_QQbar) {
        if (!((id1 <= 5 || id1 == 21) && id2 == id1) && print_bad_warnings && m_hard_scattering_warning_file)
            *m_hard_scattering_warning_file << "FC with unexpected initial state." << std::endl;
        return flavor_creat;
    }

    if (id1 == 21 && id2 == 21) {
        if (truth_id->at(GetParticleIndex(out_bar1)) == 21) return gg_gg;
        if (out_id1 > 3 && print_bad_warnings && m_hard_scattering_warning_file)
            *m_hard_scattering_warning_file << "gg->qqbar unexpected." << std::endl;
        return gg_qqbar;
    }
    if (out_id1 == 21 && out_id2 == 21) {
        if (id1 <= 3) return qqbar_gg;
    }
    if ((id1 == 21 && id2 <= 3) || (id2 == 21 && id1 <= 3)) return gq_gq;

    if ((out_id1 == 13 || out_id1 == 15) && (out_id2 == 13 || out_id2 == 15)) {
        if (id1 != id2 && print_bad_warnings && m_hard_scattering_warning_file)
            *m_hard_scattering_warning_file << "DY with unexpected incoming." << std::endl;
        return hard_scatt_drell_yan;
    }

    if ((id1 > 3 || id2 > 3 || out_id1 > 3 || out_id2 > 3) && print_bad_warnings && m_hard_scattering_warning_file)
        *m_hard_scattering_warning_file << "qq'->qq' with heavy quarks." << std::endl;
    return qqprime_qqprime;
}

// ===========================================================================
// GluonHistoryTracking
// ===========================================================================

template <class PairT, class Derived>
int PythiaTruthExtras<PairT, Derived>::GluonHistoryTracking(int gluon_bar, bool isMuon1) {
    std::vector<int> parent_bars = {gluon_bar};
    std::vector<int> parent_ids  = {21};
    bool prev_step_m23_m33 = (abs(truth_status->at(GetParticleIndex(gluon_bar))) == 23
                              || abs(truth_status->at(GetParticleIndex(gluon_bar))) == 33);

    float pt  = static_cast<float>(TruthPtAt(gluon_bar));
    float eta = static_cast<float>(TruthEtaAt(gluon_bar));
    float phi = static_cast<float>(TruthPhiAt(gluon_bar));
    Particle p {pt, eta, phi, gluon_bar, 21, truth_status->at(GetParticleIndex(gluon_bar))};
    // Note: we track in the muon's history vectors - UpdateCurParents appends
    // We do NOT push gluon's initial entry here; GluonHistoryTracking uses its own local vectors
    // and calls UpdateCurParents which appends into m1/m2_history_* for the correct muon.

    bool terminate = false;
    while (!terminate) {
        UpdateCurParents(isMuon1, parent_bars, parent_ids, false, prev_step_m23_m33);
        int pm1;
        if (isPrivate) {
            pm1 = truth_mother1->at(GetParticleIndex(parent_bars[0]));
        } else {
            int pidx = GetParticleIndex(parent_bars[0]);
            pm1 = (truth_parents && pidx < (int)truth_parents->size() && !truth_parents->at(pidx).empty())
                  ? truth_parents->at(pidx).front() : 0;
        }
        terminate = (abs(truth_status->at(GetParticleIndex(parent_bars[0]))) == 21 || pm1 <= 2);
    }

    if (abs(truth_status->at(GetParticleIndex(parent_bars[0]))) == 21) return 1; // FSR
    return 2; // ISR
}

// ===========================================================================
// FillCategoryHistograms
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::FillCategoryHistograms(TH1D* horig, TH2D* hhard_vs_orig) {
    PairT* mpair = self().mpairRef().get();
    if (!mpair || !horig || !hhard_vs_orig) return;

    horig->Fill(mpair->muon_pair_origin_category, mpair->weight);

    if (m1_hard_scatt_in_bar1 > 0) {
        if (m2_hard_scatt_in_bar1 < 0 || m1_hard_scatt_in_bar1 == m2_hard_scatt_in_bar1)
            hhard_vs_orig->Fill(mpair->muon_pair_origin_category, mpair->m1_hard_scatt_category, mpair->weight);
    } else if (m2_hard_scatt_in_bar1 > 0) {
        if (m1_hard_scatt_in_bar1 > 0 && print_bad_warnings && m_very_bad_warning_file) {
            *m_very_bad_warning_file << "Two cases should be exclusive in FillCategoryHistograms." << std::endl;
        }
        hhard_vs_orig->Fill(mpair->muon_pair_origin_category, mpair->m2_hard_scatt_category, mpair->weight);
    }
}

// ===========================================================================
// SingleMuonAncestorTracing
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::SingleMuonAncestorTracing(bool isMuon1) {
    PairT* mpair = self().mpairRef().get();
    if (!mpair) return;

    std::vector<int> parent_bars;
    std::vector<int> parent_ids;
    int first_hadron_id = 0;
    int first_hadron_barcode = 0;
    int prev_first_prt_id = -1;
    std::pair<int,int> first_hadron_bar_id;
    bool prev_step_status_m23_m33 = false;

    // Add muon itself into history
    float pt  = isMuon1 ? static_cast<float>(mpair->m1.truth_pt)  : static_cast<float>(mpair->m2.truth_pt);
    float eta = isMuon1 ? static_cast<float>(mpair->m1.truth_eta) : static_cast<float>(mpair->m2.truth_eta);
    float phi = isMuon1 ? static_cast<float>(mpair->m1.truth_phi) : static_cast<float>(mpair->m2.truth_phi);
    int barcode = isMuon1 ? mpair->m1.truth_bar : mpair->m2.truth_bar;
    int mu_idx = -1;
    try {
        mu_idx = GetParticleIndex(barcode);
    } catch (const std::exception& e) {
        std::cerr << "[SAT] Warning: muon barcode=" << barcode
                  << " cannot be resolved (" << e.what() << "). Skipping ancestor tracing for this muon." << std::endl;
        return;
    }
    int charge  = isMuon1 ? mpair->m1.truth_charge : mpair->m2.truth_charge;
    int id = -13 * charge;
    int status = truth_status->at(mu_idx);
    Particle p {pt, eta, phi, barcode, id, status};

    parent_bars.push_back(barcode);
    parent_ids.push_back(id);
    std::vector<Particle> cur_muon_profile = {p};

    if (isMuon1) { m1_history->push_back(parent_ids); m1_history_particle->push_back(cur_muon_profile); }
    else         { m2_history->push_back(parent_ids); m2_history_particle->push_back(cur_muon_profile); }

    // Trace back through leptons + c-hadrons
    // Guard against self-loop entries in truth_parents (non-private samples can have particle = its own parent)
    std::vector<int> prev_bars_snapshot;
    int trace_iter = 0;
    const int kMaxTraceIter = 500;
    while (abs(parent_ids[0]) == 13 || abs(parent_ids[0]) == 15
           || (abs(parent_ids[0]) >= 400 && abs(parent_ids[0]) < 500)
           || (abs(parent_ids[0]) >= 4000 && abs(parent_ids[0]) < 5000)) {
        if (++trace_iter > kMaxTraceIter) {
            std::cerr << "[SAT] Warning: ancestor tracing exceeded " << kMaxTraceIter
                      << " iterations (barcode=" << barcode << "). Breaking." << std::endl;
            break;
        }
        // Cycle detection: if parent set unchanged, stop
        if (parent_bars == prev_bars_snapshot) {
            std::cerr << "[SAT] Warning: detected self-loop in truth_parents at barcode="
                      << (!parent_bars.empty() ? parent_bars[0] : -1) << ". Breaking." << std::endl;
            break;
        }
        prev_bars_snapshot = parent_bars;
        if (abs(parent_ids[0]) == 15) {
            if (isMuon1) m1_from_tau = true;
            else         m2_from_tau = true;
        }
        prev_first_prt_id = parent_ids[0];
        first_hadron_bar_id = UpdateCurParents(isMuon1, parent_bars, parent_ids, true, prev_step_status_m23_m33);
        if ((isMuon1 && mpair->m1_parent_group == resonance_decay) ||
            (!isMuon1 && mpair->m2_parent_group == resonance_decay)) return;
        first_hadron_barcode = first_hadron_bar_id.first;
        first_hadron_id = first_hadron_bar_id.second;
    }

    if ((isMuon1 && mpair->m1_parent_group == resonance_decay) ||
        (!isMuon1 && mpair->m2_parent_group == resonance_decay)) return;

    bool prev_is_lepton = (abs(prev_first_prt_id) == 13 || abs(prev_first_prt_id) == 15);

    int cur_parent_group;
    if (isMuon1) {
        m1_earliest_parent_id = parent_ids[0];
        m1_youngest_non_chadron_parent_barcode = parent_bars[0];
        cur_parent_group = ParentGrouping(parent_ids, m1_c_tag, prev_is_lepton);
        mpair->m1_parent_group = cur_parent_group;
    } else {
        m2_earliest_parent_id = parent_ids[0];
        m2_youngest_non_chadron_parent_barcode = parent_bars[0];
        cur_parent_group = ParentGrouping(parent_ids, m2_c_tag, prev_is_lepton);
        mpair->m2_parent_group = cur_parent_group;
    }

    // Return if not from open HF
    if (cur_parent_group != direct_b &&
        cur_parent_group != b_to_c &&
        cur_parent_group != direct_c) return;

    // Trace through b-flavored hadrons
    if ((abs(parent_ids[0]) >= 500 && abs(parent_ids[0]) < 600) ||
        (abs(parent_ids[0]) >= 5000 && abs(parent_ids[0]) < 6000)) {
        m_n_bhadron_muons++;
        int last_bhadron_bc = parent_bars[0];
        while ((abs(parent_ids[0]) >= 500 && abs(parent_ids[0]) < 600) ||
               (abs(parent_ids[0]) >= 5000 && abs(parent_ids[0]) < 6000)) {
            last_bhadron_bc = parent_bars[0];
            first_hadron_bar_id = UpdateCurParents(isMuon1, parent_bars, parent_ids, true, prev_step_status_m23_m33);
            first_hadron_barcode = first_hadron_bar_id.first;
            first_hadron_id = first_hadron_bar_id.second;
        }
        // In overlay, some B-meson chains (Geant4 secondaries) end at a beam/proton
        // particle instead of a b-quark.  Store the B-hadron's own barcode so
        // same-b matching still works; also prevents crash in PrintHistory below.
        int store_bc = (abs(first_hadron_id) == 5) ? first_hadron_barcode : last_bhadron_bc;
        if (isMuon1) m1_eldest_bhadron_barcode = store_bc;
        else         m2_eldest_bhadron_barcode = store_bc;
    }

    // Heavy-quark level tracing
    int quark = (cur_parent_group == direct_c) ? 4 : 5;
    int quark_index = FindHeavyQuarks(parent_ids, parent_bars, quark, isMuon1, parent_bars[0], first_hadron_id);

    int prev_hq_bar = -1;
    while (quark_index >= 0) {
        prev_hq_bar = parent_bars[quark_index];
        UpdateCurParents(isMuon1, parent_bars, parent_ids, true, prev_step_status_m23_m33, quark_index);
        quark_index = FindHeavyQuarks(parent_ids, parent_bars, quark, isMuon1, prev_hq_bar);
    }

    if (prev_hq_bar < 0) {
        if (quark == 5) {
            // Truncated B-hadron chain: b-quark not found above the eldest B-meson.
            // In overlay, Geant4 secondary B-mesons end at a beam proton (bc=1) instead
            // of a b-quark, so FindHeavyQuarks returns -1.
            m_n_truncated_bhadron_muons++;
            int trunc_bc = isMuon1 ? m1_eldest_bhadron_barcode : m2_eldest_bhadron_barcode;
            if (m_truncated_bc_min < 0 || trunc_bc < m_truncated_bc_min) m_truncated_bc_min = trunc_bc;
            if (trunc_bc > m_truncated_bc_max)                            m_truncated_bc_max = trunc_bc;

            // Only call PrintHistory for standard MC chains where first_hadron_id is the
            // b-quark itself (abs==5); suppressed for truncated overlay chains to avoid crash.
            if (abs(first_hadron_id) == 5) {
                std::cout << "Previous HQ barcode is -1: heavy quark never found." << std::endl;
                PrintHistory(&std::cout, true, isMuon1);
            }
        }
        skip_event_origin_analysis = true;
        return;
    }

    if (isMuon1) {
        m1_parton_ancestor_ids  = parent_ids;
        m1_parton_ancestor_bars = parent_bars;
        GetPtEtaPhiMFromBarcode(prev_hq_bar, &mpair->m1_first_hq_ancestor_pt_eta_phi_m);
    } else {
        m2_parton_ancestor_ids  = parent_ids;
        m2_parton_ancestor_bars = parent_bars;
        GetPtEtaPhiMFromBarcode(prev_hq_bar, &mpair->m2_first_hq_ancestor_pt_eta_phi_m);
    }
}

// ===========================================================================
// HFMuonPairAnalysis
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::HFMuonPairAnalysis() {
    PairT* mpair = self().mpairRef().get();
    if (!mpair) return;

    // Check from_same_b
    if (m1_youngest_non_chadron_parent_barcode == m2_youngest_non_chadron_parent_barcode
        && m1_youngest_non_chadron_parent_barcode >= 0) {
        int prt_id = abs(truth_id->at(GetParticleIndex(m1_youngest_non_chadron_parent_barcode)));
        if ((prt_id >= 500 && prt_id < 600) || (prt_id >= 5000 && prt_id < 6000)) {
            if (m1_eldest_bhadron_barcode == 10 || m1_eldest_bhadron_barcode != m2_eldest_bhadron_barcode) {
                std::cout << "FROM SAME B EVENT WRONG!!! Printing history" << std::endl;
                PrintHistory(&std::cout, false, true);
            }
        }
    }

    if (m1_eldest_bhadron_barcode != -10 && m1_eldest_bhadron_barcode == m2_eldest_bhadron_barcode) {
        mpair->from_same_b = true;
        mpair->muon_pair_flavor_category = from_single_b;
        mpair->muon_pair_origin_category = not_both_from_open_and_different_HF_hadrons;
        return;
    }

    bool both_from_b = ((mpair->m1_parent_group == direct_b ||
                         mpair->m1_parent_group == b_to_c) &&
                        (mpair->m2_parent_group == direct_b ||
                         mpair->m2_parent_group == b_to_c));
    bool both_from_c = (mpair->m1_parent_group == direct_c &&
                        mpair->m2_parent_group == direct_c);
    bool one_b_one_c = (((mpair->m1_parent_group == direct_b ||
                          mpair->m1_parent_group == b_to_c) &&
                         mpair->m2_parent_group == direct_c) ||
                        ((mpair->m2_parent_group == direct_b ||
                          mpair->m2_parent_group == b_to_c) &&
                         mpair->m1_parent_group == direct_c));

    bool same_partonic_ancestors = (m1_parton_ancestor_bars == m2_parton_ancestor_bars);
    mpair->from_same_ancestors = same_partonic_ancestors;

    if (both_from_b)   mpair->muon_pair_flavor_category = bb;
    else if (both_from_c) mpair->muon_pair_flavor_category = cc;
    else if (one_b_one_c) mpair->muon_pair_flavor_category = one_b_one_c;

    double PI = self().pmsRef().PI;
    bool not_near = !(std::abs(mpair->truth_dphi) < PI / 2.);

    mpair->muon_pair_origin_category = others;

    if (skip_event_origin_analysis) {
        skipped_event_crossx += mpair->weight;
        return;
    }

    if (same_partonic_ancestors) {
        if (m1_parton_ancestor_bars.size() > 1) {
            // FC
            if (mpair->m1_hard_scatt_category != flavor_creat &&
                print_bad_warnings && m_very_bad_warning_file) {
                *m_very_bad_warning_file << "More than one ancestor but HS not FC." << std::endl;
                PrintHistory(m_very_bad_warning_file, false, same_partonic_ancestors);
            }
            mpair->muon_pair_origin_category = fc;
            mpair->mHard_relevant = m1_hard_scatt_Q;
        } else {
            from_same_gluon_photon_splitting_or_both_HQ_incoming = true;

            if (!m1_ancestor_is_incoming) {
                if (!both_from_b && !both_from_c && print_bad_warnings && m_very_bad_warning_file) {
                    *m_very_bad_warning_file << "Same GS but neither both from b nor both from c." << std::endl;
                    PrintHistory(m_very_bad_warning_file, false, true);
                }

                // Get barQ1, barQ2 for Qsplit computation
                int quark_type = both_from_b ? 5 : 4;
                int barQ1 = -1, barQ2 = -1;

                auto get_Q_bar = [&](std::vector<std::vector<Particle>>* hp, int& barQ) {
                    if (!hp || hp->size() < 2) return;
                    size_t n = hp->size();
                    if ((*hp)[n-2].size() == 1) {
                        barQ = (*hp)[n-2][0].barcode;
                    } else if (n >= 3 && (*hp)[n-3].size() == 1) {
                        barQ = (*hp)[n-3][0].barcode - 2;
                    }
                };
                get_Q_bar(m1_history_particle, barQ1);
                get_Q_bar(m2_history_particle, barQ2);

                if (barQ1 >= 0 && barQ2 >= 0) {
                    vQ1.SetPtEtaPhiM(TruthPtAt(barQ1), TruthEtaAt(barQ1), TruthPhiAt(barQ1), TruthMAt(barQ1));
                    vQ2.SetPtEtaPhiM(TruthPtAt(barQ2), TruthEtaAt(barQ2), TruthPhiAt(barQ2), TruthMAt(barQ2));
                }
                int barg = m1_parton_ancestor_bars[0];
                vg.SetPtEtaPhiM(TruthPtAt(barg), TruthEtaAt(barg), TruthPhiAt(barg), TruthMAt(barg));
            }

            // Determine GS/PhS origin
            if (m1_from_hard_scatt_before_gs || m2_from_hard_scatt_before_gs) {
                if (m1_from_hard_scatt_before_gs && m2_from_hard_scatt_before_gs) {
                    // same_gs_isr_both_hard_scatt — leave as "others"
                } else {
                    mpair->muon_pair_origin_category = same_gs_isr_one_hard_scatt;
                }
            } else if (m1_ancestor_is_incoming) {
                mpair->muon_pair_origin_category = same_gs_isr_zero_hard_scatt;
            } else {
                int gluon_mode = GluonHistoryTracking(m1_parton_ancestor_bars[0], true);
                mpair->m2_hard_scatt_category = mpair->m1_hard_scatt_category;
                if (gluon_mode == 1) {
                    if (m1_parton_ancestor_ids[0] == 21)
                        mpair->muon_pair_origin_category = same_gs_fsr;
                    else if (m1_parton_ancestor_ids[0] == 22)
                        mpair->muon_pair_origin_category = same_phs_fsr;
                } else {
                    mpair->muon_pair_origin_category = same_gs_isr_zero_hard_scatt;
                }
            }

            // Compute Qsplit
            if (!m1_ancestor_is_incoming) {
                int cat = mpair->muon_pair_origin_category;
                if (cat == same_gs_fsr || cat == same_phs_fsr) {
                    vg = vQ1 + vQ2;
                    mpair->Qsplit = vg.M();
                    mpair->mHard_relevant = m1_hard_scatt_Q;
                    if (h_Qsplit_gs_FSR[!mpair->truth_same_sign][not_near])
                        h_Qsplit_gs_FSR[!mpair->truth_same_sign][not_near]->Fill(std::abs(vg.M()), mpair->weight);
                    if (h_Qsplit_to_mHard_gs_FSR[!mpair->truth_same_sign][not_near] && mpair->mHard_relevant > 0)
                        h_Qsplit_to_mHard_gs_FSR[!mpair->truth_same_sign][not_near]->Fill(
                            std::abs(mpair->Qsplit) / mpair->mHard_relevant, mpair->weight);
                } else if (cat == same_gs_isr_one_hard_scatt) {
                    if (m1_hard_scatt_in_bar1 > 0) {
                        vQ1 = vg - vQ2;
                        mpair->Qsplit = vQ1.M();
                        mpair->mHard_relevant = m1_hard_scatt_Q;
                    } else {
                        vQ2 = vg - vQ1;
                        mpair->Qsplit = vQ2.M();
                        mpair->mHard_relevant = m2_hard_scatt_Q;
                    }
                    if (h_Qsplit_gs_ISR_one_hard_scatt[!mpair->truth_same_sign][not_near])
                        h_Qsplit_gs_ISR_one_hard_scatt[!mpair->truth_same_sign][not_near]->Fill(
                            std::abs(mpair->Qsplit), mpair->weight);
                    if (h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[!mpair->truth_same_sign][not_near] && mpair->mHard_relevant > 0)
                        h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[!mpair->truth_same_sign][not_near]->Fill(
                            std::abs(mpair->Qsplit) / mpair->mHard_relevant, mpair->weight);
                }
            }
        }
    } else {
        // Different ancestors
        if (m1_from_hard_scatt_before_gs && m1_hard_scatt_in_bar1 == m2_hard_scatt_in_bar1
            && mpair->m1_hard_scatt_category != hq_scatt) {
            if (self().debug_mode) std::cout << "Not same ancestors. Same HS, not Q-Q scattering." << std::endl;
        }

        auto trace_if_needed = [&](bool isMuon1) {
            bool& from_hs = isMuon1 ? m1_from_hard_scatt_before_gs : m2_from_hard_scatt_before_gs;
            std::vector<int>& anc_bars = isMuon1 ? m1_parton_ancestor_bars : m2_parton_ancestor_bars;
            std::vector<int>& anc_ids  = isMuon1 ? m1_parton_ancestor_ids  : m2_parton_ancestor_ids;
            bool& anc_inc = isMuon1 ? m1_ancestor_is_incoming : m2_ancestor_is_incoming;
            if (!from_hs && !anc_bars.empty()) {
                if ((anc_bars.size() == 1 && (anc_ids[0] == 21 || anc_ids[0] == 22)) && !anc_inc) {
                    GluonHistoryTracking(anc_bars[0], isMuon1);
                }
            }
        };
        trace_if_needed(true);
        trace_if_needed(false);

        if (m1_hard_scatt_in_bar1 == m2_hard_scatt_in_bar1 && m1_hard_scatt_in_bar1 > 0) {
            mpair->muon_pair_origin_category = diff_gs_same_hard_scatt;
            mpair->mHard_relevant = m1_hard_scatt_Q;
        }
    }

    // Print history for "others" with low minv
    if (mpair->muon_pair_origin_category == others) {
        if (print_HF_pair_origin_others_history && m_HF_pair_origin_others_category_file
            && mpair->truth_minv < low_minv_threshold) {
            *m_HF_pair_origin_others_category_file << "pair minv: " << mpair->truth_minv << std::endl;
            PrintHistory(m_HF_pair_origin_others_category_file, false, mpair->from_same_ancestors, true);
        }
    }

    // Fill category histograms
    int isign  = !mpair->truth_same_sign;
    int idphi = not_near ? 1 : 0;
    if (both_from_b) {
        FillCategoryHistograms(h_both_from_b_muon_pair_origin_categr[isign][idphi],
                               h_both_from_b_hard_scatt_categr_vs_origin_categr[isign][idphi]);
        FillCategoryHistograms(h_muon_pair_origin_categr[isign][idphi],
                               h_hard_scatt_categr_vs_origin_categr[isign][idphi]);
        if (print_prt_history && mpair->m1.ev_num < 5000)
            PrintHistory(m_b_parent_file[isign][idphi], false, same_partonic_ancestors, true);
    }
    if (one_b_one_c) {
        FillCategoryHistograms(h_one_from_b_one_from_c_muon_pair_origin_categr[isign][idphi],
                               h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[isign][idphi]);
        FillCategoryHistograms(h_muon_pair_origin_categr[isign][idphi],
                               h_hard_scatt_categr_vs_origin_categr[isign][idphi]);
    }
    if (both_from_c) {
        FillCategoryHistograms(h_both_from_c_muon_pair_origin_categr[isign][idphi],
                               h_both_from_c_hard_scatt_categr_vs_origin_categr[isign][idphi]);
        FillCategoryHistograms(h_muon_pair_origin_categr[isign][idphi],
                               h_hard_scatt_categr_vs_origin_categr[isign][idphi]);
        if (print_prt_history && mpair->m1.ev_num < 50000)
            PrintHistory(m_c_parent_file[isign][idphi], false, same_partonic_ancestors, true);
    }
}

// ===========================================================================
// MuonPairAncestorTracing
// ===========================================================================

template <class PairT, class Derived>
void PythiaTruthExtras<PairT, Derived>::MuonPairAncestorTracing() {
    MuonPairTagsReinit();
    PairT* mpair = self().mpairRef().get();
    if (!mpair) return;

    double PI = self().pmsRef().PI;
    bool not_near = !(std::abs(static_cast<double>(mpair->truth_dphi)) < PI / 2.);

    SingleMuonAncestorTracing(true);
    SingleMuonAncestorTracing(false);

    const int resonance_decay_val = static_cast<int>(resonance_decay);

    mpair->from_same_resonance = (mpair->m1_parent_group == resonance_decay_val &&
                                   mpair->m2_parent_group == resonance_decay_val &&
                                   m1_resonance_barcode == m2_resonance_barcode);
    if (mpair->from_same_resonance) {
        mpair->muon_pair_flavor_category = from_resonance;
        // update resonance crossx map
        if (!resonance_ids.empty() && m1_resonance_barcode >= 0) {
            int res_id = abs(truth_id->at(GetParticleIndex(m1_resonance_barcode))) % 10000;
            auto it = resonance_id_to_name_and_crossx_map.find(res_id);
            if (it != resonance_id_to_name_and_crossx_map.end())
                it->second.second += mpair->weight;
        }
        if (mpair->truth_minv < low_minv_threshold && print_low_minv_resonances && m_very_low_minv_resonance_file)
            PrintHistory(m_very_low_minv_resonance_file, false, true);
    }

    mpair->resonance_contaminated = (!mpair->from_same_resonance &&
                                     (mpair->m1_parent_group == resonance_decay_val ||
                                      mpair->m2_parent_group == resonance_decay_val));
    if (mpair->resonance_contaminated) mpair->muon_pair_flavor_category = resonance_contaminated;

    const int single_photon_val = static_cast<int>(single_photon);
    if (mpair->m1_parent_group == single_photon_val && mpair->m2_parent_group == single_photon_val)
        mpair->muon_pair_flavor_category = photon_to_dimuon_splitting;

    const int prt_dy = static_cast<int>(prt_drell_yan);
    if (mpair->m1_parent_group == prt_dy && mpair->m2_parent_group == prt_dy)
        mpair->muon_pair_flavor_category = drell_yan;

    mpair->m1_from_pdf = m1_ancestor_is_incoming;
    mpair->m2_from_pdf = m2_ancestor_is_incoming;

    if (mpair->m1_parent_group < 0) {
        int truth_bc_size = (truth_barcode ? (int)truth_barcode->size() : -1);
        if (print_unspecified_parent && m_unspecified_parent_file) {
            *m_unspecified_parent_file
                << "[SAT-fail] Muon1 barcode=" << mpair->m1.truth_bar
                << " not in truth_barcode table (size=" << truth_bc_size << ")."
                << " parent_group stays -10. History below (empty = SAT bailed before tracing):" << std::endl;
            PrintHistory(m_unspecified_parent_file, true, true);
        }
    }
    if (mpair->m2_parent_group < 0) {
        int truth_bc_size = (truth_barcode ? (int)truth_barcode->size() : -1);
        if (print_unspecified_parent && m_unspecified_parent_file) {
            *m_unspecified_parent_file
                << "[SAT-fail] Muon2 barcode=" << mpair->m2.truth_bar
                << " not in truth_barcode table (size=" << truth_bc_size << ")."
                << " parent_group stays -10. History below (empty = SAT bailed before tracing):" << std::endl;
            PrintHistory(m_unspecified_parent_file, true, false);
        }
    }

    int isign = !mpair->truth_same_sign;
    int idphi = not_near ? 1 : 0;
    if (isign >= 0 && isign < ParamsSet::nSigns && idphi < 2 && h_parent_groups[isign][idphi])
        h_parent_groups[isign][idphi]->Fill(mpair->m1_parent_group + 0.5, mpair->m2_parent_group + 0.5, mpair->weight);

    mpair->pair_origin_analysis_skipped = skip_event_origin_analysis;

    const int direct_b_val = static_cast<int>(direct_b);
    const int b_to_c_val   = static_cast<int>(b_to_c);
    const int direct_c_val = static_cast<int>(direct_c);
    bool m1_hf = (mpair->m1_parent_group == direct_b_val || mpair->m1_parent_group == b_to_c_val || mpair->m1_parent_group == direct_c_val);
    bool m2_hf = (mpair->m2_parent_group == direct_b_val || mpair->m2_parent_group == b_to_c_val || mpair->m2_parent_group == direct_c_val);

    if (m1_hf && m2_hf) {
        HFMuonPairAnalysis();
    } else
        mpair->muon_pair_origin_category = static_cast<int>(muon_pair_both_from_open_HF_origin_catgr::not_both_from_open_and_different_HF_hadrons);

    // Print "other flavor" with low minv
    if (mpair->muon_pair_flavor_category == other_flavor) {
        if (print_other_flavor_history && m_other_flavor_category_file && mpair->truth_minv < low_minv_threshold) {
            *m_other_flavor_category_file << "others flavor with low minv" << std::endl;
            *m_other_flavor_category_file << "pair minv: " << mpair->truth_minv << std::endl;
            PrintHistory(m_other_flavor_category_file, false, mpair->from_same_ancestors);
        }
    }

    if (self().debug_mode && mpair->m1.ev_num < 50) {
        std::cout << "[DEBUG] MuonPairAncestorTracing event " << mpair->m1.ev_num
                  << ": muon1 history" << std::endl;
        PrintHistory(&std::cout, true, true, true);
        std::cout << "[DEBUG] MuonPairAncestorTracing event " << mpair->m1.ev_num
                  << ": muon2 history" << std::endl;
        PrintHistory(&std::cout, true, false, true);
    }

    if (debug_print_history_nevents > 0 && m_debug_history_file &&
        mpair->m1.ev_num < debug_print_history_nevents) {
        *m_debug_history_file << "=== Event " << mpair->m1.ev_num << " ===" << std::endl;
        PrintHistory(m_debug_history_file, true, true,  true);
        PrintHistory(m_debug_history_file, true, false, true);
    }

    if (debug_print_bhadron_id_nevents > 0 && m_debug_bhadron_id_file &&
        mpair->m1.ev_num < debug_print_bhadron_id_nevents) {
        auto safe_id = [&](int bc) -> int {
            if (bc < 0 || !truth_id) return -999;
            int idx = GetParticleIndex(bc);
            if (idx < 0 || idx >= (int)truth_id->size()) return -999;
            return truth_id->at(idx);
        };
        *m_debug_bhadron_id_file
            << "ev=" << mpair->m1.ev_num
            << "  m1_eldest_bhadron_bc=" << m1_eldest_bhadron_barcode
            << "  m1_eldest_bhadron_id=" << safe_id(m1_eldest_bhadron_barcode)
            << "  m2_eldest_bhadron_bc=" << m2_eldest_bhadron_barcode
            << "  m2_eldest_bhadron_id=" << safe_id(m2_eldest_bhadron_barcode)
            << std::endl;
    }

}
