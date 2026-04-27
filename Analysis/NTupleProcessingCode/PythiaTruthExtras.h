#pragma once

#include <map>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include "../MuonObjectsParamsAndHelpers/struct_particle.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TInterpreter.h"
#include "TLorentzVector.h"

template <class PairT, class Derived>
class PythiaTruthExtras {
    template <class, class, class, class...> friend class PythiaAlgCoreT;

public:
    using pair_t = PairT;

    // public flags / settings - must be set before Run()
    bool isPrivate = false;

    // Centre-of-mass energy in TeV. Allowed values: 5.02, 5.36 (default).
    // Used for non-private I/O path selection only; ignored for private samples.
    double E_COM = 5.36;

    bool print_prt_history = false;
    bool print_HF_pair_origin_others_history = false;
    bool print_other_flavor_history = false;
    bool print_unspecified_parent = false;
    bool print_bad_warnings = false;
    bool print_low_minv_resonances = false;
    bool print_FE = false;

protected:
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

// --------------------- truth particle branches ---------------------------

    // always int
    std::vector<int>* truth_id       = nullptr;
    std::vector<int>* truth_barcode  = nullptr;
    std::vector<int>* truth_status   = nullptr;

    // private sample: double; non-private: float
    std::vector<double>* truth_m    = nullptr;
    std::vector<double>* truth_pt   = nullptr;
    std::vector<double>* truth_eta  = nullptr;
    std::vector<double>* truth_phi  = nullptr;
    std::vector<float>*  truth_m_f  = nullptr;
    std::vector<float>*  truth_pt_f = nullptr;
    std::vector<float>*  truth_eta_f = nullptr;
    std::vector<float>*  truth_phi_f = nullptr;

    // private: truth_mother1/2; non-private: truth_parents
    std::vector<int>*                truth_mother1  = nullptr;
    std::vector<int>*                truth_mother2  = nullptr;
    std::vector<std::vector<int>>*   truth_parents  = nullptr;

// --------------------- kinematic accessors ---------------------------

    int GetParticleIndex(int barcode) const;

    double TruthPtAt(size_t i)  const;
    double TruthEtaAt(size_t i) const;
    double TruthPhiAt(size_t i) const;
    double TruthMAt(size_t i)   const;

    mutable const std::vector<int>* barcode_lookup_source = nullptr;
    mutable size_t barcode_lookup_size = 0;
    mutable int barcode_lookup_back = -1;
    mutable std::unordered_map<int, int> barcode_to_index_cache;

// --------------------- resonance data ---------------------------

    std::vector<int> resonance_ids;
    std::map<int, std::pair<std::string, double>> resonance_id_to_name_and_crossx_map;

// --------------------- tracing state ---------------------------

    float low_minv_threshold = 0.6f;

    int m1_youngest_non_chadron_parent_barcode = -10;
    int m2_youngest_non_chadron_parent_barcode = -10;
    int m1_eldest_bhadron_barcode = -10;
    int m2_eldest_bhadron_barcode = -10;

    bool m1_ancestor_is_incoming = false;
    bool m2_ancestor_is_incoming = false;
    bool m1_from_tau = false;
    bool m2_from_tau = false;
    bool m1_c_tag = false;
    bool m2_c_tag = false;
    bool m1_osc = false;
    bool m2_osc = false;

    int m1_resonance_barcode = -10;
    int m2_resonance_barcode = -10;

    bool m1_from_hard_scatt_before_gs = false;
    bool m1_from_hard_scatt_after_gs  = false;
    bool m2_from_hard_scatt_before_gs = false;
    bool m2_from_hard_scatt_after_gs  = false;
    int  m1_hard_scatt_in_bar1 = -10;
    int  m2_hard_scatt_in_bar1 = -10;
    double m1_hard_scatt_Q = -10.;
    double m2_hard_scatt_Q = -10.;

    int m1_earliest_parent_id = 0;
    int m2_earliest_parent_id = 0;

    bool from_same_gluon_photon_splitting_or_both_HQ_incoming = false;
    bool skip_event_origin_analysis = false;

    // 4-vectors for gluon splitting kinematics
    TLorentzVector vg, vQ1, vQ2, vm1out1, vm1out2, vm2out1, vm2out2;

    // History vectors (storage keeps them alive; pointers point into storage)
    std::vector<std::vector<int>>      m1_history_storage;
    std::vector<std::vector<int>>      m2_history_storage;
    std::vector<std::vector<Particle>> m1_history_particle_storage;
    std::vector<std::vector<Particle>> m2_history_particle_storage;

    std::vector<std::vector<int>>*      m1_history          = nullptr;
    std::vector<std::vector<int>>*      m2_history          = nullptr;
    std::vector<std::vector<Particle>>* m1_history_particle = nullptr;
    std::vector<std::vector<Particle>>* m2_history_particle = nullptr;

    std::vector<int> m1_parton_ancestor_ids;
    std::vector<int> m2_parton_ancestor_ids;
    std::vector<int> m1_parton_ancestor_bars;
    std::vector<int> m2_parton_ancestor_bars;

    std::vector<int> m1_multi_hf_quark_ids;
    std::vector<int> m2_multi_hf_quark_ids;

    // Per-muon first HF hadron kinematics (stored in TruthExtras since they belong to tracing)
    std::vector<float>  m1_first_hf_hadron_prt_pt_eta_phi_m_storage;
    std::vector<float>  m2_first_hf_hadron_prt_pt_eta_phi_m_storage;
    std::vector<float>* m1_first_hf_hadron_prt_pt_eta_phi_m = nullptr;
    std::vector<float>* m2_first_hf_hadron_prt_pt_eta_phi_m = nullptr;

// --------------------- cross-section tracking ---------------------------

    int    pair_counter = 0;
    double total_crossx = 0.;
    double from_resonance_total_crossx = 0.;
    double from_same_b_total_crossx = 0.;
    double either_from_tau_total_crossx = 0.;
    double both_incoming_total_crossx = 0.;
    double both_incoming_FE_QQ_from_same_b_total_crossx = 0.;
    double FE_total_crossx = 0.;
    double from_same_gluon_spitting_total_crossx = 0.;
    double hard_QQ_scatt_total_crossx = 0.;
    double FE_from_same_b_total_crossx = 0.;
    double FE_from_same_GS_total_crossx = 0.;
    double FE_from_diff_ancestors_total_crossx = 0.;
    double FE_from_same_ancestors_not_same_b_or_gs_total_crossx = 0.;
    double skipped_event_crossx = 0.;

// --------------------- histogram & category labels ---------------------------

    std::vector<std::string> parentGroupLabels = {
        "direct b","b to c","direct c","s/light","direct photon","Drell-Yan"};
    int nParentGroups = static_cast<int>(parentGroupLabels.size());

    std::vector<std::string> HardScattCategoryLabels = {
        "gg->gg","qqbar->gg","gq->gq","flavor creation","flavor excitation","QQ'->QQ'","gg->qqbar","qq'->qq'"};
    int nHardScattCategories = static_cast<int>(HardScattCategoryLabels.size());

    std::vector<std::string> HFMuonPairOriginCategoryLabels = {
        "flavor creation","same G(/#gamma)S FSR","same GS ISR 0 hard scatt",
        "same GS ISR 1 hard scatt","same GS ISR both hard scatt","diff GS same hard scatt","others"};
    int nHFMuonPairOriginCategories = static_cast<int>(HFMuonPairOriginCategoryLabels.size());

// --------------------- output histograms ---------------------------

    TH1D* h_numMuonPairs = nullptr;

    TH2D* h_parent_groups[ParamsSet::nSigns][2] = {};
    TH2D* h_hard_scatt_categr_vs_origin_categr[ParamsSet::nSigns][2] = {};
    TH2D* h_both_from_b_hard_scatt_categr_vs_origin_categr[ParamsSet::nSigns][2] = {};
    TH2D* h_both_from_c_hard_scatt_categr_vs_origin_categr[ParamsSet::nSigns][2] = {};
    TH2D* h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[ParamsSet::nSigns][2] = {};

    TH1D* h_muon_pair_origin_categr[ParamsSet::nSigns][2] = {};
    TH1D* h_both_from_b_muon_pair_origin_categr[ParamsSet::nSigns][2] = {};
    TH1D* h_both_from_c_muon_pair_origin_categr[ParamsSet::nSigns][2] = {};
    TH1D* h_one_from_b_one_from_c_muon_pair_origin_categr[ParamsSet::nSigns][2] = {};

    TH1D* h_Qsplit_gs_ISR_one_hard_scatt[ParamsSet::nSigns][2] = {};
    TH1D* h_Qsplit_gs_FSR[ParamsSet::nSigns][2] = {};
    TH1D* h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[ParamsSet::nSigns][2] = {};
    TH1D* h_Qsplit_to_mHard_gs_FSR[ParamsSet::nSigns][2] = {};

// --------------------- output log files ---------------------------

    std::ofstream* m_crossx_summary_file = nullptr;
    std::ofstream* m_unspecified_parent_file = nullptr;
    std::ofstream* m_FE_file = nullptr;
    std::ofstream* m_very_bad_warning_file = nullptr;
    std::ofstream* m_very_low_minv_resonance_file = nullptr;
    std::ofstream* m_hard_scattering_warning_file = nullptr;
    std::ofstream* m_HF_pair_origin_others_category_file = nullptr;
    std::ofstream* m_other_flavor_category_file = nullptr;
    std::ofstream* m_b_parent_file[ParamsSet::nSigns][2] = {};
    std::ofstream* m_c_parent_file[ParamsSet::nSigns][2] = {};

// --------------------- method declarations ---------------------------

    void InitParamsExtra();
    void InitInputExtra();
    void InitTempVariablesExtra();
    void InitOutputHistsExtra();
    void InitOutputExtra();
    void HistAdjustExtra();
    void FinalizeExtra();

    void PerformTruthPairAnalysis();
    void PerPairCrossxUpdate();

    void ResonanceNameMap();
    void FillNumMuonPairsHist(int nPairsAfter, double weight) {
        if (h_numMuonPairs) h_numMuonPairs->Fill(nPairsAfter - 0.5, weight);
    }

    void MuonPairAncestorTracing();
    void MuonPairTagsReinit();
    void SingleMuonAncestorTracing(bool isMuon1);
    std::pair<int,int> UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars,
                                        std::vector<int>& cur_prt_ids, bool before_gs,
                                        bool& prev_out_hard_scatt, int hf_quark_index = -1000000);
    int FindHeavyQuarks(std::vector<int>& cur_prt_ids, std::vector<int>& cur_prt_bars,
                        int quark_type, bool isMuon1, int prev_hq_bar, int hadron_child_id = 0);
    int GluonHistoryTracking(int gluon_bar, bool isMuon1);
    int HardAnalysisCategr(int in_bar1, int in_bar2);
    int ParentGrouping(std::vector<int>& parent_ids, bool c_tag, bool prev_is_lepton);
    void GetPtEtaPhiMFromBarcode(int barcode, std::vector<float>* pt_eta_phi_m);
    void HFMuonPairAnalysis();
    void FillCategoryHistograms(TH1D* horig, TH2D* hhard_vs_orig);
    void CrossxClear();
    void WriteCrossxSummary();
    void PrintHistory(std::ostream* f, bool print_single, bool muon1_sameancestor,
                      bool print_category = false);
};
