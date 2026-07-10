#include "PythiaFullSimExtras.h"
#include "../Utilities/tchain_helpers.h"

template <class PairT, class MuonT, class Derived>
void PythiaFullSimExtras<PairT, MuonT, Derived>::InitInputExtra(){
    // Bind reco branches on every chain in evChains_kn_beam
    auto bind_reco = [&](TChain* ch) {
        if (!ch) return;
        enable_and_bind(ch, "muon_pt"           , &muon_pt);
        enable_and_bind(ch, "muon_eta"          , &muon_eta);
        enable_and_bind(ch, "muon_phi"          , &muon_phi);
        enable_and_bind(ch, "muon_quality"      , &muon_quality);
        enable_and_bind(ch, "muon_deltaP_overP" , &muon_deltaP_overP);
        enable_and_bind(ch, "muon_d0"           , &muon_d0);
        enable_and_bind(ch, "muon_z0"           , &muon_z0);
        enable_and_bind(ch, "muon_trk_pt"       , &muon_trk_pt);
        enable_and_bind(ch, "muon_trk_eta"      , &muon_trk_eta);
        enable_and_bind(ch, "muon_trk_phi"      , &muon_trk_phi);
        enable_and_bind(ch, "muon_truth_prob"   , &muon_truth_prob);
        enable_and_bind(ch, "muon_truth_barcode", &muon_truth_barcode);
        enable_and_bind(ch, "truth_muon_pt"     , &truth_muon_pt);
        enable_and_bind(ch, "truth_muon_eta"    , &truth_muon_eta);
        enable_and_bind(ch, "truth_muon_phi"    , &truth_muon_phi);
        enable_and_bind(ch, "truth_muon_ch"     , &truth_muon_ch);
    };

    for (int ikin = 0; ikin < self().nKinRanges; ikin++)
        for (int ibeam = 0; ibeam < self().nBeamTypes; ibeam++)
            bind_reco(self().evChains_kn_beam[ikin][ibeam]);

    if (self().store_mc_trigger) {
        // Trigger propagation mode: bind the trigger branches on every chain.
        // A file without them is an old trigger-off skim (e.g. a slice whose _July2026
        // grid job has not landed yet): DROP it entirely. Default-filling passmu4=false
        // for its muons would enter the denominators with an impossible numerator and
        // bias every efficiency built downstream.
        auto bind_trigger = [&](TChain* ch) -> bool {
            if (!ch) return true; // nothing to do
            // The presence check below inspects one file; the per-chain skip equals a
            // per-file skip ONLY while InitInputFullsim builds one-file chains. Guard
            // that assumption: a mixed multi-file chain would silently carry stale
            // trigger values across the file boundary (the exact bias this mode forbids).
            if (ch->GetListOfFiles()->GetEntries() != 1)
                throw std::runtime_error("store_mc_trigger: expected one file per fullsim chain; "
                                         "re-implement the trigger-branch check per file");
            if (!ch->GetBranch("muon_b_HLT_mu4_L1MU3V")) return false;
            enable_and_bind(ch, "muon_b_HLT_mu4_L1MU3V"        , &muon_b_HLT_mu4);
            enable_and_bind(ch, "dimuon_b_HLT_2mu4_L12MU3V_0_02", &dimuon_b_2mu4_mindR);
            enable_and_bind(ch, "muon_pair_muon1_index"         , &muon_pair_muon1_index);
            enable_and_bind(ch, "muon_pair_muon2_index"         , &muon_pair_muon2_index);
            enable_and_bind(ch, "b_HLT_mu4_L1MU3V"              , &b_HLT_mu4);
            enable_and_bind(ch, "b_HLT_2mu4_L12MU3V"            , &b_HLT_2mu4);
            return true;
        };

        for (int ikin = 0; ikin < self().nKinRanges; ikin++) {
            for (int ibeam = 0; ibeam < self().nBeamTypes; ibeam++) {
                TChain*& ch = self().evChains_kn_beam[ikin][ibeam];
                if (!ch) continue;
                if (!bind_trigger(ch)) {
                    const char* fname = ch->GetListOfFiles()->GetEntries() > 0
                        ? ch->GetListOfFiles()->At(0)->GetTitle() : "(unknown file)";
                    std::cerr << "==============================================================\n"
                              << "store_mc_trigger: NO trigger branches in\n  "
                              << fname << "\n"
                              << "  (old trigger-off skim) -> SKIPPING this file entirely.\n"
                              << "  Re-run when its _July2026 re-skim lands.\n"
                              << "==============================================================" << std::endl;
                    self().nentries_kn_sum[ikin] -= self().nentries_kn_beam[ikin][ibeam];
                    self().nentries_kn_beam[ikin][ibeam] = 0;
                    delete ch;
                    ch = nullptr;
                }
            }
        }
    }
}

template <class PairT, class MuonT, class Derived>
bool PythiaFullSimExtras<PairT, MuonT, Derived>::LookupPairPass2mu4(int reco_ind_a, int reco_ind_b){
    if (reco_ind_a < 0 || reco_ind_b < 0) return false; // a leg without a reco muon has no trigger info
    const int i = std::min(reco_ind_a, reco_ind_b); // skim pair block is filled with i < j
    const int j = std::max(reco_ind_a, reco_ind_b);
    for (size_t k = 0; k < muon_pair_muon1_index->size(); k++){
        if (muon_pair_muon1_index->at(k) == i && muon_pair_muon2_index->at(k) == j)
            return dimuon_b_2mu4_mindR->at(k);
    }
    throw std::runtime_error(Form(
        "LookupPairPass2mu4: reco pair (%d,%d) not found in the skim dimuon block (%zu pairs)"
        " -- raw-NTUP pair block inconsistent with muon block", i, j, muon_pair_muon1_index->size()));
}

template <class PairT, class MuonT, class Derived>
bool PythiaFullSimExtras<PairT, MuonT, Derived>::PassMuonMediumCuts(const muon_t& muon){
    if ((muon.quality&1  )==0) return false; // combined muon
    if ((muon.quality&8  )==0) return false; // Medium muon
    if ((muon.quality&32 )==0) return false; // IDCuts
    if ((muon.quality&256)==0) return false; // MuonCuts

    if (fabs(muon.eta) > 2.4) return false;
    if (muon.pt < 4) return false;

    if (muon.dP_overP > self().pmsRef().deltaP_overP_thrsh) return false;

    if (!self().disable_ip_cut) {
        double z0sinTheta = fabs(muon.z0 * sin(2.0*atan(exp(-muon.eta))));
        if (fabs(muon.d0) >= self().pmsRef().d0cut || z0sinTheta >= self().pmsRef().z0cut) return false;
    }

    if (turn_on_track_charge){
        if (muon.trk_charge != muon.charge) return false;
    }
    return true;
}

template <class PairT, class MuonT, class Derived>
void PythiaFullSimExtras<PairT, MuonT, Derived>::CheckBranchPtrsExtra(){
    if (self().debug_mode) std::cout << "Calling PythiaFullSimExtras::CheckBranchPtrsExtra" << std::endl;

    auto require = [&](auto* p, const char* name){
        if (!p) throw std::runtime_error(std::string("Null branch pointer: ") + name);
    };

    require(muon_pt,            "muon_pt");
    require(muon_eta,           "muon_eta");
    require(muon_phi,           "muon_phi");
    require(muon_quality,       "muon_quality");
    require(muon_deltaP_overP,  "muon_deltaP_overP");
    require(muon_d0,            "muon_d0");
    require(muon_z0,            "muon_z0");
    require(muon_trk_pt,        "muon_trk_pt");
    require(muon_trk_eta,       "muon_trk_eta");
    require(muon_trk_phi,       "muon_trk_phi");
    require(muon_truth_prob,    "muon_truth_prob");
    require(muon_truth_barcode, "muon_truth_barcode");
    require(truth_muon_pt,      "truth_muon_pt");
    require(truth_muon_eta,     "truth_muon_eta");
    require(truth_muon_phi,     "truth_muon_phi");
    require(truth_muon_ch,      "truth_muon_ch");
    require(self().truth_muon_barcode, "truth_muon_barcode (core)");

    if (self().store_mc_trigger){
        require(muon_b_HLT_mu4,        "muon_b_HLT_mu4_L1MU3V (store_mc_trigger)");
        require(dimuon_b_2mu4_mindR,   "dimuon_b_HLT_2mu4_L12MU3V_0_02 (store_mc_trigger)");
        require(muon_pair_muon1_index, "muon_pair_muon1_index (store_mc_trigger)");
        require(muon_pair_muon2_index, "muon_pair_muon2_index (store_mc_trigger)");
    }
}

template <class PairT, class MuonT, class Derived>
void PythiaFullSimExtras<PairT, MuonT, Derived>::ProcessEventFullsim(int ev_num){

    // ---- Build list of "real" reco muon truth barcodes (prob > threshold) ----
    std::vector<int> real_muon_truth_barcode_list;
    std::vector<int> real_muon_orig_index;
    std::vector<int> fake_muon_ind_list;

    if (muon_truth_barcode->size() != muon_truth_prob->size())
        throw std::runtime_error("muon_truth_barcode & muon_truth_prob size mismatch!");

    for (int ind = 0; ind < (int)muon_truth_barcode->size(); ind++){
        if (muon_truth_prob->at(ind) > truth_match_prob_thrsh){
            real_muon_truth_barcode_list.push_back(muon_truth_barcode->at(ind));
            real_muon_orig_index.push_back(ind);
        } else
            fake_muon_ind_list.push_back(ind);
    }

    if (muon_truth_barcode->size() != muon_pt->size())
        throw std::runtime_error("muon_truth_barcode & muon_pt size mismatch!");

    // Event weight: precomputed per beam/kn as ami_weight * isospin_ratio / N_beam
    float event_weight = static_cast<float>(self().fullsim_weight_factor);

    // ---- Build truth muon list with optional reco matching ----
    std::vector<muon_t> truth_muon_list;

    int n_pythia_truth_muons = self().GetNPythiaTruthMuons(truth_muon_pt->size());
    int n_reco = (int)muon_pt->size();

    std::vector<bool> reco_claimed(n_reco, false);

    auto fill_reco_quantities = [&](muon_t& m, int reco_ind){
        m.reco_match = true;
        m.pt      = fabs(muon_pt->at(reco_ind)) / 1000.0f;
        m.eta     = muon_eta->at(reco_ind);
        m.phi     = muon_phi->at(reco_ind);
        m.charge  = (muon_pt->at(reco_ind) > 0) ? 1 : -1;

        m.dP_overP = muon_deltaP_overP->at(reco_ind);
        m.z0       = muon_z0->at(reco_ind);
        m.d0       = muon_d0->at(reco_ind);
        m.quality  = muon_quality->at(reco_ind);

        m.trk_pt  = fabs(muon_trk_pt->at(reco_ind)) / 1000.0f;
        m.trk_eta = muon_trk_eta->at(reco_ind);
        m.trk_phi = muon_trk_phi->at(reco_ind);

        if (turn_on_track_charge)
            m.trk_charge = (muon_trk_pt->at(reco_ind) > 0) ? 1 : -1;
        else
            m.trk_charge = 0;

        m.pass_medium = PassMuonMediumCuts(m);
        m.pass_tight  = (m.pass_medium && (m.quality & 16));

        m.reco_ind = reco_ind;
        if (self().store_mc_trigger)
            // per-muon mu4 match, indexed by the raw-NTUP reco index (data-mirror,
            // DimuonDataAlgCoreT: m.passmu4 = muon_b_HLT_mu4->at(m.ind)).
            // NO mu6/mu8 OR-ing: support triggers are disabled everywhere (data D7).
            m.passmu4 = muon_b_HLT_mu4->at(reco_ind);
    };

    // ---- Truth-to-reco matching: barcode only, exclusive (one reco muon per truth) ----
    // ATLAS-standard MC matching for reconstruction efficiency (Run 2 dimuon note
    // ATL-COM-PHYS-2021-1094): a truth muon is "reconstructed" iff a reco muon carries
    // its truth barcode with truthMatchProbability > 0.5. That is the SOLE criterion.
    //
    // No geometric ΔR fallback. A ΔR<0.05-to-nearest-reco fallback (removed 2026-07-10)
    // is ALWAYS wrong here: it recovers genuine physical losses (over-estimating the
    // efficiency), it has no data analog, and — because the correction is applied to
    // data — it would inflate the corrected yield precisely for collinear (ΔR<0.05)
    // pairs. Geometric matching is a distinct alternative method (sensitive to real
    // ΔR<0.05 pairs), never a fallback bolted onto barcode matching.
    //
    // The `reco_claimed` exclusivity guard keeps a reco muon from being assigned to two
    // truth muons (possible only if two truth muons or two prob>0.5 reco muons share a
    // barcode — rare); it makes the assignment order-independent.
    std::vector<int> truth_to_reco(n_pythia_truth_muons, -1);

    for (int truth_ind = 0; truth_ind < n_pythia_truth_muons; truth_ind++){
        int truth_bar = self().truth_muon_barcode->at(truth_ind);
        for (int list_pos = 0; list_pos < (int)real_muon_truth_barcode_list.size(); list_pos++){
            if (real_muon_truth_barcode_list[list_pos] != truth_bar) continue;
            int reco_ind = real_muon_orig_index[list_pos];
            if (reco_claimed[reco_ind]) continue; // already taken by another truth muon
            truth_to_reco[truth_ind] = reco_ind;
            reco_claimed[reco_ind]   = true;
            break;
        }
    }

    for (int truth_ind = 0; truth_ind < n_pythia_truth_muons; truth_ind++){
        muon_t cur_muon;
        cur_muon.ind      = truth_ind;
        cur_muon.ev_num   = ev_num;
        cur_muon.ev_weight = event_weight;

        cur_muon.truth_pt     = fabs(truth_muon_pt->at(truth_ind)) / 1000.0f;
        cur_muon.truth_eta    = truth_muon_eta->at(truth_ind);
        cur_muon.truth_phi    = truth_muon_phi->at(truth_ind);
        cur_muon.truth_charge = truth_muon_ch->at(truth_ind);
        cur_muon.truth_bar    = self().truth_muon_barcode->at(truth_ind);

        if (truth_to_reco[truth_ind] >= 0){
            fill_reco_quantities(cur_muon, truth_to_reco[truth_ind]);
        } else {
            cur_muon.reco_match  = false;
            cur_muon.pass_medium = false;
            cur_muon.pass_tight  = false;
        }

        if constexpr (requires { self().FillMuonOverlay(cur_muon); }) {
            self().FillMuonOverlay(cur_muon);
        }

        if (self().output_single_muon_tree){
            // Nominal: truth-fiducial gate (reco-efficiency denominator = truth muons).
            // store_mc_trigger: RECO-based gate — the Step-1 trigger-efficiency
            // denominator is "offline reconstructed muon"; a truth-pT gate would sculpt
            // the reco-pT turn-on near threshold. Loose here (pT > 3, |eta| < 2.6);
            // the exact fiducial (pT > 4, |eta| < 2.4) + WP is applied downstream in RDF
            // on reco quantities. (mc_trigger_efficiency.md §3.1)
            const bool keep = self().store_mc_trigger
                ? (cur_muon.reco_match && cur_muon.pt > 3.0 && fabs(cur_muon.eta) < 2.6)
                : (cur_muon.truth_pt > 4.0 && fabs(cur_muon.truth_eta) < 2.4);
            if (keep){
                self().muon_raw_ptr = &cur_muon;
                self().FillSingleMuonTree();
            }
        }

        truth_muon_list.push_back(std::move(cur_muon));
    }

    if (self().output_single_muon_tree) return;

    // ---- Build truth pairs ----
    self().muon_pair_list_cur_event_pre_resonance_cut.clear();
    self().resonance_tagged_muon_index_list_reco.clear();
    self().resonance_tagged_muon_index_list_truth.clear();

    for (int i = 0; i < (int)truth_muon_list.size() - 1; i++){
        for (int j = i + 1; j < (int)truth_muon_list.size(); j++){
            if (self().mpairRef()) self().mpairRef()->Clear();
            else self().mpairRef() = std::make_shared<pair_t>();

            self().mpairRef()->weight  = event_weight;
            self().mpairRef()->crossx  = event_weight;

            if (self().store_mc_trigger){
                self().mpairRef()->ev_pass_mu4  = b_HLT_mu4;  // event-level decisions, no reco requirement
                self().mpairRef()->ev_pass_2mu4 = b_HLT_2mu4;
            }

            self().mpairRef()->m1 = truth_muon_list.at(i);
            self().mpairRef()->m2 = truth_muon_list.at(j);

            self().h_cutAcceptanceRef()[self().mpairRef()->m1.truth_charge != self().mpairRef()->m2.truth_charge]
                ->Fill(double(nocut) + 0.5, event_weight);

            if (!self().PassCuts()) continue; // truth fiducial cuts

            self().mpairRef()->Update(); // compute reco+truth kinematics

            if constexpr (requires { self().FillPairOverlay(); }) {
                self().FillPairOverlay();
            }

            if (self().mpairRef()->m1.reco_match && self().mpairRef()->m2.reco_match){
                ResonanceTaggingReco();
                ResonanceTaggingTruth();
                self().mpairRef()->pair_pass_medium = (self().mpairRef()->m1.pass_medium && self().mpairRef()->m2.pass_medium);
                self().mpairRef()->pair_pass_tight  = (self().mpairRef()->m1.pass_tight  && self().mpairRef()->m2.pass_tight);
                if (self().store_mc_trigger)
                    self().mpairRef()->pass2mu4 =
                        LookupPairPass2mu4(self().mpairRef()->m1.reco_ind, self().mpairRef()->m2.reco_ind);
            } else {
                self().mpairRef()->pair_pass_medium = false;
                self().mpairRef()->pair_pass_tight  = false;
            }

            self().muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(self().mpairRef()));
        }
    }

    // ---- Second loop: apply resonance veto and fill output ----
    for (int pair_ind = 0; pair_ind < (int)self().muon_pair_list_cur_event_pre_resonance_cut.size(); pair_ind++){
        self().mpairRef() = std::move(self().muon_pair_list_cur_event_pre_resonance_cut.at(pair_ind));
        if (!self().mpairRef()){
            std::cerr << "mpairRef() at second muon-pair loop NOT found!" << std::endl;
            continue;
        }

        if (self().mpairRef()->m1.reco_match && self().mpairRef()->m2.reco_match){
            self().mpairRef()->pair_pass_resonance_reco  = true;
            self().mpairRef()->pair_pass_resonance_truth = true;

            auto it1 = std::find(resonance_tagged_muon_index_list_reco.begin(),
                                 resonance_tagged_muon_index_list_reco.end(),
                                 self().mpairRef()->m1.ind);
            self().mpairRef()->pair_pass_resonance_reco &= (it1 == resonance_tagged_muon_index_list_reco.end());

            auto it2 = std::find(resonance_tagged_muon_index_list_reco.begin(),
                                 resonance_tagged_muon_index_list_reco.end(),
                                 self().mpairRef()->m2.ind);
            self().mpairRef()->pair_pass_resonance_reco &= (it2 == resonance_tagged_muon_index_list_reco.end());

            it1 = std::find(resonance_tagged_muon_index_list_truth.begin(),
                            resonance_tagged_muon_index_list_truth.end(),
                            self().mpairRef()->m1.ind);
            self().mpairRef()->pair_pass_resonance_truth &= (it1 == resonance_tagged_muon_index_list_truth.end());

            it2 = std::find(resonance_tagged_muon_index_list_truth.begin(),
                            resonance_tagged_muon_index_list_truth.end(),
                            self().mpairRef()->m2.ind);
            self().mpairRef()->pair_pass_resonance_truth &= (it2 == resonance_tagged_muon_index_list_truth.end());

            self().mpairRef()->pair_pass_medium_and_resonance =
                self().mpairRef()->pair_pass_resonance_reco && self().mpairRef()->pair_pass_medium;
            self().mpairRef()->pair_pass_tight_and_resonance =
                self().mpairRef()->pair_pass_resonance_reco && self().mpairRef()->pair_pass_tight;
        } else {
            self().mpairRef()->pair_pass_resonance_reco  = false;
            self().mpairRef()->pair_pass_resonance_truth = false;
            self().mpairRef()->pair_pass_medium                = false;
            self().mpairRef()->pair_pass_tight                 = false;
            self().mpairRef()->pair_pass_medium_and_resonance  = false;
            self().mpairRef()->pair_pass_tight_and_resonance   = false;
        }

        if (self().getPerformTruth()) self().PerformTruthPairAnalysisHook();

        self().FillMuonPairTree();
    }
}
