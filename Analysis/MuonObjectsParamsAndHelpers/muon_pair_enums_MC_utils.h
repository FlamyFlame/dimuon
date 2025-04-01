#pragma once

#include "muon_pair_enums_MC.h"


void origin_grp_map_build(std::map<int,std::string>& origin_grp_map);
void flavor_grp_map_build(std::map<int,std::string>& flavor_grp_map);

void origin_grp_map_build(std::map<int,std::string>& origin_grp_map){
    origin_grp_map[muon_pair_both_from_open_HF_origin_catgr::fc] = "_FC";
    origin_grp_map[muon_pair_both_from_open_HF_origin_catgr::same_gs_fsr] = "_gs_FSR";
    origin_grp_map[muon_pair_both_from_open_HF_origin_catgr::same_phs_fsr] = "_phs_FSR";
    origin_grp_map[muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_zero_hard_scatt] = "_GS_ISR_no_HS";
    origin_grp_map[muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_one_hard_scatt] = "_gs_ISR_one_HS";
    origin_grp_map[muon_pair_both_from_open_HF_origin_catgr::diff_gs_same_hard_scatt] = "_diff_GS_same_HS";
    origin_grp_map[muon_pair_both_from_open_HF_origin_catgr::others] = "_others";
}

void flavor_grp_map_build(std::map<int,std::string>& flavor_grp_map){
	flavor_grp_map[pair_flavor_index::from_resonance] = "_resonance";
    flavor_grp_map[pair_flavor_index::resonance_contaminated] = "_resonance_contaminated";
    flavor_grp_map[pair_flavor_index::from_single_b] = "_single_b";
    flavor_grp_map[pair_flavor_index::bb] = "_bb";
    flavor_grp_map[pair_flavor_index::cc] = "_cc";
    flavor_grp_map[pair_flavor_index::one_b_one_c] = "_one_b_one_c";
    flavor_grp_map[pair_flavor_index::photon_to_dimuon_splitting] = "_photon_splitting";
    flavor_grp_map[pair_flavor_index::drell_yan] = "_drell_yan";
    flavor_grp_map[pair_flavor_index::other_flavor] = "_other_flavors";
}
