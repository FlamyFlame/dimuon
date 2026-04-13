#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TString.h"
#include "TSystem.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../MuonObjectsParamsAndHelpers/Muon.h"
#include "../MuonObjectsParamsAndHelpers/MuonPairPowheg.h"

class PowhegSingleMuonPairMixer {
public:
    std::string mc_mode{"bb"}; // bb or cc

    Long64_t nmixed_pairs_target{5000000};
    Long64_t nmixed_pairs_generated{0};

    int file_batch{1};
    int run_year{17};
    std::string input_file_path{};
    std::string output_file_path{};
    ULong64_t rng_seed{0};

    // When true (default), apply truth_minv ∈ [1, 3] GeV fiducial filter and
    // write output to .../mixed_mass_1_3GeV/.
    // When false, accept all pairs (no minv cut) and write to .../mixed/.
    bool apply_mass_filter{true};

    explicit PowhegSingleMuonPairMixer(const std::string& mode = "bb") : mc_mode(mode) {}

    void Run(){
        ValidateMode();
        if (output_file_path.empty()) output_file_path = DefaultOutputPath();

        std::cout << "[Mixer] mc_mode = " << mc_mode << std::endl;
        std::cout << "[Mixer] file_batch = " << file_batch << std::endl;
        if (input_file_path.empty()){
            std::cout << "[Mixer] input    = combined bb+cc default files" << std::endl;
        }
        else{
            std::cout << "[Mixer] input    = " << input_file_path << std::endl;
        }
        std::cout << "[Mixer] output   = " << output_file_path << std::endl;
        std::cout << "[Mixer] target pairs = " << nmixed_pairs_target << std::endl;
        std::cout << "[Mixer] muon sampling: uniform (no pT weighting)" << std::endl;
        if (apply_mass_filter)
            std::cout << "[Mixer] pair filter: truth_minv in [1.0, 3.0] GeV" << std::endl;
        else
            std::cout << "[Mixer] pair filter: NONE (no mass window applied)" << std::endl;

        const std::size_t out_slash = output_file_path.find_last_of('/');
        if (out_slash != std::string::npos){
            const std::string out_dir = output_file_path.substr(0, out_slash);
            if (!out_dir.empty()) gSystem->mkdir(out_dir.c_str(), true);
        }

        TChain muon_chain("muon_tree");
        if (input_file_path.empty()){
            AddDefaultInputFiles(muon_chain);
        }
        else{
            if (muon_chain.Add(input_file_path.c_str(), 0) == 0){
                throw std::runtime_error("[Mixer] Cannot add input file to TChain: " + input_file_path);
            }
        }
        if (muon_chain.GetEntries() <= 0){
            throw std::runtime_error("[Mixer] No entries found in TChain 'muon_tree'.");
        }

        MuonPowhegFullSimWTruth* in_muon = nullptr;
        muon_chain.SetBranchAddress("MuonObj", &in_muon);

        const Long64_t nmuons = muon_chain.GetEntries();
        if (nmuons < 2){
            throw std::runtime_error("[Mixer] Need at least 2 muons in input tree.");
        }

        std::vector<MuonPowhegFullSimWTruth> muon_cache(static_cast<std::size_t>(nmuons));

        TFile out_file(output_file_path.c_str(), "RECREATE");
        if (!out_file.IsOpen() || out_file.IsZombie()){
            throw std::runtime_error("[Mixer] Cannot create output file: " + output_file_path);
        }

        TTree* muon_pair_tree_sign1 = new TTree("muon_pair_tree_sign1", "All mixed muon pairs, sign1");
        TTree* muon_pair_tree_sign2 = new TTree("muon_pair_tree_sign2", "All mixed muon pairs, sign2");

        MuonPairPowhegFullSimNoTruth* out_pair_ptr = nullptr;
        muon_pair_tree_sign1->Branch("MuonPairObj", &out_pair_ptr);
        muon_pair_tree_sign2->Branch("MuonPairObj", &out_pair_ptr);

        TRandom3 rng(rng_seed);

        for (Long64_t imu = 0; imu < nmuons; ++imu){
            muon_chain.GetEntry(imu);
            muon_cache[static_cast<std::size_t>(imu)] = *in_muon;
        }

        const Long64_t benchmark_step_pairs = 10000;
        const auto t_start = std::chrono::steady_clock::now();
        auto t_last_benchmark = t_start;

        const Long64_t max_attempts = std::max<Long64_t>(1000000, nmixed_pairs_target * 10000);

        auto sample_uniform_index = [&rng, nmuons]() -> Long64_t {
            return static_cast<Long64_t>(rng.Uniform(0.0, static_cast<double>(nmuons)));
        };

        // Track truth_minv range of accepted pairs for post-fill validation.
        float min_stored_minv =  1e10f;
        float max_stored_minv = -1e10f;

        auto fill_one_pair = [&](const MuonPowhegFullSimWTruth& in_m1,
                                 const MuonPowhegFullSimWTruth& in_m2){

            MuonPowhegFullSimNoTruth m1 = ConvertToNoTruth(in_m1);
            MuonPowhegFullSimNoTruth m2 = ConvertToNoTruth(in_m2);

            if (m1.ev_num == m2.ev_num && m1.ind == m2.ind) return false;

            MuonPairPowhegFullSimNoTruth pair;
            pair.Clear();

            if (in_m1.ev_weight == 0.0f || in_m2.ev_weight == 0.0f)
                throw std::runtime_error("[Mixer] ev_weight is zero on input muon. Was the single-muon tree produced with ev_weight filled?");

            const double final_weight = static_cast<double>(in_m1.ev_weight)
                * static_cast<double>(in_m2.ev_weight);
            if (!std::isfinite(final_weight) || !(final_weight > 0.0)) return false;

            pair.weight = final_weight;
            pair.crossx = final_weight;
            pair.m1 = m1;
            pair.m2 = m2;

            pair.Update();

            // NaN-safe filter: negated conjunction correctly rejects NaN
            // (NaN comparisons always return false, so !(NaN>=1 && NaN<=3) == true -> reject).
            if (apply_mass_filter && !(pair.truth_minv >= 1.0f && pair.truth_minv <= 3.0f))
                return false;

            pair.pair_pass_medium = pair.m1.pass_medium && pair.m2.pass_medium;
            pair.pair_pass_tight = pair.m1.pass_tight && pair.m2.pass_tight;
            pair.pair_pass_resonance_reco = false;
            pair.pair_pass_resonance_truth = false;
            pair.pair_pass_medium_and_resonance = false;
            pair.pair_pass_tight_and_resonance = false;

            min_stored_minv = std::min(min_stored_minv, pair.truth_minv);
            max_stored_minv = std::max(max_stored_minv, pair.truth_minv);

            out_pair_ptr = &pair;
            if (pair.truth_same_sign) muon_pair_tree_sign1->Fill();
            else muon_pair_tree_sign2->Fill();

            return true;
        };

        Long64_t attempts = 0;
        while (nmixed_pairs_generated < nmixed_pairs_target && attempts < max_attempts){
            ++attempts;

            const Long64_t idx1 = sample_uniform_index();

            const MuonPowhegFullSimWTruth& first_muon = muon_cache[static_cast<std::size_t>(idx1)];
            const int first_ev_num = first_muon.ev_num;
            const int first_ind = first_muon.ind;

            Long64_t idx2 = -1;
            bool found_distinct_second = false;
            constexpr int max_second_resample = 1000;
            for (int iresample = 0; iresample < max_second_resample; ++iresample){
                const Long64_t idx2_candidate = sample_uniform_index();
                const MuonPowhegFullSimWTruth& second_candidate = muon_cache[static_cast<std::size_t>(idx2_candidate)];

                if (second_candidate.ev_num == first_ev_num && second_candidate.ind == first_ind) continue;

                idx2 = idx2_candidate;
                found_distinct_second = true;
                break;
            }
            if (!found_distinct_second) continue;

            const MuonPowhegFullSimWTruth& second_muon = muon_cache[static_cast<std::size_t>(idx2)];
            const bool filled = fill_one_pair(first_muon, second_muon);
            if (!filled) continue;

            ++nmixed_pairs_generated;

            if (nmixed_pairs_generated % benchmark_step_pairs == 0){
                const auto t_now = std::chrono::steady_clock::now();
                const double elapsed_total_sec =
                    std::chrono::duration<double>(t_now - t_start).count();
                const double elapsed_interval_sec =
                    std::chrono::duration<double>(t_now - t_last_benchmark).count();
                const double filter_eff_pct =
                    100.0 * static_cast<double>(nmixed_pairs_generated)
                    / static_cast<double>(attempts);
                std::cout << "[Mixer] benchmark: generated="
                          << nmixed_pairs_generated
                          << " interval_pairs=" << benchmark_step_pairs
                          << " interval_time_s=" << elapsed_interval_sec
                          << " total_time_s=" << elapsed_total_sec
                          << " filter_eff=" << std::fixed << std::setprecision(3)
                          << filter_eff_pct << "%"
                          << std::endl;
                t_last_benchmark = t_now;
            }

            if (nmixed_pairs_generated % 100000 == 0){
                const double filter_eff_pct =
                    100.0 * static_cast<double>(nmixed_pairs_generated)
                    / static_cast<double>(attempts);
                std::cout << "[Mixer] attempts=" << attempts
                          << " filled_total=" << nmixed_pairs_generated
                          << " filter_eff=" << std::fixed << std::setprecision(3)
                          << filter_eff_pct << "%"
                          << std::endl;
            }
        }

        const double done_filter_eff_pct =
            (attempts > 0)
            ? 100.0 * static_cast<double>(nmixed_pairs_generated) / static_cast<double>(attempts)
            : 0.0;
        std::cout << "[Mixer] done attempts=" << attempts
                  << " filled=" << nmixed_pairs_generated
                  << " / " << nmixed_pairs_target
                  << " filter_eff=" << std::fixed << std::setprecision(3)
                  << done_filter_eff_pct << "%"
                  << std::endl;

        if (nmixed_pairs_generated < nmixed_pairs_target){
            throw std::runtime_error(
                "[Mixer] ERROR: target not reached before max attempts ("
                + std::to_string(attempts) + " attempts, "
                + std::to_string(nmixed_pairs_generated) + " / "
                + std::to_string(nmixed_pairs_target) + " pairs filled). "
                "Increase max_attempts or decrease target_pairs.");
        }

        // Post-fill validation (only when mass filter is active).
        std::cout << "[Mixer] stored truth_minv range: ["
                  << std::fixed << std::setprecision(4) << min_stored_minv
                  << ", " << max_stored_minv << "] GeV" << std::endl;
        if (apply_mass_filter && !(min_stored_minv >= 1.0f && max_stored_minv <= 3.0f)){
            throw std::runtime_error(
                "[Mixer] ERROR: stored pairs have truth_minv outside [1, 3] GeV! "
                "min=" + std::to_string(min_stored_minv) +
                " max=" + std::to_string(max_stored_minv) +
                ". The output file will NOT be written.");
        }

        out_file.cd();
        muon_pair_tree_sign1->Write();
        muon_pair_tree_sign2->Write();
        out_file.Close();
    }

private:
    void ValidateMode() const {
        if (mc_mode != "bb" && mc_mode != "cc"){
            throw std::runtime_error("[Mixer] mc_mode must be 'bb' or 'cc'.");
        }
    }

    void AddDefaultInputFiles(TChain& chain) const {
        const std::string base = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
        const std::string year = std::to_string(run_year);
        const std::string bb_file =
            base + "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp" + year
            + ".bb.Feb2026.v1._MYSTREAM/single_muon_trees_powheg_bb_fullsim_pp" + year + ".root";
        const std::string cc_file =
            base + "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp" + year
            + ".cc.Feb2026.v1._MYSTREAM/single_muon_trees_powheg_cc_fullsim_pp" + year + ".root";

        if (chain.Add(bb_file.c_str(), 0) == 0){
            throw std::runtime_error("[Mixer] Cannot add bb input file to TChain: " + bb_file);
        }
        if (chain.Add(cc_file.c_str(), 0) == 0){
            throw std::runtime_error("[Mixer] Cannot add cc input file to TChain: " + cc_file);
        }
    }

    std::string DefaultOutputPath() const {
        // Use different subdirectory depending on whether the mass filter is applied.
        const std::string subdir = apply_mass_filter ? "mixed_mass_1_3GeV" : "mixed";
        const std::string mixed_dir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/" + subdir;
        gSystem->mkdir(mixed_dir.c_str(), true);
        return mixed_dir + "/muon_pairs_powheg_bbcc_fullsim_mixed_batch" + std::to_string(file_batch) + ".root";
    }

    static MuonPowhegFullSimNoTruth ConvertToNoTruth(const MuonPowhegFullSimWTruth& in){
        MuonPowhegFullSimNoTruth out;
        out.ind = in.ind;
        out.ev_num = in.ev_num;

        out.pt = in.pt;
        out.eta = in.eta;
        out.phi = in.phi;
        out.charge = in.charge;
        out.dP_overP = in.dP_overP;
        out.z0 = in.z0;
        out.d0 = in.d0;
        out.quality = in.quality;
        out.passmu4 = in.passmu4;
        out.pass_tight = in.pass_tight;
        out.trk_pt = in.trk_pt;
        out.trk_eta = in.trk_eta;
        out.trk_phi = in.trk_phi;
        out.trk_charge = in.trk_charge;

        out.truth_pt = in.truth_pt;
        out.truth_eta = in.truth_eta;
        out.truth_phi = in.truth_phi;
        out.truth_bar = in.truth_bar;
        out.truth_charge = in.truth_charge;
        out.ev_weight = in.ev_weight;

        out.pass_medium = in.pass_medium;
        out.reco_match = in.reco_match;
        return out;
    }
};

void run_powheg_single_muon_pair_mixing(
    const std::string& mc_mode = "bb",
    Long64_t target_pairs = 5000000,
    int file_batch = 1,
    int run_year = 17,
    const std::string& input_file = "",
    const std::string& output_file = "",
    ULong64_t seed = 0,
    bool apply_mass_filter = true)
{
    PowhegSingleMuonPairMixer mixer(mc_mode);
    mixer.nmixed_pairs_target = target_pairs;
    mixer.file_batch = file_batch;
    mixer.run_year = run_year;
    mixer.rng_seed = seed;
    mixer.input_file_path = input_file;
    mixer.output_file_path = output_file;
    mixer.apply_mass_filter = apply_mass_filter;
    mixer.Run();
}
