#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TString.h"
#include "TSystem.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
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
    double truth_pt_sampling_alpha{2.5};
    ULong64_t rng_seed{0};

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
          std::cout << "[Mixer] muon sampling alpha (truth_pt^alpha) = " << truth_pt_sampling_alpha << std::endl;

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

        std::vector<double> muon_sampling_weight(static_cast<std::size_t>(nmuons), 0.0);
        std::vector<double> weight_cdf(static_cast<std::size_t>(nmuons), 0.0);
        double cumsum = 0.0;
        for (Long64_t imu = 0; imu < nmuons; ++imu){
            muon_chain.GetEntry(imu);
            muon_cache[static_cast<std::size_t>(imu)] = *in_muon;
            const double truth_pt = std::max(0.0, static_cast<double>(in_muon->truth_pt));
            double w = std::pow(truth_pt, truth_pt_sampling_alpha);
            if (!std::isfinite(w) || w < 0.0) w = 0.0;
            muon_sampling_weight[static_cast<std::size_t>(imu)] = w;
            cumsum += w;
            weight_cdf[static_cast<std::size_t>(imu)] = cumsum;
        }

        if (!(cumsum > 0.0)){
            throw std::runtime_error("[Mixer] Total muon sampling weight <= 0. Check truth_pt input.");
        }

        const Long64_t benchmark_step_pairs = 10000;
        const auto t_start = std::chrono::steady_clock::now();
        auto t_last_benchmark = t_start;

        const Long64_t max_attempts = std::max<Long64_t>(1000000, nmixed_pairs_target * 200);

        auto sample_weighted_index = [&rng, &weight_cdf](){
            const double r = rng.Uniform(0.0, weight_cdf.back());
            auto it = std::lower_bound(weight_cdf.begin(), weight_cdf.end(), r);
            if (it == weight_cdf.end()) return static_cast<Long64_t>(weight_cdf.size() - 1);
            return static_cast<Long64_t>(std::distance(weight_cdf.begin(), it));
        };

        auto fill_one_pair = [&](const MuonPowhegFullSimWTruth& in_m1,
                                 const MuonPowhegFullSimWTruth& in_m2,
                                 Long64_t idx1,
                                 Long64_t idx2){

            MuonPowhegFullSimNoTruth m1 = ConvertToNoTruth(in_m1);
            MuonPowhegFullSimNoTruth m2 = ConvertToNoTruth(in_m2);

            if (m1.ev_num == m2.ev_num && m1.ind == m2.ind) return false;

            MuonPairPowhegFullSimNoTruth pair;
            pair.Clear();
            const double w1 = muon_sampling_weight[static_cast<std::size_t>(idx1)];
            const double w2 = muon_sampling_weight[static_cast<std::size_t>(idx2)];
            if (!(w1 > 0.0) || !(w2 > 0.0) || !std::isfinite(w1) || !std::isfinite(w2)) return false;

            if (in_m1.ev_weight == 0.0f || in_m2.ev_weight == 0.0f)
                throw std::runtime_error("[Mixer] ev_weight is zero on input muon. Was the single-muon tree produced with ev_weight filled?");

            const double inv_pair_weight = 1.0 / (w1 * w2);
            if (!(inv_pair_weight > 0.0) || !std::isfinite(inv_pair_weight)) return false;

            const double final_weight = inv_pair_weight
                * static_cast<double>(in_m1.ev_weight)
                * static_cast<double>(in_m2.ev_weight);
            if (!std::isfinite(final_weight) || !(final_weight > 0.0)) return false;

            pair.weight = final_weight;
            pair.crossx = final_weight;
            pair.m1 = m1;
            pair.m2 = m2;

            pair.Update();

            pair.pair_pass_medium = pair.m1.pass_medium && pair.m2.pass_medium;
            pair.pair_pass_tight = pair.m1.pass_tight && pair.m2.pass_tight;
            pair.pair_pass_resonance_reco = false;
            pair.pair_pass_resonance_truth = false;
            pair.pair_pass_medium_and_resonance = false;
            pair.pair_pass_tight_and_resonance = false;

            out_pair_ptr = &pair;
            if (pair.truth_same_sign) muon_pair_tree_sign1->Fill();
            else muon_pair_tree_sign2->Fill();

            return true;
        };

        Long64_t attempts = 0;
        while (nmixed_pairs_generated < nmixed_pairs_target && attempts < max_attempts){
            ++attempts;

            const Long64_t idx1 = sample_weighted_index();

            const MuonPowhegFullSimWTruth& first_muon = muon_cache[static_cast<std::size_t>(idx1)];
            const int first_ev_num = first_muon.ev_num;
            const int first_ind = first_muon.ind;

            Long64_t idx2 = -1;
            bool found_distinct_second = false;
            constexpr int max_second_resample = 1000;
            for (int iresample = 0; iresample < max_second_resample; ++iresample){
                const Long64_t idx2_candidate = sample_weighted_index();
                const MuonPowhegFullSimWTruth& second_candidate = muon_cache[static_cast<std::size_t>(idx2_candidate)];

                if (second_candidate.ev_num == first_ev_num && second_candidate.ind == first_ind) continue;

                idx2 = idx2_candidate;
                found_distinct_second = true;
                break;
            }
            if (!found_distinct_second) continue;

            const MuonPowhegFullSimWTruth& second_muon = muon_cache[static_cast<std::size_t>(idx2)];
            const bool filled = fill_one_pair(first_muon, second_muon, idx1, idx2);
            if (!filled) continue;

            ++nmixed_pairs_generated;

            if (nmixed_pairs_generated % benchmark_step_pairs == 0){
                const auto t_now = std::chrono::steady_clock::now();
                const double elapsed_total_sec =
                    std::chrono::duration<double>(t_now - t_start).count();
                const double elapsed_interval_sec =
                    std::chrono::duration<double>(t_now - t_last_benchmark).count();
                std::cout << "[Mixer] benchmark: generated="
                          << nmixed_pairs_generated
                          << " interval_pairs=" << benchmark_step_pairs
                          << " interval_time_s=" << elapsed_interval_sec
                          << " total_time_s=" << elapsed_total_sec
                          << std::endl;
                t_last_benchmark = t_now;
            }

            if (nmixed_pairs_generated % 100000 == 0){
                std::cout << "[Mixer] attempts=" << attempts
                          << " filled_total=" << nmixed_pairs_generated
                          << std::endl;
            }
        }

        std::cout << "[Mixer] done attempts=" << attempts
                  << " filled=" << nmixed_pairs_generated
                  << " / " << nmixed_pairs_target << std::endl;

        if (nmixed_pairs_generated < nmixed_pairs_target){
            std::cerr << "[Mixer] WARNING: target not reached before max attempts." << std::endl;
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
        const std::string mixed_dir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/mixed";
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
    double truth_pt_sampling_alpha = 2.5)
{
    PowhegSingleMuonPairMixer mixer(mc_mode);
    mixer.nmixed_pairs_target = target_pairs;
    mixer.file_batch = file_batch;
    mixer.run_year = run_year;
    mixer.rng_seed = seed;
    mixer.input_file_path = input_file;
    mixer.output_file_path = output_file;
    mixer.truth_pt_sampling_alpha = truth_pt_sampling_alpha;
    mixer.Run();
}
