#include "TFile.h"
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

    Long64_t nmixed_pairs_target_pT_lt_20GeV{10000000};
    Long64_t nmixed_pairs_target_pT_gt_20GeV{5000000};

    Long64_t nmixed_pairs_pT_lt_20GeV{0};
    Long64_t nmixed_pairs_pT_gt_20GeV{0};

    int file_batch{1};
    int run_year{17};
    std::string input_file_path{};
    std::string output_file_path{};
    std::string high_pt_sampling_mode{"uniform"}; // uniform or pt_power
    double high_pt_truth_pt_alpha{4.0};
    ULong64_t rng_seed{0};

    explicit PowhegSingleMuonPairMixer(const std::string& mode = "bb") : mc_mode(mode) {}

    void Run(){
        ValidateMode();
        if (input_file_path.empty()) input_file_path = DefaultInputPath();
        if (output_file_path.empty()) output_file_path = DefaultOutputPath(input_file_path);

        std::cout << "[Mixer] mc_mode = " << mc_mode << std::endl;
        std::cout << "[Mixer] file_batch = " << file_batch << std::endl;
        std::cout << "[Mixer] input    = " << input_file_path << std::endl;
        std::cout << "[Mixer] output   = " << output_file_path << std::endl;
        std::cout << "[Mixer] targets  : pT<=20GeV=" << nmixed_pairs_target_pT_lt_20GeV
                  << ", pT>20GeV=" << nmixed_pairs_target_pT_gt_20GeV << std::endl;
        std::cout << "[Mixer] high-pT sampling mode = " << high_pt_sampling_mode
              << ", alpha = " << high_pt_truth_pt_alpha << std::endl;

        TFile in_file(input_file_path.c_str(), "READ");
        if (!in_file.IsOpen() || in_file.IsZombie()){
            throw std::runtime_error("[Mixer] Cannot open input file: " + input_file_path);
        }

        const std::size_t out_slash = output_file_path.find_last_of('/');
        if (out_slash != std::string::npos){
            const std::string out_dir = output_file_path.substr(0, out_slash);
            if (!out_dir.empty()) gSystem->mkdir(out_dir.c_str(), true);
        }

        TTree* muon_tree = static_cast<TTree*>(in_file.Get("muon_tree"));
        if (!muon_tree){
            throw std::runtime_error("[Mixer] Cannot find TTree 'muon_tree' in input file.");
        }

        MuonPowhegFullSimWTruth* in_muon = nullptr;
        muon_tree->SetBranchAddress("MuonObj", &in_muon);

        const Long64_t nmuons = muon_tree->GetEntries();
        if (nmuons < 2){
            throw std::runtime_error("[Mixer] Need at least 2 muons in input tree.");
        }

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

        std::vector<double> high_pt_weight_cdf;
        bool use_weighted_high_pt_sampling = (high_pt_sampling_mode == "pt_power");

        if (use_weighted_high_pt_sampling){
            high_pt_weight_cdf.resize(static_cast<std::size_t>(nmuons), 0.0);
            double cumsum = 0.0;

            for (Long64_t imu = 0; imu < nmuons; ++imu){
                muon_tree->GetEntry(imu);
                const double truth_pt = std::max(0.0, static_cast<double>(in_muon->truth_pt));
                double w = std::pow(truth_pt, high_pt_truth_pt_alpha);
                if (!std::isfinite(w) || w < 0.0) w = 0.0;
                cumsum += w;
                high_pt_weight_cdf[static_cast<std::size_t>(imu)] = cumsum;
            }

            if (!(cumsum > 0.0)){
                std::cout << "[Mixer] WARNING: high-pT weighted sampling requested but total weight <= 0. Fallback to uniform." << std::endl;
                use_weighted_high_pt_sampling = false;
                high_pt_weight_cdf.clear();
            }
        }

        const Long64_t benchmark_target_pT_lt_20GeV = 10000;
        const Long64_t benchmark_target_pT_gt_20GeV = 5000;
        bool benchmark_reported_pT_lt_20GeV = false;
        bool benchmark_reported_pT_gt_20GeV = false;
        const auto t_start = std::chrono::steady_clock::now();

        const Long64_t max_attempts_low = std::max<Long64_t>(1000000, nmixed_pairs_target_pT_lt_20GeV * 200);
        const Long64_t max_attempts_high = std::max<Long64_t>(1000000, nmixed_pairs_target_pT_gt_20GeV * 200);

        auto sample_uniform_index = [&rng, nmuons](){
            return static_cast<Long64_t>(rng.Integer(static_cast<UInt_t>(nmuons)));
        };

        auto sample_weighted_index = [&rng, &high_pt_weight_cdf, &sample_uniform_index](){
            if (high_pt_weight_cdf.empty() || !(high_pt_weight_cdf.back() > 0.0)) return sample_uniform_index();
            const double r = rng.Uniform(0.0, high_pt_weight_cdf.back());
            auto it = std::lower_bound(high_pt_weight_cdf.begin(), high_pt_weight_cdf.end(), r);
            if (it == high_pt_weight_cdf.end()) return static_cast<Long64_t>(high_pt_weight_cdf.size() - 1);
            return static_cast<Long64_t>(std::distance(high_pt_weight_cdf.begin(), it));
        };

        auto fill_one_pair_if_selected = [&](Long64_t idx1, Long64_t idx2, bool require_high_pt){
            if (idx2 == idx1) idx2 = (idx2 + 1) % nmuons;

            muon_tree->GetEntry(idx1);
            MuonPowhegFullSimNoTruth m1 = ConvertToNoTruth(*in_muon);

            muon_tree->GetEntry(idx2);
            MuonPowhegFullSimNoTruth m2 = ConvertToNoTruth(*in_muon);

            MuonPairPowhegFullSimNoTruth pair;
            pair.Clear();
            pair.weight = 1.0;
            pair.crossx = 1.0;
            pair.m1 = m1;
            pair.m2 = m2;

            pair.Update();

            pair.pair_pass_medium = pair.m1.pass_medium && pair.m2.pass_medium;
            pair.pair_pass_tight = pair.m1.pass_tight && pair.m2.pass_tight;
            pair.pair_pass_resonance_reco = false;
            pair.pair_pass_resonance_truth = false;
            pair.pair_pass_medium_and_resonance = false;
            pair.pair_pass_tight_and_resonance = false;

            const float pair_pt = pair.truth_pair_pt;
            const bool pass_category = require_high_pt ? (pair_pt > 20.0f) : (pair_pt <= 20.0f);
            if (!pass_category) return -1.0f;

            out_pair_ptr = &pair;
            if (pair.truth_same_sign) muon_pair_tree_sign1->Fill();
            else muon_pair_tree_sign2->Fill();

            return pair_pt;
        };

        Long64_t attempts_low = 0;
        Long64_t attempts_high = 0;

        while (nmixed_pairs_pT_lt_20GeV < nmixed_pairs_target_pT_lt_20GeV && attempts_low < max_attempts_low){
            ++attempts_low;

            const Long64_t idx1 = sample_uniform_index();
            Long64_t idx2 = sample_uniform_index();
            if (idx2 == idx1) idx2 = (idx2 + 1) % nmuons;

            const float pair_pt = fill_one_pair_if_selected(idx1, idx2, false);
            if (pair_pt < 0.0f) continue;

            ++nmixed_pairs_pT_lt_20GeV;

            const Long64_t nfilled = nmixed_pairs_pT_lt_20GeV + nmixed_pairs_pT_gt_20GeV;

            if (!benchmark_reported_pT_lt_20GeV
                && nmixed_pairs_pT_lt_20GeV >= benchmark_target_pT_lt_20GeV){
                const auto t_now = std::chrono::steady_clock::now();
                const double elapsed_sec =
                    std::chrono::duration<double>(t_now - t_start).count();
                std::cout << "[Mixer] benchmark: time to generate "
                          << benchmark_target_pT_lt_20GeV
                          << " pairs with pT<=20 = "
                          << elapsed_sec << " s" << std::endl;
                benchmark_reported_pT_lt_20GeV = true;
            }

            if (nfilled % 100000 == 0){
                std::cout << "[Mixer] attempts_low=" << attempts_low
                          << " attempts_high=" << attempts_high
                          << " filled_total=" << nfilled
                          << " (pT<=20=" << nmixed_pairs_pT_lt_20GeV
                          << ", pT>20=" << nmixed_pairs_pT_gt_20GeV << ")"
                          << std::endl;
            }
        }

        while (nmixed_pairs_pT_gt_20GeV < nmixed_pairs_target_pT_gt_20GeV && attempts_high < max_attempts_high){
            ++attempts_high;

            const Long64_t idx1 = use_weighted_high_pt_sampling ? sample_weighted_index() : sample_uniform_index();
            Long64_t idx2 = use_weighted_high_pt_sampling ? sample_weighted_index() : sample_uniform_index();
            if (idx2 == idx1) idx2 = (idx2 + 1) % nmuons;

            const float pair_pt = fill_one_pair_if_selected(idx1, idx2, true);
            if (pair_pt < 0.0f) continue;

            ++nmixed_pairs_pT_gt_20GeV;

            const Long64_t nfilled = nmixed_pairs_pT_lt_20GeV + nmixed_pairs_pT_gt_20GeV;

            if (!benchmark_reported_pT_gt_20GeV
                && nmixed_pairs_pT_gt_20GeV >= benchmark_target_pT_gt_20GeV){
                const auto t_now = std::chrono::steady_clock::now();
                const double elapsed_sec =
                    std::chrono::duration<double>(t_now - t_start).count();
                std::cout << "[Mixer] benchmark: time to generate "
                          << benchmark_target_pT_gt_20GeV
                          << " pairs with pT>20 = "
                          << elapsed_sec << " s" << std::endl;
                benchmark_reported_pT_gt_20GeV = true;
            }

            if (nfilled % 100000 == 0){
                std::cout << "[Mixer] attempts_low=" << attempts_low
                          << " attempts_high=" << attempts_high
                          << " filled_total=" << nfilled
                          << " (pT<=20=" << nmixed_pairs_pT_lt_20GeV
                          << ", pT>20=" << nmixed_pairs_pT_gt_20GeV << ")"
                          << std::endl;
            }
        }

        std::cout << "[Mixer] done attempts_low=" << attempts_low
              << ", attempts_high=" << attempts_high
                  << " filled pT<=20=" << nmixed_pairs_pT_lt_20GeV
                  << " / " << nmixed_pairs_target_pT_lt_20GeV
                  << ", pT>20=" << nmixed_pairs_pT_gt_20GeV
                  << " / " << nmixed_pairs_target_pT_gt_20GeV << std::endl;

        if (nmixed_pairs_pT_lt_20GeV < nmixed_pairs_target_pT_lt_20GeV
            || nmixed_pairs_pT_gt_20GeV < nmixed_pairs_target_pT_gt_20GeV){
            std::cerr << "[Mixer] WARNING: targets not reached before max attempts." << std::endl;
        }

        out_file.cd();
        muon_pair_tree_sign1->Write();
        muon_pair_tree_sign2->Write();
        out_file.Close();
        in_file.Close();
    }

private:
    void ValidateMode() const {
        if (mc_mode != "bb" && mc_mode != "cc"){
            throw std::runtime_error("[Mixer] mc_mode must be 'bb' or 'cc'.");
        }
        if (high_pt_sampling_mode != "uniform" && high_pt_sampling_mode != "pt_power"){
            throw std::runtime_error("[Mixer] high_pt_sampling_mode must be 'uniform' or 'pt_power'.");
        }
    }

    std::string DefaultInputPath() const {
        return "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/"
               "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp" + std::to_string(run_year)
               + "." + mc_mode + ".Feb2026.v1._MYSTREAM/"
               "single_muon_trees_powheg_" + mc_mode + "_fullsim_pp" + std::to_string(run_year) + "_w_truth.root";
    }

    std::string DefaultOutputPath(const std::string& input_path) const {
        const std::size_t slash = input_path.find_last_of('/');
        const std::string dir = (slash == std::string::npos) ? std::string(".") : input_path.substr(0, slash);
        const std::string mixed_dir = dir + "/mixed";
        gSystem->mkdir(mixed_dir.c_str(), true);
        return mixed_dir + "/muon_pairs_powheg_" + mc_mode + "_fullsim_mixed_batch" + std::to_string(file_batch) + ".root";
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

        out.pass_medium = in.pass_medium;
        out.reco_match = in.reco_match;
        return out;
    }
};

void run_powheg_single_muon_pair_mixing(
    const std::string& mc_mode = "bb",
    Long64_t target_lt_20 = 10000000,
    Long64_t target_gt_20 = 5000000,
    int file_batch = 1,
    int run_year = 17,
    const std::string& input_file = "",
    const std::string& output_file = "",
    ULong64_t seed = 0,
    const std::string& high_pt_sampling_mode = "uniform",
    double high_pt_truth_pt_alpha = 4.0)
{
    PowhegSingleMuonPairMixer mixer(mc_mode);
    mixer.nmixed_pairs_target_pT_lt_20GeV = target_lt_20;
    mixer.nmixed_pairs_target_pT_gt_20GeV = target_gt_20;
    mixer.file_batch = file_batch;
    mixer.run_year = run_year;
    mixer.rng_seed = seed;
    mixer.input_file_path = input_file;
    mixer.output_file_path = output_file;
    mixer.high_pt_sampling_mode = high_pt_sampling_mode;
    mixer.high_pt_truth_pt_alpha = high_pt_truth_pt_alpha;
    mixer.Run();
}
