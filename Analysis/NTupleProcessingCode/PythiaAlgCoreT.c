#include "PythiaAlgCoreT.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "Riostream.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <math.h>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <ctime>

// ---------------------------------------------------------------------------
// Params
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::InitParams_PythiaCore() {
    this->cutLabels = cutLabels_MC;
    this->numCuts = static_cast<int>(cutLabels_MC.size());
    this->isMC = true;

    isRun3 = (run_year > 18);
    if (is_fullsim || is_fullsim_overlay) {
        // Fullsim: load all kn ranges at once; AMI info read from truth directory
        nKinRanges = 6;
        kinRanges  = {8.f, 14.f, 24.f, 40.f, 70.f, 125.f, 300.f};
        // py_dir used for AMI info (same path as truth non-private)
        py_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/";
        // Fullsim NTUPs are here:
        fullsim_input_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
        outfile_name     = "muon_pairs_pythia_fullsim_pp24";
        outhistfile_name = "hists_pythia_ntuple_processing_fullsim_pp24";
        nevents.resize(nKinRanges, 0);
        nevents_accum.resize(nKinRanges, 0);
        njobs_accum.resize(nKinRanges, 0);
        njobs.resize(nKinRanges, 0);
        nentries_per_kin.resize(nKinRanges, 0);
        return;
    }

    if (getIsPrivate()) {
        py_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/";
        nKinRanges = 5;
        kin_dirs = {"k0/", "k1/", "k2/", "k3/", "k4/"};
        kinRanges = {5.f, 10.f, 25.f, 60.f, 120.f, 3200.f};
        nevents_per_file = {10, 100, 5000, 20000, 20000};
        njobs_all_files_combined = {11991, 4500, 500, 125, 125};
    } else {
        const std::string ecom_str = (std::abs(self().E_COM - 5.36) < 0.01) ? "5p36TeV" : "5p02TeV";
        const std::string ecom_subdir = (std::abs(self().E_COM - 5.36) < 0.01) ? "pythia_5p36TeV" : "pythia_5TeV";
        const std::string pythia_local_dir = "/usatlas/u/yuhanguo/dcachearea/pythia_truth_full_sample/";
        const std::string pythia_pnfs_dir = "/pnfs/usatlas.bnl.gov/users/yuhanguo/pythia_truth_full_sample/";
        py_dir = (getUseLocal() ? pythia_local_dir : pythia_pnfs_dir) + ecom_subdir + "/";
        nKinRanges = 6;
        kin_dirs = {"k0/", "k1/", "k2/", "k3/", "k4/", "k5/"};
        kinRanges = {8.f, 14.f, 24.f, 40.f, 70.f, 125.f, 300.f};
        nevents_per_file.resize(6, 0);
        njobs_all_files_combined.resize(6, 0);
        SetKnBatch(batch_num);
        std::cout << "PythiaAlgCoreT: kn_batch=" << kn_batch << ", low=" << kinRanges[kn_batch] << ", high=" << kinRanges[kn_batch + 1] << std::endl;
        if (kn_batch < 0 || kn_batch >= nKinRanges)
            throw std::runtime_error("PythiaAlgCoreT (non-private): kn_batch=" + std::to_string(kn_batch)
                + " (from constructor batch_num=" + std::to_string(batch_num)
                + ") is out of range [0," + std::to_string(nKinRanges) + ").");

        if (getUseLocal()) {
            if (std::abs(self().E_COM - 5.02) > 0.01) {
                throw std::runtime_error("PythiaAlgCoreT (non-private): useLocal=true is only supported for 5.02 TeV.");
            }
            if (kn_batch != 3) {
                throw std::runtime_error(
                    "PythiaAlgCoreT (non-private): local file does not exist for this kn batch. "
                    "Use the 40-70 GeV kn batch (kn=3) when useLocal=true.");
            }
        }

        std::cout << "PythiaAlgCoreT: useLocal=" << getUseLocal() << ", input_dir=" << py_dir << std::endl;
        outfile_name     = "muon_pairs_pythia_" + ecom_str + "_kn" + std::to_string(kn_batch);
        outhistfile_name = "hists_pythia_ntuple_processing_" + ecom_str + "_kn" + std::to_string(kn_batch);
    }
    nevents.resize(nKinRanges, 0);
    nevents_accum.resize(nKinRanges, 0);
    njobs_accum.resize(nKinRanges, 0);
    njobs.resize(nKinRanges, 0);
    nentries_per_kin.resize(nKinRanges, 0);
}

// ---------------------------------------------------------------------------
// SetInputOutputFilesFromBatch  (private only)
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::SetInputOutputFilesFromBatch_PythiaCore() {
    if (!getIsPrivate()) return;

    switch (batch_num) {
    case 1:
        new_run = false;
        batch_suffix = "_allto0318";
        outfile_name = "muon_pairs_pythia_allto0318";
        outhistfile_name = "hists_pythia_ntuple_processing_allto0318";
        job_dirs = {"0317_all_k/", "0318_all_k/", "0318_k0/"};
        kn_in_job = {{true,true,true,true,true},{true,true,true,true,true},{true,false,false,false,false}};
        nfiles_factor = {{20,10,4,1,1},{20,10,4,1,1},{40,0,0,0,0}};
        break;

    case 2:
        new_run = true;
        batch_suffix = "_after0322";
        outfile_name = "muon_pairs_pythia_after0322";
        outhistfile_name = "hists_pythia_ntuple_processing_after0322";
        job_dirs = {"0322_k0_k1/","0323_k0_k1/","0325_all_k/","0401_all_k/","0429_all_k/"};
        kn_in_job = {{true,true,false,false,false},{true,true,false,false,false},
                     {true,true,true,true,true},{true,true,true,true,true},{true,true,true,true,true}};
        nfiles_factor = {{80,20,0,0,0},{80,20,0,0,0},
                         {80,40,4,1,1},{80,40,4,1,1},{80,40,4,1,1}};
        break;

    default:
        std::cerr << "ERROR: batch_num must be 1 or 2 for private sample!" << std::endl;
        throw std::exception();
    }
}

// ---------------------------------------------------------------------------
// InputSanityCheck
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::InputSanityCheck_PythiaCore() {
    if (!getIsPrivate()) return;

    if (job_dirs.size() != kn_in_job.size() || job_dirs.size() != nfiles_factor.size()) {
        std::cout << "job_dirs, kn_in_job and nfiles_factor must have the same length!" << std::endl;
        throw std::exception();
    }
    for (auto& v : kn_in_job) {
        if ((int)v.size() != nKinRanges) {
            std::cout << "Each vector in kn_in_job must have length " << nKinRanges << "!" << std::endl;
            throw std::exception();
        }
    }
    for (auto& v : nfiles_factor) {
        if ((int)v.size() != nKinRanges) {
            std::cout << "Each vector in nfiles_factor must have length " << nKinRanges << "!" << std::endl;
            throw std::exception();
        }
    }
}

// ---------------------------------------------------------------------------
// InitInputPrivate
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::InitInputPrivate_PythiaCore() {
    evChain = new TChain("PyTree", "PyTree");
    metaChain = new TChain("meta_tree", "meta_tree");
    evChain->SetMakeClass(1);
    metaChain->SetMakeClass(1);

    for (int ikin = 0; ikin < nKinRanges; ikin++) {
        std::cout << "Loading kin range " << ikin << " into TChains." << std::endl;
        clock_t t0 = clock();

        for (int jjob = 0; jjob < (int)job_dirs.size(); jjob++) {
            if (!kn_in_job.at(jjob).at(ikin)) continue;

            for (int kbeam = 0; kbeam < nBeamTypes; kbeam++) {
                std::string job_path = py_dir + job_dirs.at(jjob) + kin_dirs.at(ikin) + beam_dirs.at(kbeam);
                if (kbeam < 0 || kbeam >= 4) throw std::out_of_range("InitInputPrivate: kbeam out of range for nfiles_base");
                int nfiles = nfiles_base[kbeam] * nfiles_factor.at(jjob).at(ikin);
                njobs.at(ikin) += nfiles;
                nevents.at(ikin) += nfiles * nevents_per_file.at(ikin);

                for (int lfile = 1; lfile <= nfiles; lfile++) {
                    std::string fpath = job_path + "pytree_" + std::to_string(lfile) + ".root";
                    std::ifstream infile(fpath.c_str());
                    if (!infile.good()) {
                        std::cout << "Warning: file not found: " << fpath << ". Skip." << std::endl;
                        njobs.at(ikin) -= 1;
                        nevents.at(ikin) -= nevents_per_file.at(ikin);
                        continue;
                    }
                    evChain->Add((fpath + "?#PyTree").c_str());
                    metaChain->Add((fpath + "?#meta_tree").c_str());
                }
            }
        }

        double dt = static_cast<double>(clock() - t0) / CLOCKS_PER_SEC;
        std::cout << "Kin range " << ikin << " loaded in " << dt << " s." << std::endl;
        std::cout << "#events in k" << ikin << ": " << nevents[ikin]
                  << ", #jobs: " << njobs[ikin] << std::endl;

        nevents_accum.at(ikin) = (ikin == 0) ? nevents.at(ikin) : nevents_accum.at(ikin-1) + nevents.at(ikin);
        njobs_accum.at(ikin)   = (ikin == 0) ? njobs.at(ikin)   : njobs_accum.at(ikin-1)   + njobs.at(ikin);
    }

    // Update nentries_per_kin from nevents
    for (int i = 0; i < nKinRanges; i++) nentries_per_kin.at(i) = nevents.at(i);

    // Meta-chain branch addresses
    metaChain->SetBranchStatus("*", 0);
    metaChain->SetBranchStatus("efficiency", 1);
    metaChain->SetBranchAddress("efficiency", &efficiency);
    if (new_run) {
        metaChain->SetBranchStatus("eventWeight", 1);
        metaChain->SetBranchAddress("eventWeight", &ev_weight);
    } else {
        metaChain->SetBranchStatus("totalSigma", 1);
        metaChain->SetBranchAddress("totalSigma", &ev_weight);
    }

    // Event-chain branch addresses (core branches)
    evChain->SetBranchStatus("*", 0);
    evChain->SetBranchStatus("QHard",  1); evChain->SetBranchAddress("QHard",  &QHard);
    evChain->SetBranchStatus("pTHat",  1); evChain->SetBranchAddress("pTHat",  &pTHat);
    evChain->SetBranchStatus("mHat",   1); evChain->SetBranchAddress("mHat",   &mHat);

    evChain->SetBranchStatus("truth_mupair_pt1",  1); evChain->SetBranchAddress("truth_mupair_pt1",  &truth_mupair_pt1);
    evChain->SetBranchStatus("truth_mupair_eta1", 1); evChain->SetBranchAddress("truth_mupair_eta1", &truth_mupair_eta1);
    evChain->SetBranchStatus("truth_mupair_phi1", 1); evChain->SetBranchAddress("truth_mupair_phi1", &truth_mupair_phi1);
    evChain->SetBranchStatus("truth_mupair_ch1",  1); evChain->SetBranchAddress("truth_mupair_ch1",  &truth_mupair_ch1);
    evChain->SetBranchStatus("truth_mupair_bar1", 1); evChain->SetBranchAddress("truth_mupair_bar1", &truth_mupair_bar1);

    evChain->SetBranchStatus("truth_mupair_pt2",  1); evChain->SetBranchAddress("truth_mupair_pt2",  &truth_mupair_pt2);
    evChain->SetBranchStatus("truth_mupair_eta2", 1); evChain->SetBranchAddress("truth_mupair_eta2", &truth_mupair_eta2);
    evChain->SetBranchStatus("truth_mupair_phi2", 1); evChain->SetBranchAddress("truth_mupair_phi2", &truth_mupair_phi2);
    evChain->SetBranchStatus("truth_mupair_ch2",  1); evChain->SetBranchAddress("truth_mupair_ch2",  &truth_mupair_ch2);
    evChain->SetBranchStatus("truth_mupair_bar2", 1); evChain->SetBranchAddress("truth_mupair_bar2", &truth_mupair_bar2);
}

// ---------------------------------------------------------------------------
// InitInputCentrProd
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::InitInputCentrProd_PythiaCore() {
    // kn_batch already validated in InitParams_PythiaCore; only load that one range.
    const std::vector<std::string> beam_names = {"pp", "pn", "np", "nn"};

    evChains_kn_beam.resize(nKinRanges);
    nentries_kn_beam.resize(nKinRanges);
    nentries_kn_sum.assign(nKinRanges, 0);
    ami_weight_kn_beam.resize(nKinRanges);

    // Initialise all rows to null/zero; populate only kn_batch
    for (int i = 0; i < nKinRanges; i++) {
        evChains_kn_beam.at(i).assign(nBeamTypes, nullptr);
        nentries_kn_beam.at(i).assign(nBeamTypes, 0);
        ami_weight_kn_beam.at(i).assign(nBeamTypes, 0.);
    }

    std::string ecom_input_tag;
    std::string ami_campaign_tag;
    std::string ami_ecom_tag;
    if (std::abs(self().E_COM - 5.36) < 0.01) {
        ecom_input_tag = "5p36TeV";
        ami_campaign_tag = "mc23";
        ami_ecom_tag = "5p36TeV";
    } else if (std::abs(self().E_COM - 5.02) < 0.01) {
        ecom_input_tag = "5TeV";
        ami_campaign_tag = "mc15";
        ami_ecom_tag = "5TeV";
    } else {
        throw std::runtime_error("InitInputCentrProd: E_COM must be 5.02 or 5.36 TeV.");
    }

    int ikin   = kn_batch;
    int kin_lo = static_cast<int>(kinRanges.at(ikin));
    int kin_hi = static_cast<int>(kinRanges.at(ikin + 1));

    for (int ibeam = 0; ibeam < nBeamTypes; ibeam++) {
        TChain* ch = new TChain("HeavyIonD3PD", "HeavyIonD3PD");
        std::string path = py_dir + "pythia_" + ecom_input_tag + "_" + beam_names.at(ibeam)
            + "_hQCD_DiMu_pTH" + std::to_string(kin_lo) + "_" + std::to_string(kin_hi) + ".TRUTH0.NTUP.root";
        std::ifstream f(path);
        if (!f.good()) {
            std::cerr << "InitInputCentrProd: WARNING - missing file (skipping): " << path << std::endl;
            delete ch;
            continue;
        }
        f.close();

        int add_ret = ch->Add(path.c_str());
        if (add_ret <= 0) {
            std::cerr << "InitInputCentrProd: WARNING - failed to add nominal file: " << path << std::endl;
            delete ch;
            continue;
        }

        ch->SetMakeClass(1);
        Long64_t nent = ch->GetEntries();
        evChains_kn_beam.at(ikin).at(ibeam) = ch;
        nentries_kn_beam.at(ikin).at(ibeam) = nent;
        nentries_kn_sum.at(ikin) += nent;
    }
    std::cout << "CentrProd kn" << ikin << " (" << kin_lo << "-" << kin_hi << " GeV): "
              << nentries_kn_sum.at(ikin) << " total entries" << std::endl;

    // AMI weights for kn_batch only
    for (int ibeam = 0; ibeam < nBeamTypes; ibeam++) {
        if (!evChains_kn_beam.at(ikin).at(ibeam)) continue;
        std::string ami_path = py_dir + "ami_info/ami_info_" + ami_campaign_tag + "_" + ami_ecom_tag + "_Py8EG_A14_" + beam_names.at(ibeam)
            + "_hQCD_DiMu_pTH" + std::to_string(kin_lo) + "_" + std::to_string(kin_hi) + ".txt";
        std::ifstream ami(ami_path);
        if (!ami.good()) {
            std::cerr << "InitInputCentrProd: WARNING - missing AMI file (skipping): " << ami_path << std::endl;
            continue;
        }
        double crossSection = 0., genFiltEff = 0.;
        std::string line;
        while (std::getline(ami, line)) {
            if (line.find("crossSection") != std::string::npos) {
                size_t c = line.find(':');
                if (c != std::string::npos) { std::istringstream(line.substr(c+1)) >> crossSection; break; }
            }
        }
        ami.close(); ami.open(ami_path);
        while (std::getline(ami, line)) {
            if (line.find("genFiltEff") != std::string::npos) {
                size_t c = line.find(':');
                if (c != std::string::npos) { std::istringstream(line.substr(c+1)) >> genFiltEff; break; }
            }
        }
        ami.close();
        ami_weight_kn_beam.at(ikin).at(ibeam) = crossSection * genFiltEff;
    }

    nominal_beam_ratio["pp"] = 4./25.;
    nominal_beam_ratio["pn"] = 6./25.;
    nominal_beam_ratio["np"] = 6./25.;
    nominal_beam_ratio["nn"] = 9./25.;

    // Bind branches on the loaded chains
    auto bind_chain = [&](TChain* ch) {
        ch->SetBranchStatus("*", 0);
        ch->SetBranchStatus("Q",   1); ch->SetBranchAddress("Q",   &Q_float);

        ch->SetBranchStatus("truth_mupair_pt1",  1); ch->SetBranchAddress("truth_mupair_pt1",  &truth_mupair_pt1_f);
        ch->SetBranchStatus("truth_mupair_eta1", 1); ch->SetBranchAddress("truth_mupair_eta1", &truth_mupair_eta1_f);
        ch->SetBranchStatus("truth_mupair_phi1", 1); ch->SetBranchAddress("truth_mupair_phi1", &truth_mupair_phi1_f);
        ch->SetBranchStatus("truth_mupair_ch1",  1); ch->SetBranchAddress("truth_mupair_ch1",  &truth_mupair_ch1);
        ch->SetBranchStatus("truth_mupair_bar1", 1); ch->SetBranchAddress("truth_mupair_bar1", &truth_mupair_bar1);

        ch->SetBranchStatus("truth_mupair_pt2",  1); ch->SetBranchAddress("truth_mupair_pt2",  &truth_mupair_pt2_f);
        ch->SetBranchStatus("truth_mupair_eta2", 1); ch->SetBranchAddress("truth_mupair_eta2", &truth_mupair_eta2_f);
        ch->SetBranchStatus("truth_mupair_phi2", 1); ch->SetBranchAddress("truth_mupair_phi2", &truth_mupair_phi2_f);
        ch->SetBranchStatus("truth_mupair_ch2",  1); ch->SetBranchAddress("truth_mupair_ch2",  &truth_mupair_ch2);
        ch->SetBranchStatus("truth_mupair_bar2", 1); ch->SetBranchAddress("truth_mupair_bar2", &truth_mupair_bar2);
    };

    for (int ibeam = 0; ibeam < nBeamTypes; ibeam++)
        if (evChains_kn_beam.at(ikin).at(ibeam))
            bind_chain(evChains_kn_beam.at(ikin).at(ibeam));
}

// ---------------------------------------------------------------------------
// InitInputFullsim_PythiaCore
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::InitInputFullsim_PythiaCore() {
    // kn_batch is unused for fullsim (all ranges processed together)
    const std::vector<std::string> beam_names = {"pp", "pn", "np", "nn"};

    evChains_kn_beam.resize(nKinRanges);
    nentries_kn_beam.resize(nKinRanges);
    nentries_kn_sum.assign(nKinRanges, 0);
    ami_weight_kn_beam.resize(nKinRanges);

    for (int i = 0; i < nKinRanges; i++) {
        evChains_kn_beam[i].assign(nBeamTypes, nullptr);
        nentries_kn_beam[i].assign(nBeamTypes, 0);
        ami_weight_kn_beam[i].assign(nBeamTypes, 0.);
    }

    for (int ikin = 0; ikin < nKinRanges; ikin++) {
        int kin_lo = static_cast<int>(kinRanges.at(ikin));
        int kin_hi = static_cast<int>(kinRanges.at(ikin + 1));

        for (int ibeam = 0; ibeam < nBeamTypes; ibeam++) {
            std::string fname = fullsim_input_dir
                + "Pythia_5p36TeV_" + beam_names.at(ibeam)
                + "_hQCD_DiMu_pTH" + std::to_string(kin_lo) + "_" + std::to_string(kin_hi)
                + ".FullSimPP24.NTUP.root";

            std::ifstream fin(fname);
            if (!fin.good()) {
                std::cout << "InitInputFullsim: missing file (skip): " << fname << std::endl;
                continue;
            }
            fin.close();

            TChain* ch = new TChain("HeavyIonD3PD", "HeavyIonD3PD");
            int add_ret = ch->Add(fname.c_str());
            if (add_ret <= 0) {
                std::cerr << "InitInputFullsim: failed to add: " << fname << std::endl;
                delete ch;
                continue;
            }
            ch->SetMakeClass(1);
            Long64_t nent = ch->GetEntries();
            evChains_kn_beam[ikin][ibeam] = ch;
            nentries_kn_beam[ikin][ibeam] = nent;
            nentries_kn_sum[ikin] += nent;

            // Disable all branches; Extras will enable what they need
            ch->SetBranchStatus("*", 0);
            // Enable and bind core branches
            ch->SetBranchStatus("truth_muon_barcode", 1);
            ch->SetBranchAddress("truth_muon_barcode", &truth_muon_barcode);

            std::cout << "Loaded " << nent << " events from " << fname << std::endl;

            // AMI weight for this kn/beam
            std::string ami_path = py_dir + "ami_info/ami_info_mc23_5p36TeV_Py8EG_A14_"
                + beam_names.at(ibeam)
                + "_hQCD_DiMu_pTH" + std::to_string(kin_lo) + "_" + std::to_string(kin_hi) + ".txt";
            std::ifstream ami(ami_path);
            if (!ami.good()) {
                std::cerr << "InitInputFullsim: WARNING - missing AMI file: " << ami_path << std::endl;
            } else {
                double crossSection = 0., genFiltEff = 0.;
                std::string line;
                while (std::getline(ami, line)) {
                    if (line.find("crossSection") != std::string::npos) {
                        size_t c = line.find(':');
                        if (c != std::string::npos) { std::istringstream(line.substr(c+1)) >> crossSection; break; }
                    }
                }
                ami.close(); ami.open(ami_path);
                while (std::getline(ami, line)) {
                    if (line.find("genFiltEff") != std::string::npos) {
                        size_t c = line.find(':');
                        if (c != std::string::npos) { std::istringstream(line.substr(c+1)) >> genFiltEff; break; }
                    }
                }
                ami.close();
                ami_weight_kn_beam[ikin][ibeam] = crossSection * genFiltEff;
                std::cout << "  AMI: crossSection=" << crossSection << " pb, genFiltEff=" << genFiltEff
                          << " -> ami_weight=" << ami_weight_kn_beam[ikin][ibeam] << " pb" << std::endl;
            }
        }
    }

    nominal_beam_ratio["pp"] = 4./25.;
    nominal_beam_ratio["pn"] = 6./25.;
    nominal_beam_ratio["np"] = 6./25.;
    nominal_beam_ratio["nn"] = 9./25.;
}

// ---------------------------------------------------------------------------
// InitInput (dispatcher)
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::InitInput_PythiaCore() {
    SetInputOutputFilesFromBatch_PythiaCore();
    InputSanityCheck_PythiaCore();
    if (is_fullsim || is_fullsim_overlay) {
        InitInputFullsim_PythiaCore();
    } else if (getIsPrivate()) {
        InitInputPrivate_PythiaCore();
    } else {
        InitInputCentrProd_PythiaCore();
    }
}

// ---------------------------------------------------------------------------
// OutputTreePathHook / OutputHistPathHook
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::OutputTreePathHook() {
    std::string apply_suffix = turn_data_resonance_cuts_on ? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";
    std::string local_suffix = getUseLocal() ? "_local_batch" : "";

    std::string output_dir;
    if (is_fullsim || is_fullsim_overlay) {
        output_dir = fullsim_input_dir;
    } else if (getIsPrivate()) {
        output_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/";
    } else {
        const std::string ecom_subdir = (std::abs(self().E_COM - 5.36) < 0.01) ? "pythia_5p36TeV" : "pythia_5TeV";
        output_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/" + ecom_subdir + "/";
    }

    this->output_file_path = output_dir + outfile_name + apply_suffix + local_suffix + ".root";
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::OutputHistPathHook() {
    std::string apply_suffix = turn_data_resonance_cuts_on ? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";
    std::string local_suffix = getUseLocal() ? "_local_batch" : "";

    std::string output_dir;
    if (is_fullsim || is_fullsim_overlay) {
        output_dir = fullsim_input_dir;
    } else if (getIsPrivate()) {
        output_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/";
    } else {
        const std::string ecom_subdir = (std::abs(self().E_COM - 5.36) < 0.01) ? "pythia_5p36TeV" : "pythia_5TeV";
        output_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/" + ecom_subdir + "/";
    }

    this->output_hist_file_path = output_dir + outhistfile_name + apply_suffix + local_suffix + ".root";
}

// ---------------------------------------------------------------------------
// InitOutputTreesExtra_PythiaCore
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::InitOutputTreesExtra_PythiaCore() {
    meta_tree_out = new TTree("meta_tree_out", "meta_tree_out");
    for (int ikin = 0; ikin < nKinRanges; ikin++)
        meta_tree_out->Branch(Form("nentries_kin%d_with_3.7GeV_cuts", ikin),
            &nentries_per_kin.at(ikin),
            Form("nentries_kin%d_with_3.7GeV_cuts/L", ikin));

    muonPairOutTreeKinRange.resize(nKinRanges);
    for (int ikin = 0; ikin < nKinRanges; ikin++) {
        muonPairOutTreeKinRange.at(ikin).resize(ParamsSet::nSigns);
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++) {
            muonPairOutTreeKinRange.at(ikin).at(ksign) = new TTree(
                Form("muon_pair_tree_kin%d_sign%u", ikin, ksign+1),
                Form("all muon pairs, kin range%u, sign%u", ikin, ksign+1));
            muonPairOutTreeKinRange.at(ikin).at(ksign)->Branch("MuonPairObj", &(this->mpair_raw_ptr));
        }
    }
}

// ---------------------------------------------------------------------------
// InitializeExtra_PythiaCore
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::InitializeExtra_PythiaCore() {
    if (getIsPrivate() && meta_tree_out) meta_tree_out->Fill();
}

// ---------------------------------------------------------------------------
// FillMuonPair_PythiaCore
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::FillMuonPair_PythiaCore(int pair_ind) {
    auto* p = mpairRef().get();
    if (!p) return;

    if (getIsPrivate()) {
        if (!truth_mupair_bar1 || !truth_mupair_bar2 || !truth_mupair_ch1 || !truth_mupair_ch2 ||
            !truth_mupair_pt1  || !truth_mupair_pt2  || !truth_mupair_eta1 || !truth_mupair_eta2 ||
            !truth_mupair_phi1 || !truth_mupair_phi2)
            throw std::runtime_error("FillMuonPair_PythiaCore (private): null branch vector pointer");
        p->m1.ind          = truth_mupair_bar1->at(pair_ind);
        p->m2.ind          = truth_mupair_bar2->at(pair_ind);
        p->m1.truth_bar    = p->m1.ind;
        p->m2.truth_bar    = p->m2.ind;
        p->m1.truth_charge = truth_mupair_ch1->at(pair_ind);
        p->m2.truth_charge = truth_mupair_ch2->at(pair_ind);
        p->m1.truth_pt     = static_cast<float>(truth_mupair_pt1->at(pair_ind));
        p->m2.truth_pt     = static_cast<float>(truth_mupair_pt2->at(pair_ind));
        p->m1.truth_eta    = static_cast<float>(truth_mupair_eta1->at(pair_ind));
        p->m2.truth_eta    = static_cast<float>(truth_mupair_eta2->at(pair_ind));
        p->m1.truth_phi    = static_cast<float>(truth_mupair_phi1->at(pair_ind));
        p->m2.truth_phi    = static_cast<float>(truth_mupair_phi2->at(pair_ind));
        p->QHard = QHard;  p->pTHat = pTHat;  p->mHat = mHat;
    } else {
        if (!truth_mupair_bar1 || !truth_mupair_bar2 || !truth_mupair_ch1 || !truth_mupair_ch2 ||
            !truth_mupair_pt1_f  || !truth_mupair_pt2_f  || !truth_mupair_eta1_f || !truth_mupair_eta2_f ||
            !truth_mupair_phi1_f || !truth_mupair_phi2_f)
            throw std::runtime_error("FillMuonPair_PythiaCore (non-private): null branch vector pointer");
        p->m1.ind          = truth_mupair_bar1->at(pair_ind);
        p->m2.ind          = truth_mupair_bar2->at(pair_ind);
        p->m1.truth_bar    = p->m1.ind;
        p->m2.truth_bar    = p->m2.ind;
        p->m1.truth_charge = truth_mupair_ch1->at(pair_ind);
        p->m2.truth_charge = truth_mupair_ch2->at(pair_ind);
        p->m1.truth_pt     = truth_mupair_pt1_f->at(pair_ind) / 1000.;  // Convert MeV to GeV
        p->m2.truth_pt     = truth_mupair_pt2_f->at(pair_ind) / 1000.;  // Convert MeV to GeV
        p->m1.truth_eta    = truth_mupair_eta1_f->at(pair_ind);
        p->m2.truth_eta    = truth_mupair_eta2_f->at(pair_ind);
        p->m1.truth_phi    = truth_mupair_phi1_f->at(pair_ind);
        p->m2.truth_phi    = truth_mupair_phi2_f->at(pair_ind);
        p->QHard = static_cast<double>(Q_float);
        p->pTHat = -1000.;  p->mHat = -1000.;
    }
    p->effcy = efficiency;
}

// ---------------------------------------------------------------------------
// PassCuts_PythiaCore
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
bool PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::PassCuts_PythiaCore() {
    auto* p = mpairRef().get();
    if (!p) return false;

    if (std::fabs(p->m1.truth_eta) > 2.4 || std::fabs(p->m2.truth_eta) > 2.4) return false;
    h_cutAcceptanceRef()[p->m1.truth_charge != p->m2.truth_charge]->Fill((int)pass_muon_eta + 0.5, p->weight);

    if (p->m1.truth_pt < 4 || p->m2.truth_pt < 4) return false;
    h_cutAcceptanceRef()[p->m1.truth_charge != p->m2.truth_charge]->Fill((int)pass_muon_pt + 0.5, p->weight);

    return true;
}

// ---------------------------------------------------------------------------
// FillMuonPairTreePythia
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::FillMuonPairTreePythia(int nkin) {
    auto* p = mpairRef().get();
    if (!p) return;
    this->mpair_raw_ptr = p;

    int nsign = p->truth_same_sign ? 0 : 1;

    // Fill global muon-pair tree (from DimuonAlgCoreT)
    if (this->muonPairOutTree[nsign])
        this->muonPairOutTree[nsign]->Fill();

    // Fill kinematic-range-binned tree
    if (nkin >= 0 && nkin < (int)muonPairOutTreeKinRange.size() &&
        nsign < (int)muonPairOutTreeKinRange.at(nkin).size() &&
        muonPairOutTreeKinRange.at(nkin).at(nsign))
        muonPairOutTreeKinRange.at(nkin).at(nsign)->Fill();
}

// ---------------------------------------------------------------------------
// ProcessDataHook (private + non-private loops)
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::ProcessDataHook() {

    // Fullsim branch: loop over all kn × beam, call ProcessEventFullsimHook per event
    if (is_fullsim || is_fullsim_overlay) {
        const std::vector<std::string> beam_names = {"pp", "pn", "np", "nn"};

        for (int ikin = 0; ikin < nKinRanges; ikin++) {
            int kin_lo = static_cast<int>(kinRanges.at(ikin));
            int kin_hi = static_cast<int>(kinRanges.at(ikin + 1));
            if (fill_kn_trees_fullsim) current_ikin = ikin;

            for (int ibeam = 0; ibeam < nBeamTypes; ibeam++) {
                TChain* ch = evChains_kn_beam.at(ikin).at(ibeam);
                Long64_t N_beam = nentries_kn_beam.at(ikin).at(ibeam);
                if (!ch || N_beam == 0) continue;

                double nom_ratio = nominal_beam_ratio.at(beam_names.at(ibeam));
                double ami_w     = ami_weight_kn_beam.at(ikin).at(ibeam);
                fullsim_weight_factor = (N_beam > 0) ? ami_w * nom_ratio / static_cast<double>(N_beam) : 0.;

                Long64_t N_proc = (this->nevents_max <= 0) ? N_beam
                                  : std::min<Long64_t>(N_beam, this->nevents_max);
                std::cout << "Fullsim pTH" << kin_lo << "_" << kin_hi
                          << " beam=" << beam_names.at(ibeam)
                          << " N=" << N_proc << "/" << N_beam
                          << " w_factor=" << fullsim_weight_factor << std::endl;

                for (Long64_t jev = 0; jev < N_proc; jev++) {
                    if (jev % 10000 == 0)
                        std::cout << "  event " << jev << " / " << N_proc << std::endl;
                    int nb = ch->GetEntry(jev);
                    if (nb <= 0) continue;
                    ProcessEventFullsimHook(static_cast<int>(jev));
                }
            }
        }
        return;
    }

    // Non-private branch: one kinematic range per job (kn_batch, validated in InitParams)
    if (!getIsPrivate()) {
        const std::vector<std::string> beam_names = {"pp", "pn", "np", "nn"};

        {
            int ikin = kn_batch;
            current_ikin = ikin;
            std::cout << "ProcessData (non-private): kn" << ikin << " ["
                      << static_cast<int>(kinRanges.at(ikin)) << "-"
                      << static_cast<int>(kinRanges.at(ikin+1)) << " GeV]" << std::endl;

            for (int ibeam = 0; ibeam < nBeamTypes; ibeam++) {
                TChain* ch = evChains_kn_beam.at(ikin).at(ibeam);
                Long64_t N_beam = nentries_kn_beam.at(ikin).at(ibeam);
                if (N_beam == 0) continue;
                double nom_ratio = nominal_beam_ratio.at(beam_names.at(ibeam));
                double w_norm = ami_weight_kn_beam.at(ikin).at(ibeam) * nom_ratio
                                / static_cast<double>(N_beam);
                efficiency = 1.;
                Long64_t N_to_process = (this->nevents_max <= 0)
                    ? N_beam : std::min<Long64_t>(N_beam, this->nevents_max);
                std::cout << "  beam " << beam_names.at(ibeam) << ": " << N_to_process
                          << " / " << N_beam << " events" << std::endl;

                for (Long64_t jevent = 0; jevent < N_to_process; jevent++) {
                    if (jevent % 10000 == 0)
                        std::cout << "  event " << jevent << " / " << N_to_process << std::endl;
                    int nb = ch->GetEntry(jevent);
                    if (nb == 0) continue;
                    if (!truth_mupair_pt1_f)
                        throw std::runtime_error("ProcessData (non-private): truth_mupair_pt1_f is null after GetEntry");

                    this->muon_pair_list_cur_event_pre_resonance_cut.clear();
                    this->resonance_tagged_muon_index_list_v2.clear();
                    this->resonance_tagged_muon_index_list.clear();
                    int NPairsAfter = 0;

                    int NPairs = static_cast<int>(truth_mupair_pt1_f->size());

                    for (int pair_ind = 0; pair_ind < NPairs; pair_ind++) {
                        mpairRef() = std::make_shared<PairT>();
                        this->FillMuonPair(pair_ind);

                        auto* p = mpairRef().get();
                        p->m1.ev_num = static_cast<int>(jevent);
                        p->m2.ev_num = static_cast<int>(jevent);
                        p->weight = w_norm;
                        p->m1.ev_weight = p->weight;
                        p->m2.ev_weight = p->weight;
                        p->crossx = p->weight * N_beam / efficiency;

                        h_cutAcceptanceRef()[p->m1.truth_charge != p->m2.truth_charge]->Fill(
                            (int)nocut + 0.5, p->weight);
                        if (!this->PassCuts()) continue;

                        mpairRef()->Update();
                        this->ResonanceTagging();
                        this->ResonanceTaggingV2();
                        this->muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(mpairRef()));
                    }

                    for (size_t pair_ind2 = 0;
                         pair_ind2 < this->muon_pair_list_cur_event_pre_resonance_cut.size();
                         pair_ind2++) {
                        mpairRef() = std::move(this->muon_pair_list_cur_event_pre_resonance_cut[pair_ind2]);
                        if (!mpairRef()) continue;

                        auto* p = mpairRef().get();
                        p->Reco_resonance_or_reso_contam_tagged_old = false;

                        auto it1 = std::find(this->resonance_tagged_muon_index_list.begin(),
                                             this->resonance_tagged_muon_index_list.end(), p->m1.ind);
                        auto it2 = std::find(this->resonance_tagged_muon_index_list.begin(),
                                             this->resonance_tagged_muon_index_list.end(), p->m2.ind);
                        if (it1 != this->resonance_tagged_muon_index_list.end() ||
                            it2 != this->resonance_tagged_muon_index_list.end()) {
                            if (turn_data_resonance_cuts_on) continue;
                            p->Reco_resonance_or_reso_contam_tagged_old = true;
                        }
                        it1 = std::find(this->resonance_tagged_muon_index_list_v2.begin(),
                                        this->resonance_tagged_muon_index_list_v2.end(), p->m1.ind);
                        it2 = std::find(this->resonance_tagged_muon_index_list_v2.begin(),
                                        this->resonance_tagged_muon_index_list_v2.end(), p->m2.ind);
                        p->Reco_resonance_or_reso_contam_tagged_new =
                            (it1 != this->resonance_tagged_muon_index_list_v2.end() ||
                             it2 != this->resonance_tagged_muon_index_list_v2.end());

                        h_cutAcceptanceRef()[p->m1.truth_charge != p->m2.truth_charge]->Fill(
                            (int)pass_resonance + 0.5, p->weight);

                        PerformTruthPairAnalysisHook();
                        FillMuonPairTreePythia(ikin);
                        (CallPerPairCrossxUpdate<Extras>(), ...);
                        NPairsAfter++;
                        (CallFillNumMuonPairsHist<Extras>(NPairsAfter, mpairRef().get()->weight), ...);
                    }
                }
            }
        }
        return;
    }

    // Private branch
    for (int ikin = 0; ikin < nKinRanges; ikin++) {
        current_ikin = ikin;
        Long64_t nevent_start = (ikin == 0) ? 0 : nevents_accum.at(ikin-1);
        int njob_start   = (ikin == 0) ? 0 : njobs_accum.at(ikin-1);
        Long64_t nevent_end = nevents_accum.at(ikin);
        if (this->nevents_max > 0)
            nevent_end = std::min<Long64_t>(nevent_end, nevent_start + this->nevents_max);
        std::cout << "Processing kin range " << ikin
                  << ", events " << nevent_start << " to " << nevent_end - 1 << std::endl;

        for (Long64_t jevent = nevent_start; jevent < nevent_end; jevent++) {
            if (jevent % 10000 == 0)
                std::cout << "  event " << jevent << " / " << nevent_end << std::endl;
            int kjob    = njob_start + static_cast<int>((jevent - nevent_start) / nevents_per_file.at(ikin));
            int nb_ev   = evChain->GetEntry(jevent);
            int nb_job  = metaChain->GetEntry(kjob);

            if (nb_ev == 0) {
                std::cout << "Error: zero bytes for event " << jevent << std::endl;
                throw std::exception();
            }
            if (nb_job == 0) {
                std::cout << "Error: zero bytes for job " << kjob << std::endl;
                throw std::exception();
            }
            if (this->nevents_max <= 0 && jevent == nevents_accum.at(ikin) - 1 && kjob != njobs_accum.at(ikin) - 1) {
                std::cout << "Error: job count mismatch." << std::endl;
                throw std::exception();
            }

            this->muon_pair_list_cur_event_pre_resonance_cut.clear();
            this->resonance_tagged_muon_index_list_v2.clear();
            this->resonance_tagged_muon_index_list.clear();
            int NPairsAfter = 0;

            if (!truth_mupair_pt1)
                throw std::runtime_error("ProcessData (private): truth_mupair_pt1 is null after GetEntry");
            int NPairs = static_cast<int>(truth_mupair_pt1->size());

            for (int pair_ind = 0; pair_ind < NPairs; pair_ind++) {
                mpairRef() = std::make_shared<PairT>();
                this->FillMuonPair(pair_ind);

                auto* p = mpairRef().get();
                p->m1.ev_num = static_cast<int>(jevent);
                p->m2.ev_num = static_cast<int>(jevent);
                p->weight = ev_weight / njobs_all_files_combined.at(ikin);
                p->m1.ev_weight = p->weight;
                p->m2.ev_weight = p->weight;
                p->crossx = p->weight * nevents.at(ikin) / efficiency;

                h_cutAcceptanceRef()[p->m1.truth_charge != p->m2.truth_charge]->Fill(
                    (int)nocut + 0.5, p->weight);
                if (!this->PassCuts()) continue;

                mpairRef()->Update();
                this->ResonanceTagging();
                this->ResonanceTaggingV2();
                this->muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(mpairRef()));
            }

            for (int pair_ind2 = 0;
                 pair_ind2 < (int)this->muon_pair_list_cur_event_pre_resonance_cut.size();
                 pair_ind2++) {
                mpairRef() = std::move(this->muon_pair_list_cur_event_pre_resonance_cut[pair_ind2]);
                if (!mpairRef()) continue;

                auto* p = mpairRef().get();
                p->Reco_resonance_or_reso_contam_tagged_old = false;

                auto it1 = std::find(this->resonance_tagged_muon_index_list.begin(),
                                     this->resonance_tagged_muon_index_list.end(), p->m1.ind);
                auto it2 = std::find(this->resonance_tagged_muon_index_list.begin(),
                                     this->resonance_tagged_muon_index_list.end(), p->m2.ind);
                if (it1 != this->resonance_tagged_muon_index_list.end() ||
                    it2 != this->resonance_tagged_muon_index_list.end()) {
                    if (turn_data_resonance_cuts_on) continue;
                    p->Reco_resonance_or_reso_contam_tagged_old = true;
                }
                it1 = std::find(this->resonance_tagged_muon_index_list_v2.begin(),
                                this->resonance_tagged_muon_index_list_v2.end(), p->m1.ind);
                it2 = std::find(this->resonance_tagged_muon_index_list_v2.begin(),
                                this->resonance_tagged_muon_index_list_v2.end(), p->m2.ind);
                p->Reco_resonance_or_reso_contam_tagged_new =
                    (it1 != this->resonance_tagged_muon_index_list_v2.end() ||
                     it2 != this->resonance_tagged_muon_index_list_v2.end());

                h_cutAcceptanceRef()[p->m1.truth_charge != p->m2.truth_charge]->Fill(
                    (int)pass_resonance + 0.5, p->weight);

                PerformTruthPairAnalysisHook();
                FillMuonPairTreePythia(ikin);
                (CallPerPairCrossxUpdate<Extras>(), ...);
                NPairsAfter++;
                (CallFillNumMuonPairsHist<Extras>(NPairsAfter, mpairRef().get()->weight), ...);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Finalize_PythiaCore
// ---------------------------------------------------------------------------

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::Finalize_PythiaCore() {
    // DimuonAlgCoreT::Finalize() writes and closes m_outfile / m_outHistFile.
    // Nothing extra for Core beyond what the base Finalize() already does.
}
