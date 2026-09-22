#include "PythiaAlgCoreT.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "Riostream.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <math.h>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <algorithm>   // std::find (AMI DSID provenance check)
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
        nKinRanges = 6;
        kinRanges  = {8.f, 14.f, 24.f, 40.f, 70.f, 125.f, 300.f};
        py_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/";
        fullsim_input_dir = fullsim_input_dir_override.empty()
            ? FullSimSampleInputDir(fullsim_sample_type, isTestSample, overlay_pbpb_year)
            : fullsim_input_dir_override;
        // AMI: the weight is that of the PYTHIA EVGEN the sample was simulated from -- never of the
        // AOD chain, never of "the sample's own ami_info/" (the 2026-09-16 version of this block
        // preferred a copy inside the sample dir, which for the r17864 overlay held the AOD chain's
        // numbers; reverted 2026-09-17, user ruling). PythiaEvgen is derived from the SAME switch
        // as the input directory and the isospin treatment (FullSimSampleType.h), so the three can
        // never disagree. `ami_info_dir_override` remains a diagnostic escape hatch only.
        ami_evgen = FullSimSampleEvgen(fullsim_sample_type, isTestSample);
        if (ami_info_dir_override.empty()) ami_info_dir_override = PythiaEvgenAmiDir(ami_evgen);
        if (expected_ami_dsids.empty())    expected_ami_dsids    = PythiaEvgenDsids(ami_evgen);
        std::cout << "PythiaAlgCoreT: fullsim sample = " << (isTestSample ? "TEST" : "FULL")
                  << ", input_dir=" << fullsim_input_dir
                  << ", evgen=" << (ami_evgen == PythiaEvgen::PDF ? "PDF (pp only, proton PDF)" : "nPDF (4 isospins, nuclear PDF)")
                  << ", ami_dir=" << ami_info_dir_override << std::endl;
        const std::string label = FullSimSampleLabel(fullsim_sample_type, overlay_pbpb_year);
        outfile_name     = "muon_pairs_pythia_fullsim_" + label;
        outhistfile_name = "hists_pythia_ntuple_processing_fullsim_" + label;
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
            // Local mode reads from dcachearea (all kn, both 5.02 and 5.36 TeV)
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
        if (only_pp_isospin && ibeam != 0) continue;
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

    // AMI weights for kn_batch only. The 5.36 TeV truth-only sample IS the nPDF evgen (e8599,
    // 802758-802781, 4 isospins, nuclear PDF), so its weights are PythiaEvgen::nPDF -- read from
    // the evgen's own ami_info_nPDF/ (FullSimSampleType.h), never from py_dir (the LGD symlink
    // farm carries no AMI files). The 5.02 TeV mc15 sample keeps its own ami_info/.
    const bool is_5p36 = (ami_ecom_tag == "5p36TeV");
    for (int ibeam = 0; ibeam < nBeamTypes; ibeam++) {
        if (!evChains_kn_beam.at(ikin).at(ibeam)) continue;
        const std::string ami_path = is_5p36
            ? PythiaEvgenAmiDir(PythiaEvgen::nPDF) + PythiaEvgenAmiFileName(PythiaEvgen::nPDF, beam_names.at(ibeam), kin_lo, kin_hi)
            : "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_" + ami_ecom_tag + "/ami_info/"
              + "ami_info_" + ami_campaign_tag + "_" + ami_ecom_tag + "_Py8EG_A14_" + beam_names.at(ibeam)
              + "_hQCD_DiMu_pTH" + std::to_string(kin_lo) + "_" + std::to_string(kin_hi) + ".txt";
        std::ifstream ami(ami_path);
        if (!ami.good()) {
            // Fatal, not a warning: a skipped AMI file left ami_weight = 0 for that beam, i.e. a
            // silently zero-weighted isospin component (docs/ami_weights.md).
            throw std::runtime_error("InitInputCentrProd: missing AMI file: " + ami_path);
        }
        double crossSection = 0., genFiltEff = 0.;
        int datasetNumber = -1;
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
        ami.close(); ami.open(ami_path);
        while (std::getline(ami, line)) {
            if (line.find("datasetNumber") != std::string::npos) {
                size_t c = line.find(':');
                if (c != std::string::npos) { std::istringstream(line.substr(c+1)) >> datasetNumber; break; }
            }
        }
        ami.close();
        if (is_5p36) {
            const std::vector<int> dsids = PythiaEvgenDsids(PythiaEvgen::nPDF);
            if (std::find(dsids.begin(), dsids.end(), datasetNumber) == dsids.end())
                throw std::runtime_error("InitInputCentrProd: AMI PROVENANCE MISMATCH: " + ami_path
                    + " has datasetNumber=" + std::to_string(datasetNumber) + ", not an nPDF evgen DSID (802758-802781)");
        }
        std::cout << "  AMI " << beam_names.at(ibeam) << " pTH" << kin_lo << "_" << kin_hi << ": DSID=" << datasetNumber
                  << " crossSection=" << crossSection << " nb, genFiltEff=" << genFiltEff
                  << " -> ami_weight=" << crossSection * genFiltEff << " nb (" << ami_path << ")" << std::endl;
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

    const std::string file_tag = FullSimSampleFileTag(fullsim_sample_type, overlay_pbpb_year);
    // Beam content follows the SIMULATED SYSTEM, not the file layout: pp-conditions
    // fullsim = pp collisions = the pp beam alone; the HIJING overlay = Pb+Pb = the 4
    // isospin beams.  Overridable per run for the two test samples (PythiaAlgCoreT.h).
    const bool four_beams = UseFourIsospinBeams();
    std::cout << "InitInputFullsim: " << (four_beams ? "4 isospin beams {pp,pn,np,nn}"
                                                     : "pp beam only (isospin weight 1)")
              << std::endl;

    for (int ikin = 0; ikin < nKinRanges; ikin++) {
        int kin_lo = static_cast<int>(kinRanges.at(ikin));
        int kin_hi = static_cast<int>(kinRanges.at(ikin + 1));

        for (int ibeam = 0; ibeam < nBeamTypes; ibeam++) {
            if (!four_beams && ibeam != 0) continue;

            const std::string fbase = fullsim_input_dir
                + "Pythia_5p36TeV_" + beam_names.at(ibeam)
                + "_hQCD_DiMu_pTH" + std::to_string(kin_lo) + "_" + std::to_string(kin_hi)
                + "." + file_tag + ".NTUP";

            // A slice is EITHER one hadded file (the test sample, and any locally downloaded
            // slice) OR the N unmerged grid outputs of the full sample, exposed as a
            // LOCALGROUPDISK symlink farm ".NTUP.partNN.root" -- the ~280 GB full sample is
            // never hadded onto GPFS (tracking doc pythia_fullsim_pp24_full_sample_skim.md,
            // Design Decision D1).
            // The two forms are MUTUALLY EXCLUSIVE: if both are present we THROW (below)
            // rather than pick one, because a stale hadded file silently shadowing a farm
            // would give a wrong N_beam with an unchanged sigma -> wrong normalization.
            const std::string fname = fbase + ".root";
            const std::string fparts = fbase + ".part*.root";

            std::ifstream fin(fname);
            const bool have_single = fin.good();
            fin.close();

            // Ambiguous input is never intentional: a stale local hadded file would silently
            // SHADOW a (possibly more complete) farm, and N_beam would then be wrong while
            // sigma stayed the same -> wrong normalization, no error. Refuse instead.
            if (have_single) {
                TChain probe("HeavyIonD3PD");
                if (probe.Add(fparts.c_str()) > 0)
                    throw std::runtime_error("InitInputFullsim: AMBIGUOUS input for pTH"
                        + std::to_string(kin_lo) + "_" + std::to_string(kin_hi) + " beam "
                        + beam_names.at(ibeam) + ": BOTH the hadded '" + fname
                        + "' and the multi-part '" + fparts + "' exist. Remove one -- silently "
                          "preferring either risks a wrong event count N_beam in the "
                          "cross-section weight sigma*eff/N_beam.");
            }

            TChain* ch = new TChain("HeavyIonD3PD", "HeavyIonD3PD");
            // TChain::Add expands the glob and returns the number of files matched.
            const int add_ret = have_single ? ch->Add(fname.c_str())
                                            : ch->Add(fparts.c_str());
            if (add_ret <= 0) {
                delete ch;
                // The pp beam is ALWAYS required. Silently dropping it would remove an entire
                // pT-hat slice from a sigma-weighted combination -- a BIAS, not just a
                // statistics loss (a partial farm is self-correcting, since N_beam is measured
                // from the files actually chained; a missing slice is not).
                if (ibeam == 0 && !allow_missing_slices)
                    throw std::runtime_error("InitInputFullsim: NO input for pTH"
                        + std::to_string(kin_lo) + "_" + std::to_string(kin_hi)
                        + " beam pp -- looked for '" + fname + "' and '" + fparts
                        + "'. Refusing to run: a missing pT-hat slice biases the "
                          "cross-section-weighted combination. (Symlink farm unreachable?)");
                if (ibeam == 0)
                    std::cout << "InitInputFullsim: WARNING - missing pT-hat slice pTH"
                              << kin_lo << "_" << kin_hi << " beam pp, SKIPPED "
                              << "(allow_missing_slices=true; diagnostic run only -- any "
                              << "cross-section-weighted result from this run is BIASED)" << std::endl;
                else
                    std::cout << "InitInputFullsim: no input (skip): " << fparts << std::endl;
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

            if (nent == 0)
                throw std::runtime_error("InitInputFullsim: 0 entries for pTH"
                    + std::to_string(kin_lo) + "_" + std::to_string(kin_hi) + " beam "
                    + beam_names.at(ibeam) + " -- input is present but unreadable. Refusing to "
                    "run: a silently-dropped pT-hat slice biases the sigma-weighted combination.");

            std::cout << "Loaded " << nent << " events from "
                      << (have_single ? fname : fparts)
                      << " (" << add_ret << " file" << (add_ret == 1 ? "" : "s") << ")" << std::endl;

            // ---- AMI weight for this kn/beam ----
            // PROVENANCE (BLOCKING): the AMI file name is keyed by BEAM+SLICE only, so it does
            // NOT distinguish productions. The pp24 TEST sample (DSIDs 802758-802781) and the
            // pp24 FULL "_pdf" sample (803015-803020) have DIFFERENT cross-sections -- e.g.
            // pTH8_14: 23.78 nb (802781) vs 36.74 nb (803020), a 55% difference, and the
            // difference VARIES BY SLICE, so it does NOT cancel in any ratio. Reading the wrong
            // production's AMI would silently corrupt every sigma-weighted quantity.
            // => the AMI dir is overridable per sample, and the DSID in the file is CHECKED.
            const std::string ami_dir = ami_info_dir_override.empty()
                ? PythiaEvgenAmiDir(ami_evgen) : ami_info_dir_override;
            std::string ami_path = ami_dir + PythiaEvgenAmiFileName(ami_evgen, beam_names.at(ibeam), kin_lo, kin_hi);
            std::ifstream ami(ami_path);
            if (!ami.good()) {
                // Fatal, not a warning: a missing AMI file left ami_weight = 0, which silently
                // gives this pT-hat slice ZERO weight in every weighted quantity.
                throw std::runtime_error("InitInputFullsim: missing AMI file: " + ami_path
                    + " -- refusing to run with a zero-weight pT-hat slice.");
            } else {
                double crossSection = 0., genFiltEff = 0.;
                int datasetNumber = -1;
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
                ami.close(); ami.open(ami_path);
                while (std::getline(ami, line)) {
                    if (line.find("datasetNumber") != std::string::npos) {
                        size_t c = line.find(':');
                        if (c != std::string::npos) { std::istringstream(line.substr(c+1)) >> datasetNumber; break; }
                    }
                }
                ami.close();

                if (!expected_ami_dsids.empty()
                    && std::find(expected_ami_dsids.begin(), expected_ami_dsids.end(), datasetNumber)
                       == expected_ami_dsids.end()) {
                    std::ostringstream oss;
                    oss << "InitInputFullsim: AMI PROVENANCE MISMATCH for pTH" << kin_lo << "_" << kin_hi
                        << " beam " << beam_names.at(ibeam) << ": " << ami_path
                        << " has datasetNumber=" << datasetNumber
                        << ", which is NOT in the expected DSID list for this sample {";
                    for (size_t q = 0; q < expected_ami_dsids.size(); ++q)
                        oss << (q ? "," : "") << expected_ami_dsids[q];
                    oss << "}. You are about to weight this sample with ANOTHER production's "
                           "cross-section. Set ami_info_dir_override to this sample's AMI directory.";
                    throw std::runtime_error(oss.str());
                }
                std::cout << "  AMI DSID=" << datasetNumber << " (" << ami_path << ")" << std::endl;
                // UNITS: ATLAS AMI `crossSection` is in **nb** (NOT pb), e.g. pp pTH8_14
                // crossSection=4.816e6 = 4.816 mb for HardQCD:All pTHat 8-14 GeV (sensible in nb,
                // absurd in pb). So ami_weight (and the per-pair `weight` derived from it below) is
                // in nb. Consumers comparing an absolute MC dsigma to pp DATA (dsigma=N/L, L in pb^-1)
                // must convert nb->pb (x1000). The weight unit CANCELS in reco-eff ratios and
                // area-normalized fit templates, so most uses are unaffected.
                ami_weight_kn_beam[ikin][ibeam] = crossSection * genFiltEff;
                std::cout << "  AMI: crossSection=" << crossSection << " nb, genFiltEff=" << genFiltEff
                          << " -> ami_weight=" << ami_weight_kn_beam[ikin][ibeam] << " nb" << std::endl;
            }
        }
    }

    // Nothing chained at all is never a legitimate diagnostic run -- with allow_missing_slices
    // it would otherwise exit 0 with empty outputs. The typical cause is a directory / file-tag
    // disagreement: fullsim_input_dir_override pointing at one conditions year while
    // overlay_pbpb_year (which names the file, FullSimSampleFileTag) says the other.
    {
        Long64_t n_total = 0;
        for (const auto& n : nentries_kn_sum) n_total += n;
        if (n_total == 0)
            throw std::runtime_error("InitInputFullsim: NO input slice found in '" + fullsim_input_dir
                + "' for file tag '" + file_tag + "' -- do the input directory and overlay_pbpb_year "
                  "(or fullsim_input_dir_override) describe the same production?");
    }

    // Isospin weight.  4:6:6:9 is the ISOSPIN CONTENT OF A Pb NUCLEUS (Z=82, N=126):
    // it exists to combine the four {pp,pn,np,nn} beams into a Pb+Pb collision.  It is
    // therefore meaningful ONLY when all four beams are read.  When the sample is read
    // as a single pp beam -- i.e. the pp-conditions fullsim, which simulates genuine pp
    // collisions and has nothing to isospin-average -- the weight is 1.
    //
    // Applying 4/25 to a single-beam sample (the previous behaviour) is a SLICE-INDEPENDENT
    // factor, and THAT -- not the absence of a weight -- is why it cancels in every ratio.
    // Be precise here, because the distinction is load-bearing:
    //   * reco-eff and detector-response hists ARE weighted. An empty weight *specifier*
    //     resolves to the `weight` COLUMN (RDFBasedHistFillingPythia.cxx:9 maps "" ->
    //     "weight"), so `make_pair(filter, "")` is weighted, not unweighted.
    //   * The MC trigger efficiency likewise carries `weight` in BOTH numerator and
    //     denominator (FillMCTrigEffHists.cxx:341, and the step-3 weight/(eps1*eps2)).
    // A factor that is the same in numerator and denominator AND the same for every pT-hat
    // slice therefore drops out of any efficiency/response ratio.
    // It does NOT cancel in an ABSOLUTE cross-section, where 4/25 = 0.16 was simply wrong.
    // NOTE the corollary: an error in the per-slice AMI cross-section is NOT slice-independent
    // and so does NOT cancel anywhere -- hence the DSID provenance guard above.
    if (UseFourIsospinBeams()) {
        nominal_beam_ratio["pp"] = 4./25.;
        nominal_beam_ratio["pn"] = 6./25.;
        nominal_beam_ratio["np"] = 6./25.;
        nominal_beam_ratio["nn"] = 9./25.;
    } else {
        nominal_beam_ratio["pp"] = 1.0;
        nominal_beam_ratio["pn"] = 0.0;
        nominal_beam_ratio["np"] = 0.0;
        nominal_beam_ratio["nn"] = 0.0;
    }
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
    std::string trig_suffix  = store_mc_trigger ? "_mc_trig" : ""; // never clobber nominal outputs

    std::string output_dir;
    if (is_fullsim || is_fullsim_overlay) {
        output_dir = fullsim_input_dir;
    } else if (getIsPrivate()) {
        output_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/";
    } else {
        const std::string ecom_subdir = (std::abs(self().E_COM - 5.36) < 0.01) ? "pythia_5p36TeV" : "pythia_5TeV";
        output_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/" + ecom_subdir + "/";
    }

    this->output_file_path = output_dir + outfile_name + apply_suffix + local_suffix + trig_suffix + this->extra_output_suffix + ".root";
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::OutputHistPathHook() {
    std::string apply_suffix = turn_data_resonance_cuts_on ? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";
    std::string local_suffix = getUseLocal() ? "_local_batch" : "";
    std::string trig_suffix  = store_mc_trigger ? "_mc_trig" : ""; // never clobber nominal outputs

    std::string output_dir;
    if (is_fullsim || is_fullsim_overlay) {
        output_dir = fullsim_input_dir;
    } else if (getIsPrivate()) {
        output_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/";
    } else {
        const std::string ecom_subdir = (std::abs(self().E_COM - 5.36) < 0.01) ? "pythia_5p36TeV" : "pythia_5TeV";
        output_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/" + ecom_subdir + "/";
    }

    this->output_hist_file_path = output_dir + outhistfile_name + apply_suffix + local_suffix + trig_suffix + this->extra_output_suffix + ".root";
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

    // FULLSIM event bookkeeping -- see the member declarations in PythiaAlgCoreT.h for why.
    // Sized ONCE here and never resized, so the addresses handed to Branch() stay valid.
    meta_nproc_kn_beam.assign(static_cast<size_t>(nKinRanges) * nBeamTypes, 0);
    meta_nbeam_kn_beam.assign(static_cast<size_t>(nKinRanges) * nBeamTypes, 0);
    for (int ikin = 0; ikin < nKinRanges; ikin++) {
        for (int ibeam = 0; ibeam < nBeamTypes; ibeam++) {
            const size_t idx = static_cast<size_t>(ikin) * nBeamTypes + ibeam;
            meta_tree_out->Branch(Form("nproc_kin%d_beam%d", ikin, ibeam),
                &meta_nproc_kn_beam.at(idx), Form("nproc_kin%d_beam%d/L", ikin, ibeam));
            meta_tree_out->Branch(Form("nbeam_kin%d_beam%d", ikin, ibeam),
                &meta_nbeam_kn_beam.at(idx), Form("nbeam_kin%d_beam%d/L", ikin, ibeam));
        }
    }

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

    // Truth analog of the data reco threshold, 4.0 -> 4.5 GeV (user decision 2026-09-08,
    // mu_pt45_gap125_pairpt9_adoption.md §3(a)).
    if (p->m1.truth_pt < 4.5 || p->m2.truth_pt < 4.5) return false;
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

    FillMuonPairTreeKinRangePythia(nkin);
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PythiaAlgCoreT<PairT, MuonT, Derived, Extras...>::FillMuonPairTreeKinRangePythia(int nkin) {
    auto* p = mpairRef().get();
    if (!p) return;
    this->mpair_raw_ptr = p;

    int nsign = p->truth_same_sign ? 0 : 1;

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

                // N_proc is hoisted ABOVE the weight on purpose (2026-09-09). The weight used to
                // divide by N_beam while this loop covers only N_proc = min(N_beam, nevents_max),
                // so a run truncated by nevents_max produced an absolute cross-section low by
                // exactly N_proc/N_beam -- undetectably, since the weight, the NTUP chain entry
                // count and AMI totalEvents all return N_beam and agree with each other.
                // Normalising by the count ACTUALLY PROCESSED is a NO-OP for every nominal run
                // (N_proc == N_beam there) and fixes precisely the truncated case.
                Long64_t N_proc = (this->nevents_max <= 0) ? N_beam
                                  : std::min<Long64_t>(N_beam, this->nevents_max);

                // fullsim_weight_factor is in **nb** (ami_w is nb; see AMI-read comment above).
                fullsim_weight_factor = (N_proc > 0) ? ami_w * nom_ratio / static_cast<double>(N_proc) : 0.;
                std::cout << "Fullsim pTH" << kin_lo << "_" << kin_hi
                          << " beam=" << beam_names.at(ibeam)
                          << " N=" << N_proc << "/" << N_beam
                          << " w_factor=" << fullsim_weight_factor << std::endl;

                // Record BOTH numbers (see PythiaAlgCoreT.h): the weight above is normalised
                // to N_beam while this loop covers N_proc, so a consumer can only detect a
                // truncated run if the file carries them both.
                {
                    const size_t idx = static_cast<size_t>(ikin) * nBeamTypes + ibeam;
                    if (idx < meta_nproc_kn_beam.size()) {
                        meta_nproc_kn_beam.at(idx) = N_proc;
                        meta_nbeam_kn_beam.at(idx) = N_beam;
                    }
                }
                if (N_proc < N_beam) {
                    meta_fullsim_truncated = true;
                    std::cout << "  ##### WARNING: TRUNCATED RUN (nevents_max) -- this chain's"
                              << " pairs are normalised to N_beam=" << N_beam
                              << " but only " << N_proc << " events were processed, so any"
                              << " ABSOLUTE cross-section from this file is low by a factor "
                              << (static_cast<double>(N_proc) / static_cast<double>(N_beam))
                              << ". Ratios (reco-eff etc.) are unaffected. #####" << std::endl;
                }

                for (Long64_t jev = 0; jev < N_proc; jev++) {
                    if (jev % 10000 == 0)
                        std::cout << "  event " << jev << " / " << N_proc << std::endl;
                    int nb = ch->GetEntry(jev);
                    if (nb <= 0) continue;
                    ProcessEventFullsimHook(static_cast<int>(jev));
                }
            }
        }

        // Fill meta_tree_out on the FULLSIM path too. It used to be Filled only under
        // getIsPrivate() (InitializeExtra_PythiaCore), so every fullsim pair file carried an
        // EMPTY meta tree and recorded nothing about how many events it had processed.
        if (meta_tree_out) meta_tree_out->Fill();
        if (meta_fullsim_truncated)
            std::cout << "##### NOTE: at least one (kn,beam) chain was truncated by"
                      << " nevents_max; meta_tree_out records N_proc vs N_beam per chain."
                      << " Do NOT use this file for an absolute cross-section. #####"
                      << std::endl;
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
                efficiency = 1.;
                // Hoisted above the weight for the same reason as the fullsim branch -- see there.
                // NO-OP when nevents_max is unset (N_to_process == N_beam).
                Long64_t N_to_process = (this->nevents_max <= 0)
                    ? N_beam : std::min<Long64_t>(N_beam, this->nevents_max);
                // w_norm (the per-pair `weight`) is in **nb** (ami_weight is nb; see AMI-read comment).
                double w_norm = (N_to_process > 0)
                    ? ami_weight_kn_beam.at(ikin).at(ibeam) * nom_ratio
                      / static_cast<double>(N_to_process)
                    : 0.;
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
                        // Must multiply back the SAME count w_norm divided by, or crossx
                        // stops being sigma*eps_filt*ratio when a run is truncated.
                        p->crossx = p->weight * N_to_process / efficiency;

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
