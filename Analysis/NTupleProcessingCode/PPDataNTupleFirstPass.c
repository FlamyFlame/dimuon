#include "PPDataNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"

PPDataNTupleFirstPass::PPDataNTupleFirstPass(){
    cutLabels.assign(cutLabelsPP.begin(), cutLabelsPP.end());
    numCuts = static_cast<int>(CutsPP::nCuts_pp_data);

    requireTight = false;
    isRun3 = true;
    file_batch = 0;

    isPbPb = false;
    PrintInstructions();
    
    std::cout << "PP Data Ntuple processing script:" << std::endl;
    std::cout << "The following public variable(s) **MUST** be set:" << std::endl;
    std::cout << "--> file_batch: int that decides which run3-file batch to process, only has effect when isRun3 is true" << std::endl;
    std::cout << "                     value = 1,2,3,4 (DEFAULT 0 --> MUST BE SET)" << std::endl;
    std::cout << std::endl;

    std::cout << "The following public variable(s) should be checked:" << std::endl;
    std::cout << "--> resonance_cut_mode: integer that determines which set of resonant cuts to apply" << std::endl;
    std::cout << "        * resonance_cut_mode = 0: NO resonance cut" << std::endl;
    std::cout << "        * resonance_cut_mode = 1: old resonance cut (default)" << std::endl;
    std::cout << "        * resonance_cut_mode = 2: new resonance cut" << std::endl;
    std::cout << "        If resonance_cut_mode value is outside {0,1,2}: assume default option" << std::endl;
    std::cout << "--> requireTight: boolean, default false - if true: require tight WP; false; require medium WP" << std::endl;
    std::cout << "--> isRun3: boolean, default true - if true: run run3 data; false: run run2 data" << std::endl;
    std::cout << "--> trigger_mode: INT, value = 1,2,3" << std::endl;
    std::cout << "                       value = 1 (DEFAULT): require single-muon trigger: tag which muon files trigger & if pair passes mu4_mu4_noL1 & 2mu4" << std::endl;
    std::cout << "                       value = 2: require mu4_mu4_noL1" << std::endl;
    std::cout << "                       value = 3: require 2mu4" << std::endl;
    std::cout << std::endl;

    std::cout << "if isRun3, output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024" << std::endl;
    std::cout << "else,      output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2" << std::endl;
    std::cout << "" << std::endl;
}

//initialize the TChain
void PPDataNTupleFirstPass::TChainFill(){

    fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain->SetMakeClass(1);
    if (!isRun3){ //run2
        fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/New_All_DataPP2017_5TeV_Dec2021.root");
    }else{ //run3
        switch (file_batch){
        case 1:
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part1.root");
            break;
        case 2:
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_1.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_2.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_3.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_4.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_5.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_6.root");
            break;
        case 3:
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_1.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_2.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_4.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_5.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_6.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_8.root");
            break;
        case 4:
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_3.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_7.root");
            break;
        default:
            throw std::runtime_error("Run3 file batch invalid/unspecified: have to be between 1 and 4. No result gets run.");
        }
    }
}

void PPDataNTupleFirstPass::InitOutput(){

    std::string tight_suffix = (requireTight)? "_tight" : "";
    
    std::string output_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    std::string run_dir = isRun3? "pp_2024/" : "pp_run2/";
    output_dir += run_dir;

    std::string run_suffix = isRun3? "_pp_2024" : "_pp_run2";

    TrigModeToSuffixMap();

    std::map<int, std::string> run3_batch_suffix_map;
    run3_batch_suffix_map[1] = "_part1";
    run3_batch_suffix_map[2] = "_part3";
    run3_batch_suffix_map[3] = "_part2-1";
    run3_batch_suffix_map[4] = "_part2-2";

    std::string run3_batch_suffix = "";
    if (isRun3) run3_batch_suffix = run3_batch_suffix_map[file_batch];

    std::string resonance_cut_suffix = "";
    switch (resonance_cut_mode){
    case 0:
      resonance_cut_suffix = "_no_res_cut";
      break;
    case 1:
      resonance_cut_suffix = "_old_res_cut";
      break;
    case 2:
      resonance_cut_suffix = "_new_res_cut";
      break;
    default:
      std::cout << "Public variable resonance_cut_mode is set to a value outside {0,1,2}: INVALID. Apply new resonance cuts by default." << std::endl;
      resonance_cut_suffix = "_new_res_cut";
    }
    std::string file_name_base = output_single_muon_tree? "single_muon_trees" : "muon_pairs";

    output_file_path = output_dir + file_name_base + run_suffix + run3_batch_suffix + trig_suffix + tight_suffix + resonance_cut_suffix + ".root";
    output_hist_file_path = output_dir + "hists_cut_acceptance" + run_suffix + run3_batch_suffix + trig_suffix + tight_suffix + resonance_cut_suffix + ".root";

    // ------------------------------------------------------------------------------

    m_outfile=new TFile(output_file_path.c_str(),"recreate");
    
    if (output_single_muon_tree){
        muonOutTree = new TTree("muon_tree","all single muons");
        muonOutTree->Branch("MuonObj",&tempmuon);
    }
    else{ // mode == 4: output muon pair trees
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
            muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair_raw_ptr);
        }
    }

    // ------------------------------------------------------------------------------

    m_outHistFile=new TFile(output_hist_file_path.c_str(),"recreate");
    DimuonAnalysisBaseClass::InitOutput();
}


void PPDataNTupleFirstPass::Run(){

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    std::cout << "Output file is " << m_outfile << std::endl;
    InitInput();
    InitOutput();
    ProcessData();
    HistAdjust();

    m_outfile->Write();
    m_outHistFile->Write();
    std::cout << "Output muon-pair trees have been written to: " << output_file_path << std::endl;
    std::cout << "Output histograms have been written to: " << output_hist_file_path << std::endl;

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << "CPU time used (in seconds) is " << cpu_time_used << std::endl;

}
