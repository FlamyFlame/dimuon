#include "PPDataNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"

PPDataNTupleFirstPass::PPDataNTupleFirstPass(){
    cutLabels.assign(cutLabelsPP.begin(), cutLabelsPP.end());
    numCuts = static_cast<int>(CutsPP::nCuts_pp_data);

    requireTight = false;
    isRun3 = true;
    isPbPb = false;

    PrintInstructions();
}

//initialize the TChain
void PPDataNTupleFirstPass::TChainFill(){

    fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain->SetMakeClass(1);
    if (!isRun3){ //run2
        std::string in_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/data_pp17_part" + std::to_string(file_batch) + ".root";
        fChain->Add(in_file.c_str());
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


void PPDataNTupleFirstPass::Run(){
    in_out_file_dir = isRun3? "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/" : "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/";
    DimuonDataAnalysisBaseClass::Run();
}
