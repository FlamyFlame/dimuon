#include "PPExtras.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"

// no default constructor: run_year & file_batch MUST be provided
PPExtras::PPExtras(int run_year_input, int file_batch_input)
    : DimuonDataAnalysisBaseClass(run_year_input, file_batch_input){
        PrintInstructions();
}

template <class Derived>
void PPExtras<Derived>::InitParamsExtra(){
    isPbPb = false;
    cutLabels.assign(cutLabelsPP.begin(), cutLabelsPP.end());
    numCuts = static_cast<int>(CutsPP::nCuts_pp_data);
}

//initialize the TChain
template <class Derived>
void PPExtras<Derived>::TChainFill(){

    fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain->SetMakeClass(1);
    if (!isRun3){ //run2
        std::string in_file = data_dir + "data_pp17_part" + std::to_string(file_batch) + ".root";
        fChain->Add(in_file.c_str());
    }else{ //run3
        switch (file_batch){
        case 1:
            fChain->Add((data_dir + "data_pp24_part1.root").c_str());
            break;
        case 2:
            fChain->Add((data_dir + "data_pp24_part3/data_pp24_part3_1.root").c_str());
            fChain->Add((data_dir + "data_pp24_part3/data_pp24_part3_2.root").c_str());
            fChain->Add((data_dir + "data_pp24_part3/data_pp24_part3_3.root").c_str());
            fChain->Add((data_dir + "data_pp24_part3/data_pp24_part3_4.root").c_str());
            fChain->Add((data_dir + "data_pp24_part3/data_pp24_part3_5.root").c_str());
            fChain->Add((data_dir + "data_pp24_part3/data_pp24_part3_6.root").c_str());
            break;
        case 3:
            fChain->Add((data_dir + "data_pp24_part2/data_pp24_part2_1.root").c_str());
            fChain->Add((data_dir + "data_pp24_part2/data_pp24_part2_2.root").c_str());
            fChain->Add((data_dir + "data_pp24_part2/data_pp24_part2_4.root").c_str());
            fChain->Add((data_dir + "data_pp24_part2/data_pp24_part2_5.root").c_str());
            fChain->Add((data_dir + "data_pp24_part2/data_pp24_part2_6.root").c_str());
            fChain->Add((data_dir + "data_pp24_part2/data_pp24_part2_8.root").c_str());
            break;
        case 4:
            fChain->Add((data_dir + "data_pp24_part2/data_pp24_part2_3.root").c_str());
            fChain->Add((data_dir + "data_pp24_part2/data_pp24_part2_7.root").c_str());
            break;
        default:
            throw std::runtime_error("Run3 file batch invalid/unspecified: have to be between 1 and 4. No result gets run.");
        }
    }
}
