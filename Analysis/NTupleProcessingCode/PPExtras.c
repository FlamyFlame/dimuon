#include "PPExtras.h"

template <class Derived>
void PPExtras<Derived>::InitParamsExtra(){
    self().isPbPb = false;
    self().cutLabels.assign(cutLabelsPP.begin(), cutLabelsPP.end());
    self().numCuts = static_cast<int>(CutsPP::nCuts_pp_data);

    if (self().run_year != 17 && self().run_year != 24) {
        std::cerr << "Error:: pp run_year must be 17 or 24" << std::endl;
        throw std::exception();
    }

    int file_batch_max = (self().run_year == 24) ? 4 : 3;
    if (self().file_batch <= 0 || self().file_batch > file_batch_max) {
        std::cerr << "Error:: pp file_batch invalid! Must be 1-4 for 2024 data, 1-3 for 2017 data" << std::endl;
        throw std::exception();
    }
}

//initialize the TChain
template <class Derived>
void PPExtras<Derived>::PerformTChainFill(){

    self().fChainRef() = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    self().fChainRef()->SetMakeClass(1);
    if (!self().isRun3){ //run2
        std::string in_file = self().data_dir + "data_pp17_part" + std::to_string(self().file_batch) + ".root";
        self().fChainRef()->Add(in_file.c_str());
    }else{ //run3
        switch (self().file_batch){
        case 1:
            self().fChainRef()->Add((self().data_dir + "data_pp24_part1.root").c_str());
            break;
        case 2:
            self().fChainRef()->Add((self().data_dir + "data_pp24_part3/data_pp24_part3_1.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part3/data_pp24_part3_2.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part3/data_pp24_part3_3.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part3/data_pp24_part3_4.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part3/data_pp24_part3_5.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part3/data_pp24_part3_6.root").c_str());
            break;
        case 3:
            self().fChainRef()->Add((self().data_dir + "data_pp24_part2/data_pp24_part2_1.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part2/data_pp24_part2_2.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part2/data_pp24_part2_4.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part2/data_pp24_part2_5.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part2/data_pp24_part2_6.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part2/data_pp24_part2_8.root").c_str());
            break;
        case 4:
            self().fChainRef()->Add((self().data_dir + "data_pp24_part2/data_pp24_part2_3.root").c_str());
            self().fChainRef()->Add((self().data_dir + "data_pp24_part2/data_pp24_part2_7.root").c_str());
            break;
        default:
            throw std::runtime_error("Run3 file batch invalid/unspecified: have to be between 1 and 4. No result gets run.");
        }
    }
}
