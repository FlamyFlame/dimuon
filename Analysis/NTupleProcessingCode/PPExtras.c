#include "PPExtras.h"

template <class Derived>
void PPExtras<Derived>::InitParamsExtra(){
    self().isPbPb = false;
    self().cutLabels.assign(cutLabelsPP.begin(), cutLabelsPP.end());
    self().numCuts = static_cast<int>(CutsPP::nCuts_pp_data);

    if (self().run_year != 24) {
        std::cerr << "Error:: pp run_year must be 24" << std::endl;
        throw std::exception();
    }

    const int file_batch_max = 12;
    if (self().file_batch <= 0 || self().file_batch > file_batch_max) {
        std::cerr << "Error:: pp file_batch invalid! Must be 1-"
                  << file_batch_max
                  << " for 20" << self().run_year << " data" << std::endl;
        throw std::exception();
    }
}

//initialize the TChain
template <class Derived>
void PPExtras<Derived>::PerformTChainFill(){

    self().fChainRef() = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    self().fChainRef()->SetMakeClass(1);

    std::string file_path = self().data_dir + "data_pp" + std::to_string(self().run_year)
                          + "_part" + std::to_string(self().file_batch) + ".root";

    if (!gSystem->AccessPathName(file_path.c_str())) {
        self().fChainRef()->Add(file_path.c_str());
    } else {
        std::cerr << "File does not exist: " << file_path << std::endl;
        throw std::exception();
    }
}
