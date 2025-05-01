#include <cmath>
#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TH2D.h>
#include <TFile.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <TChain.h>


using namespace ROOT;
using namespace ROOT::VecOps;

// -------------------------- analysis base class --------------------------
class SingleBAnalysisBase{
protected:
    std::string dataset_base_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/";
    std::string muon_pair_input_file_name;

    constexpr static double pi = 3.14159265358979323846;

    std::string signal_cuts = "minv > 1.08 && minv < 2.9 && pair_pt > 8";

    // Define binning
    std::vector<double> pT_bins_80;
    std::vector<double> pT_bins_120;
    std::vector<double> pT_bins_150;
    std::vector<double> pT_bins_200;
    std::vector<double> pT_bins_300;
    std::vector<double> pT_bins_500;
    constexpr static int n_eta_bins = 10;
    std::vector<double> eta_bins;
    constexpr static double eta_min = -2.5, eta_max = 2.5;

    virtual void Initialize();
    void fillLogBinningArray(std::vector<double>& bins, int nBins, double low, double high);

public:
    SingleBAnalysisBase(){    }
    ~SingleBAnalysisBase(){}

    void BinningPrinting();
    virtual void RunAnalysis() = 0;
};

// -------------------------- base-class member functions --------------------------

void SingleBAnalysisBase::fillLogBinningArray(std::vector<double>& bins, int nBins, double low, double high) {
    bins.clear();

    double logLow = std::log10(low);
    double logHigh = std::log10(high);
    double logStep = (logHigh - logLow) / nBins;

    for (int i = 0; i <= nBins; ++i) {
        bins.push_back(std::pow(10, logLow + i * logStep));
    }
}

void SingleBAnalysisBase::Initialize(){
    // // Create logarithmic pair-pT bin edges
    fillLogBinningArray(pT_bins_80,  12, 8.0, 80.0);  // 10 log bins from 8 to 200 GeV
    fillLogBinningArray(pT_bins_120, 14, 8.0, 120.0);  // 10 log bins from 8 to 120 GeV
    fillLogBinningArray(pT_bins_150, 30, 8.0, 150.0);  // 10 log bins from 8 to 150 GeV
    fillLogBinningArray(pT_bins_200, 18, 8.0, 200.0);  // 10 log bins from 8 to 200 GeV
    fillLogBinningArray(pT_bins_300, 30, 8.0, 300.0);  // 10 log bins from 8 to 200 GeV
    fillLogBinningArray(pT_bins_500, 30, 8.0, 500.0);  // 10 log bins from 8 to 200 GeV

    // Create eta bin edges
    eta_bins.clear();
    for (int i = 0; i <= n_eta_bins; ++i){
        eta_bins.push_back(eta_min + (eta_max - eta_min) * i / n_eta_bins);
    }
}



void SingleBAnalysisBase::BinningPrinting(){
    std::cout << "Print out bins for pair pT 8-80 GeV" << std::endl;
    for (auto bin : pT_bins_80){
        std::cout << bin << ", ";
    }
    std::cout << std::endl;

    std::cout << "Print out bins for pair pT 8-120 GeV" << std::endl;
    for (auto bin : pT_bins_120){
        std::cout << bin << ", ";
    }
    std::cout << std::endl;

    std::cout << "Print out bins for pair pT 8-150 GeV" << std::endl;
    for (auto bin : pT_bins_150){
        std::cout << bin << ", ";
    }
    std::cout << std::endl;

    std::cout << "Print out bins for pair pT 8-200 GeV" << std::endl;
    for (auto bin : pT_bins_200){
        std::cout << bin << ", ";
    }
    std::cout << std::endl;

    std::cout << "Print out bins for pair pT 8-300 GeV" << std::endl;
    for (auto bin : pT_bins_300){
        std::cout << bin << ", ";
    }
    std::cout << std::endl;

    std::cout << "Print out bins for pair pT 8-500 GeV" << std::endl;
    for (auto bin : pT_bins_500){
        std::cout << bin << ", ";
    }
    std::cout << std::endl;


}


