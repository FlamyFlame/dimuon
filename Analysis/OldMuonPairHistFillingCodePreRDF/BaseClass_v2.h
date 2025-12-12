#pragma once
#include "ParamsSet.h"
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TH2D.h>
#include <nlohmann/json.hpp>
#include <map>
#include <vector>
#include <string>
#include <array>
#include <stdexcept>
#include <fstream>
#include <iostream>


using json = nlohmann::json;

enum class DataType { Float, Int, Bool, Double, Unknown };

class MuonPairHistFillingBaseClass {
protected:

// --------------------- input files ---------------------------

    // input files
    const std::vector<std::string>& inputFiles;
    // Tree structures
    std::array<TChain*, 2> m_chains;
    std::array<TTreeReader, 2> m_readers;
    
    // Branch names
    std::vector<std::string> m_branchNames_float = {"pair_dPoverP", "pt_lead", "pair_pt", "pair_eta", "pair_y", "asym", "dpt", "deta", "etaavg", "phiavg", "dphi", "dr", "minv", "m1.pt", "m2.pt", "m1.eta", "m2.eta", "m1.phi", "m2.phi", "m1.charge", "m2.charge"};
    std::vector<std::string> m_branchNames_int;
    std::vector<std::string> m_branchNames_bool;
    std::vector<std::string> m_branchNames_double;

    // Branch maps
    std::map<std::string, TTreeReaderValue<float>*> m_branchMap_float;
    std::map<std::string, TTreeReaderValue<int>*> m_branchMap_int;
    std::map<std::string, TTreeReaderValue<bool>*> m_branchMap_bool;
    std::map<std::string, TTreeReaderValue<double>*> m_branchMap_double;

// --------------------- output file & histograms ---------------------------

    TFile *outFile = nullptr;
    std::string output_file_path;

    // Histogram storage
    std::map<std::pair<std::string, std::string>, std::vector<TH1D*>> m_hists1D;
    std::map<std::pair<std::string, std::string>, std::vector<TH2D*>> m_hists2D;

    // Configuration
    json m_config;
    std::map<std::string, std::vector<std::string>> m_cutGroups;
    std::map<std::string, std::vector<size_t>> m_groupDims;

    std::map<std::string, std::string> m_var1DMap;
    std::map<std::string, std::pair<std::string, std::string>> m_var2DMap;
    
    std::map<std::string, DataType> m_datatypeMap;

    // Event processing
    int m_maxEvents = -1;
    size_t m_processedEvents = 0;

// --------------------- class methods ---------------------------

    // Helper functions
    DataType ParseType(const std::string& typeStr);
    void LoadConfig();
    std::vector<double> GetBins(const json& binSpec, const std::string& axis = "x");
    void InitHists();
    void CleanupBranches();

    virtual void InitInput();
    virtual void InitOutput();
    virtual void InitHists();
    virtual void ProcessData();
    virtual void WriteOutput();

public:
// --------------------- public class variables ---------------------------

    std::string m_jsonConfigFilePath;

// --------------------- public class methods ---------------------------

    MuonPairHistFillingBaseClass();
    virtual ~MuonPairHistFillingBaseClass();

    void setMaxEvents(int max) { m_maxEvents = max; }
    void process();
    void ProcessMuonPair(int sign) = 0;
    size_t CalculateFlatIndex(const std::vector<size_t>& cutIndices);

protected:
    template<typename T>
    T getValue(const std::string& branchName, int sign);
};
