#include "MuonPairHistFillingBaseClass.h"
#include <iostream>

MuonPairHistFillingBaseClass::MuonPairHistFillingBaseClass() {}

MuonPairHistFillingBaseClass::~MuonPairHistFillingBaseClass() {
    cleanupBranches();
    for(auto& h : m_hists1D) for(auto* hist : h.second) delete hist;
    for(auto& h : m_hists2D) for(auto* hist : h.second) delete hist;
    delete m_chains[0];
    delete m_chains[1];
}

size_t MuonPairHistFillingBaseClass::CalculateFlatIndex(const std::vector<size_t>& cutIndices) {
    if (cutIndices.size() != m_cutOrder.size())
        throw std::runtime_error("Cut indices size mismatch");
    
    size_t index = 0;
    size_t stride = 1;
    for (int i = m_cutOrder.size()-1; i >= 0; --i) {
        if (cutIndices[i] >= m_cutDimensions[i])
            throw std::runtime_error("Cut index out of range");
        index += cutIndices[i] * stride;
        stride *= m_cutDimensions[i];
    }
    return index;
}

DataType MuonPairHistFillingBaseClass::ParseType(const std::string& typeStr) {
    if(typeStr == "float") return DataType::Float;
    if(typeStr == "int") return DataType::Int;
    if(typeStr == "bool") return DataType::Bool;
    if(typeStr == "double") return DataType::Double;
    return DataType::Unknown;
}

void MuonPairHistFillingBaseClass::LoadConfig() {
    std::ifstream f(m_jsonConfigFilePath);
    
    if (!f) {  // Check if the file exists and can be opened
        throw std::runtime_error("Error: Cannot open JSON config file: " + m_jsonConfigFilePath);
    }

    try {
        f >> m_config;  // Parse the JSON file
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error("JSON parsing error in file " + m_jsonConfigFilePath + ": " + std::string(e.what()));
    }

    // Load cut groups
    for(auto& [groupName, groupSpec] : m_config["cuts_combos"].items()) {
        m_cutGroups[groupName] = groupSpec["cuts"].get<std::vector<std::string>>();
        std::vector<size_t> dims;
        for(auto& cut : m_cutGroups[groupName]) {
            dims.push_back(m_config["cuts"][cut]["nBins"]);
        }
        m_groupDims[groupName] = dims;
    }

    // Load variables
    for(auto& var : m_config["variables1D"]) {
        std::string hist_name = var["hist_name"];
        std::string var_name = var["var_name"];
        DataType dtype = var.contains("type") ? ParseType(var["type"]) : DataType::Float;
        m_var1DMap[hist_name] = var_name;
        m_datatypeMap[var_name] = dtype;
    }

    for(auto& var : m_config["variables2D"]) {
        std::string hist_name = var["hist_name"];
        std::string xvar = var["xvar"];
        std::string yvar = var["yvar"];
        DataType xdtype = var.contains("typex") ? ParseType(var["typex"]) : DataType::Float;
        DataType ydtype = var.contains("typey") ? ParseType(var["typey"]) : DataType::Float;
        m_var2DMap[hist_name] = {xvar, yvar};
        m_datatypeMap[xvar] = xdtype;
        m_datatypeMap[yvar] = ydtype;
    }
}

void MuonPairHistFillingBaseClass::InitInput() {
    // Initialize chains
    for(int i : {0, 1}) {
        m_chains[i] = new TChain(Form("mytree_sign%d", i+1));
        for(const auto& file : inputFiles) {
            m_chains[i]->Add(file.c_str());
        }
        m_readers[i].SetTree(m_chains[i]);
    }

    // Set branches
    for(int sign : {0, 1}) {
        auto& reader = m_readers[sign];
        for (std::string branch_name : m_branchNames_float){
            m_branchMap_float[branch_name] = new TTreeReaderValue<float>(reader, branch_name.c_str());
        }
        for (std::string branch_name : m_branchNames_int){
            m_branchMap_int[branch_name] = new TTreeReaderValue<int>(reader, branch_name.c_str());
        }
        for (std::string branch_name : m_branchNames_bool){
            m_branchMap_bool[branch_name] = new TTreeReaderValue<bool>(reader, branch_name.c_str());
        }
        for (std::string branch_name : m_branchNames_double){
            m_branchMap_double[branch_name] = new TTreeReaderValue<double>(reader, branch_name.c_str());
        }
    }
}

std::vector<double> MuonPairHistFillingBaseClass::GetBins(const json& binSpec, const std::string& axis = "x") {
    if (binSpec.is_string()) {
        std::string key = binSpec.get<std::string>();
        if (!ParamsSet::log_binning_map.count(key))
            throw std::runtime_error("Unknown binning key: " + key);
        return ParamsSet::log_binning_map[key];
    } else {
        int nbins = binSpec["nbins"];
        std::string minKey = axis + "min";
        std::string maxKey = axis + "max";
        double minVal = binSpec[minKey];
        double maxVal = binSpec[maxKey];
        std::vector<double> bins;
        double step = (maxVal - minVal)/nbins;
        for (int i=0; i<=nbins; ++i)
            bins.push_back(minVal + i*step);
        return bins;
    }
}

void MuonPairHistFillingBaseClass::InitHists() {
    LoadConfig();

    // 1D Histograms
    for (const auto& var : m_config["variables1D"]) {
        const std::string name = var["name"];
        const std::string title = var.value("title", "");
        const std::string xtitle = var.value("xtitle", "");
        
        // Resolve binning
        std::vector<double> bins;
        try {
            bins = GetBins(var["binning"], "x");
        } catch (const std::exception& e) {
            std::cerr << "Error in " << name << " binning: " << e.what() << "\n";
            continue;
        }

        // Create histograms for each cut group
        for (const auto& [groupName, groupSpec] : m_config["cuts_combos"].items()) {
            const size_t totalHists = std::accumulate(
                m_groupDims[groupName].begin(), 
                m_groupDims[groupName].end(), 
                1, std::multiplies<>()
            );
            
            std::vector<TH1D*> hists;
            for (size_t i=0; i<totalHists; ++i) {
                TH1D* h = new TH1D(
                    Form("%s_%s_%zu", name.c_str(), groupSpec["suffix"].get<std::string>().c_str(), i),
                    title.c_str(),
                    bins.size()-1, bins.data()
                );
                h->SetXTitle(xtitle.c_str());
                hists.push_back(h);
            }
            m_hists1D[{name, groupName}] = hists;
        }
    }

    // 2D Histograms
    for (const auto& var : m_config["variables2D"]) {
        const std::string name = var["name"];
        const std::string title = var.value("title", "");
        const std::string xtitle = var.value("xtitle", "");
        const std::string ytitle = var.value("ytitle", "");
        
        // Resolve X and Y binnings
        std::vector<double> xBins, yBins;
        try {
            xBins = GetBins(var["binning"]["x"], "x");
            yBins = GetBins(var["binning"]["y"], "y");
        } catch (const std::exception& e) {
            std::cerr << "Error in " << name << " binning: " << e.what() << "\n";
            continue;
        }

        // Create histograms for each cut group
        for (const auto& [groupName, groupSpec] : m_config["cuts_combos"].items()) {
            const size_t totalHists = std::accumulate(
                m_groupDims[groupName].begin(),
                m_groupDims[groupName].end(),
                1, std::multiplies<>()
            );
            
            std::vector<TH2D*> hists;
            for (size_t i=0; i<totalHists; ++i) {
                TH2D* h = new TH2D(
                    Form("%s_%s_%zu", name.c_str(), groupSpec["suffix"].get<std::string>().c_str(), i),
                    title.c_str(),
                    xBins.size()-1, xBins.data(),
                    yBins.size()-1, yBins.data()
                );
                h->SetXTitle(xtitle.c_str());
                h->SetYTitle(ytitle.c_str());
                hists.push_back(h);
            }
            m_hists2D[{name, groupName}] = hists;
        }
    }
}

template<typename T>
T MuonPairHistFillingBaseClass::getValue(const std::string& branchName, int sign) {
    try {
        if constexpr(std::is_same_v<T, float>) {
            return **m_branchMap_float.at(branchName);
        }
        else if constexpr(std::is_same_v<T, int>) {
            return **m_branchMap_int.at(branchName);
        }
        // Handle other types
    }
    catch(const std::out_of_range&) {
            std::cerr << "WARNING: Missing branch " << branchName 
                      << " for type " << typeid(T).name() << "\n";
            return T{}; // Return default value
    }
}

void MuonPairHistFillingBaseClass::process() {
    for(int sign : {0, 1}) {
        m_readers[sign].Restart();
        while(m_readers[sign].Next()) {
            if(m_maxEvents > 0 && m_processedEvents >= m_maxEvents) break;
            
            try {
                ProcessMuonPair(sign);
                // auto cutIndices = CalculateFlatIndex(????????);
                // Histogram filling logic here
                m_processedEvents++;
            }
            catch(const std::exception& e) {
                std::cerr << "Error processing event: " << e.what() << std::endl;
            }
        }
    }
}

void MuonPairHistFillingBaseClass::CleanupBranches() {
    for (auto& [k,v] : m_branchMap_float) delete v;
    for (auto& [k,v] : m_branchMap_int) delete v;
    for (auto& [k,v] : m_branchMap_bool) delete v;
    for (auto& [k,v] : m_branchMap_double) delete v;
}

void MuonPairHistFillingBaseClass::FillSingleHistogram(std::string branchName, TH1* hist){
    auto it = myMap.find(branchName);
    if (it == myMap.end()) { 
        std::cerr << "Warning: variable '" << branchName << "' not found in m_datatypeMap. Return without filling histogram.\n";
        return;
    }

    auto dtype = it->second;
    switch(dtype) {
        case DataType::Float:
            return **m_branchMap_float.at(branchName);
        case DataType::Int:
            return **m_branchMap_int.at(branchName);
        case DataType::Bool:
            return **m_branchMap_bool.at(branchName);
        case DataType::Double:
            return **m_branchMap_default.at(branchName);
    }

}
