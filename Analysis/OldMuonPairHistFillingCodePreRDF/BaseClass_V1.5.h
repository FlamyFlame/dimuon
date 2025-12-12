#ifndef AnalysisBase_h
#define AnalysisBase_h

#include "ParamsSet.h"
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class AnalysisBase {
protected:
    // Data members
    TChain* m_chain1;
    TChain* m_chain2;
    std::vector<std::string> m_inputFiles;
    TFile* m_outputFile;
    std::map<std::pair<std::string, std::string>, std::vector<TH1D*>> m_map_hists1D;
    std::map<std::pair<std::string, std::string>, std::vector<TH2D*>> m_map_hists2D;
    json m_config;

    // Cut group metadata
    std::map<std::string, std::vector<std::string>> m_cutGroups;
    std::map<std::string, std::string> m_groupSuffixes;
    std::map<std::string, std::vector<size_t>> m_groupDimensions;

    // Helper: Resolve binning from ParamsSet or JSON
    std::vector<double> GetBins(const json& binSpec) {
        if (binSpec.is_string()) {
            std::string key = binSpec.get<std::string>();
            if (!ParamsSet::log_binning_map.count(key))
                throw std::runtime_error("Unknown binning key: " + key);
            return ParamsSet::log_binning_map[key];
        } else {
            // Uniform binning logic
        }
    }

    // Initialize cut groups
    void LoadCutGroups() {
        for (auto& [groupName, groupSpec] : m_config["cuts_to_hist_attributes_map"].items()) {
            m_cutGroups[groupName] = groupSpec["cuts"].get<std::vector<std::string>>();
            m_groupSuffixes[groupName] = groupSpec["suffix"].get<std::string>();
            
            // Calculate dimensions for this group
            std::vector<size_t> dims;
            for (const auto& cut : m_cutGroups[groupName]) {
                dims.push_back(m_config["cuts"][cut]["nBins"]);
            }
            m_groupDimensions[groupName] = dims;
        }
    }

public:
    virtual void InitHists(const std::string& configPath) {
        // Load JSON config
        std::ifstream configFile(configPath);
        configFile >> m_config;
        LoadCutGroups();

        // Create histograms for each variable + cut group
        for (const auto& var : m_config["variables1D"]) {
            std::string varName = var["name"];
            TH1D* prototype = CreateHist1D(var); // Implement similar to previous

            for (const auto& [groupName, _] : m_cutGroups) {
                auto& dims = m_groupDimensions[groupName];
                size_t totalHists = 1;
                for (auto d : dims) totalHists *= d;

                std::vector<TH1D*> hists;
                for (size_t i = 0; i < totalHists; ++i) {
                    TH1D* h = (TH1D*)prototype->Clone(
                        (varName + m_groupSuffixes[groupName] + "_" + std::to_string(i)).c_str()
                    );
                    hists.push_back(h);
                }
                m_map_hists1D[{varName, groupName}] = hists;
            }
            delete prototype;
        }
        // Repeat for 2D variables
    }

    void FillHist1D(const std::string& varName, const std::string& groupName, 
                   double value, const std::vector<size_t>& cutIndices) {
        auto key = std::make_pair(varName, groupName);
        if (!m_map_hists1D.count(key)) return;

        auto& dims = m_groupDimensions[groupName];
        size_t flatIdx = 0, stride = 1;
        for (int i = dims.size()-1; i >= 0; --i) {
            flatIdx += cutIndices[i] * stride;
            stride *= dims[i];
        }
        m_map_hists1D[key][flatIdx]->Fill(value);
    }

    // Similar for FillHist2D
};
#endif