#ifndef MuonPairHistFillingBaseClass_h
#define MuonPairHistFillingBaseClass_h

#include "ParamsSet.h"
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <string>
#include <map>
#include <stdexcept>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class MuonPairHistFillingBaseClass {
protected:
    TChain* m_chain1;
    TChain* m_chain2;
    std::vector<std::string> m_inputFiles;
    TFile* m_outputFile;
    std::map<std::string, std::vector<TH1D*>> m_map_hists1D;
    std::map<std::string, std::vector<TH2D*>> m_map_hists2D;
    json m_config;
    std::vector<std::string> m_cutOrder;
    std::vector<size_t> m_cutDimensions;

    // Map cut names to ParamsSet constants
    std::map<std::string, unsigned int> m_cutNameToNBins = {
        {"pt", ParamsSet::nPtBins},
        {"centrality", ParamsSet::nCtrBins},
        {"dphi", ParamsSet::ndphiRegions},
        {"deta", ParamsSet::ndetaRegions},
        {"gapcut", ParamsSet::nGapCuts},
        // Add other cuts as needed
    };

    void LoadConfig(const std::string& configPath);

    size_t CalculateFlatIndex(const std::vector<size_t>& cutIndices) {
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

    TH1D* CreateHist1D(const json& var) {
        const auto& binning = var["binning"];
        std::vector<double> bins;
        
        if (binning.is_string()) {
            std::string key = binning.get<std::string>();
            if (!ParamsSet::log_binning_map.count(key))
                throw std::runtime_error("Unknown binning key: " + key);
            bins = ParamsSet::log_binning_map[key];
        } else {
            int nbins = binning["nbins"];
            double xmin = binning["xmin"];
            double xmax = binning["xmax"];
            double step = (xmax - xmin)/nbins;
            for (int i=0; i<=nbins; ++i)
                bins.push_back(xmin + i*step);
        }
        
        return new TH1D("", var["title"].get<std::string>().c_str(),
                       bins.size()-1, bins.data());
    }

    TH2D* CreateHist2D(const json& var) {
        const auto& binning = var["binning"];
        
        // Process X axis
        std::vector<double> xbins;
        if (binning["x"].is_string()) {
            std::string key = binning["x"].get<std::string>();
            if (!ParamsSet::log_binning_map.count(key))
                throw std::runtime_error("Unknown x-binning key: " + key);
            xbins = ParamsSet::log_binning_map[key];
        } else {
            int nbins = binning["x"]["nbins"];
            double xmin = binning["x"]["xmin"];
            double xmax = binning["x"]["xmax"];
            double step = (xmax - xmin)/nbins;
            for (int i=0; i<=nbins; ++i)
                xbins.push_back(xmin + i*step);
        }

        // Process Y axis
        std::vector<double> ybins;
        if (binning["y"].is_string()) {
            std::string key = binning["y"].get<std::string>();
            if (!ParamsSet::log_binning_map.count(key))
                throw std::runtime_error("Unknown y-binning key: " + key);
            ybins = ParamsSet::log_binning_map[key];
        } else {
            int nbins = binning["y"]["nbins"];
            double ymin = binning["y"]["ymin"];
            double ymax = binning["y"]["ymax"];
            double step = (ymax - ymin)/nbins;
            for (int i=0; i<=nbins; ++i)
                ybins.push_back(ymin + i*step);
        }

        return new TH2D("", var["title"].get<std::string>().c_str(),
                       xbins.size()-1, xbins.data(),
                       ybins.size()-1, ybins.data());
    }

public:
    MuonPairHistFillingBaseClass() : m_chain1(nullptr), m_chain2(nullptr), m_outputFile(nullptr) {}
    virtual ~MuonPairHistFillingBaseClass() {
        for (auto& [name, vec] : m_map_hists1D) 
            for (auto h : vec) delete h;
        for (auto& [name, vec] : m_map_hists2D) 
            for (auto h : vec) delete h;
        if (m_outputFile) m_outputFile->Close();
    }

    virtual void InitInput(const std::vector<std::string>& inputFiles) {
        m_inputFiles = inputFiles;
        m_chain1 = new TChain("mytree_sign1");
        m_chain2 = new TChain("mytree_sign2");
        for (const auto& file : m_inputFiles) {
            m_chain1->Add(file.c_str());
            m_chain2->Add(file.c_str());
        }
    }

    virtual void InitOutput(const std::string& outputFile) {
        m_outputFile = new TFile(outputFile.c_str(), "RECREATE");
    }

    virtual void InitHists(const std::string& configPath) {
        LoadConfig(configPath);
        
        // Create 1D histograms
        size_t totalHists = 1;
        for (auto dim : m_cutDimensions) totalHists *= dim;

        for (const auto& var : m_config["variables1D"]) {
            std::string name = var["name"];
            TH1D* prototype = CreateHist1D(var);
            prototype->SetXTitle(var["xtitle"].get<std::string>().c_str());
            
            std::vector<TH1D*> hists;
            for (size_t i=0; i<totalHists; ++i) {
                TH1D* h = (TH1D*)prototype->Clone(Form("%s_%zu", name.c_str(), i));
                hists.push_back(h);
            }
            delete prototype;
            m_map_hists1D[name] = hists;
        }

        // Create 2D histograms
        for (const auto& var : m_config["variables2D"]) {
            std::string name = var["name"];
            TH2D* prototype = CreateHist2D(var);
            prototype->SetXTitle(var["xtitle"].get<std::string>().c_str());
            prototype->SetYTitle(var["ytitle"].get<std::string>().c_str());
            
            std::vector<TH2D*> hists;
            for (size_t i=0; i<totalHists; ++i) {
                TH2D* h = (TH2D*)prototype->Clone(Form("%s_%zu", name.c_str(), i));
                hists.push_back(h);
            }
            delete prototype;
            m_map_hists2D[name] = hists;
        }
    }

    virtual void ProcessData() {
        for (int sign : {1, 2}) {
            TChain* chain = (sign == 1) ? m_chain1 : m_chain2;
            Long64_t nEntries = chain->GetEntries();
            for (Long64_t i=0; i<nEntries; ++i) {
                chain->GetEntry(i);
                FillHistograms(sign);
            }
        }
    }

    virtual void WriteOutput() {
        m_outputFile->cd();
        for (auto& [name, vec] : m_map_hists1D)
            for (auto h : vec) h->Write();
        for (auto& [name, vec] : m_map_hists2D)
            for (auto h : vec) h->Write();
    }

    void FillHist1D(const std::string& name, double value, 
                   const std::vector<size_t>& cutIndices) {
        if (!m_map_hists1D.count(name)) return;
        size_t idx = CalculateFlatIndex(cutIndices);
        m_map_hists1D[name][idx]->Fill(value);
    }

    void FillHist2D(const std::string& name, double x, double y,
                   const std::vector<size_t>& cutIndices) {
        if (!m_map_hists2D.count(name)) return;
        size_t idx = CalculateFlatIndex(cutIndices);
        m_map_hists2D[name][idx]->Fill(x, y);
    }

    virtual void FillHistograms(int sign) = 0;
};

#endif