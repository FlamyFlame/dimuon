#ifndef HIST_FILLING_BASE_CLASS_H
#define HIST_FILLING_BASE_CLASS_H

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THn.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <string>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <stdexcept>
#include <memory>

#include "ParamsSet.h"   // your existing global settings

/* ----------------------- simple POD ----------------------- */
struct var1D {
    std::string name{""};
    int         nbins{0};
    double      vmin{0.};
    double      vmax{0.};
    std::string title{""};
    std::vector<double> bins{};          // optional variable binning

    bool isValid() const {
        if (nbins <= 0) return false;
        if (!bins.empty()) return (static_cast<int>(bins.size()) == nbins + 1);
        return vmin < vmax;
    }
};

/* ===========================================================
 *                    HistFillingBaseClass
 * =========================================================*/
class HistFillingBaseClass {
protected:
    /* ---- I/O ---- */
    std::string in_out_dir_;
    std::string infile_name_;
    std::string outfile_name_;
    std::unique_ptr<TFile> outFile_;
    std::unique_ptr<TFile> inFile_;
    TTree* inTree_[ParamsSet::nSigns]{nullptr};

    /* ---- run-time readers ---- */
    std::unique_ptr<TTreeReader> reader_[ParamsSet::nSigns];

    /* ---- branch containers (one per TTree) ---- */
    /*   Example: float dphi[ParamsSet::nSigns];                    */
    /*   Add all variables you need here.                           */

    /* ---- 1D-variable catalogue ---- */
    std::map<std::string,var1D> var1D_list_;
    virtual void FillVar1DList() = 0;        // child defines variables

    /* ---- filter <-> variable lists ---- */
    std::vector<std::string> filter_list_signed_;
    std::vector<std::string> filter_list_unsigned_;
    std::map<std::string,std::string> filter2suffix_;
    std::map<std::string,std::vector<std::string>>            filter2var1Ds_;
    std::map<std::string,std::vector<std::pair<std::string,std::string>>> filter2var2Ds_;
    std::map<std::string,std::vector<std::vector<std::string>>>           filter2varnDs_;
    virtual void InitFilterMaps() = 0;       // child enumerates which vars to plot

    /* ---- histogram maps ---- */
    using H1Pair   = std::pair<std::string,TH1D*>;
    using H2Pair   = std::pair<std::pair<std::string,std::string>,TH2D*>;
    using HnPair   = std::pair<std::vector<std::string>,THn*>;
    using H1PairSS = std::pair<std::string,std::pair<TH1D*,TH1D*>>;
    using H2PairSS = std::pair<std::pair<std::string,std::string>,std::pair<TH2D*,TH2D*>>;
    using HnPairSS = std::pair<std::vector<std::string>,std::pair<THn*,THn*>>;

    std::map<std::string,std::vector<H1Pair>>   u_filt_h1_;
    std::map<std::string,std::vector<H2Pair>>   u_filt_h2_;
    std::map<std::string,std::vector<HnPair>>   u_filt_hn_;

    std::map<std::string,std::vector<H1PairSS>> s_filt_h1_;
    std::map<std::string,std::vector<H2PairSS>> s_filt_h2_;
    std::map<std::string,std::vector<HnPairSS>> s_filt_hn_;

    /* ---- var-alias → branch-address map ---- */
    std::map<std::string,float*> histvar2branch_;

    /* ---------- helpers ---------- */
    static std::string DefaultSuffixMapper(const std::string& f){
        if (f=="DEFAULT") return "";
        return "_"+f;
    }

    /* ROOT object builders */
    TH1D* MakeHist1D(const var1D& var,const std::string& suffix,int hist_sign=-1);
    TH2D* MakeHist2D(const var1D& vx,const var1D& vy,const std::string& suffix,int hist_sign=-1);
    THn*  MakeHistnD(const std::vector<var1D>& vlist,const std::string& suffix,int hist_sign=-1);

    /* --------- obligatory workflow hooks --------- */
    virtual void InitInput();          // open file & set branches
    virtual void InitOutput();         // open output ROOT file
    virtual void InitHists();          // auto-define all hists
    virtual void ProcessData();        // loop over trees & call FillHistograms
    virtual void WriteOutput();        // save & close
    virtual bool PassSingleMuonGapCut(float eta,float pt,int charge);
    virtual void FillHistograms(int nsign);    // TOP entry for one sign
    virtual void FillHistogramsPerFilter(const std::string& filter,int nsign,bool filter_signed=true);

public:
    virtual ~HistFillingBaseClass() = default;
    virtual void Run();                // public entry-point
    /* ---------- user help ---------- */
    static void PrintInstructions();
};

#endif
