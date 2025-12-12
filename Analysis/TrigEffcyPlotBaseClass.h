#pragma once
// TrigEffcyPlotBaseClass.h (or .cxx if you keep it ROOT-macro style)

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <array>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "MuonObjectsParamsAndHelpers/PbPbBaseClass.h"

// ------------------------------------------------------------
// Shared tiny helpers (keep as-is from your current files)
// ------------------------------------------------------------
template<typename T>
static void SetStyle(T* g, Color_t col, Style_t m, Width_t w = 2, Style_t ls = 1)
{
    if (!g) return;
    g->SetMarkerColor(col);
    g->SetLineColor(col);
    g->SetMarkerStyle(m);
    g->SetLineWidth(w);
    g->SetLineStyle(ls);
}

static bool fileExists(const std::string& dir)
{
    return (gSystem->AccessPathName(dir.c_str()) == kFALSE);
}

static void makeDirIfNeeded(const std::string& dir)
{
    if (!fileExists(dir)) gSystem->mkdir(dir.c_str(), kTRUE /*recursive*/);
}

// ============================================================
// Base class
// ============================================================
class TrigEffcyPlotBaseClass
{
public:
    bool only_draw_single_muon_effcy = false; // skip pair observables
    bool draw_single_b_signal = false;
    bool use_sepr_for_op_and_signal = false;
    bool debug_mode = false;

    enum class SignalDrawingMode { Signed, Signal, OpAndSignal };

    using StrVec     = std::vector<std::string>;
    using XYLog      = std::map<std::string, std::pair<bool, bool>>;
    using ModeMap    = std::map<SignalDrawingMode, bool>;
    using TriggerMap = std::map<std::string, std::string>; // trigger → mapped-name
    using Rect       = std::array<double, 4>;

    virtual ~TrigEffcyPlotBaseClass()
    {
        if (fFile_mu4) { fFile_mu4->Close(); delete fFile_mu4; fFile_mu4 = nullptr; }
        if (fFile_MB)  { fFile_MB->Close();  delete fFile_MB;  fFile_MB  = nullptr; }
    }

    // Non-virtual public driver: common flow, but delegates *actual plotting*.
    void Run()
    {
        initialize();

        if (!fFile_mu4 || fFile_mu4->IsZombie()) return;
        if (isMB && (!fFile_MB || fFile_MB->IsZombie())) return;

        gStyle->SetOptStat(0);

        // derived class decides triggers (PP vs PbPb logic differs)
        trg_list = buildTriggerList();

        // Pair-observable efficiencies (skip in MB mode)
        if (!isMB && !only_draw_single_muon_effcy)
        {
            for (const auto& v : determineVar1DList()) plot1D(v);

            // 2D efficiencies
            {
                StrVec fList = filter_suffix_list.empty() ? StrVec{""} : filter_suffix_list;
                for (const auto& fs : fList)
                {
                    if (!filt_to_draw2D_map.at(fs)) continue;
                    for (const auto& v : determineVar2DList(fs)) plot2D(v, fs);
                }
            }

            // 2D profiles
            {
                StrVec fList = filter_suffix_list.empty() ? StrVec{""} : filter_suffix_list;
                for (const auto& fs : fList)
                {
                    if (!filt_to_draw2Dprofile_map.at(fs)) continue;
                    for (const auto& v : determineVar2DPList(fs)) plotProfile(v, fs);
                }
            }

            plot2Dto1DEffcyProj();
        }

        // Single-muon efficiencies (same gating logic you already use)
        if (filter_suffix_list.empty() || isMB)
        {
            for (const auto& v : vars_2D_2nd_muon) plot2D_SingleMuonEffcy(v);
            for (const auto& kv : var2Ds_for1DsingleMuonEffcy_proj) plot2Dto1DsingleMuonEffcyProj(kv);
        }
    }

protected:
    // -------------------------
    // Construction pattern (virtual class: not called as public function)
    // -------------------------
    TrigEffcyPlotBaseClass(int runYear,
                           const StrVec& vars1D_in,
                           bool draw2mu4_user,
                           const StrVec& filter_suffixes,
                           const XYLog& logopt_,
                           const StrVec& vars2D_in,
                           const StrVec& vars2DProf_in,
                           const std::map<std::string,Rect>& legSigned,
                           const std::map<std::string,Rect>& legSignal,
                           const std::map<std::string,Rect>& legOpSig,
                           const std::map<std::string,std::pair<bool,bool>>& var2DProj,
                           const std::map<std::string,std::pair<bool,bool>>& var2DSingleMuonEffcyProj)
        : runYear_(runYear),
          draw2mu4(draw2mu4_user),
          fFile_mu4(nullptr),
          fFile_MB(nullptr),
          var1Ds(vars1D_in),
          var2Ds(vars2D_in),
          var2Ds_for_profile(vars2DProf_in),
          logopt(logopt_),
          filter_suffix_list(filter_suffixes)
    {
        legendPosSigned   = legSigned;
        legendPosSignal   = legSignal;
        legendPosOpSignal = legOpSig;
        var2Ds_for1Deffcy_proj  = var2DProj;
        var2Ds_for1DsingleMuonEffcy_proj  = var2DSingleMuonEffcyProj;

        isMB = (std::find(filter_suffix_list.begin(), filter_suffix_list.end(),
                          std::string("_MB")) != filter_suffix_list.end());
    }

    // ============================================================
    // Pure-virtual hooks: you said these must differ fully
    // ============================================================
    virtual void plot1D(const std::string& var) = 0;
    virtual void plot2D(const std::string& var, const std::string& fs) = 0;
    virtual void plotProfile(const std::string& var, const std::string& fs) = 0;
    virtual void plot2Dto1DEffcyProj() = 0;
    virtual void plot2D_SingleMuonEffcy(const std::string& var) = 0;
    virtual void plot2Dto1DsingleMuonEffcyProj(const std::pair<const std::string,std::pair<bool,bool>>& kv) = 0;

    // Run-year specific setup (paths + filenames)
    virtual void configureDataFiles(int runYear) = 0;

    // Trigger list logic differs (PP has Run-2 special casing; PbPb doesn’t)
    virtual StrVec buildTriggerList() const = 0;

    // ============================================================
    // Shared data members
    // ============================================================
    int runYear_{0};
    std::string data_dir;
    std::string fname_single_mu4;
    std::string fname_MB;

    bool draw2mu4{};
    bool isMB{};

    StrVec trg_list;

    StrVec filters_weighted = {"_inv_w_by_single_mu_effcy", "_excl_inv_w_by_single_mu_effcy"};
    StrVec vars_2D_2nd_muon  = {"pt2nd_vs_q_eta2nd", "pt2nd_vs_phi2nd", "phi2nd_vs_q_eta2nd"};

    std::map<std::string, std::string> target_dimuon_trigger_to_single_muon_map =
        {{"mu4_mu4noL1", "mu4noL1"}, {"2mu4", "mu4"}};

    TFile* fFile_mu4{nullptr};
    TFile* fFile_MB{nullptr};

    // variable lists
    StrVec var1Ds;
    StrVec var2Ds;
    StrVec var2Ds_for_profile;

    // projections
    std::map<std::string, std::pair<bool,bool>> var2Ds_for1Deffcy_proj;
    std::map<std::string, std::pair<bool,bool>> var2Ds_for1DsingleMuonEffcy_proj;

    // user options
    XYLog  logopt;
    StrVec filter_suffix_list;

    // legend rectangles
    std::map<std::string,Rect> legendPosSigned;
    std::map<std::string,Rect> legendPosSignal;
    std::map<std::string,Rect> legendPosOpSignal;

    std::map<SignalDrawingMode, std::string> mode_to_png_suffix = {
        {SignalDrawingMode::Signed,      ""},
        {SignalDrawingMode::Signal,      "_w_sig_sel"},
        {SignalDrawingMode::OpAndSignal, "_all_op_and_sig_sel"}
    };

    // filter infrastructure
    mutable std::map<std::string, TriggerMap> filt_to_trig_map;
    std::map<std::string, TriggerMap>         filt_to_trig_exception_map;
    std::map<std::string, ModeMap>            filt_to_mode_map;
    std::map<std::string, ModeMap>            filt_to_mode_exception_map;
    std::map<std::string, StrVec>             filt_to_var1D_exception_map;

    std::map<std::string,bool>  filt_to_draw2D_map;
    std::map<std::string,bool>  filt_to_draw2Dprofile_map;
    std::map<std::string,StrVec> filt_to_var2D_exception_map;
    std::map<std::string,StrVec> filt_to_var2Dprofile_exception_map;

    std::map<std::string, std::string> filt_suffix_to_label_map;

    const StrVec kTriggers = {"mu4", "mu4_mu4noL1", "2mu4"};

    // ============================================================
    // Shared helpers (copied once; both your files already share these)
    // ============================================================
    virtual void initialize(){
        // derived sets filenames/dirs (pp vs pbpb)
        configureDataFiles(runYear_);

        buildFilterLabelMaps();
        buildGlobalMaps();
        consistencyChecks();
        initializeFile();
    }

    virtual void initializeFile()
    {
        fFile_mu4 = TFile::Open(fname_single_mu4.c_str(), "READ");
        if (!fFile_mu4 || fFile_mu4->IsZombie())
            std::cerr << "[TrigEffcyPlotBaseClass] ERROR: cannot open " << fname_single_mu4 << "\n";

        if (isMB)
        {
            fFile_MB = TFile::Open(fname_MB.c_str(), "READ");
            if (!fFile_MB || fFile_MB->IsZombie())
                std::cerr << "[TrigEffcyPlotBaseClass] ERROR: cannot open " << fname_MB << "\n";
        }
    }

    void buildFilterLabelMaps()
    {
        filt_suffix_to_label_map["_good_accept"]               = "good accept";
        filt_suffix_to_label_map["_inv_w_by_single_mu_effcy"]   = "1/(single muon effcy) weight";
        filt_suffix_to_label_map["_sepr"]                      = "#Delta R > 0.8";
        filt_suffix_to_label_map["_w_single_b_sig_sel"]         = "Single-b signal";
    }

    void buildGlobalMaps()
    {
        // default trigger mapping: trg -> trg + suffix
        for (const auto& fs : filter_suffix_list)
        {
            TriggerMap tmap;
            for (const auto& trg : kTriggers) tmap[trg] = trg + fs;
            tmap["mu4_mu4noL1_denom"] = "NONE";
            tmap["2mu4_denom"]        = "NONE";
            filt_to_trig_map[fs] = tmap;
        }

        // default "" entry
        {
            TriggerMap tmap;
            for (const auto& trg : kTriggers) tmap[trg] = trg;
            tmap["mu4_mu4noL1_denom"] = "NONE";
            tmap["2mu4_denom"]        = "NONE";
            filt_to_trig_map[""] = tmap;
        }

        // single-b single-muon efficiency
        {
            TriggerMap tmap;

            for (const auto& trg : kTriggers) tmap[trg] = trg + "_w_single_b_sig_sel";
            tmap["mu4_mu4noL1_denom"] = "NONE"; // no mu4_mu4noL1-specific denominator by default
            tmap["2mu4_denom"] = "NONE"; // no mu4_mu4noL1-specific denominator by default

            filt_to_trig_map["_w_single_b_sig_sel"] = tmap; // special key for default
        }

        // drawing modes default
        for (const auto& fs : filter_suffix_list)
        {
            filt_to_mode_map[fs] = {
                {SignalDrawingMode::Signed, true},
                {SignalDrawingMode::Signal, draw_single_b_signal},
                {SignalDrawingMode::OpAndSignal, draw_single_b_signal}
            };
        }
        filt_to_mode_map[""] = {
            {SignalDrawingMode::Signed, true},
            {SignalDrawingMode::Signal, draw_single_b_signal},
            {SignalDrawingMode::OpAndSignal, draw_single_b_signal}
        };

        // -------- _MB exceptions --------
        // mu4 → MB, mu4_mu4noL1 → (skip), 2mu4 → mu4
        filt_to_trig_exception_map["_MB"] = {
            {"mu4", "MB"},
            {"mu4_mu4noL1", "NONE"},
            {"2mu4", "mu4"}
        };

        // 2D flags: on only for default unless exceptions set
        filt_to_draw2D_map[""]        = true;
        filt_to_draw2Dprofile_map[""] = true;
        for (const auto& fs : filter_suffix_list)
        {
            filt_to_draw2D_map[fs]        = false;
            filt_to_draw2Dprofile_map[fs] = false;
        }

        // ---------------------------- _sepr exceptions -----------------------------------------------------
        // Only signed‑mode relevant for _sepr
        filt_to_mode_exception_map["_sepr"] = {
            {SignalDrawingMode::Signal, false},
        };

        filt_to_draw2D_map["_sepr"]              = true;
        filt_to_var2D_exception_map["_sepr"]     = {"pair_eta_vs_pair_pT", "Deta_Dphi", "eta1_eta2", "eta_avg_Deta", "eta_avg_Dphi", "minv_pair_pt_log"};

        // ---------------------------- _excl exceptions ---------------------------------------------------
        filt_to_trig_exception_map["_excl"] = {
            {"2mu4", "NONE"} // "exclusiveness" non-existent for 2mu4
        };
        // keep default modes

        // ---------------------------- _inv_w_by_single_mu_effcy exceptions ------------------------------
        filt_to_trig_exception_map["_inv_w_by_single_mu_effcy"] = {
            {"mu4", "NONE"},
            {"mu4_mu4noL1_denom", "mu4_mu4noL1_inv_w_by_single_mu_effcy_denom"}, // uses separate denominators for mu4_mu4noL1 & 2mu4 weighted efficiencies
            {"2mu4_denom", "2mu4_inv_w_by_single_mu_effcy_denom"}  // uses separate denominators for mu4_mu4noL1 & 2mu4 weighted efficiencies
        };

        filt_to_trig_exception_map["_excl_inv_w_by_single_mu_effcy"] = {
            {"mu4", "NONE"}, // means: leave mu4 unchanged
            {"mu4_mu4noL1_denom", "mu4_mu4noL1_excl_inv_w_by_single_mu_effcy_denom"}, // uses separate denominators for mu4_mu4noL1 & 2mu4 weighted efficiencies
            {"2mu4", "NONE"} // "exclusiveness" non-existent for 2mu4
        };

        filt_to_mode_exception_map["_inv_w_by_single_mu_effcy"] = {
            {SignalDrawingMode::OpAndSignal, true}
        };
        
        filt_to_var1D_exception_map["_inv_w_by_single_mu_effcy"] = {"Deta_zoomin", "Dphi_zoomin", "DR_zoomin", "DR_0_2", "minv_zoomin", "pair_pt_log", "pt2nd"};

        // filters for weighted efficiencies should not be used to draw 2D variables involving only 2nd-muon kinematics (single-muon efficiencies)
        // apply before any other 2D variable exception to be added
        for (auto fs : filters_weighted){
            filt_to_var2D_exception_map[fs] = vars_2D_2nd_muon;
        }

        // ---------------------------- good/bad accept special per‑variable list ----------------------
        // Only signed‑mode relevant for _good_accept
        filt_to_mode_exception_map["_good_accept"] = {
            {SignalDrawingMode::Signal, false},
        };

        // Only signed‑mode relevant for _bad_accept
        filt_to_mode_exception_map["_bad_accept"] = {
            {SignalDrawingMode::Signal, false},
        };

        // merge exceptions
        for (const auto& [fs, exc] : filt_to_trig_exception_map)
        {
            auto& base = filt_to_trig_map[fs];
            for (const auto& [trg, newname] : exc) base[trg] = newname;
        }
        for (const auto& [fs, exc] : filt_to_mode_exception_map)
        {
            auto& base = filt_to_mode_map[fs];
            for (const auto& [md, val] : exc) base[md] = val;
        }
    }

    void consistencyChecks() const
    {
        if (filter_suffix_list.size() > 2) {
            std::cerr << "[TrigEffPlotter] WARNING: more than two filters requested – plots may become cluttered.\n";
        }
        // Basic agreement checks across filters (vars list / drawing modes)
        if (filter_suffix_list.size() > 1) {
            // all filters that specify var‑exception must agree on list
            const auto& fs0 = filter_suffix_list.front();
            if (filt_to_var1D_exception_map.count(fs0)) {
                const auto& ref = filt_to_var1D_exception_map.at(fs0);
                for (size_t i = 1; i < filter_suffix_list.size(); ++i) {
                    const auto& fsi = filter_suffix_list[i];
                    if (!filt_to_var1D_exception_map.count(fsi) ||
                        filt_to_var1D_exception_map.at(fsi) != ref) {
                        throw std::runtime_error("Conflicting var1D exception lists between filters");
                    }
                }
            }
            // drawing‑mode agreement (signed & signal should match across filters)
            auto drawing_modes = draw_single_b_signal? std::vector<SignalDrawingMode>{SignalDrawingMode::Signed, SignalDrawingMode::Signal} : std::vector<SignalDrawingMode>{SignalDrawingMode::Signed};

            for (auto mode : drawing_modes) {
                bool val0 = filt_to_mode_map.at(fs0).at(mode);
                for (size_t i = 1; i < filter_suffix_list.size(); ++i) {
                    bool vali = filt_to_mode_map.at(filter_suffix_list[i]).at(mode);
                    if (vali != val0)
                        throw std::runtime_error("Filters disagree on required drawing modes");
                }
            }
        }
    }

    StrVec determineVar1DList() const
    {
        if (filter_suffix_list.empty()) return var1Ds;
        const auto& fs0 = filter_suffix_list.front();
        auto it = filt_to_var1D_exception_map.find(fs0);
        if (it == filt_to_var1D_exception_map.end()) return var1Ds;
        return it->second;
    }

    StrVec determineVar2DList(const std::string& fs) const
    {
        auto it = filt_to_var2D_exception_map.find(fs);
        return (it == filt_to_var2D_exception_map.end()) ? var2Ds : it->second;
    }

    StrVec determineVar2DPList(const std::string& fs) const
    {
        auto it = filt_to_var2Dprofile_exception_map.find(fs);
        return (it == filt_to_var2Dprofile_exception_map.end()) ? var2Ds_for_profile : it->second;
    }

    // =========================================================================
    //  Histogram access helper
    // =========================================================================
    template<typename T> T* getHist(const std::string& name) const
    {
        T* h = dynamic_cast<T*>(fFile_mu4->Get(name.c_str()));
        if (!h) std::cerr << "[TrigEffcyPlotBaseClass] WARNING: missing hist " << name << "\n";
        return h;
    }

    template<typename T> T* getHistMB(const std::string& name) const
    {
        T* h = dynamic_cast<T*>(fFile_MB->Get(name.c_str()));
        if (!h) std::cerr << "[TrigEffcyPlotBaseClass] WARNING: in MB file, missing hist " << name << "\n";
        return h;
    }

    // =========================================================================
    //  Name helper (after filter mapping)
    // =========================================================================
    std::string mappedTrigger(const std::string& filter_suffix, const std::string& trigger) const
    {
        if (filt_to_trig_map.find(filter_suffix) == filt_to_trig_map.end())
        {
            TriggerMap tmap;
            for (const auto& trg : kTriggers) tmap[trg] = trg + filter_suffix;
            tmap["mu4_mu4noL1_denom"] = "NONE";
            tmap["2mu4_denom"]        = "NONE";
            filt_to_trig_map[filter_suffix] = tmap;
        }
        const auto& m = filt_to_trig_map.at(filter_suffix);
        auto it = m.find(trigger);
        if (it == m.end()) return "";
        return it->second;
    }

    // =========================================================================
    // ---------- build histogram name helpers --------------------------
    // =========================================================================
    std::string hName(const std::string& var,
                      const std::string& fs,
                      const std::string& trg,
                      const std::string& suffix) const
    {
        std::string mapped = mappedTrigger(fs, trg);
        return (mapped.empty() || mapped == "NONE") ? "" :
               Form("h_%s%s_%s", var.c_str(), suffix.c_str(), mapped.c_str());
    }

    // =========================================================================
    // Ensure a TGraph drawn on a log-x pad spans the full histogram range
    // =========================================================================
    void adjustLogXRange(TGraphAsymmErrors* g, const TH1* hRef) const
    {
        if (!g || !hRef) return;
        double xmin = hRef->GetXaxis()->GetBinLowEdge(1);
        double xmax = hRef->GetXaxis()->GetBinLowEdge(hRef->GetNbinsX() + 1);
        g->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    // =========================================================================
    // Split "Y_vs_X" → {"Y","X"}; if not matching, returns {var,""}
    // =========================================================================
    std::pair<std::string,std::string> parseYvsX(const std::string& var2d) const
    {
        auto pos = var2d.find("_vs_");
        if (pos == std::string::npos) return {var2d, ""};
        return { var2d.substr(0, pos), var2d.substr(pos + 4) };
    }

    // =========================================================================
    // re-use legend placement maps with a 1D variable key (X or Y)
    // =========================================================================
    Rect legendBoxFor1DVar(const std::string& v1d,
                           SignalDrawingMode mode,
                           const StrVec& filters) const
    {
        Rect def = (mode == SignalDrawingMode::Signed)
                   ? Rect{0.43,0.12,0.80,0.35}
                   : Rect{0.50,0.14,0.77,0.30};

        if      (mode == SignalDrawingMode::Signed      && legendPosSigned.count(v1d))   return legendPosSigned.at(v1d);
        else if (mode == SignalDrawingMode::Signal      && legendPosSignal.count(v1d))   return legendPosSignal.at(v1d);
        else if (mode == SignalDrawingMode::OpAndSignal && legendPosOpSignal.count(v1d)) return legendPosOpSignal.at(v1d);

        // fallback: grow left if long filter label
        auto autoBox = [&](Rect r, const StrVec& fsList)->Rect {
            size_t L = 0;
            for (auto& fs : fsList) if (!fs.empty())
            {
                std::string lab = (filt_suffix_to_label_map.count(fs) ? filt_suffix_to_label_map.at(fs)
                                                                      : (fs[0]=='_' ? fs.substr(1) : fs));
                L = std::max(L, lab.size());
            }
            if (L > 4) { double extra = std::min(0.01*(L-4), 0.20); r[0] = std::max(0.05, r[0] - extra); }
            return r;
        };
        return autoBox(def, filters);
    }

    // You already have this (not shown in the snippet above): keep it shared here.
    // virtual/override is *not* needed; it’s just a utility.
    std::pair<bool,bool> xyFor(const std::string& var) const
    {
        auto it = logopt.find(var);
        if (it == logopt.end()) return {false,false};
        return it->second;
    }
};

// ============================================================
// Derived: PP
// ============================================================
class TrigEffPlotterPP : public TrigEffcyPlotBaseClass
{
public:
    TrigEffPlotterPP(int runYear,
                     const StrVec& vars1D_in,
                     bool draw2mu4_user = false,
                     const StrVec& filter_suffixes = {},
                     const XYLog&  logopt_ = {},
                     const StrVec& vars2D_in = {},
                     const StrVec& vars2DProf_in = {},
                     const std::map<std::string,Rect>& legSigned  = {},
                     const std::map<std::string,Rect>& legSignal  = {},
                     const std::map<std::string,Rect>& legOpSig   = {},
                     const std::map<std::string,std::pair<bool,bool>>& var2DProj = {},
                     const std::map<std::string,std::pair<bool,bool>>& var2DSingleMuonEffcyProj = {})
        : TrigEffcyPlotBaseClass(runYear, vars1D_in,
                                 draw2mu4_user,
                                 filter_suffixes, logopt_,
                                 vars2D_in, vars2DProf_in,
                                 legSigned, legSignal, legOpSig,
                                 var2DProj, var2DSingleMuonEffcyProj)
    {}

private:
    bool isRun2pp{false};

    void configureDataFiles(int runYear) override;

    StrVec buildTriggerList() const override
    {
        if (isRun2pp) return {"2mu4"};
        if (draw2mu4)  return {"2mu4", "mu4_mu4noL1"};
        return {"mu4_mu4noL1"};
    }

    // plotting functions implemented in PP file
    void plot1D(const std::string& var) override;
    void plot2D(const std::string& var, const std::string& fs) override;
    void plotProfile(const std::string& var, const std::string& fs) override;
    void plot2Dto1DEffcyProj() override;
    void plot2D_SingleMuonEffcy(const std::string& var) override;
    void plot2Dto1DsingleMuonEffcyProj(const std::pair<const std::string,std::pair<bool,bool>>& kv) override;
};

// ============================================================
// Derived: PbPb
// ============================================================
class TrigEffPlotterPbPb : public TrigEffcyPlotBaseClass, public PbPbBaseClass
{
public:
    TrigEffPlotterPbPb(int runYear,
                       const StrVec& vars1D_in,
                       const std::string ctr_binning_version_user = "default",
                       bool draw2mu4_user = false,
                       bool add_ctr_0_10_user = false,
                       const StrVec& filter_suffixes = {},
                       const XYLog&  logopt_ = {},
                       const StrVec& vars2D_in = {},
                       const StrVec& vars2DProf_in = {},
                       const std::map<std::string,Rect>& legSigned  = {},
                       const std::map<std::string,Rect>& legSignal  = {},
                       const std::map<std::string,Rect>& legOpSig   = {},
                       const std::map<std::string,std::pair<bool,bool>>& var2DProj = {},
                       const std::map<std::string,std::pair<bool,bool>>& var2DSingleMuonEffcyProj = {})
        : TrigEffcyPlotBaseClass(runYear, vars1D_in,
                                 draw2mu4_user, 
                                 filter_suffixes, logopt_,
                                 vars2D_in, vars2DProf_in,
                                 legSigned, legSignal, legOpSig,
                                 var2DProj, var2DSingleMuonEffcyProj),
          add_ctr_0_10(add_ctr_0_10_user)
    {
        ctr_binning_version = ctr_binning_version_user;
    }

private:
    bool add_ctr_0_10{};

    void configureDataFiles(int runYear) override;


    StrVec buildTriggerList() const override
    {
        if (draw2mu4) return {"2mu4", "mu4_mu4noL1"};
        return {"mu4_mu4noL1"};
    }

    // plotting functions implemented in PbPb file
    void initialize() override{
        PbPbBaseClass::InitializePbPb();
        TrigEffcyPlotBaseClass::initialize();
    }
    
    void plot1D(const std::string& var) override;
    void plot2D(const std::string& var, const std::string& fs) override;
    void plotProfile(const std::string& var, const std::string& fs) override;
    void plot2Dto1DEffcyProj() override;
    void plot2D_SingleMuonEffcy(const std::string& var) override;
    void plot2Dto1DsingleMuonEffcyProj(const std::pair<const std::string,std::pair<bool,bool>>& kv) override;
};
