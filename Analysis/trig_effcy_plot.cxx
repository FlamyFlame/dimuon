// ──────────────────────────────────────────────────────────────────────────────
//  TrigEffPlotter (updated August‑2025)
//
//  Major refactor implementing flexible filter framework & three signal drawing
//  modes requested on 2025‑08‑05.
//
//  • Replaces hard‑wired booleans (plot_sepr, plot_good_accept, …) by a vector
//    of filter suffixes (filter_suffix_list).
//  • Adds SignalDrawingMode enum {Signed, Signal, OpAndSignal} and associated
//    png‑suffix map.
//  • Introduces generic maps that describe
//      – how each filter transforms the trigger names
//      – which drawing modes are relevant for the filter
//      – per‑filter variable list exceptions
//  • plot1D() now loops over filter_suffix_list × drawing modes, building the
//    histogram names on the fly via the maps, and automatically produces the
//    requested comparison plots.
//  • The original 2‑D efficiency and profile helpers are kept (renamed
//    plot2D() / plotProfile()) so that the heavy refactor is confined to the
//    1‑D code path.
//  • Misc. cosmetic clean‑ups & better helper utilities.
//
//  NOTE:  *This single header‑implementation file is meant to be compiled
//         once in a ROOT session with:  .L trig_effcy_plot_updated.cxx+
//
//  Author: ChatGPT‑o3 (refactor for Yuhan, 2025‑08‑05)
// ──────────────────────────────────────────────────────────────────────────────
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TGraphAsymmErrors.h>

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <stdexcept>

// ──────────────────────────────────────────────────────────────────────────────
//  Small helpers
// ──────────────────────────────────────────────────────────────────────────────
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

// ──────────────────────────────────────────────────────────────────────────────
//  TrigEffPlotter class (refactored)
// ──────────────────────────────────────────────────────────────────────────────
class TrigEffPlotter
{
public:
    // ---------------- enum & typedefs ---------------------------------------
    enum class SignalDrawingMode { Signed, Signal, OpAndSignal };

    using StrVec  = std::vector<std::string>;
    using XYLog   = std::map<std::string, std::pair<bool, bool>>;
    using ModeMap = std::map<SignalDrawingMode, bool>;
    using TriggerMap = std::map<std::string, std::string>; // trigger → mapped‑name
    using Rect = std::array<double,4>;

    // ---------------- ctor ---------------------------------------------------
    TrigEffPlotter(int runYear,
                   const StrVec& vars1D_in,
                   bool          draw2mu4_user = false,
                   bool          use_sepr_for_op_and_signal_user = false,
                   bool          debug_mode_user = true,
                   const StrVec& filter_suffixes = {},
                   const XYLog&  logopt_ = {},
                   const StrVec& vars2D_in = {},
                   const StrVec& vars2DProf_in = {},
                   const std::map<std::string,Rect>& legSigned  = {},
                   const std::map<std::string,Rect>& legSignal  = {},
                   const std::map<std::string,Rect>& legOpSig   = {},
                   const std::map<std::string,std::pair<bool,bool>>& var2DProj = {})
        : draw2mu4(draw2mu4_user),
          use_sepr_for_op_and_signal (use_sepr_for_op_and_signal_user),
          debug_mode (debug_mode_user),
          fFile(nullptr),
          var1Ds(vars1D_in),
          var2Ds(vars2D_in),
          var2Ds_for_profile(vars2DProf_in),
          logopt(logopt_),
          filter_suffix_list(filter_suffixes)
    {
        legendPosSigned         = legSigned;
        legendPosSignal         = legSignal;
        legendPosOpSignal       = legOpSig;
        var2Ds_for1Deffcy_proj  = var2DProj;

        isRun2pp = (runYear % 2000 == 17);
        switch (runYear % 2000) {
        case 17:
            data_dir    = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_run2/";
            infile_name = data_dir + "histograms_real_pairs_pp_2017_single_mu4.root";
            break;
        case 24:
            data_dir    = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_2024/";
            infile_name = data_dir + "histograms_real_pairs_pp_2024_single_mu4.root";
            break;
        default:
            throw std::runtime_error("runYear must be 17 or 24 for pp data");
        }

        buildFilterLabelMaps();     // fills maps with defaults + exceptions
        buildGlobalMaps();     // fills maps with defaults + exceptions
        consistencyChecks();   // sanity on user‑provided filter list
        initializeFile();
    }

    // ---------------- public driver -----------------------------------------
    void Run()
    {
        if (!fFile || fFile->IsZombie()) return;
        gStyle->SetOptStat(0);

        if (isRun2pp){
            trg_list = {"2mu4"};
        } else{
            if (draw2mu4){
                trg_list = {"2mu4", "mu4_mu4noL1"};
            } else{
                trg_list = {"mu4_mu4noL1"};
            }
        }

        // Which 1‑D variables shall be plotted for this job?
        StrVec vars_to_use = determineVar1DList();

        for (const auto& v : vars_to_use)  plot1D(v);
        // for (const auto& v : var2Ds)       plot2D(v);
        // for (const auto& v : var2Ds_for_profile) plotProfile(v);

        // --- 2-D efficiencies --------------------------------------------------
        {
            StrVec fList = filter_suffix_list.empty() ? StrVec{""} : filter_suffix_list;
            for (const auto& fs : fList) {
                if (!filt_to_draw2D_map.at(fs)) continue;
                for (const auto& v : determineVar2DList(fs))
                    plot2D(v, fs);
            }
        }
        // --- 2-D profiles ------------------------------------------------------
        {
            StrVec fList = filter_suffix_list.empty() ? StrVec{""} : filter_suffix_list;
            for (const auto& fs : fList) {
                if (!filt_to_draw2Dprofile_map.at(fs)) continue;
                for (const auto& v : determineVar2DPList(fs))
                    plotProfile(v, fs);
            }
        }

        // --- project selected 2D histograms into 1D efficiencies ----------------
        plot2Dto1DEffcyProj(); 
    }

private:
    // =========================================================================
    //  Data members
    // =========================================================================
    bool isRun2pp{};
    std::string data_dir;
    std::string infile_name;

    bool draw2mu4{};
    bool use_sepr_for_op_and_signal{};
    bool debug_mode{};

    // trigger list
    std::vector<std::string> trg_list;
    
    // ROOT interface
    TFile* fFile;

    // variable lists
    StrVec var1Ds;
    StrVec var2Ds;
    StrVec var2Ds_for_profile;
    
    // project selected 2D hists to 1D (X/Y) and plot efficiencies like 1D
    std::map<std::string, std::pair<bool,bool>> var2Ds_for1Deffcy_proj; // key: "Y_vs_X" → {projX, projY}

    // user options
    XYLog logopt;
    StrVec filter_suffix_list; // may be empty (→ only default)

    // user-tunable legend rectangles
    std::map<std::string,Rect> legendPosSigned;
    std::map<std::string,Rect> legendPosSignal;
    std::map<std::string,Rect> legendPosOpSignal;   // for OpAndSignal mode

    // --- drawing‑mode infrastructure ----------------------------------------
    std::string op_and_sig_png_suffix = use_sepr_for_op_and_signal? "_op_sepr_and_sig_sel" : "_all_op_and_sig_sel";
    std::map<SignalDrawingMode, std::string> mode_to_png_suffix = {
        {SignalDrawingMode::Signed,        ""},
        {SignalDrawingMode::Signal,        "_w_sig_sel"},
        {SignalDrawingMode::OpAndSignal,   "_all_op_and_sig_sel"}
    };

    // --- central trigger/filter maps ----------------------------------------
    //  key = filter_suffix
    std::map<std::string, TriggerMap>             filt_to_trig_map;           // default mapping
    std::map<std::string, TriggerMap>             filt_to_trig_exception_map; // user provided exceptions
    std::map<std::string, ModeMap>                filt_to_mode_map;           // default all true
    std::map<std::string, ModeMap>                filt_to_mode_exception_map; // user provided exceptions
    std::map<std::string, StrVec>                 filt_to_var1D_exception_map;// per‑filter var list

    // ── 2-D control ──────────────────────────────────────────────────────────────
    std::map<std::string,bool>  filt_to_draw2D_map;          // draw full 2-D eff?
    std::map<std::string,bool>  filt_to_draw2Dprofile_map;   // draw profile?
    std::map<std::string,StrVec> filt_to_var2D_exception_map;
    std::map<std::string,StrVec> filt_to_var2Dprofile_exception_map;

    // --- filter suffix to label map ----------------------------------------
    std::map<std::string, std::string>            filt_suffix_to_label_map;   // map of filter suffix to label (used in legend)

    // some constants ----------------------------------------------------------
    const StrVec kTriggers = {"mu4", "mu4_mu4noL1", "2mu4"};

    // =========================================================================
    //  Initialization helpers
    // =========================================================================
    void initializeFile()
    {
        fFile = TFile::Open(infile_name.c_str(), "READ");
        if (!fFile || fFile->IsZombie())
            std::cerr << "[TrigEffPlotter] ERROR: cannot open " << infile_name << "\n";
    }

    // ------------------------------------------------------------------------
    void buildFilterLabelMaps()
    {
        filt_suffix_to_label_map["_good_accept"] = "good accept";
        filt_suffix_to_label_map["_inv_w_by_single_mu_effcy"] = "1/(single muon effcy) weight";
    }
    
    // ------------------------------------------------------------------------
    void buildGlobalMaps()
    {
        // 1)  Default trigger mapping: trigger → trigger+suffix
        for (const auto& fs : filter_suffix_list) {
            TriggerMap tmap;

            for (const auto& trg : kTriggers)
                tmap[trg] = trg + fs; // simple append rule
            tmap["mu4_mu4noL1_denom"] = "NONE"; // no mu4_mu4noL1-specific denominator by default
            tmap["2mu4_denom"] = "NONE"; // no mu4_mu4noL1-specific denominator by default
            
            filt_to_trig_map[fs] = tmap;
        }
        // The default ("no filter") entry – represented by empty string
        {
            TriggerMap tmap;

            for (const auto& trg : kTriggers) tmap[trg] = trg;
            tmap["mu4_mu4noL1_denom"] = "NONE"; // no mu4_mu4noL1-specific denominator by default
            tmap["2mu4_denom"] = "NONE"; // no mu4_mu4noL1-specific denominator by default

            filt_to_trig_map[""] = tmap; // special key for default
        }

        // 2)  Drawing‑mode map – default: signed+signal true, op_and_signal false
        for (const auto& fs : filter_suffix_list) {
            filt_to_mode_map[fs] = {{SignalDrawingMode::Signed, true},
                                    {SignalDrawingMode::Signal, true},
                                    {SignalDrawingMode::OpAndSignal, false}};
        }
        filt_to_mode_map[""] = {{SignalDrawingMode::Signed, true},
                                 {SignalDrawingMode::Signal, true},
                                 {SignalDrawingMode::OpAndSignal, false}};

        // 3)  Per‑filter *exceptions* ------------------------------------------------
        //     Hard‑code known special behaviours.  Extend here when new filter
        //     types are introduced.

        // ---------------------------- _MB exceptions -----------------------------------------------------
        // mu4 → MB, mu4_mu4noL1 → (skip), 2mu4 → mu4
        filt_to_trig_exception_map["_MB"] = {
            {"mu4", "MB"},
            {"mu4_mu4noL1", "NONE"},
            {"2mu4", "mu4"}
        };

        // 4) 2-D drawing flags – default OFF for every user filter, ON for default “”
        filt_to_draw2D_map[""]         = true;
        filt_to_draw2Dprofile_map[""]  = true;
        for (const auto& fs : filter_suffix_list) {
            filt_to_draw2D_map        [fs] = false;
            filt_to_draw2Dprofile_map [fs] = false;
        }

        // ---------------------------- _sepr exceptions -----------------------------------------------------
        // Only signed‑mode relevant for _sepr
        filt_to_mode_exception_map["_sepr"] = {
            {SignalDrawingMode::Signal, false},
        };

        filt_to_draw2D_map["_sepr"]              = true;
        filt_to_var2D_exception_map["_sepr"]     = {"pt2nd_vs_q_eta_2nd", "pair_eta_vs_pair_pT", "Deta_Dphi", "eta1_eta2", "eta_avg_Deta", "eta_avg_Dphi", "minv_pair_pt_log"};

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

        // ---------------------------- good/bad accept special per‑variable list ----------------------
        // Only signed‑mode relevant for _good_accept
        filt_to_mode_exception_map["_good_accept"] = {
            {SignalDrawingMode::Signal, false},
        };

        // Only signed‑mode relevant for _bad_accept
        filt_to_mode_exception_map["_bad_accept"] = {
            {SignalDrawingMode::Signal, false},
        };

        // ---------------------------------------------------------------------
        //  Merge exceptions into defaults
        // ---------------------------------------------------------------------
        for (const auto& [fs, exc] : filt_to_trig_exception_map) {
            auto& base = filt_to_trig_map[fs];
            for (const auto& [trg, newname] : exc) base[trg] = newname;
        }
        for (const auto& [fs, exc] : filt_to_mode_exception_map) {
            auto& base = filt_to_mode_map[fs];
            for (const auto& [md, val] : exc) base[md] = val;
        }

        // Also make sure the default "" entry exists in every map (already ok)
    }

    // ------------------------------------------------------------------------
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
            for (auto mode : {SignalDrawingMode::Signed, SignalDrawingMode::Signal}) {
                bool val0 = filt_to_mode_map.at(fs0).at(mode);
                for (size_t i = 1; i < filter_suffix_list.size(); ++i) {
                    bool vali = filt_to_mode_map.at(filter_suffix_list[i]).at(mode);
                    if (vali != val0)
                        throw std::runtime_error("Filters disagree on required drawing modes");
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    StrVec determineVar1DList() const
    {
        if (filter_suffix_list.empty()) return var1Ds;
        const auto& fs0 = filter_suffix_list.front();
        auto it = filt_to_var1D_exception_map.find(fs0);
        if (it == filt_to_var1D_exception_map.end()) return var1Ds;
        return it->second;
    }

    StrVec determineVar2DList (const std::string& fs) const
    {
        auto it = filt_to_var2D_exception_map.find(fs);
        return (it==filt_to_var2D_exception_map.end()) ? var2Ds : it->second;
    }

    StrVec determineVar2DPList (const std::string& fs) const
    {
        auto it = filt_to_var2Dprofile_exception_map.find(fs);
        return (it==filt_to_var2Dprofile_exception_map.end()) ? var2Ds_for_profile
                                                              : it->second;
    }

    // =========================================================================
    //  Histogram access helper
    // =========================================================================
    template<typename T> T* getHist(const std::string& name) const
    {
        T* h = dynamic_cast<T*>(fFile->Get(name.c_str()));
        if (!h)
            std::cerr << "[TrigEffPlotter] WARNING: missing hist " << name << "\n";
        return h;
    }

    // =========================================================================
    //  Name helper (after filter mapping)
    // =========================================================================
    std::string mappedTrigger(const std::string& filter_suffix, const std::string& trigger) const
    {
        const auto& m = filt_to_trig_map.at(filter_suffix);
        auto it = m.find(trigger);
        if (it == m.end()) return ""; // should not happen
        return it->second;
    }

    // =========================================================================
    // Clip numerator to denominator for the special single-mu inverse-weight filter
    // Also enforce the weighted-events consistency: err_num^2 ≤ err_den^2
    // We touch underflow/overflow as TEfficiency checks them too
    // =========================================================================
    std::vector<std::pair<int,double>>
    clipNumeratorIfInvW(TH1* hNum, TH1* hDen, const std::string& fs) const
    {
        if (fs != "_inv_w_by_single_mu_effcy" || !hNum || !hDen) return {};

        // Ensure Sumw2 exists so SetBinError controls ∑w²
        if (!hNum->GetSumw2N()) hNum->Sumw2();
        if (!hDen->GetSumw2N()) hDen->Sumw2();

        if (hNum->GetNbinsX() != hDen->GetNbinsX()) return {};

        std::vector<std::pair<int,double>> bins_above1;

        const int nb = hNum->GetNbinsX();
        auto clamp_one = [&](int b) {
            double n  = hNum->GetBinContent(b);
            double d  = hDen->GetBinContent(b);
            double ne = hNum->GetBinError(b);
            double de = hDen->GetBinError(b);

            if (n > d) {
                if (d > 0) bins_above1.emplace_back(b, n / d);
                hNum->SetBinContent(b, d);
            }
            // Enforce ∑w² consistency as well.
            if (ne > de) hNum->SetBinError(b, de);

            // If denominator is zero, zero numerator too to keep TEfficiency happy.
            if (d <= 0.0) {
                hNum->SetBinContent(b, 0.0);
                hNum->SetBinError(b,   0.0);
            }
        };

        // Underflow, in-range, overflow
        clamp_one(0);
        for (int b = 1; b <= nb; ++b) clamp_one(b);
        clamp_one(nb+1);

        return bins_above1;
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
    std::pair<std::string,std::string> parseYvsX(const std::string& var2d) const {
        auto pos = var2d.find("_vs_");
        if (pos == std::string::npos) return {var2d, ""};
        return { var2d.substr(0, pos), var2d.substr(pos+4) };
    }

    // =========================================================================
    // re-use legend placement maps with a 1D variable key (X or Y)
    // =========================================================================
    TrigEffPlotter::Rect legendBoxFor1DVar(const std::string& v1d,
                                           SignalDrawingMode mode,
                                           const StrVec& filters) const
    {
        Rect def = (mode==SignalDrawingMode::Signed)
                   ? Rect{0.43,0.12,0.80,0.35}
                   : Rect{0.50,0.14,0.77,0.30};

        if      (mode==SignalDrawingMode::Signed  && legendPosSigned.count(v1d))  return legendPosSigned.at(v1d);
        else if (mode==SignalDrawingMode::Signal  && legendPosSignal.count(v1d))  return legendPosSignal.at(v1d);
        else if (mode==SignalDrawingMode::OpAndSignal && legendPosOpSignal.count(v1d)) return legendPosOpSignal.at(v1d);

        // fallback: auto-resize by filter label length (same as in plot1D)
        auto autoBox = [&](Rect r, const StrVec& fsList)->Rect {
            size_t L=0;
            for (auto& fs : fsList) if (!fs.empty()) {
                std::string lab = (filt_suffix_to_label_map.count(fs) ? filt_suffix_to_label_map.at(fs)
                                                                      : (fs[0]=='_'? fs.substr(1):fs));
                L = std::max(L, lab.size());
            }
            if (L>4) { double extra = std::min(0.01*(L-4), 0.20); r[0]=std::max(0.05,r[0]-extra); }
            return r;
        };
        return autoBox(def, filters);
    }

    // =========================================================================
    //  plot1D – main refactored 1‑D efficiency routine
    // =========================================================================
    void plot1D(const std::string& var)
    {

        // ── enlarge default legend box if filter-labels are long ─────────────
        auto autoBox = [&](Rect r, const StrVec& filters) -> Rect
        {
            // longest suffix length (ignore the empty string for the default curve)
            size_t L = 0;
            for (auto& fs : filters){
                if (!fs.empty()){
                    std::string filter_label = fs;
                    if (filt_suffix_to_label_map.find(fs) != filt_suffix_to_label_map.end()){ // filter label found
                        filter_label = filt_suffix_to_label_map.at(fs);
                    }
                    L = std::max(L, filter_label.size());
                }
            }

            if (L > 4) {                                   // threshold
                double extra = std::min(0.01*(L - 4), 0.20);   // 0.01 per char, cap 0.10
                r[0] = std::max(0.05, r[0] - extra);           // move left
                // r[3] = std::min(0.95, r[3] + extra*0.1);       // comment out: do not grow talled
            }
            return r;
        };

        auto zeroUFof = [](TH1* h){
            if (!h) return;
            const int nb = h->GetNbinsX();
            h->SetBinContent(0, 0.0);   h->SetBinError(0, 0.0);
            h->SetBinContent(nb+1, 0.0); h->SetBinError(nb+1, 0.0);
        };

        // ===============================================================
        // Decide which filters will be looped over (default + comparisons)
        // ===============================================================
        StrVec filters = filter_suffix_list;
        bool plot_weighted = find(filters.begin(), filters.end(), "_inv_w_by_single_mu_effcy") != filters.end();
        if (plot_weighted){ // inverted by weight --> ploted weighted first to set the y-axis range
            filters.push_back("");
        }else{
            filters.insert(filters.begin(), ""); // no weighting, ensure default first
        }

        // For each drawing mode iterate ---------------------------------------------------
        for (auto mode : {SignalDrawingMode::Signed, SignalDrawingMode::Signal}) {
            // Is this mode requested by *any* of the (comparison) filters?
            bool modeNeeded = filter_suffix_list.empty() ? true : std::any_of(filter_suffix_list.begin(), filter_suffix_list.end(), [&](const std::string& fs){
                return filt_to_mode_map.at(fs).at(mode);
            });
            if (!modeNeeded) continue;

            // Canvas bookkeeping -------------------------------------------------------
            bool isSigned = (mode == SignalDrawingMode::Signed);
            std::string canvasName = Form("c_%s_%s", var.c_str(), (isSigned?"signed":"signal"));
            int   canvW  = isSigned ? 1200 : 700;
            int   canvH  = 500;
            TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), canvW, canvH);
            if (isSigned) c->Divide(2,1);

            // Prepare legend & colour styles ----------------------------------------
            Color_t colBase[3] = {kRed+1, kAzure+2, kGreen+2}; // 2mu4, mu4_mu4noL1, (spare)
            Color_t colFilt[3] = {kOrange+7, kAzure+1, kGreen+1}; // one-shade-lighter palette for filtered / comparison curves
            Color_t colFilt2[3] = {kMagenta, kAzure, kGreen+3};
            Style_t mkrBase[3] = {24, 20, 22};

            // Iterate over pads (sign1 / sign2 or single)
            int nPads = isSigned?2:1;
            std::string signSuffixes[2] = {"_sign1", "_sign2"};

            for (int pad = 0; pad < nPads; ++pad) {
                TPad* padPtr = isSigned? (TPad*)c->cd(pad+1) : (TPad*)c->cd();
                auto xy = xyFor(var);
                padPtr->SetLogx(xy.first); padPtr->SetLogy(xy.second);

                Rect def = {0.43,0.12,0.8,0.35};
                Rect def_excl_inv_w_by_single_mu_effcy = {0.45,0.45,0.83,0.63}; // use for "_excl_inv_w_by_single_mu_effcy"
                Rect def_default_only = {0.5,0.14,0.77,0.3};
                if (filter_suffix_list.empty()) def = def_default_only;

                if (find(filters.begin(), filters.end(), "_excl_inv_w_by_single_mu_effcy") != filters.end()) def = def_excl_inv_w_by_single_mu_effcy;

                Rect box = def;

                if      (mode==SignalDrawingMode::Signed  && legendPosSigned.count(var))
                    box = legendPosSigned.at(var);
                else if (mode==SignalDrawingMode::Signal  && legendPosSignal.count(var))
                    box = legendPosSignal.at(var);
                else if (!filter_suffix_list.empty()){ // adjust legend box size to the label length of the (first) filter
                    box = autoBox(box, filter_suffix_list);
                }

                TLegend* leg = new TLegend(box[0], box[1], box[2], box[3]);
                leg->SetBorderSize(0);

                bool firstDrawn = false;

                // loop over filters → lines (default first, then comparisons)
                int filterIdx = 0;
                for (const auto& fs : filters) {
                    if (!filt_to_mode_map.at(fs).at(mode)) continue; // skip filter that should not appear in this mode

                    // For each (valid) trigger produce TGraphAsymmErrors
                    std::vector<std::pair<std::string,int>> trigOrder; // order & colour idx

                    if (isRun2pp) trigOrder = { {"2mu4", 0} };
                    else{
                        if (draw2mu4){
                            trigOrder = { {"mu4_mu4noL1", 0}, {"2mu4", 1} }; 
                        }else{
                            trigOrder = { {"mu4_mu4noL1", 0} };
                        }
                    }

                    for (auto [trg, colourIdx] : trigOrder) {

                        // denominator hist (mu4 / trigger-specific denominator)
                        std::string denomName = (mappedTrigger(fs, trg + "_denom") != "NONE")? mappedTrigger(fs, trg + "_denom") : mappedTrigger(fs, "mu4");
                        
                        if (denomName.empty() || denomName == "NONE") {
                            throw std::runtime_error("Base trigger mu4 missing for filter " + fs);
                        }

                        std::string denomHist = Form("h_%s_%s%s", var.c_str(), denomName.c_str(),
                                                      (mode==SignalDrawingMode::Signal?"_w_single_b_sig_sel":signSuffixes[pad].c_str()));
                        
                        if (debug_mode) std::cout << "Denominator name: " << denomHist << std::endl;
                        TH1* h_denom = getHist<TH1>(denomHist);
                        if (!h_denom) continue; // can't proceed – skip filter


                        std::string mapped = mappedTrigger(fs, trg);
                        if (mapped.empty() || mapped == "NONE") continue; // skip incompatible trigger
                        std::string numHist = Form("h_%s_%s%s", var.c_str(), mapped.c_str(),
                                                   (mode==SignalDrawingMode::Signal?"_w_single_b_sig_sel":signSuffixes[pad].c_str()));
                        TH1* h_num = getHist<TH1>(numHist);
                        if (!h_num) continue;

                        auto bins_above1 = clipNumeratorIfInvW(h_num, h_denom, fs);

                        double ymax = 1.;

                        zeroUFof(h_num);
                        zeroUFof(h_denom);

                        // Build efficiency graph
                        auto g = new TGraphAsymmErrors();
                        g->BayesDivide(h_num, h_denom);

                        // if (bins_above1.size() > 0){
                        //     for (auto bpair : bins_above1){
                        //         g->SetPointY(bpair.first, bpair.second);
                        //         ymax = std::max(ymax, bpair.second);
                        //     }
                        // }

                        Color_t useCol;
                        if (!plot_weighted){
                            useCol = (filterIdx == 0) ? colBase[colourIdx]   // default
                                                    : colFilt[colourIdx];  // filtered
                            if (filterIdx == 2) useCol = colFilt2[colourIdx];
                        } else{ // plot weighted: default appears at last
                            useCol = (filterIdx == filters.size() - 1) ? colBase[colourIdx]   // filtered
                                                                       : colFilt[colourIdx];  // base
                            if (filters.size() == 3 && filterIdx == 1) useCol = colFilt2[colourIdx];
                        }

                        SetStyle(g, useCol, mkrBase[colourIdx]);
                        bool isFiltered = plot_weighted? (filterIdx < filters.size() - 1) : (filterIdx > 0);
                        if (isFiltered) {
                            // comparison filter → dashed
                            g->SetLineStyle(2);
                            g->SetMarkerStyle(mkrBase[colourIdx]+4);
                        }

                        if (!firstDrawn) {
                            g->Draw("APL");
                            // g->GetYaxis()->SetRangeUser(0.,ymax);
                            g->GetYaxis()->SetRangeUser(0.,1.);
                            g->GetXaxis()->SetTitle(h_denom->GetXaxis()->GetTitle());
                            g->GetYaxis()->SetTitle("#epsilon");
                            
                            if (xy.first) adjustLogXRange(g, h_denom);

                            firstDrawn = true;
                        } else {
                            g->Draw("PL SAME");
                        }

                        // Legend entry text
                        std::string label = trg;
                        if (!fs.empty()) {
                            std::string fs_label = (fs[0] == '_')? fs.substr(1) : fs; // default filter label: filter suffix but trimming off "_"
                            if (filt_suffix_to_label_map.find(fs) != filt_suffix_to_label_map.end()){ // filter label found
                                fs_label = filt_suffix_to_label_map.at(fs);
                            }
                            label += ", " + fs_label;
                        }
                        leg->AddEntry(g, label.c_str(), "lep");
                    } // end loop over triggers
                    ++filterIdx;
                } // end loop over filters
                leg->Draw();
            } // end loop over pads (subplots)

            // Save png -----------------------------------------------------------------
            std::string pngDir = data_dir + "trig_effcy_plots";
            for (const auto& fs : filter_suffix_list) {
                pngDir += fs;
            }
            makeDirIfNeeded(pngDir);

            std::string fn = pngDir + "/" + var + "_trig_effcy";
            for (const auto& fs : filter_suffix_list) {
                fn += fs;
            }
            fn += mode_to_png_suffix.at(mode) + ".png";
            c->SaveAs(fn.c_str());
        }

        // ================================================================
        //  Special OpAndSignal – one file per *comparison* filter
        // ================================================================
        // ----------------------------------------------------------------
        //  Op-and-signal mode:
        //    • if filter_suffix_list is empty  --> draw DEFAULT
        //    • otherwise                      --> draw one file per filter (skip default; check if OpAndSignal mode turned on)
        // ----------------------------------------------------------------
        StrVec opFilters = filter_suffix_list.empty() ? StrVec{""} : filter_suffix_list;
        for (const auto& fs : opFilters) {
            // skip comparison filters that don’t request this mode
            if (!fs.empty() &&
                !filt_to_mode_map.at(fs).at(SignalDrawingMode::OpAndSignal)) continue;


            TCanvas* c = new TCanvas(Form("c_opand_%s_%s", var.c_str(), fs.c_str()), "", 700, 500);
            auto xy = xyFor(var);
            c->SetLogx(xy.first); c->SetLogy(xy.second);

            int colourIdx = (!draw2mu4)? 1 : 0; // use Azure for 2mu4; if mu4_mu4noL1 present will be red

            // Legend
            Rect def = {0.47,0.12,0.87,0.36};
            if (legendPosOpSignal.count(var))
                def = legendPosOpSignal.at(var);
            else if (!filter_suffix_list.empty()){ // adjust legend box size to the label length of the (first) filter
                def = autoBox(def, filter_suffix_list);
            }

            TLegend* L = new TLegend(def[0],def[1],def[2],def[3]);
            L->SetBorderSize(0);

            bool firstCurve = true;

            for (auto trg : trg_list) {

                std::string denomName = (mappedTrigger(fs, trg + "_denom") != "NONE")? mappedTrigger(fs, trg + "_denom") : mappedTrigger(fs, "mu4");
                if (denomName.empty() || denomName == "NONE") continue;

                std::string denomSign2 = use_sepr_for_op_and_signal? Form("h_%s_%s_sepr_sign2", var.c_str(), denomName.c_str()) : Form("h_%s_%s_sign2", var.c_str(), denomName.c_str());
                std::string denomSel   = Form("h_%s_%s_w_single_b_sig_sel", var.c_str(), denomName.c_str());

                if (debug_mode){
                    std::cout << "Denominator name all opposite sign: " << denomSign2 << std::endl;
                    std::cout << "Denominator name signal: " << denomSel << std::endl;
                }

                TH1* h_den_sign2 = getHist<TH1>(denomSign2);
                TH1* h_den_sel   = getHist<TH1>(denomSel);
                if (!h_den_sign2 || !h_den_sel) continue;

                std::string mapped = mappedTrigger(fs, trg);
                if (mapped.empty() || mapped == "NONE") continue;

                // _sign2
                std::string numSign2 = use_sepr_for_op_and_signal? Form("h_%s_%s_sepr_sign2", var.c_str(), mapped.c_str()) : Form("h_%s_%s_sign2", var.c_str(), mapped.c_str());
                // _w_single_b_sig_sel
                std::string numSel   = Form("h_%s_%s_w_single_b_sig_sel", var.c_str(), mapped.c_str());

                TH1* h_num_s2 = getHist<TH1>(numSign2);
                TH1* h_num_sel = getHist<TH1>(numSel);
                if (!h_num_s2 || !h_num_sel) continue;

                auto bins_above1_s2 = clipNumeratorIfInvW(h_num_s2, h_den_sign2, fs);
                auto bins_above1_sel = clipNumeratorIfInvW(h_num_sel, h_den_sel, fs);

                double ymax_opandsig = 1.;

                // sign2 efficiency
                auto g_s2 = new TGraphAsymmErrors();
                g_s2->BayesDivide(h_num_s2, h_den_sign2);
                SetStyle(g_s2, colourIdx==0? kAzure+2 : kRed+1, colourIdx==0?24:20);

                // if (bins_above1_s2.size() > 0){
                //     for (auto bpair : bins_above1_s2){
                //         g_s2->SetPointY(bpair.first, bpair.second);
                //         ymax_opandsig = std::max(ymax_opandsig, bpair.second);
                //     }
                // }

                // sig‑sel efficiency
                auto g_sel = new TGraphAsymmErrors();
                g_sel->BayesDivide(h_num_sel, h_den_sel);
                SetStyle(g_sel, colourIdx==0? kAzure+1 : kRed+2, colourIdx==0?25:24, 2, 2);

                // if (bins_above1_sel.size() > 0){
                //     for (auto bpair : bins_above1_sel){
                //         g_sel->SetPointY(bpair.first, bpair.second);
                //         ymax_opandsig = std::max(ymax_opandsig, bpair.second);
                //     }
                // }

                if (firstCurve) {
                    g_s2->Draw("APL");                       // sets axes with first curve
                    // g_s2->GetYaxis()->SetRangeUser(0., ymax_opandsig);
                    g_s2->GetYaxis()->SetRangeUser(0., 1.);
                    g_s2->GetXaxis()->SetTitle(h_den_sel->GetXaxis()->GetTitle());
                    g_s2->GetYaxis()->SetTitle("#epsilon");

                    if (xy.first) adjustLogXRange(g_s2, h_den_sign2);

                    firstCurve = false;                      // ← axis is now set
                } else {
                    g_s2->Draw("PL SAME");                   // any later curve
                }
                g_sel->Draw("PL SAME");                      // dashed sig-sel curve

                colourIdx++;

                std::string op_legend_fs = trg;
                std::string sig_legend_fs = trg;

                if (!fs.empty()){
                    std::string fs_label = (fs[0] == '_')? fs.substr(1) : fs; // default filter label: filter suffix but trimming off "_"
                    if (filt_suffix_to_label_map.find(fs) != filt_suffix_to_label_map.end()){ // filter label found
                        fs_label = filt_suffix_to_label_map.at(fs);
                    }                    
                    op_legend_fs += " (all op sign, " + fs_label + ")";
                    sig_legend_fs += " (single b signal, " + fs_label + ")";
                }else{
                    op_legend_fs += " (all op sign)";
                    sig_legend_fs += " (single b signal)";
                }

                L->AddEntry(g_s2, op_legend_fs.c_str(), "lep");
                L->AddEntry(g_sel, sig_legend_fs.c_str(), "lep");

            } // finish loop over muon-pair triggers

            L->Draw();
            
            std::string outdir = data_dir + "trig_effcy_plots" + (fs.empty() ? "" : fs);
            makeDirIfNeeded(outdir + "/op_and_sig");
            std::string fn = Form("%s/op_and_sig/%s_trig_effcy%s%s.png", outdir.c_str(), var.c_str(), fs.c_str(), mode_to_png_suffix.at(SignalDrawingMode::OpAndSignal).c_str());
            c->SaveAs(fn.c_str());
        }
    }

    // =========================================================================
    //  xyFor helper (unchanged)
    // =========================================================================
    std::pair<bool,bool> xyFor(const std::string& v) const
    {
        auto it = logopt.find(v);
        return (it==logopt.end()) ? std::make_pair(false,false) : it->second;
    }

    // ======================================================================
    // 2-D efficiency canvases  (Signed + Signal modes) – one PNG per filter
    // ======================================================================
    void plot2D(const std::string& var, const std::string& fs)
    {
        const auto xy = xyFor(var);

        // ---------- build histogram name helpers --------------------------
        auto hName = [&](const std::string& trg,
                         const std::string& suffix)->std::string
        {
            std::string mapped = mappedTrigger(fs, trg);
            return (mapped.empty() || mapped=="NONE") ? "" :
                   Form("h_%s_%s%s", var.c_str(), mapped.c_str(), suffix.c_str());
        };

        // ---------- first canvas : sign-1 / sign-2 ------------------------
        std::vector<std::pair<std::string,std::string>> configs;

        if (isRun2pp) {                       // Run-2 pp → only 2mu4 / mu4
            configs = {{"2mu4","mu4"}};
        } else {                              // heavy-ion / Run-3
            if (draw2mu4){
                configs = {{"mu4_mu4noL1","mu4"},{"2mu4","mu4"}};
            } else{
                configs = {{"mu4_mu4noL1","mu4"}};
            }
        }

        int nPads = static_cast<int>(configs.size())*2;         // each config → 2 signs
        int nCols = 2;
        int nRows = (nPads/2);

        TCanvas* c = new TCanvas(Form("c2d_%s_%s",var.c_str(),fs.c_str()),"",1100,450*nRows);
        c->Divide(nCols,nRows);

        int pad=1;
        for (auto& cfg : configs)
        for (auto sign : {std::string("_sign1"),std::string("_sign2")})
        {
            auto numH = getHist<TH2>(hName(cfg.first , sign));
            auto denH = getHist<TH2>(hName(cfg.second, sign));
            if (!numH||!denH) {++pad; continue;}
            numH = (TH2*)numH->Clone();           // keep originals
            numH->Divide(denH);
            c->cd(pad++);
            gPad->SetRightMargin(0.15);
            gPad->SetLogx(xy.first); gPad->SetLogy(xy.second);
            numH->SetTitle(
                Form("%s / %s  (%s)",
                     cfg.first.c_str(),cfg.second.c_str(),
                     sign=="_sign1"?"same":"opposite"));
            numH->Draw("COLZ");
        }

        std::string outdir = data_dir + "trig_effcy_plots" + (fs.empty()?"":fs);
        makeDirIfNeeded(outdir);
        c->SaveAs(Form("%s/%s_2D_trig_effcy%s.png",
                       outdir.c_str(),var.c_str(),fs.c_str()));

        // ---------- second canvas : signal-selection ----------------------
        TCanvas* cB = new TCanvas(Form("c2ds_%s_%s",var.c_str(),fs.c_str()),"",700,450);
        cB->cd(); gPad->SetRightMargin(0.15);
        gPad->SetLogx(xy.first); gPad->SetLogy(xy.second);

        // denominator always mu4, sig-selection
        auto hDenSel = getHist<TH2>(hName("mu4","_w_single_b_sig_sel"));
        if (!hDenSel) return;

        for (auto trg : trg_list)
        {
            auto hNumSel = getHist<TH2>(hName(trg,"_w_single_b_sig_sel"));
            if (!hNumSel) continue;
            hNumSel = (TH2*)hNumSel->Clone();
            hNumSel->Divide(hDenSel);
            hNumSel->SetTitle((trg+", signal selection").c_str());
            hNumSel->Draw("COLZ SAME");           // pile in same pad
        }
        cB->SaveAs(Form("%s/%s_2D_trig_effcy%s_w_sig_sel.png",
                        outdir.c_str(),var.c_str(),fs.c_str()));
    }

    // ======================================================================
    // 2-D profile canvases  (no division) – one PNG per filter
    // ======================================================================
    void plotProfile(const std::string& var, const std::string& fs)
    {
        const auto xy = xyFor(var);

        auto hName = [&](const std::string& trg,
                         const std::string& suffix)->std::string
        {
            std::string mapped = mappedTrigger(fs, trg);
            return (mapped.empty() || mapped=="NONE") ? "" :
                   Form("h_%s_%s%s", var.c_str(), mapped.c_str(), suffix.c_str());
        };

        // -------- signed view ---------------------------------------------------
        TCanvas* c = new TCanvas(Form("cProf_%s_%s",var.c_str(),fs.c_str()),"",1100,500);
        c->Divide(2,1);

        for (int i=0;i<2;++i) {
            std::string sign = (i==0?"_sign1":"_sign2");
            auto h = getHist<TH2>(hName("mu4",sign));
            if (!h) continue;
            c->cd(i+1);
            gPad->SetRightMargin(0.15);
            gPad->SetLogx(xy.first); gPad->SetLogy(xy.second); gPad->SetLogz();
            h->SetTitle(Form("mu4 %s",sign=="_sign1"?"same":"opposite"));
            h->Draw("COLZ");
            auto pfx = h->ProfileX(Form("pfx%d",i),1,-1,"s");
            pfx->SetLineColor(kRed); pfx->SetLineWidth(2);
            pfx->Draw("E SAME");
        }

        std::string outdir = data_dir + "trig_effcy_plots" + (fs.empty()?"":fs);
        makeDirIfNeeded(outdir);
        c->SaveAs(Form("%s/%s_profile%s.png",outdir.c_str(),var.c_str(),fs.c_str()));

        // -------- signal-selection profile --------------------------------------
        auto hSel = getHist<TH2>(hName("mu4","_w_single_b_sig_sel"));
        if (!hSel) return;

        TCanvas* cS = new TCanvas(Form("cProfS_%s_%s",var.c_str(),fs.c_str()),"",550,450);
        cS->cd(); gPad->SetRightMargin(0.15);
        gPad->SetLogx(xy.first); gPad->SetLogy(xy.second); gPad->SetLogz();
        hSel->SetTitle("mu4, signal selection");
        hSel->Draw("COLZ");
        auto pS = hSel->ProfileX("pS",1,-1,"s");
        pS->SetLineColor(kRed); pS->SetLineWidth(2);
        pS->Draw("E SAME");

        cS->SaveAs(Form("%s/%s_profile%s_w_sig_sel.png",
                        outdir.c_str(),var.c_str(),fs.c_str()));
    }

    // ======================================================================
    // Project selected 2D histograms to 1D (X/Y) and plot efficiencies like plot1D()
    // ======================================================================
    void plot2Dto1DEffcyProj()
    {
        auto zeroUFof = [](TH1* h){
            if (!h) return;
            const int nb = h->GetNbinsX();
            h->SetBinContent(0, 0.0);   h->SetBinError(0, 0.0);
            h->SetBinContent(nb+1, 0.0); h->SetBinError(nb+1, 0.0);
        };

        if (var2Ds_for1Deffcy_proj.empty()) return;

        // default + comparisons (keep weighted filter behaviour consistent with plot1D)
        StrVec filters = filter_suffix_list;
        const bool hasWeighted = (std::find(filters.begin(), filters.end(),
                                            std::string("_inv_w_by_single_mu_effcy")) != filters.end());
        if (hasWeighted) filters.push_back(""); else filters.insert(filters.begin(), "");

        for (const auto& kv : var2Ds_for1Deffcy_proj)
        {
            const std::string& var2d = kv.first;
            const bool doPX = kv.second.first;   // ProjectionX → onto X (right side of "..._vs_")
            const bool doPY = kv.second.second;  // ProjectionY → onto Y (left side of "..._vs_")

            const auto parsed = parseYvsX(var2d);
            const std::string& yVar = parsed.first;
            const std::string& xVar = parsed.second;

            for (int proj = 0; proj < 2; ++proj)
            {
                const bool isProjX = (proj == 0);
                if ((isProjX && !doPX) || (!isProjX && !doPY)) continue;

                const std::string projVar = isProjX ? (xVar.empty() ? var2d : xVar)
                                                    : (yVar.empty() ? var2d : yVar);
                if (projVar.empty()) continue; // malformed name

                // ---------------- Signed / Signal (shared canvas structure) ----------------
                for (auto mode : {SignalDrawingMode::Signed, SignalDrawingMode::Signal})
                {
                    // mode gate (same rule as plot1D)
                    bool modeNeeded = filter_suffix_list.empty();
                    if (!modeNeeded) {
                        for (const auto& fs : filter_suffix_list) {
                            if (filt_to_mode_map.at(fs).at(mode)) { modeNeeded = true; break; }
                        }
                    }
                    if (!modeNeeded) continue;

                    const bool isSigned = (mode == SignalDrawingMode::Signed);
                    const int  nPads    = isSigned ? 2 : 1;
                    const char* signSuffixes[2] = {"_sign1", "_sign2"};

                    // canvas
                    std::string cname = Form("c_proj_%s_%s_%s",
                                             var2d.c_str(),
                                             isProjX ? "PX" : "PY",
                                             isSigned ? "signed" : "signal");
                    TCanvas* c = new TCanvas(cname.c_str(), "", isSigned ? 1200 : 700, 500);
                    if (isSigned) c->Divide(2,1);

                    // legend box (reuse 1D var positioning)
                    Rect box = legendBoxFor1DVar(projVar, mode, filter_suffix_list);

                    // styles: [0] → mu4_mu4noL1 (red), [1] → 2mu4 (azure)
                    Color_t colBase[2] = {kRed+1,   kAzure+2};
                    Color_t colFilt[2] = {kRed+2,   kAzure+1};
                    Style_t mkrBase[2] = {20, 24};

                    for (int pad = 0; pad < nPads; ++pad)
                    {
                        TPad* P = isSigned ? (TPad*)c->cd(pad+1) : (TPad*)c->cd();
                        (void)P;
                        auto xy = xyFor(projVar);
                        gPad->SetLogx(xy.first);
                        gPad->SetLogy(xy.second);

                        TLegend* leg = new TLegend(box[0], box[1], box[2], box[3]);
                        leg->SetBorderSize(0);

                        bool firstDrawn = false;
                        int  filterIdx  = 0;

                        for (const auto& fs : filters)
                        {
                            if (!filt_to_mode_map.at(fs).at(mode)) { ++filterIdx; continue; }

                            // triggers to overlay (skip mu4; it’s denom)
                            std::vector<std::pair<std::string,int>> trigOrder =
                                isRun2pp ? std::vector<std::pair<std::string,int>>{{"2mu4", 1}}
                                         : (draw2mu4? std::vector<std::pair<std::string,int>>{{"mu4_mu4noL1", 0}, {"2mu4", 1}}
                                                    : std::vector<std::pair<std::string,int>>{{"mu4_mu4noL1", 0}});

                            bool firstTrgInstance = true; // boolean to indicate first trigger instance: if continue without plotting, only increment filterIdx if first trigger instance

                            for (const auto& tp : trigOrder)
                            {
                                const std::string& trg = tp.first;
                                const int idx = tp.second;
                             
                                // denominator 2D → project
                                const std::string denomName = (mappedTrigger(fs, trg + "_denom") != "NONE")? mappedTrigger(fs, trg + "_denom") : mappedTrigger(fs, "mu4");
                                
                                if (denomName.empty() || denomName == "NONE") {
                                    if (firstTrgInstance){
                                        ++filterIdx;
                                        firstTrgInstance = false;
                                    }
                                    continue;
                                }

                                std::string denom2D = Form("h_%s_%s%s",
                                    var2d.c_str(), denomName.c_str(),
                                    (mode == SignalDrawingMode::Signal ? "_w_single_b_sig_sel"
                                                                       : signSuffixes[pad]));
                                if (debug_mode) std::cout << "Denominator name: " << denom2D << std::endl;
        
                                TH2* hDen2D = getHist<TH2>(denom2D);
                                if (!hDen2D || hDen2D->GetDimension() != 2) {
                                    if (firstTrgInstance){
                                        ++filterIdx;
                                        firstTrgInstance = false;
                                    }
                                    continue;
                                }

                                std::unique_ptr<TH1> hDen1D(
                                    isProjX ? hDen2D->ProjectionX(Form("px_%s_%d_den", var2d.c_str(), pad), 1, -1, "e")
                                            : hDen2D->ProjectionY(Form("py_%s_%d_den", var2d.c_str(), pad), 1, -1, "e")
                                );
                                if (!hDen1D) {
                                    if (firstTrgInstance){
                                        ++filterIdx;
                                        firstTrgInstance = false;
                                    }
                                    continue;
                                }

                                const std::string mapped = mappedTrigger(fs, trg);
                                if (mapped.empty() || mapped == "NONE") continue;

                                std::string num2D = Form("h_%s_%s%s",
                                    var2d.c_str(), mapped.c_str(),
                                    (mode == SignalDrawingMode::Signal ? "_w_single_b_sig_sel"
                                                                       : signSuffixes[pad]));
                                TH2* hNum2D = getHist<TH2>(num2D);
                                if (!hNum2D || hNum2D->GetDimension() != 2) continue;

                                std::unique_ptr<TH1> hNum1D(
                                    isProjX ? hNum2D->ProjectionX(Form("px_%s_%d_%s", var2d.c_str(), pad, trg.c_str()), 1, -1, "e")
                                            : hNum2D->ProjectionY(Form("py_%s_%d_%s", var2d.c_str(), pad, trg.c_str()), 1, -1, "e")
                                );
                                if (!hNum1D) continue;

                                // special weighting guard
                                clipNumeratorIfInvW(hNum1D.get(), hDen1D.get(), fs);

                                zeroUFof(hNum1D.get());
                                zeroUFof(hDen1D.get());

                                // graph
                                auto g = new TGraphAsymmErrors();
                                g->BayesDivide(hNum1D.get(), hDen1D.get());

                                // choose colours: default vs filtered
                                Color_t useCol = (!hasWeighted)
                                                 ? ((filterIdx == 0) ? colBase[idx] : colFilt[idx])
                                                 : ((filterIdx == (int)filters.size()-1) ? colBase[idx] : colFilt[idx]);
                                SetStyle(g, useCol, mkrBase[idx]);

                                const bool isFiltered = (!hasWeighted) ? (filterIdx > 0)
                                                                       : (filterIdx < (int)filters.size()-1);
                                if (isFiltered) {
                                    g->SetLineStyle(2);
                                    g->SetMarkerStyle(mkrBase[idx] + 4);
                                }

                                if (!firstDrawn) {
                                    g->Draw("APL");
                                    g->GetYaxis()->SetRangeUser(0., 1.);
                                    g->GetXaxis()->SetTitle(hDen1D->GetXaxis()->GetTitle());
                                    g->GetYaxis()->SetTitle("#epsilon");
                                    if (xy.first) adjustLogXRange(g, hDen1D.get());
                                    firstDrawn = true;
                                } else {
                                    g->Draw("PL SAME");
                                }

                                // legend text
                                std::string fs_label;
                                if (!fs.empty()) {
                                    fs_label = (filt_suffix_to_label_map.count(fs)
                                               ? filt_suffix_to_label_map.at(fs)
                                               : (fs[0] == '_' ? fs.substr(1) : fs));
                                }
                                std::string label = trg;
                                if (!fs_label.empty()) label += ", " + fs_label;
                                leg->AddEntry(g, label.c_str(), "lep");
                            } // loop over triggers
                            ++filterIdx;
                        } // loop over filters
                        leg->Draw();
                    } // loop over pads

                    // output (use projected 1D var name)
                    std::string pngDir = data_dir + "trig_effcy_plots";
                    for (const auto& fs : filter_suffix_list) {
                        pngDir += fs;
                    }
                    makeDirIfNeeded(pngDir);
                    std::string fn = pngDir + "/";
                    fn += projVar + "_trig_effcy";
                    for (const auto& fsu : filter_suffix_list) fn += fsu; // keep your naming convention
                    fn += mode_to_png_suffix.at(mode) + ".png";
                    c->SaveAs(fn.c_str());
                } // loop over signed & signal modes

                // ---------------- OpAndSignal (one file per filter) ----------------
                {
                    StrVec opFilters = filter_suffix_list.empty() ? StrVec{""} : filter_suffix_list;
                    for (const auto& fs : opFilters)
                    {
                        if (!fs.empty() && !filt_to_mode_map.at(fs).at(SignalDrawingMode::OpAndSignal)) continue;

                        auto xy = xyFor(projVar);
                        TCanvas* c = new TCanvas(
                            Form("c_opandsig_proj_%s_%s_%s", var2d.c_str(), (isProjX ? "PX" : "PY"), fs.c_str()),
                            "", 700, 500
                        );
                        c->SetLogx(xy.first);
                        c->SetLogy(xy.second);

                        Rect Lbox = legendBoxFor1DVar(projVar, SignalDrawingMode::OpAndSignal, filter_suffix_list);
                        TLegend* L = new TLegend(Lbox[0], Lbox[1], Lbox[2], Lbox[3]);
                        L->SetBorderSize(0);

                        bool firstCurve = true;
                        int  colourIdx  = 0;  // 0 → Azure, 1 → Red

                        for (const std::string& trg : trg_list)
                        {
                            const std::string denomName = (mappedTrigger(fs, trg + "_denom") != "NONE")? mappedTrigger(fs, trg + "_denom") : mappedTrigger(fs, "mu4");

                            if (denomName.empty() || denomName == "NONE") continue;

                            // denominators (2D → project)
                            std::string den2D_s2  = use_sepr_for_op_and_signal? Form("h_%s_%s_sepr_sign2", var2d.c_str(), denomName.c_str()) : Form("h_%s_%s_sign2", var2d.c_str(), denomName.c_str());
                            std::string den2D_sel = Form("h_%s_%s_w_single_b_sig_sel",    var2d.c_str(), denomName.c_str());
                            
                            if (debug_mode){
                                std::cout << "Denominator name all opposite sign: " << den2D_s2 << std::endl;
                                std::cout << "Denominator name signal pairs: " << den2D_sel << std::endl;
                            }

                            TH2* hDen2D_s2  = getHist<TH2>(den2D_s2);
                            TH2* hDen2D_sel = getHist<TH2>(den2D_sel);
                            if (!hDen2D_s2 || !hDen2D_sel) continue;

                            std::unique_ptr<TH1> hDen1D_s2(
                                isProjX ? hDen2D_s2->ProjectionX(Form("px_%s_den_s2",  var2d.c_str()), 1, -1, "e")
                                        : hDen2D_s2->ProjectionY(Form("py_%s_den_s2",  var2d.c_str()), 1, -1, "e")
                            );
                            std::unique_ptr<TH1> hDen1D_sel(
                                isProjX ? hDen2D_sel->ProjectionX(Form("px_%s_den_sel", var2d.c_str()), 1, -1, "e")
                                        : hDen2D_sel->ProjectionY(Form("py_%s_den_sel", var2d.c_str()), 1, -1, "e")
                            );
                            if (!hDen1D_s2 || !hDen1D_sel) continue;

                            const std::string mapped = mappedTrigger(fs, trg);
                            if (mapped.empty() || mapped == "NONE") continue;

                            std::string num2D_s2  = use_sepr_for_op_and_signal? Form("h_%s_%s_sepr_sign2", var2d.c_str(), mapped.c_str()) : Form("h_%s_%s_sign2", var2d.c_str(), mapped.c_str());
                            std::string num2D_sel = Form("h_%s_%s_w_single_b_sig_sel", var2d.c_str(), mapped.c_str());
                            TH2* hNum2D_s2  = getHist<TH2>(num2D_s2);
                            TH2* hNum2D_sel = getHist<TH2>(num2D_sel);
                            if (!hNum2D_s2 || !hNum2D_sel) continue;

                            std::unique_ptr<TH1> hNum1D_s2(
                                isProjX ? hNum2D_s2->ProjectionX(Form("px_%s_%s_s2",  var2d.c_str(), trg.c_str()), 1, -1, "e")
                                        : hNum2D_s2->ProjectionY(Form("py_%s_%s_s2",  var2d.c_str(), trg.c_str()), 1, -1, "e")
                            );
                            std::unique_ptr<TH1> hNum1D_sel(
                                isProjX ? hNum2D_sel->ProjectionX(Form("px_%s_%s_sel", var2d.c_str(), trg.c_str()), 1, -1, "e")
                                        : hNum2D_sel->ProjectionY(Form("py_%s_%s_sel", var2d.c_str(), trg.c_str()), 1, -1, "e")
                            );
                            if (!hNum1D_s2 || !hNum1D_sel) continue;

                            clipNumeratorIfInvW(hNum1D_s2.get(),  hDen1D_s2.get(),  fs);
                            clipNumeratorIfInvW(hNum1D_sel.get(), hDen1D_sel.get(), fs);

                            auto g_s2  = new TGraphAsymmErrors();
                            g_s2->BayesDivide(hNum1D_s2.get(),  hDen1D_s2.get());
                            SetStyle(g_s2,  colourIdx==0 ? kAzure+2 : kRed+1, colourIdx==0 ? 24 : 20);

                            auto g_sel = new TGraphAsymmErrors();
                            g_sel->BayesDivide(hNum1D_sel.get(), hDen1D_sel.get());
                            SetStyle(g_sel, colourIdx==0 ? kAzure+1 : kRed+2, colourIdx==0 ? 25 : 24, 2, 2);

                            if (firstCurve) {
                                g_sel->Draw("APL");  // draw dashed first so solid stays visible
                                g_sel->GetYaxis()->SetRangeUser(0., 1.);
                                g_sel->GetXaxis()->SetTitle(hDen1D_sel->GetXaxis()->GetTitle());
                                g_sel->GetYaxis()->SetTitle("#epsilon");
                                if (xy.first) adjustLogXRange(g_sel, hDen1D_sel.get());
                                firstCurve = false;
                            } else {
                                g_sel->Draw("PL SAME");
                            }
                            g_s2->Draw("PL SAME");

                            // legend text
                            std::string fs_label;
                            if (!fs.empty()) {
                                fs_label = (filt_suffix_to_label_map.count(fs)
                                           ? filt_suffix_to_label_map.at(fs)
                                           : (fs[0] == '_' ? fs.substr(1) : fs));
                            }
                            std::string labS2  = trg + " (all op sign";
                            std::string labSel = trg + " (single b signal";
                            if (!fs_label.empty()) {
                                labS2  += ", " + fs_label;
                                labSel += ", " + fs_label;
                            }
                            labS2  += ")";
                            labSel += ")";
                            L->AddEntry(g_s2,  labS2.c_str(),  "lep");
                            L->AddEntry(g_sel, labSel.c_str(), "lep");

                            ++colourIdx;
                        }
                        L->Draw();

                        std::string outdir = data_dir + "trig_effcy_plots" + (fs.empty() ? "" : fs);
                        makeDirIfNeeded(outdir + "/op_and_sig");
                        std::string fn = Form("%s/op_and_sig/%s_trig_effcy%s%s.png",
                                              outdir.c_str(), projVar.c_str(), fs.c_str(),
                                              mode_to_png_suffix.at(SignalDrawingMode::OpAndSignal).c_str());
                        c->SaveAs(fn.c_str());
                    }
                }
            }
        }
    }
};

// ──────────────────────────────────────────────────────────────────────────────
//  Steering macro (example)
// ──────────────────────────────────────────────────────────────────────────────
void trig_effcy_plot(){
    std::vector<std::string> var1Ds = {"Deta", "Deta_zoomin", "Dphi", "Dphi_zoomin", "DR", "DR_zoomin", "DR_0_2", "pt2nd", "minv_zoomin", "pair_pt_log"};
    std::vector<std::string> var2Ds = {"pt2nd_vs_q_eta_2nd", "pt2nd_vs_phi2nd", "phi2nd_vs_q_eta_2nd", "DR_zoomin_vs_pt2nd", "DR_0_2_vs_pt2nd", "pair_eta_vs_pair_pT", "Deta_Dphi", "eta1_eta2", "eta_avg_Deta", "eta_avg_Dphi", "minv_pair_pt_log",
                                       // "Deta_vs_pT_1st", "Deta_zoomin_vs_pT_1st", "Dphi_vs_pT_1st", "Dphi_zoomin_vs_pT_1st", "DR_vs_pT_1st", "DR_zoomin_vs_pT_1st", "minv_zoomin_vs_pT_1st", "pair_pt_log_vs_pT_1st", "pt2nd_vs_pT_1st", 
                                       // "Deta_vs_pair_pT", "Deta_zoomin_vs_pair_pT", "Dphi_vs_pair_pT", "Dphi_zoomin_vs_pair_pT", "DR_vs_pair_pT", "DR_zoomin_vs_pair_pT", "minv_zoomin_vs_pair_pT", "pair_pt_log_vs_pair_pT", "pt2nd_vs_pair_pT",
                                      };
    std::vector<std::string> var2DsProf = {
       "DR_zoomin_vs_pt2nd",
       // "Deta_vs_pT_1st", "Deta_zoomin_vs_pT_1st", "Dphi_vs_pT_1st", "Dphi_zoomin_vs_pT_1st", "DR_vs_pT_1st", "DR_zoomin_vs_pT_1st", "minv_zoomin_vs_pT_1st", "pair_pt_log_vs_pT_1st", "pt2nd_vs_pT_1st", 
       // "Deta_vs_pair_pT", "Deta_zoomin_vs_pair_pT", "Dphi_vs_pair_pT", "Dphi_zoomin_vs_pair_pT", "DR_vs_pair_pT", "DR_zoomin_vs_pair_pT", "minv_zoomin_vs_pair_pT", "pair_pt_log_vs_pair_pT", "pt2nd_vs_pair_pT",
    };

    std::map<std::string,std::pair<bool,bool>> logaxes = {
       {"pt2nd",{true,false}},
       {"pair_pt_log",{true,false}},
       {"pt2nd_vs_q_eta_2nd",{false,true}},
       {"pair_eta_vs_pair_pT",{true,false}},
       {"minv_pair_pt_log",{true,true}},
       {"Deta_vs_pT_1st",{true,false}},
       {"Deta_zoomin_vs_pT_1st",{true,false}},
       {"Dphi_vs_pT_1st",{true,false}},
       {"Dphi_zoomin_vs_pT_1st",{true,false}},
       {"DR_vs_pT_1st",{true,false}},
       {"DR_zoomin_vs_pT_1st",{true,false}},
       {"minv_zoomin_vs_pT_1st",{true,false}},
       {"pair_pt_log_vs_pT_1st",{true,true}},
       {"pt2nd_vs_pT_1st",{true,true}},
       {"Deta_vs_pair_pT",{true,false}},
       {"Deta_zoomin_vs_pair_pT",{true,false}},
       {"Dphi_vs_pair_pT",{true,false}},
       {"Dphi_zoomin_vs_pair_pT",{true,false}},
       {"DR_vs_pair_pT",{true,false}},
       {"DR_zoomin_vs_pair_pT",{true,false}},
       {"DR_zoomin_vs_pt2nd",{true,false}},
       {"minv_zoomin_vs_pair_pT",{true,false}},
       {"pair_pt_log_vs_pair_pT",{true,true}},
       {"pt2nd_vs_pair_p",{true,true}},
    };

    std::map<std::string,std::pair<bool,bool>> var2DProj = {
        {"pt2nd_vs_q_eta_2nd", {true,  true}},   // PX (→ q_eta_2nd) and PY (→ pt2nd)
        {"pt2nd_vs_phi2nd",    {true,  false}}  // PX only (→ phi2nd)
    };

    // TrigEffPlotter::Rect rSigned  = {0.60,0.15,0.88,0.35};
    // TrigEffPlotter::Rect rSignal  = {0.20,0.70,0.50,0.88};
    // TrigEffPlotter::Rect rOpSig   = {0.60,0.15,0.88,0.31};

    // std::map<std::string,TrigEffPlotter::Rect> legSigned  = { {"Dphi", rSigned } };
    // std::map<std::string,TrigEffPlotter::Rect> legSignal  = { {"Dphi", rSignal } };

    std::map<std::string,TrigEffPlotter::Rect> legSigned  = {};
    std::map<std::string,TrigEffPlotter::Rect> legSignal  = {};
    std::map<std::string,TrigEffPlotter::Rect> legOpSig   = {};

    std::vector<std::string> filters = {};
    // filters = {"_sepr"};
    // filters = {"_excl"};
    // filters = {"_excl", "_excl_inv_w_by_single_mu_effcy"};
    // std::map<std::string,TrigEffPlotter::Rect> legOpSig   = {{"DR_zoomin", {0.5,0.15,0.87,0.33}  }};
    // filters = {"_good_accept"};
    filters = {"_inv_w_by_single_mu_effcy"};

    TrigEffPlotter plotter(24, var1Ds, false, false, true,
                           filters, logaxes,
                           var2Ds, var2DsProf,
                           legSigned, legSignal, legOpSig,
                           var2DProj);

    plotter.Run();
}
