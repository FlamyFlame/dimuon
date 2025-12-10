// BookingHistograms.cpp — centralised histogram factory for trigger‑efficiency study
// ----------------------------------------------------------------------------------
//  * Holds three unordered_maps that own every TH1D, TH2D and the derived
//    “variable vs pT”  2‑D histograms.
//  * Exposes three public booking functions so your analysis class can just call
//       book1D();  book2D();  bookDerived2D();
//    in its constructor.
//
// 2025‑06‑05 – generated on‑demand for ChatGPT user.
// ----------------------------------------------------------------------------------
#include <unordered_map>
#include <string>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "ParamsSet.h"          // provides the global ‘pms’ object & bin counts

class TrigEffBooking
{
public:
    // three public containers (use them directly for Fill/Draw)
    std::unordered_map<std::string, TH1D*> h1d;        // all original  1‑D
    std::unordered_map<std::string, TH2D*> h2d;        // all original  2‑D
    std::unordered_map<std::string, TH2D*> h2d_vsPt;   // derived 2‑D “var vs pT”

    // -------------------------------------------------------- booking entry points
    void book1D();               // replace every  ‘new TH1D’ line
    void book2D();               // replace every  ‘new TH2D’ line (existing)
    void bookDerived2D();        // build   var×pTpair   &  var×pT1st   from 1‑D

private:
    // helpers to reduce boiler‑plate
    TH1D* make1D(const char* name, const char* title,
                 int n, double lo, double hi)
    {
        TH1D* h = new TH1D(name, title, n, lo, hi);
        h1d[name] = h;
        return h;
    }

    TH1D* make1D(const char* name, const char* title,
                 int n, const double* bins)
    {
        TH1D* h = new TH1D(name, title, n, bins);
        h1d[name] = h;
        return h;
    }

    TH2D* make2D(const char* name, const char* xtitle, const char* ytitle,
                 int nX, const double* xbins,
                 int nY, double yLo, double yHi)
    {
        std::string fullTitle = std::string(";") + xtitle + ";" + ytitle;
        TH2D* h = new TH2D(name, fullTitle.c_str(), nX, xbins, nY, yLo, yHi);
        h2d[name] = h;
        return h;
    }

    TH2D* make2D(const char* name, const char* xtitle, const char* ytitle,
                 int nX, const double* xbins,
                 int nY, const double* ybins)
    {
        std::string fullTitle = std::string(";") + xtitle + ";" + ytitle;
        TH2D* h = new TH2D(name, fullTitle.c_str(), nX, xbins, nY, ybins);
        h2d[name] = h;
        return h;
    }

    // -------------------------------------- table of 1‑D variable definitions
    struct Var1D { const char* stem; int n; double lo; double hi; };
    struct Var1DLog { const char* stem; int n; const double* bins; };
};

// ============================================================================
//   1)  Booking of all original 1‑D histograms
// ============================================================================
void TrigEffBooking::book1D()
{
    // variables with uniform bins (Δη, Δφ, ΔR, …) ---------------------------
    const Var1D uni[] = {
        {"Deta",          nDeta_bins_trig_effcy,                -4.8,              4.8},
        {"Deta_zoomin",   nDR_deta_dphi_zoomin_bins_trig_effcy, -0.8,              0.8},
        {"Dphi",          nDphi_bins_trig_effcy,                -pms.PI,           pms.PI},
        {"Dphi_zoomin",   nDR_deta_dphi_zoomin_bins_trig_effcy, -0.8,              0.8},
        {"DR",            nDR_bins_trig_effcy,                  0.0,               pms.deltaR_thrsh[2]},
        {"DR_zoomin",     nDR_deta_dphi_zoomin_bins_trig_effcy, 0.0,               0.8},
        {"minv_zoomin",   nminv_bins_trig_effcy,                0.0,               3.0}
    };

    // variables with custom bin arrays --------------------------------------
    const Var1DLog log40[] = {
        {"pt2nd",       int(pms.pT_bins_40.size()-1), pms.pT_bins_40.data()}
    };
    const Var1DLog log80[] = {
        {"pair_pt_log", int(pms.pT_bins_80.size()-1), pms.pT_bins_80.data()}
    };

    // groups (mu4, mu4_mu4noL1, …) ------------------------------------------
    const char* groups[] = {
        "mu4_w_single_b_sig_sel",
        "mu4_mu4noL1_w_single_b_sig_sel",
        "mu4_mu4noL1_excl_w_single_b_sig_sel",
        "2mu4_w_single_b_sig_sel"
    };

    // book uniform‑bin 1‑D ---------------------------------------------------
    for (auto g : groups)
        for (const auto& v : uni) {
            std::string name = std::string("h_") + v.stem + "_" + g;
            std::string title = ";" + std::string(v.stem).replace(0, 4, "#Delta") + ";"; // crude ylabel
            make1D(name.c_str(), title.c_str(), v.n, v.lo, v.hi);
        }

    // book log‑bin (40‑bin scheme) ------------------------------------------
    for (auto g : groups)
        for (const auto& v : log40) {
            std::string name = std::string("h_") + v.stem + "_" + g;
            make1D(name.c_str(), ";p_{T} [GeV];", v.n, v.bins);
        }

    // book log‑bin (80‑bin scheme) ------------------------------------------
    for (auto g : groups)
        for (const auto& v : log80) {
            std::string name = std::string("h_") + v.stem + "_" + g;
            make1D(name.c_str(), ";p_{T}^{pair} [GeV];", v.n, v.bins);
        }
}

// ============================================================================
//   2)  Booking of *existing* 2‑D histograms (Eta vs pT, Δη vs Δφ, …)
// ============================================================================
void TrigEffBooking::book2D()
{
    // For brevity we only show a couple – extend as needed -------------------
    // Example: pair η vs pair pT  (uniform Y axis)
    make2D("h_pair_eta_vs_pair_pT_mu4_w_single_b_sig_sel",
           "p_{T}^{pair} [GeV]", "#eta^{pair}",
           int(pms.pT_bins_80.size()-1), pms.pT_bins_80.data(),
           neta_bins_trig_effcy, -2.4, 2.4);

    // … repeat for every original TH2D in your list …
}

// ============================================================================
//   3)  Derived 2‑D (variable × pT1st   and   variable × pTpair)
// ============================================================================
void TrigEffBooking::bookDerived2D()
{
    const int n40 = int(pms.pT_bins_40.size()-1);
    const int n80 = int(pms.pT_bins_80.size()-1);
    const double* b40 = pms.pT_bins_40.data();
    const double* b80 = pms.pT_bins_80.data();

    // iterate over every 1‑D we have just booked ----------------------------
    for (const auto& kv : h1d) {
        const std::string& hname = kv.first;      // e.g.  h_DR_mu4_w_single_b_sig_sel
        TH1D*  h1 = kv.second;

        // keep only those matching “…_mu4_w_single_b_sig_sel”
        if (hname.find("_mu4_w_single_b_sig_sel") == std::string::npos) continue;

        // extract the variable stem  DR, Deta, …  (between h_ and _mu4_…)
        const size_t p1 = hname.find('_');            // first underscore (after h)
        const size_t p2 = hname.find("_mu4_");       // start of group suffix
        std::string varStem = hname.substr(p1+1, p2-p1-1);

        // Y‑axis spec comes from the 1‑D histogram itself --------------------
        int nY = h1->GetNbinsX();
        double yLo = h1->GetXaxis()->GetXmin();
        double yHi = h1->GetXaxis()->GetXmax();
        std::string ytitle = h1->GetXaxis()->GetTitle();

        // build & book  variable vs pTpair  ----------------------------------
        std::string namePair = "h_" + varStem + "_vs_pair_pT_mu4_w_single_b_sig_sel";
        TH2D* hPair = new TH2D(namePair.c_str(),
                               (std::string(";p_{T}^{pair} [GeV];") + ytitle).c_str(),
                               n80, b80, nY, yLo, yHi);
        h2d_vsPt[namePair] = hPair;

        // build & book  variable vs pT1st  -----------------------------------
        std::string nameLead = "h_" + varStem + "_vs_pT1st_mu4_w_single_b_sig_sel";
        TH2D* hLead = new TH2D(nameLead.c_str(),
                                (std::string(";p_{T,1st} [GeV];") + ytitle).c_str(),
                                n40, b40, nY, yLo, yHi);
        h2d_vsPt[nameLead] = hLead;
    }
}
