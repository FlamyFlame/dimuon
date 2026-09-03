// =============================================================================
// FillMCTrigEffHists.cxx
//
// MC-based trigger-efficiency histogram filling (docs/tracking/mc_trigger_efficiency.md
// Physics Procedure §3.1–§3.3). Standalone ACLiC macro over the ntuple-processing
// `_mc_trig` outputs (NEVER raw NTUPs — provenance rule).
//
//   Step 1 (§3.1): single-muon mu4 efficiency inputs from the single-muon tree
//                  (per-charge pt/eta/phi 1D + q·η×pT 2D, denom = offline reco
//                  muon w/ NO trigger requirement, num = muon's own mu4 match).
//   Step 2 (§3.2): ΔR-binned singles-efficiency inputs from the pair trees
//                  (per leg; ΔR bins 0–0.2 / 0.2–1.0 / 1.0+).
//   Step 3 (§3.3): inverse-weighted ε_ΔR inputs (do_step3=true; requires the MC
//                  turn-on fits from FitMCSinglesEffcy.cxx).
//   Step 4 (§3.4): inverse-weighted SINGLE-LEG ε_single(ΔR) inputs (do_step4=true; also
//                  requires the MC turn-on fits). Leg-level analog of Step 3 -> the ΔR
//                  correction on the PbPb union LINEAR terms. NOT needed for pp (2mu4
//                  product); pp runs only to validate the machinery on the full sample.
//
// Selection mirrors the data-side muon definition: Tight WP (analysis nominal),
// pT > 4 GeV, |η| < 2.4. Overlay restricted to 0–5% centrality (doc D2).
//
// TRUTH FIDUCIAL (round 7, user 2026-08-03 -- DEFAULT for every sample and every step):
// on top of the data-like RECO cuts, each MC muon must also satisfy
//     truth pT > 4 GeV  &&  |truth η| < 2.4
// (§3.0(c'), kTruthFiducial* below). Because the sample is truth-SEEDED (round-5 #1), every
// selected muon has a truth partner, so this is a well-defined cut on the muon itself. It
// removes muons that enter the reco fiducial region ONLY through mismeasurement (truth below
// threshold / outside acceptance) -- exactly the "bad muon" population the round-7 sanity
// check targets. It has NO data analogue: the data denominator cannot be truth-gated, so the
// Step-1 MC/data comparison acquires one more MC-only selection (documented asymmetry, §3.0).
//
// FORWARD LOW-pT VETO (round 7, user; §3.5 decision rule -- STEPS 2/3/4 ONLY): a muon is
// rejected iff it is SIMULTANEOUSLY pT < 7 GeV AND q·η < -2. The §3.5 sanity check ruled out
// "bad muons" as the cause of the saturated forward-negative MC turn-on, so those muons are
// removed from the ΔR-correlation measurement. Step 1 and the sanity check keep the full
// acceptance -- the anomaly must stay visible there. See kVetoFwdLowPt below.
// Binnings reuse the DATA conventions:
//   pt   : "pT_bins_single_muon" = pT_bins_8 + pT_bins_60 (RDFBasedHistFillingData.cxx:286-290)
//   q·η  : "eta_bins_trig_effcy" = ParamsSet::makeEtaTrigEffcyBinning(1) (ibid:294)
//   phi  : 128 uniform bins in [-pi, pi] (ibid:300-307)
//   eta  : 48 uniform bins in [-2.4, 2.4]
//   pair pT : "pair_pt_log" = pT_bins_120, 15 log bins 8–120 (var1D_pp.json:102-106)
// MC is ALWAYS weighted: ev_weight (singles) / weight (pairs).
//
// CORRECTED-MC STUDY (round-7 contract item 5, `corrected_mc = true` -- the LAST argument):
// every MC muon that FIRES the trigger carries an extra per-muon weight
//     SF(pT, q·η) = ε_data(pT, q·η) / ε_MC(pT, q·η)
// evaluated CONTINUOUSLY from the two sets of fitted turn-ons (never resample-to-nearest); the
// DENOMINATOR is untouched, so num/denom = ε_MC·<SF> ≈ ε_data. It answers: (1) does correcting
// the MC to the data single-muon efficiency give the same single-muon efficiency as using the
// data-derived ε directly (Step 1)? (2) does it change the ΔR CORRECTIONS (Steps 3/4)? In
// corrected mode Steps 3/4 divide by ε_corr from `single_mu_effcy_pT_fit_mc_corrected*.root`
// (the fit to the CORRECTED turn-ons -- §3.3/§3.4 self-consistency), so that if ε_corr were
// exactly ε_data the numerator weight SF/ε_corr would be the ORIGINAL 1/ε_MC and the ΔR
// corrections would be identical bin-by-bin. Corrected mode skips Step 2 and the L1/HLT
// numerators (no data L1-only reference exists to correct to). Nothing nominal is overwritten:
// every corrected artefact carries "_corrected" in its file name.
//
// Usage (from Analysis/RDFBasedHistFilling/):
//   root -b -l -q 'FillMCTrigEffHists.cxx+("pp")'
//   root -b -l -q 'FillMCTrigEffHists.cxx+("overlay", true)'              // Step 3
//   root -b -l -q 'FillMCTrigEffHists.cxx+("overlay", false, true, true)' // Step 4
//   root -b -l -q 'FillMCTrigEffHists.cxx+("pp_full", false, true, false, false, true)' // corrected Step 1
//
// Output: <sample dir>/mc_trig_eff_hists_<pp24|hijing_overlay_pbpb23>.root
//         (do_step3=true writes a SEPARATE ..._step3.root; do_step4=true a SEPARATE
//          ..._step4.root; neither touches the Step-1/2 file; corrected_mc=true inserts
//          "_corrected" before that suffix)
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TKey.h>
#include <TString.h>
#include <ROOT/RDataFrame.hxx>

using namespace std;

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
// Shared with plot_mc_trig_eff.cxx, which PRINTS this threshold on the sanity-check canvases.
#include "../Utilities/MCTrigEffSanityCfg.h"
#include "../Utilities/MCTrigEffPlateauWindow.h"
#include "../Utilities/MCTrigEffPairPtBinning.h"
#include "../Utilities/proj_range_to_suffix.cxx"
#include "CommonEffcyConfig.h"

namespace MCTrigEff {

// ---------- sample configuration ----------
struct SampleConfig {
    std::string dir;
    std::string label;        // pp24 | hijing_overlay_pbpb23
    std::string pair_file;
    std::string singles_file;
    bool is_overlay = false;  // overlay => 0-5% centrality restriction (doc D2) + mu4-cross-term numerator

    // ---- DATA tag-and-probe reference, for the CORRECTED-MC study only (corrected_mc=true) ----
    // Paths mirror the analysis's own data-side loader EXACTLY:
    //   pp   : RDFBasedHistFillingPP.cxx:306/323      (erf_plus_log fit dir)
    //   PbPb : RDFBasedHistFillingPbPb.cxx:678/694    (fermi_plus_log fit dir)
    // {WP} is replaced by "" (Tight, nominal) or "_medium_wp" at run time -- the data reference
    // MUST be at the SAME working point as the MC (§3.0(d)).
    std::string data_fit_file_tmpl;   // TF1s  f_pt2nd_vs_q_eta2nd<ctr>_<sign>_2mu4_sepr_py_<qeta>_divided
    std::string data_hist_file_tmpl;  // TH2Ds h_pt2nd_vs_q_eta2nd<ctr>_<sign>_2mu4_sepr_divided (gap fallback)
    std::string data_ctr;             // "" (pp) or "_ctr0_5" (overlay, doc D2)
};

// Data tag-and-probe reference paths (corrected-MC study). pp24 data is the reference for every
// pp-collision sample (pp, pp_full, noovl); PbPb23 data 0-5% for the HIJING overlay (doc D2).
namespace {
const std::string kDataPPFitTmpl =
    "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/trg_effcy_pT_fitting_to_erf_plus_log/"
    "single_mu_effcy_pT_fit{WP}.root";
const std::string kDataPPHistTmpl =
    "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
    "histograms_real_pairs_pp_2024_single_mu4_coarse_q_eta_bin_qeta_fid{WP}.root";
const std::string kDataPbPbFitTmpl =
    "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/trg_effcy_pT_fitting_to_fermi_plus_log/"
    "single_mu_effcy_pT_fit{WP}.root";
const std::string kDataPbPbHistTmpl =
    "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/"
    "histograms_real_pairs_pbpb_2023_single_mu4_coarse_q_eta_bin_qeta_fid{WP}.root";

std::string SubstWP(std::string s, const std::string& wp_suf) {
    const std::string tok = "{WP}";
    const size_t p = s.find(tok);
    if (p == std::string::npos)
        throw std::runtime_error("SubstWP: no {WP} placeholder in " + s);
    return s.replace(p, tok.size(), wp_suf);
}
} // anonymous namespace

SampleConfig GetSampleConfig(const std::string& sample) {
    SampleConfig cfg;
    if (sample == "pp") {
        cfg.dir   = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
        cfg.label = "pp24";
        cfg.pair_file    = cfg.dir + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig.root";
        cfg.singles_file = cfg.dir + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_single_muon.root";
        cfg.is_overlay = false;
    } else if (sample == "pp_full") {
        // pp24 FULL sample (the "_pdf" production, 803015-803020): identical physics to "pp"
        // (pp collisions, 2mu4 pair decision, no centrality) -- ONLY the input files differ.
        // Files end in "_full"; the "_full" label keeps its outputs (hists + plots) in their own
        // names/dirs, so they never clobber the TEST-sample trigger-efficiency results.
        cfg.dir   = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/";
        cfg.label = "pp24_full";
        cfg.pair_file    = cfg.dir + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root";
        cfg.singles_file = cfg.dir + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_single_muon_full.root";
        cfg.is_overlay = false;
    } else if (sample == "overlay") {
        cfg.dir   = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/";
        cfg.label = "hijing_overlay_pbpb23";
        cfg.pair_file    = cfg.dir + "muon_pairs_pythia_fullsim_hijing_overlay_pbpb23_no_data_resonance_cuts_mc_trig.root";
        cfg.singles_file = cfg.dir + "muon_pairs_pythia_fullsim_hijing_overlay_pbpb23_no_data_resonance_cuts_mc_trig_single_muon.root";
        cfg.is_overlay = true;
    } else if (sample == "noovl") {
        // r17663 NO-OVERLAY diagnostic (mc_trigger_efficiency.md R8, round 4): pp collisions
        // reconstructed with the overlay's PbPb23-conditions pass, but with NO overlaid event.
        // is_overlay=false: there is no centrality (no HIJING) -> no 0-5% restriction, and the
        // trigger condition is the pp-style 2mu4 pair decision, as for any quiet pp event.
        cfg.dir   = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_no_overlay_test_sample/";
        cfg.label = "r17663_no_overlay";
        cfg.pair_file    = cfg.dir + "muon_pairs_pythia_fullsim_r17663_no_overlay_no_data_resonance_cuts_mc_trig.root";
        cfg.singles_file = cfg.dir + "muon_pairs_pythia_fullsim_r17663_no_overlay_no_data_resonance_cuts_mc_trig_single_muon.root";
        cfg.is_overlay = false;
    } else {
        throw std::invalid_argument("FillMCTrigEffHists: sample must be \"pp\", \"pp_full\", "
                                    "\"overlay\" or \"noovl\", got " + sample);
    }

    // DATA tag-and-probe reference (corrected-MC study only). The overlay is compared to PbPb23
    // data restricted to 0-5% centrality (doc D2 -- the same restriction its own MC selection
    // carries); every pp-collision sample uses the pp24 data reference (identical choice to
    // plot_mc_trig_eff.cxx MakeCfg, incl. "noovl", which simulates pp collisions).
    if (cfg.is_overlay) {
        cfg.data_fit_file_tmpl  = kDataPbPbFitTmpl;
        cfg.data_hist_file_tmpl = kDataPbPbHistTmpl;
        cfg.data_ctr            = "_ctr0_5";
    } else {
        cfg.data_fit_file_tmpl  = kDataPPFitTmpl;
        cfg.data_hist_file_tmpl = kDataPPHistTmpl;
        cfg.data_ctr            = "";
    }
    return cfg;
}

// ---------- binnings (data conventions; see header comment for provenance) ----------
struct Binnings {
    std::vector<double> pt;        // pT_bins_single_muon (incl. the data's duplicated 8.0 edge)
    std::vector<double> q_eta;     // eta_bins_trig_effcy
    std::vector<double> phi;       // 128 uniform [-pi, pi]
    std::vector<double> eta;       // 48 uniform [-2.4, 2.4]
    std::vector<double> dr_zoom;   // 20 uniform [0, 1]
    std::vector<double> dr_full;   // 23 uniform [0, 5.75]
    // Step-3 pair-eta dependence (round-5 #4): coarse pair-pT (ParamsSet::pair_pt_coarse_bins,
    // the crossx binning) x coarse pair-eta (CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap).
    std::vector<double> pair_pt_coarse;
    std::vector<double> pair_eta_coarse;
};

// Contiguous (lo,hi) ranges -> bin edges. Throws if the ranges are not contiguous.
inline std::vector<double> RangesToEdges(const QEtaBinning& ranges) {
    std::vector<double> e;
    for (size_t i = 0; i < ranges.size(); ++i) {
        if (i == 0) e.push_back(ranges[i].first);
        else if (std::fabs(ranges[i].first - ranges[i - 1].second) > 1e-6)
            throw std::runtime_error("RangesToEdges: non-contiguous ranges");
        e.push_back(ranges[i].second);
    }
    return e;
}

Binnings MakeBinnings() {
    ParamsSet pms;
    Binnings b;

    // exact data pt2nd construction (RDFBasedHistFillingData.cxx:286-290): copy pT_bins_8,
    // append pT_bins_60 (this duplicates the 8.0 edge -> one zero-width bin, as in data)
    b.pt = pms.pT_bins_8;
    b.pt.insert(b.pt.end(), pms.pT_bins_60.begin(), pms.pT_bins_60.end());

    b.q_eta = ParamsSet::makeEtaTrigEffcyBinning(1);

    const int nphi = 128;
    b.phi.resize(nphi + 1);
    for (int i = 0; i <= nphi; ++i)
        b.phi[i] = -pms.PI + (2.0 * pms.PI) * (static_cast<double>(i) / nphi);

    const int neta = 48;
    b.eta.resize(neta + 1);
    for (int i = 0; i <= neta; ++i)
        b.eta[i] = -2.4 + 4.8 * (static_cast<double>(i) / neta);


    b.dr_zoom.resize(21);
    for (int i = 0; i <= 20; ++i) b.dr_zoom[i] = i * (1.0 / 20);
    b.dr_full.resize(24);
    for (int i = 0; i <= 23; ++i) b.dr_full[i] = i * (5.75 / 23);

    // round-5 #4: coarse pair-pT (crossx) x coarse pair-eta (crossx pair_eta bins)
    // Nominal 8 log bins, or the opt-in 4-bin comparison variant -- resolved in ONE place
    // (MCTrigEffPairPtBinning.h), which also supplies the matching output token, so the
    // fill and the plot stage cannot disagree about which binning is in use.
    b.pair_pt_coarse = MCTrigEffPairPt::Edges(pms);
    static const CommonEffcyConfig cfg{};
    b.pair_eta_coarse = RangesToEdges(cfg.pair_eta_proj_ranges_coarse_incl_gap); // 9 bins over [-2.4,2.4]

    return b;
}

// ---------- accumulate RDF results under target names, merge at the end ----------
template <class TH>
struct HistAccumulator {
    std::map<std::string, std::vector<ROOT::RDF::RResultPtr<TH>>> acc;
    int booking_counter = 0;

    void add(const std::string& target_name, ROOT::RDF::RResultPtr<TH> r) {
        acc[target_name].push_back(r);
    }

    // triggers the event loops; returns merged hists (caller owns; detached from any directory)
    std::map<std::string, TH*> merge() {
        std::map<std::string, TH*> out;
        for (auto& kv : acc) {
            TH* h_sum = nullptr;
            for (auto& r : kv.second) {
                if (!h_sum) {
                    h_sum = static_cast<TH*>(r->Clone(kv.first.c_str()));
                    h_sum->SetDirectory(nullptr);
                } else {
                    h_sum->Add(r.GetPtr());
                }
            }
            out[kv.first] = h_sum;
        }
        return out;
    }
};

// ---------- MC single-muon efficiency evaluation (Step 3 weights) ----------
// Mirrors the data-side EvaluateSingleMuonEffcyPtFitted (RDFBasedHistFillingData.cxx:625-655):
// TF1 per (charge, fine q·η bin), 2D-ratio fallback for gap q·η. Differences mandated by the
// task spec (pp_trig_eff_highpt_jump lesson): clamp the evaluation pt into the fit's
// [xmin, xmax], floor eps at 0.02 (counted + printed), cap at 1.0.
struct MCEffEvaluator {
    std::map<std::string, TF1*> tf1_map;   // f_mc_pt_vs_q_eta_<muplus|muminus>_<lo>_TO_<hi>
    std::map<std::string, TH2D*> ratio_map; // h_mc_pt_vs_q_eta_ratio_<muplus|muminus>
    CommonEffcyConfig cfg{};
    long long n_floor = 0;      // floor (0.02) firings
    long long n_fallback = 0;   // retired gap fallback; kept as a tripwire, must stay 0
    long long n_eval = 0;

    void LoadFits(const std::string& fit_file) {
        TFile* f = TFile::Open(fit_file.c_str(), "READ");
        if (!f || f->IsZombie())
            throw std::runtime_error("MCEffEvaluator: cannot open MC fit file " + fit_file +
                                     " (run FitMCSinglesEffcy first)");
        TIter next(f->GetListOfKeys());
        TKey* key;
        while ((key = static_cast<TKey*>(next()))) {
            const std::string cls = key->GetClassName();
            const std::string name = key->GetName();
            if (cls == "TF1" && name.rfind("f_mc_pt_vs_q_eta_", 0) == 0) {
                TF1* fn = static_cast<TF1*>(key->ReadObj());
                tf1_map[name] = fn;
            } else if (cls == "TH2D" && name.rfind("h_mc_pt_vs_q_eta_ratio_", 0) == 0) {
                TH2D* h = static_cast<TH2D*>(key->ReadObj());
                h->SetDirectory(nullptr);
                ratio_map[name] = h;
            }
        }
        // keep file open: TF1s (TFormula-based) are independent of it after ReadObj, but
        // do not delete f -> avoids any owned-object lifetime subtleties in a macro.
        std::cout << "MCEffEvaluator: loaded " << tf1_map.size() << " TF1s + "
                  << ratio_map.size() << " fallback ratio TH2Ds from " << fit_file << std::endl;
        if (tf1_map.empty())
            throw std::runtime_error("MCEffEvaluator: no f_mc_pt_vs_q_eta_* TF1s in " + fit_file);
    }

    // ROUND 8: the COARSE q·η binning is contiguous over [-2.4, 2.2) and the fiducial gap cut
    // removes q·η > 2.2, so every surviving muon lands in a fitted bin. There are no holes left
    // and therefore NO 2D fallback: an empty suffix is now a configuration error, not a gap.
    std::string FindQEtaSuffix(float q_eta) const {
        for (const auto& range : cfg.q_eta_proj_ranges_coarse_incl_gap)
            if (q_eta >= range.first && q_eta < range.second) return pairToSuffix(range);
        return "";
    }

    double Eval(float pt, float eta, int charge) {
        ++n_eval;
        const std::string chg = (charge > 0) ? "muplus" : "muminus";
        const float q_eta = charge * eta;
        const std::string q_eta_suffix = FindQEtaSuffix(q_eta);

        double val = -1.0;
        if (!q_eta_suffix.empty()) {
            auto it = tf1_map.find("f_mc_pt_vs_q_eta_" + chg + "_" + q_eta_suffix);
            if (it != tf1_map.end()) {
                // TF1 GUARD: clamp pt into the fit range (compiled/read-back TF1s return 0
                // outside [xmin,xmax] -- pp_trig_eff_highpt_jump lesson)
                const double x = std::min(std::max(static_cast<double>(pt), it->second->GetXmin()),
                                          it->second->GetXmax());
                val = it->second->Eval(x);
            }
        }
        if (val < 0.0)
            // No fitted turn-on for this muon. With the contiguous coarse binning + the gap cut
            // this is UNREACHABLE by construction, so reaching it means the sample and the fits
            // disagree (e.g. fits made with the fine binning, or the gap cut switched off).
            // THROW rather than fall back: a silent fallback here is what produced the w_trig = 0
            // pair-dropping bug on the data side (pp_trig_eff_highpt_jump.md).
            throw std::runtime_error("MCEffEvaluator: no fitted efficiency for chg=" + chg +
                                     Form(" q_eta=%.3f pt=%.2f", q_eta, pt) +
                                     " -- the fits and the selection disagree (coarse q_eta"
                                     " binning + fiducial gap cut should make this impossible)");
        if (val > 1.0) val = 1.0;    // cap: efficiency <= 1
        if (val < 0.02) { val = 0.02; ++n_floor; }  // mandated floor (counted)
        return val;
    }

    void PrintStats(const std::string& tag) const {
        std::cout << "MCEffEvaluator [" << tag << "]: " << n_eval << " evaluations, "
                  << n_fallback << " gap-q_eta 2D fallbacks (round 8: must be 0 -- the coarse "
                     "binning has no holes), "
                  << n_floor << " floor(0.02) firings ("
                  << (n_eval ? 100.0 * n_floor / n_eval : 0.0) << "%)" << std::endl;
    }
};

// ---------- DATA single-muon efficiency evaluation (CORRECTED-MC study only) ----------
// ε_data(pT, q·η) from the data tag-and-probe turn-on fits -- the SAME objects the analysis
// itself uses (RDFBasedHistFillingData.cxx:625-655 EvaluateSingleMuonEffcyPtFitted), looked up
// with the SAME per-(charge, fine q·η bin) key and the SAME unfitted-2D gap fallback:
//     TF1  f_pt2nd_vs_q_eta2nd<ctr>_<sign1|sign2>_2mu4_sepr_py_<qeta>_divided
//     TH2D h_pt2nd_vs_q_eta2nd<ctr>_<sign1|sign2>_2mu4_sepr_divided   (x = q·η, y = pT)
// with sign1 = μ⁺, sign2 = μ⁻ (data sign convention, mc_trigger_efficiency.md / _sub_mctrig_plots F3).
// Guards are MCEffEvaluator's, NOT the data class's: clamp the evaluation pT into the fit's
// [xmin, xmax] (a compiled TF1 read back from a file returns 0 outside its range -- the
// pp_trig_eff_highpt_jump lesson), cap at 1, floor at 0.02 (counted). Using the SAME guards on
// both sides is what makes SF = ε_data/ε_MC well behaved.
struct DataEffEvaluator {
    std::map<std::string, TF1*> tf1_map;
    std::map<std::string, TH2D*> ratio_map;
    std::string ctr;                   // "" or "_ctr0_5"
    CommonEffcyConfig cfg{};
    long long n_floor = 0, n_fallback = 0, n_eval = 0;

    void Load(const std::string& fit_file, const std::string& hist_file, const std::string& ctr_suffix) {
        ctr = ctr_suffix;
        TFile* ff = TFile::Open(fit_file.c_str(), "READ");
        if (!ff || ff->IsZombie())
            throw std::runtime_error("DataEffEvaluator: cannot open data fit file " + fit_file);
        TIter next(ff->GetListOfKeys());
        TKey* key;
        while ((key = static_cast<TKey*>(next())))
            if (std::string(key->GetClassName()) == "TF1")
                tf1_map[key->GetName()] = static_cast<TF1*>(key->ReadObj());

        TFile* fh = TFile::Open(hist_file.c_str(), "READ");
        if (!fh || fh->IsZombie())
            throw std::runtime_error("DataEffEvaluator: cannot open data hist file " + hist_file);
        TIter nexth(fh->GetListOfKeys());
        while ((key = static_cast<TKey*>(nexth()))) {
            const std::string nm = key->GetName();
            if (std::string(key->GetClassName()) == "TH2D" &&
                nm.rfind("h_pt2nd_vs_q_eta2nd", 0) == 0 &&
                nm.size() > 8 && nm.compare(nm.size() - 8, 8, "_divided") == 0) {
                TH2D* h = static_cast<TH2D*>(key->ReadObj());
                h->SetDirectory(nullptr);
                ratio_map[nm] = h;
            }
        }
        std::cout << "DataEffEvaluator: loaded " << tf1_map.size() << " TF1s from " << fit_file
                  << " + " << ratio_map.size() << " 2D fallbacks from " << hist_file
                  << " (ctr=\"" << ctr << "\")" << std::endl;
        if (tf1_map.empty())
            throw std::runtime_error("DataEffEvaluator: no TF1s in " + fit_file);
    }

    std::string FindQEtaSuffix(float q_eta) const {
        // ROUND 8: coarse (contiguous) binning, same as MCEffEvaluator -- see the note there.
        for (const auto& range : cfg.q_eta_proj_ranges_coarse_incl_gap)
            if (q_eta >= range.first && q_eta < range.second) return pairToSuffix(range);
        return "";
    }

    double Eval(float pt, float eta, int charge) {
        ++n_eval;
        const std::string sign = (charge > 0) ? "_sign1" : "_sign2";   // sign1 = mu+, sign2 = mu-
        const float q_eta = charge * eta;
        const std::string q_eta_suffix = FindQEtaSuffix(q_eta);

        double val = -1.0;
        if (!q_eta_suffix.empty()) {
            auto it = tf1_map.find("f_pt2nd_vs_q_eta2nd" + ctr + sign + "_2mu4_sepr_py_" +
                                   q_eta_suffix + "_divided");
            if (it != tf1_map.end()) {
                const double x = std::min(std::max(static_cast<double>(pt), it->second->GetXmin()),
                                          it->second->GetXmax());
                val = it->second->Eval(x);
            }
        }
        if (val < 0.0)
            // No 2D fallback (round 8) -- the coarse binning is contiguous and the gap cut
            // removes q*eta > 2.2, so this is unreachable unless fits and selection disagree.
            throw std::runtime_error("DataEffEvaluator: no fitted efficiency for" + sign +
                                     Form(" q_eta=%.3f pt=%.2f", q_eta, pt) +
                                     " -- fits and selection disagree");
        if (val > 1.0) val = 1.0;
        if (val < 0.02) { val = 0.02; ++n_floor; }
        return val;
    }

    void PrintStats(const std::string& tag) const {
        std::cout << "DataEffEvaluator [" << tag << "]: " << n_eval << " evaluations, "
                  << n_fallback << " gap-q_eta 2D fallbacks (round 8: must be 0 -- the coarse "
                     "binning has no holes), "
                  << n_floor << " floor(0.02) firings ("
                  << (n_eval ? 100.0 * n_floor / n_eval : 0.0) << "%)" << std::endl;
    }
};

// ---------- per-muon data/MC scale factor SF = ε_data / ε_MC (CORRECTED-MC study) ----------
// Both efficiencies are evaluated CONTINUOUSLY from the fitted TF1s at the muon's EXACT
// (pT, q·η) -- never resample-to-nearest. Because each ε is already capped at 1 and floored at
// 0.02, SF is mathematically confined to [0.02, 50]; the extra cap below is a pathology guard,
// not a physics choice, and every firing is counted and printed (a large SF can only come from
// a fit that has gone wrong, and would silently distort the corrected efficiency).
struct SFEvaluator {
    MCEffEvaluator*   mc   = nullptr;
    DataEffEvaluator* data = nullptr;
    double sf_max = 20.0, sf_min = 0.01;
    // CLOSURE mode (sf_closure): force SF ≡ 1. The corrected chain then degenerates to the
    // NOMINAL one (ε_corr = fit to the untouched turn-ons = ε_MC, numerator weight
    // SF/ε_corr = 1/ε_MC), so its Step-3/Step-4 output must reproduce the nominal output
    // bin-by-bin. That is the end-to-end closure test of the corrected-mode code path -- and it
    // is the only way to run it WITHOUT overwriting any nominal artefact (its files carry the
    // "_corrected_sfclosure" tag).
    bool force_unity = false;
    long long n_eval = 0, n_cap_hi = 0, n_cap_lo = 0;
    double sum = 0., sum2 = 0., min_seen = 1e300, max_seen = -1e300;

    double Eval(float pt, float eta, int charge) {
        const double e_mc = mc->Eval(pt, eta, charge);
        const double e_da = data->Eval(pt, eta, charge);
        double sf = force_unity ? 1.0 : e_da / e_mc;
        ++n_eval;
        min_seen = std::min(min_seen, sf);
        max_seen = std::max(max_seen, sf);
        if (sf > sf_max) { sf = sf_max; ++n_cap_hi; }
        if (sf < sf_min) { sf = sf_min; ++n_cap_lo; }
        sum += sf; sum2 += sf * sf;
        return sf;
    }

    void PrintStats(const std::string& tag) const {
        const double m = n_eval ? sum / n_eval : 0.;
        const double v = n_eval ? std::max(0., sum2 / n_eval - m * m) : 0.;
        std::cout << "SFEvaluator [" << tag << "]: " << n_eval << " evaluations, "
                  << "<SF> = " << m << " +- " << std::sqrt(v)
                  << ", range [" << min_seen << ", " << max_seen << "], "
                  << n_cap_hi << " capped at " << sf_max << " ("
                  << (n_eval ? 100.0 * n_cap_hi / n_eval : 0.0) << "%), "
                  << n_cap_lo << " floored at " << sf_min << " ("
                  << (n_eval ? 100.0 * n_cap_lo / n_eval : 0.0) << "%)" << std::endl;
    }
};

// ---------- per-leg aliases on a pair tree ----------
// Pair structs have no dictionary: use leaf-style columns. Dotted names need Alias before JIT.
ROOT::RDF::RNode AliasLeg(ROOT::RDF::RNode node, int leg, const std::string& wp_col) {
    const std::string m = (leg == 1) ? "m1." : "m2.";  // this leg
    const std::string o = (leg == 1) ? "m2." : "m1.";  // partner leg
    return node.Alias("lg_pt",      m + "pt")
               .Alias("lg_eta",     m + "eta")
               .Alias("lg_charge",  m + "charge")
               .Alias("lg_wp",      m + wp_col)
               .Alias("lg_passmu4", m + "passmu4")   // full mu4 chain (L1 && HLT)
               .Alias("lg_pass_l1", m + "pass_l1")   // L1_MU3V RoI match (round-5 #3)
               .Alias("lg_truth_pt",  m + "truth_pt")   // truth fiducial (round 7)
               .Alias("lg_truth_eta", m + "truth_eta")
               .Alias("ot_pt",      o + "pt")
               .Alias("ot_eta",     o + "eta")
               .Alias("ot_wp",      o + wp_col)
               .Alias("ot_charge",  o + "charge")     // Step-4 leg-leg covariance term
               .Alias("ot_passmu4", o + "passmu4")    // Step-4 leg-leg covariance term
               .Alias("ot_truth_pt",  o + "truth_pt")
               .Alias("ot_truth_eta", o + "truth_eta");
}

} // namespace MCTrigEff

// =============================================================================
void FillMCTrigEffHists(const std::string& sample = "pp", bool do_step3 = false,
                        bool use_tight_wp = true, bool do_step4 = false,
                        bool do_sanity = false, bool corrected_mc = false,
                        bool sf_closure = false, bool book_kn_stats = false) {
    using namespace MCTrigEff;

    if ((int)do_step3 + (int)do_step4 + (int)do_sanity > 1)
        throw std::invalid_argument("FillMCTrigEffHists: do_step3, do_step4 and do_sanity are "
                                    "mutually exclusive (each writes its own output file)");
    if (corrected_mc && do_sanity)
        throw std::invalid_argument("FillMCTrigEffHists: corrected_mc is not defined for the "
                                    "sanity mode (§3.5 compares MC variants, not MC to data)");
    if (sf_closure && !corrected_mc)
        throw std::invalid_argument("FillMCTrigEffHists: sf_closure is a mode OF the corrected-MC "
                                    "path (SF forced to 1); it needs corrected_mc = true");
    if (book_kn_stats && !do_step3)
        throw std::invalid_argument("FillMCTrigEffHists: book_kn_stats (per-pTHat-slice "
                                    "statistics) is a Step-3 addition; it needs do_step3 = true");

    const SampleConfig cfg = GetSampleConfig(sample);
    const Binnings bins = MakeBinnings();
    // WP config (registry: Analysis/docs/muon_wp_registry.md): TIGHT nominal, Medium
    // reachable for the WP systematic. Medium outputs carry the _medium_wp suffix.
    const std::string wp_col = use_tight_wp ? "pass_tight" : "pass_medium";
    const std::string wp_suf = use_tight_wp ? "" : "_medium_wp";

    // CORRECTED-MC study (mc_trigger_efficiency.md round-7 contract item 5). Every MC muon that
    // FIRES the trigger carries the extra per-muon weight SF = ε_data/ε_MC(pT, q·η); the
    // DENOMINATOR is untouched. Then ε_corr = Σ_fired w·SF / Σ_all w = ε_MC·<SF> ≈ ε_data, i.e.
    // the corrected MC reproduces the DATA single-muon efficiency. It answers two questions:
    //   (1) do "data-driven ε + MC ΔR ratios" and "MC corrected to data" give the same
    //       single-muon efficiency?  (Step 1 with corrected numerator, vs data.)
    //   (2) does correcting the MC change the ΔR CORRECTIONS?  (Steps 3/4 with numerator weight
    //       w·SF/ε_corr, ε_corr from the fit to the CORRECTED turn-ons -- self-consistency, the
    //       §3.3/§3.4 requirement that the inverse weight uses the efficiency of the sample being
    //       inverse-weighted.) If ε_corr were exactly ε_data then SF/ε_corr = 1/ε_MC, the ORIGINAL
    //       weight, and the ΔR corrections would be identical bin-by-bin -- the validity statement
    //       that the ΔR correction is insensitive to the single-muon normalisation and therefore
    //       transfers from the MC world to the data world.
    // Every corrected artefact carries "_corrected" in its file name; NOTHING nominal is touched.
    const std::string corr_suf = corrected_mc ? (sf_closure ? "_corrected_sfclosure" : "_corrected")
                                              : "";

    std::cout << "FillMCTrigEffHists: sample=" << sample << " (" << cfg.label << ")"
              << ", do_step3=" << do_step3 << ", do_step4=" << do_step4
              << ", WP=" << (use_tight_wp ? "tight" : "medium")
              << ", corrected_mc=" << corrected_mc
              << std::endl;
    std::cout << "  pair file:    " << cfg.pair_file << std::endl;
    std::cout << "  singles file: " << cfg.singles_file << std::endl;

    // TRUTH FIDUCIAL (round 7, DEFAULT -- see header). Applied on TOP of the data-like reco
    // cuts, to the muon itself and, for pairs, to BOTH legs (the pair is two muons that each
    // satisfy the analysis muon definition). Column names differ per tree: bare on the
    // single-muon tree, lg_/ot_ aliases on the pair trees, m1_/m2_ aliases in Step 3.
    const std::string kTruthFidSingle = "truth_pt > 4 && fabs(truth_eta) < 2.4";
    const std::string kTruthFidLeg    = "lg_truth_pt > 4 && fabs(lg_truth_eta) < 2.4 && "
                                        "ot_truth_pt > 4 && fabs(ot_truth_eta) < 2.4";
    const std::string kTruthFidPair   = "m1_truth_pt > 4 && fabs(m1_truth_eta) < 2.4 && "
                                        "m2_truth_pt > 4 && fabs(m2_truth_eta) < 2.4";

    // Truth-reco pT-match threshold for the round-7 SANITY CHECK (do_sanity only; it is NOT
    // part of the nominal selection). Value + justification: mc_trigger_efficiency.md §3.5.
    // Defined in ../Utilities/MCTrigEffSanityCfg.h because the plot macro draws the NUMBER on
    // the canvas -- one definition, so the drawn value cannot drift from the applied cut.
    const double kPtMatchThr = kSanityPtMatchThr;

    // FORWARD LOW-pT VETO (round 7, user; §3.5 decision rule). The §3.5 sanity check RULED OUT
    // "bad muons" as the cause of the saturated forward-negative MC turn-on: in q·η ∈ (−2.4,−2.0)
    // the MC efficiency is flat at ~0.90 from the very first pT bin (data rises 0.45 → 0.92), and
    // requiring one reconstructed vertex and |ΔpT|/pT^truth < 0.10 moves it by ≤ 0.006 -- while the
    // MIRROR bin q·η ∈ (2.0,2.2) in the same sample shows a perfectly normal turn-on 0.26 → 0.95.
    // The anomaly is therefore a real property of the r16578 trigger configuration (R8/R10), so
    // the affected muons are REMOVED from the ΔR-correlation measurement.
    //
    // The veto is ASYMMETRIC and applies ONLY to the low-pT forward-negative corner: a muon is
    // rejected iff it is SIMULTANEOUSLY pT < 7 GeV AND q·η < −2. High-pT forward muons are kept
    // (the saturation is a turn-on-region effect), and the positive-q·η side is untouched.
    //
    // SCOPE: Steps 2, 3 and 4 (the ΔR-correlation measurement) ONLY. Step 1 keeps the full
    // acceptance -- it is the ε_MC(pT,q·η) map plus the data/MC validation, and the anomalous
    // region must stay visible there. The sanity check (do_sanity) is likewise unvetoed: vetoing
    // it would erase the very effect it exists to display. A PAIR is kept only if BOTH legs pass.
    // (Kept byte-identical to the nominal macro: the corrected-MC study is only meaningful if the
    //  corrected and the original ΔR corrections are measured on EXACTLY the same sample.)
    const bool kVetoFwdLowPt = true;
    const std::string kFwdVetoSingle = "(pt > 7 || charge * eta > -2)";
    const std::string kFwdVetoLeg    = "(lg_pt > 7 || lg_charge * lg_eta > -2) && "
                                       "(ot_pt > 7 || ot_charge * ot_eta > -2)";
    const std::string kFwdVetoPair   = "(m1_pt > 7 || m1_charge * m1_eta > -2) && "
                                       "(m2_pt > 7 || m2_charge * m2_eta > -2)";
    (void)kFwdVetoSingle;   // defined for symmetry with the nominal macro; Step 1 is UNvetoed

    // ---- q*eta FIDUCIAL GAP CUT (round 8, user; ParamsSet::single_mu_fiducial_gap_cuts) -------
    // Reject a muon whose q*eta falls in a detector-gap window; a PAIR needs BOTH legs to pass.
    // UNLIKE the forward veto above this applies to EVERY step INCLUDING Step 1 and the sanity
    // check: the gap muons are being removed from the analysis altogether, so the eps_MC(pT,q*eta)
    // map itself must be measured without them. The windows are built from ParamsSet so the
    // numbers are never retyped into a JIT string.
    // The matching cut on the DATA side is applied to the PROBE only (user decision) --
    // see ParamsSet.h. The two are consistent: eps^nc is a per-muon efficiency, and after this
    // cut it is only ever evaluated for muons outside the gaps.
    // NOMINAL: ON. The one exception is the FORWARD-EDGE DECISION STUDY
    // (plot_forward_qeta_edge_scan.cxx), which has to compare candidate upper edges 2.20 / 2.25 /
    // 2.30 / 2.40 and therefore needs the q*eta > 2.2 muons this cut removes. Set the environment
    // variable MCTRIGEFF_NO_GAPCUT=1 for that one pass; it writes to a DISTINCT `_nogapcut`
    // output so it can never be mistaken for, or overwrite, the nominal.
    const bool kApplyGapCut = (gSystem->Getenv("MCTRIGEFF_NO_GAPCUT") == nullptr);
    if (!kApplyGapCut)
        std::cout << "\n  ##### MCTRIGEFF_NO_GAPCUT set: fiducial gap cut DISABLED, writing a "
                     "_nogapcut output (forward-edge study only) #####\n" << std::endl;
    const std::string kGapSingle = ParamsSet::FiducialGapCutExpr("charge * eta");
    const std::string kGapLeg    = ParamsSet::FiducialGapCutExpr("lg_charge * lg_eta") + " && "
                                 + ParamsSet::FiducialGapCutExpr("ot_charge * ot_eta");
    const std::string kGapPair   = ParamsSet::FiducialGapCutExpr("m1_charge * m1_eta") + " && "
                                 + ParamsSet::FiducialGapCutExpr("m2_charge * m2_eta");

    // common selection = data-side muon definition (nominal WP + fiducial) + truth fiducial
    const std::string sel_single = wp_col + " && pt > 4 && fabs(eta) < 2.4 && " + kTruthFidSingle
                                 + (kApplyGapCut ? " && " + kGapSingle : std::string());
    // overlay: 0-5% centrality only (doc D2; test sample is b=0-5 fm)
    const std::string sel_single_full = cfg.is_overlay
        ? sel_single + " && ev_centrality >= 0 && ev_centrality < 5"
        : sel_single;

    // Steps 2 and 4 (both legs must pass the analysis muon definition + the forward low-pT veto)
    const std::string sel_pair_legs =
        "lg_wp && lg_pt > 4 && fabs(lg_eta) < 2.4 && "
        "ot_wp && ot_pt > 4 && fabs(ot_eta) < 2.4 && " + kTruthFidLeg +
        (kVetoFwdLowPt ? " && " + kFwdVetoLeg : std::string()) +
        (kApplyGapCut  ? " && " + kGapLeg    : std::string());
    const std::string sel_pair_full = cfg.is_overlay
        ? sel_pair_legs + " && avg_centrality >= 0 && avg_centrality < 5"
        : sel_pair_legs;

    const std::vector<std::pair<std::string, std::string>> charges = {
        {"muplus", "lg_charge > 0"}, {"muminus", "lg_charge < 0"}
    };
    const std::vector<std::pair<std::string, std::string>> dr_bins = {
        {"dr0_0_2",   "dr < 0.2"},
        {"dr0_2_1_0", "dr >= 0.2 && dr < 1.0"},
        {"dr1_0_inf", "dr >= 1.0"}
    };
    // SS (sign1) + OS (sign2) are summed into the SIGN-INTEGRATED histograms.
    // ROUND 9 (user request): Steps 3 and 4 ALSO book a per-sign copy of every histogram, so the
    // same-sign and opposite-sign dR corrections are measured and fitted independently and the
    // charge dependence is CHECKED rather than assumed. The sign-integrated histograms are
    // unchanged -- the per-sign ones are additional, never a replacement.
    //
    // Why the sign-integrated series is still nominal, per step:
    //  * STEP 4 (single-leg correction eps_single(dR)): charge blindness HOLDS as far as
    //    measured -- the trigger response is a per-muon detector property, blind to the pair
    //    charge product, and R25 finds the two signs agreeing at every dR (f(0) = 1.047 same
    //    sign vs 1.148 opposite sign). Summing them is then the statistically optimal choice and
    //    the sign-integrated correction is nominal on physics grounds.
    //  * STEP 3 (dR correction to the factorised eps1*eps2): charge blindness is REFUTED. R25
    //    (measured this same round, docs/tracking/mc_trigger_efficiency.md) finds
    //    same-sign/opposite-sign = 0.427 at dR = 0.025, and the same split survives in the RAW
    //    un-inverse-weighted joint probability numraw/denom -- so it is not an artefact of the
    //    1/(eps1 eps2) weight. The sign-integrated Step-3 series therefore remains nominal only
    //    PENDING AN OPEN USER DECISION on whether to adopt the per-sign corrections; it is NOT
    //    justified by charge blindness.
    // Tree -> sign convention (DimuonAlgCoreT.c:113, MuonPairMC.h:47): sign1 = SAME sign,
    // sign2 = OPPOSITE sign, split on the TRUTH charges.
    const std::vector<std::string> pair_trees = {"muon_pair_tree_sign1", "muon_pair_tree_sign2"};
    auto sign_prefix = [](const std::string& tree) -> std::string {
        if (tree == "muon_pair_tree_sign1") return "ss_";   // same sign
        if (tree == "muon_pair_tree_sign2") return "os_";   // opposite sign
        throw std::runtime_error("FillMCTrigEffHists: unknown pair tree '" + tree + "'");
    };

    HistAccumulator<TH1D> acc1D;
    HistAccumulator<TH2D> acc2D;
    HistAccumulator<TH3D> acc3D;   // Step-3/4 pair-eta dependence (round-5 #4 / round-6 §3.4)

    MCEffEvaluator* evaluator = nullptr;  // Step-3/4 only (inverse-weight ε source)

    // ---- corrected-MC machinery (heap: must outlive the LAZY RDF loops) ----
    // sf_eval : SF = ε_data/ε_MC, applied to every FIRED muon (Steps 1, 3, 4)
    // mc_eval_nom : the NOMINAL MC ε -- it is the Bernoulli probability that governs the
    //               conditional error terms (SF is an analysis weight, it does NOT change the
    //               probability that the trigger fired), so errB/covQ use it in BOTH modes.
    // `evaluator`  : the ε that the numerator is inverse-weighted BY -- nominal ε_MC in nominal
    //               mode, ε_corr (fit to the corrected turn-ons) in corrected mode.
    SFEvaluator*      sf_eval     = nullptr;
    MCEffEvaluator*   mc_eval_nom = nullptr;
    DataEffEvaluator* data_eval   = nullptr;
    if (corrected_mc) {
        mc_eval_nom = new MCEffEvaluator();
        mc_eval_nom->LoadFits(cfg.dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root");
        data_eval = new DataEffEvaluator();
        data_eval->Load(SubstWP(cfg.data_fit_file_tmpl, wp_suf),
                        SubstWP(cfg.data_hist_file_tmpl, wp_suf), cfg.data_ctr);
        sf_eval = new SFEvaluator();
        sf_eval->mc = mc_eval_nom;
        sf_eval->data = data_eval;
        sf_eval->force_unity = sf_closure;
    }

    // keep dataframes alive until merge
    std::vector<std::unique_ptr<ROOT::RDataFrame>> rdf_store;

    int booking_id = 0;
    auto uniq = [&booking_id](const std::string& base) {
        return base + "__b" + std::to_string(booking_id++);
    };

    if (do_sanity) {
        // =====================================================================
        // (S) SANITY CHECK (§3.5, round 7): is the MC >> data single-muon efficiency in the
        //     forward q·η bins caused by "bad" muons/events rather than by the simulation?
        //     Two extra requirements are imposed on the Step-1 sample, separately and together:
        //       (1) ONE-VERTEX events (n_vtx == 1): removes any residual pile-up muon that the
        //           d0 / z0 cuts did not kill. The count is of TRACK-BEARING vertices, so the
        //           skim's dummy beamspot vertex is excluded (NTP: MuonFullsimExtra::n_vtx).
        //           n_vtx < 0 means the skim carried no vertex branch -> the variant is invalid
        //           and is left EMPTY rather than silently passing everything.
        //       (2) TRUTH-RECO pT MATCH |truth_pT - pT| / truth_pT < kPtMatchThr: removes badly
        //           mismeasured muons while keeping essentially all well-measured ones.
        //     Everything is booked in ONE file so the four variants are guaranteed to come from
        //     the same event loop; the comparison overlay is variant-vs-variant (NOT vs data --
        //     data has neither truth nor an equivalent vertex requirement).
        // =====================================================================
        rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>("muon_tree", cfg.singles_file));
        ROOT::RDF::RNode ds = *rdf_store.back();
        {   // n_vtx is REQUIRED here. Without it the `n_vtx == 1` filter cannot be built, and
            // silently dropping the requirement would turn the sanity check into a no-op that
            // looks like it passed. Fail loudly instead.
            const auto cols = ds.GetColumnNames();
            if (std::find(cols.begin(), cols.end(), "n_vtx") == cols.end())
                throw std::runtime_error(
                    "FillMCTrigEffHists(do_sanity): the single-muon tree " + cfg.singles_file +
                    " has no n_vtx branch. Re-run the store_mc_trigger single-muon NTP with the "
                    "round-7 vertex propagation (PythiaFullSimExtras: vtx_ntrk -> n_vtx).");
        }
        ds = ds.Define("q_eta", "(float)(charge * eta)")
               .Define("w", "(double)ev_weight")
               .Define("one", "0.5")   // dummy fill value for the pass-fraction counters
               .Define("pt_match", "fabs(truth_pt - pt) / truth_pt")
               .Filter(sel_single_full, "singles selection");

        const std::string cut_vtx = "n_vtx == 1";
        const std::string cut_ptm = "pt_match < " + std::to_string(kPtMatchThr);
        const std::vector<std::pair<std::string, std::string>> variants = {
            {"orig", "true"},
            {"vtx",  cut_vtx},
            {"ptm",  cut_ptm},
            {"both", cut_vtx + " && " + cut_ptm}
        };

        for (const auto& [vname, vcut] : variants) {
            auto dv = ds.Filter(vcut, "variant " + vname);
            for (const auto& [chg, chg_cut] : charges) {
                auto dc = dv.Filter(chg == "muplus" ? "charge > 0" : "charge < 0", chg);
                auto book = [&](ROOT::RDF::RNode node, const std::string& nd) {
                    const std::string s = nd + "_" + chg + "_" + vname;
                    acc1D.add("h_sanity_pt_" + s,
                        node.Histo1D({uniq("h_sanity_pt_" + s).c_str(), ";p_{T} [GeV];entries",
                                      static_cast<int>(bins.pt.size()) - 1, bins.pt.data()}, "pt", "w"));
                    acc1D.add("h_sanity_q_eta_" + s,
                        node.Histo1D({uniq("h_sanity_q_eta_" + s).c_str(), ";q#eta;entries",
                                      static_cast<int>(bins.q_eta.size()) - 1, bins.q_eta.data()},
                                     "q_eta", "w"));
                    acc2D.add("h_sanity_pt_vs_q_eta_" + s,
                        node.Histo2D({uniq("h_sanity_pt_vs_q_eta_" + s).c_str(), ";q#eta;p_{T} [GeV]",
                                      static_cast<int>(bins.q_eta.size()) - 1, bins.q_eta.data(),
                                      static_cast<int>(bins.pt.size()) - 1, bins.pt.data()},
                                     "q_eta", "pt", "w"));
                };
                book(dc, "denom");                                  // NO trigger requirement (§4)
                book(dc.Filter("passmu4", chg + " mu4 " + vname), "num");
            }
            // pass fractions (weighted and raw) -- reported for the contract's "percentage of
            // MC muons passing each requirement"
            acc1D.add("h_sanity_count_" + vname,
                ds.Filter(vcut, "count " + vname)
                  .Histo1D({uniq("h_sanity_count_" + vname).c_str(), ";;entries", 1, 0., 1.},
                           "one", "w"));
            acc1D.add("h_sanity_rawcount_" + vname,
                ds.Filter(vcut, "rawcount " + vname)
                  .Histo1D({uniq("h_sanity_rawcount_" + vname).c_str(), ";;entries", 1, 0., 1.},
                           "one"));
        }
    } else if (do_step4) {
        // =====================================================================
        // (D) Step 4 (§3.4): single-leg ΔR correction ε_single(ΔR) via inverse weighting.
        //     Leg-level analog of Step 3: numerator = the leg's OWN mu4 match weighted
        //     1/ε_MC(pT,q·η); denominator = ALL selected legs (MC weight, no trigger req §4).
        //     Ratio vs ΔR = ε_single(ΔR); plateau-normalized downstream -> ε_ΔR^single(ΔR),
        //     which dresses the PbPb union LINEAR terms (§2). NOT needed for pp (2mu4 product
        //     absorbs the single-leg ΔR into ε_ΔR^2mu4) -- pp runs only to validate the
        //     machinery on the FULL sample. Both legs of a pair are probes (role-swap) and
        //     SS+OS are summed, exactly as Step 2/3. Selection = both legs pass the analysis
        //     muon definition (partner is the "other reco muon" that defines ΔR, §3.2/§3.4).
        // =====================================================================
        // Inverse-weight ε source. NOMINAL: the §3.1 MC turn-on fits. CORRECTED: the fits to the
        // CORRECTED turn-ons (§3.3/§3.4 require the inverse weight to use the efficiency OF THE
        // SAMPLE BEING INVERSE-WEIGHTED, and the sample being weighted here is the corrected one).
        evaluator = new MCEffEvaluator();  // heap: must outlive the lazy RDF loops
        evaluator->LoadFits(cfg.dir + "single_mu_effcy_pT_fit_mc" + corr_suf + wp_suf + ".root");

        for (const auto& tree : pair_trees) {
            rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree, cfg.pair_file));
            for (int leg = 1; leg <= 2; ++leg) {
                ROOT::RDF::RNode dl = AliasLeg(*rdf_store.back(), leg, wp_col);
                dl = dl.Filter(sel_pair_full, tree + Form(" step4 leg%d selection", leg));

                // Sign-integrated + per-sign copies, exactly as Step 3 (round 9). The per-sign
                // series exist so the charge-blindness of the single-leg trigger response can be
                // checked; the sign-integrated series stays the nominal one.
                const std::string sp4 = sign_prefix(tree);
                const std::vector<std::string> name_prefixes4 = {std::string("h_mc_single_dr_"),
                                                                 "h_mc_single_dr_" + sp4};
                auto book_step4 = [&](ROOT::RDF::RNode node, const std::string& nd,
                                      const std::string& wcol) {
                  for (const auto& hp : name_prefixes4) {
                    acc1D.add(hp + "zoom_" + nd,
                        node.Histo1D({uniq(hp + "zoom_" + nd).c_str(), ";#DeltaR;entries",
                                      static_cast<int>(bins.dr_zoom.size()) - 1, bins.dr_zoom.data()},
                                     "dr", wcol));
                    acc1D.add(hp + "full_" + nd,
                        node.Histo1D({uniq(hp + "full_" + nd).c_str(), ";#DeltaR;entries",
                                      static_cast<int>(bins.dr_full.size()) - 1, bins.dr_full.data()},
                                     "dr", wcol));
                    // pair-pT x pair-eta breakdown (plateau-stability systematic, mirrors Step 3 #4)
                    acc3D.add(hp + "zoom_vs_pt_eta_" + nd,
                        node.Histo3D({uniq(hp + "zoom_vs_pt_eta_" + nd).c_str(),
                                      ";#DeltaR;p_{T}^{pair} [GeV];#eta^{pair}",
                                      static_cast<int>(bins.dr_zoom.size()) - 1, bins.dr_zoom.data(),
                                      static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                      static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                     "dr", "pair_pt", "pair_eta", wcol));
                    acc3D.add(hp + "full_vs_pt_eta_" + nd,
                        node.Histo3D({uniq(hp + "full_vs_pt_eta_" + nd).c_str(),
                                      ";#DeltaR;p_{T}^{pair} [GeV];#eta^{pair}",
                                      static_cast<int>(bins.dr_full.size()) - 1, bins.dr_full.data(),
                                      static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                      static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                     "dr", "pair_pt", "pair_eta", wcol));
                  }
                };

                // Per-leg efficiency columns. `eps_*` = the ε the numerator is DIVIDED BY (nominal
                // ε_MC, or ε_corr in corrected mode); `sf_*` = the data/MC scale factor carried by
                // every FIRED muon (identically 1 in nominal mode); `epsn_*` = the NOMINAL MC ε,
                // which is the Bernoulli probability the conditional error terms need in BOTH
                // modes (SF re-weights, it does not change whether the trigger fired).
                auto add_leg_eps = [&](ROOT::RDF::RNode n, const std::string& tag,
                                       const std::string& cpt, const std::string& ceta,
                                       const std::string& cq) {
                    n = n.Define("eps_" + tag,
                                 [ev = evaluator](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                                 {cpt, ceta, cq});
                    if (corrected_mc) {
                        n = n.Define("sf_" + tag,
                                     [ev = sf_eval](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                                     {cpt, ceta, cq})
                             .Define("epsn_" + tag,
                                     [ev = mc_eval_nom](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                                     {cpt, ceta, cq});
                    } else {
                        n = n.Define("sf_" + tag, "1.0").Define("epsn_" + tag, "eps_" + tag);
                    }
                    return n;
                };

                // denominator: ALL selected legs, no trigger requirement (§4), MC weight
                book_step4(dl, "denom", "weight");

                // numerator: leg's own mu4 match, weight a = MC weight * SF / ε(pT_leg, q·η_leg)
                // (§3.4; SF ≡ 1 nominal, so this is the unchanged `weight / eps_lg`)
                ROOT::RDF::RNode dn =
                    add_leg_eps(dl.Filter("lg_passmu4", tree + Form(" step4 leg%d mu4", leg)),
                                "lg", "lg_pt", "lg_eta", "lg_charge")
                        .Define("w_inv_single", "weight * sf_lg / eps_lg");
                book_step4(dn, "num", "w_inv_single");

                // ---- CONDITIONAL (binomial-correct) ERROR TERMS (round 7) -------------------
                // The numerator is a re-weighted SUBSET of the denominator, so propagating
                // e_R = R sqrt((e_N/N)^2 + (e_D/D)^2) -- which assumes independence -- is wrong
                // and over-states the error by sqrt((1+eps)/(1-eps)) (2-3x here). Conditioning on
                // the MC sample (D fixed; only the Bernoulli trigger decisions t_i fluctuate),
                //     N = sum_i t_i a_i ,  a_i = w_i/eps_i ,  P(t_i=1) = p_i = eps_i * R
                //     Var(N) = sum_i a_i^2 p_i(1-p_i)  +  cross terms between the two legs of a pair
                // The diagonal part is estimated from the FIRED legs as
                //     sum_fired a_i^2 (1-p_i) = A - R*B ,   A = sum w^2/eps^2 , B = sum w^2/eps
                // Both legs of a pair land in the SAME dR bin and their trigger decisions are
                // correlated (that correlation is exactly what Step 3 measures), so the cross term
                //     2 sum_pairs a_1 a_2 (p_12 - p_1 p_2) = covP - R^2 * covQ
                //     covP = sum_{both fired} 2 w^2/(eps_1 eps_2)   (estimates 2 sum a1 a2 p_12)
                //     covQ = sum_{all pairs}  2 w^2                 (2 sum a1a2 p1p2 = 2 R^2 sum w^2)
                // is booked once per pair (leg == 1). Then Var = A - R*B + covP - R^2*covQ and
                // e_R = sqrt(Var)/D, evaluated bin-by-bin in the plot macro.
                //
                // CORRECTED-MC generalization (round-7 contract item 5). With the numerator weight
                // a_i = w_i*SF_i/eps_corr,i and the trigger probability STILL governed by the
                // nominal MC efficiency, p_i = epsn_i*R, the same algebra gives
                //     A = sum_fired a_i^2 ,   B = sum_fired a_i^2 * epsn_i
                //     covP = sum_{both fired} 2 a_1 a_2 ,  covQ = sum_all 2 a_1 a_2 epsn_1 epsn_2
                // which reduce EXACTLY to the nominal expressions when SF = 1 and eps = epsn.
                auto dnerr = dn.Define("w_errA_single", "w_inv_single*w_inv_single")
                               .Define("w_errB_single", "w_inv_single*w_inv_single*epsn_lg");
                book_step4(dnerr, "errA", "w_errA_single");
                book_step4(dnerr, "errB", "w_errB_single");

                if (leg == 1) {   // pair-level terms: book ONCE per pair, not once per leg
                    ROOT::RDF::RNode dc =
                        add_leg_eps(add_leg_eps(dl, "l1", "lg_pt", "lg_eta", "lg_charge"),
                                    "l2", "ot_pt", "ot_eta", "ot_charge")
                            .Define("a_l1", "weight * sf_l1 / eps_l1")
                            .Define("a_l2", "weight * sf_l2 / eps_l2");
                    book_step4(dc.Filter("lg_passmu4 && ot_passmu4", tree + " step4 both legs mu4")
                                 .Define("w_covP", "2.0*a_l1*a_l2"),
                               "covP", "w_covP");
                    book_step4(dc.Define("w_covQ", "2.0*a_l1*a_l2*epsn_l1*epsn_l2"),
                               "covQ", "w_covQ");
                }
            }
        }
    } else if (!do_step3) {
        // =====================================================================
        // (A) Step 1 (§3.1): singles tree, per charge, denom/num
        // =====================================================================
        rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>("muon_tree", cfg.singles_file));
        ROOT::RDF::RNode ds = *rdf_store.back();
        ds = ds.Define("q_eta", "(float)(charge * eta)")
               .Define("w", "(double)ev_weight")
               .Filter(sel_single_full, "singles selection");

        for (const auto& [chg, chg_cut] : charges) {
            auto dc = ds.Filter(chg == "muplus" ? "charge > 0" : "charge < 0", chg);
            auto book_singles = [&](ROOT::RDF::RNode node, const std::string& nd,
                                    const std::string& wcol) {
                acc1D.add("h_mc_pt_" + nd + "_" + chg,
                    node.Histo1D({uniq("h_mc_pt_" + nd + "_" + chg).c_str(), ";p_{T} [GeV];entries",
                                  static_cast<int>(bins.pt.size()) - 1, bins.pt.data()}, "pt", wcol));
                acc1D.add("h_mc_eta_" + nd + "_" + chg,
                    node.Histo1D({uniq("h_mc_eta_" + nd + "_" + chg).c_str(), ";#eta;entries",
                                  static_cast<int>(bins.eta.size()) - 1, bins.eta.data()}, "eta", wcol));
                acc1D.add("h_mc_phi_" + nd + "_" + chg,
                    node.Histo1D({uniq("h_mc_phi_" + nd + "_" + chg).c_str(), ";#phi;entries",
                                  static_cast<int>(bins.phi.size()) - 1, bins.phi.data()}, "phi", wcol));
                acc2D.add("h_mc_pt_vs_q_eta_" + nd + "_" + chg,
                    node.Histo2D({uniq("h_mc_pt_vs_q_eta_" + nd + "_" + chg).c_str(),
                                  ";q#eta;p_{T} [GeV]",
                                  static_cast<int>(bins.q_eta.size()) - 1, bins.q_eta.data(),
                                  static_cast<int>(bins.pt.size()) - 1, bins.pt.data()},
                                 "q_eta", "pt", wcol));
            };
            book_singles(dc, "denom", "w");                          // NO trigger requirement (§4)
            if (!corrected_mc) {
                book_singles(dc.Filter("passmu4", chg + " mu4"), "num", "w");   // full mu4 chain
                // L1/HLT split reference (round-5 #3), for the Step-2 inclusive-singles line:
                book_singles(dc.Filter("pass_l1", chg + " L1"), "numl1", "w");  // L1_MU3V RoI
                book_singles(dc.Filter("passmu4 && pass_l1", chg + " HLT|L1"), "numhlt", "w");
            } else {
                // CORRECTED MC: the numerator (and ONLY the numerator) carries the per-muon
                // SF = ε_data/ε_MC, so num/denom = ε_MC·<SF> ≈ ε_data. No L1-only / HLT|L1
                // corrected variants: SF is defined for the FULL mu4 chain (the data tag-and-probe
                // measures the chain), and there is no data L1-only reference to correct to.
                auto dnum = dc.Filter("passmu4", chg + " mu4")
                              .Define("sf", [ev = sf_eval](float pt, float eta, int q)
                                              { return ev->Eval(pt, eta, q); },
                                      {"pt", "eta", "charge"})
                              .Define("epsn", [ev = mc_eval_nom](float pt, float eta, int q)
                                                { return ev->Eval(pt, eta, q); },
                                      {"pt", "eta", "charge"})
                              .Define("w_sf", "w * sf");
                book_singles(dnum, "num", "w_sf");
                // CONDITIONAL (binomial-correct) ERROR TERMS for the CORRECTED Step-1 efficiency.
                // The corrected numerator is a RE-WEIGHTED SUBSET of the denominator, and -- unlike
                // the nominal one -- it can EXCEED the denominator in a bin where SF > 1, so a
                // Bayesian/binomial divide is not even defined (TGraphAsymmErrors::BayesDivide
                // rejects the pair and returns an EMPTY graph, which the fitter would then "fit"
                // to its initial parameters: a silent, badly wrong turn-on). With
                //     N = sum_i t_i a_i ,  a_i = w_i*SF_i ,  P(t_i = 1) = p_i = eps_MC,i
                //     Var(N) = sum a_i^2 p_i (1-p_i)  ->  estimated from the FIRED muons as
                //     A - B ,  A = sum_fired (w*SF)^2 ,  B = sum_fired (w*SF)^2 * eps_MC
                // and e_eff = sqrt(A-B)/D. (Unweighted limit w=SF=1, eps=eff: A-B = D*eff*(1-eff),
                // the binomial variance, as it must be.) No R factor here -- unlike Steps 3/4 the
                // per-muon trigger probability IS eps_MC, with no dR-correlation factor on top.
                auto dnerr = dnum.Define("w_errA1", "w_sf * w_sf")
                                 .Define("w_errB1", "w_sf * w_sf * epsn");
                book_singles(dnerr, "errA", "w_errA1");
                book_singles(dnerr, "errB", "w_errB1");
            }
        }

        // =====================================================================
        // (B) Step 2 (§3.2): pair trees, per leg / charge / dR bin, denom/num
        // SKIPPED in corrected mode: Step 2 is the ΔR-binned factorization cross-check, it is not
        // part of the corrected-MC study, and filling it here would put unweighted Step-2 hists
        // into a file whose name says "corrected".
        // =====================================================================
        const std::vector<std::string> step2_trees =
            corrected_mc ? std::vector<std::string>{} : pair_trees;
        for (const auto& tree : step2_trees) {
            rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree, cfg.pair_file));
            for (int leg = 1; leg <= 2; ++leg) {
                ROOT::RDF::RNode dl = AliasLeg(*rdf_store.back(), leg, wp_col);
                dl = dl.Define("lg_q_eta", "(float)(lg_charge * lg_eta)")
                       .Filter(sel_pair_full, tree + Form(" leg%d selection", leg));

                for (const auto& [chg, chg_cut] : charges) {
                    auto dcc = dl.Filter(chg_cut, chg);
                    for (const auto& [drb, dr_cut] : dr_bins) {
                        auto dd = dcc.Filter(dr_cut, drb);
                        auto book_leg = [&](ROOT::RDF::RNode node, const std::string& nd) {
                            const std::string suff = nd + "_" + chg + "_" + drb;
                            acc1D.add("h_mc_pair_pt_" + suff,
                                node.Histo1D({uniq("h_mc_pair_pt_" + suff).c_str(), ";p_{T} [GeV];entries",
                                              static_cast<int>(bins.pt.size()) - 1, bins.pt.data()},
                                             "lg_pt", "weight"));
                            acc1D.add("h_mc_pair_q_eta_" + suff,
                                node.Histo1D({uniq("h_mc_pair_q_eta_" + suff).c_str(), ";q#eta;entries",
                                              static_cast<int>(bins.q_eta.size()) - 1, bins.q_eta.data()},
                                             "lg_q_eta", "weight"));
                            acc2D.add("h_mc_pair_pt_vs_q_eta_" + suff,
                                node.Histo2D({uniq("h_mc_pair_pt_vs_q_eta_" + suff).c_str(),
                                              ";q#eta;p_{T} [GeV]",
                                              static_cast<int>(bins.q_eta.size()) - 1, bins.q_eta.data(),
                                              static_cast<int>(bins.pt.size()) - 1, bins.pt.data()},
                                             "lg_q_eta", "lg_pt", "weight"));
                        };
                        // Step-2 L1/HLT split (round-5 #3). Three per-leg efficiencies:
                        //   full mu4 chain : num / denom          (current)
                        //   L1            : numl1 / denom         = P[L1 RoI | offline]
                        //   HLT | L1      : numhlt / numl1        = P[full chain | L1 RoI]
                        // so eff(chain) = eff(L1) * eff(HLT|L1). numhlt = chain && L1 keeps it
                        // <= numl1 even if the ΔR windows let a chain match miss its L1 RoI.
                        // pass_l1 is 0 on pre-reskim NTUPs -> the L1/HLT hists are empty then.
                        book_leg(dd, "denom");                                       // NO trigger req (§4)
                        book_leg(dd.Filter("lg_passmu4", drb + " leg mu4"), "num");  // full mu4 chain
                        book_leg(dd.Filter("lg_pass_l1", drb + " leg L1"), "numl1"); // L1_MU3V RoI
                        book_leg(dd.Filter("lg_passmu4 && lg_pass_l1", drb + " leg HLT|L1"), "numhlt");
                    }
                }
            }
        }
    } else {
        // =====================================================================
        // (C) Step 3 (§3.3): inverse-weighted ε_ΔR inputs
        // =====================================================================
        // Inverse-weight ε source: nominal ε_MC, or ε_corr (fits to the CORRECTED turn-ons) in
        // corrected mode -- §3.3 requires the ε of the sample being inverse-weighted.
        evaluator = new MCEffEvaluator();  // heap: must outlive the lazy RDF loops
        evaluator->LoadFits(cfg.dir + "single_mu_effcy_pT_fit_mc" + corr_suf + wp_suf + ".root");

        // pp: pair fires 2mu4; overlay (PbPb cross term): both legs mu4-matched
        const std::string trig_cond = cfg.is_overlay ? "m1_passmu4 && m2_passmu4" : "pass2mu4";

        for (const auto& tree : pair_trees) {
            rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree, cfg.pair_file));
            ROOT::RDF::RNode dp = *rdf_store.back();
            dp = dp.Alias("m1_pt", "m1.pt").Alias("m1_eta", "m1.eta").Alias("m1_charge", "m1.charge")
                   .Alias("m1_wp", "m1." + wp_col).Alias("m1_passmu4", "m1.passmu4")
                   .Alias("m1_truth_pt", "m1.truth_pt").Alias("m1_truth_eta", "m1.truth_eta")
                   .Alias("m2_pt", "m2.pt").Alias("m2_eta", "m2.eta").Alias("m2_charge", "m2.charge")
                   .Alias("m2_wp", "m2." + wp_col).Alias("m2_passmu4", "m2.passmu4")
                   .Alias("m2_truth_pt", "m2.truth_pt").Alias("m2_truth_eta", "m2.truth_eta");

            std::string sel = "m1_wp && m1_pt > 4 && fabs(m1_eta) < 2.4 && "
                              "m2_wp && m2_pt > 4 && fabs(m2_eta) < 2.4 && " + kTruthFidPair;
            if (kVetoFwdLowPt) sel += " && " + kFwdVetoPair;   // Step 3 (round-7 forward veto)
            // Step 3 builds its OWN selection string and does NOT go through sel_pair_full, so the
            // gap cut has to be repeated here -- the one place it is easy to leave out.
            if (kApplyGapCut)  sel += " && " + kGapPair;
            if (cfg.is_overlay) sel += " && avg_centrality >= 0 && avg_centrality < 5";
            dp = dp.Filter(sel, tree + " step3 selection");

            const std::string sp = sign_prefix(tree);
            // Every Step-3 histogram is booked TWICE: once sign-integrated ("h_mc_dr_...", both
            // trees merged by HistAccumulator) and once with the per-sign prefix
            // ("h_mc_dr_ss_..." / "h_mc_dr_os_...", one tree each). Downstream stages select a
            // series purely by that prefix, so no consumer has to know about the trees.
            const std::vector<std::string> name_prefixes = {std::string("h_mc_dr_"),
                                                            "h_mc_dr_" + sp};
            auto book_step3 = [&](ROOT::RDF::RNode node, const std::string& nd, const std::string& wcol) {
              for (const auto& hp : name_prefixes) {
                acc1D.add(hp + "zoom_" + nd,
                    node.Histo1D({uniq(hp + "zoom_" + nd).c_str(), ";#DeltaR;entries",
                                  static_cast<int>(bins.dr_zoom.size()) - 1, bins.dr_zoom.data()},
                                 "dr", wcol));
                acc1D.add(hp + "full_" + nd,
                    node.Histo1D({uniq(hp + "full_" + nd).c_str(), ";#DeltaR;entries",
                                  static_cast<int>(bins.dr_full.size()) - 1, bins.dr_full.data()},
                                 "dr", wcol));
                // NOTE (round 7): the separate dR x FINE-pair-pT 2D that used to live here was
                // REMOVED. It was binned on pT_bins_120 (15 log bins 8-120) and the Step-3
                // pair-pT slices panel grouped it as 8-13.8/13.8-23.6/23.6-40.6/40.6-120 --
                // a SECOND, inconsistent pair-pT binning alongside the coarse
                // canonical ParamsSet::pair_pt_coarse_bins used by the 3D below (and by
                // Step 4, and by crossx). Every pair-pT view now projects the SAME 3D, so 1D/2D/3D
                // cannot disagree. Single source of truth: ParamsSet.h (see CLAUDE.md).
                // round-5 #4: dR x coarse pair-pT x coarse pair-eta, for the pair-eta dependence
                // of eps_dR. Zoom and full dR ranges; projected per (pair pT, pair eta) cell.
                acc3D.add(hp + "zoom_vs_pt_eta_" + nd,
                    node.Histo3D({uniq(hp + "zoom_vs_pt_eta_" + nd).c_str(),
                                  ";#DeltaR;p_{T}^{pair} [GeV];#eta^{pair}",
                                  static_cast<int>(bins.dr_zoom.size()) - 1, bins.dr_zoom.data(),
                                  static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                  static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                 "dr", "pair_pt", "pair_eta", wcol));
                acc3D.add(hp + "full_vs_pt_eta_" + nd,
                    node.Histo3D({uniq(hp + "full_vs_pt_eta_" + nd).c_str(),
                                  ";#DeltaR;p_{T}^{pair} [GeV];#eta^{pair}",
                                  static_cast<int>(bins.dr_full.size()) - 1, bins.dr_full.data(),
                                  static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                  static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                 "dr", "pair_pt", "pair_eta", wcol));
              }
            };

            // denominator: ALL selected pairs, no trigger requirement (§4), weight = MC weight
            book_step3(dp, "denom", "weight");

            // ---- STATISTICS BOOKKEEPING (round 9, user request) -------------------------
            // How much sample the (pair pT, pair eta) cells of the dR correction actually have,
            // per pair charge combination. Booked HERE, on the Step-3 denominator node, because
            // that node is exactly one entry per selected pair under exactly the selection the
            // correction is measured with -- a separate macro re-deriving the selection would be
            // free to drift from it. Three quantities per cell:
            //   count  = raw number of muon pairs (unweighted; the statistical sample size)
            //   sumw   = sum of the per-pair MC weight = sigma_slice * eps_filt * r_isospin / N_slice,
            //            i.e. the cross section of that cell, in **nb** (ParamsSet / ami_weights.md;
            //            AMI crossSection is nb, NOT pb -- x1000 to compare with pp data in pb).
            //   sumw2  = sum of weight^2, so the statistical error on sumw is sqrt(sumw2).
            {
                auto dstat = dp.Define("w2_pair", "weight * weight");
                acc2D.add("h_mc_paircount_vs_pt_eta_" + std::string(sp, 0, 2),
                    dstat.Histo2D({uniq("h_mc_paircount_vs_pt_eta_" + sp).c_str(),
                                   ";p_{T}^{pair} [GeV];#eta^{pair};muon pairs",
                                   static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                   static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                  "pair_pt", "pair_eta"));
                acc2D.add("h_mc_pairsumw_vs_pt_eta_" + std::string(sp, 0, 2),
                    dstat.Histo2D({uniq("h_mc_pairsumw_vs_pt_eta_" + sp).c_str(),
                                   ";p_{T}^{pair} [GeV];#eta^{pair};#sigma [nb]",
                                   static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                   static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                  "pair_pt", "pair_eta", "weight"));
                acc2D.add("h_mc_pairsumw2_vs_pt_eta_" + std::string(sp, 0, 2),
                    dstat.Histo2D({uniq("h_mc_pairsumw2_vs_pt_eta_" + sp).c_str(),
                                   ";p_{T}^{pair} [GeV];#eta^{pair};#Sigma w^{2} [nb^{2}]",
                                   static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                   static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                  "pair_pt", "pair_eta", "w2_pair"));
            }

            // numerator: trigger condition, weight = MC weight * SF1*SF2 / (eps1 * eps2)
            //   nominal   : SF ≡ 1, eps = the §3.1 MC fits           -> weight / (eps1 eps2) (§3.3)
            //   corrected : SF = ε_data/ε_MC per muon, eps = ε_corr  -> the corrected-MC study.
            // `epsn*` is the NOMINAL MC ε (the Bernoulli probability) used by the error terms.
            auto add_pair_eps = [&](ROOT::RDF::RNode n, const std::string& tag,
                                    const std::string& cpt, const std::string& ceta,
                                    const std::string& cq) {
                n = n.Define("eps" + tag,
                             [ev = evaluator](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                             {cpt, ceta, cq});
                if (corrected_mc) {
                    n = n.Define("sf" + tag,
                                 [ev = sf_eval](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                                 {cpt, ceta, cq})
                         .Define("epsn" + tag,
                                 [ev = mc_eval_nom](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                                 {cpt, ceta, cq});
                } else {
                    n = n.Define("sf" + tag, "1.0").Define("epsn" + tag, "eps" + tag);
                }
                return n;
            };
            ROOT::RDF::RNode dn =
                add_pair_eps(add_pair_eps(dp.Filter(trig_cond, tree + " step3 trigger"),
                                          "1", "m1_pt", "m1_eta", "m1_charge"),
                             "2", "m2_pt", "m2_eta", "m2_charge")
                    .Define("w_inv", "weight * sf1 * sf2 / (eps1 * eps2)");
            book_step3(dn, "num", "w_inv");

            // RAW joint trigger probability (round 9): the SAME trigger-passing node, weighted by
            // the plain MC weight instead of 1/(eps1 eps2). numraw/denom is therefore
            // P(both mu fire | dR) with NO eps_MC anywhere in it.
            // Why it exists: R25 found the same-sign and opposite-sign eps_dR differing by a factor
            // 2.3 in the first dR bin. Two explanations compete -- a real L1 close-by-RoI effect,
            // or an eps_MC mis-parameterisation that SQUARES in the 1/(eps1 eps2) weight for a
            // close SAME-sign pair (q*eta_1 ~ q*eta_2) while partly CANCELLING for a close
            // OPPOSITE-sign pair (q*eta_1 ~ -q*eta_2). The inverse-weighted numerator cannot tell
            // them apart; this one can, because the second explanation lives entirely in the weight.
            book_step3(dn, "numraw", "weight");

            // CONDITIONAL (binomial-correct) ERROR TERMS (round 7) -- see the Step-4 block for
            // the derivation. Here each entry is one PAIR (a single Bernoulli trial with
            // p = eps_1 eps_2 R), so only the diagonal part exists: Var = A - R*B, with
            // A = sum_fired w^2/(eps1 eps2)^2 and B = sum_fired w^2/(eps1 eps2), and
            // e_R = sqrt(Var)/D. (Pairs sharing a leg -- events with >2 selected muons -- are a
            // residual, sub-leading correlation that is NOT modelled; see the tracking doc.)
            // (CORRECTED mode: A = sum_fired a^2, B = sum_fired a^2 * epsn1*epsn2 with
            //  a = w*SF1*SF2/(eps1 eps2) -- reduces EXACTLY to the nominal forms when SF = 1.)
            auto dnerr = dn.Define("w_errA", "w_inv*w_inv")
                           .Define("w_errB", "w_inv*w_inv*epsn1*epsn2");
            book_step3(dnerr, "errA", "w_errA");
            book_step3(dnerr, "errB", "w_errB");
        }

        // ---- STATISTICS BOOKKEEPING, per-pTHat-slice (2026-09-03, user request) --------------
        // Same three quantities (count/sumw/sumw2) as the STATISTICS BOOKKEEPING block above,
        // booked IDENTICALLY (same Step-3 selection, same axes), but read from the per-kn trees
        // `muon_pair_tree_kin{N}_sign{M}` that the SAME mc_trig pair file already carries
        // (PythiaAlgCoreT::fill_kn_trees_fullsim -- FillMuonPairTreePythia fills BOTH the merged
        // and the per-kn tree for every pair, so this is the ntuple-processing OUTPUT, not a
        // re-derivation from raw NTUPs). This isolates ONE pT-hat slice's contribution to each
        // (pair pT, pair eta) cell so it can be scaled to a PROJECTED sample size without
        // touching the other slices, e.g. "what would this cell look like with N events in kn4
        // instead of the N_beam actually produced". Opt-in (book_kn_stats, default false):
        // byte-unchanged for every existing call site; only the pp_full rerun that asks for it
        // gets these extra keys, additively, in the SAME step3 output file.
        if (book_kn_stats) {
            for (int ikin : {4, 5}) {  // the two highest pT-hat slices (user's stats request)
                for (int ksign = 1; ksign <= 2; ++ksign) {
                    const std::string ktree = "muon_pair_tree_kin" + std::to_string(ikin)
                                             + "_sign" + std::to_string(ksign);
                    rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(ktree, cfg.pair_file));
                    ROOT::RDF::RNode dk = *rdf_store.back();
                    dk = dk.Alias("m1_pt", "m1.pt").Alias("m1_eta", "m1.eta").Alias("m1_charge", "m1.charge")
                           .Alias("m1_wp", "m1." + wp_col).Alias("m1_passmu4", "m1.passmu4")
                           .Alias("m1_truth_pt", "m1.truth_pt").Alias("m1_truth_eta", "m1.truth_eta")
                           .Alias("m2_pt", "m2.pt").Alias("m2_eta", "m2.eta").Alias("m2_charge", "m2.charge")
                           .Alias("m2_wp", "m2." + wp_col).Alias("m2_passmu4", "m2.passmu4")
                           .Alias("m2_truth_pt", "m2.truth_pt").Alias("m2_truth_eta", "m2.truth_eta");
                    std::string ksel = "m1_wp && m1_pt > 4 && fabs(m1_eta) < 2.4 && "
                                       "m2_wp && m2_pt > 4 && fabs(m2_eta) < 2.4 && " + kTruthFidPair;
                    if (kVetoFwdLowPt) ksel += " && " + kFwdVetoPair;
                    if (kApplyGapCut)  ksel += " && " + kGapPair;
                    if (cfg.is_overlay) ksel += " && avg_centrality >= 0 && avg_centrality < 5";
                    dk = dk.Filter(ksel, ktree + " step3 selection (kn-split statistics)")
                           .Define("w2_pair_kn", "weight * weight");
                    const std::string ktag = std::string(ksign == 1 ? "ss" : "os")
                                            + "_kn" + std::to_string(ikin);
                    acc2D.add("h_mc_paircount_vs_pt_eta_" + ktag,
                        dk.Histo2D({uniq("h_mc_paircount_vs_pt_eta_" + ktag).c_str(),
                                    ";p_{T}^{pair} [GeV];#eta^{pair};muon pairs",
                                    static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                    static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                   "pair_pt", "pair_eta"));
                    acc2D.add("h_mc_pairsumw_vs_pt_eta_" + ktag,
                        dk.Histo2D({uniq("h_mc_pairsumw_vs_pt_eta_" + ktag).c_str(),
                                    ";p_{T}^{pair} [GeV];#eta^{pair};#sigma [nb]",
                                    static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                    static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                   "pair_pt", "pair_eta", "weight"));
                    acc2D.add("h_mc_pairsumw2_vs_pt_eta_" + ktag,
                        dk.Histo2D({uniq("h_mc_pairsumw2_vs_pt_eta_" + ktag).c_str(),
                                    ";p_{T}^{pair} [GeV];#eta^{pair};#Sigma w^{2} [nb^{2}]",
                                    static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                    static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                   "pair_pt", "pair_eta", "w2_pair_kn"));
                }
            }
        }
    }

    // ---------- trigger loops + merge ----------
    auto hists1D = acc1D.merge();
    auto hists2D = acc2D.merge();
    auto hists3D = acc3D.merge();

    if (evaluator) evaluator->PrintStats(cfg.label + (corrected_mc ? " eps_corr" : " eps_MC"));
    // Guard statistics of the corrected-MC study: how often the ε floors / the q·η-gap 2D
    // fallback / the SF pathology cap fired. These are the ONLY sources (besides fit quality)
    // of a residual difference between the corrected and the original ΔR corrections, so they
    // must be reported, not silently absorbed.
    if (mc_eval_nom) mc_eval_nom->PrintStats(cfg.label + " eps_MC(nominal)");
    if (data_eval)   data_eval->PrintStats(cfg.label + " eps_data");
    if (sf_eval)     sf_eval->PrintStats(cfg.label);

    // ---------- write ----------
    const std::string gap_tag = kApplyGapCut ? "" : "_nogapcut";
    const std::string ptbin_tag = MCTrigEffPairPt::FileSuffix();
    const std::string out_name = cfg.dir + "mc_trig_eff_hists_" + cfg.label + wp_suf + corr_suf +
                                 gap_tag + ptbin_tag +
                                 (do_sanity ? "_sanity.root" : do_step4 ? "_step4.root"
                                  : do_step3 ? "_step3.root" : ".root");
    TFile fout(out_name.c_str(), "RECREATE");
    if (fout.IsZombie()) throw std::runtime_error("FillMCTrigEffHists: cannot open output " + out_name);
    for (auto& kv : hists1D) kv.second->Write(kv.first.c_str());
    for (auto& kv : hists2D) kv.second->Write(kv.first.c_str());
    for (auto& kv : hists3D) kv.second->Write(kv.first.c_str());
    fout.Close();
    std::cout << "FillMCTrigEffHists: wrote " << hists1D.size() << " TH1D + "
              << hists2D.size() << " TH2D + " << hists3D.size() << " TH3D to " << out_name << std::endl;

    // ---------- sanity printout ----------
    auto PrintDrRatio = [&](const std::string& base, const std::string& tag) {
        // eps(dR) = num/denom: large-dR plateau (unweighted avg over the published window,
        // MCTrigEffPlateauWindow.h) + small-dR values from the zoom hist.
        TH1D* hn = hists1D.at("h_mc_" + base + "_dr_full_num");
        TH1D* hd = hists1D.at("h_mc_" + base + "_dr_full_denom");
        double sn = 0, sd = 0;
        for (int i = 1; i <= hd->GetNbinsX(); ++i) {
            const double c = hd->GetBinCenter(i);
            if (c >= MCTrigEffPlateau::kLo && c <= MCTrigEffPlateau::kHi) {
                sn += hn->GetBinContent(i); sd += hd->GetBinContent(i);
            }
        }
        std::cout << "\n===== " << tag << " sanity: eps_dR = num/denom, sample=" << cfg.label
                  << " =====" << std::endl;
        std::cout << Form("  large-dR average (dR in [%g,%g], weighted): ",
                          MCTrigEffPlateau::kLo, MCTrigEffPlateau::kHi)
                  << (sd > 0 ? sn / sd : -1) << std::endl;
        TH1D* hzn = hists1D.at("h_mc_" + base + "_dr_zoom_num");
        TH1D* hzd = hists1D.at("h_mc_" + base + "_dr_zoom_denom");
        for (int i = 1; i <= hzd->GetNbinsX(); ++i) {
            const double d = hzd->GetBinContent(i);
            std::cout << "  dR [" << hzd->GetBinLowEdge(i) << ", " << hzd->GetBinLowEdge(i + 1)
                      << "): eps_dR = " << (d > 0 ? hzn->GetBinContent(i) / d : -1)
                      << "  (denom w = " << d << ")" << std::endl;
        }
    };

    if (do_sanity) {
        // Sanity check: report the pass fraction of each requirement (contract item 2).
        auto frac = [&](const std::string& v, const char* base) {
            return hists1D.at(std::string(base) + v)->Integral();
        };
        const double w_all = frac("orig", "h_sanity_count_"), n_all = frac("orig", "h_sanity_rawcount_");
        std::cout << "\n===== Step-1 SANITY CHECK pass fractions, sample=" << cfg.label
                  << " (WP=" << (use_tight_wp ? "tight" : "medium") << ") =====" << std::endl;
        std::cout << "  requirement                          weighted%    raw%      raw N" << std::endl;
        for (const char* v : {"orig", "vtx", "ptm", "both"}) {
            const double wv = frac(v, "h_sanity_count_"), nv = frac(v, "h_sanity_rawcount_");
            std::cout << "  " << std::setw(34) << std::left
                      << (std::string(v) == "orig" ? "baseline (round-7 selection)"
                        : std::string(v) == "vtx"  ? "+ n_vtx == 1"
                        : std::string(v) == "ptm"  ? "+ |dpT|/truth_pT < thr"
                                                   : "+ both")
                      << "  " << (w_all > 0 ? 100.0 * wv / w_all : -1.0)
                      << "     " << (n_all > 0 ? 100.0 * nv / n_all : -1.0)
                      << "     " << nv << std::endl;
        }
        std::cout << "  (threshold = " << kPtMatchThr << ")" << std::endl;
        auto eff = [&](const std::string& v, const std::string& chg) {
            TH1D* hn = hists1D.at("h_sanity_pt_num_" + chg + "_" + v);
            TH1D* hd = hists1D.at("h_sanity_pt_denom_" + chg + "_" + v);
            return hd->Integral() > 0 ? hn->Integral() / hd->Integral() : -1.0;
        };
        for (const char* v : {"orig", "vtx", "ptm", "both"})
            std::cout << "  eps(mu4) " << std::setw(6) << std::left << v
                      << " : mu+ " << eff(v, "muplus") << " , mu- " << eff(v, "muminus") << std::endl;
    } else if (do_step4) {
        // Step-4: single-leg ε_single(dR). plateau ~1 (both samples, up to fit offset) validates
        // the machinery; the small-dR RISE is R4's saturation (pp ~1.20, overlay ~1.34 vs plateau).
        PrintDrRatio("single", "Step-4");
    } else if (!do_step3) {
        std::cout << "\n===== Step-1 sanity: weighted P(mu4 | selection), sample=" << cfg.label
                  << " =====" << std::endl;
        for (const auto& chg : {std::string("muplus"), std::string("muminus")}) {
            TH1D* hn = hists1D.at("h_mc_pt_num_" + chg);
            TH1D* hd = hists1D.at("h_mc_pt_denom_" + chg);
            const double eff_int = hd->Integral() > 0 ? hn->Integral() / hd->Integral() : -1;
            auto eff_range = [&](double lo, double hi) {
                const int b1 = hd->FindBin(lo + 1e-6), b2 = hd->FindBin(hi - 1e-6);
                const double d = hd->Integral(b1, b2);
                return d > 0 ? hn->Integral(b1, b2) / d : -1.0;
            };
            std::cout << "  " << chg << ": integrated=" << eff_int
                      << " | pt 4-5: "  << eff_range(4, 5)
                      << " | pt 6-8: "  << eff_range(6, 8)
                      << " | pt 20-60: " << eff_range(20, 60) << std::endl;
        }
    } else {
        // eps_dR ratio diagnostics: plateau at large dR (avg over the published window,
        // MCTrigEffPlateauWindow.h -- same window as plot_mc_trig_eff.cxx) + small-dR values
        TH1D* hn = hists1D.at("h_mc_dr_full_num");
        TH1D* hd = hists1D.at("h_mc_dr_full_denom");
        double sn = 0, sd = 0;
        for (int i = 1; i <= hd->GetNbinsX(); ++i) {
            const double c = hd->GetBinCenter(i);
            if (c >= MCTrigEffPlateau::kLo && c <= MCTrigEffPlateau::kHi) {
                sn += hn->GetBinContent(i); sd += hd->GetBinContent(i);
            }
        }
        std::cout << "\n===== Step-3 sanity: eps_dR = num/denom, sample=" << cfg.label
                  << " =====" << std::endl;
        std::cout << Form("  large-dR average (dR in [%g,%g], weighted): ",
                          MCTrigEffPlateau::kLo, MCTrigEffPlateau::kHi)
                  << (sd > 0 ? sn / sd : -1) << std::endl;
        TH1D* hzn = hists1D.at("h_mc_dr_zoom_num");
        TH1D* hzd = hists1D.at("h_mc_dr_zoom_denom");
        for (int i = 1; i <= hzd->GetNbinsX(); ++i) {
            const double d = hzd->GetBinContent(i);
            std::cout << "  dR [" << hzd->GetBinLowEdge(i) << ", " << hzd->GetBinLowEdge(i + 1)
                      << "): eps_dR = " << (d > 0 ? hzn->GetBinContent(i) / d : -1)
                      << "  (denom w = " << d << ")" << std::endl;
        }
    }
}
