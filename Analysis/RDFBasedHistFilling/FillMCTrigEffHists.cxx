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
// Binnings reuse the DATA conventions:
//   pt   : "pT_bins_single_muon" = pT_bins_8 + pT_bins_60 (RDFBasedHistFillingData.cxx:286-290)
//   q·η  : "eta_bins_trig_effcy" = ParamsSet::makeEtaTrigEffcyBinning(1) (ibid:294)
//   phi  : 128 uniform bins in [-pi, pi] (ibid:300-307)
//   eta  : 48 uniform bins in [-2.4, 2.4]
//   pair pT : "pair_pt_log" = pT_bins_120, 15 log bins 8–120 (var1D_pp.json:102-106)
// MC is ALWAYS weighted: ev_weight (singles) / weight (pairs).
//
// Usage (from Analysis/RDFBasedHistFilling/):
//   root -b -l -q 'FillMCTrigEffHists.cxx+("pp")'
//   root -b -l -q 'FillMCTrigEffHists.cxx+("overlay", true)'              // Step 3
//   root -b -l -q 'FillMCTrigEffHists.cxx+("overlay", false, true, true)' // Step 4
//
// Output: <sample dir>/mc_trig_eff_hists_<pp24|hijing_overlay_pbpb23>.root
//         (do_step3=true writes a SEPARATE ..._step3.root; do_step4=true a SEPARATE
//          ..._step4.root; neither touches the Step-1/2 file)
// =============================================================================

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TKey.h>
#include <TString.h>
#include <ROOT/RDataFrame.hxx>

using namespace std;

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
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
};

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
    return cfg;
}

// ---------- binnings (data conventions; see header comment for provenance) ----------
struct Binnings {
    std::vector<double> pt;        // pT_bins_single_muon (incl. the data's duplicated 8.0 edge)
    std::vector<double> q_eta;     // eta_bins_trig_effcy
    std::vector<double> phi;       // 128 uniform [-pi, pi]
    std::vector<double> eta;       // 48 uniform [-2.4, 2.4]
    std::vector<double> pair_pt;   // pT_bins_120 (15 log bins 8-120)
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

    b.pair_pt = pms.pT_bins_120;

    b.dr_zoom.resize(21);
    for (int i = 0; i <= 20; ++i) b.dr_zoom[i] = i * (1.0 / 20);
    b.dr_full.resize(24);
    for (int i = 0; i <= 23; ++i) b.dr_full[i] = i * (5.75 / 23);

    // round-5 #4: coarse pair-pT (crossx) x coarse pair-eta (crossx pair_eta bins)
    b.pair_pt_coarse = pms.pair_pt_coarse_bins;                 // {8,15,27,50,150}
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
    long long n_fallback = 0;   // gap-q·η fallback lookups
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

    std::string FindQEtaSuffix(float q_eta) const {
        for (const auto& range : cfg.q_eta_proj_ranges_fine_excl_gap)
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
        if (val < 0.0) {
            // gap q·η (or missing TF1): unfitted 2D ratio fallback
            ++n_fallback;
            auto it2d = ratio_map.find("h_mc_pt_vs_q_eta_ratio_" + chg);
            if (it2d != ratio_map.end()) {
                TH2D* h = it2d->second;
                const double x = std::min(std::max(static_cast<double>(q_eta),
                                                   h->GetXaxis()->GetXmin() + 1e-6),
                                          h->GetXaxis()->GetXmax() - 1e-6);
                const double y = std::min(std::max(static_cast<double>(pt),
                                                   h->GetYaxis()->GetXmin() + 1e-6),
                                          h->GetYaxis()->GetXmax() - 1e-6);
                val = h->GetBinContent(h->FindBin(x, y));
            }
        }
        if (val < 0.0)
            // neither TF1 nor fallback TH2D provided a value: a configuration error,
            // not a low-efficiency leg -- must not be silently absorbed into the floor
            throw std::runtime_error("MCEffEvaluator: no efficiency source for chg=" + chg +
                                     Form(" q_eta=%.3f pt=%.2f", q_eta, pt));
        if (val > 1.0) val = 1.0;    // cap: efficiency <= 1
        if (val < 0.02) { val = 0.02; ++n_floor; }  // mandated floor (counted)
        return val;
    }

    void PrintStats(const std::string& tag) const {
        std::cout << "MCEffEvaluator [" << tag << "]: " << n_eval << " evaluations, "
                  << n_fallback << " gap-q_eta 2D fallbacks ("
                  << (n_eval ? 100.0 * n_fallback / n_eval : 0.0) << "%), "
                  << n_floor << " floor(0.02) firings ("
                  << (n_eval ? 100.0 * n_floor / n_eval : 0.0) << "%)" << std::endl;
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
               .Alias("ot_pt",      o + "pt")
               .Alias("ot_eta",     o + "eta")
               .Alias("ot_wp",      o + wp_col);
}

} // namespace MCTrigEff

// =============================================================================
void FillMCTrigEffHists(const std::string& sample = "pp", bool do_step3 = false,
                        bool use_tight_wp = true, bool do_step4 = false) {
    using namespace MCTrigEff;

    if (do_step3 && do_step4)
        throw std::invalid_argument("FillMCTrigEffHists: do_step3 and do_step4 are mutually "
                                    "exclusive (each writes its own output file)");

    const SampleConfig cfg = GetSampleConfig(sample);
    const Binnings bins = MakeBinnings();
    // WP config (registry: Analysis/docs/muon_wp_registry.md): TIGHT nominal, Medium
    // reachable for the WP systematic. Medium outputs carry the _medium_wp suffix.
    const std::string wp_col = use_tight_wp ? "pass_tight" : "pass_medium";
    const std::string wp_suf = use_tight_wp ? "" : "_medium_wp";

    std::cout << "FillMCTrigEffHists: sample=" << sample << " (" << cfg.label << ")"
              << ", do_step3=" << do_step3 << ", do_step4=" << do_step4
              << ", WP=" << (use_tight_wp ? "tight" : "medium")
              << std::endl;
    std::cout << "  pair file:    " << cfg.pair_file << std::endl;
    std::cout << "  singles file: " << cfg.singles_file << std::endl;

    // common selection = data-side muon definition (nominal WP + fiducial)
    const std::string sel_single = wp_col + " && pt > 4 && fabs(eta) < 2.4";
    // overlay: 0-5% centrality only (doc D2; test sample is b=0-5 fm)
    const std::string sel_single_full = cfg.is_overlay
        ? sel_single + " && ev_centrality >= 0 && ev_centrality < 5"
        : sel_single;

    const std::string sel_pair_legs =
        "lg_wp && lg_pt > 4 && fabs(lg_eta) < 2.4 && "
        "ot_wp && ot_pt > 4 && fabs(ot_eta) < 2.4";
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
    // SS (sign1) + OS (sign2) are summed into the same histograms: the trigger response
    // is a per-muon detector property, blind to the pair charge product.
    const std::vector<std::string> pair_trees = {"muon_pair_tree_sign1", "muon_pair_tree_sign2"};

    HistAccumulator<TH1D> acc1D;
    HistAccumulator<TH2D> acc2D;
    HistAccumulator<TH3D> acc3D;   // Step-3/4 pair-eta dependence (round-5 #4 / round-6 §3.4)

    MCEffEvaluator* evaluator = nullptr;  // Step-3/4 only (inverse-weight ε source)

    // keep dataframes alive until merge
    std::vector<std::unique_ptr<ROOT::RDataFrame>> rdf_store;

    int booking_id = 0;
    auto uniq = [&booking_id](const std::string& base) {
        return base + "__b" + std::to_string(booking_id++);
    };

    if (do_step4) {
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
        evaluator = new MCEffEvaluator();  // heap: must outlive the lazy RDF loops
        evaluator->LoadFits(cfg.dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root");

        for (const auto& tree : pair_trees) {
            rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree, cfg.pair_file));
            for (int leg = 1; leg <= 2; ++leg) {
                ROOT::RDF::RNode dl = AliasLeg(*rdf_store.back(), leg, wp_col);
                dl = dl.Filter(sel_pair_full, tree + Form(" step4 leg%d selection", leg));

                auto book_step4 = [&](ROOT::RDF::RNode node, const std::string& nd,
                                      const std::string& wcol) {
                    acc1D.add("h_mc_single_dr_zoom_" + nd,
                        node.Histo1D({uniq("h_mc_single_dr_zoom_" + nd).c_str(), ";#DeltaR;entries",
                                      static_cast<int>(bins.dr_zoom.size()) - 1, bins.dr_zoom.data()},
                                     "dr", wcol));
                    acc1D.add("h_mc_single_dr_full_" + nd,
                        node.Histo1D({uniq("h_mc_single_dr_full_" + nd).c_str(), ";#DeltaR;entries",
                                      static_cast<int>(bins.dr_full.size()) - 1, bins.dr_full.data()},
                                     "dr", wcol));
                    // pair-pT x pair-eta breakdown (plateau-stability systematic, mirrors Step 3 #4)
                    acc3D.add("h_mc_single_dr_zoom_vs_pt_eta_" + nd,
                        node.Histo3D({uniq("h_mc_single_dr_zoom_vs_pt_eta_" + nd).c_str(),
                                      ";#DeltaR;p_{T}^{pair} [GeV];#eta^{pair}",
                                      static_cast<int>(bins.dr_zoom.size()) - 1, bins.dr_zoom.data(),
                                      static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                      static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                     "dr", "pair_pt", "pair_eta", wcol));
                    acc3D.add("h_mc_single_dr_full_vs_pt_eta_" + nd,
                        node.Histo3D({uniq("h_mc_single_dr_full_vs_pt_eta_" + nd).c_str(),
                                      ";#DeltaR;p_{T}^{pair} [GeV];#eta^{pair}",
                                      static_cast<int>(bins.dr_full.size()) - 1, bins.dr_full.data(),
                                      static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                      static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                     "dr", "pair_pt", "pair_eta", wcol));
                };

                // denominator: ALL selected legs, no trigger requirement (§4), MC weight
                book_step4(dl, "denom", "weight");

                // numerator: leg's own mu4 match, weight = MC weight / ε_MC(pT_leg, q·η_leg) (§3.4)
                auto dn = dl.Filter("lg_passmu4", tree + Form(" step4 leg%d mu4", leg))
                            .Define("eps_lg",
                                    [ev = evaluator](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                                    {"lg_pt", "lg_eta", "lg_charge"})
                            .Define("w_inv_single", "weight / eps_lg");
                book_step4(dn, "num", "w_inv_single");
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
            auto book_singles = [&](ROOT::RDF::RNode node, const std::string& nd) {
                acc1D.add("h_mc_pt_" + nd + "_" + chg,
                    node.Histo1D({uniq("h_mc_pt_" + nd + "_" + chg).c_str(), ";p_{T} [GeV];entries",
                                  static_cast<int>(bins.pt.size()) - 1, bins.pt.data()}, "pt", "w"));
                acc1D.add("h_mc_eta_" + nd + "_" + chg,
                    node.Histo1D({uniq("h_mc_eta_" + nd + "_" + chg).c_str(), ";#eta;entries",
                                  static_cast<int>(bins.eta.size()) - 1, bins.eta.data()}, "eta", "w"));
                acc1D.add("h_mc_phi_" + nd + "_" + chg,
                    node.Histo1D({uniq("h_mc_phi_" + nd + "_" + chg).c_str(), ";#phi;entries",
                                  static_cast<int>(bins.phi.size()) - 1, bins.phi.data()}, "phi", "w"));
                acc2D.add("h_mc_pt_vs_q_eta_" + nd + "_" + chg,
                    node.Histo2D({uniq("h_mc_pt_vs_q_eta_" + nd + "_" + chg).c_str(),
                                  ";q#eta;p_{T} [GeV]",
                                  static_cast<int>(bins.q_eta.size()) - 1, bins.q_eta.data(),
                                  static_cast<int>(bins.pt.size()) - 1, bins.pt.data()},
                                 "q_eta", "pt", "w"));
            };
            book_singles(dc, "denom");                                       // NO trigger requirement (§4)
            book_singles(dc.Filter("passmu4", chg + " mu4"), "num");         // full mu4 chain
            // L1/HLT split reference (round-5 #3), for the Step-2 inclusive-singles line:
            book_singles(dc.Filter("pass_l1", chg + " L1"), "numl1");        // L1_MU3V RoI
            book_singles(dc.Filter("passmu4 && pass_l1", chg + " HLT|L1"), "numhlt");
        }

        // =====================================================================
        // (B) Step 2 (§3.2): pair trees, per leg / charge / dR bin, denom/num
        // =====================================================================
        for (const auto& tree : pair_trees) {
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
        evaluator = new MCEffEvaluator();  // heap: must outlive the lazy RDF loops
        evaluator->LoadFits(cfg.dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root");

        // pp: pair fires 2mu4; overlay (PbPb cross term): both legs mu4-matched
        const std::string trig_cond = cfg.is_overlay ? "m1_passmu4 && m2_passmu4" : "pass2mu4";

        for (const auto& tree : pair_trees) {
            rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree, cfg.pair_file));
            ROOT::RDF::RNode dp = *rdf_store.back();
            dp = dp.Alias("m1_pt", "m1.pt").Alias("m1_eta", "m1.eta").Alias("m1_charge", "m1.charge")
                   .Alias("m1_wp", "m1." + wp_col).Alias("m1_passmu4", "m1.passmu4")
                   .Alias("m2_pt", "m2.pt").Alias("m2_eta", "m2.eta").Alias("m2_charge", "m2.charge")
                   .Alias("m2_wp", "m2." + wp_col).Alias("m2_passmu4", "m2.passmu4");

            std::string sel = "m1_wp && m1_pt > 4 && fabs(m1_eta) < 2.4 && "
                              "m2_wp && m2_pt > 4 && fabs(m2_eta) < 2.4";
            if (cfg.is_overlay) sel += " && avg_centrality >= 0 && avg_centrality < 5";
            dp = dp.Filter(sel, tree + " step3 selection");

            auto book_step3 = [&](ROOT::RDF::RNode node, const std::string& nd, const std::string& wcol) {
                acc1D.add("h_mc_dr_zoom_" + nd,
                    node.Histo1D({uniq("h_mc_dr_zoom_" + nd).c_str(), ";#DeltaR;entries",
                                  static_cast<int>(bins.dr_zoom.size()) - 1, bins.dr_zoom.data()},
                                 "dr", wcol));
                acc1D.add("h_mc_dr_full_" + nd,
                    node.Histo1D({uniq("h_mc_dr_full_" + nd).c_str(), ";#DeltaR;entries",
                                  static_cast<int>(bins.dr_full.size()) - 1, bins.dr_full.data()},
                                 "dr", wcol));
                acc2D.add("h_mc_dr_zoom_vs_pair_pt_" + nd,
                    node.Histo2D({uniq("h_mc_dr_zoom_vs_pair_pt_" + nd).c_str(),
                                  ";#DeltaR;p_{T}^{pair} [GeV]",
                                  static_cast<int>(bins.dr_zoom.size()) - 1, bins.dr_zoom.data(),
                                  static_cast<int>(bins.pair_pt.size()) - 1, bins.pair_pt.data()},
                                 "dr", "pair_pt", wcol));
                // round-5 #4: dR x coarse pair-pT x coarse pair-eta, for the pair-eta dependence
                // of eps_dR. Zoom and full dR ranges; projected per (pair pT, pair eta) cell.
                acc3D.add("h_mc_dr_zoom_vs_pt_eta_" + nd,
                    node.Histo3D({uniq("h_mc_dr_zoom_vs_pt_eta_" + nd).c_str(),
                                  ";#DeltaR;p_{T}^{pair} [GeV];#eta^{pair}",
                                  static_cast<int>(bins.dr_zoom.size()) - 1, bins.dr_zoom.data(),
                                  static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                  static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                 "dr", "pair_pt", "pair_eta", wcol));
                acc3D.add("h_mc_dr_full_vs_pt_eta_" + nd,
                    node.Histo3D({uniq("h_mc_dr_full_vs_pt_eta_" + nd).c_str(),
                                  ";#DeltaR;p_{T}^{pair} [GeV];#eta^{pair}",
                                  static_cast<int>(bins.dr_full.size()) - 1, bins.dr_full.data(),
                                  static_cast<int>(bins.pair_pt_coarse.size()) - 1, bins.pair_pt_coarse.data(),
                                  static_cast<int>(bins.pair_eta_coarse.size()) - 1, bins.pair_eta_coarse.data()},
                                 "dr", "pair_pt", "pair_eta", wcol));
            };

            // denominator: ALL selected pairs, no trigger requirement (§4), weight = MC weight
            book_step3(dp, "denom", "weight");

            // numerator: trigger condition, weight = MC weight / (eps1 * eps2), eps = MC fits (§3.3)
            auto dn = dp.Filter(trig_cond, tree + " step3 trigger")
                        .Define("eps1", [ev = evaluator](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                                {"m1_pt", "m1_eta", "m1_charge"})
                        .Define("eps2", [ev = evaluator](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                                {"m2_pt", "m2_eta", "m2_charge"})
                        .Define("w_inv", "weight / (eps1 * eps2)");
            book_step3(dn, "num", "w_inv");
        }
    }

    // ---------- trigger loops + merge ----------
    auto hists1D = acc1D.merge();
    auto hists2D = acc2D.merge();
    auto hists3D = acc3D.merge();

    if (evaluator) evaluator->PrintStats(cfg.label);

    // ---------- write ----------
    const std::string out_name = cfg.dir + "mc_trig_eff_hists_" + cfg.label + wp_suf +
                                 (do_step4 ? "_step4.root" : do_step3 ? "_step3.root" : ".root");
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
        // eps(dR) = num/denom: large-dR plateau (weighted avg over [1,4], the published window)
        // + small-dR values from the zoom hist.
        TH1D* hn = hists1D.at("h_mc_" + base + "_dr_full_num");
        TH1D* hd = hists1D.at("h_mc_" + base + "_dr_full_denom");
        double sn = 0, sd = 0;
        for (int i = 1; i <= hd->GetNbinsX(); ++i) {
            const double c = hd->GetBinCenter(i);
            if (c >= 1.0 && c <= 4.0) { sn += hn->GetBinContent(i); sd += hd->GetBinContent(i); }
        }
        std::cout << "\n===== " << tag << " sanity: eps_dR = num/denom, sample=" << cfg.label
                  << " =====" << std::endl;
        std::cout << "  large-dR average (dR in [1,4], weighted): "
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

    if (do_step4) {
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
        // eps_dR ratio diagnostics: plateau at large dR (avg over [1,4] -- same window as the published plateau in plot_mc_trig_eff.cxx) + small-dR values
        TH1D* hn = hists1D.at("h_mc_dr_full_num");
        TH1D* hd = hists1D.at("h_mc_dr_full_denom");
        double sn = 0, sd = 0;
        for (int i = 1; i <= hd->GetNbinsX(); ++i) {
            const double c = hd->GetBinCenter(i);
            if (c >= 1.0 && c <= 4.0) { sn += hn->GetBinContent(i); sd += hd->GetBinContent(i); }
        }
        std::cout << "\n===== Step-3 sanity: eps_dR = num/denom, sample=" << cfg.label
                  << " =====" << std::endl;
        std::cout << "  large-dR average (dR in [1,4], weighted): "
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
