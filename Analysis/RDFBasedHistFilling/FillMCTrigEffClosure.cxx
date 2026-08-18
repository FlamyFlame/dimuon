// =============================================================================
// FillMCTrigEffClosure.cxx
//
// MC CLOSURE of the pp24 2mu4 trigger correction
// (docs/tracking/mc_trig_eff_closure.md -- Physics Procedure §2, §3.1-§3.3;
//  executes mc_trigger_efficiency.md Remaining Work 8).
//
// THE TEST. MC is the only sample carrying BOTH an unbiased denominator (every event stored,
// `StoreAllEvents`) and the per-pair trigger decision. So the correction chain can be closed on
// itself: take the pairs that PASS 2mu4, weight each by the inverse of the per-pair trigger
// probability THE ANALYSIS ACTUALLY APPLIES, and check that the result reproduces the ALL-PAIRS
// (no trigger requirement) spectrum -- differentially, not just inclusively:
//
//       Sum_{pairs passing 2mu4}  w_MC / [ eps^nc_data(1) * eps^nc_data(2) * eps_dR(dR) ]
//   C = -------------------------------------------------------------------------------  = 1 ?
//                            Sum_{all pairs}  w_MC
//
// binned in (pair pT, pair eta). eps^nc_data is the DATA tag-and-probe single-muon turn-on and
// eps_dR the MC-derived dR correction -- i.e. the mixed data+MC configuration the analysis uses,
// which is the physics test. Weighting with eps_MC instead would close by construction (it is the
// very eps the Step-3 inverse weighting divided by) and would test only the code.
//
// THE eps_MC DIAGNOSTIC (not a plotted series). The same numerator is booked a second time with
// eps_MC in place of eps^nc_data. That variant closes at 1 by construction if -- and only if --
// the dR correction and its fit are self-consistent, because eps_MC is the very efficiency the
// Step-3 inverse weighting divided by. So the two together SEPARATE the two things a single
// closure number confounds:
//     closure(eps_MC)                       -> does the dR correction + its fit close?
//     closure(eps_data) / closure(eps_MC)   -> <r_1 * r_2>, r_j = eps_MC(leg j)/eps_data(leg j),
//                                              weighted by w/(eps_MC1 eps_MC2 eps_dR): the known
//                                              MC/data single-muon over-efficiency
//                                              (mc_trigger_efficiency.md R3/R8)
// It is the mean of the PRODUCT of the two legs' ratios, NOT <r^2>: reading a per-leg value as
// sqrt(<r_1 r_2>) is only valid where both legs sample the same q*eta regime. For an
// OPPOSITE-SIGN pair the two legs sit in MIRRORED q*eta bins (q*eta_1 ~ -q*eta_2), and the
// forward-negative MC turn-on is anomalous while its mirror is not (R8/R10/R14) -- so the
// square-root reading holds in the barrel and NOT in the endcap panels.
// Without it a reader of the figure cannot tell a broken dR correction from the MC trigger
// simulation being more efficient than the detector. It is written to the output file and printed
// per pair-eta bin; the FIGURE carries only the three series the analysis actually applies.
//
// SCOPE: pp / 2mu4 / OPPOSITE SIGN only (user). Pb+Pb needs the mu4 UNION weight and its Step-4
// single-leg correction, and its full overlay production is still in flight. Same sign is excluded
// because its no-plateau-correction fits fail in many cells and the analysis is opposite-sign.
//
// TWO SAMPLE VERSIONS, filled in ONE event loop:
//   "all_os"  every opposite-sign pair passing the MC trigger-efficiency pair selection
//   "signal"  the same, plus the DATA-LIKE single-b RECO signal cuts
// and TWO dR-fit forms per version (`expo` = nominal, `polyu_fixedRp` = backup).
//
// ERRORS. The numerator is a re-weighted SUBSET of the denominator, so the closure ratio needs the
// CONDITIONAL (binomial-correct) error, not TH1::Divide's independent propagation (which is too
// long by ~sqrt((1+eps)/(1-eps)) -- mc_trigger_efficiency.md R12). With a_i = w_i/p_i, p_i =
// eps1*eps2*eps_dR the PREDICTED per-pair trigger probability and pi_i the TRUE one,
//     Var(N) = Sum_all a^2 pi (1 - pi)  ->  estimated from the FIRED pairs as Sum_fired a^2 (1-pi)
// and pi_i is estimated to first order as R*p_i (R = the bin's closure ratio), giving
//     Var = A - R*B ,  A = Sum_fired a^2 ,  B = Sum_fired a^2 p ,  sigma_C = sqrt(Var)/D.
// A and B are booked HERE; the division and the R factor are applied at the plot stage by
// SetConditionalRatioErrors (dr_correction_ratio.h -- the repo's ONE implementation, which also
// carries the k=n boundary fallback). Note R != 1 by a wide margin in this measurement, so
// dropping it (Var = A - B) is not a small approximation: since B <= A it shifts the variance by
// (R-1)*B, OVERSTATING sigma where R > 1 (up to 1.32x measured) and UNDERSTATING it where R < 1
// (up to 2.46x) -- the second direction being the dangerous one, and common in the statistics-poor
// high-pair-pT cells.
//
// PROVENANCE: reads the ntuple-processing output (`*_mc_trig*.root`) ONLY, never raw NTUPs, and
// builds its selection from Utilities/MCTrigEffPairSelection.h so it cannot drift from the sample
// the dR correction was measured on.
//
// Usage (from Analysis/RDFBasedHistFilling/):
//   root -l -b -q 'FillMCTrigEffClosure.cxx+("pp_full", true)'                     // Tight
//   root -l -b -q 'FillMCTrigEffClosure.cxx+("pp_full", false)'                    // Medium
//   MCTRIGEFF_PAIRPT_4BIN=1 root -l -b -q 'FillMCTrigEffClosure.cxx+("pp_full", true)'
//
// Output: <sample dir>/mc_trig_eff_closure_<label><wp><ptbin>.root
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>

using namespace std;

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../Utilities/SingleMuEffEvaluator.h"
#include "../Utilities/MCTrigEffPairPtBinning.h"
#include "../Utilities/MCTrigEffPairSelection.h"
#include "../plotting_codes/trig_effcy/mc_based/dr_correction_apply.h"
#include "CommonEffcyConfig.h"

namespace MCTrigEffClosure {

// The dR-fit forms the closure is produced for. `expo` is the NOMINAL Step-3 form and
// `polyu_fixedRp` the backup; both fail in a minority of cells (mc_trigger_efficiency.md R26,
// OPEN), which is exactly why the closure is drawn for both rather than for one.
const std::vector<std::string> kMethods = {"expo", "polyu_fixedRp"};

// Sample versions. "signal" is a strict SUBSET of "all_os".
const std::vector<std::string> kVersions = {"all_os", "signal"};

// The dR-correction series: the closure sample is opposite sign, so the correction that belongs on
// it is the opposite-sign one (doc D1). The token is the repo's (mc_trigger_efficiency.md
// DrCorrSignText): "os" = opposite sign = muon_pair_tree_sign2.
const char* kSignSeries = "os";
const char* kPairTree   = "muon_pair_tree_sign2";   // opposite sign

// The PLATEAU MODE of the dR-correction fits this closure consumes: the UN-MERGED
// no-plateau-correction variant (doc 3.2). It is named here and passed EXPLICITLY to both the
// evaluator and the provenance stamp below, because `DrCorrectionEvaluator::Load` now takes the
// mode as a defaulted argument (2026-08-17): relying on that default while the stamp spells the
// token out would make the stamp lie the moment the default moved. The pp24 CROSSX application
// deliberately uses a DIFFERENT mode (DrCorrCrossxMode() = "nocorr_ptmerge", the last two pair-pT
// bins merged); this closure stays on the un-merged cells it was measured against.
const char* kDrPlateauMode = "nocorr";

// DATA tag-and-probe reference for pp (RDFBasedHistFillingPP.cxx:306/323 -- the erf_plus_log fit
// directory). {WP} is "" (Tight) or "_medium_wp": the data reference MUST be at the same working
// point as the MC (mc_trigger_efficiency.md §3.0(d)).
const std::string kDataPPFitTmpl =
    "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/trg_effcy_pT_fitting_to_erf_plus_log/"
    "single_mu_effcy_pT_fit{WP}.root";

std::string SubstWP(std::string s, const std::string& wp_suf)
{
    const std::string tok = "{WP}";
    const size_t p = s.find(tok);
    if (p == std::string::npos) throw std::runtime_error("SubstWP: no {WP} in " + s);
    return s.replace(p, tok.size(), wp_suf);
}

// Contiguous (lo,hi) ranges -> bin edges (same helper as FillMCTrigEffHists.cxx).
inline std::vector<double> RangesToEdges(const QEtaBinning& ranges)
{
    std::vector<double> e;
    for (size_t i = 0; i < ranges.size(); ++i) {
        if (i == 0) e.push_back(ranges[i].first);
        else if (std::fabs(ranges[i].first - ranges[i - 1].second) > 1e-6)
            throw std::runtime_error("RangesToEdges: non-contiguous ranges");
        e.push_back(ranges[i].second);
    }
    return e;
}

}  // namespace MCTrigEffClosure

// =============================================================================
void FillMCTrigEffClosure(const std::string& sample = "pp_full", bool use_tight_wp = true)
{
    using namespace MCTrigEffClosure;

    // ---------------------------------------------------------------- configuration
    // Sample identity comes from the SAME table the whole dR-correction chain reads, so the
    // closure cannot end up pairing one sample's tree with another sample's fits.
    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    if (cfg.key != "pp" && cfg.key != "pp_full")
        throw std::invalid_argument(
            "FillMCTrigEffClosure: pp only. The Pb+Pb weight is the mu4 UNION "
            "eps_dR^single*(eps1+eps2) - eps1*eps2*eps_dR^cross, not the 2mu4 product this macro "
            "closes, and the full HIJING overlay production is not available (doc Scope).");

    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_col  = use_tight_wp ? "pass_tight" : "pass_medium";
    const std::string wp_text = use_tight_wp ? "Tight" : "Medium";

    // The pair file is not in DrCorrSample (that table serves the fit chain), so it is built here
    // from the same directory + the ntuple-processing naming.
    const std::string pair_file = cfg.mc_dir
        + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig"
        + (cfg.key == "pp_full" ? "_full" : "") + ".root";

    std::cout << "\n================ FillMCTrigEffClosure: " << sample << " / " << wp_text
              << " muons / opposite sign ================\n"
              << "  pair file : " << pair_file << "\n"
              << "  pair pT   : " << MCTrigEffPairPt::Describe() << std::endl;

    // ---------------------------------------------------------------- binnings (never retyped)
    ParamsSet pms;
    static const CommonEffcyConfig ecfg{};
    const std::vector<double> pt_edges  = MCTrigEffPairPt::Edges(pms);
    const std::vector<double> eta_edges = RangesToEdges(ecfg.pair_eta_proj_ranges_coarse_incl_gap);
    const double pt_lo = pt_edges.front(), pt_hi = pt_edges.back();

    // ---------------------------------------------------------------- efficiency lookups
    // Heap: the lambdas below are captured by LAZY RDF nodes and must outlive the booking scope.
    auto* eps_data = new SingleMuEffEvaluator();
    eps_data->Load(SingleMuEffEvaluator::Src::kDataTagAndProbe,
                   SubstWP(kDataPPFitTmpl, wp_suf), "");   // pp: no centrality token
    // eps_MC, for the diagnostic variant only (see the header comment). SAME struct, hence the
    // same clamp/cap/floor guards -- a guard difference between the two would show up in their
    // ratio as if it were physics.
    auto* eps_mc = new SingleMuEffEvaluator();
    eps_mc->Load(SingleMuEffEvaluator::Src::kMCDirect,
                 cfg.mc_dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root");

    std::map<std::string, DrCorrectionEvaluator*> eps_dr;
    for (const auto& m : kMethods) {
        auto* e = new DrCorrectionEvaluator();
        e->Load(cfg, use_tight_wp, m, kSignSeries, kDrPlateauMode);
        // The correction cells and the histogram axes MUST be the same binning, or a pair would be
        // corrected by another cell's curve while every histogram still fills (CLAUDE.md §Binnings).
        if (e->h_fit_ok->GetNbinsX() != (int)pt_edges.size() - 1 ||
            e->h_fit_ok->GetNbinsY() != (int)eta_edges.size() - 1)
            throw std::runtime_error("FillMCTrigEffClosure: the " + m + " fit file has " +
                std::to_string(e->h_fit_ok->GetNbinsX()) + "x" +
                std::to_string(e->h_fit_ok->GetNbinsY()) + " cells but this run's binning is " +
                std::to_string(pt_edges.size() - 1) + "x" + std::to_string(eta_edges.size() - 1) +
                " -- 4-bin/8-bin mismatch (MCTRIGEFF_PAIRPT_4BIN)?");
        for (size_t i = 0; i < pt_edges.size(); ++i)
            if (std::fabs(e->h_fit_ok->GetXaxis()->GetBinLowEdge(i + 1) - pt_edges[i]) > 1e-6)
                throw std::runtime_error("FillMCTrigEffClosure: pair-pT edges differ between the "
                                         + m + " fit file and ParamsSet -- stale fit file?");
        for (size_t i = 0; i < eta_edges.size(); ++i)
            if (std::fabs(e->h_fit_ok->GetYaxis()->GetBinLowEdge(i + 1) - eta_edges[i]) > 1e-6)
                throw std::runtime_error("FillMCTrigEffClosure: pair-eta edges differ between the "
                                         + m + " fit file and CommonEffcyConfig -- stale fit file?");
        eps_dr[m] = e;
    }

    // ---------------------------------------------------------------- the pair sample
    ROOT::RDataFrame df(kPairTree, pair_file);
    ROOT::RDF::RNode d = df;
    d = d.Alias("m1_pt", "m1.pt").Alias("m1_eta", "m1.eta").Alias("m1_charge", "m1.charge")
         .Alias("m1_wp", "m1." + wp_col)
         .Alias("m1_truth_pt", "m1.truth_pt").Alias("m1_truth_eta", "m1.truth_eta")
         .Alias("m2_pt", "m2.pt").Alias("m2_eta", "m2.eta").Alias("m2_charge", "m2.charge")
         .Alias("m2_wp", "m2." + wp_col)
         .Alias("m2_truth_pt", "m2.truth_pt").Alias("m2_truth_eta", "m2.truth_eta");

    // EXACTLY the sample the dR correction was measured on (doc §3.1 / D3).
    const std::string base_sel = MCTrigEffPairSel::Step3PairSelection(true);
    d = d.Filter(base_sel, "MC trigger-efficiency pair selection");

    // Outside the cell grid no correction exists, so those pairs cannot be part of a closure test
    // of it. The fraction dropped is reported below.
    //
    // The pair-eta bound is NOT redundant with the per-leg |eta| < 2.4: pair eta is the
    // pseudorapidity of the SUM 4-vector, which is not bounded by the two legs' pseudorapidities
    // and does exceed 2.4 for a small number of forward collinear pairs (measured: 16 of 1.36 M
    // triggered pairs). Without this cut those pairs sit in the histogram overflow -- excluded
    // from every drawn panel -- while still being counted as "no correction available", which
    // reads like a lookup failure when it is simply a pair outside the measured region.
    auto d_all_pt = d;                       // kept only to count what the cell window removes
    d = d.Filter(Form("pair_pt >= %.10g && pair_pt < %.10g && "
                      "pair_eta >= %.10g && pair_eta < %.10g",
                      pt_lo, pt_hi, eta_edges.front(), eta_edges.back()),
                 "pair pT and pair eta inside the correction cells");

    std::map<std::string, ROOT::RDF::RNode> versions;
    versions.emplace("all_os", d);
    versions.emplace("signal", d.Filter(MCTrigEffPairSel::SingleBSignalCutsReco(),
                                        "data-like single-b signal cuts"));

    // ---------------------------------------------------------------- the per-pair weight
    // Evaluated ONCE, on the triggered node, before the version split ("signal" is a subset of
    // "all_os", so splitting first would evaluate every efficiency twice).
    ROOT::RDF::RNode dn = d.Filter("pass2mu4", "pair passes 2mu4");
    dn = dn.Define("eps1", [ev = eps_data](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                   {"m1_pt", "m1_eta", "m1_charge"})
           .Define("eps2", [ev = eps_data](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                   {"m2_pt", "m2_eta", "m2_charge"})
           .Define("epsmc1", [ev = eps_mc](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                   {"m1_pt", "m1_eta", "m1_charge"})
           .Define("epsmc2", [ev = eps_mc](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                   {"m2_pt", "m2_eta", "m2_charge"});
    for (const auto& m : kMethods) {
        dn = dn.Define("edr_" + m,
                       [ev = eps_dr.at(m)](float dr, float ppt, float peta) {
                           return ev->Eval(dr, ppt, peta);
                       },
                       {"dr", "pair_pt", "pair_eta"})
               // p = the per-pair trigger probability the analysis predicts; w = its inverse,
               // times the MC weight. NOTHING else enters -- no reco efficiency, no acceptance.
               // p is consumed as a Bernoulli probability by the conditional error
               // (B = sum a^2 p). eps1, eps2 <= 1 by construction but eps_dR is only capped at
               // kMaxCorr, so p CAN exceed 1 where the dR correction is large; counted here
               // rather than clamped, because clamping p would silently alter the weight too.
               .Define("p_" + m,  "eps1 * eps2 * edr_" + m)
               .Define("pgt1_" + m, "p_" + m + " > 1.0 ? 1.0 : 0.0")
               .Define("w_" + m,  "weight / p_" + m)
               .Define("wA_" + m, "w_" + m + " * w_" + m)
               .Define("wB_" + m, "w_" + m + " * w_" + m + " * p_" + m)
               // eps_MC diagnostic variant (header comment): identical in every respect except
               // which single-muon efficiency goes into the weight.
               .Define("pmc_" + m, "epsmc1 * epsmc2 * edr_" + m)
               .Define("wmc_" + m, "weight / pmc_" + m);
    }
    std::map<std::string, ROOT::RDF::RNode> trig_versions;
    trig_versions.emplace("all_os", dn);
    trig_versions.emplace("signal", dn.Filter(MCTrigEffPairSel::SingleBSignalCutsReco(),
                                              "data-like single-b signal cuts (triggered)"));

    // ---------------------------------------------------------------- booking
    const int npt = static_cast<int>(pt_edges.size()) - 1;
    const int net = static_cast<int>(eta_edges.size()) - 1;
    const std::string ttl = ";p_{T}^{pair} [GeV];#eta^{pair};#Sigma w";
    auto model = [&](const std::string& n) {
        return ROOT::RDF::TH2DModel(n.c_str(), ttl.c_str(), npt, pt_edges.data(),
                                    net, eta_edges.data());
    };

    std::map<std::string, ROOT::RDF::RResultPtr<TH2D>> books;
    for (const auto& v : kVersions) {
        const std::string pre = "h_closure_" + v + "_";
        books.emplace(pre + "den",
                      versions.at(v).Histo2D(model(pre + "den"), "pair_pt", "pair_eta", "weight"));
        // Uncorrected triggered yield: not a plotted series, but it is what makes the closure
        // interpretable (how big a correction is being asked for in this bin).
        books.emplace(pre + "numraw",
                      trig_versions.at(v).Histo2D(model(pre + "numraw"), "pair_pt", "pair_eta",
                                                  "weight"));
        for (const auto& m : kMethods) {
            books.emplace(pre + "num_" + m,
                          trig_versions.at(v).Histo2D(model(pre + "num_" + m),
                                                      "pair_pt", "pair_eta", "w_" + m));
            books.emplace(pre + "numA_" + m,
                          trig_versions.at(v).Histo2D(model(pre + "numA_" + m),
                                                      "pair_pt", "pair_eta", "wA_" + m));
            books.emplace(pre + "numB_" + m,
                          trig_versions.at(v).Histo2D(model(pre + "numB_" + m),
                                                      "pair_pt", "pair_eta", "wB_" + m));
            books.emplace(pre + "nummc_" + m,
                          trig_versions.at(v).Histo2D(model(pre + "nummc_" + m),
                                                      "pair_pt", "pair_eta", "wmc_" + m));
        }
    }
    auto n_sel     = d_all_pt.Count();
    auto n_in_cells = d.Count();
    std::map<std::string, ROOT::RDF::RResultPtr<double>> n_pgt1;
    for (const auto& m : kMethods)
        n_pgt1.emplace(m, dn.Sum<double>("pgt1_" + m));

    // ---------------------------------------------------------------- run + write
    const std::string out_name = cfg.mc_dir + "mc_trig_eff_closure_" + cfg.mc_label + wp_suf
                               + MCTrigEffPairPt::FileSuffix() + ".root";
    TFile fout(out_name.c_str(), "RECREATE");
    if (fout.IsZombie()) throw std::runtime_error("FillMCTrigEffClosure: cannot open " + out_name);
    for (auto& kv : books) kv.second->Write(kv.first.c_str());

    // Provenance travels WITH the histograms: which eps went into the weight, which fit series and
    // plateau mode, and the selection strings -- so a consumer never has to guess.
    // A consumer that opens ONLY this file must be able to see which fit files it was built from
    // (identity AND mtime -- a concurrent re-fit would otherwise be invisible) and how many cells
    // fell back to the TEMPORARY raw-bin placeholder.
    // GetPathInfo returns non-zero and leaves mt/sz UNINITIALIZED when the file is not there --
    // writing that garbage into the provenance would defeat the point of stamping it at all.
    auto stamp = [](const std::string& f) {
        Long_t id = 0, sz = 0, fl = 0, mt = 0;
        if (gSystem->GetPathInfo(f.c_str(), &id, &sz, &fl, &mt) != 0)
            return std::string(Form("%s (MISSING)", f.c_str()));
        return std::string(Form("%s (mtime %ld, %ld B)", f.c_str(), mt, sz));
    };
    std::string dr_prov;
    for (const auto& m : kMethods) {
        dr_prov += Form(" | %s: %s; cells fitted %d, RAW-BIN PLACEHOLDER %d, no correction %d"
                        " (of the rejected, %d had fit_ok=1 but an unusable baseline C)",
                        m.c_str(),
                        stamp(DrCorrFitFile(cfg, use_tight_wp, 3, m, kSignSeries, kDrPlateauMode)).c_str(),
                        eps_dr.at(m)->n_cells_fitted, eps_dr.at(m)->n_cells_raw,
                        eps_dr.at(m)->n_cells_dead, eps_dr.at(m)->n_cells_bad_C);
    }
    // EVERY efficiency input gets path+mtime+size, not just the dR fits: a concurrent session
    // rewrites the single-muon turn-ons too, and a consumer holding only this file has no other
    // way to tell which version it was weighted by.
    const std::string eps_prov =
        " | eps^nc_data file: " + stamp(SubstWP(kDataPPFitTmpl, wp_suf))
      + " | eps_MC file (diagnostic numerator): "
      + stamp(cfg.mc_dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root");
    TNamed("provenance",
           Form("MC closure of the pp 2mu4 trigger correction (docs/tracking/mc_trig_eff_closure.md)"
                " | sample=%s label=%s WP=%s | tree=%s (opposite sign)"
                " | weight = w_MC / [eps^nc_data(1) * eps^nc_data(2) * eps_dR(dR)]"
                " | eps_dR = Step-3 %s fit, no plateau correction,"
                " divided by its own fitted baseline C, applied for dR < %g only"
                " | pair pT binning: %s | base selection: %s | signal cuts: %s%s%s",
                sample.c_str(), cfg.mc_label.c_str(), wp_text.c_str(), kPairTree,
                kSignSeries,
                DrCorrectionEvaluator::kDrMax, MCTrigEffPairPt::Describe().c_str(),
                base_sel.c_str(), MCTrigEffPairSel::SingleBSignalCutsReco().c_str(),
                eps_prov.c_str(), dr_prov.c_str()))
        .Write();
    fout.Close();

    std::cout << "\nFillMCTrigEffClosure: wrote " << books.size() << " TH2D to " << out_name
              << std::endl;

    // ---------------------------------------------------------------- guard statistics
    eps_data->PrintStats(cfg.mc_label);
    eps_mc->PrintStats(cfg.mc_label);
    for (const auto& m : kMethods)
        std::cout << "  predicted per-pair probability p = eps1*eps2*eps_dR > 1 in "
                  << (long long)*n_pgt1.at(m) << " triggered pairs (" << m
                  << ") -- p is used as a Bernoulli probability by the conditional error"
                  << std::endl;
    for (const auto& m : kMethods) eps_dr.at(m)->PrintStats();
    std::cout << "  selected pairs: " << *n_sel << " ; inside the correction cells (pair pT ["
              << pt_lo << ", " << pt_hi << "), pair eta [" << eta_edges.front() << ", "
              << eta_edges.back() << ")): " << *n_in_cells << " ("
              << (*n_sel ? 100.0 * (*n_in_cells) / (*n_sel) : 0.0) << "%)" << std::endl;

    // ---------------------------------------------------------------- inclusive closure
    // The INCLUSIVE number is a weak test (it can close while the differential one does not --
    // that is the whole reason the plot is differential), but a wildly-off inclusive value means
    // something is broken before any plot is looked at.
    std::cout << "\n===== inclusive closure (integral over all cells), sample=" << cfg.mc_label
              << ", WP=" << wp_text << " =====" << std::endl;
    for (const auto& v : kVersions) {
        const std::string pre = "h_closure_" + v + "_";
        const double den = books.at(pre + "den")->Integral();
        const double raw = books.at(pre + "numraw")->Integral();
        std::cout << "  " << std::setw(7) << std::left << v
                  << " : uncorrected 2mu4 / all = " << (den > 0 ? raw / den : -1.);
        for (const auto& m : kMethods)
            std::cout << " | corrected(" << m << ")/all = "
                      << (den > 0 ? books.at(pre + "num_" + m)->Integral() / den : -1.)
                      << " [eps_MC variant " << (den > 0 ? books.at(pre + "nummc_" + m)->Integral() / den : -1.)
                      << "]";
        std::cout << std::endl;
    }

    // The eps_MC diagnostic, per pair-eta bin (header comment): which of the two effects a
    // non-unit closure comes from. The nominal (`eps^nc data`) column is the figure; the eps_MC
    // column is the dR correction alone; their ratio is <r_1*r_2> (see the header).
    std::cout << "\n===== closure per pair-eta bin, nominal method \"" << kMethods.front()
              << "\" =====\n"
              << "  (the ratio column is <r_1*r_2>, r_j = eps_MC/eps_data of leg j -- NOT <r^2>:\n"
                 "   sqrt() gives a per-leg value only where both legs sample the same q*eta regime)\n"
              << "  " << std::setw(16) << std::left << "eta^pair"
              << std::setw(14) << "eps^nc data" << std::setw(14) << "eps_MC"
              << "ratio = <r_1*r_2>" << std::endl;
    {
        const std::string pre = "h_closure_all_os_";
        auto* hd = books.at(pre + "den").GetPtr();
        auto* hn = books.at(pre + "num_" + kMethods.front()).GetPtr();
        auto* hm = books.at(pre + "nummc_" + kMethods.front()).GetPtr();
        for (int iz = 1; iz <= net; ++iz) {
            const double D = hd->Integral(1, npt, iz, iz);
            const double N = hn->Integral(1, npt, iz, iz);
            const double M = hm->Integral(1, npt, iz, iz);
            std::cout << "  " << std::setw(16) << std::left
                      << Form("[%.1f, %.1f)", eta_edges[iz - 1], eta_edges[iz])
                      << std::setw(14) << (D > 0 ? N / D : -1.)
                      << std::setw(14) << (D > 0 ? M / D : -1.)
                      << (M > 0 ? N / M : -1.) << std::endl;
        }
    }
    std::cout << "done." << std::endl;
}
