#ifndef SINGLE_MU_EFF_EVALUATOR_H
#define SINGLE_MU_EFF_EVALUATOR_H

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include <TF1.h>
#include <TFile.h>
#include <TKey.h>
#include <TString.h>

#include "../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "proj_range_to_suffix.cxx"

// =================================================================================================
// The single-muon mu4 efficiency eps(pT, q*eta), from EITHER of the two fitted turn-on sets the
// trigger chain produces. One struct, because the two differ ONLY in the lookup key -- the
// evaluation and all three guards below must be IDENTICAL, or a ratio of the two (which is what
// the MC closure and the corrected-MC study both form) would carry a guard difference instead of
// a physics difference.
//
//   kDataTagAndProbe  eps^nc, the DATA tag-and-probe turn-on the ANALYSIS applies in its per-pair
//                     trigger weight (mu4_trig_effcy_implementation.md;
//                     RDFBasedHistFillingData.cxx EvaluateSingleMuonEffcyPtFitted). Key:
//                     f_pt2nd_vs_q_eta2nd<ctr>_<sign1|sign2>_2mu4_sepr_py_<qeta>_divided
//                     (sign1 = mu+, sign2 = mu-, the data sign convention).
//   kMCDirect         eps_MC, the MC direct conditional probability P[mu4 | reco muon]
//                     (mc_trigger_efficiency.md §3.1; FitMCSinglesEffcy.cxx). Key:
//                     f_mc_pt_vs_q_eta_<muplus|muminus>_<qeta>
//
// THREE GUARDS, and each one is a scar:
//   * CLAMP pT into the fit's [xmin, xmax]. A TF1 read back from a file returns 0 outside its
//     stored range; that produced eps -> floor -> a x13-88 blow-up of the pp 2mu4 product weight
//     above ~58 GeV (docs/tracking/pp_trig_eff_highpt_jump.md).
//   * CAP at 1. An efficiency cannot exceed 1.
//   * FLOOR at 0.02, COUNTED. An unfloored eps -> 0 makes 1/eps diverge; a SILENT floor hides how
//     often it happened.
// The function is evaluated CONTINUOUSLY at the muon's exact pT -- never resample-to-nearest.
//
// NO 2D FALLBACK. Since 2026-08-04 the coarse q*eta binning is contiguous over [-2.4, 2.3) and the
// fiducial gap cut removes q*eta > 2.3, so every surviving muon lands in a fitted bin. Reaching
// the "no fitted turn-on" branch means the fits and the selection disagree -> THROW. A silent
// fallback there is what produced the w_trig = 0 pair-dropping bug on the data side.
//
// ⚠ MIRROR NOTICE. RDFBasedHistFilling/FillMCTrigEffHists.cxx carries the same lookup inline as
// `DataEffEvaluator` (~lines 398-478); this header was extracted from it verbatim. That file was
// deliberately not migrated onto this header at extraction time (a concurrent session was
// re-running it). Migrate it the next time it is touched.
// =================================================================================================
struct SingleMuEffEvaluator {
    enum class Src { kDataTagAndProbe, kMCDirect };

    Src src = Src::kDataTagAndProbe;
    std::string name;                // what the printout calls this efficiency
    std::map<std::string, TF1*> tf1_map;
    std::string ctr;                 // "" (pp) or "_ctr0_5" (PbPb 0-5%, doc D2); data key only
    CommonEffcyConfig cfg{};
    // ATOMIC (2026-09-17): Eval() runs inside RDF lambdas under ImplicitMT (the crossx SF path);
    // a plain ++ lost ~0.4 % of the counts and made the "0 floored / 0 capped" guards untrustworthy.
    std::atomic<long long> n_floor{0}, n_cap{0}, n_eval{0};

    void Load(Src source, const std::string& fit_file, const std::string& ctr_suffix = "")
    {
        src  = source;
        ctr  = ctr_suffix;
        name = (src == Src::kDataTagAndProbe) ? "eps^nc (data tag-and-probe)" : "eps_MC";
        TFile* ff = TFile::Open(fit_file.c_str(), "READ");
        if (!ff || ff->IsZombie())
            throw std::runtime_error("SingleMuEffEvaluator: cannot open " + fit_file);
        TIter next(ff->GetListOfKeys());
        TKey* key;
        while ((key = static_cast<TKey*>(next())))
            if (std::string(key->GetClassName()) == "TF1")
                tf1_map[key->GetName()] = static_cast<TF1*>(key->ReadObj());
        if (tf1_map.empty())
            throw std::runtime_error("SingleMuEffEvaluator: no TF1s in " + fit_file);
        std::cout << "SingleMuEffEvaluator [" << name << "]: loaded " << tf1_map.size()
                  << " TF1s from " << fit_file << " (ctr=\"" << ctr << "\")" << std::endl;
    }

    std::string FindQEtaSuffix(float q_eta) const
    {
        for (const auto& range : cfg.q_eta_proj_ranges_coarse_incl_gap)
            if (q_eta >= range.first && q_eta < range.second) return pairToSuffix(range);
        return "";
    }

    // The ONE place the two sets differ.
    std::string Key(int charge, const std::string& q_eta_suffix) const
    {
        if (src == Src::kDataTagAndProbe)
            return "f_pt2nd_vs_q_eta2nd" + ctr + (charge > 0 ? "_sign1" : "_sign2")
                 + "_2mu4_sepr_py_" + q_eta_suffix + "_divided";
        return std::string("f_mc_pt_vs_q_eta_") + (charge > 0 ? "muplus" : "muminus")
             + "_" + q_eta_suffix;
    }

    double Eval(float pt, float eta, int charge)
    {
        ++n_eval;
        const float q_eta = charge * eta;
        const std::string q_eta_suffix = FindQEtaSuffix(q_eta);

        double val = -1.0;
        if (!q_eta_suffix.empty()) {
            auto it = tf1_map.find(Key(charge, q_eta_suffix));
            if (it != tf1_map.end()) {
                const double x = std::min(std::max(static_cast<double>(pt), it->second->GetXmin()),
                                          it->second->GetXmax());
                val = it->second->Eval(x);
            }
        }
        if (val < 0.0)
            throw std::runtime_error("SingleMuEffEvaluator [" + name + "]: no fitted turn-on for "
                                     + Form("charge=%+d q_eta=%.3f pt=%.2f", charge, q_eta, pt) +
                                     " -- the fits and the selection disagree (the coarse"
                                     " q_eta binning + the fiducial gap cut should make this"
                                     " impossible)");
        if (val > 1.0)  { val = 1.0;  ++n_cap; }
        if (val < 0.02) { val = 0.02; ++n_floor; }
        return val;
    }

    void PrintStats(const std::string& tag) const
    {
        const long long ne = n_eval, nc = n_cap, nf = n_floor;
        std::cout << "SingleMuEffEvaluator [" << name << " / " << tag << "]: " << ne
                  << " evaluations, "
                  << nc << " capped at 1 (" << (ne ? 100.0 * nc / ne : 0.0) << "%), "
                  << nf << " floored at 0.02 ("
                  << (ne ? 100.0 * nf / ne : 0.0) << "%)" << std::endl;
    }
};

#endif  // SINGLE_MU_EFF_EVALUATOR_H
