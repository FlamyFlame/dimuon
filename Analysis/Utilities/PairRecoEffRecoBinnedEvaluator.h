#ifndef PAIR_RECO_EFF_RECO_BINNED_EVALUATOR_H
#define PAIR_RECO_EFF_RECO_BINNED_EVALUATOR_H

#include <atomic>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TParameter.h>

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../RDFBasedHistFilling/CommonEffcyConfig.h"

// =================================================================================================
// PairRecoEffRecoBinnedEvaluator -- the RECO-BINNED pair reconstruction efficiency of the 3D
// bin-by-bin correction in (pair pT, eta^pair, m_mumu)
// (docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md, Physics Procedure §2).
//
//   eps^reco_j = N[R_j] / N[T_j]          (= 1 / C_j, the bin-by-bin correction factor)
//
//   T_j : truth OS pair (any origin) in the truth signal region (truth mass window, truth pair pT;
//         NTuple truth fiducial muon pT > 4.5, |eta| < 2.4), TRUTH (pT, eta^pair, m) in cell j.
//         NO truth gap cut.
//   R_j : truth-matched OS pair, both muons reco-matched, passing the WP and the DATA reco signal
//         selection incl. the gap cuts, RECO (pT, eta^pair, m) in cell j. Its truth may be ANYWHERE,
//         including outside T (out-of-T feed-in; user decision 2026-10-11 (b): the data reco yield
//         contains it, so it sits in the numerator and N[R_j] / eps^reco_j returns N[T_j]).
//
// Numerator and denominator are binned in DIFFERENT kinematics (reco vs truth) on the SAME
// canonical axes, so this is a ratio of yields, not a binomial efficiency: eps^reco_j > 1 is
// possible (net feed-in) and is NOT capped (Physics Procedure §4.2; C_j < 1 is allowed).
//
// EVALUATED AT THE DATA PAIR'S RECO (pair pT, eta^pair, m_mumu) -- the binning it is defined in, so
// unlike PairRecoEffEvaluator (truth-binned, evaluated at reco as an approximation) this IS the
// bin-by-bin correction, exact up to the MC truth-shape (prior) dependence of C_j.
//
// GUARDS (mirroring PairRecoEffEvaluator, each counted, never silent):
//   * definition marker (kDefinitionKey/Value) -- a product built under another definition refuses
//     to load;
//   * all three axes are CHECKED against the canonical source (ParamsSet::pair_pt_coarse_bins,
//     CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, ParamsSet::SignalMinvRecoEffcyBins());
//     a mismatch THROWS;
//   * pT / eta^pair / m are CLAMPED into the measured range (counted). With the crossx signal cuts
//     only reco pT above the top pair-pT edge can occur;
//   * a cell without a measurement (content <= 0: N[T_j] = 0 or N[R_j] = 0) falls back to the
//     m-INTEGRATED (pT, eta^pair) value, then to the INCLUSIVE value; -1 only if even that is missing;
//   * FLOORED at kMinEff (counted). NO cap at 1.
// =================================================================================================
namespace PairRecoEffRecoBinnedDefinition {
inline const char* Key() { return "pair_reco_eff_reco_binned_definition"; }
inline const char* Value() {
    return "eps^reco_j = N[R_j]/N[T_j], j = (pair pT, eta^pair, m_mumu) cell on the canonical axes; "
           "T_j = truth OS pair (any origin), truth signal cuts (mass window, pair pT; no gap cut), TRUTH "
           "kinematics in j; R_j = both reco-matched + WP + reco signal cuts incl. gap cuts, RECO kinematics "
           "in j, truth anywhere (out-of-T feed-in INCLUDED, decision (b) 2026-10-11); no cap at 1; "
           "docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md";
}
// DATA-CORRECTED variant: every R_j pair weighted by the product of its two single-muon data/MC
// efficiency scale factors (m_i.eff_sf_<wp>), eps^{reco,SF}_j = sum_{R_j} w SF_1 SF_2 / N[T_j]
// (docs/tracking/crossx_raa_reco_binned_propagation_20261013.md §2). A distinct marker, so the
// MC-only and the SF-corrected products can never be swapped silently.
inline const char* ValueSF() {
    return "eps^{reco,SF}_j = sum_{R_j} w*SF_1*SF_2 / N[T_j] (SF_i = m_i.eff_sf_<wp>, skim MuonEfficiencyScaleFactors; "
           "pT < 5 GeV (no SF) -> MC-weighted mean filled SF of WP muons in [5,6) GeV at the muon eta (0.1 lookup), D5), "
           "j = (pair pT, eta^pair, m_mumu) cell on the canonical axes; T_j, R_j as the MC-only "
           "product (feed-in included, no cap); docs/tracking/crossx_raa_reco_binned_propagation_20261013.md";
}
}  // namespace PairRecoEffRecoBinnedDefinition

struct PairRecoEffRecoBinnedEvaluator {

    static constexpr double kMinEff = 0.05;   // same floor as PairRecoEffEvaluator

    std::unique_ptr<TH3D> h_eff;      // eps^reco(pair pT, eta^pair, m_mumu)
    std::unique_ptr<TH2D> h_eff2d;    // m-integrated fallback (ratio of summed N[R] / N[T])
    double eff_incl = -1.0;           // inclusive fallback
    std::string wp, path;
    // ATOMIC: Eval() runs inside RDF lambdas under ImplicitMT; only the census is written.
    std::atomic<long long> n_eval{0}, n_clamp_pt{0}, n_clamp_eta{0}, n_clamp_m{0};
    std::atomic<long long> n_3d{0}, n_2d{0}, n_incl{0}, n_none{0}, n_floor{0}, n_above1{0};

    bool with_sf = false;   // true = the SF-corrected product (ValueSF marker)

    void Load(const std::string& file_path, bool use_tight, bool use_sf = false)
    {
        wp   = use_tight ? "tight" : "medium";
        path = file_path;
        with_sf = use_sf;
        const std::string want = use_sf ? PairRecoEffRecoBinnedDefinition::ValueSF()
                                        : PairRecoEffRecoBinnedDefinition::Value();
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie())
            throw std::runtime_error("PairRecoEffRecoBinnedEvaluator: cannot open " + path + " -- build it with "
                "plotting_codes/reco_effcy/reco_binned/build_pp24_fullsim_pair_reco_eff_reco_binned.C");
        auto* def = dynamic_cast<TNamed*>(f->Get(PairRecoEffRecoBinnedDefinition::Key()));
        if (!def || std::string(def->GetTitle()) != want)
            throw std::runtime_error("PairRecoEffRecoBinnedEvaluator: " + path + " carries "
                + (def ? "a DIFFERENT" : "NO") + " '" + PairRecoEffRecoBinnedDefinition::Key()
                + "' marker -- rebuild it");

        const std::string key = "h_eps_reco_binned_" + wp;
        auto* h = dynamic_cast<TH3D*>(f->Get(key.c_str()));
        if (!h) throw std::runtime_error("PairRecoEffRecoBinnedEvaluator: no " + key + " in " + path);
        h_eff.reset(static_cast<TH3D*>(h->Clone("h_eps_reco_binned_eval")));
        h_eff->SetDirectory(nullptr);

        auto* h2 = dynamic_cast<TH2D*>(f->Get((key + "_minv_integrated").c_str()));
        if (!h2) throw std::runtime_error("PairRecoEffRecoBinnedEvaluator: no " + key
                                          + "_minv_integrated in " + path + " (fallback is REQUIRED)");
        h_eff2d.reset(static_cast<TH2D*>(h2->Clone("h_eps_reco_binned_eval2d")));
        h_eff2d->SetDirectory(nullptr);

        auto* p = dynamic_cast<TParameter<double>*>(f->Get(("eps_reco_binned_inclusive_" + wp).c_str()));
        if (!p) throw std::runtime_error("PairRecoEffRecoBinnedEvaluator: no eps_reco_binned_inclusive_" + wp
                                         + " in " + path);
        eff_incl = p->GetVal();

        CheckCanonicalBinning();

        int n_empty = 0, n_gt1 = 0;
        for (int ix = 1; ix <= h_eff->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= h_eff->GetNbinsY(); ++iy)
                for (int iz = 1; iz <= h_eff->GetNbinsZ(); ++iz) {
                    const double v = h_eff->GetBinContent(ix, iy, iz);
                    if (v <= 0.) ++n_empty;
                    if (v > 1.) ++n_gt1;
                }
        std::cout << "PairRecoEffRecoBinnedEvaluator [" << wp << " WP" << (with_sf ? ", SF-corrected" : ", MC only")
                  << "]: loaded eps^reco("
                  << h_eff->GetNbinsX() << " pT x " << h_eff->GetNbinsY() << " eta^pair x "
                  << h_eff->GetNbinsZ() << " m) = "
                  << h_eff->GetNbinsX() * h_eff->GetNbinsY() * h_eff->GetNbinsZ() << " cells, "
                  << n_empty << " without measurement (fallback), " << n_gt1
                  << " with eps > 1 (not capped), inclusive " << eff_incl << ", from " << path << std::endl;
    }

    // RECO kinematics of the pair. Returns eps^reco, or -1 if no level carries a measurement.
    double Eval(double pair_pt, double pair_eta, double minv)
    {
        ++n_eval;
        const TAxis *ax = h_eff->GetXaxis(), *ay = h_eff->GetYaxis(), *az = h_eff->GetZaxis();
        double x = pair_pt, y = pair_eta, z = minv;
        if (x < ax->GetXmin() || x >= ax->GetXmax()) { ++n_clamp_pt;  x = Clamp(x, ax); }
        if (y < ay->GetXmin() || y >= ay->GetXmax()) { ++n_clamp_eta; y = Clamp(y, ay); }
        if (z < az->GetXmin() || z >= az->GetXmax()) { ++n_clamp_m;   z = Clamp(z, az); }
        const int ix = ax->FindFixBin(x), iy = ay->FindFixBin(y), iz = az->FindFixBin(z);
        double v = h_eff->GetBinContent(ix, iy, iz);
        if (v > 0.) { ++n_3d; }
        else {
            v = h_eff2d->GetBinContent(ix, iy);
            if (v > 0.) { ++n_2d; }
            else {
                v = eff_incl;
                if (v > 0.) { ++n_incl; }
                else        { ++n_none; return -1.0; }
            }
        }
        if (v < kMinEff) { ++n_floor; return kMinEff; }
        if (v > 1.0) ++n_above1;   // counted, NOT capped
        return v;
    }

    void PrintStats() const
    {
        const long long N = n_eval;
        auto pct = [&](long long n) { return N ? 100.0 * n / N : 0.0; };
        std::cout << "PairRecoEffRecoBinnedEvaluator [" << wp << " WP]: " << N << " evaluations -- clamped: pT "
                  << n_clamp_pt << " (" << pct(n_clamp_pt) << "%), eta^pair " << n_clamp_eta << " ("
                  << pct(n_clamp_eta) << "%), m " << n_clamp_m << " (" << pct(n_clamp_m) << "%); 3D "
                  << n_3d << " (" << pct(n_3d) << "%), m-integrated fallback " << n_2d << " (" << pct(n_2d)
                  << "%), inclusive fallback " << n_incl << " (" << pct(n_incl) << "%), none " << n_none
                  << " (" << pct(n_none) << "%); floored at " << kMinEff << ": " << n_floor << " ("
                  << pct(n_floor) << "%); eps > 1 (uncapped): " << n_above1 << " (" << pct(n_above1) << "%)"
                  << std::endl;
    }

    // Canonical axes as edge vectors, shared with the builder / plotter.
    static std::vector<double> PtEdges()  { static const ParamsSet p{}; return p.pair_pt_coarse_bins; }
    static std::vector<double> EtaEdges()
    {
        static const CommonEffcyConfig c{};
        std::vector<double> e{c.pair_eta_proj_ranges_coarse_incl_gap.front().first};
        for (const auto& r : c.pair_eta_proj_ranges_coarse_incl_gap) e.push_back(r.second);
        return e;
    }
    static std::vector<double> MinvEdges() { return ParamsSet::SignalMinvRecoEffcyBins(); }

    static void CheckAxis(const TAxis* ax, const std::vector<double>& edges, const std::string& what,
                          const std::string& where)
    {
        if (ax->GetNbins() != (int)edges.size() - 1)
            throw std::runtime_error("PairRecoEffRecoBinned: " + what + " axis of " + where + " has "
                + std::to_string(ax->GetNbins()) + " bins, canonical " + std::to_string(edges.size() - 1));
        for (int i = 1; i <= ax->GetNbins() + 1; ++i)
            if (std::fabs(ax->GetBinLowEdge(i) - edges[i - 1]) > 1e-5)
                throw std::runtime_error("PairRecoEffRecoBinned: " + what + " edge " + std::to_string(i - 1)
                    + " of " + where + " is " + std::to_string(ax->GetBinLowEdge(i)) + ", canonical "
                    + std::to_string(edges[i - 1]));
    }

private:
    static double Clamp(double v, const TAxis* a)
    {
        const double lo = a->GetBinCenter(1), hi = a->GetBinCenter(a->GetNbins());
        return v < lo ? lo : (v > hi ? hi : v);
    }
    void CheckCanonicalBinning() const
    {
        CheckAxis(h_eff->GetXaxis(), PtEdges(), "pair-pT", path);
        CheckAxis(h_eff->GetYaxis(), EtaEdges(), "pair-eta", path);
        CheckAxis(h_eff->GetZaxis(), MinvEdges(), "m_mumu", path);
        CheckAxis(h_eff2d->GetXaxis(), PtEdges(), "pair-pT (2D)", path);
        CheckAxis(h_eff2d->GetYaxis(), EtaEdges(), "pair-eta (2D)", path);
    }
};

#endif  // PAIR_RECO_EFF_RECO_BINNED_EVALUATOR_H
