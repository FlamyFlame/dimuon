#ifndef PBPB_RECO_EFF_SCALED_PLACEHOLDER_EVALUATOR_H
#define PBPB_RECO_EFF_SCALED_PLACEHOLDER_EVALUATOR_H

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

#include "PairRecoEffRecoBinnedEvaluator.h"

// =================================================================================================
// PbPbRecoEffScaledPlaceholderEvaluator -- the Pb+Pb pair reco efficiency PLACEHOLDER, anchored to the pp24
// Run-3 correction (docs/tracking/crossx_raa_reco_binned_propagation_20261013.md, Physics Procedure §2, D4):
//
//   eps_PbPb(pair) = C_Run2(p, c) * eps^{reco,SF}_j
//   C_Run2(p, c)   = E_Run2(p, c) / E_Run2(p, 60-80 %),
//   E_Run2(p, c)   = <eps_Run2(pT_1, q*eta_1; c) * eps_Run2(pT_2, q*eta_2; c)> over the pp24 fullsim truth
//                    signal pairs T in pair-pT bin p (truth kinematics, MC weight)
//   eps^{reco,SF}_j: the pp24 3D bin-by-bin, SF-corrected reco efficiency of the reco (pT, eta^pair, m) cell j
//
// Run 2 supplies ONLY the centrality dependence (relative to its most peripheral interval, per pair-pT bin);
// level, pair-pT, eta^pair (gap acceptance) and m dependence and the data/MC SFs are pp24's. User decision D4
// 2026-10-13 ("anchor to pp24"): a Run-2-vs-Run-3 efficiency LEVEL difference must not enter R_AA.
// Stored as h_E_run2_<wp> (= C_Run2) and h_G_<wp> (= eps^{reco,SF}).
//
// Centrality c is binned on the Run-2 placeholder's own intervals {0,10,20,30,40,50,60,80} (the
// intervals EvaluateSingleMuonRecoEffPlaceholder maps onto; c >= 80 -> last interval, as there).
// Evaluated at the data pair's RECO kinematics. Guards: definition marker; canonical (pT, eta^pair, m)
// axes; clamping counted; empty G cell -> m-integrated G -> inclusive G; product FLOORED at 0.05
// (the floor the Run-2 placeholder product had); no cap at 1 (none needed: counted).
// =================================================================================================
namespace PbPbRecoEffScaledPlaceholderDefinition {
inline const char* Key() { return "pbpb_reco_eff_scaled_placeholder_definition"; }
inline const char* Value() {
    return "eps_PbPb = C_Run2(pT bin, centrality interval) * eps^{reco,SF}_j; C_Run2 = E_Run2(p,c)/E_Run2(p,60-80%), "
           "E_Run2 = <eps_Run2(1)*eps_Run2(2)> over pp24 fullsim truth signal pairs (truth kinematics) per pT bin; "
           "eps^{reco,SF}_j = pp24 fullsim 3D reco-binned SF-corrected (j = reco (pT, eta^pair, m)); anchored to pp24 "
           "(D4); docs/tracking/crossx_raa_reco_binned_propagation_20261013.md";
}
inline std::vector<double> CentralityEdges() { return {0., 10., 20., 30., 40., 50., 60., 80.}; }
// Representative centrality of each interval, passed to EvaluateSingleMuonRecoEffPlaceholder.
inline std::vector<int> CentralityRepresentatives() { return {5, 15, 25, 35, 45, 55, 70}; }
}  // namespace PbPbRecoEffScaledPlaceholderDefinition

struct PbPbRecoEffScaledPlaceholderEvaluator {

    static constexpr double kMinEff = 0.05;

    std::unique_ptr<TH2D> h_E;      // C_Run2(pT, centrality interval) = E_Run2 / E_Run2(60-80 %)
    std::unique_ptr<TH3D> h_G;      // eps^{reco,SF}(pT, eta^pair, m), pp24
    std::unique_ptr<TH2D> h_G2d;    // m-integrated fallback
    double G_incl = -1.0;
    std::string wp, path;
    std::atomic<long long> n_eval{0}, n_clamp{0}, n_3d{0}, n_2d{0}, n_incl{0}, n_none{0}, n_floor{0}, n_above1{0};

    void Load(const std::string& file_path, bool use_tight)
    {
        namespace D = PbPbRecoEffScaledPlaceholderDefinition;
        wp = use_tight ? "tight" : "medium";
        path = file_path;
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie())
            throw std::runtime_error("PbPbRecoEffScaledPlaceholderEvaluator: cannot open " + path + " -- build it with "
                "plotting_codes/reco_effcy/reco_binned/build_pbpb_reco_eff_scaled_placeholder.C");
        auto* def = dynamic_cast<TNamed*>(f->Get(D::Key()));
        if (!def || std::string(def->GetTitle()) != D::Value())
            throw std::runtime_error("PbPbRecoEffScaledPlaceholderEvaluator: " + path + " carries "
                + (def ? "a DIFFERENT" : "NO") + " definition marker -- rebuild it");
        auto get = [&](const std::string& n) {
            TObject* o = f->Get(n.c_str());
            if (!o) throw std::runtime_error("PbPbRecoEffScaledPlaceholderEvaluator: no " + n + " in " + path);
            return o;
        };
        h_E.reset(static_cast<TH2D*>(dynamic_cast<TH2D&>(*get("h_E_run2_" + wp)).Clone("h_E_run2_eval")));
        h_G.reset(static_cast<TH3D*>(dynamic_cast<TH3D&>(*get("h_G_" + wp)).Clone("h_G_eval")));
        h_G2d.reset(static_cast<TH2D*>(dynamic_cast<TH2D&>(*get("h_G_" + wp + "_minv_integrated")).Clone("h_G2d_eval")));
        h_E->SetDirectory(nullptr); h_G->SetDirectory(nullptr); h_G2d->SetDirectory(nullptr);
        G_incl = dynamic_cast<TParameter<double>&>(*get("G_inclusive_" + wp)).GetVal();

        using EV = PairRecoEffRecoBinnedEvaluator;
        EV::CheckAxis(h_G->GetXaxis(), EV::PtEdges(), "pair-pT (G)", path);
        EV::CheckAxis(h_G->GetYaxis(), EV::EtaEdges(), "pair-eta (G)", path);
        EV::CheckAxis(h_G->GetZaxis(), EV::MinvEdges(), "m_mumu (G)", path);
        EV::CheckAxis(h_G2d->GetXaxis(), EV::PtEdges(), "pair-pT (G 2D)", path);
        EV::CheckAxis(h_G2d->GetYaxis(), EV::EtaEdges(), "pair-eta (G 2D)", path);
        EV::CheckAxis(h_E->GetXaxis(), EV::PtEdges(), "pair-pT (E)", path);
        EV::CheckAxis(h_E->GetYaxis(), D::CentralityEdges(), "centrality (E)", path);
        for (int i = 1; i <= h_E->GetNbinsX(); ++i)
            for (int k = 1; k <= h_E->GetNbinsY(); ++k)
                if (!(h_E->GetBinContent(i, k) > 0.))
                    throw std::runtime_error("PbPbRecoEffScaledPlaceholderEvaluator: empty E_Run2 cell in " + path);
        std::cout << "PbPbRecoEffScaledPlaceholderEvaluator [" << wp << " WP]: loaded C_Run2 = E_Run2/E_Run2(60-80%) ("
                  << h_E->GetNbinsX() << " pT x " << h_E->GetNbinsY() << " centrality) and eps^{reco,SF} ("
                  << h_G->GetNbinsX() * h_G->GetNbinsY() * h_G->GetNbinsZ() << " cells, inclusive " << G_incl
                  << ") from " << path << std::endl;
    }

    // RECO pair kinematics + event centrality (>= 0; the crossx filters 0 <= c < 80 upstream).
    double Eval(double pair_pt, double pair_eta, double minv, double centrality)
    {
        ++n_eval;
        const TAxis *ax = h_G->GetXaxis(), *ay = h_G->GetYaxis(), *az = h_G->GetZaxis();
        double x = Clamp(pair_pt, ax), y = Clamp(pair_eta, ay), z = Clamp(minv, az);
        if (x != pair_pt || y != pair_eta || z != minv) ++n_clamp;
        const int ix = ax->FindFixBin(x), iy = ay->FindFixBin(y), iz = az->FindFixBin(z);
        const TAxis* ac = h_E->GetYaxis();
        const int ic = std::min(std::max(ac->FindFixBin(centrality), 1), ac->GetNbins());
        const double E = h_E->GetBinContent(ix, ic);
        double G = h_G->GetBinContent(ix, iy, iz);
        if (G > 0.) ++n_3d;
        else if ((G = h_G2d->GetBinContent(ix, iy)) > 0.) ++n_2d;
        else if ((G = G_incl) > 0.) ++n_incl;
        else { ++n_none; return -1.0; }
        double v = E * G;
        if (v < kMinEff) { ++n_floor; v = kMinEff; }
        if (v > 1.0) ++n_above1;
        return v;
    }

    void PrintStats() const
    {
        const long long N = n_eval;
        auto pct = [&](long long n) { return N ? 100.0 * n / N : 0.0; };
        std::cout << "PbPbRecoEffScaledPlaceholderEvaluator [" << wp << " WP]: " << N << " evaluations -- clamped "
                  << n_clamp << " (" << pct(n_clamp) << "%); G 3D " << n_3d << " (" << pct(n_3d) << "%), m-integrated "
                  << n_2d << " (" << pct(n_2d) << "%), inclusive " << n_incl << " (" << pct(n_incl) << "%), none "
                  << n_none << "; floored " << n_floor << " (" << pct(n_floor) << "%); eps > 1: " << n_above1
                  << " (" << pct(n_above1) << "%)" << std::endl;
    }

private:
    static double Clamp(double v, const TAxis* a)
    {
        if (v >= a->GetXmin() && v < a->GetXmax()) return v;
        const double lo = a->GetBinCenter(1), hi = a->GetBinCenter(a->GetNbins());
        return v < lo ? lo : (v > hi ? hi : v);
    }
};

#endif  // PBPB_RECO_EFF_SCALED_PLACEHOLDER_EVALUATOR_H
