#ifndef PAIR_RECO_EFF_EVALUATOR_H
#define PAIR_RECO_EFF_EVALUATOR_H

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
#include <TParameter.h>

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../RDFBasedHistFilling/CommonEffcyConfig.h"

// =================================================================================================
// PairRecoEffEvaluator -- the single-b PAIR reconstruction efficiency
// eps_reco(pair pT, pair eta, dR), measured in the pp24-condition Pythia fullsim FULL sample and
// applied per pair in the cross-section.
//
// WHY A PAIR EFFICIENCY, NOT eps_1 * eps_2. The two signal muons come from ONE b-hadron chain, so
// they are close in angle: they share ID hits and compete in the ambiguity resolution, and the
// reconstruction of one is not independent of the other. The dR axis is exactly what carries that
// correlation (docs/analysis_overview.md 4b). This REPLACES the Run-2 single-muon product
// placeholder (docs/tracking/reco_eff_placeholder_run2.md), which had no dR dependence at all.
//
// DEFINITION of the histogram this reads (built by
// plotting_codes/reco_effcy/build_pp24_fullsim_pair_reco_eff.C from the RDF hist-filling output):
//
//   eps_reco(cell) =  N[ single-b OS pair, both muons truth-matched to reco, pair passes the WP,
//                        and the RECO pair passes the signal region incl. the fiducial gap cut ]
//                  /  N[ single-b OS pair whose TRUTH pair passes the signal region
//                        incl. the fiducial gap cut on TRUTH q*eta ]
//
// binned in TRUTH kinematics and weighted by the MC event weight. Because the gap cut sits in BOTH
// legs, eps_reco is a FIDUCIAL efficiency and does NOT contain the truth-level gap acceptance
// eps_acc (muon_gap_cuts_acceptance.md F12) -- that stays a separate, not-yet-applied factor.
//
// EVALUATED AT RECO KINEMATICS. The data pair's own (pair pT, pair eta, dR) select the cell. This
// is the standard pre-unfolding approximation and is what the placeholder it replaces did; it
// becomes exact once detector-response unfolding lands, at which point the correct ordering is
// unfold first, then apply eps_reco in TRUTH bins (RDFBasedHistFilling/CorrectionStages.h).
//
// BINNING is NOT chosen here -- it is read off the histogram's own axes and CHECKED against
// ParamsSet::pair_pt_coarse_bins and
// CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap (.claude/CLAUDE.md Binnings). A
// mismatch THROWS: a pair corrected by another cell's efficiency is a silent failure -- every
// histogram still fills and every number still comes out.
//
// GUARDS, each one deliberate:
//   * pair pT / pair eta / dR are CLAMPED into the measured range. A pair above the top pair-pT
//     edge or beyond the top dR edge is corrected by the edge cell rather than dropped; the count
//     is reported, never silent.
//   * an EMPTY 3D cell does NOT return "no correction". Inside the signal region the three
//     variables are strongly correlated (dR ~< 2 m_uu / pT^pair), so ~60 % of the 3D cells carry
//     no measurement; leaving those pairs uncorrected would be a silent, ONE-SIDED bias. The
//     evaluator falls back, per pair and COUNTED:
//         3D cell -> the dR-INTEGRATED efficiency of the same (pair pT, pair eta) cell
//                 -> the INCLUSIVE efficiency of the signal region.
//     Only if even the inclusive number is missing does it return -1 ("no correction").
//   * the efficiency is FLOORED at kMinEff before inversion, so a statistically-starved cell
//     cannot blow a single pair's weight up without bound. Counted.
// =================================================================================================
struct PairRecoEffEvaluator {

    static constexpr double kMinEff = 0.05;   // same floor the Run-2 placeholder used

    std::unique_ptr<TH3D> h_eff;      // eps_reco(pair pT, pair eta, dR)
    std::unique_ptr<TH2D> h_eff2d;    // dR-integrated fallback
    double eff_incl = -1.0;           // inclusive fallback
    std::string wp, path;
    // ATOMIC: Eval() runs inside RDF lambdas under ImplicitMT. The lookups themselves are const
    // reads, so a race here would corrupt only the census -- which is precisely the thing that
    // tells us how many pairs took a fallback or were extrapolated past the measured domain.
    std::atomic<long long> n_eval{0}, n_clamp_pt{0}, n_clamp_eta{0}, n_clamp_dr{0};
    std::atomic<long long> n_3d{0}, n_2d{0}, n_incl{0}, n_none{0}, n_floor{0};

    // `use_tight` selects the working point; it MUST match the WP the crossx filters on.
    void Load(const std::string& file_path, bool use_tight)
    {
        wp   = use_tight ? "tight" : "medium";
        path = file_path;
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie())
            throw std::runtime_error("PairRecoEffEvaluator: cannot open " + path +
                                     " -- build it with "
                                     "plotting_codes/reco_effcy/build_pp24_fullsim_pair_reco_eff.C");
        const std::string key = "h_pair_reco_eff_" + wp;
        auto* h = dynamic_cast<TH3D*>(f->Get(key.c_str()));
        if (!h)
            throw std::runtime_error("PairRecoEffEvaluator: no " + key + " in " + path);
        h_eff.reset(static_cast<TH3D*>(h->Clone("h_pair_reco_eff_clone")));
        h_eff->SetDirectory(nullptr);

        auto* h2 = dynamic_cast<TH2D*>(f->Get((key + "_dr_integrated").c_str()));
        if (!h2)
            throw std::runtime_error("PairRecoEffEvaluator: no " + key + "_dr_integrated in "
                                     + path + " -- the dR-integrated fallback is REQUIRED, not "
                                     "optional (about 60 % of the 3D cells carry no measurement)");
        h_eff2d.reset(static_cast<TH2D*>(h2->Clone("h_pair_reco_eff2d_clone")));
        h_eff2d->SetDirectory(nullptr);

        auto* p_incl = dynamic_cast<TParameter<double>*>(
            f->Get(("pair_reco_eff_inclusive_" + wp).c_str()));
        if (!p_incl)
            throw std::runtime_error("PairRecoEffEvaluator: no pair_reco_eff_inclusive_" + wp
                                     + " in " + path);
        eff_incl = p_incl->GetVal();
        f->Close();

        CheckCanonicalBinning();

        int n_empty_cells = 0, n_cells = 0;
        for (int ix = 1; ix <= h_eff->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= h_eff->GetNbinsY(); ++iy)
                for (int iz = 1; iz <= h_eff->GetNbinsZ(); ++iz) {
                    ++n_cells;
                    if (h_eff->GetBinContent(ix, iy, iz) <= 0.) ++n_empty_cells;
                }
        std::cout << "PairRecoEffEvaluator [" << wp << " WP]: loaded eps_reco("
                  << h_eff->GetNbinsX() << " pair pT x " << h_eff->GetNbinsY() << " pair eta x "
                  << h_eff->GetNbinsZ() << " dR) = " << n_cells << " cells, " << n_empty_cells
                  << " empty (no correction there), from " << path << std::endl;
    }

    // Returns eps_reco, or -1 only if no level of the map carries a measurement at all
    // (caller -> w_reco = 1).
    double Eval(double pair_pt, double pair_eta, double dr)
    {
        ++n_eval;
        const TAxis* ax = h_eff->GetXaxis();
        const TAxis* ay = h_eff->GetYaxis();
        const TAxis* az = h_eff->GetZaxis();

        double x = pair_pt, y = pair_eta, z = dr;
        if (x < ax->GetXmin() || x >= ax->GetXmax()) { ++n_clamp_pt;  x = Clamp(x, ax); }
        if (y < ay->GetXmin() || y >= ay->GetXmax()) { ++n_clamp_eta; y = Clamp(y, ay); }
        if (z < az->GetXmin() || z >= az->GetXmax()) { ++n_clamp_dr;  z = Clamp(z, az); }

        const int ix = ax->FindBin(x), iy = ay->FindBin(y), iz = az->FindBin(z);
        double v = h_eff->GetBinContent(ix, iy, iz);
        if (v > 0.)      { ++n_3d; }
        else {
            v = h_eff2d->GetBinContent(ix, iy);
            if (v > 0.)  { ++n_2d; }
            else {
                v = eff_incl;
                if (v > 0.) { ++n_incl; }
                else        { ++n_none; return -1.0; }
            }
        }
        if (v < kMinEff) { ++n_floor; return kMinEff; }
        return (v > 1.0) ? 1.0 : v;
    }

    void PrintStats() const
    {
        const long long N = n_eval;
        auto pct = [&](long long n) { return N ? 100.0 * n / N : 0.0; };
        std::cout << "PairRecoEffEvaluator [" << wp << " WP]: " << N << " evaluations -- "
                  << n_clamp_pt  << " clamped in pair pT ("  << pct(n_clamp_pt)  << "%), "
                  << n_clamp_eta << " in pair eta ("         << pct(n_clamp_eta) << "%), "
                  << n_clamp_dr  << " in dR ("               << pct(n_clamp_dr)  << "%); "
                  << n_3d   << " from the 3D map ("      << pct(n_3d)   << "%), "
                  << n_2d   << " from the dR-INTEGRATED fallback (" << pct(n_2d)   << "%), "
                  << n_incl << " from the INCLUSIVE fallback ("     << pct(n_incl) << "%), "
                  << n_none << " with no efficiency at all ("       << pct(n_none)
                  << "%, no reco correction); " << n_floor << " floored at " << kMinEff
                  << " (" << pct(n_floor) << "%)" << std::endl;
    }

private:
    static double Clamp(double v, const TAxis* a)
    {
        const double lo = a->GetBinCenter(1), hi = a->GetBinCenter(a->GetNbins());
        return v < lo ? lo : (v > hi ? hi : v);
    }

    void CheckCanonicalBinning() const
    {
        static const ParamsSet pms{};
        static const CommonEffcyConfig cfg{};

        std::vector<double> eta;
        eta.push_back(cfg.pair_eta_proj_ranges_coarse_incl_gap.front().first);
        for (const auto& r : cfg.pair_eta_proj_ranges_coarse_incl_gap) eta.push_back(r.second);

        auto check = [this](const TAxis* ax, const std::vector<double>& edges, const char* what) {
            if (ax->GetNbins() != (int)edges.size() - 1)
                throw std::runtime_error("PairRecoEffEvaluator: " + std::string(what) + " axis of "
                    + path + " has " + std::to_string(ax->GetNbins()) + " bins but the canonical "
                    "binning has " + std::to_string(edges.size() - 1) + " -- stale efficiency file");
            for (size_t i = 0; i < edges.size(); ++i) {
                const double got = (i + 1 <= (size_t)ax->GetNbins())
                                 ? ax->GetBinLowEdge(i + 1) : ax->GetBinUpEdge(ax->GetNbins());
                if (std::fabs(got - edges[i]) > 1e-6)
                    throw std::runtime_error("PairRecoEffEvaluator: " + std::string(what)
                        + " edge " + std::to_string(i) + " is " + std::to_string(got) + " in "
                        + path + " but " + std::to_string(edges[i])
                        + " canonically -- stale efficiency file");
            }
        };
        check(h_eff->GetXaxis(), pms.pair_pt_coarse_bins, "pair-pT");
        check(h_eff->GetYaxis(), eta,                     "pair-eta");
        // The dR axis has no canonical vector of its own -- it is defined by the producer
        // (RDFBasedHistFillingPythia.h dr_bins_edges_for_pair_reco_eff) and simply read here.
    }
};

#endif  // PAIR_RECO_EFF_EVALUATOR_H
