// =================================================================================================
// FillMCTrigEffPairEff.cxx
//
// THE SINGLE-VALUE PAIR 2mu4 EFFICIENCY, mass-windowed
// (docs/tracking/mc_trigeff_single_value_pair_eff.md -- Physics Procedure §2, §PP-1, §PP-2, §PP-3)
//
// WHAT IS MEASURED. One number per cell = (pair pT coarse bin) x (|eta^pair| group) x (pair sign),
// on the pairs of one dimuon MASS WINDOW, from three sums over the SAME pairs:
//
//     S0 = sum_{all pairs}         w
//     S1 = sum_{pairs firing 2mu4} w
//     S2 = sum_{pairs firing 2mu4} w / [eps_MC(1) eps_MC(2)]
//
//     eps_2mu4^pair = S1 / S0     PURE        -- replaces the WHOLE per-pair trigger weight
//     K             = S2 / S0     CALIBRATED  -- MULTIPLIES the two single-muon efficiencies,
//                                                i.e. the single-number analogue of eps_dR
//
// Both are defined so the corresponding inverse-weighted estimator is UNBIASED when integrated
// over the cell: applying w/eps^pair to the firing pairs returns S0 exactly, and so does applying
// w/[eps(1) eps(2) K], because S2 is the MC estimate of sum_all w * P_pair/(eps1 eps2) and the
// unbiased single factor is the 1/(eps1 eps2)-weighted mean of P_pair/(eps1 eps2).
//
// WHY THIS EXISTS. In the top pair-pT cells the pp24 fullsim sample cannot constrain the SHAPE of
// eps_dR(dR): fits are rejected, one cell is corrected by the raw-bin placeholder at 2.157 where
// the physics requires <= 1 (mc_trigger_efficiency.md R26/R32, mc_trig_eff_closure.md R5). At
// fixed pair pT and pair eta the dimuon MASS essentially fixes the opening angle (dR ~ 2m/pT), so
// a mass window plus a (pT, |eta|) cell already pins down the pair kinematics the trigger responds
// to, and the single ratio measures in one number what the factorized product needs a fit for.
// The cost, stated rather than hidden: the correction is exact only INTEGRATED over its cell and
// inherits the MC pair spectrum inside it.
//
// WHY TWO MASS WINDOWS. The number is valid only for the mass mixture it was measured on. `sig`
// (1.08-2.9 GeV) is the signal window; `wide` (1-4 GeV) is the template-fit window containing the
// phi / J/psi / psi(2S) region. Their DIFFERENCE is the measurement of the trigger efficiency's
// mass dependence, and it decides whether one number may serve the wider fit.
//
// WHY THE SIGNS ARE SEPARATE. The dR correction is genuinely charge-dependent at small dR
// (same-sign/opposite-sign = 0.43 at dR = 0.025, an L1 close-by-RoI signature -- parent R25), and
// the background estimate is OS - SS.
//
// THE SAMPLE IS THE dR CORRECTION's OWN (MCTrigEffPairSel::Step3PairSelection), so the two
// procedures are measured on identical pairs and the closure comparison is a comparison of
// PROCEDURES and not of samples. Nothing is re-derived from raw NTUPs (.claude/CLAUDE.md
// §NTuple-Processing Provenance).
//
// Usage (from Analysis/RDFBasedHistFilling/):
//   root -l -b -q 'FillMCTrigEffPairEff.cxx+("pp_full", true)'    // Tight  (nominal)
//   root -l -b -q 'FillMCTrigEffPairEff.cxx+("pp_full", false)'   // Medium (WP systematic)
// =================================================================================================

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../Utilities/MCTrigEffPairSelection.h"
#include "../Utilities/PairTrigEffEvaluator.h"
#include "../Utilities/SingleMuEffEvaluator.h"
#include "../plotting_codes/trig_effcy/mc_based/dr_correction_ratio.h"
#include "../plotting_codes/trig_effcy/mc_based/dr_correction_sample_cfg.h"
#include "CommonEffcyConfig.h"

namespace {

// A file-level stamp, identical in spirit to FillMCTrigEffClosure's: a consumer holding only the
// output must be able to see WHICH eps_MC it was built from, mtime included -- a concurrent re-fit
// of the single-muon turn-ons is otherwise invisible.
std::string Stamp(const std::string& f)
{
    Long_t id = 0, sz = 0, fl = 0, mt = 0;
    if (gSystem->GetPathInfo(f.c_str(), &id, &sz, &fl, &mt) != 0)
        return std::string(Form("%s (MISSING)", f.c_str()));
    return std::string(Form("%s (mtime %ld, %ld B)", f.c_str(), mt, sz));
}

}  // namespace

// =================================================================================================
void FillMCTrigEffPairEff(const std::string& sample = "pp_full", bool use_tight_wp = true)
{
    // The `sig` window must still BE the signal region's mass window; a drift there would measure
    // the efficiency on a different sample from the one it is applied to.
    PairTrigEff::CheckSignalWindowMirror();

    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    if (cfg.key != "pp" && cfg.key != "pp_full")
        throw std::invalid_argument(
            "FillMCTrigEffPairEff: pp only. The Pb+Pb per-pair weight is the mu4 UNION "
            "eps_dR^single (eps1+eps2) - eps1 eps2 eps_dR^cross, a different formula, and the "
            "single-value replacement has not been defined for it (doc Scope).");

    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_col  = use_tight_wp ? "pass_tight" : "pass_medium";
    const std::string wp_text = use_tight_wp ? "Tight" : "Medium";

    // Same pair file the closure reads, built from the same directory + ntuple-processing naming.
    const std::string pair_file = cfg.mc_dir
        + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig"
        + (cfg.key == "pp_full" ? "_full" : "") + ".root";

    // ------------------------------------------------------------------ the cells (never retyped)
    const std::vector<double> pt_edges  = PairTrigEff::PairPtEdges();
    const DrAxisGroups        eta_grp   = PairTrigEff::AbsEtaGroups();
    const std::vector<double> eta_edges = eta_grp.edges;
    const int npt  = static_cast<int>(pt_edges.size())  - 1;
    const int neta = static_cast<int>(eta_edges.size()) - 1;
    const int first_delivered = PairTrigEff::FirstDeliveredPtBin();

    std::cout << "\n================ FillMCTrigEffPairEff: " << sample << " / " << wp_text
              << " muons ================\n"
              << "  pair file    : " << pair_file << "\n"
              << "  pair pT cells: " << npt << " (ParamsSet::pair_pt_coarse_bins)\n"
              << "  |eta^pair|   : " << DrGroupsDescribe(eta_grp, "|eta^pair|", "") << "\n"
              << "  delivered    : pair-pT bins " << first_delivered << ".." << npt << "  = ["
              << pt_edges[first_delivered - 1] << ", " << pt_edges[npt] << ") GeV  (bin "
              << first_delivered << " is the CONTROL region where the dR procedure also works)"
              << std::endl;

    // ------------------------------------------------------------------ eps_MC (for K only)
    // The SAME single-muon efficiency the Step-3 inverse weighting divides by, so K is the
    // single-number analogue of eps_dR and not a differently-normalized object. Heap-allocated:
    // the lambdas below are captured by LAZY RDF nodes and must outlive this scope.
    const std::string eps_mc_file = cfg.mc_dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root";
    auto* eps_mc = new SingleMuEffEvaluator();
    eps_mc->Load(SingleMuEffEvaluator::Src::kMCDirect, eps_mc_file);

    // ------------------------------------------------------------------ book, per sign x window
    std::map<std::string, ROOT::RDF::RResultPtr<TH2D>> books;
    std::map<std::string, ROOT::RDF::RResultPtr<ULong64_t>> counts;
    auto model = [&](const std::string& n, const char* zt) {
        return ROOT::RDF::TH2DModel(n.c_str(),
                                    Form(";p_{T}^{pair} [GeV];|#eta^{pair}|;%s", zt),
                                    npt, pt_edges.data(), neta, eta_edges.data());
    };

    // Built ONCE, before the sign loop: they are the same for both signs, and assigning them
    // inside the loop made the provenance stamp record whichever sign happened to run last.
    const std::string base_sel_text = MCTrigEffPairSel::Step3PairSelection(true);
    const std::string cell_sel_text =
        Form("pair_pt >= %.10g && pair_pt < %.10g && abs_pair_eta < %.10g",
             pt_edges.front(), pt_edges.back(), eta_edges.back());

    for (const auto& S : PairTrigEff::Signs()) {
        ROOT::RDataFrame df(S.tree, pair_file);
        ROOT::RDF::RNode d = df;
        d = d.Alias("m1_pt", "m1.pt").Alias("m1_eta", "m1.eta").Alias("m1_charge", "m1.charge")
             .Alias("m1_wp", "m1." + wp_col)
             .Alias("m1_truth_pt", "m1.truth_pt").Alias("m1_truth_eta", "m1.truth_eta")
             .Alias("m2_pt", "m2.pt").Alias("m2_eta", "m2.eta").Alias("m2_charge", "m2.charge")
             .Alias("m2_wp", "m2." + wp_col)
             .Alias("m2_truth_pt", "m2.truth_pt").Alias("m2_truth_eta", "m2.truth_eta");

        // EXACTLY the sample the dR correction is measured on -- that is what makes the closure a
        // comparison of two PROCEDURES rather than of two samples.
        d = d.Filter(base_sel_text, "MC trigger-efficiency pair selection");

        // Inside the cell grid. The |eta^pair| axis is the sign-independent FOLD, so the cut is on
        // |pair_eta|; applying the folded bounds to the SIGNED branch would silently drop every
        // negative-eta pair (the bug mc_trigeff_dr_binning_approaches.md D11 records).
        d = d.Define("abs_pair_eta", "fabs(pair_eta)");
        d = d.Filter(cell_sel_text, "pair inside the cell grid");

        // The trigger decision and the two single-muon efficiencies, evaluated ONCE on the firing
        // node before the mass-window split (every window is a subset of the same node).
        ROOT::RDF::RNode dn = d.Filter("pass2mu4", "pair fires 2mu4");
        dn = dn.Define("epsmc1", [ev = eps_mc](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                       {"m1_pt", "m1_eta", "m1_charge"})
               .Define("epsmc2", [ev = eps_mc](float pt, float eta, int q) { return ev->Eval(pt, eta, q); },
                       {"m2_pt", "m2_eta", "m2_charge"})
               .Define("epsprod", "epsmc1 * epsmc2")
               // The two numerator weights, and the A/B sums their CONDITIONAL (binomial-correct)
               // errors need. Numerator and denominator are not independent -- the numerator is a
               // re-weighted SUBSET of the denominator -- and TH1::Divide's independent-error form
               // over-states these bars by 1.5-3.2x (parent R12). Convention exactly as in
               // FillMCTrigEffClosure: Var = A - R*B.
               //   eps = S1/S0 : a = w,            A = sum a^2,  B = sum a^2       (p ~ R in-cell)
               //   K   = S2/S0 : a = w/(e1 e2),    A = sum a^2,  B = sum a^2 e1 e2 (p = e1 e2 K)
               .Define("w_k",   "weight / epsprod")
               .Define("A_eps", "weight * weight")
               .Define("B_eps", "weight * weight")
               .Define("A_k",   "w_k * w_k")
               .Define("B_k",   "w_k * w_k * epsprod");

        for (const auto& W : PairTrigEff::Windows()) {
            const std::string mass = Form("minv > %.10g && minv < %.10g", W.lo, W.hi);
            auto d_all  = d .Filter(mass, "mass window " + W.token);
            auto d_pass = dn.Filter(mass, "mass window " + W.token + " (firing)");
            const std::string sw = S.token + "_" + W.token;
            auto N = [&](const char* q) { return PairTrigEff::HistName(q, S.token, W.token); };

            books.emplace(N("den"),  d_all .Histo2D(model(N("den"),  "#Sigma w"),
                                                    "pair_pt", "abs_pair_eta", "weight"));
            books.emplace(N("num"),  d_pass.Histo2D(model(N("num"),  "#Sigma w"),
                                                    "pair_pt", "abs_pair_eta", "weight"));
            books.emplace(N("numk"), d_pass.Histo2D(model(N("numk"), "#Sigma w/#varepsilon_{1}#varepsilon_{2}"),
                                                    "pair_pt", "abs_pair_eta", "w_k"));
            books.emplace(N("Aeps"), d_pass.Histo2D(model(N("Aeps"), "#Sigma a^{2}"),
                                                    "pair_pt", "abs_pair_eta", "A_eps"));
            books.emplace(N("Beps"), d_pass.Histo2D(model(N("Beps"), "#Sigma a^{2}p"),
                                                    "pair_pt", "abs_pair_eta", "B_eps"));
            books.emplace(N("Ak"),   d_pass.Histo2D(model(N("Ak"),   "#Sigma a^{2}"),
                                                    "pair_pt", "abs_pair_eta", "A_k"));
            books.emplace(N("Bk"),   d_pass.Histo2D(model(N("Bk"),   "#Sigma a^{2}p"),
                                                    "pair_pt", "abs_pair_eta", "B_k"));
            // RAW (UNWEIGHTED) counts -- the statistical reach the weighted numbers hide. Never
            // omitted: a cell can carry a smooth-looking efficiency and a handful of pairs.
            books.emplace(N("nraw"),     d_all .Histo2D(model(N("nraw"),     "pairs"),
                                                        "pair_pt", "abs_pair_eta"));
            books.emplace(N("nrawpass"), d_pass.Histo2D(model(N("nrawpass"), "pairs"),
                                                        "pair_pt", "abs_pair_eta"));
            counts.emplace(sw + "_all",  d_all .Count());
            counts.emplace(sw + "_pass", d_pass.Count());
        }
    }

    // ------------------------------------------------------------------ the cells, both modes
    // `nomerge` is the canonical 8-cell pair-pT axis. `ptmerge` combines the last two cells into
    // one [72.08, 150) GeV cell -- the pair-pT analogue of the dR correction's own `nocorr_ptmerge`
    // (mc_trigger_efficiency.md R32), and the same trade: it buys statistics where the sample runs
    // out, at the cost of describing a MIXTURE of two cells. It is NOT a new binning: the merged
    // cell is the two source cells' num / den / A / B SUMMED BEFORE the ratio, which is exactly
    // what filling a 7-bin axis would have given (.claude/CLAUDE.md Binnings item 4), and it is
    // opt-in and suffixed so it cannot be mistaken for the un-merged measurement.
    std::vector<std::unique_ptr<TH2D>> owned;
    std::map<std::string, TH2D*> H;          // every delivered histogram, keyed by its final name
    // A COLLIDING NAME MUST THROW. std::map::emplace is a no-op on a duplicate key, so a collision
    // would silently keep the stale histogram, drop the new one from the file, and make every later
    // H.at(name) return the wrong object -- with no diagnostic at all.
    auto add = [&H](const std::string& name, TH2D* h) {
        if (!H.emplace(name, h).second)
            throw std::runtime_error("FillMCTrigEffPairEff: duplicate histogram name '" + name
                                     + "' -- two (quantity, sign, window, cell mode) combinations "
                                       "map to one name");
    };
    for (auto& kv : books) add(kv.first, kv.second.GetPtr());

    const std::vector<double> mpt = PairTrigEff::PairPtEdges("ptmerge");
    const int nmpt = (int)mpt.size() - 1;
    auto merge_last_two_pt = [&](const TH2D* h, const std::string& name) {
        auto m = std::unique_ptr<TH2D>(new TH2D(name.c_str(), h->GetTitle(),
                                                nmpt, mpt.data(), neta, eta_edges.data()));
        m->SetDirectory(nullptr);
        m->Sumw2();
        for (int iy = 1; iy <= neta; ++iy) {
            for (int ix = 1; ix <= nmpt - 1; ++ix) {           // the untouched low cells
                m->SetBinContent(ix, iy, h->GetBinContent(ix, iy));
                m->SetBinError  (ix, iy, h->GetBinError  (ix, iy));
            }
            // the merged top cell: contents ADD, errors add in quadrature (these are all SUMS --
            // of weights, of squared weights -- never ratios, so summing them is exact.)
            const double c = h->GetBinContent(npt - 1, iy) + h->GetBinContent(npt, iy);
            const double e1 = h->GetBinError(npt - 1, iy), e2 = h->GetBinError(npt, iy);
            m->SetBinContent(nmpt, iy, c);
            m->SetBinError  (nmpt, iy, std::sqrt(e1 * e1 + e2 * e2));
        }
        TH2D* raw = m.get();
        owned.push_back(std::move(m));
        return raw;
    };

    // ------------------------------------------------------------------ form the two ratios
    // Built here rather than by every consumer, so the conditional error is applied in ONE place
    // and a reader of the file gets the number AND its bar without re-deriving either. The SAME
    // routine serves both cell modes -- the merged inputs are ordinary sums, so nothing about the
    // ratio or its error changes.
    auto form_ratio = [&](const std::string& q, const std::string& sign, const std::string& win,
                          const std::string& mode, const char* num_key, const char* A_key,
                          const char* B_key) {
        const std::string name = PairTrigEff::HistName(q, sign, win, mode);
        TH2D* den = H.at(PairTrigEff::HistName("den",  sign, win, mode));
        TH2D* num = H.at(PairTrigEff::HistName(num_key, sign, win, mode));
        TH2D* A   = H.at(PairTrigEff::HistName(A_key,   sign, win, mode));
        TH2D* B   = H.at(PairTrigEff::HistName(B_key,   sign, win, mode));
        auto r = std::unique_ptr<TH2D>(static_cast<TH2D*>(num->Clone(name.c_str())));
        r->SetDirectory(nullptr);
        r->Divide(den);
        // SetConditionalRatioErrors is written for TH1D (one row at a time); the cells here are
        // independent of one another, so it is applied row by row on projections and the errors
        // copied back. Same single implementation, no second copy of the formula.
        for (int iy = 1; iy <= neta; ++iy) {
            const std::string tag = "_" + name + "_row" + std::to_string(iy);
            std::unique_ptr<TH1D> rr(r  ->ProjectionX(("r" + tag).c_str(), iy, iy, "e"));
            std::unique_ptr<TH1D> dd(den->ProjectionX(("d" + tag).c_str(), iy, iy, "e"));
            std::unique_ptr<TH1D> aa(A  ->ProjectionX(("a" + tag).c_str(), iy, iy, "e"));
            std::unique_ptr<TH1D> bb(B  ->ProjectionX(("b" + tag).c_str(), iy, iy, "e"));
            for (auto* h : {rr.get(), dd.get(), aa.get(), bb.get()}) h->SetDirectory(nullptr);
            SetConditionalRatioErrors(rr.get(), dd.get(), aa.get(), bb.get());
            for (int ix = 1; ix <= r->GetNbinsX(); ++ix)
                r->SetBinError(ix, iy, rr->GetBinError(ix));
        }
        r->SetTitle(Form(";p_{T}^{pair} [GeV];|#eta^{pair}|;%s",
                         q == "eps" ? "#varepsilon_{2#mu4}^{pair}" : "K"));
        TH2D* raw = r.get();
        owned.push_back(std::move(r));
        add(name, raw);
    };

    for (const auto& S : PairTrigEff::Signs())
        for (const auto& W : PairTrigEff::Windows()) {
            // the merged SUMS first -- the ratios of both modes are then formed the same way
            for (const char* q : {"den", "num", "numk", "Aeps", "Beps", "Ak", "Bk",
                                  "nraw", "nrawpass"}) {
                const std::string src = PairTrigEff::HistName(q, S.token, W.token, "nomerge");
                const std::string dst = PairTrigEff::HistName(q, S.token, W.token, "ptmerge");
                add(dst, merge_last_two_pt(H.at(src), dst));
            }
            for (const auto& M : PairTrigEff::CellModes()) {
                form_ratio("eps", S.token, W.token, M.token, "num",  "Aeps", "Beps");
                form_ratio("k",   S.token, W.token, M.token, "numk", "Ak",   "Bk");
            }
        }


    // ------------------------------------------------------------------ write
    const std::string out_name = PairTrigEff::FileName(cfg.mc_dir, cfg.mc_label, wp_suf);
    TFile fout(out_name.c_str(), "RECREATE");
    if (fout.IsZombie()) throw std::runtime_error("FillMCTrigEffPairEff: cannot open " + out_name);
    for (auto& kv : H) kv.second->Write(kv.first.c_str());

    std::string wins;
    for (const auto& W : PairTrigEff::Windows())
        wins += Form("%s=[%g,%g] ", W.token.c_str(), W.lo, W.hi);
    TNamed("provenance",
           Form("SINGLE-VALUE pair 2mu4 efficiency "
                "(docs/tracking/mc_trigeff_single_value_pair_eff.md) | sample=%s label=%s WP=%s"
                " | eps^pair = sum_pass w / sum_all w  (PURE: replaces the whole per-pair trigger"
                " weight); K = sum_pass w/(eps_MC1 eps_MC2) / sum_all w  (CALIBRATED: multiplies"
                " the two single-muon efficiencies, the single-number analogue of eps_dR)"
                " | eps_MC: %s | cells: %d pair pT (ParamsSet::pair_pt_coarse_bins) x %d |eta^pair|"
                " groups (%s) x signs os,ss | mass windows: %s| CELL MODES: `` = the canonical %d"
                " pair-pT cells, `_ptmerge` = the last two combined into [%g, %g) GeV (a projection"
                " of the same filled sums, NOT a new binning) | DELIVERED from %g GeV up"
                " | pair file: %s | base selection: %s | cell selection: %s"
                " | delivery gate: >= %d raw pairs and a value in (0, %g]"
                " | NOT wired into the cross-section",
                sample.c_str(), cfg.mc_label.c_str(), wp_text.c_str(), Stamp(eps_mc_file).c_str(),
                npt, neta, DrGroupsDescribe(eta_grp, "|eta^pair|", "").c_str(), wins.c_str(),
                npt, mpt[nmpt - 1], mpt[nmpt], PairTrigEff::FirstDeliveredPtEdge(),
                Stamp(pair_file).c_str(), base_sel_text.c_str(), cell_sel_text.c_str(),
                PairTrigEff::MinCellPairs(), PairTrigEff::MaxDeliveredValue()))
        .Write();
    fout.Close();
    std::cout << "\nFillMCTrigEffPairEff: wrote " << H.size() << " TH2D to " << out_name
              << std::endl;

    // ------------------------------------------------------------------ the tables
    eps_mc->PrintStats(cfg.mc_label);
    for (const auto& S : PairTrigEff::Signs())
        for (const auto& W : PairTrigEff::Windows())
            std::cout << "  " << S.token << " / " << W.token << ": selected pairs "
                      << *counts.at(S.token + "_" + W.token + "_all") << ", firing 2mu4 "
                      << *counts.at(S.token + "_" + W.token + "_pass") << std::endl;

    for (const auto& S : PairTrigEff::Signs()) {
        std::cout << "\n===== single-value pair 2mu4 efficiency, " << S.text << ", " << wp_text
                  << " =====\n"
                  << "  (delivered cells are pair-pT bins " << first_delivered << ".." << npt
                  << "; the lower bins are the sanity check against the dR procedure)\n  "
                  << std::setw(22) << std::left << "pair pT [GeV]"
                  << std::setw(16) << "|eta^pair|";
        for (const auto& W : PairTrigEff::Windows())
            std::cout << std::setw(34) << ("eps^pair " + W.token)
                      << std::setw(22) << ("K " + W.token);
        std::cout << "raw pairs (all / 2mu4), sig" << std::endl;

        for (int ix = 1; ix <= npt; ++ix)
            for (int iy = 1; iy <= neta; ++iy) {
                std::cout << (ix >= first_delivered ? "  " : "  (") << std::setw(ix >= first_delivered ? 22 : 21)
                          << std::left << Form("[%.2f, %.2f)", pt_edges[ix - 1], pt_edges[ix])
                          << std::setw(16) << Form("[%.1f, %.1f)", eta_edges[iy - 1], eta_edges[iy]);
                for (const auto& W : PairTrigEff::Windows()) {
                    const TH2D* e = H.at(PairTrigEff::HistName("eps", S.token, W.token));
                    const TH2D* k = H.at(PairTrigEff::HistName("k",   S.token, W.token));
                    std::cout << std::setw(34) << Form("%.4f +- %.4f", e->GetBinContent(ix, iy),
                                                       e->GetBinError(ix, iy))
                              << std::setw(22) << Form("%.4f +- %.4f", k->GetBinContent(ix, iy),
                                                       k->GetBinError(ix, iy));
                }
                const TH2D* na = H.at(PairTrigEff::HistName("nraw",     S.token, "sig"));
                const TH2D* np = H.at(PairTrigEff::HistName("nrawpass", S.token, "sig"));
                std::cout << Form("%.0f / %.0f", na->GetBinContent(ix, iy), np->GetBinContent(ix, iy))
                          << (ix >= first_delivered ? "" : ")") << std::endl;
            }

        // THE MASS-DEPENDENCE DIAGNOSTIC. `wide` is a yield-weighted mixture of `sig` and the rest
        // of 1-4 GeV, so the efficiency of that REST is exact arithmetic on the two windows --
        // no third fill, and no new mass binning invented (.claude/CLAUDE.md §Binnings).
        std::cout << "\n  mass dependence, " << S.text
                  << ": eps^pair inside the signal window vs OUTSIDE it within 1-4 GeV\n  "
                  << std::setw(22) << std::left << "pair pT [GeV]" << std::setw(16) << "|eta^pair|"
                  << std::setw(14) << "eps(sig)" << std::setw(14) << "eps(outside)"
                  << std::setw(14) << "eps(wide)" << "raw pairs outside" << std::endl;
        for (int ix = first_delivered; ix <= npt; ++ix)
            for (int iy = 1; iy <= neta; ++iy) {
                const TH2D* ds = H.at(PairTrigEff::HistName("den", S.token, "sig"));
                const TH2D* dw = H.at(PairTrigEff::HistName("den", S.token, "wide"));
                const TH2D* ns = H.at(PairTrigEff::HistName("num", S.token, "sig"));
                const TH2D* nw = H.at(PairTrigEff::HistName("num", S.token, "wide"));
                const TH2D* ra = H.at(PairTrigEff::HistName("nraw", S.token, "sig"));
                const TH2D* rw = H.at(PairTrigEff::HistName("nraw", S.token, "wide"));
                const double dD = dw->GetBinContent(ix, iy) - ds->GetBinContent(ix, iy);
                const double dN = nw->GetBinContent(ix, iy) - ns->GetBinContent(ix, iy);
                std::cout << "  " << std::setw(22) << std::left
                          << Form("[%.2f, %.2f)", pt_edges[ix - 1], pt_edges[ix])
                          << std::setw(16) << Form("[%.1f, %.1f)", eta_edges[iy - 1], eta_edges[iy])
                          << std::setw(14)
                          << Form("%.4f", H.at(PairTrigEff::HistName("eps", S.token, "sig"))
                                              ->GetBinContent(ix, iy))
                          << std::setw(14) << (dD > 0 ? Form("%.4f", dN / dD) : "-")
                          << std::setw(14)
                          << Form("%.4f", H.at(PairTrigEff::HistName("eps", S.token, "wide"))
                                              ->GetBinContent(ix, iy))
                          << Form("%.0f", rw->GetBinContent(ix, iy) - ra->GetBinContent(ix, iy))
                          << std::endl;
            }
    }
    // ------------------------------------------------------------------ the pT-merged variant
    // Only the SIGNAL window is tabulated (the request's scope); the merged `wide` cells are
    // written to the file all the same, because the merge is a projection of sums that are already
    // there and costs nothing to produce.
    for (const auto& S : PairTrigEff::Signs()) {
        std::cout << "\n===== single-value pair 2mu4 efficiency, LAST TWO pair-pT CELLS MERGED, "
                  << S.text << ", " << wp_text << ", signal window =====\n  "
                  << std::setw(22) << std::left << "pair pT [GeV]" << std::setw(16) << "|eta^pair|"
                  << std::setw(26) << "eps^pair sig" << std::setw(26) << "K sig"
                  << "raw pairs (all / 2mu4)" << std::endl;
        const TH2D* e  = H.at(PairTrigEff::HistName("eps",      S.token, "sig", "ptmerge"));
        const TH2D* k  = H.at(PairTrigEff::HistName("k",        S.token, "sig", "ptmerge"));
        const TH2D* na = H.at(PairTrigEff::HistName("nraw",     S.token, "sig", "ptmerge"));
        const TH2D* np = H.at(PairTrigEff::HistName("nrawpass", S.token, "sig", "ptmerge"));
        for (int ix = PairTrigEff::FirstDeliveredPtBin("ptmerge"); ix <= nmpt; ++ix)
            for (int iy = 1; iy <= neta; ++iy)
                std::cout << "  " << std::setw(22) << std::left
                          << Form("[%.2f, %.2f)", mpt[ix - 1], mpt[ix])
                          << std::setw(16) << Form("[%.1f, %.1f)", eta_edges[iy - 1], eta_edges[iy])
                          << std::setw(26) << Form("%.4f +- %.4f", e->GetBinContent(ix, iy),
                                                   e->GetBinError(ix, iy))
                          << std::setw(26) << Form("%.4f +- %.4f", k->GetBinContent(ix, iy),
                                                   k->GetBinError(ix, iy))
                          << Form("%.0f / %.0f", na->GetBinContent(ix, iy),
                                  np->GetBinContent(ix, iy)) << std::endl;
    }

    // ------------------------------------------------------------------ the delivery gate
    // Which delivered cells PairTrigEffEvaluator will refuse, and why. Nothing is removed from the
    // file -- every cell above is measured, written and printed; this states which of them a
    // consumer must not apply, so a 4-pair fluctuation cannot reach the cross-section through the
    // same API as a 4830-pair measurement (PairTrigEff::MinCellPairs / MaxDeliveredValue).
    std::cout << "\n===== delivery gate: >= " << PairTrigEff::MinCellPairs()
              << " raw pairs and a value in (0, " << PairTrigEff::MaxDeliveredValue()
              << "] =====" << std::endl;
    int n_gated = 0;
    for (const auto& S : PairTrigEff::Signs())
      for (const auto& W : PairTrigEff::Windows())
        for (const auto& M : PairTrigEff::CellModes()) {
            const std::vector<double>& mp = (M.token == "ptmerge") ? mpt : pt_edges;
            const int nb = (int)mp.size() - 1;
            const int lo = PairTrigEff::FirstDeliveredPtBin(M.token);
            const TH2D* nr = H.at(PairTrigEff::HistName("nraw", S.token, W.token, M.token));
            for (const auto& q : {std::string("eps"), std::string("k")}) {
                const TH2D* v = H.at(PairTrigEff::HistName(q, S.token, W.token, M.token));
                for (int ix = lo; ix <= nb; ++ix)
                    for (int iy = 1; iy <= neta; ++iy) {
                        const double n = nr->GetBinContent(ix, iy);
                        const double x = v->GetBinContent(ix, iy);
                        const char* why = (n <= 0) ? "no pairs"
                                        : (n < PairTrigEff::MinCellPairs()) ? "too few raw pairs"
                                        : (x <= 0) ? "measured zero (pairs exist, none fired)"
                                        : (x > PairTrigEff::MaxDeliveredValue())
                                              ? "value above the physical maximum" : nullptr;
                        if (!why) continue;
                        ++n_gated;
                        std::cout << "  REFUSED  " << S.token << " / " << W.token << " / "
                                  << M.token << " / " << q
                                  << Form("  pT [%.2f, %.2f)  |eta| [%.1f, %.1f)  value %.4f "
                                          "+- %.4f  raw pairs %.0f  -- %s",
                                          mp[ix - 1], mp[ix], eta_edges[iy - 1],
                                          eta_edges[iy], x, v->GetBinError(ix, iy), n, why)
                                  << std::endl;
                    }
            }
        }
    if (n_gated == 0) std::cout << "  every delivered cell passes." << std::endl;

    std::cout << "done." << std::endl;
}
