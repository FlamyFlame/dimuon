#ifndef PAIR_TRIG_EFF_EVALUATOR_H
#define PAIR_TRIG_EFF_EVALUATOR_H

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TFile.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TString.h>

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../MuonObjectsParamsAndHelpers/FullSimSampleType.h"   // FullSimMCTrigEffPairEffDir
#include "../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../plotting_codes/trig_effcy/mc_based/dr_correction_cell_groups.h"
#include "MCTrigEffPairSelection.h"

// =================================================================================================
// THE SINGLE-VALUE PAIR 2mu4 EFFICIENCY -- cell definition, on-disk naming, and the reader
// (docs/tracking/mc_trigeff_single_value_pair_eff.md, Physics Procedure §2 / §PP-1 / §PP-2).
//
// WHAT THIS OBJECT IS. One directly measured number per
//     (pair pT coarse bin) x (|eta^pair| group) x (pair sign),
// on the pairs of a dimuon MASS WINDOW:
//
//     eps_2mu4^pair(cell) = sum_{pairs firing 2mu4} w / sum_{all pairs} w          (PURE)
//     K(cell)             = sum_{pairs firing 2mu4} w/(eps_MC(1) eps_MC(2))
//                           / sum_{all pairs} w                                    (CALIBRATED)
//
// It is the ALTERNATIVE to the factorized eps(pT1,q*eta1) eps(pT2,q*eta2) eps_dR(dR) weight in the
// top pair-pT cells, where the dR-shape fit is statistics-starved (mc_trigger_efficiency.md
// R26/R32, mc_trig_eff_closure.md R5).
//
// ⚠ IT IS NOT eps_dR AND MUST NEVER BE ROUTED THROUGH DrCorrectionEvaluator. It carries NO dR
// dependence, it replaces the WHOLE per-pair weight (the pure form) rather than multiplying it,
// it is mass-windowed, and it exists only in the cells it was measured in. Every consumer of
// DrCorrectionEvaluator assumes the form f(dR)/C; this object has no such form.
//
// WHICH OF THE TWO TO APPLY. `eps` is the user's stated procedure and is what the closure figure
// draws; `K` is the one the analysis's calibration chain wants, because the cross-section takes
// its single-muon efficiencies from the DATA tag-and-probe and only the correlation from MC. The
// MC/data pair over-efficiency is measured: <r_1 r_2> = 1.27 (mc_trig_eff_closure.md R1), so
// applying the pure MC eps on data would bias those bins by ~25 % -- invisibly to an MC closure,
// which is self-contained by construction. K is measured with eps_MC in MC exactly as eps_dR is,
// and is applied on top of the data eps^nc exactly as eps_dR is. NEITHER is wired into the
// cross-section by the doc that created this header.
//
// BINNINGS ARE READ, NEVER RETYPED (.claude/CLAUDE.md §Binnings): the pair-pT axis is
// ParamsSet::pair_pt_coarse_bins and the |eta^pair| axis is the sign-independent 3-group fold
// MakeDrEtaGroups(..., true) builds from CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap
// -- the SAME fold the `nocorr_etamerge*` dR modes use, so the two groupings cannot drift.
// CheckCanonicalBinning() below re-derives both from those sources and THROWS on any mismatch, so
// a stale file cannot be read as if it described today's cells.
// =================================================================================================
namespace PairTrigEff {

// ---------------------------------------------------------------------------- the mass windows
struct MassWindow {
    std::string token;   // file/histogram token
    double      lo, hi;  // GeV
    std::string text;    // for a legend / a printout
};

// `sig`  = the single-b signal window. Combined with MCTrigEffPairSel::Step3PairSelection it IS
//          the signal region: SingleBSignalCutsReco adds only the pair-pT threshold
//          `ParamsSet::signal_pair_pt_min` (implied by the coarse pair-pT axis, which starts at
//          the same value) and the fiducial gap cut (already in the base selection).
//          The value is deliberately NOT written here: it moved 8 -> 9 GeV on 2026-09-08
//          (docs/tracking/mu_pt45_gap125_pairpt9_adoption.md D3) and a retyped copy in a comment
//          is how the two silently drift -- note CheckSignalWindowMirror below compares only the
//          `minv` half of that string, so a pair-pT mismatch would NOT be caught.
// `wide` = the template-fit window, which contains the phi, J/psi and psi(2S) mass region. The MC
//          pair file carries NO resonance veto, so this is the honest 1-4 GeV mixture the MC
//          produces (docs/tracking/mc_trigeff_single_value_pair_eff.md §PP-2).
//
// ⚠ MIRROR of the signal window in RDFBasedHistFilling/RDFBasedHistFillingPP.cxx `signal_cuts`
// and Utilities/MCTrigEffPairSelection.h SingleBSignalCutsReco(). There is no ParamsSet member to
// read it from today, so CheckSignalWindowMirror() below asserts at RUN TIME that the values here
// still appear verbatim in SingleBSignalCutsReco() -- a drift is a throw, not a silent shift.
inline const std::vector<MassWindow>& Windows()
{
    static const std::vector<MassWindow> w = {
        {"sig",  1.08, 2.9, "1.08 < m_{#mu#mu} < 2.9 GeV (signal window)"},
        {"wide", 1.00, 4.0, "1 < m_{#mu#mu} < 4 GeV (template-fit window)"},
    };
    return w;
}

inline const MassWindow& Window(const std::string& token)
{
    for (const auto& w : Windows()) if (w.token == token) return w;
    throw std::invalid_argument("PairTrigEff::Window: unknown mass-window token '" + token + "'");
}

inline void CheckSignalWindowMirror()
{
    const MassWindow& s = Window("sig");
    const std::string sel = MCTrigEffPairSel::SingleBSignalCutsReco();
    const std::string need = "minv > " + std::string(Form("%g", s.lo))
                           + " && minv < " + std::string(Form("%g", s.hi));
    if (sel.find(need) == std::string::npos)
        throw std::runtime_error(
            "PairTrigEff: the `sig` mass window (" + need + ") is no longer the signal region's "
            "mass window -- MCTrigEffPairSel::SingleBSignalCutsReco() now reads '" + sel +
            "'. The single-value efficiency would be measured on a different sample from the one "
            "it is applied to. Update PairTrigEff::Windows() deliberately, do not widen this "
            "check.");
}

// ---------------------------------------------------------------------------- the pair signs
// Repo convention (memory reference_sign_convention): muon_pair_tree_sign1 = SAME sign,
// muon_pair_tree_sign2 = OPPOSITE sign. The dR-correction chain's tokens are "ss" / "os"
// (mc_trigger_efficiency.md DrCorrSignText) and are reused here unchanged.
struct PairSign { std::string token, tree, text; };
inline const std::vector<PairSign>& Signs()
{
    static const std::vector<PairSign> s = {
        {"os", "muon_pair_tree_sign2", "opposite sign"},
        {"ss", "muon_pair_tree_sign1", "same sign"},
    };
    return s;
}
inline const PairSign& Sign(const std::string& token)
{
    for (const auto& s : Signs()) if (s.token == token) return s;
    throw std::invalid_argument("PairTrigEff::Sign: unknown pair-sign token '" + token + "'");
}

// ---------------------------------------------------------------------------- the cell modes
// A second way of buying statistics in the top cells, on the pair-pT axis this time: the
// alternative to a finer correction is a COARSER one. `ptmerge` combines the last two coarse
// pair-pT cells into a single [72.08, 150) GeV cell, exactly as the dR correction's own
// `nocorr_ptmerge` mode does (mc_trigger_efficiency.md R32).
//
// ⚠ IT IS NOT A NEW BINNING (.claude/CLAUDE.md §Binnings item 4). The FILLED histograms keep the
// canonical 8 pair-pT bins; a merged cell is the two source bins' num / den / A / B SUMMED BEFORE
// the ratio, which is numerically identical to having filled a 7-bin axis. The variant is opt-in
// and SUFFIXED, so it can neither overwrite nor be mistaken for the un-merged one.
struct CellMode {
    std::string token;    // "nomerge" / "ptmerge"
    std::string suffix;   // "" / "_ptmerge" -- empty for the default so existing names are unchanged
    std::string text;
};
inline const std::vector<CellMode>& CellModes()
{
    static const std::vector<CellMode> m = {
        {"nomerge", "",          "the canonical 8 coarse pair-pT cells"},
        {"ptmerge", "_ptmerge",  "the last two coarse pair-pT cells combined into one"},
    };
    return m;
}
inline const CellMode& Mode(const std::string& token)
{
    for (const auto& m : CellModes()) if (m.token == token) return m;
    throw std::invalid_argument("PairTrigEff::Mode: unknown cell-mode token '" + token + "'");
}

// ---------------------------------------------------------------------------- on-disk naming
// One TH2D per (quantity, sign, window). X = the canonical coarse pair-pT axis, Y = the 3
// |eta^pair| groups. Quantities:
//   eps   the pure single-value efficiency S1/S0        (bin error = conditional/binomial)
//   k     the calibrated correction       S2/S0         (bin error = conditional/binomial)
//   den   S0 = sum_all w                                (Sumw2)
//   num   S1 = sum_pass w
//   numk  S2 = sum_pass w/(eps_MC1 eps_MC2)
//   nraw      raw (UNWEIGHTED) all-pairs count   -- the statistical reach, never hidden
//   nrawpass  raw (UNWEIGHTED) firing-pair count
inline std::string HistName(const std::string& quantity, const std::string& sign,
                            const std::string& window, const std::string& mode = "nomerge")
{
    return "h_paireff_" + quantity + "_" + sign + "_" + window + Mode(mode).suffix;
}

// The deliverable file. One file per sample x working point; the sign and the window are inside.
// `sample_dir` is the sample ROOT directory (DrCorrSample::sample_dir); the file lives in its
// mc_trig_eff/pair_eff/ subtree (FullSimSampleType.h "PER-SAMPLE DIRECTORY LAYOUT").
inline std::string FileName(const std::string& sample_dir, const std::string& mc_label,
                            const std::string& wp_suffix)
{
    return FullSimMCTrigEffPairEffDir(sample_dir) + "pair_trig_eff_" + mc_label + wp_suffix + ".root";
}

// ---------------------------------------------------------------------------- the cell axes
// Both are DERIVED, never typed. The |eta| fold is built by the caller (the fill macro and the
// evaluator both hand in the edges they got from MakeDrEtaGroups) -- this header only states what
// the pair-pT axis is, because that one has a single canonical source with no grouping step.
inline std::vector<double> PairPtEdges(const std::string& mode = "nomerge")
{
    static const ParamsSet pms;
    std::vector<double> e = pms.pair_pt_coarse_bins;
    if (Mode(mode).token == "ptmerge") {
        if (e.size() < 3)
            throw std::runtime_error("PairTrigEff::PairPtEdges: cannot merge the last two pair-pT "
                                     "cells of a binning with fewer than 2");
        e.erase(e.end() - 2);         // drop the edge BETWEEN the last two cells
    }
    return e;
}

// The pair-pT edge the DELIVERED range starts at: the lower edge of the request's control cell.
// The VALUE is deliberately not written here -- it moves with ParamsSet::pair_pt_coarse_bins (it
// was 49.97 GeV on the retired 8 GeV axis and is 52.2273 since the pair-pT floor moved to 9 GeV,
// 2026-09-08), and a retyped edge in a comment is exactly what goes stale (CLAUDE.md Binnings
// rule 1). Read from the canonical vector, never typed, and independent of the cell mode:
// merging the cells ABOVE it cannot move it, which is why the delivered range is expressed as an
// edge here and looked up in whatever axis is in use, rather than as a bin index.
inline double FirstDeliveredPtEdge()
{
    const std::vector<double> e = PairPtEdges("nomerge");
    return e[e.size() - 4];
}

// Contiguous (lo,hi) ranges -> bin edges (same helper as FillMCTrigEffHists.cxx /
// FillMCTrigEffClosure.cxx; kept local so this header has no dependency on either).
inline std::vector<double> RangesToEdgesLocal(const QEtaBinning& ranges)
{
    std::vector<double> e;
    for (size_t i = 0; i < ranges.size(); ++i) {
        if (i == 0) e.push_back(ranges[i].first);
        else if (std::fabs(ranges[i].first - ranges[i - 1].second) > 1e-6)
            throw std::runtime_error("PairTrigEff::RangesToEdgesLocal: non-contiguous ranges");
        e.push_back(ranges[i].second);
    }
    return e;
}

// The |eta^pair| GROUP edges, derived from the live coarse pair-eta ranges through the SAME
// MakeDrEtaGroups fold the `nocorr_etamerge*` dR-correction modes use. Deriving instead of typing
// is what keeps the two groupings from drifting apart, and it means the outer edge follows the
// fiducial cut automatically (2.4 -> 2.2 on 2026-09-07) instead of being a second place to
// remember (.claude/CLAUDE.md §Binnings item 1).
inline DrAxisGroups AbsEtaGroups()
{
    static const CommonEffcyConfig ecfg{};
    std::vector<double> src = RangesToEdgesLocal(ecfg.pair_eta_proj_ranges_coarse_incl_gap);
    TAxis ax((int)src.size() - 1, src.data());
    DrAxisGroups G = MakeDrEtaGroups(&ax, true);

    // ⚠ GUARD, AND THE REASON IT EXISTS. `MakeDrEtaGroups(..., true)` returns the sign-independent
    // |eta| fold ONLY in the version of dr_correction_cell_groups.h that carries
    // mc_trigeff_dr_binning_approaches.md D11. An older version of that header returns the SIGNED
    // three-region grouping {-2.2, -1.0, 1.0, 2.2} for the same call, and every histogram here is
    // FILLED with |pair_eta| -- so against it the first cell would be permanently empty, |eta| in
    // [1, 2.2) would collapse into one cell, and the labels would read "-2.2..-1": a completely
    // different measurement, with no crash and no warning, because CheckCanonicalBinning() only
    // checks the file against whatever this function returns and is self-consistent either way.
    // (`/review-analysis-code` 2026-09-08 found the repository in exactly that state -- the D11
    // header is a concurrent thread's uncommitted work; see the Progress Log Stage-0 entry.)
    // A folded axis starts at 0 and is non-negative; anything else is not the axis this object is
    // defined on, and must stop the job rather than quietly re-bin it.
    if (!G.folded || G.edges.empty() || std::fabs(G.edges.front()) > 1e-9)
        throw std::runtime_error(
            "PairTrigEff::AbsEtaGroups: MakeDrEtaGroups(..., true) did not return a "
            "sign-independent |eta^pair| fold (first edge "
            + (G.edges.empty() ? std::string("<none>") : std::to_string(G.edges.front()))
            + ", folded=" + (G.folded ? "true" : "false")
            + "). The single-value pair efficiency is binned in |eta^pair| and would be silently "
              "measured on a different axis. Use the dr_correction_cell_groups.h that carries the "
              "D11 |eta| fold (mc_trigeff_dr_binning_approaches.md).");
    return G;
}

// MINIMUM CELL POPULATION for a cell to be DELIVERED through the reader. A single-value
// efficiency is one binomial number, so its whole value is its statistical precision: at
// eps ~ 0.5 the error is sqrt(eps(1-eps)/N), i.e. +-0.07 at N = 50 and +-0.35 at N = 2. Below this
// count the "measurement" is a fluctuation, and delivering it through the same API as a
// 4830-pair cell is a silent hazard -- the forward, highest-pT opposite-sign cell rests on FOUR
// pairs and reads eps = 0.807 +- 0.315, HIGHER than every other cell in its row when the physics
// requires it to be the lowest. Counted on the RAW (unweighted) all-pairs histogram, because it is
// the number of Bernoulli trials that sets the error, not the weighted yield.
// The cell is still MEASURED, WRITTEN and PRINTED -- only its automatic delivery is gated, and
// every gated cell is listed by the fill macro and in this file's provenance stamp.
inline int MinCellPairs() { return 50; }

// A delivered value must be a probability-like quantity. eps^pair = S1/S0 cannot exceed 1 by
// construction; K can, as a fluctuation, and a K > 1 would turn the per-pair weight 1/(e1 e2 K)
// into an ENHANCEMENT where the close-by-RoI physics is a loss. Such a cell is refused rather
// than capped, because at these statistics a K > 1 means the cell is unmeasured, not that the
// correction is really above 1.
inline double MaxDeliveredValue() { return 1.0; }

// The first pair-pT bin (1-based) the measurement is DELIVERED for. The request is about the three
// highest cells; the measurement itself is made in all of them (it costs nothing and the low bins
// are the sanity check against the well-determined dR correction), so the delivered range is a
// separate, explicit statement rather than a truncated axis.
inline int FirstDeliveredPtBin(const std::string& mode = "nomerge")
{
    const std::vector<double> e = PairPtEdges(mode);
    for (size_t i = 0; i + 1 < e.size(); ++i)
        if (std::fabs(e[i] - FirstDeliveredPtEdge()) < 1e-6) return (int)i + 1;
    throw std::runtime_error("PairTrigEff::FirstDeliveredPtBin: the delivered-range edge is not an "
                             "edge of the cell mode's own pair-pT axis");
}

}  // namespace PairTrigEff

// =================================================================================================
// THE READER. One instance = one (sign, mass window) of one file.
// =================================================================================================
class PairTrigEffEvaluator {
public:
    // `form` selects which of the two delivered numbers Eval() returns:
    //   kPure        eps_2mu4^pair -- replaces the WHOLE per-pair weight
    //   kCalibrated  K             -- MULTIPLIES the two single-muon efficiencies
    // There is no default: the two are applied in different places and a wrong choice is a ~25 %
    // normalization error that no MC closure can see (header block above).
    enum class ApplyForm { kPure, kCalibrated };

    void Load(const std::string& file, const std::string& sign, const std::string& window,
              ApplyForm form, const std::string& mode = "nomerge")
    {
        PairTrigEff::CheckSignalWindowMirror();
        sign_ = sign; window_ = window; form_ = form; file_ = file; mode_ = mode;
        TFile* f = TFile::Open(file.c_str(), "READ");
        if (!f || f->IsZombie())
            throw std::runtime_error("PairTrigEffEvaluator: cannot open " + file
                + " -- run FillMCTrigEffPairEff first");
        // BOTH quantities are loaded, whichever one Eval() returns. The delivery gate must be the
        // SAME predicate for the pure and the calibrated form: `den_paireff_<window>`, the coverage
        // denominator the closure books, is built once per window, so a cell that one form
        // delivered and the other refused would leave that form's numerator missing pairs its
        // denominator still counts -- a bias that reads as non-closure. K >= eps always, so only
        // the calibrated form can fail the upper bound; requiring BOTH in range makes the two
        // agree by construction. (No cell differs today: max delivered K = 0.78.)
        h_eps_.reset(Grab(f, PairTrigEff::HistName("eps", sign, window, mode)));
        h_k_  .reset(Grab(f, PairTrigEff::HistName("k",   sign, window, mode)));
        h_.reset(static_cast<TH2D*>((form == ApplyForm::kPure ? h_eps_ : h_k_)->Clone(
            (PairTrigEff::HistName(form == ApplyForm::kPure ? "eps" : "k", sign, window, mode)
             + "_sel").c_str())));
        h_->SetDirectory(nullptr);
        h_den_     .reset(Grab(f, PairTrigEff::HistName("den",      sign, window, mode)));
        h_nraw_    .reset(Grab(f, PairTrigEff::HistName("nraw",     sign, window, mode)));
        h_nrawpass_.reset(Grab(f, PairTrigEff::HistName("nrawpass", sign, window, mode)));
        auto* prov = dynamic_cast<TNamed*>(f->Get("provenance"));
        provenance_ = prov ? prov->GetTitle() : "";
        f->Close();
        CheckCanonicalBinning();
    }

    // Is there a DELIVERABLE number for this pair? A cell qualifies when it is inside the axes, at
    // or above the first delivered pair-pT bin, carries at least PairTrigEff::MinCellPairs() raw
    // pairs, and holds a value in (0, MaxDeliveredValue()]. Callers must ask BEFORE Eval(): a pair
    // outside the delivered region has no single-value correction at all, and inventing one (1, or
    // a neighbouring cell's) would be a silent extrapolation.
    //
    // The three rejection reasons are kept DISTINCT by Reject() below, because they mean different
    // things: "no pairs" is an empty cell, "too few pairs" is a cell that exists but cannot be
    // measured, and "value 0" is a cell where pairs exist and NONE fired -- a real zero efficiency,
    // not an absence. Only the first two are absences of information.
    enum class Reject { kDelivered, kOffAxis, kBelowDeliveredPt, kNoPairs, kTooFewPairs,
                        kZeroValue, kAboveMax };

    Reject Status(double pair_pt, double pair_eta) const
    {
        int ix = 0, iy = 0;
        if (!Cell(pair_pt, pair_eta, ix, iy))              return Reject::kOffAxis;
        if (ix < PairTrigEff::FirstDeliveredPtBin(mode_)) return Reject::kBelowDeliveredPt;
        // Emptiness is decided on the RAW histogram, the same one the fill macro's REFUSED list
        // uses -- the weighted denominator would be a second, near-equivalent definition of the
        // same predicate, and near-equivalent is how these drift apart.
        if (h_nraw_->GetBinContent(ix, iy) <= 0.)          return Reject::kNoPairs;
        if (h_nraw_->GetBinContent(ix, iy) < PairTrigEff::MinCellPairs())
                                                           return Reject::kTooFewPairs;
        if (h_->GetBinContent(ix, iy) <= 0.)               return Reject::kZeroValue;
        if (h_eps_->GetBinContent(ix, iy) > PairTrigEff::MaxDeliveredValue() ||
            h_k_  ->GetBinContent(ix, iy) > PairTrigEff::MaxDeliveredValue())
                                                           return Reject::kAboveMax;
        return Reject::kDelivered;
    }

    static const char* StatusText(Reject r)
    {
        switch (r) {
            case Reject::kDelivered:        return "delivered";
            case Reject::kOffAxis:          return "outside the cell axes";
            case Reject::kBelowDeliveredPt: return "below the delivered pair-pT range";
            case Reject::kNoPairs:          return "no pairs in the cell";
            case Reject::kTooFewPairs:      return "too few raw pairs to measure";
            case Reject::kZeroValue:        return "measured zero (pairs exist, none fired)";
            case Reject::kAboveMax:         return "value above the physical maximum";
        }
        return "?";
    }

    bool Covered(double pair_pt, double pair_eta) const
    {
        return Status(pair_pt, pair_eta) == Reject::kDelivered;
    }

    // The measured number for this pair's cell. Throws for an undelivered pair rather than
    // returning a neutral 1: a silently neutral correction is exactly the failure mode the
    // Covered() gate exists to prevent.
    double Eval(double pair_pt, double pair_eta) const
    {
        const Reject r = Status(pair_pt, pair_eta);
        if (r != Reject::kDelivered)
            throw std::runtime_error(Form("PairTrigEffEvaluator::Eval: no single-value pair "
                "efficiency for pair pT = %.3f GeV, eta^pair = %.3f (sign %s, window %s): %s -- "
                "call Covered() first", pair_pt, pair_eta, sign_.c_str(), window_.c_str(),
                StatusText(r)));
        int ix = 0, iy = 0; Cell(pair_pt, pair_eta, ix, iy);
        return h_->GetBinContent(ix, iy);
    }

    // Gated exactly like Eval(): an error bar quoted for a cell the reader would refuse to deliver
    // would be read as if the cell were usable.
    double Error(double pair_pt, double pair_eta) const
    {
        const Reject r = Status(pair_pt, pair_eta);
        if (r != Reject::kDelivered)
            throw std::runtime_error(Form("PairTrigEffEvaluator::Error: cell not delivered for "
                "pair pT = %.3f GeV, eta^pair = %.3f (sign %s, window %s): %s",
                pair_pt, pair_eta, sign_.c_str(), window_.c_str(), StatusText(r)));
        int ix = 0, iy = 0; Cell(pair_pt, pair_eta, ix, iy);
        return h_->GetBinError(ix, iy);
    }

    // The RAW (unweighted) populations behind a cell -- the statistical reach the weighted value
    // hides. Ungated on purpose: a consumer asking "how many pairs is this built on?" must be able
    // to ask about a cell the reader refuses to deliver.
    double NRaw(double pair_pt, double pair_eta) const
    {
        int ix = 0, iy = 0;
        return Cell(pair_pt, pair_eta, ix, iy) ? h_nraw_->GetBinContent(ix, iy) : 0.;
    }
    double NRawPass(double pair_pt, double pair_eta) const
    {
        int ix = 0, iy = 0;
        return Cell(pair_pt, pair_eta, ix, iy) ? h_nrawpass_->GetBinContent(ix, iy) : 0.;
    }

    const TH2D* Hist() const { return h_.get(); }
    const std::string& Provenance() const { return provenance_; }
    std::string Describe() const
    {
        return "single-value pair 2mu4 efficiency ["
             + std::string(form_ == ApplyForm::kPure ? "PURE eps^pair, replaces the whole weight"
                                                : "CALIBRATED K, multiplies eps(1)eps(2)")
             + "], sign " + sign_ + ", " + PairTrigEff::Window(window_).text
             + ", cells: " + PairTrigEff::Mode(mode_).text
             + Form(", delivered only for cells with >= %d raw pairs and a value in (0, %g]",
                    PairTrigEff::MinCellPairs(), PairTrigEff::MaxDeliveredValue())
             + ", from " + file_;
    }

private:
    static TH2D* Grab(TFile* f, const std::string& n)
    {
        auto* h = dynamic_cast<TH2D*>(f->Get(n.c_str()));
        if (!h) throw std::runtime_error("PairTrigEffEvaluator: missing '" + n + "' in "
                                         + f->GetName());
        h = static_cast<TH2D*>(h->Clone((n + "_loaded").c_str()));
        h->SetDirectory(nullptr);
        return h;
    }

    // 1-based cell indices, or false if the pair is off either axis. The pair-eta lookup is by
    // |eta^pair| BY CONSTRUCTION: the Y axis is the sign-independent fold, so a signed lookup
    // would send every negative-eta pair into the underflow bin -- the exact silent bug D11 of
    // mc_trigeff_dr_binning_approaches.md records for the dR correction's own folded modes.
    bool Cell(double pair_pt, double pair_eta, int& ix, int& iy) const
    {
        const TAxis* ax = h_->GetXaxis();
        const TAxis* ay = h_->GetYaxis();
        const double aeta = std::fabs(pair_eta);
        if (pair_pt < ax->GetXmin() || pair_pt >= ax->GetXmax()) return false;
        if (aeta    < ay->GetXmin() || aeta    >= ay->GetXmax()) return false;
        ix = ax->FindFixBin(pair_pt);
        iy = ay->FindFixBin(aeta);
        return ix >= 1 && ix <= ax->GetNbins() && iy >= 1 && iy <= ay->GetNbins();
    }

    // The file's axes must still BE the canonical ones. A moved edge is silent otherwise: every
    // lookup succeeds and returns a number measured in a different cell (.claude/CLAUDE.md
    // §Binnings). The pair-pT axis is compared against ParamsSet::pair_pt_coarse_bins directly;
    // the |eta| axis is compared against the live fold of
    // CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, built by the SAME MakeDrEtaGroups
    // helper the `nocorr_etamerge*` dR modes use, so the two groupings cannot drift apart.
    void CheckCanonicalBinning() const
    {
        const std::vector<double> px = PairTrigEff::PairPtEdges(mode_);
        if (h_->GetNbinsX() != (int)px.size() - 1)
            throw std::runtime_error("PairTrigEffEvaluator: " + file_ + " has "
                + std::to_string(h_->GetNbinsX()) + " pair-pT cells but cell mode '" + mode_
                + "' of ParamsSet::pair_pt_coarse_bins has " + std::to_string(px.size() - 1)
                + " -- stale file, or the wrong cell mode");
        for (size_t i = 0; i < px.size(); ++i)
            if (std::fabs(h_->GetXaxis()->GetBinLowEdge((int)i + 1) - px[i]) > 1e-6)
                throw std::runtime_error("PairTrigEffEvaluator: pair-pT edge " + std::to_string(i)
                    + " of " + file_ + " is " + std::to_string(h_->GetXaxis()->GetBinLowEdge((int)i + 1))
                    + " but cell mode '" + mode_ + "' of ParamsSet says " + std::to_string(px[i])
                    + " -- stale file");

        const std::vector<double> ay = PairTrigEff::AbsEtaGroups().edges;
        if (h_->GetNbinsY() != (int)ay.size() - 1)
            throw std::runtime_error("PairTrigEffEvaluator: " + file_ + " has "
                + std::to_string(h_->GetNbinsY()) + " |eta^pair| groups but the live fold of "
                "CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap has "
                + std::to_string(ay.size() - 1) + " -- stale file");
        for (size_t i = 0; i < ay.size(); ++i)
            if (std::fabs(h_->GetYaxis()->GetBinLowEdge((int)i + 1) - ay[i]) > 1e-6)
                throw std::runtime_error("PairTrigEffEvaluator: |eta^pair| edge "
                    + std::to_string(i) + " of " + file_ + " is "
                    + std::to_string(h_->GetYaxis()->GetBinLowEdge((int)i + 1))
                    + " but the live fold says " + std::to_string(ay[i]) + " -- stale file "
                    "(the coarse pair-eta outer edge moved 2.4 -> 2.2 on 2026-09-07)");
    }

    std::unique_ptr<TH2D> h_, h_eps_, h_k_, h_den_, h_nraw_, h_nrawpass_;
    std::string sign_, window_, mode_ = "nomerge", file_, provenance_;
    ApplyForm form_ = ApplyForm::kPure;
};

#endif  // PAIR_TRIG_EFF_EVALUATOR_H
