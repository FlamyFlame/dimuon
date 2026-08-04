// dr_correction_ratio.h
//
// The ONE implementation of the round-7 conditional (binomial-correct) error on an
// inverse-weighted efficiency ratio, plus the (pair pT, pair eta) cell projection built on top
// of it. Shared by plot_mc_trig_eff.cxx (which measures eps_dR / eps_single) and by
// fit_dr_corrections.cxx / plot_dr_correction_fits.cxx (which fit and draw them).
//
// It lives in a header because the fit is judged on chi2/ndf: a second, drifted copy of the
// error formula would silently change every chi2 and therefore every method comparison.
//
// =============================================================================
// CONDITIONAL (binomial-correct) ERROR ON AN INVERSE-WEIGHTED EFFICIENCY RATIO  (round 7)
//
// R(dR) = N/D with D = sum_all w  and  N = sum_fired w/eps  is an EFFICIENCY: the numerator
// is a re-weighted SUBSET of the denominator. `TH1::Divide` without option "B" propagates
//     e_R = R sqrt((e_N/N)^2 + (e_D/D)^2)
// which assumes num and den are INDEPENDENT. They are not, and the resulting bars are too
// long by ~sqrt((1+eps)/(1-eps)) -- 1.5x to 3.2x here, worst in the highest pair-pT bin where
// eps is largest. Symptom: chi2/ndf of a constant fit over the plateau ~0.25 instead of ~1
// (mc_trigger_efficiency.md R12).
//
// Conditioning on the MC sample (D fixed; only the Bernoulli trigger decisions fluctuate),
//     Var(N) = sum_i a_i^2 p_i (1-p_i) + 2 sum_pairs a_1 a_2 (p_12 - p_1 p_2),  a_i = w_i/eps_i
// estimated from the histograms booked in FillMCTrigEffHists.cxx as
//     Var = A - R*B  (+ covP - R^2*covQ  for Step 4, where both legs of a pair share a dR bin)
//     e_R = sqrt(Var)/D
// Unweighted limit (w=1, eps=1): A=B=N -> Var = N(1-R) -> e_R = sqrt(R(1-R)/D), as it must be.
// covP/covQ are absent for Step 3 (one entry = one pair = one Bernoulli trial) -> pass nullptr.
// =============================================================================

#ifndef DR_CORRECTION_RATIO_H
#define DR_CORRECTION_RATIO_H

#include <TH1D.h>
#include <TH3D.h>

#include <algorithm>
#include <cmath>

inline void SetConditionalRatioErrors(TH1D* r, const TH1D* den, const TH1D* A, const TH1D* B,
                                      const TH1D* covP = nullptr, const TH1D* covQ = nullptr)
{
    for (int i = 1; i <= r->GetNbinsX(); ++i) {
        const double D = den->GetBinContent(i);
        if (D <= 0.) { r->SetBinError(i, 0.); continue; }
        const double R  = r->GetBinContent(i);
        const double a  = A->GetBinContent(i), b = B->GetBinContent(i);
        const double cp = covP ? covP->GetBinContent(i) : 0.;
        const double cq = covQ ? covQ->GetBinContent(i) : 0.;
        double var = a - R * b + cp - R * R * cq;

        // BOUNDARY CASE. When every effective entry in the bin fired (p_i -> 1) the conditional
        // binomial variance genuinely vanishes -- the k=n binomial artefact -- and `var` comes
        // out as 0, slightly negative, or a catastrophic cancellation of terms many orders of
        // magnitude larger (seen in the near-empty high-pair-pT cells of the 10k-event overlay
        // TEST sample: A=1.1e-06 vs var=1e-22). A ~0 error is NOT safe: the plateau weighted
        // mean weights by 1/e^2, so such a bin would either be dropped (e=0) or completely
        // dominate the mean (e=1e-7). Detect it by comparing var with the SCALE of the terms
        // that built it, and fall back to the "1/n rule" (the 68% bound for k=n) with the
        // effective denominator count n_eff = (D/e_D)^2, e_D = sqrt(sum w^2) from Sumw2.
        const double scale = a + R * b + std::fabs(cp) + R * R * cq;
        if (var > 1e-6 * scale) { r->SetBinError(i, std::sqrt(var) / D); continue; }

        // BOTH binomial boundaries: k = n gives err ~ R/n_eff, k = 0 gives err ~ 1/n_eff.
        // Using R/n_eff alone returns EXACTLY 0 when R = 0 (denominator > 0, nothing fired),
        // and both consumers drop zero-error points -- silently deleting a genuine
        // zero-efficiency dR bin from the fit and the plot. max(R,1) covers both.
        const double eD = den->GetBinError(i);
        const double neff = (eD > 0.) ? (D / eD) * (D / eD) : 1.0;
        r->SetBinError(i, (neff > 0.) ? std::max(R, 1.0) / neff : 0.);
    }
}

// eps_dR(dR) (Step 3) or eps_single(dR) (Step 4) for ONE (pair pT, pair eta) cell of the 3D
// histograms. iy = pair-pT bin, iz = pair-eta bin; **0 means "integrate over that axis"**, so
// (iy=0, iz=0) is the fully inclusive curve. `hp`/`hq` (the Step-4 leg-leg covariance terms) may
// be null. The returned TH1D is detached from any file and owned by the caller.
inline TH1D* DrCellRatio(TH3D* hn, TH3D* hd, TH3D* ha, TH3D* hb, TH3D* hp, TH3D* hq,
                         int iy, int iz, const char* nm)
{
    const int npt = hn->GetYaxis()->GetNbins(), neta = hn->GetZaxis()->GetNbins();
    const int ylo = (iy == 0) ? 1 : iy, yhi = (iy == 0) ? npt  : iy;
    const int zlo = (iz == 0) ? 1 : iz, zhi = (iz == 0) ? neta : iz;
    TH1D* n = hn->ProjectionX(Form("%s_n", nm), ylo, yhi, zlo, zhi, "e");
    TH1D* d = hd->ProjectionX(Form("%s_d", nm), ylo, yhi, zlo, zhi, "e");
    TH1D* a = ha->ProjectionX(Form("%s_a", nm), ylo, yhi, zlo, zhi, "e");
    TH1D* b = hb->ProjectionX(Form("%s_b", nm), ylo, yhi, zlo, zhi, "e");
    TH1D* p = hp ? hp->ProjectionX(Form("%s_p", nm), ylo, yhi, zlo, zhi, "e") : nullptr;
    TH1D* q = hq ? hq->ProjectionX(Form("%s_q", nm), ylo, yhi, zlo, zhi, "e") : nullptr;
    auto* r = (TH1D*)n->Clone(nm);
    r->SetDirectory(nullptr);
    r->Divide(d);
    SetConditionalRatioErrors(r, d, a, b, p, q);
    delete n; delete d; delete a; delete b; delete p; delete q;
    return r;
}

#endif // DR_CORRECTION_RATIO_H
