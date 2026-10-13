// =================================================================================================
// syst_pp24_pair_eta_within_cell_shape.cxx
//
// Within-cell eta^pair shape systematic of the reco-binned bin-by-bin correction
// (docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md §3d; option A, user 2026-10-13).
//
// The correction uses ONE eps^reco per coarse eta^pair cell c (canonical axis), while the single-muon
// gap windows carve narrow holes into eta^pair inside some cells. The cell eps is right for data only
// to the extent MC models the reco eta^pair shape inside the cell. Variation: reweight the MC truth
// eta^pair shape inside the cell so that MC reco reproduces the (trigger-corrected) data reco shape in
// the fine bins k (ParamsSet::pair_eta_measured_bins, 0.1), and recompute the cell eps:
//
//     eps_k   = N[R_k] / N[T_k]                    (fine-bin bin-by-bin eps, m-integrated)
//     eps_c   = sum_k N[R_k] / sum_k N[T_k]        (nominal cell eps, same m-integrated scheme)
//     eps'_c  = sum_k D_k / sum_k (D_k / eps_k)    (MC truth weight w_k ∝ (D_k/eps_k)/N[T_k] ⇒ MC reco = D)
//     delta_c = eps_c / eps'_c - 1                 (= relative shift of the corrected cell yield)
//
// D_k = data, crossx selection, trigger-weighted (h_pair_eta_fine_trig); R_k, T_k = MC closure tree
// histograms (reco / truth eta^pair). Systematic = |delta_c| (symmetric). Stat error of delta from toys:
// data = Poisson on the raw pair count n_k of the bin (h_pair_eta_fine_counts), D_k = n_k x the bin's mean
// trigger weight; eps_k Gaussian with a binomial width from N_eff of T_k (approximation: R_k is not a strict
// subset of T_k because of bin migration, reco vs truth eta^pair), and eps_c is RECOMPUTED per toy as sum T_k eps_k / sum T_k from the
// same toy eps_k (eps_c and eps_k are fully correlated; holding eps_c fixed inflates the error).
// delta uses the m-integrated cell eps (pT-summed for merged groups), an approximation to the shift of the
// applied 3D (m-resolved) correction. D_k contains the data background (fake/hadronic pairs) that the
// truth-matched MC R_k does not, so part of a data/MC shape difference may be background composition.
// High pT: for p_T^pair bins with index >= kMergeFromPtBin the fine gap-hole bins have N[R_k] ~ 0
// (eps_k undefined / toys unstable), so delta is taken from the merged range of those bins and
// assigned to each of them.
//
// Inputs: the eta^pair smoothness ROOT files (no refill). Outputs (SmoothnessPlotDir()/<wp>/):
//   within_cell_shape_syst_<wp>WP.png, within_cell_shape_syst_<wp>WP.txt, within_cell_shape_syst_<wp>WP.root
// Usage: root -l -b -q 'syst_pp24_pair_eta_within_cell_shape.cxx+(true)'
// =================================================================================================
#include <cmath>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "pair_eta_smoothness_common.h"

namespace {

const int kMergeFromPtBin = 4;   // canonical pT bin index 4 (= PtEdges()[4]) and above: merged (see header)
const int kNToys = 2000;

TH1D* Sum(TFile* f, const std::string& base, int lo, int hi)
{
    TH1D* s = nullptr;
    for (int p = lo; p <= hi; ++p) {
        auto* h = dynamic_cast<TH1D*>(f->Get((base + "_pt" + std::to_string(p)).c_str()));
        if (!h) throw std::runtime_error("missing " + base + "_pt" + std::to_string(p) + " in " + f->GetName());
        if (!s) { s = static_cast<TH1D*>(h->Clone((base + "_s" + std::to_string(lo) + std::to_string(hi)).c_str())); s->SetDirectory(nullptr); }
        else s->Add(h);
    }
    return s;
}

void CheckFine(const TH1D* h, const std::vector<double>& fine)
{
    if (h->GetNbinsX() != (int)fine.size() - 1) throw std::runtime_error(std::string(h->GetName()) + ": not the fine binning");
    for (int i = 1; i <= h->GetNbinsX() + 1; ++i)
        if (std::fabs(h->GetXaxis()->GetBinLowEdge(i) - fine[i - 1]) > 1e-6)
            throw std::runtime_error(std::string(h->GetName()) + ": fine edge mismatch");
}

struct Delta { double val = NAN, stat = NAN, eps_c = NAN; int toys_used = 0; };

// delta for the fine bins [f, l] of one coarse cell
Delta CellDelta(const TH1D* D, const TH1D* N, const TH1D* R, const TH1D* T, int f, int l, TRandom3& rng)
{
    struct Bin { double d, n, t, eps, seps; };
    std::vector<Bin> b;
    double sR = 0, sT = 0;
    for (int k = f; k <= l; ++k) {
        const double r = R->GetBinContent(k), t = T->GetBinContent(k), et = T->GetBinError(k);
        sR += r; sT += t;
        const double eps = t > 0 ? r / t : 0;
        const double neff = et > 0 ? t * t / (et * et) : 0;
        const double seps = neff > 0 ? std::sqrt(std::max(eps * std::fabs(1 - eps), 0.0) / neff) : 0;
        b.push_back({D->GetBinContent(k), N->GetBinContent(k), t, eps, seps});
    }
    Delta out;
    if (sT <= 0 || sR <= 0) return out;
    out.eps_c = sR / sT;
    auto eval = [&](bool toy) -> double {
        double sD = 0, sDe = 0, sTe = 0, sTt = 0;
        bool undefined = false;
        for (const auto& x : b) {
            const double d = toy ? (x.n > 0 ? rng.Poisson(x.n) * x.d / x.n : 0.0) : x.d;
            const double e = toy ? rng.Gaus(x.eps, x.seps) : x.eps;
            sTe += x.t * e; sTt += x.t;
            if (d <= 0) continue;
            if (e <= 0) undefined = true;   // data where MC predicts no reco pair: eps' undefined
            else { sD += d; sDe += d / e; }
        }
        if (undefined || sD <= 0 || sTt <= 0) return NAN;
        return (sTe / sTt) / (sD / sDe) - 1;
    };
    out.val = eval(false);
    double s1 = 0, s2 = 0;
    for (int i = 0; i < kNToys; ++i) {
        const double z = eval(true);
        if (!std::isfinite(z)) continue;
        s1 += z; s2 += z * z; ++out.toys_used;
    }
    if (out.toys_used > 1) out.stat = std::sqrt(std::max(0.0, s2 / out.toys_used - std::pow(s1 / out.toys_used, 2)));
    return out;
}

}  // namespace

void syst_pp24_pair_eta_within_cell_shape(bool tight_WP = true)
{
    gROOT->SetBatch(kTRUE);
    const std::string wp = tight_WP ? "tight" : "medium";
    const std::vector<double> fine = PairEtaFineSmoothnessBins();
    const std::vector<double> eta = PairRecoEffRecoBinnedEvaluator::EtaEdges();
    const std::vector<double> pt = PairRecoEffRecoBinnedEvaluator::PtEdges();
    const int npt = pt.size() - 1, neta = eta.size() - 1;
    if (kMergeFromPtBin < 1 || kMergeFromPtBin >= npt) throw std::runtime_error("bad kMergeFromPtBin");

    std::unique_ptr<TFile> fd(TFile::Open(SmoothnessDataFile(tight_WP).c_str()));
    std::unique_ptr<TFile> fm(TFile::Open(SmoothnessMCFile().c_str()));
    if (!fd || fd->IsZombie() || !fm || fm->IsZombie()) throw std::runtime_error("cannot open smoothness inputs");

    // groups: each canonical pT bin below the merge, the merged high-pT range, all pT
    struct Group { int lo, hi; std::string label; };
    std::vector<Group> groups;
    for (int p = 0; p < kMergeFromPtBin; ++p)
        groups.push_back({p, p, Form("%.4g < p_{T}^{pair} < %.4g GeV", pt[p], pt[p + 1])});
    groups.push_back({kMergeFromPtBin, npt - 1, Form("%.4g < p_{T}^{pair} < %.4g GeV", pt[kMergeFromPtBin], pt[npt])});
    groups.push_back({0, npt - 1, Form("%.4g < p_{T}^{pair} < %.4g GeV (all)", pt[0], pt[npt])});

    TRandom3 rng(20261013);
    std::vector<std::vector<Delta>> res(groups.size(), std::vector<Delta>(neta));
    for (size_t g = 0; g < groups.size(); ++g) {
        TH1D* D = Sum(fd.get(), "h_pair_eta_fine_trig", groups[g].lo, groups[g].hi);
        TH1D* N = Sum(fd.get(), "h_pair_eta_fine_counts", groups[g].lo, groups[g].hi);
        TH1D* R = Sum(fm.get(), "h_mcR_fine_raw_" + wp, groups[g].lo, groups[g].hi);
        TH1D* T = Sum(fm.get(), "h_mcT_fine", groups[g].lo, groups[g].hi);
        CheckFine(D, fine); CheckFine(N, fine); CheckFine(R, fine); CheckFine(T, fine);
        for (int c = 0; c < neta; ++c) {
            int f = -1, l = -1;
            for (int k = 1; k <= D->GetNbinsX(); ++k) {
                const double x = D->GetBinCenter(k);
                if (x > eta[c] && x < eta[c + 1]) { if (f < 0) f = k; l = k; }
            }
            res[g][c] = CellDelta(D, N, R, T, f, l, rng);
        }
    }

    // per canonical (pT, eta) cell: assigned delta (merged above kMergeFromPtBin)
    const std::string odir = SmoothnessPlotDir() + wp + "/";
    gSystem->mkdir(odir.c_str(), kTRUE);
    TFile fo((odir + "within_cell_shape_syst_" + wp + "WP.root").c_str(), "RECREATE");
    TH2D h2("h_within_cell_shape_delta", ";p_{T}^{pair} [GeV];#eta^{pair};#delta", npt, pt.data(), neta, eta.data());
    TH2D h2s("h_within_cell_shape_syst", ";p_{T}^{pair} [GeV];#eta^{pair};|#delta|", npt, pt.data(), neta, eta.data());
    for (int p = 0; p < npt; ++p)
        for (int c = 0; c < neta; ++c) {
            const Delta& d = res[p < kMergeFromPtBin ? p : kMergeFromPtBin][c];
            h2.SetBinContent(p + 1, c + 1, d.val); h2.SetBinError(p + 1, c + 1, d.stat);
            h2s.SetBinContent(p + 1, c + 1, std::fabs(d.val)); h2s.SetBinError(p + 1, c + 1, d.stat);
        }
    h2.Write(); h2s.Write();
    TNamed("definition", "delta = eps_c/eps'_c - 1, eps'_c = sum D_k / sum D_k/eps_k (fine 0.1 eta^pair, m-integrated); "
           "syst = |delta|; pT bins >= kMergeFromPtBin share the merged-range delta; "
           "docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md 3d").Write();
    fo.Close();

    std::ofstream txt(odir + "within_cell_shape_syst_" + wp + "WP.txt");
    txt << "# within-cell eta^pair shape systematic, " << wp << " WP; delta = eps_c/eps'_c - 1 (corrected-yield shift)\n"
        << "# rows: pT group; columns: coarse eta^pair cells; entry delta (toy stat) [toys used of " << kNToys << "]\n";
    txt << "# eta cells:";
    for (int c = 0; c < neta; ++c) txt << Form(" [%g,%g]", eta[c], eta[c + 1]);
    txt << "\n";
    for (size_t g = 0; g < groups.size(); ++g) {
        txt << Form("pT %6.2f-%6.2f:", pt[groups[g].lo], pt[groups[g].hi + 1]);
        for (int c = 0; c < neta; ++c)
            txt << Form(" %+7.4f(%.4f)[%d]", res[g][c].val, res[g][c].stat, res[g][c].toys_used);
        txt << "\n";
    }
    txt.close();

    // figure: delta vs eta^pair, one panel per group (nrows >= ncols)
    const int np = groups.size();
    const int ncol = (np <= 3) ? np : 2, nrow = (np + ncol - 1) / ncol;
    gStyle->SetOptStat(0);
    TCanvas cv("cv", "", 520 * ncol, 360 * nrow + 150);
    TPad* head = new TPad("head", "", 0, 1 - 150.0 / (360.0 * nrow + 150), 1, 1);
    TPad* body = new TPad("body", "", 0, 0, 1, 1 - 150.0 / (360.0 * nrow + 150));
    head->Draw(); body->Draw();
    head->cd();
    {
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.165);
        t.DrawLatex(0.02, 0.76, "#it{pp} 2024 data + Pythia8 full simulation, #sqrt{s} = 5.36 TeV, opposite-sign pairs, signal region");
        t.DrawLatex(0.02, 0.44, "#delta = #varepsilon_{cell} / #varepsilon'_{cell} #minus 1,   #varepsilon'_{cell} = #Sigma_{k} D_{k} / #Sigma_{k} (D_{k}/#varepsilon_{k}),   k: 0.1-wide #eta^{pair} bins in the cell");
        t.DrawLatex(0.02, 0.12, "D_{k}: trigger-corrected data;   #varepsilon_{k} = N^{MC}_{reco,k} / N^{MC}_{truth,k};   #varepsilon_{cell} = #Sigma_{k} N^{MC}_{reco,k} / #Sigma_{k} N^{MC}_{truth,k}");
        t.SetTextAlign(31); t.DrawLatex(0.98, 0.12, tight_WP ? "Tight muons" : "Medium muons");
    }
    body->cd();
    body->Divide(ncol, nrow);
    std::vector<TH1D*> keep;
    for (int g = 0; g < np; ++g) {
        body->cd(g + 1);
        gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.10); gPad->SetTicks(1, 1);
        TH1D* h = new TH1D(Form("hd%d", g), "", neta, eta.data());
        keep.push_back(h);
        double lo = -0.02, hi = 0.02;
        for (int c = 0; c < neta; ++c) {
            const Delta& d = res[g][c];
            if (!std::isfinite(d.val)) continue;
            h->SetBinContent(c + 1, d.val); h->SetBinError(c + 1, std::isfinite(d.stat) ? d.stat : 0);
            lo = std::min(lo, d.val - h->GetBinError(c + 1)); hi = std::max(hi, d.val + h->GetBinError(c + 1));
        }
        const double pad = 0.15 * (hi - lo);
        h->GetYaxis()->SetRangeUser(lo - pad, hi + pad);
        h->GetXaxis()->SetTitle("#eta^{pair}"); h->GetYaxis()->SetTitle("#delta (corrected-yield shift)");
        h->GetXaxis()->SetTitleSize(0.060); h->GetXaxis()->SetLabelSize(0.050);
        h->GetYaxis()->SetTitleSize(0.060); h->GetYaxis()->SetLabelSize(0.050); h->GetYaxis()->SetTitleOffset(1.30);
        h->SetMarkerStyle(20); h->SetMarkerSize(0.9); h->SetLineColor(kAzure + 2); h->SetMarkerColor(kAzure + 2);
        h->Draw("PE");
        TLine* z = new TLine(eta.front(), 0, eta.back(), 0); z->SetLineStyle(2); z->Draw();
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.058); t.SetTextAlign(31);
        t.DrawLatex(0.95, 0.92, groups[g].label.c_str());
    }
    cv.SaveAs((odir + "within_cell_shape_syst_" + wp + "WP.png").c_str());
}
