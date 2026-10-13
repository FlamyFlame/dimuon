#ifndef RECO_BINNED_CELLS_H
#define RECO_BINNED_CELLS_H
// =================================================================================================
// reco_binned_cells.h -- ratios over SETS of (pair pT, eta^pair, m_mumu) cells from the bin-migration
// file (docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md §2, §3a). Shared by the builder and
// the plotter so the 3D product and every plotted projection come from ONE formula.
//
// Cell flattening of the migration file (fill_pp24_fullsim_pair_bin_migration.C):
//   k = (i_pt * N_eta + i_eta) * N_m + i_m, 0-based, histogram bin k + 1.
//
// For a set S of cells (one 3D cell, a (pT, eta) column summed over m, a projection, everything):
//   N[R_S] = sum_{j in S} R_j            (reco in S, truth anywhere -- incl. out-of-T feed-in)
//   N[T_S] = sum_{i in S} T_i
//   A      = sum_{i in S, j in S} M_ij   (truth in S AND reco in S)
//   B = N[R_S] - A,  D = N[T_S] - A      (disjoint from A and from each other)
//   eps^reco_S = (A + B) / (A + D), linear propagation with Sum w^2 of each piece:
//     var = ((D - B)/(A+D)^2)^2 s2A + (1/(A+D))^2 s2B + ((A+B)/(A+D)^2)^2 s2D,
//     s2B = Sum w^2(R_S) - Sum w^2(A), s2D = Sum w^2(T_S) - Sum w^2(A).
//   Ratio of SUMMED counts -- never a projection of a ratio.
// Truth-binned comparison series: eps^truth_S = N[T_S and R](truth) / N[T_S] (TR subset of T),
//   weighted-binomial var = [(1 - 2r) s2_TR + r^2 s2_T] / T^2.
// =================================================================================================
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

struct MigrationCells {
    int npt = 0, neta = 0, nm = 0;
    const TH1D* T = nullptr;    // h_cell_minv_T          (truth cell)
    const TH1D* R = nullptr;    // h_cell_minv_R_<wp>     (reco cell, truth anywhere)
    const TH1D* TR = nullptr;   // h_cell_minv_TR_<wp>    (truth cell, reco selected anywhere)
    const TH2D* M = nullptr;    // h_cell_minv_M_<wp>     (X truth cell, Y reco cell)
    int Cell(int ipt, int ieta, int im) const { return (ipt * neta + ieta) * nm + im; }   // 0-based
};

struct CellRatio {
    double eps = 0, err = 0;       // eps^reco_S
    double nT = 0, nR = 0, nA = 0; // N[T_S], N[R_S], A
    double s2T = 0, s2R = 0, s2A = 0;
    double eps_truth = 0, err_truth = 0, nTR = 0;   // truth-binned eps
    bool measured = false;         // N[T_S] > 0 and N[R_S] > 0
    bool truth_measured = false;   // N[T_S] > 0
    double NeffT() const { return s2T > 0 ? nT * nT / s2T : 0.; }
};

inline CellRatio RatioOverCells(const MigrationCells& C, const std::vector<int>& cells)
{
    CellRatio r;
    double s2TR = 0;
    for (int k : cells) {
        r.nT  += C.T->GetBinContent(k + 1);   r.s2T  += std::pow(C.T->GetBinError(k + 1), 2);
        r.nR  += C.R->GetBinContent(k + 1);   r.s2R  += std::pow(C.R->GetBinError(k + 1), 2);
        r.nTR += C.TR->GetBinContent(k + 1);  s2TR   += std::pow(C.TR->GetBinError(k + 1), 2);
        for (int l : cells) {
            r.nA  += C.M->GetBinContent(k + 1, l + 1);
            r.s2A += std::pow(C.M->GetBinError(k + 1, l + 1), 2);
        }
    }
    if (r.nT > 0) {
        r.truth_measured = true;
        r.eps_truth = r.nTR / r.nT;
        const double v = ((1 - 2 * r.eps_truth) * s2TR + r.eps_truth * r.eps_truth * r.s2T) / (r.nT * r.nT);
        r.err_truth = std::sqrt(std::max(0.0, v));
    }
    if (r.nT > 0 && r.nR > 0) {
        r.measured = true;
        const double A = r.nA, B = r.nR - r.nA, D = r.nT - r.nA;
        const double s2B = std::max(0.0, r.s2R - r.s2A), s2D = std::max(0.0, r.s2T - r.s2A);
        const double den = A + D;
        r.eps = (A + B) / den;
        const double dA = (D - B) / (den * den), dB = 1.0 / den, dD = -(A + B) / (den * den);
        r.err = std::sqrt(dA * dA * r.s2A + dB * dB * s2B + dD * dD * s2D);
    }
    return r;
}

// Cells selected by inclusive 0-based index ranges on each axis.
inline std::vector<int> CellsIn(const MigrationCells& C, int pt_lo, int pt_hi, int eta_lo, int eta_hi,
                                int m_lo, int m_hi)
{
    std::vector<int> v;
    for (int i = pt_lo; i <= pt_hi; ++i)
        for (int j = eta_lo; j <= eta_hi; ++j)
            for (int k = m_lo; k <= m_hi; ++k) v.push_back(C.Cell(i, j, k));
    return v;
}

inline MigrationCells LoadMigrationCells(TFile* f, const std::string& wp, int npt, int neta, int nm)
{
    MigrationCells C;
    C.npt = npt; C.neta = neta; C.nm = nm;
    auto get1 = [&](const std::string& n) {
        auto* h = dynamic_cast<TH1D*>(f->Get(n.c_str()));
        if (!h) throw std::runtime_error("reco_binned_cells: missing " + n + " in " + f->GetName());
        if (h->GetNbinsX() != npt * neta * nm)
            throw std::runtime_error("reco_binned_cells: " + n + " has " + std::to_string(h->GetNbinsX())
                                     + " cells, canonical axes give " + std::to_string(npt * neta * nm));
        return h;
    };
    C.T  = get1("h_cell_minv_T");
    C.R  = get1("h_cell_minv_R_" + wp);
    C.TR = get1("h_cell_minv_TR_" + wp);
    C.M  = dynamic_cast<TH2D*>(f->Get(("h_cell_minv_M_" + wp).c_str()));
    if (!C.M || C.M->GetNbinsX() != npt * neta * nm || C.M->GetNbinsY() != npt * neta * nm)
        throw std::runtime_error("reco_binned_cells: missing/inconsistent h_cell_minv_M_" + wp);
    return C;
}

#endif
