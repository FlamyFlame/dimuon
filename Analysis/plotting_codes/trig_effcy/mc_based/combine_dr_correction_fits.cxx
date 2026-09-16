// combine_dr_correction_fits.cxx
//
// DELIVERABLE-PACKAGING stage of the DeltaR-correction chain (user, 2026-09-03; parent
// docs/tracking/mc_trigeff_dr_binning_approaches.md D11). Not a measurement or a fit -- it opens
// the per-METHOD fit files fit_dr_corrections.cxx already wrote (one per `expo` / `polyu_fixedRp`
// / `interp`) for ONE (sample, WP, step, sign, plateau mode) and re-packages them into a SINGLE
// ROOT file that carries, for every (pair pT, pair eta) fit cell:
//   * the RAW measured eps_dR points (the inverse-weighting result, before any fit) -- IDENTICAL
//     across the three method files by construction (same source histograms, same window,
//     same DrGroupCellRatio projection; only the SUBSEQUENT fit differs), so it is taken from one
//     file and cross-checked point-by-point against the other two, THROWING on any disagreement;
//   * the fitted TF1 for `expo` and `polyu_fixedRp` -- persisted with its exact formula string and
//     current (fitted) parameter values, exactly as fit_dr_corrections.cxx wrote it;
//   * the interpolation knots TGraph for `interp` (no TF1: the method has no closed-form formula);
//   * the per-method fit_ok map, whose axes ARE the canonical (pair pT, pair eta) cell binning of
//     this plateau mode (folded |eta| axis for the two eta-merge modes) -- so a consumer can read
//     the cell edges directly off any of the three maps without re-deriving the grouping.
// A cell absent from a method's file (too few measured points to attempt that method's fit) is
// simply absent here too -- this macro invents nothing.
//
// Nothing here changes any measurement: every object it writes is copied, not recomputed.

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "dr_correction_sample_cfg.h"

namespace {

// A TF1 built with NUMERIC-indexed parameters ("[0]", "[1]", ...) and later renamed with
// SetParName is exactly what fit_dr_corrections.cxx does -- and reading such a TF1 back and then
// WRITING IT AGAIN in the SAME process (TF1::Clone + Write, or a plain re-Write of the object
// TFile::Get() returned) corrupts its persisted parameter-name bookkeeping in this ROOT build
// (6.34.04): the re-serialized object comes back with the wrong parameter COUNT and
// GetParName()/GetExpFormula() failing, even though GetParameter(i)/Eval() still return the
// right numbers. Verified empirically 2026-09-03 (scratch tests): a FRESH construction (`new
// TF1` with the same formula string, values copied by GetParameter(i)/GetParError(i), names set
// via SetParName(i, ...), one Write()) reproduces fit_dr_corrections.cxx's own working
// first-generation write and survives a clean round trip. So this macro NEVER Clones/re-Writes a
// TF1 object it read from a fit file; it always rebuilds one from scratch, using the SAME
// formula string and parameter names fit_dr_corrections.cxx used for that method's NOCORR fit
// (MakeMethodCfg / the per-method SetParNames block there) -- every plateau mode this macro
// handles is a nocorr* one (checked below), so only the nocorr formulas are needed here.
struct NocorrMethodFormula {
    std::string              formula;
    std::vector<std::string> parNames;
    int                       fixedParIndex = -1;   // -1 = none; else that parameter is FIXED
};

const std::map<std::string, NocorrMethodFormula>& NocorrMethodFormulas()
{
    static const std::map<std::string, NocorrMethodFormula> m = {
        {"expo", {"[3]+[0]*TMath::Exp(-TMath::Power(x/[1],[2]))",
                  {"A", "#lambda", "p", "C"}, -1}},
        // p0 is A = f(0)-C, NOT a2: the Step-3 restriction of 2026-09-08 required the
        // constrained combination a2+a3+a4 to BE a parameter (fit_dr_corrections.cxx
        // MakeMethodCfg). Same function family, a3/a4 unchanged.
        {"polyu_fixedRp", {"[4]+TMath::Power(TMath::Max(0.,1.-x/[3]),2)*([0]+[1]*(TMath::Max(0.,1.-x/[3])"
                           "-1)+[2]*(TMath::Power(TMath::Max(0.,1.-x/[3]),2)-1))",
                           {"A", "a_{3}", "a_{4}", "R_{p}", "C"}, 3}},
    };
    return m;
}

// Rebuilds a NEW TF1 with the same formula, parameter values/errors and (for polyu_fixedRp) the
// same fixed parameter as `src` -- see the header comment on why this replaces Clone().
TF1* RebuildNocorrTF1(const TF1* src, const std::string& method, const std::string& newName)
{
    const auto it = NocorrMethodFormulas().find(method);
    if (it == NocorrMethodFormulas().end())
        throw std::runtime_error("combine_dr_correction_fits: no known NOCORR formula for "
                                 "method '" + method + "' -- add it to NocorrMethodFormulas()");
    const auto& mf = it->second;
    if ((int)mf.parNames.size() != src->GetNpar())
        throw std::runtime_error("combine_dr_correction_fits: method '" + method + "' TF1 '" +
                                 src->GetName() + "' has " + std::to_string(src->GetNpar()) +
                                 " parameters but the known NOCORR formula for it expects " +
                                 std::to_string(mf.parNames.size()) +
                                 " -- MakeMethodCfg (fit_dr_corrections.cxx) must have changed");
    auto* nf = new TF1(newName.c_str(), mf.formula.c_str(), src->GetXmin(), src->GetXmax());
    for (int i = 0; i < src->GetNpar(); ++i) {
        nf->SetParameter(i, src->GetParameter(i));
        nf->SetParError(i, src->GetParError(i));
        nf->SetParName(i, mf.parNames[i].c_str());
    }
    if (mf.fixedParIndex >= 0) nf->FixParameter(mf.fixedParIndex, src->GetParameter(mf.fixedParIndex));
    // The formula above is a COPY of the producer's, and the parameter-count check alone cannot
    // catch it drifting: a re-parametrized shape (as `polyu_fixedRp` was on 2026-09-08) keeps its
    // parameter count while changing what p0 MEANS, and the rebuilt TF1 would then be silently
    // wrong. So compare the rebuild against the source it copied, point by point, over the whole
    // stored range -- the source TF1 read back from the fit file evaluates correctly even in the
    // ROOT build whose name bookkeeping is broken (see the header comment).
    for (int k = 0; k <= 100; ++k) {
        const double x  = src->GetXmin() + (src->GetXmax() - src->GetXmin()) * k / 100.0;
        const double a  = src->Eval(x), b = nf->Eval(x);
        if (std::fabs(a - b) > 1e-9 * std::max(1e-3, std::fabs(a)))
            throw std::runtime_error(std::string("combine_dr_correction_fits: the rebuilt TF1 for "
                "method '") + method + "' does not reproduce the fitted one it copied (at dR=" +
                std::to_string(x) + ": " + std::to_string(a) + " vs " + std::to_string(b) +
                ") -- EITHER NocorrMethodFormulas() has drifted from fit_dr_corrections.cxx"
                "::MakeMethodCfg, OR this fit file was written under an older parametrization"
                " (polyu_fixedRp was re-written in A = f(0)-C on 2026-09-08, same parameter"
                " COUNT, different meaning of p0) -- re-run fit_dr_corrections for it");
    }
    return nf;
}

std::string CellName(const std::string& base, int step, int iy, int iz)
{
    if (iy == 0 && iz == 0) return Form("%s_step%d_incl", base.c_str(), step);
    return Form("%s_step%d_pt%d_eta%d", base.c_str(), step, iy, iz);
}

// Point-by-point comparison of two TGraphErrors that are supposed to be the SAME raw measurement,
// re-derived independently in two method files. Any disagreement means the "raw points are
// method-independent" assumption this macro relies on is wrong, and it must not silently pick one.
void CheckGraphsMatch(const TGraphErrors* ref, const TGraphErrors* other,
                      const std::string& refLabel, const std::string& otherLabel,
                      const std::string& cellTag)
{
    if (!other) return;   // that method simply has no graph for this cell -- not a mismatch
    if (ref->GetN() != other->GetN())
        throw std::runtime_error("combine_dr_correction_fits: raw-graph point count differs for "
                                 "cell " + cellTag + " between " + refLabel + " (" +
                                 std::to_string(ref->GetN()) + ") and " + otherLabel + " (" +
                                 std::to_string(other->GetN()) + ")");
    for (int i = 0; i < ref->GetN(); ++i) {
        double xr, yr, xo, yo;
        ref->GetPoint(i, xr, yr);
        other->GetPoint(i, xo, yo);
        const double er = ref->GetErrorY(i), eo = other->GetErrorY(i);
        if (std::fabs(xr - xo) > 1e-9 || std::fabs(yr - yo) > 1e-9 * std::max(1.0, std::fabs(yr))
            || std::fabs(er - eo) > 1e-9 * std::max(1.0, er))
            throw std::runtime_error("combine_dr_correction_fits: raw-graph point " +
                                     std::to_string(i) + " differs for cell " + cellTag +
                                     " between " + refLabel + " and " + otherLabel +
                                     " -- the 'raw points are method-independent' assumption "
                                     "does not hold, refusing to silently pick one");
    }
}

}  // namespace

// sample/use_tight_wp/step/sign/plateau_mode : identical meaning to fit_dr_corrections.cxx.
// `methods` : which per-method fit files to pull from; default is the repo's standard trio.
void combine_dr_correction_fits(const std::string& sample = "pp_full", bool use_tight_wp = true,
                                int step = 3, const std::string& sign = "os",
                                const std::string& plateau_mode = "nocorr_etamerge",
                                std::vector<std::string> methods =
                                    {"expo", "polyu_fixedRp", "interp"},
                                int overlay_year = 24)
{
    gROOT->SetBatch(kTRUE);

    if (!DrCorrModeMergeEta(plateau_mode))
        throw std::runtime_error("combine_dr_correction_fits: written for the pair-eta-merged "
                                 "modes ('nocorr_etamerge' / 'nocorr_etamerge_ptmerge'), got '" +
                                 plateau_mode + "'");

    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp, overlay_year);

    // ---- open every method's fit file, keep them open for the whole run -----------------------
    struct MethodFile {
        std::string method;
        TFile*      f = nullptr;
        TH2D*       fit_ok = nullptr;
    };
    std::vector<MethodFile> mf;
    for (const auto& m : methods) {
        const std::string path = DrCorrFitFile(cfg, use_tight_wp, step, m, sign, plateau_mode);
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie())
            throw std::runtime_error("combine_dr_correction_fits: cannot open " + path +
                                     " -- run fit_dr_corrections for method '" + m + "' first");
        auto* ok = dynamic_cast<TH2D*>(f->Get(Form("h_step%d_fit_ok", step)));
        if (!ok)
            throw std::runtime_error("combine_dr_correction_fits: no h_step" +
                                     std::to_string(step) + "_fit_ok in " + path);
        mf.push_back({m, f, ok});
        std::cout << "  opened [" << m << "] " << path << std::endl;
    }

    const int npt  = mf[0].fit_ok->GetNbinsX();
    const int neta = mf[0].fit_ok->GetNbinsY();
    for (size_t i = 1; i < mf.size(); ++i)
        if (mf[i].fit_ok->GetNbinsX() != npt || mf[i].fit_ok->GetNbinsY() != neta)
            throw std::runtime_error("combine_dr_correction_fits: method '" + mf[i].method +
                                     "' has a " + std::to_string(mf[i].fit_ok->GetNbinsX()) + "x" +
                                     std::to_string(mf[i].fit_ok->GetNbinsY()) +
                                     " cell grid but method '" + mf[0].method + "' has " +
                                     std::to_string(npt) + "x" + std::to_string(neta) +
                                     " -- the fit files disagree about the mode's binning");

    const std::string out_dir = cfg.out_base + "step" + std::to_string(step) + "_dr_fit/"
                              + DrCorrPlateauModeDir(plateau_mode);
    gSystem->mkdir(out_dir.c_str(), true);
    const std::string wp_tag = use_tight_wp ? "tight" : "medium";
    const std::string out_path = out_dir + "dr_correction_data_and_fits_" + wp_tag + "_" +
                                 (sign.empty() ? "signintgr" : sign) + ".root";
    TFile* fout = TFile::Open(out_path.c_str(), "RECREATE");
    fout->cd();

    int n_cells_written = 0, n_graph_mismatch_checks = 0;
    for (int iy = 0; iy <= npt; ++iy) {
        for (int iz = 0; iz <= neta; ++iz) {
            const bool inclusive = (iy == 0 && iz == 0);
            if (!inclusive && (iy == 0 || iz == 0)) continue;   // only the full inclusive cell

            const std::string cellTag = inclusive ? "incl" : Form("pt%d_eta%d", iy, iz);
            const std::string gname = CellName("g", step, iy, iz);

            // ---- the raw points: read from the FIRST method that has them, cross-checked ------
            TGraphErrors* gref = nullptr;
            std::string   gref_method;
            for (auto& m : mf) {
                auto* g = dynamic_cast<TGraphErrors*>(m.f->Get(gname.c_str()));
                if (g && !gref) { gref = g; gref_method = m.method; }
            }
            if (!gref) continue;   // no method attempted a fit here (too few measured points)

            for (auto& m : mf) {
                if (m.method == gref_method) continue;
                auto* g = dynamic_cast<TGraphErrors*>(m.f->Get(gname.c_str()));
                CheckGraphsMatch(gref, g, gref_method, m.method, cellTag);
                ++n_graph_mismatch_checks;
            }
            fout->cd();
            auto* gwrite = (TGraphErrors*)gref->Clone(("g_" + cellTag).c_str());
            gwrite->SetTitle(Form("raw #varepsilon_{#DeltaR} (inverse weighting); #DeltaR; %s",
                                  cfg.eps_dr_text.c_str()));
            gwrite->Write();
            delete gwrite;

            // ---- the fitted objects, one per method ------------------------------------------
            for (auto& m : mf) {
                const std::string fname = CellName("f", step, iy, iz);
                auto* f = dynamic_cast<TF1*>(m.f->Get(fname.c_str()));
                if (f) {
                    TF1* fwrite = RebuildNocorrTF1(f, m.method, "f_" + m.method + "_" + cellTag);
                    fout->cd();
                    fwrite->Write();
                    delete fwrite;
                    continue;
                }
                const std::string kname = CellName("gknots", step, iy, iz);
                auto* gk = dynamic_cast<TGraph*>(m.f->Get(kname.c_str()));
                if (gk) {
                    fout->cd();
                    auto* gkwrite = (TGraph*)gk->Clone(("gknots_" + m.method + "_" + cellTag).c_str());
                    gkwrite->Write();
                    delete gkwrite;
                }
            }
            ++n_cells_written;
        }
    }

    // ---- the per-method fit_ok maps: also the canonical (pair pT, pair eta) axes of this mode --
    fout->cd();
    for (auto& m : mf) {
        auto* h = (TH2D*)m.fit_ok->Clone(("h_fit_ok_" + m.method).c_str());
        h->SetDirectory(fout);
        h->Write();
    }

    std::string methods_csv;
    for (size_t i = 0; i < mf.size(); ++i) methods_csv += (i ? "," : "") + mf[i].method;
    TNamed("provenance",
          Form("sample=%s (%s); WP=%s; step=%d; sign=%s; plateau_mode=%s; methods=%s; "
               "cells=%dx%d (pair pT x pair eta); cells_written=%d; raw-graph cross-checks=%d",
               sample.c_str(), cfg.mc_label.c_str(), wp_tag.c_str(), step,
               (sign.empty() ? "sign-integrated" : sign.c_str()), plateau_mode.c_str(),
               methods_csv.c_str(), npt, neta, n_cells_written, n_graph_mismatch_checks))
        .Write();

    fout->Write();
    fout->Close();
    for (auto& m : mf) m.f->Close();

    std::cout << "combine_dr_correction_fits: wrote " << n_cells_written << " cells ("
              << npt << "x" << neta << " grid + inclusive) to " << out_path
              << " -- " << n_graph_mismatch_checks
              << " raw-graph cross-method checks, all matched" << std::endl;
}
