// =================================================================================================
// plot_forward_qeta_edge_scan.cxx
//
// DECISION PLOT: where should the forward edge of the q*eta fiducial gap cut sit?
//
// The forward acceptance collapses one-sidedly -- muon_gap_cuts_acceptance.md F6 measures the pp
// yield at q*eta = 2.38-2.40 as only 8.8% of its 2.20-2.22 value, while the NEGATIVE side still
// holds ~71% out to -2.4. The question this figure answered was whether the forward edge had to
// sit at 2.20 or could be loosened to 2.25 / 2.30, recovering 1.4% (pp) / 2.2% (PbPb) of muons.
// It was loosened to 2.30 on 2026-08-04 -- and the user moved it BACK to 2.20 on 2026-09-07:
// `ParamsSet::single_mu_fiducial_gap_cuts` now rejects q*eta in [2.20, 2.40], so the edge again
// coincides with the historical per-muon `q*eta < 2.2` signal cut. This macro is the supporting
// evidence for that scan, kept live so the choice can be re-checked; it does not itself set the
// edge, which is read from ParamsSet.
//
// WHAT THIS DRAWS: the single-muon mu4 efficiency vs pT in the MOST POSITIVE q*eta bin, with the
// bin's UPPER edge scanned over 2.20 / 2.25 / 2.30 / 2.40 (lower edge fixed at 2.0). If the 2.40
// curve is visibly degraded relative to 2.20 but 2.25 / 2.30 are not, the cut can be loosened.
// (2.40 is known to be bad -- it is drawn as the reference for "how bad".)
//
// LAYOUT
//   pp   (unchanged): ONE canvas, 4 panels -- mu+ LEFT, mu- RIGHT, data TOP, MC BOTTOM.
//   PbPb (2026-08-13, user): DATA ONLY, TWO canvases, each 3 rows x 2 columns.
//        The HIJING overlay was dropped from the PbPb figure: it is still a TEST sample, so its
//        statistics in a single forward q*eta slice are too thin to say anything about where the
//        edge belongs, and a half-empty MC row only invited over-reading. The two data canvases
//        use the axis that PbPb actually has and pp does not -- year and centrality:
//        (1) `..._by_year.png`      -- centrality-integrated; mu+ LEFT, mu- RIGHT;
//                                      rows = PbPb 2023 / 2024 / 2025.
//                                      Answers: is the forward edge behaving the same in all
//                                      three running periods, per charge?
//        (2) `..._centrality.png`   -- PbPb 23+24+25 summed AND mu+ + mu- summed, one panel per
//                                      centrality bin of `ParamsSet::ctrbins`
//                                      (0-5 / 5-10 / 10-20 / 20-30 / 30-50 / 50-80 %) = exactly
//                                      the 3 x 2 grid. Answers: does occupancy move the forward
//                                      edge, i.e. is the loosened window safe in central events?
//        Centrality binning is READ FROM `ParamsSet::ctrbins`, never retyped (CLAUDE.md
//        BLOCKING binning rule); the histogram key `_ctr<lo>_<hi>_` and the panel label are both
//        built from that one vector, so they cannot drift apart.
//
//        TWO q*eta GRIDS ARE IN PLAY, and the two canvases handle it differently ON PURPOSE.
//        The q*eta binning of the tag-and-probe histograms is registered PER CENTRALITY BIN in
//        `hist_binning_map`, and the peripheral bins are deliberately coarser: 184 bins for
//        0-5 / 5-10 / 10-20 / 20-30 %, 102 for 30-50 %, 61 for 50-80 %.
//          - `_by_year` reads only the centrality-INTEGRATED histograms, which are on the same
//            184-bin grid in all three years, so it scans the canonical candidate ladder.
//          - `_centrality` spans all three grids at once, so it scans only the boundaries COMMON
//            to all six panels (`CommonQEtaBoundaries`) -- otherwise the same coloured curve
//            would silently cover a different q*eta window in different panels, which is exactly
//            what CLAUDE.md §Binnings forbids. Every legend entry is written from the edge the
//            axis actually delivered, so no label can disagree with the data behind it.
//
// IMPORTANT -- INPUTS MUST HAVE NO FORWARD GAP CUT. The nominal chain now rejects q*eta > 2.3
// outright, so its outputs contain nothing to compare. This macro therefore reads files produced
// WITHOUT the gap cut (the un-suffixed `_fine_q_eta_bin` tag-and-probe outputs; the gap-cut
// variant carries `_qeta_fid`), and refuses to run if the region above 2.2 is empty (which is
// exactly what a gap-cut input would look like) rather than silently drawing identical curves.
//
// The efficiency is formed by projecting the 2D (q*eta, pT) numerator and denominator over each
// candidate q*eta range and dividing -- so the scan is independent of the fitted q*eta binning.
//   data: h_pt2nd_vs_q_eta2nd[_ctr<lo>_<hi>]_<sign1|sign2>_2mu4_sepr   (numerator, probe fired)
//         h_pt2nd_vs_q_eta2nd[_ctr<lo>_<hi>]_<sign1|sign2>_mu4_sepr    (denominator, tag fired)
//         i.e. the tag-and-probe pair {"_2mu4","_mu4"} of RDFBasedHistFillingData::trigs_pair.
//   MC  : h_mc_pt_vs_q_eta_num_<muplus|muminus>        /  h_mc_pt_vs_q_eta_denom_<chg>   (pp only)
// Both are (x = q*eta, y = pT).
//
// Errors: TGraphAsymmErrors::BayesDivide, matching how the nominal turn-on points are built.
//
// Run (from this directory):
//   root -l -b -q 'plot_forward_qeta_edge_scan.cxx+("pp")'
//   root -l -b -q 'plot_forward_qeta_edge_scan.cxx+("pbpb")'
// =================================================================================================

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "dr_correction_sample_cfg.h"
#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"

namespace {

// The candidate upper edges. Lower edge fixed at 2.0 (the last coarse q*eta bin's lower edge).
// These are REQUESTS: the q*eta axis has 0.02-wide bins above 2.20, so 2.20 / 2.30 / 2.40 land
// on bin boundaries but 2.25 does not (it sits inside [2.24, 2.26)) and resolves to 2.26. Every
// label is drawn from the achieved edge, not from this list -- see ResolveQEtaRange.
const double kQEtaLo = 2.0;
const std::vector<double> kUpperEdges = {2.20, 2.25, 2.30, 2.40};
const std::vector<Color_t> kEdgeCol   = {kBlack, kBlue + 1, kGreen + 2, kRed + 1};
const std::vector<Style_t> kEdgeMark  = {20, 21, 22, 23};

// The three PbPb running periods, in the order they are drawn as rows.
const std::vector<int> kPbPbYears = {2023, 2024, 2025};

template <typename T>
T* GetObj(TFile* f, const std::string& name) {
    auto* o = dynamic_cast<T*>(f->Get(name.c_str()));
    if (!o) throw std::runtime_error("plot_forward_qeta_edge_scan: missing '" + name +
                                     "' in " + f->GetName());
    return o;
}

// Detached clone, so the histogram survives closing its file and never collides on name.
TH2D* CloneDetached(TH2D* h, const std::string& name) {
    auto* c = static_cast<TH2D*>(h->Clone(name.c_str()));
    c->SetDirectory(nullptr);
    return c;
}

// Resolve the half-open q*eta window [lo, hi) onto the histogram's bin grid, and report the
// upper edge that was ACTUALLY achieved.
//
// This matters, and it bit this macro: the fine q*eta axis is NOT uniform -- it is 0.10 wide up
// to 2.20 and 0.02 wide above it (edges ... 2.10, 2.20, 2.22, 2.24, 2.26 ...). So 2.20, 2.30 and
// 2.40 are real bin boundaries but **2.25 is not**: it falls inside [2.24, 2.26), and a
// projection asked to stop "below 2.25" necessarily stops at 2.26. Labelling that curve "2.25"
// puts a number on the canvas that the drawn data does not have (CLAUDE.md §Binnings rule 5:
// if the label and the binning disagree, one of them is a bug). We therefore always return the
// achieved edge and label the curve with THAT -- the text is derived from the axis, so it cannot
// drift away from it.
void ResolveQEtaRange(const TAxis* ax, double lo, double hi, int& b1, int& b2,
                      double& lo_achieved, double& hi_achieved) {
    // half-open [lo, hi): FindBin(hi) would include the bin containing hi, so step back.
    b1 = ax->FindBin(lo + 1e-6);
    b2 = ax->FindBin(hi - 1e-6);
    if (b2 < b1) b2 = b1;
    // BOTH edges come back from the axis. 2.0 happens to be a boundary on all three grids today,
    // but half-applying the principle is how the 2.25 defect survived in the first place.
    lo_achieved = ax->GetBinLowEdge(b1);
    hi_achieved = ax->GetBinUpEdge(b2);
}

// eff(pT) for q*eta in [lo, hi), by projecting both 2D hists onto pT over that x-range.
// `hi_achieved` returns the bin-grid upper edge actually used (see ResolveQEtaRange).
TGraphAsymmErrors* EffInQEtaRange(TH2D* num, TH2D* den, double lo, double hi,
                                  const std::string& tag, double& n_den_entries,
                                  double& lo_achieved, double& hi_achieved) {
    int b1 = 0, b2 = 0;
    ResolveQEtaRange(den->GetXaxis(), lo, hi, b1, b2, lo_achieved, hi_achieved);
    TH1D* n = num->ProjectionY((tag + "_n").c_str(), b1, b2, "e");
    TH1D* d = den->ProjectionY((tag + "_d").c_str(), b1, b2, "e");
    n_den_entries = d->Integral();
    auto* g = new TGraphAsymmErrors();
    g->BayesDivide(n, d);
    g->SetName((tag + "_g").c_str());
    delete n; delete d;
    return g;
}

// The q*eta bin boundaries in [lo_scan, hi_scan] that are shared by EVERY histogram in `dens`.
//
// Needed because the per-centrality q*eta binning is NOT uniform across centrality: it is
// registered per bin in `hist_binning_map` (binName + "_ctr<lo>_<hi>") and the peripheral bins
// are deliberately coarser -- 184 bins for 0-5 / 5-10 / 10-20 / 20-30 %, 102 for 30-50 %, 61 for
// 50-80 %. A single candidate edge therefore does NOT mean the same window in every panel: 2.30
// exists on the 184-bin grid but on the 61-bin grid the neighbouring boundaries are 2.28 and
// 2.36. Drawing "the same" curve in six panels while it silently spans different q*eta in each
// is precisely the failure CLAUDE.md §Binnings exists to prevent, so the centrality figure
// scans the boundaries COMMON to all six panels ({2.20, 2.28, 2.36, 2.40} for the current
// registration) -- identical physical windows everywhere, and derived from the axes rather than
// retyped, so it follows the registration if it ever changes.
std::vector<double> CommonQEtaBoundaries(const std::vector<TH2D*>& dens,
                                         double lo_scan, double hi_scan) {
    const double tol = 1e-4;
    std::vector<double> common;
    const TAxis* a0 = dens.front()->GetXaxis();
    for (int i = 1; i <= a0->GetNbins() + 1; ++i) {
        const double e = a0->GetBinLowEdge(i);
        if (e < lo_scan - tol || e > hi_scan + tol) continue;
        bool in_all = true;
        for (size_t k = 1; k < dens.size() && in_all; ++k) {
            const TAxis* a = dens[k]->GetXaxis();
            bool found = false;
            for (int j = 1; j <= a->GetNbins() + 1 && !found; ++j)
                found = std::fabs(a->GetBinLowEdge(j) - e) < tol;
            in_all = found;
        }
        if (in_all) common.push_back(e);
    }
    if (common.size() < 2)
        throw std::runtime_error("plot_forward_qeta_edge_scan: fewer than 2 q*eta boundaries are "
                                 "common to all centrality panels -- cannot scan the edge on a "
                                 "single grid.");
    return common;
}

void DrawPanel(TH2D* num, TH2D* den, const std::string& panel_title, const std::string& tag,
               bool draw_legend, const std::vector<double>& upper_edges) {
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    // pT axis: range taken straight from the input histogram, which already carries the CANONICAL
    // single-muon trigger-efficiency pT binning (41 log bins, 4-60 GeV). A sanity / special-
    // selection plot must never invent its own binning or range -- it exists to be compared
    // against the regular efficiency plots, so it has to share their axis. Log binning => log
    // axis (memory `feedback_log_scale_plots`).
    const TAxis* pax = den->GetYaxis();
    gPad->SetLogx();
    auto* fr = gPad->DrawFrame(pax->GetXmin(), 0.0, pax->GetXmax(), 1.15);
    fr->GetXaxis()->SetTitle("p_{T} [GeV]");
    fr->GetYaxis()->SetTitle("#varepsilon(mu4)");
    fr->GetXaxis()->SetTitleSize(0.05);
    fr->GetYaxis()->SetTitleSize(0.05);
    // 4-60 GeV spans well under two decades, so the default log-x labelling draws a single "10".
    fr->GetXaxis()->SetMoreLogLabels();
    fr->GetXaxis()->SetNoExponent();

    auto* leg = new TLegend(0.45, 0.18, 0.93, 0.44);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);

    for (size_t i = 0; i < upper_edges.size(); ++i) {
        double nden = 0, lo = 0, hi = 0;
        auto* g = EffInQEtaRange(num, den, kQEtaLo, upper_edges[i],
                                 tag + "_" + std::to_string(i), nden, lo, hi);
        g->SetMarkerStyle(kEdgeMark[i % kEdgeMark.size()]);
        g->SetMarkerColor(kEdgeCol[i % kEdgeCol.size()]);
        g->SetLineColor(kEdgeCol[i % kEdgeCol.size()]);
        g->SetLineWidth(2);
        g->Draw("P same");
        // Label the edges the axis actually gave us, never the ones we asked for.
        leg->AddEntry(g, Form("%.2f < q#eta < %.2f", lo, hi), "lp");
    }
    if (draw_legend) leg->Draw();

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.048);
    t.DrawLatex(0.17, 0.90, panel_title.c_str());
}

// Guard: an input that already has the forward gap cut applied has NOTHING above q*eta = 2.2,
// which would silently render identical curves and invite exactly the wrong conclusion.
void AssertForwardRegionPopulated(TH2D* den, const std::string& what) {
    const TAxis* ax = den->GetXaxis();
    const double above = den->Integral(ax->FindBin(2.2 + 1e-6), ax->FindBin(2.4 - 1e-6),
                                       1, den->GetNbinsY());
    if (above <= 0.)
        throw std::runtime_error(
            "plot_forward_qeta_edge_scan: " + what + " has NO entries with q*eta > 2.2. That is "
            "what a gap-cut input looks like, and every candidate edge would draw the same curve. "
            "Re-run the producer with the forward fiducial gap cut DISABLED "
            "(RDFBasedHistFillingData::apply_fiducial_gap_cut = false).");
}

// Plateau efficiency (pT > 8 GeV) and probe count per candidate edge -- the numbers behind the
// figure, so the decision can be read off values rather than eyeballed off a curve.
void PrintPlateauRow(const std::string& label, TH2D* num, TH2D* den,
                     const std::vector<double>& upper_edges) {
    std::cout << "    " << label << ":";
    for (double hi_req : upper_edges) {
        int b1 = 0, b2 = 0;
        double lo = 0, hi = 0;
        ResolveQEtaRange(den->GetXaxis(), kQEtaLo, hi_req, b1, b2, lo, hi);
        TH1D* n = num->ProjectionY("tn", b1, b2, "e");
        TH1D* d = den->ProjectionY("td", b1, b2, "e");
        const int p1 = d->FindBin(8.0 + 1e-6);
        const double N = n->Integral(p1, d->GetNbinsX());
        const double D = d->Integral(p1, d->GetNbinsX());
        std::cout << Form("  [%.2f,%.2f) %.4f (D=%.4g)", lo, hi, D > 0 ? N / D : 0., D);
        delete n; delete d;
    }
    std::cout << "\n";
}

// ---------------------------------------------------------------------------------------------
// PbPb data plumbing
// ---------------------------------------------------------------------------------------------

std::string PbPbDataFile(int year, const std::string& wp) {
    return std::string(gSystem->Getenv("HOME")) +
           "/usatlasdata/dimuon_data/pbpb_" + std::to_string(year) +
           "/histograms_real_pairs_pbpb_" + std::to_string(year) +
           "_single_mu4_fine_q_eta_bin" + wp + ".root";
}

// `_ctr<lo>_<hi>` token as written by RDFBasedHistFillingPbPb, built from ParamsSet::ctrbins.
std::string CtrToken(size_t ibin) {
    return "_ctr" + std::to_string(ParamsSet::ctrbins[ibin]) + "_" +
           std::to_string(ParamsSet::ctrbins[ibin + 1]);
}
std::string CtrLabel(size_t ibin) {
    return Form("%d-%d%%", ParamsSet::ctrbins[ibin], ParamsSet::ctrbins[ibin + 1]);
}

// Tag-and-probe numerator / denominator name for one (centrality token, charge) cell.
// Empty `ctr` = centrality-integrated. sign1 = mu+, sign2 = mu-.
std::string TnPName(const std::string& ctr, const std::string& sign, bool numerator) {
    return "h_pt2nd_vs_q_eta2nd" + ctr + "_" + sign + (numerator ? "_2mu4_sepr" : "_mu4_sepr");
}

const std::string kPbPbHeader = "Pb+Pb  #sqrt{s_{NN}} = 5.36 TeV,  ";
const std::string kOutName    = "single_mu_eff_forward_qeta_edge_scan";

// The forward edge the analysis actually adopted, read from the cut vector itself so the canvas
// can never quote a value the code no longer applies. Drawn on the figure because it is the one
// number a reader needs to connect this scan to the decision it supports -- the adopted edge is
// not necessarily one of the drawn curves (in the centrality figure the common q*eta grid
// brackets it between two curves), so without it the figure cannot be mapped onto the cut.
double AdoptedForwardEdge() {
    double lo = 0.;
    for (const auto& w : ParamsSet::single_mu_fiducial_gap_cuts)
        if (w.first > lo) lo = w.first;   // the most forward window's LOWER edge = the cut
    return lo;
}
std::string AdoptedForwardEdgeLabel() {
    return Form(",  adopted edge q#eta = %.2f", AdoptedForwardEdge());
}

std::string PbPbOutDir(bool use_tight_wp) {
    return std::string(gSystem->Getenv("HOME")) + "/usatlasdata/dimuon_data/plots/"
           "pbpb_trigger_efficiency/mc_based" +
           std::string(use_tight_wp ? "" : "_medium") + "/forward_qeta_edge_scan/";
}

void DrawHeader(const std::string& text) {
    TLatex head;
    head.SetNDC();
    head.SetTextFont(42);
    head.SetTextSize(0.018);
    head.DrawLatex(0.03, 0.986, text.c_str());
}

// ---------------------------------------------------------------------------------------------
// PbPb canvas 1: centrality-integrated, mu+ / mu- columns, one row per running period.
// ---------------------------------------------------------------------------------------------
void PlotPbPbByYear(bool use_tight_wp, const std::string& wp, const std::string& wpt) {
    TCanvas c("c_fwd_edge_year", "", 1400, 1650);
    c.Divide(2, static_cast<int>(kPbPbYears.size()));

    std::vector<std::pair<std::string, std::pair<TH2D*, TH2D*>>> rows;  // for the number table

    for (size_t iy = 0; iy < kPbPbYears.size(); ++iy) {
        const int year = kPbPbYears[iy];
        const std::string fname = PbPbDataFile(year, wp);
        std::unique_ptr<TFile> f(TFile::Open(fname.c_str(), "READ"));
        if (!f || f->IsZombie()) throw std::runtime_error("cannot open data file " + fname);

        for (size_t is = 0; is < 2; ++is) {
            const std::string sign  = (is == 0) ? "sign1" : "sign2";
            const std::string chg   = (is == 0) ? "#mu^{+}" : "#mu^{-}";
            const std::string tag   = Form("y%d_%s", year, sign.c_str());
            auto* num = CloneDetached(GetObj<TH2D>(f.get(), TnPName("", sign, true)),  tag + "_n2");
            auto* den = CloneDetached(GetObj<TH2D>(f.get(), TnPName("", sign, false)), tag + "_d2");
            AssertForwardRegionPopulated(den, Form("PbPb %d %s", year, chg.c_str()));

            c.cd(static_cast<int>(2 * iy + is + 1));
            DrawPanel(num, den, Form("Pb+Pb %d,  %s", year, chg.c_str()), tag,
                      /*draw_legend=*/(iy == 0 && is == 0), kUpperEdges);
            rows.emplace_back(Form("PbPb %d %s", year, (is == 0) ? "mu+" : "mu-"),
                              std::make_pair(num, den));
        }
        f->Close();
    }

    c.cd(0);
    // NOT "0-80%": the centrality-INTEGRATED histograms carry no centrality filter at all, so
    // they also hold the events beyond 80% and those with no centrality assignment (~0.3%).
    DrawHeader(kPbPbHeader + wpt + " muons,  forward q#eta edge scan,  "
               "data, all centralities" + AdoptedForwardEdgeLabel());

    const std::string out_dir = PbPbOutDir(use_tight_wp);
    gSystem->mkdir(out_dir.c_str(), kTRUE);
    const std::string png = out_dir + kOutName + "_by_year.png";
    c.SaveAs(png.c_str());
    std::cout << "  wrote " << png << std::endl;

    std::cout << "\n  [by year] plateau eff (pT > 8 GeV) and probe count per candidate edge:\n";
    for (auto& r : rows) PrintPlateauRow(r.first, r.second.first, r.second.second, kUpperEdges);
}

// ---------------------------------------------------------------------------------------------
// PbPb canvas 2: 2023+2024+2025 summed and mu+ + mu- summed, one panel per centrality bin.
// ---------------------------------------------------------------------------------------------
void PlotPbPbByCentrality(bool use_tight_wp, const std::string& wp, const std::string& wpt) {
    const size_t nctr = ParamsSet::ctrbins.size() - 1;

    // Sum over the three years AND over both charges, per centrality bin.
    std::vector<TH2D*> num(nctr, nullptr), den(nctr, nullptr);
    for (int year : kPbPbYears) {
        const std::string fname = PbPbDataFile(year, wp);
        std::unique_ptr<TFile> f(TFile::Open(fname.c_str(), "READ"));
        if (!f || f->IsZombie()) throw std::runtime_error("cannot open data file " + fname);

        for (size_t ic = 0; ic < nctr; ++ic) {
            const std::string ctr = CtrToken(ic);
            for (const std::string& sign : {std::string("sign1"), std::string("sign2")}) {
                auto* n = GetObj<TH2D>(f.get(), TnPName(ctr, sign, true));
                auto* d = GetObj<TH2D>(f.get(), TnPName(ctr, sign, false));
                if (!num[ic]) {
                    num[ic] = CloneDetached(n, Form("sum_n_c%zu", ic));
                    den[ic] = CloneDetached(d, Form("sum_d_c%zu", ic));
                } else {
                    // Add() refuses (and only prints) when the axes differ, which would silently
                    // drop a whole year from the sum. Make that fatal instead.
                    if (!num[ic]->Add(n) || !den[ic]->Add(d))
                        throw std::runtime_error(
                            "plot_forward_qeta_edge_scan: cannot sum year " + std::to_string(year) +
                            " into centrality bin " + CtrLabel(ic) + " -- the q*eta/pT binning "
                            "differs between years, so the combined figure would silently omit it.");
                }
            }
        }
        f->Close();
    }

    // Guard on the summed sample: with all years and both charges in, an empty forward region
    // can only mean gap-cut inputs.
    AssertForwardRegionPopulated(den[0], "PbPb 23+24+25, both charges, most central bin");

    // The six panels do NOT share one q*eta grid (184 / 102 / 61 bins, per-centrality by
    // registration), so scan only the boundaries all six have in common -- otherwise the same
    // coloured curve would span a different q*eta window in different panels.
    const std::vector<double> edges = CommonQEtaBoundaries(den, kUpperEdges.front(),
                                                           kUpperEdges.back());
    std::cout << "  [centrality] q*eta boundaries common to all " << nctr << " panels:";
    for (double e : edges) std::cout << Form(" %.2f", e);
    std::cout << "\n";

    // nrows >= ncols, nrows ~ sqrt(N) (memory `feedback_subplot_layout`); the canonical
    // ParamsSet::ctrbins has 6 bins, i.e. exactly the 3 x 2 grid the user asked for.
    const int ncols = 2;
    const int nrows = static_cast<int>((nctr + ncols - 1) / ncols);
    TCanvas c("c_fwd_edge_ctr", "", 700 * ncols, 550 * nrows);
    c.Divide(ncols, nrows);

    for (size_t ic = 0; ic < nctr; ++ic) {
        c.cd(static_cast<int>(ic) + 1);
        DrawPanel(num[ic], den[ic], "Centrality " + CtrLabel(ic), Form("ctr%zu", ic),
                  /*draw_legend=*/(ic == 0), edges);
    }

    c.cd(0);
    DrawHeader(kPbPbHeader + wpt + " muons,  forward q#eta edge scan,  "
               "data 2023+2024+2025,  #mu^{+} + #mu^{-}" + AdoptedForwardEdgeLabel());

    const std::string out_dir = PbPbOutDir(use_tight_wp);
    gSystem->mkdir(out_dir.c_str(), kTRUE);
    const std::string png = out_dir + kOutName + "_centrality.png";
    c.SaveAs(png.c_str());
    std::cout << "  wrote " << png << std::endl;

    std::cout << "\n  [centrality] plateau eff (pT > 8 GeV) and probe count per candidate edge:\n";
    for (size_t ic = 0; ic < nctr; ++ic)
        PrintPlateauRow("centrality " + CtrLabel(ic), num[ic], den[ic], edges);
}

// ---------------------------------------------------------------------------------------------
// pp canvas (unchanged): data / MC rows, mu+ / mu- columns.
// ---------------------------------------------------------------------------------------------
void PlotPP(bool use_tight_wp, const std::string& wp, const std::string& wpt) {
    const std::string home = std::string(gSystem->Getenv("HOME"));
    const std::string data_file = home + "/usatlasdata/dimuon_data/pp_2024/"
        "histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin" + wp + ".root";
    const std::string mc_file = home + "/usatlasdata/pythia_fullsim_full_sample/"
        "mc_trig_eff_hists_pp24_full" + wp + "_nogapcut.root";

    TFile* fd = TFile::Open(data_file.c_str(), "READ");
    if (!fd || fd->IsZombie())
        throw std::runtime_error("cannot open data file " + data_file);
    TFile* fm = TFile::Open(mc_file.c_str(), "READ");
    if (!fm || fm->IsZombie())
        throw std::runtime_error(
            "cannot open MC file " + mc_file +
            " -- produce it by running FillMCTrigEffHists with the gap cut DISABLED");

    // data: sign1 = mu+, sign2 = mu-   |   MC: muplus / muminus
    TH2D* dnum_p = GetObj<TH2D>(fd, TnPName("", "sign1", true));
    TH2D* dden_p = GetObj<TH2D>(fd, TnPName("", "sign1", false));
    TH2D* dnum_m = GetObj<TH2D>(fd, TnPName("", "sign2", true));
    TH2D* dden_m = GetObj<TH2D>(fd, TnPName("", "sign2", false));
    TH2D* mnum_p = GetObj<TH2D>(fm, "h_mc_pt_vs_q_eta_num_muplus");
    TH2D* mden_p = GetObj<TH2D>(fm, "h_mc_pt_vs_q_eta_denom_muplus");
    TH2D* mnum_m = GetObj<TH2D>(fm, "h_mc_pt_vs_q_eta_num_muminus");
    TH2D* mden_m = GetObj<TH2D>(fm, "h_mc_pt_vs_q_eta_denom_muminus");

    AssertForwardRegionPopulated(dden_p, "data mu+");
    AssertForwardRegionPopulated(mden_p, "MC mu+");

    TCanvas c("c_fwd_edge", "", 1400, 1100);
    c.Divide(2, 2);
    c.cd(1); DrawPanel(dnum_p, dden_p, "data, #mu^{+}", "d_p", true,  kUpperEdges);
    c.cd(2); DrawPanel(dnum_m, dden_m, "data, #mu^{-}", "d_m", false, kUpperEdges);
    c.cd(3); DrawPanel(mnum_p, mden_p, "MC, #mu^{+}",   "m_p", false, kUpperEdges);
    c.cd(4); DrawPanel(mnum_m, mden_m, "MC, #mu^{-}",   "m_m", false, kUpperEdges);

    c.cd(0);
    DrawHeader(Form("pp  #sqrt{s} = 5.36 TeV,  %s muons,  forward q#eta edge scan", wpt.c_str())
               + AdoptedForwardEdgeLabel());

    const std::string out_dir = std::string(gSystem->Getenv("HOME")) +
        "/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based" +
        std::string(use_tight_wp ? "" : "_medium") + "/forward_qeta_edge_scan/";
    gSystem->mkdir(out_dir.c_str(), kTRUE);
    const std::string png = out_dir + kOutName + ".png";
    c.SaveAs(png.c_str());
    std::cout << "  wrote " << png << std::endl;

    std::cout << "\n  plateau eff (pT > 8 GeV) and probe count per candidate edge:\n";
    PrintPlateauRow("data mu+", dnum_p, dden_p, kUpperEdges);
    PrintPlateauRow("data mu-", dnum_m, dden_m, kUpperEdges);
    PrintPlateauRow("MC   mu+", mnum_p, mden_p, kUpperEdges);
    PrintPlateauRow("MC   mu-", mnum_m, mden_m, kUpperEdges);
    fd->Close();
    fm->Close();
}

}  // namespace

void plot_forward_qeta_edge_scan(const std::string& sample = "pp", bool use_tight_wp = true) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const bool is_pbpb = (sample == "pbpb");
    if (!is_pbpb && sample != "pp")
        throw std::runtime_error("plot_forward_qeta_edge_scan: sample must be 'pp' or 'pbpb'");

    // Working point exposed on every plot set (memory `feedback_plots_wp_config_var`); nominal
    // Tight is un-suffixed, Medium routes to its own file and its own plot subtree.
    const std::string wp  = use_tight_wp ? "" : "_medium_wp";
    const std::string wpt = use_tight_wp ? "Tight" : "Medium";

    if (is_pbpb) {
        PlotPbPbByYear(use_tight_wp, wp, wpt);
        PlotPbPbByCentrality(use_tight_wp, wp, wpt);
    } else {
        PlotPP(use_tight_wp, wp, wpt);
    }
}
