// plot_muon_q_eta_pt_dependence.cxx
//
// p_T dependence of the single-muon q*eta SHAPE, in pp24 and in Pb+Pb centrality slices,
// with the adopted fiducial gap windows overlaid.
// See Analysis/docs/tracking/muon_gap_cuts_acceptance.md (Step 8).
//
// WHY THIS FIGURE EXISTS: the legacy gap cut was p_T-DEPENDENT -- PassSingleMuonGapCut
// (RDFBasedHistFillingData.cxx:141) applies ParamsSet::charge_eta_gap_cuts only below
// 6 GeV -- whereas the adopted fiducial cut (ParamsSet::single_mu_fiducial_gap_cuts) is
// p_T-INDEPENDENT by design, so that the acceptance factors as a pure function of q*eta.
// That choice is only justified if the q*eta SHAPE is p_T-stable where the windows sit.
// This figure is the direct test: if the dips move, deepen or vanish with p_T, the
// p_T-independent cut is mis-specified.
//
// INPUT (NTuple-processing OUTPUT -- never the raw NTUPs):
//   single_muon_trees_<sample>_part<N>_<trig>_mindR_0_02.root : TTree "muon_tree".
//   Same population as plot_muon_q_eta_spectrum.cxx: muons belonging to a pair that
//   survives the full nominal pair chain, de-duplicated by muon index. NO trigger
//   efficiency correction is applied -- the low-p_T slices in particular are shaped by
//   the mu4 turn-on, which is exactly why they are compared as SHAPES and not yields.
//
// BINNINGS -- both read from their single sources of truth, nothing retyped:
//   q*eta   : ParamsSet::makeEtaTrigEffcyBinning(1), the same canonical fine axis every
//             other panel in muon_gap_cuts/ uses (1D/2D/3D views of one quantity share
//             one binning). Variable width => everything is a DENSITY, Scale(1.,"width").
//             CAVEAT for the reader: over |q*eta| in [1.4,2.2] the canonical step is 0.10,
//             coarser than the real structure there, so the apparent step up at
//             |q*eta| ~ 1.4 is the binning change, not detector structure, and a feature
//             position cannot be read off this figure to better than 0.1 in that region.
//             The canonical binning is NOT changed to suit a plot.
//   centrality : ParamsSet::ctrbins = {0,5,10,20,30,50,80} with 0-5 and 5-10 merged into
//             0-10, identical to the sibling macro's Plot C.
//   p_T     : {4.5, 5, 6, +inf} GeV. These three slices are the USER'S EXPLICIT CHOICE
//             for this figure (2026-09-03) and are deliberately NOT
//             ParamsSet::single_mu_pt_coarse_bins (READ it from ParamsSet; its values are not
//             repeated here, because a retyped copy is what goes stale -- this comment used to
//             claim {4,8,14,25,100} and the canonical vector now starts at 4.5): the question
//             here was the behaviour just above the 4 GeV threshold, where the canonical coarse
//             binning has a single first bin. They are a plot-local diagnostic axis and
//             feed no physics result, no fit and no correction.
//
//             RESOLVED 2026-09-09 (user): the [4, 4.5) slice is DROPPED. It would have been empty
//             by construction after the threshold moved to 4.5 GeV (an NTuple-stage cut), while
//             still drawing a legend entry. The original question -- what the [4,4.5) muons do --
//             is the one the 4.5 GeV adoption answered (muon_gap_cuts_acceptance.md F14).
//
// NORMALISATION: each curve is a unit-area PDF over the plotted q*eta range,
//   Scale(1.,"width") then divide by Integral("width"). Shapes only -- the p_T slices
//   have wildly different yields and the figure is about shape.
//
// WORKING POINT: use_tight_wp (default TRUE) selects quality&16 (Tight), false quality&8
//   (Medium). Per-muon, with the same caveat as the sibling macro: the nominal analysis
//   requires the WP on BOTH muons of the pair, but the single-muon tree carries no partner
//   information. That is the right population for a SINGLE-MUON fiducial cut.
//
// Usage:
//   root -l -b -q 'plot_muon_q_eta_pt_dependence.cxx+()'        // Tight
//   root -l -b -q 'plot_muon_q_eta_pt_dependence.cxx+(false)'   // Medium

#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

// Tight lands in new_gap_cuts/, Medium in new_gap_cuts/medium/ -- the same split the four
// sibling figures use, so the two working points never share a directory.
const char* kOutDirBase =
    "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/muon_gap_cuts/"
    "new_gap_cuts";
const char* kDataBase = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";

// Set once from the canonical q*eta binning; never retyped.
double kAxisLo = 0., kAxisHi = 0.;

// The adopted fiducial windows, READ from ParamsSet::single_mu_fiducial_gap_cuts.
std::vector<std::pair<double, double>> g_windows;

void BuildWindows() {
    g_windows.clear();
    for (const auto& g : ParamsSet::single_mu_fiducial_gap_cuts)
        g_windows.emplace_back(g.first, g.second);
    std::sort(g_windows.begin(), g_windows.end());
}

// p_T-INDEPENDENT by design (see the header note).
bool PassFiducialCut(double eta, int charge) {
    const double x = charge * eta;
    for (const auto& w : g_windows)
        if (x > w.first && x < w.second) return false;
    return true;
}

// ---- the plot-local p_T slices (user-specified; see header) ----
// The [4.0, 4.5) slice was DROPPED on 2026-09-09 (user decision). The muon p_T threshold moved
// 4 -> 4.5 GeV and it is an NTUPLE-STAGE cut, so the single-muon tree now contains nothing below
// 4.5 GeV: that slice would have been empty by construction while still drawing a legend entry.
// The figure's question -- the q*eta shape just above the threshold -- is unchanged; the
// threshold is simply 4.5 now.
const std::vector<double> kPtEdges = {4.5, 5.0, 6.0,
                                      std::numeric_limits<double>::infinity()};
const int kNPt = static_cast<int>(kPtEdges.size()) - 1;

std::string PtLabel(int i) {
    if (i == kNPt - 1) return Form("p_{T} > %.0f GeV", kPtEdges[i]);
    return Form("%.1f < p_{T} < %.1f GeV", kPtEdges[i], kPtEdges[i + 1]);
}

// Index of the p_T slice a muon falls in, or -1 if below threshold.
int PtSlice(double pt) {
    if (pt < kPtEdges.front()) return -1;  // cannot happen (tree cut), guarded anyway
    for (int i = kNPt - 1; i >= 0; --i)
        if (pt >= kPtEdges[i]) return i;
    return -1;
}

struct Sample {
    std::string dir, tag, trig;
    int nparts;
};

void AddParts(TChain& ch, const Sample& s) {
    for (int p = 1; p <= s.nparts; ++p)
        ch.Add((std::string(kDataBase) + s.dir + "single_muon_trees_" + s.tag + "_part" +
                std::to_string(p) + "_" + s.trig + "_mindR_0_02.root")
                   .c_str());
}

// Fill h[panel][ptslice] from one sample. pp fills panel 0 only; a PbPb sample fills the
// centrality panel its event belongs to. Muons with ev_centrality < 0 (the GetCentrality
// >85% sentinel) or >= 80 are dropped from the histograms AND from the counters, so every
// printed number describes exactly the plotted population.
void FillSample(const Sample& s, bool use_tight_wp, bool is_pbpb, int pp_panel,
                const std::vector<int>& ctr_lo, const std::vector<int>& ctr_hi,
                std::vector<std::vector<TH1D*>>& h,
                std::vector<std::vector<Long64_t>>& n_all,
                std::vector<std::vector<Long64_t>>& n_cut) {
    TChain ch("muon_tree");
    AddParts(ch, s);

    // SetMakeClass(1): the MuonObj branch is split, so leaves are addressable by name
    // without a compiled dictionary for MuonPbPb / MuonPP. The branch TYPES below must
    // match the tree exactly (Float_t eta/pt, Int_t charge/quality/ev_centrality) --
    // a mismatch leaves the variable untouched at its initial value and silently
    // produces an all-zero spectrum.
    ch.SetMakeClass(1);
    Float_t eta = 0.f, pt = 0.f;
    Int_t charge = 0, quality = 0, ev_centrality = -1;
    ch.SetBranchStatus("*", 0);
    for (const char* b : {"eta", "pt", "charge", "quality"}) ch.SetBranchStatus(b, 1);
    ch.SetBranchAddress("eta", &eta);
    ch.SetBranchAddress("pt", &pt);
    ch.SetBranchAddress("charge", &charge);
    ch.SetBranchAddress("quality", &quality);
    if (is_pbpb) {
        ch.SetBranchStatus("ev_centrality", 1);
        ch.SetBranchAddress("ev_centrality", &ev_centrality);
    }

    const int wp_bit = use_tight_wp ? 16 : 8;
    const Long64_t n_tree = ch.GetEntries();
    Long64_t n_used = 0;

    for (Long64_t i = 0; i < n_tree; ++i) {
        ch.GetEntry(i);
        if ((quality & wp_bit) == 0) continue;
        if (is_pbpb && (ev_centrality < 0 || ev_centrality >= ctr_hi.back())) continue;

        const int ip = PtSlice(pt);
        if (ip < 0) continue;

        int panel = pp_panel;
        if (is_pbpb) {
            panel = -1;
            for (size_t k = 0; k < ctr_lo.size(); ++k)
                if (ev_centrality >= ctr_lo[k] && ev_centrality < ctr_hi[k])
                    panel = static_cast<int>(k) + 1;  // panel 0 is pp
            if (panel < 0) continue;
        }

        h[panel][ip]->Fill(charge * eta);
        ++n_all[panel][ip];
        if (!PassFiducialCut(eta, charge)) ++n_cut[panel][ip];
        ++n_used;
    }
    printf("  %-12s : %10lld muons in tree, %10lld in the plotted population (%s%s)\n",
           s.tag.c_str(), n_tree, n_used, use_tight_wp ? "Tight" : "Medium",
           is_pbpb ? ", 0-80%" : "");
}

// ------------------------------ drawing helpers ------------------------------

void StyleShape(TH1D* h, int color) {
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetStats(0);
    h->GetXaxis()->SetTitle("q #times #eta^{#mu}");
    h->GetYaxis()->SetTitle("(1/N) dN_{#mu} / d(q#times#eta)");
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.70);
}

// Common LINEAR-y range over a GROUP of histograms. The grouping matters and is deliberate
// (see the call site): the five Pb+Pb panels share one range so the centrality dependence
// can be read across them, while pp gets its own -- a single range over all six would be set
// by pp's peak and would squeeze the Pb+Pb panels toward the bottom of their frames.
// Floor pinned at 0 -- a linear axis must start at zero or the apparent gap depth is
// visually distorted -- and ceiling at kHeadroom above the group peak. Unlike a log floor,
// nothing is ever pushed off scale.
const double kHeadroom = 1.15;
void SetCommonLinRange(const std::vector<TH1D*>& hs) {
    double hi = 0.;
    for (TH1D* h : hs)
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            hi = std::max(hi, h->GetBinContent(i));
    for (TH1D* h : hs) {
        h->SetMinimum(0.);
        h->SetMaximum(hi * kHeadroom);
    }
}

// Shaded bands for the adopted fiducial windows, drawn under the curves, which are then
// redrawn on top.
void DrawGapBands(const std::vector<TH1D*>& redraw) {
    const double y0 = gPad->GetUymin();
    const double y1 = gPad->GetUymax();
    for (const auto& w : g_windows) {
        TBox* b = new TBox(w.first, y0, w.second, y1);
        b->SetFillColorAlpha(kGreen - 7, 0.45);
        b->SetLineColor(kGreen - 7);
        b->Draw("same");
    }
    for (TH1D* h : redraw) h->Draw("hist same");
    gPad->RedrawAxis();
}

void Header(const std::string& text, const std::string& sub) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    const double frame_top = 1.0 - gPad->GetTopMargin();
    t.SetTextSize(0.042);
    t.DrawLatex(gPad->GetLeftMargin(), frame_top + (sub.empty() ? 0.025 : 0.067),
                text.c_str());
    if (!sub.empty()) {
        t.SetTextSize(0.035);
        t.DrawLatex(gPad->GetLeftMargin(), frame_top + 0.020, sub.c_str());
    }
}

// Both keys are drawn ONCE across the top of the canvas rather than inside the frames.
//
// This is not a style choice, it is the fix for a measured defect (2026-09-03 review; the
// same defect F13 fixed in the sibling macro). An in-frame key must be opaque, because it
// sits over the shaded gap bands -- so wherever it lands it ERASES what passes underneath.
// A 4-entry box needs ~0.29 of the frame height while the shared Pb+Pb axis leaves only
// ~0.15 above the peak, so there is no free corner for auto-placement to find: the box at
// bottom left covered the eta ~ 0 crack minimum in all six panels, making the crack read
// as ~0.44 of plateau where it is really ~0.25 -- i.e. it hid the very feature the middle
// fiducial window is drawn around. The keys are identical in all six panels, so the canvas
// band costs nothing and is the only placement that cannot cover data.
void CanvasKey(TCanvas& c, double y1, double y2, double text_size, int ncol,
               const std::vector<std::pair<TObject*, std::string>>& entries,
               const char* opt) {
    c.cd();
    TLegend* l = new TLegend(0.06, y1, 0.98, y2);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextSize(text_size);
    l->SetNColumns(ncol);
    // The default margin reserves a QUARTER of each column for the swatch; across a
    // full-canvas box that is a slab wider than the text it labels.
    l->SetMargin(0.06);
    for (const auto& e : entries) l->AddEntry(e.first, e.second.c_str(), opt);
    l->Draw();
}

// The key for the green bands. Values are FORMATTED from ParamsSet, never retyped, so the
// text cannot drift from the bands.
void CanvasGapKey(TCanvas& c, double y1, double y2, double text_size) {
    c.cd();
    TLegend* l = new TLegend(0.06, y1, 0.98, y2);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextSize(text_size);
    // The default margin reserves a QUARTER of the box width for the colour swatch; across
    // a full-canvas box that is a green slab wider than the text it labels.
    l->SetMargin(0.03);
    TBox* b = new TBox();
    b->SetFillColorAlpha(kGreen - 7, 0.45);
    std::string txt = "rejected at all p_{T}, q#times#eta in:";
    for (size_t i = 0; i < g_windows.size(); ++i)
        txt += Form("%s (%.2f, %.2f)", i ? "," : "", g_windows[i].first,
                    g_windows[i].second);
    l->AddEntry(b, txt.c_str(), "f");
    l->Draw();
}

TPad* MakePad(TCanvas& c, const std::string& name, double x1, double y1, double x2,
              double y2) {
    c.cd();  // pads belong to the CANVAS, not to whatever sub-pad is current
    TPad* p = new TPad(name.c_str(), "", x1, y1, x2, y2);
    p->SetLeftMargin(0.20);
    p->SetRightMargin(0.04);
    p->SetTopMargin(0.14);
    p->SetBottomMargin(0.14);
    p->SetTickx(1);
    p->SetTicky(1);
    p->Draw();
    return p;
}

}  // namespace

void plot_muon_q_eta_pt_dependence(bool use_tight_wp = true) {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    const std::string wp_tag = use_tight_wp ? "" : "_medium_wp";
    const std::string wp_lbl = use_tight_wp ? "Tight" : "Medium";

    BuildWindows();
    const std::string outdir =
        std::string(kOutDirBase) + (use_tight_wp ? "" : "/medium");
    gSystem->mkdir(outdir.c_str(), kTRUE);

    printf("fiducial windows from ParamsSet::single_mu_fiducial_gap_cuts =");
    for (const auto& w : g_windows) printf(" (%.2f, %.2f)", w.first, w.second);
    printf("\n");

    // ---- canonical binnings, read from their single sources of truth ----
    const std::vector<double> qeta_bins = ParamsSet::makeEtaTrigEffcyBinning(1);
    const int nb = static_cast<int>(qeta_bins.size()) - 1;
    kAxisLo = qeta_bins.front();
    kAxisHi = qeta_bins.back();
    printf("q*eta binning: %d bins over [%.2f, %.2f] "
           "(ParamsSet::makeEtaTrigEffcyBinning(1))\n",
           nb, kAxisLo, kAxisHi);

    // ParamsSet::ctrbins = {0,5,10,20,30,50,80}; merge the FIRST TWO bins into 0-10,
    // identically to plot_muon_q_eta_spectrum.cxx.
    std::vector<int> ctr_lo, ctr_hi;
    for (size_t i = 1; i + 1 < ParamsSet::ctrbins.size(); ++i) {
        ctr_lo.push_back(i == 1 ? ParamsSet::ctrbins.front() : ParamsSet::ctrbins[i]);
        ctr_hi.push_back(ParamsSet::ctrbins[i + 1]);
    }
    const int nctr = static_cast<int>(ctr_lo.size());
    const int npanel = nctr + 1;  // panel 0 = pp24
    printf("centrality bins (ParamsSet::ctrbins, 0-5 + 5-10 merged):");
    for (int k = 0; k < nctr; ++k) printf(" %d-%d%%", ctr_lo[k], ctr_hi[k]);
    printf("\np_T slices (plot-local, user-specified):");
    for (int i = 0; i < kNPt; ++i)
        printf(" [%.1f,%s)", kPtEdges[i],
               i == kNPt - 1 ? "inf" : Form("%.1f", kPtEdges[i + 1]));
    printf("\n");

    std::vector<std::vector<TH1D*>> h(npanel, std::vector<TH1D*>(kNPt, nullptr));
    std::vector<std::vector<Long64_t>> n_all(npanel, std::vector<Long64_t>(kNPt, 0));
    std::vector<std::vector<Long64_t>> n_cut(npanel, std::vector<Long64_t>(kNPt, 0));
    for (int p = 0; p < npanel; ++p)
        for (int i = 0; i < kNPt; ++i) {
            h[p][i] = new TH1D(Form("h_p%d_pt%d", p, i), "", nb, qeta_bins.data());
            h[p][i]->Sumw2();
        }

    printf("\n=== pp24 ===\n");
    FillSample({"pp_2024/", "pp_2024", "2mu4", 12}, use_tight_wp, false, 0, ctr_lo, ctr_hi,
               h, n_all, n_cut);

    printf("\n=== PbPb 2023 + 2024 + 2025 (combined) ===\n");
    const std::vector<Sample> pbpb = {{"pbpb_2023/", "pbpb_2023", "single_mu4", 4},
                                      {"pbpb_2024/", "pbpb_2024", "single_mu4", 2},
                                      {"pbpb_2025/", "pbpb_2025", "single_mu4", 6}};
    for (const auto& s : pbpb)
        FillSample(s, use_tight_wp, true, -1, ctr_lo, ctr_hi, h, n_all, n_cut);

    // ---- unit-area PDFs on a variable-width axis ----
    // Scale(1.,"width") first makes the content a density; Integral("width") of that is
    // the original entry count, so the second scale normalises the area to 1.
    for (int p = 0; p < npanel; ++p)
        for (int i = 0; i < kNPt; ++i) {
            h[p][i]->Scale(1., "width");
            const double integ = h[p][i]->Integral("width");
            if (integ > 0) h[p][i]->Scale(1. / integ);
        }

    // ---- panel headers ----
    std::vector<std::string> head(npanel), sub(npanel);
    // The trigger belongs on EVERY panel, not just pp: the two samples are taken with
    // different triggers (pp 2mu4, both legs; Pb+Pb single mu4) and the low-p_T slices are
    // shaped by that turn-on, which is deliberately not corrected here. Naming it on pp
    // alone invites the reader to assume one trigger across the canvas. It goes on the
    // sub-line so neither line overflows a half-width pad.
    head[0] = "pp #sqrt{s} = 5.36 TeV, 2024";
    sub[0] = "HLT_2mu4, " + wp_lbl + " WP (per muon)";
    for (int k = 0; k < nctr; ++k) {
        head[k + 1] = Form("Pb+Pb #sqrt{s_{NN}} = 5.36 TeV, 23+24+25, %d-%d%%", ctr_lo[k],
                           ctr_hi[k]);
        sub[k + 1] = "HLT_mu4, " + wp_lbl + " WP (per muon)";
    }

    // ---- draw: 3 rows x 2 cols (nrows >= ncols for N = 6) ----
    // Colours avoid green, which is the gap-band colour.
    const int cols[4] = {kBlack, kRed + 1, kAzure + 2, kMagenta + 1};
    for (int p = 0; p < npanel; ++p)
        for (int i = 0; i < kNPt; ++i) StyleShape(h[p][i], cols[i]);
    // Two range groups (see SetCommonLinRange): pp alone, and the five Pb+Pb panels
    // together so their centrality dependence is readable on one axis.
    SetCommonLinRange(h[0]);
    std::vector<TH1D*> pb_h;
    for (int p = 1; p < npanel; ++p)
        pb_h.insert(pb_h.end(), h[p].begin(), h[p].end());
    SetCommonLinRange(pb_h);

    TCanvas c("c_qeta_pt", "", 1400, 1760);
    const int ncol = 2, nrow = 3;
    const double kPanelTop = 0.942;  // band above the panels for the two canvas-level keys
    for (int p = 0; p < npanel; ++p) {
        const int col = p % ncol, row = p / ncol;
        MakePad(c, Form("pad%d", p), col / double(ncol),
                kPanelTop * (1.0 - (row + 1) / double(nrow)), (col + 1) / double(ncol),
                kPanelTop * (1.0 - row / double(nrow)))
            ->cd();
        for (int i = 0; i < kNPt; ++i) h[p][i]->Draw(i == 0 ? "hist" : "hist same");
        gPad->Update();
        DrawGapBands(h[p]);
        Header(head[p], sub[p]);
        // NO in-frame legend -- see the comment on CanvasKey. The p_T key is drawn once
        // across the top of the canvas below.
    }

    // p_T key first (nearest the panels, four columns in the muon colours), gap key above.
    std::vector<std::pair<TObject*, std::string>> pt_entries;
    for (int i = 0; i < kNPt; ++i) pt_entries.emplace_back(h[0][i], PtLabel(i));
    CanvasKey(c, kPanelTop + 0.003, kPanelTop + 0.027, 0.0130, kNPt, pt_entries, "l");
    CanvasGapKey(c, kPanelTop + 0.031, kPanelTop + 0.055, 0.0130);
    const std::string out = outdir + "/muon_q_eta_pt_dependence" + wp_tag + ".png";
    c.SaveAs(out.c_str());

    // ---------------- numeric summary ----------------
    // The cost of the p_T-INDEPENDENT fiducial cut, resolved in p_T: if these columns
    // drift strongly with p_T the single-window-set choice is carrying a p_T dependence.
    printf("\n=== fraction of single muons rejected by the fiducial cut, by p_T slice ===\n");
    printf("%-22s", "panel");
    for (int i = 0; i < kNPt; ++i)
        printf("%14s", i == kNPt - 1 ? Form(">%.0f GeV", kPtEdges[i])
                                     : Form("%.1f-%.1f GeV", kPtEdges[i], kPtEdges[i + 1]));
    printf("%14s\n", "all pT");
    for (int p = 0; p < npanel; ++p) {
        const std::string name =
            p == 0 ? "pp24" : Form("PbPb %d-%d%%", ctr_lo[p - 1], ctr_hi[p - 1]);
        printf("%-22s", name.c_str());
        Long64_t tot = 0, tot_cut = 0;
        for (int i = 0; i < kNPt; ++i) {
            printf("%14s", n_all[p][i] > 0
                               ? Form("%.4f", double(n_cut[p][i]) / n_all[p][i])
                               : "-");
            tot += n_all[p][i];
            tot_cut += n_cut[p][i];
        }
        printf("%14s\n", tot > 0 ? Form("%.4f", double(tot_cut) / tot) : "-");
    }

    printf("\nPNG written to %s\n", out.c_str());
}
