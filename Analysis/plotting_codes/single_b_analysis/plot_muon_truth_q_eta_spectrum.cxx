// plot_muon_truth_q_eta_spectrum.cxx
//
// TRUTH single-muon q*eta spectrum, dsigma_mu/d(q*eta), with NO reco match, NO working point
// and NO trigger requirement -- the truth fiducial and nothing else. It is the reference
// the data spectra of plot_muon_q_eta_spectrum.cxx cannot provide: a q*eta distribution
// carrying neither the trigger nor the reconstruction efficiency, and the raw input to the
// gap-cut acceptance eps_acc.
// See Analysis/docs/tracking/muon_gap_cuts_acceptance.md (sub-step 6).
//
// INPUT (NTuple-processing OUTPUT -- never the raw NTUPs): TTree "muon_tree" written by
//   PythiaFullSimExtras.c:395-412 in its NOMINAL mode, whose gate is exactly
//       truth_pt > 4.0 && fabs(truth_eta) < 2.4
//   with NO reco match required (it is the single-muon reco-efficiency denominator).
//   The store_mc_trigger variant of the same tree is NOT interchangeable: its gate is a
//   RECO one (reco_match && pt > 3 && |eta| < 2.6), so it would reintroduce exactly the
//   reconstruction efficiency this figure exists to remove. The truth fiducial is
//   RE-APPLIED in the fill loop below and the number of entries failing it is reported --
//   a non-zero count means the wrong file was given.
//   The muon loop is seeded from the Pythia-only truth block (PythiaTruthExtras.h:68-88),
//   so HIJING muons never enter, in the overlay either.
//
//   pp24  : Pythia fullsim FULL "_pdf" production (DSIDs 803015-803020), pp beam only.
//   PbPb  : Pythia + HIJING overlay, PbPb23 conditions, TEST sample -- ~150x fewer muons
//           than pp24. Its statistics are stated on the figure.
//
// WEIGHT: every histogram is filled with ev_weight = ami_w * nom_ratio / N_beam
//   (PythiaAlgCoreT.c:810-812), i.e. the AMI sigma * eps_filt of the pT-hat slice times
//   the isospin beam ratio. MC is never plotted unweighted: the slice weights set the
//   pT-hat mixture, and getting them wrong is the silent, non-cancelling failure of
//   Analysis/docs/ami_weights.md. Because that weight is sigma * eps_filt / N_events in nb,
//   the plotted density is a CROSS SECTION dsigma_mu/d(q*eta) in nb, not a raw muon count;
//   the raw counts behind it are printed and quoted on the figure.
//
// WORKING POINT: deliberately ABSENT. A truth muon has no reconstruction quality, so the
//   Tight/Medium config var that every other plot set carries (docs/muon_wp_registry.md)
//   has no meaning here; adding one would suggest a choice that does not exist.
//
// BINNING: q*eta uses ParamsSet::makeEtaTrigEffcyBinning(1), the SAME canonical fine axis
//   as the data figures, so truth and data panels are cell-for-cell comparable. Not a new
//   binning, and nothing is retyped. Bins are variable-width => everything is a DENSITY
//   dN/d(q*eta) via Scale(1., "width").
//
// GAP WINDOWS: the settled fiducial definition, READ from
//   ParamsSet::single_mu_fiducial_gap_cuts -- re-tuning the cut means editing ParamsSet.h
//   alone and re-running. The truth-level fraction each window removes is the raw
//   acceptance cost eps_acc has to compensate, and is printed per window and in total.
//
// Usage:
//   root -l -b -q 'plot_muon_truth_q_eta_spectrum.cxx+()'

#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace {

const char* kOutDir =
    "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/muon_gap_cuts/"
    "truth_q_eta";

// The truth fiducial, restated from PythiaFullSimExtras.c:405 (the gate that filled the
// tree) and identical to FillMCTrigEffHists.cxx:603 kTruthFidSingle.
const double kTruthPtMin = 4.0;
const double kTruthEtaMax = 2.4;

// Split point of the third panel. Kept at the value used by the DATA figures so the two
// sets can be laid side by side; at truth level it carries no cut meaning (the fiducial
// windows below are pT-independent by construction).
const double kPtSplit = 6.0;

// Plot axis limits, set once from the canonical binning. Never retyped.
double kAxisLo = 0., kAxisHi = 0.;

// ---------------------------------------------------------------------------
// The fiducial gap cut, READ from ParamsSet::single_mu_fiducial_gap_cuts.
// pT-INDEPENDENT by design, so the acceptance factorises as a pure function of q*eta.
// ---------------------------------------------------------------------------
std::vector<std::pair<double, double>> g_windows;

void BuildWindows() {
    g_windows.clear();
    for (const auto& w : ParamsSet::single_mu_fiducial_gap_cuts)
        g_windows.emplace_back(w.first, w.second);
    std::sort(g_windows.begin(), g_windows.end());
}

// Index of the window a muon falls in, or -1 if it passes. Returning the index (rather
// than a bool) is what lets the per-window acceptance cost be counted in one pass.
int FailingWindow(double q_eta) {
    for (size_t i = 0; i < g_windows.size(); ++i)
        if (q_eta > g_windows[i].first && q_eta < g_windows[i].second) return int(i);
    return -1;
}

// ---------------------------------------------------------------------------

struct Sample {
    std::string file;
    std::string headline;  // drawn across the top of the canvas
    std::string fname;     // output PNG stem
};

// Weighted and raw yields, and what the fiducial windows remove from them. Both are kept
// because they answer different questions: the weighted fraction is the acceptance cost
// eps_acc must compensate, the raw one is the statistical precision behind it.
struct FillCounts {
    Long64_t n_raw = 0;
    double n_wgt = 0.;
    std::vector<Long64_t> cut_raw;
    std::vector<double> cut_wgt;
    Long64_t n_outside_fiducial = 0;  // provenance check: must be 0 on a nominal file
};

FillCounts FillSample(const Sample& s, TH1D* h_all, TH1D* h_pos, TH1D* h_neg, TH1D* h_lopt,
                      TH1D* h_hipt, std::string* centrality_note) {
    FillCounts fc;
    fc.cut_raw.assign(g_windows.size(), 0);
    fc.cut_wgt.assign(g_windows.size(), 0.);

    TFile* f = TFile::Open(s.file.c_str(), "READ");
    if (!f || f->IsZombie()) {
        printf("  ERROR: cannot open %s\n", s.file.c_str());
        return fc;
    }
    TTree* t = dynamic_cast<TTree*>(f->Get("muon_tree"));
    if (!t) {
        printf("  ERROR: no TTree 'muon_tree' in %s\n", s.file.c_str());
        f->Close();
        return fc;
    }

    // SetMakeClass(1): the MuonObj branch is split, so leaves are addressable by name
    // without a compiled dictionary for the MuonPythiaFullSim* classes.
    t->SetMakeClass(1);
    Float_t truth_pt = 0.f, truth_eta = 0.f, ev_weight = 0.f;
    Int_t truth_charge = 0, ev_centrality = -1;
    t->SetBranchStatus("*", 0);
    for (const char* b : {"truth_pt", "truth_eta", "truth_charge", "ev_weight"})
        t->SetBranchStatus(b, 1);
    t->SetBranchAddress("truth_pt", &truth_pt);
    t->SetBranchAddress("truth_eta", &truth_eta);
    t->SetBranchAddress("truth_charge", &truth_charge);
    t->SetBranchAddress("ev_weight", &ev_weight);
    // Overlay only. Not a selection -- a truth spectrum has no centrality dependence to
    // cut on -- but the centrality range the sample covers belongs on the figure.
    const bool has_ctr = t->GetBranch("ev_centrality") != nullptr;
    if (has_ctr) {
        t->SetBranchStatus("ev_centrality", 1);
        t->SetBranchAddress("ev_centrality", &ev_centrality);
    }
    int ctr_min = 1000, ctr_max = -1000;

    const Long64_t n_tree = t->GetEntries();
    for (Long64_t i = 0; i < n_tree; ++i) {
        t->GetEntry(i);
        // Re-applied, not assumed: see the INPUT note at the top.
        if (!(truth_pt > kTruthPtMin && std::fabs(truth_eta) < kTruthEtaMax)) {
            ++fc.n_outside_fiducial;
            continue;
        }
        if (has_ctr) {
            ctr_min = std::min(ctr_min, int(ev_centrality));
            ctr_max = std::max(ctr_max, int(ev_centrality));
        }

        const double q_eta = truth_charge * truth_eta;
        const double w = ev_weight;
        ++fc.n_raw;
        fc.n_wgt += w;
        const int iw = FailingWindow(q_eta);
        if (iw >= 0) {
            ++fc.cut_raw[iw];
            fc.cut_wgt[iw] += w;
        }

        h_all->Fill(q_eta, w);
        (truth_charge > 0 ? h_pos : h_neg)->Fill(q_eta, w);
        (truth_pt < kPtSplit ? h_lopt : h_hipt)->Fill(q_eta, w);
    }

    if (centrality_note) {
        // The MEASURED span of the sample, not an assumed range: the overlay test sample
        // is not guaranteed to cover 0-80%, and stating a range it does not cover would
        // misdescribe the population.
        *centrality_note = (has_ctr && ctr_max >= ctr_min)
                               ? std::string(Form("centrality %d-%d%%", ctr_min, ctr_max))
                               : std::string();
    }

    printf("  %-40s : %10lld entries, %10lld in the truth fiducial "
           "(sum of weights %.4g), %lld outside it\n",
           s.fname.c_str(), n_tree, fc.n_raw, fc.n_wgt, fc.n_outside_fiducial);
    if (fc.n_outside_fiducial > 0)
        printf("  WARNING: %lld entries fail truth_pt > %.1f && |truth_eta| < %.1f -- this "
               "file was not written with the nominal truth gate\n",
               fc.n_outside_fiducial, kTruthPtMin, kTruthEtaMax);

    f->Close();
    return fc;
}

// ------------------------------ drawing helpers ------------------------------

void StyleSpectrum(TH1D* h, int color) {
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetStats(0);
    h->GetXaxis()->SetTitle("q #times #eta^{#mu}_{truth}");
    h->GetYaxis()->SetTitle("d#sigma_{#mu} / d(q#times#eta)  [nb]");
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.70);
}

// Legends must be opaque: they sit over the shaded gap bands.
void StyleLegend(TLegend& l, double text_size) {
    l.SetBorderSize(0);
    l.SetFillColor(kWhite);
    l.SetFillStyle(1001);
    l.SetTextSize(text_size);
}

// Height a legend needs so its entries do not print on top of each other: ROOT divides the
// box evenly over the entries, so a box sized independently of the entry count overlaps its
// own text once the count grows.
double LegendHeight(int n_entries, double text_size) {
    return n_entries * 1.8 * text_size + 0.02;
}

// Where to put a legend so it hides no curve. These legends are opaque (they sit over the
// gap bands) and drawn last, so a fixed corner erases whatever passes underneath it -- and
// the corner that is empty for one sample is crossed by the next one's acceptance edge.
// Scores the four in-frame corners by how much of the drawn curve falls inside the box and
// returns the emptiest. MUST be called after gPad->Update(): it reads the pad user range.
// Scoring uses consecutive-bin SEGMENTS, not bin contents, because a step histogram draws
// the vertical connector between bins and a steep plunge can cross a box with no single bin
// content inside it. Kept identical to plot_muon_q_eta_spectrum.cxx so the truth and data
// figures place their legends by the same rule.
struct LegBox {
    double x1, y1, x2, y2;
};

LegBox AutoLegendBox(const std::vector<TH1D*>& hs, double w, double h) {
    const double lm = gPad->GetLeftMargin(), rm = gPad->GetRightMargin();
    const double bm = gPad->GetBottomMargin(), tm = gPad->GetTopMargin();
    const double inset = 0.02;
    const double ux0 = gPad->GetUxmin(), ux1 = gPad->GetUxmax();
    const double uy0 = gPad->GetUymin(), uy1 = gPad->GetUymax();  // log10 units if logy
    const bool logy = gPad->GetLogy();

    // Conventional order, so a tie keeps the familiar top-right placement.
    const std::vector<std::pair<double, double>> cands = {
        {1 - rm - inset - w, 1 - tm - inset - h},  // top right
        {lm + inset, 1 - tm - inset - h},          // top left
        {1 - rm - inset - w, bm + inset},          // bottom right
        {lm + inset, bm + inset}};                 // bottom left

    auto to_x = [&](double xn) { return ux0 + (xn - lm) / (1 - lm - rm) * (ux1 - ux0); };
    auto to_y = [&](double yn) {
        const double v = uy0 + (yn - bm) / (1 - bm - tm) * (uy1 - uy0);
        return logy ? std::pow(10., v) : v;
    };

    size_t best = 0;
    int best_score = -1;
    for (size_t i = 0; i < cands.size(); ++i) {
        const double bx0 = to_x(cands[i].first), bx1 = to_x(cands[i].first + w);
        const double by0 = to_y(cands[i].second), by1 = to_y(cands[i].second + h);
        int hits = 0;
        for (const TH1D* hh : hs)
            for (int b = 1; b < hh->GetNbinsX(); ++b) {
                const double xa = hh->GetBinCenter(b), xb = hh->GetBinCenter(b + 1);
                if (xb < bx0 || xa > bx1) continue;
                const double ca = hh->GetBinContent(b), cb = hh->GetBinContent(b + 1);
                if (std::max(ca, cb) >= by0 && std::min(ca, cb) <= by1) ++hits;
            }
        if (best_score < 0 || hits < best_score) {
            best_score = hits;
            best = i;
        }
    }
    return {cands[best].first, cands[best].second, cands[best].first + w,
            cands[best].second + h};
}

// Shaded bands for the fiducial gap windows, drawn under the curves, which are then
// redrawn on top.
void DrawGapMarkers(const std::vector<TH1D*>& redraw, const char* opt = "hist e same") {
    const double ymin = gPad->GetUymin(), ymax = gPad->GetUymax();
    const double y0 = gPad->GetLogy() ? std::pow(10., ymin) : ymin;
    const double y1 = gPad->GetLogy() ? std::pow(10., ymax) : ymax;
    for (const auto& w : g_windows) {
        TBox* b = new TBox(w.first, y0, w.second, y1);
        b->SetFillColorAlpha(kGreen - 7, 0.45);
        b->SetLineColor(kGreen - 7);
        b->Draw("same");
    }
    for (TH1D* h : redraw) h->Draw(opt);
    gPad->RedrawAxis();
}

// The key for the gap-band overlay. Values are FORMATTED from the constants, never
// retyped, so the text cannot drift from the bands.
void AddGapKeyEntries(TLegend& l, bool compact = false) {
    TBox* b = new TBox();
    b->SetFillColorAlpha(kGreen - 7, 0.45);
    if (compact) {
        std::string txt = "rejected by the fiducial gap cut, q#times#eta in:";
        for (size_t i = 0; i < g_windows.size(); ++i)
            txt += Form("%s (%.2f, %.2f)", i ? "," : "", g_windows[i].first,
                        g_windows[i].second);
        l.AddEntry(b, txt.c_str(), "f");
    } else {
        l.AddEntry(b, "rejected by the fiducial gap cut:", "f");
        for (const auto& w : g_windows)
            l.AddEntry((TObject*)nullptr,
                       Form("   %.2f < q#times#eta < %.2f", w.first, w.second), "");
    }
}

// The gap key drawn ONCE across the top of a canvas. In a 1/3-width panel an in-frame box
// carrying the three window lines has nowhere to sit that is not on top of the spectrum, and
// the overlay is the same in every panel anyway.
void CanvasGapKey(TCanvas& c, double y1, double y2, double text_size) {
    c.cd();
    TLegend* l = new TLegend(0.06, y1, 0.98, y2);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextSize(text_size);
    // Default margin reserves a QUARTER of the box width for the colour swatch; across a
    // full-canvas box that is a green slab wider than the text it labels.
    l->SetMargin(0.03);
    AddGapKeyEntries(*l, /*compact=*/true);
    l->Draw();
}

// Common LINEAR-y range for a set of overlaid histograms. Floor pinned at 0 -- a linear
// axis must start at zero or the apparent depth of the gap structure this figure exists
// to show is visually distorted -- and ceiling at kHeadroom above the peak so the curve
// never touches the frame. Unlike a log floor, nothing is ever pushed off scale.
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

void Header(const std::string& text) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.042);
    t.DrawLatex(gPad->GetLeftMargin(), 1.0 - gPad->GetTopMargin() + 0.025, text.c_str());
}

// Sample headline drawn ONCE across the top of the canvas: a per-panel pad is 1/3 of the
// canvas wide and would clip it, and the sample is common to all panels anyway.
void CanvasHeadline(TCanvas& c, const std::string& text, const std::string& sub,
                    double y_main = 0.955, double y_sub = 0.915, double size_main = 0.038,
                    double size_sub = 0.032) {
    c.cd();
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(size_main);
    t.DrawLatex(0.065, y_main, text.c_str());
    if (!sub.empty()) {
        t.SetTextSize(size_sub);
        t.DrawLatex(0.065, y_sub, sub.c_str());
    }
}

// A main+ratio pad pair inside the rectangle [x1,x2] x [y1,y2] of the current canvas.
std::pair<TPad*, TPad*> MakeRatioPads(TCanvas& c, const std::string& name, double x1,
                                      double y1, double x2, double y2) {
    c.cd();  // pads belong to the CANVAS, not to whatever sub-pad is current
    const double ysplit = y1 + 0.34 * (y2 - y1);
    TPad* p_main = new TPad((name + "_main").c_str(), "", x1, ysplit, x2, y2);
    TPad* p_rat = new TPad((name + "_rat").c_str(), "", x1, y1, x2, ysplit);
    p_main->SetLeftMargin(0.20);
    p_main->SetRightMargin(0.04);
    p_main->SetTopMargin(0.14);
    p_main->SetBottomMargin(0.02);
    p_main->SetTickx(1);
    p_main->SetTicky(1);
    p_rat->SetLeftMargin(0.20);
    p_rat->SetRightMargin(0.04);
    p_rat->SetTopMargin(0.03);
    p_rat->SetBottomMargin(0.34);
    p_rat->SetTickx(1);
    p_rat->SetTicky(1);
    p_main->Draw();
    p_rat->Draw();
    return {p_main, p_rat};
}

TPad* MakeSinglePad(TCanvas& c, const std::string& name, double x1, double y1, double x2,
                    double y2) {
    c.cd();
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

// Build num/den. The two histograms are always DISJOINT samples here (mu+ vs mu-, low-pT
// vs high-pT, pp vs Pb+Pb), so the default independent error propagation of TH1::Divide is
// correct -- unlike a subset/superset ratio, which needs the conditional (binomial) form.
TH1D* MakeRatio(const TH1D* num, const TH1D* den, const std::string& name, int color,
                const std::string& ytitle) {
    TH1D* r = (TH1D*)num->Clone(name.c_str());
    r->Divide(den);
    r->SetStats(0);
    r->SetLineColor(color);
    r->SetMarkerColor(color);
    r->SetLineWidth(2);
    r->GetYaxis()->SetTitle(ytitle.c_str());
    r->GetXaxis()->SetTitle("q #times #eta^{#mu}_{truth}");
    // The ratio pad is ~1/3 the height of the main pad, so text/ticks need scaling up.
    r->GetXaxis()->SetTitleSize(0.115);
    r->GetXaxis()->SetLabelSize(0.100);
    r->GetXaxis()->SetTitleOffset(1.20);
    r->GetYaxis()->SetTitleSize(0.095);
    r->GetYaxis()->SetLabelSize(0.090);
    r->GetYaxis()->SetTitleOffset(0.80);
    r->GetYaxis()->SetNdivisions(505);
    return r;
}

// Y-range covering every plotted ratio point INCLUDING its error bar, so nothing is
// cropped, and always containing the unity reference the dashed line marks.
void SetRatioRange(const std::vector<TH1D*>& rs) {
    double lo = 1e300, hi = -1e300;
    for (TH1D* r : rs)
        for (int i = 1; i <= r->GetNbinsX(); ++i) {
            if (r->GetBinContent(i) == 0 && r->GetBinError(i) == 0) continue;
            lo = std::min(lo, r->GetBinContent(i) - r->GetBinError(i));
            hi = std::max(hi, r->GetBinContent(i) + r->GetBinError(i));
        }
    if (lo > hi) return;
    lo = std::min(lo, 1.0);
    hi = std::max(hi, 1.0);
    const double pad = 0.08 * (hi - lo) + 1e-6;
    for (TH1D* r : rs) {
        r->SetMinimum(lo - pad);
        r->SetMaximum(hi + pad);
    }
}

void DrawRatioUnityLine() {
    TLine* l = new TLine(kAxisLo, 1., kAxisHi, 1.);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw("same");
}

// Local minima of the truth spectrum. At truth level there is no detector, so any dip here
// is generator kinematics -- reported precisely because it must NOT be mistaken for the
// detector structure the data spectra show at the same q*eta.
void ReportDips(const TH1D* h, const char* label) {
    printf("\n--- %s : local minima of the TRUTH dsigma/d(q*eta) (bins below 60%% of the "
           "local 21-bin median) ---\n",
           label);
    const int n = h->GetNbinsX();
    bool any = false;
    for (int i = 1; i <= n; ++i) {
        std::vector<double> loc;
        for (int j = std::max(1, i - 10); j <= std::min(n, i + 10); ++j)
            loc.push_back(h->GetBinContent(j));
        std::sort(loc.begin(), loc.end());
        const double med = loc[loc.size() / 2];
        const double c = h->GetBinContent(i);
        if (med > 0 && c < 0.60 * med) {
            any = true;
            printf("   q*eta in [%+.3f, %+.3f]  density = %10.4g   local median = "
                   "%10.4g   ratio = %.2f\n",
                   h->GetBinLowEdge(i), h->GetBinLowEdge(i + 1), c, med, c / med);
        }
    }
    if (!any) printf("   none -- the truth spectrum is smooth over the whole axis\n");
}

void ReportAcceptance(const char* label, const FillCounts& fc) {
    printf("\n%s -- fraction of TRUTH muons removed by the fiducial gap windows\n", label);
    printf("   (weighted = the acceptance cost eps_acc must compensate; "
           "raw = the statistics behind it)\n");
    double sum_wgt = 0.;
    Long64_t sum_raw = 0;
    for (size_t i = 0; i < g_windows.size(); ++i) {
        sum_wgt += fc.cut_wgt[i];
        sum_raw += fc.cut_raw[i];
        printf("   window [%+.2f, %+.2f] : weighted %.4f   raw %.4f  (%lld / %lld)\n",
               g_windows[i].first, g_windows[i].second,
               fc.n_wgt > 0 ? fc.cut_wgt[i] / fc.n_wgt : 0.,
               fc.n_raw > 0 ? double(fc.cut_raw[i]) / fc.n_raw : 0., fc.cut_raw[i],
               fc.n_raw);
    }
    printf("   ALL WINDOWS         : weighted %.4f   raw %.4f  (%lld / %lld)\n",
           fc.n_wgt > 0 ? sum_wgt / fc.n_wgt : 0.,
           fc.n_raw > 0 ? double(sum_raw) / fc.n_raw : 0., sum_raw, fc.n_raw);
    printf("   => single-muon acceptance of the fiducial cut, eps_acc = %.4f (weighted)\n",
           fc.n_wgt > 0 ? 1. - sum_wgt / fc.n_wgt : 0.);
}

}  // namespace

void plot_muon_truth_q_eta_spectrum() {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    BuildWindows();
    gSystem->mkdir(kOutDir, kTRUE);

    printf("fiducial gap windows from ParamsSet::single_mu_fiducial_gap_cuts =");
    for (const auto& w : g_windows) printf(" (%.2f, %.2f)", w.first, w.second);
    printf("\ntruth fiducial: truth_pt > %.1f GeV && |truth_eta| < %.1f "
           "(no reco match, no working point, no trigger)\n",
           kTruthPtMin, kTruthEtaMax);

    // ---- canonical binning, read from its single source of truth ----
    const std::vector<double> qeta_bins = ParamsSet::makeEtaTrigEffcyBinning(1);
    const int nb = static_cast<int>(qeta_bins.size()) - 1;
    kAxisLo = qeta_bins.front();  // never retyped
    kAxisHi = qeta_bins.back();
    printf("q*eta binning: %d bins over [%.2f, %.2f] "
           "(ParamsSet::makeEtaTrigEffcyBinning(1))\n",
           nb, qeta_bins.front(), qeta_bins.back());

    auto book = [&](const char* name) {
        TH1D* h = new TH1D(name, "", nb, qeta_bins.data());
        h->Sumw2();
        return h;
    };

    const std::string kMCBase = "/usatlas/u/yuhanguo/usatlasdata/";
    const std::vector<Sample> samples = {
        {kMCBase + "pythia_fullsim_full_sample/"
                   "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_"
                   "single_muon_full.root",
         "Pythia8 hard-QCD dimuon, pp #sqrt{s} = 5.36 TeV, full simulation",
         "muon_truth_q_eta_spectrum_pp24_fullsim"},
        {kMCBase + "pythia_fullsim_hijing_overlay_test_sample/"
                   "muon_pairs_pythia_fullsim_hijing_overlay_pbpb23_no_data_resonance_"
                   "cuts_single_muon.root",
         "Pythia8 hard-QCD dimuon + HIJING overlay, Pb+Pb #sqrt{s_{NN}} = 5.36 TeV",
         "muon_truth_q_eta_spectrum_pbpb_overlay"}};

    struct Filled {
        TH1D *all, *pos, *neg, *lopt, *hipt;
        FillCounts fc;
        std::string sub, headline, fname;
    };
    std::vector<Filled> res;

    for (size_t i = 0; i < samples.size(); ++i) {
        printf("\n=== %s ===\n", samples[i].fname.c_str());
        Filled r;
        r.all = book(Form("h_all_%zu", i));
        r.pos = book(Form("h_pos_%zu", i));
        r.neg = book(Form("h_neg_%zu", i));
        r.lopt = book(Form("h_lopt_%zu", i));
        r.hipt = book(Form("h_hipt_%zu", i));
        std::string ctr;
        r.fc = FillSample(samples[i], r.all, r.pos, r.neg, r.lopt, r.hipt, &ctr);
        r.headline = samples[i].headline;
        r.fname = samples[i].fname;
        r.sub = Form("truth muons, p_{T} > %.0f GeV, |#eta| < %.1f; "
                     "no reco match / working point / trigger%s%s  (%lld muons)",
                     kTruthPtMin, kTruthEtaMax, ctr.empty() ? "" : "; ",
                     ctr.c_str(), r.fc.n_raw);
        for (TH1D* h : {r.all, r.pos, r.neg, r.lopt, r.hipt}) h->Scale(1., "width");
        res.push_back(r);
    }

    // ================== Plots A / B : one canvas per sample ==================
    // Three columns (N <= 3 => one row): the spectrum with the fiducial windows overlaid;
    // the same split by truth charge; and split at truth pT = 6 GeV. Columns 2 and 3 carry
    // a ratio pad, because the ratio is exactly what those comparisons are about.
    for (const auto& ps : res) {
        if (ps.fc.n_raw == 0) continue;
        TCanvas c(("c_" + ps.fname).c_str(), "", 2100, 840);

        // --- column 1: the spectrum + fiducial windows ---
        MakeSinglePad(c, ps.fname + "_p1", 0.000, 0.0, 0.334, 0.875)->cd();
        StyleSpectrum(ps.all, kBlack);
        SetCommonLinRange({ps.all});
        ps.all->Draw("hist e");
        gPad->Update();
        DrawGapMarkers({ps.all});
        Header("all truth muons; fiducial gap cut overlaid");
        // No in-frame legend: one curve, named by the header, and the window key is drawn
        // once across the canvas below. A box carrying the key would span the frame.

        // --- column 2: by truth charge, with the mu+/mu- ratio ---
        {
            auto pads = MakeRatioPads(c, ps.fname + "_p2", 0.334, 0.0, 0.667, 0.875);
            pads.first->cd();
            StyleSpectrum(ps.pos, kRed + 1);
            StyleSpectrum(ps.neg, kAzure + 2);
            for (TH1D* h : {ps.pos, ps.neg}) {
                h->GetXaxis()->SetLabelSize(0.);
                h->GetXaxis()->SetTitleSize(0.);
            }
            SetCommonLinRange({ps.pos, ps.neg});
            ps.pos->Draw("hist e");
            ps.neg->Draw("hist e same");
            gPad->Update();
            DrawGapMarkers({ps.pos, ps.neg});
            Header("by truth charge");
            {
                const LegBox lb = AutoLegendBox({ps.pos, ps.neg}, 0.23,
                                                LegendHeight(2, 0.050));
                TLegend l(lb.x1, lb.y1, lb.x2, lb.y2);
                StyleLegend(l, 0.050);
                l.AddEntry(ps.pos, "#mu^{+}", "l");
                l.AddEntry(ps.neg, "#mu^{-}", "l");
                l.DrawClone();
            }
            pads.second->cd();
            TH1D* r = MakeRatio(ps.pos, ps.neg, ps.fname + "_r2", kBlack,
                                "#mu^{+} / #mu^{-}");
            SetRatioRange({r});
            r->Draw("hist e");
            gPad->Update();
            DrawGapMarkers({r});
            DrawRatioUnityLine();
        }

        // --- column 3: split at truth pT = 6 GeV, with the ratio ---
        {
            auto pads = MakeRatioPads(c, ps.fname + "_p3", 0.667, 0.0, 1.000, 0.875);
            pads.first->cd();
            StyleSpectrum(ps.lopt, kMagenta + 1);
            StyleSpectrum(ps.hipt, kGreen + 2);
            for (TH1D* h : {ps.lopt, ps.hipt}) {
                h->GetXaxis()->SetLabelSize(0.);
                h->GetXaxis()->SetTitleSize(0.);
            }
            SetCommonLinRange({ps.lopt, ps.hipt});
            ps.hipt->Draw("hist e");
            ps.lopt->Draw("hist e same");
            gPad->Update();
            DrawGapMarkers({ps.hipt, ps.lopt});
            Header(Form("split at truth p_{T} = %.0f GeV", kPtSplit));
            {
                const LegBox lb = AutoLegendBox({ps.lopt, ps.hipt}, 0.32,
                                                LegendHeight(2, 0.045));
                TLegend l(lb.x1, lb.y1, lb.x2, lb.y2);
                StyleLegend(l, 0.045);
                l.AddEntry(ps.lopt, Form("p_{T} < %.0f GeV", kPtSplit), "l");
                l.AddEntry(ps.hipt, Form("p_{T} #geq %.0f GeV", kPtSplit), "l");
                l.DrawClone();
            }
            pads.second->cd();
            TH1D* r = MakeRatio(ps.lopt, ps.hipt, ps.fname + "_r3", kBlack,
                                Form("#frac{p_{T} < %.0f GeV}{p_{T} #geq %.0f GeV}",
                                     kPtSplit, kPtSplit));
            SetRatioRange({r});
            r->Draw("hist e");
            gPad->Update();
            DrawGapMarkers({r});
            DrawRatioUnityLine();
        }

        CanvasGapKey(c, 0.878, 0.906, 0.020);
        CanvasHeadline(c, ps.headline, ps.sub);
        c.SaveAs((std::string(kOutDir) + "/" + ps.fname + ".png").c_str());
    }

    // ============ Plot C : pp vs Pb+Pb-overlay truth shape, with ratio ============
    // Unit area, because the two samples carry different cross-section normalisations
    // (different beam mixtures and, for the overlay, a test-sample AMI set) -- only the
    // SHAPE of the truth q*eta distribution is comparable between them.
    if (res.size() == 2 && res[0].fc.n_raw > 0 && res[1].fc.n_raw > 0) {
        TCanvas c("cC", "", 950, 850);
        auto pads = MakeRatioPads(c, "cC", 0.0, 0.0, 1.0, 0.875);
        TH1D* pp_s = (TH1D*)res[0].all->Clone("h_pp_truth_shape");
        TH1D* pb_s = (TH1D*)res[1].all->Clone("h_pb_truth_shape");
        for (TH1D* h : {pp_s, pb_s}) {
            const double integ = h->Integral("width");
            if (integ > 0) h->Scale(1. / integ);
        }
        StyleSpectrum(pp_s, kBlack);
        StyleSpectrum(pb_s, kRed + 1);
        // AFTER StyleSpectrum, which would otherwise restore the x labels onto the main
        // pad (they belong to the ratio pad below).
        for (TH1D* h : {pp_s, pb_s}) {
            h->GetYaxis()->SetTitle("normalised d#sigma_{#mu} / d(q#times#eta)");
            h->GetXaxis()->SetLabelSize(0.);
            h->GetXaxis()->SetTitleSize(0.);
        }
        pads.first->cd();
        SetCommonLinRange({pp_s, pb_s});
        pp_s->Draw("hist e");
        pb_s->Draw("hist e same");
        gPad->Update();
        DrawGapMarkers({pp_s, pb_s});
        Header("truth q#times#eta shape: pp vs Pb+Pb");
        {
            const LegBox lb = AutoLegendBox({pp_s, pb_s}, 0.49, LegendHeight(2, 0.036));
            TLegend l(lb.x1, lb.y1, lb.x2, lb.y2);
            StyleLegend(l, 0.036);
            l.AddEntry(pp_s, Form("pp full simulation (%lld #mu)", res[0].fc.n_raw), "l");
            l.AddEntry(pb_s, Form("Pb+Pb HIJING overlay (%lld #mu)", res[1].fc.n_raw),
                       "l");
            l.DrawClone();
        }
        pads.second->cd();
        TH1D* r = MakeRatio(pp_s, pb_s, "h_truth_pppb_ratio", kBlack, "pp / Pb+Pb");
        SetRatioRange({r});
        r->Draw("hist e");
        gPad->Update();
        DrawGapMarkers({r});
        DrawRatioUnityLine();
        CanvasGapKey(c, 0.885, 0.912, 0.018);
        // The fiducial IS what defines this population and sets the shape being compared, so
        // it belongs on the canvas, formatted from the constants rather than retyped.
        CanvasHeadline(
            c, "Pythia8 hard-QCD dimuon, full simulation, #sqrt{s} = #sqrt{s_{NN}} = 5.36 TeV",
            Form("truth muons, p_{T} > %.0f GeV, |#eta| < %.1f; "
                 "no reco match / working point / trigger; unit area",
                 kTruthPtMin, kTruthEtaMax),
            0.972, 0.938, 0.030, 0.026);
        c.SaveAs((std::string(kOutDir) + "/muon_truth_q_eta_spectrum_pp_vs_pbpb_shape.png")
                     .c_str());
    }

    // ---------------- numeric summary ----------------
    for (size_t i = 0; i < res.size(); ++i) {
        if (res[i].fc.n_raw == 0) continue;
        ReportDips(res[i].all, res[i].fname.c_str());
        ReportAcceptance(res[i].fname.c_str(), res[i].fc);
    }

    printf("\nPNGs written to %s\n", kOutDir);
}
