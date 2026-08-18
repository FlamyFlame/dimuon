// =================================================================================================
// plot_mc_data_pair_pt_in_eta.cxx
//
// pp24 data vs Pythia truth: the single-b dimuon pair-pT cross-section over 8-150 GeV, both
// pair-eta-integrated and split into the 9 canonical pair-eta bins. The MC counterpart of
// plots/single_b_analysis/pp24/pp24_crossx_pair_pt_in_eta_subplots.png (user, 2026-08-18).
//
// WHY THIS IS A SEPARATE MACRO FROM plot_mc_data_compr.cxx. That one overlays 1D shape
// distributions (dR, dphi, ...) taken from the Pythia PRIVATE sample. This one is a
// CROSS-SECTION comparison and must use the AMI-normalized FULL sample, because:
//   * the private sample's pair-pT axis is 30 uniform bins 0-30 GeV -- it does not reach 150, and
//     2.8 % of its weight is already in the pT overflow;
//   * its pair-eta axis has 24 uniform bins, and 4 of the 9 canonical pair-eta boundaries
//     (+-0.5, +-1.5) fall at BIN CENTRES there, so the panels could not be cut at the right
//     edges without straddling and double-counting;
//   * its selection is `from_same_b` only -- no minv window, no pair-pT threshold -- i.e. not the
//     signal region at all;
//   * and its weight never touches AMI (PythiaAlgCoreT.c:983, eventWeight / njobs), so it has no
//     cross-section normalization to convert.
//
// THE PAIR USED HERE HAS BIT-IDENTICAL AXES (verified edge by edge: all 16 pT edges from
// ParamsSet::pT_bins_150 and all 45 pair-eta edges), so both sides are cut at exactly the same
// panel boundaries and no binning mismatch is introduced (.claude/CLAUDE.md Binnings):
//
//   data : dimuon_data/pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root
//          h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts
//          = (1/L) * sum 1/(eps_trig * eps_reco), L = 400.412 pb^-1, single-b signal region
//            (m in (1.08,2.9), pair pT > 8, BOTH muons passing the fiducial gap cut), Tight WP.
//   MC   : pythia_truth_full_sample/pythia_5p36TeV/histograms_pythia_5p36TeV_no_data_resonance_cuts.root
//          h2d_sig_accept_num_pt_150_eta
//          = truth from_same_b in the truth signal region, AMI-weighted.
//
// NORMALIZATION: the AMI `crossSection` is in **nb** (PythiaAlgCoreT.c:549-557), so the MC needs
// exactly one factor, nb -> pb = 1000, and nothing else. The data is already in pb.
//
// KNOWN, DELIBERATE MISMATCHES -- state them wherever this plot is used, they are NOT bugs:
//   1. FIDUCIAL REGION. The data uses the 3 gap windows; the Pythia TRUTH acceptance was
//      deliberately left on the one-sided `q*eta < 2.2`
//      (docs/tracking/pp24_crossx_rerun_2026_08.md, "Pb+Pb / truth acceptance not migrated").
//      So the MC keeps muons the data cuts, and cuts 2.2 < q*eta < 2.3 that the data keeps.
//   2. SIGNAL DEFINITION. MC is pure truth single-b; the data is raw OS with NO same-sign
//      subtraction and no template-fit background subtraction, so it still contains
//      gluon-splitting and combinatorial background.
//   3. ISOSPIN. The full sample is the 4-beam 4:6:6:9 Pb-averaged NN cross-section, compared
//      against pp data.
//   4. The data is reco-level, corrected but NOT unfolded; the MC is truth-level.
//
// Usage:  root -l -b -q 'plot_mc_data_pair_pt_in_eta.cxx+()'
// =================================================================================================

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>

#include "../helper_functions.c"
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"

namespace {

const char* kDataFile =
    "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
    "histograms_real_pairs_pp_2024_2mu4_nominal.root";
const char* kDataHist = "h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts";

const char* kMcFile =
    "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/"
    "histograms_pythia_5p36TeV_no_data_resonance_cuts.root";
const char* kMcHist = "h2d_sig_accept_num_pt_150_eta";

// AMI crossSection is in nb; the data is normalized by a luminosity in pb^-1.
const double kNbToPb = 1.0e3;

const char* kOutDir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/";

TH2D* Fetch(const char* file, const char* hist)
{
    TFile* f = TFile::Open(file, "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error(std::string("cannot open ") + file);
    auto* h = dynamic_cast<TH2D*>(f->Get(hist));
    if (!h)
        throw std::runtime_error(std::string("no ") + hist + " in " + file);
    auto* c = static_cast<TH2D*>(h->Clone((std::string(hist) + "_clone").c_str()));
    c->SetDirectory(nullptr);
    f->Close();
    return c;
}

// The two axes MUST agree bin for bin: a panel cut at a different edge on each side would compare
// different eta ranges while every histogram still fills.
void AssertSameAxes(const TH2D* a, const TH2D* b)
{
    auto cmp = [](const TAxis* x, const TAxis* y, const char* what) {
        if (x->GetNbins() != y->GetNbins())
            throw std::runtime_error(std::string("axis bin count differs on ") + what);
        for (int i = 1; i <= x->GetNbins() + 1; ++i)
            if (std::fabs(x->GetBinLowEdge(i) - y->GetBinLowEdge(i)) > 1e-6)
                throw std::runtime_error(std::string("axis edge ") + std::to_string(i)
                                         + " differs on " + what);
    };
    cmp(a->GetXaxis(), b->GetXaxis(), "pair pT");
    cmp(a->GetYaxis(), b->GetYaxis(), "pair eta");
}

TH1D* Project(TH2D* h, double eta_lo, double eta_hi, const char* name, double scale,
              Color_t col, Style_t mstyle)
{
    const int y1 = h->GetYaxis()->FindBin(eta_lo + 1e-6);
    const int y2 = h->GetYaxis()->FindBin(eta_hi - 1e-6);
    TH1D* p = h->ProjectionX(name, y1, y2, "e");
    p->SetDirectory(nullptr);
    p->Scale(scale, "width");          // -> dsigma/dpT [pb/GeV]
    p->SetLineWidth(2);
    p->SetLineColor(col);
    p->SetMarkerColor(col);
    p->SetMarkerStyle(mstyle);
    p->SetMarkerSize(0.9);
    p->SetStats(0);
    p->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
    p->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb GeV^{-1}]");
    return p;
}

void StyleFrame(TH1D* h, double lo, double hi)
{
    h->GetYaxis()->SetRangeUser(lo, hi);
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.25);
    // On a log x axis the "10^{2}" label is tall; at the default offset the axis TITLE lands on
    // top of it. Push the title clear of the labels.
    h->GetXaxis()->SetTitleOffset(1.15);
}

}  // namespace

void plot_mc_data_pair_pt_in_eta()
{
    gStyle->SetOptStat(0);

    std::unique_ptr<TH2D> h_data(Fetch(kDataFile, kDataHist));
    std::unique_ptr<TH2D> h_mc  (Fetch(kMcFile,   kMcHist));
    AssertSameAxes(h_data.get(), h_mc.get());

    static const CommonEffcyConfig cfg{};
    const auto& eta_bins = cfg.pair_eta_proj_ranges_coarse_incl_gap;

    printf("data integral = %.6g pb ; MC integral = %.6g nb -> %.6g pb ; MC/data = %.4g\n",
           h_data->Integral(), h_mc->Integral(), h_mc->Integral() * kNbToPb,
           h_mc->Integral() * kNbToPb / h_data->Integral());

    // ---------------------------------------------------------------- pair-eta integrated
    {
        TH1D* d = Project(h_data.get(), -2.4, 2.4, "hpt_int_data", 1.0,      kGreen + 2, 20);
        TH1D* m = Project(h_mc.get(),   -2.4, 2.4, "hpt_int_mc",   kNbToPb,  kBlue,      21);

        double lo = 1e300, hi = -1e300;
        for (TH1D* h : {d, m})
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                const double v = h->GetBinContent(b);
                if (v > 0.) { lo = std::min(lo, v); hi = std::max(hi, v); }
            }
        StyleFrame(d, lo * 0.3, hi * 3.0);

        TCanvas c("c_int", "pair pT, pair-eta integrated", 700, 600);
        c.SetLogx(); c.SetLogy();
        c.SetLeftMargin(0.16); c.SetBottomMargin(0.15);
        d->Draw("E");
        m->Draw("E,same");
        TLegend l(0.20, 0.20, 0.55, 0.36);
        l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextFont(42); l.SetTextSize(0.038);
        l.AddEntry(d, "pp data 2024", "lp");
        l.AddEntry(m, "Pythia", "lp");
        l.Draw();
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.038);
        t.DrawLatex(0.20, 0.86, "|#eta^{pair}| < 2.4");
        c.SaveAs((std::string(kOutDir) + "pair_pt_mc_data_compr.png").c_str());
    }

    // ---------------------------------------------------------------- 9 pair-eta panels
    {
        int nrow = 1, ncol = 1;
        DetermineSubplotGrid(static_cast<int>(eta_bins.size()), nrow, ncol);

        std::vector<TH1D*> ds, ms;
        double lo = 1e300, hi = -1e300;
        for (size_t i = 0; i < eta_bins.size(); ++i) {
            TH1D* d = Project(h_data.get(), eta_bins[i].first, eta_bins[i].second,
                              ("hd" + std::to_string(i)).c_str(), 1.0,     kGreen + 2, 20);
            TH1D* m = Project(h_mc.get(),   eta_bins[i].first, eta_bins[i].second,
                              ("hm" + std::to_string(i)).c_str(), kNbToPb, kBlue,      21);
            ds.push_back(d); ms.push_back(m);
            for (TH1D* h : {d, m})
                for (int b = 1; b <= h->GetNbinsX(); ++b) {
                    const double v = h->GetBinContent(b);
                    if (v > 0.) { lo = std::min(lo, v); hi = std::max(hi, v); }
                }
        }

        // ONE common log-y range across all panels, so the nine are directly comparable by eye.
        TCanvas c("c_eta", "pair pT in pair-eta bins", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);
        for (size_t i = 0; i < eta_bins.size(); ++i) {
            c.cd(static_cast<int>(i) + 1);
            gPad->SetLogx(); gPad->SetLogy();
            gPad->SetLeftMargin(0.18); gPad->SetBottomMargin(0.17);
            StyleFrame(ds[i], lo * 0.3, hi * 3.0);
            ds[i]->Draw("E");
            ms[i]->Draw("E,same");

            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.05);
            t.DrawLatex(0.22, 0.87, Form("#eta^{pair} #in [%.1f, %.1f]",
                                         eta_bins[i].first, eta_bins[i].second));
            if (i == 0) {
                TLegend* l = new TLegend(0.22, 0.20, 0.62, 0.36);
                l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(42); l->SetTextSize(0.05);
                l->AddEntry(ds[i], "pp data 2024", "lp");
                l->AddEntry(ms[i], "Pythia", "lp");
                l->Draw();
            }
        }
        c.SaveAs((std::string(kOutDir) + "pair_pt_in_eta_subplots_mc_data_compr.png").c_str());
    }
}
