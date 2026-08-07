// =================================================================================================
// plot_forward_qeta_edge_scan.cxx
//
// DECISION PLOT: where should the forward edge of the q*eta fiducial gap cut sit?
//
// `ParamsSet::single_mu_fiducial_gap_cuts` currently rejects q*eta in (2.20, 2.40). The forward
// acceptance collapses one-sidedly -- muon_gap_cuts_acceptance.md F6 measures the pp yield at
// q*eta = 2.38-2.40 as only 8.8% of its 2.20-2.22 value, while the NEGATIVE side still holds ~71%
// out to -2.4. The open question is whether 2.20 is really necessary or whether the cut can be
// loosened to 2.25 or 2.30, which would recover 1.4% (pp) / 2.2% (PbPb) of muons.
//
// WHAT THIS DRAWS: the single-muon mu4 efficiency vs pT in the MOST POSITIVE q*eta bin, with the
// bin's UPPER edge scanned over 2.20 / 2.25 / 2.30 / 2.40 (lower edge fixed at 2.0). If the 2.40
// curve is visibly degraded relative to 2.20 but 2.25 / 2.30 are not, the cut can be loosened.
// (2.40 is known to be bad -- it is drawn as the reference for "how bad".)
//
// LAYOUT (user): 4 panels -- mu+ LEFT, mu- RIGHT, data TOP, MC BOTTOM.
//
// IMPORTANT -- INPUTS MUST HAVE NO FORWARD GAP CUT. The nominal chain now rejects q*eta > 2.2
// outright, so its outputs contain nothing to compare. This macro therefore reads files produced
// WITHOUT the gap cut, and refuses to run if the region above 2.2 is empty (which is exactly what
// a gap-cut input would look like) rather than silently drawing four identical curves.
//
// The efficiency is formed by projecting the 2D (q*eta, pT) numerator and denominator over each
// candidate q*eta range and dividing -- so the scan is independent of the fitted q*eta binning.
//   data: h_pt2nd_vs_q_eta2nd_<sign1|sign2>_2mu4_sepr  /  ..._<sign>_mu4_sepr   (tag-and-probe)
//   MC  : h_mc_pt_vs_q_eta_num_<muplus|muminus>        /  h_mc_pt_vs_q_eta_denom_<chg>
// Both are (x = q*eta, y = pT).
//
// Errors: TGraphAsymmErrors::BayesDivide, matching how the nominal turn-on points are built.
//
// Run (from this directory):
//   root -l -b -q 'plot_forward_qeta_edge_scan.cxx+'
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
#include <stdexcept>
#include <string>
#include <vector>

#include "dr_correction_sample_cfg.h"

namespace {

// The candidate upper edges. Lower edge fixed at 2.0 (the last coarse q*eta bin's lower edge).
const double kQEtaLo = 2.0;
const std::vector<double> kUpperEdges = {2.20, 2.25, 2.30, 2.40};
const std::vector<Color_t> kEdgeCol   = {kBlack, kBlue + 1, kGreen + 2, kRed + 1};
const std::vector<Style_t> kEdgeMark  = {20, 21, 22, 23};

template <typename T>
T* GetObj(TFile* f, const std::string& name) {
    auto* o = dynamic_cast<T*>(f->Get(name.c_str()));
    if (!o) throw std::runtime_error("plot_forward_qeta_edge_scan: missing '" + name +
                                     "' in " + f->GetName());
    return o;
}

// eff(pT) for q*eta in [lo, hi), by projecting both 2D hists onto pT over that x-range.
TGraphAsymmErrors* EffInQEtaRange(TH2D* num, TH2D* den, double lo, double hi,
                                  const std::string& tag, double& n_den_entries) {
    const TAxis* ax = den->GetXaxis();
    // half-open [lo, hi): FindBin(hi) would include the bin containing hi, so step back when hi
    // lands exactly on an edge (it does, by construction of the q*eta axis).
    const int b1 = ax->FindBin(lo + 1e-6);
    int b2 = ax->FindBin(hi - 1e-6);
    if (b2 < b1) b2 = b1;
    TH1D* n = num->ProjectionY((tag + "_n").c_str(), b1, b2, "e");
    TH1D* d = den->ProjectionY((tag + "_d").c_str(), b1, b2, "e");
    n_den_entries = d->Integral();
    auto* g = new TGraphAsymmErrors();
    g->BayesDivide(n, d);
    g->SetName((tag + "_g").c_str());
    delete n; delete d;
    return g;
}

void DrawPanel(TH2D* num, TH2D* den, const std::string& panel_title, const std::string& tag,
               bool draw_legend) {
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

    auto* leg = new TLegend(0.45, 0.18, 0.93, 0.44);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);

    for (size_t i = 0; i < kUpperEdges.size(); ++i) {
        double nden = 0;
        auto* g = EffInQEtaRange(num, den, kQEtaLo, kUpperEdges[i],
                                 tag + "_" + std::to_string(i), nden);
        g->SetMarkerStyle(kEdgeMark[i]);
        g->SetMarkerColor(kEdgeCol[i]);
        g->SetLineColor(kEdgeCol[i]);
        g->SetLineWidth(2);
        g->Draw("P same");
        leg->AddEntry(g, Form("%.1f < q#eta < %.2f", kQEtaLo, kUpperEdges[i]), "lp");
    }
    if (draw_legend) leg->Draw();

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.048);
    t.DrawLatex(0.17, 0.90, panel_title.c_str());
}

// Guard: an input that already has the forward gap cut applied has NOTHING above q*eta = 2.2,
// which would silently render four identical curves and invite exactly the wrong conclusion.
void AssertForwardRegionPopulated(TH2D* den, const std::string& what) {
    const TAxis* ax = den->GetXaxis();
    const double above = den->Integral(ax->FindBin(2.2 + 1e-6), ax->FindBin(2.4 - 1e-6),
                                       1, den->GetNbinsY());
    if (above <= 0.)
        throw std::runtime_error(
            "plot_forward_qeta_edge_scan: " + what + " has NO entries with q*eta > 2.2. That is "
            "what a gap-cut input looks like, and every candidate edge would draw the same curve. "
            "Re-run the producer with the forward fiducial gap cut DISABLED.");
}

}  // namespace

void plot_forward_qeta_edge_scan(const std::string& sample = "pp", bool use_tight_wp = true) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const bool is_pbpb = (sample == "pbpb");
    if (!is_pbpb && sample != "pp")
        throw std::runtime_error("plot_forward_qeta_edge_scan: sample must be 'pp' or 'pbpb'");

    const std::string wp  = use_tight_wp ? "" : "_medium_wp";
    const std::string wpt = use_tight_wp ? "Tight" : "Medium";
    const std::string home = std::string(gSystem->Getenv("HOME"));

    // NO-GAP-CUT inputs (see the header note). The data files are the pre-round-8 fine-q*eta
    // outputs, which predate the fiducial gap cut and so still contain q*eta > 2.2.
    const std::string data_file = is_pbpb
        ? home + "/usatlasdata/dimuon_data/pbpb_2023/"
                 "histograms_real_pairs_pbpb_2023_single_mu4_fine_q_eta_bin" + wp + ".root"
        : home + "/usatlasdata/dimuon_data/pp_2024/"
                 "histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin" + wp + ".root";
    const std::string mc_file = is_pbpb
        ? home + "/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/"
                 "mc_trig_eff_hists_hijing_overlay_pbpb23" + wp + "_nogapcut.root"
        : home + "/usatlasdata/pythia_fullsim_full_sample/"
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
    TH2D* dnum_p = GetObj<TH2D>(fd, "h_pt2nd_vs_q_eta2nd_sign1_2mu4_sepr");
    TH2D* dden_p = GetObj<TH2D>(fd, "h_pt2nd_vs_q_eta2nd_sign1_mu4_sepr");
    TH2D* dnum_m = GetObj<TH2D>(fd, "h_pt2nd_vs_q_eta2nd_sign2_2mu4_sepr");
    TH2D* dden_m = GetObj<TH2D>(fd, "h_pt2nd_vs_q_eta2nd_sign2_mu4_sepr");
    TH2D* mnum_p = GetObj<TH2D>(fm, "h_mc_pt_vs_q_eta_num_muplus");
    TH2D* mden_p = GetObj<TH2D>(fm, "h_mc_pt_vs_q_eta_denom_muplus");
    TH2D* mnum_m = GetObj<TH2D>(fm, "h_mc_pt_vs_q_eta_num_muminus");
    TH2D* mden_m = GetObj<TH2D>(fm, "h_mc_pt_vs_q_eta_denom_muminus");

    AssertForwardRegionPopulated(dden_p, "data mu+");
    AssertForwardRegionPopulated(mden_p, "MC mu+");

    TCanvas c("c_fwd_edge", "", 1400, 1100);
    c.Divide(2, 2);
    c.cd(1); DrawPanel(dnum_p, dden_p, "data, #mu^{+}", "d_p", true);
    c.cd(2); DrawPanel(dnum_m, dden_m, "data, #mu^{-}", "d_m", false);
    c.cd(3); DrawPanel(mnum_p, mden_p, "MC, #mu^{+}",   "m_p", false);
    c.cd(4); DrawPanel(mnum_m, mden_m, "MC, #mu^{-}",   "m_m", false);

    c.cd(0);
    TLatex head;
    head.SetNDC();
    head.SetTextFont(42);
    head.SetTextSize(0.022);
    head.DrawLatex(0.03, 0.982,
                   is_pbpb
                     ? Form("Pb+Pb  #sqrt{s_{NN}} = 5.36 TeV,  %s muons,  forward q#eta edge scan",
                            wpt.c_str())
                     : Form("pp  #sqrt{s} = 5.36 TeV,  %s muons,  forward q#eta edge scan",
                            wpt.c_str()));

    const std::string out_dir = home + "/usatlasdata/dimuon_data/plots/"
                              + std::string(is_pbpb ? "pbpb" : "pp") + "_trigger_efficiency/"
                                "mc_based" + std::string(use_tight_wp ? "" : "_medium") + "/forward_qeta_edge_scan/"
                              ;   // WP is in the top-level tree, not a subdirectory
    gSystem->mkdir(out_dir.c_str(), kTRUE);
    const std::string png = out_dir + "single_mu_eff_forward_qeta_edge_scan.png";
    c.SaveAs(png.c_str());
    std::cout << "  wrote " << png << std::endl;

    // Print the numbers behind the figure: the plateau efficiency and the probe count per edge,
    // so the decision can be read off values rather than eyeballed off a curve.
    std::cout << "\n  plateau eff (pT > 8 GeV) and probe count per candidate edge:\n";
    struct Src { const char* tag; TH2D* n; TH2D* d; };
    for (const Src& s : {Src{"data mu+", dnum_p, dden_p}, Src{"data mu-", dnum_m, dden_m},
                         Src{"MC   mu+", mnum_p, mden_p}, Src{"MC   mu-", mnum_m, mden_m}}) {
        std::cout << "    " << s.tag << ":";
        for (double hi : kUpperEdges) {
            const TAxis* ax = s.d->GetXaxis();
            const int b1 = ax->FindBin(kQEtaLo + 1e-6), b2 = ax->FindBin(hi - 1e-6);
            TH1D* n = s.n->ProjectionY("tn", b1, b2, "e");
            TH1D* d = s.d->ProjectionY("td", b1, b2, "e");
            const int p1 = d->FindBin(8.0 + 1e-6);
            const double N = n->Integral(p1, d->GetNbinsX());
            const double D = d->Integral(p1, d->GetNbinsX());
            std::cout << Form("  [2.0,%.2f) %.4f (D=%.3g)", hi, D > 0 ? N / D : 0., D);
            delete n; delete d;
        }
        std::cout << "\n";
    }
    fd->Close();
    fm->Close();
}
