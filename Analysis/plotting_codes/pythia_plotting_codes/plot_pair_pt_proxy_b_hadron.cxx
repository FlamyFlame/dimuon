// How good a proxy is the dimuon pair pT for the parent b-hadron pT? (truth level, single-b OS pairs)
//
// Pair population: opposite-sign truth pair tree, `from_same_b` (both muons from one b-hadron decay chain).
// b-hadron pT: pT of the LAST b hadron in the chain of the leading muon (m1; m1 is the leading muon in
// both trees). For the same-b pairs the last b hadron of the subleading muon is the same particle in
// all but ~1e-4 of the pairs.
//
// "Truth signal selection" = SignalAcceptanceDefinition::NumExpr(): from_same_b, truth pair pT > 9 GeV,
// |truth eta^pair| < 2.2, truth mass window. NO detector-gap cuts. The single-muon truth fiducial
// (pT > 4.5 GeV, |eta| < 2.4) is applied at the NTuple stage in BOTH trees, so "from_same_b only" below
// means from_same_b + that muon fiducial.
//
// Samples (both are the NTuple-processing output, tree muon_pair_tree_sign2, per-pair `weight`):
//   central truth : CENTRAL Pythia8 5.36 TeV truth-only evgen (4 isospin beams 4:6:6:9, nuclear PDF)
//   pp24 fullsim  : CENTRAL Pythia8 5.36 TeV full-simulation FULL production (pp beam only, proton PDF);
//                   only the TRUTH columns are used.
// UNITS: `weight` carries the AMI cross-section in nb (docs/ami_weights.md) => d sigma/dR is in nb.
//
// Outputs (PNG only) in OutputDir():
//   pair_pt_over_b_hadron_pt.png            R = pT/pT(b) for pair, leading muon, subleading muon  (central, signal)
//   pair_pt_vs_b_hadron_pt.png              2D (b-hadron pT, pair pT) (central, signal)
//   sanity_sample_pair_pt_over_b_hadron_pt.png   pp24 fullsim truth vs central truth (unit area)
//   sanity_cuts_pair_pt_over_b_hadron_pt.png     from_same_b only vs truth signal selection (unit area)

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cstdio>
#include <string>

#include "../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../Utilities/SignalAcceptanceDefinition.h"

namespace {

const char* kTruthFile =
    "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/"
    "muon_pairs_pythia_5p36TeV_no_data_resonance_cuts.root";

std::string FullsimFile() {
    return FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/false)
         + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root";
}
std::string OutputDir() {
    return "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/plots/truth_level_studies/pair_pt_proxy_b_hadron/";
}

const int    kNbinsR = 50;
const double kRmin = 0., kRmax = 1.;

// pairs with `from_same_b` and a valid last-b-hadron pT; R columns
ROOT::RDF::RNode Prepare(const char* file, const std::string& selection) {
    ROOT::RDataFrame df("muon_pair_tree_sign2", file);
    return ROOT::RDF::RNode(df.Filter(selection)
        .Define("b_pt",
            [](const std::vector<float>& v){ return v.empty() ? 0.f : v[0] / 1000.f; },
            {"m1_last_b_hadron_prt_pt_eta_phi_m"})
        .Filter("b_pt > 0")
        .Define("R_pair", "truth_pair_pt / b_pt")
        .Define("R_m1",   "m1.truth_pt / b_pt")
        .Define("R_m2",   "m2.truth_pt / b_pt")
        .Define("wR_pair", "weight * R_pair")
        .Define("wR_m1",   "weight * R_m1")
        .Define("wR_m2",   "weight * R_m2"));
}

struct Dist {
    ROOT::RDF::RResultPtr<TH1D> h;
    ROOT::RDF::RResultPtr<double> sumw, sumwR;
    ROOT::RDF::RResultPtr<ULong64_t> n;
    TH1D* hc = nullptr;
    double mean = 0, norm_all = 0;
};
Dist BookR(ROOT::RDF::RNode df, const std::string& var, const std::string& name) {
    Dist d;
    d.h     = df.Histo1D({name.c_str(), "", kNbinsR, kRmin, kRmax}, var, "weight");
    d.sumw  = df.Sum("weight");
    d.sumwR = df.Sum("w" + var);
    d.n     = df.Count();
    return d;
}
void Finish(Dist& d) {
    d.hc = (TH1D*) d.h->Clone();
    d.hc->SetDirectory(nullptr);
    d.norm_all = *d.sumw;
    d.mean = *d.sumwR / *d.sumw;
}

void Style(TH1D* h, int color, int marker) {
    h->SetLineColor(color); h->SetMarkerColor(color);
    h->SetMarkerStyle(marker); h->SetMarkerSize(0.8); h->SetLineWidth(1);
}
// (1/N) dN/dR
void UnitArea(TH1D* h) { h->Scale(1. / h->Integral()); h->Scale(1., "width"); }

void Header(const char* l1, const char* l2, double x, double y) {
    TLatex t; t.SetNDC(); t.SetTextFont(43); t.SetTextSize(15);
    t.DrawLatex(x, y, l1);
    if (l2) t.DrawLatex(x, y - 0.045, l2);
}
std::string Mean(const char* name, double m) { return Form("%s, #LTR#GT = %.3f", name, m); }

}  // namespace

void plot_pair_pt_proxy_b_hadron() {
    gSystem->mkdir(OutputDir().c_str(), true);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const std::string sel_signal = SignalAcceptanceDefinition::NumExpr();   // incl. from_same_b
    const std::string sel_sameb  = "from_same_b";

    auto df_c_sig = Prepare(kTruthFile, sel_signal);
    auto df_c_sb  = Prepare(kTruthFile, sel_sameb);
    const std::string fs_file = FullsimFile();
    auto df_f_sig = Prepare(fs_file.c_str(), sel_signal);

    // ---- book everything, single event loop per input tree ----
    Dist c_pair = BookR(df_c_sig, "R_pair", "c_pair");
    Dist c_m1   = BookR(df_c_sig, "R_m1",   "c_m1");
    Dist c_m2   = BookR(df_c_sig, "R_m2",   "c_m2");
    Dist f_pair = BookR(df_f_sig, "R_pair", "f_pair");
    Dist s_pair = BookR(df_c_sb,  "R_pair", "s_pair");   // from_same_b only (central)

    // 2D: b-hadron pT (x) vs pair pT (y) on the canonical fine pair-pT axis (ParamsSet::pT_bins_150)
    ParamsSet pms;
    const std::vector<double> edges = pms.pT_bins_150;
    const int nb = (int) edges.size() - 1;
    auto h2 = df_c_sig.Histo2D({"h2_pair_vs_b", "", nb, edges.data(), nb, edges.data()}, "b_pt", "truth_pair_pt", "weight");
    auto w2_all = df_c_sig.Sum("weight");

    Finish(c_pair); Finish(c_m1); Finish(c_m2); Finish(f_pair); Finish(s_pair);
    TH2D* h2c = (TH2D*) h2->Clone("h2c");
    h2c->SetDirectory(nullptr);

    printf("[proxy] central truth, signal selection: pairs=%llu  sum(w)=%.6g nb\n", (ULong64_t)*c_pair.n, c_pair.norm_all);
    printf("[proxy]   <R_pair>=%.4f  <R_lead>=%.4f  <R_sublead>=%.4f\n", c_pair.mean, c_m1.mean, c_m2.mean);
    printf("[proxy] pp24 fullsim truth, signal sel.: pairs=%llu  <R_pair>=%.4f\n", (ULong64_t)*f_pair.n, f_pair.mean);
    printf("[proxy] central truth, from_same_b only: pairs=%llu  <R_pair>=%.4f\n", (ULong64_t)*s_pair.n, s_pair.mean);
    for (auto* d : {&c_pair, &c_m1, &c_m2, &f_pair, &s_pair})
        printf("[proxy]   weight in R overflow: %.3e (of %.6g)\n", d->hc->GetBinContent(kNbinsR + 1), d->norm_all);
    printf("[proxy] 2D: weight inside the pT_bins_150 box = %.4f of the total\n", h2c->Integral() / *w2_all);

    // ===================== 1. R for pair / leading / subleading =====================
    {
        TH1D* hp = (TH1D*) c_pair.hc->Clone("hp"); TH1D* h1 = (TH1D*) c_m1.hc->Clone("h1"); TH1D* h2_ = (TH1D*) c_m2.hc->Clone("h2_");
        for (TH1D* h : {hp, h1, h2_}) h->Scale(1., "width");   // d sigma / dR [nb]
        Style(hp, kRed, 20); Style(h1, kBlack, 20); Style(h2_, kBlue, 20);
        const double ymax = std::max({hp->GetMaximum(), h1->GetMaximum(), h2_->GetMaximum()});

        TCanvas c("c1", "c1", 700, 600);
        c.SetLeftMargin(0.17); c.SetBottomMargin(0.13); c.SetRightMargin(0.05); c.SetTopMargin(0.13);
        hp->GetXaxis()->SetTitle("R = p_{T} / p_{T}^{b hadron}"); hp->GetXaxis()->SetTitleOffset(1.2);
        hp->GetYaxis()->SetTitle("d#sigma / dR [nb]");
        hp->GetYaxis()->SetRangeUser(0., ymax * 1.4);
        hp->Draw("E"); h1->Draw("E same"); h2_->Draw("E same");
        TLegend leg(0.40, 0.62, 0.94, 0.88);
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(43); leg.SetTextSize(16);
        leg.AddEntry(hp, Mean("pair p_{T}", c_pair.mean).c_str(), "ep");
        leg.AddEntry(h1, Mean("leading #mu p_{T}", c_m1.mean).c_str(), "ep");
        leg.AddEntry(h2_, Mean("subleading #mu p_{T}", c_m2.mean).c_str(), "ep");
        leg.Draw();
        Header("Pythia 8, #sqrt{s} = 5.36 TeV, truth level", "Single-b opposite-sign #mu#mu, signal selection", 0.17, 0.945);
        c.SaveAs((OutputDir() + "pair_pt_over_b_hadron_pt.png").c_str());
    }

    // ===================== 2. 2D pair pT vs b-hadron pT =====================
    {
        TH2D* h = (TH2D*) h2c->Clone("h2draw");
        h->Scale(1., "width");   // d^2 sigma / (d p_T^b d p_T^pair) [nb/GeV^2]
        double zmin = 1e300;     // log-z floor = smallest populated cell, so no filled cell is drawn blank
        for (int i = 1; i <= h->GetNbinsX(); i++) for (int j = 1; j <= h->GetNbinsY(); j++)
            if (h->GetBinContent(i, j) > 0) zmin = std::min(zmin, h->GetBinContent(i, j));
        h->SetMinimum(0.5 * zmin);
        TCanvas c("c2", "c2", 760, 650);
        // log x/y: pT_bins_150 is log-binned. log z: d^2sigma/dpT dpT spans ~8 decades (steeply falling pT-hat
        // spectrum); linear z would show only the lowest-pT corner.
        c.SetLogx(); c.SetLogy(); c.SetLogz();
        c.SetLeftMargin(0.16); c.SetBottomMargin(0.13); c.SetRightMargin(0.17); c.SetTopMargin(0.11);
        h->GetXaxis()->SetTitle("p_{T}^{b hadron} [GeV]");
        h->GetYaxis()->SetTitle("p_{T}^{pair} [GeV]");
        h->GetZaxis()->SetTitle("d^{2}#sigma / (dp_{T}^{b hadron} dp_{T}^{pair}) [nb/GeV^{2}]");
        h->GetZaxis()->SetTitleOffset(1.5);
        h->GetXaxis()->SetMoreLogLabels(); h->GetYaxis()->SetTitleOffset(1.4);
        h->Draw("COLZ");
        TLine diag(edges.front(), edges.front(), edges.back(), edges.back());
        diag.SetLineStyle(2); diag.SetLineWidth(2); diag.Draw();
        Header("Pythia 8, #sqrt{s} = 5.36 TeV, truth level", "Single-b opposite-sign #mu#mu, signal selection", 0.16, 0.945);
        c.SaveAs((OutputDir() + "pair_pt_vs_b_hadron_pt.png").c_str());
    }

    // ===================== sanity plots: unit-area overlay + ratio pad =====================
    auto sanity = [&](TH1D* a, const std::string& la, TH1D* b, const std::string& lb,
                      const char* ratio_title, const char* l1, const char* l2, const char* out, double rlo, double rhi) {
        TH1D* ha = (TH1D*) a->Clone(); TH1D* hb = (TH1D*) b->Clone();
        UnitArea(ha); UnitArea(hb);
        Style(ha, kBlack, 20); Style(hb, kRed, 24);
        TH1D* hr = (TH1D*) hb->Clone(); hr->Divide(ha);
        TCanvas c("cs", "cs", 700, 700);
        TPad p1("p1", "", 0, 0.30, 1, 1), p2("p2", "", 0, 0, 1, 0.30);
        p1.SetBottomMargin(0.02); p1.SetLeftMargin(0.17); p1.SetRightMargin(0.05); p1.SetTopMargin(0.12);
        p2.SetTopMargin(0.03); p2.SetBottomMargin(0.36); p2.SetLeftMargin(0.17); p2.SetRightMargin(0.05);
        p1.Draw(); p2.Draw();
        p1.cd();
        ha->GetYaxis()->SetTitle("(1/N) dN / dR");
        ha->GetYaxis()->SetRangeUser(0., std::max(ha->GetMaximum(), hb->GetMaximum()) * 1.45);
        ha->GetXaxis()->SetLabelSize(0); ha->GetXaxis()->SetTitleSize(0);
        ha->Draw("E"); hb->Draw("E same");
        TLegend leg(0.19, 0.64, 0.94, 0.88);
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(43); leg.SetTextSize(15);
        leg.AddEntry(ha, la.c_str(), "ep"); leg.AddEntry(hb, lb.c_str(), "ep"); leg.Draw();
        Header(l1, l2, 0.17, 0.945);
        p2.cd();
        hr->GetXaxis()->SetTitle("R = p_{T}^{pair} / p_{T}^{b hadron}");
        hr->GetYaxis()->SetTitle(ratio_title);
        hr->GetYaxis()->SetRangeUser(rlo, rhi);
        hr->GetYaxis()->SetNdivisions(505);
        hr->GetXaxis()->SetTitleSize(0.13); hr->GetXaxis()->SetLabelSize(0.11); hr->GetXaxis()->SetTitleOffset(1.1);
        hr->GetYaxis()->SetTitleSize(0.085); hr->GetYaxis()->SetLabelSize(0.10); hr->GetYaxis()->SetTitleOffset(0.75);
        hr->Draw("E");
        TLine one(kRmin, 1., kRmax, 1.); one.SetLineStyle(2); one.Draw();
        c.SaveAs((OutputDir() + out).c_str());
    };

    // sample: pp24 fullsim (pp, proton PDF) truth vs central (4 isospin, nuclear PDF) truth; both signal selection
    sanity(c_pair.hc, Mean("4 isospin beams, nuclear PDF", c_pair.mean),
           f_pair.hc, Mean("pp beam, proton PDF", f_pair.mean),
           "pp / 4 isospin",
           "Pythia 8, #sqrt{s} = 5.36 TeV, truth level",
           "Single-b opposite-sign #mu#mu, signal selection",
           "sanity_sample_pair_pt_over_b_hadron_pt.png", 0.8, 1.2);

    // selection: central truth, signal selection vs from_same_b only
    sanity(s_pair.hc, Mean("same b hadron only", s_pair.mean),
           c_pair.hc, Mean("same b hadron + signal selection", c_pair.mean),
           "with / without cuts",
           "Pythia 8, #sqrt{s} = 5.36 TeV, truth level, 4 isospin beams, nuclear PDF",
           "Single-b opposite-sign #mu#mu",
           "sanity_cuts_pair_pt_over_b_hadron_pt.png", 0.0, 1.6);
}
