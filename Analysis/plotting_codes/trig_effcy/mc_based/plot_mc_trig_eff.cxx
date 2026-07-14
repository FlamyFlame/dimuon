// plot_mc_trig_eff.cxx
// Step-1/2/3 MC trigger-efficiency plots (mc_trigger_efficiency.md §3.1-§3.3).
//
//   Step 1 (§3.1): MC direct P[mu4-matched | offline reco muon] vs data tag-and-probe
//                  P[2mu4 | mu4-tag, DeltaR>0.8] (the data eps^nc estimate), per charge.
//   Step 2 (§3.2): DeltaR-binned MC singles efficiency (factorization cross-check),
//                  3 DeltaR bins overlaid + Step-1 inclusive reference.
//   Step 3 (§3.3): eps_dR(dR) = inverse-weighted num / unweighted denom (TH1::Divide with
//                  error propagation -- weights > 1 make Bayes invalid, same convention as
//                  the data cross-term), plateau = weighted mean over dR in [1,3].
//
// One macro for both samples:
//   plot_mc_trig_eff("pp")      : Pythia8 pp24 fullsim  vs pp24 data       -> eps_dR^2mu4
//   plot_mc_trig_eff("overlay") : HIJING overlay PbPb23 vs PbPb23 data 0-5% (D2) -> eps_dR^cross
//
// WP config (registry: Analysis/docs/muon_wp_registry.md): use_tight_wp default TRUE (Tight
// nominal, unsuffixed inputs); false selects the _medium_wp MC input files + label text.
//
// Data sign convention (evidence in docs/tracking/_sub_mctrig_plots.md F3):
//   sign1 = mu+ (RDFBasedHistFillingPP.cxx:168 Filter("charge2nd > 0");
//   SingleMuEffcyPtTurnOnFitter.cxx:153-154), sign2 = mu-.
//
// Compile/run (ACLiC, from this directory):
//   root -l -b -q 'plot_mc_trig_eff.cxx+("pp")'
//   root -l -b -q 'plot_mc_trig_eff.cxx+("overlay")'

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------- helpers

template <typename T>
T* GetObj(TFile* f, const std::string& name)
{
    T* obj = dynamic_cast<T*>(f->Get(name.c_str()));
    if (!obj)
        throw std::runtime_error("plot_mc_trig_eff: missing object '" + name +
                                 "' in file " + f->GetName());
    return obj;
}

TFile* OpenFile(const std::string& path)
{
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("plot_mc_trig_eff: cannot open " + path);
    return f;
}

// Bayesian efficiency graph (efficiencies bounded in [0,1]; num subset of denom).
TGraphAsymmErrors* BayesEff(TH1* num, TH1* denom)
{
    auto* g = new TGraphAsymmErrors();
    g->Divide(num, denom, "cl=0.683 b(1,1) mode");
    return g;
}

void StyleGraph(TGraphAsymmErrors* g, Color_t col, Style_t marker,
                double msize = 0.9, Width_t lw = 2, Style_t ls = 1)
{
    g->SetMarkerColor(col);
    g->SetLineColor(col);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(msize);
    g->SetLineWidth(lw);
    g->SetLineStyle(ls);
}

TH1* DrawEffFrame(double xlo, double xhi, const std::string& xtitle,
                  double ylo = 0.0, double yhi = 1.1,
                  const std::string& ytitle = "efficiency")
{
    TH1* frame = gPad->DrawFrame(xlo, ylo, xhi, yhi);
    frame->GetXaxis()->SetTitle(xtitle.c_str());
    frame->GetYaxis()->SetTitle(ytitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    return frame;
}

void DrawUnityLine(double xlo, double xhi)
{
    auto* l = new TLine(xlo, 1.0, xhi, 1.0);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw("same");
}

// SF-ratio graph from TH1::Divide output, dropping bins with no denominator content
// or zero width (the duplicated-8.0-edge bin of the data pt binning; BayesDivide skips
// those automatically, TGraphAsymmErrors(hist) does not).
TGraphAsymmErrors* DivideGraphClean(TH1D* num, TH1D* den)
{
    auto* h = static_cast<TH1D*>(num->Clone());
    h->Divide(num, den, 1.0, 1.0);
    auto* g = new TGraphAsymmErrors();
    int ip = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        if (den->GetBinContent(b) <= 0 || h->GetXaxis()->GetBinWidth(b) <= 0) continue;
        g->SetPoint(ip, h->GetXaxis()->GetBinCenter(b), h->GetBinContent(b));
        g->SetPointError(ip, h->GetXaxis()->GetBinWidth(b)/2, h->GetXaxis()->GetBinWidth(b)/2,
                         h->GetBinError(b), h->GetBinError(b));
        ++ip;
    }
    delete h;
    return g;
}

void SaveCanvas(TCanvas& c, const std::string& path)
{
    c.SaveAs(path.c_str());
    std::cout << "  wrote " << path << std::endl;
}

// ---------------------------------------------------------------- config

struct SampleCfg {
    std::string mc_dir;
    std::string mc_label;      // file label
    std::string data_file;
    std::string ctr;           // data ctr token ("" or "_ctr0_5")
    std::string out_base;      // plots root for this sample
    std::string sample_text;   // canvas headline (without WP)
    std::string data_text;     // data legend entry
    std::string eps_dr_text;   // eps_DR symbol for step 3
};

SampleCfg MakeCfg(const std::string& sample)
{
    SampleCfg c;
    if (sample == "pp") {
        c.mc_dir      = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
        c.mc_label    = "pp24";
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                        "histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin.root";
        c.ctr         = "";
        c.out_base    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "pp_trigger_efficiency/mc_based/";
        c.sample_text = "Pythia8 pp24 fullsim";
        c.data_text   = "pp24 data";
        c.eps_dr_text = "#varepsilon_{#DeltaR}^{2mu4}";
    } else if (sample == "overlay") {
        c.mc_dir      = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/";
        c.mc_label    = "hijing_overlay_pbpb23";
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/"
                        "histograms_real_pairs_pbpb_2023_single_mu4_fine_q_eta_bin.root";
        c.ctr         = "_ctr0_5";   // D2: overlay compares ONLY to PbPb23 data 0-5%
        c.out_base    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "pbpb_trigger_efficiency/mc_based/";
        c.sample_text = "HIJING overlay Pb+Pb23 cond., 0-5%";
        c.data_text   = "Pb+Pb23 data 0-5%";
        c.eps_dr_text = "#varepsilon_{#DeltaR}^{cross}";
    } else {
        throw std::runtime_error("plot_mc_trig_eff: sample must be 'pp' or 'overlay', got " + sample);
    }
    return c;
}

// charge <-> data-sign mapping (F3): sign1 = mu+, sign2 = mu-
const std::vector<std::string> kCharges     = {"muplus", "muminus"};
const std::vector<std::string> kDataSigns   = {"sign1", "sign2"};
const std::vector<std::string> kChargeTex   = {"#mu^{+}", "#mu^{-}"};

// fine q.eta bins (CommonEffcyConfig.h q_eta_proj_ranges_fine_excl_gap)
const std::vector<std::string> kQEtaSuffix = {
    "minus2_40_TO_minus2_00", "minus2_00_TO_minus1_60", "minus1_60_TO_minus1_30",
    "minus0_90_TO_minus0_50", "minus0_50_TO_minus0_10", "0_10_TO_0_50",
    "0_50_TO_1_00", "1_30_TO_1_60", "1_60_TO_2_00", "2_00_TO_2_20"};
const std::vector<std::pair<double,double>> kQEtaRange = {
    {-2.4,-2.0},{-2.0,-1.6},{-1.6,-1.3},{-0.9,-0.5},{-0.5,-0.1},
    { 0.1, 0.5},{ 0.5, 1.0},{ 1.3, 1.6},{ 1.6, 2.0},{ 2.0, 2.2}};

// DeltaR bins of Step 2
const std::vector<std::string> kDrSuffix = {"dr0_0_2", "dr0_2_1_0", "dr1_0_inf"};
const std::vector<std::string> kDrTex    = {"#DeltaR < 0.2", "0.2 #leq #DeltaR < 1.0",
                                            "#DeltaR #geq 1.0"};
const std::vector<Color_t>     kDrColor  = {kRed + 1, kBlue + 1, kGreen + 2};
const std::vector<Style_t>     kDrMarker = {20, 21, 22};

const Color_t kMCColor   = kRed + 1;
const Color_t kDataColor = kBlack;

void DrawHeadline(const std::string& text, double x = 0.12, double y = 0.955,
                  double size = 0.038)
{
    TLatex tl;
    tl.SetNDC();
    tl.SetTextSize(size);
    tl.SetTextFont(42);
    tl.DrawLatex(x, y, text.c_str());
}

// weighted mean of ratio-hist bins with centers in [xlo,xhi]; returns {mean, err}
std::pair<double,double> PlateauWeightedMean(TH1* ratio, double xlo, double xhi)
{
    double sumw = 0., sumwv = 0.;
    for (int i = 1; i <= ratio->GetNbinsX(); ++i) {
        const double xc = ratio->GetBinCenter(i);
        if (xc < xlo || xc > xhi) continue;
        const double v = ratio->GetBinContent(i);
        const double e = ratio->GetBinError(i);
        if (e <= 0. || v == 0.) continue;
        const double w = 1. / (e * e);
        sumw  += w;
        sumwv += w * v;
    }
    if (sumw <= 0.) throw std::runtime_error("PlateauWeightedMean: no usable bins in window");
    return {sumwv / sumw, std::sqrt(1. / sumw)};
}

} // namespace

// ================================================================= main

void plot_mc_trig_eff(const std::string& sample = "pp", bool use_tight_wp = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gErrorIgnoreLevel = kWarning;

    const SampleCfg cfg = MakeCfg(sample);

    // WP config (registry: Analysis/docs/muon_wp_registry.md): Tight nominal unsuffixed;
    // Medium MC inputs carry _medium_wp. Data tag-and-probe files exist for Tight only.
    const std::string wp_suf  = use_tight_wp ? "" : "_medium_wp";
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";
    const std::string headline = cfg.sample_text + ", " + wp_text;

    TFile* fmc   = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf + ".root");
    TFile* fmc3  = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf + "_step3.root");
    TFile* ffit  = OpenFile(cfg.mc_dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root");
    TFile* fdata = OpenFile(cfg.data_file);

    const std::string dir1 = cfg.out_base + "step1_singles_data_mc/";
    const std::string dir2 = cfg.out_base + "step2_dr_binned_singles/";
    const std::string dir3 = cfg.out_base + "step3_dr_correction/";
    for (const auto& d : {dir1, dir2, dir3}) gSystem->mkdir(d.c_str(), kTRUE);

    const std::string mc_leg   = "MC direct P(mu4 | reco #mu)";
    const std::string data_leg = cfg.data_text + " T&P P(2mu4 | mu4 tag, #DeltaR>0.8)";

    // ================================================================
    // Step 1 (§3.1): data vs MC singles efficiency
    // ================================================================
    std::cout << "\n===== Step 1 (" << sample << ", " << wp_text << ") =====\n";

    struct Var { std::string mc, data, xtitle; double xlo, xhi; bool logx; };
    const std::vector<Var> vars = {
        {"pt",  "pt2nd",  "p_{T} [GeV]", 4.0, 60.0, true},
        {"eta", "eta2nd", "#eta",       -2.4,  2.4, false},
        {"phi", "phi2nd", "#phi",       -M_PI, M_PI, false}};

    for (const auto& v : vars) {
        TCanvas c(("c_step1_" + v.mc).c_str(), "", 1400, 600);
        c.Divide(2, 1);
        for (int ic = 0; ic < 2; ++ic) {
            c.cd(ic + 1);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            if (v.logx) gPad->SetLogx();

            TH1D* mnum = GetObj<TH1D>(fmc, "h_mc_" + v.mc + "_num_"   + kCharges[ic]);
            TH1D* mden = GetObj<TH1D>(fmc, "h_mc_" + v.mc + "_denom_" + kCharges[ic]);
            TH1D* dnum = GetObj<TH1D>(fdata, "h_" + v.data + cfg.ctr + "_" + kDataSigns[ic] + "_2mu4_sepr");
            TH1D* dden = GetObj<TH1D>(fdata, "h_" + v.data + cfg.ctr + "_" + kDataSigns[ic] + "_mu4_sepr");

            auto* gmc = BayesEff(mnum, mden);
            auto* gda = BayesEff(dnum, dden);
            StyleGraph(gmc, kMCColor, 21);
            StyleGraph(gda, kDataColor, 20);

            DrawEffFrame(v.xlo, v.xhi, v.xtitle);
            DrawUnityLine(v.xlo, v.xhi);
            gda->Draw("PZ same");
            gmc->Draw("PZ same");

            DrawHeadline(headline + ", " + kChargeTex[ic]);
            // wide box + smaller text: the overlay data label ("... T&P P(2mu4 | mu4 tag,
            // DR>0.8)") is long and must not clip at the pad edge (plot review iter 1)
            auto* leg = new TLegend(0.18, 0.16, 0.93, 0.32);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.027);
            leg->AddEntry(gmc, mc_leg.c_str(), "lp");
            leg->AddEntry(gda, data_leg.c_str(), "lp");
            leg->Draw();

            // record MC/data ratio magnitudes (turn-on + plateau) from the pt hists
            if (v.mc == "pt") {
                std::cout << "  " << kCharges[ic] << " MC/data eff ratio vs pT:\n";
                for (double pt : {4.3, 5.5, 7.0, 10.0, 20.0, 40.0}) {
                    const int bm = mden->FindBin(pt);
                    const int bd = dden->FindBin(pt);
                    const double em = mden->GetBinContent(bm) > 0
                        ? mnum->GetBinContent(bm) / mden->GetBinContent(bm) : 0.;
                    const double ed = dden->GetBinContent(bd) > 0
                        ? dnum->GetBinContent(bd) / dden->GetBinContent(bd) : 0.;
                    printf("    pT=%5.1f  MC=%.4f  data=%.4f  MC/data=%.3f\n",
                           pt, em, ed, ed > 0 ? em / ed : 0.);
                }
            }
        }
        SaveCanvas(c, dir1 + "step1_eff_" + v.mc + wp_suf + ".png");
    }

    // --- pT in q.eta bins: per charge, 3x4 grid (10 bins + legend pad) ---
    for (int ic = 0; ic < 2; ++ic) {
        TCanvas c(("c_step1_qeta_" + kCharges[ic]).c_str(), "", 1500, 1800);
        c.Divide(3, 4);  // ncols=3, nrows=4 (nrows >= ncols)
        for (size_t iq = 0; iq < kQEtaSuffix.size(); ++iq) {
            c.cd(static_cast<int>(iq) + 1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.12);
            gPad->SetLogx();

            auto* gmc = GetObj<TGraphAsymmErrors>(ffit,
                "g_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq]);
            auto* fmy = GetObj<TF1>(ffit,
                "f_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq]);
            auto* gda = GetObj<TGraphAsymmErrors>(fdata,
                "g_pt2nd_vs_q_eta2nd" + cfg.ctr + "_" + kDataSigns[ic] +
                "_2mu4_sepr_py_" + kQEtaSuffix[iq] + "_divided");

            gmc = (TGraphAsymmErrors*)gmc->Clone();
            gda = (TGraphAsymmErrors*)gda->Clone();
            StyleGraph(gmc, kMCColor, 21, 0.7);
            StyleGraph(gda, kDataColor, 20, 0.7);

            DrawEffFrame(4.0, 60.0, "p_{T} [GeV]");
            DrawUnityLine(4.0, 60.0);
            auto* fdraw = (TF1*)fmy->Clone();
            fdraw->SetLineColor(kMCColor);
            fdraw->SetLineWidth(1);
            fdraw->Draw("same");
            gda->Draw("PZ same");
            gmc->Draw("PZ same");

            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize(0.055);
            tl.SetTextFont(42);
            tl.DrawLatex(0.35, 0.24, Form("%.1f < q#upoint#eta < %.1f",
                                          kQEtaRange[iq].first, kQEtaRange[iq].second));
        }
        // legend / label pad
        c.cd(11);
        DrawHeadline(headline + ", " + kChargeTex[ic], 0.02, 0.88, 0.048);
        auto* leg = new TLegend(0.02, 0.45, 0.98, 0.80);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.05);
        auto* gm = new TGraphAsymmErrors(); StyleGraph(gm, kMCColor, 21);
        auto* gd = new TGraphAsymmErrors(); StyleGraph(gd, kDataColor, 20);
        leg->AddEntry(gm, (mc_leg + " (+ fit)").c_str(), "lp");
        leg->AddEntry(gd, (cfg.data_text + " tag&probe").c_str(), "lp");
        leg->Draw();
        TLatex note;
        note.SetNDC();
        note.SetTextSize(0.045);
        note.SetTextFont(42);
        note.DrawLatex(0.02, 0.30, "Data: T&P P(2mu4 | mu4 tag, #DeltaR>0.8 pairs);");
        note.DrawLatex(0.02, 0.22, "MC: direct conditional, no T&P");
        SaveCanvas(c, dir1 + "step1_eff_pt_in_q_eta_bins_" + kCharges[ic] + wp_suf + ".png");
    }

    // ================================================================
    // Step 1-SF variant: data vs MC vs MC x SF (reco/ID WP scale factors)
    // SF from the skim MuonEfficiencyScaleFactors tools, WP-matched; the SF
    // maps start at pT = 5 GeV, below which the NTP sets SF = 1 (uncorrected).
    // ================================================================
    if (fmc->Get(("h_mc_pt_num_sf_muplus")) != nullptr) {
        const std::string dir1sf = cfg.out_base + "step1_singles_data_mc_sf/";
        gSystem->mkdir(dir1sf.c_str(), kTRUE);
        const Color_t kSFColor = kBlue + 1;
        const std::string sf_leg = "MC #times SF (reco/ID " +
            std::string(use_tight_wp ? "Tight" : "Medium") + "-WP scale factor)";
        const char* sf_note = "SF maps start at p_{T} = 5 GeV; below, SF = 1 (uncorrected)";

        std::cout << "\n===== Step 1-SF (" << sample << ", " << wp_text << ") =====\n";

        for (const auto& v : vars) {
            TCanvas c(("c_step1sf_" + v.mc).c_str(), "", 1400, 600);
            c.Divide(2, 1);
            for (int ic = 0; ic < 2; ++ic) {
                c.cd(ic + 1);
                gPad->SetLeftMargin(0.12);
                gPad->SetBottomMargin(0.12);
                if (v.logx) gPad->SetLogx();

                TH1D* mnum  = GetObj<TH1D>(fmc, "h_mc_" + v.mc + "_num_"    + kCharges[ic]);
                TH1D* msfnum= GetObj<TH1D>(fmc, "h_mc_" + v.mc + "_num_sf_" + kCharges[ic]);
                TH1D* mden  = GetObj<TH1D>(fmc, "h_mc_" + v.mc + "_denom_"  + kCharges[ic]);
                TH1D* dnum  = GetObj<TH1D>(fdata, "h_" + v.data + cfg.ctr + "_" + kDataSigns[ic] + "_2mu4_sepr");
                TH1D* dden  = GetObj<TH1D>(fdata, "h_" + v.data + cfg.ctr + "_" + kDataSigns[ic] + "_mu4_sepr");

                auto* gmc = BayesEff(mnum, mden);
                auto* gda = BayesEff(dnum, dden);
                // SF-weighted numerator (weights < 1): TH1::Divide with error propagation,
                // not Bayes (same reasoning as the step-3 inverse-weighted ratio)
                auto* gsf = DivideGraphClean(msfnum, mden);
                StyleGraph(gmc, kMCColor, 21);
                StyleGraph(gsf, kSFColor, 22);
                StyleGraph(gda, kDataColor, 20);

                DrawEffFrame(v.xlo, v.xhi, v.xtitle);
                DrawUnityLine(v.xlo, v.xhi);
                gda->Draw("PZ same");
                gmc->Draw("PZ same");
                gsf->Draw("PZ same");

                DrawHeadline(headline + ", " + kChargeTex[ic]);
                auto* leg = new TLegend(0.18, 0.14, 0.93, 0.34);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.025);
                leg->AddEntry(gmc, mc_leg.c_str(), "lp");
                leg->AddEntry(gsf, sf_leg.c_str(), "lp");
                leg->AddEntry(gda, data_leg.c_str(), "lp");
                leg->Draw();
                TLatex tn; tn.SetNDC(); tn.SetTextSize(0.022); tn.SetTextFont(42);
                tn.DrawLatex(0.18, 0.355, sf_note);

                if (v.mc == "pt") {
                    std::cout << "  " << kCharges[ic] << " (MC x SF)/data eff ratio vs pT:\n";
                    for (double pt : {4.3, 5.5, 7.0, 10.0, 20.0, 40.0}) {
                        const int bm = mden->FindBin(pt);
                        const int bd = dden->FindBin(pt);
                        const double es = mden->GetBinContent(bm) > 0
                            ? msfnum->GetBinContent(bm) / mden->GetBinContent(bm) : 0.;
                        const double ed = dden->GetBinContent(bd) > 0
                            ? dnum->GetBinContent(bd) / dden->GetBinContent(bd) : 0.;
                        printf("    pT=%5.1f  MCxSF=%.4f  data=%.4f  ratio=%.3f\n",
                               pt, es, ed, ed > 0 ? es / ed : 0.);
                    }
                }
            }
            SaveCanvas(c, dir1sf + "step1_eff_" + v.mc + "_sf" + wp_suf + ".png");
        }

        // --- pT in q.eta bins, SF variant ---
        for (int ic = 0; ic < 2; ++ic) {
            TH2D* h2num  = GetObj<TH2D>(fmc, "h_mc_pt_vs_q_eta_num_"    + kCharges[ic]);
            TH2D* h2sf   = GetObj<TH2D>(fmc, "h_mc_pt_vs_q_eta_num_sf_" + kCharges[ic]);
            TH2D* h2den  = GetObj<TH2D>(fmc, "h_mc_pt_vs_q_eta_denom_"  + kCharges[ic]);
            TCanvas c(("c_step1sf_qeta_" + kCharges[ic]).c_str(), "", 1500, 1800);
            c.Divide(3, 4);
            for (size_t iq = 0; iq < kQEtaSuffix.size(); ++iq) {
                c.cd(static_cast<int>(iq) + 1);
                gPad->SetLeftMargin(0.13);
                gPad->SetBottomMargin(0.12);
                gPad->SetLogx();

                const int b1 = h2den->GetXaxis()->FindBin(kQEtaRange[iq].first  + 1e-6);
                const int b2 = h2den->GetXaxis()->FindBin(kQEtaRange[iq].second - 1e-6);
                TH1D* pn  = h2num->ProjectionY(Form("pn_%s_%zu",  kCharges[ic].c_str(), iq), b1, b2, "e");
                TH1D* ps  = h2sf ->ProjectionY(Form("ps_%s_%zu",  kCharges[ic].c_str(), iq), b1, b2, "e");
                TH1D* pd  = h2den->ProjectionY(Form("pd_%s_%zu",  kCharges[ic].c_str(), iq), b1, b2, "e");

                auto* gmc = BayesEff(pn, pd);
                auto* gsf = DivideGraphClean(ps, pd);
                auto* gda = GetObj<TGraphAsymmErrors>(fdata,
                    "g_pt2nd_vs_q_eta2nd" + cfg.ctr + "_" + kDataSigns[ic] +
                    "_2mu4_sepr_py_" + kQEtaSuffix[iq] + "_divided");
                gda = (TGraphAsymmErrors*)gda->Clone();
                StyleGraph(gmc, kMCColor, 21, 0.7);
                StyleGraph(gsf, kSFColor, 22, 0.7);
                StyleGraph(gda, kDataColor, 20, 0.7);

                DrawEffFrame(4.0, 60.0, "p_{T} [GeV]");
                DrawUnityLine(4.0, 60.0);
                gda->Draw("PZ same");
                gmc->Draw("PZ same");
                gsf->Draw("PZ same");

                TLatex tl;
                tl.SetNDC(); tl.SetTextSize(0.055); tl.SetTextFont(42);
                tl.DrawLatex(0.35, 0.24, Form("%.1f < q#upoint#eta < %.1f",
                                              kQEtaRange[iq].first, kQEtaRange[iq].second));
            }
            c.cd(11);
            DrawHeadline(headline + ", " + kChargeTex[ic], 0.02, 0.88, 0.048);
            auto* leg = new TLegend(0.02, 0.42, 0.98, 0.80);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.045);
            auto* gm = new TGraphAsymmErrors(); StyleGraph(gm, kMCColor, 21);
            auto* gs = new TGraphAsymmErrors(); StyleGraph(gs, kSFColor, 22);
            auto* gd = new TGraphAsymmErrors(); StyleGraph(gd, kDataColor, 20);
            leg->AddEntry(gm, mc_leg.c_str(), "lp");
            leg->AddEntry(gs, sf_leg.c_str(), "lp");
            leg->AddEntry(gd, (cfg.data_text + " tag&probe").c_str(), "lp");
            leg->Draw();
            TLatex note;
            note.SetNDC(); note.SetTextSize(0.04); note.SetTextFont(42);
            note.DrawLatex(0.02, 0.30, sf_note);
            SaveCanvas(c, dir1sf + "step1_eff_pt_in_q_eta_bins_" + kCharges[ic] + "_sf" + wp_suf + ".png");
        }
    } else {
        std::cout << "plot_mc_trig_eff: no h_mc_*_num_sf_* hists in MC file - skipping Step 1-SF"
                  << " (rerun FillMCTrigEffHists on SF-enabled NTP output)" << std::endl;
    }

    // ================================================================
    // Step 2 (§3.2): DeltaR-binned MC singles efficiency
    // ================================================================
    std::cout << "\n===== Step 2 (" << sample << ", " << wp_text << ") =====\n";

    const std::string leg_incl = "inclusive singles (Step 1)";

    // --- 1D pt and 1D q.eta, 2-pad canvases (mu+ | mu-) -------------
    struct PairVar { std::string tag, xtitle; double xlo, xhi; bool logx; };
    const std::vector<PairVar> pvars = {
        {"pt",    "p_{T} [GeV]",   4.0, 60.0, true},
        {"q_eta", "q#upoint#eta", -2.4,  2.4, false}};

    for (const auto& v : pvars) {
        TCanvas c(("c_step2_" + v.tag).c_str(), "", 1400, 600);
        c.Divide(2, 1);
        for (int ic = 0; ic < 2; ++ic) {
            c.cd(ic + 1);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            if (v.logx) gPad->SetLogx();

            DrawEffFrame(v.xlo, v.xhi, v.xtitle);
            DrawUnityLine(v.xlo, v.xhi);

            // Step-1 inclusive singles reference (thin black). No 1D q.eta singles hist
            // exists -> project the singles 2D (x = q.eta) over all pT.
            TH1* rnum = nullptr;
            TH1* rden = nullptr;
            if (v.tag == "pt") {
                rnum = GetObj<TH1D>(fmc, "h_mc_pt_num_"   + kCharges[ic]);
                rden = GetObj<TH1D>(fmc, "h_mc_pt_denom_" + kCharges[ic]);
            } else {
                TH2D* h2n = GetObj<TH2D>(fmc, "h_mc_pt_vs_q_eta_num_"   + kCharges[ic]);
                TH2D* h2d = GetObj<TH2D>(fmc, "h_mc_pt_vs_q_eta_denom_" + kCharges[ic]);
                rnum = h2n->ProjectionX(Form("px_num_%s_%d", v.tag.c_str(), ic));
                rden = h2d->ProjectionX(Form("px_den_%s_%d", v.tag.c_str(), ic));
            }
            auto* gref = BayesEff(rnum, rden);
            StyleGraph(gref, kBlack, 1, 0.4, 1);
            gref->Draw("LX same");   // thin black reference line, no error bars

            auto* leg = new TLegend(0.38, 0.16, 0.92, 0.38);
            leg->SetBorderSize(0);
            // semi-opaque backing: this legend can sit inside a dense error-bar cloud
            // on the low-stats overlay q.eta panel (plot review iter 1, INFO)
            leg->SetFillColorAlpha(kWhite, 0.75);
            leg->SetFillStyle(1001);
            leg->SetTextSize(0.034);
            for (size_t id = 0; id < kDrSuffix.size(); ++id) {
                TH1D* n = GetObj<TH1D>(fmc, "h_mc_pair_" + v.tag + "_num_"   +
                                            kCharges[ic] + "_" + kDrSuffix[id]);
                TH1D* d = GetObj<TH1D>(fmc, "h_mc_pair_" + v.tag + "_denom_" +
                                            kCharges[ic] + "_" + kDrSuffix[id]);
                auto* g = BayesEff(n, d);
                StyleGraph(g, kDrColor[id], kDrMarker[id], 0.8);
                g->Draw("PZ same");
                leg->AddEntry(g, kDrTex[id].c_str(), "lp");
            }
            leg->AddEntry(gref, leg_incl.c_str(), "l");
            leg->Draw();
            DrawHeadline(headline + ", " + kChargeTex[ic]);
        }
        SaveCanvas(c, dir2 + "step2_eff_" + v.tag + "_dr_bins" + wp_suf + ".png");
    }

    // --- pT in q.eta bins per charge, 3 DeltaR lines per pad --------
    for (int ic = 0; ic < 2; ++ic) {
        TCanvas c(("c_step2_qeta_" + kCharges[ic]).c_str(), "", 1500, 1800);
        c.Divide(3, 4);
        for (size_t iq = 0; iq < kQEtaSuffix.size(); ++iq) {
            c.cd(static_cast<int>(iq) + 1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.12);
            gPad->SetLogx();
            DrawEffFrame(4.0, 60.0, "p_{T} [GeV]");
            DrawUnityLine(4.0, 60.0);
            for (size_t id = 0; id < kDrSuffix.size(); ++id) {
                TH2D* h2n = GetObj<TH2D>(fmc, "h_mc_pair_pt_vs_q_eta_num_" +
                                              kCharges[ic] + "_" + kDrSuffix[id]);
                TH2D* h2d = GetObj<TH2D>(fmc, "h_mc_pair_pt_vs_q_eta_denom_" +
                                              kCharges[ic] + "_" + kDrSuffix[id]);
                const int blo = h2n->GetXaxis()->FindBin(kQEtaRange[iq].first  + 1e-6);
                const int bhi = h2n->GetXaxis()->FindBin(kQEtaRange[iq].second - 1e-6);
                TH1D* n = h2n->ProjectionY(Form("py_n_%d_%zu_%zu", ic, iq, id), blo, bhi);
                TH1D* d = h2d->ProjectionY(Form("py_d_%d_%zu_%zu", ic, iq, id), blo, bhi);
                auto* g = BayesEff(n, d);
                StyleGraph(g, kDrColor[id], kDrMarker[id], 0.7);
                g->Draw("PZ same");
            }
            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize(0.055);
            tl.SetTextFont(42);
            tl.DrawLatex(0.35, 0.24, Form("%.1f < q#upoint#eta < %.1f",
                                          kQEtaRange[iq].first, kQEtaRange[iq].second));
        }
        c.cd(11);
        DrawHeadline(headline + ", " + kChargeTex[ic], 0.02, 0.88, 0.048); // 0.048: long overlay headline + charge must fit (review iter 1)
        auto* leg = new TLegend(0.05, 0.35, 0.95, 0.80);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.06);
        for (size_t id = 0; id < kDrSuffix.size(); ++id) {
            auto* g = new TGraphAsymmErrors();
            StyleGraph(g, kDrColor[id], kDrMarker[id]);
            leg->AddEntry(g, kDrTex[id].c_str(), "lp");
        }
        leg->Draw();
        SaveCanvas(c, dir2 + "step2_eff_pt_in_q_eta_bins_" + kCharges[ic] +
                          "_dr_bins" + wp_suf + ".png");
    }

    // ================================================================
    // Step 3 (§3.3): eps_dR(dR) = inverse-weighted num / denom
    // ================================================================
    std::cout << "\n===== Step 3 (" << sample << ", " << wp_text << ") =====\n";

    // ratios with TH1::Divide error propagation (inverse weights > 1 -> Bayes invalid)
    auto MakeRatio = [&](const std::string& tag) -> TH1D* {
        TH1D* num = GetObj<TH1D>(fmc3, "h_mc_dr_" + tag + "_num");
        TH1D* den = GetObj<TH1D>(fmc3, "h_mc_dr_" + tag + "_denom");
        auto* r = (TH1D*)num->Clone(("r_dr_" + tag).c_str());
        r->SetDirectory(nullptr);
        r->Divide(den);
        return r;
    };
    TH1D* r_zoom = MakeRatio("zoom");
    TH1D* r_full = MakeRatio("full");

    const auto plateau = PlateauWeightedMean(r_full, 1.0, 3.0);
    printf("  %s plateau (weighted mean, dR in [1,3]): %.4f +- %.4f\n",
           sample.c_str(), plateau.first, plateau.second);

    auto DrawStep3 = [&](TH1D* r, double xlo, double xhi, const std::string& png,
                         bool plateau_in_range) {
        TCanvas c(("c_" + png).c_str(), "", 900, 700);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        double ymax = 0.;
        for (int i = 1; i <= r->GetNbinsX(); ++i)
            ymax = std::max(ymax, r->GetBinContent(i) + r->GetBinError(i));
        ymax = std::max(1.15, 1.15 * ymax);   // y linear, auto but include 1
        DrawEffFrame(xlo, xhi, "#DeltaR", 0.0, ymax, cfg.eps_dr_text);
        DrawUnityLine(xlo, xhi);

        // fitted plateau line: solid over the fit window [1,3] where in range,
        // dotted across the pad otherwise (zoom canvas)
        if (plateau_in_range) {
            auto* lp = new TLine(1.0, plateau.first, 3.0, plateau.first);
            lp->SetLineColor(kBlue + 1);
            lp->SetLineWidth(3);
            lp->Draw("same");
        } else {
            auto* lp = new TLine(xlo, plateau.first, xhi, plateau.first);
            lp->SetLineColor(kBlue + 1);
            lp->SetLineWidth(2);
            lp->SetLineStyle(3);
            lp->Draw("same");
        }

        r->SetMarkerStyle(20);
        r->SetMarkerColor(kMCColor);
        r->SetLineColor(kMCColor);
        r->SetLineWidth(2);
        r->Draw("E1 same");

        DrawHeadline(headline);
        TLatex tl;
        tl.SetNDC();
        tl.SetTextFont(42);
        tl.SetTextSize(0.035);
        tl.DrawLatex(0.40, 0.86, Form("plateau #LT#DeltaR#in[1,3]#GT = %.3f #pm %.3f",
                                      plateau.first, plateau.second));
        tl.DrawLatex(0.40, 0.80, (cfg.eps_dr_text +
            " = P(trig | #DeltaR) / (#varepsilon_{1}#varepsilon_{2}), MC #varepsilon in weights").c_str());
        SaveCanvas(c, dir3 + png + wp_suf + ".png");
    };
    DrawStep3(r_zoom, 0.0, 1.0,  "step3_eps_dr_zoom", false);
    DrawStep3(r_full, 0.0, 5.75, "step3_eps_dr_full", true);

    // --- pair-pT-binned zoom ratio (4 slices of the pT_bins_120 axis) ---
    {
        TH2D* h2n = GetObj<TH2D>(fmc3, "h_mc_dr_zoom_vs_pair_pt_num");
        TH2D* h2d = GetObj<TH2D>(fmc3, "h_mc_dr_zoom_vs_pair_pt_denom");
        // slice edges aligned to the pair-pT axis bin edges (F2)
        const std::vector<std::pair<int,int>> ybins = {{1,3},{4,6},{7,9},{10,15}};
        const std::vector<Color_t> scol   = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 2};
        const std::vector<Style_t> smark  = {20, 21, 22, 23};

        TCanvas c("c_step3_ptslices", "", 900, 700);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        std::vector<TH1D*> ratios;
        std::vector<std::string> slabels;
        double ymax = 0.;
        for (size_t is = 0; is < ybins.size(); ++is) {
            TH1D* n = h2n->ProjectionX(Form("s3_n_%zu", is), ybins[is].first, ybins[is].second);
            TH1D* d = h2d->ProjectionX(Form("s3_d_%zu", is), ybins[is].first, ybins[is].second);
            auto* r = (TH1D*)n->Clone(Form("s3_r_%zu", is));
            r->SetDirectory(nullptr);
            r->Divide(d);
            ratios.push_back(r);
            const double plo = h2n->GetYaxis()->GetBinLowEdge(ybins[is].first);
            const double phi = h2n->GetYaxis()->GetBinUpEdge(ybins[is].second);
            slabels.push_back(Form("%.1f < p_{T}^{pair} < %.1f GeV", plo, phi));
            for (int i = 1; i <= r->GetNbinsX(); ++i)
                ymax = std::max(ymax, r->GetBinContent(i) + r->GetBinError(i));
        }
        ymax = std::max(1.15, 1.15 * ymax);
        DrawEffFrame(0.0, 1.0, "#DeltaR", 0.0, ymax, cfg.eps_dr_text);
        DrawUnityLine(0.0, 1.0);

        auto* leg = new TLegend(0.45, 0.68, 0.92, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        for (size_t is = 0; is < ratios.size(); ++is) {
            ratios[is]->SetMarkerStyle(smark[is]);
            ratios[is]->SetMarkerColor(scol[is]);
            ratios[is]->SetLineColor(scol[is]);
            ratios[is]->SetLineWidth(2);
            ratios[is]->Draw("E1 same");
            leg->AddEntry(ratios[is], slabels[is].c_str(), "lp");
        }
        leg->Draw();
        DrawHeadline(headline);
        SaveCanvas(c, dir3 + "step3_eps_dr_zoom_pair_pt_slices" + wp_suf + ".png");
    }

    fmc->Close(); fmc3->Close(); ffit->Close(); fdata->Close();
    std::cout << "\nplot_mc_trig_eff(" << sample << ") done.\n";
}
