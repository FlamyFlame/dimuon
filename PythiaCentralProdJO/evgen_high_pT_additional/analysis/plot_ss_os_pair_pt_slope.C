// Helper plot (tracking doc additional_mc_request_jo_test.md, follow-up 4): does the same-sign (SS)
// low-mass pair-pT spectrum fall more steeply than the opposite-sign (OS) one, per pT-hat slice?
// dsigma/dpT^mumu for 45-150 GeV, pTH70_125 and pTH125_300, OS and SS (4 curves) + SS/OS ratio.
//
// Inputs = ntuple-processing output (PythiaTruthAnalysis external raw-counts mode, run_truth_processing.sh):
//   OS pairs from the sign-INDEPENDENT test samples (filter: any pair, m<10, pT^mumu>45)  -> unbiased above 45 GeV
//   SS pairs from the SAME-SIGN test samples        (filter: SS pair, m<10, pT^mumu>20)  -> unbiased above 20 GeV
// (any pair passing the cuts below also passes its sample's filter, whose legs are pT>3.7, |eta|<2.5).
// Pair = truth pair, 0 < m < 10 GeV, truth fiducial of both legs (MCTrigEffPairSel::TruthFiducial: pT>4.5,
// |eta|<2.4). Weight per pair = sigma * eps_filt / N_events, sigma and eps_filt from THIS sample's own grid
// log.generate files (D7), N_events from the processing meta tree. Units nb.
// Axis: the pair-pT edges of plot_test_sample_sanity.C above 45 GeV (13 log-uniform bins 45-150).
//
// Usage: root -l -b -q 'plot_ss_os_pair_pt_slope.C+'
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/Utilities/MCTrigEffPairSelection.h"
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/Utilities/MCRequestLooseMassFilter.h"
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/Utilities/PairTrigEffEvaluator.h"

namespace {
const std::string kBase = "/usatlas/u/yuhanguo/usatlasdata/evgen_high_pT_additional_test";

std::vector<std::string> FindFiles(const std::string& dir, const std::string& must) {
    std::vector<std::string> out;
    FILE* p = gSystem->OpenPipe(("find " + dir + " -type f -name '*" + must + "*' | sort").c_str(), "r");
    char buf[4096];
    while (p && fgets(buf, sizeof buf, p)) { std::string s(buf); s.erase(s.find_last_not_of("\n") + 1); out.push_back(s); }
    if (p) gSystem->ClosePipe(p);
    return out;
}
// sigma*eps_filt [nb] of one config: mean sigma over its grid jobs x sum(acc)/sum(gen) (as plot_test_sample_sanity.C)
double SigmaEps(const std::string& dir) {
    std::regex rs(R"(MetaData: cross-section \(nb\)=\s*([0-9.eE+-]+))");
    std::regex rf(R"(EvgenFilterSeq\s+INFO Filter Efficiency = \S+ \[(\d+) / (\d+)\])");
    double sum_sigma = 0; long acc = 0, gen = 0; int nj = 0;
    for (auto& f : FindFiles(dir, "log.generate")) {
        std::ifstream in(f); std::string l; bool gs = false, gf = false;
        while (std::getline(in, l)) {
            std::smatch m;
            if (std::regex_search(l, m, rs)) { sum_sigma += std::stod(m[1]); gs = true; }
            if (std::regex_search(l, m, rf)) { acc += std::stol(m[1]); gen += std::stol(m[2]); gf = true; }
        }
        if (!(gs && gf)) throw std::runtime_error("incomplete log " + f);
        ++nj;
    }
    if (nj == 0 || gen == 0) throw std::runtime_error("no log.generate under " + dir);
    return sum_sigma / nj * double(acc) / gen;
}
Long64_t NProc(const std::string& file, int kn) {
    TFile f(file.c_str(), "READ");
    TTree* t = (TTree*)f.Get("meta_tree_out");
    const std::string b = "nproc_kin" + std::to_string(kn) + "_beam0";
    if (!t || !t->GetBranch(b.c_str())) throw std::runtime_error("no " + b + " in " + file);
    Long64_t n = 0; t->SetBranchAddress(b.c_str(), &n); t->GetEntry(0);
    if (n <= 0) throw std::runtime_error(b + " <= 0 in " + file);
    return n;
}
// histogram -> graph with x shifted by a factor (keeps overlapping curves visible); x errors = bin edges;
// empty bins skipped
TGraphAsymmErrors* Shifted(const TH1D* h, double fx) {
    auto* g = new TGraphAsymmErrors();
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        if (h->GetBinContent(b) <= 0) continue;
        const double lo = h->GetXaxis()->GetBinLowEdge(b), hi = h->GetXaxis()->GetBinUpEdge(b), x = std::sqrt(lo * hi) * fx;
        const int n = g->GetN();
        g->SetPoint(n, x, h->GetBinContent(b));
        g->SetPointError(n, x - lo, hi - x, h->GetBinError(b), h->GetBinError(b));
    }
    g->SetLineColor(h->GetLineColor()); g->SetMarkerColor(h->GetMarkerColor());
    g->SetMarkerStyle(h->GetMarkerStyle()); g->SetMarkerSize(h->GetMarkerSize());
    return g;
}
std::vector<double> LogEdges(double lo, double hi, int n) {
    std::vector<double> e(n + 1);
    for (int i = 0; i <= n; ++i) e[i] = lo * std::pow(hi / lo, double(i) / n);
    return e;
}
} // namespace

void plot_ss_os_pair_pt_slope(const char* outdir = "/usatlas/u/yuhanguo/usatlasdata/evgen_high_pT_additional_test/plots") {
    gStyle->SetOptStat(0);
    const std::vector<double> edges = LogEdges(45., 150., 13);   // = plot_test_sample_sanity.C edges above 45 GeV
    struct Curve { std::string sign, slice, ntup_dir, label; int kn, tree_sign; };
    const std::vector<Curve> cv = {
        {"OS", "70-125",  "run_sign_indep_pT70_125",  "5p36_SI_m10_pT45", 4, 2},
        {"SS", "70-125",  "run_same_sign_pT70_125",   "5p36_SS_m10_pT20", 4, 1},
        {"OS", "125-300", "run_sign_indep_pT125_300", "5p36_SI_m10_pT45", 5, 2},
        {"SS", "125-300", "run_same_sign_pT125_300",  "5p36_SS_m10_pT20", 5, 1},
    };
    const std::string sel = MCTrigEffPairSel::TruthFiducial() + " && truth_minv > 0 && truth_minv < "
                          + std::to_string(MCRequest::kLooseMassMax);
    // stdout summary: SS/OS sigma ratio integrated over the 3 delivered trig-eff pair-pT cells (canonical edges)
    const std::vector<double> pe = PairTrigEff::PairPtEdges();
    const int first = PairTrigEff::FirstDeliveredPtBin();
    const std::vector<double> cells(pe.begin() + first - 1, pe.end());
    std::vector<TH1D*> h, hc;
    for (const auto& c : cv) {
        const std::string f = kBase + "/truth_processing/muon_pairs_pythia_truth_" + c.label + "_kn" + std::to_string(c.kn)
                            + "_no_data_resonance_cuts.root";
        const double w = SigmaEps(kBase + "/" + c.ntup_dir) / NProc(f, c.kn);
        ROOT::RDataFrame d("muon_pair_tree_kin" + std::to_string(c.kn) + "_sign" + std::to_string(c.tree_sign), f);
        auto hh = d.Alias("m1_truth_pt", "m1.truth_pt").Alias("m1_truth_eta", "m1.truth_eta")
                   .Alias("m2_truth_pt", "m2.truth_pt").Alias("m2_truth_eta", "m2.truth_eta")
                   .Filter(sel)
                   .Histo1D({(c.sign + c.slice).c_str(), "", int(edges.size()) - 1, edges.data()}, "truth_pair_pt");
        auto hcell = d.Alias("m1_truth_pt", "m1.truth_pt").Alias("m1_truth_eta", "m1.truth_eta")
                      .Alias("m2_truth_pt", "m2.truth_pt").Alias("m2_truth_eta", "m2.truth_eta")
                      .Filter(sel).Histo1D({(c.sign + c.slice + "c").c_str(), "", int(cells.size()) - 1, cells.data()}, "truth_pair_pt");
        TH1D* xc = (TH1D*)hcell->Clone((c.sign + "_" + c.slice + "_cells").c_str()); xc->SetDirectory(nullptr); xc->Sumw2(); xc->Scale(w);
        hc.push_back(xc);
        TH1D* x = (TH1D*)hh->Clone((c.sign + "_" + c.slice).c_str()); x->SetDirectory(nullptr);
        std::cout << c.sign << " pTH" << c.slice << ": " << x->Integral() << " pairs in 45-150 GeV, weight " << w << " nb/pair" << std::endl;
        x->Sumw2(); x->Scale(w, "width");
        h.push_back(x);
    }

    for (int k = 0; k < 2; ++k) {
        std::cout << "SS/OS sigma ratio, pTH" << cv[2 * k].slice << ":";
        for (int b = 1; b <= hc[2 * k]->GetNbinsX(); ++b) {
            const double ss = hc[2 * k + 1]->GetBinContent(b), os = hc[2 * k]->GetBinContent(b);
            const double e = (ss > 0 && os > 0) ? ss / os * std::hypot(hc[2 * k + 1]->GetBinError(b) / ss, hc[2 * k]->GetBinError(b) / os) : 0;
            printf("  [%.1f,%.1f) %.4f +- %.4f", cells[b - 1], cells[b], os > 0 ? ss / os : 0., e);
        }
        std::cout << std::endl;
    }
    TCanvas cnv("cnv", "", 800, 850);
    TPad top("top", "", 0, 0.33, 1, 1), bot("bot", "", 0, 0, 1, 0.33);
    top.SetBottomMargin(0.02); top.SetLeftMargin(0.14); bot.SetTopMargin(0.03); bot.SetBottomMargin(0.30); bot.SetLeftMargin(0.14);
    top.Draw(); bot.Draw();
    top.cd(); gPad->SetLogx(); gPad->SetLogy();   // log x: log-binned pair pT; log y: power-law spectrum
    const int col[2] = {kAzure + 2, kOrange + 7};
    double ymax = 0, ymin = 1e30;
    for (auto* x : h) for (int b = 1; b <= x->GetNbinsX(); ++b) if (x->GetBinContent(b) > 0) {
        ymax = std::max(ymax, x->GetBinContent(b)); ymin = std::min(ymin, x->GetBinContent(b)); }
    TLegend lg(0.50, 0.62, 0.90, 0.89); lg.SetBorderSize(0);
    for (size_t i = 0; i < h.size(); ++i) {
        const int k = i / 2; const bool ss = cv[i].sign == "SS";
        h[i]->SetLineColor(col[k]); h[i]->SetMarkerColor(col[k]); h[i]->SetMarkerStyle(ss ? 24 : 20); h[i]->SetMarkerSize(0.9);
        h[i]->SetMaximum(ymax * 5); h[i]->SetMinimum(ymin * 0.3);
        h[i]->GetYaxis()->SetTitle("d#sigma/dp_{T}^{#mu#mu} [nb/GeV] (pairs)"); h[i]->GetYaxis()->SetTitleOffset(1.3);
        h[i]->GetXaxis()->SetLabelSize(0);
        if (i == 0) h[i]->Draw("AXIS");
        TGraphAsymmErrors* g = Shifted(h[i], k ? 1.012 : 0.988);   // 70-125 slightly left, 125-300 right
        g->Draw("P SAME");
        lg.AddEntry(g, ((ss ? "same sign, " : "opposite sign, ") + std::string("#hat{p}_{T} ") + cv[i].slice + " GeV").c_str(), "lep");
    }
    lg.Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.035);
    t.DrawLatex(0.16, 0.92, "Pythia8 pp, #sqrt{s} = 5.36 TeV, truth, m_{#mu#mu} < 10 GeV, p_{T}^{#mu} > 4.5 GeV, |#eta^{#mu}| < 2.4");

    bot.cd(); gPad->SetLogx();   // ratio vs the same log-binned pair pT
    std::vector<TH1D*> r;
    for (int k = 0; k < 2; ++k) {
        TH1D* x = (TH1D*)h[2 * k + 1]->Clone(("ratio" + std::to_string(k)).c_str());
        x->Divide(h[2 * k]);
        x->SetMarkerStyle(24); x->SetLineColor(col[k]); x->SetMarkerColor(col[k]);
        r.push_back(x);
    }
    double rmax = 0;
    for (auto* x : r) for (int b = 1; b <= x->GetNbinsX(); ++b) if (x->GetBinContent(b) > 0) rmax = std::max(rmax, x->GetBinContent(b) + x->GetBinError(b));
    r[0]->SetMinimum(0); r[0]->SetMaximum(rmax * 1.15);
    r[0]->GetYaxis()->SetTitle("SS / OS"); r[0]->GetYaxis()->SetNdivisions(505);
    r[0]->GetYaxis()->SetTitleSize(0.09); r[0]->GetYaxis()->SetLabelSize(0.08); r[0]->GetYaxis()->SetTitleOffset(0.6);
    r[0]->GetXaxis()->SetTitle("p_{T}^{#mu#mu} [GeV]"); r[0]->GetXaxis()->SetTitleSize(0.10); r[0]->GetXaxis()->SetLabelSize(0.09);
    r[0]->GetXaxis()->SetMoreLogLabels(); r[0]->GetXaxis()->SetNoExponent(); r[0]->GetXaxis()->SetLabelOffset(0.01);
    r[0]->Draw("AXIS");
    for (int k = 0; k < 2; ++k) Shifted(r[k], k ? 1.012 : 0.988)->Draw("P SAME");
    cnv.SaveAs((std::string(outdir) + "/ss_os_pair_pt_45_150.png").c_str());
}
