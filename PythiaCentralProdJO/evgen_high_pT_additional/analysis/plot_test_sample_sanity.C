// Sanity check + plots of the grid test samples of the additional high-pT / same-sign request.
// Truth-only NTUPs (SkimCode TruthTrigRates.py). Definitions are the GENERATOR FILTER's
// (run/Py8Control_hQCD_dimu_pdf_dilepton_filters.py): muon = stable truth muon with
// pT > 3.7 GeV, |eta| < 2.5; pair = any two such muons. This validates the filter; it is not an
// analysis-level selection (tracking doc additional_mc_request_jo_test.md, Physics Procedure §3c).
//
// Per-slice weight (combined + overlay plots) = sigma * eps_filt / N_accepted, all three read from
// THIS sample's own grid log.generate files (MetaData cross-section, EvgenFilterSeq counts) — no
// other production's numbers are used. Units: nb.
//
// Usage: root -l -b -q 'plot_test_sample_sanity.C+("<base dir>", "<out dir>")'
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <regex>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <algorithm>

namespace {
const double kMuPt = 3.7, kMuEta = 2.5, kMassMax = 10.;   // filter cuts [GeV], JO control file

struct Cfg {
    std::string dir, group, slice;   // group: "SI" / "SS"
    double pair_pt_min;               // pair-pT filter cut [GeV]; < 0 = none
    bool same_sign;
};
// the 6 JOs in run/ (PairPtMin / PairSameSign per JO)
const std::vector<Cfg> kCfgs = {
    {"run_sign_indep_pT70_125",  "SI", "70-125",  45., false},
    {"run_sign_indep_pT125_300", "SI", "125-300", 45., false},
    {"run_same_sign_pT24_40",    "SS", "24-40",   -1., true},
    {"run_same_sign_pT40_70",    "SS", "40-70",   20., true},
    {"run_same_sign_pT70_125",   "SS", "70-125",  20., true},
    {"run_same_sign_pT125_300",  "SS", "125-300", 20., true},
};

std::vector<std::string> FindFiles(const std::string& dir, const std::string& must, const std::string& mustnot = "") {
    std::vector<std::string> out;
    std::string cmd = "find " + dir + " -type f -name '*" + must + "*'";
    FILE* p = gSystem->OpenPipe(cmd.c_str(), "r");
    char buf[4096];
    while (p && fgets(buf, sizeof buf, p)) {
        std::string s(buf); s.erase(s.find_last_not_of("\n") + 1);
        if (!mustnot.empty() && s.find(mustnot) != std::string::npos) continue;
        out.push_back(s);
    }
    if (p) gSystem->ClosePipe(p);
    std::sort(out.begin(), out.end());
    return out;
}

struct GenInfo { double sigma_nb = 0; long acc = 0, gen = 0; int njobs = 0; };

// sigma: mean over jobs of "MetaData: cross-section (nb)= X"; eps = sum(acc)/sum(gen) from
// "Filter Efficiency = e [acc / gen]" (EvgenFilterSeq, unweighted line).
GenInfo ReadLogs(const std::string& dir) {
    GenInfo g; double sum_sigma = 0;
    std::regex rs(R"(MetaData: cross-section \(nb\)=\s*([0-9.eE+-]+))");
    std::regex rf(R"(EvgenFilterSeq\s+INFO Filter Efficiency = \S+ \[(\d+) / (\d+)\])");
    for (auto& f : FindFiles(dir, "log.generate")) {
        std::ifstream in(f); std::string l; bool got_s = false, got_f = false;
        while (std::getline(in, l)) {
            std::smatch m;
            if (std::regex_search(l, m, rs)) { sum_sigma += std::stod(m[1]); got_s = true; }
            if (std::regex_search(l, m, rf)) { g.acc += std::stol(m[1]); g.gen += std::stol(m[2]); got_f = true; }
        }
        if (got_s && got_f) g.njobs++;
        else std::cerr << "WARNING: incomplete log " << f << std::endl;
    }
    if (g.njobs) g.sigma_nb = sum_sigma / g.njobs;
    return g;
}

std::vector<double> LogEdges(double lo, double hi, int n) {
    std::vector<double> e(n + 1);
    for (int i = 0; i <= n; ++i) e[i] = lo * std::pow(hi / lo, double(i) / n);
    return e;
}
} // namespace

void plot_test_sample_sanity(const char* base = "/usatlas/u/yuhanguo/usatlasdata/evgen_high_pT_additional_test",
                             const char* outdir = "/usatlas/u/yuhanguo/usatlasdata/evgen_high_pT_additional_test/plots") {
    gStyle->SetOptStat(0);
    gSystem->mkdir(outdir, kTRUE);
    const auto pt_edges = LogEdges(9., 150., 30);   // user: pair pT 9-150 GeV, log scale

    struct Hs { TH1D *nmu, *pt, *m, *sign; };
    std::map<std::string, Hs> raw, wgt;
    std::ofstream tab(std::string(outdir) + "/sanity_table.csv");
    tab << "# filter check: muon = stable truth mu pT>3.7 GeV |eta|<2.5; event must have >=2 muons and a pair with "
           "0<m<10 GeV, pT(mumu) > cut, same sign if required\n"
        << "config,njobs_logs,sigma_nb,accepted,generated,eps_filt,ntup_events,events_failing_check,"
           "events_nmu_lt2,weight_per_event_nb\n";

    for (const auto& c : kCfgs) {
        const std::string d = std::string(base) + "/" + c.dir;
        GenInfo g = ReadLogs(d);
        auto files = FindFiles(d, ".root", "log");
        TChain ch("HeavyIonD3PD");
        for (auto& f : files) ch.Add(f.c_str());
        const long n = ch.GetEntries();
        if (n == 0 || g.acc == 0) { std::cerr << "SKIP " << c.dir << ": no NTUP entries or no logs" << std::endl; continue; }
        const double w = g.sigma_nb * (double(g.acc) / g.gen) / n;   // nb per event

        std::vector<float> *mu_pt = nullptr, *mu_eta = nullptr;
        std::vector<float> *p_pt1 = nullptr, *p_eta1 = nullptr, *p_pt2 = nullptr, *p_eta2 = nullptr, *p_m = nullptr, *p_pt = nullptr;
        std::vector<int> *p_ch1 = nullptr, *p_ch2 = nullptr;
        ch.SetBranchAddress("truth_muon_pt", &mu_pt);   ch.SetBranchAddress("truth_muon_eta", &mu_eta);
        ch.SetBranchAddress("truth_mupair_pt1", &p_pt1); ch.SetBranchAddress("truth_mupair_eta1", &p_eta1);
        ch.SetBranchAddress("truth_mupair_ch1", &p_ch1); ch.SetBranchAddress("truth_mupair_pt2", &p_pt2);
        ch.SetBranchAddress("truth_mupair_eta2", &p_eta2); ch.SetBranchAddress("truth_mupair_ch2", &p_ch2);
        ch.SetBranchAddress("truth_mupair_m", &p_m);    ch.SetBranchAddress("truth_mupair_pt", &p_pt);

        for (auto* hs : {&raw, &wgt}) {
            const std::string t = (hs == &raw ? "raw_" : "w_") + c.dir;
            (*hs)[c.dir] = {new TH1D((t + "_nmu").c_str(), "", 8, -0.5, 7.5),
                            new TH1D((t + "_pt").c_str(), "", 30, pt_edges.data()),
                            new TH1D((t + "_m").c_str(), "", 40, 0., 40.),
                            new TH1D((t + "_sign").c_str(), "", 2, 0., 2.)};
            for (TH1D* h : {(*hs)[c.dir].nmu, (*hs)[c.dir].pt, (*hs)[c.dir].m, (*hs)[c.dir].sign}) h->Sumw2();
        }
        long nfail = 0, nlt2 = 0;
        for (long i = 0; i < n; ++i) {
            ch.GetEntry(i);
            int nmu = 0;
            for (size_t k = 0; k < mu_pt->size(); ++k)
                if (std::fabs(mu_pt->at(k)) * 1e-3 > kMuPt && std::fabs(mu_eta->at(k)) < kMuEta) ++nmu;
            bool pass = false;
            for (size_t k = 0; k < p_m->size(); ++k) {
                if (std::fabs(p_pt1->at(k)) * 1e-3 <= kMuPt || std::fabs(p_pt2->at(k)) * 1e-3 <= kMuPt ||
                    std::fabs(p_eta1->at(k)) >= kMuEta || std::fabs(p_eta2->at(k)) >= kMuEta) continue;
                const bool ss = p_ch1->at(k) * p_ch2->at(k) > 0;
                const double m = p_m->at(k) * 1e-3, pt = p_pt->at(k) * 1e-3;
                if (m > 0 && m < kMassMax && pt > c.pair_pt_min && (ss || !c.same_sign)) pass = true;
                for (auto* hs : {&raw, &wgt}) {
                    const double ww = (hs == &raw ? 1. : w);
                    (*hs)[c.dir].pt->Fill(pt, ww);   // outside the drawn ranges -> under/overflow, not drawn
                    (*hs)[c.dir].m->Fill(m, ww);
                    (*hs)[c.dir].sign->Fill(ss ? 0.5 : 1.5, ww);
                }
            }
            raw[c.dir].nmu->Fill(std::min(nmu, 7)); wgt[c.dir].nmu->Fill(std::min(nmu, 7), w);
            if (nmu < 2) ++nlt2;
            if (nmu < 2 || !pass) ++nfail;
        }
        tab << c.dir << "," << g.njobs << "," << g.sigma_nb << "," << g.acc << "," << g.gen << ","
            << double(g.acc) / g.gen << "," << n << "," << nfail << "," << nlt2 << "," << w << "\n";
        std::cout << c.dir << ": " << n << " events, failing filter check " << nfail
                  << ", eps_filt = " << g.acc << "/" << g.gen << " = " << double(g.acc) / g.gen << std::endl;
    }
    tab.close();

    const int colors[] = {kAzure + 2, kOrange + 7, kGreen + 2, kMagenta + 1};
    for (std::string grp : {"SI", "SS"}) {
        std::vector<Cfg> cs;
        for (auto& c : kCfgs) if (c.group == grp && raw.count(c.dir)) cs.push_back(c);
        if (cs.empty()) continue;
        const std::string glabel = grp == "SI" ? "sign-indep. filter" : "same-sign filter";

        // (1) pT-hat combined, sigma*eps-weighted: 2x2 panels
        TCanvas c1(("c_" + grp).c_str(), "", 1200, 1000); c1.Divide(2, 2);
        TH1D* sum[4] = {nullptr, nullptr, nullptr, nullptr};
        for (auto& c : cs) {
            TH1D* hs[4] = {wgt[c.dir].nmu, wgt[c.dir].pt, wgt[c.dir].m, wgt[c.dir].sign};
            for (int j = 0; j < 4; ++j) {
                if (!sum[j]) sum[j] = (TH1D*)hs[j]->Clone((grp + "_sum" + std::to_string(j)).c_str());
                else sum[j]->Add(hs[j]);
            }
        }
        sum[1]->Scale(1., "width"); sum[2]->Scale(1., "width");
        const char* xt[4] = {"N_{#mu} per event (p_{T}>3.7 GeV, |#eta|<2.5)", "p_{T}^{#mu#mu} [GeV]", "m_{#mu#mu} [GeV]", ""};
        const char* yt[4] = {"#sigma #times #varepsilon_{filt} per bin [nb]", "d#sigma/dp_{T}^{#mu#mu} [nb/GeV] (pairs)",
                             "d#sigma/dm_{#mu#mu} [nb/GeV] (pairs)", "#sigma #times #varepsilon_{filt} [nb] (pairs)"};
        for (int j = 0; j < 4; ++j) {
            c1.cd(j + 1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13);
            sum[j]->SetLineColor(kBlack); sum[j]->SetMarkerStyle(20); sum[j]->SetMarkerSize(0.7);
            sum[j]->GetXaxis()->SetTitle(xt[j]); sum[j]->GetYaxis()->SetTitle(yt[j]);
            sum[j]->GetYaxis()->SetTitleOffset(1.5);
            if (j == 1) { gPad->SetLogx(); gPad->SetLogy(); sum[j]->GetXaxis()->SetMoreLogLabels(); sum[j]->GetXaxis()->SetNoExponent(); sum[j]->GetXaxis()->SetTitleOffset(1.2); }  // pair pT: log-binned axis, power-law spectrum
            if (j == 3) { sum[j]->GetXaxis()->SetBinLabel(1, "same sign"); sum[j]->GetXaxis()->SetBinLabel(2, "opposite sign");
                          sum[j]->GetXaxis()->SetLabelSize(0.06); }
            if (j != 1) sum[j]->SetMinimum(0);
            sum[j]->Draw(j == 3 ? "HIST E" : "E1");
            if (j == 3) sum[j]->Draw("HIST SAME");
            TLatex t; t.SetNDC(); t.SetTextSize(0.045);
            t.DrawLatex(0.17, 0.93, ("Pythia8 pp #sqrt{s}=5.36 TeV, " + glabel).c_str());
        }
        // pair-pT filter cut (per slice; SS 24-40 has none)
        bool common_cut = true;   // drawn only when every slice of the panel has the same cut
        for (auto& c : cs) common_cut &= (c.pair_pt_min > 0 && c.pair_pt_min == cs.front().pair_pt_min);
        if (common_cut) { c1.cd(2); double ylo = 1e30; for (int b = 1; b <= sum[1]->GetNbinsX(); ++b) if (sum[1]->GetBinContent(b) > 0) ylo = std::min(ylo, sum[1]->GetBinContent(b));
                          TLine l(cs.front().pair_pt_min, 0.5 * ylo,
                                             cs.front().pair_pt_min, sum[1]->GetBinContent(sum[1]->GetMaximumBin()));
                          l.SetLineStyle(2); l.SetLineColor(kRed); l.DrawClone(); }
        // filter cuts drawn on the pair panels
        c1.cd(3); { TLine l(kMassMax, 0, kMassMax, sum[2]->GetBinContent(sum[2]->GetMaximumBin())); l.SetLineStyle(2); l.SetLineColor(kRed); l.DrawClone(); }
        c1.SaveAs((std::string(outdir) + "/" + grp + "_combined_sanity.png").c_str());

        // (2) pair pT per pT-hat slice overlaid, sigma*eps-weighted
        TCanvas c2(("c2_" + grp).c_str(), "", 800, 650);
        c2.SetLogx(); c2.SetLogy(); c2.SetLeftMargin(0.14); c2.SetBottomMargin(0.12);   // log x: log-binned; log y: power-law spectrum
        TLegend lg(0.62, 0.68, 0.9, 0.88); lg.SetBorderSize(0);
        lg.SetHeader("#hat{p}_{T} slice [GeV]");
        double ymax = 0, ymin = 1e30; std::vector<TH1D*> hh;
        for (size_t k = 0; k < cs.size(); ++k) {
            TH1D* h = (TH1D*)wgt[cs[k].dir].pt->Clone((grp + "_ov" + std::to_string(k)).c_str());
            h->Scale(1., "width"); hh.push_back(h);
            ymax = std::max(ymax, h->GetBinContent(h->GetMaximumBin()));
            for (int b = 1; b <= h->GetNbinsX(); ++b) if (h->GetBinContent(b) > 0) ymin = std::min(ymin, h->GetBinContent(b));
        }
        for (size_t k = 0; k < hh.size(); ++k) {
            hh[k]->SetLineColor(colors[k]); hh[k]->SetMarkerColor(colors[k]); hh[k]->SetMarkerStyle(20 + k); hh[k]->SetMarkerSize(0.8);
            hh[k]->SetMaximum(ymax * 5); hh[k]->SetMinimum(ymin * 0.3);
            hh[k]->GetXaxis()->SetTitle("p_{T}^{#mu#mu} [GeV]"); hh[k]->GetXaxis()->SetMoreLogLabels(); hh[k]->GetXaxis()->SetNoExponent(); hh[k]->GetXaxis()->SetTitleOffset(1.2); hh[k]->GetYaxis()->SetTitle("d#sigma/dp_{T}^{#mu#mu} [nb/GeV] (pairs)");
            hh[k]->Draw(k ? "E1 SAME" : "E1");
            std::string lab = cs[k].slice + (cs[k].pair_pt_min < 0 ? " (no p_{T}^{#mu#mu} cut)" : "");
            lg.AddEntry(hh[k], lab.c_str(), "lep");
        }
        lg.Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.035);
        t.DrawLatex(0.16, 0.92, ("Pythia8 pp #sqrt{s}=5.36 TeV, " + glabel).c_str());
        c2.SaveAs((std::string(outdir) + "/" + grp + "_pair_pt_by_pthat.png").c_str());
    }
}
