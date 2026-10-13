// Fill stage of the Pb+Pb event-selection run / LB / ZDC-timing diagnostics (evidence for the
// DQ experts -- nothing here changes a cut or a run list).  One call per merged skim part; the
// per-part outputs are hadd-ed and drawn by plot_pbpb_evsel_run_lb_diag.cxx.
//
// Event reading and cuts 1-5 mirror the analysis application PbPbExtras::PassEventSel exactly
// (same cuts file, same per-run preamp tree, same strict/non-strict inequalities, same units):
//   cut 1  zdc_tot = (E_A+E_C)/1000 TeV  >  g_ZDC_FCal_cut(FCal)        -> fail
//   cut 2  |t_A| >= cut2_ns || |t_C| >= cut2_ns                         -> fail
//   cut 3  preamp_A >= cut_A && preamp_C >= cut_C (per-run if present)  -> fail
//   cut 4  N_trk^HItight/N_trk^total < g_ntrk_frac_cut_lo(N_trk^total) -> fail (N_trk^total>0)
//   cut 5  N_trk^HItight outside [g_ntrk_fcal_cut_lo, _hi](FCal)       -> fail
// N_trk^total = trk_numqual[0], N_trk^HItight = trk_numqual[3] (both pT > 400 MeV);
// FCal = (FCal_Et_P + FCal_Et_N) x 1e-6 TeV; ZDC index [1]=A, [0]=C.  HLT_mu4_L1MU3V required.
// Vertex: vtx_z[0] = PV0, the vertex the skim's |z| < MaxZvtx (250 mm) requirement acts on.
//
// Contents of the output file:
//   h_frac_c13_<run>      N_trk^HItight/N_trk^total vs N_trk^total, after cuts 1-3 (= the
//                         sample of event_sel_cut4_nTrk_frac_standalone), axes as that figure
//   h_ntrkfcal_c14_<run>  N_trk^HItight vs FCal, after cuts 1-4 (= cut-5 standalone sample/axes)
//   h_ntrkfcal_nocut_<run>, h_ratio_nocut_<run>  no cut (as plot_pbpb_ntrk_lowband_runs)
//   h_cnt_<run>           counters: 1 HLT, 2 pass 1-3, 3 fail 4 (after 1-3), 4 pass 1-4, 5 fail 5 (after 1-4)
//   h3_zdct_vtxz          t_A x t_C x z(PV0), HLT events, no cut (= ZDC_t_A_vs_C sample/axes)
//   t_c4fail              every event passing cuts 1-3 and failing cut 4 (run, lbn, bcid, ...)
//   t_lbscan              every HLT event of the runs in lbscan_runs, with per-cut pass bits
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TTree.h"

#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

void fill_pbpb_evsel_run_lb_diag(int run_year, int part, std::vector<int> lbscan_runs, std::string out_file,
                                Long64_t nmax = -1)
{
    const int yr = run_year % 2000;
    const std::string yrs  = std::to_string(yr);
    const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";

    // ---- cuts (as PbPbExtras::InitEventSel) ----
    const std::string cuts_path = PbPbEvSelCutsPath(run_year);
    TFile* fc = TFile::Open(cuts_path.c_str(), "READ");
    if (!fc || fc->IsZombie()) throw std::runtime_error("cannot open " + cuts_path);
    auto get_graph = [&](const char* k) {
        TGraph* g = dynamic_cast<TGraph*>(fc->Get(k));
        if (!g) throw std::runtime_error(std::string("cuts file lacks ") + k);
        return g;
    };
    auto get_param = [&](const char* k) {
        auto* p = dynamic_cast<TParameter<double>*>(fc->Get(k));
        if (!p) throw std::runtime_error(std::string("cuts file lacks ") + k);
        return p->GetVal();
    };
    TGraph* g1   = get_graph(PbPbEvSelKey::kZDCFCalCut);
    TGraph* g4   = get_graph(PbPbEvSelKey::kNTrkFracCutLo);
    TGraph* g5lo = get_graph(PbPbEvSelKey::kNTrkFCalCutLo);
    TGraph* g5hi = get_graph(PbPbEvSelKey::kNTrkFCalCutHi);
    const double cut2_ns = get_param(PbPbEvSelKey::kZDCTimeCutNs);
    const float  cut3_A0 = (float)get_param(PbPbEvSelKey::kPreampACutADC);
    const float  cut3_C0 = (float)get_param(PbPbEvSelKey::kPreampCCutADC);
    std::map<int, std::pair<float, float>> cut3_run;
    if (TTree* t = dynamic_cast<TTree*>(fc->Get(PbPbEvSelKey::kPreampPerRunTree))) {
        Int_t r = 0; Double_t a = 0., c = 0.;
        t->SetBranchAddress("run_number", &r); t->SetBranchAddress("cut_A_ADC", &a); t->SetBranchAddress("cut_C_ADC", &c);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); cut3_run[r] = {(float)a, (float)c}; }
    }
    std::cout << "cuts: " << cuts_path << "  cut2 " << cut2_ns << " ns, cut3 scalar " << cut3_A0 << "/" << cut3_C0
              << ", per-run preamp cuts for " << cut3_run.size() << " runs" << std::endl;
    const std::unordered_set<int> bad_runs = PbPbBadRuns(run_year);

    // ---- input ----
    const std::string fin = base + "pbpb_20" + yrs + "/data_pbpb" + yrs + "_part" + std::to_string(part) + ".root";
    if (gSystem->AccessPathName(fin.c_str())) throw std::runtime_error("missing " + fin);
    TChain chain("HeavyIonD3PD");
    chain.Add(fin.c_str());
    chain.SetMakeClass(1);
    chain.SetBranchStatus("*", 0);
    UInt_t run = 0, lbn = 0, bcid = 0;
    Int_t b_HLT = 0;
    Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f, zdc_E[2] = {}, zdc_t[2] = {}, preamp[2][4] = {};
    std::vector<int>* trk_numqual = nullptr;
    std::vector<float>* vtx_z = nullptr;
    for (const char* b : {"b_HLT_mu4_L1MU3V", "RunNumber", "lbn", "bcid", "FCal_Et_P", "FCal_Et_N", "zdc_ZdcEnergy",
                          "zdc_ZdcTime", "zdc_ZdcModulePreSampleAmp", "trk_numqual", "vtx_z"})
        chain.SetBranchStatus(b, 1);
    chain.SetBranchAddress("b_HLT_mu4_L1MU3V", &b_HLT);
    chain.SetBranchAddress("RunNumber", &run);
    chain.SetBranchAddress("lbn", &lbn);
    chain.SetBranchAddress("bcid", &bcid);
    chain.SetBranchAddress("FCal_Et_P", &FCal_Et_P);
    chain.SetBranchAddress("FCal_Et_N", &FCal_Et_N);
    chain.SetBranchAddress("zdc_ZdcEnergy", zdc_E);
    chain.SetBranchAddress("zdc_ZdcTime", zdc_t);
    chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp", preamp);
    chain.SetBranchAddress("trk_numqual", &trk_numqual);
    chain.SetBranchAddress("vtx_z", &vtx_z);

    // ---- output ----
    TFile fo(out_file.c_str(), "RECREATE");
    TH3F* h3 = new TH3F("h3_zdct_vtxz", ";ZDC t^{A} [ns];ZDC t^{C} [ns];z_{vtx}(PV0) [mm]",
                        100, -10., 10., 100, -10., 10., 64, -320., 320.);
    struct RunH { TH2D *frac_c13, *nf_c14, *nf_nocut; TH1D *ratio, *cnt; };
    std::map<int, RunH> rh;
    auto get_run = [&](int r) -> RunH& {
        auto it = rh.find(r);
        if (it != rh.end()) return it->second;
        RunH h;
        h.frac_c13 = new TH2D(Form("h_frac_c13_%d", r), ";N_{trk}^{total} (p_{T}>400 MeV);N_{trk}^{HItight} / N_{trk}^{total}",
                              100, 0., 5000., 110, 0., 1.1);
        h.nf_c14   = new TH2D(Form("h_ntrkfcal_c14_%d", r), ";FCal E_{T}^{A+C} [TeV];N_{trk}^{HItight}",
                              120, -0.5, 5.5, 100, 0., 3500.);
        h.nf_nocut = new TH2D(Form("h_ntrkfcal_nocut_%d", r), ";FCal E_{T}^{A+C} [TeV];N_{trk}^{HItight}",
                              120, -0.5, 5.5, 100, 0., 3500.);
        h.ratio    = new TH1D(Form("h_ratio_nocut_%d", r), ";N_{trk}^{HItight} / FCal E_{T}^{A+C} [TeV^{-1}];events",
                              100, 0., 1000.);
        h.cnt      = new TH1D(Form("h_cnt_%d", r), "", 5, 0.5, 5.5);
        return rh.emplace(r, h).first->second;
    };

    Int_t o_run, o_lbn, o_bcid, o_ntot, o_ntight, o_pass;
    Float_t o_fcal, o_tA, o_tC, o_zdc, o_pA, o_pC, o_vz;
    auto book_tree = [&](const char* name) {
        TTree* t = new TTree(name, "");
        t->Branch("run", &o_run); t->Branch("lbn", &o_lbn); t->Branch("bcid", &o_bcid);
        t->Branch("ntrk_total", &o_ntot); t->Branch("ntrk_tight", &o_ntight);
        t->Branch("fcal", &o_fcal); t->Branch("tA", &o_tA); t->Branch("tC", &o_tC); t->Branch("zdc_tot", &o_zdc);
        t->Branch("preamp_A", &o_pA); t->Branch("preamp_C", &o_pC); t->Branch("vtx_z0", &o_vz);
        t->Branch("pass_bits", &o_pass);  // bit (k-1) set = passes cut k ALONE (k = 1..5)
        return t;
    };
    TTree* t_c4 = book_tree("t_c4fail");
    TTree* t_lb = book_tree("t_lbscan");
    const std::set<int> lbset(lbscan_runs.begin(), lbscan_runs.end());

    const Long64_t n = nmax > 0 ? std::min(nmax, chain.GetEntries()) : chain.GetEntries();
    std::cout << fin << ": " << n << " events" << std::endl;
    for (Long64_t i = 0; i < n; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;
        const int r = (int)run;
        if (bad_runs.count(r)) continue;
        if (!vtx_z || vtx_z->empty()) throw std::runtime_error("event without vtx_z -- skim guarantees PV0");
        const float fcal = (FCal_Et_P + FCal_Et_N) * 1e-6f;
        const float zdc  = (zdc_E[0] + zdc_E[1]) / 1000.f;
        const float tA = zdc_t[1], tC = zdc_t[0];
        float pA = 0.f, pC = 0.f;
        for (int k = 0; k < 4; ++k) { pA += preamp[1][k]; pC += preamp[0][k]; }
        int ntot = 0, ntight = 0;
        const bool has_trk = trk_numqual && (int)trk_numqual->size() >= 4;
        if (has_trk) { ntot = (*trk_numqual)[0]; ntight = (*trk_numqual)[3]; }

        float c3A = cut3_A0, c3C = cut3_C0;
        if (auto it = cut3_run.find(r); it != cut3_run.end()) { c3A = it->second.first; c3C = it->second.second; }
        const bool p1 = !(zdc > (float)PbPbEvSelEvalCut(g1, fcal));
        const bool p2 = std::abs(tA) < (float)cut2_ns && std::abs(tC) < (float)cut2_ns;
        const bool p3 = !(pA >= c3A && pC >= c3C);
        const double frac = ntot > 0 ? (double)ntight / ntot : -1.;
        const bool p4 = has_trk && !(ntot > 0 && frac < PbPbEvSelEvalCut(g4, ntot));
        const bool p5 = has_trk && !((double)ntight < PbPbEvSelEvalCut(g5lo, fcal) || (double)ntight > PbPbEvSelEvalCut(g5hi, fcal));
        const bool p13 = p1 && p2 && p3;

        h3->Fill(tA, tC, vtx_z->at(0));
        RunH& h = get_run(r);
        h.cnt->Fill(1);
        if (has_trk) {
            h.nf_nocut->Fill(fcal, ntight);
            if (fcal > 1.0f) h.ratio->Fill(ntight / fcal);
        }
        if (p13 && has_trk) {
            h.cnt->Fill(2);
            if (ntot > 0) h.frac_c13->Fill(ntot, frac);
            if (!p4) h.cnt->Fill(3);
            else { h.cnt->Fill(4); h.nf_c14->Fill(fcal, ntight); if (!p5) h.cnt->Fill(5); }
        }

        const bool in_c4 = p13 && has_trk && !p4;
        const bool in_lb = lbset.count(r) > 0;
        if (in_c4 || in_lb) {
            o_run = r; o_lbn = (int)lbn; o_bcid = (int)bcid; o_ntot = ntot; o_ntight = ntight;
            o_fcal = fcal; o_tA = tA; o_tC = tC; o_zdc = zdc; o_pA = pA; o_pC = pC; o_vz = vtx_z->at(0);
            o_pass = (p1 << 0) | (p2 << 1) | (p3 << 2) | (p4 << 3) | (p5 << 4);
            if (in_c4) t_c4->Fill();
            if (in_lb) t_lb->Fill();
        }
        if (i % 10000000 == 0 && i) std::cout << "  " << i << " / " << n << std::endl;
    }
    fo.Write();
    fo.Close();
    fc->Close();
    std::cout << "[DONE] " << out_file << std::endl;
}
