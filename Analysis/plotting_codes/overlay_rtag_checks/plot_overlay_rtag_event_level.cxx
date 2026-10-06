// Event-level sanity plots for the HIJING-overlay Pb+Pb-2024-conditions test sample
// (pTH125_300 only, vtx z = -71.2 mm, HIJING b 0-5 fm) against the Pb+Pb-2023-conditions test
// sample (r17618), SAME pT-hat slice so the two event populations match.
// Tracking doc: docs/tracking/hijing_overlay_pbpb24_test_sample_skim.md (R1, R4, R6-R12).
//
// R-TAG UNDER TEST: set by NEW_RTAG below; currently **r17913** (AOD merge tid52733789, skim
// jediTaskID 52785715, downloaded 2026-09-29). r17913 is "clone r17864 with additional preExec
// parameters" -- the digitisation re-configuration answering R9: it adds
// Digitization.PU.Initial/FinalBunchCrossing = 0 (no out-of-time crossings at all),
// LAr.ROD.NumberOfCollisions = 0, LAr.ROD.nSamples = 4 and Digitization.HighGainFCal = True.
// Its two predecessors under r17864 -- the first production (tid52480369) and its bit-identical
// re-run (tid52629267), both with the broken calorimeter of R1/R4/R8/R9 -- are archived in
// broken_calo_tid52480369/ and broken_calo_r17864_rerun_tid52629267/ inside the sample dir.
// The NTUP PATH does NOT carry the r-tag: it is keyed on the production CONFIGURATION
// (vtxz-71_2mm_b0_5fm, FullSimSampleType.h), which is identical across all three, which is
// exactly why only one of them may sit under the canonical name at a time.
//
// Reads the raw skim NTUPs (HeavyIonD3PD) directly: FCal_Et(_P/_N), L1TE, bcid, centrality,
// trk_numqual and vtx_z / truth_vtx_z are per-event skim quantities with no selection, weight or
// truth-navigation procedure behind them, so nothing from the ntuple-processing chain is needed or
// bypassed. `bcid` in particular: the only processing output carrying it is the DATA muon-pair tree
// (DimuonDataAlgCoreT.c), which is per-pair and post-selection and so cannot give the unbiased
// per-event bunch population these checks need.
//
//   root -l -b -q 'plot_overlay_rtag_event_level.cxx+("fcal")'   FCal sum-ET, both samples
//   root -l -b -q 'plot_overlay_rtag_event_level.cxx+("vtx")'    reco vs truth vertex z
//   root -l -b -q 'plot_overlay_rtag_event_level.cxx+("ntrk")'   ID track multiplicities (trk_numqual), both samples
//   root -l -b -q 'plot_overlay_rtag_event_level.cxx+("l1te")'   L1 total transverse energy (L1TE), both samples
//   root -l -b -q 'plot_overlay_rtag_event_level.cxx+("centr")'  centrality population table, both samples (no plot)
//   root -l -b -q 'plot_overlay_rtag_event_level.cxx+("bcid")'   filled-bunch pattern and <FCal> vs bcid, both samples
//   root -l -b -q 'plot_overlay_rtag_event_level.cxx+("trainpos")' <FCal> and <L1TE> vs position in the bunch train
//
// Every event is unweighted: both samples are ONE pT-hat slice (no slice mixing), and the
// quantities are underlying-event / beam-spot properties independent of the hard scatter.
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>
#include <TProfile.h>
#include <map>
#include <algorithm>
#include <vector>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <cstdio>
#include "../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"

namespace {
    // The r-tag under test.  Legend labels, the output directory and every output file name
    // derive from these two constants; the NTUP path does not (configuration-keyed, see above).
    const char* NEW_RTAG = "r17913";
    const char* OLD_RTAG = "r17618";
    // Resolved at first use, like OUTDIR(): a namespace-scope std::string built by concatenation
    // breaks cling's later RDF JIT symbol resolution (same trap as noted for OUTDIR below).
    const std::string& NEW_LABEL() { static const std::string s = std::string("2024 conditions (") + NEW_RTAG + ")"; return s; }
    const std::string& OLD_LABEL() { static const std::string s = std::string("2023 conditions (") + OLD_RTAG + ")"; return s; }
    // pTH125_300 NTUP of each conditions year, named by FullSimSampleType.h (dir + file tag,
    // the latter carrying the production configuration vtxz/b since 2026-09-17).
    std::string ntup(int pbpb_year) {
        return FullSimSampleInputDir(FullSimSampleType::hijing, true, pbpb_year)
             + "Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300." + FullSimSampleFileTag(FullSimSampleType::hijing, pbpb_year)
             + ".NTUP.root";
    }
    // Resolved at first use, not at library load: a namespace-scope std::string initialised
    // from the header's inline functions breaks cling's later RDF JIT symbol resolution.
    const std::string& NEW_FILE() { static const std::string f = ntup(24); return f; }
    const std::string& OLD_FILE() { static const std::string f = ntup(23); return f; }
    const std::string& OUTDIR() {
        static const std::string d = FullSimPlotsDir(FullSimSampleInputDir(FullSimSampleType::hijing, true, 24))
                                   + NEW_RTAG + "_rtag_sanity/";
        return d;
    }
    // Common output-file name: <base>_<new r-tag>_vs_<old r-tag>_pTH125_300.<ext>
    std::string outfile(const char* base, const char* ext) {
        return OUTDIR() + base + "_" + NEW_RTAG + "_vs_" + OLD_RTAG + "_pTH125_300." + ext;
    }
    const int kNew = kRed + 1, kOld = kBlue + 1;

    void style() {
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gSystem->mkdir(OUTDIR().c_str(), true);
    }
    void axes(TH1* h, const char* xt, const char* yt) {
        h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
        h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
        h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetTitleOffset(1.3);
    }
    // FCal_Et is stored in MeV; TeV on the plots.
    ROOT::RDF::RNode fcal_node(ROOT::RDataFrame& df) {
        return df.Define("fcal_tev",   "FCal_Et   / 1e6f")
                 .Define("fcal_a_tev", "FCal_Et_P / 1e6f")   // A side = eta > 0
                 .Define("fcal_c_tev", "FCal_Et_N / 1e6f");  // C side = eta < 0
    }
}

void plot_fcal() {
    style();
    ROOT::RDataFrame dnew("HeavyIonD3PD", NEW_FILE()), dold("HeavyIonD3PD", OLD_FILE());
    auto nn = fcal_node(dnew), no = fcal_node(dold);
    const ROOT::RDF::TH1DModel mtot("", "", 110, -3.0, 8.0), mside("", "", 90, -2.0, 4.0);
    auto h_new = nn.Histo1D(mtot, "fcal_tev"),    h_old = no.Histo1D(mtot, "fcal_tev");
    auto a_new = nn.Histo1D(mside, "fcal_a_tev"), c_new = nn.Histo1D(mside, "fcal_c_tev");
    auto a_old = no.Histo1D(mside, "fcal_a_tev"), c_old = no.Histo1D(mside, "fcal_c_tev");
    const double n_new = *nn.Count(), n_old = *no.Count();
    const double neg_new = *nn.Filter("fcal_tev < 0").Count(), neg_old = *no.Filter("fcal_tev < 0").Count();

    TCanvas c("c", "", 1600, 700);
    c.Divide(2, 1);
    // (a) total FCal sum-ET, per-event fraction so the two 10k samples overlay directly
    c.cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    h_new->Scale(1.0 / n_new); h_old->Scale(1.0 / n_old);
    axes(h_new.GetPtr(), "FCal #Sigma E_{T} [TeV]", "fraction of events / 0.1 TeV");
    h_new->SetLineColor(kNew); h_new->SetLineWidth(2);
    h_old->SetLineColor(kOld); h_old->SetLineWidth(2);
    const double ymax = 1.75 * std::max(h_new->GetBinContent(h_new->GetMaximumBin()),
                                        h_old->GetBinContent(h_old->GetMaximumBin()));
    h_new->SetMaximum(ymax); h_new->SetMinimum(0);
    h_new->Draw("hist"); h_old->Draw("hist same");
    TLine zero(0, 0, 0, ymax); zero.SetLineStyle(2); zero.SetLineColor(kGray + 2); zero.Draw();
    TLegend l1(0.32, 0.62, 0.86, 0.89); l1.SetBorderSize(0); l1.SetFillStyle(0); l1.SetTextSize(0.031);
    l1.SetHeader("HIJING overlay, b < 5 fm, #hat{p}_{T} 125-300 GeV");
    l1.AddEntry(h_new.GetPtr(), Form("%s", NEW_LABEL().c_str()), "l");
    l1.AddEntry((TObject*)nullptr, Form("median %.2f TeV, %.0f%% negative", [&]{ double q, p = 0.5; h_new->GetQuantiles(1, &q, &p); return q; }(), 100.0 * neg_new / n_new), "");
    l1.AddEntry(h_old.GetPtr(), Form("%s", OLD_LABEL().c_str()), "l");
    l1.AddEntry((TObject*)nullptr, Form("median %.2f TeV, %.0f%% negative", [&]{ double q, p = 0.5; h_old->GetQuantiles(1, &q, &p); return q; }(), 100.0 * neg_old / n_old), "");
    l1.Draw();
    // (b) per side, A (eta>0) solid vs C (eta<0) dashed -- a timing/readout problem need not be symmetric
    c.cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    for (auto* h : {a_new.GetPtr(), c_new.GetPtr()}) { h->Scale(1.0 / n_new); h->SetLineColor(kNew); h->SetLineWidth(2); }
    for (auto* h : {a_old.GetPtr(), c_old.GetPtr()}) { h->Scale(1.0 / n_old); h->SetLineColor(kOld); h->SetLineWidth(2); }
    c_new->SetLineStyle(2); c_old->SetLineStyle(2);
    axes(a_new.GetPtr(), "FCal #Sigma E_{T} per side [TeV]", "fraction of events / 0.067 TeV");
    double ymax2 = 0; for (auto* h : {a_new.GetPtr(), c_new.GetPtr(), a_old.GetPtr(), c_old.GetPtr()}) ymax2 = std::max(ymax2, h->GetBinContent(h->GetMaximumBin()));
    a_new->SetMaximum(1.7 * ymax2); a_new->SetMinimum(0);
    a_new->Draw("hist"); c_new->Draw("hist same"); a_old->Draw("hist same"); c_old->Draw("hist same");
    TLine zero2(0, 0, 0, 1.7 * ymax2); zero2.SetLineStyle(2); zero2.SetLineColor(kGray + 2); zero2.Draw();
    // the header defines solid = A and dashed = C, so the per-entry ", side A/C" suffixes are dropped:
    // with them the two longest entries ran into the right frame edge and over the axis ticks
    TLegend l2(0.30, 0.64, 0.89, 0.89); l2.SetBorderSize(0); l2.SetFillStyle(0); l2.SetTextSize(0.032);
    l2.SetHeader("side A: #eta > 0 (solid), side C: #eta < 0 (dashed)");
    l2.AddEntry(a_new.GetPtr(), NEW_LABEL().c_str(), "l");
    l2.AddEntry(a_old.GetPtr(), OLD_LABEL().c_str(), "l");
    l2.Draw();
    c.SaveAs(outfile("fcal_sum_et", "png").c_str());
    std::printf("FCal sum-ET: new %.0f events (%.1f%% < 0), old %.0f events (%.1f%% < 0)\n",
                n_new, 100.0 * neg_new / n_new, n_old, 100.0 * neg_old / n_old);
}

void plot_vtx() {
    style();
    ROOT::RDataFrame dnew("HeavyIonD3PD", NEW_FILE()), dold("HeavyIonD3PD", OLD_FILE());
    // vtx_z[0] = the reconstructed primary vertex; truth_vtx_z[0] = the generated interaction
    // point (every stored truth vertex sits at the signal PV, checked 2026-09-16). The generated
    // vertex is UNSMEARED: truth z is exactly -71.2000 mm in every event of this configuration
    // (the 0.01 mm sigma_z in the r-tag postExec is the reconstruction beam-spot width, not a
    // truth smearing), so a 2D
    // reco-vs-truth view carries no information beyond the residual.
    auto def = [](ROOT::RDataFrame& d) {
        return d.Filter("vtx_z.size() > 0 && truth_vtx_z.size() > 0")
                .Define("z_reco",  "vtx_z[0]").Define("z_truth", "truth_vtx_z[0]")
                .Define("dz_um",   "(vtx_z[0] - truth_vtx_z[0]) * 1000.f");     // mm -> um
    };
    auto nn = def(dnew), no = def(dold);
    // edges offset by half a bin so the truth value -71.2000 mm is bin-CENTRED, not on an edge
    auto h_reco  = nn.Histo1D({"", "", 40, -71.2975, -71.0975}, "z_reco");   // data span -71.25..-71.16 mm
    auto h_truth = nn.Histo1D({"", "", 40, -71.2975, -71.0975}, "z_truth");
    auto r_new   = nn.Histo1D({"", "", 80, -80., 80.}, "dz_um");
    auto r_old   = no.Histo1D({"", "", 80, -80., 80.}, "dz_um");

    TCanvas c("c", "", 1600, 700);
    c.Divide(2, 1);
    // log y: the truth point is a delta (all 1e4 events in one 5 um bin) next to a reco
    // distribution falling over ~3 decades to single events; linear would hide the reco tails.
    c.cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogy();
    axes(h_reco.GetPtr(), "vertex z [mm]", "events / 5 #mum");
    h_reco->GetXaxis()->SetNdivisions(505);
    h_reco->SetLineColor(kNew); h_reco->SetLineWidth(2);
    h_truth->SetLineColor(kBlack); h_truth->SetLineWidth(2); h_truth->SetLineStyle(2);
    h_reco->SetMaximum(30 * std::max(h_reco->GetBinContent(h_reco->GetMaximumBin()), h_truth->GetBinContent(h_truth->GetMaximumBin())));
    h_reco->SetMinimum(0.5);
    h_reco->Draw("hist"); h_truth->Draw("hist same");
    TLegend l1(0.16, 0.72, 0.92, 0.90); l1.SetBorderSize(0); l1.SetFillStyle(0); l1.SetTextSize(0.035);
    l1.SetHeader(Form("Pb+Pb %s, #hat{p}_{T} 125-300 GeV", NEW_LABEL().c_str()));
    l1.AddEntry(h_reco.GetPtr(),  "reconstructed primary vertex", "l");
    l1.AddEntry(h_truth.GetPtr(), "truth interaction point", "l");
    l1.Draw();
    c.cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    axes(r_new.GetPtr(), "z_{reco} - z_{truth} [#mum]", "events / 2 #mum");
    r_new->SetLineColor(kNew); r_new->SetLineWidth(2);
    r_old->SetLineColor(kOld); r_old->SetLineWidth(2);
    r_new->SetMaximum(1.7 * std::max(r_new->GetBinContent(r_new->GetMaximumBin()), r_old->GetBinContent(r_old->GetMaximumBin())));
    r_new->SetMinimum(0);
    r_new->Draw("hist"); r_old->Draw("hist same");
    TLegend l3(0.16, 0.64, 0.92, 0.90); l3.SetBorderSize(0); l3.SetFillStyle(0); l3.SetTextSize(0.034);
    l3.AddEntry(r_new.GetPtr(), Form("%s, z_{vtx} = -71.2 mm", NEW_LABEL().c_str()), "l");
    l3.AddEntry((TObject*)nullptr, Form("mean %+.1f #mum, RMS %.1f #mum", r_new->GetMean(), r_new->GetRMS()), "");
    l3.AddEntry(r_old.GetPtr(), Form("%s, z_{vtx} = -3.3 mm", OLD_LABEL().c_str()), "l");
    l3.AddEntry((TObject*)nullptr, Form("mean %+.1f #mum, RMS %.1f #mum", r_old->GetMean(), r_old->GetRMS()), "");
    l3.Draw();
    c.SaveAs((OUTDIR() + "vertex_z_reco_vs_truth_" + NEW_RTAG + ".png").c_str());
    std::printf("vertex z: new reco mean %.4f mm, truth mean %.4f mm; residual mean %+.2f um RMS %.2f um (old: %+.2f / %.2f um)\n",
                h_reco->GetMean(), h_truth->GetMean(), r_new->GetMean(), r_new->GetRMS(), r_old->GetMean(), r_old->GetRMS());
}

// ID track multiplicity per event, two skim counters (TrigRates.cxx m_trk_numqual): the
// question is whether the 2024-conditions calorimeter problem (R1) also reaches the
// tracker. No centrality cut (an FCal-derived centrality is meaningless while the FCal is
// wrong, and this check must not depend on it either way). Unweighted, per-event
// fraction, so the two 10k samples overlay directly.
void plot_ntrk() {
    style();
    ROOT::RDataFrame dnew("HeavyIonD3PD", NEW_FILE()), dold("HeavyIonD3PD", OLD_FILE());
    // Two panels: [0] p_T > 400 MeV, no quality requirement; [3] p_T > 400 MeV + HITight.
    // The untagged counter reaches ~12k tracks in the 2024 sample (250 / bin), HITight ~4k (80/bin).
    //
    // WHY ONLY TWO. m_trk_numqual has 8 slots: [0-3] are {no quality, PP_MIN_BIAS, HI_LOOSE,
    // HI_TIGHT} WITH p_T > 400 MeV, and [4-7] are the SAME four WITHOUT the p_T cut (added
    // Dec 2023 for exactly that purpose, TrigRates.cxx:1838). They carry no extra information
    // here: [4]==[0], [5]==[1], [6]==[2], [7]==[3] in 10 000 of 10 000 events in BOTH samples,
    // because the loop runs over InDetTrackParticles (TrigRates.cxx:85) and the HI track
    // reconstruction of these AODs stores nothing below 400 MeV. The p_T threshold is therefore
    // upstream of the skim: it cannot be lifted by re-skimming, only by a different reconstruction,
    // and there is no fallback either (the per-track track_pt branch is not stored in these
    // overlay NTUPs -- StoreBasic is off). Plotting [4] and [7] as separate panels only invited
    // the reading "the p_T cut removes nothing", which is backwards.
    struct Panel { int idx; const char* title; int nbins; double xmax; double width; };
    const Panel panels[2] = {
        {0, "p_{T} > 400 MeV",          56, 14000., 250.},
        {3, "p_{T} > 400 MeV, HITight", 50,  4000.,  80.},
    };
    const double n_new = *dnew.Count(), n_old = *dold.Count();
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hn, ho;
    for (const auto& p : panels) {
        const std::string col = "ntrk_" + std::to_string(p.idx);
        const std::string expr = "(double)trk_numqual[" + std::to_string(p.idx) + "]";
        const ROOT::RDF::TH1DModel m("", "", p.nbins, 0., p.xmax);
        hn.push_back(dnew.Define(col, expr).Histo1D(m, col));
        ho.push_back(dold.Define(col, expr).Histo1D(m, col));
    }
    TCanvas c("c", "", 1600, 700);
    c.Divide(2, 1);
    for (int i = 0; i < 2; ++i) {
        c.cd(i + 1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        TH1D *n = hn[i].GetPtr(), *o = ho[i].GetPtr();
        n->Scale(1.0 / n_new); o->Scale(1.0 / n_old);
        axes(n, Form("N_{trk} (%s)", panels[i].title), Form("fraction of events / %.0f tracks", panels[i].width));
        n->GetXaxis()->SetNdivisions(507);
        n->SetLineColor(kNew); n->SetLineWidth(2);
        o->SetLineColor(kOld); o->SetLineWidth(2);
        const double ymax = 1.6 * std::max(n->GetBinContent(n->GetMaximumBin()), o->GetBinContent(o->GetMaximumBin()));
        n->SetMaximum(ymax); n->SetMinimum(0);
        n->Draw("hist"); o->Draw("hist same");
        TLegend l(0.40, 0.62, 0.90, 0.89); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.030);
        l.SetHeader("HIJING overlay, b < 5 fm, #hat{p}_{T} 125-300 GeV");
        l.AddEntry(n, NEW_LABEL().c_str(), "l");
        l.AddEntry((TObject*)nullptr, Form("mean %.0f, RMS %.0f", n->GetMean(), n->GetRMS()), "");
        l.AddEntry(o, OLD_LABEL().c_str(), "l");
        l.AddEntry((TObject*)nullptr, Form("mean %.0f, RMS %.0f", o->GetMean(), o->GetRMS()), "");
        l.DrawClone();
        std::printf("trk_numqual[%d] (%s): new mean %.0f RMS %.0f | old mean %.0f RMS %.0f | ratio of means %.3f\n",
                    panels[i].idx, panels[i].title, n->GetMean(), n->GetRMS(), o->GetMean(), o->GetRMS(), n->GetMean() / o->GetMean());
    }
    c.SaveAs(outfile("trk_multiplicity", "png").c_str());
}

// L1 total transverse energy (LVL1EnergySumRoI energyT, stored in GeV by the skim) -- the
// L1Calo view of the same calorimeter whose offline sums were broken in r17864 (R1). No
// centrality cut, unweighted, per-event fraction.
void plot_l1te() {
    style();
    ROOT::RDataFrame dnew("HeavyIonD3PD", NEW_FILE()), dold("HeavyIonD3PD", OLD_FILE());
    auto nn = dnew.Define("l1te_tev", "L1TE / 1e3f"), no = dold.Define("l1te_tev", "L1TE / 1e3f");
    // the 2023-conditions sample reaches ~25 TeV (b < 5 fm), the 2024 one sits at ~3 TeV
    const ROOT::RDF::TH1DModel m("", "", 100, 0., 25.);
    auto h_new = nn.Histo1D(m, "l1te_tev"), h_old = no.Histo1D(m, "l1te_tev");
    const double n_new = *nn.Count(), n_old = *no.Count();
    TCanvas c("c", "", 800, 700);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    h_new->Scale(1.0 / n_new); h_old->Scale(1.0 / n_old);
    axes(h_new.GetPtr(), "L1 #Sigma E_{T} [TeV]", "fraction of events / 0.25 TeV");
    h_new->SetLineColor(kNew); h_new->SetLineWidth(2);
    h_old->SetLineColor(kOld); h_old->SetLineWidth(2);
    const double ymax = 1.6 * std::max(h_new->GetBinContent(h_new->GetMaximumBin()), h_old->GetBinContent(h_old->GetMaximumBin()));
    h_new->SetMaximum(ymax); h_new->SetMinimum(0);
    h_new->Draw("hist"); h_old->Draw("hist same");
    auto median = [](TH1D* h) { double q, p = 0.5; h->GetQuantiles(1, &q, &p); return q; };
    TLegend l(0.36, 0.62, 0.90, 0.89); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.032);
    l.SetHeader("HIJING overlay, b < 5 fm, #hat{p}_{T} 125-300 GeV");
    l.AddEntry(h_new.GetPtr(), NEW_LABEL().c_str(), "l");
    l.AddEntry((TObject*)nullptr, Form("median %.2f TeV", median(h_new.GetPtr())), "");
    l.AddEntry(h_old.GetPtr(), OLD_LABEL().c_str(), "l");
    l.AddEntry((TObject*)nullptr, Form("median %.2f TeV", median(h_old.GetPtr())), "");
    l.Draw();
    c.SaveAs(outfile("l1te", "png").c_str());
    std::printf("L1TE: new median %.2f TeV (mean %.2f), old median %.2f TeV (mean %.2f); overflow >25 TeV new %.0f old %.0f events\n",
                median(h_new.GetPtr()), h_new->GetMean(), median(h_old.GetPtr()), h_old->GetMean(),
                h_new->GetBinContent(101) * n_new, h_old->GetBinContent(101) * n_old);
}


// Centrality population table (check-list item 5). `centrality` is an INT percentile index
// 0..84 from the 2023 FCal calibration (SkimCode Module_EventShape.cxx: 0-based index into
// FCal_ET_Bins_PbPb2023, applied to ALL Run-3 years), -1 = below the 85 % edge, which includes
// every negative-FCal event. No plot: the doc records this as a table.
void plot_centr() {
    style();
    ROOT::RDataFrame dnew("HeavyIonD3PD", NEW_FILE()), dold("HeavyIonD3PD", OLD_FILE());
    auto cn = *dnew.Take<int>("centrality"), co = *dold.Take<int>("centrality");
    struct Bin { const char* name; int lo, hi; };                 // [lo, hi) in percent
    const Bin bins[] = {{"0-10 %", 0, 10}, {"10-20 %", 10, 20}, {"20-30 %", 20, 30},
                        {"30-40 %", 30, 40}, {"40-50 %", 40, 50}, {"50-80 %", 50, 80},
                        {"80-100 %", 80, 101}};
    auto frac = [](const std::vector<int>& v, int lo, int hi) {
        return 100.0 * std::count_if(v.begin(), v.end(), [&](int c){ return c >= lo && c < hi; }) / v.size();
    };
    auto undef = [](const std::vector<int>& v) {
        return 100.0 * std::count_if(v.begin(), v.end(), [](int c){ return c < 0; }) / v.size();
    };
    const std::string out = outfile("centrality_population", "txt");
    FILE* f = std::fopen(out.c_str(), "w");
    auto emit = [&](const char* fmt, auto... a) { std::printf(fmt, a...); std::fprintf(f, fmt, a...); };
    emit("centrality population, pT-hat 125-300 GeV, %zu vs %zu events (2023 FCal calibration)\n",
         cn.size(), co.size());
    emit("%-12s %14s %14s\n", "bin", (std::string("pbpb24 ") + NEW_RTAG).c_str(),
         (std::string("pbpb23 ") + OLD_RTAG).c_str());
    for (const auto& b : bins) emit("%-12s %13.2f%% %13.2f%%\n", b.name, frac(cn, b.lo, b.hi), frac(co, b.lo, b.hi));
    emit("%-12s %13.2f%% %13.2f%%\n", "undefined", undef(cn), undef(co));
    emit("%-12s %13.2f%% %13.2f%%\n", "(0-5 %)", frac(cn, 0, 5), frac(co, 0, 5));
    emit("%-12s %13.2f%% %13.2f%%\n", "(5-10 %)", frac(cn, 5, 10), frac(co, 5, 10));
    std::fclose(f);
    std::printf("wrote %s\n", out.c_str());
}

// Bunch-crossing structure. The ONE stated r17835 -> r17864 change was a DIGITISATION bunch-structure
// option (Digitization.PU.BunchStructureConfig = BunchStructureHeavyIon2022), which sets where in the
// bunch train the signal sits; sampling the bipolar LAr pulse off-peak is the standard way to get
// energies that are far too low AND negative. Two views: (a) which bcids are filled at all -- does the
// 2024 production use a different bunch pattern than the 2023 one; (b) <FCal Sum E_T> vs bcid -- is the
// deficit bcid-dependent (a timing explanation) or common to every crossing.
// Measured filled pattern, IDENTICAL in both samples (1010 of 3564 crossings): 882 gaps of 2 bcid and
// 108 of 4, forming 18 trains of 7x8 bunches -- 8 bunches at 2-bcid spacing, then one crossing skipped,
// so a sub-burst repeats every 18 bcid -- plus 2 ISOLATED single bunches (826, 2611).
void plot_bcid() {
    style();
    ROOT::RDataFrame dnew("HeavyIonD3PD", NEW_FILE()), dold("HeavyIonD3PD", OLD_FILE());
    auto nn = dnew.Define("bc", "(double)bcid").Define("fcal_tev", "FCal_Et / 1e6f");
    auto no = dold.Define("bc", "(double)bcid").Define("fcal_tev", "FCal_Et / 1e6f");
    // (a) zoom on the first 400 crossings, ONE bin per bcid, so individual filled bunches are visible
    auto z_new = nn.Histo1D({"", "", 400, 0.5, 400.5}, "bc");
    auto z_old = no.Histo1D({"", "", 400, 0.5, 400.5}, "bc");
    // (b) profile over the full LHC orbit, 12 bcid per bin (~10 events/bcid -> ~120 per bin)
    auto p_new = nn.Profile1D({"", "", 297, 0., 3564.}, "bc", "fcal_tev");
    auto p_old = no.Profile1D({"", "", 297, 0., 3564.}, "bc", "fcal_tev");

    TCanvas c("c", "", 1600, 700);
    c.Divide(2, 1);
    c.cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    axes(z_new.GetPtr(), "bunch crossing ID (first 400)", "events / bcid");
    z_new->SetLineColor(kNew); z_new->SetLineWidth(2);
    z_old->SetLineColor(kOld); z_old->SetLineWidth(2);
    z_new->SetMaximum(1.7 * std::max(z_new->GetBinContent(z_new->GetMaximumBin()),
                                     z_old->GetBinContent(z_old->GetMaximumBin())));
    z_new->SetMinimum(0);
    z_new->Draw("hist"); z_old->Draw("hist same");
    TLegend l1(0.34, 0.70, 0.90, 0.89); l1.SetBorderSize(0); l1.SetFillStyle(0); l1.SetTextSize(0.032);
    l1.SetHeader("filled-bunch pattern");
    l1.AddEntry(z_new.GetPtr(), NEW_LABEL().c_str(), "l");
    l1.AddEntry(z_old.GetPtr(), OLD_LABEL().c_str(), "l");
    l1.Draw();
    c.cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    axes(p_new.GetPtr(), "bunch crossing ID", "#LT FCal #Sigma E_{T} #GT [TeV]");
    p_new->SetLineColor(kNew); p_new->SetMarkerColor(kNew); p_new->SetMarkerStyle(20); p_new->SetMarkerSize(0.6);
    p_old->SetLineColor(kOld); p_old->SetMarkerColor(kOld); p_old->SetMarkerStyle(24); p_old->SetMarkerSize(0.6);
    p_new->SetMinimum(-2.0); p_new->SetMaximum(8.0);
    p_new->Draw("p"); p_old->Draw("p same");
    TLine zero(0., 0., 3564., 0.); zero.SetLineStyle(2); zero.SetLineColor(kGray + 2); zero.Draw();
    TLegend l2(0.34, 0.70, 0.90, 0.89); l2.SetBorderSize(0); l2.SetFillStyle(0); l2.SetTextSize(0.032);
    l2.SetHeader("12 bcid per point");
    l2.AddEntry(p_new.GetPtr(), NEW_LABEL().c_str(), "p");
    l2.AddEntry(p_old.GetPtr(), OLD_LABEL().c_str(), "p");
    l2.Draw();
    c.SaveAs(outfile("bcid_structure", "png").c_str());

    // numbers for the doc: spread of the profile tells flat-vs-bcid-dependent
    auto spread = [](TProfile* p) {
        double lo = 1e9, hi = -1e9, sum = 0; int n = 0;
        for (int i = 1; i <= p->GetNbinsX(); ++i) {
            if (p->GetBinEntries(i) < 20) continue;
            const double v = p->GetBinContent(i);
            lo = std::min(lo, v); hi = std::max(hi, v); sum += v; ++n;
        }
        std::printf("      bins with >=20 events: %d, <FCal> range %.2f .. %.2f TeV, mean %.2f TeV\n", n, lo, hi, sum / n);
    };
    auto filled = [](ROOT::RDF::RNode d) {
        auto h = d.Histo1D({"", "", 3564, 0.5, 3564.5}, "bc");
        int n = 0;
        for (int i = 1; i <= h->GetNbinsX(); ++i) if (h->GetBinContent(i) > 0) ++n;
        return n;
    };
    std::printf("bcid: distinct filled bcids (of 3564) -- pbpb24 %s %d, pbpb23 %s %d\n",
                NEW_RTAG, filled(nn), OLD_RTAG, filled(no));
    std::printf("   pbpb24 %s:", NEW_RTAG); spread(p_new.GetPtr());
    std::printf("   pbpb23 %s:", OLD_RTAG); spread(p_old.GetPtr());
}


// Position WITHIN the bunch train. This is the mode that localised the r17864 failure (R9) and
// therefore the mode that decides whether a re-configured r-tag has fixed it. The filled
// pattern (identical in both samples) is 18 trains of 7x8 bunches plus 2 ISOLATED single bunches;
// a train boundary is a bcid DIFFERENCE of more than 5 (i.e. more than 4 empty crossings), and an
// isolated bunch is a train of length 1.
// In the 2024 sample the ONLY crossings that reproduce the 2023 reference are the two isolated ones.
// Every crossing inside a train is wrong, INCLUDING the train head (FCal 15 % high, L1TE 36 % low),
// and beyond the head <FCal> is not a monotonic decay but PERIODIC with the 8-bunch sub-burst:
// it falls to about -0.75 TeV and recovers to about +2.8 TeV once per sub-burst.
// That a bunch with no predecessors is correct while every bunch with predecessors is not points at
// the treatment of the bunch-train signal history in digitisation (see the tracking doc, R9).
namespace {
    // An isolated bunch has NO preceding bunch at all, so it is the control: whatever a train does to
    // its members cannot act on it. It must NOT share a position bin with the train heads -- they are
    // different populations (the heads are 15 % high in FCal and 36 % low in L1TE, the isolated
    // bunches agree with the 2023 reference in both).
    constexpr int kIsolated = -2;
    // bcid -> position in its train (or kIsolated), built from the filled pattern of the sample itself
    std::map<unsigned int, int> train_positions(ROOT::RDataFrame& df) {
        auto bc = *df.Take<unsigned int>("bcid");
        std::vector<unsigned int> f(bc.begin(), bc.end());
        std::sort(f.begin(), f.end());
        f.erase(std::unique(f.begin(), f.end()), f.end());
        std::map<unsigned int, int> pos;
        size_t i = 0;
        while (i < f.size()) {
            size_t j = i;                                    // [i, j] = one train (gap > 5 ends it)
            while (j + 1 < f.size() && f[j+1] - f[j] <= 5) ++j;   // bcid diff > 5 ends the train
            const bool isolated = (j == i);
            for (size_t k = i; k <= j; ++k) pos[f[k]] = isolated ? kIsolated : (int)(k - i);
            i = j + 1;
        }
        return pos;
    }
}

void plot_trainpos() {
    style();
    ROOT::RDataFrame dnew("HeavyIonD3PD", NEW_FILE()), dold("HeavyIonD3PD", OLD_FILE());
    auto pn = train_positions(dnew), po = train_positions(dold);
    auto with_pos = [](ROOT::RDataFrame& d, const std::map<unsigned int,int>& m) {
        return d.Define("tpos", [m](unsigned int b) { auto it = m.find(b); return it == m.end() ? -1 : it->second; }, {"bcid"})
                .Define("fcal_tev", "FCal_Et / 1e6f")
                .Define("l1te_tev", "L1TE / 1e3f");
    };
    auto nn = with_pos(dnew, pn), no = with_pos(dold, po);
    auto in_train = [](ROOT::RDF::RNode d) { return d.Filter("tpos >= 0"); };
    auto isolated = [](ROOT::RDF::RNode d) { return d.Filter("tpos == -2"); };
    const ROOT::RDF::TProfile1DModel m("", "", 56, -0.5, 55.5);
    auto f_new = in_train(nn).Profile1D(m, "tpos", "fcal_tev"), f_old = in_train(no).Profile1D(m, "tpos", "fcal_tev");
    auto l_new = in_train(nn).Profile1D(m, "tpos", "l1te_tev"), l_old = in_train(no).Profile1D(m, "tpos", "l1te_tev");
    // the isolated bunches: the only crossings in the sample with no preceding bunch
    auto iso_f_new = isolated(nn).Mean("fcal_tev"), iso_l_new = isolated(nn).Mean("l1te_tev");
    auto iso_f_old = isolated(no).Mean("fcal_tev"), iso_l_old = isolated(no).Mean("l1te_tev");
    auto iso_n_new = isolated(nn).Count(), iso_n_old = isolated(no).Count();
    // beyond the train head the modulation repeats with the 8-bunch sub-burst
    auto fold = [](ROOT::RDF::RNode d) {
        return d.Filter("tpos >= 8").Define("phase", "tpos % 8")
                .Profile1D({"", "", 8, -0.5, 7.5}, "phase", "fcal_tev");
    };
    auto fold_new = fold(nn), fold_old = fold(no);

    TCanvas c("c", "", 1600, 700);
    c.Divide(2, 1);
    // The data fill the frame, so the legend gets its OWN reserved strip above it (a transparent
    // legend inside the frame would be painted over by the points, and these entries are long).
    auto draw = [&](int pad, TProfile* pnew, TProfile* pold, double iso_new, const char* yt,
                    const char* unit, double lo, double hi) {
        c.cd(pad); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.24);
        axes(pnew, "position in bunch train", yt);
        pnew->SetLineColor(kNew); pnew->SetMarkerColor(kNew); pnew->SetMarkerStyle(20); pnew->SetMarkerSize(0.8);
        pold->SetLineColor(kOld); pold->SetMarkerColor(kOld); pold->SetMarkerStyle(24); pold->SetMarkerSize(0.8);
        pnew->SetMinimum(lo); pnew->SetMaximum(hi);
        pnew->Draw("p"); pold->Draw("p same");
        TLine* z = new TLine(-0.5, 0., 55.5, 0.); z->SetLineStyle(2); z->SetLineColor(kGray + 2); z->Draw();
        TLine* iso = new TLine(-0.5, iso_new, 55.5, iso_new);
        iso->SetLineStyle(7); iso->SetLineColor(kNew); iso->SetLineWidth(2); iso->Draw();
        TLegend* l = new TLegend(0.13, 0.765, 0.97, 0.995);
        l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.030);
        l->SetHeader("18 trains of 7#times8 bunches + 2 isolated bunches");
        l->AddEntry(pnew, NEW_LABEL().c_str(), "p");
        l->AddEntry(iso, Form("isolated bunches, no train: %.2f %s", iso_new, unit), "l");
        l->AddEntry(pold, OLD_LABEL().c_str(), "p");
        l->Draw();
    };
    draw(1, f_new.GetPtr(), f_old.GetPtr(), *iso_f_new, "#LT FCal #Sigma E_{T} #GT [TeV]", "TeV", -2.0, 6.5);
    draw(2, l_new.GetPtr(), l_old.GetPtr(), *iso_l_new, "#LT L1 #Sigma E_{T} #GT [TeV]",   "TeV", -2.0, 19.0);
    c.SaveAs(outfile("train_position", "png").c_str());

    std::printf("position in train (IN-TRAIN bunches only):  <FCal> 2024 / 2023 [TeV]   <L1TE> 2024 / 2023 [TeV]\n");
    for (int i = 1; i <= 10; ++i)
        std::printf("   pos %2d (N=%5.0f): %+7.3f / %7.3f      %7.3f / %7.3f\n",
                    i-1, f_new->GetBinEntries(i), f_new->GetBinContent(i), f_old->GetBinContent(i),
                    l_new->GetBinContent(i), l_old->GetBinContent(i));
    std::printf("   isolated bunches (no train): 2024 <FCal> %.3f, <L1TE> %.3f (N=%lld);"
                " 2023 <FCal> %.3f, <L1TE> %.3f (N=%lld)\n",
                *iso_f_new, *iso_l_new, (long long)*iso_n_new, *iso_f_old, *iso_l_old, (long long)*iso_n_old);
    std::printf("   positions 8-55 folded modulo 8 (the sub-burst period), <FCal> 2024 / 2023 [TeV]:\n");
    for (int i = 1; i <= 8; ++i)
        std::printf("      phase %d: %+7.3f +- %.3f / %7.3f\n", i-1,
                    fold_new->GetBinContent(i), fold_new->GetBinError(i), fold_old->GetBinContent(i));
}

void plot_overlay_rtag_event_level(const char* what = "fcal") {
    const std::string w = what;
    if      (w == "fcal") plot_fcal();
    else if (w == "vtx")  plot_vtx();
    else if (w == "ntrk") plot_ntrk();
    else if (w == "l1te") plot_l1te();
    else if (w == "centr") plot_centr();
    else if (w == "bcid") plot_bcid();
    else if (w == "trainpos") plot_trainpos();
    else std::printf("unknown mode '%s' (fcal|vtx|ntrk|l1te|centr|bcid|trainpos)\n", what);
}
