// =================================================================================================
// plot_mc_data_compr_signal.cxx -- the SIGNAL family of the pp24 MC-vs-data comparison.
//
// Data and Pythia carry the SAME cuts: the single-b signal region (m_uu in (1.08, 2.9) GeV,
// pair pT > ParamsSet::signal_pair_pt_min (9 GeV since 2026-09-08), the fiducial gap cut on BOTH muons), taken in TRUTH quantities on the MC side.
// Every efficiency correction on the data side is therefore applied inside the region in which it
// was measured, which is what makes this family -- and not the generic one -- a genuine (fiducial,
// reco-level, background-unsubtracted) cross-section comparison.
//
//   DATA:   d(sigma)/dX = (1/L_int) * SUM_{OS pairs in the signal region} 1/(eps_trig * eps_reco)
//           The histograms READ HERE already carry the 1/L_int factor (`*_dsigma`), so
//           `data_already_dsigma` is set and the luminosity factor is NOT applied a second time.
//   PYTHIA: d(sigma)/dX from the pp24-condition FULLSIM FULL sample -- pp beam only, isospin
//           weight 1, AMI-weighted -- times nb -> pb = 1000 and nothing else.
//
// The two pads use DIFFERENT MC selections, and that is physics, not an oversight:
//   OS pad: `_single_b_pass_signal_truth_gapcut` (truth OS pairs from the SAME b -- the single-b
//           signal the measurement is after) AND `_op_pass_signal_truth_gapcut` (ALL truth OS pairs), the
//           like-for-like partner of the data OS yield, which still contains gluon-splitting and
//           combinatorial background. Drawing both makes that one remaining deliberate mismatch
//           visible instead of implicit.
//   SS pad: `_ss_pass_signal_truth_gapcut` (ALL truth same-sign pairs). `from_same_b` has no same-sign
//           counterpart by construction, and the data SS pad is the combinatorial-background
//           estimate, whose MC partner is the inclusive SS yield.
// Each pad's legend names its own MC selection.
//
// RATIO PAD. Every panel carries an MC/data ratio pad. On the OS pad BOTH MC curves get a ratio
// marker -- single-b in black (the headline number) and all-OS in blue -- because the gap between
// those two ratios IS the background the data curve still contains, and showing only one of them
// would hide the size of the one deliberate mismatch this family has. On the SS pad the ratio is
// the pad's own MC curve over the data. Errors: TH1::Divide without "B" (independent samples).
//
// NOT in this family: unity-normalized variants, jacobian-corrected variants, POWHEG.
// pair_pt and the 9-panel pair_pt_in_eta view come from plot_mc_data_pair_pt_in_eta.cxx, which
// projects the SAME 2D histogram for both, so the two views cannot drift apart.
//
// Output: /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/signal/
// Usage:  root -l -b -q 'plot_mc_data_compr_signal.cxx+()'            (Tight, nominal)
//         root -l -b -q 'plot_mc_data_compr_signal.cxx+("medium")'    (Medium WP systematic)
// =================================================================================================

#include "PlotMCDataComprBaseClass.c"
#include <TLatex.h>

class PlotMCDataComprSignalKinematics : public PlotMCDataComprBaseClass{
protected:
    // One drawable series: a histogram key, its legend text and its colour.
    struct Series {
        std::string key;
        std::string legend;
        Color_t     color;
        Style_t     mstyle;
        TH1D*       h = nullptr;
    };

public:
    // The MC OS-inclusive curve. On by default: it is the like-for-like partner of the data OS
    // yield, and without it the reader cannot tell how much of the data/MC gap is background.
    bool draw_mc_all_os = true;

    PlotMCDataComprSignalKinematics(std::string kin_input){
        kin = kin_input;
        kin_title = Observable(kin_input).x_title;
        output_subdir = "signal";
        data_already_dsigma = true;   // the `*_dsigma` histograms already carry 1/L_int
        logy = true;                  // 1D distributions, per the repo standard
    }

    ~PlotMCDataComprSignalKinematics(){}

    void Run() override;
};


void PlotMCDataComprSignalKinematics::Run(){
    initialize();

    const McDataObservable& obs = Observable(kin);
    const std::string ytitle = YTitle(obs, /*jacobian=*/false);

    // 700 px tall, not 500: the bottom 30 % is now the ratio pad, and the main pad must not end
    // up shorter than it was before the ratio pad existed.
    TCanvas* c = new TCanvas("c_sig","c_sig",1200,700);
    c->Divide(s_nSigns,1);

    std::vector<TH1D*> owned;

    for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){
        TPad* pad_main = nullptr;
        TPad* pad_ratio = nullptr;
        McDataComprRatio::SplitPadForRatio(c->cd(ksign + 1), pad_main, pad_ratio);

        // --- the series of this pad -----------------------------------------------------------
        // COLOUR CONVENTION of the single-b (signal) family, user instruction 2026-08-25:
        //   pp data          -> RED
        //   Pythia (primary) -> BLACK   (single-b on the OS pad, all-SS on the SS pad -- the
        //                                one curve that is the data's counterpart on that pad)
        //   Pythia, all OS   -> BLUE    (the secondary reference curve, deliberately NOT black
        //                                so it can never be mistaken for the single-b signal)
        // These are LOCAL to the signal family; the generic family keeps `colors[]`.
        std::vector<Series> series;
        series.push_back({SignalDataKey(obs, ksign), "pp data 2024",
                          kSignalDataColor, 20, nullptr});
        const int i_data = 0;

        if (ksign == 0){
            series.push_back({SignalMcKey(obs, McSignalCatAllSS()), "Pythia, all SS",
                              kSignalMcColor, 21, nullptr});
        } else {
            series.push_back({SignalMcKey(obs, McSignalCatSingleB()), "Pythia, single-b",
                              kSignalMcColor, 21, nullptr});
            if (draw_mc_all_os)
                series.push_back({SignalMcKey(obs, McSignalCatAllOS()), "Pythia, all OS",
                                  kSignalMcRefColor, 22, nullptr});
        }

        // --- fetch ----------------------------------------------------------------------------
        for (size_t i = 0; i < series.size(); ++i){
            const int idt = (int)i == i_data ? (int)DataType::pp_2024_2mu4 : (int)DataType::pythia;
            try {
                series[i].h = GetHist1D(idt, series[i].key);
                owned.push_back(series[i].h);
            } catch (const std::exception& e){
                std::cerr << "[SKIP] " << e.what() << std::endl;
            }
        }

        // Data and MC are booked from the SAME registered binning, so a divergence is a producer
        // bug and must stop the plot rather than silently compare different cells.
        for (size_t i = 1; i < series.size(); ++i)
            if (series[i_data].h && series[i].h)
                AssertSameAxis1D(series[i_data].h, series[i].h,
                                 "signal " + kin + " (" + series[i_data].key + " vs "
                                 + series[i].key + ")");

        // --- normalize & style ------------------------------------------------------------------
        for (size_t i = 0; i < series.size(); ++i){
            if (!series[i].h) continue;
            const float norm = ((int)i == i_data) ? norm_factor[DataType::pp_2024_2mu4]
                                                  : norm_factor[DataType::pythia];
            hist_helper(series[i].h, norm, false, signTitles[ksign].c_str(), ytitle);
            series[i].h->SetLineColor(series[i].color);
            series[i].h->SetMarkerColor(series[i].color);
            series[i].h->SetMarkerStyle(series[i].mstyle);
            series[i].h->GetXaxis()->SetTitle(obs.x_title.c_str());
        }

        int first = -1;
        for (size_t i = 0; i < series.size(); ++i) if (series[i].h) { first = (int)i; break; }
        if (first < 0){
            std::cerr << "[WARN] No histograms available for " << kin << " sign " << ksign << std::endl;
            continue;
        }

        // --- ratios, built BEFORE the main pad hides its x axis -----------------------------------
        // Order matters: McDataComprRatio::MakeRatio CLONES the numerator, so a ratio built after
        // HideXAxis() would inherit the suppressed x labels and title and the ratio pad would come
        // out with no x axis at all. (It silently did on the generic family, where the first drawn
        // series is the numerator of the first ratio.)
        // Every MC curve on this pad is divided by the pad's own DATA curve -- always the data, so
        // the OS pad's two markers are directly comparable and their vertical separation is the
        // background content of the data OS yield.
        std::vector<TH1*> ratios;
        if (series[i_data].h){
            for (size_t i = 0; i < series.size(); ++i){
                if ((int)i == i_data || !series[i].h) continue;
                TH1D* r = McDataComprRatio::MakeRatio(
                    series[i].h, series[i_data].h,
                    Form("r_%s_%u_%zu", kin.c_str(), ksign, i));
                if (r){ ratios.push_back(r); owned.push_back(r); }
            }
        }

        // ================================ MAIN PAD ==============================================
        pad_main->cd();
        gPad->SetLogx(obs.logx);
        gPad->SetLogy(logy);

        // --- common y range across the drawn series ---------------------------------------------
        double ymax = 0., ymin = 1e300;
        for (const auto& s : series){
            if (!s.h) continue;
            for (int b = 1; b <= s.h->GetNbinsX(); ++b){
                const double y = s.h->GetBinContent(b);
                if (y > ymax) ymax = y;
                if (y > 0. && y < ymin) ymin = y;
            }
        }
        if (ymin > 1e299) ymin = (ymax > 0. ? ymax * 1e-6 : 1e-6);
        // HEADROOM. The top strip is left deliberately empty so ROOT's auto-placement has a gap
        // to put the legend in; with the old x3 the box was pushed onto the top frame line and
        // onto the outermost markers.
        // Room at BOTH ends: the top strip is for the auto-placed legend, and a bottom gap is
        // what stops ROOT squeezing the box onto the lowest curve when the top happens to be taken.
        const double ylo = ymin * 0.25;
        const double yhi = logy ? ymax * 6.0 : ymax * 1.35;
        if (logy) series[first].h->GetYaxis()->SetRangeUser(ylo, yhi);
        else      series[first].h->GetYaxis()->SetRangeUser(0., yhi);
        if (logy) McDataComprRatio::ApplyLogYLabelPolicy(series[first].h, ylo, yhi);
        McDataComprRatio::HideXAxis(series[first].h);   // the ratio pad owns the x axis

        series[first].h->Draw("E");
        for (size_t i = 0; i < series.size(); ++i)
            if (series[i].h && (int)i != first) series[i].h->Draw("E,same");

        // The sign goes in the legend HEADER; each entry names its own selection, so the two pads
        // cannot be confused for the same MC object.
        TLegend* l = new TLegend();
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(42);
        l->SetMargin(0.2);
        l->SetTextSize(0.055);
        l->SetHeader(signTitles[ksign].c_str());
        for (const auto& s : series) if (s.h) l->AddEntry(s.h, s.legend.c_str(), "lp");
        l->Draw();
        McDataComprRatio::FixLegendTopOverlap(l, pad_main);

        // ================================ RATIO PAD =============================================
        pad_ratio->cd();
        gPad->SetLogx(obs.logx);

        if (!ratios.empty()){
            double rlo = 0., rhi = 2.; bool rlog = false;
            McDataComprRatio::AutoRatioRange(ratios, rlo, rhi, rlog);
            printf("[ratio] signal %s %s : y-range %.4g .. %.4g%s\n",
                   kin.c_str(), signTitles[ksign].c_str(), rlo, rhi, rlog ? " (log)" : "");
            gPad->SetLogy(rlog);
            // hist_helper set font 43 (PIXEL sizes) on these histograms, so the ratio-pad text is
            // already the same physical size as the main pad's: rel_scale = 1.
            McDataComprRatio::StyleRatioFrame(ratios[0], obs.x_title, rlo, rhi, rlog, 1.);
            ratios[0]->Draw("E");
            for (size_t i = 1; i < ratios.size(); ++i) ratios[i]->Draw("E,same");
            McDataComprRatio::DrawUnityLine(ratios[0]);
            ratios[0]->Draw("E,same");   // markers on top of the reference line
        }
    }

    c->SaveAs(OutputPath(kin + "_mc_data_compr.png").c_str());
    c->Close();
    delete c;
    for (TH1D* h : owned) delete h;
}


// `wp` = "tight" (nominal) or "medium" (WP systematic). Repo rule: every plot set exposes a
// Medium/Tight config var and defaults to Tight. It switches the DATA input file only -- the MC
// histograms here are truth quantities and have no reconstruction working point. See
// McDataComprConfig::MuonWP.
void plot_mc_data_compr_signal(const char* wp = "tight"){
    const McDataComprConfig::MuonWP muon_wp = McDataComprConfig::ParseWP(wp);
    // The five 1D signal-region observables. pair_pt and pair_pt_in_eta_subplots come from
    // plot_mc_data_pair_pt_in_eta.cxx (both projected from the SAME 2D histogram).
    for (const char* kin : {"DR_zoomin", "Deta_zoomin", "Dphi_zoomin", "minv_zoomin", "pair_eta"}){
        PlotMCDataComprSignalKinematics p(kin);
        p.muon_wp = muon_wp;
        p.Run();
    }
}
