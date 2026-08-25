// =================================================================================================
// plot_mc_data_compr.cxx -- the GENERIC family of the pp24 MC-vs-data comparison.
//
// No signal-region cut on either side: all OS (resp. SS) pairs, fiducial gap cut only. Data,
// Pythia (pp24-condition fullsim FULL sample, truth quantities) and POWHEG (truth, bb only)
// therefore describe the same region, and the panels are read for SHAPES.
//
// The absolute level here is NOT a cross-section: the data-side efficiencies were measured only
// inside the single-b signal region and are extrapolated (eps_reco clamped at the edge cells,
// eps_dR identically 1 outside its grid), the MC carries no data resonance veto, and -- the
// largest effect of all -- the MC is TRUTH with no reconstruction or trigger requirement, so at
// small dR it is a different object entirely (MC/data ~ 1e3 in the first DR_zoomin bin). Those
// warnings live in plots/mc_data_compr/generic/README.md -- deliberately NOT on the plots, which
// are held to publication standard.
//
// RATIO PAD. Every panel carries an MC/data ratio pad, with the Pythia and (where it exists)
// POWHEG curves each divided by the data. POWHEG is in the ratio because its 2026-08-25 rework
// put it on the SAME registered axes as the data and Pythia -- verified edge by edge -- so it now
// passes the same axis assertion as every other curve. Errors: TH1::Divide without "B"
// (independent samples). The ratio pad switches to a log scale by itself where the span demands
// it, which on DR_zoomin it does.
//
// Output: /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/generic/
// Usage:  root -l -b -q 'plot_mc_data_compr.cxx+()'           (Tight, nominal)
//         root -l -b -q 'plot_mc_data_compr.cxx+("medium")'   (Medium WP systematic)
// =================================================================================================

#include "PlotMCDataComprBaseClass.c"
#include <TLatex.h>

class PlotMCDataComprSingleKinematics : public PlotMCDataComprBaseClass{
protected:
    // One drawable series. Same structure as the signal family's, so the ratio-pad code below is
    // the same code -- the two families used to diverge here for no reason.
    struct Series {
        int         idt;        // DataType, for the normalization and the legend title
        std::string key;
        TH1D*       h = nullptr;
    };

public:
    bool logx = false;
    // These are 1D DISTRIBUTIONS, not efficiencies, so the repo rule is a log y axis. It is not
    // cosmetic here: on the JACOBIAN-CORRECTED dR panels the Pythia dR -> 0 divergence sets the
    // linear scale and pushes data flat onto zero, so the comparison the plot exists to make is
    // unreadable. Off by default so the pre-existing linear panels are unchanged unless a caller
    // asks. `logy` is inherited from PlotMCDataComprBaseClass (set_legend_position needs it).
    bool jacobian_correct = false;

    // Neither the Pythia fullsim nor the POWHEG truth file holds a 1/dR-weighted histogram, so on
    // a jacobian panel the ONLY 1/dR-weighted curve is the data one. Drawing an unweighted MC
    // curve beside it is the very defect (F1) this rework exists to remove, so both MC samples
    // are suppressed there and the panel is data-only. Set from jacobian_correct in Run().
    bool draw_mc = true;

    float powheg_scale = 1.;
    float pythia_scale = 1.;

    PlotMCDataComprSingleKinematics(std::string kin_input){
                kin = kin_input;
                kin_title = Observable(kin_input).x_title;
                output_subdir = "generic";
    }

    ~PlotMCDataComprSingleKinematics(){}

    void Run() override;
};


void PlotMCDataComprSingleKinematics::Run(){
    initialize();

    const McDataObservable& obs = Observable(kin);
    if (jacobian_correct) draw_mc = false;

    // 700 px tall, not 500: the bottom 30 % is the ratio pad.
    TCanvas* c = new TCanvas("c","c",1200,700);
    c->Divide(s_nSigns,1);

    // The 1/dR-weighted panels are a DIFFERENT quantity from the unweighted ones beside them, and
    // the y title has to say so -- the two used to be byte-identical, which made the two PNGs
    // indistinguishable from the image alone.
    const std::string ytitle = YTitle(obs, jacobian_correct);

    std::vector<TH1D*> owned;

    for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){
        // A data-only panel (the 1/dR-weighted ones) has nothing to take a ratio OF, so it keeps
        // the full cell height instead of being given an empty ratio pad.
        TVirtualPad* cell = c->cd(ksign + 1);
        TPad* pad_main = nullptr;
        TPad* pad_ratio = nullptr;
        if (draw_mc) {
            McDataComprRatio::SplitPadForRatio(cell, pad_main, pad_ratio);
        } else {
            // The un-split cell needs the same margins SplitPadForRatio gives its sub-pads:
            // hist_helper puts the y title at a font-43 offset that does not fit in ROOT's
            // default 10 % left margin, and the axis labels get clipped with it.
            cell->SetLeftMargin(0.23);
            cell->SetRightMargin(0.04);
            cell->SetTopMargin(0.06);
            cell->SetBottomMargin(0.135);
        }

        // --- the series of this pad, in draw order --------------------------------------------
        // 2026-08-17: the DATA histograms are no longer taken with the legacy `_wgapcut` suffix.
        // That suffix applies `PassSingleMuonGapCut` -- |eta| < eta_gap_cut1 = 0.135 at all pT
        // plus charge_eta_gap_cuts only below 6 GeV -- which is a DIFFERENT and wider object than
        // the analysis fiducial cut. Since the detector-gap fiducial cut
        // (ParamsSet::single_mu_fiducial_gap_cuts, both muons) moved into the base pp24
        // selection, the PLAIN histograms already carry exactly the analysis fiducial region.
        // (user decision 2026-08-17; docs/tracking/pp24_crossx_rerun_2026_08.md)
        std::vector<Series> series;
        if (draw_mc && !obs.powheg_var.empty())
            series.push_back({DataType::powheg, PowhegKey(obs, ksign), nullptr});
        if (draw_mc)
            series.push_back({DataType::pythia, GenericMcKey(obs, ksign, jacobian_correct), nullptr});
        series.push_back({DataType::pp_2024_2mu4, GenericDataKey(obs, ksign, jacobian_correct), nullptr});
        const int i_data = (int)series.size() - 1;

        for (auto& s : series){
            try {
                s.h = GetHist1D(s.idt, s.key);
                owned.push_back(s.h);
            } catch (const std::exception& e) {
                std::cerr << "[SKIP] " << e.what() << std::endl;
            }
        }

        // Data, Pythia AND POWHEG are booked from the SAME registered binning; a divergence is a
        // producer bug, not something to plot silently. POWHEG was exempt from this check until
        // its 2026-08-25 rework, when its dR axis was 100 bins against the data's 40.
        for (int i = 0; i < i_data; ++i)
            if (series[i_data].h && series[i].h)
                AssertSameAxis1D(series[i_data].h, series[i].h,
                                 "generic " + kin + " (" + series[i_data].key + " vs "
                                 + series[i].key + ")");

        // --- normalize & style ------------------------------------------------------------------
        for (auto& s : series){
            if (!s.h) continue;
            if (s.idt == DataType::powheg && powheg_scale != 1.) s.h->Scale(powheg_scale);
            if (s.idt == DataType::pythia && pythia_scale != 1.) s.h->Scale(pythia_scale);
            // norm_factor[pythia] is the nb -> pb conversion of the AMI weight and NOTHING else;
            // norm_factor[powheg] is exactly 1, its weight already being a pb cross-section.
            hist_helper(s.h, norm_factor[s.idt], false, signTitles[ksign].c_str(), ytitle);
            s.h->SetMarkerColor(colors[s.idt]);
            s.h->SetLineColor(colors[s.idt]);
            s.h->GetXaxis()->SetTitle(obs.x_title.c_str());
        }

        int first = -1;
        for (size_t i = 0; i < series.size(); ++i) if (series[i].h) { first = (int)i; break; }
        if (first < 0){
            std::cerr << "[WARN] No histograms available for sign " << ksign << std::endl;
            continue;
        }

        // --- ratios, built BEFORE the main pad hides its x axis -----------------------------------
        // Order matters: McDataComprRatio::MakeRatio CLONES the numerator, so a ratio built after
        // HideXAxis() inherits the suppressed x labels and title and the ratio pad comes out with
        // no x axis at all. Here the first drawn series IS the numerator of the first ratio, so
        // getting this order wrong is not hypothetical.
        std::vector<TH1*> ratios;
        if (series[i_data].h){
            for (int i = 0; i < i_data; ++i){
                if (!series[i].h) continue;
                TH1D* r = McDataComprRatio::MakeRatio(
                    series[i].h, series[i_data].h,
                    Form("r_%s%s_%u_%d", kin.c_str(),
                         jacobian_correct ? "_jac" : "", ksign, i));
                if (r){ ratios.push_back(r); owned.push_back(r); }
            }
        }

        // ================================ MAIN PAD ==============================================
        if (pad_main) pad_main->cd(); else cell->cd();
        gPad->SetLogx(logx || obs.logx);
        gPad->SetLogy(logy);

        double ymax = 0., ymin = 1e300;
        for (const auto& s : series){
            if (!s.h) continue;
            for (int b = 1; b <= s.h->GetNbinsX(); ++b){
                const double y = s.h->GetBinContent(b);
                if (y > ymax) ymax = y;
                // On a log axis the floor cannot be 0. Use the smallest positive bin content
                // across the drawn series, so no curve is silently clipped off the bottom.
                if (y > 0. && y < ymin) ymin = y;
            }
        }
        if (ymin > 1e299) ymin = (ymax > 0. ? ymax * 1e-6 : 1e-6);
        // HEADROOM for the auto-placed legend (see the signal family for the same rule).
        // Room at BOTH ends: the top strip is for the auto-placed legend, and a bottom gap is
        // what stops ROOT squeezing the box onto the lowest curve when the top happens to be
        // taken (generic/Deta_zoomin SS).
        const double ylo = logy ? ymin * 0.25 : 0.;
        const double yhi = logy ? ymax * 6.0 : ymax * 1.35;
        series[first].h->GetYaxis()->SetRangeUser(ylo, yhi);
        if (logy) McDataComprRatio::ApplyLogYLabelPolicy(series[first].h, ylo, yhi);
        // The ratio pad owns the x axis -- but only when there IS one.
        if (pad_ratio) McDataComprRatio::HideXAxis(series[first].h);

        series[first].h->Draw("E");
        for (size_t i = 0; i < series.size(); ++i)
            if (series[i].h && (int)i != first) series[i].h->Draw("E,same");

        // Legend placement. An UNSET position (s_NULL) used to build a TLegend from garbage
        // coordinates, which is how the entries ended up sitting on the markers. A
        // default-constructed TLegend is auto-placed into a gap the pad finds itself, which is
        // the only thing that works across these very differently shaped distributions (dR peaks
        // at the top, dphi is U-shaped); FixLegendTopOverlap then keeps it off the frame line.
        const std::array<float,4>& lp =
            (ksign == 0) ? legend_position_same_sign : legend_position_opp_sign;
        TLegend* l = (lp[0] < -999.) ? new TLegend() : new TLegend(lp[0], lp[1], lp[2], lp[3]);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(42);
        l->SetMargin(0.2);
        l->SetTextColor(1);
        l->SetTextSize(0.055);
        // The sign goes in the legend HEADER, not as an extra entry and not as a free-floating
        // label. As an entry it grew the box downward onto the markers; as a fixed-coordinate
        // label it could collide with the auto-placed legend. As the header it travels WITH the
        // box, so whatever gap ROOT finds holds both.
        l->SetHeader(signTitles[ksign].c_str());
        for (const auto& s : series) if (s.h) l->AddEntry(s.h, dtTitles[s.idt].c_str(), "lp");
        l->Draw();
        McDataComprRatio::FixLegendTopOverlap(l, gPad);

        // ================================ RATIO PAD =============================================
        if (pad_ratio && !ratios.empty()){
            pad_ratio->cd();
            gPad->SetLogx(logx || obs.logx);
            double rlo = 0., rhi = 2.; bool rlog = false;
            McDataComprRatio::AutoRatioRange(ratios, rlo, rhi, rlog);
            printf("[ratio] generic %s%s %s : y-range %.4g .. %.4g%s\n",
                   kin.c_str(), jacobian_correct ? "_jacobian_corrected" : "",
                   signTitles[ksign].c_str(), rlo, rhi, rlog ? " (log)" : "");
            gPad->SetLogy(rlog);
            // hist_helper set font 43 (PIXEL sizes), so the ratio-pad text is already the same
            // physical size as the main pad's: rel_scale = 1.
            McDataComprRatio::StyleRatioFrame(ratios[0], obs.x_title, rlo, rhi, rlog, 1.);
            ratios[0]->Draw("E");
            for (size_t i = 1; i < ratios.size(); ++i) ratios[i]->Draw("E,same");
            McDataComprRatio::DrawUnityLine(ratios[0]);
            ratios[0]->Draw("E,same");
        }
    }

    std::string jac_save_suffix = jacobian_correct ? "_jacobian_corrected" : "";

    if (powheg_scale != 1. || pythia_scale != 1.){
        c->SaveAs(OutputPath(Form("%s_mc_data_compr%s_powheg_%.2f_pythia_%.2f.png",
                        kin.c_str(), jac_save_suffix.c_str(), powheg_scale, pythia_scale)).c_str());
    } else {
        c->SaveAs(OutputPath(Form("%s_mc_data_compr%s.png",
                        kin.c_str(), jac_save_suffix.c_str())).c_str());
    }

    c->Close();
    delete c;
    for (TH1D* h : owned) delete h;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

class PlotMCDataComprClass{
protected:
    ParamsSet pms;
    McDataComprConfig::MuonWP muon_wp = McDataComprConfig::MuonWP::Tight;

    // LOG y on every one of these: they are 1D DISTRIBUTIONS, not efficiencies, and the repo
    // standard is a log axis for those. It is not cosmetic here -- each panel spans 1-6 decades,
    // and on a linear axis the tallest series sets the scale and flattens the others onto zero
    // (worst on the jacobian-corrected dR panels, where the data vanished entirely).
    void One(const std::string& kin, bool jacobian = false){
        PlotMCDataComprSingleKinematics p(kin);
        p.jacobian_correct = jacobian;
        p.logy = true;
        p.muon_wp = muon_wp;
        p.Run();
    }

public:
    PlotMCDataComprClass(McDataComprConfig::MuonWP wp = McDataComprConfig::MuonWP::Tight)
        : muon_wp(wp) {}
    ~PlotMCDataComprClass(){}

    void plot_mc_data_compr_generic();
};


void PlotMCDataComprClass::plot_mc_data_compr_generic(){
    One("DR");
    One("DR", /*jacobian=*/true);
    One("DR_zoomin");
    One("DR_zoomin", /*jacobian=*/true);
    One("Dphi");
    One("Dphi_zoomin");
    One("Deta_zoomin");
    One("minv_zoomin");
}


// `wp` = "tight" (nominal) or "medium" (WP systematic). See McDataComprConfig::MuonWP.
void plot_mc_data_compr(const char* wp = "tight"){
    PlotMCDataComprClass myobj(McDataComprConfig::ParseWP(wp));
    myobj.plot_mc_data_compr_generic();
}
