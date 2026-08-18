#include "PlotMCDataComprBaseClass.c"
#include <TLatex.h>

class PlotMCDataComprSingleKinematics : public PlotMCDataComprBaseClass{
protected:
    TCanvas* c;
    TH1D* h[s_nDtTypes][s_nSigns];
    bool h_exists[s_nDtTypes][s_nSigns];

public:
    bool norm_unity = false;
    bool logx = false;
    // These are 1D DISTRIBUTIONS, not efficiencies, so the repo rule is a log y axis. It is not
    // cosmetic here: on the JACOBIAN-CORRECTED dR panels the Pythia dR -> 0 divergence sets the
    // linear scale and pushes data and POWHEG flat onto zero, so the comparison the plot exists to
    // make is unreadable. Off by default so the pre-existing linear panels are unchanged unless a
    // caller asks.
    // `logy` is inherited from PlotMCDataComprBaseClass (set_legend_position needs it).
    bool jacobian_correct = false;
    float powheg_scale = 1.;
    float pythia_scale = 1.;

    PlotMCDataComprSingleKinematics(std::string kin_input, std::string kin_title_input){
                kin = kin_input;
                kin_title = kin_title_input;
    }

    PlotMCDataComprSingleKinematics(std::string kin_input, bool norm_unity_input, std::string kin_title_input): norm_unity(norm_unity_input){
                kin = kin_input;
                kin_title = kin_title_input;
    }

    ~PlotMCDataComprSingleKinematics(){}

    void Run() override;
};


void PlotMCDataComprSingleKinematics::Run(){
    initialize();

    TCanvas* c = new TCanvas("c","c",1200,500);
    c->Divide(s_nSigns,1);

    for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){
        c->cd(ksign + 1);
        gPad->SetLeftMargin(0.23);
        gPad->SetBottomMargin(0.135);
        gPad->SetLogx(logx);
        gPad->SetLogy(logy);

        // Legend placement. An UNSET position (s_NULL) used to build a TLegend from garbage
        // coordinates, which is how the entries ended up sitting on the markers. Fall back to a
        // corner the data does not reach: with a log y axis the curves fill the upper half, so
        // the legend goes low; on a linear axis the upper right is empty.
        const std::array<float,4>& lp =
            (ksign == 0) ? legend_position_same_sign : legend_position_opp_sign;
        // A still-unset position means "let ROOT place it": a default-constructed TLegend is
        // auto-placed into a gap the pad finds itself, which is the only thing that works across
        // these very differently shaped distributions (dR peaks at the top, dphi is U-shaped).
        TLegend* l = (lp[0] < -999.) ? new TLegend() : new TLegend(lp[0], lp[1], lp[2], lp[3]);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(42);
        l->SetMargin(0.2);
        l->SetTextColor(1);
        l->SetTextSize(0.045);

        for (unsigned int idt = 0; idt < s_nDtTypes; idt++){
            // 2026-08-17: the DATA histograms are no longer taken with the legacy `_wgapcut`
            // suffix. That suffix applies `PassSingleMuonGapCut` -- |eta| < eta_gap_cut1 = 0.135
            // at all pT plus charge_eta_gap_cuts only below 6 GeV -- which is a DIFFERENT and
            // wider object than the analysis fiducial cut. Since the detector-gap fiducial cut
            // (ParamsSet::single_mu_fiducial_gap_cuts, both muons) moved into the base pp24
            // selection, the PLAIN histograms already carry exactly the analysis fiducial region,
            // the same one the cross-section measures; `_wgapcut` on top would be a selection
            // used nowhere else. (user decision 2026-08-17;
            // docs/tracking/pp24_crossx_rerun_2026_08.md)
            bool use_gapcut = false;
            std::string hist_name = BuildHistName(idt, kin, ksign, jacobian_correct, use_gapcut);

            try {
                h[idt][ksign] = GetHist1D(idt, hist_name);
                h_exists[idt][ksign] = true;
            } catch (const std::exception& e) {
                std::cerr << "[SKIP] " << e.what() << std::endl;
                h[idt][ksign] = nullptr;
                h_exists[idt][ksign] = false;
            }

            if (h_exists[idt][ksign]) {
                h[idt][ksign]->SetMarkerColor(colors[idt]);
                h[idt][ksign]->SetLineColor(colors[idt]);
            }
        }

        // POWHEG: add cc to bb
        if (h_exists[DataType::powheg_bb][ksign] && h_exists[DataType::powheg_cc][ksign]) {
            h[DataType::powheg_bb][ksign]->Add(h[DataType::powheg_cc][ksign]);
        }

        // Scale and draw
        if (powheg_scale != 1. && h_exists[0][ksign]) h[0][ksign]->Scale(powheg_scale);

        std::string ytitle = Form("#frac{d#sigma}{d %s} [pb]", kin_title.c_str());

        // Find first existing histogram to draw first
        int first_drawn = -1;

        // Draw POWHEG (bb+cc) if exists
        if (h_exists[0][ksign]) {
            hist_helper(h[0][ksign], norm_factor[0], norm_unity, signTitles[ksign].c_str(), ytitle);
            l->AddEntry(h[0][ksign], dtTitles[0].c_str(), "lp");
            first_drawn = 0;
        }

        // Pythia scale to bring the truth Pythia dsigma onto the data [pb] scale.
        // NOTE (units, 2026-06-30): this 1e6 is the product of TWO factors:
        //   (i)  NB_TO_PB = 1e3 — the stored MC `weight` is in nb (ATLAS AMI `crossSection`
        //        is nb, NOT pb — see PythiaAlgCoreT.c AMI-read comment); data here is pb.
        //   (ii) ~1e3 — a SEPARATE, pre-existing under-normalization of the TRUTH COMBINED
        //        sample (`histograms_pythia_combined*`): its per-pTHat-slice combine left the
        //        absolute scale ~1e3 low (stored weights span 0..3.5e6). This is immaterial to
        //        the analysis (templates are area-normalized, k=G_SS/G_OS is a ratio) and is
        //        tracked separately (low_mass_dimuon_template_fit.md). It must be fixed at the
        //        truth-combine step before this hand factor can be reduced to the clean 1e3.
        // Until the truth combine is fixed, DO NOT change this value — it keeps the known-good
        // data-vs-Pythia overlay on scale.
        static const double NB_TO_PB             = 1.0e3;
        static const double TRUTH_COMBINE_RENORM = 1.0e3;  // pre-existing truth-combine under-norm (see above)
        if (h_exists[2][ksign]) {
            h[2][ksign]->Scale(NB_TO_PB * TRUTH_COMBINE_RENORM);
            if (pythia_scale != 1.) h[2][ksign]->Scale(pythia_scale);
            hist_helper(h[2][ksign], norm_factor[2], norm_unity, signTitles[ksign].c_str(), ytitle);
            l->AddEntry(h[2][ksign], dtTitles[2].c_str(), "lp");
            if (first_drawn < 0) first_drawn = 2;
        }

        // Data
        if (h_exists[3][ksign]) {
            hist_helper(h[3][ksign], norm_factor[3], norm_unity, signTitles[ksign].c_str(), ytitle);
            l->AddEntry(h[3][ksign], dtTitles[3].c_str(), "lp");
            if (first_drawn < 0) first_drawn = 3;
        }

        if (first_drawn < 0) {
            std::cerr << "[WARN] No histograms available for sign " << ksign << std::endl;
            continue;
        }

        // Y-axis range
        float ylim = 0;
        for (int idt = 0; idt < s_nDtTypes; idt++){
            if (idt == 1) continue; // powheg_cc added to bb
            if (h_exists[idt][ksign] && h[idt][ksign]->GetMaximum() > ylim)
                ylim = h[idt][ksign]->GetMaximum();
        }
        ylim *= 1.1;
        if (logy) {
            // On a log axis the floor cannot be 0. Use the smallest positive bin content across
            // the drawn series, so no curve is silently clipped off the bottom.
            double ymin = 1e300;
            for (int idt = 0; idt < s_nDtTypes; idt++){
                if (idt == 1 || !h_exists[idt][ksign]) continue;
                for (int b = 1; b <= h[idt][ksign]->GetNbinsX(); ++b) {
                    const double y = h[idt][ksign]->GetBinContent(b);
                    if (y > 0. && y < ymin) ymin = y;
                }
            }
            if (ymin > 1e299) ymin = ylim * 1e-6;
            h[first_drawn][ksign]->GetYaxis()->SetRangeUser(ymin * 0.5, ylim * 2.0);
        } else {
            h[first_drawn][ksign]->GetYaxis()->SetRangeUser(0, ylim);
        }

        // Draw: first one with "E", rest with "E,same"
        h[first_drawn][ksign]->Draw("E");
        for (int idt = 0; idt < s_nDtTypes; idt++){
            if (idt == 1) continue; // powheg_cc
            if (idt == first_drawn) continue;
            if (h_exists[idt][ksign]) h[idt][ksign]->Draw("E,same");
        }

        // The sign goes in the legend HEADER, not as an extra entry and not as a free-floating
        // label. As an entry it grew the box downward onto the markers; as a fixed-coordinate
        // label it could collide with the auto-placed legend. As the header it travels WITH the
        // box, so whatever gap ROOT finds holds both.
        l->SetHeader(signTitles[ksign].c_str());
        l->Draw();
    }

    std::string unity_suffix = norm_unity ? "_unity" : "";
    std::string jac_save_suffix = jacobian_correct ? "_jacobian_corrected" : "";

    if (powheg_scale != 1. || pythia_scale != 1.){
        c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/%s_mc_data_compr%s_powheg_%.2f_pythia_%.2f.png",
                        kin.c_str(), jac_save_suffix.c_str(), powheg_scale, pythia_scale));
    } else {
        c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/%s_mc_data_compr%s%s.png",
                        kin.c_str(), jac_save_suffix.c_str(), unity_suffix.c_str()));
    }

    c->Close();
    delete c;

    for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){
        for (unsigned int idt = 0; idt < s_nDtTypes; idt++){
            if (h_exists[idt][ksign]) delete h[idt][ksign];
        }
    }
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

class PlotMCDataComprClass{
protected:
    ParamsSet pms;
public:
    PlotMCDataComprClass(){}
    ~PlotMCDataComprClass(){}

    void plot_mc_data_compr_1D();
    void plot_mc_data_compr_DR_zoomin();
    void plot_mc_data_compr_1D_normalize_to_unity();
};


void PlotMCDataComprClass::plot_mc_data_compr_1D(){
    // LOG y on every one of these: they are 1D DISTRIBUTIONS, not efficiencies, and the repo
    // standard is a log axis for those. It is not cosmetic here -- each panel spans 1-6 decades,
    // and on a linear axis the tallest series sets the scale and flattens the others onto zero
    // (worst on the jacobian-corrected dR panels, where data and POWHEG vanished entirely).
    // With logy set, the legend default moves to the lower left, clear of the curves
    // (PlotMCDataComprBaseClass::set_legend_position).
    PlotMCDataComprSingleKinematics p_DR("DR", "#Delta R");
    p_DR.logy = true;
    p_DR.Run();

    PlotMCDataComprSingleKinematics p_DR_jac("DR", "#Delta R");
    p_DR_jac.jacobian_correct = true;
    p_DR_jac.logy = true;
    p_DR_jac.Run();

    PlotMCDataComprSingleKinematics p_Dphi("Dphi", "#Delta #phi");
    p_Dphi.logy = true;
    p_Dphi.Run();

    // Added 2026-08-18 (user): these three exist on BOTH sides under the same name, so the only
    // thing that was missing was the data histogram -- now filled by the extended generic 1D list
    // in RDFBasedHistFillingData::BuildFilterToVarListMapDataCommon. POWHEG carries none of them
    // (its combined file has only DR with mechanism suffixes), so it is skipped with a [SKIP]
    // line, exactly as it already is for Dphi.
    PlotMCDataComprSingleKinematics p_Dphi_zoomin("Dphi_zoomin", "#Delta #phi");
    p_Dphi_zoomin.logy = true;
    p_Dphi_zoomin.Run();

    PlotMCDataComprSingleKinematics p_Deta_zoomin("Deta_zoomin", "#Delta #eta");
    p_Deta_zoomin.logy = true;
    p_Deta_zoomin.Run();

    PlotMCDataComprSingleKinematics p_minv_zoomin("minv_zoomin", "m_{#mu#mu} [GeV]");
    p_minv_zoomin.logy = true;
    p_minv_zoomin.Run();
}

void PlotMCDataComprClass::plot_mc_data_compr_DR_zoomin(){
    PlotMCDataComprSingleKinematics p_DR_zoomin("DR_zoomin", "#Delta R");
    p_DR_zoomin.logy = true;
    p_DR_zoomin.Run();

    PlotMCDataComprSingleKinematics p_DR_zoomin_jac("DR_zoomin", "#Delta R");
    p_DR_zoomin_jac.jacobian_correct = true;
    p_DR_zoomin_jac.logy = true;
    p_DR_zoomin_jac.Run();
}

void PlotMCDataComprClass::plot_mc_data_compr_1D_normalize_to_unity(){
    PlotMCDataComprSingleKinematics p_DR_unity("DR", true, "#Delta R");
    p_DR_unity.logy = true;
    p_DR_unity.Run();

    PlotMCDataComprSingleKinematics p_Dphi_unity("Dphi", true, "#Delta #phi");
    p_Dphi_unity.logy = true;
    p_Dphi_unity.Run();
}


void plot_mc_data_compr(){
    PlotMCDataComprClass myobj;
    myobj.plot_mc_data_compr_1D();
    myobj.plot_mc_data_compr_DR_zoomin();
    myobj.plot_mc_data_compr_1D_normalize_to_unity();
}
