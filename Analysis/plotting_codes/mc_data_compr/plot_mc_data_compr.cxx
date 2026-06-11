#include "PlotMCDataComprBaseClass.c"

class PlotMCDataComprSingleKinematics : public PlotMCDataComprBaseClass{
protected:
    TCanvas* c;
    TH1D* h[s_nDtTypes][s_nSigns];
    bool h_exists[s_nDtTypes][s_nSigns];

public:
    bool norm_unity = false;
    bool logx = false;
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

        TLegend* l;
        if (ksign == 0) l = new TLegend(legend_position_same_sign[0], legend_position_same_sign[1], legend_position_same_sign[2], legend_position_same_sign[3]);
        else            l = new TLegend(legend_position_opp_sign[0], legend_position_opp_sign[1], legend_position_opp_sign[2], legend_position_opp_sign[3]);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(42);
        l->SetMargin(0.2);
        l->SetTextColor(1);

        for (unsigned int idt = 0; idt < s_nDtTypes; idt++){
            bool use_gapcut = is_data[idt];
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

        // Pythia scale
        if (h_exists[2][ksign]) {
            h[2][ksign]->Scale(pow(10,6));
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
        h[first_drawn][ksign]->GetYaxis()->SetRangeUser(0, ylim);

        // Draw: first one with "E", rest with "E,same"
        h[first_drawn][ksign]->Draw("E");
        for (int idt = 0; idt < s_nDtTypes; idt++){
            if (idt == 1) continue; // powheg_cc
            if (idt == first_drawn) continue;
            if (h_exists[idt][ksign]) h[idt][ksign]->Draw("E,same");
        }

        l->AddEntry("", (signTitles[ksign]).c_str(), "");
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
    PlotMCDataComprSingleKinematics p_DR("DR", "#Delta R");
    p_DR.legend_position_same_sign = {0.62,0.6,0.95,0.89};
    p_DR.Run();

    PlotMCDataComprSingleKinematics p_DR_jac("DR", "#Delta R");
    p_DR_jac.jacobian_correct = true;
    p_DR_jac.Run();

    PlotMCDataComprSingleKinematics p_Dphi("Dphi", "#Delta #phi");
    p_Dphi.Run();
}

void PlotMCDataComprClass::plot_mc_data_compr_DR_zoomin(){
    PlotMCDataComprSingleKinematics p_DR_zoomin("DR_zoomin", "#Delta R");
    p_DR_zoomin.Run();

    PlotMCDataComprSingleKinematics p_DR_zoomin_jac("DR_zoomin", "#Delta R");
    p_DR_zoomin_jac.jacobian_correct = true;
    p_DR_zoomin_jac.Run();
}

void PlotMCDataComprClass::plot_mc_data_compr_1D_normalize_to_unity(){
    PlotMCDataComprSingleKinematics p_DR_unity("DR", true, "#Delta R");
    p_DR_unity.Run();

    PlotMCDataComprSingleKinematics p_Dphi_unity("Dphi", true, "#Delta #phi");
    p_Dphi_unity.Run();
}


void plot_mc_data_compr(){
    PlotMCDataComprClass myobj;
    myobj.plot_mc_data_compr_1D();
    myobj.plot_mc_data_compr_DR_zoomin();
    myobj.plot_mc_data_compr_1D_normalize_to_unity();
}
