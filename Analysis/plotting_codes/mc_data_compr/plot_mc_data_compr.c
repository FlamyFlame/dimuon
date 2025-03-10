#include "PlotMCDataComprBaseClass.c"


class PlotMCDataComprSingleKinematics : public PlotMCDataComprBaseClass{
protected:
    TCanvas* c;
    TH1D* h[s_nDtTypes][s_nSigns][s_nDphis];

public:

    // ------------------- public class variables -------------------

    bool projx_2d = false;
    bool projy_2d = false;

    bool norm_unity = false;
    bool logx=false;
    bool jacobian_correct = false;
    float powheg_scale = 1.;
    float pythia_scale = 1.;
    bool near_away_sepr = false; // plot and save near region (dPhi < pi/2) and away region (dPhi >= pi/2) in separate png files

    bool ratio_plot_mode = false;
    bool ratio_plot_draw_errorbar = true;

    // ------------------- public class methods -------------------

    // constructor for "normal" 1D histogram
    PlotMCDataComprSingleKinematics(std::string kin_input, std::string kin_title_input){
                kin = kin_input;
                kin_title = kin_title_input;
    }

    // convenient constructor for 1D histogram requiring normalization to unity
    PlotMCDataComprSingleKinematics(std::string kin_input, bool norm_unity_input, std::string kin_title_input): norm_unity (norm_unity_input){
                kin = kin_input;
                kin_title = kin_title_input;
    }

    // constructor for 2D histogram projecting onto an 1D histogram (only projx OR projy can be set true)
    PlotMCDataComprSingleKinematics(std::string kin_input, bool projx_2d_input, bool projy_2d_input, std::string kin1d_input, std::string kin_title_input): 
                projx_2d (projx_2d_input), projy_2d (projy_2d_input){
                kin = kin_input;
                kin1d = kin1d_input;
                kin_title = kin_title_input;
    }

    // convenient constructor for 2D histogram projecting onto an 1D histogram requiring normalization to unity (only projx OR projy can be set true)
    PlotMCDataComprSingleKinematics(std::string kin_input, bool projx_2d_input, bool projy_2d_input, bool norm_unity_input, std::string kin1d_input, std::string kin_title_input): 
                projx_2d (projx_2d_input), projy_2d (projy_2d_input), norm_unity (norm_unity_input){
                kin = kin_input;
                kin1d = kin1d_input;
                kin_title = kin_title_input;
    }
    ~PlotMCDataComprSingleKinematics(){}

    void Run() override;
};

// ------------------------------------------------------------------------------------------------------------------------------------------------------


// void plot_mc_data_compr_single_kinematic(std::string kin, bool projx_2d, bool projy_2d, bool norm_unity, std::string kin1d, std::string kin_title, bool logx=false, bool jacobian_correct=false, float powheg_scale = 1., float pythia_scale = 1.){
void PlotMCDataComprSingleKinematics::Run(){

    initialize();

    if (kin1d == "") kin1d = kin;

    if (projx_2d && projy_2d){
        std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
        throw std::exception();
    }

    if (norm_unity && (powheg_scale != 1. || pythia_scale != 1.)){
        std::cout << "Cannot both normalize to unity and scale powheg/pythia distributions." << std::endl;
        throw std::exception();
    }

    if (powheg_scale <= 0 || pythia_scale <= 0){
        std::cout << "Powheg & pythia scales must both be positive." << std::endl;
        throw std::exception();
    }
    
    TCanvas* c = new TCanvas("c","c",1200,500);
    c->Divide(s_nSigns,1);

    for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){
        c->cd(ksign + 1);

        gPad->SetLeftMargin(0.23);
        gPad->SetBottomMargin(0.135);
        // gPad->SetTopMargin(0.18);
        gPad->SetLogx(logx);


        // TLegend* l = new TLegend(0.24,0.6,0.5,0.89);
        TLegend* l;

        if (ksign == 0) l = new TLegend(legend_position_same_sign[0], legend_position_same_sign[1], legend_position_same_sign[2], legend_position_same_sign[3]);
        else            l = new TLegend(legend_position_opp_sign[0], legend_position_opp_sign[1], legend_position_opp_sign[2], legend_position_opp_sign[3]);

        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(42);
        // l->SetTextSize(20);
        l->SetMargin(0.2);
        l->SetTextColor(1);

        for (unsigned int jdphi = 0; jdphi < s_nSigns; jdphi++){
            for (unsigned int idt = 0; idt < s_nDtTypes; idt++){
                std::string hist_gapcut_postfix = (is_data[idt])? "_gapcut1" : "";
                std::string hist_name = "h_" + kin + dphis[jdphi] + signs[ksign] + hist_gapcut_postfix;
                
                if (projx_2d){
                    TH2D* h2d = (TH2D*) f[idt]->Get(hist_name.c_str());
                    if (norm_unity){
                        h[idt][ksign][jdphi] = (TH1D*) h2d->ProjectionX(Form("hx_unity_%d%d%d",idt+1,ksign+1,jdphi+1));
                    }else{
                        h[idt][ksign][jdphi] = (TH1D*) h2d->ProjectionX(Form("hx_%d%d%d",idt+1,ksign+1,jdphi+1));
                    }
                }else if (projy_2d){
                    TH2D* h2d = (TH2D*) f[idt]->Get(hist_name.c_str());
                    if (norm_unity){
                        h[idt][ksign][jdphi] = (TH1D*) h2d->ProjectionY(Form("hy_unity_%d%d%d",idt+1,ksign+1,jdphi+1));
                    }else{
                        h[idt][ksign][jdphi] = (TH1D*) h2d->ProjectionY(Form("hy_%d%d%d",idt+1,ksign+1,jdphi+1));
                    }
                }else{ // 1D
                    h[idt][ksign][jdphi] = (TH1D*) f[idt]->Get(hist_name.c_str());
                }
                if (!h[idt][ksign][jdphi]){
                    std::cerr << "Hist name: " << hist_name << "; file: " << dt_paths[idt] + fnames[idt] << std::endl;
                    throw std::runtime_error("Histogram does not exist (pointer is s_NULLptr)!");
                }

                h[idt][ksign][jdphi]->SetMarkerColor(colors[idt]);
                h[idt][ksign][jdphi]->SetLineColor(colors[idt]);
            }

            h[DataType::powheg_bb][ksign][jdphi]->Add(h[DataType::powheg_cc][ksign][jdphi]); // powheg - add cc to bb for all sign & dphi regions
        }

        for (unsigned int idt = 0; idt < s_nDtTypes; idt++){
            h[idt][ksign][0]->Add(h[idt][ksign][1]); // add up the two dphi regions for all data types
        }

        if (powheg_scale != 1.) h[0][ksign][0]->Scale(powheg_scale);

        std::string ytitle = ratio_plot_mode? Form("#frac{d#sigma}{d %s} / #frac{d#sigma_{pp2017}}{d %s}", kin_title.c_str(), kin_title.c_str()) : Form("#frac{d#sigma}{d %s} [pb]", kin_title.c_str());
        hist_helper(h[0][ksign][0], norm_factor[0], norm_unity, signTitles[ksign].c_str(), ytitle);
        l->AddEntry(h[0][ksign][0], dtTitles[0].c_str(), "lp");

        h[2][ksign][0]->Scale(pow(10,6));
        if (pythia_scale != 1.) h[2][ksign][0]->Scale(pythia_scale);

        for (unsigned int idt = 2; idt < s_nDtTypes; idt++){
            hist_helper(h[idt][ksign][0], norm_factor[idt], norm_unity, signTitles[ksign].c_str(), ytitle);
            l->AddEntry(h[idt][ksign][0], dtTitles[idt].c_str(), "lp");
        }

        if (ratio_plot_mode){
            try{
                // set correct ratio-plot bin content for all but the reference histogram
                for (unsigned int idt = 0; idt < s_nDtTypes; idt++){
                    if (idt == DataType::pp_2017) continue; // this is required for the datasets following the reference histogram to have correct bin content for the ratio plot

                    // find ratio of #bins(reference histogram) to #bins(current histogram) - currently only allow ratio to be 1 or 2

                    int nbins_refr_to_current_ratio;
                    if (h[idt][ksign][0]->GetNbinsX() == h[DataType::pp_2017][ksign][0]->GetNbinsX()){
                        nbins_refr_to_current_ratio = 1;
                    } else if (h[DataType::pp_2017][ksign][0]->GetNbinsX() == h[idt][ksign][0]->GetNbinsX() * 2){
                        nbins_refr_to_current_ratio = 2;
                    }else{
                        std::cerr << "Required to perform ratio plot for kinematic " << kin1d << " but histogram has different #bins from the reference histogram!" << std::endl;
                        std::cerr << "Current dt index: " << idt << "; #bins for kinematic " << kin1d << " in current dt is: " << h[idt][ksign][0]->GetNbinsX() << std::endl;
                        std::cerr << "#bins for kinematic " << kin1d << "in pp2017 (reference) data is: " << h[DataType::pp_2017][ksign][0]->GetNbinsX() << std::endl;
                        throw std::runtime_error("Required to perform ratio plot but histogram has different #bins from the reference histogram!");
                    }

                    for (int bin = 1; bin <= h[idt][ksign][0]->GetNbinsX(); bin++){ // loop over the bins
                        double bin_content_current = h[idt][ksign][0]->GetBinContent(bin);
                        double bin_err_current = h[idt][ksign][0]->GetBinError(bin);

                        double bin_content_refr = (nbins_refr_to_current_ratio == 1)? h[DataType::pp_2017][ksign][0]->GetBinContent(bin) : h[DataType::pp_2017][ksign][0]->GetBinContent(2*bin-1) + h[DataType::pp_2017][ksign][0]->GetBinContent(2*bin);
                        double bin_err_refr = (nbins_refr_to_current_ratio == 1)? h[DataType::pp_2017][ksign][0]->GetBinError(bin) : h[DataType::pp_2017][ksign][0]->GetBinError(2*bin-1) + h[DataType::pp_2017][ksign][0]->GetBinContent(2*bin);

                        if (bin_content_refr == 0){ // for current bin, #events in reference histogram is zero
                            h[idt][ksign][0]->SetBinContent(bin, s_NULL); // ratio is invalid - set to be s_NULL
                            h[idt][ksign][0]->SetBinError(bin, 0); // error of ratio is invalid - set to be zero (to avoid plot distortion)
                            continue; // skip to next bin
                        }

                        // set bin content
                        double ratio = bin_content_current * 1. / bin_content_refr;
                        h[idt][ksign][0]->SetBinContent(bin, ratio);

                        // set bin error
                        if (bin_content_current == 0){
                            h[idt][ksign][0]->SetBinError(bin, 0);
                        }else{ // both the current histogram & the reference histogram have nonzero bin content --> propogate error
                            if (ratio_plot_draw_errorbar){
                                double error_ratio = ratio * sqrt( pow(bin_err_current / bin_content_current, 2) + pow(bin_err_refr / bin_content_refr, 2) );
                                h[idt][ksign][0]->SetBinError(bin, error_ratio);                                
                            }else{
                                h[idt][ksign][0]->SetBinError(bin, 0);
                            }
                        }
                    }
                }

                // set correct ratio-plot bin content for the reference histogram
                for (int bin = 1; bin <= h[DataType::pp_2017][ksign][0]->GetNbinsX(); bin++){
                    h[DataType::pp_2017][ksign][0]->SetBinContent(bin, 1.);
                    double bin_content_pp17 = h[DataType::pp_2017][ksign][0]->GetBinContent(bin);
                    double bin_err_pp17 = h[DataType::pp_2017][ksign][0]->GetBinError(bin);
                    double relative_err_pp17 = (bin_content_pp17 == 0)? 0 : bin_err_pp17 * 1. / bin_content_pp17;
                    // h[DataType::pp_2017][ksign][0]->SetBinError(bin, sqrt(2) * relative_err_pp17); 
                    h[DataType::pp_2017][ksign][0]->SetBinError(bin, 0.); // ratio-plot error bar for pp17 looks suspicious --> set error to be zero for now
                }
            }catch (const std::runtime_error& err){
                // std::cout << "WARNING: Runtime error caught when attempting ratio plot for the kinematic " << kin1d << ": " << err.what() << std::endl;
                std::cout << "WARNING: Runtime error caught when attempting ratio plot for the kinematic " << kin1d << std::endl;
                std::cout << "Return without drawing ratio plot!" << std::endl;
                return;
            }
        }

        float ylim = (h[0][ksign][0]->GetMaximum() > h[2][ksign][0]->GetMaximum())? h[0][ksign][0]->GetMaximum() : h[2][ksign][0]->GetMaximum();
        ylim = (ylim > h[3][ksign][0]->GetMaximum())? ylim : h[3][ksign][0]->GetMaximum();
        ylim = (ylim > h[4][ksign][0]->GetMaximum())? ylim : h[4][ksign][0]->GetMaximum();
        ylim *= 1.1;
        
        if (ratio_plot_mode)    h[0][ksign][0]->GetYaxis()->SetRangeUser(-0.1,ylim);
        else                    h[0][ksign][0]->GetYaxis()->SetRangeUser(0,ylim);

        h[0][ksign][0]->Draw("E");

        for (unsigned int idt = 2; idt < s_nDtTypes; idt++){
            h[idt][ksign][0]->Draw("E,same");
        }

        l->AddEntry("",(signTitles[ksign]).c_str(),"");
        if (ratio_plot_mode){
            if (ratio_plot_draw_errorbar)   l->AddEntry("","Ratio plot","");
            else                            l->AddEntry("","Ratio plot; error bar not drawn","");
        }
        l->Draw();
    }

    std::string unity_suffix = norm_unity? "_unity" : "";
    std::string ratio_suffix = ratio_plot_mode? "_ratio_to_pp17" : "";
    std::string ratio_errorbar_suffix = ratio_plot_draw_errorbar? "" : "_noerrorbar";

    if (powheg_scale != 1. || pythia_scale != 1.){ // can't coexist with normalizing-to-unity requirement
        c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/%s_mc_data_compr%s%s%_powheg_%.2f_pythia_%.2f.png", kin1d.c_str(), ratio_suffix.c_str(), ratio_errorbar_suffix.c_str(), powheg_scale, pythia_scale));
    }else{
        c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/%s_mc_data_compr%s%s%%s.png", kin1d.c_str(), ratio_suffix.c_str(), ratio_errorbar_suffix.c_str(), unity_suffix.c_str()));
    }

    c->Close();
    delete c;

    for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){
        for (unsigned int jdphi = 0; jdphi < s_nSigns; jdphi++){
            for (unsigned int idt = 0; idt < s_nDtTypes; idt++){
                delete h[idt][ksign][jdphi];
            }
        }
    }
}


void plot_mc_data_compr_1D(){
    // ----------------- 1D histogram plotting with #events / luminosity-based normalization -----------------
    PlotMCDataComprSingleKinematics p_DR("DR", "#Delta R");
    p_DR.legend_position_same_sign = {0.62,0.6,0.95,0.89};
    p_DR.Run();

    PlotMCDataComprSingleKinematics p_DR_jacobian_corr("DR", "#Delta R");
    p_DR_jacobian_corr.jacobian_correct = true;
    p_DR_jacobian_corr.Run();

    PlotMCDataComprSingleKinematics p_Dphi("Dphi", "#Delta #phi");
    p_Dphi.Run();

    PlotMCDataComprSingleKinematics p_pt_asym("pt_asym", "A");
    p_pt_asym.Run();

    PlotMCDataComprSingleKinematics p_pair_pt_ptlead_ratio("pair_pt_ptlead_ratio", "#frac{p_{T}^{pair}}{p_{T}^{lead}}");
    p_pair_pt_ptlead_ratio.Run();
}

void plot_mc_data_compr_1D_ratio(){
    // ----------------- 1D histogram plotting with #events / luminosity-based normalization -----------------
    PlotMCDataComprSingleKinematics p_DR_ratio("DR", "#Delta R");
    p_DR_ratio.ratio_plot_mode = true;
    p_DR_ratio.Run();

    PlotMCDataComprSingleKinematics p_Dphi_ratio("Dphi", "#Delta #phi");
    p_Dphi_ratio.ratio_plot_mode = true;
    p_Dphi_ratio.Run();

    PlotMCDataComprSingleKinematics p_pt_asym_ratio("pt_asym", "A");
    p_pt_asym_ratio.ratio_plot_mode = true;
    p_pt_asym_ratio.Run();

    PlotMCDataComprSingleKinematics p_pair_pt_ptlead_ratio_ratio("pair_pt_ptlead_ratio", "#frac{p_{T}^{pair}}{p_{T}^{lead}}");
    p_pair_pt_ptlead_ratio_ratio.ratio_plot_mode = true;
    p_pair_pt_ptlead_ratio_ratio.Run();
}

void plot_mc_data_compr_1D_ratio_noerrorbar(){
    // ----------------- 1D histogram plotting with #events / luminosity-based normalization -----------------
    PlotMCDataComprSingleKinematics p_DR_ratio("DR", "#Delta R");
    p_DR_ratio.ratio_plot_mode = true;
    p_DR_ratio.ratio_plot_draw_errorbar = false;
    p_DR_ratio.Run();

    PlotMCDataComprSingleKinematics p_DR_jacobian_corr_ratio("DR", "#Delta R");
    p_DR_jacobian_corr_ratio.ratio_plot_mode = true;
    p_DR_jacobian_corr_ratio.jacobian_correct = true;
    p_DR_jacobian_corr_ratio.ratio_plot_draw_errorbar = false;
    p_DR_jacobian_corr_ratio.Run();

    PlotMCDataComprSingleKinematics p_Dphi_ratio("Dphi", "#Delta #phi");
    p_Dphi_ratio.ratio_plot_mode = true;
    p_Dphi_ratio.ratio_plot_draw_errorbar = false;
    p_Dphi_ratio.Run();

    PlotMCDataComprSingleKinematics p_pt_asym_ratio("pt_asym", "A");
    p_pt_asym_ratio.ratio_plot_mode = true;
    p_pt_asym_ratio.ratio_plot_draw_errorbar = false;
    p_pt_asym_ratio.Run();

    PlotMCDataComprSingleKinematics p_pair_pt_ptlead_ratio_ratio("pair_pt_ptlead_ratio", "#frac{p_{T}^{pair}}{p_{T}^{lead}}");
    p_pair_pt_ptlead_ratio_ratio.ratio_plot_mode = true;
    p_pair_pt_ptlead_ratio_ratio.ratio_plot_draw_errorbar = false;
    p_pair_pt_ptlead_ratio_ratio.Run();
}

void plot_mc_data_compr_DR_zoomin(){
    PlotMCDataComprSingleKinematics p_DR_zoomin("DR_zoomin", "#Delta R");
    // p_DR_zoomin.legend_position_same_sign = {0.62,0.6,0.95,0.89};
    p_DR_zoomin.Run();

    PlotMCDataComprSingleKinematics p_DR_zoomin_ratio("DR_zoomin", "#Delta R");
    p_DR_zoomin_ratio.ratio_plot_mode = true;
    p_DR_zoomin_ratio.Run();

    PlotMCDataComprSingleKinematics p_DR_zoomin_jacobian_corr("DR_zoomin", "#Delta R");
    p_DR_zoomin_jacobian_corr.jacobian_correct = true;
    p_DR_zoomin_jacobian_corr.Run();

}


void plot_mc_data_compr_2D_proj(){
    // ----------------- 2D histograms projecting onto 1D -----------------
    PlotMCDataComprSingleKinematics p_ptlead_pair_pt_x("ptlead_pair_pt", true, false, false, "pair_pt", "p_{T}^{pair}");
    p_ptlead_pair_pt_x.Run();

    PlotMCDataComprSingleKinematics p_ptlead_pair_pt_y("ptlead_pair_pt", false, true, false, "ptlead", "p_{T}^{lead}");
    p_ptlead_pair_pt_y.Run();

    PlotMCDataComprSingleKinematics p_minv_pair_pt_y("minv_pair_pt", false, true, false, "minv","m_{#mu#mu}");
    p_minv_pair_pt_y.Run();

    PlotMCDataComprSingleKinematics p_minv_pair_pt_zoomin_y("minv_pair_pt_zoomin", false, true, false, "minv_zoomin","m_{#mu#mu}");
    p_minv_pair_pt_zoomin_y.Run();

    PlotMCDataComprSingleKinematics p_Deta_Dphi_y("Deta_Dphi", false, true, false, "Deta", "#Delta #eta");
    p_Deta_Dphi_y.legend_position_same_sign = {0.62,0.6,0.95,0.89};
    p_Deta_Dphi_y.legend_position_opp_sign = {0.6,0.6,0.93,0.89};
    p_Deta_Dphi_y.Run();
}

void plot_mc_data_compr_2D_proj_ratio(){
    // ----------------- 2D histograms projecting onto 1D -----------------
    PlotMCDataComprSingleKinematics p_ptlead_pair_pt_x_ratio("ptlead_pair_pt", true, false, false, "pair_pt", "p_{T}^{pair}");
    p_ptlead_pair_pt_x_ratio.ratio_plot_mode = true;
    p_ptlead_pair_pt_x_ratio.Run();

    PlotMCDataComprSingleKinematics p_ptlead_pair_pt_y_ratio("ptlead_pair_pt", false, true, false, "ptlead", "p_{T}^{lead}");
    p_ptlead_pair_pt_y_ratio.ratio_plot_mode = true;
    p_ptlead_pair_pt_y_ratio.Run();

    PlotMCDataComprSingleKinematics p_minv_pair_pt_y_ratio("minv_pair_pt", false, true, false, "minv","m_{#mu#mu}");
    p_minv_pair_pt_y_ratio.ratio_plot_mode = true;
    p_minv_pair_pt_y_ratio.Run();

    PlotMCDataComprSingleKinematics p_minv_pair_pt_zoomin_y_ratio("minv_pair_pt_zoomin", false, true, false, "minv_zoomin","m_{#mu#mu}");
    p_minv_pair_pt_zoomin_y_ratio.ratio_plot_mode = true;
    p_minv_pair_pt_zoomin_y_ratio.Run();

    PlotMCDataComprSingleKinematics p_Deta_Dphi_y_ratio("Deta_Dphi", false, true, false, "Deta", "#Delta #eta");
    p_Deta_Dphi_y_ratio.ratio_plot_mode = true;
    p_Deta_Dphi_y_ratio.Run();
}

void plot_mc_data_compr_2D_proj_ratio_noerrorbar(){
    // ----------------- 2D histograms projecting onto 1D -----------------
    PlotMCDataComprSingleKinematics p_ptlead_pair_pt_x_ratio("ptlead_pair_pt", true, false, false, "pair_pt", "p_{T}^{pair}");
    p_ptlead_pair_pt_x_ratio.ratio_plot_mode = true;
    p_ptlead_pair_pt_x_ratio.ratio_plot_draw_errorbar = false;
    p_ptlead_pair_pt_x_ratio.Run();

    PlotMCDataComprSingleKinematics p_ptlead_pair_pt_y_ratio("ptlead_pair_pt", false, true, false, "ptlead", "p_{T}^{lead}");
    p_ptlead_pair_pt_y_ratio.ratio_plot_mode = true;
    p_ptlead_pair_pt_y_ratio.ratio_plot_draw_errorbar = false;
    p_ptlead_pair_pt_y_ratio.Run();

    PlotMCDataComprSingleKinematics p_minv_pair_pt_y_ratio("minv_pair_pt", false, true, false, "minv","m_{#mu#mu}");
    p_minv_pair_pt_y_ratio.ratio_plot_mode = true;
    p_minv_pair_pt_y_ratio.ratio_plot_draw_errorbar = false;
    p_minv_pair_pt_y_ratio.Run();

    PlotMCDataComprSingleKinematics p_minv_pair_pt_zoomin_y_ratio("minv_pair_pt_zoomin", false, true, false, "minv_zoomin","m_{#mu#mu}");
    p_minv_pair_pt_zoomin_y_ratio.ratio_plot_mode = true;
    p_minv_pair_pt_zoomin_y_ratio.ratio_plot_draw_errorbar = false;
    p_minv_pair_pt_zoomin_y_ratio.Run();

    PlotMCDataComprSingleKinematics p_Deta_Dphi_y_ratio("Deta_Dphi", false, true, false, "Deta", "#Delta #eta");
    p_Deta_Dphi_y_ratio.ratio_plot_mode = true;
    p_Deta_Dphi_y_ratio.ratio_plot_draw_errorbar = false;
    p_Deta_Dphi_y_ratio.Run();
}

void plot_mc_data_compr_1D_normalize_to_unity(){
    // ----------------- 1D histogram plotting normalized to unity -----------------
    PlotMCDataComprSingleKinematics p_DR_unity("DR", true, "#Delta R");
    p_DR_unity.Run();

    PlotMCDataComprSingleKinematics p_Dphi_unity("Dphi", true, "#Delta #phi");
    p_Dphi_unity.Run();

    PlotMCDataComprSingleKinematics p_pt_asym_unity("pt_asym", true, "A");
    p_pt_asym_unity.Run();

    PlotMCDataComprSingleKinematics p_pair_pt_ptlead_ratio_unity("pair_pt_ptlead_ratio", true, "#frac{p_{T}^{pair}}{p_{T}^{lead}}");
    p_pair_pt_ptlead_ratio_unity.Run();
}

void plot_mc_data_compr_2D_proj_normalize_to_unity(){
    // ----------------- 2D histograms projecting onto 1D normalized to unity -----------------

    PlotMCDataComprSingleKinematics p_ptlead_pair_pt_norm_x("ptlead_pair_pt", true, false, true, "pair_pt", "p_{T}^{pair}");
    p_ptlead_pair_pt_norm_x.Run();

    PlotMCDataComprSingleKinematics p_ptlead_pair_pt_norm_y("ptlead_pair_pt", false, true, true, "ptlead", "p_{T}^{lead}");
    p_ptlead_pair_pt_norm_y.Run();

    PlotMCDataComprSingleKinematics p_minv_pair_pt_norm_y("minv_pair_pt", false, true, true, "minv","m_{#mu#mu}");
    p_minv_pair_pt_norm_y.Run();

    PlotMCDataComprSingleKinematics p_minv_pair_pt_zoomin_norm_y("minv_pair_pt_zoomin", false, true, true, "minv_zoomin","m_{#mu#mu}");
    p_minv_pair_pt_zoomin_norm_y.Run();

    PlotMCDataComprSingleKinematics p_Deta_Dphi_norm_y("Deta_Dphi", false, true, true, "Deta", "#Delta #eta");
    p_Deta_Dphi_norm_y.Run();
}

void plot_mc_data_compr(){
    // plot_mc_data_compr_1D();
    // plot_mc_data_compr_2D_proj();
    // plot_mc_data_compr_1D_ratio();
    // plot_mc_data_compr_2D_proj_ratio();
    // plot_mc_data_compr_1D_ratio_noerrorbar();
    // plot_mc_data_compr_2D_proj_ratio_noerrorbar();
    plot_mc_data_compr_DR_zoomin();
    // plot_mc_data_compr_1D_normalize_to_unity();
    // plot_mc_data_compr_2D_proj_normalize_to_unity();

}


