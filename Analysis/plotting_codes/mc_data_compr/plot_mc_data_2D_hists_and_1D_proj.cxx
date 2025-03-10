#include "PlotMCDataComprBaseClass.c"

class PlotMCDataSingle2DHistogram : public PlotMCDataComprBaseClass{
protected:
    TH2D* h2d_list[s_nDtTypes][s_nSigns][s_nDphis];
    std::vector<TH1D*> projY [s_nDtTypes][s_nSigns];

    std::string dt_suffix[s_nDtTypes] = {"_", "_powheg", "_pythia", "_pp2017", "_pp2024"};
    bool draw_2d_hist = true; // can ONLY plot exclusively 2D histogram or projected intervals (ensure correct scaling)

    void initialize() override;
public:

    // ------------------- public class variables -------------------

    std::vector<bool> turn_on_plotting_vec;
    std::vector<std::pair<int, int>> proj_intervals; // intervals for projectiony (use transpose if need projectionx)
    std::vector<std::string> proj_intvl_labels;
    std::vector<Color_t> proj_intvl_colors;

    std::string kiny="";
    std::string kinx="";
    std::string kiny_title=""; // used for interval projection
    bool intvls_same_canvas = true;
    bool transpose=false;
    bool logx=false;
    bool logy=false;
    bool debug_mode = false;

    // ------------------- public class methods -------------------
    PlotMCDataSingle2DHistogram(){}
    PlotMCDataSingle2DHistogram(std::string kin_in){
        kin = kin_in;
    }
    PlotMCDataSingle2DHistogram(std::string kin_in, std::vector<bool> turn_on_plotting_vec_in, std::vector<std::pair<int, int>> proj_intervals_in, std::string kiny_in, std::string kinx_in, std::string kiny_title_in):
        kiny (kiny_in), kinx (kinx_in), kiny_title (kiny_title_in){
            kin = kin_in;
            turn_on_plotting_vec = turn_on_plotting_vec_in;
            proj_intervals = proj_intervals_in;
        }

    ~PlotMCDataSingle2DHistogram(){}

    void Run() override;
};


void PlotMCDataSingle2DHistogram::initialize(){
    PlotMCDataComprBaseClass::initialize(); // call parent initialize function
    if (!proj_intervals.empty()){ // plot projected intervals: do NOT draw 2D histograms
        draw_2d_hist = false;
    }
}

void PlotMCDataSingle2DHistogram::Run(){

    initialize();
  
    // for each data type
    for (unsigned int idt = 0; idt < s_nDtTypes ; idt++){

        if (!turn_on_plotting_vec[idt]) continue; // if do NOT turn on plotting for the data type for a given observable: skip to next data type

        // ----------------------------------------------------------------------
        // retrieve all the 2D histograms without drawing
        // ----------------------------------------------------------------------

        for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){
            for (unsigned int jdphi = 0; jdphi < s_nSigns; jdphi++){
                
                std::string hist_gapcut_postfix = (is_data[idt])? "_gapcut1" : "";
                std::string hist_name = "h_" + kin + dphis[jdphi] + signs[ksign] + hist_gapcut_postfix;
                
                h2d_list[idt][ksign][jdphi] = (TH2D*) f[idt]->Get(hist_name.c_str());
                
                if (!h2d_list[idt][ksign][jdphi]){
                    std::cerr << "Hist name: " << hist_name << "; file: " << dt_paths[idt] + fnames[idt] << std::endl;
                    throw std::runtime_error("Histogram does not exist (pointer is nullptr)!");
                }
            }
        
            h2d_list[idt][ksign][0]->Add(h2d_list[idt][ksign][1]); // add up the two dphi regions for all data types
        
            if (idt == DataType::powheg_cc){
                h2d_list[idt][ksign][0]->Add(h2d_list[DataType::powheg_bb][ksign][0]); // powheg - add cc to bb for all sign & dphi regions
            }else if (idt == DataType::pythia){
                h2d_list[idt][ksign][0]->Scale(pow(10,6));
            }
        }

        // ----------------------------------------------------------------------
        // plot 2D histograms + save as png file
        // ----------------------------------------------------------------------

        if (draw_2d_hist){
            TCanvas* c = new TCanvas("c","c",1200,500);
            c->Divide(s_nSigns,1);
        
            for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){
                c->cd(ksign + 1);
          
                gPad->SetLeftMargin(0.23);
                gPad->SetBottomMargin(0.135);
                // gPad->SetTopMargin(0.18);
                gPad->SetLogx(logx);
                if (idt != DataType::powheg_bb){
                    hist_helper(h2d_list[idt][ksign][0], norm_factor[idt], false, signTitles[ksign].c_str());
                    h2d_list[idt][ksign][0]->Draw("colz");
                }
            }
        
            c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/hist_2D/%s.png", kin.c_str()));
            c->Close();
            delete c;
        }
        
        // ----------------------------------------------------------------------
        // Do interval projections if proj_intervals is not empty
        // ----------------------------------------------------------------------

        if (!proj_intervals.empty() && idt != DataType::powheg_bb) { // powheg bb already added to cc
            if (intvls_same_canvas){
                // ----------------------------------------------------------------------
                // if plot the intervals on the same canvas (default)
                // ----------------------------------------------------------------------
                TCanvas *projCanvas = new TCanvas(("projCanvas_same_canvas" + dt_suffix[idt]).c_str(), ("projCanvas_same_canvas" + dt_suffix[idt]).c_str(), 1200, 600);
                projCanvas->Divide(s_nSigns, 1);

                for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){ // different subplots (signs)
                    projCanvas->cd(ksign + 1);

                    TLegend* l;

                    if (ksign == 0) l = new TLegend(legend_position_same_sign[0], legend_position_same_sign[1], legend_position_same_sign[2], legend_position_same_sign[3]);
                    else            l = new TLegend(legend_position_opp_sign[0], legend_position_opp_sign[1], legend_position_opp_sign[2], legend_position_opp_sign[3]);

                    l->SetBorderSize(0);
                    l->SetFillStyle(0);
                    l->SetTextFont(42);
                    // l->SetTextSize(20);
                    l->SetMargin(0.2);
                    l->SetTextColor(1);

                    float ylim = 0;

                    for (size_t intvl = 0; intvl < proj_intervals.size(); ++intvl) { // different lines (proj_intervals)
                        int binX1 = proj_intervals[intvl].first;
                        int binX2 = proj_intervals[intvl].second;
                        projY[idt][ksign].push_back(h2d_list[idt][ksign][0]->ProjectionY((kiny + signs[ksign] + dt_suffix[idt] + "_" + std::to_string(intvl)).c_str(), binX1, binX2));
                                                
                        if (kiny_title == "") kiny_title = kiny;
                        std::string ytitle = Form("#frac{d#sigma}{d %s} [pb]", kiny_title.c_str());
                        // hist_helper(projY[idt][ksign][intvl], norm_factor[idt], false, signTitles[ksign].c_str(), ytitle.c_str());
                        hist_helper(projY[idt][ksign][intvl], 1., true, signTitles[ksign].c_str(), ytitle.c_str()); // normalize to unity
                        if (debug_mode) cout << "DEBUG: #Entries in projected histogram = " << projY[idt][ksign][intvl]->GetEntries() << std::endl;

                        if (proj_intvl_colors.size() > intvl){
                            projY[idt][ksign][intvl]->SetMarkerColor(proj_intvl_colors[intvl]);
                            projY[idt][ksign][intvl]->SetLineColor(proj_intvl_colors[intvl]);
                        }

                        std::string intvl_label = (proj_intvl_labels.size() > intvl)? proj_intvl_labels.at(intvl) : "Bins: " + std::to_string(binX1) + "-" + std::to_string(binX2);
                        if (debug_mode) cout << "legend same sign x start, y start, x end, y end: " << legend_position_same_sign[0] << ", " << legend_position_same_sign[1] << ", " << legend_position_same_sign[2] << ", " << legend_position_same_sign[3] << std::endl;
                        if (debug_mode) cout << "interval label: " << intvl_label << std::endl;
                        l->AddEntry(projY[idt][ksign][intvl], intvl_label.c_str(), "lp");

                        if (projY[idt][ksign][intvl]->GetMaximum() > ylim) ylim = projY[idt][ksign][intvl]->GetMaximum();
                    }

                    ylim *= 1.1;
                    projY[idt][ksign][0]->GetYaxis()->SetRangeUser(0,ylim);

                    for (size_t intvl = 0; intvl < proj_intervals.size(); ++intvl) { // different lines (proj_intervals)
                        std::string draw_opt = (intvl == 0)? "" : "same";
                        projY[idt][ksign][intvl]->Draw(draw_opt.c_str());
                    }

                    l->Draw();
                }

                std::string proj_filename;
                if (kiny != "" && kinx != "")   proj_filename = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/hist_2D_interval_projections/" + kiny + dt_suffix[idt] + "_" + kinx + "_intvls" + ".png";
                else              proj_filename = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/hist_2D_interval_projections/" + kin + "_projY_" + dt_suffix[idt] + "_intvls" + ".png";

                projCanvas->SaveAs(proj_filename.c_str());
        
                // Clean up
                projCanvas->Close();
                delete projCanvas;

                for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){
                    for (size_t intvl = 0; intvl < proj_intervals.size(); ++intvl) {
                        delete projY[idt][ksign][intvl];
                    }
                    projY[idt][ksign].clear();
                }
            }else{
                // ----------------------------------------------------------------------
                // if plot the intervals on the same canvas (default)
                // ----------------------------------------------------------------------
                for (size_t intvl = 0; intvl < proj_intervals.size(); ++intvl) {
                    int binX1 = proj_intervals[intvl].first;
                    int binX2 = proj_intervals[intvl].second;
        
                    // Project Y axis for both histograms
                    // Create canvas for projection
                    TCanvas *projCanvas = new TCanvas(("projCanvas_" + dt_suffix[idt] + std::to_string(intvl)).c_str(), ("projCanvas_" + dt_suffix[idt] + std::to_string(intvl)).c_str(), 1200, 600);
                    projCanvas->Divide(s_nSigns, 1);
        
                    
                    for (unsigned int ksign = 0; ksign < s_nSigns; ksign++){    
                        projCanvas->cd(ksign + 1);

                        projY[idt][ksign].push_back(h2d_list[idt][ksign][0]->ProjectionY((kiny + signs[ksign] + dt_suffix[idt] + "_" + std::to_string(intvl)).c_str(), binX1, binX2));
                        hist_helper(projY[idt][ksign][intvl], norm_factor[idt], false, signTitles[ksign].c_str());
                        if (debug_mode) cout << "DEBUG: #Entries in projected histogram = " << projY[idt][ksign][intvl]->GetEntries() << std::endl;
                        
                        projY[idt][ksign][intvl]->Draw();
        
                        // delete projY[idt][ksign][intvl];
                    }
        
                    // Save the projections
                    std::string proj_filename;
                    if (kiny != "" && kinx != "")   proj_filename = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/hist_2D_interval_projections/" + kiny + dt_suffix[idt] + "_" + kinx + "_intvl" + std::to_string(intvl) + ".png";
                    else              proj_filename = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/hist_2D_interval_projections/" + kin + "_projY_" + dt_suffix[idt] + "_intvl" + std::to_string(binX1) + "-" + std::to_string(binX2) + ".png";
                    projCanvas->SaveAs(proj_filename.c_str());
        
                    // Clean up
                    projCanvas->Close();
                    delete projCanvas;
                    for (unsigned int ksign = 0; ksign < s_nSigns; ksign++) delete projY[idt][ksign][intvl];
                }
                for (unsigned int ksign = 0; ksign < s_nSigns; ksign++) projY[idt][ksign].clear();
            }
        }
    }
}

void plot_mc_data_2D_hists_and_1D_proj(){
    std::vector<std::string> pair_pt_intvls_labels = {"p_{T}^{pair}: 0-4GeV", "p_{T}^{pair}: 4-8GeV", "p_{T}^{pair}: 8-12GeV", "p_{T}^{pair}: 12-16GeV", "p_{T}^{pair}: 16-20GeV"};
    std::vector<Color_t> pair_pt_intvl_colors = {kMagenta, kRed, kBlue, kBlack, kGreen+2};
    std::vector<std::string> pair_pt_intvls_labels_no_4_to_8GeV = {"p_{T}^{pair}: 8-12GeV", "p_{T}^{pair}: 12-16GeV", "p_{T}^{pair}: 16-20GeV"};
    std::vector<Color_t> pair_pt_intvl_colors_no_4_to_8GeV = {kBlue, kBlack, kGreen+2};
    // ------------------------- data -------------------------

    std::vector<std::pair<int, int>> pair_pt_intvls_data = {};
    // pair_pt_intvls_data.push_back(std::make_pair(1,20));
    // pair_pt_intvls_data.push_back(std::make_pair(21,40));
    pair_pt_intvls_data.push_back(std::make_pair(41,60));
    pair_pt_intvls_data.push_back(std::make_pair(61,80));
    pair_pt_intvls_data.push_back(std::make_pair(81,100));

    std::vector<bool> minv_proj_turn_on_vec_data {false,false,false,true,true};
    PlotMCDataSingle2DHistogram minv_distr_pairpt_intvls_data("minv_pair_pt_zoomin", minv_proj_turn_on_vec_data, pair_pt_intvls_data, "minv", "pair_pt", "m_{#mu#mu}");
    minv_distr_pairpt_intvls_data.proj_intvl_labels = pair_pt_intvls_labels_no_4_to_8GeV;
    minv_distr_pairpt_intvls_data.proj_intvl_colors = pair_pt_intvl_colors_no_4_to_8GeV;
    minv_distr_pairpt_intvls_data.debug_mode = true;
    minv_distr_pairpt_intvls_data.Run();

    // ------------------------- MC -------------------------

    std::vector<std::pair<int, int>> pair_pt_intvls_mc = {};
    // pair_pt_intvls_mc.push_back(std::make_pair(1,10))
    // pair_pt_intvls_mc.push_back(std::make_pair(11,20));
    pair_pt_intvls_mc.push_back(std::make_pair(21,30));
    pair_pt_intvls_mc.push_back(std::make_pair(31,40));
    pair_pt_intvls_mc.push_back(std::make_pair(41,50));
    std::vector<bool> minv_proj_turn_on_vec_mc {true,true,true,false,false};

    PlotMCDataSingle2DHistogram minv_distr_pairpt_intvls_mc("minv_pair_pt_zoomin", minv_proj_turn_on_vec_mc, pair_pt_intvls_mc, "minv", "pair_pt", "m_{#mu#mu}");
    minv_distr_pairpt_intvls_mc.proj_intvl_labels = pair_pt_intvls_labels_no_4_to_8GeV;
    minv_distr_pairpt_intvls_mc.proj_intvl_colors = pair_pt_intvl_colors_no_4_to_8GeV;
    minv_distr_pairpt_intvls_mc.debug_mode = true;
    minv_distr_pairpt_intvls_mc.Run();
}
