#include "plot_mc_data_compr.c"

TH2D* h2d[nDtTypes][nSigns][nDphis];
std::vector<TH1D*> projY [nDtTypes][nSigns];

std::string dt_suffix[nDtTypes] = {"_", "_powheg", "_pythia", "_pp2017", "_pp2024"};


void plot_mc_data_single_2D_histogram(std::string kin, std::vector<bool> turn_on_plotting, const std::vector<std::pair<int, int>>& intervals, std::string kiny="", std::string kin_title="", bool intvls_same_canvas = false, bool transpose=false, bool logx=false, bool logy=false){

    initialize();
  
    for (unsigned int idt = 0; idt < nDtTypes ; idt++){
  
        // plot 2D histograms
        TCanvas* c = new TCanvas("c","c",1200,500);
        c->Divide(nSigns,1);
    
        for (unsigned int ksign = 0; ksign < nSigns; ksign++){
            c->cd(ksign + 1);
      
            gPad->SetLeftMargin(0.23);
            gPad->SetBottomMargin(0.135);
            // gPad->SetTopMargin(0.18);
            gPad->SetLogx(logx);
      
            for (unsigned int jdphi = 0; jdphi < nSigns; jdphi++){
                if (!turn_on_plotting[idt]) continue;
                std::string hist_gapcut_postfix = (is_data[idt])? "_gapcut1" : "";
                std::string hist_name = "h_" + kin + dphis[jdphi] + signs[ksign] + hist_gapcut_postfix;
                
                h2d[idt][ksign][jdphi] = (TH2D*) f[idt]->Get(hist_name.c_str());
                
                if (!h2d[idt][ksign][jdphi]){
                    std::cerr << "Hist name: " << hist_name << "; file: " << dt_paths[idt] + fnames[idt] << std::endl;
                    throw std::runtime_error("Histogram does not exist (pointer is nullptr)!");
                }
            }
        
            h2d[idt][ksign][0]->Add(h2d[idt][ksign][1]); // add up the two dphi regions for all data types
        
            if (idt == 1){
                h2d[idt][ksign][0]->Add(h2d[0][ksign][0]); // powheg - add cc to bb for all sign & dphi regions
            }else if (idt == 2){
                h2d[idt][ksign][0]->Scale(pow(10,6));
            }

            if (idt != 0){
                if (kin_title == "") kin_title = kiny;
                hist_helper(h2d[idt][ksign][0], norm_factor[idt], false, signTitles[ksign].c_str());
                h2d[idt][ksign][0]->Draw("colz");
            }
        }
    
        c->SaveAs(Form("/usatlas/u/yuhanguo/workarea/dimuon_codes/plots/hist_2D/%s.png", kin.c_str()));
        c->Close();
        delete c;
    
    
        // Process projectionY if there are intervals
        if (!intervals.empty() && idt != 0) {
            for (size_t intvl = 0; intvl < intervals.size(); ++intvl) {
                int binX1 = intervals[intvl].first;
                int binX2 = intervals[intvl].second;
    
                // Project Y axis for both histograms
                // Create canvas for projection
                TCanvas *projCanvas = new TCanvas(("projCanvas_" + dt_suffix[idt] + std::to_string(intvl)).c_str(), ("projCanvas_" + dt_suffix[idt] + std::to_string(intvl)).c_str(), 1200, 600);
                projCanvas->Divide(nSigns, 1);
    
                
                for (unsigned int ksign = 0; ksign < nSigns; ksign++){    
                    projCanvas->cd(ksign + 1);

                    projY[idt][ksign].push_back(h2d[idt][ksign][0]->ProjectionY((kiny + signs[ksign] + dt_suffix[idt] + "_" + std::to_string(intvl)).c_str(), binX1, binX2));
                    cout << projY[idt][ksign][intvl]->GetEntries() << std::endl;
                    hist_helper(projY[idt][ksign][intvl], norm_factor[idt], false, signTitles[ksign].c_str());
                    
                    projY[idt][ksign][intvl]->Draw();
    
                    // delete projY[idt][ksign][intvl];
                }
    
                // Save the projections
                std::string proj_filename;
                if (kiny != "")   proj_filename = "/usatlas/u/yuhanguo/workarea/dimuon_codes/plots/hist_2D/" + kiny + dt_suffix[idt] + "_intvl" + std::to_string(intvl) + ".png";
                else              proj_filename = "/usatlas/u/yuhanguo/workarea/dimuon_codes/plots/hist_2D/" + kin + "_projY_" + dt_suffix[idt] + "_intvl" + std::to_string(intvl) + ".png";
                projCanvas->SaveAs(proj_filename.c_str());
    
                // Clean up
                projCanvas->Close();
                delete projCanvas;
                for (unsigned int ksign = 0; ksign < nSigns; ksign++) delete projY[idt][ksign][intvl];
            }
        }
    }
}

void plot_mc_data_2D_hists_and_1D_proj(){
    // std::vector<std::pair<int, int>> pair_pt_intvls_data = {};
    // pair_pt_intvls_data.push_back(std::make_pair(1,20));
    // pair_pt_intvls_data.push_back(std::make_pair(21,40));
    // pair_pt_intvls_data.push_back(std::make_pair(41,60));
    // pair_pt_intvls_data.push_back(std::make_pair(61,80));
    // pair_pt_intvls_data.push_back(std::make_pair(81,100));

    // std::vector<bool> minv_proj_turn_on_vec {false,false,false,true,true};
    // plot_mc_data_single_2D_histogram("minv_pair_pt_zoomin", minv_proj_turn_on_vec, pair_pt_intvls_data, "minv", "m_{#mu#mu} versus p_{T}^{pair}");
    
    std::vector<std::pair<int, int>> pair_pt_intvls_mc = {};
    pair_pt_intvls_mc.push_back(std::make_pair(1,10));
    pair_pt_intvls_mc.push_back(std::make_pair(11,20));
    pair_pt_intvls_mc.push_back(std::make_pair(21,30));
    pair_pt_intvls_mc.push_back(std::make_pair(31,40));
    pair_pt_intvls_mc.push_back(std::make_pair(41,50));
    std::vector<bool> minv_proj_turn_on_vec {true,true,true,false,false};
    plot_mc_data_single_2D_histogram("minv_pair_pt_zoomin", minv_proj_turn_on_vec, pair_pt_intvls_mc, "minv", "m_{#mu#mu} versus p_{T}^{pair}");
}
