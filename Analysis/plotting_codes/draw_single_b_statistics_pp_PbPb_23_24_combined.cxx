
// function to draw single-b statistics as function of pT
// for Pb+Pb 2023 (mu4_mu4noL1) & 2024 (mu4) combined

void draw_single_b_statistics_pp_PbPb_23_24_combined() {

    std::string dataset_base_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/";
    std::string f_pp_24_single_b = "pp_2024/pp_2024_single_b_ana_hists_mu4_mu4noL1.root";
    std::string f_PbPb_23_single_b = "pbpb_2023/pbpb_2023_single_b_ana_hists_mu4_mu4noL1_no_res_cut.root";
    std::string f_PbPb_24_single_b = "pbpb_2024/pbpb_2024_single_b_ana_hists_single_mu4.root";

    // --- settings ---
    bool combine_0_5_AND_5_10 = true; // if true, combine 0-5 and 5-10 statistics
    bool subtract_pp_same_sign = true; // if true, subtract same-sign pairs from opposite-sign pairs for pp24 data
    bool plot_single_pad = true; // if true, plot single pad; if false, plot 2 pads

    std::vector<int> ctr_bin_edges = {0,5,10,20,30,50,80}; // if change ctr binning --> adjust this vector for histogram retrieving updates
    int ctr_0_5_ind = 0; // needed if combine_0_5_AND_5_10
    int pT_max = 150; // if change pT ranges --> adjust this number for histogram retrieving updates

    std::vector<Color_t> colors = {kRed, kGreen+2, kBlue, kOrange+7, kMagenta, kCyan+2, kGray}; // vector of colors

    std::map<int, double> pT_max_to_legend_lowxedge_map_single_pad;
    std::map<int, double> pT_max_to_legend_lowyedge_map_single_pad;

    std::map<int, double> pT_max_to_legend_lowxedge_map_pad1;
    std::map<int, double> pT_max_to_legend_lowyedge_map_pad1;
    std::map<int, double> pT_max_to_legend_lowxedge_map_pad2;
    std::map<int, double> pT_max_to_legend_lowyedge_map_pad2;

    pT_max_to_legend_lowxedge_map_single_pad[200] = 0.35;
    pT_max_to_legend_lowyedge_map_single_pad[200] = 0.45;
    pT_max_to_legend_lowxedge_map_single_pad[150] = 0.35;
    pT_max_to_legend_lowyedge_map_single_pad[150] = 0.45;

    pT_max_to_legend_lowxedge_map_pad1[200] = 0.35;
    pT_max_to_legend_lowyedge_map_pad1[200] = 0.45;
    pT_max_to_legend_lowxedge_map_pad2[200] = 0.33;
    pT_max_to_legend_lowyedge_map_pad2[200] = 0.45;

    pT_max_to_legend_lowxedge_map_pad1[150] = 0.35;
    pT_max_to_legend_lowyedge_map_pad1[150] = 0.53;
    pT_max_to_legend_lowxedge_map_pad2[150] = 0.33;
    pT_max_to_legend_lowyedge_map_pad2[150] = 0.45;

    // --- open files ---
    TFile *fpp24 = TFile::Open((dataset_base_dir + f_pp_24_single_b).c_str());
    TFile *fPbPb23 = TFile::Open((dataset_base_dir + f_PbPb_23_single_b).c_str());
    TFile *fPbPb24 = TFile::Open((dataset_base_dir + f_PbPb_24_single_b).c_str());
    if (!fpp24 || fpp24->IsZombie() || !fPbPb23 || fPbPb23->IsZombie() || !fPbPb24 || fPbPb24->IsZombie()) {
        std::cerr << "Error opening input files!" << std::endl;
        return;
    }

    // --- retrieve pp histograms ---

    std::string h_name_pp_ss =  "h_pT_" + to_string(pT_max) + "_ss";
    std::string h_name_pp_op =  "h_pT_" + to_string(pT_max) + "_op";
    
    TH1D* h_pp24_ss = (TH1D*)fpp24->Get(h_name_pp_ss.c_str());
    TH1D* h_pp24 = (TH1D*)fpp24->Get(h_name_pp_op.c_str());

    if (!h_pp24 && !h_pp24_ss){
        std::cerr << "Missing pp 24 histogram: " << h_name_pp_op << std::endl;
        return;
    }

    if (subtract_pp_same_sign) h_pp24->Add(h_pp24_ss, -1);

    h_pp24->SetLineColor(kBlack);
    h_pp24->SetLineStyle(kDashed);
    h_pp24->SetMarkerColor(kBlack);
    h_pp24->SetLineWidth(2);

    // --- vectors for centrality-binned Pb+Pb histograms ---
    std::vector<TH1D*> v_PbPb_23_ss, v_PbPb_23_op, v_PbPb_24_ss, v_PbPb_24_op, v_stats_PbPb_23_24_comb;

    // --- retrieve histograms in each Pb+Pb centrality bin, perform bkg subtraction & add 23, 24 statistics ---
    for (int i = 0; i < ctr_bin_edges.size()-1; ++i) {
        std::string h_name_ss =  "hss_pT_" + to_string(pT_max) + "_counts_ctr" + to_string(ctr_bin_edges.at(i)) + "_" + to_string(ctr_bin_edges.at(i+1)) + "_w_signal_cuts";
        std::string h_name_op =  "hop_pT_" + to_string(pT_max) + "_counts_ctr" + to_string(ctr_bin_edges.at(i)) + "_" + to_string(ctr_bin_edges.at(i+1)) + "_w_signal_cuts";

        v_PbPb_23_ss.push_back((TH1D*)fPbPb23->Get(h_name_ss.c_str()));
        v_PbPb_23_op.push_back((TH1D*)fPbPb23->Get(h_name_op.c_str()));
        v_PbPb_24_ss.push_back((TH1D*)fPbPb24->Get(h_name_ss.c_str()));
        v_PbPb_24_op.push_back((TH1D*)fPbPb24->Get(h_name_op.c_str()));

        if (!v_PbPb_23_ss.back() || !v_PbPb_23_op.back()) {
            std::cerr << "Missing Pb+Pb 23 histogram: " << h_name_ss << " or " << h_name_op << std::endl;
            return;
        }

        if (!v_PbPb_24_ss.back() || !v_PbPb_24_op.back()) {
            std::cerr << "Missing Pb+Pb 24 histogram: " << h_name_ss << " or " << h_name_op << std::endl;
            return;
        }

        // subtract same-sign pairs (with the same minv & pair pT cuts for single-b signal pairs) as estimated combinatoric background
        // add 2023 + 2024 Pb+Pb statistics in each centrality bin
        TH1D *h_stats_23_24_comb = (TH1D*)v_PbPb_23_op.at(i)->Clone(Form("h%d_%d_stats_23_24_comb", ctr_bin_edges.at(i), ctr_bin_edges.at(i+1)));
        h_stats_23_24_comb->Add(v_PbPb_23_ss.at(i), -1);
        h_stats_23_24_comb->Add(v_PbPb_24_op.at(i), 1);
        h_stats_23_24_comb->Add(v_PbPb_24_ss.at(i), -1);
        h_stats_23_24_comb->SetLineColor(colors.at(i));
        h_stats_23_24_comb->SetMarkerColor(i + 1);
        h_stats_23_24_comb->SetLineWidth(2);
        v_stats_PbPb_23_24_comb.push_back(h_stats_23_24_comb);
    }

    // --- draw canvas ---
    TCanvas *c;
    if (plot_single_pad){
        c = new TCanvas("c", "Comparison Results", 700, 500);
        c->cd(1);
        gPad->SetLogy();

        TLegend *leg = new TLegend(pT_max_to_legend_lowxedge_map_single_pad[pT_max], pT_max_to_legend_lowyedge_map_single_pad[pT_max], 0.9, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);

        h_pp24->SetStats(0);
        h_pp24->SetTitle("");
        h_pp24->GetXaxis()->SetTitle("p_{T} [GeV]");
        h_pp24->GetYaxis()->SetTitle("Counts per bin");
        h_pp24->Draw("ELP");

        for (int i = 0; i < v_stats_PbPb_23_24_comb.size(); ++i) {

            if (combine_0_5_AND_5_10){
                if (i == ctr_0_5_ind) v_stats_PbPb_23_24_comb.at(i)->Add(v_stats_PbPb_23_24_comb.at(i+1)); // add 5-10 centrality bin to 0-5 centrality bin
                else if (i == ctr_0_5_ind + 1) continue; // do not plot for centrality 5-10%
            }
            
            v_stats_PbPb_23_24_comb.at(i)->Draw("ELP SAME");

            if (combine_0_5_AND_5_10 && i == ctr_0_5_ind)
                leg->AddEntry(v_stats_PbPb_23_24_comb.at(i), "Centrality: 0-10%", "lp");
            else
                leg->AddEntry(v_stats_PbPb_23_24_comb.at(i), Form("Centrality: %d-%d%%", ctr_bin_edges.at(i), ctr_bin_edges.at(i+1)), "lp");
        }
        leg->AddEntry("", "pp (mu4_mu4noL1)", "");
        leg->AddEntry("", "Pb+Pb, 2023 (mu4_mu4noL1), 2024 (mu4) combined", "");
        leg->AddEntry("", "(opposite sign - same sign) pairs", "");
        leg->AddEntry("", "medium muon WP", "");
        leg->AddEntry("", "m_{#mu#mu} #in (1.08, 2.9) GeV, p_{T}^{pair} > 8 GeV", "");
        leg->Draw();
        if (subtract_pp_same_sign)
            c->SaveAs(Form("%ssingle_b_stats_pp_PbPb_23_24_combined_pT%d_subtr_pp_same_sign.png", dataset_base_dir.c_str(), pT_max));
        else
            c->SaveAs(Form("%ssingle_b_stats_pp_PbPb_23_24_combined_pT%d.png", dataset_base_dir.c_str(), pT_max));

        
    } else{
        c = new TCanvas("c", "Comparison Results", 1400, 500);
        c->Divide(2, 1);

        // --- pad 1: h1-h3 ---
        c->cd(1);
        gPad->SetLogy();

        TLegend *leg1 = new TLegend(pT_max_to_legend_lowxedge_map_pad1[pT_max], pT_max_to_legend_lowyedge_map_pad1[pT_max], 0.9, 0.88);
        leg1->SetFillStyle(0);
        leg1->SetBorderSize(0);

        h_pp24->SetStats(0);
        h_pp24->SetTitle("");
        h_pp24->GetXaxis()->SetTitle("p_{T} [GeV]");
        h_pp24->GetYaxis()->SetTitle("Counts per bin");
        h_pp24->Draw("ELP");

        for (int i = 0; i < v_stats_PbPb_23_24_comb.size()/2; ++i) {

            if (combine_0_5_AND_5_10){
                if (i == ctr_0_5_ind) v_stats_PbPb_23_24_comb.at(i)->Add(v_stats_PbPb_23_24_comb.at(i+1)); // add 5-10 centrality bin to 0-5 centrality bin
                else if (i == ctr_0_5_ind + 1) continue; // do not plot for centrality 5-10%
            }
            
            v_stats_PbPb_23_24_comb.at(i)->Draw("ELP SAME");

            if (combine_0_5_AND_5_10 && i == ctr_0_5_ind)
                leg1->AddEntry(v_stats_PbPb_23_24_comb.at(i), "Centrality: 0-10%", "lp");
            else
                leg1->AddEntry(v_stats_PbPb_23_24_comb.at(i), Form("Centrality: %d-%d%%", ctr_bin_edges.at(i), ctr_bin_edges.at(i+1)), "lp");
        }
        leg1->AddEntry("", "pp (mu4_mu4noL1)", "");
        leg1->AddEntry("", "medium muon WP", "");
        leg1->AddEntry("", "m_{#mu#mu} #in (1.08, 2.9) GeV, p_{T}^{pair} > 8 GeV", "");
        leg1->Draw();

        // --- pad 2: h4-h6 ---
        c->cd(2);
        gPad->SetLogy();

        TLegend *leg2 = new TLegend(pT_max_to_legend_lowxedge_map_pad2[pT_max], pT_max_to_legend_lowyedge_map_pad2[pT_max], 0.88, 0.88);
        leg2->SetFillStyle(0);
        leg2->SetBorderSize(0);

        h_pp24->SetStats(0);
        h_pp24->SetTitle("");
        h_pp24->GetXaxis()->SetTitle("p_{T} [GeV]");
        h_pp24->GetYaxis()->SetTitle("Counts per bin");
        h_pp24->Draw("ELP");

        for (int i = v_stats_PbPb_23_24_comb.size()/2; i < v_stats_PbPb_23_24_comb.size(); ++i) {
            v_stats_PbPb_23_24_comb.at(i)->Draw("ELP SAME");
            leg2->AddEntry(v_stats_PbPb_23_24_comb.at(i), Form("Centrality: %d-%d%%", ctr_bin_edges.at(i), ctr_bin_edges.at(i+1)), "lp");
        }
        leg2->AddEntry("", "Pb+Pb, 2023 (mu4_mu4noL1), 2024 (mu4) combined", "");
        leg2->AddEntry("", "(opposite sign - same sign) pairs", "");
        leg2->Draw();

        if (subtract_pp_same_sign)
            c->SaveAs(Form("%ssingle_b_stats_pp_PbPb_23_24_combined_pT%d_double_panel_subtr_pp_same_sign.png", dataset_base_dir.c_str(), pT_max));
        else
            c->SaveAs(Form("%ssingle_b_stats_pp_PbPb_23_24_combined_pT%d_double_panel.png", dataset_base_dir.c_str(), pT_max));
    }

}
