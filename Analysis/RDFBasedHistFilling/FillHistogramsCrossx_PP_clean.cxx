// FillHistogramsCrossx: Single B → J/psi crossx measurements (trigger_mode==2)
// ============================================================================

void RDFBasedHistFillingPP::FillHistogramsCrossx(){
    std::cout << "[PP] FillHistogramsCrossx: Building opposite-sign dataframe with signal cuts" << std::endl;
    
    // Signal cuts: minv > 1.08 && minv < 2.9 && pair_pt > 8 && pair_eta < 2.2 && dr > 0.05
    const std::string signal_cuts = "minv > 1.08 && minv < 2.9 && pair_pt > 8 && pair_eta < 2.2 && dr > 0.05";
    
    // Get opposite-sign dataframe (use df_op_weighted if available, else df_op)
    ROOT::RDF::RNode df_op_rdf;
    if (df_map.find("df_op_weighted") != df_map.end()) {
        df_op_rdf = df_map.at("df_op_weighted");
    } else {
        df_op_rdf = map_at_checked(df_map, "df_op", "FillHistogramsCrossx: df_op");
    }
    
    // Apply signal cuts
    ROOT::RDF::RNode df_crossx = df_op_rdf.Filter(signal_cuts);
    
    std::cout << "[PP] Filling 2D crossx histograms..." << std::endl;
    
    // 2D: pair_pt vs pair_eta
    std::string h_name = "h2d_crossx_pair_pt_pair_eta";
    if (hist2D_map.find(h_name) == hist2D_map.end()) {
        hist2D_map[h_name] = new TH2D(h_name.c_str(), "crossx: pair pT vs pair eta", 
                                       50, 8, 80, 44, -2.4, 2.4);
    }
    df_crossx.Foreach([this, h_name](float pt, float eta) {
        map_at_checked(hist2D_map, h_name, "FillHistogramsCrossx")->Fill(pt, eta);
    }, {"pair_pt", "pair_eta"});
    
    // 2D: pair_pt vs minv
    h_name = "h2d_crossx_pair_pt_minv";
    if (hist2D_map.find(h_name) == hist2D_map.end()) {
        hist2D_map[h_name] = new TH2D(h_name.c_str(), "crossx: pair pT vs minv", 
                                       50, 8, 80, 50, 1, 3);
    }
    df_crossx.Foreach([this, h_name](float pt, float minv) {
        map_at_checked(hist2D_map, h_name, "FillHistogramsCrossx")->Fill(pt, minv);
    }, {"pair_pt", "minv"});
    
    // 2D: pair_pt vs dr
    h_name = "h2d_crossx_pair_pt_dr";
    if (hist2D_map.find(h_name) == hist2D_map.end()) {
        hist2D_map[h_name] = new TH2D(h_name.c_str(), "crossx: pair pT vs dR", 
                                       50, 8, 80, 50, 0.05, 3);
    }
    df_crossx.Foreach([this, h_name](float pt, float dr) {
        map_at_checked(hist2D_map, h_name, "FillHistogramsCrossx")->Fill(pt, dr);
    }, {"pair_pt", "dr"});
    
    std::cout << "[PP] FillHistogramsCrossx completed" << std::endl;
}
