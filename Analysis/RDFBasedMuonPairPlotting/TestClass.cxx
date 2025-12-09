void TestClass::Run(){
    Initialize();
    ProcessData();    
    Finalize();
}

void TestClass::Initialize(){
    input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024_combined_single_mu4.root");
    output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_single_mu4_new_RDF.root";
    BuildFilterToVarListMap();
    BuildHistBinningMap();
    ReadVar1DJson();
}

void TestClass::ProcessData(){
    rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree_op, input_files));
    ROOT::RDF::RNode node_op = *(rdf_store.back());
    df_map["df_op"] = node_op;
    FillHistogramsSingleDataFrame("_op", node_op);
}

void TestClass::BuildFilterToVarListMap(){
        df_filter_to_var1D_list_map["_op"]         = {"Dphi", "DR", "DR_zoomin", "minv_zoomin"};
        df_filter_to_var2D_list_map["_op"]         = {{"pair_pT_log", "minv_log"}};
        df_filter_to_var3D_list_map["_op"]         = {{"pair_pT_log", "minv_log", "DR"}};
}

void TestClass::BuildHistBinningMap(){

    // ------- pT binnings -------
    hist_binning_map["pT_bins_80"] = pms.pT_bins_80;

    // ------- pT binning for single-muon trigger efficieny -------

    std::vector<double> pT_bins_single_muon (pms.pT_bins_8); // make a copy of a suitable set of single-muon pT bins (adjustable) --> use the copy for histogram settings

    pT_bins_single_muon.insert(pT_bins_single_muon.end(), pms.pT_bins_60.begin(), pms.pT_bins_60.end());

    hist_binning_map["pT_bins_single_muon"] = pT_bins_single_muon;

    // ------- eta binning for single-muon trigger efficiency -------

    std::vector<double> eta_bins_trig_effcy = ParamsSet::makeEtaTrigEffcyBinning();
    int neta_bins_trig_effcy = eta_bins_trig_effcy.size() - 1;
    hist_binning_map["eta_bins_trig_effcy"] = eta_bins_trig_effcy;
    
    // ------- eta binning for single-muon trigger efficiency -------

    // Build uniform phi edges so we can use the (xbins, ybins, zbins) TH3D ctor
    int nphi_bins_trig_effcy = 128; // phi 2nd muon
    
    std::vector<double> phi2nd_bins(nphi_bins_trig_effcy + 1);
    for (int i = 0; i <= nphi_bins_trig_effcy; ++i) {
        phi2nd_bins[i] = -pms.PI + (2.0 * pms.PI) * (static_cast<double>(i) / nphi_bins_trig_effcy);
    }

    hist_binning_map["phi2nd_bins"] = phi2nd_bins;

    // ------- minv log binning -------

    static const int nminv_bins_log = 40;
    float minv_logpow[ParamsSet::nSigns];
    minv_logpow[0] = 0.062;
    minv_logpow[1] = 0.045;
    float minv_max = 60;
    std::vector<double>minv_bins_log[ParamsSet::nSigns];

    std::string dphi_regions[2] = {"near", "away"};

    for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
        for(int iminv = 0; iminv <= nminv_bins_log; iminv++){
            minv_bins_log[isign].push_back(minv_max * pow(10.0, (static_cast<float>(iminv - nminv_bins_log)) * minv_logpow[isign]));
        }
    }

    hist_binning_map["minv_log_bins_ss"] = minv_bins_log[0];
    hist_binning_map["minv_log_bins_op"] = minv_bins_log[1];
}

// ---------- WRITE OUTPUT & FINALIZE ----------

#include <algorithm> // std::find

// 1) For maps like std::map<std::string, ROOT::RDF::RResultPtr<THnD>>
template <typename H>
void write_hist_map_rresult(
    std::map<std::string, ROOT::RDF::RResultPtr<H>>& m,
    const std::vector<std::string>& hists_to_not_write)
{
    for (auto& kv : m) {
        const auto& name = kv.first;
        if (std::find(hists_to_not_write.begin(), hists_to_not_write.end(), name)
            != hists_to_not_write.end()) continue;

        kv.second->Write();  // RResultPtr<THnD>
    }
}

// 2) For maps like std::map<std::string, std::vector<THnD>>
template <typename H>
void write_hist_map_vector(
    std::map<std::string, std::vector<H>>& m,
    const std::vector<std::string>& hists_to_not_write)
{
    for (auto& kv : m) {
        const auto& name = kv.first;
        if (std::find(hists_to_not_write.begin(), hists_to_not_write.end(), name)
            != hists_to_not_write.end()) continue;

        for (auto& h : kv.second) {
            h.Write();        // THnD in a vector
        }
    }
}

void TestClass::Finalize(){
    TFile * m_outfile=new TFile(output_file.c_str(), file_mode.c_str());
    write_hist_map_rresult(hist1d_rresultptr_map, hists_to_not_write);
    write_hist_map_rresult(hist2d_rresultptr_map, hists_to_not_write);
    write_hist_map_rresult(hist3d_rresultptr_map, hists_to_not_write);
    delete m_outfile;
}


// ---------- HELPER FUNCTIONS ----------

void TestClass::FillHistogramsSingleDataFrame(const std::string& filter,
                                            ROOT::RDF::RNode df,
                                            bool hists_not_write = false,
                                            std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0}) {
    static const std::vector<std::string> empty1D;
    static const std::vector<std::array<std::string, 2>> empty2D;
    static const std::vector<std::array<std::string, 3>> empty3D;

    auto it1D = df_filter_to_var1D_list_map.find(filter);
    auto it2D = df_filter_to_var2D_list_map.find(filter);
    auto it3D = df_filter_to_var3D_list_map.find(filter);

    const auto& vars1D = (it1D != df_filter_to_var1D_list_map.end()) ? it1D->second : empty1D;
    const auto& vars2D = (it2D != df_filter_to_var2D_list_map.end()) ? it2D->second : empty2D;
    const auto& vars3D = (it3D != df_filter_to_var3D_list_map.end()) ? it3D->second : empty3D;

    FillHistogramsSingleDataFrame(filter, df, "", vars1D, vars2D, vars3D, hists_not_write, hists_1_2_3D_not_write);
}

void TestClass::FillHistogramsSingleDataFrame(const std::string& filter,
                                            const std::string& weight,
                                            ROOT::RDF::RNode df,
                                            bool weight_before_filter = false,
                                            bool hists_not_write = false,
                                            std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0}) {
    static const std::vector<std::string> empty1D;
    static const std::vector<std::array<std::string, 2>> empty2D;
    static const std::vector<std::array<std::string, 3>> empty3D;

    std::pair<std::string, std::string> key(filter, weight);

    auto it1D = df_filter_and_weight_to_var1D_list_map.find(key);
    auto it2D = df_filter_and_weight_to_var2D_list_map.find(key);
    auto it3D = df_filter_and_weight_to_var3D_list_map.find(key);

    const auto& vars1D = (it1D != df_filter_and_weight_to_var1D_list_map.end()) ? it1D->second : empty1D;
    const auto& vars2D = (it2D != df_filter_and_weight_to_var2D_list_map.end()) ? it2D->second : empty2D;
    const auto& vars3D = (it3D != df_filter_and_weight_to_var3D_list_map.end()) ? it3D->second : empty3D;

    // resolve weight column
    std::string wCol;
    auto itW = weight_specifier_to_column_map.find(weight);
    if (itW != weight_specifier_to_column_map.end()) {
        wCol = itW->second;
    } else {
        wCol = weight.substr(1); // trim off the "_" at the beginning
    }

    std::string suffix = (weight_before_filter)? weight + filter : filter + weight;
    FillHistogramsSingleDataFrame(suffix, df, wCol, vars1D, vars2D, vars3D, hists_not_write, hists_1_2_3D_not_write);
}

void TestClass::FillHistogramsSingleDataFrame(const std::string& suffix, // filter or filter & weight concatenated with custom order
                                            ROOT::RDF::RNode df,
                                            const std::string& weight_col, // weight column, "" if unweighted
                                            const std::vector<std::string>& vars1D,
                                            const std::vector<std::array<std::string, 2>>& vars2D, 
                                            const std::vector<std::array<std::string, 3>>& vars3D,
                                            bool hists_not_write = false,
                                            std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0}) {


    if (vars1D.empty() && vars2D.empty() && vars3D.empty()) {
        std::cerr << "Warning: FillHistogramsSingleDataFrame:: no variables configured "
                  << "for suffix '" << suffix << "'\n";
        return;
    }

    // 1D
    for (const auto& vName : vars1D) {
        var1D* v = nullptr;
        try {
            v = Var1DSearch(vName);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            throw; // as you requested: error if not found
        }

        AxisInfo b = GetAxisInfo(*v, suffix);
        std::string hname = "h_" + v->name + suffix;
        std::string htitle = ";" + v->title + ";";

        if (hists_not_write || hists_1_2_3D_not_write.at(0)){
            hists_to_not_write.push_back(hname);
        }
        
        try {
            if (b.bin_edges) {
                ROOT::RDF::TH1DModel model(hname.c_str(), htitle.c_str(),
                                           b.nbins, b.bin_edges);
                hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var, weight_col);
            } else {
                ROOT::RDF::TH1DModel model(hname.c_str(), htitle.c_str(),
                                           b.nbins, b.min, b.max);
                hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var, weight_col);
            }
        } catch (const std::exception& e) {
            std::cerr << "FillHistogramsSingleDataFrame (1D) exception for suffix '"
                      << suffix << "', column '" << v->var
                      << "': " << e.what() << "\n";
            continue;
        }
    }

    // 2D
    for (const auto& arr : vars2D) {
        const std::string& xName = arr[0];
        const std::string& yName = arr[1];

        var1D *vx = nullptr, *vy = nullptr;
        try {
            vx = Var1DSearch(xName);
            vy = Var1DSearch(yName);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            throw;
        }

        AxisInfo bx = GetAxisInfo(*vx, suffix);
        AxisInfo by = GetAxisInfo(*vy, suffix);

        std::string hname  = "h_" + vy->name + "_vs_" + vx->name + suffix;
        ROOT::RDF::TH2DModel model = MakeTH2DModel(hname, vx->title, vy->title, bx, by);

        if (hists_not_write || hists_1_2_3D_not_write.at(1)){
            hists_to_not_write.push_back(hname);
        }

        try {
            hist2d_rresultptr_map[hname] = df.Histo2D(model, vx->var, vy->var, weight_col);
        } catch (const std::exception& e) {
            std::cerr << "FillHistogramsSingleDataFrame (2D) exception for suffix '"
                      << suffix << "', columns ('" << vx->var << "', '"
                      << vy->var << "'): " << e.what() << "\n";
            continue;
        }
    }

    // 3D
    for (const auto& arr : vars3D) {
        const std::string& xName = arr[0];
        const std::string& yName = arr[1];
        const std::string& zName = arr[2];

        var1D *vx = nullptr, *vy = nullptr, *vz = nullptr;
        try {
            vx = Var1DSearch(xName);
            vy = Var1DSearch(yName);
            vz = Var1DSearch(zName);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            throw;
        }

        AxisInfo bx = GetAxisInfo(*vx, suffix);
        AxisInfo by = GetAxisInfo(*vy, suffix);
        AxisInfo bz = GetAxisInfo(*vz, suffix);

        std::string hname  = "h_" + vz->name + "_vs_" + vy->name + "_vs_" + vx->name + suffix;
        ROOT::RDF::TH3DModel model = MakeTH3DModel(hname, vx->title, vy->title, vz->title,
                                                   bx, by, bz);
        
        if (hists_not_write || hists_1_2_3D_not_write.at(2)){
            hists_to_not_write.push_back(hname);
        }

        try {
            hist3d_rresultptr_map[hname] = df.Histo3D(model, vx->var, vy->var, vz->var, weight_col);
        } catch (const std::exception& e) {
            std::cerr << "FillHistogramsSingleDataFrame (3D) exception for suffix '"
                      << suffix << "', columns ('" << vx->var << "', '"
                      << vy->var << "', '" << vz->var << "'): " << e.what() << "\n";
            continue;
        }
    }
}

var1D* TestClass::Var1DSearch(const std::string& var1DName) const {
    auto it = var1D_dict.find(var1DName);
    if (it == var1D_dict.end() || !(it->second)) {
        std::ostringstream oss;
        oss << "Var1DSearch: variable '" << var1DName << "' not found in var1D_dict";
        throw std::runtime_error(oss.str());
    }
    return it->second;
}

TestClass::AxisInfo TestClass::GetAxisInfo(const var1D& v, const std::string& filter) const {
    bool hasSS = (filter.find("ss") != std::string::npos);
    bool hasOP = (filter.find("op") != std::string::npos);

    if (hasSS && hasOP) {
        std::ostringstream oss;
        oss << "GetAxisInfo: filter '" << filter
            << "' contains both 'ss' and 'op'; ambiguous for variable '" << v.name << "'";
        throw std::runtime_error(oss.str());
    }

    AxisInfo b;

    // default: use global nbins/bins
    int nbins = v.nbins;
    const std::vector<double>* binsVec = nullptr;
    if (!v.bins.empty()) binsVec = &v.bins;
    double vmin = v.min;
    double vmax = v.max;

    // ss/op override
    if (hasSS && v.nbins_ss > 0 && !v.bins_ss.empty()) {
        nbins = v.nbins_ss;
        binsVec = &v.bins_ss;
    } else if (hasOP && v.nbins_op > 0 && !v.bins_op.empty()) {
        nbins = v.nbins_op;
        binsVec = &v.bins_op;
    }

    b.nbins = nbins;

    if (binsVec && !binsVec->empty()) {
        b.bin_edges = binsVec->data();
        b.min  = 0.;
        b.max  = 0.;
    } else {
        b.bin_edges = nullptr;
        b.min  = vmin;
        b.max  = vmax;
    }

    return b;
}

ROOT::RDF::TH2DModel TestClass::MakeTH2DModel(const std::string& hname,
                                            const std::string& xtitle,
                                            const std::string& ytitle,
                                            const AxisInfo& bx,
                                            const AxisInfo& by) const {
    std::string htitle = ";" + xtitle + ";" + ytitle;

    if (bx.bin_edges && by.bin_edges) {
        return ROOT::RDF::TH2DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.bin_edges,
                                    by.nbins, by.bin_edges);
    } else if (bx.bin_edges && !by.bin_edges) {
        return ROOT::RDF::TH2DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.bin_edges,
                                    by.nbins, by.min, by.max);
    } else if (!bx.bin_edges && by.bin_edges) {
        return ROOT::RDF::TH2DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.bin_edges);
    } else {
        return ROOT::RDF::TH2DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.min, by.max);
    }
}

ROOT::RDF::TH3DModel TestClass::MakeTH3DModel(const std::string& hname,
                                            const std::string& xtitle,
                                            const std::string& ytitle,
                                            const std::string& ztitle,
                                            const AxisInfo& bx,
                                            const AxisInfo& by,
                                            const AxisInfo& bz) const {
    std::string htitle = ";" + xtitle + ";" + ytitle + ";" + ztitle;

    const bool ex = bx.bin_edges != nullptr;
    const bool ey = by.bin_edges != nullptr;
    const bool ez = bz.bin_edges != nullptr;

    if (ex && ey && ez) {
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.bin_edges,
                                    by.nbins, by.bin_edges,
                                    bz.nbins, bz.bin_edges);
    } else if (ex && ey && !ez) {
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.bin_edges,
                                    by.nbins, by.bin_edges,
                                    bz.nbins, bz.min, bz.max);
    } else if (ex && !ey && ez) {
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.bin_edges,
                                    by.nbins, by.min, by.max,
                                    bz.nbins, bz.bin_edges);
    } else if (!ex && ey && ez) {
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.bin_edges,
                                    bz.nbins, bz.bin_edges);
    } else if (ex && !ey && !ez) {
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.bin_edges,
                                    by.nbins, by.min, by.max,
                                    bz.nbins, bz.min, bz.max);
    } else if (!ex && ey && !ez) {
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.bin_edges,
                                    bz.nbins, bz.min, bz.max);
    } else if (!ex && !ey && ez) {
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.min, by.max,
                                    bz.nbins, bz.bin_edges);
    } else { // all uniform
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.min, by.max,
                                    bz.nbins, bz.min, bz.max);
    }
}


void TestClass::ReadVar1DJson() {
    // Open file
    std::ifstream in(infile_var1D_json);
    if (!in) {
        throw std::runtime_error(
            "ReadVar1DJson: cannot open JSON file '" + infile_var1D_json + "'");
    }

    // Parse JSON
    json j;
    try {
        in >> j;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            "ReadVar1DJson: failed to parse JSON file '" + infile_var1D_json +
            "': " + e.what());
    }

    // Check that "variables1D" exists and is an array
    if (!j.contains("variables1D") || !j["variables1D"].is_array()) {
        throw std::runtime_error(
            "ReadVar1DJson: mandatory top-level field 'variables1D' "
            "is missing or is not an array");
    }

    // Loop over entries
    for (const auto& jv : j["variables1D"]) {
        // hist_name (mandatory, non-empty string)
        std::string histName;
        if (jv.contains("hist_name") && jv["hist_name"].is_string()) {
            histName = jv["hist_name"].get<std::string>();
        }
        if (histName.empty()) {
            ThrowMissingField("hist_name", /*histName*/ "");
        }

        // Create var1D object
        auto* v = new var1D{};
        v->name = histName;

        // var_name: optional, default to name
        if (jv.contains("var_name") && jv["var_name"].is_string()) {
            std::string tmp = jv["var_name"].get<std::string>();
            if (!tmp.empty()) {
                v->var = tmp;
            }
        }
        if (v->var.empty()) {
            v->var = v->name;
        }

        // title: optional, default to name
        if (jv.contains("title") && jv["title"].is_string()) {
            std::string tmp = jv["title"].get<std::string>();
            if (!tmp.empty()) {
                v->title = tmp;
            }
        }
        if (v->title.empty()) {
            v->title = v->name;
        }

        // binning: mandatory
        if (!jv.contains("binning") || jv["binning"].is_null()) {
            delete v;
            ThrowMissingField("binning", histName);
        }

        const auto& jb = jv["binning"];

        if (jb.is_string()) {
            // binning by name: look up in hist_binning_map
            std::string binName = jb.get<std::string>();
            if (binName.empty()) {
                delete v;
                ThrowMissingField("binning", histName);
            }

            auto it = hist_binning_map.find(binName);
            auto it_ss = hist_binning_map.find(binName + "_ss");
            auto it_op = hist_binning_map.find(binName + "_op");
            if (it == hist_binning_map.end() && (it_ss == hist_binning_map.end() || it_op == hist_binning_map.end())) {
                delete v;
                std::ostringstream oss;
                oss << "ReadVar1DJson: binning name '" << binName
                    << "' (referenced by hist_name='" << histName
                    << "') not found in hist_binning_map";
                throw std::runtime_error(oss.str());
            }

            if (it_ss != hist_binning_map.end() && it_op != hist_binning_map.end()){
                v->bins_ss = it_ss->second;
                v->nbins_ss = static_cast<int>(v->bins_ss.size()) - 1;
                v->bins_op = it_op->second;
                v->nbins_op = static_cast<int>(v->bins_op.size()) - 1;
            } else{
                v->bins = it->second;
                v->nbins = static_cast<int>(v->bins.size()) - 1;                
            }

        } else if (jb.is_object()) {
            // binning as {nbins, min, max}
            if (!jb.contains("nbins") || !jb.contains("min") || !jb.contains("max")) {
                delete v;
                ThrowMissingField("binning.nbins/min/max", histName);
            }

            if (!jb["nbins"].is_number_integer() ||
                !jb["min"].is_number() ||
                !jb["max"].is_number()) {
                delete v;
                std::ostringstream oss;
                oss << "ReadVar1DJson: binning object for hist_name='"
                    << histName
                    << "' must have integer 'nbins' and numeric 'min', 'max'";
                throw std::runtime_error(oss.str());
            }

            v->nbins = jb["nbins"].get<int>();
            v->vmin  = jb["min"].get<double>();
            v->vmax  = jb["max"].get<double>();

        } else {
            delete v;
            std::ostringstream oss;
            oss << "ReadVar1DJson: 'binning' for hist_name='" << histName
                << "' must be either a string or an object {nbins,min,max}";
            throw std::runtime_error(oss.str());
        }

        // Final validation using your isValid()
        if (!v->isValid()) {
            std::ostringstream oss;
            oss << "ReadVar1DJson: invalid binning specification for hist_name='"
                << histName << "'";
            delete v;
            throw std::runtime_error(oss.str());
        }

        // All good, store pointer in dictionary
        var1D_dict[v->name] = v;
    }
}

void TestClass::ThrowMissingField(const std::string& field,
                                const std::string& histName) {
    std::ostringstream oss;
    oss << "ReadVar1DJson: mandatory field '" << field
        << "' is missing or empty";
    if (!histName.empty()) {
        oss << " (hist_name='" << histName << "')";
    }
    throw std::runtime_error(oss.str());
}

void TestClass::PrintVar1DList() const {
    std::cout << "===== var1D_dict contents (" 
              << var1D_dict.size() << " entries) =====\n";

    for (auto pair : var1D_dict){
        const var1D* v = pair.second;
        if (!v) continue;

        std::cout << "  name  = " << v->name  << "\n"
                  << "  var   = " << v->var   << "\n"
                  << "  title = " << v->title << "\n";

        if (!v->bins_ss.empty()) {
            std::cout << "Same-sign variable-binning (" << v->nbins_ss
                      << " bins): [";
            for (size_t j = 0; j < v->bins_ss.size(); ++j) {
                std::cout << v->bins_ss[j];
                if (j + 1 < v->bins_ss.size()) std::cout << ", ";
            }
            std::cout << "]\n";

            std::cout << "Opposite-sign variable-binning (" << v->nbins_op
                      << " bins): [";
            for (size_t j = 0; j < v->bins_op.size(); ++j) {
                std::cout << v->bins_op[j];
                if (j + 1 < v->bins_op.size()) std::cout << ", ";
            }
            std::cout << "]\n";
        } else if (!v->bins.empty()) {
            std::cout << "  variable-binning (" << v->nbins
                      << " bins): [";
            for (size_t j = 0; j < v->bins.size(); ++j) {
                std::cout << v->bins[j];
                if (j + 1 < v->bins.size()) std::cout << ", ";
            }
            std::cout << "]\n";
        } else {
            std::cout << "  nbins = " << v->nbins << "\n"
                      << "  vmin  = " << v->vmin  << "\n"
                      << "  vmax  = " << v->vmax  << "\n";
        }

        std::cout << "  isValid = " << (v->isValid() ? "true" : "false") << "\n";
        std::cout << "--------------------------------------------\n";
    }
}
