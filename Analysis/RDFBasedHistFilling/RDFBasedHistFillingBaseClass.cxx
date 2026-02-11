#include "RDFBasedHistFillingBaseClass.h"
#include <ROOT/RDFHelpers.hxx>     // ROOT::RDF::RunGraphs
#include <ROOT/RResultHandle.hxx>  // ROOT::RDF::RResultHandle
#include <map>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <algorithm>

void RDFBasedHistFillingBaseClass::Run(){

	std::cout << "start the run" << std::endl;
	auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);

	Initialize();

    ProcessData();

    // info summary for the current for loop (input file)
    std::cout << "#event loop runs per data frame is: " << df_map.at("df_ss").GetNRuns() << std::endl;
    std::cout << "#slots per data frame is: " << df_map.at("df_ss").GetNSlots() << std::endl;
    std::cout << "# 1D histograms: " << hist1D_map.size() << std::endl;
    std::cout << "# 2D histograms: " << hist2D_map.size() << std::endl;
    std::cout << "# 3D histograms: " << hist2D_map.size() << std::endl;

	Finalize();
}

// ------- PROCESS DATA -------
void RDFBasedHistFillingBaseClass::ProcessData(){    
    ROOT::EnableImplicitMT();
    
    CreateRDFs();
    FillHistograms();

    HistPostProcess();
}

// ------- Create RDataFrame and keep ownership -------
void RDFBasedHistFillingBaseClass::CreateRDFs(){
    rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree_ss, input_files));
    ROOT::RDF::RNode node_ss = *(rdf_store.back());
    df_map.emplace("df_ss", node_ss);

    rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree_op, input_files));
    ROOT::RDF::RNode node_op = *(rdf_store.back());
    df_map.emplace("df_op", node_op);
}


void RDFBasedHistFillingBaseClass::Initialize(){
	std::cout << "Calling Initialize." << std::endl;
    InitializeBaseCommon();
    InitAnalysisSettingsHook(); // hook to set datatype-specific analysis settings
    SetIOPathsHook();
    InitOutput();
	BuildFilterToVarListMap();
	BuildHistBinningMap();
	ReadVar1DJson();
}

void RDFBasedHistFillingBaseClass::InitializeBaseCommon(){
    hist1d_rresultptr_map.clear();
    hist2d_rresultptr_map.clear();
    hist3d_rresultptr_map.clear();
    
    TH1::SetDefaultSumw2(kTRUE); // turn on Sumw2 for all histograms
}

void RDFBasedHistFillingBaseClass::InitOutput(){
    std::string file_mode = "recreate";
    m_outfile = new TFile(output_file.c_str(), file_mode.c_str());
}

// ---------- ----------
void RDFBasedHistFillingBaseClass::HistPostProcess(){
    HistPostProcessBaseCommon();
    HistPostProcessExtra();
}

void RDFBasedHistFillingBaseClass::HistPostProcessBaseCommon(){
    if (debug_mode) std::cout << "Calling HistPostProcess" << std::endl;

    // 1) Force completion of all booked actions (1D/2D/3D) before GetPtr()
    std::vector<ROOT::RDF::RResultHandle> handles;
    handles.reserve(hist1d_rresultptr_map.size()
                  + hist2d_rresultptr_map.size()
                  + hist3d_rresultptr_map.size());

    for (auto &kv : hist1d_rresultptr_map) handles.emplace_back(ROOT::RDF::RResultHandle{kv.second});
    for (auto &kv : hist2d_rresultptr_map) handles.emplace_back(ROOT::RDF::RResultHandle{kv.second});
    for (auto &kv : hist3d_rresultptr_map) handles.emplace_back(ROOT::RDF::RResultHandle{kv.second});

    ROOT::RDF::RunGraphs(std::move(handles));

    // 2) Convert to raw histogram pointers (RResultPtr must outlive your use of these!)
    hist1D_map.clear();
    hist2D_map.clear();
    hist3D_map.clear();

    for (auto &kv : hist1d_rresultptr_map) hist1D_map.emplace(kv.first, kv.second.GetPtr());
    for (auto &kv : hist2d_rresultptr_map) hist2D_map.emplace(kv.first, kv.second.GetPtr());
    for (auto &kv : hist3d_rresultptr_map) hist3D_map.emplace(kv.first, kv.second.GetPtr());
}

//--------- ---------
void RDFBasedHistFillingBaseClass::BuildFilterToVarListMap(){
    BuildFilterToVarListMapBaseCommon();
    BuildFilterToVarListMapExtra();
}

void RDFBasedHistFillingBaseClass::BuildFilterToVarListMapBaseCommon(){
    // Add common filter to var maps here
    df_filter_and_weight_to_var1D_list_map[{"_single_b", "_crossx"}] = {"pair_pt", "pair_eta"}; // crossx as a weighting
}

// ---------- ----------
void RDFBasedHistFillingBaseClass::BuildHistBinningMap(){
    BuildHistBinningMapBaseCommon();
    BuildHistBinningMapExtra();
}

void RDFBasedHistFillingBaseClass::BuildHistBinningMapBaseCommon(){

    // ------- pT binnings -------
    hist_binning_map["pT_bins_80"] = pms.pT_bins_80;

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

void RDFBasedHistFillingBaseClass::Finalize(){
    WriteOutput();
    Cleanup();
}

void RDFBasedHistFillingBaseClass::WriteOutput(){
    WriteOutputBaseCommon();
    WriteOutputExtra();
}

void RDFBasedHistFillingBaseClass::WriteOutputBaseCommon(){
    // ----- write output -----

    // std::vector<THnD> maps
    TrigEffcyUtils::write_hist_map_vector(hist1D_map, hists_to_not_write);
    TrigEffcyUtils::write_hist_map_vector(hist2D_map, hists_to_not_write);
    TrigEffcyUtils::write_hist_map_vector(hist3D_map, hists_to_not_write);
}

void RDFBasedHistFillingBaseClass::Cleanup(){
    CleanupBaseCommon();
    CleanupExtra();
}

void RDFBasedHistFillingBaseClass::CleanupBaseCommon(){
    // ----- delete pointers -----
    m_outfile->Close();
    delete m_outfile;

    for (auto kv : var1D_dict){
        delete kv.second;
    }
}


// ---------- THROW MISSING FIELD IN JSON READING HELPER ----------
void RDFBasedHistFillingBaseClass::ThrowMissingField(const std::string& field,
                                const std::string& histName) {
    std::ostringstream oss;
    oss << "ReadVar1DJson: mandatory field '" << field
        << "' is missing or empty";
    if (!histName.empty()) {
        oss << " (hist_name='" << histName << "')";
    }
    throw std::runtime_error(oss.str());
}

// ---------- VAR1D JSON FILE READING ----------
void RDFBasedHistFillingBaseClass::ReadVar1DJson() {
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

// ---------- VAR1D LIST PRINTINT FOR DEBUG ----------
void RDFBasedHistFillingBaseClass::PrintVar1DList() const {
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

// ---------- HELPER FUNCTIONS ----------

void RDFBasedHistFillingBaseClass::FillHistogramsSingleDataFrame(const std::string& filter,
                                            ROOT::RDF::RNode df,
                                            bool hists_write = true,
                                            std::array<bool, 3> hists_1_2_3D_write = {1,1,1}) {
    static const std::vector<std::string> empty1D;
    static const std::vector<std::array<std::string, 2>> empty2D;
    static const std::vector<std::array<std::string, 3>> empty3D;

    auto it1D = df_filter_to_var1D_list_map.find(filter);
    auto it2D = df_filter_to_var2D_list_map.find(filter);
    auto it3D = df_filter_to_var3D_list_map.find(filter);

    const auto& vars1D = (it1D != df_filter_to_var1D_list_map.end()) ? it1D->second : empty1D;
    const auto& vars2D = (it2D != df_filter_to_var2D_list_map.end()) ? it2D->second : empty2D;
    const auto& vars3D = (it3D != df_filter_to_var3D_list_map.end()) ? it3D->second : empty3D;

    FillHistogramsSingleDataFrame(filter, df, "", vars1D, vars2D, vars3D, hists_write, hists_1_2_3D_write);
}

void RDFBasedHistFillingBaseClass::FillHistogramsSingleDataFrame(const std::string& filter,
                                            const std::string& weight,
                                            ROOT::RDF::RNode df,
                                            bool weight_before_filter = false,
                                            bool hists_write = true,
                                            std::array<bool, 3> hists_1_2_3D_write = {1,1,1}) {
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
    FillHistogramsSingleDataFrame(suffix, df, wCol, vars1D, vars2D, vars3D, hists_write, hists_1_2_3D_write);
}

void RDFBasedHistFillingBaseClass::FillHistogramsSingleDataFrame(const std::string& suffix, // filter or filter & weight concatenated with custom order
                                            ROOT::RDF::RNode df,
                                            const std::string& weight_col, // weight column, "" if unweighted
                                            const std::vector<std::string>& vars1D,
                                            const std::vector<std::array<std::string, 2>>& vars2D, 
                                            const std::vector<std::array<std::string, 3>>& vars3D,
                                            bool hists_write = true,
                                            std::array<bool, 3> hists_1_2_3D_write = {1,1,1}) {


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

        if (!hists_write || !hists_1_2_3D_write.at(0)){
            hists_to_not_write.push_back(hname);
        }
        
        try {
            if (b.bin_edges) {
                ROOT::RDF::TH1DModel model(hname.c_str(), htitle.c_str(),
                                           b.nbins, b.bin_edges);
                if (weight_col.empty()) hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var);
                else                    hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var, weight_col);
            } else {
                ROOT::RDF::TH1DModel model(hname.c_str(), htitle.c_str(),
                                           b.nbins, b.min, b.max);
                if (weight_col.empty()) hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var);
                else                    hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var, weight_col);
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

        if (!hists_write || !hists_1_2_3D_write.at(1)){
            hists_to_not_write.push_back(hname);
        }

        try {
            if (weight_col.empty()) hist2d_rresultptr_map[hname] = df.Histo2D(model, vx->var, vy->var);
            else                    hist2d_rresultptr_map[hname] = df.Histo2D(model, vx->var, vy->var, weight_col);
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
        
        if (!hists_write || !hists_1_2_3D_write.at(2)){
            hists_to_not_write.push_back(hname);
        }

        try {
            if (weight_col.empty()) hist3d_rresultptr_map[hname] = df.Histo3D(model, vx->var, vy->var, vz->var);
            else                    hist3d_rresultptr_map[hname] = df.Histo3D(model, vx->var, vy->var, vz->var, weight_col);

        } catch (const std::exception& e) {
            std::cerr << "FillHistogramsSingleDataFrame (3D) exception for suffix '"
                      << suffix << "', columns ('" << vx->var << "', '"
                      << vy->var << "', '" << vz->var << "'): " << e.what() << "\n";
            continue;
        }
    }
}

var1D* RDFBasedHistFillingBaseClass::Var1DSearch(const std::string& var1DName) const {
    auto it = var1D_dict.find(var1DName);
    if (it == var1D_dict.end() || !(it->second)) {
        std::ostringstream oss;
        oss << "Var1DSearch: variable '" << var1DName << "' not found in var1D_dict";
        throw std::runtime_error(oss.str());
    }
    return it->second;
}

AxisInfo RDFBasedHistFillingBaseClass::GetAxisInfo(const var1D& v, const std::string& filter) const {
    bool hasSS = (filter.find("_ss") != std::string::npos);
    bool hasOP = (filter.find("_op") != std::string::npos);

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
    double vmin = v.vmin;
    double vmax = v.vmax;

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

ROOT::RDF::TH2DModel RDFBasedHistFillingBaseClass::MakeTH2DModel(const std::string& hname,
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

ROOT::RDF::TH3DModel RDFBasedHistFillingBaseClass::MakeTH3DModel(const std::string& hname,
                                              const std::string& xtitle,
                                              const std::string& ytitle,
                                              const std::string& ztitle,
                                              const AxisInfo& bx,
                                              const AxisInfo& by,
                                              const AxisInfo& bz) const
{
    std::string htitle = ";" + xtitle + ";" + ytitle + ";" + ztitle;

    const bool ex = (bx.bin_edges != nullptr);
    const bool ey = (by.bin_edges != nullptr);
    const bool ez = (bz.bin_edges != nullptr);

    // 1) NONE has bin_edges → use uniform on all
    if (!ex && !ey && !ez) {
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.min, by.max,
                                    bz.nbins, bz.min, bz.max);
    }

    // 2) At least one has bin_edges → all-variable constructor must be used
    // Define storage for any axis that lacks explicit bin_edges
    std::vector<double> xStore, yStore, zStore;

    auto makeEdges = [](const AxisInfo& a, std::vector<double>& store) -> const double* {
        if (a.bin_edges) return a.bin_edges;
        store.resize(a.nbins + 1);
        double step = (a.max - a.min) / a.nbins;
        for (int i = 0; i <= a.nbins; ++i)
            store[i] = a.min + i * step;
        return store.data();
    };

    const double* xEdges = makeEdges(bx, xStore);
    const double* yEdges = makeEdges(by, yStore);
    const double* zEdges = makeEdges(bz, zStore);

    return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                bx.nbins, xEdges,
                                by.nbins, yEdges,
                                bz.nbins, zEdges);
}

