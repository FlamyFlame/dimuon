namespace HistFillUtils {

    TGraphAsymmErrors* proj_divide_and_write (TH2* hNum2D, TH2* hDen2D, std::map<std::string, TGraphAsymmErrors*>* graph_map = nullptr);
    TGraphAsymmErrors* proj_divide_and_write (TH2* hNum2D, TH2* hDen2D, bool projy, int firstbin, int lastbin, std::string proj_range_str, std::map<std::string, TGraphAsymmErrors*>* graph_map = nullptr);
    TH1D* proj_and_write (TH2* h2D, bool projy, int firstbin = 1, int lastbin = -1, std::string proj_range_str = "", std::map<std::string, TH1D*>* h1D_map = nullptr);

    TGraphAsymmErrors* divide_and_write (TH1* hNum, TH1* hDen, std::map<std::string, TGraphAsymmErrors*>* graph_map = nullptr);

    void flatten_levels(const std::vector<std::vector<std::string>>& levels,
                        std::vector<std::string>& flattened);

    void write_post_sum_levels(const std::vector<std::vector<std::string>>& levels_pre_sum,
                               const std::vector<int>& levels_to_be_summed,
                               std::vector<std::vector<std::string>>& levels_post_sum);

    template <class TH, class Var, class MakeBaseName>
    void SumTrigEffHistsGeneric(
        const std::vector<Var>& vars,
        const std::vector<std::string>& filters_post_sum,
        const std::vector<std::string>& filters_to_be_summed,
        std::map<std::string, TH*>& hist_map,
        MakeBaseName makeBaseName);

    template <typename H>
    void write_hist_map_vector(
        std::map<std::string, H*>& m,
        const std::vector<std::string>& hists_to_not_write);
}


// helper for projection, making & writing of TEfficiency graphs
TGraphAsymmErrors* HistFillUtils::proj_divide_and_write (TH2* hNum2D, TH2* hDen2D, std::map<std::string, TGraphAsymmErrors*>* graph_map = nullptr) {
    return proj_divide_and_write(hNum2D, hDen2D, true, 1, -1, "", graph_map);
}

TGraphAsymmErrors* HistFillUtils::proj_divide_and_write (TH2* hNum2D, TH2* hDen2D, bool projy, int firstbin, int lastbin, std::string proj_range_str, std::map<std::string, TGraphAsymmErrors*>* graph_map = nullptr) {
    if (!hNum2D || !hDen2D) return (TGraphAsymmErrors*)nullptr;

    // suffix that captures projection axis & range
    std::string proj_suffix = projy? "_py" : "_px";
    if (proj_range_str != "") proj_suffix += "_" + proj_range_str;

    std::unique_ptr<TH1> hNum1D(
        projy ? hNum2D->ProjectionY(Form("%s%s", hNum2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
              : hNum2D->ProjectionX(Form("%s%s", hNum2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
    );

    std::unique_ptr<TH1> hDen1D(
        projy ? hDen2D->ProjectionY(Form("%s%s", hDen2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
              : hDen2D->ProjectionX(Form("%s%s", hDen2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
    );

    auto g = new TGraphAsymmErrors();
    g->BayesDivide(hNum1D.get(), hDen1D.get());

    // Name: "h_<...>" → "g_<...>" from the *numerator* histo name
    std::string n = hNum1D->GetName() ? hNum1D->GetName() : "graph";
    
    n += "_divided";
    if (n.rfind("h_", 0) == 0) n.replace(0, 2, "g_"); else n = "g_" + n;
    g->SetName(n.c_str());

    if (graph_map == nullptr)   g->Write(g->GetName(), TObject::kOverwrite);
    else{
        (*graph_map)[g->GetName()] = g;
    }

    return g;
}

TH1D* HistFillUtils::proj_and_write (TH2* h2D, bool projy, int firstbin = 1, int lastbin = -1, std::string proj_range_str = "", std::map<std::string, TH1D*>* h1D_map = nullptr){
    if (!h2D) return (TH1D*)nullptr;

    // suffix that captures projection axis & range
    std::string proj_suffix = projy? "_py" : "_px";
    if (proj_range_str != "") proj_suffix += "_" + proj_range_str;

    TH1D* h1D(
        projy ? h2D->ProjectionY(Form("%s%s", h2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
              : h2D->ProjectionX(Form("%s%s", h2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
    );
    if (!h1D) return nullptr;
    
    if (h1D_map == nullptr)     h1D->Write(h1D->GetName(), TObject::kOverwrite);
    else                        (*h1D_map)[h1D->GetName()] = h1D;

    return h1D;
}

TGraphAsymmErrors* HistFillUtils::divide_and_write (TH1* hNum, TH1* hDen, std::map<std::string, TGraphAsymmErrors*>* graph_map = nullptr) {
    if (!hNum || !hDen) return (TGraphAsymmErrors*)nullptr;

    auto g = new TGraphAsymmErrors();
    g->BayesDivide(hNum, hDen);

    // Name: "h_<...>" → "g_<...>" from the *numerator* histo name
    std::string n = hNum->GetName() ? hNum->GetName() : "graph";
    
    n += "_divided";
    if (n.rfind("h_", 0) == 0) n.replace(0, 2, "g_"); else n = "g_" + n;
    g->SetName(n.c_str());

    if (graph_map == nullptr)   g->Write(g->GetName(), TObject::kOverwrite);
    else{
        (*graph_map)[g->GetName()] = g;
    }

    return g;
}

void HistFillUtils::flatten_levels(
    const std::vector<std::vector<std::string>>& levels,
    std::vector<std::string>& flattened)
{
    std::vector<int> level_sizes;
    int total_size = 1;

    for (const auto& vlevel : levels) {
        level_sizes.push_back(static_cast<int>(vlevel.size()));
        if (!vlevel.empty()) {
            total_size *= static_cast<int>(vlevel.size());
        }
    }

    flattened.clear();
    flattened.assign(total_size, "");

    for (int flattened_ind = 0; flattened_ind < total_size; ++flattened_ind) {
        int remaining_ind = flattened_ind;

        for (int level_ind = 0; level_ind < static_cast<int>(levels.size()); ++level_ind) {
            if (level_sizes.at(level_ind) == 0) continue;

            int dim_remaining_levels = 1;
            for (int level_ind_remain = level_ind + 1;
                 level_ind_remain < static_cast<int>(levels.size());
                 ++level_ind_remain)
            {
                if (level_sizes.at(level_ind_remain) != 0) {
                    dim_remaining_levels *= level_sizes.at(level_ind_remain);
                }
            }

            int cur_level_index = remaining_ind / dim_remaining_levels;
            remaining_ind       = remaining_ind % dim_remaining_levels;

            flattened.at(flattened_ind) += levels.at(level_ind).at(cur_level_index);
        }
    }
}

void HistFillUtils::write_post_sum_levels(
    const std::vector<std::vector<std::string>>& levels_pre_sum,
    const std::vector<int>& levels_to_be_summed,
    std::vector<std::vector<std::string>>& levels_post_sum)
{
    levels_post_sum.clear();
    for (int level_ind = 0; level_ind < static_cast<int>(levels_pre_sum.size()); ++level_ind) {
        // level is NOT contracted → copy it into post_sum
        if (std::find(levels_to_be_summed.begin(),
                      levels_to_be_summed.end(),
                      level_ind) == levels_to_be_summed.end())
        {
            levels_post_sum.push_back(levels_pre_sum.at(level_ind));
        }
    }
}

//--------- TEMPLATE FUNCTION FOR SUMMING TRIGGER EFFICIENCY HISTOGRAMS WITH AN ARBITRARY DATATYPE (DIMENSION) ---------
template <class TH, class Var, class MakeBaseName>
void HistFillUtils::SumTrigEffHistsGeneric(
    const std::vector<Var>& vars,
    const std::vector<std::string>& filters_post_sum,
    const std::vector<std::string>& filters_to_be_summed,
    std::map<std::string, TH*>& hist_map,
    MakeBaseName makeBaseName)        // only depends on Var
{
    for (const auto& var : vars) {
        const std::string base = makeBaseName(var);

        for (const auto& post_sum_filter : filters_post_sum) {

            std::string hname_post_sum = base + post_sum_filter;

            TH*  h_post_sum     = nullptr;
            bool first_instance = true;

            for (const auto& to_sum_filter : filters_to_be_summed) {
                std::string hname_pre_sum = base + to_sum_filter + post_sum_filter;

                auto it = hist_map.find(hname_pre_sum);
                if (it == hist_map.end()) {
                    std::cerr << "SumTrigEffHistsGeneric: histogram "
                              << hname_pre_sum << " not found!\n";
                    continue;
                }

                TH* h_pre_sum = it->second;

                if (first_instance) {
                    h_post_sum = static_cast<TH*>(
                        h_pre_sum->Clone(hname_post_sum.c_str())
                    );
                    first_instance = false;
                } else {
                    h_post_sum->Add(h_pre_sum);
                }
            }

            if (h_post_sum != nullptr) {
                hist_map.emplace(hname_post_sum, h_post_sum);
            }
        }
    }
}

template <typename H>
void HistFillUtils::write_hist_map_vector(
    std::map<std::string, H*>& m,
    const std::vector<std::string>& hists_to_not_write)
{
    for (auto& kv : m) {
        const auto& name = kv.first;
        if (std::find(hists_to_not_write.begin(), hists_to_not_write.end(), name)
            != hists_to_not_write.end()) continue;

        kv.second->Write();        // THnD in a vector
    }
}
