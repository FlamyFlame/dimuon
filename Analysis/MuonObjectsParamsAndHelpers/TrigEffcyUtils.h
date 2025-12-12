namespace TrigEffcyUtils {

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
}

void TrigEffcyUtils::flatten_levels(
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

void TrigEffcyUtils::write_post_sum_levels(
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
void TrigEffcyUtils::SumTrigEffHistsGeneric(
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
