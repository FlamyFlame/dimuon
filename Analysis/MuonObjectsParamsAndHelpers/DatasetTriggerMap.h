#pragma once
#include <map>
#include <string>
#include <stdexcept>
#include <utility>

// DatasetTriggerMap: fixed mapping from (run_year, data_type) to trigger.
//
// Two separate maps:
//   FileSuffix: the trigger string as it appears in ROOT file names
//               (e.g. PbPb 2024 uses "single_mu4" in filenames)
//   Label:      the human-readable trigger string for plot legends
//               (e.g. PbPb 2024 displays as "mu4")
//
// data_type strings: "PbPb", "pp", "pp_2mu4"

class DatasetTriggerMap {
public:
    // Map to file-name suffix (matches RDFBasedHistFillingData trig_suffix, no leading _)
    static const std::map<std::pair<int,std::string>, std::string>& FileSuffixMap() {
        static const std::map<std::pair<int,std::string>, std::string> m = {
            {{23, "PbPb"},    "single_mu4" },
            {{24, "PbPb"},    "single_mu4" },
            {{25, "PbPb"},    "single_mu4" },
            {{24, "pp"},      "2mu4"       },
            {{24, "pp_2mu4"}, "2mu4"       },
        };
        return m;
    }

    // Map to display label for legends (shorter, physicist-friendly)
    static const std::map<std::pair<int,std::string>, std::string>& LabelMap() {
        static const std::map<std::pair<int,std::string>, std::string> m = {
            {{23, "PbPb"},    "mu4"        },
            {{24, "PbPb"},    "mu4"        },
            {{25, "PbPb"},    "mu4"        },
            {{24, "pp"},      "2mu4"       },
            {{24, "pp_2mu4"}, "2mu4"       },
        };
        return m;
    }

    // Returns file-name trigger suffix. Throws if not found.
    // run_year: 2-digit or 4-digit accepted.
    static const std::string& GetTrigger(int run_year, const std::string& data_type) {
        run_year %= 2000;
        const auto& m = FileSuffixMap();
        auto it = m.find({run_year, data_type});
        if (it == m.end()) {
            throw std::runtime_error(
                "DatasetTriggerMap::GetTrigger: no entry for (run_year=" +
                std::to_string(run_year) + ", data_type=\"" + data_type + "\"). "
                "Valid keys: (23,PbPb),(24,PbPb),(25,PbPb),(24,pp),(24,pp_2mu4)."
            );
        }
        return it->second;
    }

    // Returns display label for legends. Throws if not found.
    static const std::string& GetTriggerLabel(int run_year, const std::string& data_type) {
        run_year %= 2000;
        const auto& m = LabelMap();
        auto it = m.find({run_year, data_type});
        if (it == m.end()) {
            throw std::runtime_error(
                "DatasetTriggerMap::GetTriggerLabel: no entry for (run_year=" +
                std::to_string(run_year) + ", data_type=\"" + data_type + "\"). "
                "Valid keys: (23,PbPb),(24,PbPb),(25,PbPb),(24,pp),(24,pp_2mu4)."
            );
        }
        return it->second;
    }

    // Convenience: return "DataType YYYY, trigger_label" string for plot legends.
    static std::string GetLabel(int run_year, const std::string& data_type) {
        run_year %= 2000;
        const std::string& trig = GetTriggerLabel(run_year, data_type);
        std::string dtype_display = (data_type.rfind("pp", 0) == 0) ? "pp" : data_type;
        return dtype_display + " 20" + std::to_string(run_year) + ", " + trig;
    }
};
