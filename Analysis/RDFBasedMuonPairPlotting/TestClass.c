#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>
#include <string>

#include <nlohmann/json.hpp>
using nlohmann::json;

// assuming var1D is defined as you wrote
// and that this is inside some class that has:
//   std::map<std::string, std::vector<double>> hist_binning_map;
//   std::string infile_var1D_json;
//   std::vector<var1D*> var1D_dict;

struct var1D {
    std::string var{""}; // e.g, dR / minv / m1.pt
    std::string name{""}; // e.g, dR_zoomin / minv_log / pt_lead
    std::string title{""}; // e.g, #Delta R / m_{#mu#mu}
    int         nbins{0};
    int         nbins_ss{0};
    int         nbins_op{0};
    double      vmin{0.};
    double      vmax{0.};
    std::vector<double> bins{};          // optional variable binning
    std::vector<double> bins_ss{};          // optional variable binning
    std::vector<double> bins_op{};          // optional variable binning

    bool isValid() const {
        if (nbins <= 0) return false;
        if (!bins.empty()) return (static_cast<int>(bins.size()) == nbins + 1);
        return vmin < vmax;
    }
};


class TestClass {
private:
    // members
    std::map<std::string, std::vector<double>> hist_binning_map;
    std::string infile_var1D_json = "var1D.json";
    std::vector<var1D*> var1D_dict;

    void ReadVar1DJson();

private:
    static void BuildHistBinningMap();
    static void throwMissingField(const std::string& field,
                                  const std::string& histName);
public:
    void Run();
};

void TestClass::Run(){
    BuildHistBinningMap();
    ReadVar1DJson();
    PrintVar1DList();
}


void TestClass::BuildHistBinningMap(){
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

    hist1d_rresultptrs["minv_log"] = minv_bins_log[0];
}


void TestClass::throwMissingField(const std::string& field,
                                const std::string& histName) {
    std::ostringstream oss;
    oss << "ReadVar1DJson: mandatory field '" << field
        << "' is missing or empty";
    if (!histName.empty()) {
        oss << " (hist_name='" << histName << "')";
    }
    throw std::runtime_error(oss.str());
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
            throwMissingField("hist_name", /*histName*/ "");
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
            throwMissingField("binning", histName);
        }

        const auto& jb = jv["binning"];

        if (jb.is_string()) {
            // binning by name: look up in hist_binning_map
            std::string binName = jb.get<std::string>();
            if (binName.empty()) {
                delete v;
                throwMissingField("binning", histName);
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
                throwMissingField("binning.nbins/min/max", histName);
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
