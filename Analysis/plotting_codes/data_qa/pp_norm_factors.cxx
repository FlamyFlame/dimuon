// Abnormal-run investigation (docs/tracking/data_per_run_signal_yield_qa.md): pp24 per-run N_sig/L factored as
//   N_sig/L = (N_HLT_2mu4/L) x (N_sig/N_HLT_2mu4)
//   N_HLT_2mu4 = events with the event-level HLT_2mu4_L12MU3V decision in the raw skim NTUP (hlt_events_per_run.cxx)
//   N_sig      = signal-region OS Tight distinct events (per_run_signal_yield_table.cxx, <tag>_per_run.csv)
//   L          = Prescale-Corrected HLT_2mu4_L12MU3V lumi [pb^-1] (unprescaled; the table the N/L plot uses)
// pp has no event selection.  User choice 2026-10-12: the analysis trigger 2mu4 (unprescaled), not the prescaled mu4.
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "TStyle.h"
#include "TSystem.h"
#include "plotting_codes/data_qa/per_run_factor_plot.h"

namespace ppnf {
inline std::map<int, std::vector<double>> ReadCsv(const std::string& path) {   // run -> remaining columns
    std::map<int, std::vector<double>> m; std::ifstream f(path); std::string line; std::getline(f, line);
    while (std::getline(f, line)) {
        std::stringstream ss(line); std::string c; std::vector<std::string> v;
        while (std::getline(ss, c, ',')) v.push_back(c);
        if (v.empty() || v[0].find_first_not_of("0123456789 ") != std::string::npos) continue;   // skip Total row
        std::vector<double> d; for (size_t i = 1; i < v.size(); ++i) d.push_back(std::stod(v[i]));
        m[std::stoi(v[0])] = d;
    }
    if (m.empty()) throw std::runtime_error("empty csv " + path);
    return m;
}
}

void pp_norm_factors(const std::string& hlt_csv, const std::string& sig_csv, const std::string& lumi_csv, const std::string& outdir) {
    using namespace ppnf;
    gStyle->SetOptStat(0);
    gSystem->mkdir(outdir.c_str(), true);
    const auto hlt = ReadCsv(hlt_csv), sig = ReadCsv(sig_csv), lumi = ReadCsv(lumi_csv);
    std::vector<int> runs; std::vector<double> nh, ns, L;
    double SH = 0, SS = 0, SL = 0;
    std::ofstream o(outdir + "/pp24_run_norm_factors.csv");
    o << "run,L_pb,n_hlt_2mu4,n_sig,hlt_per_L,sig_per_hlt\n";
    for (auto& kv : lumi) {
        const int r = kv.first; const double l = kv.second.at(5);     // Prescale Corrected
        if (l <= 0) continue;
        if (!hlt.count(r) || !sig.count(r)) throw std::runtime_error("run " + std::to_string(r) + " missing in hlt or signal csv");
        const double h = hlt.at(r).at(0), s = sig.at(r).at(1);        // n_hlt ; n_events
        runs.push_back(r); nh.push_back(h); ns.push_back(s); L.push_back(l);
        SH += h; SS += s; SL += l;
        o << r << "," << l << "," << h << "," << s << "," << h / l << "," << s / h << "\n";
    }
    for (auto& kv : hlt) if (!lumi.count(kv.first)) throw std::runtime_error("run " + std::to_string(kv.first) + " not in lumi table");
    const std::string lab = "pp 2024, HLT_2mu4_L12MU3V, #sqrt{s} = 5.36 TeV";
    DrawPerRunFactor(runs, nh, L, SH / SL, "N_{events}^{HLT_2mu4} / L", " pb", lab, outdir + "/pp24_hlt_events_per_lumi_per_run.png");
    DrawPerRunFactor(runs, ns, nh, SS / SH, "N_{events}^{signal} / N_{events}^{HLT_2mu4}", "", lab, outdir + "/pp24_signal_per_hlt_event_per_run.png");
}
