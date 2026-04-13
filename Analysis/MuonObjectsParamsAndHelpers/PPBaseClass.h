#pragma once
#include <map>
#include <string>
#include <stdexcept>
#include <utility>

// PPBaseClass: single source of truth for pp integrated-luminosity normalization factors.
// Maps (run_year % 2000, trigger_string) --> 1/L_int [pb^{-1}], used to normalize
// event counts to differential cross-sections.
//
// Trigger strings match the trig_suffix convention used in RDFBasedHistFillingData
// (without the leading underscore): "2mu4", "mu4_mu4noL1".

class PPBaseClass {
public:
    // Map: {run_year (2-digit), trigger} --> crossx_factor = 1/L_int [pb^{-1}]
    static const std::map<std::pair<int,std::string>, double>& CrossxFactorMap() {
        static const std::map<std::pair<int,std::string>, double> m = {
            {{17, "2mu4"},          1./256.8  },   // pp Run 2 (2017), 2mu4,        L = 256.8  pb^{-1}
            {{24, "mu4_mu4noL1"},   1./113.999},   // pp 2024, mu4_mu4noL1,         L = 113.999 pb^{-1}
            {{24, "2mu4"},          1./410.815},   // pp 2024, 2mu4,                L = 410.815 pb^{-1}
        };
        return m;
    }

    // Lookup helper: throws std::runtime_error if (run_year, trig) not in map.
    // run_year: 2-digit (e.g. 17, 24); trig: e.g. "2mu4", "mu4_mu4noL1"
    static double GetCrossxFactor(int run_year, const std::string& trig) {
        run_year %= 2000; // accept 4-digit or 2-digit year
        const auto& m = CrossxFactorMap();
        auto it = m.find({run_year, trig});
        if (it == m.end()) {
            throw std::runtime_error(
                "PPBaseClass::GetCrossxFactor: no entry for (run_year=" +
                std::to_string(run_year) + ", trig=\"" + trig + "\"). "
                "Valid combinations: (17,2mu4), (24,mu4_mu4noL1), (24,2mu4)."
            );
        }
        return it->second;
    }
};
