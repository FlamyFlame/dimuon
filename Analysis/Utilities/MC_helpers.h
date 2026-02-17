#pragma once

double SumMetaNentriesBeforeFilter(const std::vector<std::string>& files,
                                  const char* metaTreeName = "meta_tree",
                                  const char* brName       = "nentries_before_cuts")
{
    double sum = 0.0;

    for (const auto& fn : files) {
        TFile f(fn.c_str(), "READ");
        if (f.IsZombie()) throw std::runtime_error("Cannot open file: " + fn);

        auto* t = dynamic_cast<TTree*>(f.Get(metaTreeName));
        if (!t) throw std::runtime_error("Missing meta tree '" + std::string(metaTreeName) +
                                         "' in file: " + fn);

        Long64_t v = 0;
        if (t->SetBranchAddress(brName, &v) < 0)
            throw std::runtime_error("Missing branch '" + std::string(brName) +
                                     "' in meta tree for file: " + fn);

        const auto n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            sum += static_cast<double>(v);
        }
    }
    return sum;
}
