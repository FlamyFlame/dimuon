// For each pair pT bin, list pTHat contributors with >10% cross-section share,
// ranked by cross-section, with their crossx% and full-sample error fraction%.
// Outputs a Markdown file. Runs for truth and reco pair pT, both binnings.

#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>

void write_table(std::ofstream& out,
                 const std::vector<TH1D*>& hists,
                 const std::vector<double>& edges,
                 const std::array<std::string, 6>& kn_names,
                 const std::array<double, 6>& scale_factors,
                 const std::string& section_title) {

    const int nkn   = 6;
    const int nbins  = hists[0]->GetNbinsX();

    out << "## " << section_title << "\n\n";

    for (int ib = 1; ib <= nbins; ib++) {
        double lo = edges[ib - 1], hi = edges[ib];

        double tot_crossx = 0.;
        for (int ikn = 0; ikn < nkn; ikn++)
            tot_crossx += hists[ikn]->GetBinContent(ib);
        if (tot_crossx <= 0.) continue;

        double tot_err2 = 0.;
        for (int ikn = 0; ikn < nkn; ikn++) {
            double e = hists[ikn]->GetBinError(ib) / std::sqrt(scale_factors[ikn]);
            tot_err2 += e * e;
        }
        const double tot_err = std::sqrt(tot_err2);

        struct Entry { int ikn; double xfrac; double efrac; };
        std::vector<Entry> entries;
        for (int ikn = 0; ikn < nkn; ikn++) {
            double c = hists[ikn]->GetBinContent(ib);
            double xfrac = c / tot_crossx;
            if (xfrac < 0.10) continue;
            double e = hists[ikn]->GetBinError(ib) / std::sqrt(scale_factors[ikn]);
            entries.push_back({ikn, xfrac, tot_err > 0 ? e / tot_err : 0.});
        }
        if (entries.empty()) continue;
        std::sort(entries.begin(), entries.end(),
                  [](const Entry& a, const Entry& b){ return a.xfrac > b.xfrac; });

        std::ostringstream line;
        line << std::fixed << std::setprecision(1);
        line << "* **" << lo << " – " << hi << " GeV**: ";
        for (size_t i = 0; i < entries.size(); i++) {
            const auto& e = entries[i];
            if (i > 0) line << ", ";
            line << "pTHat " << kn_names[e.ikn] << " GeV"
                 << " (crossx: " << e.xfrac * 100. << "%"
                 << ", err: "    << e.efrac * 100. << "%)";
        }
        out << line.str() << "\n";
    }
    out << "\n";
}

void make_kn_contributor_table() {

    const std::string input_file =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/"
        "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts.root";
    const std::string output_md =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/plots/"
        "kn_contributor_table.md";

    const int nkn = 6;
    const std::array<std::string, nkn> kn_names = {
        "8-14", "14-24", "24-40", "40-70", "70-125", "125-300"
    };
    const std::array<double, nkn> scale_factors = {51., 90., 60., 30., 8., 1.};

    // Two binnings
    const std::array<int,    2> nbins_arr = {20, 25};
    const std::array<double, 2> xmax_arr  = {120., 150.};

    // Two variables
    const std::array<std::string, 2> var_names   = {"truth_pair_pt", "pair_pt"};
    const std::array<std::string, 2> filter_strs = {
        "from_same_b", "from_same_b && pair_pass_medium"
    };
    const std::array<std::string, 2> var_labels  = {
        "truth pair pT", "reco pair pT"
    };

    std::ofstream out(output_md);
    out << "# pTHat contributor table — Pythia fullsim pp24, single-b (from_same_b)\n";
    out << "Full-sample scale factors: "
        << "kn0(8-14)=x51, kn1(14-24)=x90, kn2(24-40)=x60, "
        << "kn3(40-70)=x30, kn4(70-125)=x8, kn5(125-300)=x1\n";
    out << "Only contributors with >10% of total cross-section shown, ranked by crossx.\n";
    out << "err = σ_err,kn / σ_err,total (quadrature fraction).\n\n";

    for (int iv = 0; iv < 2; iv++) {
        for (int ib2 = 0; ib2 < 2; ib2++) {
            const int    nbins = nbins_arr[ib2];
            const double xmin  = 8., xmax = xmax_arr[ib2];

            std::vector<double> edges(nbins + 1);
            const double lmin = std::log(xmin), lmax = std::log(xmax);
            for (int i = 0; i <= nbins; i++)
                edges[i] = std::exp(lmin + i * (lmax - lmin) / nbins);

            std::vector<TH1D*> hists(nkn);
            for (int ikn = 0; ikn < nkn; ikn++) {
                const std::string tree =
                    "muon_pair_tree_kin" + std::to_string(ikn) + "_sign2";
                ROOT::RDataFrame df(tree, input_file);
                auto hptr = df.Filter(filter_strs[iv])
                              .Histo1D(ROOT::RDF::TH1DModel{
                                  ("htab_" + std::to_string(iv) + "_"
                                   + std::to_string(ib2) + "_kn"
                                   + std::to_string(ikn)).c_str(),
                                  "", nbins, edges.data()
                              }, var_names[iv], "weight");
                hists[ikn] = (TH1D*)hptr->Clone();
                hists[ikn]->SetDirectory(nullptr);
                hists[ikn]->Scale(1., "width");
            }

            const std::string title = var_labels[iv]
                + ", " + std::to_string(nbins) + " bins, 8–"
                + std::to_string((int)xmax) + " GeV";
            write_table(out, hists, edges, kn_names, scale_factors, title);

            for (auto* h : hists) delete h;
        }
    }

    printf("Written: %s\n", output_md.c_str());
}
