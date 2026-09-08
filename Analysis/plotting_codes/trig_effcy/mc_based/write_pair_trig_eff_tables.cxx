// =================================================================================================
// write_pair_trig_eff_tables.cxx
//
// THE SINGLE-VALUE PAIR 2mu4 EFFICIENCY AS CSV
// (docs/tracking/mc_trigeff_single_value_pair_eff.md §2 / §PP-1; user request 2026-09-08)
//
// Pure formatting: every number is READ from `pair_trig_eff_<label><wp>.root` exactly as
// FillMCTrigEffPairEff wrote it. Nothing is recomputed here, so a value in a CSV is the same
// object as the value in that file and in the doc's R1 table.
//
// THREE FILES, all with pair-pT ROWS and |eta^pair| COLUMNS on the canonical axes:
//   single_value_pair_eff_opposite_sign.csv   eps^pair and K, both mass windows, + delivery status
//   single_value_pair_eff_same_sign.csv       the same for same-sign pairs
//   single_value_pair_stats_same_sign.csv     the same-sign STATISTICS in the SIGNAL window
//                                             (raw, unweighted pair counts -- the statistical
//                                             reach the weighted efficiencies hide)
//
// THE DELIVERY STATUS COLUMN IS THE READER'S OWN PREDICATE, not a copy of it: the status is asked
// of PairTrigEffEvaluator at each cell's centre, so a cell marked `delivered` in a CSV is exactly a
// cell a consumer can Eval(). A second implementation of the gate is how the two would drift apart.
//
// Usage (from Analysis/plotting_codes/trig_effcy/mc_based/):
//   root -l -b -q 'write_pair_trig_eff_tables.cxx+("pp_full", true)'    // Tight  (nominal)
//   root -l -b -q 'write_pair_trig_eff_tables.cxx+("pp_full", false)'   // Medium (WP systematic)
// =================================================================================================

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TSystem.h>

#include "dr_correction_sample_cfg.h"
#include "../../../Utilities/PairTrigEffEvaluator.h"

namespace {

// The |eta^pair| column header for group g, written from the axis the numbers were binned on.
std::string EtaCol(const std::vector<double>& e, int g)
{
    return Form("|eta|=%.1f-%.1f", e[g], e[g + 1]);
}

TH2D* Get(TFile* f, const std::string& n)
{
    auto* h = dynamic_cast<TH2D*>(f->Get(n.c_str()));
    if (!h) throw std::runtime_error("write_pair_trig_eff_tables: missing '" + n + "' in "
                                     + f->GetName());
    return h;
}

}  // namespace

// =================================================================================================
void write_pair_trig_eff_tables(const std::string& sample = "pp_full", bool use_tight_wp = true)
{
    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_text = use_tight_wp ? "Tight" : "Medium";

    const std::string in_path = PairTrigEff::FileName(cfg.mc_dir, cfg.mc_label, wp_suf);
    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie())
        throw std::runtime_error("write_pair_trig_eff_tables: cannot open " + in_path
                                 + " -- run FillMCTrigEffPairEff first");
    auto* prov = dynamic_cast<TNamed*>(fin->Get("provenance"));

    const std::vector<double> pt  = PairTrigEff::PairPtEdges();
    const std::vector<double> eta = PairTrigEff::AbsEtaGroups().edges;
    const int npt = (int)pt.size() - 1, neta = (int)eta.size() - 1;
    const int first_delivered = PairTrigEff::FirstDeliveredPtBin();

    const std::string outdir = cfg.out_base + "single_value_pair_eff_tables/";
    gSystem->mkdir(outdir.c_str(), kTRUE);

    // ---------------------------------------------------------------- the two efficiency tables
    for (const auto& S : PairTrigEff::Signs()) {
        // One evaluator per (window, form) so the status column is the reader's own answer.
        std::map<std::string, PairTrigEffEvaluator> ev;
        for (const auto& W : PairTrigEff::Windows()) {
            ev[W.token].Load(in_path, S.token, W.token, PairTrigEffEvaluator::ApplyForm::kPure);
        }

        const std::string path = outdir + "single_value_pair_eff_"
                               + (S.token == "os" ? "opposite_sign" : "same_sign") + wp_suf + ".csv";
        std::ofstream out(path);
        if (!out) throw std::runtime_error("write_pair_trig_eff_tables: cannot write " + path);
        out << std::setprecision(6) << std::fixed;

        out << "# single-value pair 2mu4 efficiency, " << cfg.sample_text << ", " << wp_text
            << " muons, " << S.text << " pairs\n"
            << "# eps = sum_{pairs firing 2mu4} w / sum_{all pairs} w   (PURE: replaces the whole "
               "per-pair trigger weight)\n"
            << "# K   = sum_{pairs firing 2mu4} w/(eps_MC1 eps_MC2) / sum_{all pairs} w   "
               "(CALIBRATED: multiplies the two single-muon efficiencies)\n"
            << "# mass windows: sig = " << PairTrigEff::Window("sig").lo << "-"
            << PairTrigEff::Window("sig").hi << " GeV (the single-b signal window); wide = "
            << PairTrigEff::Window("wide").lo << "-" << PairTrigEff::Window("wide").hi
            << " GeV (the template-fit window)\n"
            << "# status: `delivered` = PairTrigEffEvaluator will return this cell; anything else "
               "is the reason it refuses (>= " << PairTrigEff::MinCellPairs()
            << " raw pairs and a value in (0, " << PairTrigEff::MaxDeliveredValue()
            << "] are required).\n"
            << "# the measurement is made in all " << npt << " pair-pT bins; only bins "
            << first_delivered << "-" << npt << " are DELIVERED (the request's three highest cells, "
            << "of which the first is the control region).\n"
            << "# source: " << in_path << "\n";
        if (prov) out << "# provenance: " << prov->GetTitle() << "\n";

        out << "pair_pt_lo_GeV,pair_pt_hi_GeV";
        for (int g = 0; g < neta; ++g) {
            const std::string c = EtaCol(eta, g);
            for (const auto& W : PairTrigEff::Windows())
                out << ",eps_" << W.token << "[" << c << "],eps_" << W.token << "_err[" << c << "]"
                    << ",K_" << W.token << "[" << c << "],K_" << W.token << "_err[" << c << "]";
            for (const auto& W : PairTrigEff::Windows())
                out << ",status_" << W.token << "[" << c << "]";
        }
        out << "\n";

        for (int ix = 1; ix <= npt; ++ix) {
            out << pt[ix - 1] << "," << pt[ix];
            const double pt_c = 0.5 * (pt[ix - 1] + pt[ix]);
            for (int iy = 1; iy <= neta; ++iy) {
                const double eta_c = 0.5 * (eta[iy - 1] + eta[iy]);
                for (const auto& W : PairTrigEff::Windows()) {
                    const TH2D* he = Get(fin, PairTrigEff::HistName("eps", S.token, W.token));
                    const TH2D* hk = Get(fin, PairTrigEff::HistName("k",   S.token, W.token));
                    out << "," << he->GetBinContent(ix, iy) << "," << he->GetBinError(ix, iy)
                        << "," << hk->GetBinContent(ix, iy) << "," << hk->GetBinError(ix, iy);
                }
                for (const auto& W : PairTrigEff::Windows())
                    out << "," << PairTrigEffEvaluator::StatusText(
                                      ev.at(W.token).Status(pt_c, eta_c));
            }
            out << "\n";
        }
        out.close();
        std::cout << "  wrote " << path << std::endl;
    }

    // ---------------------------------------------------------------- the same-sign statistics
    // RAW, UNWEIGHTED counts in the SIGNAL mass window -- what the same-sign efficiency is built on,
    // and therefore what decides whether it can be measured at all. Written for same sign because
    // that is the sparse one: the background estimate is OS - SS, so the same-sign reach is the
    // binding constraint on the whole procedure at high pair pT.
    {
        const std::string win = "sig";
        const TH2D* na = Get(fin, PairTrigEff::HistName("nraw",     "ss", win));
        const TH2D* np = Get(fin, PairTrigEff::HistName("nrawpass", "ss", win));
        const std::string path = outdir + "single_value_pair_stats_same_sign" + wp_suf + ".csv";
        std::ofstream out(path);
        if (!out) throw std::runtime_error("write_pair_trig_eff_tables: cannot write " + path);

        out << "# same-sign pair STATISTICS in the single-b signal mass window "
            << PairTrigEff::Window(win).lo << " < m_mumu < " << PairTrigEff::Window(win).hi
            << " GeV, " << cfg.sample_text << ", " << wp_text << " muons\n"
            << "# RAW, UNWEIGHTED pair counts on the same cells as the efficiency tables.\n"
            << "# n_all  = pairs passing the MC trigger-efficiency pair selection in the cell\n"
            << "# n_2mu4 = of those, the ones that fired 2mu4   (eps^pair is their WEIGHTED ratio)\n"
            << "# source: " << in_path << "\n"
            << "pair_pt_lo_GeV,pair_pt_hi_GeV";
        for (int g = 0; g < neta; ++g)
            out << ",n_all[" << EtaCol(eta, g) << "],n_2mu4[" << EtaCol(eta, g) << "]";
        out << ",n_all[all |eta|],n_2mu4[all |eta|]\n";

        double tot_a = 0., tot_p = 0.;
        for (int ix = 1; ix <= npt; ++ix) {
            out << pt[ix - 1] << "," << pt[ix];
            double ra = 0., rp = 0.;
            for (int iy = 1; iy <= neta; ++iy) {
                const double a = na->GetBinContent(ix, iy), p = np->GetBinContent(ix, iy);
                out << "," << (long long)a << "," << (long long)p;
                ra += a; rp += p;
            }
            out << "," << (long long)ra << "," << (long long)rp << "\n";
            tot_a += ra; tot_p += rp;
        }
        out << pt.front() << "," << pt.back();
        for (int iy = 1; iy <= neta; ++iy)
            out << "," << (long long)na->Integral(1, npt, iy, iy)
                << "," << (long long)np->Integral(1, npt, iy, iy);
        out << "," << (long long)tot_a << "," << (long long)tot_p << "\n";
        out.close();
        std::cout << "  wrote " << path << "  (last row = all pair-pT bins)" << std::endl;
    }

    fin->Close();
    std::cout << "done." << std::endl;
}
