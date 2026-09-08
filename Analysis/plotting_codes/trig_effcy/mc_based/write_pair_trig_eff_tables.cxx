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

    const std::vector<double> eta = PairTrigEff::AbsEtaGroups().edges;
    const int neta = (int)eta.size() - 1;

    const std::string outdir = cfg.out_base + "single_value_pair_eff_tables/";
    gSystem->mkdir(outdir.c_str(), kTRUE);

    // ---------------------------------------------------------------- the efficiency tables
    // One CSV per (sign, cell mode). The merged mode is a DIFFERENT pair-pT axis -- seven cells,
    // the last one [72.08, 150) GeV -- so it cannot share a table with the un-merged one without
    // either padding a row or repeating a number as if it were two measurements.
    for (const auto& M : PairTrigEff::CellModes())
    for (const auto& S : PairTrigEff::Signs()) {
        const std::vector<double> pt = PairTrigEff::PairPtEdges(M.token);
        const int npt = (int)pt.size() - 1;
        const int first_delivered = PairTrigEff::FirstDeliveredPtBin(M.token);
        // One evaluator per (window, form) so the status column is the reader's own answer.
        std::map<std::string, PairTrigEffEvaluator> ev;
        for (const auto& W : PairTrigEff::Windows()) {
            ev[W.token].Load(in_path, S.token, W.token, PairTrigEffEvaluator::ApplyForm::kPure,
                             M.token);
        }

        const std::string path = outdir + "single_value_pair_eff_"
                               + (S.token == "os" ? "opposite_sign" : "same_sign")
                               + M.suffix + wp_suf + ".csv";
        std::ofstream out(path);
        if (!out) throw std::runtime_error("write_pair_trig_eff_tables: cannot write " + path);
        out << std::setprecision(6) << std::fixed;

        out << "# single-value pair 2mu4 efficiency, " << cfg.sample_text << ", " << wp_text
            << " muons, " << S.text << " pairs\n"
            << "# eps = sum_{pairs firing 2mu4} w / sum_{all pairs} w   (PURE: replaces the whole "
               "per-pair trigger weight)\n"
            << "# K   = sum_{pairs firing 2mu4} w/(eps_MC1 eps_MC2) / sum_{all pairs} w   "
               "(CALIBRATED: multiplies the two single-muon efficiencies)\n"
            << "# cells: " << M.text << " x the " << neta << " |eta^pair| groups\n"
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
                    const TH2D* he = Get(fin, PairTrigEff::HistName("eps", S.token, W.token, M.token));
                    const TH2D* hk = Get(fin, PairTrigEff::HistName("k",   S.token, W.token, M.token));
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
    for (const auto& M : PairTrigEff::CellModes()) {
        const std::string win = "sig";
        const std::vector<double> pt = PairTrigEff::PairPtEdges(M.token);
        const int npt = (int)pt.size() - 1;
        const TH2D* na = Get(fin, PairTrigEff::HistName("nraw",     "ss", win, M.token));
        const TH2D* np = Get(fin, PairTrigEff::HistName("nrawpass", "ss", win, M.token));
        const std::string path = outdir + "single_value_pair_stats_same_sign" + M.suffix
                               + wp_suf + ".csv";
        std::ofstream out(path);
        if (!out) throw std::runtime_error("write_pair_trig_eff_tables: cannot write " + path);

        out << "# same-sign pair STATISTICS in the single-b signal mass window "
            << PairTrigEff::Window(win).lo << " < m_mumu < " << PairTrigEff::Window(win).hi
            << " GeV, " << cfg.sample_text << ", " << wp_text << " muons\n"
            << "# RAW, UNWEIGHTED pair counts on the same cells as the efficiency tables.\n"
            << "# cells: " << M.text << "\n"
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

    // ------------------------------------------------- the compact pT-merged tables (user request)
    // FOUR files beside the figure they belong to, each a strict matrix of the DELIVERED merged
    // cells: 2 pair-pT rows x 3 |eta^pair| columns. Quantities that need more than one number per
    // cell (the error, the second form, the delivery status) are stacked as further blocks of the
    // SAME shape rather than widened into extra columns, so every block can be read as the matrix
    // it is. Signal mass window only -- that is what the pt_merge comparison is about.
    {
        const std::string win  = "sig";
        const std::string mode = "ptmerge";
        const std::vector<double> pt = PairTrigEff::PairPtEdges(mode);
        const int npt = (int)pt.size() - 1;
        const int lo  = PairTrigEff::FirstDeliveredPtBin(mode);
        const std::string dir = cfg.out_base + "closure/single_value_highpt_comparison/"
                                               "pt_merge_compr/";
        gSystem->mkdir(dir.c_str(), kTRUE);

        auto row_label = [&](int ix) { return Form("%.2f-%.2f", pt[ix - 1], pt[ix]); };
        auto header = [&](std::ofstream& o) {
            o << "pair_pt_GeV";
            for (int g = 0; g < neta; ++g) o << "," << EtaCol(eta, g);
            o << "\n";
        };
        auto block = [&](std::ofstream& o, const char* title, const TH2D* h, bool err, int prec) {
            o << "# " << title << "\n";
            header(o);
            for (int ix = lo; ix <= npt; ++ix) {
                o << row_label(ix);
                for (int iy = 1; iy <= neta; ++iy)
                    o << "," << Form("%.*f", prec,
                                     err ? h->GetBinError(ix, iy) : h->GetBinContent(ix, iy));
                o << "\n";
            }
            o << "\n";
        };

        for (const auto& S : PairTrigEff::Signs()) {
            const std::string tag = (S.token == "os" ? "opposite_sign" : "same_sign");

            // ---- efficiencies
            PairTrigEffEvaluator ev;
            ev.Load(in_path, S.token, win, PairTrigEffEvaluator::ApplyForm::kPure, mode);
            const std::string pe = dir + "pt_merge_pair_eff_" + tag + wp_suf + ".csv";
            std::ofstream oe(pe);
            if (!oe) throw std::runtime_error("write_pair_trig_eff_tables: cannot write " + pe);
            oe << "# single-value pair 2mu4 efficiency, LAST TWO pair-pT CELLS COMBINED\n"
               << "# " << cfg.sample_text << ", " << wp_text << " muons, " << S.text << " pairs\n"
               << "# mass window: " << PairTrigEff::Window(win).lo << " < m_mumu < "
               << PairTrigEff::Window(win).hi << " GeV (the single-b signal window)\n"
               << "# eps = sum_{firing 2mu4} w / sum_{all} w      (PURE: replaces the whole "
                  "per-pair trigger weight)\n"
               << "# K   = sum_{firing 2mu4} w/(eps_MC1 eps_MC2) / sum_{all} w   (CALIBRATED: "
                  "multiplies the two single-muon efficiencies)\n"
               << "# errors are conditional/binomial; `status` is PairTrigEffEvaluator's own "
                  "answer (>= " << PairTrigEff::MinCellPairs() << " raw pairs and a value in (0, "
               << PairTrigEff::MaxDeliveredValue() << "] are required to deliver a cell)\n"
               << "# source: " << in_path << "\n\n";
            block(oe, "eps",     Get(fin, PairTrigEff::HistName("eps", S.token, win, mode)), false, 6);
            block(oe, "eps_err", Get(fin, PairTrigEff::HistName("eps", S.token, win, mode)), true,  6);
            block(oe, "K",       Get(fin, PairTrigEff::HistName("k",   S.token, win, mode)), false, 6);
            block(oe, "K_err",   Get(fin, PairTrigEff::HistName("k",   S.token, win, mode)), true,  6);
            oe << "# status\n";
            header(oe);
            for (int ix = lo; ix <= npt; ++ix) {
                oe << row_label(ix);
                for (int iy = 1; iy <= neta; ++iy)
                    oe << "," << PairTrigEffEvaluator::StatusText(
                                    ev.Status(0.5 * (pt[ix - 1] + pt[ix]),
                                              0.5 * (eta[iy - 1] + eta[iy])));
                oe << "\n";
            }
            oe.close();
            std::cout << "  wrote " << pe << std::endl;

            // ---- statistics
            const std::string ps = dir + "pt_merge_pair_stats_" + tag + wp_suf + ".csv";
            std::ofstream os(ps);
            if (!os) throw std::runtime_error("write_pair_trig_eff_tables: cannot write " + ps);
            os << "# " << S.text << " pair STATISTICS on the same cells as the efficiency table\n"
               << "# " << cfg.sample_text << ", " << wp_text << " muons, mass window "
               << PairTrigEff::Window(win).lo << " < m_mumu < " << PairTrigEff::Window(win).hi
               << " GeV\n"
               << "# RAW, UNWEIGHTED counts -- the number of Bernoulli trials the efficiency and "
                  "its error are built on\n"
               << "# n_all  = pairs passing the MC trigger-efficiency pair selection in the cell\n"
               << "# n_2mu4 = of those, the ones that fired 2mu4\n"
               << "# source: " << in_path << "\n\n";
            block(os, "n_all",  Get(fin, PairTrigEff::HistName("nraw",     S.token, win, mode)),
                  false, 0);
            block(os, "n_2mu4", Get(fin, PairTrigEff::HistName("nrawpass", S.token, win, mode)),
                  false, 0);
            os.close();
            std::cout << "  wrote " << ps << std::endl;
        }
    }

    fin->Close();
    std::cout << "done." << std::endl;
}
