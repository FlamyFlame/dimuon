// =================================================================================================
// build_pp24_fullsim_pair_reco_eff_reco_binned.C
//
// The RECO-BINNED pair reco efficiency of the 3D bin-by-bin correction in (pair pT, eta^pair,
// m_mumu), docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md Physics Procedure §2, §3a:
//
//   eps^reco_j = N[R_j] / N[T_j]   (numerator binned in RECO, denominator in TRUTH, same axes;
//                                   N[R_j] = sum_i M_ij: truth in T, out-of-T feed-in EXCLUDED
//                                   (template-fit fake, user decision D7 2026-10-13); NOT capped)
//
// Input : reco_eff/pair_bin_migration_pp24_full.root (fill_pp24_fullsim_pair_bin_migration.C), built
//         from the ntuple-processing output with the RDF class's expression builders. Its cell
//         index histograms are re-verified here: cell count vs the canonical axes, the 1D axis
//         histograms' edges vs the canonical source, the definition TNamed, and N[T], N[T and R]
//         re-closed cell by cell against the RDF-class 3D legs of the CURRENT fullsim hist file.
// Output: reco_eff/pair_reco_eff_reco_binned_pp24_full.root, per WP (tight, medium):
//   h_eps_reco_binned_<wp>                    TH3D eps^reco (pT, eta^pair, m), errors per reco_binned_cells.h
//   h_eps_reco_binned_<wp>_minv_integrated    TH2D eps^reco (pT, eta^pair) = sum_m R / sum_m T (A over the column)
//   eps_reco_binned_inclusive_<wp>            TParameter<double> (+ _err)
//   h_NR_<wp>, h_NT (TH3D, Sum w^2 kept), h_NA_<wp> (TH3D diagonal M_jj), h_NF_<wp> (out-of-T feed-in, info)
//   + definition marker (PairRecoEffRecoBinnedDefinition) and the input path.
//
// SF mode (with_sf = true; docs/tracking/crossx_raa_reco_binned_propagation_20261013.md §2, §3b):
//   input  reco_eff/pair_reco_sf_cells_pp24_full.root (fill_pp24_fullsim_pair_reco_sf_cells.C; its unweighted
//          cells were verified identical to the migration file there), cell set "sf_<wp>" = R, TR, M weighted
//          by w*SF_1*SF_2, T unweighted; output reco_eff/pair_reco_eff_reco_binned_sf_pp24_full.root with the
//          same keys and the ValueSF marker. Value = N^SF[R]/N[T] (exact); error = the MC-only A/B/D RELATIVE
//          error carried over (see the RatioOverCells lambda below).
//   BLOCKED since D7 (2026-10-13): the SF product is the crossx input and still carries the out-of-T feed-in
//   (ValueSF); rebuilding it with the D7 numerator is the crossx propagation step (user: not before the pbpb26
//   ntuple rerun), which must update ValueSF together with it. with_sf = true THROWS until then.
//
// Usage (from Analysis/plotting_codes/reco_effcy/reco_binned/):
//   root -l -b -q 'build_pp24_fullsim_pair_reco_eff_reco_binned.C+'            (MC only)
//   root -l -b -q 'build_pp24_fullsim_pair_reco_eff_reco_binned.C+(true)'      (SF-corrected)
// =================================================================================================
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TSystem.h>

#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../Utilities/PairRecoEffRecoBinnedEvaluator.h"
#include "reco_binned_cells.h"

namespace {

void CompareCells(TFile* fh, const TH1D* mine, const std::string& name3d)
{
    auto* h = dynamic_cast<TH3D*>(fh->Get(name3d.c_str()));
    if (!h) throw std::runtime_error("closure: missing " + name3d);
    const int nx = h->GetNbinsX(), ny = h->GetNbinsY(), nz = h->GetNbinsZ();
    if (nx * ny * nz != mine->GetNbinsX()) throw std::runtime_error("closure: cell count differs for " + name3d);
    double max_rel = 0;
    for (int ix = 1; ix <= nx; ++ix)
        for (int iy = 1; iy <= ny; ++iy)
            for (int iz = 1; iz <= nz; ++iz) {
                const int k = ((ix - 1) * ny + (iy - 1)) * nz + (iz - 1);
                const double a = h->GetBinContent(ix, iy, iz), b = mine->GetBinContent(k + 1);
                if (a == 0 && b == 0) continue;
                max_rel = std::max(max_rel, std::fabs(a - b) / std::max(std::fabs(a), std::fabs(b)));
            }
    std::cout << "[closure] " << mine->GetName() << " vs " << name3d << ": max |rel diff| = " << max_rel << "\n";
    if (max_rel > 1e-6) throw std::runtime_error("closure FAILED: " + std::string(mine->GetName()));
}

}  // namespace

void build_pp24_fullsim_pair_reco_eff_reco_binned(bool with_sf = false)
{
    if (with_sf)
        throw std::runtime_error("build_pp24_fullsim_pair_reco_eff_reco_binned: the SF product (crossx input) is NOT rebuilt "
                                 "under the D7 numerator before the crossx propagation (user 2026-10-13); update "
                                 "PairRecoEffRecoBinnedDefinition::ValueSF together with it -- "
                                 "docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md D7");
    using EV = PairRecoEffRecoBinnedEvaluator;
    const std::vector<double> ept = EV::PtEdges(), eeta = EV::EtaEdges(), em = EV::MinvEdges();
    const int npt = ept.size() - 1, neta = eeta.size() - 1, nm = em.size() - 1;

    const std::string sample_dir = FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/false);
    const std::string reco_dir   = FullSimRecoEffDir(sample_dir);
    const std::string in_path    = reco_dir + (with_sf ? "pair_reco_sf_cells_pp24_full.root"
                                                           : "pair_bin_migration_pp24_full.root");
    const std::string hist_path  = sample_dir + "histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root";
    const std::string ntup_path  = sample_dir + (with_sf ? "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root"
                                                         : "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root");
    const std::string out_path   = reco_dir + (with_sf ? "pair_reco_eff_reco_binned_sf_pp24_full.root"
                                                       : "pair_reco_eff_reco_binned_pp24_full.root");
    const std::string set_prefix = with_sf ? "sf_" : "";   // cell set read by LoadMigrationCells

    std::unique_ptr<TFile> fin(TFile::Open(in_path.c_str(), "READ"));
    if (!fin || fin->IsZombie()) throw std::runtime_error("cannot open " + in_path);

    // ---- provenance checks of the migration file ----
    auto* tin = dynamic_cast<TNamed*>(fin->Get("input"));
    if (!tin || std::string(tin->GetTitle()) != ntup_path)
        throw std::runtime_error("migration file input is not the current ntuple: "
                                 + std::string(tin ? tin->GetTitle() : "none"));
    auto* tdef = dynamic_cast<TNamed*>(fin->Get("definition"));
    if (!with_sf && (!tdef || std::string(tdef->GetTitle()).find("R = pair_pass_<wp> && pass_signal_reco") == std::string::npos
              || std::string(tdef->GetTitle()).find("category = all truth OS pairs (_op)") == std::string::npos))
        throw std::runtime_error("migration file definition TNamed missing/unexpected");
    if (with_sf && (!tdef || std::string(tdef->GetTitle()).find("_sf_<wp> = weight*SF_1*SF_2") == std::string::npos))
        throw std::runtime_error("SF cells file definition TNamed missing/unexpected");
    {
        FileStat_t s_in, s_nt;
        gSystem->GetPathInfo(in_path.c_str(), s_in);
        gSystem->GetPathInfo(ntup_path.c_str(), s_nt);
        if (s_in.fMtime < s_nt.fMtime) throw std::runtime_error("cell file is OLDER than the ntuple -- refill it");
        std::cout << "[provenance] migration file mtime " << s_in.fMtime << " > ntuple mtime " << s_nt.fMtime << "\n";
    }
    auto axis_of = [&](const std::string& n) {
        auto* h = dynamic_cast<TH1D*>(fin->Get(n.c_str()));
        if (!h) throw std::runtime_error("missing " + n);
        return h->GetXaxis();
    };
    EV::CheckAxis(axis_of("h_pair_pt_T"),  ept,  "pair-pT",  in_path);
    if (!with_sf) {   // the SF cells file carries only the pT axis histogram (cells are re-verified there)
        EV::CheckAxis(axis_of("h_pair_eta_T"), eeta, "pair-eta", in_path);
        EV::CheckAxis(axis_of("h_minv_T"),     em,   "m_mumu",   in_path);
    }
    {   // re-close against the CURRENT RDF-class hist file
        std::unique_ptr<TFile> fh(TFile::Open(hist_path.c_str(), "READ"));
        if (!fh || fh->IsZombie()) throw std::runtime_error("cannot open " + hist_path);
        const std::string b3 = "h_truth_minv_effcy_vs_truth_pair_eta_coarse_vs_truth_pair_pt_coarse_op";
        CompareCells(fh.get(), (TH1D*)fin->Get("h_cell_minv_T"), b3 + "_pass_signal_truth");
        for (const std::string wp : {"tight", "medium"})   // unweighted TR (present in both input files)
            CompareCells(fh.get(), (TH1D*)fin->Get(("h_cell_minv_TR_" + wp).c_str()),
                         b3 + "_pass_" + wp + "_and_signal_truth_and_reco");
    }

    TFile fout(out_path.c_str(), "RECREATE");
    auto make3 = [&](const std::string& n, const std::string& t) {
        auto* h = new TH3D(n.c_str(), t.c_str(), npt, ept.data(), neta, eeta.data(), nm, em.data());
        h->SetDirectory(nullptr);
        h->Sumw2();
        return h;
    };
    auto make2 = [&](const std::string& n, const std::string& t) {
        auto* h = new TH2D(n.c_str(), t.c_str(), npt, ept.data(), neta, eeta.data());
        h->SetDirectory(nullptr);
        h->Sumw2();
        return h;
    };
    const std::string ax3 = ";reco/truth p_{T}^{pair} [GeV];#eta^{pair};m_{#mu#mu} [GeV]";

    bool wrote_T = false;
    for (const std::string wp : {"tight", "medium"}) {
        const MigrationCells C  = LoadMigrationCells(fin.get(), set_prefix + wp, npt, neta, nm);
        const MigrationCells Cu = LoadMigrationCells(fin.get(), wp, npt, neta, nm);   // unweighted (== C if !with_sf)
        // Value: (A+B)/(A+D) == N[R_S]/N[T_S] algebraically, so the SF-weighted R over the unweighted T is exact.
        // Error: the A/B/D propagation assumes ONE weight per pair in numerator and denominator, which SF-on-R
        // breaks; the SF product is ~constant inside a cell (rms 0.025 per muon), so the MC-only RELATIVE error
        // is carried over: err_SF = eps_SF * err_MC / eps_MC.
        auto RatioOverCells = [&](const MigrationCells&, const std::vector<int>& cells) {
            CellRatio r = ::RatioOverCells(C, cells);
            if (with_sf && r.measured) {
                const CellRatio u = ::RatioOverCells(Cu, cells);
                r.err = u.measured && u.eps > 0 ? r.eps * u.err / u.eps : 0.0;
            }
            return r;
        };
        TH3D* heps = make3("h_eps_reco_binned_" + wp, "#varepsilon^{reco}_{j} = N[R_{j}]/N[T_{j}]" + ax3);
        TH3D* hNR  = make3("h_NR_" + wp, "N[R_{j}] (reco cell, truth in T)" + ax3);
        TH3D* hNF  = make3("h_NF_" + wp, "N[F_{j}] (reco cell, truth outside T: template-fit fake)" + ax3);
        TH3D* hNA  = make3("h_NA_" + wp, "M_{jj} = N[T_{j} #wedge R_{j}]" + ax3);
        TH3D* hNT  = make3("h_NT", "N[T_{j}] (truth cell)" + ax3);
        TH2D* h2   = make2("h_eps_reco_binned_" + wp + "_minv_integrated", ";p_{T}^{pair} [GeV];#eta^{pair}");
        int n_gt1 = 0, n_meas = 0, n_unmeas_T = 0, n_unmeas_R = 0;
        double max_eps = 0, min_eps = 1e9;
        for (int i = 0; i < npt; ++i)
            for (int j = 0; j < neta; ++j) {
                for (int k = 0; k < nm; ++k) {
                    const int c = C.Cell(i, j, k);
                    const CellRatio r = RatioOverCells(C, {c});
                    hNR->SetBinContent(i + 1, j + 1, k + 1, r.nR); hNR->SetBinError(i + 1, j + 1, k + 1, std::sqrt(r.s2R));
                    hNT->SetBinContent(i + 1, j + 1, k + 1, r.nT); hNT->SetBinError(i + 1, j + 1, k + 1, std::sqrt(r.s2T));
                    hNA->SetBinContent(i + 1, j + 1, k + 1, r.nA); hNA->SetBinError(i + 1, j + 1, k + 1, std::sqrt(r.s2A));
                    hNF->SetBinContent(i + 1, j + 1, k + 1, r.nF);
                    if (!r.measured) {
                        if (r.nT <= 0) ++n_unmeas_T; else ++n_unmeas_R;
                        continue;
                    }
                    ++n_meas;
                    heps->SetBinContent(i + 1, j + 1, k + 1, r.eps);
                    heps->SetBinError(i + 1, j + 1, k + 1, r.err);
                    if (r.eps > 1) ++n_gt1;
                    max_eps = std::max(max_eps, r.eps); min_eps = std::min(min_eps, r.eps);
                }
                const CellRatio r2 = RatioOverCells(C, CellsIn(C, i, i, j, j, 0, nm - 1));
                if (r2.measured) { h2->SetBinContent(i + 1, j + 1, r2.eps); h2->SetBinError(i + 1, j + 1, r2.err); }
            }
        std::vector<int> all;
        for (int c = 0; c < npt * neta * nm; ++c) all.push_back(c);
        const CellRatio ri = RatioOverCells(C, all);
        // Reco pairs outside the axes (only reco pT > top edge), truth in T: the reco-index underflow of M.
        const int nc = npt * neta * nm;
        const double r_out = C.M->Integral(1, nc, 0, 0), r_tot = C.M->Integral(1, nc, 0, nc);
        int n2_gt1 = 0;
        for (int i = 1; i <= npt; ++i) for (int j = 1; j <= neta; ++j) if (h2->GetBinContent(i, j) > 1) ++n2_gt1;

        fout.cd();
        heps->Write(); h2->Write(); hNR->Write(); hNA->Write(); hNF->Write();
        if (!wrote_T) { hNT->Write(); wrote_T = true; }
        TParameter<double>(("eps_reco_binned_inclusive_" + wp).c_str(), ri.eps).Write();
        TParameter<double>(("eps_reco_binned_inclusive_" + wp + "_err").c_str(), ri.err).Write();
        TParameter<double>(("frac_R_outside_axes_" + wp).c_str(), r_out / r_tot).Write();

        std::cout << Form("[%s] cells measured %d/%d (N[T]=0: %d, N[R]=0 with N[T]>0: %d); eps>1: %d 3D cells, %d 2D cells; "
                          "3D eps range [%.4f, %.4f]; inclusive eps^reco = %.5f +- %.5f (N[R]=%.5g, N[T]=%.5g, A=%.5g, "
                          "out-of-T feed-in N[F]=%.5g = %.3f %% of R+F); R outside axes: %.3g %% of R\n",
                          wp.c_str(), n_meas, npt * neta * nm, n_unmeas_T, n_unmeas_R, n_gt1, n2_gt1, min_eps, max_eps,
                          ri.eps, ri.err, ri.nR, ri.nT, ri.nA, ri.nF, 100.0 * ri.nF / (ri.nR + ri.nF), 100.0 * r_out / r_tot);
    }
    TNamed(PairRecoEffRecoBinnedDefinition::Key(),
           with_sf ? PairRecoEffRecoBinnedDefinition::ValueSF() : PairRecoEffRecoBinnedDefinition::Value()).Write();
    TNamed("input", in_path.c_str()).Write();
    fout.Close();

    // Read-back through the evaluator: marker + canonical axes + a lookup.
    PairRecoEffRecoBinnedEvaluator ev;
    ev.Load(out_path, true, with_sf);
    std::cout << "[readback] eps(20 GeV, 0.0, 1.5 GeV) = " << ev.Eval(20.0, 0.0, 1.5) << "\n[done] " << out_path << "\n";
}
