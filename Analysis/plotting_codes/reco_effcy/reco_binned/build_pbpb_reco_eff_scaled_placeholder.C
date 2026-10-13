// =================================================================================================
// build_pbpb_reco_eff_scaled_placeholder.C
//
// Builds the Pb+Pb SCALED reco-efficiency placeholder (docs/tracking/crossx_raa_reco_binned_propagation_20261013.md,
// Physics Procedure §2, §3c, D4; evaluator Utilities/PbPbRecoEffScaledPlaceholderEvaluator.h):
//
//   eps_PbPb(j, c) = C_Run2(p, c) * eps^{reco,SF}_j
//   C_Run2(p, c)   = E_Run2(p, c) / E_Run2(p, 60-80 %)          (stored as h_E_run2_<wp>)
//   E_Run2(p, c)   = sum_{T, pT bin p} w * eps_Run2(pT_1, q*eta_1; c) * eps_Run2(pT_2, q*eta_2; c) / sum_{T, p} w
//   G_j            = eps^{reco,SF}_j                             (stored as h_G_<wp>)
//
// D4 (user 2026-10-13, "anchor to pp24"): Run 2 supplies ONLY the centrality dependence relative to its most
// peripheral interval; level, pair-pT, eta^pair (gap acceptance) and m dependence and the SFs all come from the
// pp24 eps^{reco,SF} -- so 60-80 % Pb+Pb = pp24 exactly and no Run-2-vs-Run-3 level factor enters R_AA.
//
//   T            : pp24 fullsim truth signal pairs (truth mass window, truth pair pT, truth (pT, eta^pair, m)
//                  inside the canonical axes) -- the SAME definition as fill_pp24_fullsim_pair_bin_migration.C;
//                  binned in TRUTH pair pT; eps_Run2 at the TRUTH muon (pT, q*eta) = the Run-2 single-muon
//                  MC fits via RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder (the function the
//                  crossx used), centrality = the interval's representative value.
//   eps^{reco,SF}: reco_eff/pair_reco_eff_reco_binned_sf_pp24_full.root (3D, m-integrated 2D, inclusive)
//   eps^{nogap}_p: written for reference only (h_eps_nogap_<wp>; the D1 normalisation it served is superseded by D4)
//
// Output: reco_eff/pbpb_reco_eff_scaled_placeholder.root (in the pp24 fullsim reco_eff dir, next to its inputs):
//   h_E_run2_<wp> (TH2D C_Run2, pT x centrality interval), h_E_run2_abs_<wp> (absolute E_Run2, reference),
//   h_G_<wp> (TH3D), h_G_<wp>_minv_integrated (TH2D), G_inclusive_<wp>, h_eps_nogap_<wp> (TH1D), definition marker.
//
// Usage (needs the compiled data class for the Run-2 lookup):
//   cd Analysis/plotting_codes/reco_effcy/reco_binned
//   root -l -b -e '.L ../../../RDFBasedHistFilling/RDFBasedHistFillingPbPb.cxx+' -q build_pbpb_reco_eff_scaled_placeholder.C
// =================================================================================================
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../Utilities/PbPbRecoEffScaledPlaceholderEvaluator.h"

// The Run-2 lookup the crossx uses is a protected static of the data class: re-export it, do not copy it.
struct RunTwoRecoEffPlaceholder : RDFBasedHistFillingData {
    using RDFBasedHistFillingData::OpenRecoEffPlaceholderFile;
    using RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder;
};

void build_pbpb_reco_eff_scaled_placeholder()
{
    using EV = PairRecoEffRecoBinnedEvaluator;
    namespace D = PbPbRecoEffScaledPlaceholderDefinition;
    const std::vector<double> ept = EV::PtEdges(), eeta = EV::EtaEdges(), em = EV::MinvEdges();
    const std::vector<double> ec = D::CentralityEdges();
    const std::vector<int> crep = D::CentralityRepresentatives();
    const int npt = ept.size() - 1, neta = eeta.size() - 1, nm = em.size() - 1, nc = ec.size() - 1;
    if ((int)crep.size() != nc) throw std::runtime_error("centrality representatives / edges mismatch");

    const std::string sample_dir = FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/false);
    const std::string reco_dir   = FullSimRecoEffDir(sample_dir);
    const std::string ntup       = sample_dir + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root";
    const std::string sf_prod    = reco_dir + "pair_reco_eff_reco_binned_sf_pp24_full.root";
    const std::string sf_cells   = reco_dir + "pair_reco_sf_cells_pp24_full.root";
    const std::string out_path   = reco_dir + "pbpb_reco_eff_scaled_placeholder.root";

    // ---- E_Run2: eta^pair-averaged Run-2 pair placeholder over T ----
    TFile out(out_path.c_str(), "RECREATE");
    for (const std::string wp : {"tight", "medium"}) {
        RunTwoRecoEffPlaceholder::OpenRecoEffPlaceholderFile(wp == "tight");   // sets the WP of the key lookup
        ROOT::RDataFrame df0("muon_pair_tree_sign2", ntup);
        TAxis apt(npt, ept.data()), aeta(neta, eeta.data()), am(nm, em.data());
        auto in_ax = [](const TAxis& a, float v) { const int b = a.FindFixBin(v); return b >= 1 && b <= a.GetNbins(); };
        const std::string sig_truth = ParamsSet::SignalMinvCutExpr("truth_minv") + " && "
                                    + ParamsSet::SignalPairPtCutExpr("truth_pair_pt");
        auto dT = df0.Filter(sig_truth)
                     .Filter([&](float pt, float eta, float m) { return in_ax(apt, pt) && in_ax(aeta, eta) && in_ax(am, m); },
                             {"truth_pair_pt", "truth_pair_eta", "truth_minv"});
        auto hw = dT.Histo1D({"hw", "", npt, ept.data()}, "truth_pair_pt", "weight");
        std::vector<ROOT::RDF::RResultPtr<TH1D>> he;
        std::vector<ROOT::RDF::RResultPtr<ULong64_t>> nbad;
        for (int k = 0; k < nc; ++k) {
            const int c = crep[k];
            auto dk = dT.Define("eps_run2_pair", [c](float pt1, float eta1, int q1, float pt2, float eta2, int q2) {
                            const float e1 = RunTwoRecoEffPlaceholder::EvaluateSingleMuonRecoEffPlaceholder(c, pt1, q1 * eta1);
                            const float e2 = RunTwoRecoEffPlaceholder::EvaluateSingleMuonRecoEffPlaceholder(c, pt2, q2 * eta2);
                            return (e1 > 0 && e2 > 0) ? (double)e1 * e2 : -1.0;
                        }, {"m1.truth_pt", "m1.truth_eta", "m1.truth_charge", "m2.truth_pt", "m2.truth_eta", "m2.truth_charge"})
                        .Define("w_eps", "weight * eps_run2_pair");
            nbad.push_back(dk.Filter("eps_run2_pair < 0").Count());
            he.push_back(dk.Histo1D({Form("he%d", k), "", npt, ept.data()}, "truth_pair_pt", "w_eps"));
        }
        TH2D hEabs(("h_E_run2_abs_" + wp).c_str(),
                "E_{Run2} = #LT#varepsilon_{Run2}(1)#varepsilon_{Run2}(2)#GT_{T};p_{T}^{pair} [GeV];centrality [%]",
                npt, ept.data(), nc, ec.data());
        for (int k = 0; k < nc; ++k) {
            if (*nbad[k] != 0)
                throw std::runtime_error("Run-2 placeholder lookup failed for " + std::to_string(*nbad[k])
                                         + " T pairs (centrality interval " + std::to_string(k) + ")");
            for (int i = 1; i <= npt; ++i) {
                const double w = hw->GetBinContent(i);
                if (!(w > 0)) throw std::runtime_error("empty T pT bin " + std::to_string(i));
                hEabs.SetBinContent(i, k + 1, he[k]->GetBinContent(i) / w);
            }
        }
        // C_Run2 = E_Run2 / E_Run2(most peripheral interval), per pT bin (D4).
        TH2D hE(("h_E_run2_" + wp).c_str(),
                "C_{Run2} = E_{Run2}(p_{T},c)/E_{Run2}(p_{T},60-80%);p_{T}^{pair} [GeV];centrality [%]",
                npt, ept.data(), nc, ec.data());
        for (int i = 1; i <= npt; ++i)
            for (int k = 1; k <= nc; ++k)
                hE.SetBinContent(i, k, hEabs.GetBinContent(i, k) / hEabs.GetBinContent(i, nc));

        // ---- G = eps^{reco,SF}_j / eps^{nogap}_p ----
        std::unique_ptr<TFile> fp(TFile::Open(sf_prod.c_str(), "READ"));
        std::unique_ptr<TFile> fc(TFile::Open(sf_cells.c_str(), "READ"));
        if (!fp || fp->IsZombie() || !fc || fc->IsZombie()) throw std::runtime_error("cannot open SF product / cells");
        auto* def = dynamic_cast<TNamed*>(fp->Get(PairRecoEffRecoBinnedDefinition::Key()));
        if (!def || std::string(def->GetTitle()) != PairRecoEffRecoBinnedDefinition::ValueSF())
            throw std::runtime_error(sf_prod + " is not the SF-corrected product");
        auto* e3 = dynamic_cast<TH3D*>(fp->Get(("h_eps_reco_binned_" + wp).c_str()));
        auto* e2 = dynamic_cast<TH2D*>(fp->Get(("h_eps_reco_binned_" + wp + "_minv_integrated").c_str()));
        auto* ei = dynamic_cast<TParameter<double>*>(fp->Get(("eps_reco_binned_inclusive_" + wp).c_str()));
        auto* hRn = dynamic_cast<TH1D*>(fc->Get(("h_pair_pt_Rnogap_" + wp).c_str()));
        auto* hTp = dynamic_cast<TH1D*>(fc->Get("h_pair_pt_T"));
        if (!e3 || !e2 || !ei || !hRn || !hTp) throw std::runtime_error("missing G inputs");
        EV::CheckAxis(e3->GetXaxis(), ept, "pair-pT", sf_prod);
        EV::CheckAxis(e3->GetYaxis(), eeta, "pair-eta", sf_prod);
        EV::CheckAxis(e3->GetZaxis(), em, "m_mumu", sf_prod);
        EV::CheckAxis(hRn->GetXaxis(), ept, "pair-pT", sf_cells);
        EV::CheckAxis(hTp->GetXaxis(), ept, "pair-pT", sf_cells);
        // T per pT bin from the cells file must equal the T used for E (same tree, same definition).
        for (int i = 1; i <= npt; ++i)
            if (std::fabs(hTp->GetBinContent(i) - hw->GetBinContent(i)) > 1e-6 * hw->GetBinContent(i))
                throw std::runtime_error("T per pT bin differs between the E pass and the SF cells file");

        TH1D hng(("h_eps_nogap_" + wp).c_str(), "#varepsilon^{nogap};p_{T}^{pair} [GeV]", npt, ept.data());
        TH3D hG(("h_G_" + wp).c_str(), "G = #varepsilon^{reco,SF}_{j};p_{T}^{pair} [GeV];#eta^{pair};m_{#mu#mu} [GeV]",
                npt, ept.data(), neta, eeta.data(), nm, em.data());
        TH2D hG2((("h_G_" + wp) + "_minv_integrated").c_str(), ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ept.data(), neta, eeta.data());
        for (int i = 1; i <= npt; ++i) {
            const double ng = hRn->GetBinContent(i) / hTp->GetBinContent(i);
            hng.SetBinContent(i, ng);
            for (int j = 1; j <= neta; ++j) {
                if (e2->GetBinContent(i, j) > 0) hG2.SetBinContent(i, j, e2->GetBinContent(i, j));
                for (int k = 1; k <= nm; ++k)
                    if (e3->GetBinContent(i, j, k) > 0) hG.SetBinContent(i, j, k, e3->GetBinContent(i, j, k));
            }
        }
        const double ng_incl = hRn->Integral() / hTp->Integral();
        out.cd();
        hE.Write(); hEabs.Write(); hG.Write(); hG2.Write(); hng.Write();
        TParameter<double>(("G_inclusive_" + wp).c_str(), ei->GetVal()).Write();

        std::cout << "[" << wp << "] eps^nogap inclusive (reference) " << ng_incl << ", eps^{reco,SF} inclusive "
                  << ei->GetVal() << "\n  pT bin: eps^nogap | E_Run2 (abs) per centrality | C_Run2 per centrality\n";
        for (int i = 1; i <= npt; ++i) {
            std::cout << Form("  [%6.2f,%6.2f): %.4f |", ept[i - 1], ept[i], hng.GetBinContent(i));
            for (int k = 1; k <= nc; ++k) std::cout << Form(" %.4f", hEabs.GetBinContent(i, k));
            std::cout << " |";
            for (int k = 1; k <= nc; ++k) std::cout << Form(" %.4f", hE.GetBinContent(i, k));
            std::cout << "\n";
        }
        std::cout << "  G (m-integrated) per pT bin over eta^pair cells:\n";
        for (int i = 1; i <= npt; ++i) {
            std::cout << Form("  [%6.2f,%6.2f):", ept[i - 1], ept[i]);
            for (int j = 1; j <= neta; ++j) std::cout << Form(" %.3f", hG2.GetBinContent(i, j));
            std::cout << "\n";
        }
    }
    out.cd();
    TNamed(D::Key(), D::Value()).Write();
    TNamed("inputs", (ntup + " | " + sf_prod + " | " + sf_cells).c_str()).Write();
    out.Close();

    PbPbRecoEffScaledPlaceholderEvaluator ev;
    ev.Load(out_path, true);
    std::cout << "[readback] eps_PbPb(20 GeV, 0.0, 1.5 GeV, 5 %) = " << ev.Eval(20., 0., 1.5, 5.)
              << ", (20, 0, 1.5, 70 %) = " << ev.Eval(20., 0., 1.5, 70.) << "\n[done] " << out_path << "\n";
}
