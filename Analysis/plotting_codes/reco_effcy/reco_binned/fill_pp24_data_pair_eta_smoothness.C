// =================================================================================================
// fill_pp24_data_pair_eta_smoothness.C
//
// pp24 DATA eta^pair smoothness test (docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md §3c):
// the crossx selection and trigger weight EXACTLY, reco weight = 1 / eps^reco evaluated at the pair's
// reco (pT, eta^pair, m_mumu) on the CANONICAL coarse axes, eta^pair filled in the FINE plot-only
// binning (pair_eta_smoothness_common.h).
//
// MECHANISM: a subclass of RDFBasedHistFillingPP. It reuses the class's own initialization, input
// selection (SetIOPathsHook: V1 _mindR_0_02 2mu4 ntuple), dataframe construction, WP filter
// (ApplyMuonWorkingPointFilter), trigger-efficiency inputs (OpenEffcyPtFitFile,
// OpenPairEfficiencyInputs) and weight columns (AddPairEfficiencyWeightColumns: w_trig = the hybrid
// pair trigger weight, w_reco = the CURRENT applied truth-binned eps(pT, eta, dR)), and the
// lumi factor pp_crossx_lumi_factor. It overrides ONLY
//   * SetIOPathsHook  -> after the parent's, redirect output_file to THIS test's file (the base
//                        InitOutput RECREATEs output_file, so this is what keeps the crossx hist file
//                        untouched) and point the var1D json to its absolute path;
//   * FillHistograms  -> book this test's histograms instead of the crossx/generic families.
// Configuration = run_crossx_hist_filling_pp24.sh: RDFBasedHistFillingPP(24), trigger_mode = 3,
// every other flag default (isTight = true nominal; isTight = false -> Medium, as the crossx).
// The selection string is built from the SAME ParamsSet builders, in the same order, as
// RDFBasedHistFillingPP::FillHistogramsCrossx `signal_cuts`.
//
// Histograms (weights include weight * pp_crossx_lumi_factor except `counts`):
//   h_pair_eta_fine_<w>[_ptK]  w in {counts, raw, trig, trig_reco_rb, trig_reco_tb}; K = coarse pT bin
//   raw = no efficiency; trig = w_trig; trig_reco_rb = w_trig / eps^reco (reco-binned, NEW);
//   trig_reco_tb = w_trig * w_reco (w_reco of AddPairEfficiencyWeightColumns: the truth-binned eps(pT, eta, dR)
//                  when this macro was run (2026-10-11); since a2121af (2026-10-13) w_reco is the SF-corrected
//                  reco-binned eps -- a rerun would mislabel this series).
//
// Usage (from Analysis/plotting_codes/reco_effcy/reco_binned/):
//   root -l -b -q 'fill_pp24_data_pair_eta_smoothness.C+(true)'    // Tight
//   root -l -b -q 'fill_pp24_data_pair_eta_smoothness.C+(false)'   // Medium
// =================================================================================================
#include "../../../RDFBasedHistFilling/RDFBasedHistFillingPP.cxx"
#include "pair_eta_smoothness_common.h"

static PairRecoEffRecoBinnedEvaluator* s_pair_reco_eff_rb = nullptr;

class PairEtaSmoothnessPP : public RDFBasedHistFillingPP {
public:
    explicit PairEtaSmoothnessPP(const std::string& out) : RDFBasedHistFillingPP(24), my_out(out) {}

protected:
    std::string my_out;

    void SetIOPathsHook() override
    {
        RDFBasedHistFillingPP::SetIOPathsHook();
        output_file = my_out;   // NEVER the crossx hist file
        infile_var1D_json = "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/RDFBasedHistFilling/var1D_pp.json";
        gSystem->mkdir(gSystem->DirName(my_out.c_str()), kTRUE);
        std::cout << "[smoothness] input " << input_files.front() << " -> output " << output_file << std::endl;
    }

    void FillHistograms() override
    {
        if (trigger_effcy_calc) throw std::runtime_error("smoothness: must run in the crossx mode (trigger_mode 3)");
        if (pp_crossx_lumi_factor < 0.) throw std::runtime_error("smoothness: pp_crossx_lumi_factor not set");
        OpenEffcyPtFitFile();
        OpenPairEfficiencyInputs();
        ApplyMuonWorkingPointFilter();
        if (!s_pair_reco_eff_rb) {
            s_pair_reco_eff_rb = new PairRecoEffRecoBinnedEvaluator();
            s_pair_reco_eff_rb->Load(RecoBinnedProductFile(), isTight);
        }

        // == RDFBasedHistFillingPP::FillHistogramsCrossx `signal_cuts` (same builders, same order)
        const std::string signal_cuts =
            ParamsSet::SignalMinvCutExpr("minv") + " && " + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
            + ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && "
            + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta") + " && "
            + ParamsSet::PairFiducialEtaCutExpr("pair_eta");
        std::cout << "[smoothness] signal_cuts = " << signal_cuts << std::endl;

        const double lumi_factor = pp_crossx_lumi_factor;
        ROOT::RDF::RNode df = AddPairEfficiencyWeightColumns(
                map_at_checked(df_map, "df_op", "smoothness: df_op").Filter(signal_cuts, "crossx signal_cuts"))
            .Define("effcy_reco_rb", [](float pt, float eta, float m) {
                return s_pair_reco_eff_rb->Eval(pt, eta, m);
            }, {"pair_pt", "pair_eta", "minv"})
            .Define("w_reco_rb", "effcy_reco_rb > 0 ? 1.0 / effcy_reco_rb : 1.0")   // same convention as w_reco
            .Define("crossx_weight", [lumi_factor](double w) { return w * lumi_factor; }, {"weight"})
            .Define("cw_raw", "crossx_weight")
            .Define("cw_trig", "crossx_weight * w_trig")
            .Define("cw_trig_reco_rb", "crossx_weight * w_trig * w_reco_rb")
            .Define("cw_trig_reco_tb", "crossx_weight * w_trig * w_reco");

        const std::vector<double> fine = PairEtaFineSmoothnessBins();
        const std::vector<double> ept = PairRecoEffRecoBinnedEvaluator::PtEdges();
        const std::vector<std::pair<std::string, std::string>> weights = {
            {"counts", ""}, {"raw", "cw_raw"}, {"trig", "cw_trig"},
            {"trig_reco_rb", "cw_trig_reco_rb"}, {"trig_reco_tb", "cw_trig_reco_tb"}};
        auto book = [&](ROOT::RDF::RNode d, const std::string& tag) {
            for (const auto& w : weights) {
                const std::string n = "h_pair_eta_fine_" + w.first + tag;
                ROOT::RDF::TH1DModel m(n.c_str(), ";#eta^{pair};", (int)fine.size() - 1, fine.data());
                hist1d_rresultptr_map[n] = w.second.empty() ? d.Histo1D(m, "pair_eta") : d.Histo1D(m, "pair_eta", w.second);
            }
        };
        book(df, "");
        for (size_t k = 0; k + 1 < ept.size(); ++k) {
            // pT bins half-open [lo, hi) as TAxis; the last bin is closed at the top by the
            // evaluator's clamp only for eps -- here reco pT above the top edge is OUTSIDE every
            // panel (reported by the integrated-vs-panel-sum check in the plotter).
            const std::string cut = Form("pair_pt >= %.10g && pair_pt < %.10g", ept[k], ept[k + 1]);
            book(df.Filter(cut), Form("_pt%zu", k));
        }
        std::cout << "[smoothness] booked " << hist1d_rresultptr_map.size() << " histograms" << std::endl;
    }

    void HistPostProcessDataExtra() override
    {
        PrintPairEfficiencyStats();
        if (s_pair_reco_eff_rb) s_pair_reco_eff_rb->PrintStats();
    }

    void WriteOutputExtra() override
    {
        RDFBasedHistFillingPP::WriteOutputExtra();
        m_outfile->cd();
        TNamed("input", input_files.front().c_str()).Write();
        TNamed("reco_binned_product", RecoBinnedProductFile().c_str()).Write();
        TNamed("working_point", isTight ? "tight" : "medium").Write();
        TParameter<double>("pp_crossx_lumi_factor", pp_crossx_lumi_factor).Write();
        TNamed("definition",
               "pp24 2mu4 data, RDFBasedHistFillingPP crossx selection (df_op, WP filter, signal_cuts), weight*lumi; "
               "w_trig = AddPairEfficiencyWeightColumns; trig_reco_rb = /eps^reco(reco pT, eta^pair, m) coarse axes "
               "(PairRecoEffRecoBinnedEvaluator); trig_reco_tb = *w_reco (applied truth-binned eps(pT,eta,dR)); "
               "eta^pair in pair_eta_fine_smoothness (=ParamsSet::pair_eta_measured_bins); "
               "docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md 3c").Write();
    }
};

void fill_pp24_data_pair_eta_smoothness(bool tight_WP = true)
{
    const std::string out = SmoothnessDataFile(tight_WP);
    if (!gSystem->AccessPathName(out.c_str()))
        std::cout << "[smoothness] NOTE: overwriting this test's own previous output " << out << std::endl;
    PairEtaSmoothnessPP pp(out);
    pp.trigger_mode = 3;          // 2mu4 crossx, as run_crossx_hist_filling_pp24.sh
    pp.isTight = tight_WP;
    pp.Run();
    std::cout << "[smoothness] done -> " << out << std::endl;
}
