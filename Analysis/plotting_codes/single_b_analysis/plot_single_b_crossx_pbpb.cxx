#include "SingleBCrossxPlotterBase.cxx"

class SingleBCrossxPlotterPbPb : public SingleBCrossxPlotterBase {
private:
    int run_year;
    std::vector<std::string> ctr_candidates{"ctr0_5", "ctr5_10", "ctr10_20", "ctr20_30", "ctr30_50", "ctr50_80", "ctr50_100"};

    std::vector<std::string> DetectAvailableCtrBins() {
        std::vector<std::string> bins;
        for (const auto& ctr : ctr_candidates) {
            const std::string hname = "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr;
            TObject* obj = fin->Get(hname.c_str());
            if (obj && dynamic_cast<TH2D*>(obj)) bins.push_back(ctr);
        }
        if (bins.empty()) {
            throw std::runtime_error("No centrality-sliced crossx histograms found in input file.");
        }
        return bins;
    }

    static std::string CtrLabelFromSuffix(const std::string& ctr) {
        std::string s = ctr;
        if (s.rfind("ctr", 0) == 0) s = s.substr(3);
        std::replace(s.begin(), s.end(), '_', '-');
        return s;
    }

public:
    SingleBCrossxPlotterPbPb(const std::string& in_path, int run_year_in)
        : SingleBCrossxPlotterBase(in_path, "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis"), run_year(run_year_in) {}

    void Run() override {
        if (!Init()) return;

        const auto ctr_bins = DetectAvailableCtrBins();
        for (const auto& ctr : ctr_bins) {
            const std::string tag = "pbpb" + std::to_string(run_year) + "_" + ctr;

            Save2DColz("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr, tag + "_pair_pt_pair_eta_colz.png");
            Save2DColz("h2d_crossx_pair_pt_minv_w_signal_cuts_" + ctr, tag + "_pair_pt_minv_colz.png");
            Save2DColz("h2d_crossx_pair_pt_dr_w_signal_cuts_" + ctr, tag + "_pair_pt_dr_colz.png");

            DrawPairPtByEtaWithDrLines(
                "h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts_" + ctr,
                "Pb+Pb 20" + std::to_string(run_year) + ", medium WP, mu4_mu4noL1, " + CtrLabelFromSuffix(ctr) + "%",
                tag + "_pair_pt_in_eta_subplots_dr_lines.png");
        }
    }
};

void plot_single_b_crossx_pbpb(
    int run_year = 24,
    const std::string& input_file = "")
{
    std::string in_path = input_file;
    if (in_path.empty()) {
        const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + std::to_string(run_year)
            + "/histograms_real_pairs_pbpb_20" + std::to_string(run_year);
        const std::vector<std::string> candidates = {
            base + "_mu4_mu4noL1_coarse_q_eta_bin.root",
            base + "_mu4_mu4noL1_mindR_noCut.root",
            base + "_mu4_mu4noL1.root",
            base + "_single_mu4.root",
            base + "_single_mu4_coarse_q_eta_bin.root"
        };

        for (const auto& p : candidates) {
            if (!gSystem->AccessPathName(p.c_str())) {
                in_path = p;
                break;
            }
        }
        if (in_path.empty()) {
            in_path = candidates.front();
        }
    }

    SingleBCrossxPlotterPbPb pl(in_path, run_year);
    pl.Run();
}
