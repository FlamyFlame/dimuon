#include "SingleBCrossxPlotterBase.cxx"

// ======================== PbPb combined-year plotter ========================
// Sums histograms from all available PbPb years. This is the only supported
// mode — per-year plots are not produced.

class SingleBCrossxPlotterPbPbCombined : public SingleBCrossxPlotterBase {
    std::vector<std::pair<int,std::string>> year_paths_;  // (2-digit year, file path)
    std::vector<TFile*>  files_;
    std::map<std::string,TH1*> hist_cache_;  // owned combined histograms

    std::string label_line1_;  // e.g. "Pb+Pb 2023, 2024 combined"
    std::string label_line2_;  // "medium WP"
    std::string label_line3_;  // e.g. "mu4 (2023), mu4 (2024)"

    std::vector<std::string> ctr_candidates_{"ctr0_5","ctr5_10","ctr10_20","ctr20_30","ctr30_50","ctr50_80","ctr50_100"};

    void BuildLabels() {
        std::string yrs_str;
        label_line3_ = "";
        for (size_t i = 0; i < year_paths_.size(); ++i) {
            int yr = year_paths_.at(i).first;
            if (i > 0) { yrs_str  += ", "; label_line3_ += ", "; }
            yrs_str  += "20" + std::to_string(yr);
            label_line3_ += DatasetTriggerMap::GetTriggerLabel(yr, "PbPb")
                            + " (20" + std::to_string(yr) + ")";
        }
        label_line1_ = "Pb+Pb " + yrs_str + " combined";
        label_line2_ = "medium WP";
    }

    std::string OutDirName() const {
        std::string s = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pbpb";
        for (auto& [yr, path] : year_paths_) s += "_" + std::to_string(yr);
        s += "_combined";
        return s;
    }

    bool InitCombined() {
        BuildLabels();
        output_dir = OutDirName();
        gSystem->mkdir(output_dir.c_str(), true);
        for (auto& [yr, path] : year_paths_) {
            TFile* f = TFile::Open(path.c_str(), "READ");
            if (!f || f->IsZombie())
                throw std::runtime_error("Cannot open: " + path);
            files_.push_back(f);
        }
        gStyle->SetOptStat(0);
        return true;
    }

    std::vector<std::string> DetectAvailableCtrBins() {
        if (files_.empty()) return {};
        std::vector<std::string> bins;
        for (const auto& ctr : ctr_candidates_) {
            const std::string hname = "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr;
            if (files_.front()->Get(hname.c_str())) bins.push_back(ctr);
        }
        if (bins.empty()) throw std::runtime_error("No centrality crossx histograms found.");
        return bins;
    }

    bool HasCountsHists() {
        if (files_.empty()) return false;
        const std::string test = "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_"
                                 + ctr_candidates_.front() + "_counts";
        return files_.front()->Get(test.c_str()) != nullptr;
    }

    static std::string CtrLabelFromSuffix(const std::string& ctr) {
        std::string s = ctr;
        if (s.rfind("ctr", 0) == 0) s = s.substr(3);
        std::replace(s.begin(), s.end(), '_', '-');
        return s;
    }

protected:
    TObject* GetHistObject(const std::string& name) override {
        auto it = hist_cache_.find(name);
        if (it != hist_cache_.end()) return it->second;

        TH1* combined = nullptr;
        for (TFile* f : files_) {
            TH1* h = dynamic_cast<TH1*>(f->Get(name.c_str()));
            if (!h) continue;
            if (!combined) {
                combined = dynamic_cast<TH1*>(h->Clone(name.c_str()));
                combined->SetDirectory(nullptr);
            } else {
                combined->Add(h);
            }
        }
        if (!combined) return nullptr;
        hist_cache_[name] = combined;
        return combined;
    }

public:
    bool use_pt_bins_150 = false;

    SingleBCrossxPlotterPbPbCombined(const std::vector<std::pair<int,std::string>>& year_paths)
        : SingleBCrossxPlotterBase("", ""), year_paths_(year_paths) {}

    ~SingleBCrossxPlotterPbPbCombined() override {
        for (auto& [name, h] : hist_cache_) delete h;
        for (TFile* f : files_) { f->Close(); delete f; }
    }

    void Run() override {
        InitCombined();

        const auto ctr_bins   = DetectAvailableCtrBins();
        const bool has_counts = HasCountsHists();

        const std::string base_out   = output_dir;
        const std::string counts_dir = base_out + "/counts";
        const std::string taa_dir    = base_out + "/TAA_weighted";
        gSystem->mkdir(taa_dir.c_str(), true);
        if (has_counts) gSystem->mkdir(counts_dir.c_str(), true);

        for (const auto& ctr : ctr_bins) {
            const std::string tag     = "pbpb_combined_" + ctr;
            const std::string ctr_pct = CtrLabelFromSuffix(ctr) + "%";
            const std::string l1 = label_line1_ + ", " + ctr_pct;

            if (has_counts) {
                output_dir = counts_dir;
                Save2DColz("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr + "_counts",
                           tag + "_pair_pt_pair_eta.png",
                           "d^{2}N_{events}/dp_{T}d#eta [GeV^{-1}]");
                Save2DColz("h2d_crossx_pair_pt_minv_w_signal_cuts_" + ctr + "_counts",
                           tag + "_pair_pt_minv.png",
                           "d^{2}N_{events}/dp_{T}dm_{#mu#mu} [GeV^{-1} GeV^{-1}]");
                Save2DColz("h2d_crossx_pair_pt_dr_w_signal_cuts_" + ctr + "_counts",
                           tag + "_pair_pt_dr.png",
                           "d^{2}N_{events}/dp_{T}d#DeltaR [GeV^{-1}]");
                DrawPairPtByEtaWithDrLines(
                    "h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts_" + ctr + "_counts",
                    l1, label_line3_,
                    tag + "_pair_pt_in_eta_subplots_dr_lines.png",
                    "dN_{events}/dp_{T} [GeV^{-1}]");
                DrawPairPtByEta(
                    "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr + "_counts",
                    l1, label_line3_,
                    tag + "_pair_pt_in_eta_subplots.png",
                    "dN_{events}/dp_{T} [GeV^{-1}]");
                DrawPairPtByEta(
                    "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr + "_counts",
                    l1, label_line3_,
                    tag + "_pair_pt_in_eta_subplots_nondifferential.png",
                    "N_{events}", false);
            }

            output_dir = taa_dir;
            Save2DColz("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr,
                       tag + "_pair_pt_pair_eta.png",
                       "d^{2}#sigma/dp_{T}d#eta [pb GeV^{-1}]");
            Save2DColz("h2d_crossx_pair_pt_minv_w_signal_cuts_" + ctr,
                       tag + "_pair_pt_minv.png",
                       "d^{2}#sigma/dp_{T}dm_{#mu#mu} [pb GeV^{-1} GeV^{-1}]");
            Save2DColz("h2d_crossx_pair_pt_dr_w_signal_cuts_" + ctr,
                       tag + "_pair_pt_dr.png",
                       "d^{2}#sigma/dp_{T}d#DeltaR [pb GeV^{-1}]");
            DrawPairPtByEtaWithDrLines(
                "h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts_" + ctr,
                l1, label_line3_,
                tag + "_pair_pt_in_eta_subplots_dr_lines.png",
                "d#sigma/dp_{T} [pb GeV^{-1}]");
            DrawPairPtByEta(
                "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr,
                l1, label_line3_,
                tag + "_pair_pt_in_eta_subplots.png",
                "d#sigma/dp_{T} [pb GeV^{-1}]");
        }

        if (use_pt_bins_150) {
            const std::string pt150_dir    = base_out + "_pt_150";
            const std::string pt150_cntdir = pt150_dir + "/counts";
            const std::string pt150_taadir = pt150_dir + "/TAA_weighted";
            gSystem->mkdir(pt150_taadir.c_str(), true);
            if (has_counts) gSystem->mkdir(pt150_cntdir.c_str(), true);

            for (const auto& ctr : ctr_bins) {
                const std::string tag     = "pbpb_combined_" + ctr;
                const std::string ctr_pct = CtrLabelFromSuffix(ctr) + "%";
                const std::string l1 = label_line1_ + ", " + ctr_pct;

                if (has_counts) {
                    output_dir = pt150_cntdir;
                    DrawPairPtByEtaWithDrLines(
                        "h3d_crossx_dr_vs_pair_eta_vs_pt_150_w_signal_cuts_" + ctr + "_counts",
                        l1, label_line3_,
                        tag + "_pair_pt_in_eta_subplots_dr_lines.png",
                        "dN_{events}/dp_{T} [GeV^{-1}]");
                    DrawPairPtByEta(
                        "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pt_150_" + ctr + "_counts",
                        l1, label_line3_,
                        tag + "_pair_pt_in_eta_subplots.png",
                        "dN_{events}/dp_{T} [GeV^{-1}]");
                }
                output_dir = pt150_taadir;
                DrawPairPtByEtaWithDrLines(
                    "h3d_crossx_dr_vs_pair_eta_vs_pt_150_w_signal_cuts_" + ctr,
                    l1, label_line3_,
                    tag + "_pair_pt_in_eta_subplots_dr_lines.png",
                    "d#sigma/dp_{T} [pb GeV^{-1}]");
                DrawPairPtByEta(
                    "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pt_150_" + ctr,
                    l1, label_line3_,
                    tag + "_pair_pt_in_eta_subplots.png",
                    "d#sigma/dp_{T} [pb GeV^{-1}]");
            }
            output_dir = base_out;
        }
    }
};

void plot_single_b_crossx_pbpb(bool use_pt_bins_150 = false)
{
    static const std::string plots_base =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/";

    std::vector<std::pair<int,std::string>> year_paths;
    for (int yr : {23, 24, 25}) {
        const std::string& trig = DatasetTriggerMap::GetTrigger(yr, "PbPb");
        const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20"
                                 + std::to_string(yr) + "/histograms_real_pairs_pbpb_20"
                                 + std::to_string(yr);
        std::string path = base + "_" + trig + "_coarse_q_eta_bin.root";
        if (gSystem->AccessPathName(path.c_str())) {
            path = base + "_" + trig + ".root";
            if (gSystem->AccessPathName(path.c_str())) {
                std::cout << "[INFO] No input file found for PbPb 20" << yr << " — skipping." << std::endl;
                continue;
            }
        }
        std::cout << "[INFO] Found PbPb 20" << yr << ": " << path << std::endl;
        year_paths.push_back({yr, path});
    }
    if (year_paths.empty())
        throw std::runtime_error("plot_single_b_crossx_pbpb: no PbPb input files found.");

    SingleBCrossxPlotterPbPbCombined pl(year_paths);
    pl.use_pt_bins_150 = use_pt_bins_150;
    pl.Run();
}
