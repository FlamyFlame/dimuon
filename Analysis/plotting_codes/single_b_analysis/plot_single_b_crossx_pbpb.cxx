#include "SingleBCrossxPlotterBase.cxx"
#include "../../Utilities/PbPbSampledLumi.h"

// ======================== PbPb combined-year plotter ========================
// Sums histograms from all available PbPb years. This is the only supported
// mode — per-year plots are not produced.

class SingleBCrossxPlotterPbPbCombined : public SingleBCrossxPlotterBase {
    std::vector<std::pair<int,std::string>> year_paths_;  // (2-digit year, file path)
    std::vector<TFile*>  files_;
    std::map<std::string,TH1*> hist_cache_;  // owned combined histograms
    // Which YEARS actually contributed to each combined histogram. A year whose histogram is
    // absent is skipped by GetHistObject, and before 2026-09-09 that happened in silence while
    // the canvas label still claimed all three years -- a mislabelled final-results figure.
    std::map<std::string,std::vector<int>> combined_years_;

    std::string label_line1_;  // e.g. "Pb+Pb 2023, 2024 combined"
    std::string label_line2_;  // "tight WP"
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
        label_line2_ = "tight WP";
    }

    // The drawn label must describe the years that ACTUALLY contributed to THIS histogram, not the
    // years the job was asked to combine. GetCombined() already records them per histogram name and
    // warns when a year is missing -- but the label was still built once from every entry in
    // year_paths_, so the warning said "the figure label must not claim the missing year(s)" and
    // then the figure claimed them anyway. This reads the recorded list back.
    std::string YearsLabelFor(const std::string& hname) const {
        auto it = combined_years_.find(hname);
        if (it == combined_years_.end() || it->second.size() == year_paths_.size()) return label_line1_;
        std::string yrs;
        for (size_t i = 0; i < it->second.size(); ++i) {
            if (i > 0) yrs += ", ";
            yrs += "20" + std::to_string(it->second.at(i));
        }
        return "Pb+Pb " + yrs + " combined";
    }

    // Fetch-then-label. The fetch is what POPULATES combined_years_, so it has to happen before the
    // label is formed; doing it here rather than relying on argument evaluation order keeps the two
    // correct whichever order the compiler picks (the result is cached, so this costs nothing).
    std::string L1For(const std::string& hname, const std::string& ctr_pct) {
        (void)GetHistObject(hname);
        return YearsLabelFor(hname) + ", " + ctr_pct;
    }

    std::string OutDirName() const {
        std::string s = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pbpb";
        for (auto& [yr, path] : year_paths_) s += "_" + std::to_string(yr);
        s += "_combined";
        // "" for the nominal 9 -> 150 GeV view, "_pt_120" for the opt-in alternative.
        s += PtAxisDirSuffix();
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
            // Probe the family actually being drawn, not the other one: a centrality bin whose
            // selected-axis histogram is absent must not be reported as available.
            const std::string hname =
                PtAxisHist("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr);
            // Probe EVERY year, not just files_.front(). Probing only the first year silently
            // dropped a centrality bin that the first file happens to lack but the others have.
            bool any = false, all = true;
            for (auto* f : files_) { if (f->Get(hname.c_str())) any = true; else all = false; }
            if (any) {
                bins.push_back(ctr);
                if (!all)
                    std::cerr << "[WARN] centrality " << ctr << ": histogram " << hname
                              << " is missing from at least one year -- the combined figure for"
                              << " this bin averages only the years that have it." << std::endl;
            }
        }
        if (bins.empty()) throw std::runtime_error("No centrality crossx histograms found.");
        return bins;
    }

    bool HasCountsHists() {
        if (files_.empty()) return false;
        const std::string test = PtAxisHist("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_"
                                            + ctr_candidates_.front() + "_counts");
        // Probe EVERY year: keyed on files_.front() alone, one year lacking the counts family
        // silently suppressed the ENTIRE counts figure set with no message.
        for (auto* f : files_) if (f->Get(test.c_str())) return true;
        return false;
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

        // Year combination (HF R_AA note HION-2019-58 §4.1 Eq.3):
        //  - cross-section histos: LUMINOSITY-WEIGHTED AVERAGE  Sum(L_y h_y)/Sum(L_y).
        //    Each h_y is already crossx_factor(~1/L_y)-weighted, so this yields the
        //    correct combined ΣN/(f·σ·T_AA·ΣL), NOT a naive sum (which over-counts ~×years).
        //  - "_counts" histos (raw event counts, weight=1): SIMPLE SUM (total counts).
        const bool is_counts = (name.find("_counts") != std::string::npos);
        TH1* combined = nullptr;
        double sumL = 0.;
        std::vector<int> contributing_years;
        std::vector<int> missing_years;
        for (size_t i = 0; i < files_.size(); ++i) {
            TH1* h = dynamic_cast<TH1*>(files_[i]->Get(name.c_str()));
            if (!h) { missing_years.push_back(year_paths_[i].first); continue; }
            contributing_years.push_back(year_paths_[i].first);
            const double w = is_counts ? 1.0 : PbPbMu4SampledLumiNb(year_paths_[i].first);
            if (!combined) {
                combined = dynamic_cast<TH1*>(h->Clone(name.c_str()));
                combined->SetDirectory(nullptr);
                combined->Scale(w);
            } else {
                combined->Add(h, w);
            }
            sumL += w;
        }
        if (!combined) return nullptr;
        // A year missing this histogram used to be dropped in SILENCE while the canvas label
        // still read "Pb+Pb 2023, 2024, 2025 combined" -- a mislabelled final-results figure,
        // with sumL normalising over the subset so the number looked entirely reasonable.
        if (!missing_years.empty()) {
            std::cerr << "[WARN] " << name << ": missing from year(s)";
            for (int y : missing_years) std::cerr << " 20" << y;
            std::cerr << "; the combined result uses only";
            for (int y : contributing_years) std::cerr << " 20" << y;
            std::cerr << ". The figure label must not claim the missing year(s)." << std::endl;
        }
        combined_years_[name] = contributing_years;
        if (!is_counts && sumL > 0.) combined->Scale(1.0 / sumL);
        hist_cache_[name] = combined;
        return combined;
    }

public:
    SingleBCrossxPlotterPbPbCombined(const std::vector<std::pair<int,std::string>>& year_paths)
        : SingleBCrossxPlotterBase("", ""), year_paths_(year_paths) {}

    ~SingleBCrossxPlotterPbPbCombined() override {
        for (auto& [name, h] : hist_cache_) delete h;
        for (TFile* f : files_) { f->Close(); delete f; }
    }

    // One invocation draws exactly ONE pair-pT axis into its own directory tree (both the
    // TAA_weighted and the counts variant): the nominal 9 -> 150 GeV view into
    // pbpb_<years>_combined/, the opt-in 9 -> 120 GeV view into pbpb_<years>_combined_pt_120/.
    // Every histogram is named ONCE in its canonical unsuffixed form and mapped onto the selected
    // family by PtAxisHist(), so no figure in a directory can end up on the other axis.
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

            if (has_counts) {
                output_dir = counts_dir;
                Save2DColz(PtAxisHist("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr + "_counts"),
                           tag + "_pair_pt_pair_eta.png",
                           "d^{2}N_{events}/dp_{T}d#eta [GeV^{-1}]");
                Save2DColz(PtAxisHist("h2d_crossx_pair_pt_minv_w_signal_cuts_" + ctr + "_counts"),
                           tag + "_pair_pt_minv.png",
                           "d^{2}N_{events}/dp_{T}dm_{#mu#mu} [GeV^{-1} GeV^{-1}]");
                Save2DColz(PtAxisHist("h2d_crossx_pair_pt_dr_w_signal_cuts_" + ctr + "_counts"),
                           tag + "_pair_pt_dr.png",
                           "d^{2}N_{events}/dp_{T}d#DeltaR [GeV^{-1}]");
                {
                    const std::string hn = PtAxisHist("h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts_" + ctr + "_counts");
                    DrawPairPtByEtaWithDrLines(
                        hn, L1For(hn, ctr_pct), label_line3_,
                        tag + "_pair_pt_in_eta_subplots_dr_lines.png",
                        "dN_{events}/dp_{T} [GeV^{-1}]");
                }
                {
                    const std::string hn = PtAxisHist("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr + "_counts");
                    DrawPairPtByEta(
                        hn, L1For(hn, ctr_pct), label_line3_,
                        tag + "_pair_pt_in_eta_subplots.png",
                        "dN_{events}/dp_{T} [GeV^{-1}]");
                }
                {
                    const std::string hn = PtAxisHist("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr + "_counts");
                    DrawPairPtByEta(
                        hn, L1For(hn, ctr_pct), label_line3_,
                        tag + "_pair_pt_in_eta_subplots_nondifferential.png",
                        "N_{events}", false);
                }
            }

            output_dir = taa_dir;
            Save2DColz(PtAxisHist("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr),
                       tag + "_pair_pt_pair_eta.png",
                       "#frac{1}{#LTT_{AA}#GT N_{evt}} #frac{d^{2}n_{AA}}{dp_{T}d#eta} [pb GeV^{-1}]");
            Save2DColz(PtAxisHist("h2d_crossx_pair_pt_minv_w_signal_cuts_" + ctr),
                       tag + "_pair_pt_minv.png",
                       "#frac{1}{#LTT_{AA}#GT N_{evt}} #frac{d^{2}n_{AA}}{dp_{T}dm_{#mu#mu}} [pb GeV^{-1} GeV^{-1}]");
            Save2DColz(PtAxisHist("h2d_crossx_pair_pt_dr_w_signal_cuts_" + ctr),
                       tag + "_pair_pt_dr.png",
                       "#frac{1}{#LTT_{AA}#GT N_{evt}} #frac{d^{2}n_{AA}}{dp_{T}d#DeltaR} [pb GeV^{-1}]");
            {
                const std::string hn = PtAxisHist("h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts_" + ctr);
                DrawPairPtByEtaWithDrLines(
                    hn, L1For(hn, ctr_pct), label_line3_,
                    tag + "_pair_pt_in_eta_subplots_dr_lines.png",
                    "#frac{1}{#LTT_{AA}#GT N_{evt}} #frac{dn_{AA}}{dp_{T}} [pb GeV^{-1}]");
            }
            {
                const std::string hn = PtAxisHist("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr);
                DrawPairPtByEta(
                    hn, L1For(hn, ctr_pct), label_line3_,
                    tag + "_pair_pt_in_eta_subplots.png",
                    "#frac{1}{#LTT_{AA}#GT N_{evt}} #frac{dn_{AA}}{dp_{T}} [pb GeV^{-1}]");
            }
        }

        output_dir = base_out;
    }
};

// `also_pt_120 = true` additionally refreshes the OPT-IN 9 -> 120 GeV alternative view in
// pbpb_<years>_combined_pt_120/, in the SAME invocation as the nominal 9 -> 150 GeV one --
// refreshing only one of the two is how they drifted apart for months (2026-06-19 to 2026-08-04,
// back when the 150 view was the opt-in one).
void plot_single_b_crossx_pbpb(bool also_pt_120 = false)
{
    static const std::string plots_base =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/";

    std::vector<std::pair<int,std::string>> year_paths;
    for (int yr : {23, 24, 25}) {
        const std::string& trig = DatasetTriggerMap::GetTrigger(yr, "PbPb");
        const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20"
                                 + std::to_string(yr) + "/histograms_real_pairs_pbpb_20"
                                 + std::to_string(yr);
        std::vector<std::string> candidates = {
            base + "_" + trig + "_no_trg_plots_nominal.root",
            base + "_" + trig + "_no_trg_plots_coarse_q_eta_bin.root",
            base + "_" + trig + "_coarse_q_eta_bin.root",
            base + "_" + trig + "_no_trg_plots_fine_q_eta_bin.root",
            base + "_" + trig + "_fine_q_eta_bin.root",
        };
        std::string path;
        for (const auto& c : candidates) {
            if (!gSystem->AccessPathName(c.c_str())) { path = c; break; }
        }
        if (path.empty()) {
            std::cout << "[INFO] No input file found for PbPb 20" << yr << " — skipping." << std::endl;
            continue;
        }
        std::cout << "[INFO] Found PbPb 20" << yr << ": " << path << std::endl;
        year_paths.push_back({yr, path});
    }
    if (year_paths.empty())
        throw std::runtime_error("plot_single_b_crossx_pbpb: no PbPb input files found.");

    SingleBCrossxPlotterPbPbCombined pl(year_paths);  // DEFAULT: pT_bins_150 (9 -> 150 GeV)
    pl.Run();

    if (also_pt_120) {
        SingleBCrossxPlotterPbPbCombined pl120(year_paths);
        pl120.use_pt_bins_120 = true;                 // OPT-IN: pT_bins_120 (9 -> 120 GeV)
        pl120.Run();
    }
}
