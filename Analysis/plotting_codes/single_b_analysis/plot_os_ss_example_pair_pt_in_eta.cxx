// Illustration of what the opposite-sign, same-sign and (opposite - same)-sign pair spectra look
// like in three representative (sub-)datasets: pp 2024, Pb+Pb 2023-26 combined 0-10 % and
// Pb+Pb 2023-26 combined 50-80 %. One pair_pt_in_eta_subplots figure per dataset, three lines per
// eta^pair panel:  (1) opposite sign   (2) same sign   (3) opposite minus same sign.
//
// Both signs go through the IDENTICAL selection (every reco-level signal cut, Tight WP, trigger,
// centrality) and the IDENTICAL weight chain; only the pair sign differs. The histograms are the
// ones the producers already book for the R_AA subtraction / the OS-SS crossx pass:
//   pp 2024  : h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts  / h2d_ss_<same>   (RDFBasedHistFillingPP)
//   Pb+Pb    : h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_<ctr> / h2d_ss_<same>
//              (RDFBasedHistFillingPbPb ss_subtr_crossx pass, *_ss_subtr_crossx.root)
// Axes: pair pT = pT_bins_150 (as booked), eta^pair panels = pair_eta_proj_ranges_coarse_incl_gap
// (read from CommonEffcyConfig by the base class). Nothing retyped, nothing rebinned.
//
// Pb+Pb 0-10 % is NOT a filled centrality bin: it is built here from the filled 0-5 % and 5-10 %.
// Each filled bin i carries Y_i = n_i * factor_i with factor_i = 1/(f_i sigma_PbPb <T_AA>_i L)
// (PbPbBaseClass::make_crossx_factors_pbpb_<yr>), i.e. Y_i = n_i / D_i with D_i the bin's
// "T_AA-weighted event normalisation". The merged bin has D = sum_i D_i  (<T_AA>_{0-10} =
// sum_i f_i <T_AA>_i / sum_i f_i), so  Y_{0-10} = sum_i Y_i D_i / sum_i D_i,  D_i = 1/factor_i.
// This is exact (the histograms are linear in the weights) and year-independent in its ratio.
// The factors are READ from PbPbBaseClass, not retyped.
//
// Year combination: luminosity-weighted average sum_y L_y h_y / sum_y L_y, exactly as
// SingleBCrossxPlotterPbPbCombined::GetCombined does for the weighted histograms.
//
// Output: .../plots/single_b_analysis/os_ss_examples/<dataset>_pair_pt_in_eta_subplots_os_ss.png
// Run from the old-GRL run-only worktree for the frozen-input Pb+Pb 2025/26 deliverable (its
// PbPbSampledLumi.h carries the old 2026 lumi that matches those histograms).

#include "SingleBCrossxPlotterBase.cxx"
#include "../../Utilities/PbPbSampledLumi.h"
#include "../../MuonObjectsParamsAndHelpers/PbPbBaseClass.h"

namespace OsSsExample {
// Dummy Derived: PbPbBaseClass is a CRTP template; only its static crossx-factor tables are used.
struct FactorsProbe {};
using PbPbFactors = PbPbBaseClass<FactorsProbe>;

inline std::vector<double> CrossxFactors(int yr) {
    switch (yr) {
        case 23: return PbPbFactors::make_crossx_factors_pbpb_2023();
        case 24: return PbPbFactors::make_crossx_factors_pbpb_2024();
        case 25: return PbPbFactors::make_crossx_factors_pbpb_2025();
        case 26: return PbPbFactors::make_crossx_factors_pbpb_2026();
    }
    throw std::runtime_error("OsSsExample::CrossxFactors: no Pb+Pb crossx factors for year "
                             + std::to_string(yr));
}

// Index of a filled centrality bin in the crossx-factor vector (same order as ctr_bin_edges
// {0,5,10,20,30,50,80}).
inline int CtrFactorIndex(const std::string& ctr) {
    static const std::vector<std::string> order{"ctr0_5","ctr5_10","ctr10_20","ctr20_30","ctr30_50","ctr50_80"};
    for (size_t i = 0; i < order.size(); ++i) if (order[i] == ctr) return static_cast<int>(i);
    throw std::runtime_error("OsSsExample::CtrFactorIndex: unknown centrality bin " + ctr);
}

struct OsSsPair { std::unique_ptr<TH2D> os, ss; };

inline TH2D* GetTH2Checked(TFile* f, const std::string& name) {
    TH2D* h = dynamic_cast<TH2D*>(f->Get(name.c_str()));
    if (!h) throw std::runtime_error("OsSsExample: missing TH2D '" + name + "' in " + f->GetName());
    return h;
}

// Pb+Pb: opposite- and same-sign (pair pT x eta^pair) TAA-weighted yields for the union of the
// given filled centrality bins, lumi-combined over `years`.
inline OsSsPair LoadPbPb(const std::vector<int>& years, const std::vector<std::string>& ctrs) {
    OsSsPair out;
    double sumL = 0.;
    for (int yr : years) {
        const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + std::to_string(yr)
            + "/histograms_real_pairs_pbpb_20" + std::to_string(yr) + "_"
            + DatasetTriggerMap::GetTrigger(yr, "PbPb") + "_no_trg_plots_nominal_ss_subtr_crossx.root";
        std::unique_ptr<TFile> f(TFile::Open(base.c_str(), "READ"));
        if (!f || f->IsZombie()) throw std::runtime_error("OsSsExample: cannot open " + base);
        const std::vector<double> fac = CrossxFactors(yr);
        double sumD = 0.;
        for (const auto& c : ctrs) sumD += 1. / fac.at(CtrFactorIndex(c));
        const double lumi = PbPbMu4SampledLumiNb(yr);
        for (const auto& c : ctrs) {
            const double wD = (1. / fac.at(CtrFactorIndex(c))) / sumD;   // bin's share of D
            const std::string nm_os = "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + c;
            const std::string nm_ss = "h2d_ss_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + c;
            TH2D* hos = GetTH2Checked(f.get(), nm_os);
            TH2D* hss = GetTH2Checked(f.get(), nm_ss);
            if (!out.os) {
                out.os.reset(dynamic_cast<TH2D*>(hos->Clone(("os_" + std::to_string(std::rand())).c_str())));
                out.ss.reset(dynamic_cast<TH2D*>(hss->Clone(("ss_" + std::to_string(std::rand())).c_str())));
                out.os->SetDirectory(nullptr); out.ss->SetDirectory(nullptr);
                out.os->Reset(); out.ss->Reset();   // keeps Sumw2 structure; zeroes content
            }
            out.os->Add(hos, wD * lumi);
            out.ss->Add(hss, wD * lumi);
        }
        sumL += lumi;
    }
    if (!out.os || sumL <= 0.) throw std::runtime_error("OsSsExample::LoadPbPb: nothing loaded");
    out.os->Scale(1. / sumL);
    out.ss->Scale(1. / sumL);
    return out;
}

inline OsSsPair LoadPP24() {
    const std::string path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_"
        + DatasetTriggerMap::GetTrigger(24, "pp") + "_nominal.root";
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) throw std::runtime_error("OsSsExample: cannot open " + path);
    OsSsPair out;
    out.os.reset(dynamic_cast<TH2D*>(GetTH2Checked(f.get(), "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts")->Clone("pp_os")));
    out.ss.reset(dynamic_cast<TH2D*>(GetTH2Checked(f.get(), "h2d_ss_crossx_pair_pt_pair_eta_binned_w_signal_cuts")->Clone("pp_ss")));
    out.os->SetDirectory(nullptr); out.ss->SetDirectory(nullptr);
    return out;
}
}  // namespace OsSsExample

class OsSsExamplePlotter : public SingleBCrossxPlotterBase {
public:
    OsSsExamplePlotter()
        : SingleBCrossxPlotterBase("",
              "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/os_ss_examples")
    {
        subplot_ceil_pad = 30.0;   // header + 3-entry legend sit above the highest points
        gSystem->mkdir(output_dir.c_str(), true);
        gStyle->SetOptStat(0);
    }

    void Run() override {}   // unused: figures are drawn through Draw() from the macro entry point

    // One figure: eta^pair panels, three lines each (opposite sign, same sign, opposite - same sign).
    void Draw(const OsSsExample::OsSsPair& p, const std::string& line1, const std::string& line2,
              const std::string& png_name, const std::string& y_title)
    {
        TH2D* os2 = p.os.get();
        TH2D* ss2 = p.ss.get();
        std::unique_ptr<TH2D> df2(dynamic_cast<TH2D*>(os2->Clone("os_minus_ss_2d")));
        df2->SetDirectory(nullptr);
        df2->Add(ss2, -1.0);   // Sumw2 present on both -> errors add in quadrature

        int nrow = 1, ncol = 1;
        DetermineSubplotGrid(static_cast<int>(q_eta_bins.size()), nrow, ncol);
        TCanvas c("c_os_ss", "os ss pair_pt by eta", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        // line order = drawing order = legend order
        const std::array<int, 3> col{{kBlack, kRed + 1, kBlue + 1}};
        const std::array<int, 3> mst{{20, 25, 22}};
        const std::array<const char*, 3> lab{{"opposite sign", "same sign", "opposite #minus same sign"}};
        std::vector<std::array<TH1D*, 3>> lines(q_eta_bins.size());
        std::vector<TH1*> flat;
        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            const auto& eb = q_eta_bins.at(ieta);
            const TH2D* src[3] = {os2, ss2, df2.get()};
            for (int il = 0; il < 3; ++il) {
                const auto yb = PanelEtaBins(src[il]->GetYaxis(), eb, "OsSsExamplePlotter (pT-in-eta panels)");
                TH1D* hp = const_cast<TH2D*>(src[il])->ProjectionX(
                    ("hp" + std::to_string(ieta) + "_" + std::to_string(il) + "_" + std::to_string(std::rand())).c_str(),
                    yb.first, yb.second, "e");
                hp->SetDirectory(nullptr);
                hp->Scale(1.0, "width");   // differential: sum(w)/bin width
                hp->SetLineWidth(2);
                hp->SetLineColor(col[il]); hp->SetMarkerColor(col[il]);
                hp->SetMarkerStyle(mst[il]); hp->SetMarkerSize(0.9);
                hp->SetTitle("");
                hp->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
                hp->GetYaxis()->SetTitle(y_title.c_str());
                hp->GetXaxis()->SetTitleSize(0.06); hp->GetYaxis()->SetTitleSize(0.06);
                hp->GetXaxis()->SetLabelSize(0.05); hp->GetYaxis()->SetLabelSize(0.05);
                hp->GetYaxis()->SetTitleOffset(1.45);
                lines[ieta][il] = hp;
                flat.push_back(hp);
            }
            // Cells where (opposite - same) <= 0 cannot be drawn on the log axis: report them.
            int n_nonpos = 0, n_os = 0;
            for (int b = 1; b <= lines[ieta][2]->GetNbinsX(); ++b) {
                if (lines[ieta][0]->GetBinContent(b) <= 0.) continue;
                ++n_os;
                if (lines[ieta][2]->GetBinContent(b) <= 0.) ++n_nonpos;
            }
            std::cout << "[INFO] " << png_name << " eta^pair [" << eb.first << "," << eb.second << "]: "
                      << n_nonpos << " of " << n_os << " pT bins with opposite sign filled have"
                      << " (opposite - same) <= 0 and are not drawn." << std::endl;
        }
        // ONE log-y range for the whole figure, from the positive points of all three lines.
        ApplyCommonLogYRange(flat, subplot_ceil_pad);
        // One more decade below the lowest point: the 3-entry key sits in the bottom-left corner and
        // must not cover the (same-sign) points. Only adds empty space, never hides a point.
        for (TH1* h : flat) h->SetMinimum(h->GetMinimum() * 0.1);

        std::vector<TLegend*> legs;
        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            c.cd(static_cast<int>(ieta) + 1);
            gPad->SetLogx(); gPad->SetLogy();
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.13);
            lines[ieta][0]->Draw("E1");
            lines[ieta][1]->Draw("E1 SAME");
            lines[ieta][2]->Draw("E1 SAME");

            const auto& eb = q_eta_bins.at(ieta);
            TLegend* info = new TLegend(0.38, 0.70, 1.0 - gPad->GetRightMargin(), 0.90);
            info->SetBorderSize(0); info->SetFillStyle(0);
            info->SetTextSize(0.045); info->SetTextAlign(32); info->SetMargin(0.01);
            info->AddEntry((TObject*)0, line1.c_str(), "");
            info->AddEntry((TObject*)0, line2.c_str(), "");
            info->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eb.first, eb.second), "");
            info->Draw();
            TLegend* key = new TLegend(0.17, 0.14, 0.60, 0.33);   // bottom-left: empty for a falling spectrum
            key->SetBorderSize(0); key->SetFillColor(kWhite); key->SetFillStyle(1001);  // opaque: lower error bars of low-yield points run to the axis bottom and would cut through the text
            key->SetTextSize(0.045); key->SetTextAlign(12); key->SetMargin(0.12);
            for (int il = 0; il < 3; ++il) key->AddEntry(lines[ieta][il], lab[il], "lep");
            key->Draw();
            legs.push_back(info); legs.push_back(key);
        }
        const std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;
        for (auto& a : lines) for (TH1D* h : a) delete h;
        for (TLegend* l : legs) delete l;
    }
};

void plot_os_ss_example_pair_pt_in_eta()
{
    OsSsExamplePlotter pl;
    const std::string pbpb_ytitle = "#frac{1}{#LTT_{AA}#GT N_{evt}} #frac{dn_{AA}}{dp_{T}} [pb GeV^{-1}]";
    const std::vector<int> yrs{23, 24, 25, 26};
    std::string yrs_str;
    for (size_t i = 0; i < yrs.size(); ++i) yrs_str += (i ? ", " : "") + std::string("20") + std::to_string(yrs[i]);
    const std::string pbpb_l2 = DatasetTriggerMap::GetTriggerLabel(yrs.front(), "PbPb") + ", tight WP";

    pl.Draw(OsSsExample::LoadPP24(), "PP 2024, tight WP",
            DatasetTriggerMap::GetTriggerLabel(24, "pp"),
            "pp24_pair_pt_in_eta_subplots_os_ss.png", "d#sigma/dp_{T} [pb GeV^{-1}]");
    pl.Draw(OsSsExample::LoadPbPb(yrs, {"ctr0_5", "ctr5_10"}), "Pb+Pb " + yrs_str + ", 0-10%",
            pbpb_l2, "pbpb_23_24_25_26_combined_ctr0_10_pair_pt_in_eta_subplots_os_ss.png", pbpb_ytitle);
    pl.Draw(OsSsExample::LoadPbPb(yrs, {"ctr50_80"}), "Pb+Pb " + yrs_str + ", 50-80%",
            pbpb_l2, "pbpb_23_24_25_26_combined_ctr50_80_pair_pt_in_eta_subplots_os_ss.png", pbpb_ytitle);
}
