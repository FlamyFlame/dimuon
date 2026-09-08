#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../helper_functions.c"
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../Utilities/CommonLogYRange.h"
#include "../../Utilities/PairEtaPanelBins.h"
#include "../../Utilities/proj_range_to_suffix.cxx"
#include "../../MuonObjectsParamsAndHelpers/DatasetTriggerMap.h"
#include "../../MuonObjectsParamsAndHelpers/PPBaseClass.h"

class SingleBCrossxPlotterBase {
protected:
    std::string input_file_path;
    std::string output_dir;
    TFile* fin{nullptr};

    CommonEffcyConfig cfg;
    std::vector<std::pair<float, float>> q_eta_bins;
    // READ from RDFBasedHistFillingPythia::dr_bins_edges_for_reco_effcy, the single source of the
    // reco-efficiency dR binning, instead of retyping the edges (.claude/CLAUDE.md §Binnings).
    // The values happened to agree, which is exactly how two binnings coexist unnoticed.
    std::vector<std::pair<float, float>> dr_bins = [] {
        std::vector<std::pair<float, float>> v;
        const auto& e = CommonEffcyConfig{}.dr_bins_edges_for_reco_effcy;
        for (std::size_t i = 0; i + 1 < e.size(); ++i) v.emplace_back(e[i], e[i + 1]);
        return v;
    }();

    // Panel/axis alignment guard: hoisted 2026-09-07 into Utilities/PairEtaPanelBins.h, because
    // five macros OUTSIDE this class project the same panels off UNGUARDED axes and a guard only
    // some call sites use is not a guard. See that header for the full rationale (F18).
    static std::pair<int, int> PanelEtaBins(const TAxis* ax, const std::pair<float, float>& panel,
                                            const std::string& context)
    { return PairEtaPanels::Bins(ax, panel, context); }

    std::array<int, 6> line_colors{{kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2}};
    std::array<int, 6> marker_styles{{20, 21, 22, 33, 34, 29}};

public:
    SingleBCrossxPlotterBase(const std::string& in_path, const std::string& out_dir)
        : input_file_path(in_path), output_dir(out_dir), q_eta_bins(cfg.pair_eta_proj_ranges_coarse_incl_gap) {}

    virtual ~SingleBCrossxPlotterBase() {
        if (fin) {
            fin->Close();
            delete fin;
            fin = nullptr;
        }
    }

    bool Init() {
        // Check input file exists
        if (gSystem->AccessPathName(input_file_path.c_str())) {
            throw std::runtime_error("Input file does not exist: " + input_file_path);
        }

        gSystem->mkdir(output_dir.c_str(), true);
        fin = TFile::Open(input_file_path.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            throw std::runtime_error("Cannot open input file: " + input_file_path);
        }
        gStyle->SetOptStat(0);
        return true;
    }

    bool CheckHistogramExists(const std::string& hname, const std::string& htype = "TH1") {
        TObject* obj = GetHistObject(hname);
        if (!obj) {
            throw std::runtime_error("Missing histogram: " + hname);
        }
        if ((htype == "TH2D" && !dynamic_cast<TH2D*>(obj)) ||
            (htype == "TH3D" && !dynamic_cast<TH3D*>(obj))) {
            throw std::runtime_error("Histogram " + hname + " is not of type " + htype);
        }
        return true;
    }

    bool CheckHistogramNonEmpty(const std::string& hname) {
        TH1* h = dynamic_cast<TH1*>(GetHistObject(hname));
        if (!h) {
            throw std::runtime_error("Cannot retrieve histogram to check non-empty: " + hname);
        }
        if (h->GetEntries() == 0) {
            throw std::runtime_error("Histogram is empty (Entries=0): " + hname);
        }
        return true;
    }

    virtual void Run() = 0;

protected:
    virtual TObject* GetHistObject(const std::string& name) {
        if (!fin) throw std::runtime_error("Input file not opened, cannot retrieve: " + name);
        return fin->Get(name.c_str());
    }

    void Save2DColz(const std::string& hname, const std::string& png_name,
                    const std::string& z_title = "") {
        CheckHistogramExists(hname, "TH2D");
        CheckHistogramNonEmpty(hname);

        TH2D* h = dynamic_cast<TH2D*>(GetHistObject(hname));
        if (!h) {
            throw std::runtime_error("Failed to retrieve TH2D: " + hname);
        }
        // Clone so we don't modify the in-file object
        TH2D* hplot = dynamic_cast<TH2D*>(h->Clone((hname + "_plot").c_str()));
        hplot->SetDirectory(nullptr);
        // Patch axis titles to physics notation
        {
            auto fix = [](TAxis* ax) {
                const std::string t = ax->GetTitle();
                if (t == "p_{T} [GeV]")  ax->SetTitle("p_{T}^{pair} [GeV]");
                else if (t == "#eta")     ax->SetTitle("#eta^{pair}");
            };
            fix(hplot->GetXaxis());
            fix(hplot->GetYaxis());
        }
        // Make differential: divide by bin widths of both axes
        hplot->Scale(1., "width");

        TCanvas c("c2d", "c2d", 900, 750);
        c.cd();
        c.SetRightMargin(0.17);
        c.SetLeftMargin(0.12);
        c.SetBottomMargin(0.12);

        const std::string xt = hplot->GetXaxis()->GetTitle();
        const std::string yt = hplot->GetYaxis()->GetTitle();
        const bool x_is_pair_pt = xt.find("p_{T}") != std::string::npos;
        const bool y_is_pair_pt = yt.find("p_{T}") != std::string::npos;
        if (x_is_pair_pt) c.SetLogx();
        if (y_is_pair_pt) c.SetLogy();
        if (x_is_pair_pt || y_is_pair_pt) c.SetLogz();

        hplot->GetXaxis()->SetTitleSize(0.05);
        hplot->GetYaxis()->SetTitleSize(0.05);
        hplot->GetZaxis()->SetTitleSize(0.045);
        hplot->GetXaxis()->SetLabelSize(0.04);
        hplot->GetYaxis()->SetLabelSize(0.04);
        hplot->GetZaxis()->SetLabelSize(0.035);
        hplot->GetYaxis()->SetTitleOffset(1.35);

        hplot->GetZaxis()->SetTitle(z_title.empty() ? "d^{2}N/d[X]d[Y]" : z_title.c_str());
        hplot->GetZaxis()->SetTitleOffset(1.35);
        hplot->Draw("colz");
        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;
        delete hplot;
    }

    void DrawPairPtByEtaWithDrLines(
        const std::string& h3_name,
        const std::string& data_info_line1,
        const std::string& data_info_line2,
        const std::string& png_name,
        const std::string& y_title = "d#sigma/dp_{T} [pb GeV^{-1}]")
    {
        CheckHistogramExists(h3_name, "TH3D");
        CheckHistogramNonEmpty(h3_name);

        TH3D* h3 = dynamic_cast<TH3D*>(GetHistObject(h3_name));
        if (!h3) {
            throw std::runtime_error("Failed to retrieve TH3D: " + h3_name);
        }

        int nrow = 1;
        int ncol = 1;
        DetermineSubplotGrid(static_cast<int>(q_eta_bins.size()), nrow, ncol);

        TCanvas c("cpt", "pair_pt by eta with dr lines", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        std::vector<std::vector<TH1D*>> all_lines(q_eta_bins.size());
        std::vector<TLegend*> all_legends;
        all_legends.reserve(q_eta_bins.size());

        // PASS 1 — build every projection of every panel WITHOUT drawing, so the
        // common log-y range can be derived from all of them before the first
        // frame is painted (see Utilities/CommonLogYRange.h).
        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            const auto& eta_bin = q_eta_bins.at(ieta);
            const auto ybins = PanelEtaBins(h3->GetYaxis(), eta_bin,
                                            "SingleBCrossxPlotterBase (3D dR-in-eta panels)");
            const int y1 = ybins.first, y2 = ybins.second;

            std::vector<TH1D*> lines;
            lines.reserve(dr_bins.size());

            for (size_t idr = 0; idr < dr_bins.size(); ++idr) {
                const auto& dr_bin = dr_bins.at(idr);
                const int z1 = h3->GetZaxis()->FindBin(dr_bin.first + 1e-6);
                const int z2 = h3->GetZaxis()->FindBin(dr_bin.second - 1e-6);

                const std::string proj_name = "hpt_eta" + std::to_string(ieta) + "_dr" + std::to_string(idr) + "_" + std::to_string(std::rand());
                TH1D* hp = h3->ProjectionX(proj_name.c_str(), y1, y2, z1, z2, "e");
                hp->SetDirectory(nullptr);
                hp->Scale(1.0, "width");
                hp->SetLineWidth(2);
                hp->SetLineColor(line_colors.at(idr % line_colors.size()));
                hp->SetMarkerColor(line_colors.at(idr % line_colors.size()));
                hp->SetMarkerStyle(marker_styles.at(idr % marker_styles.size()));
                hp->SetMarkerSize(0.9);
                hp->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
                hp->GetYaxis()->SetTitle(y_title.c_str());
                hp->GetXaxis()->SetTitleSize(0.06);
                hp->GetYaxis()->SetTitleSize(0.06);
                hp->GetXaxis()->SetLabelSize(0.05);
                hp->GetYaxis()->SetLabelSize(0.05);
                hp->GetYaxis()->SetTitleOffset(1.45);
                hp->SetTitle("");

                lines.push_back(hp);
            }

            all_lines.at(ieta) = std::move(lines);
        }

        // ONE log-y scale for the whole PNG, low enough to keep every non-empty
        // point of every panel AND every dR curve inside the frame. The scan must
        // cover all curves, not just the first-drawn one: the pad's frame is
        // defined by the first histogram, so a lower point on an overlaid dR
        // curve would otherwise be silently drawn off-frame.
        {
            std::vector<TH1*> flat;
            for (const auto& lines : all_lines)
                for (TH1D* h : lines) flat.push_back(h);
            ApplyCommonLogYRange(flat);
        }

        // PASS 2 — draw.
        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            auto& lines = all_lines.at(ieta);
            if (lines.empty()) continue;

            c.cd(static_cast<int>(ieta) + 1);
            gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.13);

            const auto& eta_bin = q_eta_bins.at(ieta);

            lines.at(0)->Draw("E1");
            for (size_t il = 1; il < lines.size(); ++il) {
                lines.at(il)->Draw("E1 SAME");
            }

            // Info text (no symbols): right-aligned, minimal margin
            TLegend* leg_info = new TLegend(0.38, 0.72, 0.93, 0.90);
            leg_info->SetBorderSize(0);
            leg_info->SetFillStyle(0);
            leg_info->SetTextSize(0.045);
            leg_info->SetTextAlign(32);
            leg_info->SetMargin(0.01);
            leg_info->AddEntry((TObject*)0, data_info_line1.c_str(), "");
            leg_info->AddEntry((TObject*)0, data_info_line2.c_str(), "");
            leg_info->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eta_bin.first, eta_bin.second), "");
            leg_info->Draw();
            all_legends.push_back(leg_info);

            // dR lines: narrow legend, left-aligned text, symbol just left of text
            TLegend* leg_dr = new TLegend(0.70, 0.48, 0.93, 0.72);
            leg_dr->SetBorderSize(0);
            leg_dr->SetFillStyle(0);
            leg_dr->SetTextSize(0.045);
            leg_dr->SetTextAlign(12);
            leg_dr->SetMargin(0.22);
            for (size_t idr = 0; idr < dr_bins.size(); ++idr) {
                const auto& dr_bin = dr_bins.at(idr);
                leg_dr->AddEntry(lines.at(idr), Form("#DeltaR #in [%.1f, %.1f]", dr_bin.first, dr_bin.second), "lep");
            }
            leg_dr->Draw();
            all_legends.push_back(leg_dr);
        }

        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;

        for (auto& lines : all_lines) {
            for (TH1D* h : lines) delete h;
        }
        for (TLegend* leg : all_legends) delete leg;
    }

    // dR-integrated pair-pT vs pair-eta: one subplot per eta bin, single line per subplot.
    // h2_name: TH2D with pair_pt on X, pair_eta on Y.
    // differential=true: divide by bin width (d/dp_T); false: raw counts per eta bin.
    void DrawPairPtByEta(
        const std::string& h2_name,
        const std::string& data_info_line1,
        const std::string& data_info_line2,
        const std::string& png_name,
        const std::string& y_title = "d#sigma/dp_{T} [pb GeV^{-1}]",
        bool differential = true)
    {
        CheckHistogramExists(h2_name, "TH2D");
        CheckHistogramNonEmpty(h2_name);

        TH2D* h2 = dynamic_cast<TH2D*>(GetHistObject(h2_name));
        if (!h2) {
            throw std::runtime_error("Failed to retrieve TH2D: " + h2_name);
        }

        int nrow = 1, ncol = 1;
        DetermineSubplotGrid(static_cast<int>(q_eta_bins.size()), nrow, ncol);

        TCanvas c("cpt_eta", "pair_pt by eta dR-integrated", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        std::vector<TH1D*> all_hists;
        all_hists.reserve(q_eta_bins.size());
        std::vector<TLegend*> all_legends;
        all_legends.reserve(q_eta_bins.size());

        // PASS 1 — build every panel's projection WITHOUT drawing, so the common
        // log-y range can be derived from all panels (Utilities/CommonLogYRange.h).
        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            const auto& eta_bin = q_eta_bins.at(ieta);
            const auto ybins = PanelEtaBins(h2->GetYaxis(), eta_bin,
                                            "SingleBCrossxPlotterBase (2D pT-in-eta panels)");
            const int y1 = ybins.first, y2 = ybins.second;

            const std::string proj_name = "hpt_eta" + std::to_string(ieta) + "_"
                                          + std::to_string(std::rand());
            TH1D* hp = h2->ProjectionX(proj_name.c_str(), y1, y2, "e");
            hp->SetDirectory(nullptr);
            if (differential) hp->Scale(1.0, "width");
            hp->SetLineWidth(2);
            hp->SetLineColor(kBlack);
            hp->SetMarkerColor(kBlack);
            hp->SetMarkerStyle(20);
            hp->SetMarkerSize(0.9);
            hp->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            hp->GetYaxis()->SetTitle(y_title.c_str());
            hp->GetXaxis()->SetTitleSize(0.06);
            hp->GetYaxis()->SetTitleSize(0.06);
            hp->GetXaxis()->SetLabelSize(0.05);
            hp->GetYaxis()->SetLabelSize(0.05);
            hp->GetYaxis()->SetTitleOffset(1.45);
            hp->SetTitle("");

            all_hists.push_back(hp);
        }

        // ONE log-y scale for the whole PNG, low enough to keep every non-empty
        // point of every panel inside the frame.
        ApplyCommonLogYRange(std::vector<TH1*>(all_hists.begin(), all_hists.end()));

        // PASS 2 — draw.
        for (size_t ieta = 0; ieta < all_hists.size(); ++ieta) {
            c.cd(static_cast<int>(ieta) + 1);
            gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.13);

            const auto& eta_bin = q_eta_bins.at(ieta);
            all_hists.at(ieta)->Draw("E1");

            // x2 is the FRAME edge, not 0.93: the pad keeps ROOT's default right margin 0.1, so
            // a right-aligned (SetTextAlign(32)) legend ending at 0.93 pushes its last character
            // 0.03 NDC OUTSIDE the frame and the frame line is drawn through it. Found in review
            // 2026-09-03 on both this figure and the counts figure that reuses this method.
            TLegend* leg = new TLegend(0.60, 0.70, 1.0 - gPad->GetRightMargin(), 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.045);
            leg->SetTextAlign(32);
            leg->SetMargin(0.06);
            leg->AddEntry((TObject*)0, data_info_line1.c_str(), "");
            leg->AddEntry((TObject*)0, data_info_line2.c_str(), "");
            leg->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eta_bin.first, eta_bin.second), "");
            leg->Draw();
            all_legends.push_back(leg);
        }

        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;

        for (TH1D* h : all_hists) delete h;
        for (TLegend* leg : all_legends) delete leg;
    }

    // =============================================================================================
    // TEMPORARY DIAGNOSTIC helpers — muon reconstructed pT > 4.5 GeV cut study
    // (docs/tracking/muon_pt45_cut_diagnostic.md). Three generic two-series overlay Draw methods,
    // modeled on DrawPairPtByEtaWithDrLines (N-series-from-TH3D pattern) and DrawPairPtByEta
    // (panel-projection logic), so no panel/binning code is duplicated. NOT tied to the muon-pT
    // study specifically -- they take two histogram NAMES and two labels, so any two-variant
    // comparison on the SAME pair-pT x pair-eta axes can reuse them. If the 4.5 GeV cut is dropped
    // as a diagnostic dead end, these three methods and their two driver macros can be deleted
    // without touching anything else in this file.
    // =============================================================================================

    // Panel-per-pair-eta (9 canonical panels), two-series overlay: pair-pT dependence in each
    // pair-eta bin. h2_name_a/b must share the SAME pair-pT x pair-eta axes.
    void DrawPairPtByEtaTwoSeries(
        const std::string& h2_name_a, const std::string& h2_name_b,
        const std::string& label_a, const std::string& label_b,
        const std::string& data_info_line1, const std::string& data_info_line2,
        const std::string& png_name,
        const std::string& y_title = "N_{pairs}")
    {
        CheckHistogramExists(h2_name_a, "TH2D");
        CheckHistogramExists(h2_name_b, "TH2D");
        CheckHistogramNonEmpty(h2_name_a);
        CheckHistogramNonEmpty(h2_name_b);

        TH2D* ha = dynamic_cast<TH2D*>(GetHistObject(h2_name_a));
        TH2D* hb = dynamic_cast<TH2D*>(GetHistObject(h2_name_b));
        if (!ha || !hb) {
            throw std::runtime_error("DrawPairPtByEtaTwoSeries: failed to retrieve TH2D(s): "
                                     + h2_name_a + ", " + h2_name_b);
        }

        int nrow = 1, ncol = 1;
        DetermineSubplotGrid(static_cast<int>(q_eta_bins.size()), nrow, ncol);

        TCanvas c("cpt_eta_2s", "pair_pt by eta, two series", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        // Legend-icon PROXIES, styled like the two series but carrying no real bin content/error.
        // NOT the real per-panel histograms: TLegend's "lep" icon for a referenced TH1 with a
        // real (tiny, sqrt(N)) error on a log-y pad was found to render as a full-frame-height
        // spike instead of a small icon tick (review 2026-09-06 self-check on real PbPb output,
        // reproduced with independently-verified clean bin content/error -- a TLegend/log-pad
        // rendering quirk, not a data bug). A 1-bin, zero-content, zero-error dummy sidesteps it.
        TH1D proxy_a("proxy_pteta2s_a", "", 1, 0, 1);
        TH1D proxy_b("proxy_pteta2s_b", "", 1, 0, 1);
        proxy_a.SetLineWidth(2); proxy_a.SetLineColor(line_colors.at(0));
        proxy_a.SetMarkerColor(line_colors.at(0)); proxy_a.SetMarkerStyle(marker_styles.at(0));
        proxy_b.SetLineWidth(2); proxy_b.SetLineColor(line_colors.at(1));
        proxy_b.SetMarkerColor(line_colors.at(1)); proxy_b.SetMarkerStyle(marker_styles.at(1));

        std::vector<std::array<TH1D*, 2>> all_lines(q_eta_bins.size());
        std::vector<TLegend*> all_legends;

        // PASS 1 — build every panel's two projections WITHOUT drawing, so the common log-y
        // range can be derived from all of them first (Utilities/CommonLogYRange.h).
        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            const auto& eta_bin = q_eta_bins.at(ieta);
            const auto ybins_a = PanelEtaBins(ha->GetYaxis(), eta_bin,
                                              "SingleBCrossxPlotterBase (2-sample overlay, A)");
            const auto ybins_b = PanelEtaBins(hb->GetYaxis(), eta_bin,
                                              "SingleBCrossxPlotterBase (2-sample overlay, B)");
            const int y1a = ybins_a.first, y2a = ybins_a.second;
            const int y1b = ybins_b.first, y2b = ybins_b.second;

            const std::string nA = "hpt2s_a_eta" + std::to_string(ieta) + "_" + std::to_string(std::rand());
            const std::string nB = "hpt2s_b_eta" + std::to_string(ieta) + "_" + std::to_string(std::rand());
            TH1D* hpa = ha->ProjectionX(nA.c_str(), y1a, y2a, "e");
            TH1D* hpb = hb->ProjectionX(nB.c_str(), y1b, y2b, "e");
            hpa->SetDirectory(nullptr);
            hpb->SetDirectory(nullptr);

            auto style = [&](TH1D* h, int idx) {
                h->SetLineWidth(2);
                h->SetLineColor(line_colors.at(idx));
                h->SetMarkerColor(line_colors.at(idx));
                h->SetMarkerStyle(marker_styles.at(idx));
                h->SetMarkerSize(0.9);
                h->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
                h->GetYaxis()->SetTitle(y_title.c_str());
                h->GetXaxis()->SetTitleSize(0.06);
                h->GetYaxis()->SetTitleSize(0.06);
                h->GetXaxis()->SetLabelSize(0.05);
                h->GetYaxis()->SetLabelSize(0.05);
                h->GetYaxis()->SetTitleOffset(1.45);
                h->SetTitle("");
            };
            style(hpa, 0);
            style(hpb, 1);

            all_lines.at(ieta) = {hpa, hpb};
        }

        {
            std::vector<TH1*> flat;
            for (auto& pr : all_lines) { flat.push_back(pr[0]); flat.push_back(pr[1]); }
            ApplyCommonLogYRange(flat);
        }

        // PASS 2 — draw.
        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            c.cd(static_cast<int>(ieta) + 1);
            gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.13);

            const auto& eta_bin = q_eta_bins.at(ieta);
            all_lines.at(ieta)[0]->Draw("E1");
            all_lines.at(ieta)[1]->Draw("E1 SAME");

            // Info text (no symbols). x2 is the FRAME edge (1 - right margin), not 0.93 -- a
            // right-aligned legend ending at a hardcoded 0.93 pushes past the frame when the pad's
            // default right margin (0.1) is in effect (found in review 2026-09-03).
            // Geometry mirrors DrawPairPtByEtaWithDrLines' two-box layout above (info box wider
            // and higher, symbol box narrower/further right so it clears the falling spectrum's
            // tail): found overlapping data at the wider/lower placement tried first (review
            // 2026-09-06 self-check on real pp24 output).
            TLegend* leg_info = new TLegend(0.38, 0.72, 1.0 - gPad->GetRightMargin(), 0.90);
            leg_info->SetBorderSize(0);
            leg_info->SetFillStyle(0);
            leg_info->SetTextSize(0.042);
            leg_info->SetTextAlign(32);
            leg_info->SetMargin(0.01);
            leg_info->AddEntry((TObject*)0, data_info_line1.c_str(), "");
            leg_info->AddEntry((TObject*)0, data_info_line2.c_str(), "");
            leg_info->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eta_bin.first, eta_bin.second), "");
            leg_info->Draw();
            all_legends.push_back(leg_info);

            // Series legend: narrow, left-aligned text, symbol just left of text (matches the
            // dR-lines legend pattern), positioned further right/lower so it sits over the
            // spectrum's low-value tail rather than its descending shoulder.
            TLegend* leg_series = new TLegend(0.70, 0.58, 1.0 - gPad->GetRightMargin(), 0.72);
            leg_series->SetBorderSize(0);
            leg_series->SetFillStyle(0);
            leg_series->SetTextSize(0.045);
            leg_series->SetTextAlign(12);
            leg_series->SetMargin(0.22);
            leg_series->AddEntry(&proxy_a, label_a.c_str(), "p");
            leg_series->AddEntry(&proxy_b, label_b.c_str(), "p");
            leg_series->Draw();
            all_legends.push_back(leg_series);
        }

        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;

        for (auto& pr : all_lines) { delete pr[0]; delete pr[1]; }
        for (TLegend* leg : all_legends) delete leg;
    }

    // Pair-eta dependence, pair-pT integrated: ONE panel, X = the 9 canonical pair-eta panels
    // (variable-width bins built from q_eta_bins, which is contiguous by construction), two-series
    // overlay. Row sums of h2_name_a/b over the full pair-pT axis, per panel.
    void DrawPairEtaIntegratedTwoSeries(
        const std::string& h2_name_a, const std::string& h2_name_b,
        const std::string& label_a, const std::string& label_b,
        const std::string& data_info_line1, const std::string& data_info_line2,
        const std::string& png_name,
        const std::string& y_title = "N_{pairs}")
    {
        CheckHistogramExists(h2_name_a, "TH2D");
        CheckHistogramExists(h2_name_b, "TH2D");
        CheckHistogramNonEmpty(h2_name_a);
        CheckHistogramNonEmpty(h2_name_b);

        TH2D* ha = dynamic_cast<TH2D*>(GetHistObject(h2_name_a));
        TH2D* hb = dynamic_cast<TH2D*>(GetHistObject(h2_name_b));
        if (!ha || !hb) {
            throw std::runtime_error("DrawPairEtaIntegratedTwoSeries: failed to retrieve TH2D(s): "
                                     + h2_name_a + ", " + h2_name_b);
        }

        // Variable bin edges from the 9 canonical panels -- contiguous by construction
        // (pair_eta_proj_ranges_coarse_incl_gap), never retyped.
        std::vector<double> edges;
        edges.reserve(q_eta_bins.size() + 1);
        edges.push_back(q_eta_bins.front().first);
        for (const auto& b : q_eta_bins) edges.push_back(b.second);

        auto build = [&](TH2D* h2, const std::string& nm) {
            TH1D* h1 = new TH1D(nm.c_str(), "", (int)q_eta_bins.size(), edges.data());
            h1->Sumw2();
            for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
                const auto& eta_bin = q_eta_bins.at(ieta);
                const auto ybins = PanelEtaBins(h2->GetYaxis(), eta_bin,
                                                "SingleBCrossxPlotterBase (eta-panel ratios)");
                const int y1 = ybins.first, y2 = ybins.second;
                double err = 0.;
                const double n = h2->IntegralAndError(1, h2->GetNbinsX(), y1, y2, err);
                h1->SetBinContent((int)ieta + 1, n);
                h1->SetBinError((int)ieta + 1, err);
            }
            h1->SetDirectory(nullptr);
            return h1;
        };

        TH1D* hpa = build(ha, "heta2s_a_" + std::to_string(std::rand()));
        TH1D* hpb = build(hb, "heta2s_b_" + std::to_string(std::rand()));

        auto style = [&](TH1D* h, int idx) {
            h->SetLineWidth(2);
            h->SetLineColor(line_colors.at(idx));
            h->SetMarkerColor(line_colors.at(idx));
            h->SetMarkerStyle(marker_styles.at(idx));
            h->SetMarkerSize(1.0);
            h->GetXaxis()->SetTitle("#eta^{pair}");
            h->GetYaxis()->SetTitle(y_title.c_str());
            h->GetYaxis()->SetNoExponent(kFALSE);
            h->SetTitle("");
        };
        style(hpa, 0);
        style(hpb, 1);

        // Scientific notation on the y-axis (e.g. "1 #times 10^{4}"): TAxis has no per-label
        // format hook (unlike TGraph/TF1), so this is ROOT's actual mechanism -- a shared
        // "#times 10^{p}" header with reduced tick values -- forced on for BOTH plots via
        // TGaxis::SetMaxDigits (global/static, so save+restore around this Draw call only).
        // Without this, ROOT's default threshold (5 digits) fires for pp24's 6-digit range but
        // NOT PbPb's 5-digit range, leaving one plot with the header and the other with bare
        // long integers -- an inconsistent look between the two datasets this forces to agree.
        const int prev_max_digits = TGaxis::GetMaxDigits();
        TGaxis::SetMaxDigits(3);

        TCanvas c("ceta_2s", "pair eta integrated, two series", 700, 550);
        c.SetLeftMargin(0.17);
        c.SetBottomMargin(0.12);

        // LINEAR y (2026-09-07, corrected from log-y): pair_eta is not a momentum-like
        // observable and its axis is not log-binned, so per the linear-y-default rule
        // (.claude/commands/review-plot.md R4 / feedback_log_scale_plots memory) this plot
        // takes a linear y-scale by default -- there is no power-law or log-binning
        // justification for log here. This quantity's shape vs pair-eta is NOT monotonic (a
        // broad hill with a central dip at the eta~0 gap cut) and its scale differs between
        // datasets (pp24 ~5x, PbPb ~10x less), so no FIXED corner is safe for both: a
        // bottom-left box that cleared pp24 overlapped PbPb's central bump (review 2026-09-06
        // self-check on real output, when this was still log-y). Keep the same fix -- headroom
        // ABOVE every point regardless of shape or dataset, then place both legends there (the
        // muon_gap_cuts_acceptance.md F13/F14 fix for the same class of problem on q*eta plots)
        // -- just computed directly for a linear axis instead of via the log-only
        // ApplyCommonLogYRange (Utilities/CommonLogYRange.h; its floor/ceil padding and
        // positive-bin-only scan are log-scale-specific and do not apply here).
        double eta2s_ymax = 0.0;
        for (const TH1D* h : {hpa, hpb}) {
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                eta2s_ymax = std::max(eta2s_ymax, h->GetBinContent(b));
            }
        }
        // 1.4x (not 1.35x): a /review-plot pass (2026-09-07) pixel-measured the legend's top
        // row bisected by the frame's top border at 1.35x with the legend box unmoved -- widen
        // the headroom AND move the legend down (below) to guarantee clearance from both the
        // frame border above and the tallest point below.
        for (TH1D* h : {hpa, hpb}) {
            h->SetMinimum(0.0);
            h->SetMaximum(eta2s_ymax * 1.4);
        }

        hpa->Draw("E1");
        hpb->Draw("E1 SAME");

        // Legend-icon PROXY objects (not hpa/hpb) -- see DrawPairPtByEtaTwoSeries for why: a
        // referenced real histogram's "lep" icon can render as a full-frame spike rather than a
        // small tick. ONE combined TLegend (info text + series icons), not two separate boxes:
        // two separate TLegend objects on this log-y pad -- one text-only, one symbol-bearing --
        // was found to reproduce the same spike even with proxy icons (review 2026-09-06
        // self-check on real PbPb output; isolated with a minimal standalone repro), while the
        // panel method's two-box layout above does not show it. Root cause not fully understood
        // (a ROOT TLegend/log-pad rendering interaction), so the single-box layout the panel
        // method's own precedent (DrawPairPtByEta) already uses is the safer, tested pattern here.
        TH1D proxy_a("proxy_eta2s_a", "", 1, 0, 1);
        TH1D proxy_b("proxy_eta2s_b", "", 1, 0, 1);
        proxy_a.SetLineWidth(2); proxy_a.SetLineColor(line_colors.at(0));
        proxy_a.SetMarkerColor(line_colors.at(0)); proxy_a.SetMarkerStyle(marker_styles.at(0));
        proxy_b.SetLineWidth(2); proxy_b.SetLineColor(line_colors.at(1));
        proxy_b.SetMarkerColor(line_colors.at(1)); proxy_b.SetMarkerStyle(marker_styles.at(1));

        // Top bound 0.90 = default frame top (ROOT top margin 0.1, not overridden on this
        // canvas), matching the safe bound already used elsewhere in this file (e.g. the
        // panel-method legends above); bottom bound 0.73 clears the tallest point, which the
        // 1.4x headroom above places at NDC ~0.68 regardless of dataset shape or scale.
        TLegend leg(0.15, 0.73, 0.97, 0.90);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetNColumns(2);
        leg.SetTextSize(0.032);
        leg.SetTextAlign(12);
        leg.SetMargin(0.18);
        leg.AddEntry((TObject*)0, data_info_line1.c_str(), "");
        leg.AddEntry(&proxy_a, label_a.c_str(), "p");
        leg.AddEntry((TObject*)0, data_info_line2.c_str(), "");
        leg.AddEntry(&proxy_b, label_b.c_str(), "p");
        leg.Draw();

        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;

        TGaxis::SetMaxDigits(prev_max_digits);  // restore -- static/global setting

        delete hpa;
        delete hpb;
    }

    // Pair-pT dependence, pair-eta integrated: ONE panel, X = the pT axis of h2_name_a/b
    // (ParamsSet::pT_bins_120 by construction of the input histograms), two-series overlay.
    void DrawPairPtIntegratedTwoSeries(
        const std::string& h2_name_a, const std::string& h2_name_b,
        const std::string& label_a, const std::string& label_b,
        const std::string& data_info_line1, const std::string& data_info_line2,
        const std::string& png_name,
        const std::string& y_title = "N_{pairs}")
    {
        CheckHistogramExists(h2_name_a, "TH2D");
        CheckHistogramExists(h2_name_b, "TH2D");
        CheckHistogramNonEmpty(h2_name_a);
        CheckHistogramNonEmpty(h2_name_b);

        TH2D* ha = dynamic_cast<TH2D*>(GetHistObject(h2_name_a));
        TH2D* hb = dynamic_cast<TH2D*>(GetHistObject(h2_name_b));
        if (!ha || !hb) {
            throw std::runtime_error("DrawPairPtIntegratedTwoSeries: failed to retrieve TH2D(s): "
                                     + h2_name_a + ", " + h2_name_b);
        }

        const std::string nA = "hpt1d_a_" + std::to_string(std::rand());
        const std::string nB = "hpt1d_b_" + std::to_string(std::rand());
        TH1D* hpa = ha->ProjectionX(nA.c_str(), 1, ha->GetNbinsY(), "e");
        TH1D* hpb = hb->ProjectionX(nB.c_str(), 1, hb->GetNbinsY(), "e");
        hpa->SetDirectory(nullptr);
        hpb->SetDirectory(nullptr);

        auto style = [&](TH1D* h, int idx) {
            h->SetLineWidth(2);
            h->SetLineColor(line_colors.at(idx));
            h->SetMarkerColor(line_colors.at(idx));
            h->SetMarkerStyle(marker_styles.at(idx));
            h->SetMarkerSize(1.0);
            h->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            h->GetYaxis()->SetTitle(y_title.c_str());
            h->SetTitle("");
        };
        style(hpa, 0);
        style(hpb, 1);

        TCanvas c("cpt1d_2s", "pair pT integrated, two series", 700, 550);
        c.SetLogx();
        c.SetLogy();
        c.SetLeftMargin(0.13);
        c.SetBottomMargin(0.12);

        // Extra ceil_pad headroom (see DrawPairEtaIntegratedTwoSeries for the same mechanism):
        // the falling spectrum's top-left is occupied by the peak, so the legend needs both
        // genuine top clearance AND a box wide enough for its longest label -- a box that fits
        // "PP 2024, tight WP" clipped past the frame for the longer "Pb+Pb 2023+2024+2025
        // combined, tight WP" info line (review 2026-09-06 self-check on real PbPb output).
        ApplyCommonLogYRange({hpa, hpb}, /*ceil_pad=*/3.0);

        hpa->Draw("E1");
        hpb->Draw("E1 SAME");

        TLegend leg_info(0.40, 0.80, 0.97, 0.92);
        leg_info.SetBorderSize(0);
        leg_info.SetFillStyle(0);
        leg_info.SetTextSize(0.032);
        leg_info.AddEntry((TObject*)0, data_info_line1.c_str(), "");
        leg_info.AddEntry((TObject*)0, data_info_line2.c_str(), "");
        leg_info.Draw();

        // Legend-icon PROXY objects (not hpa/hpb) -- see DrawPairPtByEtaTwoSeries for why: a
        // referenced real histogram's "lep" icon on a log-y pad was found to render as a
        // full-frame spike rather than a small tick.
        TH1D proxy_a("proxy_pt1d2s_a", "", 1, 0, 1);
        TH1D proxy_b("proxy_pt1d2s_b", "", 1, 0, 1);
        proxy_a.SetLineWidth(2); proxy_a.SetLineColor(line_colors.at(0));
        proxy_a.SetMarkerColor(line_colors.at(0)); proxy_a.SetMarkerStyle(marker_styles.at(0));
        proxy_b.SetLineWidth(2); proxy_b.SetLineColor(line_colors.at(1));
        proxy_b.SetMarkerColor(line_colors.at(1)); proxy_b.SetMarkerStyle(marker_styles.at(1));

        TLegend leg_series(0.40, 0.65, 0.97, 0.79);
        leg_series.SetBorderSize(0);
        leg_series.SetFillStyle(0);
        leg_series.SetTextSize(0.036);
        leg_series.SetTextAlign(12);
        leg_series.SetMargin(0.15);
        leg_series.AddEntry(&proxy_a, label_a.c_str(), "p");
        leg_series.AddEntry(&proxy_b, label_b.c_str(), "p");
        leg_series.Draw();

        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;

        delete hpa;
        delete hpb;
    }
};
