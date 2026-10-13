#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../helper_functions.c"
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../Utilities/PairEtaPanelBins.h"
#include "../../Utilities/SignalAcceptanceDefinition.h"
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"

class SignalAcceptancePlotter {
protected:
    std::string input_file_path;
    std::string output_dir;
    TFile* fin{nullptr};

    CommonEffcyConfig cfg;
    ParamsSet pms;
    std::vector<std::pair<float, float>> q_eta_bins;

public:
    // WHICH pair-pT axis (user decision 2026-09-08; mirrors SingleBCrossxPlotterBase):
    //   false (DEFAULT) -> the UNSUFFIXED "h2d_sig_accept_*_pt_eta" family, booked on
    //                      ParamsSet::pT_bins_150 (16 log bins 9 -> 150 GeV), into the unsuffixed
    //                      output directory. The signal acceptance multiplies the cross-section,
    //                      so it MUST be binned on the cross-section's nominal axis.
    //   true  (OPT-IN)  -> the "h2d_sig_accept_*_pt_120_eta" family, booked on
    //                      ParamsSet::pT_bins_120 (16 log bins 9 -> 120 GeV), into "_pt_120".
    // Producers: RDFBasedHistFillingPythiaTruth.cxx / RDFBasedHistFillingPowhegTruth.cxx.
    bool use_pt_bins_120 = false;

    std::string PtAxisDirSuffix() const { return use_pt_bins_120 ? "_pt_120" : ""; }

    // WHICH acceptance factor (TARGET definition, docs/tracking/signal_acceptance_mass_window.md
    // PD4/D7): "" = the full A_sig numerator h2d_sig_accept_num_*; "mupt" / "mueta" / "mu" / "mass"
    // = the single-factor numerators h2d_sig_accept_num_<tag>_* (target-mode files only). The
    // denominator is always h2d_sig_accept_denom_*.
    std::string num_tag;
    // The definition marker the input MUST carry (Value() = interim, TargetValue() = target).
    std::string expected_marker = SignalAcceptanceDefinition::Value();
    std::string y_title = "#alpha_{signal}";

    std::string AxisSuffix() const { return use_pt_bins_120 ? "pt_120_eta" : "pt_eta"; }
    std::string NumName() const { return "h2d_sig_accept_num_" + (num_tag.empty() ? std::string() : num_tag + "_") + AxisSuffix(); }
    std::string DenName() const { return "h2d_sig_accept_denom_" + AxisSuffix(); }

    SignalAcceptancePlotter(const std::string& in_path, const std::string& out_dir)
        : input_file_path(in_path), output_dir(out_dir),
          q_eta_bins(cfg.pair_eta_proj_ranges_coarse_incl_gap) {}

    virtual ~SignalAcceptancePlotter() {
        if (fin) { fin->Close(); delete fin; fin = nullptr; }
    }

    bool Init() {
        if (gSystem->AccessPathName(input_file_path.c_str())) {
            throw std::runtime_error("Input file not found: " + input_file_path);
        }
        gSystem->mkdir(output_dir.c_str(), true);
        fin = TFile::Open(input_file_path.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            throw std::runtime_error("Cannot open input file: " + input_file_path);
        }
        // A_sig definition guard (Utilities/SignalAcceptanceDefinition.h): a file filled before
        // 2026-10-05 holds the SAME histogram names with the old definition (gap cuts in the
        // numerator, no measured-range cut in the denominator) -- refuse it rather than mislabel.
        TNamed* marker = dynamic_cast<TNamed*>(fin->Get(SignalAcceptanceDefinition::Key()));
        if (!marker || std::string(marker->GetTitle()) != expected_marker) {
            throw std::runtime_error(std::string("SignalAcceptancePlotter: ") + input_file_path
                + " lacks the expected '" + SignalAcceptanceDefinition::Key()
                + "' marker (interim vs target definition, or stale) -- refill it with the current truth RDF producer.");
        }
        gStyle->SetOptStat(0);
        return true;
    }

    // Restrict a pair-eta axis to the measured range |eta^pair| < pair_eta_fiducial_max BY BIN
    // INDEX: SetRangeUser(+-2.2f) keeps one extra 0.1 bin per side (2.2f = 2.2000000477).
    static void SetMeasuredEtaRange(TAxis* ax) {
        const double half = 0.5 * ax->GetBinWidth(1);
        const int j1 = ax->FindFixBin(-ParamsSet::pair_eta_fiducial_max + half);
        const int j2 = ax->FindFixBin( ParamsSet::pair_eta_fiducial_max - half);
        ax->SetRange(j1, j2);
    }

    // 2D colz plot of the acceptance ratio histogram.
    void Draw2DColz(const std::string& png_name) {
        TH2D* hplot = nullptr;
        if (num_tag.empty()) {
            const std::string hname = use_pt_bins_120 ? "h2d_sig_accept_pt_120_eta" : "h2d_sig_accept_pt_eta";
            TH2D* h = dynamic_cast<TH2D*>(fin->Get(hname.c_str()));
            if (!h) throw std::runtime_error("Missing histogram: " + hname);
            if (h->GetEntries() == 0) throw std::runtime_error(hname + " is empty");
            hplot = dynamic_cast<TH2D*>(h->Clone((hname + "_plot").c_str()));
        } else {   // single-factor numerators carry no stored ratio: num / den here
            TH2D* hn = dynamic_cast<TH2D*>(fin->Get(NumName().c_str()));
            TH2D* hd = dynamic_cast<TH2D*>(fin->Get(DenName().c_str()));
            if (!hn) throw std::runtime_error("Missing histogram: " + NumName());
            if (!hd) throw std::runtime_error("Missing histogram: " + DenName());
            if (hd->GetEntries() == 0) throw std::runtime_error(DenName() + " is empty");
            hplot = dynamic_cast<TH2D*>(hn->Clone((NumName() + "_ratio_plot").c_str()));
            hplot->Divide(hd);
        }
        hplot->SetDirectory(nullptr);

        TCanvas c("c2d", "", 900, 750);
        c.SetRightMargin(0.17);
        c.SetLeftMargin(0.12);
        c.SetBottomMargin(0.12);
        c.SetLogx();

        hplot->GetXaxis()->SetTitleSize(0.05);
        hplot->GetYaxis()->SetTitleSize(0.05);
        hplot->GetZaxis()->SetTitleSize(0.045);
        hplot->GetXaxis()->SetLabelSize(0.04);
        hplot->GetYaxis()->SetLabelSize(0.04);
        hplot->GetZaxis()->SetLabelSize(0.035);
        hplot->GetYaxis()->SetTitleOffset(1.35);
        hplot->GetZaxis()->SetTitleOffset(1.35);
        hplot->GetZaxis()->SetTitle(y_title.c_str());
        hplot->SetMinimum(0.);
        hplot->SetMaximum(1.);
        // Measured range only: |eta^pair| >= pair_eta_fiducial_max is cut from num AND den.
        SetMeasuredEtaRange(hplot->GetYaxis());
        hplot->Draw("colz");

        const std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;
        delete hplot;
    }

    // 1D acceptance vs truth pair pT, one subplot per coarse eta bin.
    // Projects num and denom 2D histograms per eta bin, divides bin-by-bin.
    void DrawAcceptancePtByEta(const std::string& png_name,
                                const std::string& mc_label_line1,
                                const std::string& mc_label_line2) {
        const std::string num_name = NumName();
        const std::string den_name = DenName();
        TH2D* hnum = dynamic_cast<TH2D*>(fin->Get(num_name.c_str()));
        TH2D* hden = dynamic_cast<TH2D*>(fin->Get(den_name.c_str()));
        if (!hnum) throw std::runtime_error("Missing histogram: " + num_name);
        if (!hden) throw std::runtime_error("Missing histogram: " + den_name);
        if (hden->GetEntries() == 0) throw std::runtime_error(den_name + " is empty");

        int nrow = 1, ncol = 1;
        DetermineSubplotGrid(static_cast<int>(q_eta_bins.size()), nrow, ncol);

        TCanvas c("cacc", "signal acceptance by eta", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        std::vector<TH1D*> all_ratios;
        std::vector<TLegend*> all_legs;
        all_legs.reserve(q_eta_bins.size());

        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            c.cd(static_cast<int>(ieta) + 1);
            gPad->SetLogx();
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.13);

            const auto& eta_bin = q_eta_bins.at(ieta);
            const auto yb = PairEtaPanels::Bins(hnum->GetYaxis(), eta_bin,
                                               "SignalAcceptancePlotter");
            const int y1 = yb.first, y2 = yb.second;

            const std::string tag = std::to_string(ieta);
            TH1D* hn = hnum->ProjectionX(("hn_acc_eta" + tag).c_str(), y1, y2, "e");
            TH1D* hd = hden->ProjectionX(("hd_acc_eta" + tag).c_str(), y1, y2, "e");
            hn->SetDirectory(nullptr);
            hd->SetDirectory(nullptr);

            TH1D* hratio = static_cast<TH1D*>(hn->Clone(("hratio_acc_eta" + tag).c_str()));
            hratio->SetDirectory(nullptr);
            hratio->Divide(hn, hd, 1., 1., "B");   // numerator is a subset of the denominator

            hratio->SetLineWidth(2);
            hratio->SetLineColor(kBlack);
            hratio->SetMarkerColor(kBlack);
            hratio->SetMarkerStyle(20);
            hratio->SetMarkerSize(0.9);
            hratio->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            hratio->GetYaxis()->SetTitle(y_title.c_str());
            hratio->GetXaxis()->SetTitleSize(0.06);
            hratio->GetYaxis()->SetTitleSize(0.06);
            hratio->GetXaxis()->SetLabelSize(0.05);
            hratio->GetYaxis()->SetLabelSize(0.05);
            hratio->GetYaxis()->SetTitleOffset(1.45);
            hratio->SetTitle("");
            hratio->SetMinimum(0.);
            hratio->SetMaximum(1.3);
            hratio->Draw("E1");

            TLegend* leg = new TLegend(0.18, 0.65, 0.62, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.045);
            leg->SetTextAlign(12);
            leg->SetMargin(0.06);
            leg->AddEntry((TObject*)0, mc_label_line1.c_str(), "");
            leg->AddEntry((TObject*)0, mc_label_line2.c_str(), "");
            leg->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eta_bin.first, eta_bin.second), "");
            leg->Draw();
            all_legs.push_back(leg);

            all_ratios.push_back(hratio);
            delete hn;
            delete hd;
        }

        const std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;

        for (TH1D* h : all_ratios) delete h;
        for (TLegend* l : all_legs) delete l;
    }

    // 1D acceptance vs truth eta^pair, one subplot per CANONICAL coarse pair-pT bin
    // (ParamsSet::pair_pt_coarse_bins). Projects the SAME num/denom 2D histograms as the other
    // views (CLAUDE.md Binnings rule 2); the coarse edges must be edges of the fine pT axis
    // (true for pT_bins_150, 2 fine bins per coarse bin) -- checked, never assumed.
    void DrawAcceptanceEtaByPt(const std::string& png_name,
                               const std::string& mc_label_line1,
                               const std::string& mc_label_line2) {
        const std::string num_name = NumName();
        const std::string den_name = DenName();
        TH2D* hnum = dynamic_cast<TH2D*>(fin->Get(num_name.c_str()));
        TH2D* hden = dynamic_cast<TH2D*>(fin->Get(den_name.c_str()));
        if (!hnum) throw std::runtime_error("Missing histogram: " + num_name);
        if (!hden) throw std::runtime_error("Missing histogram: " + den_name);
        if (hden->GetEntries() == 0) throw std::runtime_error(den_name + " is empty");

        const std::vector<double>& pt_edges = pms.pair_pt_coarse_bins;
        if (static_cast<int>(pt_edges.size()) != ParamsSet::N_COARSE_PAIR_PT_BINS + 1)
            throw std::runtime_error("pair_pt_coarse_bins size does not match N_COARSE_PAIR_PT_BINS");
        const TAxis* xax = hnum->GetXaxis();
        const int npanel = ParamsSet::N_COARSE_PAIR_PT_BINS;

        int nrow = 1, ncol = 1;
        DetermineSubplotGrid(npanel, nrow, ncol);
        TCanvas c("caccEta", "signal acceptance by pair pT", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        std::vector<TH1D*> all_ratios;
        std::vector<TLegend*> all_legs;
        for (int ipt = 0; ipt < npanel; ++ipt) {
            const double lo = pt_edges.at(ipt), hi = pt_edges.at(ipt + 1);
            const int x1 = xax->FindBin(lo + 1e-6 * lo);
            const int x2 = xax->FindBin(hi - 1e-6 * hi);
            if (std::fabs(xax->GetBinLowEdge(x1) - lo) > 1e-6 * lo || std::fabs(xax->GetBinUpEdge(x2) - hi) > 1e-6 * hi)
                throw std::runtime_error(Form("SignalAcceptancePlotter: coarse pair-pT edges [%g, %g] are not "
                                              "bin edges of %s's pT axis", lo, hi, num_name.c_str()));

            c.cd(ipt + 1);
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.13);

            const std::string tag = std::to_string(ipt);
            TH1D* hn = hnum->ProjectionY(("hn_accpt" + tag).c_str(), x1, x2, "e");
            TH1D* hd = hden->ProjectionY(("hd_accpt" + tag).c_str(), x1, x2, "e");
            hn->SetDirectory(nullptr);
            hd->SetDirectory(nullptr);
            TH1D* hratio = static_cast<TH1D*>(hn->Clone(("hratio_accpt" + tag).c_str()));
            hratio->SetDirectory(nullptr);
            hratio->Divide(hn, hd, 1., 1., "B");   // numerator is a subset of the denominator

            hratio->SetLineWidth(2);
            hratio->SetLineColor(kBlack);
            hratio->SetMarkerColor(kBlack);
            hratio->SetMarkerStyle(20);
            hratio->SetMarkerSize(0.9);
            hratio->GetXaxis()->SetTitle("#eta^{pair}");
            hratio->GetYaxis()->SetTitle(y_title.c_str());
            hratio->GetXaxis()->SetTitleSize(0.06);
            hratio->GetYaxis()->SetTitleSize(0.06);
            hratio->GetXaxis()->SetLabelSize(0.05);
            hratio->GetYaxis()->SetLabelSize(0.05);
            hratio->GetYaxis()->SetTitleOffset(1.45);
            // Measured range only: |eta^pair| >= pair_eta_fiducial_max is cut from num AND den.
            SetMeasuredEtaRange(hratio->GetXaxis());
            hratio->SetTitle("");
            hratio->SetMinimum(0.);
            hratio->SetMaximum(1.3);
            hratio->Draw("E1");

            TLegend* leg = new TLegend(0.18, 0.65, 0.62, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.045);
            leg->SetTextAlign(12);
            leg->SetMargin(0.06);
            leg->AddEntry((TObject*)0, mc_label_line1.c_str(), "");
            leg->AddEntry((TObject*)0, mc_label_line2.c_str(), "");
            leg->AddEntry((TObject*)0, Form("p_{T}^{pair} #in [%.3g, %.3g] GeV", lo, hi), "");
            leg->Draw();
            all_legs.push_back(leg);
            all_ratios.push_back(hratio);
            delete hn;
            delete hd;
        }

        const std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;
        for (TH1D* h : all_ratios) delete h;
        for (TLegend* l : all_legs) delete l;
    }

    virtual void Run() = 0;
};
