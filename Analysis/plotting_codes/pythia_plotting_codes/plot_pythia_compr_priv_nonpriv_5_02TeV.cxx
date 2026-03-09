#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TAxis.h>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../helper_functions.c"

namespace {

std::string trim_title_before_bracket(const std::string& s) {
    size_t pos1 = s.find(" [");
    size_t pos2 = s.find("[");
    size_t pos = std::string::npos;

    if (pos1 != std::string::npos) pos = pos1;
    if (pos2 != std::string::npos) pos = (pos == std::string::npos ? pos2 : std::min(pos, pos2));
    if (pos == std::string::npos) return s;
    return s.substr(0, pos);
}

std::string make_default_ytitle_from_xtitle(const std::string& xtitle) {
    const std::string prefix = "FULL: ";
    if (xtitle.rfind(prefix, 0) == 0) {
        return xtitle.substr(prefix.size());
    }
    return "d#sigma/d" + trim_title_before_bracket(xtitle) + " [nb]";
}

TFile* open_first_existing(const std::vector<std::string>& candidates, std::string& chosen_path) {
    for (const auto& p : candidates) {
        TFile* f = TFile::Open(p.c_str(), "READ");
        if (f && !f->IsZombie()) {
            chosen_path = p;
            return f;
        }
        if (f) delete f;
    }
    return nullptr;
}

TH1D* clone_hist_from_file(TFile* f, const std::string& hname, const std::string& new_name) {
    if (!f) return nullptr;
    TH1* h = dynamic_cast<TH1*>(f->Get(hname.c_str()));
    if (!h) return nullptr;
    TH1D* out = dynamic_cast<TH1D*>(h->Clone(new_name.c_str()));
    if (!out) return nullptr;
    out->SetDirectory(nullptr);
    return out;
}

TH1D* clone_proj_from_2d(TFile* f,
                         const std::string& h2name,
                         bool project_x,
                         const std::string& new_name) {
    if (!f) return nullptr;
    TH2D* h2 = dynamic_cast<TH2D*>(f->Get(h2name.c_str()));
    if (!h2) return nullptr;

    TH1D* h1 = project_x
        ? dynamic_cast<TH1D*>(h2->ProjectionX(new_name.c_str()))
        : dynamic_cast<TH1D*>(h2->ProjectionY(new_name.c_str()));

    if (!h1) return nullptr;
    h1->SetDirectory(nullptr);
    return h1;
}

TH1D* fetch_hist_with_fallback(TFile* f,
                               const std::string& var1d,
                               const std::string& sign_suffix,
                               const std::string& sample_tag) {
    const std::string h1_name = "h_" + var1d + sign_suffix;
    if (TH1D* h1 = clone_hist_from_file(f, h1_name, h1_name + "_" + sample_tag)) {
        return h1;
    }

    // Fallbacks used by the categorized plotting flow when 1D is not directly present.
    if (var1d == "pair_pt") {
        return clone_proj_from_2d(f,
                                  "h_pt_lead_vs_pair_pt" + sign_suffix,
                                  true,
                                  "h_" + var1d + sign_suffix + "_" + sample_tag + "_projx");
    }

    if (var1d == "minv_zoomin") {
        return clone_proj_from_2d(f,
                                  "h_minv_zoomin_vs_pair_pt" + sign_suffix,
                                  false,
                                  "h_" + var1d + sign_suffix + "_" + sample_tag + "_projy");
    }

    if (var1d == "Deta_zoomin") {
        TH1D* h = clone_proj_from_2d(f,
                                     "h_Deta_vs_Dphi" + sign_suffix,
                                     false,
                                     "h_" + var1d + sign_suffix + "_" + sample_tag + "_projy");
        if (h) return h;
        return clone_proj_from_2d(f,
                                  "h_Deta_zoomin_vs_Dphi_zoomin" + sign_suffix,
                                  false,
                                  "h_" + var1d + sign_suffix + "_" + sample_tag + "_projy_zoom");
    }

    if (var1d == "Dphi_zoomin") {
        TH1D* h = clone_proj_from_2d(f,
                                     "h_Deta_vs_Dphi" + sign_suffix,
                                     true,
                                     "h_" + var1d + sign_suffix + "_" + sample_tag + "_projx");
        if (h) return h;
        return clone_proj_from_2d(f,
                                  "h_Deta_zoomin_vs_Dphi_zoomin" + sign_suffix,
                                  true,
                                  "h_" + var1d + sign_suffix + "_" + sample_tag + "_projx_zoom");
    }

    return nullptr;
}

TH1D* fetch_hist_or_throw(TFile* f,
                          const std::string& var1d,
                          const std::string& sign_suffix,
                          const std::string& sample_tag,
                          const std::string& file_path) {
    TH1D* h = fetch_hist_with_fallback(f, var1d, sign_suffix, sample_tag);
    if (h) return h;

    std::ostringstream oss;
    oss << "Missing histogram for var='" << var1d
        << "', sign='" << sign_suffix
        << "', sample='" << sample_tag
        << "' in file: " << file_path
        << ". Tried 1D key h_" << var1d << sign_suffix
        << " and 2D projection fallbacks where applicable.";
    throw std::runtime_error(oss.str());
}

void style_hist(TH1D* h, Color_t color, float norm) {
    if (!h) return;
    hist_helper(h, norm, false, "");
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->SetStats(0);
}

} // namespace

void plot_pythia_compr_priv_nonpriv_5p02TeV() {
    gROOT->SetBatch(kTRUE);

    const std::vector<std::string> vars = {
        "DR", "DR_zoomin", "Deta_zoomin", "Dphi_zoomin", "pair_pt", "minv_zoomin"
    };

    const std::map<int, std::string> sign_suffix = {
        {0, "_sign1"},
        {1, "_sign2"}
    };

    const std::map<int, std::string> panel_title = {
        {0, "same-sign"},
        {1, "opposite-sign"}
    };

    std::string private_path;
    std::string nonprivate_path;

    std::vector<std::string> private_candidates = {
        "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/histograms_pythia_combined_no_data_resonance_cuts.root",
        "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/histograms_pythia_5p36TeV_no_data_resonance_cuts.root",
        "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/histograms_pythia_5p36TeV_no_data_resonance_cuts.root"
    };

    std::vector<std::string> nonprivate_candidates = {
        "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5TeV/histograms_pythia_5p02TeV_no_data_resonance_cuts.root",
        "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5TeV/histograms_pythia_5p36TeV_no_data_resonance_cuts.root"
    };

    std::unique_ptr<TFile> f_priv(open_first_existing(private_candidates, private_path));
    std::unique_ptr<TFile> f_nonpriv(open_first_existing(nonprivate_candidates, nonprivate_path));

    if (!f_priv) {
        std::ostringstream oss;
        oss << "Invalid private input file. None of these files can be opened:";
        for (const auto& p : private_candidates) oss << "\n  - " << p;
        throw std::runtime_error(oss.str());
    }
    if (!f_nonpriv) {
        std::ostringstream oss;
        oss << "Invalid non-private input file. None of these files can be opened:";
        for (const auto& p : nonprivate_candidates) oss << "\n  - " << p;
        throw std::runtime_error(oss.str());
    }

    std::cout << "[plot_pythia_compr_priv_nonpriv_5p02TeV] private file: " << private_path << std::endl;
    std::cout << "[plot_pythia_compr_priv_nonpriv_5p02TeV] non-private file: " << nonprivate_path << std::endl;

    const std::string outdir = "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5TeV/plots";
    gSystem->mkdir(outdir.c_str(), kTRUE);

    for (const auto& var : vars) {
        std::unique_ptr<TCanvas> c(new TCanvas(("c_" + var).c_str(), ("c_" + var).c_str(), 1500, 600));
        c->Divide(2, 1);

        // Keep histograms alive until after SaveAs; otherwise pads can appear empty.
        std::vector<std::unique_ptr<TH1D>> keep_alive;
        keep_alive.reserve(4);

        for (int isign = 0; isign < 2; ++isign) {
            c->cd(isign + 1);
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.14);

            const std::string ss = sign_suffix.at(isign);

            TH1D* h_priv_raw = fetch_hist_or_throw(f_priv.get(), var, ss, "priv", private_path);
            TH1D* h_nonpriv_raw = fetch_hist_or_throw(f_nonpriv.get(), var, ss, "nonpriv", nonprivate_path);

            keep_alive.emplace_back(h_priv_raw);
            keep_alive.emplace_back(h_nonpriv_raw);

            TH1D* h_priv = keep_alive[keep_alive.size() - 2].get();
            TH1D* h_nonpriv = keep_alive[keep_alive.size() - 1].get();

            style_hist(h_priv, kBlack, 1.0e6f);
            style_hist(h_nonpriv, kRed, 1000.0f);

            std::string xtitle = h_priv->GetXaxis()->GetTitle();
            if (xtitle.empty()) xtitle = h_nonpriv->GetXaxis()->GetTitle();
            const std::string ytitle = make_default_ytitle_from_xtitle(xtitle);

            h_priv->GetXaxis()->SetTitle(xtitle.c_str());
            h_nonpriv->GetXaxis()->SetTitle(xtitle.c_str());
            h_priv->GetYaxis()->SetTitle(ytitle.c_str());
            h_nonpriv->GetYaxis()->SetTitle(ytitle.c_str());

            const double ymax = std::max(h_priv->GetMaximum(), h_nonpriv->GetMaximum());
            h_priv->GetYaxis()->SetRangeUser(0.0, ymax > 0.0 ? ymax * 1.25 : 1.0);

            h_priv->Draw("hist e");
            h_nonpriv->Draw("hist e same");

            TLegend* leg = new TLegend(0.5, 0.69, 0.89, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(43);
            leg->SetTextSize(20);
            leg->AddEntry(h_priv, "private", "l");
            leg->AddEntry((TObject*)nullptr, "nCTEQ15npFullNuc_208_82", "");
            leg->AddEntry(h_nonpriv, "non-private (5.02 TeV)", "l");
            leg->AddEntry((TObject*)nullptr, "nNNPDF30_nlo_as_0118_A208_Z82", "");
            leg->Draw();

            h_priv->SetTitle((var + " " + panel_title.at(isign)).c_str());
        }

        const std::string outname = outdir + "/pythia_" + var + "_priv_nonpriv_5.02TeV_compr.png";
        c->SaveAs(outname.c_str());
    }
}

void plot_pythia_compr_priv_nonpriv_5_02TeV() {
    plot_pythia_compr_priv_nonpriv_5p02TeV();
}
