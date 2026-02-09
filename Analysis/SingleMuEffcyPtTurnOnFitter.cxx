#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cstdio>

#include "MuonObjectsParamsAndHelpers/PbPbBaseClass.h"
#include "MuonObjectsParamsAndHelpers/proj_range_to_suffix.cxx"

// ============================================================
// Base class: all shared logic
// ============================================================
class SingleMuEffcyPtTurnOnFitterBase {
protected:
    // ---- I/O that differs by dataset ----
    std::string data_dir;
    std::string fitting_outdir;
    std::string infile_name;
    std::string outfile_name;

    // ---- ROOT I/O ----
    TFile* f = nullptr;
    TFile* fout_pT_fit = nullptr;
    FILE* f_txtout = nullptr;

    // ---- reference histogram ----
    std::string h2d_ref_name;
    TH1D* h_pt_ref = nullptr; // reference histogram for setting the log axis correctly for the TGraph

    // ---- common maps ----
    std::map<std::string, std::string> trg_maps;
    std::map<std::string, std::string> musign_maps;

    std::vector<std::string> musigns = {"sign1", "sign2"};

    std::vector<std::string> q_eta_bins_for_pT_trg_effcy_graphs = {
        "minus2_40_TO_minus2_00", "minus2_00_TO_minus1_60", "minus1_60_TO_minus1_30", "minus0_90_TO_minus0_50",
        "minus0_50_TO_minus0_10", "0_10_TO_0_50", "0_50_TO_1_00", "1_30_TO_1_60", "1_60_TO_2_00", "2_00_TO_2_20"
    };

    static bool fileExists(const std::string& dir){
        return (gSystem->AccessPathName(dir.c_str()) == kFALSE);
    }
    static void makeDirIfNeeded(const std::string& dir){
        if (!fileExists(dir)) gSystem->mkdir(dir.c_str(), kTRUE /*recursive*/);
    }

public:
    bool debug_mode = false;

    enum FittingMode {
        erf_plus_log,
        fermi_plus_log,
        erf_plus_linear,
        fermi_plus_linear,
        erf
    };
    int fitting_mode = erf_plus_log;

    explicit SingleMuEffcyPtTurnOnFitterBase(bool debug_mode_input=false)
      : debug_mode(debug_mode_input) {}

    virtual ~SingleMuEffcyPtTurnOnFitterBase() {
        if (f_txtout) { fclose(f_txtout); f_txtout = nullptr; }
        if (fout_pT_fit) { fout_pT_fit->Close(); delete fout_pT_fit; fout_pT_fit = nullptr; }
        if (f) { f->Close(); delete f; f = nullptr; }
    }

    // ---- Run is shared except for centrality list + dataset naming ----
    void Run() {
        configureIO();     // dataset-specific (pp vs pbpb)
        initialize();      // shared

        // shared loops; centrality list comes from derived
        const auto ctrs = centralitySuffixes(); // default: {""} for pp

        for (const auto& [single_mu_trg_name, /*dimu_trg_name*/ _] : trg_maps) {
            for (const auto& musign : musigns) {
                for (const auto& ctr : ctrs) {
                    drawOneCanvas(single_mu_trg_name, musign, ctr);
                }
            }
        }

        // close happens in destructor too, but keeping your explicit style is fine:
        if (f_txtout) { fclose(f_txtout); f_txtout = nullptr; }
        if (f) { f->Close(); }
        if (fout_pT_fit) { fout_pT_fit->Close(); }
    }

protected:
    // dataset-specific knobs:
    virtual void configureIO() = 0;  // must set data_dir + infile_name
    virtual std::vector<std::string> centralitySuffixes() const { return {""}; } // pp default

    // =========================================================================
    // Ensure a TGraph drawn on a log-x pad spans the full histogram range
    // =========================================================================
    void adjustLogXRange(TGraphAsymmErrors* g, const TH1* hRef) const
    {
        if (!g || !hRef) return;
        double xmin = hRef->GetXaxis()->GetBinLowEdge(1);
        double xmax = hRef->GetXaxis()->GetBinLowEdge(hRef->GetNbinsX() + 1);
        g->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    virtual void initialize() {
        // fitting outdir (same as your two versions) :contentReference[oaicite:5]{index=5} :contentReference[oaicite:6]{index=6}
        switch (fitting_mode){
        case erf_plus_log:       fitting_outdir = "trg_effcy_pT_fitting_to_erf_plus_log/"; break;
        case fermi_plus_log:     fitting_outdir = "trg_effcy_pT_fitting_to_fermi_plus_log/"; break;
        case erf:                fitting_outdir = "trg_effcy_pT_fitting_to_erf/"; break;
        case erf_plus_linear:    fitting_outdir = "trg_effcy_pT_fitting_to_erf_plus_linear/"; break;
        case fermi_plus_linear:  fitting_outdir = "trg_effcy_pT_fitting_to_fermi_plus_linear/"; break;
        default:
            std::cerr << "Fitting mode invalid; defaulting to erf_plus_log\n";
            fitting_outdir = "trg_effcy_pT_fitting_to_erf_plus_log/";
        }

        // outfile always under data_dir + fitting_outdir (same as your code) :contentReference[oaicite:7]{index=7} :contentReference[oaicite:8]{index=8}
        outfile_name = data_dir + fitting_outdir + "single_mu_effcy_pT_fit.root";
        makeDirIfNeeded(data_dir + fitting_outdir);

        // text output
        f_txtout = fopen((data_dir + fitting_outdir + "fit_results.txt").c_str(), "w");
        if (!f_txtout) perror("Error opening fit_results.txt");

        // trigger maps (same as your code) :contentReference[oaicite:9]{index=9} :contentReference[oaicite:10]{index=10}
        trg_maps.clear();
        trg_maps["mu4noL1"] = "mu4_mu4noL1";
        trg_maps["mu4"]     = "2mu4";
        // trg_maps["mu4_AND_mu4noL1"]     = "2mu4_AND_mu4_mu4noL1";

        musign_maps.clear();
        musign_maps["sign1"] = "mu+";
        musign_maps["sign2"] = "mu-";

        // open input/output ROOT files
        f = TFile::Open(infile_name.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Error: cannot open file " << infile_name << "\n";
            return;
        }

        fout_pT_fit = TFile::Open(outfile_name.c_str(), "recreate");
        if (!fout_pT_fit || fout_pT_fit->IsZombie()) {
            std::cerr << "Error: cannot open file " << outfile_name << "\n";
            return;
        }

        initializeHistRef();
    }

    virtual void initializeHistRef() {
        TH2D* h2d_ref = dynamic_cast<TH2D*>(f->Get(h2d_ref_name.c_str()));
        if (!h2d_ref) {
            std::cerr << "Warning: reference TH2D histogram " << h2d_ref_name << " not found in file " << infile_name << std::endl;
            return;
        }

        h_pt_ref = h2d_ref->ProjectionY(Form("%s_py", h2d_ref->GetName()));
    }

private:
    TF1* fitTurnOnErfPlusLog(TGraphAsymmErrors* g){
        std::string fname = g->GetName() ? g->GetName() : "graph";
        if (fname.rfind("g_", 0) == 0) fname.replace(0, 2, "f_"); else fname = "f_" + fname;

        double pT_min = 4;
        double pT_max = 60;

        TF1* fTurnOn = new TF1(fname.c_str(),
                               [this, pT_min](double* x, double* p){
                                    double pT = x[0];
                                    double mean = p[0];
                                    double sigma = p[1];
                                    double norm = p[2];
                                    double corrCoef = p[3];

                                    if (fitting_mode == erf_plus_log)
                                        return norm * 0.5 * (1.0 + TMath::Erf((pT - mean)/(TMath::Sqrt2()*sigma)))
                                               * (1 + corrCoef * TMath::Log(1 + (pT - pT_min)/pT_min));
                                    else
                                        return norm * 0.5 * (1.0 + TMath::Erf((pT - mean)/(TMath::Sqrt2()*sigma)))
                                               * (1 + corrCoef * (pT - pT_min)/pT_min);
                               },
                               pT_min, pT_max, 4);

        fTurnOn->SetParNames("mean", "sigma", "plateau", "corrCoef");
        fTurnOn->SetParameters(4.0, 2.0, 1.0, 0.);
        fTurnOn->SetParLimits(0, 0, 10);
        fTurnOn->SetParLimits(1, 0.01, 50);
        fTurnOn->SetParLimits(2, 0.5, 1.2);
        fTurnOn->SetParLimits(3, 0, (fitting_mode == erf_plus_log ? 0.1 : 0.2));

        g->Fit(fTurnOn, "QR");

        if (f_txtout){
            fprintf(f_txtout, "=== Fit results ===\n");
            fprintf(f_txtout, "Mean   = %.3f\n", fTurnOn->GetParameter(0));
            fprintf(f_txtout, "Sigma  = %.3f\n", fTurnOn->GetParameter(1));
            fprintf(f_txtout, "Plateau= %.3f\n", fTurnOn->GetParameter(2));
            fprintf(f_txtout, "Correction Coefficient= %.3f\n", fTurnOn->GetParameter(3));
        }
        return fTurnOn;
    }

    TF1* fitTurnOnFermiPlusLog(TGraphAsymmErrors* g){
        std::string fname = g->GetName() ? g->GetName() : "graph";
        if (fname.rfind("g_", 0) == 0) fname.replace(0, 2, "f_"); else fname = "f_" + fname;

        double pT_min = 4;
        double pT_max = 60;

        TF1* fTurnOn = new TF1(fname.c_str(),
                               [this, pT_min](double* x, double* p){
                                    double pT = x[0];
                                    double normFermi = p[0];
                                    double pT0 = p[1];
                                    double Delta = p[2];
                                    double corrCoef = p[3];

                                    if (fitting_mode == fermi_plus_log)
                                        return normFermi / (1 + TMath::Exp((pT0 - pT)/ Delta))
                                               * (1 + corrCoef * TMath::Log(1 + (pT - pT_min)/pT_min));
                                    else
                                        return normFermi / (1 + TMath::Exp((pT0 - pT)/ Delta))
                                               * (1 + corrCoef * ((pT - pT_min)/pT_min));
                               },
                               pT_min, pT_max, 4);

        fTurnOn->SetParNames("normFermi", "pT0", "Delta", "corrCoef");
        fTurnOn->SetParameters(0.9, 4.0, 1.5, 0.);
        fTurnOn->SetParLimits(0, 0.6, 1.0);
        fTurnOn->SetParLimits(1, 2.5, 5.5);
        fTurnOn->SetParLimits(2, 0, 10);
        fTurnOn->SetParLimits(3, 0, (fitting_mode == fermi_plus_log ? 0.1 : 0.2));

        g->Fit(fTurnOn, "QR");

        if (f_txtout){
            fprintf(f_txtout, "=== Fit results ===\n");
            fprintf(f_txtout, "normFermi   = %.3f\n", fTurnOn->GetParameter(0));
            fprintf(f_txtout, "pT0  = %.3f\n", fTurnOn->GetParameter(1));
            fprintf(f_txtout, "Delta  = %.3f\n", fTurnOn->GetParameter(2));
            fprintf(f_txtout, "Correction Coefficient  = %.3f\n", fTurnOn->GetParameter(3));
        }
        return fTurnOn;
    }

    TF1* fitTurnOnErf(TGraphAsymmErrors* g){
        std::string fname = g->GetName() ? g->GetName() : "graph";
        if (fname.rfind("g_", 0) == 0) fname.replace(0, 2, "f_"); else fname = "f_" + fname;

        double pT_min = 4;
        double pT_max = 60;

        TF1* fTurnOn = new TF1(fname.c_str(),
                               [](double* x, double* p){
                                   double pT = x[0];
                                   double mean = p[0];
                                   double sigma = p[1];
                                   double norm = p[2];
                                   return norm * 0.5 * (1.0 + TMath::Erf((pT - mean)/(TMath::Sqrt2()*sigma)));
                               },
                               pT_min - 0.1, pT_max, 3);

        fTurnOn->SetParNames("mean", "sigma", "plateau");
        fTurnOn->SetParameters(4.0, 2.0, 1.0);
        fTurnOn->SetParLimits(0, 0, 10);
        fTurnOn->SetParLimits(1, 0.01, 50);
        fTurnOn->SetParLimits(2, 0.5, 1.2);

        g->Fit(fTurnOn, "QR");

        if (f_txtout){
            fprintf(f_txtout, "=== Fit results ===\n");
            fprintf(f_txtout, "Mean   = %.3f\n", fTurnOn->GetParameter(0));
            fprintf(f_txtout, "Sigma  = %.3f\n", fTurnOn->GetParameter(1));
            fprintf(f_txtout, "Plateau= %.3f\n", fTurnOn->GetParameter(2));
        }
        return fTurnOn;
    }

    TF1* fitTurnOn(TGraphAsymmErrors* g){
        switch(fitting_mode){
        case erf_plus_log:      return fitTurnOnErfPlusLog(g);
        case fermi_plus_log:    return fitTurnOnFermiPlusLog(g);
        case erf:               return fitTurnOnErf(g);
        case erf_plus_linear:   return fitTurnOnErfPlusLog(g);
        case fermi_plus_linear: return fitTurnOnFermiPlusLog(g);
        default:
            std::cerr << "Fitting mode invalid; defaulting to erf_plus_log\n";
            return fitTurnOnErfPlusLog(g);
        }
    }

    // shared canvas drawing; ctr is "" for pp, and "_ctrX_Y" for PbPb :contentReference[oaicite:11]{index=11} :contentReference[oaicite:12]{index=12}
    void drawOneCanvas(const std::string& trg, const std::string& musign, const std::string& ctr)
    {
        if (debug_mode)
            std::cout << "drawOneCanvas trg=" << trg << " musign=" << musign_maps[musign] << " ctr=" << ctr << "\n";

        gStyle->SetOptStat(0);

        std::string cname = ctr.empty()
            ? Form("c_%s_%s", trg.c_str(), musign.c_str())
            : Form("c_%s%s_%s", trg.c_str(), ctr.c_str(), musign.c_str());

        TCanvas* c = new TCanvas(cname.c_str(), "Trigger Turn-on Curves", 1500, 1000);
        c->Divide(4,3);
        c->SetGrid();

        int idx = 0;
        for (const auto& q_eta_bin : q_eta_bins_for_pT_trg_effcy_graphs){
            c->cd(idx + 1);
            gPad->SetLogx();

            TLegend* leg = new TLegend(0.35, 0.25, 0.88, 0.5);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            // graph name differs only by ctr insertion (exactly like your two files) :contentReference[oaicite:13]{index=13} :contentReference[oaicite:14]{index=14}
            std::string gname;
            if (ctr.empty()){
                gname = "g_pt2nd_vs_q_eta2nd_" + musign + "_" + trg_maps.at(trg) + "_sepr_py_" + q_eta_bin + "_divided";
            } else {
                gname = "g_pt2nd_vs_q_eta2nd";
                gname += ctr + "_" + musign + "_" + trg_maps.at(trg) + "_sepr_py_" + q_eta_bin + "_divided";
            }

            TGraphAsymmErrors* g = dynamic_cast<TGraphAsymmErrors*>(f->Get(gname.c_str()));
            if (!g) {
                std::cerr << "Warning: graph " << gname << " not found in file.\n";
                ++idx;
                continue;
            }

            TF1* fit = fitTurnOn(g);
            if (!fit){
                std::cerr << "Warning: fit for graph " << gname << " is null\n";
                ++idx;
                continue;
            }

            // write fit function into output root file
            fit->Write();

            g->SetMarkerColor(kBlack);
            g->SetLineColor(kBlack);
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.9);

            fit->SetLineColor(kRed);
            fit->SetLineWidth(2);

            g->GetXaxis()->SetTitle("p_{T} [GeV]");
            g->GetYaxis()->SetTitle("#epsilon");
            g->GetYaxis()->SetRangeUser(0, 1.1);
            g->Draw("AP");
            if (h_pt_ref) adjustLogXRange(g, h_pt_ref);
            fit->Draw("SAME");

            std::string musign_label = (musign == "sign1")? "#mu^{+}" : "#mu^{-}";
            auto q_eta_pair = hProjNameToPair(gname);
            std::string q_eta_label = pairToLegendLabel(q_eta_pair);

            leg->AddEntry(g, (trg + ", " + musign_label).c_str(), "lp");
            leg->AddEntry("", q_eta_label.c_str(), "");
            leg->Draw("SAME");

            ++idx;
        }

        // filename differs only by ctr insertion (exactly like your two files) :contentReference[oaicite:15]{index=15} :contentReference[oaicite:16]{index=16}
        if (ctr.empty()){
            c->SaveAs(Form("%s%strg_effcy_pT_fitting_%s_%s.png",
                           data_dir.c_str(), fitting_outdir.c_str(),
                           trg_maps.at(trg).c_str(), musign_maps.at(musign).c_str()));
        } else {
            c->SaveAs(Form("%s%strg_effcy_pT_fitting_%s%s_%s.png",
                           data_dir.c_str(), fitting_outdir.c_str(),
                           trg_maps.at(trg).c_str(), ctr.c_str(), musign_maps.at(musign).c_str()));
        }
    }
};

// ============================================================
// Derived: PP
// ============================================================
class SingleMuEffcyPtTurnOnFitterPP : public SingleMuEffcyPtTurnOnFitterBase {
public:
    using SingleMuEffcyPtTurnOnFitterBase::SingleMuEffcyPtTurnOnFitterBase;

protected:
    void configureIO() override {
        data_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_2024/";
        infile_name = data_dir + "histograms_real_pairs_pp_2024_single_mu4.root";
        h2d_ref_name = "h_pt2nd_vs_q_eta2nd_sign1_mu4";
    }
};

// ============================================================
// Derived: PbPb
// ============================================================
class SingleMuEffcyPtTurnOnFitterPbPb : public SingleMuEffcyPtTurnOnFitterBase, public PbPbBaseClass {
public:
    SingleMuEffcyPtTurnOnFitterPbPb(
                                    const std::string ctr_binning_version_input = "default",
                                    bool debug_mode_input=false)
                    : SingleMuEffcyPtTurnOnFitterBase(debug_mode_input){
                        ctr_binning_version = ctr_binning_version_input;
                    }


protected:
    void initialize() override {
        PbPbBaseClass::InitializePbPb();
        SingleMuEffcyPtTurnOnFitterBase::initialize();
    }

    void configureIO() override {
        data_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pbpb_2023/";
        infile_name = data_dir + "histograms_real_pairs_pbpb_2023_single_mu4.root";
        h2d_ref_name = "h_pt2nd_vs_q_eta2nd_ctr0_5_sign1_mu4";
    }

    std::vector<std::string> centralitySuffixes() const override {
        return ctr_bins;
    }
};

// ------------------------------------------------------------
// ROOT entry points (like your originals)
// ------------------------------------------------------------
void single_muon_trig_effcy_pT_fitting() {
    auto* fitter = new SingleMuEffcyPtTurnOnFitterPP();
    fitter->fitting_mode = SingleMuEffcyPtTurnOnFitterBase::erf_plus_log;
    fitter->Run();
    delete fitter;
}

void single_muon_trig_effcy_pT_fitting_PbPb() {
    auto* fitter = new SingleMuEffcyPtTurnOnFitterPbPb();
    fitter->fitting_mode = SingleMuEffcyPtTurnOnFitterBase::fermi_plus_log;
    fitter->Run();
    delete fitter;
}
