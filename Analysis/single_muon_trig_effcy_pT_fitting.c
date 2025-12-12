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
#include "MuonObjectsParamsAndHelpers/proj_range_to_suffix.cxx"
#include <fstream>
#include <iostream>
#include <cstdio>  // for printf, fprintf, fopen, fclose


class SingleMuEffcyPtTurnOnFitter {
private:
    std::string data_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pp_2024/";
    std::string fitting_outdir;
    std::string infile_name;
    std::string outfile_name;
        
    TFile* f = nullptr;
    TFile* fout_pT_fit = nullptr;

    FILE* f_txtout;

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

    void initialize(){
        switch (fitting_mode){
        case erf_plus_log:
            fitting_outdir = "trg_effcy_pT_fitting_to_erf_plus_log/";
            break;
        case fermi_plus_log:
            fitting_outdir = "trg_effcy_pT_fitting_to_fermi_plus_log/";
            break;
        case erf:
            fitting_outdir = "trg_effcy_pT_fitting_to_erf/";
            break;
        case erf_plus_linear:
            fitting_outdir = "trg_effcy_pT_fitting_to_erf_plus_linear/";
            break;
        case fermi_plus_linear:
            fitting_outdir = "trg_effcy_pT_fitting_to_fermi_plus_linear/";
            break;
        default:
            std::cerr << "Fitting mode must be between 1 and 3! Set to default value 1: erf_plus_log" << std::endl;
            fitting_outdir = "trg_effcy_pT_fitting_to_fermi_plus_log/";
        }

        infile_name = data_dir + "histograms_real_pairs_pp_2024_single_mu4_new_RDF.root";
        outfile_name = data_dir + fitting_outdir + "single_mu_effcy_pT_fit.root";
        
        makeDirIfNeeded(data_dir + fitting_outdir);

        f_txtout = fopen((data_dir + fitting_outdir + "fit_results.txt").c_str(), "w");  // open file for writing
        if (!f_txtout) {
            perror("Error opening file");
        }

        trg_maps["mu4noL1"] = "mu4_mu4noL1";
        trg_maps["mu4"] = "2mu4";

        musign_maps["sign1"] = "mu+";
        musign_maps["sign2"] = "mu-";

        f = TFile::Open(infile_name.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Error: cannot open file " << infile_name << std::endl;
            return;
        }

        fout_pT_fit = TFile::Open(outfile_name.c_str(), "recreate");
        if (!fout_pT_fit || fout_pT_fit->IsZombie()) {
            std::cerr << "Error: cannot open file " << outfile_name << std::endl;
            return;
        }
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
    
    // ---- Constructor ----
    SingleMuEffcyPtTurnOnFitter(bool debug_mode_input = false) : debug_mode(debug_mode_input) {}

    // ---- Destructor ----
    ~SingleMuEffcyPtTurnOnFitter() {
        if (f && f->IsOpen()) f->Close();
    }

    // ---- Run fitting & plotting ----
    void Run() {
        initialize();

        for (const auto& [single_mu_trg_name, dimu_trg_name] : trg_maps) {
            for (auto musign : musigns){
                drawOneCanvas(single_mu_trg_name, musign);
            }
        }

        fclose(f_txtout);  // always close
        f->Close();
        fout_pT_fit->Close();
    }

private:
    // ---- Helper: fit one graph ----
    TF1* fitTurnOnErfPlusLog(TGraphAsymmErrors* g){
        std::string fname = g->GetName() ? g->GetName() : "graph";
        if (fname.rfind("g_", 0) == 0) fname.replace(0, 2, "f_"); else fname = "f_" + fname;

        double pT_min = 4;
        double pT_max = 60;

        TF1* fTurnOn = new TF1(fname.c_str(),
                               [&](double* x, double* p){
                                    double pT_min = 4;
                                    double pT = x[0];
                                    double mean = p[0];
                                    double sigma = p[1];
                                    double norm = p[2]; // max plateau (≤1)
                                    double corrCoef = p[3]; // max plateau (≤1)

                                    const double eps = 1e-6;

                                    if (fitting_mode == erf_plus_log) // log
                                        return norm * 0.5 * (1.0 + TMath::Erf((pT - mean)/(TMath::Sqrt2()*sigma))) * (1 + corrCoef * TMath::Log(1 + (pT - pT_min)/pT_min));
                                        // return norm * 0.5 * (1.0 + TMath::Erf((pT - mean)/(TMath::Sqrt2()*sigma))) * (1 + corrCoef * TMath::Log(TMath::Max((pT - pT_min)/pT_min, eps)));
                                    else // linear
                                        return norm * 0.5 * (1.0 + TMath::Erf((pT - mean)/(TMath::Sqrt2()*sigma))) * (1 + corrCoef * (pT - pT_min)/pT_min);
                               },
                               pT_min, pT_max, 4); // range & #parameters

        fTurnOn->SetParNames("mean", "sigma", "plateau", "corrCoef");
        fTurnOn->SetParameters(4.0, 2.0, 1.0, 0.); // initial guesses

        // Reasonable limits to help convergence
        fTurnOn->SetParLimits(0, 0, 10);
        fTurnOn->SetParLimits(1, 0.01, 50);
        fTurnOn->SetParLimits(2, 0.5, 1.2);
        if (fitting_mode == erf_plus_log)   fTurnOn->SetParLimits(3, 0, 0.1); // corrCoef limit for log correction
        else                                fTurnOn->SetParLimits(3, 0, 0.2); // corrCoef limit for linear correction

        g->Fit(fTurnOn, "QR"); // "R" = respect function range

        // Print fit parameters
        if (f_txtout){
            fprintf(f_txtout, "=== Fit results ===\n");
            fprintf(f_txtout, "Mean   = %.3f\n", fTurnOn->GetParameter(0));
            fprintf(f_txtout, "Sigma  = %.3f\n", fTurnOn->GetParameter(1));
            fprintf(f_txtout, "Plateau= %.3f\n", fTurnOn->GetParameter(2));
            fprintf(f_txtout, "Correction Coefficient= %.3f\n",   fTurnOn->GetParameter(3));
        }

        return fTurnOn;
    }

    TF1* fitTurnOnFermiPlusLog(TGraphAsymmErrors* g){
        std::string fname = g->GetName() ? g->GetName() : "graph";
        if (fname.rfind("g_", 0) == 0) fname.replace(0, 2, "f_"); else fname = "f_" + fname;

        double pT_min = 4;
        double pT_max = 60;

        TF1* fTurnOn = new TF1(fname.c_str(),
                               [&](double* x, double* p){
                                    double pT_min = 4;
                                    double pT = x[0];
                                    double normFermi = p[0];
                                    double pT0 = p[1];
                                    double Delta = p[2]; // max plateau (≤1)
                                    double corrCoef = p[3]; // max plateau (≤1)
                                     
                                    const double eps = 1e-6;
    
                                    if (fitting_mode == fermi_plus_log) // log
                                        return normFermi / (1 + TMath::Exp((pT0 - pT)/ Delta)) * (1 + corrCoef * TMath::Log(1 + (pT - pT_min)/pT_min));
                                        // return normFermi / (1 + TMath::Exp((pT0 - pT)/ Delta)) * (1 + corrCoef * TMath::Log(TMath::Max((pT - pT_min)/pT_min, eps)));
                                    else // linear
                                        return normFermi / (1 + TMath::Exp((pT0 - pT)/ Delta)) * (1 + corrCoef * ((pT - pT_min)/pT_min));
                               },
                               pT_min, pT_max, 4); // range & #parameters

        fTurnOn->SetParNames("normFermi", "pT0", "Delta", "corrCoef");
        fTurnOn->SetParameters(0.9, 4.0, 1.5, 0.); // initial guesses

        // Reasonable limits to help convergence
        fTurnOn->SetParLimits(0, 0.6, 1); // norm Fermi: pT --> infty, no correction: Fermi --> normFermi
        fTurnOn->SetParLimits(1, 2.5, 5.5); // pT0 (turning point) < 6GeV
        fTurnOn->SetParLimits(2, 0, 10); // Delta > 0
        if (fitting_mode == fermi_plus_log) fTurnOn->SetParLimits(3, 0, 0.1); // corrCoef limit for log correction
        else                                fTurnOn->SetParLimits(3, 0, 0.2); // corrCoef limit for linear correction

        g->Fit(fTurnOn, "QR"); // "R" = respect function range

        // Print fit parameters
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
                                   double norm = p[2]; // max plateau (≤1)
                                   return norm * 0.5 * (1.0 + TMath::Erf((pT - mean)/(TMath::Sqrt2()*sigma)));
                               },
                               pT_min - 0.1, pT_max, 3); // range & #parameters

        fTurnOn->SetParNames("mean", "sigma", "plateau");
        fTurnOn->SetParameters(4.0, 2.0, 1.0); // initial guesses

        // Reasonable limits to help convergence
        fTurnOn->SetParLimits(0, 0, 10);
        fTurnOn->SetParLimits(1, 0.01, 50);
        fTurnOn->SetParLimits(2, 0.5, 1.2);

        g->Fit(fTurnOn, "QR"); // "R" = respect function range

        // Print fit parameters
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
        case erf_plus_log:
            return fitTurnOnErfPlusLog(g);
        case fermi_plus_log:
            return fitTurnOnFermiPlusLog(g);
        case erf:
            return fitTurnOnErf(g);
        case erf_plus_linear:
            return fitTurnOnErfPlusLog(g);
        case fermi_plus_linear:
            return fitTurnOnFermiPlusLog(g);
        default:
            std::cerr << "Fitting mode invalid: must be between 1 and 3! Set to default value 1: erf_plus_log" << std::endl;
            return fitTurnOnErfPlusLog(g);
        }

        return nullptr;
    }

    // ---- Helper: draw all on one canvas ----
    void drawOneCanvas(std::string trg, std::string musign) // draw one canvas for one (trigger, muon sign) pair
    {
        if (debug_mode) std::cout << "Calling function drawOneCanvas for trigger " << trg << ", muon sign " << musign_maps[musign] << std::endl;
        gStyle->SetOptStat(0);

        TCanvas* c = new TCanvas(Form("c_%s_%s", trg.c_str(), musign.c_str()), "Trigger Turn-on Curves", 1500, 1000);
        c->Divide(4,3);
        c->SetGrid();


        int idx = 0;
        for (auto q_eta_bin : q_eta_bins_for_pT_trg_effcy_graphs){
            c->cd(idx + 1);
            gPad->SetLogx();

            TLegend* leg = new TLegend(0.35, 0.25, 0.88, 0.5);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            std::string gname = "g_pt2nd_vs_q_eta2nd_" + musign + "_" + trg_maps[trg] + "_sepr_py_" + q_eta_bin + "_divided";
            

            if (debug_mode) std::cout << "Graph name: " << gname << std::endl;

            TGraphAsymmErrors* g = dynamic_cast<TGraphAsymmErrors*>(f->Get(gname.c_str()));
            if (!g) {
                std::cerr << "Warning: graph " << gname << " not found in file." << std::endl;
                continue;
            }
            
            TF1* fit = fitTurnOn(g);

            if (debug_mode) std::cout << "test fitting drawOneCanvas: " << fit->Eval(8) << std::endl;

            if (!fit){
                std::cerr << "Warning: fitted function for graph " << gname << " does not exist!" << std::endl;
                return;
            }

            fit->Write();

            g->SetMarkerColor(kBlack);
            g->SetLineColor(kBlack);
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.9);

            g->GetXaxis()->SetTitle("p_{T}");
            g->GetYaxis()->SetTitle("#epsilon");

            fit->SetLineColor(kRed);
            fit->SetLineWidth(2);

            g->GetXaxis()->SetTitle("p_{T} [GeV]");
            g->GetYaxis()->SetTitle("Trigger Efficiency");
            g->GetYaxis()->SetRangeUser(0, 1.1);
            g->Draw("AP");

            fit->Draw("SAME");

            std::string musign_label = (musign == "sign1")? "#mu^{+}" : "#mu^{-}";
            auto q_eta_pair = hProjNameToPair(gname);
            std::string q_eta_label = pairToLegendLabel(q_eta_pair);
            leg->AddEntry(g, (trg + ", " + musign_label).c_str(), "lp");
            leg->AddEntry("", q_eta_label.c_str(), "");
            
            leg->Draw("SAME");

            ++idx;
        } // loop over q * eta bins (subplots)

        c->SaveAs(Form("%s%strg_effcy_pT_fitting_%s_%s.png", data_dir.c_str(), fitting_outdir.c_str(), trg_maps.at(trg).c_str(), musign_maps.at(musign).c_str()));
    } // end of function drawOneCanvas
};

void single_muon_trig_effcy_pT_fitting (){
    SingleMuEffcyPtTurnOnFitter* fitter = new SingleMuEffcyPtTurnOnFitter();
    fitter->fitting_mode = SingleMuEffcyPtTurnOnFitter::erf_plus_log;
    fitter->Run();
}
