#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <cmath>
#include <stdexcept>
#include "TMath.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cstdio>

#include "MuonObjectsParamsAndHelpers/PbPbBaseClass.h"
#include "Utilities/proj_range_to_suffix.cxx"
#include "RDFBasedHistFilling/CommonEffcyConfig.h"

// ============================================================
// Base class: all shared logic
// ============================================================
class SingleMuEffcyPtTurnOnFitterBase {
protected:
    // ---- I/O that differs by dataset ----
    std::string data_dir;
    std::string fitting_outdir;
    std::string plot_outdir;   // PNG output directory (under plots/)
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

    // DERIVED from CommonEffcyConfig -- never retyped. This list used to be a hand-maintained
    // copy of q_eta_proj_ranges_fine_excl_gap that had to be edited in lockstep with three other
    // copies. Since round 8 the NOMINAL binning is the CONTIGUOUS COARSE one (gaps INCLUDED), so
    // every muon has a fitted turn-on and the unfitted-2D gap fallback is gone.
    std::vector<std::string> q_eta_bins_for_pT_trg_effcy_graphs = [] {
        static const CommonEffcyConfig cfg{};
        std::vector<std::string> v;
        for (const auto& r : cfg.q_eta_proj_ranges_coarse_incl_gap) v.push_back(pairToSuffix(r));
        return v;
    }();

    static bool fileExists(const std::string& dir){
        return (gSystem->AccessPathName(dir.c_str()) == kFALSE);
    }
    static void makeDirIfNeeded(const std::string& dir){
        if (!fileExists(dir)) gSystem->mkdir(dir.c_str(), kTRUE /*recursive*/);
    }

public:
    bool debug_mode = false;
    bool isBNL = true;

    // Muon working-point selector for the trigger turn-on. Default = TIGHT (nominal WP as of
    // 2026-07-07; docs/tracking/tight_wp_default_change.md). The WP is NOT applied in this fitter:
    // it is applied UPSTREAM in the RDF hist-filling that produces the input q*eta trigger-eff
    // graphs (probe/pair tight selection). This suffix is inserted into BOTH the input-graph
    // filename and the output fit filename + PNG dir, so it MUST match the suffix the hist-filling
    // writes on the graph file. Medium (legacy, unsuffixed graph/fit files) = "".
    std::string wp_suffix = "";

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

        outfile_name = data_dir + fitting_outdir + "single_mu_effcy_pT_fit" + wp_suffix + ".root";
        makeDirIfNeeded(data_dir + fitting_outdir);
        if (!plot_outdir.empty()) makeDirIfNeeded(plot_outdir);

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

        // Build from a TFormula STRING (not a lambda) so the analytic formula
        // PERSISTS when the TF1 is written to a ROOT file. A lambda/compiled TF1
        // can only be saved as a sampled grid (fSave) over [pT_min,pT_max], and on
        // read-back Eval interpolates that grid and returns 0 outside the range —
        // the trigger-eff high-pT bug. The 4.0 pivot in the log/linear term is
        // pT_min (kept numeric to match the original).
        std::string erf_core = "[2]*0.5*(1.0+TMath::Erf((x-[0])/(sqrt(2.0)*[1])))";
        std::string formula  = (fitting_mode == erf_plus_log)
            ? erf_core + "*(1.0+[3]*TMath::Log(1.0+(x-4.0)/4.0))"
            : erf_core + "*(1.0+[3]*((x-4.0)/4.0))";
        TF1* fTurnOn = new TF1(fname.c_str(), formula.c_str(), pT_min, pT_max);

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

        // TFormula STRING (not a lambda) so the analytic formula PERSISTS on Write
        // (a lambda TF1 saves only a sampled grid → Eval returns 0 outside [4,60]).
        // The 4.0 pivot in the log/linear term is pT_min (numeric, matches original).
        std::string fermi_core = "[0]/(1.0+TMath::Exp(([1]-x)/[2]))";
        std::string formula    = (fitting_mode == fermi_plus_log)
            ? fermi_core + "*(1.0+[3]*TMath::Log(1.0+(x-4.0)/4.0))"
            : fermi_core + "*(1.0+[3]*((x-4.0)/4.0))";
        TF1* fTurnOn = new TF1(fname.c_str(), formula.c_str(), pT_min, pT_max);

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

        // TFormula STRING (not a lambda) so the analytic formula PERSISTS on Write.
        TF1* fTurnOn = new TF1(fname.c_str(),
                               "[2]*0.5*(1.0+TMath::Erf((x-[0])/(sqrt(2.0)*[1])))",
                               pT_min - 0.1, pT_max);

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

        TCanvas* c = new TCanvas(cname.c_str(), "Trigger Turn-on Curves", 1500, 1300);
        // pad grid derived from the bin count (10 coarse bins + legend since round 8);
        // a hardcoded 4x3 silently drew into the wrong pads when the count changed
        {
            const int nb = static_cast<int>(q_eta_bins_for_pT_trg_effcy_graphs.size());
            const int ncols = 3, nrows = (nb + 1 + ncols - 1) / ncols;  // +1 legend pad
            c->Divide(ncols, nrows);
        }
        c->SetGrid();

        int idx = 0;
        for (const auto& q_eta_bin : q_eta_bins_for_pT_trg_effcy_graphs){
            c->cd(idx + 1);
            gPad->SetLogx();

            // The legend lives in the reserved strip ABOVE the frame. Inside the frame there is
            // no safe corner: the turn-on sweeps the lower left and plateaus across the upper
            // right, and on the peripheral Pb+Pb panels the old box was drawn straight through
            // both the points and the fit.
            TLegend* leg = new TLegend(0.56, 0.735, 0.98, 0.795);
            leg->SetNColumns(2);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.055);

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

            // Axis text sized for a 3-column pad grid on a 1500x1000 canvas: at the ROOT
            // defaults the titles and labels come out at ~3 pt in the saved PNG, which is
            // unreadable once the figure is placed in a document. SetMoreLogLabels() adds the
            // intermediate decade ticks -- without it a 4-60 GeV log axis carries a single "10".
            gPad->SetLeftMargin(0.17);
            gPad->SetBottomMargin(0.18);
            gPad->SetTopMargin(0.30);   // reserved strip: equation, parameters, bin label, legend
            g->GetXaxis()->SetTitle("p_{T} [GeV]");
            g->GetYaxis()->SetTitle("#varepsilon_{single}");
            g->GetYaxis()->SetRangeUser(0, 1.1);
            g->GetXaxis()->SetTitleSize(0.085);  g->GetXaxis()->SetLabelSize(0.072);
            g->GetYaxis()->SetTitleSize(0.085);  g->GetYaxis()->SetLabelSize(0.072);
            g->GetYaxis()->SetTitleOffset(0.90);
            g->GetXaxis()->SetTitleOffset(0.95);
            g->GetXaxis()->SetMoreLogLabels();
            g->GetXaxis()->SetNoExponent();
            g->Draw("AP");
            // Draw over the FITTED range, taken from the TF1 itself so the number is never
            // retyped. adjustLogXRange stretches the axis to the reference histogram's first bin
            // (0.4 GeV), a decade below the lowest measured point, which squeezed the whole
            // turn-on into the right half of the frame.
            g->GetXaxis()->SetLimits(fit->GetXmin(), fit->GetXmax());
            fit->Draw("SAME");

            std::string musign_label = (musign == "sign1")? "#mu^{+}" : "#mu^{-}";
            auto q_eta_pair = hProjNameToPair(gname);
            std::string q_eta_label = pairToLegendLabel(q_eta_pair);

            leg->AddEntry(g, (trg_maps.at(trg) + ", " + musign_label).c_str(), "lp");
            leg->AddEntry(fit, "fit", "l");
            leg->Draw("SAME");

            // The fitted parameters, with the equation above them (atlas-plotting.md P2: a
            // drawn fit without its exact form and its parameter values fails). Written in the
            // upper-left, which the turn-on leaves empty at low pT.
            TLatex ft;
            ft.SetNDC();
            ft.SetTextFont(42);
            ft.SetTextSize(0.056);
            ft.DrawLatex(0.19, 0.750, q_eta_label.c_str());
            ft.SetTextSize(0.058);
            // Written FLAT, not as nested fractions: a superscript inside a fraction denominator
            // inside an outer fraction renders as a smudge at this text size, which defeats the
            // point of drawing the equation at all.
            if (fitting_mode != erf_plus_log && fitting_mode != fermi_plus_log)
                throw std::runtime_error("SingleMuEffcyPtTurnOnFitter: the drawn equation is only "
                                         "written for erf_plus_log and fermi_plus_log; add the "
                                         "form before enabling another fitting_mode");
            ft.DrawLatex(0.19, 0.945, (fitting_mode == erf_plus_log)
                ? "#varepsilon = 0.5 P [1 + erf((p_{T}-m)/(#sqrt{2}s))] "
                  "[1 + c ln(1+(p_{T}-4)/4)]"
                : "#varepsilon = P / [1 + exp((x_{0}-p_{T})/w)] #times "
                  "[1 + c ln(1+(p_{T}-4)/4)]");
            // EVERY free parameter of the drawn equation is printed: the logarithmic term adds
            // up to ~20 % by 60 GeV, so a reader given only P cannot reconstruct the curve they
            // are looking at.
            const double chi2ndf = fit->GetNDF() > 0 ? fit->GetChisquare() / fit->GetNDF() : 0.0;

            // A parameter that has run onto one of its own fit limits is a CONSTRAINT, not a
            // measurement, and must be identifiable as such. The test is relative to the width of
            // the allowed interval: an absolute tolerance flags a value pinned at 0.10000000 and
            // misses its neighbour at 0.09999849, which is the same boundary. Applied to every
            // printed parameter, not only to c.
            auto at_limit = [&fit](int ipar) {
                double lo = 0., hi = 0.;
                fit->GetParLimits(ipar, lo, hi);
                if (hi <= lo) return "";                      // unbounded: nothing to report
                const double v = fit->GetParameter(ipar);
                const double tol = 1e-4 * (hi - lo);
                return (std::fabs(v - lo) < tol || std::fabs(v - hi) < tol) ? "*" : "";
            };
            const int i_plateau = (fitting_mode == erf_plus_log) ? 2 : 0;
            const int i_shape1  = (fitting_mode == erf_plus_log) ? 0 : 1;
            const int i_shape2  = (fitting_mode == erf_plus_log) ? 1 : 2;
            const bool any_at_limit = *at_limit(i_plateau) || *at_limit(i_shape1)
                                   || *at_limit(i_shape2)  || *at_limit(3);

            ft.SetTextSize(0.048);
            ft.DrawLatex(0.19, 0.885,
                         Form((fitting_mode == erf_plus_log)
                                  ? "P = %.3f%s, m = %.2f%s GeV, s = %.2f%s GeV"
                                  : "P = %.3f%s, x_{0} = %.2f%s GeV, w = %.2f%s GeV",
                              fit->GetParameter(i_plateau), at_limit(i_plateau),
                              fit->GetParameter(i_shape1),  at_limit(i_shape1),
                              fit->GetParameter(i_shape2),  at_limit(i_shape2)));
            ft.DrawLatex(0.19, 0.825,
                         Form("c = %.4f%s,  #chi^{2}/ndf = %.2f%s",
                              fit->GetParameter(3), at_limit(3), chi2ndf,
                              any_at_limit ? "    * at a fit limit" : ""));

            ++idx;
        }

        std::string png_dir = plot_outdir.empty() ? (data_dir + fitting_outdir) : (plot_outdir + "/");
        if (ctr.empty()){
            c->SaveAs(Form("%strg_effcy_pT_fitting_%s_%s.png",
                           png_dir.c_str(),
                           trg_maps.at(trg).c_str(), musign_maps.at(musign).c_str()));
        } else {
            c->SaveAs(Form("%strg_effcy_pT_fitting_%s%s_%s.png",
                           png_dir.c_str(),
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
        std::string base = isBNL ? "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/"
                                 : "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/";
        data_dir = base + "pp_2024/";
        infile_name = data_dir + "histograms_real_pairs_pp_2024_single_mu4_coarse_q_eta_bin_qeta_fid" + wp_suffix + ".root";
        h2d_ref_name = "h_pt2nd_vs_q_eta2nd_sign1_mu4";
        plot_outdir = base + "plots/pp_trigger_efficiency/mu4/pT_fitting/pp24" + wp_suffix;
    }
};

// ============================================================
// Derived: PbPb
// ============================================================
class SingleMuEffcyPtTurnOnFitterPbPb : public SingleMuEffcyPtTurnOnFitterBase, public PbPbBaseClass<SingleMuEffcyPtTurnOnFitterPbPb> {
public:
    int run_year = 25;
    int RunYear() const { return run_year; }

    SingleMuEffcyPtTurnOnFitterPbPb(
                                    int run_year_input = 25,
                                    const std::string ctr_binning_version_input = "default",
                                    bool debug_mode_input=false)
                    : SingleMuEffcyPtTurnOnFitterBase(debug_mode_input),
                      run_year(run_year_input){
                        ctr_binning_version = ctr_binning_version_input;
                    }


protected:
    void initialize() override {
        PbPbBaseClass<SingleMuEffcyPtTurnOnFitterPbPb>::InitializePbPb();
        SingleMuEffcyPtTurnOnFitterBase::initialize();
    }

    void configureIO() override {
        std::string yr = std::to_string(run_year);
        std::string base = isBNL ? "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/"
                                 : "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/";
        data_dir = base + "pbpb_20" + yr + "/";
        infile_name = data_dir + "histograms_real_pairs_pbpb_20" + yr + "_single_mu4_coarse_q_eta_bin_qeta_fid" + wp_suffix + ".root";
        h2d_ref_name = "h_pt2nd_vs_q_eta2nd_ctr0_5_sign1_mu4";
        plot_outdir = base + "plots/pbpb_trigger_efficiency/mu4/pT_fitting/pbpb" + yr + wp_suffix;
    }

    std::vector<std::string> centralitySuffixes() const override {
        return ctr_bins;
    }
};

// ------------------------------------------------------------
// ROOT entry points (like your originals)
// ------------------------------------------------------------
// wp_suffix: "" (nominal = TIGHT, default — unsuffixed filenames, matches the crossx convention) or
// "_medium_wp" (Medium systematic). Must match the suffix the RDF
// hist-filling writes on the input q*eta trigger-eff graph file.
void single_muon_trig_effcy_pT_fitting(const std::string& wp_suffix = "") {
    auto* fitter = new SingleMuEffcyPtTurnOnFitterPP();
    fitter->wp_suffix = wp_suffix;
    fitter->fitting_mode = SingleMuEffcyPtTurnOnFitterBase::erf_plus_log;
    fitter->Run();
    delete fitter;
}

void single_muon_trig_effcy_pT_fitting_PbPb(int year = 25, const std::string& wp_suffix = "") {
    auto* fitter = new SingleMuEffcyPtTurnOnFitterPbPb(year);
    fitter->wp_suffix = wp_suffix;
    fitter->fitting_mode = SingleMuEffcyPtTurnOnFitterBase::fermi_plus_log;
    fitter->Run();
    delete fitter;
}
