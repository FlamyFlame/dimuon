#include"TrigEffPlotterBaseClass.h"

// =========================================================================
// configure input data files
// =========================================================================
void TrigEffPlotterPbPb::configureDataFiles(int runYear)
{
    switch (runYear % 2000)
    {
    case 23:
        data_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pbpb_2023/";
        fname_single_mu4 = data_dir + "histograms_real_pairs_pbpb_2023_single_mu4.root";
        fname_MB         = data_dir + "histograms_real_pairs_pbpb_2023_MB.root";
        break;
    case 24:
        data_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pbpb_2024/";
        fname_single_mu4 = data_dir + "histograms_real_pairs_pbpb_2024_single_mu4.root"; // (fix your extra .root.root)
        fname_MB         = data_dir + "histograms_real_pairs_pbpb_2024_MB.root";
        break;
    default:
        throw std::runtime_error("runYear must be 23 or 24 for PbPb");
    }
}

// =========================================================================
//  plot1D – main refactored 1‑D efficiency routine
// =========================================================================
void TrigEffPlotterPbPb::plot1D(const std::string& var)
{

    // ── enlarge default legend box if filter-labels are long ─────────────
    auto autoBox = [&](Rect r, const StrVec& filters) -> Rect
    {
        // longest suffix length (ignore the empty string for the default curve)
        size_t L = 0;
        for (auto& fs : filters){
            if (!fs.empty()){
                std::string filter_label = fs;
                if (filt_suffix_to_label_map.find(fs) != filt_suffix_to_label_map.end()){ // filter label found
                    filter_label = filt_suffix_to_label_map.at(fs);
                }
                L = std::max(L, filter_label.size());
            }
        }

        if (L > 4) {                                   // threshold
            double extra = std::min(0.01*(L - 4), 0.20);   // 0.01 per char, cap 0.10
            r[0] = std::max(0.05, r[0] - extra);           // move left
            // r[3] = std::min(0.95, r[3] + extra*0.1);       // comment out: do not grow talled
        }
        return r;
    };

    // ===============================================================
    // Decide which filters will be looped over (default + comparisons)
    // ===============================================================
    StrVec filters = filter_suffix_list;
    bool plot_weighted = find(filters.begin(), filters.end(), "_inv_w_by_single_mu_effcy") != filters.end();
    if (plot_weighted){ // inverted by weight --> ploted weighted first to set the y-axis range
        filters.push_back("");
    }else{
        filters.insert(filters.begin(), ""); // no weighting, ensure default first
    }

    // For each drawing mode iterate ---------------------------------------------------
    
    auto drawing_modes = draw_single_b_signal? std::vector<SignalDrawingMode>{SignalDrawingMode::Signed, SignalDrawingMode::Signal} : std::vector<SignalDrawingMode>{SignalDrawingMode::Signed};

    for (auto mode : drawing_modes) {
        // Is this mode requested by *any* of the (comparison) filters?
        bool modeNeeded = filter_suffix_list.empty() ? true : std::any_of(filter_suffix_list.begin(), filter_suffix_list.end(), [&](const std::string& fs){
            return filt_to_mode_map.at(fs).at(mode);
        });
        if (!modeNeeded) continue;

        // Canvas bookkeeping -------------------------------------------------------
        bool isSigned = (mode == SignalDrawingMode::Signed);
        std::string canvasName = Form("c_%s_%s", var.c_str(), (isSigned?"signed":"signal"));
        int   canvW  = isSigned ? 1200 : 700;
        int   canvH  = 500;
        TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), canvW, canvH);
        if (isSigned) c->Divide(2,1);

        // Prepare legend & colour styles ----------------------------------------
        Color_t colBase[3] = {kRed+1, kAzure+2, kGreen+2}; // 2mu4, mu4_mu4noL1, (spare)
        Color_t colFilt[3] = {kOrange+7, kAzure+1, kGreen+1}; // one-shade-lighter palette for filtered / comparison curves
        Color_t colFilt2[3] = {kMagenta, kAzure, kGreen+3};
        Style_t mkrBase[3] = {24, 20, 22};

        // Iterate over pads (sign1 / sign2 or single)
        int nPads = isSigned?2:1;
        std::string signSuffixes[2] = {"_sign1", "_sign2"};

        for (int pad = 0; pad < nPads; ++pad) {
            TPad* padPtr = isSigned? (TPad*)c->cd(pad+1) : (TPad*)c->cd();
            auto xy = xyFor(var);
            padPtr->SetLogx(xy.first); padPtr->SetLogy(xy.second);

            Rect def = {0.43,0.12,0.8,0.35};
            Rect def_excl_inv_w_by_single_mu_effcy = {0.45,0.45,0.83,0.63}; // use for "_excl_inv_w_by_single_mu_effcy"
            Rect def_default_only = {0.5,0.14,0.77,0.3};
            if (filter_suffix_list.empty()) def = def_default_only;

            if (find(filters.begin(), filters.end(), "_excl_inv_w_by_single_mu_effcy") != filters.end()) def = def_excl_inv_w_by_single_mu_effcy;

            Rect box = def;

            if      (mode==SignalDrawingMode::Signed  && legendPosSigned.count(var))
                box = legendPosSigned.at(var);
            else if (mode==SignalDrawingMode::Signal  && legendPosSignal.count(var))
                box = legendPosSignal.at(var);
            else if (!filter_suffix_list.empty()){ // adjust legend box size to the label length of the (first) filter
                box = autoBox(box, filter_suffix_list);
            }

            TLegend* leg = new TLegend(box[0], box[1], box[2], box[3]);
            leg->SetBorderSize(0);

            bool firstDrawn = false;

            // loop over filters → lines (default first, then comparisons)
            int filterIdx = 0;
            for (const auto& fs : filters) {
                if (!filt_to_mode_map.at(fs).at(mode)) continue; // skip filter that should not appear in this mode

                // For each (valid) trigger produce TGraphAsymmErrors
                std::vector<std::pair<std::string,int>> trigOrder; // order & colour idx

                if (draw2mu4){
                    trigOrder = { {"mu4_mu4noL1", 0}, {"2mu4", 1} }; 
                }else{
                    trigOrder = { {"mu4_mu4noL1", 0} };
                }

                for (auto [trg, colourIdx] : trigOrder) {

                    // denominator hist (mu4 / trigger-specific denominator)
                    std::string denomName = (mappedTrigger(fs, trg + "_denom") != "NONE")? mappedTrigger(fs, trg + "_denom") : mappedTrigger(fs, "mu4");
                    
                    if (denomName.empty() || denomName == "NONE") {
                        throw std::runtime_error("Base trigger mu4 missing for filter " + fs);
                    }

                    std::string denomHist = Form("h_%s%s_%s", var.c_str(),
                                                  (mode==SignalDrawingMode::Signal?"_w_single_b_sig_sel":signSuffixes[pad].c_str()),
                                                  denomName.c_str());
                    
                    if (debug_mode) std::cout << "Denominator name: " << denomHist << std::endl;
                    TH1* h_denom = getHist<TH1>(denomHist);
                    if (!h_denom) continue; // can't proceed – skip filter


                    std::string mapped = mappedTrigger(fs, trg);
                    if (mapped.empty() || mapped == "NONE") continue; // skip incompatible trigger

                    std::string numHist = Form("h_%s%s_%s", var.c_str(),
                                                  (mode==SignalDrawingMode::Signal?"_w_single_b_sig_sel":signSuffixes[pad].c_str()),
                                                  mapped.c_str());

                    TH1* h_num = getHist<TH1>(numHist);
                    if (!h_num) continue;

                    Color_t useCol;
                    if (!plot_weighted){
                        useCol = (filterIdx == 0) ? colBase[colourIdx]   // default
                                                : colFilt[colourIdx];  // filtered
                        if (filterIdx == 2) useCol = colFilt2[colourIdx];
                    } else{ // plot weighted: default appears at last
                        useCol = (filterIdx == filters.size() - 1) ? colBase[colourIdx]   // filtered
                                                                   : colFilt[colourIdx];  // base
                        if (filters.size() == 3 && filterIdx == 1) useCol = colFilt2[colourIdx];
                    }

                    TH1* h_ratio = nullptr;
                    TGraphAsymmErrors* g = nullptr;

                    if (std::find(filters_weighted.begin(), filters_weighted.end(), fs) != filters_weighted.end()){ // filtered ratio is weighted binomial --> use TH1 instead of TGraphAsymmErrors
                        h_ratio = (TH1*)h_num->Clone(Form("%s_divided", h_num->GetName()));
                        h_ratio->Divide(h_denom);

                        SetStyle(h_ratio, useCol, mkrBase[colourIdx]);
                        bool isFiltered = plot_weighted? (filterIdx < filters.size() - 1) : (filterIdx > 0);
                        if (isFiltered) {
                            // comparison filter → dashed
                            h_ratio->SetLineStyle(2);
                            h_ratio->SetMarkerStyle(mkrBase[colourIdx]+4);
                        }

                        if (!firstDrawn) {
                            h_ratio->Draw("E1");
                            h_ratio->GetYaxis()->SetRangeUser(0,1.05);
                            h_ratio->GetXaxis()->SetTitle(h_denom->GetXaxis()->GetTitle());
                            h_ratio->GetYaxis()->SetTitle("#epsilon");
                            
                            gPad->RedrawAxis();
                            gPad->Update(); // ensure pad limits are known

                            double xmin = gPad->GetUxmin();
                            double xmax = gPad->GetUxmax();
                            
                            TLine *line = new TLine(xmin, 1.0, xmax, 1.0);
                            line->SetLineColor(kBlack);
                            line->SetLineStyle(2); // dashed
                            line->SetLineWidth(2);
                            line->Draw("same");

                            firstDrawn = true;
                        } else {
                            h_ratio->Draw("E1 SAME");
                        }
                    } else{ // unweighted --> use TGraphAsymmErrors
                        g = new TGraphAsymmErrors();
                        g->BayesDivide(h_num, h_denom);

                        SetStyle(g, useCol, mkrBase[colourIdx]);
                        bool isFiltered = plot_weighted? (filterIdx < filters.size() - 1) : (filterIdx > 0);
                        if (isFiltered) {
                            // comparison filter → dashed
                            g->SetLineStyle(2);
                            g->SetMarkerStyle(mkrBase[colourIdx]+4);
                        }

                        if (!firstDrawn) {
                            g->Draw("APL");
                            g->GetYaxis()->SetRangeUser(0,1.);
                            g->GetXaxis()->SetTitle(h_denom->GetXaxis()->GetTitle());
                            g->GetYaxis()->SetTitle("#epsilon");
                            
                            if (xy.first) adjustLogXRange(g, h_denom);

                            firstDrawn = true;
                        } else {
                            g->Draw("PL SAME");
                        }
                    }

                    // Legend entry text
                    std::string label = trg;
                    if (!fs.empty()) {
                        std::string fs_label = (fs[0] == '_')? fs.substr(1) : fs; // default filter label: filter suffix but trimming off "_"
                        if (filt_suffix_to_label_map.find(fs) != filt_suffix_to_label_map.end()){ // filter label found
                            fs_label = filt_suffix_to_label_map.at(fs);
                        }
                        label += ", " + fs_label;
                    }
                    if (std::find(filters_weighted.begin(), filters_weighted.end(), fs) != filters_weighted.end()){ // filtered ratio is weighted binomial --> use TH1 instead of TGraphAsymmErrors
                        leg->AddEntry(h_ratio, label.c_str(), "lep");
                    } else{
                        leg->AddEntry(g, label.c_str(), "lep");
                    }
                } // end loop over triggers
                ++filterIdx;
            } // end loop over filters
            leg->Draw();
        } // end loop over pads (subplots)

        // Save png -----------------------------------------------------------------
        std::string pngDir = data_dir + "trig_effcy_plots";
        for (const auto& fs : filter_suffix_list) {
            pngDir += fs;
        }
        if (draw2mu4) pngDir += "_w_2mu4";
        makeDirIfNeeded(pngDir);

        std::string fn = pngDir + "/" + var + "_trig_effcy";
        for (const auto& fs : filter_suffix_list) {
            fn += fs;
        }
        fn += mode_to_png_suffix.at(mode) + ".png";
        c->SaveAs(fn.c_str());
    }

    // ================================================================
    //  Special OpAndSignal – one file per *comparison* filter
    // ================================================================
    // ----------------------------------------------------------------
    //  Op-and-signal mode:
    //    • if filter_suffix_list is empty  --> draw DEFAULT
    //    • otherwise                      --> draw one file per filter (skip default; check if OpAndSignal mode turned on)
    // ----------------------------------------------------------------
    StrVec opFilters = filter_suffix_list.empty() ? StrVec{""} : filter_suffix_list;
    for (const auto& fs : opFilters) {
        // skip comparison filters that don’t request this mode
        if (!fs.empty() &&
            !filt_to_mode_map.at(fs).at(SignalDrawingMode::OpAndSignal)) continue;


        TCanvas* c = new TCanvas(Form("c_opand_%s_%s", var.c_str(), fs.c_str()), "", 700, 500);
        auto xy = xyFor(var);
        c->SetLogx(xy.first); c->SetLogy(xy.second);

        int colourIdx = (!draw2mu4)? 1 : 0; // use Azure for 2mu4; if mu4_mu4noL1 present will be red

        // Legend
        Rect def = {0.47,0.12,0.87,0.36};
        if (legendPosOpSignal.count(var))
            def = legendPosOpSignal.at(var);
        else if (!filter_suffix_list.empty()){ // adjust legend box size to the label length of the (first) filter
            def = autoBox(def, filter_suffix_list);
        }

        TLegend* L = new TLegend(def[0],def[1],def[2],def[3]);
        L->SetBorderSize(0);

        bool firstCurve = true;

        for (auto trg : trg_list) {

            std::string denomName = (mappedTrigger(fs, trg + "_denom") != "NONE")? mappedTrigger(fs, trg + "_denom") : mappedTrigger(fs, "mu4");
            if (denomName.empty() || denomName == "NONE") continue;

            std::string denomSign2 = use_sepr_for_op_and_signal? Form("h_%s_sign2_%s_sepr", var.c_str(), denomName.c_str()) : Form("h_%s_sign2_%s", var.c_str(), denomName.c_str());
            std::string denomSel = Form("h_%s_single_b_%s", var.c_str(), denomName.c_str());

            if (debug_mode){
                std::cout << "Denominator name all opposite sign: " << denomSign2 << std::endl;
                std::cout << "Denominator name signal: " << denomSel << std::endl;
            }

            TH1* h_den_sign2 = getHist<TH1>(denomSign2);
            TH1* h_den_sel   = getHist<TH1>(denomSel);
            if (!h_den_sign2 || !h_den_sel) continue;

            std::string mapped = mappedTrigger(fs, trg);
            if (mapped.empty() || mapped == "NONE") continue;

            std::string numSign2 = use_sepr_for_op_and_signal? Form("h_%s_sign2_%s_sepr", var.c_str(), mapped.c_str()) : Form("h_%s_sign2_%s", var.c_str(), mapped.c_str());
            std::string numSel = Form("h_%s_single_b_%s", var.c_str(), mapped.c_str());

            TH1* h_num_s2 = getHist<TH1>(numSign2);
            TH1* h_num_sel = getHist<TH1>(numSel);
            if (!h_num_s2 || !h_num_sel) continue;

            TH1* h_ratio_s2 = nullptr;
            TH1* h_ratio_sel = nullptr;
            TGraphAsymmErrors* g_s2 = nullptr;
            TGraphAsymmErrors* g_sel = nullptr;

            // sign2 efficiency
            if (std::find(filters_weighted.begin(), filters_weighted.end(), fs) != filters_weighted.end()){ // filtered ratio is weighted binomial --> use TH1 instead of TGraphAsymmErrors
                
                h_ratio_s2 = (TH1*)h_num_s2->Clone(Form("%s_divided", h_num_s2->GetName()));
                h_ratio_s2->Divide(h_den_sign2);

                SetStyle(h_ratio_s2, colourIdx==0? kAzure+2 : kRed+1, colourIdx==0?24:20);

                // sig‑sel efficiency
                h_ratio_sel = (TH1*)h_num_sel->Clone(Form("%s_divided", h_num_sel->GetName()));
                h_ratio_sel->Divide(h_den_sel);

                SetStyle(h_ratio_sel, colourIdx==0? kAzure+1 : kRed+2, colourIdx==0?25:24, 2, 2);
                
                if (firstCurve) {
                    h_ratio_s2->Draw("E1");                       // sets axes with first curve
                    h_ratio_s2->GetYaxis()->SetRangeUser(0., 1.05);
                    h_ratio_s2->GetXaxis()->SetTitle(h_den_sel->GetXaxis()->GetTitle());
                    h_ratio_s2->GetYaxis()->SetTitle("#epsilon");
                    
                    gPad->RedrawAxis();
                    gPad->Update(); // ensure pad limits are known

                    double xmin = gPad->GetUxmin();
                    double xmax = gPad->GetUxmax();
                    
                    TLine *line = new TLine(xmin, 1.0, xmax, 1.0);
                    line->SetLineColor(kBlack);
                    line->SetLineStyle(2); // dashed
                    line->SetLineWidth(2);
                    line->Draw("same");

                    firstCurve = false;                      // ← axis is now set
                } else {
                    h_ratio_s2->Draw("E1 SAME");                   // any later curve
                }
                h_ratio_sel->Draw("PL SAME");                      // dashed sig-sel curve
            } else{
                g_s2 = new TGraphAsymmErrors();
                g_s2->BayesDivide(h_num_s2, h_den_sign2);
                SetStyle(g_s2, colourIdx==0? kAzure+2 : kRed+1, colourIdx==0?24:20);

                // sig‑sel efficiency
                g_sel = new TGraphAsymmErrors();
                g_sel->BayesDivide(h_num_sel, h_den_sel);
                SetStyle(g_sel, colourIdx==0? kAzure+1 : kRed+2, colourIdx==0?25:24, 2, 2);

                if (firstCurve) {
                    g_s2->Draw("APL");                       // sets axes with first curve
                    g_s2->GetYaxis()->SetRangeUser(0., 1.);
                    g_s2->GetXaxis()->SetTitle(h_den_sel->GetXaxis()->GetTitle());
                    g_s2->GetYaxis()->SetTitle("#epsilon");

                    if (xy.first) adjustLogXRange(g_s2, h_den_sign2);

                    firstCurve = false;                      // ← axis is now set
                } else {
                    g_s2->Draw("PL SAME");                   // any later curve
                }
                g_sel->Draw("PL SAME");                      // dashed sig-sel curve
            }

            colourIdx++;

            std::string op_legend_fs = trg;
            std::string sig_legend_fs = trg;

            if (!fs.empty()){
                std::string fs_label = (fs[0] == '_')? fs.substr(1) : fs; // default filter label: filter suffix but trimming off "_"
                if (filt_suffix_to_label_map.find(fs) != filt_suffix_to_label_map.end()){ // filter label found
                    fs_label = filt_suffix_to_label_map.at(fs);
                }                    
                op_legend_fs += " (all op sign, " + fs_label + ")";
                sig_legend_fs += " (single b signal, " + fs_label + ")";
            }else{
                op_legend_fs += " (all op sign)";
                sig_legend_fs += " (single b signal)";
            }

            if (std::find(filters_weighted.begin(), filters_weighted.end(), fs) != filters_weighted.end()){ // filtered ratio is weighted binomial --> use TH1 instead of TGraphAsymmErrors
                L->AddEntry(h_ratio_s2, op_legend_fs.c_str(), "lep");
                L->AddEntry(h_ratio_sel, sig_legend_fs.c_str(), "lep");
            } else{
                L->AddEntry(g_s2, op_legend_fs.c_str(), "lep");
                L->AddEntry(g_sel, sig_legend_fs.c_str(), "lep");                    
            }
        } // finish loop over muon-pair triggers

        L->Draw();
        
        std::string outdir = data_dir + "trig_effcy_plots" + (fs.empty() ? "" : fs);
        if (draw2mu4) outdir += "_w_2mu4";
        makeDirIfNeeded(outdir + "/op_and_sig");
        std::string fn = Form("%s/op_and_sig/%s_trig_effcy%s%s.png", outdir.c_str(), var.c_str(), fs.c_str(), mode_to_png_suffix.at(SignalDrawingMode::OpAndSignal).c_str());
        c->SaveAs(fn.c_str());
    }
} // end function plot1D

#include <TLatex.h>
#include <TLine.h>

// ----------------------------------------------------------------------
// Plot single-muon efficiency vs centrality, for fixed pT bins (mid-rapidity)
// One PNG per trigger-pair, canvas 2x2, one pad per pT bin.
// ----------------------------------------------------------------------
void TrigEffPlotterPbPb::plotSingleMuEffCtrDep()
{
    // user-specified config
    const std::vector<std::pair<std::string,std::string>> trig_pairs = {
        {"mu4_mu4noL1","mu4"},
        {"2mu4","mu4"}
    };

    const std::vector<std::string> pT_bins_for_ctr_dep = {
        "_pt4_6", "_pt6_8", "_pt8_12", "_pt12_30"
    };

    const std::vector<std::string> pT_labels_for_ctr_dep = {
        "p_{T}: [4, 6) GeV",
        "p_{T}: [6, 8) GeV",
        "p_{T}: [8, 12) GeV",
        "p_{T}: [12, 30) GeV"
    };

    // output dir
    std::string outdir = data_dir + "trig_effcy_plots/single_muon_effcy/ctr_dep";
    makeDirIfNeeded(outdir);  // helper already in your framework :contentReference[oaicite:1]{index=1}

    // simple style palette (reuse your SetStyle helper) :contentReference[oaicite:2]{index=2}
    const Color_t col = kAzure + 2;
    const Style_t mkr = 20;

    for (const auto& tp : trig_pairs)
    {
        const std::string& num_trg = tp.first;
        const std::string& den_trg = tp.second;

        // canvas: 2x2
        std::string cname = Form("c_singleMuEffCtr_%s_over_%s", num_trg.c_str(), den_trg.c_str());
        TCanvas* c = new TCanvas(cname.c_str(), cname.c_str(), 1100, 900);
        c->Divide(2, 2);

        for (int i = 0; i < (int)pT_bins_for_ctr_dep.size(); ++i)
        {
            c->cd(i+1);
            gPad->SetGrid();

            const std::string& pt_bin = pT_bins_for_ctr_dep[i];

            // names exactly as you specified
            const std::string hnum_name =
                "h_ctr" + pt_bin + "_mid_rapidity_" + num_trg + "_sepr";
            const std::string hden_name =
                "h_ctr" + pt_bin + "_mid_rapidity_" + den_trg + "_sepr";

            TH1D* hnum = getHist<TH1D>(hnum_name);
            TH1D* hden = getHist<TH1D>(hden_name);

            if (!hnum || !hden) {
                // getHist already warns; just leave pad blank
                continue;
            }

            // graph efficiency
            auto* g = new TGraphAsymmErrors();
            g->BayesDivide(hnum, hden);

            SetStyle(g, col, mkr);  // your shared helper :contentReference[oaicite:3]{index=3}
            // g->GetYaxis()->SetRangeUser(0.0, 1.00);
            g->GetYaxis()->SetTitle("#epsilon");
            g->GetXaxis()->SetTitle(hden->GetXaxis()->GetTitle()); // should be centrality axis title

            g->Draw("APL");

            // horizontal line at 1
            const double xmin = gPad->GetUxmin();
            const double xmax = gPad->GetUxmax();
            TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");

            // pad label: 2 lines, as requested
            TLatex latex;
            latex.SetNDC(true);
            latex.SetTextSize(0.04);
            latex.DrawLatex(0.16, 0.86, pT_labels_for_ctr_dep[i].c_str());
            latex.DrawLatex(0.16, 0.80, "-0.5 <= #eta < 0.5");
        }

        // save png
        std::string outpng = Form("%s/singleMuEffCtrDep_%s_over_%s.png",
                                  outdir.c_str(), num_trg.c_str(), den_trg.c_str());
        c->SaveAs(outpng.c_str());
    }
} // end function plotSingleMuEffCtrDep

// ======================================================================
// 2-D efficiency canvases  (Signed + Signal modes) – one PNG per filter
// ======================================================================
void TrigEffPlotterPbPb::plot2D(const std::string& var, const std::string& fs)
{
    const auto xy = xyFor(var);

    // ---------- first canvas : sign-1 / sign-2 ------------------------
    std::vector<std::pair<std::string,std::string>> configs;

    if (draw2mu4){
        configs = {{"mu4_mu4noL1","mu4"},{"2mu4","mu4"}};
    } else{
        configs = {{"mu4_mu4noL1","mu4"}};
    }

    int nPads = static_cast<int>(configs.size())*2;         // each config → 2 signs
    int nCols = 2;
    int nRows = (nPads/2);

    TCanvas* c = new TCanvas(Form("c2d_%s_%s",var.c_str(),fs.c_str()),"",1100,450*nRows);
    c->Divide(nCols,nRows);

    int pad=1;
    for (auto& cfg : configs)
    for (auto sign : {std::string("_sign1"),std::string("_sign2")})
    {
        auto numH = getHist<TH2>(hName(var, fs, cfg.first , sign));
        auto denH = getHist<TH2>(hName(var, fs, cfg.second, sign));
        if (!numH||!denH) {++pad; continue;}
        numH = (TH2*)numH->Clone();           // keep originals
        numH->Divide(denH);
        c->cd(pad++);
        gPad->SetRightMargin(0.15);
        gPad->SetLogx(xy.first); gPad->SetLogy(xy.second);
        numH->SetTitle(
            Form("%s / %s  (%s)",
                 cfg.first.c_str(),cfg.second.c_str(),
                 sign=="_sign1"?"same":"opposite"));
        numH->Draw("COLZ");
    }

    std::string outdir = data_dir + "trig_effcy_plots" + (fs.empty()?"":fs);
    if (draw2mu4) outdir += "_w_2mu4";
    
    makeDirIfNeeded(outdir);
    c->SaveAs(Form("%s/%s_2D_trig_effcy%s.png",
                   outdir.c_str(),var.c_str(),fs.c_str()));

    // ---------- second canvas : signal-selection ----------------------
    if (!draw_single_b_signal) return;
    TCanvas* cB = new TCanvas(Form("c2ds_%s_%s",var.c_str(),fs.c_str()),"",700,450);
    cB->cd(); gPad->SetRightMargin(0.15);
    gPad->SetLogx(xy.first); gPad->SetLogy(xy.second);

    // denominator always mu4, sig-selection
    auto hDenSel = getHist<TH2>(hName(var, fs,"mu4","_w_single_b_sig_sel"));
    if (!hDenSel) return;

    for (auto trg : trg_list)
    {
        auto hNumSel = getHist<TH2>(hName(var, fs,trg,"_w_single_b_sig_sel"));
        if (!hNumSel) continue;
        hNumSel = (TH2*)hNumSel->Clone();
        hNumSel->Divide(hDenSel);
        hNumSel->SetTitle((trg+", signal selection").c_str());
        hNumSel->Draw("COLZ SAME");           // pile in same pad
    }
    cB->SaveAs(Form("%s/%s_2D_trig_effcy%s_w_sig_sel.png",
                    outdir.c_str(),var.c_str(),fs.c_str()));
}

// ======================================================================
// 2-D single-muon efficiency canvases - signed by 2nd-muon charge, plot for separated & signal selection
// ======================================================================
// for 2D single-muon effcy plots, always plot 2mu4 --> can cut later
void TrigEffPlotterPbPb::plot2D_SingleMuonEffcy(const std::string& var)
{
    std::cout << "Call function plot2D_SingleMuonEffcy" << std::endl;
    
    const auto xy = xyFor(var);

    // ---------- first canvas : sign-1 / sign-2 ------------------------
    std::vector<std::pair<std::string,std::string>> configs;

    configs = {{"mu4_mu4noL1","mu4"}};

    int nPads = static_cast<int>(configs.size())*2;         // each config → 2 signs
    int nCols = 2;
    int nRows = (nPads/2);

    std::vector<std::string> filters_single_muon_effcy = draw_single_b_signal? std::vector<std::string>{"_sepr", "_w_single_b_sig_sel"} : std::vector<std::string>{"_sepr"};

    for (auto fs : filters_single_muon_effcy){
        for (auto& ctr : ctr_bins){
            TCanvas* c = new TCanvas(Form("c2d_%s%s_%s",var.c_str(),ctr.c_str(),fs.c_str()),"",1100,450*nRows);
            c->Divide(nCols,nRows);

            int pad=1;
            for (auto& cfg : configs){
                for (auto sign : {std::string("_sign1"),std::string("_sign2")})
                {
                    auto numH = getHist<TH2>(hName(var, fs, cfg.first , ctr + sign));
                    auto denH = getHist<TH2>(hName(var, fs, cfg.second, ctr + sign));
                    if (!numH||!denH) {++pad; continue;}

                    numH = (TH2*)numH->Clone();           // keep originals
                    numH->Divide(denH);
                    c->cd(pad++);
                    gPad->SetRightMargin(0.15);
                    gPad->SetLogx(xy.first); gPad->SetLogy(xy.second);
                    std::string single_mu_trg = (cfg.first == "mu4_mu4noL1")? "mu4noL1" : "mu4";
                    numH->SetTitle(
                        Form("%s, %s",
                             single_mu_trg.c_str(),
                             sign=="_sign1"?"#mu^{+}":"#mu^{-}"));
                    numH->Draw("COLZ");
                }
            }
            std::string outdir = data_dir + "trig_effcy_plots/single_muon_effcy";
            makeDirIfNeeded(outdir);
            c->SaveAs(Form("%s/%s_2D_trig_effcy%s%s.png",
                           outdir.c_str(),var.c_str(),ctr.c_str(),fs.c_str()));
        }
    }
}

// ======================================================================
// 2-D profile canvases  (no division) – one PNG per filter
// ======================================================================
void TrigEffPlotterPbPb::plotProfile(const std::string& var, const std::string& fs)
{
    const auto xy = xyFor(var);

    // -------- signed view ---------------------------------------------------
    TCanvas* c = new TCanvas(Form("cProf_%s_%s",var.c_str(),fs.c_str()),"",1100,500);
    c->Divide(2,1);

    for (int i=0;i<2;++i) {
        std::string sign = (i==0?"_sign1":"_sign2");
        auto h = getHist<TH2>(hName(var, fs,"mu4",sign));
        if (!h) continue;
        c->cd(i+1);
        gPad->SetRightMargin(0.15);
        gPad->SetLogx(xy.first); gPad->SetLogy(xy.second); gPad->SetLogz();
        h->SetTitle(Form("mu4 %s",sign=="_sign1"?"same":"opposite"));
        h->Draw("COLZ");
        auto pfx = h->ProfileX(Form("pfx%d",i),1,-1,"s");
        pfx->SetLineColor(kRed); pfx->SetLineWidth(2);
        pfx->Draw("E SAME");
    }

    std::string outdir = data_dir + "trig_effcy_plots" + (fs.empty()?"":fs);
    if (draw2mu4) outdir += "_w_2mu4";
    
    makeDirIfNeeded(outdir);
    c->SaveAs(Form("%s/%s_profile%s.png",outdir.c_str(),var.c_str(),fs.c_str()));

    // -------- signal-selection profile --------------------------------------
    if (!draw_single_b_signal) return;
    auto hSel = getHist<TH2>(hName(var, fs,"mu4","_w_single_b_sig_sel"));
    if (!hSel) return;

    TCanvas* cS = new TCanvas(Form("cProfS_%s_%s",var.c_str(),fs.c_str()),"",550,450);
    cS->cd(); gPad->SetRightMargin(0.15);
    gPad->SetLogx(xy.first); gPad->SetLogy(xy.second); gPad->SetLogz();
    hSel->SetTitle("mu4, signal selection");
    hSel->Draw("COLZ");
    auto pS = hSel->ProfileX("pS",1,-1,"s");
    pS->SetLineColor(kRed); pS->SetLineWidth(2);
    pS->Draw("E SAME");

    cS->SaveAs(Form("%s/%s_profile%s_w_sig_sel.png",
                    outdir.c_str(),var.c_str(),fs.c_str()));
}

// ======================================================================
// Project selected 2D histograms (pair observables) to 1D (X/Y) and plot efficiencies like plot1D()
// ======================================================================
void TrigEffPlotterPbPb::plot2Dto1DEffcyProj()
{
    if (var2Ds_for1Deffcy_proj.empty()) return;

    // default + comparisons (keep weighted filter behaviour consistent with plot1D)
    StrVec filters = filter_suffix_list;
    const bool hasWeighted = (std::find(filters.begin(), filters.end(),
                                        std::string("_inv_w_by_single_mu_effcy")) != filters.end());
    if (hasWeighted) filters.push_back(""); else filters.insert(filters.begin(), "");

    for (const auto& kv : var2Ds_for1Deffcy_proj)
    {
        const std::string& var2d = kv.first;
        const bool doPX = kv.second.first;   // ProjectionX → onto X (right side of "..._vs_")
        const bool doPY = kv.second.second;  // ProjectionY → onto Y (left side of "..._vs_")

        const auto parsed = parseYvsX(var2d);
        const std::string& yVar = parsed.first;
        const std::string& xVar = parsed.second;

        for (int proj = 0; proj < 2; ++proj)
        {
            const bool isProjX = (proj == 0);
            if ((isProjX && !doPX) || (!isProjX && !doPY)) continue;

            const std::string projVar = isProjX ? (xVar.empty() ? var2d : xVar)
                                                : (yVar.empty() ? var2d : yVar);
            if (projVar.empty()) continue; // malformed name

            // ---------------- Signed / Signal (shared canvas structure) ----------------
            auto drawing_modes = draw_single_b_signal? std::vector<SignalDrawingMode>{SignalDrawingMode::Signed, SignalDrawingMode::Signal} : std::vector<SignalDrawingMode>{SignalDrawingMode::Signed};
            for (auto mode : drawing_modes)
            {
                // mode gate (same rule as plot1D)
                bool modeNeeded = filter_suffix_list.empty();
                if (!modeNeeded) {
                    for (const auto& fs : filter_suffix_list) {
                        if (filt_to_mode_map.at(fs).at(mode)) { modeNeeded = true; break; }
                    }
                }
                if (!modeNeeded) continue;

                const bool isSigned = (mode == SignalDrawingMode::Signed);
                const int  nPads    = isSigned ? 2 : 1;
                const char* signSuffixes[2] = {"_sign1", "_sign2"};

                // canvas
                std::string cname = Form("c_proj_%s_%s_%s",
                                         var2d.c_str(),
                                         isProjX ? "PX" : "PY",
                                         isSigned ? "signed" : "signal");
                TCanvas* c = new TCanvas(cname.c_str(), "", isSigned ? 1200 : 700, 500);
                if (isSigned) c->Divide(2,1);

                // legend box (reuse 1D var positioning)
                Rect box = legendBoxFor1DVar(projVar, mode, filter_suffix_list);

                // styles: [0] → mu4_mu4noL1 (red), [1] → 2mu4 (azure)
                Color_t colBase[2] = {kRed+1,   kAzure+2};
                Color_t colFilt[2] = {kRed+2,   kAzure+1};
                Style_t mkrBase[2] = {20, 24};

                for (int pad = 0; pad < nPads; ++pad)
                {
                    TPad* P = isSigned ? (TPad*)c->cd(pad+1) : (TPad*)c->cd();
                    (void)P;
                    auto xy = xyFor(projVar);
                    gPad->SetLogx(xy.first);
                    gPad->SetLogy(xy.second);

                    TLegend* leg = new TLegend(box[0], box[1], box[2], box[3]);
                    leg->SetBorderSize(0);

                    bool firstDrawn = false;
                    int  filterIdx  = 0;

                    for (const auto& fs : filters)
                    {
                        if (!filt_to_mode_map.at(fs).at(mode)) { ++filterIdx; continue; }

                        // triggers to overlay (skip mu4; it’s denom)
                        std::vector<std::pair<std::string,int>> trigOrder =
                                    draw2mu4? std::vector<std::pair<std::string,int>>{{"mu4_mu4noL1", 0}, {"2mu4", 1}}
                                            : std::vector<std::pair<std::string,int>>{{"mu4_mu4noL1", 0}};

                        bool firstTrgInstance = true; // boolean to indicate first trigger instance: if continue without plotting, only increment filterIdx if first trigger instance

                        for (const auto& tp : trigOrder)
                        {
                            const std::string& trg = tp.first;
                            const int idx = tp.second;
                         
                            // denominator 2D → project
                            const std::string denomName = (mappedTrigger(fs, trg + "_denom") != "NONE")? mappedTrigger(fs, trg + "_denom") : mappedTrigger(fs, "mu4");
                            
                            if (denomName.empty() || denomName == "NONE") {
                                if (firstTrgInstance){
                                    ++filterIdx;
                                    firstTrgInstance = false;
                                }
                                continue;
                            }

                            std::string denom2D = Form("h_%s%s_%s",
                                                var2d.c_str(),
                                                (mode == SignalDrawingMode::Signal ? "_w_single_b_sig_sel"
                                                                                   : signSuffixes[pad]),
                                                denomName.c_str());

                            if (debug_mode) std::cout << "Denominator name: " << denom2D << std::endl;
    
                            TH2* hDen2D = getHist<TH2>(denom2D);
                            if (!hDen2D || hDen2D->GetDimension() != 2) {
                                if (firstTrgInstance){
                                    ++filterIdx;
                                    firstTrgInstance = false;
                                }
                                continue;
                            }

                            std::unique_ptr<TH1> hDen1D(
                                isProjX ? hDen2D->ProjectionX(Form("px_%s_%d_den", var2d.c_str(), pad), 1, -1, "e")
                                        : hDen2D->ProjectionY(Form("py_%s_%d_den", var2d.c_str(), pad), 1, -1, "e")
                            );
                            if (!hDen1D) {
                                if (firstTrgInstance){
                                    ++filterIdx;
                                    firstTrgInstance = false;
                                }
                                continue;
                            }

                            const std::string mapped = mappedTrigger(fs, trg);
                            if (mapped.empty() || mapped == "NONE") continue;

                            std::string num2D = Form("h_%s%s_%s",
                                                var2d.c_str(),
                                                (mode == SignalDrawingMode::Signal ? "_w_single_b_sig_sel"
                                                                                   : signSuffixes[pad]),
                                                mapped.c_str());

                            TH2* hNum2D = getHist<TH2>(num2D);
                            if (!hNum2D || hNum2D->GetDimension() != 2) continue;

                            std::unique_ptr<TH1> hNum1D(
                                isProjX ? hNum2D->ProjectionX(Form("px_%s_%d_%s", var2d.c_str(), pad, trg.c_str()), 1, -1, "e")
                                        : hNum2D->ProjectionY(Form("py_%s_%d_%s", var2d.c_str(), pad, trg.c_str()), 1, -1, "e")
                            );
                            if (!hNum1D) continue;

                            // choose colours: default vs filtered
                            Color_t useCol = (!hasWeighted)
                                             ? ((filterIdx == 0) ? colBase[idx] : colFilt[idx])
                                             : ((filterIdx == (int)filters.size()-1) ? colBase[idx] : colFilt[idx]);

                            const bool isFiltered = (!hasWeighted) ? (filterIdx > 0)
                                                                   : (filterIdx < (int)filters.size()-1);

                            TH1* h_ratio;
                            TGraphAsymmErrors* g;

                            if (std::find(filters_weighted.begin(), filters_weighted.end(), fs) != filters_weighted.end()){ // filtered ratio is weighted binomial --> use TH1 instead of TGraphAsymmErrors
                                h_ratio = (TH1*)hNum1D->Clone(Form("%s_divided", hNum1D->GetName()));
                                h_ratio->Divide(hDen1D.get());

                                SetStyle(h_ratio, useCol, mkrBase[idx]);
                                
                                if (isFiltered) {
                                    h_ratio->SetLineStyle(2);
                                    h_ratio->SetMarkerStyle(mkrBase[idx] + 4);
                                }
                                
                                if (!firstDrawn) {
                                    h_ratio->Draw("E1");
                                    h_ratio->GetYaxis()->SetRangeUser(0., 1.05);
                                    h_ratio->GetXaxis()->SetTitle(hDen1D->GetXaxis()->GetTitle());
                                    h_ratio->GetYaxis()->SetTitle("#epsilon");
                                    
                                    gPad->RedrawAxis();
                                    gPad->Update(); // ensure pad limits are known

                                    double xmin = gPad->GetUxmin();
                                    double xmax = gPad->GetUxmax();
                                    
                                    TLine *line = new TLine(xmin, 1.0, xmax, 1.0);
                                    line->SetLineColor(kBlack);
                                    line->SetLineStyle(2); // dashed
                                    line->SetLineWidth(2);
                                    line->Draw("same");

                                    firstDrawn = true;
                                } else {
                                    h_ratio->Draw("E1 SAME");
                                }
                            } else{
                                g = new TGraphAsymmErrors();
                                g->BayesDivide(hNum1D.get(), hDen1D.get());

                                SetStyle(g, useCol, mkrBase[idx]);
                                
                                if (isFiltered) {
                                    g->SetLineStyle(2);
                                    g->SetMarkerStyle(mkrBase[idx] + 4);
                                }
                                
                                if (!firstDrawn) {
                                    g->Draw("APL");
                                    g->GetYaxis()->SetRangeUser(0., 1.);
                                    g->GetXaxis()->SetTitle(hDen1D->GetXaxis()->GetTitle());
                                    g->GetYaxis()->SetTitle("#epsilon");
                                    if (xy.first) adjustLogXRange(g, hDen1D.get());
                                    firstDrawn = true;
                                } else {
                                    g->Draw("PL SAME");
                                }
                            }

                            // legend text
                            std::string fs_label;
                            if (!fs.empty()) {
                                fs_label = (filt_suffix_to_label_map.count(fs)
                                           ? filt_suffix_to_label_map.at(fs)
                                           : (fs[0] == '_' ? fs.substr(1) : fs));
                            }
                            std::string label = trg;
                            if (!fs_label.empty()) label += ", " + fs_label;
                            if (std::find(filters_weighted.begin(), filters_weighted.end(), fs) != filters_weighted.end()){ // filtered ratio is weighted binomial --> use TH1 instead of TGraphAsymmErrors
                                leg->AddEntry(h_ratio, label.c_str(), "lep");
                            } else{
                                leg->AddEntry(g, label.c_str(), "lep");
                            }
                        } // loop over triggers
                        ++filterIdx;
                    } // loop over filters
                    leg->Draw();
                } // loop over pads

                // output (use projected 1D var name)
                std::string pngDir = data_dir + "trig_effcy_plots";
                for (const auto& fs : filter_suffix_list) {
                    pngDir += fs;
                }
                if (draw2mu4) pngDir += "_w_2mu4";

                makeDirIfNeeded(pngDir);
                std::string fn = pngDir + "/";
                fn += projVar + "_trig_effcy";
                for (const auto& fsu : filter_suffix_list) fn += fsu; // keep your naming convention
                fn += mode_to_png_suffix.at(mode) + ".png";
                c->SaveAs(fn.c_str());
            } // loop over signed & signal modes

            // ---------------- OpAndSignal (one file per filter) ----------------
            {
                StrVec opFilters = filter_suffix_list.empty() ? StrVec{""} : filter_suffix_list;
                for (const auto& fs : opFilters)
                {
                    if (!fs.empty() && !filt_to_mode_map.at(fs).at(SignalDrawingMode::OpAndSignal)) continue;

                    auto xy = xyFor(projVar);
                    TCanvas* c = new TCanvas(
                        Form("c_opandsig_proj_%s_%s_%s", var2d.c_str(), (isProjX ? "PX" : "PY"), fs.c_str()),
                        "", 700, 500
                    );
                    c->SetLogx(xy.first);
                    c->SetLogy(xy.second);

                    Rect Lbox = legendBoxFor1DVar(projVar, SignalDrawingMode::OpAndSignal, filter_suffix_list);
                    TLegend* L = new TLegend(Lbox[0], Lbox[1], Lbox[2], Lbox[3]);
                    L->SetBorderSize(0);

                    bool firstCurve = true;
                    int  colourIdx  = 0;  // 0 → Azure, 1 → Red

                    for (const std::string& trg : trg_list)
                    {
                        const std::string denomName = (mappedTrigger(fs, trg + "_denom") != "NONE")? mappedTrigger(fs, trg + "_denom") : mappedTrigger(fs, "mu4");

                        if (denomName.empty() || denomName == "NONE") continue;

                        // denominators (2D → project)
                        std::string den2D_s2 = use_sepr_for_op_and_signal? Form("h_%s_sign2_%s_sepr", var2d.c_str(), denomName.c_str()) : Form("h_%s_sign2_%s", var2d.c_str(), denomName.c_str());
                        std::string den2D_sel = Form("h_%s_single_b_%s",    var2d.c_str(), denomName.c_str());
                        
                        if (debug_mode){
                            std::cout << "Denominator name all opposite sign: " << den2D_s2 << std::endl;
                            std::cout << "Denominator name signal pairs: " << den2D_sel << std::endl;
                        }

                        TH2* hDen2D_s2  = getHist<TH2>(den2D_s2);
                        TH2* hDen2D_sel = getHist<TH2>(den2D_sel);
                        if (!hDen2D_s2 || !hDen2D_sel) continue;

                        std::unique_ptr<TH1> hDen1D_s2(
                            isProjX ? hDen2D_s2->ProjectionX(Form("px_%s_den_s2",  var2d.c_str()), 1, -1, "e")
                                    : hDen2D_s2->ProjectionY(Form("py_%s_den_s2",  var2d.c_str()), 1, -1, "e")
                        );
                        std::unique_ptr<TH1> hDen1D_sel(
                            isProjX ? hDen2D_sel->ProjectionX(Form("px_%s_den_sel", var2d.c_str()), 1, -1, "e")
                                    : hDen2D_sel->ProjectionY(Form("py_%s_den_sel", var2d.c_str()), 1, -1, "e")
                        );
                        if (!hDen1D_s2 || !hDen1D_sel) continue;

                        const std::string mapped = mappedTrigger(fs, trg);
                        if (mapped.empty() || mapped == "NONE") continue;

                        std::string num2D_s2 = use_sepr_for_op_and_signal? Form("h_%s_sign2_%s_sepr", var2d.c_str(), mapped.c_str()) : Form("h_%s_sign2_%s", var2d.c_str(), mapped.c_str());
                        std::string num2D_sel = Form("h_%s_single_b_%s",    var2d.c_str(), mapped.c_str());

                        TH2* hNum2D_s2  = getHist<TH2>(num2D_s2);
                        TH2* hNum2D_sel = getHist<TH2>(num2D_sel);
                        if (!hNum2D_s2 || !hNum2D_sel) continue;

                        std::unique_ptr<TH1> hNum1D_s2(
                            isProjX ? hNum2D_s2->ProjectionX(Form("px_%s_%s_s2",  var2d.c_str(), trg.c_str()), 1, -1, "e")
                                    : hNum2D_s2->ProjectionY(Form("py_%s_%s_s2",  var2d.c_str(), trg.c_str()), 1, -1, "e")
                        );
                        std::unique_ptr<TH1> hNum1D_sel(
                            isProjX ? hNum2D_sel->ProjectionX(Form("px_%s_%s_sel", var2d.c_str(), trg.c_str()), 1, -1, "e")
                                    : hNum2D_sel->ProjectionY(Form("py_%s_%s_sel", var2d.c_str(), trg.c_str()), 1, -1, "e")
                        );
                        if (!hNum1D_s2 || !hNum1D_sel) continue;

                        TH1* h_ratio_s2 = nullptr;
                        TH1* h_ratio_sel = nullptr;
                        TGraphAsymmErrors* g_s2 = nullptr;
                        TGraphAsymmErrors* g_sel = nullptr;

                        if (std::find(filters_weighted.begin(), filters_weighted.end(), fs) != filters_weighted.end()){ // filtered ratio is weighted binomial --> use TH1 instead of TGraphAsymmErrors
                            
                            h_ratio_s2 = (TH1*)hNum1D_s2->Clone(Form("%s_divided", hNum1D_s2->GetName()));
                            h_ratio_s2->Divide(hDen1D_s2.get());

                            SetStyle(h_ratio_s2,  colourIdx==0 ? kAzure+2 : kRed+1, colourIdx==0 ? 24 : 20);

                            h_ratio_sel = (TH1*)hNum1D_sel->Clone(Form("%s_divided", hNum1D_sel->GetName()));
                            h_ratio_sel->Divide(hDen1D_sel.get());

                            SetStyle(h_ratio_sel, colourIdx==0 ? kAzure+1 : kRed+2, colourIdx==0 ? 25 : 24, 2, 2);

                            if (firstCurve) {
                                h_ratio_sel->Draw("E1");  // draw dashed first so solid stays visible
                                h_ratio_sel->GetYaxis()->SetRangeUser(0., 1.05);
                                h_ratio_sel->GetXaxis()->SetTitle(hDen1D_sel->GetXaxis()->GetTitle());
                                h_ratio_sel->GetYaxis()->SetTitle("#epsilon");

                                gPad->RedrawAxis();
                                gPad->Update(); // ensure pad limits are known

                                double xmin = gPad->GetUxmin();
                                double xmax = gPad->GetUxmax();
                                
                                TLine *line = new TLine(xmin, 1.0, xmax, 1.0);
                                line->SetLineColor(kBlack);
                                line->SetLineStyle(2); // dashed
                                line->SetLineWidth(2);
                                line->Draw("same");

                                firstCurve = false;
                            } else {
                                h_ratio_sel->Draw("E1 SAME");
                            }
                            h_ratio_s2->Draw("PL SAME");
                        } else{
                            g_s2  = new TGraphAsymmErrors();
                            g_s2->BayesDivide(hNum1D_s2.get(),  hDen1D_s2.get());
                            SetStyle(g_s2,  colourIdx==0 ? kAzure+2 : kRed+1, colourIdx==0 ? 24 : 20);

                            g_sel = new TGraphAsymmErrors();
                            g_sel->BayesDivide(hNum1D_sel.get(), hDen1D_sel.get());
                            SetStyle(g_sel, colourIdx==0 ? kAzure+1 : kRed+2, colourIdx==0 ? 25 : 24, 2, 2);

                            if (firstCurve) {
                                g_sel->Draw("APL");  // draw dashed first so solid stays visible
                                g_sel->GetYaxis()->SetRangeUser(0., 1.);
                                g_sel->GetXaxis()->SetTitle(hDen1D_sel->GetXaxis()->GetTitle());
                                g_sel->GetYaxis()->SetTitle("#epsilon");
                                if (xy.first) adjustLogXRange(g_sel, hDen1D_sel.get());
                                firstCurve = false;
                            } else {
                                g_sel->Draw("PL SAME");
                            }
                            g_s2->Draw("PL SAME");
                        }
                        

                        // legend text
                        std::string fs_label;
                        if (!fs.empty()) {
                            fs_label = (filt_suffix_to_label_map.count(fs)
                                       ? filt_suffix_to_label_map.at(fs)
                                       : (fs[0] == '_' ? fs.substr(1) : fs));
                        }
                        std::string labS2  = trg + " (all op sign";
                        std::string labSel = trg + " (single b signal";
                        if (!fs_label.empty()) {
                            labS2  += ", " + fs_label;
                            labSel += ", " + fs_label;
                        }
                        labS2  += ")";
                        labSel += ")";
                        if (std::find(filters_weighted.begin(), filters_weighted.end(), fs) != filters_weighted.end()){ // filtered ratio is weighted binomial --> use TH1 instead of TGraphAsymmErrors
                            L->AddEntry(h_ratio_s2,  labS2.c_str(),  "lep");
                            L->AddEntry(h_ratio_sel, labSel.c_str(), "lep");
                        } else{                                
                            L->AddEntry(g_s2,  labS2.c_str(),  "lep");
                            L->AddEntry(g_sel, labSel.c_str(), "lep");
                        }

                        ++colourIdx;
                    }
                    L->Draw();

                    std::string outdir = data_dir + "trig_effcy_plots" + (fs.empty() ? "" : fs);
                    if (draw2mu4) outdir += "_w_2mu4";
                    
                    makeDirIfNeeded(outdir + "/op_and_sig");
                    std::string fn = Form("%s/op_and_sig/%s_trig_effcy%s%s.png",
                                          outdir.c_str(), projVar.c_str(), fs.c_str(),
                                          mode_to_png_suffix.at(SignalDrawingMode::OpAndSignal).c_str());
                    c->SaveAs(fn.c_str());
                }
            }
        }
    }
}

// ======================================================================
// Project 2nd-muon-kinematic 2D histograms to 1D (X/Y) and plot efficiencies
// Gives signed (by 2nd-muon charge) single-muon trigger efficiencies
// ======================================================================
void TrigEffPlotterPbPb::plot2Dto1DsingleMuonEffcyProj (const std::pair<const std::string, std::pair<bool,bool>>& kv){

    std::cout << "Call function plot2Dto1DsingleMuonEffcyProj" << std::endl;
    const std::string& var2d = kv.first;
    const bool doPX = kv.second.first;   // ProjectionX → onto X (right side of "..._vs_")
    const bool doPY = kv.second.second;  // ProjectionY → onto Y (left side of "..._vs_")

    const auto parsed = parseYvsX(var2d);
    const std::string& yVar = parsed.first;
    const std::string& xVar = parsed.second;


    std::vector<std::pair<std::string,std::string>> trig_pairs;

    if (isMB) { // Run-2 pp or MB mode → only 2mu4 / mu4
        trig_pairs = {{"2mu4","mu4"}};
    } else { // draw for both mu4noL1 and mu4 single-muon efficiency
        trig_pairs = {{"mu4_mu4noL1","mu4"},{"2mu4","mu4"}};
    }

    std::vector<std::string> filters_single_muon_effcy;
    if (isMB)   filters_single_muon_effcy = {"_sepr", "_MB"};
    else if (draw_single_b_signal)  filters_single_muon_effcy = {"_sepr", "", "_w_single_b_sig_sel"};
    else                            filters_single_muon_effcy = {"_sepr", ""};
    

    std::map<std::string, std::string> filter_to_label_map_single_muon_effcy;
    filter_to_label_map_single_muon_effcy["_sepr"] = "#Delta R > 0.8";
    filter_to_label_map_single_muon_effcy[""] = "no selection";
    filter_to_label_map_single_muon_effcy["_w_single_b_sig_sel"] = "Single-b signal";

    std::map<std::string, std::string> filter_to_label_map_single_muon_effcy_MB;
    filter_to_label_map_single_muon_effcy_MB["_sepr"] = "P[2mu4 | mu4 && #Delta R > 0.8]";
    filter_to_label_map_single_muon_effcy_MB["_MB"] = "P[mu4 | MB]";


    for (auto& trg_pair : trig_pairs){ // loop over configurations (single-muon triggers): mu4noL1 or mu4

        // --------------------- Step 1: retrieve 2D histogram & perform x/y projections ------------------------------------------
        std::map<std::pair<std::string, std::string>, std::array<std::array<std::array<TH1D*, 2>, 2>, 2>> h1D_proj;
        
        for (auto& ctr : ctr_bins){
            // map of ({filter, ctr}, numerator / denominator trigger, +/- 2nd-muon sign, x/y) to the projected TH1D*

            for (auto fs : filters_single_muon_effcy){ // loop over filters
                for (int trg_ind = 0; trg_ind < 2; ++trg_ind){ // trigger index (numerator / denominator)
                    std::string trg = (trg_ind == 0)? trg_pair.first : trg_pair.second;
                    
                    for (int sign_ind = 0; sign_ind < 2; ++sign_ind){ // charge sign index (+/-)
                        
                        std::string sign = "_sign" + std::to_string(sign_ind + 1);                    
                        std::string hname = hName(var2d, fs, trg , ctr + sign);

                        TH2* h2D;

                        if (fs != "_MB"){
                            h2D = getHist<TH2>(hname);
                        } else{
                            h2D = getHistMB<TH2>(hname);
                        }

                        if (!h2D) continue; // can't proceed with projection: skip current histogram

                        for (int proj_ind = 0; proj_ind < 2; ++proj_ind) { // x/y axis index
                            
                            const bool isProjX = (proj_ind == 0);
                            if ((isProjX && !doPX) || (!isProjX && !doPY)) continue;

                            const std::string projVar = isProjX ? (xVar.empty() ? var2d : xVar)
                                                                : (yVar.empty() ? var2d : yVar);
                            if (projVar.empty()) continue; // malformed name

                            h1D_proj[{ctr, fs}][trg_ind][sign_ind][proj_ind] = 
                                isProjX ? h2D->ProjectionX(Form("px_%s", hname.c_str()), 1, -1, "e")
                                        : h2D->ProjectionY(Form("py_%s", hname.c_str()), 1, -1, "e");
                        }            
                    } // end loop over signs
                } // end loop over triggers (numerator / denominator)
            } // end loop over filters
        } // end loop over ctrs

        // --------------------- Step 2: drawing & saving canvases ------------------------------------------
        Color_t colors[3] = {kRed+1, kAzure+2, kGreen+2};
        Style_t mkrs[3] = {24, 20, 22};

        std::vector<Color_t> ctr_colors = {kBlack, kBlue, kMagenta, kGreen+2, kOrange + 2, kRed, kGray, kBrown};

        for (int proj_ind = 0; proj_ind < 2; ++proj_ind) // x, y projection --> fixes the variable of interest
        {
            const bool isProjX = (proj_ind == 0);
            if ((isProjX && !doPX) || (!isProjX && !doPY)) continue;

            const std::string projVar = isProjX ? (xVar.empty() ? var2d : xVar)
                                                : (yVar.empty() ? var2d : yVar);
            if (projVar.empty()) continue; // malformed name

            auto xy = xyFor(projVar);

            for (auto& ctr : ctr_bins){
                // ----------- if not MB mode: draw first canvas: +/- comparison for well-separated pairs ------------------------------------------
                if (!isMB){
                    // canvas
                    std::string cname = Form("c_%s%s_%s_charge_sign_compr",
                                             projVar.c_str(), ctr.c_str(), trg_pair.first.c_str());
                    TCanvas* c = new TCanvas(cname.c_str(), "", 700, 500);

                    Rect box = legendBoxFor1DVar(projVar, SignalDrawingMode::OpAndSignal, filter_suffix_list);

                    c->cd();
                    gPad->SetLogx(xy.first);
                    gPad->SetLogy(xy.second);

                    TLegend* leg = new TLegend(box[0], box[1], box[2], box[3]);
                    leg->SetBorderSize(0);

                    for (int sign_ind = 0; sign_ind < 2; ++sign_ind){ // charge sign index (+/-)
                        
                        // draw the trigger-efficiency graphs
                        auto g = new TGraphAsymmErrors();

                        // h1D_proj[{ctr,fs}][trg_ind][sign_ind][proj_ind]
                        if (!h1D_proj[{ctr, "_sepr"}][0][sign_ind][proj_ind] || !h1D_proj[{ctr, "_sepr"}][1][sign_ind][proj_ind]){
                            std::cerr << "The well-separated single-muon trigger-efficiency histograms cannot be found for the variable " << projVar << std::endl;
                            continue;
                        }
                        
                        g->BayesDivide(h1D_proj[{ctr, "_sepr"}][0][sign_ind][proj_ind], h1D_proj[{ctr, "_sepr"}][1][sign_ind][proj_ind]);
                        
                        SetStyle(g, colors[sign_ind], mkrs[sign_ind]);

                        if (sign_ind == 0) {
                            g->Draw("APL");
                            g->GetYaxis()->SetRangeUser(0., 1.);
                            g->GetXaxis()->SetTitle(h1D_proj[{ctr, "_sepr"}][1][sign_ind][proj_ind]->GetXaxis()->GetTitle());
                            g->GetYaxis()->SetTitle("#epsilon");
                            if (xy.first) adjustLogXRange(g, h1D_proj[{ctr, "_sepr"}][1][sign_ind][proj_ind]);
                        } else {
                            g->Draw("PL SAME");
                        }
                        
                        // legend entry
                        std::string sign_label = (sign_ind == 0)? "#mu^{+}" : "#mu^{-}";
                        leg->AddEntry(g, sign_label.c_str(), "lep");
                    } // end loop over charge signs (lines)

                    leg->Draw();

                    std::string outdir = data_dir + "trig_effcy_plots";
                    makeDirIfNeeded(outdir + "/single_muon_effcy");
                    std::string fn = Form("%s/single_muon_effcy/%s_trig_effcy%s_%s_charge_sign_compr.png",
                                          outdir.c_str(), projVar.c_str(), ctr.c_str(), target_dimuon_trigger_to_single_muon_map[trg_pair.first].c_str());
                    c->SaveAs(fn.c_str());

                } // end if !MB statement & first canvas drawing

                // ----------- draw second canvas: for +/- 2nd muon, draw EITHER (!isMB) comparison between well-separated, no-selection & signal pairs ------------------------------------------
                // ------------------------------------------------------ OR (isMB) comparison between well-separated, no-selection & signal pairs ------------------------------------------
                {
                    // canvas
                    std::string cname = Form("c_%s%s_%s_signed_compr",
                                             projVar.c_str(), ctr.c_str(), trg_pair.first.c_str());
                    TCanvas* c = new TCanvas(cname.c_str(), "", 1200, 500);
                    c->Divide(2,1);

                    Rect box = legendBoxFor1DVar(projVar, SignalDrawingMode::Signed, filter_suffix_list);

                    for (int sign_ind = 0; sign_ind < 2; ++sign_ind){ // charge sign index (+/-)
                        c->cd(sign_ind + 1);
                        gPad->SetLogx(xy.first);
                        gPad->SetLogy(xy.second);

                        TLegend* leg = new TLegend(box[0], box[1], box[2], box[3]);
                        leg->SetBorderSize(0);

                        for (int filter_ind = 0; filter_ind < filters_single_muon_effcy.size(); filter_ind++){
                            std::string filter = filters_single_muon_effcy.at(filter_ind);

                            // draw the trigger-efficiency graphs
                            auto g = new TGraphAsymmErrors();

                            // h1D_proj[{ctr,fs}][trg_ind][sign_ind][proj_ind]
                            if (!h1D_proj[{ctr,filter}][0][sign_ind][proj_ind] || !h1D_proj[{ctr,filter}][1][sign_ind][proj_ind]){
                                std::cerr << "The well-separated single-muon trigger-efficiency histograms cannot be found for the variable " << projVar << std::endl;
                                continue;
                            }
                            
                            g->BayesDivide(h1D_proj[{ctr,filter}][0][sign_ind][proj_ind], h1D_proj[{ctr,filter}][1][sign_ind][proj_ind]);
                            
                            SetStyle(g, colors[filter_ind], mkrs[filter_ind]);

                            if (filter_ind == 0) {
                                g->Draw("APL");
                                g->GetYaxis()->SetRangeUser(0., 1.);
                                g->GetXaxis()->SetTitle(h1D_proj[{ctr,filter}][1][sign_ind][proj_ind]->GetXaxis()->GetTitle());
                                g->GetYaxis()->SetTitle("#epsilon");
                                if (xy.first) adjustLogXRange(g, h1D_proj[{ctr,filter}][1][sign_ind][proj_ind]);
                            } else {
                                g->Draw("PL SAME");
                            }
                            
                            // legend entry
                            std::string label = (isMB)? filter_to_label_map_single_muon_effcy_MB[filter] : filter_to_label_map_single_muon_effcy[filter];
                            leg->AddEntry(g, label.c_str(), "lep");
                        } // end loop over filters (lines)

                        leg->Draw();
                    } // end loop over 2nd-muon signs (pads)

                    std::string compr_type = (isMB)? "MB_2mu4" : "selection";
                    std::string outdir = data_dir + "trig_effcy_plots";
                    makeDirIfNeeded(outdir + "/single_muon_effcy");
                    std::string fn = Form("%s/single_muon_effcy/%s_trig_effcy%s_%s_%s_compr.png",
                                          outdir.c_str(), projVar.c_str(), ctr.c_str(), target_dimuon_trigger_to_single_muon_map[trg_pair.first].c_str(), compr_type.c_str());
                    c->SaveAs(fn.c_str());
                } // end second canvas drawing
            } // end loop over ctr bins

            // --------------------- Third canvas: ctr comparison (sum over mu+, mu-) ------------------------------------------
            // canvas
            std::string cname = Form("c_%s_%s_ctr_compr",
                                     projVar.c_str(), trg_pair.first.c_str());
            TCanvas* c = new TCanvas(cname.c_str(), "", 700, 500);

            Rect box = legendBoxFor1DVar(projVar, SignalDrawingMode::OpAndSignal, filter_suffix_list);

            c->cd();
            gPad->SetLogx(xy.first);
            gPad->SetLogy(xy.second);

            TLegend* leg = new TLegend(box[0], box[1], box[2], box[3]);
            leg->SetBorderSize(0);

            for (int ctr_ind = 0; ctr_ind < ctr_bins.size(); ++ctr_ind){
                
                auto ctr = ctr_bins.at(ctr_ind);

                // h1D_proj[{ctr, fs}][trg_ind][sign_ind][proj_ind]
                for (int sign_ind = 0; sign_ind < 2; ++sign_ind){ // charge sign index (+/-)
                    if (!h1D_proj[{ctr, "_sepr"}][0][sign_ind][proj_ind] || !h1D_proj[{ctr, "_sepr"}][1][sign_ind][proj_ind]){
                        std::cerr << "The well-separated single-muon trigger-efficiency histograms cannot be found for the variable " << projVar << ", " << ctr << std::endl;
                        return;
                    }
                }

                // h1D_proj[{ctr, fs}][trg_ind][sign_ind][proj_ind]
                TH1D* hnum_muplus = h1D_proj[{ctr, "_sepr"}][0][0][proj_ind];
                TH1D* hnum_sumsign = (TH1D*) hnum_muplus->Clone(Form("%s_added", hnum_muplus->GetName()));
                hnum_sumsign->Add(h1D_proj[{ctr, "_sepr"}][0][1][proj_ind]);

                TH1D* hdom_muplus = h1D_proj[{ctr, "_sepr"}][1][0][proj_ind];
                TH1D* hdom_sumsign = (TH1D*) hdom_muplus->Clone(Form("%s_added", hdom_muplus->GetName()));
                hdom_sumsign->Add(h1D_proj[{ctr, "_sepr"}][1][1][proj_ind]);
                
                // draw the trigger-efficiency graphs
                auto g = new TGraphAsymmErrors();
                
                g->BayesDivide(hnum_sumsign, hdom_sumsign);
                
                SetStyle(g, ctr_colors.at(ctr_ind), 20);

                if (ctr_ind == 0) {
                    g->Draw("APL");
                    g->GetYaxis()->SetRangeUser(0., 1.05);
                    g->GetXaxis()->SetTitle("Centrality [%]");
                    g->GetYaxis()->SetTitle("#epsilon");
                } else {
                    g->Draw("PL SAME");
                }
                
                // legend entry 
                std::string ctr_label = std::to_string(ctr_bin_edges.at(ctr_ind)) + "-" + std::to_string(ctr_bin_edges.at(ctr_ind+1)) + "%";
                leg->AddEntry(g, ctr_label.c_str(), "lep");
            } // end loop over ctrs (lines)

            leg->Draw();

            std::string outdir = data_dir + "trig_effcy_plots";
            makeDirIfNeeded(outdir + "/single_muon_effcy/ctr_dep");
            std::string fn = Form("%s/single_muon_effcy/ctr_dep/%s_trig_effcy_%s_ctr_compr.png",
                                  outdir.c_str(), projVar.c_str(), target_dimuon_trigger_to_single_muon_map[trg_pair.first].c_str());
            c->SaveAs(fn.c_str());
        } // end loop over x/y projection (determines 1D variable to plot efficiencies for)
    } // end loop over (target, supporting) trigger pairs
} // end function plot2Dto1DsingleMuonEffcyProj

