#pragma once
#include "PlotCommonConfig.h"
// ------------------------------------------------------------
// Shared helpers
// ------------------------------------------------------------
template<typename T>
static void SetStyle(T* g, Color_t col, Style_t m, Width_t w = 2, Style_t ls = 1)
{
    if (!g) return;
    g->SetMarkerColor(col);
    g->SetLineColor(col);
    g->SetMarkerStyle(m);
    g->SetLineWidth(w);
    g->SetLineStyle(ls);
}


static bool fileExists(const std::string& path)
{
    // AccessPathName returns kTRUE if it CANNOT be accessed
    return !gSystem->AccessPathName(path.c_str());
}

static void makeDirIfNeeded(const std::string& dir)
{
    if (!fileExists(dir)) gSystem->mkdir(dir.c_str(), kTRUE /*recursive*/);
}

void adjustLogXRange(TGraphAsymmErrors* g, const TH1* hRef)
{
    if (!g || !hRef) return;
    double xmin = hRef->GetXaxis()->GetBinLowEdge(1);
    double xmax = hRef->GetXaxis()->GetBinLowEdge(hRef->GetNbinsX() + 1);
    g->GetXaxis()->SetRangeUser(xmin, xmax);
}

std::string StripTruthPrefix(const std::string& var)
{
    if (var.rfind("truth_", 0) == 0)
        return var.substr(6);

    return var;
}


bool LogAx(const std::string& var, const PlotCommonConfig& cfg){ // whether to set axis on log scale
    std::string var_stripped = StripTruthPrefix(var);
    bool log_ax = (cfg.log_map.find(var_stripped) != cfg.log_map.end())? cfg.log_map.at(var_stripped) : false;
    return log_ax;
}
