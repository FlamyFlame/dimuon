#ifndef ParamsSet_h
#define ParamsSet_h

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <utility> 
#include <functional>

struct AxisInfo{ // a simplified version of axis information - for constructing TH1DModel
    int        	nbins		= 1;         ///<  Number of bins
    double     	min			= 0.;          ///<  Low edge of first bin
    double     	max			= 1.;          ///<  Upper edge of last bin
    const double*     bin_edges		= nullptr;         ///<  Bin edges array in X
};

class ParamsSet{
public:
	static const unsigned int nCtrBins=6; // number of coarse bins; each can be studied by itself
   	static const unsigned int nPtBins=10;
   	static const unsigned int nPythiaKinRanges=5;
   	static const unsigned int nSigns=2;

   	static const unsigned int ndRselcs=3;
   	static const unsigned int ndphiselcs=3;
   	static const unsigned int ndphiRegions=3;
   	static const unsigned int ndetaRegions=2;

   	static const unsigned int nEffCorr=2; // 0: no efficiency correction; 1: having effiency correction
   	static const unsigned int nPhotoProdCuts=2; //0: no gap cut; 1: having gap cut
   	static const unsigned int nGapCuts=2; //0: no gap cut; 1: having gap cut
  	static const unsigned int CtrStep = 5;
  	static const unsigned int nCtrIntvls = 20; // number of small intervals; the intervals not to be studied by themselves, but to allow easy combinations

	static std::vector<std::function<bool(float)>> dphi_cut_funcs;
	static std::vector<std::function<bool(float)>> deta_cut_funcs;

  	// float pTbins[nPtBins] = {4.,5.,6.,7.,8.,9.,10.,12.,15.,20.};
	// int scaleFactorCtrs[nCtrBins] = {1,1,2,3,3};
	// std::vector<float> ctrbins = {0, 5, 10, 20, 30, 50, 80};
    std::vector<double> pT_bins_40;
    std::vector<double> pT_bins_8;
    std::vector<double> pT_bins_60;
    std::vector<double> pT_bins_80;
	
	static std::vector<float> pTbins;
	static std::vector<double> pTHatbins_pythia;
    static std::vector<int> ctrbins;

    std::vector<std::string> pt_bin_labels;
    std::vector<std::string> pt_bin_exprs;
    std::vector<std::string> pt_titles = {"#bar{p}_T 4-5 GeV", "#bar{p}_T 5-6 GeV", "#bar{p}_T 6-7 GeV", "#bar{p}_T 7-8 GeV", "#bar{p}_T 8-9 GeV", "#bar{p}_T 9-10 GeV", "#bar{p}_T 10-12 GeV", "#bar{p}_T 12-15 GeV", "#bar{p}_T 15-20 GeV", "#bar{p}_T > 20 GeV"};
    
    std::vector<std::string> ctr_bin_labels;
    std::vector<std::string> ctr_bin_exprs;
	std::vector<std::string> ctr_titles = {"Centrality 0-5", "Centrality 5-10", "Centrality 10-20", "Centrality 20-30", "Centrality 30-50", "Centrality 50-80"};
	std::vector<std::string> ctrNpp_titles = {"Centrality 0-5", "Centrality 5-10", "Centrality 10-20", "Centrality 20-30", "Centrality 30-50", "Centrality 50-80", "pp"};


	std::string signs[nSigns] = {"ss","op"};
	std::vector<std::string> sign_labels = {"_ss","_op"};
	std::vector<std::string> sign_titles = {"same sign", "opposite sign"};

	std::vector<std::string> deta_cut_labels = {"_deta_lt0_8", "_deta_gt0_8"};
	std::vector<std::string> deta_titles = {"#Delta #eta < 0.8","#Delta #eta #geq 0.8"};

	std::vector<std::string> dphi_cut_labels = {"_near", "_away", "_dphi_lt0_6"};
	std::vector<std::string> dphi_titles = {"#Delta #phi < #pi/2","#Delta #phi #geq #pi/2", "#Delta #phi < 0.6"};

	std::vector<std::string> gapcut_labels = {"_nogapcut", "_wgapcut"};
	std::vector<std::string> gapcut_titles = {"no gap cut", "with gap cut"};

	std::vector<std::string> photocut_labels = {"_nophotocut", "_wphotocut"};
	std::vector<std::string> photocut_titles = {"no photoproduction cut", "with photoproduction cut"};


   	float deltaP_overP_thrsh; // single-muon dP/P cut
   	float deltaP_overP_max;
   	float deltaP_overP_step;
   	int deltaP_overP_nbins;
   	
   	float d0cut; // single-muon |d0| cut (unit: mm)
   	float z0cut; // single-muon |z0 sin(theta)| cut (unit: mm)

   	static float deltaR_thrsh[ndRselcs];
   	static float deltaR_thrsh_zoomin;
   	float 	deltaR_step = 0.01;
   	int deltaR_nbins[ndRselcs];

   	// float eta_gap_cut1 = 0.1;
   	static constexpr float eta_gap_cut1 = 0.135;
   	// float eta_gap_cut2[2] = {1.05,1.29};

	// // no log
   	// float minv_min[nSigns] = {0.2,1.06};
   	// float minv_max[ndRselcs] = {10,15,100};
   	// int minv_nbins[ndRselcs] = {100,150,1000};

   	float minv_max[ndRselcs] = {10,15,60};
   	// static constexpr int minv_nbins[ndRselcs] = {25,40,120};
   	static constexpr int minv_nbins[ndRselcs] = {20,30,120};
   		// the C++ standard does not specifiy how floating point should be implemented and is left to the processor. 
   		// To get around this and other limitations constexpr was introduced.
   	// float minv_logpow[nSigns][ndRselcs];
   	// std::vector<float> minv_bins[nSigns][ndRselcs];

   	std::vector<std::array<float,2>> minv_cuts;
   	std::vector<std::array<float,2>> minv_cuts_v2;
   	static std::vector<std::array<float,2>> charge_eta_gap_cuts;

   	float minv_upper = 60;

  	static double PI;

  	static AxisInfo ptlead_axis;
	static AxisInfo dr_axis;
	static AxisInfo dphi_axis;


  	static const int npt_bins = 50;
  	float ptlogpow = 0.02194;
  	float ptmax = 50;
  	double pTBins[npt_bins+1];
	
  	static const int npairPT_bins = 40;
  	float pairPTlogpow[nSigns][ndRselcs];
  	float pairPTmax = 80;
  	double pairPTBins[nSigns][ndRselcs][npairPT_bins+1];

  	static const int n_hq_pt_bins = 40;
  	float hq_ptlogpow = 0.034;
  	float hq_ptmax = 110;
  	double hq_pTBins[n_hq_pt_bins+1];

  	static const int n_hq_minv_bins = 40;
  	float hq_minvlogpow = 0.049;
  	float hq_minvmax = 220;
  	double hq_minvBins[n_hq_minv_bins+1];
  	
  	ParamsSet();
  	~ParamsSet(){}
    static std::vector<double> makeEtaTrigEffcyBinning(int rebin_factor = 1);
  	void fillLogBinningArray(std::vector<double>& bins, int nBins, double low, double high);
  	template <typename T>
  	std::string write_single_bin_expr (std::string kinvar, T a, T b);
  	template <typename T>
	std::string write_single_bin_label (std::string kinvar, T a, T b);
	template <typename T>
	void write_cut_label_strs (std::string kinvar_expr, std::string kinvar_label, const std::vector<T> & var_bin_bdrys, std::vector<std::string> & var_bin_exprs, std::vector<std::string> & var_bin_labels, bool open_ended);
};

void ParamsSet::fillLogBinningArray(std::vector<double>& bins, int nBins, double low, double high) {
    bins.clear();

    double logLow = std::log10(low);
    double logHigh = std::log10(high);
    double logStep = (logHigh - logLow) / nBins;

    for (int i = 0; i <= nBins; ++i) {
        bins.push_back(std::pow(10, logLow + i * logStep));
    }
}

#include <vector>
#include <algorithm>
#include <cmath>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

std::vector<double> ParamsSet::makeEtaTrigEffcyBinning(int rebin_factor)
{
    const double minEdge = -2.4;
    const double maxEdge =  2.4;

    // Base steps (before rebinning)
    const double ultraFineBase = 0.01;  // adjust to 0.005 if you really want that
    const double fineBase      = 0.02;
    const double coarseStep    = 0.10;

    const double ultraFineStep = ultraFineBase * rebin_factor;
    const double fineStep      = fineBase      * rebin_factor;

    // Fine ranges (to be merged)
    std::vector<std::pair<double,double>> fineRanges = {
        {-1.3,  1.0},
        {-0.8,  1.4},
        { 2.2,  2.4}
    };

    // Ultra-fine range
    const std::pair<double,double> ultraFineRange = {-0.2, 0.2};

    // --- merge overlapping fine ranges ---
    std::sort(fineRanges.begin(), fineRanges.end());
    std::vector<std::pair<double,double>> merged;
    for (auto &r : fineRanges) {
        if (merged.empty() || r.first > merged.back().second) {
            merged.push_back(r);
        } else {
            merged.back().second = std::max(merged.back().second, r.second);
        }
    }

    std::vector<double> edges;
    edges.reserve(200); // just to avoid reallocations
    const double eps = 1e-10;

    auto push_edge = [&](double x) {
        // round to kill FP noise
        x = std::round(x * 1e12) / 1e12;
        if (edges.empty() || std::fabs(edges.back() - x) > eps) {
            edges.push_back(x);
        }
    };

    // Add segment [start, end] with given step
    auto add_segment = [&](double start, double end, double step) {
        if (end <= start + eps) {
            push_edge(start);
            return;
        }
        push_edge(start);
        double x = start;
        while (x + step < end - eps) {
            x += step;
            push_edge(x);
        }
        push_edge(end);
    };

    double x = minEdge;

    // Walk over all merged fine ranges in order
    for (const auto &fr : merged) {
        const double a = fr.first;
        const double b = fr.second;

        // 1) Coarse region up to the start of this fine range
        if (x < a - eps) {
            add_segment(x, a, coarseStep);
            x = a;
        }

        // 2) Inside the fine range [a, b]
        double uf_lo = std::max(a, ultraFineRange.first);
        double uf_hi = std::min(b, ultraFineRange.second);

        // 2a) Fine region before ultra-fine part
        if (x < uf_lo - eps) {
            add_segment(x, uf_lo, fineStep);
            x = uf_lo;
        }

        // 2b) Ultra-fine region (overlap)
        if (uf_lo < uf_hi - eps) {
            add_segment(x, uf_hi, ultraFineStep);
            x = uf_hi;
        }

        // 2c) Fine region after ultra-fine part within this fine range
        if (x < b - eps) {
            add_segment(x, b, fineStep);
            x = b;
        }
    }

    // 3) After the last fine range: coarse step to maxEdge
    if (x < maxEdge - eps) {
        add_segment(x, maxEdge, coarseStep);
    }

    return edges;
}


double ParamsSet::PI = acos(-1.0);

float ParamsSet::deltaR_thrsh[ndRselcs] = {0.8,1.2,5.75};
float ParamsSet::deltaR_thrsh_zoomin = 0.8;

std::vector<std::array<float,2>> ParamsSet::charge_eta_gap_cuts = {{0.56,0.67}, {1.064,1.29}, {-1.29,-1.12}};
std::vector<float> ParamsSet::pTbins = {4.,5.,6.,7.,8.,9.,10.,12.,15.,20.};
std::vector<int> ParamsSet::ctrbins = {0, 5, 10, 20, 30, 50, 80};
// std::vector<double> ParamsSet::pTHatbins_pythia = {5, 10, 25, 60, 120, 3200};
std::vector<double> ParamsSet::pTHatbins_pythia = {4.5, 10, 25, 60, 120, 3200}; // lower to 4.5 since a few events have pTHat < 5GeV (leakage)

AxisInfo ParamsSet::ptlead_axis = {100, 0, 30., nullptr};
AxisInfo ParamsSet::dr_axis = {100, 0, ParamsSet::deltaR_thrsh[2], nullptr};
AxisInfo ParamsSet::dphi_axis = {32, -ParamsSet::PI/2., ParamsSet::PI * 3. / 2., nullptr};

std::vector<std::function<bool(float)>> ParamsSet::dphi_cut_funcs {[](float x){return abs(x) < 1.5708;}, [](float x){return abs(x) >= 1.5708;}, [](float x){return abs(x) < 0.6;}};
std::vector<std::function<bool(float)>> ParamsSet::deta_cut_funcs {[](float x){return abs(x) < 0.8;}, [](float x){return abs(x) >= 0.8;}};
const unsigned int ParamsSet::nCtrBins; // number of coarse bins; each can be studied by itself
const unsigned int ParamsSet::nPtBins;
const unsigned int ParamsSet::nPythiaKinRanges;
const unsigned int ParamsSet::nSigns;
const unsigned int ParamsSet::ndRselcs;
const unsigned int ParamsSet::ndphiselcs;
const unsigned int ParamsSet::ndphiRegions;
const unsigned int ParamsSet::ndetaRegions;
const unsigned int ParamsSet::nPhotoProdCuts; //0: no gap cut; 1: having gap cut
const unsigned int ParamsSet::nGapCuts; //0: no gap cut; 1: having gap cut
const unsigned int ParamsSet::CtrStep;
const unsigned int ParamsSet::nCtrIntvls; // number of small intervals; the intervals not to be studied by themselves, but to allow easy combinations

template <typename T>
std::string ParamsSet::write_single_bin_expr (std::string kinvar, T a, T b){
    return kinvar + " >= " + std::to_string(a) + " && " + kinvar + " < " + std::to_string(b);
}

template <typename T>
std::string ParamsSet::write_single_bin_label (std::string kinvar, T a, T b){
    return "_" + kinvar + std::to_string(static_cast<int>(a)) + "to" + std::to_string(static_cast<int>(b));
}

template <typename T>
void ParamsSet::write_cut_label_strs (std::string kinvar_expr, std::string kinvar_label, const std::vector<T> & var_bin_bdrys, std::vector<std::string> & var_bin_exprs, std::vector<std::string> & var_bin_labels, bool open_ended){
    var_bin_exprs.clear();
    var_bin_labels.clear();
    for (typename std::vector<T>::const_iterator it = var_bin_bdrys.begin() + 1; it < var_bin_bdrys.end(); it++){
        var_bin_exprs.push_back(write_single_bin_expr(kinvar_expr, *(it-1), *it));
        var_bin_labels.push_back(write_single_bin_label(kinvar_label, *(it-1), *it));
    }

    if (open_ended){
        var_bin_exprs.push_back(kinvar_expr + " >= " + std::to_string(var_bin_bdrys[var_bin_bdrys.size()-1]));
        var_bin_labels.push_back("_" + kinvar_label + "_above" + std::to_string(static_cast<int>(var_bin_bdrys[var_bin_bdrys.size()-1])));
    }
}

ParamsSet::ParamsSet(){
	deltaP_overP_thrsh = 0.12;
 	deltaP_overP_max = 0.12 * sqrt(2);
 	deltaP_overP_step = 0.002;
 	deltaP_overP_nbins = static_cast<int>(deltaP_overP_max/deltaP_overP_step);

	d0cut = 2; // 2mm
	z0cut = 2; // 2mm

// ------------------------------------ pt & ctr labels & expressions ------------------------------------
	    
    // write_cut_label_strs("(m1.pt + m2.pt) / 2.", "ptavg", pTbins, pt_bin_exprs, pt_bin_labels, true);
    write_cut_label_strs("pt_avg", "ptavg", pTbins, pt_bin_exprs, pt_bin_labels, true);
    write_cut_label_strs("avg_centrality", "ctr", ctrbins, ctr_bin_exprs, ctr_bin_labels, false);

  	for (unsigned int idr = 0; idr < ndRselcs; idr++){
    	deltaR_nbins[idr] = static_cast<int>(deltaR_thrsh[idr]/deltaR_step);
  	}

	for(int i = 0; i <= npt_bins; i++){
    	pTBins[i] = ptmax * pow(10.0, (static_cast<float>(i-npt_bins))*ptlogpow);
  	}

	for(int i = 0; i <= n_hq_pt_bins; i++){
    	hq_pTBins[i] = hq_ptmax * pow(10.0, (static_cast<float>(i-n_hq_pt_bins)) * hq_ptlogpow);
  	}

	for(int i = 0; i <= n_hq_minv_bins; i++){
    	hq_minvBins[i] = hq_minvmax * pow(10.0, (static_cast<float>(i-n_hq_minv_bins)) * hq_minvlogpow);
  	}

	pairPTlogpow[0][0] = 0.0198;
	pairPTlogpow[0][1] = 0.0198;
	// pairPTlogpow[0][2] = 0.052;
	pairPTlogpow[0][2] = 0.078;
	pairPTlogpow[1][0] = 0.0198;
	pairPTlogpow[1][1] = 0.0198;
	pairPTlogpow[1][2] = 0.068;


	// ={{0.0198,0.0198,0.052},{0.0198,0.0198,0.082}};
  	for (int isign = 0; isign < nSigns; isign++){
  		for (int idr = 0; idr < ndRselcs; idr++){
  			for(int ipt = 0; ipt <= npairPT_bins; ipt++){
    			pairPTBins[isign][idr][ipt] = pairPTmax * pow(10.0, (static_cast<float>(ipt - npairPT_bins)) * pairPTlogpow[isign][idr]);
  			}
  		}
  	}

    fillLogBinningArray(pT_bins_40,  18, 4.0, 40.0);  // 18 log bins from 4 to 40 GeV
    fillLogBinningArray(pT_bins_8 ,  20, 4.0, 8.0 );  // 15 log bins from 4 to 	8 GeV
    fillLogBinningArray(pT_bins_60,  20, 8.0, 60.0);  // 20 log bins from 8 to 60 GeV
    fillLogBinningArray(pT_bins_80,  12, 8.0, 80.0);  // 12 log bins from 8 to 200 GeV


  	// old set of minv cut
  	minv_cuts.push_back({0,1.06});
   	minv_cuts.push_back({2.9,3.3});
   	minv_cuts.push_back({3.55,3.8});
   	// minv_cuts.push_back({9.08,9.8});
   	minv_cuts.push_back({9.08,10.5}); // previously 9 - 9.8

   	// new set of minv cut
  	minv_cuts_v2.push_back({0.,0.6});
  	minv_cuts_v2.push_back({0.72,0.85});
  	minv_cuts_v2.push_back({0.94,1.06});
   	minv_cuts_v2.push_back({2.9,3.3});
   	minv_cuts_v2.push_back({3.55,3.8});
   	// minv_cuts_v2.push_back({9.08,9.8});
   	minv_cuts_v2.push_back({9.08,10.5}); // previously 9 - 9.8
}

#endif
