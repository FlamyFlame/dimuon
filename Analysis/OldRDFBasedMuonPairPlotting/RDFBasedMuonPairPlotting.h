#ifndef __RDFBasedMuonPairPlotting_h__
#define __RDFBasedMuonPairPlotting_h__

#include <map>
#include <vector>
#include <string>
#include <utility> 
#include <numeric>
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "time.h"
#include <ROOT/RLogger.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>


struct Hist1DInfo{
	std::string var_expr;
	std::string h1d_name_prefix;
	std::string hist1d_titles;
	AxisInfo xaxis;
};


struct Hist2DInfo{

std::string xvar_expr;
std::string yvar_expr;
std::string h2d_name_prefix;
std::string hist2d_titles;
AxisInfo xaxis;
AxisInfo yaxis;
};



class RDFBasedMuonPairPlotting{
protected:
	ParamsSet pms;

	std::string_view tree_ss = "muon_pair_tree_sign1";
	std::string_view tree_op = "muon_pair_tree_sign2";

	std::vector<Hist1DInfo*> hist1Dinfo_list;
		
	std::vector<Hist2DInfo*> hist2Dinfo_list;

	std::vector<std::vector<TH1D>>* hist1Ds; // first "vector" for the different batches/input files (size=1 for data)
	std::vector<std::vector<TH2D>>* hist2Ds; // second "vector" for the different variables to be plotted & the different bins

	std::vector<ROOT::RDF::RResultPtr< ::TH1D>>* hist1d_rresultptrs_cur_file;
	std::vector<ROOT::RDF::RResultPtr< ::TH2D>>* hist2d_rresultptrs_cur_file;

	std::vector<std::string> input_files {};
	int nInputFiles;
	std::string output_file;
	
    int nDR_bins =                          80;
    int nDphi_bins =                        64;
    int neta_bins =                         50;
    int nDeta_bins =                        100;
    int npair_y_bins =                      45;
    int npt_asym_bins =                     50;
    int npair_pt_ptlead_ratio_bins =        50;
    int nminv_bins_linear = 50;
    int npair_pT_bins_linear = 50;
    int npT_lead_bins_linear = 50;
    static const int nminv_bins_log = 40;
    float minv_logpow[ParamsSet::nSigns] = {0.062, 0.045};
	double minv_bins_log[ParamsSet::nSigns][nminv_bins_log+1];

public:
	RDFBasedMuonPairPlotting(){

	    float minv_max = 60;
	    for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
	        for(int iminv = 0; iminv <= nminv_bins_log; iminv++){
	            minv_bins_log[ksign][iminv] = minv_max * pow(10.0, (static_cast<float>(iminv - nminv_bins_log)) * minv_logpow[ksign]);
	        }
	    }

		auto title_str = [](std::string kin){return ";" + kin + ";#frac{d#sigma}{d" + kin + "}";};
		auto title_str_w_unit = [](std::string kin, std::string unit){return ";" + kin + " [" + unit + "];#frac{d#sigma}{d" + kin + "}";};

		std::cout << "constructor for base class called" << std::endl;
		std::cout << "The base class constructor initiates and collects (into a vector) 1D & 2D histograms one needs to make for all data samples." << std::endl;
		std::cout << "If needed, one can add additional 1D & 2D histograms of interest to each data sample in the corresponding child class constructors." << std::endl;

		hist1Dinfo_list.push_back(new Hist1DInfo({"m1.pt", "h_pt_lead", title_str_w_unit("p_{T}^{lead}", "GeV"), ParamsSet::ptlead_axis}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"dr", "h_dr", title_str("#Delta R"), ParamsSet::dr_axis}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"dphi_winded", "h_dphi", title_str("#Delta #phi"), ParamsSet::dphi_axis}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"pair_y", "h_pair_y", title_str("y_{pair}"), {npair_y_bins,-3,3,nullptr}}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"asym", "h_pt_asym", "A #coloneqq #frac{p_{T1} - p_{T2}}{p_{T1} + p_{T2}};#frac{d#sigma}{dA}", {npt_asym_bins,0,1.,nullptr}}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"pair_pt_ptlead_ratio", "h_pair_pt_ptlead_ratio", title_str("#frac{p_{T}^{pair}}{p_{T}^{lead}}"), {50,0,2.,nullptr}}));

		hist2Dinfo_list.push_back(new Hist2DInfo({"dphi", "etaavg", "h_eta_avg_Dphi", 	";#Delta#phi;#bar{#eta}", 		{nDphi_bins,-pms.PI,pms.PI,nullptr},{neta_bins,-2.4,2.4,nullptr}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"dphi", "deta", "h_Deta_Dphi", 		";#Delta#phi;#Delta#eta", 		{nDphi_bins,-pms.PI,pms.PI,nullptr},{nDeta_bins,-4.8,4.8,nullptr}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"m1.eta", "m2.eta", "h_eta2_eta1", 		";#eta_{sublead};#eta_{lead}",	{neta_bins,-2.4,2.4,nullptr}, 		{neta_bins,-2.4,2.4,nullptr}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"deta", "etaavg", "h_eta_avg_Deta", 	";#Delta#eta;#bar{#eta}",		{nDeta_bins,-4.8,4.8,nullptr},		{neta_bins,-2.4,2.4,nullptr}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"m2.pt", "m1.pt", "h_pt1_pt2", 			";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",		{npT_lead_bins_linear,0,15,nullptr},{npT_lead_bins_linear,0,15,nullptr}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"pair_pt", "m1.pt", "h_ptlead_pair_pt", 			";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",	{npair_pT_bins_linear,0,30,nullptr},{npT_lead_bins_linear,0,30,nullptr}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"pair_pt", "minv", "h_minv_pair_pt", 			";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",		{npair_pT_bins_linear,0,30,nullptr},{nminv_bins_linear,0,30,nullptr}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"pair_pt", "m1.pt", "h_ptlead_pair_pt_zoomin", 	";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",	{npair_pT_bins_linear,0,20,nullptr},{npT_lead_bins_linear,0,15,nullptr}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"pair_pt", "minv", "h_minv_pair_pt_zoomin", 		";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",		{npair_pT_bins_linear,0,20,nullptr},{nminv_bins_linear,5,15,nullptr}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"pair_pt", "m1.pt", "h_ptlead_pair_pt_log", 		";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",	{pms.npairPT_bins,0,0,pms.pairPTBins[1][2]},{pms.npt_bins,0,0,pms.pTBins}}));
		hist2Dinfo_list.push_back(new Hist2DInfo({"pair_pt", "minv", "h_minv_pair_pt_log", 		";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",		{pms.npairPT_bins,0,0,pms.pairPTBins[1][2]},{nminv_bins_log,0,0,minv_bins_log[1]}}));
		// hist2Dinfo_list.push_back(new Hist2DInfo({"pair_pt", "m1.pt", "h_ptlead_pair_pt_log", 		";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",	{pms.npairPT_bins,0,0,pms.pairPTBins[ksign][2]},{pms.npt_bins,0,0,pms.pTBins}}));
		// hist2Dinfo_list.push_back(new Hist2DInfo({"pair_pt", "minv", "h_minv_pair_pt_log", 		";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",		{pms.npairPT_bins,0,0,pms.pairPTBins[ksign][2]},{nminv_bins_log,0,0,minv_bins_log[ksign]}}));
	}


	~RDFBasedMuonPairPlotting(){
		std::cout << "destructor for base class called" << std::endl;
	}
	static bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
	static bool BothMuonPassGapCuts(float m1eta, float m1pt, int m1charge, float m2eta, float m2pt, int m2charge);
	static bool PhotoProductionCut(float x, float y);

	static ROOT::VecOps::RVec<float> HadronMuonDphiCalc(float m1_phi, float m2_phi, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m);
	static ROOT::VecOps::RVec<float> CHadronMuonDphiCalc(bool both_from_c, float m1_phi, float m2_phi, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m);
	static ROOT::VecOps::RVec<float> MuonHQJetPtRatioCalc(float m1_pt, float m2_pt, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m);
	static ROOT::VecOps::RVec<float> MuonCJetPtRatioCalc(bool both_from_c, float m1_pt, float m2_pt, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m);

	virtual void Initialize() = 0;
	virtual void Hist1DMaking(const Hist1DInfo * h1d, bool same_sign, int ifile, ROOT::RDF::RNode df) = 0;
	virtual void Hist2DMaking(const Hist2DInfo * h2d, bool same_sign, int ifile, ROOT::RDF::RNode df) = 0;

	void FillSingleHist1D(const Hist1DInfo * h1d, std::string hist_name_base, ROOT::RDF::RNode df_filtered, bool gapphoto);
	// void FillSingleHist1D(const Hist1DInfo * h1d, std::string hist_name_base, ROOT::RDF::RNode df_filtered, bool gapphoto, bool isRVec = false);
	void FillSingleHist2D(const Hist2DInfo * h2d, std::string hist_name_base, ROOT::RDF::RNode df_filtered, bool gapphoto);

	void Finalize();
	void Run();
	void RunTest();

	// template <typename T, typename Integer>
	// void WeightCalc(T & weight_factor);  // calculting meta-data-dependent weight factors needed for combining multiple files

	// template <typename T, typename Integer>
	// void WeightCalc(std::vector<T> & weight_factors);  // calculting meta-data-dependent weight factors needed for combining multiple files

	// template <typename T, typename Integer>
	// void WeightCalc(std::vector<std::vector<T>> & weight_factors); // calculting meta-data-dependent weight factors needed for combining multiple files

	template<typename T>
	static bool ApplySingleBinSelection(T value, int bin_num, const std::vector<T>& bin_bdrys, bool last_bin_open_ended = false){
		if (last_bin_open_ended){
			if (bin_num != bin_bdrys.size() - 1){
				std::cout << "The bin num being looked at must be (vector size - 1)!" << std::endl;
				throw std::exception();
			}
			return (value > bin_bdrys[bin_num]);
		}else{
			if (bin_num > bin_bdrys.size() - 2 || bin_num < 0){
				std::cout << "The bin num must be between 0 and (vector size - 2)!" << std::endl;
				std::cout << "The passed-in bin num is: " << bin_num << std::endl;
				std::cout << "vector size is: " << bin_bdrys.size() << std::endl;
				throw std::exception();
			}
			return (value >= bin_bdrys[bin_num] && value < bin_bdrys[bin_num+1]);
		}
	}
};



// ----------------------------------------------------------------------------------------------------------------
class RDFBasedMuonPairPlottingPbPb : public RDFBasedMuonPairPlotting{
public:

	enum DataMode{
		real,
		scram
	};

	int data_mode = DataMode::real;

	RDFBasedMuonPairPlottingPbPb(){
		std::cout << "constructor for pbpb called" << std::endl; 
	}
	~RDFBasedMuonPairPlottingPbPb(){
		std::cout << "destructor for pbpb called" << std::endl; 
	}

	void Initialize();
	void Hist1DMaking(const Hist1DInfo * h1d, bool same_sign, int ifile, ROOT::RDF::RNode df);
	void Hist2DMaking(const Hist2DInfo * h2d, bool same_sign, int ifile, ROOT::RDF::RNode df);
};


// ----------------------------------------------------------------------------------------------------------------
class RDFBasedMuonPairPlottingPP : public RDFBasedMuonPairPlotting{
public:

	enum DataMode{
		real,
		scram
	};

	int data_mode = DataMode::real;

	RDFBasedMuonPairPlottingPP(){
		std::cout << "constructor for pp called" << std::endl; 
	}
	~RDFBasedMuonPairPlottingPP(){
		std::cout << "destructor for pp called" << std::endl;
	}
	
	void Initialize();
	void Hist1DMaking(const Hist1DInfo * h1d, bool same_sign, int ifile, ROOT::RDF::RNode df);
	void Hist2DMaking(const Hist2DInfo * h2d, bool same_sign, int ifile, ROOT::RDF::RNode df);
};


// ----------------------------------------------------------------------------------------------------------------
class RDFBasedMuonPairPlottingPowheg : public RDFBasedMuonPairPlotting{
private:
	double powheg_weight_factor = 0.;
public:

	enum PowhegMode{
		bb,
		cc
	};

	int powheg_mode;
	RDFBasedMuonPairPlottingPowheg(){
		std::cout << "constructor for powheg called" << std::endl; 
		auto title_str = [](std::string kin){return ";" + kin + ";#frac{d#sigma}{d" + kin + "}";};
		// auto title_str_w_unit = [](std::string kin, std::string unit){return ";" + kin + " [" + unit + "];#frac{d#sigma}{d" + kin + "}";};

		hist1Dinfo_list.push_back(new Hist1DInfo({"b_hadron_muon_dphi", "h_b_hadron_muon_dphi", title_str("#Delta #phi_{#mu, b hadron}"), {32, -1.,1., nullptr}}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"c_hadron_muon_dphi", "h_c_hadron_muon_dphi", title_str("#Delta #phi_{#mu, c hadron}"), {32, -1.,1., nullptr}}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"muon_b_jet_pt_ratio", "h_muon_b_jet_pt_ratio", title_str("#frac{p_{T}^{muon}}{p_{T}^{b}}"), 40,0,1.,nullptr}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"muon_c_jet_pt_ratio", "h_muon_c_jet_pt_ratio", title_str("#frac{p_{T}^{muon}}{p_{T}^{c}}"), 40,0,1.,nullptr}));
	}
	~RDFBasedMuonPairPlottingPowheg(){
		std::cout << "destructor for powheg called" << std::endl; 
	}

	template <typename Integer>
	void WeightCalc();

	void Initialize();
	void Hist1DMaking(const Hist1DInfo * h1d, bool same_sign, int ifile, ROOT::RDF::RNode df);
	void Hist2DMaking(const Hist2DInfo * h2d, bool same_sign, int ifile, ROOT::RDF::RNode df);
};


// ----------------------------------------------------------------------------------------------------------------
class RDFBasedMuonPairPlottingPythia : public RDFBasedMuonPairPlotting{
private:
	std::vector<std::vector<double>> pythia_weight_factors;
public:
	
	RDFBasedMuonPairPlottingPythia(){
		std::cout << "constructor for pythia called" << std::endl;
		auto title_str = [](std::string kin){return ";" + kin + ";#frac{d#sigma}{d" + kin + "}";};
		// auto title_str_w_unit = [](std::string kin, std::string unit){return ";" + kin + " [" + unit + "];#frac{d#sigma}{d" + kin + "}";};

		hist1Dinfo_list.push_back(new Hist1DInfo({"b_hadron_muon_dphi", "h_b_hadron_muon_dphi", title_str("#Delta #phi_{#mu, b hadron}"), ParamsSet::dphi_axis}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"c_hadron_muon_dphi", "h_c_hadron_muon_dphi", title_str("#Delta #phi_{#mu, c hadron}"), ParamsSet::dphi_axis}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"muon_b_jet_pt_ratio", "h_muon_b_jet_pt_ratio", title_str("#frac{p_{T}^{muon}}{p_{T}^{b}}"), 40,0,1.,nullptr}));
		hist1Dinfo_list.push_back(new Hist1DInfo({"muon_c_jet_pt_ratio", "h_muon_c_jet_pt_ratio", title_str("#frac{p_{T}^{muon}}{p_{T}^{c}}"), 40,0,1.,nullptr}));
	}
	~RDFBasedMuonPairPlottingPythia(){
		std::cout << "destructor for pythia called" << std::endl; 
	}

	template <typename Integer>
	void WeightCalc();

	void Initialize();
	void Hist1DMaking(const Hist1DInfo * h1d, bool same_sign, int ifile, ROOT::RDF::RNode df);
	void Hist2DMaking(const Hist2DInfo * h2d, bool same_sign, int ifile, ROOT::RDF::RNode df);

	template <typename WeightType, typename PTHatType>
	static WeightType PTHatBinnedWeightCorrection(WeightType weight_raw, PTHatType pTHat, int ifile, const std::vector<std::vector<WeightType>> & weight_factors);

};


#endif