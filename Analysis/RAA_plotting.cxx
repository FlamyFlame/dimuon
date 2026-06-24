#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include "AtlasStyle.C"   // ATLAS frame style (SetAtlasStyle); header AtlasStyle.h
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "MuonObjectsParamsAndHelpers/DatasetTriggerMap.h"
#include "Utilities/PbPbSampledLumi.h"

class RAAPlotting{
private:

	std::string pp_infile;
	std::string pbpb_infile;
	std::vector<std::string> pbpb_infiles; // if non-empty: combine these PbPb-year RDF files (luminosity-weighted, per HF R_AA note Eq.3)
	std::vector<double>      pbpb_lumis;   // per-year HLT_mu4 sampled lumi [nb^-1], parallel to pbpb_infiles (weights for the combine)
	// Cluster data area. Inputs are the RDFBasedHistFilling crossx outputs
	// (histograms_real_pairs_*), now reco+trig efficiency corrected.
	std::string base_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
	std::string out_dir;
	std::string run_year_trigger_suffix;
	std::string legend_pbpb_label;
	std::string legend_pp_label;

	// x: pair pt; y: pair eta; z: centrality
    std::string h3d_name_ss = "h3d_ss_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt";
    std::string h3d_name_op = "h3d_op_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt";

    std::string h2d_name    = "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts";    // pp OS
    std::string h2d_name_ss = "h2d_ss_crossx_pair_pt_pair_eta_binned_w_signal_cuts"; // pp SS (combinatorial)

	int nctr_bins = 6;
	std::vector<std::vector<int>> ctr_rebins = {{1},{2},{3},{4},{5},{6}};
	std::vector<Color_t> ctr_colors = {kBlack, kBlue, kRed, kGreen+2, kMagenta, kOrange+2};
	std::vector<std::string> ctr_labels = {"Centrality: 0-5%", "Centrality: 5-10%", "Centrality: 10-20%", "Centrality: 20-30%", "Centrality: 30-50%", "Centrality: 50-80%"};
	
	// RDF crossx pT axis has 15 bins, edges {8,10,11,14,16,20,24,28,34,41,49,58,70,84,100,120}.
	// Group into 3 physics ranges covering ALL 15 bins (legacy used 12 bins / pT_bins_80).
	std::vector<std::vector<int>> pair_pt_rebins = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
	std::vector<Color_t> pair_pt_colors = {kBlue, kRed, kGreen+2};
	std::vector<std::string> pair_pt_labels = {"p_{T}^{pair} 8-20GeV", "p_{T}^{pair} 20-49GeV", "p_{T}^{pair} 49-120GeV"};

	// the rebins, colors & labels for the variable showing up as different lines
	std::vector<std::vector<int>> line_var_rebins;
	std::vector<Color_t> line_var_colors;
	std::vector<std::string> line_var_labels;

    TFile* infile_pp = nullptr;
    TFile* infile_pbpb = nullptr;

	TH2D* h2d_crossx_pp = nullptr;
	TH3D* h3d_crossx_pbpb_ss = nullptr;
	TH3D* h3d_crossx_pbpb_op = nullptr;

	TH1D* hcrossx_pp_proj;
	std::vector<double> pp_crossx_intgr_in_pair_pT_bin;
	std::vector<TH1D*> hcrossx_pbpb_proj_list;
	std::vector<TH1D*> h_RAA_list;

	std::string obs_suffix;
	std::string ytitle;

	int mode;
	int run_year_trigger_mode;
	int nBinsThresh = 5;
	int nSubplots = 1;

	void InputOutputPrepare();
	void ModePrepare();
	void HistRetrieve();
	void HistProject();
public:
	bool use_same_y_axis_range_23_24 = false;
	vector<vector<double>> ymax_pt_mode_23_24_common = {{0.92, 1.2}, {2.4, 3}}; // top dim: trigger; 2nd dim: ctr-based subplot

	RAAPlotting(int mode_in, int run_year_trigger_mode_in): mode(mode_in), run_year_trigger_mode(run_year_trigger_mode_in){}
	~RAAPlotting();
	void RunPlotting();
	void PrintInstructions();
};

RAAPlotting::~RAAPlotting(){
	// infile_pbpb stays null in the combined (multi-file) path — year files are
	// opened and closed inside HistRetrieve. Guard both.
	if (infile_pp)   { infile_pp->Close();   delete infile_pp; }
	if (infile_pbpb) { infile_pbpb->Close(); delete infile_pbpb; }
}

void RAAPlotting::InputOutputPrepare(){
	switch(run_year_trigger_mode){
	case 1: // Run2 PbPb + pp data
		pp_infile = base_dir + "pp_run2/pp_run2_single_b_ana_hists.root";
		pbpb_infile = base_dir + "pbpb_run2/pbpb_run2_single_b_ana_hists.root";
		run_year_trigger_suffix = "_pbpbrun2_pp17_2mu4";
		out_dir = base_dir + "pbpb_run2/plots/";
		legend_pbpb_label = "PbPb Run2 data, mu4_mu4noL1";
		legend_pp_label   = "pp 2017 data, 2mu4";
		break;
	case 2: // 2023 PbPb, 2024 pp mu4mu4noL1
		pp_infile = base_dir + "pp_2024/pp_2024_single_b_ana_hists_mu4_mu4noL1.root";
		pbpb_infile = base_dir + "pbpb_2023/pbpb_2023_single_b_ana_hists.root";
		run_year_trigger_suffix = "_pbpb23_pp24_mu4mu4noL1";
		out_dir = base_dir + "pbpb_2023/plots/";
		legend_pbpb_label = "PbPb 2023 data, " + DatasetTriggerMap::GetTriggerLabel(23, "PbPb");
		legend_pp_label   = "pp 2024 data, "   + DatasetTriggerMap::GetTriggerLabel(24, "pp");
		break;
	case 3: // 2023 PbPb, 2024 pp 2mu4
		pp_infile = base_dir + "pp_2024/pp_2024_single_b_ana_hists_2mu4.root";
		pbpb_infile = base_dir + "pbpb_2023/pbpb_2023_single_b_ana_hists.root";
		run_year_trigger_suffix = "_pbpb23_pp24_2mu4";
		out_dir = base_dir + "pbpb_2023/plots/";
		legend_pbpb_label = "PbPb 2023 data, " + DatasetTriggerMap::GetTriggerLabel(23, "PbPb");
		legend_pp_label   = "pp 2024 data, "   + DatasetTriggerMap::GetTriggerLabel(24, "pp_2mu4");
		break;
	case 4: // 2024 PbPb, 2024 pp mu4mu4noL1
		pp_infile = base_dir + "pp_2024/pp_2024_single_b_ana_hists_mu4_mu4noL1.root";
		pbpb_infile = base_dir + "pbpb_2024/pbpb_2024_single_b_ana_hists.root";
		run_year_trigger_suffix = "_pbpb24_pp24_mu4mu4noL1";
		out_dir = base_dir + "pbpb_2024/plots/";
		legend_pbpb_label = "PbPb 2024 data, " + DatasetTriggerMap::GetTriggerLabel(24, "PbPb");
		legend_pp_label   = "pp 2024 data, "   + DatasetTriggerMap::GetTriggerLabel(24, "pp");
		break;
	case 5: // 2024 PbPb, 2024 pp 2mu4
		pp_infile = base_dir + "pp_2024/pp_2024_single_b_ana_hists_2mu4.root";
		pbpb_infile = base_dir + "pbpb_2024/pbpb_2024_single_b_ana_hists.root";
		run_year_trigger_suffix = "_pbpb24_pp24_2mu4";
		out_dir = base_dir + "pbpb_2024/plots/";
		legend_pbpb_label = "PbPb 2024 data, " + DatasetTriggerMap::GetTriggerLabel(24, "PbPb");
		legend_pp_label   = "pp 2024 data, "   + DatasetTriggerMap::GetTriggerLabel(24, "pp_2mu4");
		break;
	case 6: // Run 3: PbPb 2023+2024+2025 combined (single_mu4) + pp 2024 (2mu4), from RDF crossx outputs
		pp_infile   = base_dir + "pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root";
		pbpb_infiles = {
			base_dir + "pbpb_2023/histograms_real_pairs_pbpb_2023_single_mu4_no_trg_plots_nominal.root",
			base_dir + "pbpb_2024/histograms_real_pairs_pbpb_2024_single_mu4_no_trg_plots_nominal.root",
			base_dir + "pbpb_2025/histograms_real_pairs_pbpb_2025_single_mu4_no_trg_plots_nominal.root",
		};
		pbpb_lumis = { PbPbMu4SampledLumiNb(23), PbPbMu4SampledLumiNb(24), PbPbMu4SampledLumiNb(25) };
		run_year_trigger_suffix = "_pbpb23_24_25_combined_pp24_2mu4";
		out_dir = base_dir + "plots/single_b_analysis/RAA/";
		legend_pbpb_label = "Pb+Pb 2023+2024+2025, HLT_mu4";
		legend_pp_label   = "pp 2024, HLT_2mu4";
		break;

	default:
		std::cout << "Run year + trigger mode must be in range 1-6! Either not set or out of range!" << std::endl;
		throw std::exception();
	}
}


void RAAPlotting::ModePrepare(){
	if (mode == 1 || mode == 2){ // var1 (rebinned) = centrality
		line_var_rebins = ctr_rebins;
		line_var_colors = ctr_colors;
		line_var_labels = ctr_labels;
	}else{ // mode == 3: var1 (rebinned) = pair pT
		line_var_rebins = pair_pt_rebins;
		line_var_colors = pair_pt_colors;
		line_var_labels = pair_pt_labels;
	}

	switch (mode){
	case 1:
		obs_suffix = "pair_pt";
		ytitle = "RAA";
		break;
	case 2:
		obs_suffix = "pair_eta";
		ytitle = "RAA";
		break;
	case 3:
		obs_suffix = "ctr";
		ytitle = "RAA";
		break;
	default:
		throw std::runtime_error("INVALID MODE: mode must be 1, 2, or 3!");
	}
}

void RAAPlotting::HistRetrieve(){

    // ---- pp cross-section: OS - SS (combinatorial subtraction, analysis_overview §4a) ----
    infile_pp = TFile::Open(pp_infile.c_str());
    if (!infile_pp || infile_pp->IsZombie()) {
        std::cerr << "Error opening pp file: " << pp_infile << std::endl;
        throw std::exception();
    }
    TH2D* h2d_pp_os = dynamic_cast<TH2D*>(infile_pp->Get(h2d_name.c_str()));
    if (!h2d_pp_os) { std::cerr << "Error getting pp TH2D " << h2d_name << std::endl; throw std::exception(); }
    h2d_crossx_pp = dynamic_cast<TH2D*>(h2d_pp_os->Clone("h2d_crossx_pp_osmss"));
    h2d_crossx_pp->SetDirectory(nullptr);
    TH2D* h2d_pp_ss = dynamic_cast<TH2D*>(infile_pp->Get(h2d_name_ss.c_str()));
    if (h2d_pp_ss) h2d_crossx_pp->Add(h2d_pp_ss, -1.);
    else std::cerr << "WARNING: pp SS hist " << h2d_name_ss
                   << " missing -> OS-only pp (no combinatorial subtraction)" << std::endl;

    // ---- PbPb yield: LUMINOSITY-WEIGHTED combine over year files, then OS - SS ----
    // combined = Sum_y(L_y * h_y) / Sum_y(L_y)  (HF R_AA note HION-2019-58 §4.1 Eq.3).
    // Each h_y is already crossx_factor-weighted (~1/L_y), so L_y*h_y recovers the
    // year-independent (counts * K) and the result is the correct ΣN/(...·ΣL).
    std::vector<std::string> pbpb_files = pbpb_infiles.empty()
        ? std::vector<std::string>{pbpb_infile} : pbpb_infiles;
    bool any_pbpb_ss = false;
    double sumL = 0.;
    for (size_t iy = 0; iy < pbpb_files.size(); ++iy){
        const std::string& fpath = pbpb_files[iy];
        const double L = (iy < pbpb_lumis.size()) ? pbpb_lumis[iy] : 1.0; // legacy single-file -> weight 1
        TFile* f = TFile::Open(fpath.c_str());
        if (!f || f->IsZombie()) { std::cerr << "Error opening PbPb file: " << fpath << std::endl; throw std::exception(); }
        TH3D* op = dynamic_cast<TH3D*>(f->Get(h3d_name_op.c_str()));
        TH3D* ss = dynamic_cast<TH3D*>(f->Get(h3d_name_ss.c_str()));
        if (!op) { std::cerr << "Error getting TH3D " << h3d_name_op << " in " << fpath << std::endl; throw std::exception(); }
        if (!h3d_crossx_pbpb_op){
            h3d_crossx_pbpb_op = dynamic_cast<TH3D*>(op->Clone("h3d_op_comb"));
            h3d_crossx_pbpb_op->SetDirectory(nullptr);
            h3d_crossx_pbpb_op->Scale(L);
        } else h3d_crossx_pbpb_op->Add(op, L);
        if (ss){
            any_pbpb_ss = true;
            if (!h3d_crossx_pbpb_ss){
                h3d_crossx_pbpb_ss = dynamic_cast<TH3D*>(ss->Clone("h3d_ss_comb"));
                h3d_crossx_pbpb_ss->SetDirectory(nullptr);
                h3d_crossx_pbpb_ss->Scale(L);
            } else h3d_crossx_pbpb_ss->Add(ss, L);
        }
        sumL += L;
        f->Close(); delete f;
    }
    if (sumL > 0.){
        h3d_crossx_pbpb_op->Scale(1.0 / sumL);
        if (h3d_crossx_pbpb_ss) h3d_crossx_pbpb_ss->Scale(1.0 / sumL);
    }
    if (any_pbpb_ss && h3d_crossx_pbpb_ss) h3d_crossx_pbpb_op->Add(h3d_crossx_pbpb_ss, -1.);
    else std::cerr << "WARNING: PbPb SS hist " << h3d_name_ss
                   << " missing -> OS-only PbPb (no combinatorial subtraction)" << std::endl;
}


void RAAPlotting::HistProject(){
    if (mode == 1 || mode == 3){
		hcrossx_pp_proj = h2d_crossx_pp->ProjectionX("h_pp_crossx_pair_pt");
    }else{ // mode == 2: RAA vs pair eta
		hcrossx_pp_proj = h2d_crossx_pp->ProjectionY("h_pp_crossx_pair_eta");
    }
	
	hcrossx_pp_proj->Scale(1.,"width");		    	

	pp_crossx_intgr_in_pair_pT_bin.clear();

    for (int ibin = 0; ibin < line_var_rebins.size(); ibin++){
    	auto line_bins_cur_rebin = line_var_rebins.at(ibin);
    	int bin_num_first = line_bins_cur_rebin.at(0);
    	int bin_num_last  = line_bins_cur_rebin.at(line_bins_cur_rebin.size()-1);
    	TH1D* h_proj;
    	switch (mode){
    	case 1:
    		// project x-axis (pair pT) over rebinned z bin (ctr); sum over y axis (pair eta)
    		h_proj = h3d_crossx_pbpb_op->ProjectionX(Form("h_proj_bin%d", ibin),0,-1,bin_num_first,bin_num_last);
	    	h_proj->Scale(1.,"width");
    		hcrossx_pbpb_proj_list.push_back(h_proj);
    		break;
    	case 2:
    		// project y-axis (pair eta) over rebinned z bin (ctr); sum over x axis (pair pT)
    		h_proj = h3d_crossx_pbpb_op->ProjectionY(Form("h_proj_bin%d", ibin),0,-1,bin_num_first,bin_num_last);
	    	h_proj->Scale(1.,"width"); 
    		hcrossx_pbpb_proj_list.push_back(h_proj);
    		break;
    	case 3:
    		// project z-axis (ctr) over rebinned x bin (pair pT); sum over x axis (pair eta)
    		h_proj = h3d_crossx_pbpb_op->ProjectionZ(Form("h_proj_bin%d", ibin),bin_num_first,bin_num_last,0,-1);
    		// divide by the pair-pT rebinned bin width
	    	// h_proj->Scale(1./(hcrossx_pp_proj->GetBinLowEdge(bin_num_last + 1) - hcrossx_pp_proj->GetBinLowEdge(bin_num_first)));
    		hcrossx_pbpb_proj_list.push_back(h_proj);

    		// integrate over pp crossx in pp_crossx_intgr_in_pair_pT_bin
			pp_crossx_intgr_in_pair_pT_bin.push_back(hcrossx_pp_proj->Integral(bin_num_first,bin_num_last,"width"));
    		break;
    	}
    }
}


void RAAPlotting::PrintInstructions(){
	std::cout << "Class RAAPlotting: runs RAA plotting with one of the 3 modes" << std::endl; 
	std::cout << "Mode 1: plots RAA vs pair pT in given centrality bins (integrate over pair eta)" << std::endl;
	std::cout << "Mode 2: plots RAA vs pair eta in given centrality bins (integrate over pair pT)" << std::endl;
	std::cout << "Mode 3: plots RAA vs centrality in given pair pT bins (integrate over pair eta)" << std::endl;
	std::cout << "For a given run year: if run year < 2020, use run2 pp + PbPb data" << std::endl;
	std::cout << "else, check that (run year) % 2000 == 23 or 24, and uses corresponding data" << std::endl;
}


void RAAPlotting::RunPlotting(){
	SetAtlasStyle();   // apply ATLAS frame style before any object is created
	InputOutputPrepare();
	gSystem->mkdir(out_dir.c_str(), true); // ensure the (cluster) output dir exists
	ModePrepare();
	HistRetrieve();
	HistProject();

	Int_t canvas_width = (line_var_rebins.size() > nBinsThresh)? 1200 : 700;
    TCanvas* c_raa = new TCanvas(Form("c_raa_%s", obs_suffix.c_str()), Form("c_raa_%s", obs_suffix.c_str()), canvas_width, 500);

    if (line_var_rebins.size() > nBinsThresh){
    	nSubplots = 2;
	    c_raa->Divide(2,1);
	}

	cout << "pp_crossx_intgr_in_pair_pT_bin content: " << std::endl;
	for (auto pp_intgr : pp_crossx_intgr_in_pair_pT_bin){
	 	cout << pp_intgr << ", ";
	}
	cout << endl;

	for (int subplot = 1; subplot <= nSubplots; subplot++){
	    c_raa->cd(subplot);

	    // Upper-centre, compact: with the y-headroom below, the legend clears the
	    // data bulk everywhere; placing it in x∈[0.30,0.70] also avoids the tall
	    // low-stat error bars at the extreme eta acceptance edges (far left/right).
	    TLegend* l = new TLegend(0.30,0.64,0.70,0.90);

		l->SetBorderSize(0);
		l->SetFillStyle(0);
		l->SetTextFont(42);
		l->SetTextSize(0.030);
		l->SetMargin(0.2);
		l->SetTextColor(1);

		double ylim = 0.;

	    for (int iline = (subplot-1) * line_var_rebins.size() / nSubplots; iline < subplot * line_var_rebins.size() / nSubplots; iline++){
			TH1D* hcrossx_pbpb_centr_cur_ctr = hcrossx_pbpb_proj_list.at(iline);
	    	TH1D* h_RAA_cur_bin = (TH1D*) hcrossx_pbpb_centr_cur_ctr->Clone("h_RAA_cur_bin");
			if (mode == 1 || mode == 2){
				// Index-wise ratio R_AA = PbPb / pp. Bin counts match by construction
				// (same pT/eta axis); done explicitly (not TH1::Divide) because the
				// PbPb 3D eta axis is an explicit edge array while the pp 2D eta axis
				// is uniform — numerically identical edges but different internal
				// representation, which makes TH1::Divide warn. Independent ratio error.
				for (int ib = 1; ib <= h_RAA_cur_bin->GetNbinsX(); ++ib){
					double a = hcrossx_pbpb_centr_cur_ctr->GetBinContent(ib), ea = hcrossx_pbpb_centr_cur_ctr->GetBinError(ib);
					double b = hcrossx_pp_proj->GetBinContent(ib),           eb = hcrossx_pp_proj->GetBinError(ib);
					double r = (b != 0.) ? a / b : 0.;
					double er = (a != 0. && b != 0.) ? r * std::sqrt((ea/a)*(ea/a) + (eb/b)*(eb/b)) : 0.;
					h_RAA_cur_bin->SetBinContent(ib, r);
					h_RAA_cur_bin->SetBinError(ib, er);
				}
			}else{
				// mode 3 (R_AA vs centrality): PbPb yield (ProjectionZ over a pair-pT
				// group) divided by the pp cross-section integrated over the SAME
				// pair-pT group (a per-line scalar). Iterate REAL bins 1..N (bin 0 is
				// underflow); the denominator is constant, so the error scales by 1/pp_int.
				const double pp_int = (iline < (int)pp_crossx_intgr_in_pair_pT_bin.size())
					? pp_crossx_intgr_in_pair_pT_bin.at(iline) : 0.;
				for (int ib = 1; ib <= h_RAA_cur_bin->GetNbinsX(); ++ib){
					double a  = h_RAA_cur_bin->GetBinContent(ib);
					double ea = h_RAA_cur_bin->GetBinError(ib);
					h_RAA_cur_bin->SetBinContent(ib, (pp_int != 0.) ? a / pp_int : 0.);
					h_RAA_cur_bin->SetBinError(ib,   (pp_int != 0.) ? ea / pp_int : 0.);
				}
			}
			h_RAA_cur_bin->SetLineColor(line_var_colors.at(iline));
			h_RAA_cur_bin->SetLineWidth(2);
			h_RAA_cur_bin->SetStats(0);

			// assume the last bin (largest pair pT has the largest statistical errorbar)
			// double max_plus_error = h_RAA_cur_bin->GetMaximum() + h_RAA_cur_bin->GetBinError(h_RAA_cur_bin->GetNbinsX());
			double max_plus_error = h_RAA_cur_bin->GetMaximum();
			if (max_plus_error> ylim) ylim = max_plus_error;
			h_RAA_cur_bin->SetTitle("");
			h_RAA_list.push_back(h_RAA_cur_bin);
			l->AddEntry(h_RAA_cur_bin, line_var_labels.at(iline).c_str(), "lp");
		}
		l->AddEntry("", legend_pbpb_label.c_str(), "");
		l->AddEntry("", legend_pp_label.c_str(), "");
		l->AddEntry("","m_{#mu#mu}#in(1.08,2.9), p_{T}^{pair}>8 GeV","");
		l->AddEntry("","OS #minus SS subtracted","");
		l->AddEntry("","reco-eff + T_{AA}: PLACEHOLDERS","");
	    

		double ymax = ylim * 1.7; // headroom so the upper-right legend/info block clears the data (esp. the ~flat eta distribution)
		if (mode == 1 && use_same_y_axis_range_23_24 && run_year_trigger_mode > 1){ // pt mode; use same range for 23 + 24; running 23/24 plotting
			ymax = ymax_pt_mode_23_24_common.at(run_year_trigger_mode % 2).at(subplot-1);
		}
	    TH1D* h_first = h_RAA_list.at((subplot-1) * line_var_rebins.size() / nSubplots);
	    h_first->GetYaxis()->SetRangeUser(0,ymax);
	    h_first->GetYaxis()->SetTitle("R_{AA}");
	    h_first->GetXaxis()->SetTitle( (mode==1) ? "p_{T}^{pair} [GeV]" : (mode==2) ? "#eta^{pair}" : "Centrality [%]" );
	    // mode 2: restrict the displayed eta to the well-populated region; the
	    // |eta| ~ 2.3-2.4 acceptance-edge bins are very low-stat (giant error bars)
	    // and not meaningful for R_AA (and would otherwise spike into the legend).
	    if (mode == 2) h_first->GetXaxis()->SetRangeUser(-2.3, 2.3);
	    
	    for (int iline = (subplot-1) * line_var_rebins.size() / nSubplots; iline < subplot * line_var_rebins.size() / nSubplots; iline++){
	    	TH1D* h_RAA_cur_bin = h_RAA_list.at(iline);
	    	std::string draw_options = (iline == 0)? "" : "same";
	    	h_RAA_cur_bin->Draw(draw_options.c_str());
		}
		l->Draw("same");

		// ATLAS-convention label for talks/slides: "ATLAS" (font 72, bold-italic)
		// + status word (font 42), stacked in the free upper-left strip (legend
		// starts at x=0.30). "Internal" is correct for unapproved results; switch
		// to "Preliminary" only if the result is approved for public showing.
		// Physics caveats remain in the legend (reco-eff + T_AA placeholders, OS-SS).
		TLatex* atlas = new TLatex();
		atlas->SetNDC();
		atlas->SetTextFont(72);           // ATLAS bold-italic
		atlas->SetTextColor(kBlack);
		atlas->SetTextSize(0.052);
		atlas->DrawLatex(0.185, 0.875, "ATLAS");
		TLatex* lab = new TLatex();
		lab->SetNDC();
		lab->SetTextFont(42);
		lab->SetTextColor(kBlack);
		lab->SetTextSize(0.044);
		lab->DrawLatex(0.185, 0.815, "Internal");
	}

	std::string use_same_y_axis_range_23_24_suffix = (mode == 1 && use_same_y_axis_range_23_24 && run_year_trigger_mode > 1)? "_same_23_24_ymax" : "";
	std::string outfile_path = out_dir + "raa_" + obs_suffix + run_year_trigger_suffix + use_same_y_axis_range_23_24_suffix + ".png";
	c_raa->SaveAs(outfile_path.c_str());
	c_raa->Close();
	delete c_raa;
}


void RAA_plotting_single_run_year_trigger_mode(int run_year_trigger_mode, bool use_same_y_axis_range_23_24_arg = false){
	// std::vector raa_plotting_list;
	for (int mode = 1; mode <= 3; mode++){
		auto raa_plotting_local = RAAPlotting(mode, run_year_trigger_mode);
		raa_plotting_local.use_same_y_axis_range_23_24 = use_same_y_axis_range_23_24_arg;
		raa_plotting_local.RunPlotting();
	} 
}

void RAA_plotting(){
	// -------------------- Run 3: PbPb 23+24+25 combined + pp24 2mu4, from RDF crossx --------------------
	// (reco-eff PLACEHOLDER + OS-SS combinatorial subtraction; T_AA 2023 placeholder)
	RAA_plotting_single_run_year_trigger_mode(6);

	// -------------------- legacy cases (SingleBAnalysis-era inputs; not maintained) --------------------
	// RAA_plotting_single_run_year_trigger_mode(2,true); // 2023 PbPb, 2024 pp mu4mu4noL1
	// RAA_plotting_single_run_year_trigger_mode(3,true); // 2023 PbPb, 2024 pp 2mu4
	// RAA_plotting_single_run_year_trigger_mode(4,true); // 2024 PbPb, 2024 pp mu4mu4noL1
	// RAA_plotting_single_run_year_trigger_mode(5,true); // 2024 PbPb, 2024 pp 2mu4
}