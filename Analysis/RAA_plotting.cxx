#include "MuonObjectsParamsAndHelpers/DatasetTriggerMap.h"

class RAAPlotting{
private:

	std::string pp_infile;
	std::string pbpb_infile;
	std::string base_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/";
	std::string out_dir;
	std::string run_year_trigger_suffix;
	std::string legend_pbpb_label;
	std::string legend_pp_label;

	// x: pair pt; y: pair eta; z: centrality
    std::string h3d_name_ss = "h3d_ss_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt";
    std::string h3d_name_op = "h3d_op_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt";
    
    std::string h2d_name = "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts";

	int nctr_bins = 6;
	std::vector<std::vector<int>> ctr_rebins = {{1},{2},{3},{4},{5},{6}};
	std::vector<Color_t> ctr_colors = {kBlack, kBlue, kRed, kGreen+2, kMagenta, kOrange+2};
	std::vector<std::string> ctr_labels = {"Centrality: 0-5%", "Centrality: 5-10%", "_Centrality: 10-20%", "Centrality: 20-30%", "Centrality: 30-50%", "Centrality: 50-80%"};
	
	std::vector<std::vector<int>> pair_pt_rebins = {{1,2,3,4,5},{6,7,8,9},{10,11,12}};
	std::vector<Color_t> pair_pt_colors = {kBlue, kRed, kGreen+2};
	std::vector<std::string> pair_pt_labels = {"p_{T}^{pair} 8-21GeV", "p_{T}^{pair} 21-45GeV", "p_{T}^{pair} 45-80GeV"};

	// the rebins, colors & labels for the variable showing up as different lines
	std::vector<std::vector<int>> line_var_rebins;
	std::vector<Color_t> line_var_colors;
	std::vector<std::string> line_var_labels;

    TFile* infile_pp = nullptr;
    TFile* infile_pbpb = nullptr;

	TH2D* h2d_crossx_pp;
	TH3D* h3d_crossx_pbpb_ss;
	TH3D* h3d_crossx_pbpb_op;

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
	infile_pp->Close();
	delete infile_pp;
	infile_pbpb->Close();
	delete infile_pbpb;
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
	default:
		std::cout << "Run year + trigger mode must be in range 1-5! Either not set or out of range!" << std::endl;
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

    infile_pp = TFile::Open(pp_infile.c_str());
    if (!infile_pp || infile_pp->IsZombie()) {
        std::cerr << "Error opening file: " << pp_infile << std::endl;
        throw std::exception();
    }
    infile_pbpb = TFile::Open(pbpb_infile.c_str());
    if (!infile_pbpb || infile_pbpb->IsZombie()) {
        std::cerr << "Error opening file: " << pbpb_infile << std::endl;
        throw std::exception();
    }

	h2d_crossx_pp = dynamic_cast<TH2D*>(infile_pp->Get(h2d_name.c_str()));
	h3d_crossx_pbpb_ss = dynamic_cast<TH3D*>(infile_pbpb->Get(h3d_name_ss.c_str()));
	h3d_crossx_pbpb_op = dynamic_cast<TH3D*>(infile_pbpb->Get(h3d_name_op.c_str()));


    if (!h3d_crossx_pbpb_ss) {
        std::cerr << "Error getting TH3D " << h3d_name_ss << std::endl;
        throw std::exception();
    }
    if (!h3d_crossx_pbpb_op) {
        std::cerr << "Error getting TH3D " << h3d_name_op << std::endl;
        throw std::exception();
    }

    h3d_crossx_pbpb_op->Add(h3d_crossx_pbpb_ss, -1.);

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
	InputOutputPrepare();
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

	    TLegend* l = new TLegend(0.15,0.4,0.6,0.89);

		l->SetBorderSize(0);
		l->SetFillStyle(0);
		l->SetTextFont(42);
		l->SetMargin(0.2);
		l->SetTextColor(1);

		double ylim = 0.;

	    for (int iline = (subplot-1) * line_var_rebins.size() / nSubplots; iline < subplot * line_var_rebins.size() / nSubplots; iline++){
			TH1D* hcrossx_pbpb_centr_cur_ctr = hcrossx_pbpb_proj_list.at(iline);
	    	TH1D* h_RAA_cur_bin = (TH1D*) hcrossx_pbpb_centr_cur_ctr->Clone("h_RAA_cur_bin");
			if (mode == 1 || mode == 2){
				h_RAA_cur_bin->Divide(hcrossx_pp_proj);
			}else{
				cout << "before + after dividing by pp:" << endl;
				for (int ibin = 0; ibin < h_RAA_cur_bin->GetNbinsX(); ibin++){
					try{
						cout << "bin " << ibin << ": " << h_RAA_cur_bin->GetBinContent(ibin) << ", ";
						h_RAA_cur_bin->SetBinContent(ibin, h_RAA_cur_bin->GetBinContent(ibin) / pp_crossx_intgr_in_pair_pT_bin.at(iline));
						cout << h_RAA_cur_bin->GetBinContent(ibin) << endl;
					}catch (const std::out_of_range& e) {
					    std::cerr << "Out of Range error when setting RAA vs centrality bin content" << std::endl;
					    std::cerr << "Occurs at the " << ibin << "-th rebinned pair-pT bin" << std::endl;
					    std::cerr << "Results in incorrect RAA vs centrality plot!" << std::endl;
					}
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
		l->AddEntry("","m_{#mu#mu} #in (1.08, 2.9) GeV, p_{T}^{pair} > 8 GeV","");
		l->AddEntry("","opposite sign","");
	    

		double ymax = ylim * 1.1;
		if (mode == 1 && use_same_y_axis_range_23_24 && run_year_trigger_mode > 1){ // pt mode; use same range for 23 + 24; running 23/24 plotting
			ymax = ymax_pt_mode_23_24_common.at(run_year_trigger_mode % 2).at(subplot-1);
		}
	    h_RAA_list.at((subplot-1) * line_var_rebins.size() / nSubplots)->GetYaxis()->SetRangeUser(0,ymax);
	    
	    for (int iline = (subplot-1) * line_var_rebins.size() / nSubplots; iline < subplot * line_var_rebins.size() / nSubplots; iline++){
	    	TH1D* h_RAA_cur_bin = h_RAA_list.at(iline);
	    	std::string draw_options = (iline == 0)? "" : "same";
	    	h_RAA_cur_bin->Draw(draw_options.c_str());
		}
		l->Draw("same");
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
	// -------------------- Run2 --------------------
	// RAA_plotting_single_run_year_trigger_mode(1);
	
	// -------------------- 2023 PbPb, 2024 pp mu4mu4noL1 --------------------
	RAA_plotting_single_run_year_trigger_mode(2,true);

	// -------------------- 2023 PbPb, 2024 pp 2mu4 --------------------
	RAA_plotting_single_run_year_trigger_mode(3,true);

	// -------------------- 2024 PbPb, 2024 pp mu4mu4noL1 --------------------
	RAA_plotting_single_run_year_trigger_mode(4,true);

	// -------------------- 2024 PbPb, 2024 pp 2mu4 --------------------
	RAA_plotting_single_run_year_trigger_mode(5,true);



}