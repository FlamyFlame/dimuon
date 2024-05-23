#include "RDFBasedMuonPairPlotting.h"



template<typename T>
bool RDFBasedMuonPairPlotting::ApplySingleBinSelection(T value, int bin_num, const std::vector<T>& bin_bdrys, bool last_bin_open_ended = false);

bool RDFBasedMuonPairPlotting::PassSingleMuonGapCut(float meta, float mpt, int mcharge){
    // parameters: eta and pT of the same muon
    if (fabs(meta) < ParamsSet::eta_gap_cut1) return false;
    if (mpt < 6){
        for (array<float,2> charge_eta_gap_cut : ParamsSet::charge_eta_gap_cuts){
            if (mcharge * meta > charge_eta_gap_cut[0] && mcharge * meta < charge_eta_gap_cut[1]) return false;
        }
    }
    // if (mpt < 6 && fabs(meta) > pms.eta_gap_cut2[0] && fabs(meta) < pms.eta_gap_cut2[1]) return false;
    return true;
}

bool RDFBasedMuonPairPlotting::BothMuonPassGapCuts(float m1eta, float m1pt, int m1charge, float m2eta, float m2pt, int m2charge){
	return PassSingleMuonGapCut(m1eta, m1pt, m1charge) && PassSingleMuonGapCut(m2eta, m2pt, m2charge);
}

bool RDFBasedMuonPairPlotting::PhotoProductionCut(float asym, float acop){
	// pass the cut if NOT from photo production
	return !(asym < 0.05 && acop < 0.01);
}

ROOT::VecOps::RVec<float> RDFBasedMuonPairPlotting::HadronMuonDphiCalc(float m1_phi, float m2_phi, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m){
	// function to be applied to either both-b case or both-c case
	ROOT::VecOps::RVec<float> dphi_vec = {};
	if (hadron1_pt_eta_phi_m.size() == 0 || hadron2_pt_eta_phi_m.size() == 0){
		return dphi_vec; // return empty vector if either of the two muons aren't originally from b
	}
	float dphi1 = m1_phi - hadron1_pt_eta_phi_m[2];
	dphi1 = atan2(sin(dphi1),cos(dphi1));
	float dphi2 = m2_phi - hadron2_pt_eta_phi_m[2];
	dphi2 = atan2(sin(dphi2),cos(dphi2));
	
	dphi_vec.push_back(dphi1);
	dphi_vec.push_back(dphi2);
	return dphi_vec;
}

ROOT::VecOps::RVec<float> RDFBasedMuonPairPlotting::CHadronMuonDphiCalc(bool both_from_c, float m1_phi, float m2_phi, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m){
	ROOT::VecOps::RVec<float> c_dphi_vec = {};
	if (!both_from_c){
		return c_dphi_vec; // return empty if not both from c
	}
	return HadronMuonDphiCalc(m1_phi, m2_phi, hadron1_pt_eta_phi_m, hadron2_pt_eta_phi_m);
}

ROOT::VecOps::RVec<float> RDFBasedMuonPairPlotting::MuonHQJetPtRatioCalc(float m1_pt, float m2_pt, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m){
	// function to be applied to either both-b case or both-c case
	ROOT::VecOps::RVec<float> pt_ratio_vec = {};
	if (hadron1_pt_eta_phi_m.size() == 0 || hadron2_pt_eta_phi_m.size() == 0){
		return pt_ratio_vec; // return empty vector if either of the two muons aren't originally from b
	}
	float pt_ratio1 = m1_pt / hadron1_pt_eta_phi_m[0];
	float pt_ratio2 = m2_pt / hadron2_pt_eta_phi_m[0];
	
	pt_ratio_vec.push_back(pt_ratio1);
	pt_ratio_vec.push_back(pt_ratio2);
	return pt_ratio_vec;
}

ROOT::VecOps::RVec<float> RDFBasedMuonPairPlotting::MuonCJetPtRatioCalc(bool both_from_c, float m1_pt, float m2_pt, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m){
	ROOT::VecOps::RVec<float> c_pt_ratio_vec = {};
	if (!both_from_c){
		return c_pt_ratio_vec;
	}
	return MuonHQJetPtRatioCalc(m1_pt, m2_pt, hadron1_pt_eta_phi_m, hadron2_pt_eta_phi_m);
}

void RDFBasedMuonPairPlotting::FillSingleHist1D(const Hist1DInfo * h1d, std::string hist_name_base, ROOT::RDF::RNode df_filtered, bool gapphoto){
// void RDFBasedMuonPairPlotting::FillSingleHist1D(const Hist1DInfo * h1d, std::string hist_name_base, ROOT::RDF::RNode df_filtered, bool gapphoto, bool isRVec){
	std::string hist_name;
	if (gapphoto){
		if (h1d->xaxis.bin_edges == nullptr){
			hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[0];
			hist1d_rresultptrs_cur_file->push_back(df_filtered.Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
			// hist1d_rresultptrs_cur_file->push_back(df_filtered.Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
			hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[1];
			hist1d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
			// hist1d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
			hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[0];
			hist1d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
			// hist1d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
			hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[1];
			hist1d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
			// hist1d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
		}else{
			std::cout << "AT THIS STAGE - SHOULD BE NULLPTR!!! WHAT'S WRONG?" << std::endl;
			std::cout << h1d->xaxis.bin_edges[0]<< ", " << h1d->xaxis.bin_edges[1] << std::endl;
			std::cout << "*********************************************" << std::endl;
			hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[0];
			hist1d_rresultptrs_cur_file->push_back(df_filtered.Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.bin_edges}, h1d->var_expr, "weight"));
			hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[1];
			hist1d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.bin_edges}, h1d->var_expr, "weight"));
			hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[0];
			hist1d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.bin_edges}, h1d->var_expr, "weight"));
			hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[1];
			hist1d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.bin_edges}, h1d->var_expr, "weight"));
		}
	}else{
		hist_name = hist_name_base;

		if (h1d->xaxis.bin_edges == nullptr){
			// hist1d_rresultptrs_cur_file->push_back(df_filtered.Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
			hist1d_rresultptrs_cur_file->push_back(df_filtered.Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.min, h1d->xaxis.max}, h1d->var_expr, "weight"));
		}else{
			std::cout << "AT THIS STAGE - SHOULD BE NULLPTR!!! WHAT'S WRONG?" << std::endl;
			std::cout << h1d->xaxis.bin_edges[0]<< ", " << h1d->xaxis.bin_edges[1] << std::endl;
			std::cout << "*********************************************" << std::endl;

			hist1d_rresultptrs_cur_file->push_back(df_filtered.Histo1D({hist_name.c_str(), h1d->hist1d_titles.c_str(), h1d->xaxis.nbins, h1d->xaxis.bin_edges}, h1d->var_expr, "weight"));
		}
	}
}


void RDFBasedMuonPairPlotting::FillSingleHist2D(const Hist2DInfo * h2d, std::string hist_name_base, ROOT::RDF::RNode df_filtered, bool gapphoto){
	std::string hist_name;
	if (gapphoto){
		if (h2d->xaxis.bin_edges == nullptr){
			if (h2d->yaxis.bin_edges == nullptr){
				hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[0];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[1];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[0];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[1];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
			}else{
				hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[0];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[1];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[0];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[1];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
			}
		}else{
			if (h2d->yaxis.bin_edges == nullptr){
				hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[0];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[1];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[0];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[1];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
			}else{
				hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[0];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[0] + pms.gapcut_labels[1];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[0];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
				hist_name = hist_name_base + pms.photocut_labels[1] + pms.gapcut_labels[1];
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Filter(RDFBasedMuonPairPlotting::PhotoProductionCut, {"asym", "acop"}).Filter(RDFBasedMuonPairPlotting::BothMuonPassGapCuts, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"}).Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
			}
		}
	}else{
		hist_name = hist_name_base;

		if (h2d->xaxis.bin_edges == nullptr){
			if (h2d->yaxis.bin_edges == nullptr){
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
			}else{
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.min, h2d->xaxis.max, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
			}
		}else{
			if (h2d->yaxis.bin_edges == nullptr){
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.min, h2d->yaxis.max}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
			}else{
				hist2d_rresultptrs_cur_file->push_back(df_filtered.Histo2D<float, float, double>({hist_name.c_str(), h2d->hist2d_titles.c_str(), h2d->xaxis.nbins, h2d->xaxis.bin_edges, h2d->yaxis.nbins, h2d->yaxis.bin_edges}, h2d->xvar_expr, h2d->yvar_expr, "weight"));
			}
		}
	}
}


void RDFBasedMuonPairPlotting::Finalize(){
	for (auto ptr : hist1Dinfo_list){
		delete ptr;
	}
	for (auto ptr : hist2Dinfo_list){
		delete ptr;
	}
	delete hist1Ds;
	delete hist2Ds;
	delete hist1d_rresultptrs_cur_file;
	delete hist2d_rresultptrs_cur_file;
}

void RDFBasedMuonPairPlotting::RunTest(){

	std::cout << "start the run" << std::endl;
	auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);

	Initialize();
}

void RDFBasedMuonPairPlotting::Run(){

	std::cout << "start the run" << std::endl;
	auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);

	Initialize();
	// hist1Ds = new std::vector<std::vector<ROOT::RDF::RResultPtr< ::TH1D>>>;
	hist1Ds = new std::vector<std::vector<TH1D>>;
	hist2Ds = new std::vector<std::vector<TH2D>>;

	hist1d_rresultptrs_cur_file = new std::vector<ROOT::RDF::RResultPtr< ::TH1D>>;
	hist2d_rresultptrs_cur_file = new std::vector<ROOT::RDF::RResultPtr< ::TH2D>>;


	ROOT::EnableImplicitMT();
	// clock_t start, end;
	// double cpu_time_used;
	// start = clock();

	for (int ifile = 0; ifile < nInputFiles; ifile++){
		std::string input_file = input_files[ifile];
		ROOT::RDataFrame df_ss(tree_ss, input_file);
		ROOT::RDataFrame df_op(tree_op, input_file);
		
		// define columns needed for binning
		auto df_ss_updated = df_ss.Define("pt_avg",[](float x, float y){return static_cast<Float_t>((x + y)/2.);},{"m1.pt","m2.pt"}).Define("dphi_winded", [pi = ParamsSet::PI](float x){return (x > -pi / 2.)? static_cast<Float_t>(x) : static_cast<Float_t>(x + static_cast<Float_t>(pi) * 2.);}, {"dphi"}).Define("pair_pt_ptlead_ratio", [](float x, float y){return static_cast<Float_t>(x / y);}, {"pair_pt", "m1.pt"});
		auto df_op_updated = df_op.Define("pt_avg",[](float x, float y){return static_cast<Float_t>((x + y)/2.);},{"m1.pt","m2.pt"}).Define("dphi_winded", [pi = ParamsSet::PI](float x){return (x > -pi / 2.)? static_cast<Float_t>(x) : static_cast<Float_t>(x + static_cast<Float_t>(pi) * 2.);}, {"dphi"}).Define("pair_pt_ptlead_ratio", [](float x, float y){return static_cast<Float_t>(x / y);}, {"pair_pt", "m1.pt"});
				
		hist1d_rresultptrs_cur_file->clear();
		hist2d_rresultptrs_cur_file->clear();
		std::vector<TH1D> hist1Ds_cur_file {};
		std::vector<TH2D> hist2Ds_cur_file {};
		
		for (auto it = hist1Dinfo_list.begin(); it < hist1Dinfo_list.end(); it++){
			Hist1DMaking(*it, true, ifile, df_ss_updated);
			Hist1DMaking(*it, false, ifile, df_op_updated);
		}
		
		for (auto it = hist2Dinfo_list.begin(); it < hist2Dinfo_list.end(); it++){
			Hist2DMaking(*it, true, ifile, df_ss_updated);
			Hist2DMaking(*it, false, ifile, df_op_updated);
		}

		for (auto ptr : *hist1d_rresultptrs_cur_file){
			hist1Ds_cur_file.push_back(*ptr); // take action (run event loop) upon dereferencing the first pointer (only once for the entire function call)
		}

		for (auto ptr : *hist2d_rresultptrs_cur_file){
			hist2Ds_cur_file.push_back(*ptr); // take action (run event loop) upon dereferencing the first pointer (only once for the entire function call)
		}

		hist1Ds->push_back(hist1Ds_cur_file);
		hist2Ds->push_back(hist2Ds_cur_file);

		// info summary for the current for loop (input file)
		std::cout << "#event loop runs per data frame is: " << df_ss_updated.GetNRuns() << std::endl;
		std::cout << "#slots per data frame is: " << df_ss_updated.GetNSlots() << std::endl;
		std::cout << "# 1D histograms for the current input file: " << hist1Ds_cur_file.size() << std::endl;
		std::cout << "# 2D histograms for the current input file: " << hist2Ds_cur_file.size() << std::endl;

	}

	if (hist1Ds->size() != nInputFiles || hist2Ds->size() != nInputFiles){
		std::cout << "Size of hist1Ds and hist2Ds must equal #input files!" << std::endl;
		throw std::exception();
	}
	for (int ihist = 0; ihist < (*hist1Ds)[0].size(); ihist++){
		for (int ifile = 1; ifile < nInputFiles; ifile++)
		{
			(*hist1Ds)[0][ihist].Add(&(*hist1Ds)[ifile][ihist]);
		}
	}
	for (int ihist = 0; ihist < (*hist2Ds)[0].size(); ihist++){
		for (int ifile = 1; ifile < nInputFiles; ifile++)
		{
			(*hist2Ds)[0][ihist].Add(&(*hist2Ds)[ifile][ihist]);
		}
	}
	
	TFile * m_outfile=new TFile(output_file.c_str(), "recreate");

	// for (int ihist = 0; ihist < (*hist1Ds)[0].size(); ihist++){
	// 	(*hist1Ds)[0][ihist].Write();
	// }

	for (auto h : (*hist1Ds)[0]){
		h.Write();
	}
	for (auto h : (*hist2Ds)[0]){
		h.Write();
	}

	delete m_outfile;
	Finalize();
}

//______________________________________________________________________________
//                     RDFBasedMuonPairPlottingPbPb methods
//______________________________________________________________________________

void RDFBasedMuonPairPlottingPbPb::Initialize(){
	if (data_mode == DataMode::real) {
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_small_ctr_intvls.root");
		output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/hists_pbpb_data.root";
	}	
	else{
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/scrambled_muon_pairs_small_ctr_intvls.root");
		output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/hists_pbpb_scrambled.root";
		tree_ss = "scrambled_muon_pair_tree_sign1";
		tree_op = "scrambled_muon_pair_tree_sign2";
	} 								
	
	
	nInputFiles = input_files.size();
}



void RDFBasedMuonPairPlottingPbPb::Hist1DMaking(const Hist1DInfo * h1d, bool same_sign, int ifile, ROOT::RDF::RNode df){

	if (h1d->xaxis.nbins <= 0){
		std::cout << "Number of bins must be positive!" << std::endl;
		throw std::exception();
	}

	std::string hist_name_base;
	std::string hist_name;

	for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){ // pt bins are open-ended
		for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
			for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){

				ROOT::RDF::RNode df_filtered = df.Filter([ipt, &vec = std::as_const(ParamsSet::pTbins)](float x){return ApplySingleBinSelection<float>(x, ipt, vec, ipt == vec.size()-1);}, {"pt_avg"});
				df_filtered = df_filtered.Filter(ParamsSet::dphi_cut_funcs[idphi], {"dphi"});
				df_filtered = df_filtered.Filter(ParamsSet::deta_cut_funcs[ideta],{"deta"});

				for (unsigned int ictr = 0; ictr < ParamsSet::nCtrBins - 1; ictr++){ // ctr bins are close-ended
					hist_name_base = h1d->h1d_name_prefix + pms.pt_bin_labels[ipt] + pms.ctr_bin_labels[ictr] + pms.sign_labels[!same_sign] + pms.dphi_cut_labels[idphi] + pms.deta_cut_labels[ideta];
					df_filtered = df_filtered.Filter([ictr, &vec = std::as_const(ParamsSet::ctrbins)](int x){return ApplySingleBinSelection<int>(x, ictr, vec);}, {"avg_centrality"});
					
					FillSingleHist1D(h1d, hist_name_base, df_filtered, true);
				}
			}
		}
	}
}


void RDFBasedMuonPairPlottingPbPb::Hist2DMaking(const Hist2DInfo * h2d, bool same_sign, int ifile, ROOT::RDF::RNode df){

	if (h2d->xaxis.nbins <= 0 || h2d->yaxis.nbins <= 0){
		std::cout << "Number of bins must be positive!" << std::endl;
		throw std::exception();
	}

	std::string hist_name_base;


	for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){ // pt bins are open-ended
		for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
			for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){

				ROOT::RDF::RNode df_filtered = df.Filter([ipt, &vec = std::as_const(ParamsSet::pTbins)](float x){return ApplySingleBinSelection<float>(x, ipt, vec, ipt == vec.size()-1);}, {"pt_avg"});
				df_filtered = df_filtered.Filter(ParamsSet::dphi_cut_funcs[idphi], {"dphi"});
				df_filtered = df_filtered.Filter(ParamsSet::deta_cut_funcs[ideta],{"deta"});

				for (unsigned int ictr = 0; ictr < ParamsSet::nCtrBins - 1; ictr++){ // ctr bins are close-ended
					hist_name_base = h2d->h2d_name_prefix + pms.pt_bin_labels[ipt] + pms.ctr_bin_labels[ictr] + pms.sign_labels[!same_sign] + pms.dphi_cut_labels[idphi] + pms.deta_cut_labels[ideta];
					df_filtered = df_filtered.Filter([ictr, &vec = std::as_const(ParamsSet::ctrbins)](int x){return ApplySingleBinSelection<int>(x, ictr, vec);}, {"avg_centrality"});
					
					FillSingleHist2D(h2d, hist_name_base, df_filtered, true);
				}
			}
		}
	}
}


//______________________________________________________________________________
//                     RDFBasedMuonPairPlottingPP methods
//______________________________________________________________________________

void RDFBasedMuonPairPlottingPP::Initialize(){
	std::cout << "start Initialization." << std::endl;
	input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp.root");
	output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/hists_pp_data.root";
	nInputFiles = input_files.size();
}

void RDFBasedMuonPairPlottingPP::Hist1DMaking(const Hist1DInfo * h1d, bool same_sign, int ifile, ROOT::RDF::RNode df){

	// std::cout << "start hist1d making." << std::endl;
	if (h1d->xaxis.nbins <= 0){
		std::cout << "Number of bins must be positive!" << std::endl;
		throw std::exception();
	}

	std::string hist_name_base;
	std::string hist_name;
	for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){ // pt bins are open-ended
		for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
			for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){

				ROOT::RDF::RNode df_filtered = df.Filter([ipt, &vec = std::as_const(ParamsSet::pTbins)](float x){return ApplySingleBinSelection<float>(x, ipt, vec, ipt == vec.size()-1);}, {"pt_avg"});
				df_filtered = df_filtered.Filter(ParamsSet::dphi_cut_funcs[idphi], {"dphi"});
				df_filtered = df_filtered.Filter(ParamsSet::deta_cut_funcs[ideta],{"deta"});

				hist_name_base = h1d->h1d_name_prefix + pms.pt_bin_labels[ipt] + pms.sign_labels[!same_sign] + pms.dphi_cut_labels[idphi] + pms.deta_cut_labels[ideta];
				
				FillSingleHist1D(h1d, hist_name_base, df_filtered, true);
			}
		}
	}
}


void RDFBasedMuonPairPlottingPP::Hist2DMaking(const Hist2DInfo * h2d, bool same_sign, int ifile, ROOT::RDF::RNode df){

	if (h2d->xaxis.nbins <= 0 || h2d->yaxis.nbins <= 0){
		std::cout << "Number of bins must be positive!" << std::endl;
		throw std::exception();
	}

	std::string hist_name_base;


	for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){ // pt bins are open-ended
		for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
			for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){

				ROOT::RDF::RNode df_filtered = df.Filter([ipt, &vec = std::as_const(ParamsSet::pTbins)](float x){return ApplySingleBinSelection<float>(x, ipt, vec, ipt == vec.size()-1);}, {"pt_avg"});
				df_filtered = df_filtered.Filter(ParamsSet::dphi_cut_funcs[idphi], {"dphi"});
				df_filtered = df_filtered.Filter(ParamsSet::deta_cut_funcs[ideta],{"deta"});

				hist_name_base = h2d->h2d_name_prefix + pms.pt_bin_labels[ipt] + pms.sign_labels[!same_sign] + pms.dphi_cut_labels[idphi] + pms.deta_cut_labels[ideta];
				
				FillSingleHist2D(h2d, hist_name_base, df_filtered, true);
			}
		}
	}
}


//______________________________________________________________________________
//                     RDFBasedMuonPairPlottingPowheg methods
//______________________________________________________________________________

void RDFBasedMuonPairPlottingPowheg::Initialize(){
	std::cout << "The current powheg mode is " << powheg_mode << std::endl;
	switch (powheg_mode){
	case PowhegMode::bb:
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/muon_pairs_mc_truth_bb_1-5.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/muon_pairs_mc_truth_bb_6-10.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/muon_pairs_mc_truth_bb_11-15.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/muon_pairs_mc_truth_bb_16-20.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/muon_pairs_mc_truth_bb_21-25.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/muon_pairs_mc_truth_bb_26-30.root");
		output_file = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/hists_powheg_bb.root";
		nInputFiles = input_files.size();
		break;
	case PowhegMode::cc:
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/muon_pairs_mc_truth_cc_1-5.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/muon_pairs_mc_truth_cc_6-10.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/muon_pairs_mc_truth_cc_11-15.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/muon_pairs_mc_truth_cc_16-20.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/muon_pairs_mc_truth_cc_21-25.root");
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/muon_pairs_mc_truth_cc_26-30.root");
		output_file = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/hists_powheg_cc.root";
		nInputFiles = input_files.size();
		break;
	default:
		std::cout << "Powheg sample/mode must be bb or cc!" << std::endl;
		throw std::exception();
	}

	WeightCalc<Long64_t>();
}


template <typename Integer>
void RDFBasedMuonPairPlottingPowheg::WeightCalc(){
	
	Integer sum_nevents = 0;
	for (int ifile = 0; ifile < nInputFiles; ifile++){
		ROOT::RDataFrame df_meta("meta_tree", input_files[ifile]);
		sum_nevents += df_meta.Take<Integer>("nentries_before_cuts")->at(0);
	}

	powheg_weight_factor = static_cast<double>(1.) / sum_nevents;
}

void RDFBasedMuonPairPlottingPowheg::Hist1DMaking(const Hist1DInfo * h1d, bool same_sign, int ifile, ROOT::RDF::RNode df){

	ROOT::RDF::RNode df_mc_updated = df.Define("b_hadron_muon_dphi", HadronMuonDphiCalc, {"m1.phi", "m2.phi", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"});
	df_mc_updated = df_mc_updated.Define("c_hadron_muon_dphi", CHadronMuonDphiCalc, {"both_from_c", "m1.phi", "m2.phi", "m1_last_hf_hadron_prt_pt_eta_phi_m", "m2_last_hf_hadron_prt_pt_eta_phi_m"});
	df_mc_updated = df_mc_updated.Define("muon_b_jet_pt_ratio", MuonHQJetPtRatioCalc, {"m1.pt", "m2.pt", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"});
	df_mc_updated = df_mc_updated.Define("muon_c_jet_pt_ratio", MuonCJetPtRatioCalc, {"both_from_c", "m1.pt", "m2.pt", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"});
	if (h1d->xaxis.nbins <= 0){
		std::cout << "Number of bins must be positive!" << std::endl;
		throw std::exception();
	}

	std::string hist_name_base;
	// std::string hist_name;
	for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){ // pt bins are open-ended
		for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
			for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){

				ROOT::RDF::RNode df_filtered = df_mc_updated.Filter([ipt, &vec = std::as_const(ParamsSet::pTbins)](float x){return ApplySingleBinSelection<float>(x, ipt, vec, ipt == vec.size()-1);}, {"pt_avg"});
				df_filtered = df_filtered.Filter(ParamsSet::dphi_cut_funcs[idphi], {"dphi"});
				df_filtered = df_filtered.Filter(ParamsSet::deta_cut_funcs[ideta],{"deta"});

				df_filtered = df_filtered.Redefine("weight", [x = powheg_weight_factor](double w){return w * x;}, {"weight"});

				hist_name_base = h1d->h1d_name_prefix + pms.pt_bin_labels[ipt] + pms.sign_labels[!same_sign] + pms.dphi_cut_labels[idphi] + pms.deta_cut_labels[ideta];
				
				FillSingleHist1D(h1d, hist_name_base, df_filtered, true);
			}
		}
	}
}


void RDFBasedMuonPairPlottingPowheg::Hist2DMaking(const Hist2DInfo * h2d, bool same_sign, int ifile, ROOT::RDF::RNode df){

	if (h2d->xaxis.nbins <= 0 || h2d->yaxis.nbins <= 0){
		std::cout << "Number of bins must be positive!" << std::endl;
		throw std::exception();
	}

	std::string hist_name_base;


	for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){ // pt bins are open-ended
		for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
			for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){

				ROOT::RDF::RNode df_filtered = df.Filter([ipt, &vec = std::as_const(ParamsSet::pTbins)](float x){return ApplySingleBinSelection<float>(x, ipt, vec, ipt == vec.size()-1);}, {"pt_avg"});
				df_filtered = df_filtered.Filter(ParamsSet::dphi_cut_funcs[idphi], {"dphi"});
				df_filtered = df_filtered.Filter(ParamsSet::deta_cut_funcs[ideta],{"deta"});

				df_filtered = df_filtered.Redefine("weight", [x = powheg_weight_factor](double w){return w * x;}, {"weight"});
				
				hist_name_base = h2d->h2d_name_prefix + pms.pt_bin_labels[ipt] + pms.sign_labels[!same_sign] + pms.dphi_cut_labels[idphi] + pms.deta_cut_labels[ideta];

				FillSingleHist2D(h2d, hist_name_base, df_filtered, true);
			}
		}
	}
}

//______________________________________________________________________________
//                     RDFBasedMuonPairPlottingPythia methods
//______________________________________________________________________________


void RDFBasedMuonPairPlottingPythia::Initialize(){
	input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_allto0318.root");
	input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_after0322.root");
	output_file = "/usatlas/u/yuhanguo/usatlasdata/pythia/hists_pythia.root";
	nInputFiles = input_files.size();
	WeightCalc<Long64_t>();
}

template <typename Integer>
void RDFBasedMuonPairPlottingPythia::WeightCalc(){
	if (pythia_weight_factors.size() != 0){
		std::cout << "vector must be empty!" << std::endl;
		throw std::exception();
	}

	int nKinRanges = 5;
	std::vector<std::vector<Integer>> nevents_before_cuts = {{20000, 99910}, {50000, 400000}, {1000000, 1500000}, {1000000, 1500000}, {1000000, 1500000}};
		// much easier to have an array of size nKinRanges * nInputFiles then the other way around

    if (nevents_before_cuts.size() != nKinRanges){
        std::cout << "The vectors nevents_before_cuts must have size equal to the number of kinematic ranges" << std::endl;
        throw std::exception();
    }else{
        for (auto nevents_cur_kin_range : nevents_before_cuts){
            if (nevents_cur_kin_range.size() != nInputFiles){
                std::cout << "All vectors in nevents_before_cuts must have the same size: ";
                std::cout << "equal to the number of batches/files to be added." << std::endl;
                throw std::exception();
            }
        }
    }

    for (int ikin = 0; ikin < nKinRanges; ikin++){
    	std::vector<Integer> nevents_cur_kin_range = nevents_before_cuts[ikin];
        Integer sum_nevents_cur_kin_range = std::accumulate(nevents_cur_kin_range.begin(),nevents_cur_kin_range.end(),0);
        std::vector<double> weight_factors_cur_kin_range {};
        for (int jfile = 0; jfile < nInputFiles; jfile++){
            weight_factors_cur_kin_range.push_back(static_cast<double>(nevents_cur_kin_range[jfile]) / sum_nevents_cur_kin_range);
        }
        pythia_weight_factors.push_back(weight_factors_cur_kin_range);
    }

	for (auto v : pythia_weight_factors){
		cout << "pythia_weight_factors size: " << v.size() << endl << "elements: ";
		for (auto vv : v){
			cout << vv << " ";
		}
		cout << endl;
	}
}

template <typename WeightType, typename PTHatType>
WeightType RDFBasedMuonPairPlottingPythia::PTHatBinnedWeightCorrection(WeightType weight_raw, PTHatType pTHat, int ifile, const std::vector<std::vector<WeightType>> & weight_factors){
	for (typename std::vector<PTHatType>::iterator it = ParamsSet::pTHatbins_pythia.begin(); it < ParamsSet::pTHatbins_pythia.end()-1; it++){
		if (pTHat >= *it && pTHat < *(it+1)){
			return weight_raw * weight_factors[it - ParamsSet::pTHatbins_pythia.begin()][ifile];
		}
	}
	std::cout << "Input pTHat value " << pTHat << " is not in the generated range!" << std::endl;
	throw std::exception();
	return -1;
}

void RDFBasedMuonPairPlottingPythia::Hist1DMaking(const Hist1DInfo * h1d, bool same_sign, int ifile, ROOT::RDF::RNode df){

	ROOT::RDF::RNode df_mc_updated = df.Define("b_hadron_muon_dphi", HadronMuonDphiCalc, {"m1.phi", "m2.phi", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"});
	df_mc_updated = df_mc_updated.Define("c_hadron_muon_dphi", CHadronMuonDphiCalc, {"both_from_c", "m1.phi", "m2.phi", "m1_last_hf_hadron_prt_pt_eta_phi_m", "m2_last_hf_hadron_prt_pt_eta_phi_m"});
	df_mc_updated = df_mc_updated.Define("muon_b_jet_pt_ratio", MuonHQJetPtRatioCalc, {"m1.pt", "m2.pt", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"});
	df_mc_updated = df_mc_updated.Define("muon_c_jet_pt_ratio", MuonCJetPtRatioCalc, {"both_from_c", "m1.pt", "m2.pt", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"});
	
	if (h1d->xaxis.nbins <= 0){
		std::cout << "Number of bins must be positive!" << std::endl;
		throw std::exception();
	}

	std::string hist_name_base;
	// std::string hist_name;
	for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){ // pt bins are open-ended
		for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
			for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){

				ROOT::RDF::RNode df_filtered = df_mc_updated.Filter([ipt, &vec = std::as_const(ParamsSet::pTbins)](float x){return ApplySingleBinSelection<float>(x, ipt, vec, ipt == vec.size()-1);}, {"pt_avg"});
				df_filtered = df_filtered.Filter(ParamsSet::dphi_cut_funcs[idphi], {"dphi"});
				df_filtered = df_filtered.Filter(ParamsSet::deta_cut_funcs[ideta],{"deta"});

				df_filtered = df_filtered.Redefine("weight", [ifile, &vec = std::as_const(pythia_weight_factors)](double w, double pthat){return RDFBasedMuonPairPlottingPythia::PTHatBinnedWeightCorrection<double, double>(w, pthat, ifile, vec);}, {"weight", "pTHat"});
				
				hist_name_base = h1d->h1d_name_prefix + pms.pt_bin_labels[ipt] + pms.sign_labels[!same_sign] + pms.dphi_cut_labels[idphi] + pms.deta_cut_labels[ideta];

				FillSingleHist1D(h1d, hist_name_base, df_filtered, true);
			}
		}
	}
}


void RDFBasedMuonPairPlottingPythia::Hist2DMaking(const Hist2DInfo * h2d, bool same_sign, int ifile, ROOT::RDF::RNode df){

	if (h2d->xaxis.nbins <= 0 || h2d->yaxis.nbins <= 0){
		std::cout << "Number of bins must be positive!" << std::endl;
		throw std::exception();
	}

	std::string hist_name_base;


	for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){ // pt bins are open-ended
		for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
			for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){

				ROOT::RDF::RNode df_filtered = df.Filter([ipt, &vec = std::as_const(ParamsSet::pTbins)](float x){return ApplySingleBinSelection<float>(x, ipt, vec, ipt == vec.size()-1);}, {"pt_avg"});
				df_filtered = df_filtered.Filter(ParamsSet::dphi_cut_funcs[idphi], {"dphi"});
				df_filtered = df_filtered.Filter(ParamsSet::deta_cut_funcs[ideta],{"deta"});

				df_filtered = df_filtered.Redefine("weight", [ifile, &vec = std::as_const(pythia_weight_factors)](double w, double pthat){return RDFBasedMuonPairPlottingPythia::PTHatBinnedWeightCorrection<double, double>(w, pthat, ifile, vec);}, {"weight", "pTHat"});
				
				hist_name_base = h2d->h2d_name_prefix + pms.pt_bin_labels[ipt] + pms.sign_labels[!same_sign] + pms.dphi_cut_labels[idphi] + pms.deta_cut_labels[ideta];

				FillSingleHist2D(h2d, hist_name_base, df_filtered, true);
			}
		}
	}
}
