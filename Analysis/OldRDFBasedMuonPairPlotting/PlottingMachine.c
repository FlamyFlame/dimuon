#include <TROOT.h>
// #include <TChain.h>
#include <TFile.h>
#include  <stdlib.h>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include<algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "ParamsSet.h"


// still need to implement:
// 1D / projection from 2D + 2D + rebinning + scaling/normalization + legend / subplot title + png naming








std::vector<std::string> PlottingMachine::dataset_types = {"pbpb_real", "pbpb_scr", "pp", "powheg_bb", "powheg_cc", "pythia"};
// std::vector<std::string> PlottingMachine::category_names = {"ptavg", "ctr", "sign", "deta", "dphi", "photoprodcut", "gapcut", "flavor", "origin", "data"};
std::vector<std::string> PlottingMachine::category_names = {"ptavg", "ctr", "sign", "deta", "dphi", "photoprodcut", "gapcut", "data"};
std::vector<std::string> PlottingMachine::levels = {"canvas", "row", "col", "line"};

std::map<std::string, DataSet*> 	PlottingMachine::dataset_map;




void PlottingMachine::DatasetTypenameToStructMapImpl(){
	dataset_map["pbpb_real"] = new DataSet({{"ptavg", "ctr", "sign", "deta", "dphi", "photoprodcut", "gapcut"}, "data/hists_pbpb_data.root" nullptr});
	dataset_map["pbpb_scr"] = new DataSet({{"ptavg", "ctr", "sign", "deta", "dphi", "photoprodcut", "gapcut"}, "data/hists_pbpb_scrambled.root" nullptr});
	dataset_map["pp"] = new DataSet({{"ptavg", "sign", "deta", "dphi", "photoprodcut", "gapcut"}, "data/hists_pp_data.root" nullptr});
	// dataset_map["powheg_bb"] = new DataSet({{"ptavg", "sign", "deta", "dphi", "flavor", "origin"}, "data/hists_powheg_bb.root" nullptr});
	dataset_map["powheg_bb"] = new DataSet({{"ptavg", "sign", "deta", "dphi"}, "data/hists_powheg_bb.root" nullptr});
	// dataset_map["powheg_cc"] = new DataSet({{"ptavg", "sign", "deta", "dphi", "flavor", "origin"}, "data/hists_powheg_cc.root" nullptr});
	dataset_map["powheg_cc"] = new DataSet({{"ptavg", "sign", "deta", "dphi"}, "data/hists_powheg_cc.root" nullptr});
	// dataset_map["pythia"] = new DataSet({{"ptavg", "sign", "deta", "dphi", "flavor", "origin"}, "data/hists_pythia.root" nullptr});
	dataset_map["pythia"] = new DataSet({{"ptavg", "sign", "deta", "dphi"}, "data/hists_pythia.root" nullptr});
}

void PlottingMachine::DefaultCategoryMapImpl(){
	category_map["ptavg"] 			= new KinematicCategory({{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {}, ParamsSet::pt_bin_labels});
	category_map["ctr"]				= new KinematicCategory({{0, 1, 2, 3, 4, 5}, {}, ParamsSet::ctr_bin_labels});
	category_map["sign"]			= new KinematicCategory({{0, 1}, {}, ParamsSet::sign_labels});
	category_map["deta"]			= new KinematicCategory({{0, 1}, {}, ParamsSet::deta_cut_labels});
	category_map["dphi"]			= new KinematicCategory({{0, 1}, {}, ParamsSet::dphi_cut_labels});
	category_map["photoprodcut"]	= new KinematicCategory({{0}, {}, ParamsSet::photocut_labels});
	category_map["gapcut"]			= new KinematicCategory({{0}, {}, ParamsSet::gapcut_labels});
	// category_map["flavor"]		= // NOT implemented yet
	// category_map["origin"]		= // NOT implemented yet

}

PlottingMachine::Initialize(){
	std::cout << "Require the following inputs from child class implementation:" << std::endl;
	std::cout << "(1) observable to be plotted" << std::endl;
	std::cout << "(2) categories to turn splitting on" << std::endl;
	std::cout << "(3) input datasets' typenames (with an additional vector of reclustering if needed) " << std::endl;
	std::cout << "(4) map from the levels (canvas, row, column, lines) to categories" << std::endl;
	std::cout << "(5) observables to be plotted (along with the basic plot settings), combined into a vector of instances of the Hist1D structure" << std::endl;
	std::cout << "This code works only if the binning for each category shared among the datasets are identical." << std::endl;
	std::cout << "Also, it cannot handle yet reclustering of kinematic binning." << std::endl;

	// implement upon construction the maps essential to ensure correct automation 
	DatasetTypenameToStructMapImpl();
	DefaultCategoryMapImpl();
}


PlottingMachine::PlottingMachine(){
	Initialize();
}


void PlottingMachine::LevelDimCalc(){

	// calculate level "dimension"
	int data_dim = (input_dataset_indices_reclustered.size() == 0)? input_dataset_typenames.size() : input_dataset_indices_reclustered.size();
	
	for (std::string s : levels){
		if (level_to_category_map[s] == ""){
			level_dim_map[s] = 1;
		} else if (level_to_category_map[s] == "data"){
			level_dim_map[s] = data_dim;
		} else{
			std:vector<int> bins_on = category_map[level_to_category_map[s]]->bins_turned_on;
			level_dim_map[s] = bins_on.size();
		}
	}

	// calculate level multiplicative scale
	for (std::string s : {"canvas", "row", "col", "line"}){
    	if (level_to_category_map[s] == ""){
    		level_to_multpl_scale_map[s] = 0;
    	} else if (level_to_category_map[s] == "data"){
            level_to_multpl_scale_map[s] = prod_nbins_final[categories_to_turn_on.size()];
        } else{
            int mapped_index = std::find(categories_to_turn_on.begin(), categories_to_turn_on.end(), level_to_category_map[s]) - categories_to_turn_on.begin();
            if (mapped_index == categories_to_turn_on.size()){
                std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
                throw std::exception();
            }
            level_to_multpl_scale_map[s] = prod_nbins_final[mapped_index];
        }
    }
}


void PlottingMachine::GetHistsPerDataset(std::string dataset){
    std::vector<string> categories_cur_dataset = dataset_map[dataset]->existing_categories; // get the strings for the categories in the data set for future use

    std::vector<int> nbins = {}; // to record, for each category in the current data set, # bins turned on
    std::vector<int> prod_nbins; // to record the product of #bins for all previous categories --> needed for flattening/expanding

    for (int i = 0; i < categories_cur_dataset.size(); i++){
        nbins.push_back(category_map[categories_cur_dataset[i]]->bins_turned_on.size());
        prod_nbins.push_back((i == 0)? 1 : prod_nbins[i-1] * nbins[i-1]);
    }
    int total_prod_nbins = std::accumulate(nbins.begin(), nbins.end(), 1, std::multiplies<int>());

    std::vector<std::vector<TH1D*>> * flattened_raw_hist_list_list_cur_dataset = new std::vector<std::vector<TH1D*>>; // to store all histograms (un-summed) in the current dataset
    for (int iobserv = 0; iobserv < hist1d_info_list.size(); iobserv++){
    	flattened_raw_hist_list_list_cur_dataset->push_back({});
    }

    // ----------------------------------------------------------------------------------------------------------------
    // loop over all flattened (1D) indices & fill in the 1D vector of histograms
    for (int iraw = 0; iraw < total_prod_nbins; iraw++){

        // transform the flattened 1D vector index back to the multi-dimensional vector indices (dimension = #categories existing for the current dataset)
        std::string hname_raw_postfix = ""; // to store the histogram name
        if (iraw < 30) cout << "1D index: " << iraw << ", multi-d indices: ";
        
        for (int icatgry = 0; icatgry < nbins.size(); icatgry++){ // loop over the categories
            int multi_d_index = iraw / prod_nbins[icatgry] % nbins[icatgry];
            
            if (iraw < 30) cout << multi_d_index << " ";
            Category* cur_catgry = category_map[categories_cur_dataset[icatgry]];
            hname_raw_postfix += cur_catgry->bin_labels.at(cur_catgry->bins_turned_on[multi_d_index]); // refering to bins_turned_on is necessary since bin_labels contain labels for ALL bins, but only a subset is turned on
        }

        if (iraw < 30) cout << endl;

        for (int iobserv = 0; iobserv < hist1d_info_list.size(); iobserv++){
        	std::string observable_name = hist1d_info_list[iobserv]->kin;
	        if (projx_2d){
	            TH2D* h2d = (TH2D*) dataset_map[dataset]->f->Get(("h_" + hist1d_info_list[iobserv]->kin + hname_raw_postfix).c_str());
	            flattened_raw_hist_list_list_cur_dataset[iobserv]->push_back((TH1D*) h2d->ProjectionX(("hx_" + hist1d_info_list[iobserv]->kin + hname_raw_postfix).c_str()));
	        } else if (projy_2d){
	            TH2D* h2d = (TH2D*) dataset_map[dataset]->f->Get(("h_" + hist1d_info_list[iobserv]->kin + hname_raw_postfix).c_str());
	            flattened_raw_hist_list_list_cur_dataset[iobserv]->push_back((TH1D*) h2d->ProjectionY(("hy_" + hist1d_info_list[iobserv]->kin + hname_raw_postfix).c_str()));
	        } else{
	            flattened_raw_hist_list_list_cur_dataset[iobserv]->push_back((TH1D*) dataset_map[dataset]->f->Get(("h_" + hist1d_info_list[iobserv]->kin + hname_raw_postfix).c_str()));
	        }
	    }
    }

    // ----------------------------------------------------------------------------------------------------------------
    // finding positions of the turned-on categories among all categories in the current dataset
    std::vector<int> flattened_indices_cur_hist = {}; // to store the indices of all sub-histograms that contribute to each histogram
    std::vector<int> indices_in_category_list = {};

    // find the indices of the turned-on categories
    for (int icatgry = 0; icatgry < categories_to_turn_on.size(); icatgry++){
        int ind = std::find(categories_cur_dataset.begin(), categories_cur_dataset.end(), categories_to_turn_on[icatgry]) - categories_cur_dataset.begin();
        if (ind == categories_cur_dataset.size()){ // not found
            std::cout << "The current category (" << categories_to_turn_on[icatgry] << ") turned on is NOT found in the current dataset (" << dataset << ")." << std::endl;
            throw std::exception();
        }
        
        indices_in_category_list.push_back(ind);
    }

    // ----------------------------------------------------------------------------------------------------------------
    // loop over the flattened-1D-vector indices of the TARGET (FINAL) histograms
    // translate each index into corresponding n-D vector indices (n = # categories turned on)
    // also find the positions i1, ..., in of the n turned-on categories in the vector of the m categories in the current dataset
    // loop over the raw histogram indices, change each one into corresponding m-D vector indice s(m = # categories in the current dataset)
    // record the flattened indices of the raw histogram to sum over for each final histogram (sharing the same bin_i1, ..., bin_in)
    for (ind_summed = 0; ind_summed < total_nhists_final_per_dataset; ind_summed++){
        std::vector<int> indices_raw_to_sum = {};
        for (int ind_raw = 0; ind_raw < total_prod_nbins; ind_raw){
            bool include_cur_raw = true;
            for (int icatgry = 0; icatgry < categories_to_turn_on.size(); icatgry++){
                int index_in_category_list = indices_in_category_list[icatgry];
                int target_ind = ind_summed / prod_nbins_final[icatgry] % nbins_final[icatgry]; // the target bin-index of the i-th turned-on category in the current final histogram to produce
                int actual_ind = ind_raw / prod_nbins[index_in_category_list] % nbins[index_in_category_list]; // the actual bin-index of the i-th turned-on category in the current raw histogram being looked at

                if (target_ind != actual_ind){ // if not match: the current raw histogram is not to include in the sum
                    include_cur_raw = false;
                    break; // no need to look at the rest of the categories
                }
            }
            
            if (include_cur_raw){
                indices_raw_to_sum.push_back(ind_raw);
            }
        }
        
        // sanity check - the number of histograms to sum should equal the product of #bins of the summed-over dimensions
        int total_nhists_final_per_dataset = std::accumulate(nbins_final.begin(), nbins_final.begin() + categories_to_turn_on.size(), 1, std::multiplies<int>());
        if (indices_raw_to_sum.size() != total_prod_nbins / total_nhists_final_per_dataset){
            std::cout << "The number of raw histograms to be summed for each result histogram should be " << total_prod_nbins / total_nhists_final_per_dataset << std::endl;
            std::cout << "The calculated # raw histograms to include is " << indices_raw_to_sum.size() << std::endl;
            throw std::exception();
        }

		for (int iobserv = 0; iobserv < hist1d_info_list.size(); iobserv++){
	        TH1D* cur_hist_summed = (*flattened_raw_hist_list_list_cur_dataset)[iobserv][0];
	        for (int ind = 1; ind < indices_raw_to_sum.size(); ind++){
	            cur_hist_summed->Add((*flattened_raw_hist_list_list_cur_dataset)[iobserv][ind]); // perform the sum
	            delete (*flattened_raw_hist_list_list_cur_dataset)[iobserv][ind]; // delete unneeded pointers
	        }

	        flattened_final_hist_list_list[iobserv]->push_back(cur_hist_summed);
	    }
    }

    delete flattened_raw_hist_list_list_cur_dataset; // delete the pointer to the flattened vector of raw histograms
}


void PlottingMachine::DatasetReclusteringHandling(){
    int total_nhists_final_per_dataset = std::accumulate(nbins_final.begin(), nbins_final.begin() + categories_to_turn_on.size(), 1, std::multiplies<int>());

    if (input_dataset_indices_reclustered.size() != 0){ // if reclustering required
        std::vector<std::vector<TH1D*>> * flattened_final_hist_list_list_reclustered = new std::vector<std::vector<TH1D*>>; // dimension: # datafiles
        for (cluster : input_dataset_indices_reclustered){
            if (cluster.size() == 0){
                std::cout << "All input-dataset clusters must have nonzero size!" << std::endl;
                throw std::exception();
            }
            if (cluster.size() > 1){
            	for (int iobserv = 0; iobserv < hist1d_info_list.size(); iobserv++){
	                for (int ihist = 0; ihist < total_nhists_final_per_dataset; ihist++){
	                    for (int ids = 1; ids < cluster.size(); ids++){
	                        (*flattened_final_hist_list_list[iobserv])[cluster[0] * total_nhists_final_per_dataset + ihist]->Add((*flattened_final_hist_list_list[iobserv])[cluster[ids] * total_nhists_final_per_dataset + ihist]);
	                        delete (*flattened_final_hist_list_list[iobserv])[cluster[ids] * total_nhists_final_per_dataset + ihist];
	                    }
	                }
	            }
            }
            for (int iobserv = 0; iobserv < hist1d_info_list.size(); iobserv++){
	            flattened_final_hist_list_list_reclustered[iobserv]->insert(flattened_final_hist_list_list_reclustered[iobserv]->end(), flattened_final_hist_list_list[iobserv]->begin() + cluster[0] * total_nhists_final_per_dataset, flattened_final_hist_list_list[iobserv]->begin() + (cluster[0] + 1) * total_nhists_final_per_dataset);
    		}
        }

        delete flattened_final_hist_list_list; // delete the old pointer
        flattened_final_hist_list_list = flattened_final_hist_list_list_reclustered; // change the address to that of the reclustered list
    }
}


void PlottingMachine::GetHists(){

    // -------------------------- in preparation for flattening & expanding --------------------------

    flattened_final_hist_list_list = new std::vector<std::vector<TH1D*>>;
    for (int iobserv = 0; iobserv < hist1d_info_list.size(); iobserv++){
    	flattened_final_hist_list_list->push_back({});
    }


    for (int icatgry = 0; icatgry < categories_to_turn_on.size(); icatgry++){

        int nbins_cur_catgry = category_map[categories_to_turn_on[icatgry]]->bins_turned_on.size();
        nbins_final.push_back(nbins_cur_catgry);
        prod_nbins_final.push_back((icatgry == 0)? 1 : prod_nbins_final[icatgry-1] * nbins_final[icatgry-1]);
    }

    nbins_final.push_back(input_dataset_typenames.size());
    prod_nbins_final.push_back(total_nhists_final_per_dataset);


    for (std::string dataset : input_dataset_typenames){ // loop over the input data sets
        GetHistsPerDataset(dataset);
    }

    DatasetReclusteringHandling();
}


void PlottingMachine::Plotting(){
    std::vector<int> level_index = {};
    int nkinmetic_catgries = 0;

    for (std::string s : {"canvas", "row", "col", "line"}){
    	if (level_to_category_map[s] == ""){
    		level_to_multpl_scale_map[s] = 0;
    	} else if (level_to_category_map[s] == "data"){
            level_to_multpl_scale_map[s] = prod_nbins_final[categories_to_turn_on.size()];
        } else{
            int mapped_index = std::find(categories_to_turn_on.begin(), categories_to_turn_on.end(), level_to_category_map[s]) - categories_to_turn_on.begin();
            if (mapped_index == categories_to_turn_on.size()){
                std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
                throw std::exception();
            }
            level_to_multpl_scale_map[s] = prod_nbins_final[mapped_index];
        }
    }

    for (int icanvas = 0; icanvas < level_dim_map["canvas"]; icanvas++){
		// tcanvas_list->push_back(new TCanvas)
    	TCanvas* c = new TCanvas(Form("c%d",icanvas), Form("c%d",icanvas), 750 * level_dim_map["col"], 500 * level_dim_map["row"]);
		c->Divide(level_dim_map["row"], level_dim_map["col"]);

		for (int irow = 0; irow < level_dim_map["row"]; irow++){
			for (int icol = 0; icol < level_dim_map["col"]; icol++){

                c->cd(irow * level_dim_map["col"] + icol + 1);
                gPad->SetLeftMargin(0.2);
                gPad->SetBottomMargin(0.135);
                // gPad->SetTopMargin(0.18);
                gPad->SetLogx(logx);

                // TLegend* l = new TLegend(0.24,0.6,0.5,0.89);
                TLegend* l;
                if (powheg_scale == 1. && pythia_scale == 1.){
                    if (irow == 0)  l = new TLegend(0.73,0.6,0.89,0.89);
                    else            l = new TLegend(0.7,0.6,0.91,0.89);        
                }else{
                    if (irow == 0)  l = new TLegend(0.68,0.6,0.89,0.89);
                    else            l = new TLegend(0.65,0.6,0.91,0.89);        
                }

                l->SetBorderSize(0);
                l->SetFillStyle(0);
                l->SetTextFont(42);
                // l->SetTextSize(20);
                l->SetMargin(0.2);
                l->SetTextColor(1);

				for (int iline = 0; iline < level_dim_map["line"]; iline++){
					/// 重要！！！ -1 ！！！！！
					int index_flattened = icanvas * level_to_multpl_scale_map["canvas"] + irow * level_to_multpl_scale_map["row"] + icol * level_to_multpl_scale_map["col"] + iline * level_to_multpl_scale_map["line"];
					flattened_final_hist_list_list[index_flattened]->Draw();
				}
			}
		}

		c->SaveAs(Form("%s%s_mc_data_compr_near_away.png", plot_output_path.c_str(), kin1d.c_str()));
		c->Clear();
		delete c;
	}

}


void PlottingMachine::Finalize(){
	// delete the pointers associated with maps
	for (std::string ds : input_dataset_typenames){
		delete dataset_map[ds]->f;
	}

	for (std::string ds : dataset_types){
		delete dataset_map[ds];
	}

	for (std::string catgry : category_names){
		delete category_map[catgry];
	}

}


void PlottingMachine::Run(){

	InputInit();
	CategoryStructsAdjust();


	if (observable.size() == 0 || input_dataset_typenames.size() == 0 || level_to_category_map.size() == 0){
		// categories_to_turn_on.size() == 0 is allowed: in this case, we make a single canvas with a single subplot and a single line
		std::cout << "Input insufficient! Must provide (1) observable to be plotted, ";
		std::cout << "(2) categories to turn splitting on, (3) input datasets' typenames, ";
		std::cout << "(4) map from the levels (canvas, row, column, lines) to categories" << std::endl;
		throw std::exception();
	}

	//------------------------------------------------------------------------------------------------------------------------------


	for (std::string ds : input_dataset_typenames){
		dataset_map[ds]->f = new TFile*(dataset_map[ds]->filename, "read");
	}

	// calculate the total number of final (i.e, with "invisible" binning dimensions summed over) histograms for each dataset
	total_nhists_final_per_dataset = 1; // this number is identical for all datasets and only depend on the categories turned on
	auto prod_increment = [&map = this->category_map, &total_prod = this->total_nhists_final_per_dataset](string catgry){total_prod *= (map[catgry]->bins_turned_on).size();};
	for_each(categories_to_turn_on.begin(), categories_to_turn_on.end(), prod_increment);

	LevelDimCalc();

	GetHists();
	Plotting();

	Finalize();
}





