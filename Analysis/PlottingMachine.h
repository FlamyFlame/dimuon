#ifndef __PlottingMachine_h__
#define __PlottingMachine_h__

struct DataSet{
	std::vector<std::string> existing_categories;
	std::string filename;
	TFile* f;
	float scale = 1.;
};

struct KinematicCategory{
	// bool if_exist = false; // if exists: e.g, false for ctr for all datasets but the PbPb ones (real / scrambled) 
	std::vector<int> bins_turned_on;
	std::vector<std::vector<int>> bins_turned_on_reclustered = {};
	std::vector<std::string> bin_labels; // expect to be defined in ParamsSet.h
	std::vector<std::string> bin_titles; // expect to be defined in ParamsSet.h
};

struct Hist1D{
	std::string kin;
	bool projx_2d;
	bool projy_2d;
	bool staggered;
	bool norm_unity;
	std::string kin1d;
	std::string kin_title;
	bool logx = false;
	bool logy = false;
};

struct Hist2D{
	std::string kin2d;
	// std::string kinx_title = "";
	// std::string kiny_title = "";
	bool logx = false;
	bool logy = false;
	bool norm_unity = false;
};

class PlottingMachine{
public:
	
	// ------------------------------------------------------------------------------------------------------------------------------
	// ----------------------------------- input required ---------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------------------------------------------------
	std::string observable = "";
	std::string plot_output_path = "plots/"; // e.g, "plots/mc_data_compr/"
	std::vector<std::string> 		categories_to_turn_on = {};
	std::vector<std::string> 		input_dataset_typenames = {};
	std::vector<float> 				input_dataset_scales = {};
	std::vector<std::vector<int>> 	input_dataset_indices_reclustered = {}; // e.g, {{"pp"}, {"powheg_bb", "powheg_cc"}, {"pythia"}}

	std::map<std::string, std::string> level_to_category_map;
	std::map<std::string, int> level_to_multpl_scale_map;

	std::vector<Hist1D*> hist1d_info_list;

	// ------------------------------------------------------------------------------------------------------------------------------
	// ----------------------------------- helper member variables and maps to be implemented automatically  ------------------------
	// ----------------------------------- if needed, can be modified in child classes ----------------------------------------------
	// ------------------------------------------------------------------------------------------------------------------------------

	static std::vector<string> dataset_types;
	static std::vector<string> category_names;
	static std::vector<string> levels;
	std::map<std::string, int> level_dim_map;

    std::vector<int> nbins_final;
    std::vector<int> prod_nbins_final;
	int total_nhists_final_per_dataset;

    std::vector<std::vector<TH1D*>> * flattened_final_hist_list_list = nullptr;

	std::vector<std::vector<std::vector<std::vector<TH1D*>>>> hist1d_list;
	std::vector<std::vector<std::vector<std::vector<TH2D*>>>> hist2d_list;
	// std::vector<std::vector<std::vector<std::vector<TH1*>>>> hist_list;

	static std::map<std::string, DataSet*> dataset_map;
	std::map<std::string, KinematicCategory*> category_map;

	std::map<std::pair<std::string, std::string>, int> category_index_in_dataset_map;
	// Remark: we do NOT need a category_exist_map; if category_split_on_map contains a certain <dataset, category_name> pair as a valid key, it already shows that the category exists for this dataset
	
	// ------------------------------------------------------------------------------------------------------------------------------
	// ----------------------------------- member functions -------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------------------------------------------------

	PlottingMachine();
	~PlottingMachine(){}

	void Initialize();
	static void DatasetTypenameToStructMapImpl();
	void DefaultCategoryMapImpl();

	virtual void InputInit() = 0;
	virtual void CategoryStructsAdjust() = 0;

	void LevelDimCalc();
	void GetHists();
	void GetHistsPerDataset(std::string dataset);
	void DatasetReclusteringHandling();
	void Plotting();
	void Finalize();
	void Run();

};


#endif