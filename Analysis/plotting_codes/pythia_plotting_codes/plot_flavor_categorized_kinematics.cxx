#include "PythonCategorizedPlottingBaseClass.cxx"


class PythonFlavorCategorizedPlotting : public PythonCategorizedPlottingBaseClass{
protected:
    void initialize();
    void fill_line_map();
    void fill_thstack_order_map();
public:
    bool turn_single_b_resonance_on = true;

    PythonFlavorCategorizedPlotting(std::string kin_in, bool projx_2d_in, bool projy_2d_in, bool staggered_in, bool norm_unity_in, std::string kin1d_in, std::string kin_title_in)
            : PythonCategorizedPlottingBaseClass(kin_in, projx_2d_in, projy_2d_in, staggered_in, norm_unity_in, kin1d_in, kin_title_in){}
    PythonFlavorCategorizedPlotting(std::string kin_in, bool projx_2d_in, bool projy_2d_in, bool staggered_in, bool norm_unity_in, std::string kin1d_in, std::string kin_title_in, std::vector<std::array<float,2>> cuts_in)
            : PythonCategorizedPlottingBaseClass(kin_in, projx_2d_in, projy_2d_in, staggered_in, norm_unity_in, kin1d_in, kin_title_in, cuts_in){}
    
    ~PythonFlavorCategorizedPlotting(){}
};

void PythonFlavorCategorizedPlotting::initialize(){
    PythonCategorizedPlottingBaseClass::initialize();
    subdir_name = "flavor_categoried/";
    if (!turn_single_b_resonance_on){
        optional_suffix = "_no_single_b_no_res";
    }
}

void PythonFlavorCategorizedPlotting::fill_line_map(){
    line_map[pair_flavor_index::from_resonance] = new line({"_resonance", "resonances", kCyan+1, kCyan+1});
    line_map[pair_flavor_index::resonance_contaminated] = new line({"_resonance_contaminated", "resonance contaminated", kPink-4, kPink-4});
    line_map[pair_flavor_index::from_single_b] = new line({"_single_b", "single b", kOrange, kYellow});
    line_map[pair_flavor_index::bb] = new line({"_bb", "bb", kBlue, kBlue});
    line_map[pair_flavor_index::cc] = new line({"_cc", "cc", kRed, kRed});
    line_map[pair_flavor_index::one_b_one_c] = new line({"_one_b_one_c", "one b one c", kGreen+2, kGreen+2});
    line_map[pair_flavor_index::photon_to_dimuon_splitting] = new line({"_photon_splitting", "photon splitting", kViolet, kViolet});
    line_map[pair_flavor_index::other_flavor] = new line({"_other_flavors", "other flavors", kGray, kGray});
}

void PythonFlavorCategorizedPlotting::fill_thstack_order_map(){
    thstack_order_map[sign::same_sign] = {  
        pair_flavor_index::other_flavor, 
        pair_flavor_index::one_b_one_c,
        pair_flavor_index::cc, 
        pair_flavor_index::resonance_contaminated, 
        pair_flavor_index::bb
    };

    thstack_order_map[sign::op_sign] = {
        pair_flavor_index::other_flavor,
        // pair_flavor_index::drell_yan,
        pair_flavor_index::resonance_contaminated, 
        pair_flavor_index::one_b_one_c,
        pair_flavor_index::cc, 
        pair_flavor_index::bb,
    };

    if (turn_single_b_resonance_on){
        thstack_order_map[sign::same_sign].push_back(pair_flavor_index::from_single_b);
        thstack_order_map[sign::same_sign].push_back(pair_flavor_index::from_resonance);

        thstack_order_map[sign::op_sign].push_back(pair_flavor_index::from_single_b);
        thstack_order_map[sign::op_sign].push_back(pair_flavor_index::from_resonance);
        thstack_order_map[sign::op_sign].push_back(pair_flavor_index::photon_to_dimuon_splitting);
    }
}

void plot_flavor_categorized_kinematics(){
    ParamsSet pms;
    bool with_data_resonance_cuts = false;

    std::vector<PythonFlavorCategorizedPlotting*> flavor_plotting_list = {};
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR", false, false, false, false, "DR", "#Delta R"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR", false, false, true, false, "DR", "#Delta R")); // accumulative
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR", false, false, false, true, "DR", "#Delta R")); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR_zoomin", false, false, false, false, "DR_zoomin", "#Delta R"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR_zoomin", false, false, true, false, "DR_zoomin", "#Delta R")); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR_zoomin", false, false, false, true, "DR_zoomin", "#Delta R")); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Deta_zoomin", false, false, false, false, "Deta_zoomin", "#Delta #eta"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Deta_zoomin", false, false, true, false, "Deta_zoomin", "#Delta #eta")); // accumulative
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Dphi_zoomin", false, false, false, false, "Dphi_zoomin", "#Delta #phi"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Dphi_zoomin", false, false, true, false, "Dphi_zoomin", "#Delta #phi")); // accumulative
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR_jacobian_corrected", false, false, false, false, "DR_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR_jacobian_corrected", false, false, true, false, "DR_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}")); // accumulative
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR_jacobian_corrected", false, false, false, true, "DR_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}")); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR_zoomin_jacobian_corrected", false, false, false, false, "DR_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR_zoomin_jacobian_corrected", false, false, true, false, "DR_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}")); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("DR_zoomin_jacobian_corrected", false, false, false, true, "DR_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}")); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("pt_asym", false, false, false, false, "pt_asym", "A"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("pt_asym", false, false, true, false, "pt_asym", "A")); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("pt_asym", false, false, false, true, "pt_asym", "A")); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("psrapidity_ordered_pt_asym", false, false, false, false, "psrapidity_ordered_pt_asym", "#Tilde{A}"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("psrapidity_ordered_pt_asym", false, false, true, false, "psrapidity_ordered_pt_asym", "#Tilde{A}")); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("psrapidity_ordered_pt_asym", false, false, false, true, "psrapidity_ordered_pt_asym", "#Tilde{A}")); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("pair_pt_ptlead_ratio", false, false, false, false, "pair_pt_ptlead_ratio", "p_{T}^{pair} / p_{T}^{lead}"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("pair_pt_ptlead_ratio", false, false, true, false, "pair_pt_ptlead_ratio", "p_{T}^{pair} / p_{T}^{lead}")); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("pair_pt_ptlead_ratio", false, false, false, true, "pair_pt_ptlead_ratio", "p_{T}^{pair} / p_{T}^{lead}")); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("ptlead_pair_pt", true, false, false, false, "pair_pt", "p_{T}^{pair}"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("ptlead_pair_pt", true, false, true, false, "pair_pt", "p_{T}^{pair}")); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("ptlead_pair_pt", true, false, false, true, "pair_pt", "p_{T}^{pair}")); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("ptlead_pair_pt", false, true, false, false, "ptlead", "p_{T}^{lead}"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("ptlead_pair_pt", false, true, true, false, "ptlead", "p_{T}^{lead}")); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("ptlead_pair_pt", false, true, false, true, "ptlead", "p_{T}^{lead}")); // norm to unity

    std::vector<std::array<float,2>> minv_cuts;
    if (with_data_resonance_cuts) minv_cuts = pms.minv_cuts;

    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt", false, true, false, false, "minv", "m_{#mu#mu}", minv_cuts));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt", false, true, true, false, "minv", "m_{#mu#mu}", minv_cuts)); // accumulative
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt", false, true, false, true, "minv", "m_{#mu#mu}", minv_cuts)); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_zoomin", false, true, false, false, "minv_zoomin", "m_{#mu#mu}", minv_cuts));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_zoomin", false, true, true, false, "minv_zoomin", "m_{#mu#mu}", minv_cuts)); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_zoomin", false, true, false, true, "minv_zoomin", "m_{#mu#mu}", minv_cuts)); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_jacobian_corrected", false, true, false, false, "minv_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_jacobian_corrected", false, true, true, false, "minv_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts)); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_jacobian_corrected", false, true, false, true, "minv_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts)); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_zoomin_jacobian_corrected", false, true, false, false, "minv_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_zoomin_jacobian_corrected", false, true, true, false, "minv_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts)); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_zoomin_jacobian_corrected", false, true, false, true, "minv_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts)); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Deta_Dphi", true, false, false, false, "Dphi", "#Delta #phi"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Deta_Dphi", true, false, true, false, "Dphi", "#Delta #phi")); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Deta_Dphi", true, false, false, true, "Dphi", "#Delta #phi")); // norm to unity
    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Deta_Dphi", false, true, false, false, "Deta", "#Delta #eta"));
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Deta_Dphi", false, true, true, false, "Deta", "#Delta #eta")); // accumulative
    // flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("Deta_Dphi", false, true, false, true, "Deta", "#Delta #eta")); // norm to unity

    for (auto& flavor_plot : flavor_plotting_list){
        flavor_plot->with_data_resonance_cuts = with_data_resonance_cuts;
        flavor_plot->Run();
        delete flavor_plot;
    }

    PythonFlavorCategorizedPlotting DR_no_single_b_no_res ("DR", false, false, false, false, "DR", "#Delta R");
    DR_no_single_b_no_res.turn_single_b_resonance_on = false;
    DR_no_single_b_no_res.Run();

    PythonFlavorCategorizedPlotting DR_jacobian_corrected_no_single_b_no_res ("DR_jacobian_corrected", false, false, false, false, "DR_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}");
    DR_jacobian_corrected_no_single_b_no_res.turn_single_b_resonance_on = false;
    DR_jacobian_corrected_no_single_b_no_res.Run();

    PythonFlavorCategorizedPlotting minv_no_single_b_no_res ("minv_pair_pt", false, true, false, false, "minv", "m_{#mu#mu}", minv_cuts);
    minv_no_single_b_no_res.turn_single_b_resonance_on = false;
    minv_no_single_b_no_res.Run();

}
