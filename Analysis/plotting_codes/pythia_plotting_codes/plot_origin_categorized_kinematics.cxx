#include "PythonCategorizedPlottingBaseClass.cxx"


class PythonOriginCategorizedPlotting : public PythonCategorizedPlottingBaseClass{
protected:
    void initialize();
    void fill_line_map();
    void fill_thstack_order_map();
public:

    bool turn_single_b_resonance_on = true;

    PythonOriginCategorizedPlotting(std::string kin_in, bool projx_2d_in, bool projy_2d_in, bool staggered_in, bool norm_unity_in, std::string kin1d_in, std::string kin_title_in)
            : PythonCategorizedPlottingBaseClass(kin_in, projx_2d_in, projy_2d_in, staggered_in, norm_unity_in, kin1d_in, kin_title_in){}
    PythonOriginCategorizedPlotting(std::string kin_in, bool projx_2d_in, bool projy_2d_in, bool staggered_in, bool norm_unity_in, std::string kin1d_in, std::string kin_title_in, std::vector<std::array<float,2>> cuts_in)
            : PythonCategorizedPlottingBaseClass(kin_in, projx_2d_in, projy_2d_in, staggered_in, norm_unity_in, kin1d_in, kin_title_in, cuts_in){}
    
    ~PythonOriginCategorizedPlotting(){}
};

void PythonOriginCategorizedPlotting::initialize(){
    PythonCategorizedPlottingBaseClass::initialize();
    subdir_name = "origin_categoried/";
    if (!turn_single_b_resonance_on){
        optional_suffix = "_no_single_b_no_res";
    }
}

void PythonOriginCategorizedPlotting::fill_line_map(){
line_map[muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_zero_hard_scatt] = new line({"_GS_ISR_no_HS", "GS in ISR no HS", kOrange, kOrange});
line_map[muon_pair_both_from_open_HF_origin_catgr::diff_gs_same_hard_scatt]   = new line({"_diff_GS_same_HS", "diff GS same HS", kGray, kGray});
line_map[muon_pair_both_from_open_HF_origin_catgr::others]            = new line({"_others", "others", kViolet+1, kViolet+1});
line_map[muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_one_hard_scatt]  = new line({"_gs_ISR_one_HS", "flavor excitation (FE)", kBlue, kBlue});
line_map[muon_pair_both_from_open_HF_origin_catgr::fc]              = new line({"_FC", "flavor creation (FC)", kRed, kRed});
line_map[muon_pair_both_from_open_HF_origin_catgr::same_gs_fsr]         = new line({"_gs_FSR", "GS in FSR", kGreen + 2, kGreen});
line_map[muon_pair_both_from_open_HF_origin_catgr::nOrigins + pair_flavor_index::from_single_b]   = new line({"_single_b", "single b", kBlack, kYellow});
line_map[muon_pair_both_from_open_HF_origin_catgr::nOrigins + pair_flavor_index::from_resonance]  = new line({"_resonance", "resonances", kCyan+1, kCyan+1});}

void PythonOriginCategorizedPlotting::fill_thstack_order_map(){
    thstack_order_map[sign::same_sign] = {
        muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_zero_hard_scatt, 
        muon_pair_both_from_open_HF_origin_catgr::diff_gs_same_hard_scatt, 
        muon_pair_both_from_open_HF_origin_catgr::others, 
        muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_one_hard_scatt, 
        muon_pair_both_from_open_HF_origin_catgr::fc, 
        muon_pair_both_from_open_HF_origin_catgr::same_gs_fsr
    };
    thstack_order_map[sign::op_sign] = {  
        muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_zero_hard_scatt, 
        muon_pair_both_from_open_HF_origin_catgr::diff_gs_same_hard_scatt, 
        muon_pair_both_from_open_HF_origin_catgr::others, 
        muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_one_hard_scatt, 
        muon_pair_both_from_open_HF_origin_catgr::fc, 
        muon_pair_both_from_open_HF_origin_catgr::same_gs_fsr
    };
    if (turn_single_b_resonance_on){
        thstack_order_map[sign::same_sign].push_back(muon_pair_both_from_open_HF_origin_catgr::nOrigins + pair_flavor_index::from_single_b);
        thstack_order_map[sign::same_sign].push_back(muon_pair_both_from_open_HF_origin_catgr::nOrigins + pair_flavor_index::from_resonance);
        thstack_order_map[sign::op_sign].push_back(muon_pair_both_from_open_HF_origin_catgr::nOrigins + pair_flavor_index::from_single_b);
        thstack_order_map[sign::op_sign].push_back(muon_pair_both_from_open_HF_origin_catgr::nOrigins + pair_flavor_index::from_resonance);
    }
}


void plot_origin_categorized_kinematics(){
    ParamsSet pms;
    bool with_data_resonance_cuts = false;

    std::vector<PythonOriginCategorizedPlotting*> origin_plotting_list = {};
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR", false, false, false, false, "DR", "#Delta R"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR", false, false, true, false, "DR", "#Delta R")); // accumulative
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR", false, false, false, true, "DR", "#Delta R")); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR_zoomin", false, false, false, false, "DR_zoomin", "#Delta R"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR_zoomin", false, false, true, false, "DR_zoomin", "#Delta R")); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR_zoomin", false, false, false, true, "DR_zoomin", "#Delta R")); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Deta_zoomin", false, false, false, false, "Deta_zoomin", "#Delta #eta"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Deta_zoomin", false, false, true, false, "Deta_zoomin", "#Delta #eta")); // accumulative
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Dphi_zoomin", false, false, false, false, "Dphi_zoomin", "#Delta #phi"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Dphi_zoomin", false, false, true, false, "Dphi_zoomin", "#Delta #phi")); // accumulative
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR_jacobian_corrected", false, false, false, false, "DR_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR_jacobian_corrected", false, false, true, false, "DR_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}")); // accumulative
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR_jacobian_corrected", false, false, false, true, "DR_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}")); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR_zoomin_jacobian_corrected", false, false, false, false, "DR_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR_zoomin_jacobian_corrected", false, false, true, false, "DR_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}")); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("DR_zoomin_jacobian_corrected", false, false, false, true, "DR_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}")); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("pt_asym", false, false, false, false, "pt_asym", "A"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("pt_asym", false, false, true, false, "pt_asym", "A")); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("pt_asym", false, false, false, true, "pt_asym", "A")); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("psrapidity_ordered_pt_asym", false, false, false, false, "psrapidity_ordered_pt_asym", "#Tilde{A}"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("psrapidity_ordered_pt_asym", false, false, true, false, "psrapidity_ordered_pt_asym", "#Tilde{A}")); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("psrapidity_ordered_pt_asym", false, false, false, true, "psrapidity_ordered_pt_asym", "#Tilde{A}")); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("pair_pt_ptlead_ratio", false, false, false, false, "pair_pt_ptlead_ratio", "p_{T}^{pair} / p_{T}^{lead}"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("pair_pt_ptlead_ratio", false, false, true, false, "pair_pt_ptlead_ratio", "p_{T}^{pair} / p_{T}^{lead}")); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("pair_pt_ptlead_ratio", false, false, false, true, "pair_pt_ptlead_ratio", "p_{T}^{pair} / p_{T}^{lead}")); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("ptlead_pair_pt", true, false, false, false, "pair_pt", "p_{T}^{pair}"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("ptlead_pair_pt", true, false, true, false, "pair_pt", "p_{T}^{pair}")); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("ptlead_pair_pt", true, false, false, true, "pair_pt", "p_{T}^{pair}")); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("ptlead_pair_pt", false, true, false, false, "ptlead", "p_{T}^{lead}"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("ptlead_pair_pt", false, true, true, false, "ptlead", "p_{T}^{lead}")); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("ptlead_pair_pt", false, true, false, true, "ptlead", "p_{T}^{lead}")); // norm to unity

    std::vector<std::array<float,2>> minv_cuts;
    if (with_data_resonance_cuts) minv_cuts = pms.minv_cuts;

    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt", false, true, false, false, "minv", "m_{#mu#mu}", minv_cuts));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt", false, true, true, false, "minv", "m_{#mu#mu}", minv_cuts)); // accumulative
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt", false, true, false, true, "minv", "m_{#mu#mu}", minv_cuts)); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt_zoomin", false, true, false, false, "minv_zoomin", "m_{#mu#mu}", minv_cuts));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt_zoomin", false, true, true, false, "minv_zoomin", "m_{#mu#mu}", minv_cuts)); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt_zoomin", false, true, false, true, "minv_zoomin", "m_{#mu#mu}", minv_cuts)); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt_jacobian_corrected", false, true, false, false, "minv_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt_jacobian_corrected", false, true, true, false, "minv_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts)); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt_jacobian_corrected", false, true, false, true, "minv_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts)); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt_zoomin_jacobian_corrected", false, true, false, false, "minv_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt_zoomin_jacobian_corrected", false, true, true, false, "minv_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts)); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("minv_pair_pt_zoomin_jacobian_corrected", false, true, false, true, "minv_zoomin_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{dm_{#mu#mu}}", minv_cuts)); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Deta_Dphi", true, false, false, false, "Dphi", "#Delta #phi"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Deta_Dphi", true, false, true, false, "Dphi", "#Delta #phi")); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Deta_Dphi", true, false, false, true, "Dphi", "#Delta #phi")); // norm to unity
    
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Deta_Dphi", false, true, false, false, "Deta", "#Delta #eta"));
    origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Deta_Dphi", false, true, true, false, "Deta", "#Delta #eta")); // accumulative
    // origin_plotting_list.push_back(new PythonOriginCategorizedPlotting("Deta_Dphi", false, true, false, true, "Deta", "#Delta #eta")); // norm to unity

    for (auto& origin_plot : origin_plotting_list){
        origin_plot->with_data_resonance_cuts = with_data_resonance_cuts;
        origin_plot->Run();
        delete origin_plot;
    }

    PythonOriginCategorizedPlotting DR_no_single_b_no_res ("DR", false, false, false, false, "DR", "#Delta R");
    DR_no_single_b_no_res.turn_single_b_resonance_on = false;
    DR_no_single_b_no_res.Run();

    PythonOriginCategorizedPlotting DR_jacobian_corrected_no_single_b_no_res ("DR_jacobian_corrected", false, false, false, false, "DR_jacobian_corrected", "FULL: #frac{1}{#Delta R} #frac{d#sigma}{d#Delta R}");
    DR_jacobian_corrected_no_single_b_no_res.turn_single_b_resonance_on = false;
    DR_jacobian_corrected_no_single_b_no_res.Run();

    PythonOriginCategorizedPlotting minv_no_single_b_no_res ("minv_pair_pt", false, true, false, false, "minv", "m_{#mu#mu}", minv_cuts);
    minv_no_single_b_no_res.turn_single_b_resonance_on = false;
    minv_no_single_b_no_res.Run();

}
