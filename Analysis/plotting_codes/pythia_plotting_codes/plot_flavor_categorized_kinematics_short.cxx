#include "plot_flavor_categorized_kinematics.cxx"


void plot_flavor_categorized_kinematics_short(){
    ParamsSet pms;
    bool with_data_resonance_cuts = false;

    std::vector<PythonFlavorCategorizedPlotting*> flavor_plotting_list = {};
    
    std::vector<std::array<float,2>> minv_cuts;
    if (with_data_resonance_cuts) minv_cuts = pms.minv_cuts;

    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt", false, true, true, false, "minv", "m_{#mu#mu}", minv_cuts)); // accumulative    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_10GeV_pair_pt", false, true, true, false, "minv_10GeV", "m_{#mu#mu}", minv_cuts)); // accumulative    
    flavor_plotting_list.push_back(new PythonFlavorCategorizedPlotting("minv_pair_pt_zoomin", false, true, true, false, "minv_zoomin", "m_{#mu#mu}", minv_cuts)); // accumulative
    
    for (auto& flavor_plot : flavor_plotting_list){
        flavor_plot->with_data_resonance_cuts = with_data_resonance_cuts;
        flavor_plot->Run();
        delete flavor_plot;
    }
}
