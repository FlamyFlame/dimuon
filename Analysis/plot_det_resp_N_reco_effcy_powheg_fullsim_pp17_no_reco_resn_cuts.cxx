#include "RecoEffyRetRespPlotter.cxx"

void plot_det_resp_N_reco_effcy_powheg_fullsim_pp17_no_reco_resn_cuts(){ // detector response & reco efficiencies
    DetRespPlotter pl(17, false, true);
    pl.Run();
}
