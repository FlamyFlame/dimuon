#include "PowhegFullsimRecoEffPlotter.cxx"

void plot_det_resp_N_reco_effcy_powheg_fullsim_pp17_require_signal_cuts(){ // detector response & reco efficiencies
    PowhegFullsimRecoEffPlotter pl(17, true, true);  // NOMINAL WP = Tight (pass false for Medium)
    pl.Run();
}
