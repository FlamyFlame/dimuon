#include "RecoEffyRetRespPlotter.cxx"
#include <TROOT.h>

void plot_det_resp_N_reco_effcy_powheg_fullsim_pp17(){ // detector response & reco efficiencies
    gROOT->SetBatch(kTRUE);
    DetRespPlotter pl;
    pl.Run();
}
