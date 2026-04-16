#include "PowhegFullsimDetRespPlotterSingleMuon.cxx"
#include <TROOT.h>

void plot_det_resp_N_reco_effcy_powheg_fullsim_pp17_single_muon(){ // detector response & reco efficiencies, single muon
    gROOT->SetBatch(kTRUE);
    PowhegFullsimDetRespPlotterSingleMuon pl(17, false); // run_year=17, tight_WP=false (medium)
    pl.Run();
}
