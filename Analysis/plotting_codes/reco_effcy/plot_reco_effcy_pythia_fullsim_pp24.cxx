#include "PythiaFullsimRecoEffPlotter.cxx"
#include <TROOT.h>

void plot_reco_effcy_pythia_fullsim_pp24(){
    gROOT->SetBatch(kTRUE);
    PythiaFullsimRecoEffPlotter pl;
    pl.Run();
}
