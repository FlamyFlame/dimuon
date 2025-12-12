#include "SingleMuEffcyPtTurnOnFitter.cxx"

void single_muon_trig_effcy_pT_fitting_pp (){
    SingleMuEffcyPtTurnOnFitterPP* fitter = new SingleMuEffcyPtTurnOnFitterPP();
    fitter->fitting_mode = SingleMuEffcyPtTurnOnFitterPP::erf_plus_log;
    fitter->Run();
}
