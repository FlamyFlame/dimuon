#include "SingleMuEffcyPtTurnOnFitter.cxx"

void single_muon_trig_effcy_pT_fitting_pbpb (){
    SingleMuEffcyPtTurnOnFitterPbPb* fitter = new SingleMuEffcyPtTurnOnFitterPbPb("include_upc");
    fitter->fitting_mode = SingleMuEffcyPtTurnOnFitterPbPb::fermi_plus_log;
    fitter->Run();
}
