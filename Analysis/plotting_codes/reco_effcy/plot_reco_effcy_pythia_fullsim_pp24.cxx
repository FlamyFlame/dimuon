#include "PythiaFullsimRecoEffPlotter.cxx"
#include <TROOT.h>

// pp24 Pythia-fullsim reco-efficiency + detector-response plots.
//   is_test_sample : true = the small TEST sample (4 isospin beams); false = the FULL production
//                    (pp beam only, isospin weight 1). It selects the input DIR *and* the "_full"
//                    hist-file suffix, and the plots land in that sample's own plots/ dir, so the
//                    two sets never clobber each other.
//   tight_WP       : NOMINAL muon WP = Tight.
void plot_reco_effcy_pythia_fullsim_pp24(bool is_test_sample = true, bool tight_WP = true){
    gROOT->SetBatch(kTRUE);
    PythiaFullsimRecoEffPlotter pl(tight_WP);
    pl.is_test_sample = is_test_sample;
    pl.Run();
}
