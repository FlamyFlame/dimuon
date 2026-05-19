#include "RDFBasedHistFillingData.h"

void run_pipeline3_pbpb(int year = 25) {
    gROOT->SetBatch(kTRUE);

    RDFBasedHistFillingPbPb pbpb(year, false);
    pbpb.trigger_mode = 1;
    pbpb.hist_filling_cycle = RDFBasedHistFillingData::inv_weight_by_single_mu_effcy;
    pbpb.useCoarseQEtaBin = true;
    pbpb.mindR_trig = 0.02;
    pbpb.Run();
}
