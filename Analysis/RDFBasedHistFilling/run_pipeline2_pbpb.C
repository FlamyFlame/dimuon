#include "RDFBasedHistFillingData.h"

void run_pipeline2_pbpb(int year = 25, bool coarse = false) {
    gROOT->SetBatch(kTRUE);

    RDFBasedHistFillingPbPb pbpb(year, false);
    pbpb.trigger_mode = 1;
    pbpb.hist_filling_cycle = RDFBasedHistFillingData::generic;
    pbpb.output_generic_hists = false;
    pbpb.useCoarseQEtaBin = coarse;
    pbpb.mindR_trig = 0.02;
    pbpb.Run();
}
