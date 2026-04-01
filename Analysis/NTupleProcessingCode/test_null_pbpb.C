{
    gROOT->ProcessLine(".L DataAnalysisClasses.h");
    PbPbAnalysis* pbpb_23 = new PbPbAnalysis(23, 1);
    pbpb_23->nevents_max = 200;
    pbpb_23->trigger_mode = 1;
    pbpb_23->Run();
    delete pbpb_23;
}
