{
    gROOT->ProcessLine(".L DataAnalysisClasses.h");
    PPAnalysis* pp_24 = new PPAnalysis(24, 1);
    pp_24->nevents_max = 200;
    pp_24->trigger_mode = 1;
    pp_24->Run();
    delete pp_24;
}
