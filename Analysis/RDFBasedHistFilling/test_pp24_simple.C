void test_pp24_crossx()
{
    gROOT->SetBatch(kTRUE);
    
    std::cout << "\n===== Compiling classes =====" << std::endl;
    gROOT->LoadMacro("RDFBasedHistFillingPP.cxx++O");
    
    std::cout << "\n===== Creating pp24 filler =====" << std::endl;
    RDFBasedHistFillingPP filler(24, false);
    
    // Disable mindR requirement (input files don't have _mindR_0_02 suffix)
    filler.mindR_trig = -1;
    
    std::cout << "Setting trigger_mode to 2..." << std::endl;
    filler.trigger_mode = 2;
    
    std::cout << "Running FillHistograms..." << std::endl;
    filler.FillHistograms();
    
    std::cout << "\n[SUCCESS] pp24 crossx test completed!" << std::endl;
}
