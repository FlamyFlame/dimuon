// Test crossx histogram filling for all three datasets
void test_crossx_all() {
    gROOT->SetBatch(kTRUE);
    
    // Load compiled libraries
    gROOT->LoadMacro("RDFBasedHistFillingPP.cxx++O");
    gROOT->LoadMacro("RDFBasedHistFillingPbPb.cxx++O");
    
    // Test pp24
    {
        std::cout << "\n===== Testing pp24 with trigger_mode==2 =====" << std::endl;
        RDFBasedHistFillingPP pp24(24, false);
        pp24.trigger_mode = 2;
        pp24.FillHistograms();
        std::cout << "[SUCCESS] pp24 crossx test completed!" << std::endl;
    }
    
    // Test pbpb23
    {
        std::cout << "\n===== Testing pbpb23 with trigger_mode==2 =====" << std::endl;
        RDFBasedHistFillingPbPb pbpb23(23, false);
        pbpb23.trigger_mode = 2;
        pbpb23.FillHistograms();
        std::cout << "[SUCCESS] pbpb23 crossx test completed!" << std::endl;
    }
    
    // Test pbpb24
    {
        std::cout << "\n===== Testing pbpb24 with trigger_mode==2 =====" << std::endl;
        RDFBasedHistFillingPbPb pbpb24(24, false);
        pbpb24.trigger_mode = 2;
        pbpb24.FillHistograms();
        std::cout << "[SUCCESS] pbpb24 crossx test completed!" << std::endl;
    }
    
    std::cout << "\n===== All tests completed =====" << std::endl;
}
