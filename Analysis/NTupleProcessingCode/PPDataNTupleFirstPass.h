#pragma once

#include "DimuonDataAnalysisBaseClass.c"

class PPDataNTupleFirstPass : public DimuonDataAnalysisBaseClass{
    // Read through the N-tuple, apply appropriate cuts
    // Then fill in histograms and/or output trees
    // mode = 3: kill resonances & output single-muon information into a TTree
    // mode = 4: kill resonances & output muon-pair information into a TTree
    // NOW ONLY HAVE MODE = 3, 4


private:
      
// --------------------- class methods ---------------------------
  
    void TChainFill() override;    

public :
    PPDataNTupleFirstPass();
    ~PPDataNTupleFirstPass(){}
    void Run() override;
};


