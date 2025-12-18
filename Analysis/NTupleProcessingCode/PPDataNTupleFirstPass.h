#pragma once

#include "DimuonDataAnalysisBaseClass.c"

class PPDataNTupleFirstPass : public DimuonDataAnalysisBaseClass{

private:
      
// --------------------- class methods ---------------------------
  
    void TChainFill() override;    

public :
    explicit PPDataNTupleFirstPass(int run_year_input, int file_batch_input);
    ~PPDataNTupleFirstPass(){}
};


