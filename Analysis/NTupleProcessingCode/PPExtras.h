#pragma once

#include "DimuonDataAnalysisBaseClass.c"

class PPExtras : public DimuonDataAnalysisBaseClass{

private:
      
// --------------------- class methods ---------------------------
  
    void TChainFill();
    void InitParamsExtra();

public :
    explicit PPExtras(int run_year_input, int file_batch_input);
    ~PPExtras(){}
};


class PPAnalysis
  : public DimuonDataAlgCoreT<
        MuonPairPP, MuonPP,
        PPAnalysis,
        PPExtras<PPAnalysis>
    >
  , public PPExtras<PPAnalysis>
{};
