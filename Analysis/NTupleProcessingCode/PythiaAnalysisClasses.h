#ifndef __PythiaAnalysisClasses_h__
#define __PythiaAnalysisClasses_h__

#include "../MuonObjectsParamsAndHelpers/Muon.h"
#include "../MuonObjectsParamsAndHelpers/MuonPairPythia.h"
#include "PythiaAlgCoreT.c"
#include "PythiaTruthExtras.c"

class PythiaTruthAnalysis
  : public PythiaAlgCoreT<
        MuonPairPythiaTruth, MuonPythiaTruth,
        PythiaTruthAnalysis,
        PythiaTruthExtras<MuonPairPythiaTruth, PythiaTruthAnalysis>
    >
  , public PythiaTruthExtras<MuonPairPythiaTruth, PythiaTruthAnalysis>
{
public:
    PythiaTruthAnalysis(int batch_num_input, bool is_private, double e_com = 5.36, bool use_local = false)
        : PythiaAlgCoreT(batch_num_input, use_local)
    {
        if (std::abs(e_com - 5.36) > 0.01 && std::abs(e_com - 5.02) > 0.01)
            throw std::runtime_error("PythiaTruthAnalysis: E_COM must be 5.02 or 5.36, got "
                + std::to_string(e_com));
        this->isPrivate = is_private;
        this->E_COM = e_com;
    }
};

#endif
