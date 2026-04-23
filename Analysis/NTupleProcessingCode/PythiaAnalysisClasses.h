#ifndef __PythiaAnalysisClasses_h__
#define __PythiaAnalysisClasses_h__

#include "../MuonObjectsParamsAndHelpers/Muon.h"
#include "../MuonObjectsParamsAndHelpers/MuonPairPythia.h"
#include "PythiaAlgCoreT.c"
#include "PythiaTruthExtras.c"
#include "PythiaFullSimExtras.c"
#include "PythiaFullSimOverlayExtras.c"

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

// Full-sim pp24 analysis: truth ancestry tracing + reco muon matching
class PythiaFullSimAnalysis
  : public PythiaAlgCoreT<
        MuonPairPythiaFullSimWTruth, MuonPythiaFullSimWTruth,
        PythiaFullSimAnalysis,
        PythiaFullSimExtras<MuonPairPythiaFullSimWTruth, MuonPythiaFullSimWTruth, PythiaFullSimAnalysis>,
        PythiaTruthExtras<MuonPairPythiaFullSimWTruth, PythiaFullSimAnalysis>
    >
  , public PythiaFullSimExtras<MuonPairPythiaFullSimWTruth, MuonPythiaFullSimWTruth, PythiaFullSimAnalysis>
  , public PythiaTruthExtras<MuonPairPythiaFullSimWTruth, PythiaFullSimAnalysis>
{
public:
    PythiaFullSimAnalysis(int batch_num_input = 0, bool use_local = false)
        : PythiaAlgCoreT(batch_num_input, use_local)
    {
        this->isPrivate = false;
        this->E_COM = 5.36;
        this->run_year = 24;
        this->fullsim_sample_type = FullSimSampleType::pp;
    }
};

// Placeholder: full-sim PbPb overlay analysis (to be implemented when overlay samples available)
class PythiaFullSimOverlayAnalysis
  : public PythiaAlgCoreT<
        MuonPairPythiaFullSimOverlayWTruth, MuonPythiaFullSimOverlayWTruth,
        PythiaFullSimOverlayAnalysis,
        PythiaFullSimExtras<MuonPairPythiaFullSimOverlayWTruth, MuonPythiaFullSimOverlayWTruth, PythiaFullSimOverlayAnalysis>,
        PythiaFullSimOverlayExtras<PythiaFullSimOverlayAnalysis>,
        PythiaTruthExtras<MuonPairPythiaFullSimOverlayWTruth, PythiaFullSimOverlayAnalysis>
    >
  , public PythiaFullSimExtras<MuonPairPythiaFullSimOverlayWTruth, MuonPythiaFullSimOverlayWTruth, PythiaFullSimOverlayAnalysis>
  , public PythiaFullSimOverlayExtras<PythiaFullSimOverlayAnalysis>
  , public PythiaTruthExtras<MuonPairPythiaFullSimOverlayWTruth, PythiaFullSimOverlayAnalysis>
{
public:
    PythiaFullSimOverlayAnalysis(FullSimSampleType sample_type = FullSimSampleType::hijing,
                                 int batch_num_input = 0, bool use_local = false)
        : PythiaAlgCoreT(batch_num_input, use_local)
    {
        this->isPrivate = false;
        this->E_COM = 5.36;
        this->run_year = 24;
        this->fullsim_sample_type = sample_type;
    }
};

#endif
