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
    PythiaTruthAnalysis(int batch_num_input, bool is_private, double e_com = 5.36,
                        bool use_local = false, bool pp_only = false)
        : PythiaAlgCoreT(batch_num_input, use_local)
    {
        if (std::abs(e_com - 5.36) > 0.01 && std::abs(e_com - 5.02) > 0.01)
            throw std::runtime_error("PythiaTruthAnalysis: E_COM must be 5.02 or 5.36, got "
                + std::to_string(e_com));
        this->isPrivate = is_private;
        this->E_COM = e_com;
        if (pp_only) {
            this->only_pp_isospin = true;
            this->extra_output_suffix = "_pp_only";
        }
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
    // Sample + isospin are BOTH driven by `isTestSample` (default false = the FULL production;
    // see PythiaAlgCoreT.h and FullSimSampleType.h). Do NOT force the isospin here -- that would
    // override the isTestSample-derived default and silently give the 4-beam TEST sample a
    // pp-only weight.
    //   isTestSample=false -> pp full sample, pp beam only, isospin weight 1   (the physics case)
    //   isTestSample=true  -> pp24 test sample, 4 beams, Pb ratio 4:6:6:9      (produced by mistake)
    // `pp_only` = force the pp beam alone regardless (the legacy pp-only cross-check on the
    // 4-beam TEST sample); it also tags the output "_pp_only".
    // `sample_type` is pp for the nominal pp24 fullsim. The ONLY other admissible value is
    // FullSimSampleType::noovl (the r17663 no-overlay diagnostic): it has no overlaid event,
    // so it needs the pp class (no overlay Extras, no centrality) while carrying its own
    // input dir / file tag / label. Overlay samples must use PythiaFullSimOverlayAnalysis.
    PythiaFullSimAnalysis(int batch_num_input = 0, bool use_local = false, bool pp_only = false,
                          FullSimSampleType sample_type = FullSimSampleType::pp)
        : PythiaAlgCoreT(batch_num_input, use_local)
    {
        if (FullSimSampleIsOverlay(sample_type))
            throw std::runtime_error("PythiaFullSimAnalysis: sample_type is an OVERLAY sample "
                "-- use PythiaFullSimOverlayAnalysis (it needs the overlay Extras: centrality, "
                "FCal, HIJING truth).");
        this->isPrivate = false;
        this->E_COM = 5.36;
        this->run_year = 24;
        this->fullsim_sample_type = sample_type;
        if (pp_only) {
            this->setIsospinBeams(false);   // escape hatch: pp beam alone, weight 1
            this->extra_output_suffix = "_pp_only";
        }
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
    // Beam content: the HIJING overlay simulates Pb+Pb, whose nucleons are a p/n mix, so
    // the 4 isospin beams {pp,pn,np,nn} combined with the Pb ratio 4:6:6:9 are the
    // DEFAULT (inherited from FullSimSampleIsOverlay).  The overlay TEST sample on disk
    // has ONLY the pp beam, so runs over it must opt out with setIsospinBeams(false);
    // the 4-beam overlay FULL sample now in production uses the default.
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
