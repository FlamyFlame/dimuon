#pragma once

#include "RDFBasedHistFillingBaseClass.cxx"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include <cmath>
#include <map>
#include <string>
#include <vector>

class RDFBasedHistFillingPythia : public virtual RDFBasedHistFillingBaseClass{
protected:
    bool isPrivate{false};
    bool with_data_resonance_cuts{false};
    double E_COM{5.36};

    bool is_fullsim{false};
    bool is_fullsim_overlay{false};

    void            InitializePythiaImpl();
    void            InitializePythiaCommon();
    virtual void    InitializePythiaExtra(){}

    virtual void    InitAnalysisSettingsHook() override{ return InitializePythiaImpl(); }
    virtual void    SetIOPathsHook() override;

    void            CreateBaseRDFsPythiaImpl();
    void            CreateBaseRDFsPythiaCommon();
    virtual void    CreateBaseRDFsPythiaExtra(){}
    virtual void    CreateBaseRDFsExtra() override{ return CreateBaseRDFsPythiaImpl(); }

    virtual void    FillHistograms() override;
    virtual void    FillHistogramsTruth(){}
    virtual void    FillHistogramsFullSim(){} // Placeholder only for now

public:
    explicit RDFBasedHistFillingPythia(bool isPrivate_input = false,
                                       double ecom_input = 5.36,
                                       bool with_data_resonance_cuts_input = false)
    : isPrivate(isPrivate_input),
      with_data_resonance_cuts(with_data_resonance_cuts_input),
      E_COM(ecom_input)
    {}

    ~RDFBasedHistFillingPythia(){}

    void SetIsPrivate(bool v){ isPrivate = v; }
    bool GetIsPrivate() const { return isPrivate; }

    void SetECom(double v){ E_COM = v; }
    double GetECom() const { return E_COM; }

    void SetWithDataResonanceCuts(bool v){ with_data_resonance_cuts = v; }
    bool GetWithDataResonanceCuts() const { return with_data_resonance_cuts; }
};

class RDFBasedHistFillingPythiaTruth : public virtual RDFBasedHistFillingPythia{
protected:
    std::vector<std::string>                vars1D_general;
    std::vector<std::string>                vars1D_general_over_dr;
    std::vector<std::string>                vars1D_general_over_pair_pt;
    std::vector<std::array<std::string, 2>> vars2D_general;
    std::vector<std::array<std::string, 2>> vars2D_general_over_dr;

    std::vector<std::string>                vars1D_flavor_origin;
    std::vector<std::string>                vars1D_flavor_origin_over_dr;
    std::vector<std::string>                vars1D_flavor_origin_over_pair_pt;
    std::vector<std::array<std::string, 2>> vars2D_flavor_origin;
    std::vector<std::array<std::string, 2>> vars2D_flavor_origin_over_dr;

    std::map<int, std::string> flavor_suffix_map;
    std::map<int, std::string> origin_suffix_map;

    virtual void    InitializePythiaExtra() override;
    virtual void    BuildHistBinningMapExtra() override;
    virtual void    CreateBaseRDFsPythiaExtra() override;
    virtual void    FillHistogramsTruth() override;
    virtual void    HistPostProcessExtra() override;
    virtual void    WriteOutputExtra() override;

    void FillHistogramsGeneral();
    void FillHistogramsFlavorBinned();
    void FillHistogramsOriginBinned();
    void FillHistogramsResonanceStudy();
    void FillHistogramsCrossxAndSpecialEta();
    void FillHistogramsSignalAcceptance();
    void ValidatePythiaTruthSchemaAndCoverage();
    void BuildAndStoreNearAwaySummedHistograms();
    void MarkNearAwayHistogramsForNominalExclusion();
    bool IsNearAwayDividedHistName(const std::string& hname) const;
    std::string NearAwayToSummedName(const std::string& hname) const;

public:
    explicit RDFBasedHistFillingPythiaTruth(bool isPrivate_input = false,
                                            double ecom_input = 5.36,
                                            bool with_data_resonance_cuts_input = false)
    : RDFBasedHistFillingPythia(isPrivate_input, ecom_input, with_data_resonance_cuts_input)
    {}

    ~RDFBasedHistFillingPythiaTruth(){}
};
