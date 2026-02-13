#include "RDFBasedHistFillingBaseClass.cxx"

class RDFBasedHistFillingPowheg : public virtual RDFBasedHistFillingBaseClass{
protected:
    int     run_year;
    bool    isRun3; // auto-determined from run_year
    bool    is_fullsim{false}; // if is fullsim, load reco quantities
    bool    is_fullsim_overlay{false}; // if is fullsim overlay, load centrality info
    bool    perform_truth{false};

    // ----- initialize powheg-specific settings -----
    virtual void    InitAnalysisSettingsHook() override{
        return InitializePowhegImpl();
    }

    virtual void    BuildFilterToVarListMapExtra() override{ 
        return BuildFilterToVarListMapPowhegImpl();
    }

    virtual void    FillHistograms() override;
    virtual void    FillHistogramsTruth();
    virtual void    FillHistogramsFullSim();
    virtual void    FillHistogramsFullSimTruthRecoCompr();
    virtual void    FillHistogramsFullSimRecoEffcies();


public:

};

class RDFBasedHistFillingPowhegTruth : public RDFBasedHistFillingPowheg{
protected:
    virtual void        SetIOPathsHook() override;

public:

};

class RDFBasedHistFillingPowhegFullsim : public RDFBasedHistFillingPowheg{
protected:
    virtual void        SetIOPathsHook() override;

public:

};