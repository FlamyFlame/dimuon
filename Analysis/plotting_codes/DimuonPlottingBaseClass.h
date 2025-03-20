#ifndef __DimuonPlottingBaseClass_h__
#define __DimuonPlottingBaseClass_h__

class DimuonPlottingBaseClass{
protected:

    static const int s_nSigns = 2;
    static const int s_nDphis = 2;

    static constexpr float s_NULL = -1000.;
    
    enum Sign{
        same_sign,
        opposite_sign
    };

    std::vector<std::string> single_b_analysis_key_observables = {}; // key observables for single-b analysis
    std::vector<std::string> gs_hf_analysis_key_observables = {}; // key observables for gluon-splitting analysis


// --------------------------------------------------------------------------------

    std::array<std::string,s_nSigns> signs = {"_sign1", "_sign2"};
    std::array<std::string,s_nSigns> signTitles = {"same sign", "opposite sign"};

    std::array<std::string,s_nDphis> dphis = {"_near", "_away"};
    std::array<std::string,s_nDphis> dphiTitles = {"#Delta #phi < #pi/2","#Delta #phi #geq #pi/2"};


public:
    std::string kin;
    std::string kin_title;


    virtual void Run() = 0;

};


#endif