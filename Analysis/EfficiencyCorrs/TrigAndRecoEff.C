#include "Bins.h"
#include "NoOverlayMC.C"
#include "TrigEffMu4NoL1.C"


vector<float>   *muon_eff_corr_medium_data=nullptr; // need to set branch address + read these branches
vector<float>   *muon_eff_SF_medium       =nullptr;
vector<float>   *muon_eff_corr_tight_data =nullptr;
vector<float>   *muon_eff_SF_tight        =nullptr;


int m_data_type    =Bins::DataType  ::PbPbAll; // see the enums in Bins.h - set to reflect the data type etc. I'm using
int m_quality_cut  =Bins::QualityCut::MEDIUM ;
const int m_eff_cor=Bins::EffCor    ::COR    ;

//in main code
/*
    const float pt1           = fabs(muon_pair_muon1_pt     ->at(i))/1000;
    const float pt2           = fabs(muon_pair_muon2_pt     ->at(i))/1000;
    const float eta1          =      muon_pair_muon1_eta    ->at(i);
    const float eta2          =      muon_pair_muon2_eta    ->at(i);
    const float phi1          =      muon_pair_muon1_phi    ->at(i);
    const float phi2          =      muon_pair_muon2_phi    ->at(i);
    const int   charge1       =     (muon_pair_muon1_pt     ->at(i)>0)?1:-1;
    const int   charge2       =     (muon_pair_muon2_pt     ->at(i)>0)?1:-1;

    const int   quality1      =      muon_pair_muon1_quality->at(i);
    const int   quality2      =      muon_pair_muon2_quality->at(i);
    const int   index1        =      muon_pair_muon1_index  ->at(i);
    const int   index2        =      muon_pair_muon2_index  ->at(i);

    bool trig_match_2mu4       =true/false;
    bool trig_match_mu4_mu4noL1=true/fasle;

    float eff1=1.0,eff2=1.0,eff_total=1.0;
    TrigAndRecoEff(eff1,eff2,eff_total, //Pass by reference to set the efficiencies
                  eta1,pt1 ,charge1  ,
                  eta2,pt2 ,charge2  ,
                  centrality,
                  index1,index2,
                  trig_match_2mu4,trig_match_mu4_mu4noL1);

    Should fill some MONITORING histograms for the efficiency corrections
    Plot eff_total in bins of (as function of) physical variables, e.g, centrality + ptavg
*/



void TrigAndRecoEff(float &eff1, float &eff2, float& eff_total, //calculated efficiency corrections
                    float eta1 , float pt1  , int charge1     , //information for muon-1: pt in GeV, charge=+-1
                    float eta2 , float pt2  , int charge2     , //information for muon-2: pt in GeV, charge=+-1
                    int centrality_percentile,                  //centrality percentile in 1% bin  
                    int m1, int m2,                             //index for the muons in the vector in the TTree entry
                    bool trig_match_2mu4, bool trig_match_mu4_mu4noL1 //The trigger whose eff we want, only one of these should be true
                    )
{
  
  
  static NoOverlayMC    *m_NoOverlayMC   =new NoOverlayMC   (true,m_quality_cut);//Object for Pb+Pb Reco efficiency
  static TrigEffMu4NoL1 *m_TrigEffMu4NoL1=new TrigEffMu4NoL1(2   ,m_quality_cut);//Object for Pb+Pb Trigger efficiency


  //-----------------------------------------------------
  //single muon reconstruction efficiencies

  //For Pb+Pb we use the HIJING Overlay
  if(m_data_type==Bins::DataType::PbPbAll  ||
     m_data_type==Bins::DataType::PbPb2015 ||
     m_data_type==Bins::DataType::PbPb2018)
  {
    //MC Scale factors for Reco Eff
    double mc_scale_factor1=1;
    double mc_scale_factor2=1;
    switch(m_quality_cut){
      case Bins::QualityCut::MEDIUM:
        //mc_scale_factor1=muon_eff_SF_medium->at(m1);
        //mc_scale_factor2=muon_eff_SF_medium->at(m2);
        break;
      case Bins::QualityCut::TIGHT:
        //mc_scale_factor1=muon_eff_SF_tight->at(m1);
        //mc_scale_factor2=muon_eff_SF_tight->at(m2);
        break;
      default: 
        Common::Exception(__LINE__,__FILE__);
    }
    eff1=m_NoOverlayMC->GetReconstructionEfficiency(eta1,pt1,charge1,centrality_percentile)*mc_scale_factor1;
    eff2=m_NoOverlayMC->GetReconstructionEfficiency(eta2,pt2,charge2,centrality_percentile)*mc_scale_factor2;
  }
  //For pp reconstruction efficiencies we use the ones from MCP 
  //http://atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData/MuonEfficiencyCorrections/180808_SummerUpdate/
  else if(m_data_type==Bins::DataType::pp2017){
    eff1=Bins::EffMCP2017(pt1,eta1,m_quality_cut);
    eff2=Bins::EffMCP2017(pt2,eta2,m_quality_cut);
  }
  else Common::Exception(__LINE__,__FILE__,"unimplemented data_type");
  //-----------------------------------------------------



  //-----------------------------------------------------
  //systematic uncertainties for reco eff 
  if(m_eff_cor==Bins::EffCor::CorrUp ||  //vary reco eff up
     m_eff_cor==Bins::EffCor::CorrDown){ //vary reco eff down

    float sign=(m_eff_cor==Bins::EffCor::CorrUp)? 1:-1;

    if(m_data_type==Bins::DataType::PbPbAll){
      eff1 +=eff1*sign*Bins::EffMCP2016(pt1,eta1,m_quality_cut,1);
      eff2 +=eff2*sign*Bins::EffMCP2016(pt2,eta2,m_quality_cut,1);
    }
    else if(m_data_type==Bins::DataType::pp2017){
      float syst1=Bins::EffMCP2017(pt1,eta1,m_quality_cut,1);
      float syst2=Bins::EffMCP2017(pt2,eta2,m_quality_cut,1);
      //TODO Require atleast 1% uncertainty : DONE
      syst1=(syst1>0.01)? syst1:0.01;
      syst2=(syst2>0.01)? syst2:0.01;
      eff1 +=eff1*sign*syst1;
      eff2 +=eff2*sign*syst2;
    }
    else Common::Exception(__LINE__,__FILE__,"unimplemented data_type");
  }
  //-----------------------------------------------------


  //-----------------------------------------------------
  //Based on comments from Cesare
  if(eff1>1.0) eff1=1.0;
  if(eff2>1.0) eff2=1.0;
  //-----------------------------------------------------
  


  //-----------------------------------------------------
  //single-muon trigger efficiencies
  float trig_eff1=1.0,trig_eff1a=0.0,trig_eff11a=0.0;
  float trig_eff2=1.0,trig_eff2a=0.0,trig_eff22a=0.0;
  if(m_data_type==Bins::DataType::PbPbAll  ||
     m_data_type==Bins::DataType::PbPb2015 ||
     m_data_type==Bins::DataType::PbPb2018
    ){
    trig_eff1 =m_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Mu4    ,eta1,pt1,charge1,centrality_percentile);
    trig_eff2 =m_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Mu4    ,eta2,pt2,charge2,centrality_percentile);
    trig_eff1a=m_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Mu4NoL1,eta1,pt1,charge1,centrality_percentile);
    trig_eff2a=m_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Mu4NoL1,eta2,pt2,charge2,centrality_percentile);
    //single muon trigger efficiencies for muon to fire both mu4 and mu4noL1
    trig_eff11a=m_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Both  ,eta1,pt1,charge1,centrality_percentile);
    trig_eff22a=m_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Both  ,eta2,pt2,charge2,centrality_percentile);

    if(m_eff_cor==Bins::EffCor::TrigCorrUp  || m_eff_cor==Bins::EffCor::TrigCorrDown){
      double sign=(m_eff_cor==Bins::EffCor::TrigCorrDown)? -1:1;
      trig_eff1   *= (1+sign*m_TrigEffMu4NoL1->CentMultiplierSystematic(TrigEffMu4NoL1::Mu4    ,pt1,centrality_percentile));
      trig_eff2   *= (1+sign*m_TrigEffMu4NoL1->CentMultiplierSystematic(TrigEffMu4NoL1::Mu4    ,pt2,centrality_percentile));
      trig_eff1a  *= (1+sign*m_TrigEffMu4NoL1->CentMultiplierSystematic(TrigEffMu4NoL1::Mu4NoL1,pt1,centrality_percentile));
      trig_eff2a  *= (1+sign*m_TrigEffMu4NoL1->CentMultiplierSystematic(TrigEffMu4NoL1::Mu4NoL1,pt2,centrality_percentile));
      trig_eff11a *= (1+sign*m_TrigEffMu4NoL1->CentMultiplierSystematic(TrigEffMu4NoL1::Both   ,pt1,centrality_percentile));
      trig_eff22a *= (1+sign*m_TrigEffMu4NoL1->CentMultiplierSystematic(TrigEffMu4NoL1::Both   ,pt2,centrality_percentile));
    }
  }
  //for pp efficiencies we use the ones from Qipeng
  else if(m_data_type==Bins::DataType::pp2017){
    static TH2D  *hist  =nullptr;
    static TH2D  *histSF=nullptr;
    static TH2D  *histSF2=nullptr;
    if(!hist){
      TFile *_file=gFile;
      TFile *file0  = new TFile("EffFiles/TriggerEfficiency_pp.root","read");
      if(file0->IsZombie()) Common::Exception(__LINE__,__FILE__);
      if(_file) _file->cd();//switch back to previous
      hist   =(TH2D*)file0->Get("trigeff2017_mu4_MC_StatError");
      histSF =(TH2D*)file0->Get("trigscalefactor2017_mu4_StatError");
      histSF2=(TH2D*)file0->Get("trigscalefactor2017_mu4_SystError");
    }
    //These are the Medium muon TrigEffs
    {
      float pt1_=(pt1<29)? pt1:29;//stay within the limits of the histogram
      float pt2_=(pt2<29)? pt2:29;

      int ipt1 =hist->GetYaxis()->FindBin(pt1_);
      int ipt2 =hist->GetYaxis()->FindBin(pt2_);
      int ieta1=hist->GetXaxis()->FindBin(eta1*charge1);
      int ieta2=hist->GetXaxis()->FindBin(eta2*charge2);

      int ipt1a =histSF->GetYaxis()->FindBin(pt1_);
      int ipt2a =histSF->GetYaxis()->FindBin(pt2_);
      int ieta1a=histSF->GetXaxis()->FindBin(fabs(eta1));
      int ieta2a=histSF->GetXaxis()->FindBin(fabs(eta2));

      trig_eff1 =hist  ->GetBinContent(ieta1 ,ipt1 );
      trig_eff2 =hist  ->GetBinContent(ieta2 ,ipt2 );

      double trig_syst=0;
      if(m_eff_cor==Bins::EffCor::TrigCorrUp  ) trig_syst=+1;
      if(m_eff_cor==Bins::EffCor::TrigCorrDown) trig_syst=-1;

      trig_eff1*=(histSF->GetBinContent(ieta1a,ipt1a) + trig_syst*histSF2->GetBinError(ieta1a,ipt1a));
      trig_eff2*=(histSF->GetBinContent(ieta2a,ipt2a) + trig_syst*histSF2->GetBinError(ieta2a,ipt2a));
    }
    //For Tight muons we dont have Trig-effs from Qipeng,
    //We use the ratio of Tight/Medium Effs that we have in Pb+Pb
    //  to scale the medium muon effs
    if(m_quality_cut==Bins::QualityCut::TIGHT){
      TFile *_file=gFile;
      static TrigEffMu4NoL1* l_TrigEffMu4NoL1=new TrigEffMu4NoL1(2   ,Bins::QualityCut::MEDIUM);//For medium muons
      if(_file) _file->cd();
      float n1=m_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Mu4    ,eta1,pt1,charge1,centrality_percentile);
      float n2=m_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Mu4    ,eta2,pt2,charge2,centrality_percentile);
      float d1=l_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Mu4    ,eta1,pt1,charge1,centrality_percentile);
      float d2=l_TrigEffMu4NoL1->GetTriggerEfficiency(TrigEffMu4NoL1::Mu4    ,eta2,pt2,charge2,centrality_percentile);
      trig_eff1*=n1/d1;
      trig_eff2*=n2/d2;
    }
  }
  else Common::Exception(__LINE__,__FILE__,"unimplemented data_type");
  //-----------------------------------------------------



  //-----------------------------------------------------
  //dimuon trigger efficiency
  //For pairs coming from the mu4_mu4_noL1 trigger
  float net_trig_eff=(trig_eff1*trig_eff2a) + 
                     (trig_eff2*trig_eff1a) - 
                     (trig_eff11a*trig_eff22a);
  //If mu4_mu4_noL1 did not fire then the pair comes from 2mu4, so eff=e1*e2
  if(!trig_match_mu4_mu4noL1) net_trig_eff=trig_eff1*trig_eff2;
  if(net_trig_eff>1) net_trig_eff=1.0;

  if(m_eff_cor==Bins::EffCor::NOCOR){ //Disable eff corrections in this case
    eff1=1.0;
    eff2=1.0;
    eff_total=1.0;
  }
  else{
    eff_total=1.0/(eff1*eff2*net_trig_eff);//pair efficiency
    eff1=1.0/(eff1*net_trig_eff);//eff for muon-1 only
    eff2=1.0/(eff2*net_trig_eff);//eff for muon-2 only
  }
  //-----------------------------------------------------
}






