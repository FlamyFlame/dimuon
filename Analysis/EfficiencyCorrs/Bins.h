#ifndef __BINS_H__
#define __BINS_H__
#include "common.C"

namespace Bins{
const double SmallNumber=1e-6;

enum PlotMode{
  IntNote  =0,
  ConfNote =1,
  EbeChecks=2,
  ConfNote2=3,
  PaperPRL =4,
  CheckV3  =5,
};
//bool Mode=IntNote;

int GetCentBin(int icentrality);
int Initialize_CentAdd();
int GetCentIndex(double cent_low,double cent_high);
std::string LabelCent(int icent,int option=0,float X=0.44, float Y=0.92, int SIZE=22, int color=1);


bool PassQualityCut(const int quality1,const int quality2, const int i_quality_cut);
bool PassPtCut(const float pt1, const float pt2, const int i_pt);
  
bool PassMomImbCut(const float momimb_rad, const int i_momimb_cut, const float i_momimb1, const float i_momimb2);
bool PassMomImbSigCut(const float momimbsig_rad, const int i_momimbsig_cut, const float i_momimbsig1, const float i_momimbsig2, const float rms);

std::string LabelQuality(const int i_quality_cut,int option=0,float X=0.44,float Y=0.92,int SIZE=22,int color=1);

std::string LabelCharge(int ich,int option=0,float X=0.44, float Y=0.92, int SIZE=22, int color=1);
void LabelATLAS (int l_data_type, float X,float Y,int SIZE,float dX1,float dY);
void LabelATLAS2(                 float X,float Y,int SIZE,float dX1,float dX2,float dY);


const double PbPbLumi     = 1.944*1e6;//lumi in mb^-1 //1943.895=(504.125+1439.77);//lumi in mub-1, computed by lumicalc using the HLT_2mu4 trigger for 2015 + 2018
const double PbPbLumi2015 = 0.504*1e6;
const double PbPbLumi2018 = 1.440*1e6;
const double PbPbCS       = 7.66*1e3;//CS in mb
const double ppLumi       = 258 *1e9;//HLT_2mu4 lumi from lumicalc in mb-1

//-----------------------------------------------------------------------------------------------------
enum{
  NCENT=18,
  NCENT_ORIGINAL=11,
  NCENT_ADD     = 7,
};
//Binning
int CENT_LO[NCENT]={0, 5,10,15,20,30,40,50,60,80, 90};
int CENT_HI[NCENT]={5,10,15,20,30,40,50,60,80,90,100};
int cent_add_lo[NCENT_ADD]  ={0};
int cent_add_up[NCENT_ADD]  ={0};

int GetCentBin(int icentrality)
{
  for(int i=0;i<NCENT_ORIGINAL;i++){
    if(icentrality<CENT_HI[i]) return i;
  }
  std::cout<<"Could not determine centrality bin for centrality="<<icentrality<<std::endl;
  Common::Exception(__LINE__, __FILE__);
  return -1;
}

int Initialize_CentAdd()
{
  //std::vector<std::pair<double,double>> new_cent_bins={{0,100},{0,10},{10,20},{0,20},{20,40},{40,60},{40,80}};
  std::vector<std::pair<double,double>> new_cent_bins={{0,100},{0,10},{10,20},{20,40},{40,60},{0,20},{80,100}};

  if(new_cent_bins.size()!=NCENT_ADD){
     std::cout<<"Initialize_CentAdd()::new_cent_bins.size()!=NCENT_ADD "
              <<new_cent_bins.size()<<"  "<<NCENT_ADD<<std::endl;
     throw std::exception();
  }

  int ibin=0;
  for(auto new_bin:new_cent_bins){
    int low=-1,high=-1;
    for(int icent=0;icent<NCENT_ORIGINAL;icent++){
      if(fabs(CENT_LO[icent]-new_bin.first )<0.001 ) low  =icent;
      if(fabs(CENT_HI[icent]-new_bin.second)<0.001 ) high =icent+1;
    }
    if(low==-1 || high==-1 || low>=high){
      std::cout<<"Initialize_CentAdd():: Problem adding new bin"<<std::endl;
      throw std::exception();
    }

    cent_add_lo[ibin]=low;
    cent_add_up[ibin]=high;
    CENT_LO[NCENT_ORIGINAL+ibin]=new_bin.first;
    CENT_HI[NCENT_ORIGINAL+ibin]=new_bin.second;
    std::cout<<"To get the centrality "
             <<CENT_LO[NCENT_ORIGINAL+ibin]<<"--"<<CENT_HI[NCENT_ORIGINAL+ibin]
             <<" we will add the bins ["<<low<<"--"<<high<<")"
             <<std::endl;
    ibin++;
  }
  std::cout<<"Initialize_CentAdd() finished"<<std::endl;
  return 1;
}

int GetCentIndex(double cent_low,double cent_high)
{
  for(int index=0;index<NCENT;index++){
    if(fabs(CENT_LO[index]-cent_low)<0.001 && fabs(CENT_HI[index]-cent_high)<0.001 ) return index;
  }
  std::cout<<"This Centrality doesnot exist "<<cent_low<<","<<cent_high<<std::endl;
  throw std::exception();
}


std::string LabelCent(int icent,int option,float X,float Y,int SIZE,int color)
{
  char label[600];
  sprintf(label,"%d-%d%%",CENT_LO[icent],CENT_HI[icent]);
  std::string ret=label;
  if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);
  return ret;
}
//-----------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------
enum PairSignCut{
  NCH =3,
    SAME_SIGN    =0,
    OPPOSITE_SIGN=1,
    COMBINED_SIGN=2,
};
std::string LabelCharge(int ich,int option,float X,float Y,int SIZE,int color){
   std::string ret;
   if     (ich==PairSignCut::SAME_SIGN    ) ret="Same-sign pairs";
   else if(ich==PairSignCut::OPPOSITE_SIGN) ret="Opp-sign pairs";
   else if(ich==PairSignCut::COMBINED_SIGN) ret="All pairs";
   else Common::Exception(__LINE__,__FILE__," unknown ich");

   if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);
   return ret;
}
//-----------------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------------
enum QualityCut{
  MEDIUM=0,
  TIGHT =1,
};
std::map<int,std::string> QualityCutLabel={{QualityCut::MEDIUM,"Medium-#it{#mu}"},
                                           {QualityCut::TIGHT ,"Tight-#it{#mu}" }};
bool PassQualityCut(const int quality1,const int quality2, const int i_quality_cut){
  if ((quality1&quality2&1)==0) return false;
  if (i_quality_cut==QualityCut::MEDIUM) {if((quality1&quality2&8 )!=8 ) return false; else return true;}
  if (i_quality_cut==QualityCut::TIGHT ) {if((quality1&quality2&16)!=16) return false; else return true;}
  return false;
}
std::string LabelQuality(const int i_quality_cut,int option,float X,float Y,int SIZE,int color){
  std::string ret;
  if      (i_quality_cut==QualityCut::MEDIUM) ret="Medium muons";
  else if (i_quality_cut==QualityCut::TIGHT ) ret="Tight muons";
  else Common::Exception(__LINE__,__FILE__," unknown quality_cut");
  if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);
  return ret;
}
//-----------------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------------
enum PTCut{
  PT4PLUS_4PLUS =0,
  PT4to4p5_4PLUS=1,
  PT4p5to5_4PLUS=2,
  PT5to5p5_4PLUS=3,
  PT5p5to6_4PLUS=4,

  PT4to5_4PLUS  =7,
  PT5to6_4PLUS  =8,
  PT6to8_4PLUS  =5,
  PT8to10_4PLUS =6,

  //Soumya
  PT5PLUS_5PLUS =20,
  PT4to5_4to5   =21,

  PTBar_4to5    =22,
  PTBar_GT5     =23,
  PTBar_GT6     =24,
  PTBar_GT7     =25,
  PTBar_GT6p5   =26,

  PTBar_5to6    =27,
  PTBar_4to6    =28,
  PTBar_6to8    =29,
  PTBar_8to10   =30,
  PTBar_6to10   =31,
  PTBar_GT4     =32,
  PTBar_GT4p5   =33,
};

bool PassPtCut(const float pt1, const float pt2, const int i_pt){
  float ptbar=(pt1+pt2)/2.0;
  if(pt1<4 || pt2<4) return false; //always require pt>4

  if     (i_pt==PTCut:: PT4PLUS_4PLUS){if(( 4    <pt1)      && (4 <pt2))                            return true; else return false;}
  else if(i_pt==PTCut::PT4to4p5_4PLUS){if((((4   <pt1)&&(pt1< 4.5)) && (4 <pt2)) || (((4   <pt2)&&(pt2< 4.5)) && (4 <pt1))) return true; else return false;}
  else if(i_pt==PTCut::PT4p5to5_4PLUS){if((((4.5 <pt1)&&(pt1< 5  )) && (4 <pt2)) || (((4.5 <pt2)&&(pt2< 5  )) && (4 <pt1))) return true; else return false;}
  else if(i_pt==PTCut::PT5to5p5_4PLUS){if((((5   <pt1)&&(pt1< 5.5)) && (4 <pt2)) || (((5   <pt2)&&(pt2< 5.5)) && (4 <pt1))) return true; else return false;}
  else if(i_pt==PTCut::PT5p5to6_4PLUS){if((((5.5 <pt1)&&(pt1< 6  )) && (4 <pt2)) || (((5.5 <pt2)&&(pt2< 6  )) && (4 <pt1))) return true; else return false;}

  else if(i_pt==PTCut::  PT6to8_4PLUS){if((((6   <pt1)&&(pt1< 8  )) && (4 <pt2)) || (((6   <pt2)&&(pt2< 8  )) && (4 <pt1))) return true; else return false;}
  else if(i_pt==PTCut:: PT8to10_4PLUS){if((((8   <pt1)&&(pt1< 10 )) && (4 <pt2)) || (((8   <pt2)&&(pt2< 10 )) && (4 <pt1))) return true; else return false;}
  else if(i_pt==PTCut::  PT4to5_4PLUS){if((((4   <pt1)&&(pt1<   5)) && (4 <pt2)) || (((4   <pt2)&&(pt2<   5)) && (4 <pt1))) return true; else return false;}
  else if(i_pt==PTCut::  PT5to6_4PLUS){if((((5   <pt1)&&(pt1<   6)) && (4 <pt2)) || (((5   <pt2)&&(pt2<   6)) && (4 <pt1))) return true; else return false;}

  else if(i_pt==PTCut:: PT5PLUS_5PLUS){if((  5    <pt1)             && (5 <pt2))                   return true; else return false;}
  else if(i_pt==PTCut:: PT4to5_4to5  ){if((  4<pt1)&&(pt1<5)        && (4 <pt2)&&(pt2<5))          return true; else return false;}

  else if(i_pt==PTCut::PTBar_4to5  ){if(ptbar>4 && ptbar<5  ) return true; else return false;}
  else if(i_pt==PTCut::PTBar_GT5   ){if(           ptbar>5  ) return true; else return false;}
  else if(i_pt==PTCut::PTBar_GT6   ){if(           ptbar>6  ) return true; else return false;}
  else if(i_pt==PTCut::PTBar_GT7   ){if(           ptbar>7  ) return true; else return false;}
  else if(i_pt==PTCut::PTBar_GT6p5 ){if(           ptbar>6.5) return true; else return false;}
  else if(i_pt==PTCut::PTBar_GT4   ){if(           ptbar>4  ) return true; else return false;}
  else if(i_pt==PTCut::PTBar_GT4p5 ){if(           ptbar>4.5) return true; else return false;}

  else if(i_pt==PTCut::PTBar_5to6  ){if(ptbar>5 && ptbar<6  ) return true; else return false;}
  else if(i_pt==PTCut::PTBar_4to6  ){if(ptbar>4 && ptbar<6  ) return true; else return false;}
  else if(i_pt==PTCut::PTBar_6to8  ){if(ptbar>6 && ptbar<8  ) return true; else return false;}
  else if(i_pt==PTCut::PTBar_8to10 ){if(ptbar>8 && ptbar<10 ) return true; else return false;}
  else if(i_pt==PTCut::PTBar_6to10 ){if(ptbar>6 && ptbar<10 ) return true; else return false;}
  else Common::Exception(__LINE__,__FILE__," unknown ipt");

  return false;
}
std::string LabelPt(int ipt,int option,float X=.2,float Y=.8,int SIZE=20,int color=1){
   std::string ret;
   if     (ipt==PTCut:: PT4PLUS_4PLUS) ret="#it{p}_{T}^{a,b}>4 GeV";
   else if(ipt==PTCut::PT4to4p5_4PLUS) ret="#it{p}_{T}^{a}>4 GeV,4<#it{p}_{T}^{b}<4.5 GeV";
   else if(ipt==PTCut::PT4p5to5_4PLUS) ret="#it{p}_{T}^{a}>4 GeV,4.5<#it{p}_{T}^{b}<5 GeV";
   else if(ipt==PTCut::PT5to5p5_4PLUS) ret="#it{p}_{T}^{a}>4 GeV,5<#it{p}_{T}^{b}<5.5 GeV";
   else if(ipt==PTCut::PT5p5to6_4PLUS) ret="#it{p}_{T}^{a}>4 GeV,5.5<#it{p}_{T}^{b}<6 GeV";
   else if(ipt==PTCut::  PT6to8_4PLUS) ret="#it{p}_{T}^{a}>4 GeV,6<#it{p}_{T}^{b}<8 GeV";
   else if(ipt==PTCut:: PT8to10_4PLUS) ret="#it{p}_{T}^{a}>4 GeV,8<#it{p}_{T}^{b}<10 GeV";
   else if(ipt==PTCut::  PT4to5_4PLUS) ret="#it{p}_{T}^{a}>4 GeV,4<#it{p}_{T}^{b}<5 GeV";
   else if(ipt==PTCut::  PT5to6_4PLUS) ret="#it{p}_{T}^{a}>4 GeV,5<#it{p}_{T}^{b}<6 GeV";

   else if(ipt==PTCut:: PT5PLUS_5PLUS) ret="#it{p}_{T}^{a,b}>5 GeV";
   else if(ipt==PTCut:: PT4to5_4to5  ) ret="4<#it{p}_{T}^{a,b}<5 GeV";

   else if(ipt==PTCut:: PTBar_4to5 ) ret="4<#bar{#it{p}}_{T}<5 GeV";
   else if(ipt==PTCut:: PTBar_GT5  ) ret="#bar{#it{p}}_{T}>5 GeV"  ;
   else if(ipt==PTCut:: PTBar_GT6  ) ret="#bar{#it{p}}_{T}>6 GeV"  ;
   else if(ipt==PTCut:: PTBar_GT7  ) ret="#bar{#it{p}}_{T}>7 GeV"  ;
   else if(ipt==PTCut:: PTBar_GT6p5) ret="#bar{#it{p}}_{T}>6.5 GeV";
   else if(ipt==PTCut:: PTBar_GT4  ) ret="#bar{#it{p}}_{T}>4 GeV";
   else if(ipt==PTCut:: PTBar_GT4p5) ret="#bar{#it{p}}_{T}>4.5 GeV";

   else if(ipt==PTCut:: PTBar_5to6 ) ret="5<#bar{#it{p}}_{T}<6 GeV";
   else if(ipt==PTCut:: PTBar_4to6 ) ret="4<#bar{#it{p}}_{T}<6 GeV";
   else if(ipt==PTCut:: PTBar_6to8 ) ret="6<#bar{#it{p}}_{T}<8 GeV";
   else if(ipt==PTCut:: PTBar_8to10) ret="8<#bar{#it{p}}_{T}<10 GeV";
   else if(ipt==PTCut:: PTBar_6to10) ret="6<#bar{#it{p}}_{T}<10 GeV";
   else Common::Exception(__LINE__,__FILE__," unknown ipt");

   if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);
   return ret;
}
//-----------------------------------------------------------------------------------------------------

bool TrueSignal(float recphi, float receta, float sigphi, float sigeta){
  if (recphi<-Common::PI) recphi=recphi+2*Common::PI;
  if (recphi> Common::PI) recphi=recphi-2*Common::PI;
  if (sigphi<-Common::PI) sigphi=sigphi+2*Common::PI;
  if (sigphi> Common::PI) sigphi=sigphi-2*Common::PI;

  if ((sqrt(pow((sigeta-receta),2)+pow((sigphi-recphi),2)))<0.02) return true; else return false;
  
  return false;
}
//-----------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------
/*
enum MomImbCut{
  RADCUT0p00= 0,
  RADCUT0p01= 1,
  RADCUT0p02= 2,
  RADCUT0p03= 3,
  RADCUT0p04= 4,
  RADCUT0p05= 5,
  RADCUT0p06= 6,
  RADCUT0p07= 7,
  RADCUT0p08= 8,
  RADCUT0p09= 9,
  RADCUT0p10=10,
  RADCUT0p11=11,
  RADCUT0p12=12,
  RADCUT0p13=13,
  RADCUT0p14=14,
  RADCUT0p15=15,
  RADCUT5p00=16,
  SQRCUTn0p3to0p0  =17,
  SQRCUT0p0to0p15  =18,
  SQRCUT0p15to0p3  =19,
  SQRCUT0p3plus    =20,
  SQRTCUTn0p3to0p15=21,
  SQRTCUTn0p3to0p0and0p15plus=22
  
};

bool PassMomImbCut(const float momimb_rad, const int i_momimb_cut, const float momimb1, const float momimb2){
  if     (i_momimb_cut==MomImbCut::RADCUT0p00){if(momimb_rad<0.00) return true; else return false;}//No signal
  else if(i_momimb_cut==MomImbCut::RADCUT0p01){if(momimb_rad<0.01) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p02){if(momimb_rad<0.02) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p03){if(momimb_rad<0.03) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p04){if(momimb_rad<0.04) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p05){if(momimb_rad<0.05) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p06){if(momimb_rad<0.06) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p07){if(momimb_rad<0.07) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p08){if(momimb_rad<0.08) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p09){if(momimb_rad<0.09) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p10){if(momimb_rad<0.10) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p11){if(momimb_rad<0.11) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p12){if(momimb_rad<0.12) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p13){if(momimb_rad<0.13) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p14){if(momimb_rad<0.14) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT0p15){if(momimb_rad<0.15) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::RADCUT5p00){if(momimb_rad<5.00) return true; else return false;}//All signal
  else if(i_momimb_cut==MomImbCut::  SQRCUTn0p3to0p0){if(((-0.3<momimb1)&&(momimb1< 0.0))&&((-0.3<momimb2)&&(momimb2< 0.0))) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::  SQRCUT0p0to0p15){if((( 0.0<momimb1)&&(momimb1<0.15))&&(( 0.0<momimb2)&&(momimb2<0.15))) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::  SQRCUT0p15to0p3){if(((0.15<momimb1)&&(momimb1< 0.3))&&((0.15<momimb2)&&(momimb2< 0.3))) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::    SQRCUT0p3plus){if((  0.3<momimb1)&&(0.3<momimb2))                                     return true; else return false;}
  else if(i_momimb_cut==MomImbCut::SQRTCUTn0p3to0p15){if(((-0.3<momimb1)&&(momimb1<0.15))&&((-0.3<momimb2)&&(momimb2<0.15))) return true; else return false;}
  else if(i_momimb_cut==MomImbCut::SQRTCUTn0p3to0p0and0p15plus){
    if((((-0.3<momimb1)&&(momimb1<0.0))&&(0.15<momimb2))||
       (((-0.3<momimb2)&&(momimb2<0.0))&&(0.15<momimb1))){
      return true;
    }
    else return false;}

  return false;
}

std::string LabelMomImbCut(int imomimb_cut,int option,float X,float Y,int SIZE,int color){
  std::string ret;
  if     (imomimb_cut==MomImbCut::SQRCUTn0p3to0p0)   ret="-0.3 <#Deltap/p<0.0 ";
  else if(imomimb_cut==MomImbCut::SQRCUT0p0to0p15)   ret=" 0.0 <#Deltap/p<0.15";
  else if(imomimb_cut==MomImbCut::SQRCUT0p15to0p3)   ret=" 0.15<#Deltap/p<0.3 ";
  else if(imomimb_cut==MomImbCut::SQRCUT0p3plus)     ret=" 0.3 <#Deltap/p"     ;
  else if(imomimb_cut==MomImbCut::SQRTCUTn0p3to0p15) ret="-0.3 <#Deltap/p<0.15";
  else if(imomimb_cut==MomImbCut::SQRTCUTn0p3to0p0and0p15plus) ret="-0.3<#Deltap_{a}/p_{a}<0.0, 0.15<#Deltap_{b}/p_{b} [GeV]";
  else if(imomimb_cut==MomImbCut::RADCUT0p12                 ) ret="R(#Deltap/p)<0.12";
  else Common::Exception(__LINE__,__FILE__," unknown momimbcut");

  if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);  return ret;
}
*/
//-----------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------
enum MomImbSigCut{
  RADSIGCUT0to2p5RMS   = 0,
  RADSIGCUT2p5to4p5RMS = 1,
  RADSIGCUTI4p5t7RMS   = 2,
  RADSIGCUTA3p5RMS     = 3,
  RADSIGCUTA7p0RMS     = 4,
  SQRSIGCUT1p5RMS      = 5,
  SQRSIGCUT2p0RMS      = 6,
  SQRSIGCUT2p5RMS      = 7,
  SQRSIGCUT3p0RMS      = 8,
  SQRSIGCUT3p5RMS      = 9,
  NOSIGCUT             =10,
  RADSIGCUT0to1p5RMS   =11,
  RADSIGCUT0to1p0RMS   =12,
};

bool PassMomImbSigCut(const float momimbsig_rad, 
                      const int i_momimbsig_cut, 
                      const float i_momimbsig1, 
                      const float i_momimbsig2, 
                      const float rms){
  if     (i_momimbsig_cut==MomImbSigCut::  RADSIGCUT0to2p5RMS){if(momimbsig_rad<2.5*rms) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut:: RADSIGCUT2p5to4p5RMS){if((momimbsig_rad>2.5*rms)&&(momimbsig_rad<4.5*rms)) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::RADSIGCUTI4p5t7RMS){if((momimbsig_rad>4.5*rms)&&(momimbsig_rad<7.0*rms)) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::  RADSIGCUTA3p5RMS){if(momimbsig_rad>2.5*rms) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::  RADSIGCUTA7p0RMS){if(momimbsig_rad>4.5*rms) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::   SQRSIGCUT1p5RMS){if(((-1.5*rms<i_momimbsig1)&&(i_momimbsig1<1.5*rms))&&((-1.5*rms<i_momimbsig2)&&(i_momimbsig2<1.5*rms))) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::   SQRSIGCUT2p0RMS){if(((-2.0*rms<i_momimbsig1)&&(i_momimbsig1<2.0*rms))&&((-2.0*rms<i_momimbsig2)&&(i_momimbsig2<2.0*rms))) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::   SQRSIGCUT2p5RMS){if(((-2.5*rms<i_momimbsig1)&&(i_momimbsig1<2.5*rms))&&((-2.5*rms<i_momimbsig2)&&(i_momimbsig2<2.5*rms))) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::   SQRSIGCUT3p0RMS){if(((-3.0*rms<i_momimbsig1)&&(i_momimbsig1<3.0*rms))&&((-3.0*rms<i_momimbsig2)&&(i_momimbsig2<3.0*rms))) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::   SQRSIGCUT3p5RMS){if(((-5.0*rms<i_momimbsig1)&&(i_momimbsig1<2.5*rms))&&((-5.0*rms<i_momimbsig2)&&(i_momimbsig2<2.5*rms))) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::          NOSIGCUT) return true;
  else if(i_momimbsig_cut==MomImbSigCut::  RADSIGCUT0to1p5RMS){if(momimbsig_rad<1.5*rms) return true; else return false;}
  else if(i_momimbsig_cut==MomImbSigCut::  RADSIGCUT0to1p0RMS){if(momimbsig_rad<1.0*rms) return true; else return false;}
  else Common::Exception(__LINE__,__FILE__," unknown momimbsigcut");
  

  return false;
}

std::string LabelMomImbSigCut(int imomimbsig_cut,int option,float X,float Y,int SIZE,int color){
  std::string ret;
  if     (imomimbsig_cut==MomImbSigCut::RADSIGCUT0to2p5RMS  )   ret="#Deltap/p pair sig<2.5";
  else if(imomimbsig_cut==MomImbSigCut::RADSIGCUT2p5to4p5RMS)   ret="2.5<#Deltap/p pair sig<4.5";
  else if(imomimbsig_cut==MomImbSigCut::RADSIGCUTI4p5t7RMS  )   ret="4.5<#Deltap/p pair sig<7.0";
  else if(imomimbsig_cut==MomImbSigCut::RADSIGCUTA3p5RMS    )   ret="#Deltap/p pair sig>3.5";
  else if(imomimbsig_cut==MomImbSigCut::RADSIGCUTA7p0RMS    )   ret="#Deltap/p pair sig>7.0";
  else if(imomimbsig_cut==MomImbSigCut::NOSIGCUT            )   ret="No #Deltap/p pair sig cuts";
  else if(imomimbsig_cut==MomImbSigCut::RADSIGCUT0to1p5RMS  )   ret="#Deltap/p pair sig<1.5";
  else if(imomimbsig_cut==MomImbSigCut::RADSIGCUT0to1p0RMS  )   ret="#Deltap/p pair sig<1.0";
  else Common::Exception(__LINE__,__FILE__," unknown momimbsigcut");

  if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);
  return ret;
}
//-----------------------------------------------------------------------------------------------------  



//-----------------------------------------------------------------------------------------------------
enum DataType{
  PbPbAll  =0,
  PbPb2015 =1,
  PbPb2018 =2,
  pp2017   =3,
  HIJINGsig=4,
  ppBG     =5,

  //Only used for Trigger Effs
  pp2017_13TeV    =10,
  pp2017_13TeV_MB =11,
  pp2017_MB       =12,
  pp2017_13TeV_HMT=13,
  pp2017_HMT      =14,
};
//-----------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------  
enum EffCor{
  NOCOR    =0,
  COR      =1,
  CorrUp   =2,//Systematic variation upward
  CorrDown =3,//Systematic variation downward
  InbuiltPP=4,//efficiencies obtained from Qipeng

  TrigCorrUp   =5,//Systematic variation upward for TrigEff
  TrigCorrDown =6,//Systematic variation downward for TrigEff
};
std::string LabelEff(int i_eff_cor,int option,float X,float Y,int SIZE,int color){
   std::string ret;
   if     (i_eff_cor==EffCor::NOCOR    ) ret="No Eff Corr";
   else if(i_eff_cor==EffCor::COR      ) ret="Eff Corrected";
   else if(i_eff_cor==EffCor::CorrUp   ) ret="Eff Corrected (high)";
   else if(i_eff_cor==EffCor::CorrDown ) ret="Eff Corrected (Low)";
   else if(i_eff_cor==EffCor::InbuiltPP) ret="Eff Corrected (Tool)";
   else if(i_eff_cor==EffCor::TrigCorrUp   ) ret="TrigEff (high)";
   else if(i_eff_cor==EffCor::TrigCorrDown ) ret="TrigEff (Low)";
   else Common::Exception(__LINE__,__FILE__," unknown itrig");

   if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);
   return ret;
}


//https://atlas-groupdata.web.cern.ch/atlas-groupdata/MuonEfficiencyCorrections/210222_Precision_r21/
float EffMCP2017(float pt, float eta, int quality, int get_syst_error=0){
  static TH2D  *Eff_2017         [2]= {nullptr};
  static TH2D  *Eff_SYS_1UP_2017 [2]= {nullptr};
  static TH2D  *Eff_SYS_1DN_2017 [2]= {nullptr};
  if(!Eff_2017[quality]){
    TFile *oldFile=gFile;

    std::string                          eff_file_name="EffFiles/210222_Precision_r21/Reco_Medium_JPsi.root";
    if(quality==Bins::QualityCut::TIGHT) eff_file_name="EffFiles/210222_Precision_r21/Reco_Tight_JPsi.root" ;

    TFile *file0 = new TFile(eff_file_name.c_str(),"read");
    if(file0->IsZombie()) Common::Exception(__LINE__,__FILE__);

    if(oldFile) oldFile->cd();

    Eff_2017        [quality]=(TH2D*) file0->Get("Eff_2017"        );
    Eff_SYS_1UP_2017[quality]=(TH2D*) file0->Get("Eff_SYS_1UP_2017");
    Eff_SYS_1DN_2017[quality]=(TH2D*) file0->Get("Eff_SYS_1DN_2017");
  }

  float pt_ =(pt<19)? pt:19;//to stay within the limits of the histogram

  //return syst error
  if(get_syst_error==1){
    int  ipt   =Eff_SYS_1UP_2017[quality]->GetXaxis()->FindBin(pt_);
    int  ieta  =Eff_SYS_1UP_2017[quality]->GetYaxis()->FindBin(eta);
    float err1 =fabs(Eff_SYS_1UP_2017[quality]->GetBinContent(ipt,ieta));
    float err2 =fabs(Eff_SYS_1DN_2017[quality]->GetBinContent(ipt,ieta));

    return (err1>err2)? err1:err2;
  }

  
  int  ipt  =Eff_2017[quality]->GetXaxis()->FindBin(pt_);
  int  ieta =Eff_2017[quality]->GetYaxis()->FindBin(eta);
  float eff =Eff_2017[quality]->GetBinContent(ipt,ieta);

  return eff;
}


//http://atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData/MuonEfficiencyCorrections/180808_SummerUpdate/
float EffMCP2016(float pt, float eta, int quality, int get_syst_error=0){
   static TH2F* Eff_2016    [2]={nullptr};
   static TH2F* Eff_sys_2016[2]={nullptr};
   if(!Eff_2016[quality]){
     TFile *oldFile=gFile;

     std::string                          eff_file_name="EffFiles/180808_SummerUpdate/Reco_Medium_JPsi.root";
     if(quality==Bins::QualityCut::TIGHT) eff_file_name="EffFiles/180808_SummerUpdate/Reco_Tight_JPsi.root" ;
     
     TFile* InFile=new TFile(eff_file_name.c_str(),"read");
     if(InFile->IsZombie()) Common::Exception(__LINE__,__FILE__);

     if(oldFile) oldFile->cd();

     Eff_2016    [quality]=(TH2F*) InFile->Get("Eff_2016");
     Eff_sys_2016[quality]=(TH2F*) InFile->Get("Eff_sys_2016");
   }

   float pt_=(pt<19)? pt:19;//to stay within the limits of the histogram

   //return syst error
   if(get_syst_error==1){
     int ieta=Eff_sys_2016[quality]->GetXaxis()->FindFixBin(eta);
     int ipt =Eff_sys_2016[quality]->GetYaxis()->FindFixBin(pt_);
     return   fabs(Eff_sys_2016[quality]->GetBinContent(ieta,ipt));
   }

   //Note this histograms x-axis is eta (not q*eta)
   int ieta=Eff_2016[quality]->GetXaxis()->FindFixBin(eta);
   int ipt =Eff_2016[quality]->GetYaxis()->FindFixBin(pt_);
   return   fabs(Eff_2016[quality]->GetBinContent(ieta,ipt));
}
//-----------------------------------------------------------------------------------------------------  


//-----------------------------------------------------------------------------------------------------  
enum Trigger{
  AllTrigs=0,
  Only2Mu4=1,
};
std::string LabelTrigger(int itrig,int option,float X,float Y,int SIZE,int color){
   std::string ret;
   if     (itrig==Trigger::AllTrigs ) ret="All Trigs";
   else if(itrig==Trigger::Only2Mu4 ) ret="2mu4";
   else Common::Exception(__LINE__,__FILE__," unknown itrig");

   if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);
   return ret;
}
//-----------------------------------------------------------------------------------------------------  

//-----------------------------------------------------------------------------------------------------  
enum DetaCut{
  Deta0p8 =0,
  Deta0p9 =1,
  Deta1p0 =2,
  Deta1p2 =3,
  Deta1p4 =4,
  Deta1p6 =5,
  Deta1p8 =6,
  Deta2p0 =7,
  Deta2p2 =8,
};
std::map<int,float> DetaCutVals={
  {Deta0p8,0.8},
  {Deta0p9,0.9},
  {Deta1p0,1.0},
  {Deta1p2,1.2},
  {Deta1p4,1.4},
  {Deta1p6,1.6},
  {Deta1p8,1.8},
  {Deta2p0,2.0},
  {Deta2p2,2.2},
  };
bool PassDetaCut(int icut,float deta_val){
  if(DetaCutVals.count(icut)!=1) Common::Exception(__LINE__,__FILE__," unknown DetaCut");
  return (fabs(deta_val)>DetaCutVals[icut]);
}
std::string LabelDetaCut(int ideta,int option,float X=0.2,float Y=0.88,int SIZE=20,int color=1){
   std::string ret;
   if     (ideta==DetaCut::Deta0p8 ) ret="|#Delta#eta|>0.8";
   else if(ideta==DetaCut::Deta0p9 ) ret="|#Delta#eta|>0.9";
   else if(ideta==DetaCut::Deta1p0 ) ret="|#Delta#eta|>1.0";
   else if(ideta==DetaCut::Deta1p2 ) ret="|#Delta#eta|>1.2";
   else if(ideta==DetaCut::Deta1p4 ) ret="|#Delta#eta|>1.4";
   else if(ideta==DetaCut::Deta1p6 ) ret="|#Delta#eta|>1.6";
   else if(ideta==DetaCut::Deta1p8 ) ret="|#Delta#eta|>1.8";
   else if(ideta==DetaCut::Deta2p0 ) ret="|#Delta#eta|>2.0";
   else if(ideta==DetaCut::Deta2p2 ) ret="|#Delta#eta|>2.2";

   else Common::Exception(__LINE__,__FILE__," unknown DetaCut");

   if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);
   return ret;
}
//-----------------------------------------------------------------------------------------------------  

//-----------------------------------------------------------------------------------------------------
enum PDF{
  NTYPE=4,
    nNNPDF20 =0,
    nCTEQ15  =1,
    NNPDF30  =2,
    NoDYCorr =3,//no corrections applied, returns empty histograms 
};
map<int,string> LabelPDF={
  {PDF::NoDYCorr,"NoDYCorr"},
  {PDF::nNNPDF20,"nNNPDF20"},
  {PDF::nCTEQ15 ,"nCTEQ15" },
  {PDF::NNPDF30 ,"NNPDF30" },
};
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------  
enum OtherFlags{
  None               = 0,
  
  NoMinvCuts         = 1,
  NoMinvEff          = 2,
  MinEffFromEventMix = 4,
};
//-----------------------------------------------------------------------------------------------------  

//-----------------------------------------------------------------------------------------------------
struct Cuts{
  int m_quality_cut   =QualityCut  ::MEDIUM       ;
  int m_pt_cut        =PTCut       ::PT4PLUS_4PLUS;
  int m_momimbsig_cut =MomImbSigCut::NOSIGCUT     ;
  int m_data_type     =DataType    ::PbPbAll      ;
  int m_eff_cor       =EffCor      ::COR          ;
  int m_trig          =Trigger     ::AllTrigs     ;
  int m_deta          =DetaCut     ::Deta0p8      ;
  int m_pdf           =PDF         ::nNNPDF20     ;
  int m_other_flags   =OtherFlags  ::None         ;
};
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
const std::vector<std::pair<float,float>> MinvCuts={{0.6,1.05},{2.9,3.3},{3.6,3.8},{9.2,10.4},{70.0,110.0}};
const double AcopCut =0.01;
const double AsymmCut=0.08;
//-----------------------------------------------------------------------------------------------------





int LabelATLASColor=1;
void LabelATLAS(int l_data_type,float X,float Y,int SIZE,float dX1,float dY){
    Common::myText2(X     , Y    ,LabelATLASColor, "ATLAS"   , SIZE, 73);
    Common::myText2(X+dX1 , Y    ,LabelATLASColor, Common::Internal, SIZE, 43);
    std::string label= "Pb+Pb 5.02 TeV, 1.94 nb^{-1}";
    if(l_data_type==DataType::PbPb2018) label="Pb+Pb Data 2018 #sqrt{#it{s}_{NN}}=5.02 TeV, 1.44 nb^{-1}";
    if(l_data_type==DataType::PbPb2015) label="Pb+Pb Data 2015 #sqrt{#it{s}_{NN}}=5.02 TeV, 0.44 nb^{-1}";
    if(l_data_type==DataType::pp2017  ) label="#it{pp} 5.02 TeV, 0.26 fb^{-1}";
    Common::myText2(X     , Y-dY ,LabelATLASColor,label, SIZE, 43);
}


void LabelATLAS2(float X,float Y,int SIZE,float dX1,float dX2,float dY){
    Common::myText2(X     , Y    ,LabelATLASColor, "ATLAS"   , SIZE, 73);
    Common::myText2(X+dX1 , Y    ,LabelATLASColor, Common::Internal, SIZE, 43);
    Common::myText2(X     , Y-dY ,LabelATLASColor, "Pb+Pb", SIZE, 43);
    Common::myText2(X+dX2 , Y-dY ,LabelATLASColor, "#sqrt{#it{s}_{NN}}=5.02 TeV", SIZE, 43);
}
  


//-----------------------------------------------------------------------------------------------------
  enum MCType{
    NumMCType=4,
      PYTHIA   =0,
      POWHEG_BB=1,
      POWHEG_CC=2,
      POWHEG   =3,
  };
  map<int,string> LabelMC={
    {MCType::PYTHIA   ,"PYTHIA"   },
    {MCType::POWHEG_BB,"POWHEG_BB"},
    {MCType::POWHEG_CC,"POWHEG_CC"},
    {MCType::POWHEG   ,"POWHEG"   },
  };
  map<int,string> LabelMCLegend={
    {MCType::PYTHIA   ,"PYTHIA"   },
    {MCType::POWHEG_BB,"POWHEG b#bar{b}"},
    {MCType::POWHEG_CC,"POWHEG c#bar{c}"},
    {MCType::POWHEG   ,"POWHEG b#bar{b}+c#bar{c}"},
  };
  std::string LabelMC2 (int imc,int option=0,float X=0.44,float Y=0.92,int SIZE=22,int color=1){
    std::string ret=LabelMC[imc];
    if(option==0) Common::myText2(X,Y,color,ret,SIZE,43);
    return ret;
  }
  void LabelMC3(int l_mc_type,float X,float Y,int SIZE,float dX1,float dY){
      Common::myText2(X     , Y    ,LabelATLASColor, "ATLAS"   , SIZE, 73);
      Common::myText2(X+dX1 , Y    ,LabelATLASColor, Common::Internal+" Simulation", SIZE, 43);
      std::string label;
      if     (l_mc_type==MCType::PYTHIA   ) label="PYTHIA #sqrt{#it{s}}=5.02 TeV";
      else if(l_mc_type==MCType::POWHEG_BB) label="POWHEG #sqrt{#it{s}}=5.02 TeV, b#bar{b}#rightarrow#it{#mu#mu}";
      else if(l_mc_type==MCType::POWHEG_CC) label="POWHEG #sqrt{#it{s}}=5.02 TeV, c#bar{c}#rightarrow#it{#mu#mu}";
      else if(l_mc_type==MCType::POWHEG   ) label="POWHEG #sqrt{#it{s}}=5.02 TeV, b#bar{b}#rightarrow#it{#mu#mu} + c#bar{c}#rightarrow#it{#mu#mu}";
      else Common::Exception(__LINE__,__FILE__," Unknown MC");
      Common::myText2(X     , Y-dY ,LabelATLASColor,label, SIZE, 43);
  }
//-----------------------------------------------------------------------------------------------------


int is_initialized=Initialize_CentAdd();

std::string SimpleFileName(const int l_quality_cut){
  std::string name="";
  if      (l_quality_cut==Bins::QualityCut::MEDIUM) name+="";
  else if (l_quality_cut==Bins::QualityCut::TIGHT ) name+="_Tight" ;
  else Common::Exception(__LINE__,__FILE__," Unknown quality_cut");

  return name;
}

std::string FileName(const int l_quality_cut  , 
                     const int l_pt_cut       , 
                     const int l_momimbsig_cut, 
                     const int l_data_type    , 
                     const int l_eff_cor      ,
                     const int l_trig         ,
                     const int l_deta         ,
                     const int l_pdf          ,
                     const int l_other_flag    
                    )
{
   cout<<"qual="<<l_quality_cut<<" pt="<<l_pt_cut<<" momimb="<<l_momimbsig_cut<<" data="<<l_data_type<<" eff="<<l_eff_cor<<" trig="<<l_trig<<" deta="<<l_deta<<" pdf="<<l_pdf<<" other="<<l_other_flag<<endl;
   std::string name="";
   if      (l_quality_cut==Bins::QualityCut::MEDIUM) name+="";
   else if (l_quality_cut==Bins::QualityCut::TIGHT ) name+="_Tight" ;
   else Common::Exception(__LINE__,__FILE__," Unknown quality_cut");

   if     (l_pt_cut==Bins::PTCut:: PT4PLUS_4PLUS) name+="_pT_GT4";
   else if(l_pt_cut==Bins::PTCut::PT4to4p5_4PLUS) name+="_pT4to4p5_4plus";
   else if(l_pt_cut==Bins::PTCut::PT4p5to5_4PLUS) name+="_pT4p5to5_4plus";
   else if(l_pt_cut==Bins::PTCut::PT5to5p5_4PLUS) name+="_pT5to5p5_4plus";
   else if(l_pt_cut==Bins::PTCut::PT5p5to6_4PLUS) name+="_pT5p5to6_4plus";
   else if(l_pt_cut==Bins::PTCut::  PT6to8_4PLUS) name+=  "_pT6to8_4plus";
   else if(l_pt_cut==Bins::PTCut:: PT8to10_4PLUS) name+= "_pT8to10_4plus";
   else if(l_pt_cut==Bins::PTCut::  PT4to5_4PLUS) name+=  "_pT4to5_4plus";
   else if(l_pt_cut==Bins::PTCut::  PT5to6_4PLUS) name+=  "_pT5to6_4plus";

   else if(l_pt_cut==Bins::PTCut:: PT5PLUS_5PLUS) name+="_pT5plus_5plus";
   else if(l_pt_cut==Bins::PTCut:: PT4to5_4to5  ) name+="_pT4to5_4to5"  ;

   else if(l_pt_cut==Bins::PTCut:: PTBar_4to5  ) name+="_pTBar_4to5"  ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_GT5   ) name+="_pTBar_GT5"   ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_GT6   ) name+="_pTBar_GT6"   ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_GT7   ) name+="_pTBar_GT7"   ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_GT6p5 ) name+="_pTBar_GT6p5" ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_GT4   ) name+="_pTBar_GT4"   ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_GT4p5 ) name+="_pTBar_GT4p5" ;

   else if(l_pt_cut==Bins::PTCut:: PTBar_5to6  ) name+="_pTBar_5to6"  ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_4to6  ) name+="_pTBar_4to6"  ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_6to8  ) name+="_pTBar_6to8"  ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_8to10 ) name+="_pTBar_8to10" ;
   else if(l_pt_cut==Bins::PTCut:: PTBar_6to10 ) name+="_pTBar_6to10" ;

   else Common::Exception(__LINE__, __FILE__, " Unknown pt_cut "+std::to_string(l_pt_cut));
   /*
   if     (l_momimb_cut==Bins::MomImbCut::RADCUT0p00) name+="_RadCut0p00";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p01) name+="_RadCut0p01";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p02) name+="_RadCut0p02";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p03) name+="_RadCut0p03";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p04) name+="_RadCut0p04";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p05) name+="_RadCut0p05";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p06) name+="_RadCut0p06";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p07) name+="_RadCut0p07";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p08) name+="_RadCut0p08";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p09) name+="_RadCut0p09";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p10) name+="_RadCut0p10";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p11) name+="_RadCut0p11";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p12) name+="";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p13) name+="_RadCut0p13";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p14) name+="_RadCut0p14";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT0p15) name+="_RadCut0p15";
   else if(l_momimb_cut==Bins::MomImbCut::RADCUT5p00) name+="_RadCut5p00";
   else if(l_momimb_cut==Bins::MomImbCut::SQRCUTn0p3to0p0)   name+="_SQRCUTn0p3to0p0";
   else if(l_momimb_cut==Bins::MomImbCut::SQRCUT0p0to0p15)   name+="_SQRCUT0p0to0p15";
   else if(l_momimb_cut==Bins::MomImbCut::SQRCUT0p15to0p3)   name+="_SQRCUT0p15to0p3";
   else if(l_momimb_cut==Bins::MomImbCut::SQRCUT0p3plus)     name+="_SQRCUT0p3plus";
   else if(l_momimb_cut==Bins::MomImbCut::SQRTCUTn0p3to0p15) name+="_SQRTCUTn0p3to0p15";
   else if(l_momimb_cut==Bins::MomImbCut::SQRTCUTn0p3to0p0and0p15plus) name+="_SQRTCUTn0p3to0p0and0p15plus";
   else Common::Exception(__LINE__, __FILE__, " Unknown momimb_cut");
   */

   if     (l_momimbsig_cut==Bins::MomImbSigCut::RADSIGCUT0to2p5RMS  ) name+="_RADSIGCUT0to2p5RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::RADSIGCUT2p5to4p5RMS) name+="_RADSIGCUT2p5to4p5RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::RADSIGCUTI4p5t7RMS  ) name+="_RADSIGCUTI4p5t7RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::RADSIGCUTA3p5RMS    ) name+="_RADSIGCUTA3p5RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::RADSIGCUTA7p0RMS    ) name+="_RADSIGCUTA7p0RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::SQRSIGCUT1p5RMS     ) name+="_SqrSigCut1p5RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::SQRSIGCUT2p0RMS     ) name+="_SqrSigCut2p0RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::SQRSIGCUT2p5RMS     ) name+="_SqrSigCut2p5RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::SQRSIGCUT3p0RMS     ) name+="_SqrSigCut3p0RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::SQRSIGCUT3p5RMS     ) name+="_SqrSigCut3p5RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::NOSIGCUT            ) name+="";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::RADSIGCUT0to1p5RMS  ) name+="_RADSIGCUT0to1p5RMS";
   else if(l_momimbsig_cut==Bins::MomImbSigCut::RADSIGCUT0to1p0RMS  ) name+="_RADSIGCUT0to1p0RMS";
   else Common::Exception(__LINE__, __FILE__, " Unknown momimbsig_cut");
   
   if     (l_data_type==Bins::DataType::PbPbAll ) name+="";
   else if(l_data_type==Bins::DataType::PbPb2018) name+="_PbPb2018";
   else if(l_data_type==Bins::DataType::PbPb2015) name+="_PbPb2015";
   else if(l_data_type==Bins::DataType::pp2017  ) name+="_pp2017";
   else Common::Exception(__LINE__, __FILE__, " Unknown data_type");

   if      (l_eff_cor==Bins::EffCor::NOCOR    ) name+="";
   else if (l_eff_cor==Bins::EffCor::COR      ) name+="_EffCor";
   else if (l_eff_cor==Bins::EffCor::CorrUp   ) name+="_EffCorUp";
   else if (l_eff_cor==Bins::EffCor::CorrDown ) name+="_EffCorDown";
   else if (l_eff_cor==Bins::EffCor::InbuiltPP) name+="_ToolPP";
   else if (l_eff_cor==Bins::EffCor::TrigCorrUp  ) name+="_TrigCorrUp";
   else if (l_eff_cor==Bins::EffCor::TrigCorrDown) name+="_TrigCorrDown";
   else Common::Exception(__LINE__,__FILE__," Unknown Efficiency Correction");

   if     (l_trig==Trigger::AllTrigs ) name+="";
   else if(l_trig==Trigger::Only2Mu4 ) name+="_2mu4";
   else Common::Exception(__LINE__,__FILE__," unknown itrig");

   if     (l_deta==DetaCut::Deta0p8 ) name+="";
   else if(l_deta==DetaCut::Deta0p9 ) name+="_Deta0p9";
   else if(l_deta==DetaCut::Deta1p0 ) name+="_Deta1p0";
   else if(l_deta==DetaCut::Deta1p2 ) name+="_Deta1p2";
   else if(l_deta==DetaCut::Deta1p4 ) name+="_Deta1p4";
   else if(l_deta==DetaCut::Deta1p6 ) name+="_Deta1p6";
   else if(l_deta==DetaCut::Deta1p8 ) name+="_Deta1p8";
   else if(l_deta==DetaCut::Deta2p0 ) name+="_Deta2p0";
   else if(l_deta==DetaCut::Deta2p2 ) name+="_Deta2p2";
   else Common::Exception(__LINE__,__FILE__," unknown deta");

   if     (l_pdf==PDF::NoDYCorr) name+="_NoDY";
   else if(l_pdf==PDF::nNNPDF20) name+="";//"_nNNPDF20";
   else if(l_pdf==PDF::nCTEQ15 ) name+="_nCTEQ15" ;
   else if(l_pdf==PDF::NNPDF30 ) name+="_NNPDF30" ;
   else Common::Exception(__LINE__,__FILE__," unknown pdf:"+to_string(l_pdf));

   if(l_other_flag!=0){
     name+="_other";
     name+=std::to_string(l_other_flag);
   }
   return name;
}







std::string replaceChar(string str, char ch1, char ch2) {
  //char str[600];
  //sprintf(str,"%.1f",range);
  //string str = to_string(range);
  for (unsigned int i = 0; i < str.length(); ++i) {
    if (str[i] == ch1)
      str[i] = ch2;
  }
  return str;
}


//https://twiki.cern.ch/twiki/pub/AtlasProtected/HeavyIonAnalysis2015/centrality_values_Gv32_update_26Jul19.txt
std::pair<double,double> GetTAB(int icent){
  //5% wide bins
  if     (CENT_LO[icent]== 0 && CENT_HI[icent]==  5) return std::make_pair(26.0271   ,0.119849      );
  else if(CENT_LO[icent]== 5 && CENT_HI[icent]== 10) return std::make_pair(20.4019   ,0.133136      );
  else if(CENT_LO[icent]==10 && CENT_HI[icent]== 15) return std::make_pair(16.0997   ,0.133496      );
  else if(CENT_LO[icent]==15 && CENT_HI[icent]== 20) return std::make_pair(12.6593   ,0.133066      );
  //10% wide bins
  else if(CENT_LO[icent]== 0 && CENT_HI[icent]== 10) return std::make_pair(23.2145    ,0.124284     );
  else if(CENT_LO[icent]==10 && CENT_HI[icent]== 20) return std::make_pair(14.3795    ,0.132406     );
  else if(CENT_LO[icent]==20 && CENT_HI[icent]== 30) return std::make_pair( 8.76883   ,0.13243      );
  else if(CENT_LO[icent]==30 && CENT_HI[icent]== 40) return std::make_pair( 5.08891   ,0.11986      );
  else if(CENT_LO[icent]==40 && CENT_HI[icent]== 50) return std::make_pair( 2.74507   ,0.0946701    );
  else if(CENT_LO[icent]==50 && CENT_HI[icent]== 60) return std::make_pair( 1.35167   ,0.0651435    );
  //20%wide bins
  else if(CENT_LO[icent]==60 && CENT_HI[icent]== 80) return std::make_pair( 0.420011  ,0.0296767    );
  else if(CENT_LO[icent]==80 && CENT_HI[icent]==100) return std::make_pair( 0.0548188 ,0.00539919   );
  //100% wide bin
  else if(CENT_LO[icent]== 0 && CENT_HI[icent]==100) return std::make_pair( 5.64981387,0.);

  else Common::Exception(__LINE__,__FILE__);
  return std::make_pair(0.0,0.0);
}

//Multiplicity based centrality bins
float Multiplicity_Bins_HI_Tight[100]={
2925.02,  2795.89,  2688.74,  2588.16,  2493.03,  2402.52,  2315.93,  2232.9 ,  2153.02,  2076.17, 
2001.96,  1930.21,  1860.89,  1793.91,  1729.01,  1665.95,  1604.85,  1545.6 ,  1488   ,  1432.13, 
1377.84,  1325.06,  1273.86,  1224.18,  1175.95,  1129.16,  1083.58,  1039.34,  996.578,  955.089, 
914.719,  875.583,  837.645,  800.756,  765.063,  730.56 ,  697.111,  664.798,  633.495,  603.248, 
574.085,  545.872,  518.634,  492.399,  467.098,  442.728,  419.28 ,  396.701,  374.965,  354.134, 
334.209,  315.071,  296.762,  279.228,  262.433,  246.447,  231.194,  216.634,  202.767,  189.575, 
177.028,  165.138,  153.909,  143.261,  133.183,  123.677,  114.712,  106.247,  98.2822,  90.8317, 
83.8275,  77.2581,  71.1107,  65.3626,  60.0129,  55.0229,   50.382,  46.0484,  42.0099,  38.2649, 
 34.764,  31.5105,    28.49,  25.6756,  23.0646,  20.6384,  18.3798,  16.2808,  14.3403,  12.5676, 
 10.9803,  9.57461,  8.34568,  7.30048,  5.97227,  2.56114,      2.1,      2.1,      2.1,        0  //TODO (OK) last 4 bins are manually edited
 };
int GetCentralityHITight(int ntrk){
  for(int i=0;i<100;i++){
    if(ntrk>Multiplicity_Bins_HI_Tight[i]) return i;
  }
  return -1;
}


float FCal_ET_Bins_DATA[100]={
4.26258 ,4.08137,3.91763 ,3.7635  ,3.61844 ,3.48077  ,3.34945 ,3.22397 ,3.10407 ,2.98931 ,//0-10
2.87864 ,2.77237 ,2.66999,2.57162 ,2.47658 ,2.38468  ,2.29572 ,2.21002 ,2.12711 ,2.04651 ,//10-20
1.96859 ,1.89316 ,1.81997 ,1.74932,1.68058 ,1.61434  ,1.55005 ,1.48744 ,1.42719 ,1.36875 ,//20-30
1.31197 ,1.25693 ,1.20373 ,1.15214 ,1.10211,1.05367  ,1.0068  ,0.961609,0.917795,0.87541 ,//30-40
0.834538,0.795018,0.756791,0.719896,0.684377,0.65018 ,0.617108,0.585275,0.554569,0.525092,//40-50
0.49675 ,0.46959 ,0.443549,0.418573,0.394518,0.371561,0.349697,0.328744,0.308686,0.289595,//50-60
0.27137 ,0.25407 ,0.237615,0.22199 ,0.207148,0.193096,0.179776,0.167193,0.155307,0.14414 ,//60-70
0.133573,0.123657,0.114352,0.105619,0.097388,0.089723,0.082548,0.075838,0.06956 ,0.063719,//70-80
0.058250,0.053153,0.048423,0.044023,0.039948,0.036189,0.032726,0.029543,0.026643,0.024007,//80-90
0.020   ,   0.015,  0.010 ,   0.005,  0.0025, //90-95
0.000   ,      -1,    -2  ,      -3,      -4  //95-100
};
int GetCentrality(float FCal_Et){
  for(int i=0;i<100;i++){
    if(FCal_Et>FCal_ET_Bins_DATA[i]) return i;
  }
  return -1;
}

}
#endif

