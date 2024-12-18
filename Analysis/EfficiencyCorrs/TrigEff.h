#ifndef __TrigEff_h__
#define __TrigEff_h__

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "vector"
#include "TH1D.h"
#include "TH2D.h"
#include "Bins.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"



class TrigEff
{
private:
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   /*
   vector<float>   *vtx_z=nullptr;
   vector<float>   *vtx_x=nullptr;
   vector<float>   *vtx_y=nullptr;
   vector<int>     *vtx_ntrk=nullptr;
   */
   vector<float>   *muon_pt     =nullptr;
   vector<float>   *muon_eta    =nullptr;
   vector<int>     *muon_quality=nullptr;
   vector<float>   *muon_mspt   =nullptr;
   vector<float>   *muon_msp    =nullptr;
   vector<float>   *muon_d0     =nullptr;
   vector<float>   *muon_z0     =nullptr;

   Bool_t          b_HLT_mu4               =false;
   Bool_t          b_HLT_mu4_rerun_decision=false;
   Float_t         f_HLT_mu4_prescale      =0;
   vector<bool>    *muon_b_HLT_mu4   =nullptr;
   vector<bool>    *muon_b_HLT_mu4_V2=nullptr;


   //----------------------------------------
   Float_t         FCal_Et=0;
   Int_t           centrality=0;
   Bool_t          b_HLT_mb_sptrk_L1ZDC_A_C_VTE50=0;
   Bool_t          b_HLT_noalg_pc_L1TE50_VTE600_0ETA49=0;
   Bool_t          b_HLT_noalg_cc_L1TE600_0ETA49=0;
   Float_t         f_HLT_mb_sptrk_L1ZDC_A_C_VTE50_prescale=0;
   Float_t         f_HLT_noalg_pc_L1TE50_VTE600_0ETA49_prescale=0;
   Float_t         f_HLT_noalg_cc_L1TE600_0ETA49_prescale=0;

   Float_t         f_HLT_mu3_prescale=0;
   vector<bool>    *muon_b_HLT_mu3=nullptr;
   vector<bool>    *muon_b_HLT_mu3_V2=nullptr;
   //----------------------------------------


   //----------------------------------------
   vector<int>     *trk_numqual=nullptr;
   Bool_t          b_HLT_mb_sptrk=false;
   //PP 2017 5 TeV  
   Bool_t          b_HLT_mb_sp400_trk40_hmt_L1MBTS_1_1=false;
   Bool_t          b_HLT_mb_sp900_trk60_hmt_L1MBTS_1_1=false;
   Bool_t          b_HLT_mb_sp900_pusup400_trk60_hmt_L1TE5=false;
   Bool_t          b_HLT_mb_sp1000_pusup450_trk70_hmt_L1TE5=false;
   Bool_t          b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE10=false;
   Bool_t          b_HLT_mb_sp1800_pusup700_trk110_hmt_L1TE30=false;
   Bool_t          b_HLT_mb_sp600_trk40_hmt_L1TE5=false;
   Bool_t          b_HLT_mb_sp900_trk50_hmt_L1TE5=false;
   Bool_t          b_HLT_mb_sp900_trk60_hmt_L1TE5=false;
   Bool_t          b_HLT_mb_sp1000_trk70_hmt_L1TE5=false;
   Bool_t          b_HLT_mb_sp1400_trk90_hmt_L1TE5=false;
   Bool_t          b_HLT_mb_sp600_trk40_hmt_L1TE10=false;
   Bool_t          b_HLT_mb_sp700_trk50_hmt_L1TE10=false;
   Bool_t          b_HLT_mb_sp900_trk60_hmt_L1TE10=false;
   Bool_t          b_HLT_mb_sp1100_trk70_hmt_L1TE10=false;
   Bool_t          b_HLT_mb_sp1400_trk90_hmt_L1TE10=false;
   Bool_t          b_HLT_mb_sp700_trk50_hmt_L1TE20=false;
   Bool_t          b_HLT_mb_sp900_trk60_hmt_L1TE20=false;
   Bool_t          b_HLT_mb_sp1100_trk70_hmt_L1TE20=false;
   Bool_t          b_HLT_mb_sp1200_trk80_hmt_L1TE20=false;
   Bool_t          b_HLT_mb_sp1400_trk90_hmt_L1TE20=false;
   Bool_t          b_HLT_mb_sp1100_trk70_hmt_L1TE30=false;
   Bool_t          b_HLT_mb_sp1200_trk80_hmt_L1TE30=false;
   Bool_t          b_HLT_mb_sp1400_trk90_hmt_L1TE30=false;
   Bool_t          b_HLT_mb_sp1400_trk90_hmt_L1TE40=false;
   Bool_t          b_HLT_mb_sp1600_trk100_hmt_L1TE40=false;
   Bool_t          b_HLT_mb_sp1900_trk120_hmt_L1TE40=false;
   Bool_t          b_HLT_mb_sp1600_trk100_hmt_L1TE50=false;
   Bool_t          b_HLT_mb_sp1700_trk110_hmt_L1TE50=false;
   Bool_t          b_HLT_mb_sp1900_trk120_hmt_L1TE50=false;
   Bool_t          b_HLT_noalg_mb_L1MBTS_1=false;
   Bool_t          b_HLT_noalg_mb_L1MBTS_1_1=false;
   //PP 2017 13 TeV
   Bool_t          b_HLT_mb_sp700_pusup350_trk50_hmt_L1TE10=false;
   Bool_t          b_HLT_mb_sp1100_pusup450_trk70_hmt_L1TE20=false;
   Bool_t          b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE30=false;
   Bool_t          b_HLT_mb_sp1600_pusup600_trk100_hmt_L1TE40=false;
   Bool_t          b_HLT_mb_sp1700_pusup650_trk110_hmt_L1TE50=false;
   //----------------------------------------



   TFile *m_outfile=nullptr;
   void Init();
   void InitHists(bool read_flag);
   int  GetEtaBin(float eta,float charge);
   void ProcessSingles();
   void MakeEffFits();

   //configurables
   Bins::QualityCut m_muon_selection_cut=Bins::QualityCut::MEDIUM;
   std::vector<TObject**> m_all_objects;

   TFile *File=nullptr;

public :
   enum TrigEffBins{
     N_EFF_TYPE=2,
       PRESCALE_CORRECTED=0,
       RERUN             =1,
     NETA  =9,
     NTRIGS=2,
       HLT_MU4=0,
       HLT_MU3=1,
     N_CENT=14,
   };
   std::map<int,std::string> LabelEta={
     {0,"-2.4<q*#eta<-2.0"},
     {1,"-2.0<q*#eta<-1.5"},
     {2,"-1.5<q*#eta<-1.0"},
     {3,"-1.0<q*#eta<-0.5"},
     {4,"-0.5<q*#eta<0.5"},
     {5," 0.5<q*#eta<1.0"},
     {6," 1.0<q*#eta<1.5"},
     {7," 1.5<q*#eta<2.0"},
     {8," 2.0<q*#eta<2.4"},
   };
   std::map<int,std::string> LabelTrig={{HLT_MU4,"HLT_MU4"},{HLT_MU3,"HLT_MU3"}};
   std::map<int,std::string> LabelType={{PRESCALE_CORRECTED,"PS Weighted"},{HLT_MU3,"Rerun"}};

   void Plot();
   void PlotEff  (int option);
   void Plot2DEff(int itype,int itrig);
   std::vector<TCanvas*> m_can_vec;
   std::string m_name;
   int m_data_type=Bins::DataType::PbPbAll;

   //returns single-muon Trigger efficiency
   float GetTriggerEfficiency(float eta, float pt, float charge,int itrig, int centrality_percentile);

   TrigEff(int read_only_flag    =1,
           int muon_selection_cut=static_cast<int>(Bins::QualityCut::MEDIUM),
           int data_type         =Bins::DataType    ::PbPbAll);

   ~TrigEff(){
     for(auto *obj:m_all_objects){
       if(*obj){
         delete *obj;
         *obj=nullptr;
       }
     }
   }
   TH1D*h_FCal_Et[2]={nullptr};
//private:
   TH1D             *h_Eff_Den   [TrigEffBins::N_EFF_TYPE][TrigEffBins::NETA][TrigEffBins::NTRIGS]={{{nullptr}}};
   TH1D             *h_Eff_Num   [TrigEffBins::N_EFF_TYPE][TrigEffBins::NETA][TrigEffBins::NTRIGS]={{{nullptr}}};
   TGraphAsymmErrors*gr_Eff      [TrigEffBins::N_EFF_TYPE][TrigEffBins::NETA][TrigEffBins::NTRIGS]={{{nullptr}}};
   TF1              *tf1_eff_fit [TrigEffBins::N_EFF_TYPE][TrigEffBins::NETA][TrigEffBins::NTRIGS]={{{nullptr}}};
   TH2D             *h_Eff_Den_2D[TrigEffBins::N_EFF_TYPE]                   [TrigEffBins::NTRIGS]={{nullptr}};
   TH2D             *h_Eff_Num_2D[TrigEffBins::N_EFF_TYPE]                   [TrigEffBins::NTRIGS]={{nullptr}};
   TH1D*h_cent   [2]={nullptr};
   TH1D*h_cent_v1[2]={nullptr};
   TH1D*h_cent_v2[2]={nullptr};
   TH1D*h_ntrk      = nullptr;

   //----------------------------------------------------------
   TH1D             *h_Eff_Den_cent  [TrigEffBins::N_EFF_TYPE][TrigEffBins::NETA][TrigEffBins::NTRIGS][TrigEffBins::N_CENT]={{{{nullptr}}}};
   TH1D             *h_Eff_Num_cent  [TrigEffBins::N_EFF_TYPE][TrigEffBins::NETA][TrigEffBins::NTRIGS][TrigEffBins::N_CENT]={{{{nullptr}}}};
   //TGraphAsymmErrors*gr_Eff_cent     [TrigEffBins::N_EFF_TYPE][TrigEffBins::NETA][TrigEffBins::NTRIGS][TrigEffBins::N_CENT]={{{{nullptr}}}};
   //TF1              *tf1_eff_fit_cent[TrigEffBins::N_EFF_TYPE][TrigEffBins::NETA][TrigEffBins::NTRIGS][TrigEffBins::N_CENT]={{{{nullptr}}}};
   //eta-integrated
   TH1D             *h_Eff_Den_cent2 [TrigEffBins::N_EFF_TYPE]                   [TrigEffBins::NTRIGS][TrigEffBins::N_CENT]={{{nullptr}}};
   TH1D             *h_Eff_Num_cent2 [TrigEffBins::N_EFF_TYPE]                   [TrigEffBins::NTRIGS][TrigEffBins::N_CENT]={{{nullptr}}};
   double CENT_DN[TrigEffBins::N_CENT]={ 0,20, 40,  0,  0,10,20,30,40,50,60, 80, 60,40};
   double CENT_UP[TrigEffBins::N_CENT]={20,40,100,100, 10,20,30,40,50,60,80,100,100,60};
   bool IsInCentBin(int icent,int ibin){
      if (ibin>=N_CENT) Common::Exception(__LINE__,__FILE__);
     if(icent>=CENT_DN[ibin] && icent<CENT_UP[ibin]) return true;
     return false;
   }
   int GetCentIndexTrig(double cent_lo,double cent_hi){
     for(int i=0;i<N_CENT;i++){
       if(CENT_DN[i]==cent_lo && CENT_UP[i]==cent_hi) return i;
     }
     Common::Exception(__LINE__,__FILE__);
     return -1;
   }
   void InitHistsCentDepEffs(bool read_flag);
   void PlotEffCentDep(double ptmin, double ptmax,int option=0);
   //----------------------------------------------------------
};


void TrigEff::Init()
{
   fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
   fChain->SetMakeClass(1);

   switch(m_data_type){
     case Bins::DataType::PbPbAll:
     {
       fChain->Add("Data/ForTrigEffs/MinBias/*.root");
       //fChain->Add("Data/user.soumya.TrigRates.physics_MinBias.PbPb2018.27April2019.AthAnalysis.DimuonMatch._MYSTREAM/all*.root");
       break;
     }
     case Bins::DataType::pp2017:
     case Bins::DataType::pp2017_MB:
     case Bins::DataType::pp2017_HMT:
     {
       //fChain->Add("Data/PP2017_5TeV/TrigRates_ForTrigEff.physics_Main.PP2017_5TeV.10Feb2022.root");
       fChain->Add("Data/PP2017_5TeV/user.soumya.TrigRates_ForTrigEff.physics_Main.PP2017_5TeV.22Feb2022.v1_MYSTREAM/*");
       break;
     }
     case Bins::DataType::pp2017_13TeV:
     case Bins::DataType::pp2017_13TeV_MB:
     case Bins::DataType::pp2017_13TeV_HMT:
     {
       fChain->Add("Data/PP2017_5TeV/user.soumya.TrigRates.PPMu22017.MinBias.Feb10.2022_MYSTREAM/*");
       break;
     }
     default: Common::Exception(__LINE__,__FILE__);
   }//switch

   // Declaration of leaf types
   fChain->SetMakeClass(1);
   //fChain->SetBranchAddress("vtx_z", &vtx_z);
   //fChain->SetBranchAddress("vtx_x", &vtx_x);
   //fChain->SetBranchAddress("vtx_y", &vtx_y);
   //fChain->SetBranchAddress("vtx_ntrk", &vtx_ntrk);
   fChain->SetBranchAddress("muon_pt"     , &muon_pt     );
   fChain->SetBranchAddress("muon_eta"    , &muon_eta    );
   fChain->SetBranchAddress("muon_quality", &muon_quality);
   fChain->SetBranchAddress("muon_mspt"   , &muon_mspt   );
   fChain->SetBranchAddress("muon_msp"    , &muon_msp    );
   fChain->SetBranchAddress("muon_d0"     , &muon_d0     );
   fChain->SetBranchAddress("muon_z0"     , &muon_z0     );



   switch(m_data_type){
     case Bins::DataType::PbPbAll:
     {
       fChain->SetBranchAddress("muon_b_HLT_mu4", &muon_b_HLT_mu4);
       fChain->SetBranchAddress("muon_b_HLT_mu4_V2", &muon_b_HLT_mu4_V2);

       fChain->SetBranchAddress("b_HLT_mu4",&b_HLT_mu4);
       fChain->SetBranchAddress("b_HLT_mu4_rerun_decision", &b_HLT_mu4_rerun_decision);
       fChain->SetBranchAddress("f_HLT_mu4_prescale", &f_HLT_mu4_prescale);

       //fChain->SetBranchAddress("b_HLT_mu3",&b_HLT_mu3);
       fChain->SetBranchAddress("f_HLT_mu3_prescale", &f_HLT_mu3_prescale);
       fChain->SetBranchAddress("muon_b_HLT_mu3", &muon_b_HLT_mu3);
       fChain->SetBranchAddress("muon_b_HLT_mu3_V2", &muon_b_HLT_mu3_V2);


       fChain->SetBranchAddress("FCal_Et", &FCal_Et);
       fChain->SetBranchAddress("centrality", &centrality);
       fChain->SetBranchAddress("b_HLT_mb_sptrk_L1ZDC_A_C_VTE50", &b_HLT_mb_sptrk_L1ZDC_A_C_VTE50);
       fChain->SetBranchAddress("b_HLT_noalg_pc_L1TE50_VTE600.0ETA49", &b_HLT_noalg_pc_L1TE50_VTE600_0ETA49);
       fChain->SetBranchAddress("b_HLT_noalg_cc_L1TE600.0ETA49", &b_HLT_noalg_cc_L1TE600_0ETA49);
       fChain->SetBranchAddress("f_HLT_mb_sptrk_L1ZDC_A_C_VTE50_prescale", &f_HLT_mb_sptrk_L1ZDC_A_C_VTE50_prescale);
       fChain->SetBranchAddress("f_HLT_noalg_pc_L1TE50_VTE600.0ETA49_prescale", &f_HLT_noalg_pc_L1TE50_VTE600_0ETA49_prescale);
       fChain->SetBranchAddress("f_HLT_noalg_cc_L1TE600.0ETA49_prescale", &f_HLT_noalg_cc_L1TE600_0ETA49_prescale);


       break;
     }
     case Bins::DataType::pp2017:
     case Bins::DataType::pp2017_MB:
     case Bins::DataType::pp2017_HMT:
     {
       fChain->SetBranchAddress("muon_b_HLT_mu4", &muon_b_HLT_mu4);
       fChain->SetBranchAddress("muon_b_HLT_mu4_V3", &muon_b_HLT_mu4_V2);
       //fChain->SetBranchAddress("muon_b_HLT_mu4_V2", &muon_b_HLT_mu4_V2);
       fChain->SetBranchAddress("b_HLT_mu4",&b_HLT_mu4);
       fChain->SetBranchAddress("b_HLT_mu4_rerun_decision", &b_HLT_mu4_rerun_decision);
       fChain->SetBranchAddress("f_HLT_mu4_prescale", &f_HLT_mu4_prescale);


       fChain->SetBranchAddress("trk_numqual", &trk_numqual);
       fChain->SetBranchAddress("b_HLT_mb_sptrk", &b_HLT_mb_sptrk);
       fChain->SetBranchAddress("b_HLT_noalg_mb_L1MBTS_1", &b_HLT_noalg_mb_L1MBTS_1);
       fChain->SetBranchAddress("b_HLT_noalg_mb_L1MBTS_1_1", &b_HLT_noalg_mb_L1MBTS_1_1);

       fChain->SetBranchAddress("b_HLT_mb_sp400_trk40_hmt_L1MBTS_1_1", &b_HLT_mb_sp400_trk40_hmt_L1MBTS_1_1);
       fChain->SetBranchAddress("b_HLT_mb_sp900_trk60_hmt_L1MBTS_1_1", &b_HLT_mb_sp900_trk60_hmt_L1MBTS_1_1);
       fChain->SetBranchAddress("b_HLT_mb_sp900_pusup400_trk60_hmt_L1TE5", &b_HLT_mb_sp900_pusup400_trk60_hmt_L1TE5);
       fChain->SetBranchAddress("b_HLT_mb_sp1000_pusup450_trk70_hmt_L1TE5", &b_HLT_mb_sp1000_pusup450_trk70_hmt_L1TE5);
       fChain->SetBranchAddress("b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE10", &b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE10);
       fChain->SetBranchAddress("b_HLT_mb_sp1800_pusup700_trk110_hmt_L1TE30", &b_HLT_mb_sp1800_pusup700_trk110_hmt_L1TE30);
       fChain->SetBranchAddress("b_HLT_mb_sp600_trk40_hmt_L1TE5", &b_HLT_mb_sp600_trk40_hmt_L1TE5);
       fChain->SetBranchAddress("b_HLT_mb_sp900_trk50_hmt_L1TE5", &b_HLT_mb_sp900_trk50_hmt_L1TE5);
       fChain->SetBranchAddress("b_HLT_mb_sp900_trk60_hmt_L1TE5", &b_HLT_mb_sp900_trk60_hmt_L1TE5);
       fChain->SetBranchAddress("b_HLT_mb_sp1000_trk70_hmt_L1TE5", &b_HLT_mb_sp1000_trk70_hmt_L1TE5);
       fChain->SetBranchAddress("b_HLT_mb_sp1400_trk90_hmt_L1TE5", &b_HLT_mb_sp1400_trk90_hmt_L1TE5);
       fChain->SetBranchAddress("b_HLT_mb_sp600_trk40_hmt_L1TE10", &b_HLT_mb_sp600_trk40_hmt_L1TE10);
       fChain->SetBranchAddress("b_HLT_mb_sp700_trk50_hmt_L1TE10", &b_HLT_mb_sp700_trk50_hmt_L1TE10);
       fChain->SetBranchAddress("b_HLT_mb_sp900_trk60_hmt_L1TE10", &b_HLT_mb_sp900_trk60_hmt_L1TE10);
       fChain->SetBranchAddress("b_HLT_mb_sp1100_trk70_hmt_L1TE10", &b_HLT_mb_sp1100_trk70_hmt_L1TE10);
       fChain->SetBranchAddress("b_HLT_mb_sp1400_trk90_hmt_L1TE10", &b_HLT_mb_sp1400_trk90_hmt_L1TE10);
       fChain->SetBranchAddress("b_HLT_mb_sp700_trk50_hmt_L1TE20", &b_HLT_mb_sp700_trk50_hmt_L1TE20);
       fChain->SetBranchAddress("b_HLT_mb_sp900_trk60_hmt_L1TE20", &b_HLT_mb_sp900_trk60_hmt_L1TE20);
       fChain->SetBranchAddress("b_HLT_mb_sp1100_trk70_hmt_L1TE20", &b_HLT_mb_sp1100_trk70_hmt_L1TE20);
       fChain->SetBranchAddress("b_HLT_mb_sp1200_trk80_hmt_L1TE20", &b_HLT_mb_sp1200_trk80_hmt_L1TE20);
       fChain->SetBranchAddress("b_HLT_mb_sp1400_trk90_hmt_L1TE20", &b_HLT_mb_sp1400_trk90_hmt_L1TE20);
       fChain->SetBranchAddress("b_HLT_mb_sp1100_trk70_hmt_L1TE30", &b_HLT_mb_sp1100_trk70_hmt_L1TE30);
       fChain->SetBranchAddress("b_HLT_mb_sp1200_trk80_hmt_L1TE30", &b_HLT_mb_sp1200_trk80_hmt_L1TE30);
       fChain->SetBranchAddress("b_HLT_mb_sp1400_trk90_hmt_L1TE30", &b_HLT_mb_sp1400_trk90_hmt_L1TE30);
       fChain->SetBranchAddress("b_HLT_mb_sp1400_trk90_hmt_L1TE40", &b_HLT_mb_sp1400_trk90_hmt_L1TE40);
       fChain->SetBranchAddress("b_HLT_mb_sp1600_trk100_hmt_L1TE40", &b_HLT_mb_sp1600_trk100_hmt_L1TE40);
       fChain->SetBranchAddress("b_HLT_mb_sp1900_trk120_hmt_L1TE40", &b_HLT_mb_sp1900_trk120_hmt_L1TE40);
       fChain->SetBranchAddress("b_HLT_mb_sp1600_trk100_hmt_L1TE50", &b_HLT_mb_sp1600_trk100_hmt_L1TE50);
       fChain->SetBranchAddress("b_HLT_mb_sp1700_trk110_hmt_L1TE50", &b_HLT_mb_sp1700_trk110_hmt_L1TE50);
       fChain->SetBranchAddress("b_HLT_mb_sp1900_trk120_hmt_L1TE50", &b_HLT_mb_sp1900_trk120_hmt_L1TE50);
       break;
     }
     case Bins::DataType::pp2017_13TeV:
     case Bins::DataType::pp2017_13TeV_MB:
     case Bins::DataType::pp2017_13TeV_HMT:
     {
       //NOTE HLT_mu4_L1MU4 and not "HLT_mu4"
       fChain->SetBranchAddress("muon_b_HLT_mu4_L1MU4", &muon_b_HLT_mu4);
       fChain->SetBranchAddress("muon_b_HLT_mu4_L1MU4_V2", &muon_b_HLT_mu4_V2);
       fChain->SetBranchAddress("b_HLT_mu4_L1MU4",&b_HLT_mu4);
       fChain->SetBranchAddress("b_HLT_mu4_L1MU4_rerun_decision", &b_HLT_mu4_rerun_decision);
       fChain->SetBranchAddress("f_HLT_mu4_L1MU4_prescale", &f_HLT_mu4_prescale);

       fChain->SetBranchAddress("trk_numqual", &trk_numqual);
       fChain->SetBranchAddress("b_HLT_mb_sptrk", &b_HLT_mb_sptrk);
       fChain->SetBranchAddress("b_HLT_mb_sp700_pusup350_trk50_hmt_L1TE10"  , &b_HLT_mb_sp700_pusup350_trk50_hmt_L1TE10);
       fChain->SetBranchAddress("b_HLT_mb_sp1100_pusup450_trk70_hmt_L1TE20" , &b_HLT_mb_sp1100_pusup450_trk70_hmt_L1TE20);
       fChain->SetBranchAddress("b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE30" , &b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE30);
       fChain->SetBranchAddress("b_HLT_mb_sp1600_pusup600_trk100_hmt_L1TE40", &b_HLT_mb_sp1600_pusup600_trk100_hmt_L1TE40);
       fChain->SetBranchAddress("b_HLT_mb_sp1700_pusup650_trk110_hmt_L1TE50", &b_HLT_mb_sp1700_pusup650_trk110_hmt_L1TE50);
       break;
     }
     default: Common::Exception(__LINE__,__FILE__);
   }//switch 



   fChain->SetBranchStatus("*"           ,0);
   fChain->SetBranchStatus("muon_pt"     ,1);
   fChain->SetBranchStatus("muon_eta"    ,1);
   fChain->SetBranchStatus("muon_quality",1);
   fChain->SetBranchStatus("muon_mspt"   ,1);
   fChain->SetBranchStatus("muon_msp"    ,1);
   fChain->SetBranchStatus("muon_d0"     ,1);
   fChain->SetBranchStatus("muon_z0"     ,1);

   
   switch(m_data_type){
     case Bins::DataType::PbPbAll:
     {
       fChain->SetBranchStatus("muon_b_HLT_mu4",1);
       fChain->SetBranchStatus("muon_b_HLT_mu4_V2",1);
       fChain->SetBranchStatus("f_HLT_mu4_prescale",1);
       fChain->SetBranchStatus("b_HLT_mu4_rerun_decision",1);
       fChain->SetBranchStatus("b_HLT_mu4",1);

       fChain->SetBranchStatus("muon_b_HLT_mu3",1);
       fChain->SetBranchStatus("muon_b_HLT_mu3_V2",1);
       fChain->SetBranchStatus("f_HLT_mu3_prescale",1);


       fChain->SetBranchStatus("FCal_Et" ,1);
       fChain->SetBranchStatus("centrality",1);
       fChain->SetBranchStatus("b_HLT_mb_sptrk_L1ZDC_A_C_VTE50",1);
       fChain->SetBranchStatus("b_HLT_noalg_pc_L1TE50_VTE600.0ETA49",1);
       fChain->SetBranchStatus("b_HLT_noalg_cc_L1TE600.0ETA49",1);
       fChain->SetBranchStatus("f_HLT_mb_sptrk_L1ZDC_A_C_VTE50_prescale"     ,1);
       fChain->SetBranchStatus("f_HLT_noalg_pc_L1TE50_VTE600.0ETA49_prescale",1);
       fChain->SetBranchStatus("f_HLT_noalg_cc_L1TE600.0ETA49_prescale"      ,1);

       break;
     }
     case Bins::DataType::pp2017:
     case Bins::DataType::pp2017_MB:
     case Bins::DataType::pp2017_HMT:
     {
       fChain->SetBranchStatus("muon_b_HLT_mu4",1);
       fChain->SetBranchStatus("muon_b_HLT_mu4_V3",1);
       //fChain->SetBranchStatus("muon_b_HLT_mu4_V2",1);
       fChain->SetBranchStatus("f_HLT_mu4_prescale",1);
       fChain->SetBranchStatus("b_HLT_mu4_rerun_decision",1);
       fChain->SetBranchStatus("b_HLT_mu4",1);

       fChain->SetBranchStatus("trk_numqual", 1);
       fChain->SetBranchStatus("b_HLT_mb_sptrk",1);
       fChain->SetBranchStatus("b_HLT_noalg_mb_L1MBTS_1",1);
       fChain->SetBranchStatus("b_HLT_noalg_mb_L1MBTS_1_1",1);

       fChain->SetBranchStatus("b_HLT_mb_sp400_trk40_hmt_L1MBTS_1_1",1);
       fChain->SetBranchStatus("b_HLT_mb_sp900_trk60_hmt_L1MBTS_1_1",1);
       fChain->SetBranchStatus("b_HLT_mb_sp900_pusup400_trk60_hmt_L1TE5",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1000_pusup450_trk70_hmt_L1TE5",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE10",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1800_pusup700_trk110_hmt_L1TE30",1);
       fChain->SetBranchStatus("b_HLT_mb_sp600_trk40_hmt_L1TE5",1);
       fChain->SetBranchStatus("b_HLT_mb_sp900_trk50_hmt_L1TE5",1);
       fChain->SetBranchStatus("b_HLT_mb_sp900_trk60_hmt_L1TE5",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1000_trk70_hmt_L1TE5",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1400_trk90_hmt_L1TE5",1);
       fChain->SetBranchStatus("b_HLT_mb_sp600_trk40_hmt_L1TE10",1);
       fChain->SetBranchStatus("b_HLT_mb_sp700_trk50_hmt_L1TE10",1);
       fChain->SetBranchStatus("b_HLT_mb_sp900_trk60_hmt_L1TE10",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1100_trk70_hmt_L1TE10",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1400_trk90_hmt_L1TE10",1);
       fChain->SetBranchStatus("b_HLT_mb_sp700_trk50_hmt_L1TE20",1);
       fChain->SetBranchStatus("b_HLT_mb_sp900_trk60_hmt_L1TE20",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1100_trk70_hmt_L1TE20",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1200_trk80_hmt_L1TE20",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1400_trk90_hmt_L1TE20",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1100_trk70_hmt_L1TE30",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1200_trk80_hmt_L1TE30",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1400_trk90_hmt_L1TE30",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1400_trk90_hmt_L1TE40",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1600_trk100_hmt_L1TE40",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1900_trk120_hmt_L1TE40",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1600_trk100_hmt_L1TE50",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1700_trk110_hmt_L1TE50",1);
       fChain->SetBranchStatus("b_HLT_mb_sp1900_trk120_hmt_L1TE50",1);
       break;
     }
     case Bins::DataType::pp2017_13TeV:
     case Bins::DataType::pp2017_13TeV_MB:
     case Bins::DataType::pp2017_13TeV_HMT:
     {
       fChain->SetBranchStatus("muon_b_HLT_mu4_L1MU4",1);
       fChain->SetBranchStatus("muon_b_HLT_mu4_L1MU4_V2",1);
       fChain->SetBranchStatus("f_HLT_mu4_L1MU4_prescale",1);
       fChain->SetBranchStatus("b_HLT_mu4_L1MU4_rerun_decision",1);
       fChain->SetBranchStatus("b_HLT_mu4_L1MU4",1);

       fChain->SetBranchStatus("trk_numqual", 1);
       fChain->SetBranchStatus("b_HLT_mb_sptrk", 1);
       fChain->SetBranchStatus("b_HLT_mb_sp700_pusup350_trk50_hmt_L1TE10"  , 1);
       fChain->SetBranchStatus("b_HLT_mb_sp1100_pusup450_trk70_hmt_L1TE20" , 1);
       fChain->SetBranchStatus("b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE30" , 1);
       fChain->SetBranchStatus("b_HLT_mb_sp1600_pusup600_trk100_hmt_L1TE40", 1);
       fChain->SetBranchStatus("b_HLT_mb_sp1700_pusup650_trk110_hmt_L1TE50", 1);
       break;
     }
     default: Common::Exception(__LINE__,__FILE__);
   }//switch
}

#endif



