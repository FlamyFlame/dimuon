#ifndef __TrigEff_C__
#define __TrigEff_C__

#include "TrigEff.h"
#include "Bins.h"
#include "TStyle.h"

#define  CREATE(VAR_NAME,TYPE,SUMW2,NAME,...) \
  if(!read_flag) {\
    VAR_NAME=new TYPE((NAME).c_str(),__VA_ARGS__);\
    if(SUMW2) VAR_NAME->Sumw2();\
  }\
  else{\
    VAR_NAME=(TYPE*)File->Get((NAME).c_str());\
    Common::CheckObject2(VAR_NAME,NAME,File);\
  }\
  m_all_objects.push_back((TObject**)&VAR_NAME)

#define  CREATE_GR(VAR_NAME,NAME,...) \
  if(!read_flag) {\
    VAR_NAME=new TGraphAsymmErrors(__VA_ARGS__);\
    VAR_NAME->SetName((NAME).c_str());\
  }\
  else{\
    VAR_NAME=(TGraphAsymmErrors*)File->Get((NAME).c_str());\
    Common::CheckObject2(VAR_NAME,NAME,File);\
  }\
  m_all_objects.push_back((TObject**)&VAR_NAME)

#define READ_TF1(VAR_NAME,NAME) \
  if(read_flag) {\
    VAR_NAME=(TF1*)File->Get((NAME).c_str());\
    Common::CheckObject2(VAR_NAME,NAME,File);\
  }\
  m_all_objects.push_back((TObject**)&VAR_NAME)


void TrigEff::InitHists(bool read_flag)
{
  std::string ext;
  const int NBINSPT=20;
  double eff_pt_bins[NBINSPT+1]={3.5,3.6,3.7,3.8,3.9,4,4.1,4.25,4.5,5,5.5,6,7,8,9,10,12,14,16,18,20};
  for(int itype=0;itype<TrigEffBins::N_EFF_TYPE;itype++){
    for(int itrig=0;itrig<TrigEffBins::NTRIGS;itrig++){
      for(int ieta=0;ieta<TrigEffBins::NETA;ieta++){
        ext="_type"+std::to_string(itype)+"_eta"+std::to_string(ieta)+"_trig"+std::to_string(itrig);
        CREATE   (h_Eff_Num   [itype][ieta][itrig],TH1D,true,("h_Eff_Num"  +ext),";#it{p}_{T} [GeV];",NBINSPT,eff_pt_bins);
        CREATE   (h_Eff_Den   [itype][ieta][itrig],TH1D,true,("h_Eff_Den"  +ext),";#it{p}_{T} [GeV];",NBINSPT,eff_pt_bins);
        CREATE_GR(gr_Eff      [itype][ieta][itrig],("gr_Eff"     +ext),NBINSPT);
        READ_TF1 (tf1_eff_fit [itype][ieta][itrig],("tf1_eff_fit"+ext));
      }
      ext="_type"+std::to_string(itype)+"_trig"+std::to_string(itrig);
      CREATE(h_Eff_Num_2D[itype][itrig],TH2D,true,("h_Eff_Num_2D"+ext),";q*#eta;#it{p}_{T} [GeV];",48,-2.4,2.4,NBINSPT,eff_pt_bins);
      CREATE(h_Eff_Den_2D[itype][itrig],TH2D,true,("h_Eff_Den_2D"+ext),";q*#eta;#it{p}_{T} [GeV];",48,-2.4,2.4,NBINSPT,eff_pt_bins);
    }
  }
  for(int weight:{0,1}){
    ext="_PSweight"+std::to_string(weight);
    CREATE(h_cent   [weight],TH1D,true,("h_cent"   +ext),";Centrality [%];"      ,100 ,-0.5,99.5);
    CREATE(h_cent_v1[weight],TH1D,true,("h_cent_v1"+ext),";Centrality [%];"      ,100 ,-0.5,99.5);
    CREATE(h_cent_v2[weight],TH1D,true,("h_cent_v2"+ext),";Centrality [%];"      ,100 ,-0.5,99.5);
    CREATE(h_FCal_Et[weight],TH1D,true,("h_FCal_Et"+ext),";FCal #it{E}_{T}[TeV];",5500,-0.5,5.0 );
  }
  ext="";
  CREATE(h_ntrk ,TH1D,true,("h_ntrk"   +ext),";Nchrec;"      ,300 ,-0.5,299.5);
}

void TrigEff::InitHistsCentDepEffs(bool read_flag){
  const int NBINSPT=20;
  double eff_pt_bins[NBINSPT+1]={3.5,3.6,3.7,3.8,3.9,4,4.1,4.25,4.5,5,5.5,6,7,8,9,10,12,14,16,18,20};
  for(int itype=0;itype<TrigEffBins::N_EFF_TYPE;itype++){
    for(int itrig=0;itrig<TrigEffBins::NTRIGS;itrig++){
      for(int ieta=0;ieta<TrigEffBins::NETA;ieta++){
        for(int icent=0;icent<TrigEffBins::N_CENT;icent++){
          string ext="_type"+std::to_string(itype)+"_eta"+std::to_string(ieta)+"_trig"+std::to_string(itrig)+"_cent"+std::to_string(icent);
          CREATE   (h_Eff_Num_cent   [itype][ieta][itrig][icent],TH1D,true,("h_Eff_Num"  +ext),";#it{p}_{T} [GeV];",NBINSPT,eff_pt_bins);
          CREATE   (h_Eff_Den_cent   [itype][ieta][itrig][icent],TH1D,true,("h_Eff_Den"  +ext),";#it{p}_{T} [GeV];",NBINSPT,eff_pt_bins);
          //CREATE_GR(gr_Eff_cent      [itype][ieta][itrig][icent],("gr_Eff"     +ext),NBINSPT);
          if(ieta==0){
                 ext="_type"+std::to_string(itype)+"_trig"+std::to_string(itrig)+"_cent"+std::to_string(icent);
            CREATE   (h_Eff_Num_cent2  [itype]      [itrig][icent],TH1D,true,("h_Eff_Num"  +ext),";#it{p}_{T} [GeV];",NBINSPT,eff_pt_bins);
            CREATE   (h_Eff_Den_cent2  [itype]      [itrig][icent],TH1D,true,("h_Eff_Den"  +ext),";#it{p}_{T} [GeV];",NBINSPT,eff_pt_bins);
          }
        }
      }
    }
  }
}
#undef CREATE
#undef CREATE_GR
#undef READ_TF1

TrigEff::TrigEff(int read_flag,int muon_selection_cut,int data_type)
{

   m_muon_selection_cut=static_cast<Bins::QualityCut>(muon_selection_cut);
   m_data_type=data_type;
   switch (data_type){
     case Bins::DataType::PbPbAll:
       m_name="MuonSelectionCut"+std::to_string(muon_selection_cut);break;
     case Bins::DataType::pp2017:
       m_name="MuonSelectionCut"+std::to_string(muon_selection_cut)+"_pp2017";break;
     case Bins::DataType::pp2017_MB:
       m_name="MuonSelectionCut"+std::to_string(muon_selection_cut)+"_pp2017_MB";break;
     case Bins::DataType::pp2017_HMT:
       m_name="MuonSelectionCut"+std::to_string(muon_selection_cut)+"_pp2017_HMT";break;
     case Bins::DataType::pp2017_13TeV:
       m_name="MuonSelectionCut"+std::to_string(muon_selection_cut)+"_pp2017_13TeV";break;
     case Bins::DataType::pp2017_13TeV_MB:
       m_name="MuonSelectionCut"+std::to_string(muon_selection_cut)+"_pp2017_13TeV_MB";break;
     case Bins::DataType::pp2017_13TeV_HMT:
       m_name="MuonSelectionCut"+std::to_string(muon_selection_cut)+"_pp2017_13TeV_HMT";break;
     default : 
       Common::Exception(__LINE__,__FILE__);
   }



   //process data
   if(read_flag==0){
     Init();              //Initialize branches
     File=new TFile(("01RootFiles/TrigEff_"+m_name+".root").c_str(),"recreate");//Create output file
     InitHists(false);//Create histograms and TGraphs
     if(m_data_type==Bins::DataType::PbPbAll) InitHistsCentDepEffs (false);

     ProcessSingles();    //Prcess data to fill single-particle histograms
     MakeEffFits   ();    //Fit efficiency graphs to get fit functions

     File->Write();       //Write output file
     for(int itype=0;itype<TrigEffBins::N_EFF_TYPE;itype++){
       for(int ieta=0;ieta<TrigEffBins::NETA;ieta++){//Explicitly write TGraphs and TF1s
         for(int itrig=0;itrig<TrigEffBins::NTRIGS;itrig++){
           gr_Eff     [itype][ieta][itrig]->Write();
           tf1_eff_fit[itype][ieta][itrig]->Write();
         }
       }
     }
   }
   //read in histograms from already processed data file
   else if (read_flag==1){
     File=new TFile(("01RootFiles/TrigEff_"+m_name+".root").c_str(),"read");
     InitHists(true);
   }
   else if(read_flag==2){
     File=new TFile(("01RootFiles/TrigEff_"+m_name+".root").c_str(),"read");
     InitHists(true);
     Plot();
     Common::SaveCanvas(m_can_vec,"_"+m_name);
   }
   else if(read_flag==3 && m_data_type==Bins::DataType::PbPbAll){
     File=new TFile(("01RootFiles/TrigEff_"+m_name+".root").c_str(),"read");
     InitHists(true);
     InitHistsCentDepEffs (true);
     //PlotEffCentDep(4, 5,0);
     //PlotEffCentDep(5, 6,0);
     //PlotEffCentDep(6,20,0);
     PlotEffCentDep(4,20,0);
     Common::SaveCanvas(m_can_vec,"_"+m_name);
   }
}



void TrigEff::ProcessSingles()
{
   Long64_t nentries = fChain->GetEntries();
   int badcount_mu4[TrigEffBins::N_EFF_TYPE]={0};
   int badcount_mu3[TrigEffBins::N_EFF_TYPE]={0};
   for(Long64_t jentry=0; jentry<nentries;jentry++){
     int nb = fChain->GetEntry(jentry);
     if(jentry%100000==0) std::cout<<jentry<<" "<<nentries<<std::endl;
     int centrality_v1=0;
     int centrality_v2=0;

     switch (m_data_type){
       case Bins::DataType::PbPbAll:
       {
         if(!(b_HLT_noalg_cc_L1TE600_0ETA49 || b_HLT_noalg_pc_L1TE50_VTE600_0ETA49 || b_HLT_mb_sptrk_L1ZDC_A_C_VTE50)) continue;//require min-bias trigger

         float PS=b_HLT_mb_sptrk_L1ZDC_A_C_VTE50     *f_HLT_mb_sptrk_L1ZDC_A_C_VTE50_prescale 
                 +b_HLT_noalg_pc_L1TE50_VTE600_0ETA49*f_HLT_noalg_pc_L1TE50_VTE600_0ETA49_prescale 
                 +b_HLT_noalg_cc_L1TE600_0ETA49      *f_HLT_noalg_cc_L1TE600_0ETA49_prescale;
         centrality_v1=Bins::GetCentrality(FCal_Et       /1e6);
         centrality_v2=Bins::GetCentrality(FCal_Et*1.0015/1e6);
         h_cent   [0]->Fill(centrality      );
         h_cent   [1]->Fill(centrality   ,PS);
         h_cent_v1[0]->Fill(centrality_v1   );
         h_cent_v1[1]->Fill(centrality_v1,PS);
         h_cent_v2[0]->Fill(centrality_v2   );
         h_cent_v2[1]->Fill(centrality_v2,PS);
         h_FCal_Et[0]->Fill(FCal_Et/1e6     );
         h_FCal_Et[1]->Fill(FCal_Et/1e6  ,PS);
         break;
       }
       case Bins::DataType::pp2017:
       case Bins::DataType::pp2017_MB:
       case Bins::DataType::pp2017_HMT:
       {
          bool any_mb=
          b_HLT_mb_sptrk             ||
          b_HLT_noalg_mb_L1MBTS_1    ||
          b_HLT_noalg_mb_L1MBTS_1_1;

          bool any_hmt=
          b_HLT_mb_sp400_trk40_hmt_L1MBTS_1_1        ||
          b_HLT_mb_sp900_trk60_hmt_L1MBTS_1_1        ||
          b_HLT_mb_sp900_pusup400_trk60_hmt_L1TE5    ||
          b_HLT_mb_sp1000_pusup450_trk70_hmt_L1TE5   ||
          b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE10  ||
          b_HLT_mb_sp1800_pusup700_trk110_hmt_L1TE30 ||
          b_HLT_mb_sp600_trk40_hmt_L1TE5   ||
          b_HLT_mb_sp900_trk50_hmt_L1TE5   ||
          b_HLT_mb_sp900_trk60_hmt_L1TE5   ||
          b_HLT_mb_sp1000_trk70_hmt_L1TE5  ||
          b_HLT_mb_sp1400_trk90_hmt_L1TE5  ||
          b_HLT_mb_sp600_trk40_hmt_L1TE10  ||
          b_HLT_mb_sp700_trk50_hmt_L1TE10  ||
          b_HLT_mb_sp900_trk60_hmt_L1TE10  ||
          b_HLT_mb_sp1100_trk70_hmt_L1TE10 ||
          b_HLT_mb_sp1400_trk90_hmt_L1TE10 ||
          b_HLT_mb_sp700_trk50_hmt_L1TE20  ||
          b_HLT_mb_sp900_trk60_hmt_L1TE20  ||
          b_HLT_mb_sp1100_trk70_hmt_L1TE20 ||
          b_HLT_mb_sp1200_trk80_hmt_L1TE20 ||
          b_HLT_mb_sp1400_trk90_hmt_L1TE20 ||
          b_HLT_mb_sp1100_trk70_hmt_L1TE30 ||
          b_HLT_mb_sp1200_trk80_hmt_L1TE30 ||
          b_HLT_mb_sp1400_trk90_hmt_L1TE30 ||
          b_HLT_mb_sp1400_trk90_hmt_L1TE40 ||
          b_HLT_mb_sp1600_trk100_hmt_L1TE40||
          b_HLT_mb_sp1900_trk120_hmt_L1TE40||
          b_HLT_mb_sp1600_trk100_hmt_L1TE50||
          b_HLT_mb_sp1700_trk110_hmt_L1TE50||
          b_HLT_mb_sp1900_trk120_hmt_L1TE50;

          bool any_trig=any_mb || any_hmt;

          if(m_data_type==Bins::DataType::pp2017     && (!any_trig) ) continue;
          if(m_data_type==Bins::DataType::pp2017_MB  && (!any_mb  ) ) continue;
          if(m_data_type==Bins::DataType::pp2017_HMT && (!any_hmt ) ) continue;

          int ntrk=trk_numqual->at(1);
          h_ntrk->Fill(ntrk);
          break;
       }
       case Bins::DataType::pp2017_13TeV:
       case Bins::DataType::pp2017_13TeV_MB:
       case Bins::DataType::pp2017_13TeV_HMT:
       {
          bool any_mb =b_HLT_mb_sptrk;
          bool any_hmt=b_HLT_mb_sp700_pusup350_trk50_hmt_L1TE10   ||
                       b_HLT_mb_sp1100_pusup450_trk70_hmt_L1TE20  ||
                       b_HLT_mb_sp1400_pusup550_trk90_hmt_L1TE30  ||
                       b_HLT_mb_sp1600_pusup600_trk100_hmt_L1TE40 ||
                       b_HLT_mb_sp1700_pusup650_trk110_hmt_L1TE50;
          bool any_trig=any_mb || any_hmt;

          if(m_data_type==Bins::DataType::pp2017_13TeV     && (!any_trig)) continue;
          if(m_data_type==Bins::DataType::pp2017_13TeV_MB  && (!any_mb  )) continue;
          if(m_data_type==Bins::DataType::pp2017_13TeV_HMT && (!any_hmt )) continue;

          int ntrk=trk_numqual->at(1);
          h_ntrk->Fill(ntrk);
          break;
       }
       default : Common::Exception(__LINE__,__FILE__);
     }//switch 





     //--------------------------------------------------------------------------------------
     const int N=muon_pt->size();
     for(int imuon1=0;imuon1<N;imuon1++){
       const float pt1    =fabs(muon_pt ->at(imuon1)/1000.);
       const float eta1   =     muon_eta->at(imuon1);
       const int   charge1=    (muon_pt ->at(imuon1)>0)?1:-1;

       const int ieta1=GetEtaBin(eta1,charge1);
       if(ieta1<0) continue;
       const int quality=muon_quality->at(imuon1);
       if(!Bins::PassQualityCut(quality,quality,m_muon_selection_cut)) continue;

       if((quality&32 )==0) continue;//IDCuts
       if((quality&256)==0) continue;//MuonCuts
       const float d01 = fabs(muon_d0->at(imuon1));
       const float z01 = fabs(muon_z0->at(imuon1));
       const float z0sinTheta1 = z01*sin(2.0*atan(exp(-eta1)));
       if(d01>1.0 || z0sinTheta1>1.0) continue;


       //Fill histograms for efficiency calculation
       for(int itype=0;itype<TrigEffBins::N_EFF_TYPE;itype++){
       //for(int itype=0;itype<1;itype++){}
         //Check if muon is matched to Trigger
         for(int itrig=0;itrig<TrigEffBins::NTRIGS;itrig++){
           bool pass=false;
           float prescale=1;
           if(itrig==TrigEffBins::HLT_MU4){
             if(itype==PRESCALE_CORRECTED) {
               pass    =muon_b_HLT_mu4->at(imuon1);
               prescale=f_HLT_mu4_prescale;
             }
             else{
               pass    =muon_b_HLT_mu4_V2->at(imuon1);// && (b_HLT_mu4||b_HLT_mu4_rerun_decision);
             }
             if(pass==true && prescale<=0.0001){
               std::cout<<"PS="<<prescale<<" for trigger HLT_mu4"<<std::endl;
               badcount_mu4[itype]++;
               //continue;
             }
           }
           else if(itrig==TrigEffBins::HLT_MU3){
             if(itype==PRESCALE_CORRECTED) {
               pass    =(muon_b_HLT_mu3)? muon_b_HLT_mu3->at(imuon1):false;
               prescale=(muon_b_HLT_mu3)? f_HLT_mu3_prescale:1;
             }
             else{
               pass    =(muon_b_HLT_mu3_V2)? muon_b_HLT_mu3_V2->at(imuon1):false;
             }
             if(pass==true && prescale<=0.0001){
               std::cout<<"PS="<<prescale<<" for trigger HLT_mu3"<<std::endl;
               badcount_mu3[itype]++;
               //continue;
             }
           }
           if(pass          ) {
             h_Eff_Num   [itype][ieta1][itrig]->Fill(pt1);
             h_Eff_Num_2D[itype]       [itrig]->Fill(eta1*charge1,pt1);
           }
           if(prescale>0.0001) {
             h_Eff_Den   [itype][ieta1][itrig]->Fill(pt1             ,1.0/prescale);
             h_Eff_Den_2D[itype]       [itrig]->Fill(eta1*charge1,pt1,1.0/prescale);
           }
           //----------------------------------------------------------
           //New Centrality depedent effs
           if(m_data_type==Bins::DataType::PbPbAll){
             for(int ibin=0;ibin<N_CENT;ibin++){
               if(IsInCentBin(centrality_v1,ibin)){ 
                 if(pass          ) h_Eff_Num_cent [itype][ieta1][itrig][ibin]->Fill(pt1);
                 if(prescale>0.001) h_Eff_Den_cent [itype][ieta1][itrig][ibin]->Fill(pt1,1.0/prescale);
                 //eta integrated
                 if(pass          ) h_Eff_Num_cent2[itype]       [itrig][ibin]->Fill(pt1);
                 if(prescale>0.001) h_Eff_Den_cent2[itype]       [itrig][ibin]->Fill(pt1,1.0/prescale);
               }
             }
           }
           //----------------------------------------------------------
         }
       }
     }
     //--------------------------------------------------------------------------------------
   }
   std::cout<<"Badcounts mu3="<<badcount_mu3[0]<<" "<<badcount_mu3[1]
                     <<" mu4="<<badcount_mu4[0]<<" "<<badcount_mu4[1]<<std::endl;


   //TGraphs storing efficiency and Fake
   for(int itype=0;itype<TrigEffBins::N_EFF_TYPE;itype++){
     for(int ieta=0;ieta<TrigEffBins::NETA;ieta++){
       for(int itrig=0;itrig<TrigEffBins::NTRIGS;itrig++){
         if(itype==TrigEffBins::PRESCALE_CORRECTED){
           //For Prescale corrected case, the denomenator might be smaller than numerator
           gr_Eff [itype][ieta][itrig]->Divide(h_Eff_Num[itype][ieta][itrig],h_Eff_Den[itype][ieta][itrig],"pois");
         }
         else{
           gr_Eff [itype][ieta][itrig]->BayesDivide(h_Eff_Num[itype][ieta][itrig],h_Eff_Den[itype][ieta][itrig]);
         }
       }
     }
   }
}


void TrigEff::MakeEffFits()
{
   char name[600];
   for(int itype=0;itype<TrigEffBins::N_EFF_TYPE;itype++){
     for(int ieta=0;ieta<TrigEffBins::NETA;ieta++){
       for(int itrig=0;itrig<TrigEffBins::NTRIGS;itrig++){
         sprintf(name,"tf1_eff_fit_type%d_eta%d_trig%d",itype,ieta,itrig);
         if     (itrig==TrigEffBins::HLT_MU4){
           tf1_eff_fit[itype][ieta][itrig]= new TF1(name,"[0]/(1+exp(-(x-[1])/[2])) +[3]*(x-4)*(x>4)",4,20);
           tf1_eff_fit[itype][ieta][itrig]->SetParameter(0,0.9);
           tf1_eff_fit[itype][ieta][itrig]->SetParameter(1,2.0);
           tf1_eff_fit[itype][ieta][itrig]->SetParameter(2,0.2);
           tf1_eff_fit[itype][ieta][itrig]->SetParameter(3,0.02);
         }
         else if(itrig==TrigEffBins::HLT_MU3){
           tf1_eff_fit[itype][ieta][itrig]= new TF1(name,"[0]/(1+exp(-(x-[1])/[2])) +[3]*(x-3.5)*(x>3.5)",3.5,20);
           tf1_eff_fit[itype][ieta][itrig]->SetParameter(0,0.9);
           tf1_eff_fit[itype][ieta][itrig]->SetParameter(1,2.0);
           tf1_eff_fit[itype][ieta][itrig]->SetParameter(2,0.2);
           tf1_eff_fit[itype][ieta][itrig]->SetParameter(3,0.02);
         }
         gr_Eff     [itype][ieta][itrig]->Fit(tf1_eff_fit[itype][ieta][itrig],"RI");
         gr_Eff     [itype][ieta][itrig]->GetListOfFunctions()->Clear();
       }
     }
   }
}


int TrigEff::GetEtaBin(float eta,float charge)
{
   if(fabs(eta)>2.4) return -1;
   //if(fabs(eta)<0.1) return -1;//TODO maybe this cut should be imposed

   if(charge<0) eta=-eta;

   if(eta<-2.0) return 0;
   if(eta<-1.5) return 1;
   if(eta<-1.0) return 2;
   if(eta<-0.5) return 3;
   if(eta< 0.5) return 4;
   if(eta< 1.0) return 5;
   if(eta< 1.5) return 6;
   if(eta< 2.0) return 7;
   if(eta< 2.4) return 8;
   return -1;
}


float TrigEff::GetTriggerEfficiency(float eta, float pt, float charge,int itrig, int cent_percentile)
{
  int ieta=GetEtaBin(eta,charge);
  if(ieta<0){
    cout<<__PRETTY_FUNCTION__<<": exception at line "<<__LINE__<<":  eta="<<eta<<std::endl;
    throw std::exception();
  }
  if(!tf1_eff_fit[TrigEffBins::RERUN][ieta][itrig]){
    cout<<__PRETTY_FUNCTION__<<": exception at line "<<__LINE__<<":  function is not initialized"<<std::endl;
    throw std::exception();
  }

  float scale=1.045;//>60%
  if     (cent_percentile<10) scale=0.950;
  else if(cent_percentile<20) scale=1.000;
  else if(cent_percentile<30) scale=1.021;
  else if(cent_percentile<40) scale=1.029;
  else if(cent_percentile<50) scale=1.037;
  else if(cent_percentile<60) scale=1.041;


  //TODO : currently the fits are not good in the turn-on region
  //but the stat erros in the efficiency graphs are small
  //so we use the fits above 8 GeV and below 8 GeV we use the Tgraphs themselves
  if(pt>8) return tf1_eff_fit[TrigEffBins::RERUN][ieta][itrig]->Eval(pt)*scale;
  else     return gr_Eff     [TrigEffBins::RERUN][ieta][itrig]->Eval(pt)*scale;//TODO fix fits
}


void TrigEff::Plot(){
  PlotEff(0);
  PlotEff(1);
  PlotEff(2);
  Plot2DEff(TrigEffBins::RERUN,TrigEffBins::HLT_MU4);
  Plot2DEff(TrigEffBins::RERUN,TrigEffBins::HLT_MU3);
}


void TrigEff::PlotEff(int option)
{
  TCanvas *C1=Common::StandardCanvas9("Can_eff_trig_option"+std::to_string(option));
  m_can_vec.push_back(C1);
  for(int ieta=0;ieta<TrigEffBins::NETA;ieta++){
    C1->cd(ieta+1);

    TH1* hist=new TH1D(Common::UniqueName().c_str(),"",50,0,20);
    hist->SetLineColor(0);
    hist->Draw();

    if(option==0){//Mu4 Eff and fit only
      int itype=TrigEffBins::RERUN,itrig=TrigEffBins::HLT_MU4;
      TGraphAsymmErrors *gr   =(TGraphAsymmErrors*)gr_Eff     [itype][ieta][itrig]->Clone(Common::UniqueName().c_str());
      TF1               *func =(TF1*              )tf1_eff_fit[itype][ieta][itrig]->Clone(Common::UniqueName().c_str());
      gr  ->Draw("P");
      func->Draw("same");
      Common::format(gr ,1,24);
      Common::myMarkerText(0.4,0.30,1,24," "+LabelTrig[itrig],1.0,0.050);
    }
    else if(option==1){//Mu3 and Mu4 Eff
      int itype=TrigEffBins::RERUN,itrig1=TrigEffBins::HLT_MU3,itrig2=TrigEffBins::HLT_MU4;
      TGraphAsymmErrors *gr   =(TGraphAsymmErrors*)gr_Eff     [itype][ieta][itrig1]->Clone(Common::UniqueName().c_str());
      TGraphAsymmErrors *gr2  =(TGraphAsymmErrors*)gr_Eff     [itype][ieta][itrig2]->Clone(Common::UniqueName().c_str());
      TF1               *func =(TF1*              )tf1_eff_fit[itype][ieta][itrig1]->Clone(Common::UniqueName().c_str());
      TF1               *func2=(TF1*              )tf1_eff_fit[itype][ieta][itrig2]->Clone(Common::UniqueName().c_str());
      Common::format(gr ,1,24);
      func->SetLineColor(1);
      Common::myMarkerText(0.4,0.30,1,24," "+LabelTrig[itrig1],1.2,0.040);
      Common::format(gr2,2,25);
      func2->SetLineColor(2);
      Common::myMarkerText(0.4,0.23,2,25," "+LabelTrig[itrig2],1.2,0.040);
      gr   ->Draw("P");
      gr2  ->Draw("P");
      func ->Draw("same");
      func2->Draw("same");
    }
    else if(option==2){//Mu3 Minbias vs Prescale_Corrected
      int itype1=TrigEffBins::RERUN,itype2=TrigEffBins::PRESCALE_CORRECTED,itrig=TrigEffBins::HLT_MU3;
      TGraphAsymmErrors *gr  =(TGraphAsymmErrors*)gr_Eff[itype1][ieta][itrig]->Clone(Common::UniqueName().c_str());
      TGraphAsymmErrors *gr2 =(TGraphAsymmErrors*)gr_Eff[itype2][ieta][itrig]->Clone(Common::UniqueName().c_str());
      Common::format(gr ,1,24);
      Common::myMarkerText(0.4,0.30,1,24," "+LabelTrig[itrig]+" "+LabelType[itype1],1.2,0.040);
      Common::format(gr2,2,25);
      Common::myMarkerText(0.4,0.23,2,25," "+LabelTrig[itrig]+" "+LabelType[itype2],1.2,0.040);
      gr ->Draw("P");
      gr2->Draw("P");
    }
    gPad->SetGrid();

    Common::FormatHist(hist,Common::StandardFormat());
    Bins::LabelATLAS2(.2, .92, 15, .17, .12, .07);
    hist->GetYaxis()->SetRangeUser(0.001, 1.62);
    hist->GetXaxis()->SetRangeUser(4.0, 20);
    hist->GetYaxis()->SetTitle("Trigger Efficiency");
    hist->GetXaxis()->SetTitle("#it{p}_{T}");

    Common::myText2(.2, .78, 1, LabelEta[ieta], 15, 43);
  }
}

void TrigEff::Plot2DEff(int itype,int itrig){
  //gStyle->SetPalette(kBird);
  char name[600];
  sprintf(name,"Can_trigeff_type%d_trig%d",itype,itrig);
  TCanvas *C1=Common::StandardCanvas1(name);
  C1->SetRightMargin(0.2);
  m_can_vec.push_back(C1);

  TH2D *heff  =(TH2D*)h_Eff_Num_2D[itype][itrig]->Clone(Common::UniqueName().c_str());
  Common::FormatHist(heff,Common::StandardFormat());
  heff->Divide(h_Eff_Den_2D[itype][itrig]);
  heff->GetZaxis()->SetTitle("Trigger Efficiency");
  heff->GetYaxis()->SetRangeUser(4,20);
  heff->Draw("COLZ");

  Bins::LabelATLAS2(.2,.90,15,.17,.12,.06);
  Common::myText2(.2, .78, 1,LabelTrig[itrig], 15, 43);
  Common::myText2(.2, .72, 1,Bins::QualityCutLabel[m_muon_selection_cut], 15, 43);
}




//Compare the Trigger Effs between different datasets
void PlotEffCompare(int quality_cut=Bins::QualityCut::MEDIUM,int icase=0,int itype=TrigEff::RERUN){
  
  vector<int>  data_types={Bins::DataType::PbPbAll,Bins::DataType::pp2017,Bins::DataType::pp2017_13TeV};
  if(icase==1) data_types={Bins::DataType::pp2017,Bins::DataType::pp2017_MB,Bins::DataType::pp2017_HMT};
  if(icase==2) data_types={Bins::DataType::pp2017_13TeV,Bins::DataType::pp2017_13TeV_MB,Bins::DataType::pp2017_13TeV_HMT};

  //const int itype = TrigEff::RERUN  ;
  const int itrig = TrigEff::HLT_MU4;

  //gStyle->SetPalette(kBird);
  char name[600];
  sprintf(name,"Can_trigeff_type%d_trig%d_qual%d",itype,itrig,quality_cut);
  TCanvas *C1=Common::StandardCanvas9(name);
  //m_can_vec.push_back(C1);


  //--------------------------------------
  //Medium muon Trigeffs from Qipeng
  TFile *file0  = new TFile("EffFiles/TriggerEfficiency_pp.root","read");
  if(file0->IsZombie()) Common::Exception(__LINE__,__FILE__);
  TH2D* hist   =(TH2D*)file0->Get("trigeff2017_mu4_MC_StatError");
  TH2D* histSF =(TH2D*)file0->Get("trigscalefactor2017_mu4_StatError");
  for(int ibinx=1;ibinx<=hist->GetNbinsX();ibinx++){
    for(int ibiny=1;ibiny<=hist->GetNbinsY();ibiny++){
      float q_eta=hist->GetXaxis()->GetBinCenter(ibinx);
      float pt   =hist->GetYaxis()->GetBinCenter(ibiny);

      int _ibinx  =histSF->GetXaxis()->FindFixBin(fabs(q_eta));
      int _ibiny  =histSF->GetYaxis()->FindFixBin(pt         );
      float sf   =histSF->GetBinContent(_ibinx,_ibiny);

      float val=hist->GetBinContent(ibinx,ibiny);
      hist->SetBinContent(ibinx,ibiny,val*sf);
    }
  }
  //--------------------------------------


  
  int idraw=0;
  for(int idata:data_types){

    TrigEff *trigEffs=new TrigEff(true,quality_cut,idata);
    TH2D *hNum  =(TH2D*)trigEffs->h_Eff_Num_2D[itype][itrig]->Clone(Common::UniqueName().c_str());
    TH2D *hDen  =(TH2D*)trigEffs->h_Eff_Den_2D[itype][itrig]->Clone(Common::UniqueName().c_str());

    for(int ipad=0;ipad<9;ipad++){
      C1->cd(ipad+1);

      const double etabins[]={-2.4,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.4};
      hNum->GetXaxis()->SetRangeUser(etabins[ipad]+0.00001,etabins[ipad+1]-.00001);
      hDen->GetXaxis()->SetRangeUser(etabins[ipad]+0.00001,etabins[ipad+1]-.00001);

      TH1D* h1=hNum->ProjectionY(Common::UniqueName().c_str());
      TH1D* h2=hDen->ProjectionY(Common::UniqueName().c_str());
      h1->Divide(h2);

      int col[]={1,2,4};
      int sty[]={20,24,25};
      Common::format(h1,col[idraw],sty[idraw]);
      
      if(idraw==0){
        h1->Draw();
        h1->GetYaxis()->SetRangeUser(0,1.2);
        h1->GetXaxis()->SetRangeUser(4,10 );
        Bins::LabelATLAS2(.2,.90,15,.17,.12,.06);
        Common::myText2(.2, .78, 1,trigEffs->LabelTrig[itrig], 15, 43);
        Common::myText2(.2, .72, 1,Bins::QualityCutLabel[quality_cut], 15, 43);
        
        hist->GetXaxis()->SetRangeUser(etabins[ipad]+0.00001,etabins[ipad+1]-.00001);
        TH1D* htemp=hist->ProjectionY(Common::UniqueName().c_str());
        htemp->Scale( 0.1/(etabins[ipad+1]-etabins[ipad]) );
        Common::SetYError(htemp,0);//Remove Y-errors
        Common::format(htemp,kGreen-2,28);
        htemp->Draw("same");

      }
      else h1->Draw("same");
    }
    idraw++;
  }
}

//option==0 draw the % change in eff
//option==1 draw the absolute change in eff
void TrigEff::PlotEffCentDep(double ptmin, double ptmax,int option)
{
  int itype=TrigEffBins::RERUN;
  int itrig=TrigEffBins::HLT_MU4;
  Double_t bins[2]={ptmin,ptmax};

  char name[600];
  TGraphAsymmErrors *gr_Eff=new TGraphAsymmErrors(1);
  //const int Nbins=3;double centbins[Nbins+1]={0,20,40,100};
  const int Nbins=6;double centbins[Nbins+1]={0,10,20,30,40,60,100};

  //Eta differential
  sprintf(name,"Can_eff_cent_trig%d_ptmin%d_ptmax%d_option%d",itrig,int(ptmin*10),int(ptmax*10),option);
  TCanvas *C1=Common::StandardCanvas9(name);
  m_can_vec.push_back(C1);
  for(int ieta=0;ieta<TrigEffBins::NETA;ieta++){
    C1->cd(ieta+1);
    TH1* hist=new TH1D(Common::UniqueName().c_str(),"",Nbins,centbins);

    //the centrality averaged hist
    TH1D* hist_num=(TH1D*) h_Eff_Num  [itype][ieta][itrig]->Rebin(1,Common::UniqueName().c_str(),bins);
    TH1D* hist_den=(TH1D*) h_Eff_Den  [itype][ieta][itrig]->Rebin(1,Common::UniqueName().c_str(),bins);
    hist_num->Divide(hist_den);
    double avg_eff=hist_num->GetBinContent(1);
    delete hist_num;hist_num=nullptr;
    delete hist_den;hist_den=nullptr;

    
    for(int ibin=0;ibin<Nbins;ibin++){
      int icent=GetCentIndexTrig(centbins[ibin],centbins[ibin+1]);
      hist_num=(TH1D*) h_Eff_Num_cent  [itype][ieta][itrig][icent]->Rebin(1,Common::UniqueName().c_str(),bins);
      hist_den=(TH1D*) h_Eff_Den_cent  [itype][ieta][itrig][icent]->Rebin(1,Common::UniqueName().c_str(),bins);

      gr_Eff->BayesDivide(hist_num,hist_den);
      Double_t x=0,y=0,y_err=0;
      gr_Eff->GetPoint(0,x,y);
      y_err=gr_Eff->GetErrorYhigh(0);

      hist->SetBinContent(ibin+1,y);
      hist->SetBinError  (ibin+1,y_err);
      //cout<<ieta<<" "<<icent<<" "<<y<<" "<<y_err<<"  "<<gr_Eff->GetN()<<endl;

      delete hist_num;
      delete hist_den;
    }
    cout<<endl;
    if(option==0) hist->Scale(1.0/avg_eff);
    hist->Draw();
    
    gPad->SetGrid();

    Common::FormatHist(hist,Common::StandardFormat());
    Bins::LabelATLAS2(.2, .92, 15, .17, .12, .07);
    hist->GetYaxis()->SetRangeUser(0.001, 1.62);
    hist->GetYaxis()->SetTitle("Efficiency");
    if(option==0) hist->GetYaxis()->SetTitle("Relative Change");
    hist->GetXaxis()->SetTitle("Centrality [%]");

    Common::myText2(.2, .78, 1, LabelEta[ieta], 15, 43);
  }


  //Eta Integrated
  sprintf(name,"Can_eff_cent2_trig%d_ptmin%d_ptmax%d_option%d",itrig,int(ptmin*10),int(ptmax*10),option);
  TCanvas *C2=Common::StandardCanvas1(name);
  m_can_vec.push_back(C2);
  {
    C2->cd();
    TH1* hist=new TH1D(Common::UniqueName().c_str(),"",Nbins,centbins);

    //the centrality averaged hist
    int icent=GetCentIndexTrig(0,100);
    TH1D* hist_num=(TH1D*) h_Eff_Num_cent2  [itype][itrig][icent]->Rebin(1,Common::UniqueName().c_str(),bins);
    TH1D* hist_den=(TH1D*) h_Eff_Den_cent2  [itype][itrig][icent]->Rebin(1,Common::UniqueName().c_str(),bins);
    hist_num->Divide(hist_den);
    double avg_eff=hist_num->GetBinContent(1);
    delete hist_num;hist_num=nullptr;
    delete hist_den;hist_den=nullptr;

    
    for(int ibin=0;ibin<Nbins;ibin++){
      int icent=GetCentIndexTrig(centbins[ibin],centbins[ibin+1]);
      hist_num=(TH1D*) h_Eff_Num_cent2  [itype]      [itrig][icent]->Rebin(1,Common::UniqueName().c_str(),bins);
      hist_den=(TH1D*) h_Eff_Den_cent2  [itype]      [itrig][icent]->Rebin(1,Common::UniqueName().c_str(),bins);

      gr_Eff->BayesDivide(hist_num,hist_den);
      Double_t x=0,y=0,y_err=0;
      gr_Eff->GetPoint(0,x,y);
      y_err=gr_Eff->GetErrorYhigh(0);

      hist->SetBinContent(ibin+1,y);
      hist->SetBinError  (ibin+1,y_err);

      delete hist_num;
      delete hist_den;
    }
    cout<<endl;
    if(option==0) hist->Scale(1.0/avg_eff);
    hist->Draw();
    
    gPad->SetGrid();

    Common::FormatHist(hist,Common::StandardFormat());
    Bins::LabelATLAS2( .2, .92, 15, .17, .12, .07);
    hist->GetYaxis()->SetRangeUser(0.001, 1.62);
    hist->GetYaxis()->SetTitle("Efficiency");
    if(option==0) hist->GetYaxis()->SetTitle("Fractional change");
    hist->GetXaxis()->SetTitle("Centrality [%]");
  }
}
#endif
