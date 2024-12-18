#ifndef __NoOverlayMC_C__
#define __NoOverlayMC_C__

#include "NoOverlayMC.h"
#include "Bins.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TDatime.h"
#include "TRandom3.h"

#define  CREATE(VAR_NAME,TYPE,SUMW2,NAME,...) \
  if(!read_flag) {\
    VAR_NAME=new TYPE((NAME+ext).c_str(),__VA_ARGS__);\
    if(SUMW2) VAR_NAME->Sumw2();\
  }\
  else{\
    VAR_NAME=(TYPE*)File->Get((NAME+ext).c_str());\
    Common::CheckObject2(VAR_NAME,NAME+ext,File);\
  }\
  m_all_objects.push_back((TObject**)&VAR_NAME)

#define  CREATE_GR(VAR_NAME,NAME,...) \
  if(!read_flag) {\
    VAR_NAME=new TGraphAsymmErrors(__VA_ARGS__);\
    VAR_NAME->SetName((NAME+ext).c_str());\
  }\
  else{\
    VAR_NAME=(TGraphAsymmErrors*)File->Get((NAME+ext).c_str());\
    Common::CheckObject2(VAR_NAME,NAME+ext,File);\
  }\
  m_all_objects.push_back((TObject**)&VAR_NAME)

#define READ_TF1(VAR_NAME,NAME) \
  if(read_flag) {\
    VAR_NAME=(TF1*)File->Get((NAME+ext).c_str());\
    Common::CheckObject2(VAR_NAME,NAME+ext,File);\
  }\
  m_all_objects.push_back((TObject**)&VAR_NAME)



void NoOverlayMC::InitHists(bool read_flag)
{
  char name [1200];
  char name1[600];

	std::string ext="";
  CREATE(h_cent  ,TH1D,false,"h_cent"  ,";Centrality[%];Entries"  ,101 ,-0.5 ,100.5);
  CREATE(h_FCalEt,TH1D,false,"h_FCalEt",";FCal #it{E}_{T};Entries",1200,-0.5 ,5.5  );

  CREATE(h_Deta_Match ,TH1D,false,"h_Deta_Match" ,";#it{#eta}^{rec}-#it{#eta}^{gen};Entries"        ,100,-0.02,0.02);
  CREATE(h_Dphi_Match ,TH1D,false,"h_Dphi_Match" ,";#it{#phi}^{rec}-#it{#phi}^{gen};Entries"        ,100,-0.02,0.02);
  CREATE(h_DpT_Match  ,TH1D,false,"h_DpT_Match"  ,";#it{p}_{T}^{rec}-#it{p}_{T}^{gen} [GeV];Entries",100,-1   ,1   );
  CREATE(h_DpT_MatchV2,TH1D,false,"h_DpT_MatchV2",";#it{p}_{T}^{rec}-#it{p}_{T}^{gen} [GeV];Entries",100,-1   ,1   );
  CREATE(h_Prob_Match ,TH1D,false,"h_Prob_Match" ,";Truth Match Probability;Entries"                ,100, 0   ,1   );

  CREATE(h_count      ,TH1D,false,"h_count"       ,";Matched muons;"      ,5                     ,-0.5 ,4.5         );

  CREATE(p_D0_rec     ,TProfile,true,"p_D0"       ,";#it{p}_{T} [GeV];<d0> [mm];"           ,50,0,50);
  CREATE(p_Z0_rec     ,TProfile,true,"p_Z0"       ,";#it{p}_{T} [GeV];<z0sin(#theta)> [mm];",50,0,50);

  CREATE(p_Deta_Match ,TProfile,true,"p_Deta_Match",";#it{#eta}^{rec} ;#LT#it{#eta}^{rec}-#it{#eta}^{gen}#GT;"              ,50,-2.5,2.5,-.5,.5,"s");
  CREATE(p_DpT_Match  ,TProfile,true,"p_DpT_Match" ,";#it{p}_{T}^{rec} [GeV];#LT#it{p}_{T}^{rec}-#it{p}_{T}^{gen} [GeV]#GT;",32,   4,20 ,-2,2,"s");
  CREATE(h2_DpT_Match ,TH2D,false,"h2_DpT_Match" ,";#it{p}_{T}^{gen} [GeV];#it{p}_{T}^{rec}-#it{p}_{T}^{gen} [GeV];Entries",170,3,20,100,-1   ,1   );


  double xBins[101];
  double yBins[101];
  double xMaxExp = 0.03, yMaxExp = 0.03;
  for(int i = 0; i <= 100; i++)
    {
      xBins[i] = pow(10.0, ((float)(i-100))*xMaxExp );
      yBins[i] = pow(10.0, ((float)(i-100))*yMaxExp );
    }
  /*
  CREATE(h_deltaP_overP_dist_2D,TH2D,false,"h_deltaP_overP_dist_2D",";#it{#Deltap_{1}/p_{1}};#it{#Deltap_{2}/p_{2}}",100,-0.5,0.5,100,-0.5,0.5);
  CREATE(h_deltaP_overP_dist_1D,TH1D,false,"h_deltaP_overP_dist_1D",";#it{#Deltap/p};Entries",100,-0.5,0.5);
  //CREATE(h_d0_dist_2D          ,TH2D,false,"h_d0_dist_2D"          ,";#it{d0_{1}};#it{d0_{2}}",100,xBins,100,yBins);
  h_d0_dist_2D          = new TH2D("h_d0_dist_2D",";#it{d0_{1}};#it{d0_{2}}",100,xBins,100,yBins);
  CREATE(h_d0_dist_1D          ,TH1D,false,"h_d0_dist_1D"          ,";#it{d0};Entries"       ,200,-1,1);
  CREATE(h_d0_pair_dist_1D     ,TH1D,false,"h_d0_pair_dist_1D"     ,";#it{d0_{pair}};Entries",100, 0,1);
  */

  for (int ich:{Bins::PairSignCut::SAME_SIGN,
                Bins::PairSignCut::OPPOSITE_SIGN,
                Bins::PairSignCut::COMBINED_SIGN}){
    for(int idphi:{DPhi::NearSide,DPhi::AwaySide,DPhi::All}){
      sprintf(name1,"ch%d_dphi%d",ich,idphi);
      
      sprintf(name,"h_deltaP_overP_dist_2D_%s",name1);
      h_deltaP_overP_dist_2D[ich][idphi]= new TH2D(name,";#it{#Deltap_{1}/p_{1}};#it{#Deltap_{2}/p_{2}}",100,-0.5,0.5,100,-0.5,0.5);
      
      sprintf(name,"h_deltaP_overP_dist_1D_%s",name1);
      h_deltaP_overP_dist_1D[ich][idphi]= new TH1D(name,";#it{#Deltap/p}",100,-0.5,0.5);
      
      sprintf(name,"h_d0_dist_2D_%s",name1);
      h_d0_dist_2D[ich][idphi]          = new TH2D(name,";#it{d0_{1}};#it{d0_{2}}",100,xBins,100,yBins);
      
      sprintf(name,"h_d0_dist_1D_%s",name1);
      h_d0_dist_1D[ich][idphi]          = new TH1D(name,";#it{d0}",100,-1,1);
      
      sprintf(name,"h_d0_pair_dist_1D_%s",name1);
      h_d0_pair_dist_1D[ich][idphi]     = new TH1D(name,";#it{d0_{pair}}",50,0,1);

      sprintf(name,"h_deltaP_overPpairsig_dist_1D_%s",name1);
      h_deltaP_overPpairsig_dist_1D[ich][idphi]= new TH1D(name,";#it{#Deltap/p pair sig};Entries",100,0,20);
    }
  }

  h_deltaP_overP_dist_1D_singles= new TH1D("h_deltaP_overP_dist_1D_singles",";#it{#Deltap/p}",100,-0.5,0.5);
  h_d0_dist_1D_singles          = new TH1D("h_d0_dist_1D_singles",";#it{d0}",100,-1,1);
  
  h_pt1_minus_pt2_1D      = new TH1D("h_pt1_minus_pt2_1D",";#it{p_{T_{1}}-p_{T_{2}}};Entries",100,-25,25);
  h_pt1_minus_pt2_rand_1D = new TH1D("h_pt1_minus_pt2_rand_1D",";#it{p_{T_{1}}-p_{T_{2}}};Entries",100,-25,25);

  for (int ipt=0;ipt<4;ipt++){
    for(int ieta=0;ieta<3;ieta++){
      sprintf(name1,"pt%d_eta%d",ipt,ieta);
      
      sprintf(name ,"h_deltaP_overP_dist_1D_singles_%s",name1);
      h_deltaP_overP_dist_pt_and_eta_1D_singles[ipt][ieta] = new TH1D(name,";#it{#Deltap/p};Entries",100,-0.5,0.5);

      sprintf(name ,"h_deltaP_overP_dist_1D_%s",name1);
      h_deltaP_overP_dist_pt_and_eta_1D[ipt][ieta] = new TH1D(name,";#it{#Deltap/p};Entries",100,-0.5,0.5);
    }
  }

  double dpop_bins [    4]={0,3.5,7,50};
  double deta_bins [ 28+1]={0};
  double dphi_bins [128+1]={0};
  for(int i = 0;i<  20;i++) deta_bins[    i]=   0.05*i;
  for(int i = 0;i<   9;i++) deta_bins[20 +i]=1.0+0.5*i;
  for(int i = 0;i< 129;i++) dphi_bins[    i]=-Common::PI/2 + i*(2*Common::PI/128.0);
  
  for (int i=0;i<Bins::NCENT;i++){
    for(int ich=0;ich<Bins::PairSignCut::NCH;ich++){
      sprintf(name1,"cent_%dto%d_ch%d",Bins::CENT_LO[i], Bins::CENT_HI[i],ich);
      sprintf(name,"corr_sig_signif3D_%s",name1);

      if(!read_flag) corr_sig_signif3D[i][ich] = new TH3D(name,";#it{m}_{#mu#mu} [GeV];#Delta#eta;#Delta#phi",3,dpop_bins,28,deta_bins,128,dphi_bins);
      else           corr_sig_signif3D[i][ich] =(TH3D*)File->Get(name);
    }
  }
  
  ext="";
  CREATE(h_Spectra_Reco_EffCorrected,TH1D,false,"h_Spectra_Reco_EffCorrected",";#it{p}_{T} [GeV];Muons/GeV",30,0,30);
  for(int itype=0;itype<TruthType::NTYPE;itype++){
    ext="_type"+std::to_string(itype);
    CREATE(h_Spectra[itype],TH1D,false,"h_Spectra",";#it{p}_{T} [GeV];Muons/GeV",30,0,30);
    CREATE(h_eta    [itype],TH1D,false,"h_eta"    ,";#it{#eta};Entries",20,-2.5,2.5);
    CREATE(h_phi    [itype],TH1D,false,"h_phi"    ,";#it{#phi};Entries",32,-Common::PI,Common::PI);
  }

  const int NBINSPT=20;
  double eff_pt_bins[NBINSPT+1]={3.5,3.6,3.7,3.8,3.9,4,4.1,4.25,4.5,5,5.5,6,7,8,9,10,12,14,16,18,20};
  ext="";
  CREATE(h_Eff_Den_2D    ,TH2D,true,"h_Eff_Den_2D"    ,";q*#eta;#it{p}_{T} [GeV];",48,-2.4,2.4,NBINSPT,eff_pt_bins);
  CREATE(h_Eff_Num_2D    ,TH2D,true,"h_Eff_Num_2D"    ,";q*#eta;#it{p}_{T} [GeV];",48,-2.4,2.4,NBINSPT,eff_pt_bins);
  CREATE(h_Eff_Num_2D_rec,TH2D,true,"h_Eff_Num_2D_rec",";q*#eta;#it{p}_{T} [GeV];",48,-2.4,2.4,NBINSPT,eff_pt_bins);
  for(int ieta=0;ieta<TruthType::NETA;ieta++){
    ext="_eta"+std::to_string(ieta);
    CREATE   (h_Eff_Num    [ieta],TH1D,false,"h_Eff_Num"  ,";#it{p}_{T}^{gen} [GeV];",NBINSPT,eff_pt_bins);
    CREATE   (h_Eff_Den    [ieta],TH1D,false,"h_Eff_Den"  ,";#it{p}_{T}^{gen} [GeV];",NBINSPT,eff_pt_bins);
    CREATE   (h_Fake_Num   [ieta],TH1D,false,"h_Fake_Num" ,";#it{p}_{T}^{rec} [GeV];"  ,NBINSPT,eff_pt_bins);
    CREATE   (h_Fake_Den   [ieta],TH1D,false,"h_Fake_Den" ,";#it{p}_{T}^{rec} [GeV];"  ,NBINSPT,eff_pt_bins);
    CREATE_GR(gr_Eff       [ieta],"gr_Eff"     ,NBINSPT);
    CREATE_GR(gr_Fake      [ieta],"gr_Fake"    ,NBINSPT);
    READ_TF1 (tf1_eff_fit  [ieta],"tf1_eff_fit");
    CREATE   (h_Eff_Num_rec[ieta],TH1D,false,"h_Eff_Num_rec"  ,";#it{p}_{T}^{gen} [GeV];",NBINSPT,eff_pt_bins);
    CREATE_GR(gr_Eff_rec   [ieta],"gr_Eff_rec" ,NBINSPT);
  }
  for(int icent=0;icent<Bins::NCENT;icent++){
    ext="_cent"+std::to_string(icent);
    CREATE(h_Eff_Den_2D_cent    [icent],TH2D,true,"h_Eff_Den_2D"    ,";q*#eta;#it{p}_{T} [GeV];",24,-2.4,2.4,NBINSPT,eff_pt_bins);
    CREATE(h_Eff_Num_2D_cent    [icent],TH2D,true,"h_Eff_Num_2D"    ,";q*#eta;#it{p}_{T} [GeV];",24,-2.4,2.4,NBINSPT,eff_pt_bins);
    CREATE(h_Eff_Num_2D_cent_rec[icent],TH2D,true,"h_Eff_Num_2D_rec",";q*#eta;#it{p}_{T} [GeV];",24,-2.4,2.4,NBINSPT,eff_pt_bins);
    for(int ieta=0;ieta<TruthType::NETA;ieta++){
      ext="_cent"+std::to_string(icent)+"_eta"+std::to_string(ieta);
      CREATE   (h_Eff_Den_cent    [icent][ieta],TH1D,false,"h_Eff_Den"  ,";#it{p}_{T}^{gen} [GeV];",NBINSPT,eff_pt_bins);
      CREATE   (h_Eff_Num_cent    [icent][ieta],TH1D,false,"h_Eff_Num"  ,";#it{p}_{T}^{gen} [GeV];",NBINSPT,eff_pt_bins);
      CREATE_GR(gr_Eff_cent       [icent][ieta],"gr_Eff"     ,NBINSPT);
      READ_TF1 (tf1_eff_fit_cent  [icent][ieta],"tf1_eff_fit");
      CREATE   (h_Eff_Num_cent_rec[icent][ieta],TH1D,false,"h_Eff_Num_rec",";#it{p}_{T}^{gen} [GeV];",NBINSPT,eff_pt_bins);
      CREATE_GR(gr_Eff_cent_rec   [icent][ieta],"gr_Eff_rec" ,NBINSPT);
    }
  }
}
#undef CREATE
#undef CREATE_GR
#undef READ_TF1

NoOverlayMC::NoOverlayMC(int read_flag,int muon_selection_cut)
{
   //---------------------------------------------------------------
   std::string FileName="01RootFiles/MC_MuonSelection"+std::to_string(muon_selection_cut) +".root";//Create output file
   //---------------------------------------------------------------

   m_muon_selection_cut=static_cast<Bins::QualityCut>(muon_selection_cut);

   //process data
   if(read_flag==0){
     Init();              //Initialize branches
     File=new TFile(FileName.c_str(),"recreate");//Create output file
     InitHists(read_flag);//Create histograms and TGraphs
     std::cout<<"InitHists() Done"<<std::endl;

     ProcessSingles();    //Process data to fill single-particle histograms
     std::cout<<"ProcessSingles() Done"<<std::endl;
     MakeEffFits   ();    //Fit efficiency graphs to get fit functions
     std::cout<<"MakeEffFits() Done"<<std::endl;
     CheckRecoEff();
     std::cout<<"CheckRecoEff() Done"<<std::endl;
     ProcessPairs();
     std::cout<<"ProcessPairs() Done"<<std::endl;
     RebinCentralityAndCharge();
     std::cout<<"RebinCentralityAndCharge() Done"<<std::endl;
     
     File->Write();       //Write output file
     for(int ieta=0;ieta<TruthType::NETA;ieta++){//Explicitly write TGraphs and TF1s
       gr_Eff     [ieta]->Write();
       gr_Fake    [ieta]->Write();
       tf1_eff_fit[ieta]->Write();
       gr_Eff_rec [ieta]->Write();
       for(int icent=0;icent<Bins::NCENT;icent++){
         gr_Eff_cent     [icent][ieta]->Write();
         tf1_eff_fit_cent[icent][ieta]->Write();
         gr_Eff_cent_rec [icent][ieta]->Write();
       }
     }
   }
   //read in histograms from already processed data file
   else{
     File=new TFile(FileName.c_str(),"read");
     if(File->IsZombie()) Common::Exception(__LINE__,__FILE__," TFile "+FileName+" doesnot exist");
     InitHists(true);

     //make plots of the efficiency
     if(read_flag>=2) MakePlots(read_flag);
   }
}


void NoOverlayMC::ProcessSingles()
{
   Long64_t nentries = fChain->GetEntries();
   for(Long64_t jentry=0; jentry<nentries;jentry++){
     int nb = fChain->GetEntry(jentry);
     if(jentry%10000==0) std::cout<<"Running "<<jentry<<" out of "<<nentries<<std::endl;
     const int centrality=Bins::GetCentralityHITight(trk_numqual->at(3));//Multiplicity based centrality for HIJING overlay

     h_cent  ->Fill(centrality);
     h_FCalEt->Fill(FCal_Et/1e6);

     //--------------------------------------------------------------------------------------
     //Truth Muons
     int count=0;
     const int NTruth=2;//truth_pt->size(); Only the first two particles are the Starlight overlay particles, barcode check is done below
     for(int itruth1=0;itruth1<NTruth;itruth1++){
       if(truth_barcode->at(itruth1) >10002 || 
          truth_barcode->at(itruth1)<=10000 || 
          fabs(truth_id->at(itruth1))!=13) Common::Exception(__LINE__,__FILE__," Incorrect Barcode");

       const float pt1    =truth_pt    ->at(itruth1)/1000.;
       const float eta1   =truth_eta   ->at(itruth1);
       const float phi1   =truth_phi   ->at(itruth1);
       const float charge1=truth_charge->at(itruth1);

       const int ieta1=GetEtaBin(eta1,charge1);
       if(ieta1<0) continue;

       //Fill pT, eta and phi distributions
       h_Spectra[TRUTH]->Fill(pt1 );
       h_eta    [TRUTH]->Fill(eta1);
       h_phi    [TRUTH]->Fill(phi1);

       //Check if truth muon is matched to reconstructed muon
       bool is_matched1=false;
       float pt1_rec=-1;
       ///*
       const int imuon=truth_muon_index->at(itruth1);
       if(imuon>=0){
         if(muon_truth_prob->at(imuon)>=0.5){
           int qual=muon_quality->at(imuon);
           if(Bins::PassQualityCut(qual,qual,m_muon_selection_cut) && ((qual&32)==32) ) {
             is_matched1=true;
             pt1_rec=fabs(muon_pt->at(imuon)/1000.0);
             h_DpT_MatchV2 ->Fill(pt1_rec-pt1);
             count++;
           }
         }
       }
       // */

      /*
       const int NRec=muon_pt->size();
       for(int imuon2=0;imuon2<NRec;imuon2++){
         const float eta2 = muon_eta    ->at(imuon2);
         const float phi2 = muon_phi    ->at(imuon2);
         const int qual2  = muon_quality->at(imuon2);
         if(!Bins::PassQualityCut(qual2,qual2,m_muon_selection_cut) || ((qual2&32)!=32) ) continue;
         float deta=fabs(eta2-eta1);
         float dphi=fabs(phi2-phi1);
         if(dphi>Common::PI) dphi=2*Common::PI-dphi;
         float dr=sqrt( deta*deta + dphi*dphi );
         if(dr<0.1) {
           is_matched1=true;
           count++;
           break;
         }
       }
       // */

       //Fill histograms for efficiency calculation
       float qeta1=(charge1>0)? eta1:-eta1;
       h_Eff_Den         [ieta1]->Fill(pt1);
       h_Eff_Den_2D             ->Fill(qeta1,pt1);
       if(is_matched1){
         h_Eff_Num       [ieta1]->Fill(pt1);
         h_Eff_Num_2D           ->Fill(qeta1,pt1);
         h_Eff_Num_rec   [ieta1]->Fill(pt1_rec);
         h_Eff_Num_2D_rec       ->Fill(qeta1,pt1_rec);
       }
       int icent=Bins::GetCentBin(centrality);
       h_Eff_Den_cent   [icent][ieta1]->Fill(pt1);
       h_Eff_Den_2D_cent[icent]       ->Fill(qeta1,pt1);
       if(is_matched1){
         h_Eff_Num_cent       [icent][ieta1]->Fill(pt1);
         h_Eff_Num_2D_cent    [icent]       ->Fill(qeta1,pt1);
         h_Eff_Num_cent_rec   [icent][ieta1]->Fill(pt1_rec);
         h_Eff_Num_2D_cent_rec[icent]       ->Fill(qeta1,pt1_rec);
       }
     }
     h_count->Fill(count);
     //--------------------------------------------------------------------------------------


     //--------------------------------------------------------------------------------------
     //Reconstructed Muons
     const int NRec=muon_pt->size();
     for(int imuon1=0;imuon1<NRec;imuon1++){
       const float pt1           =fabs(muon_pt      ->at(imuon1)/1000.);
       const float eta1          = muon_eta         ->at(imuon1);
       const float phi1          = muon_phi         ->at(imuon1);
       const float deltaP_overP1 = muon_deltaP_overP->at(imuon1);
       const float d01           = muon_d0          ->at(imuon1);
       const float z01           = muon_z0          ->at(imuon1);
       const int qual1           = muon_quality     ->at(imuon1);
       const int charge1         =(muon_pt          ->at(imuon1)>0)? 1:-1;

       if(!Bins::PassQualityCut(qual1,qual1,m_muon_selection_cut)) continue;
       if((qual1&32 )==0) continue;//IDCuts

       int ieta1=GetEtaBin(eta1,charge1);
       if(ieta1<0) continue;

       //Fill pT, eta and phi distributions
       h_Spectra[RECO]->Fill(pt1);
       h_eta    [RECO]->Fill(eta1);
       h_phi    [RECO]->Fill(phi1);

       //pt and eta binning for dp/p
       int pt1_position=0; int eta1_position=0;
       if (4<=pt1&&pt1<5) pt1_position=0;
       if (5<=pt1&&pt1<6) pt1_position=1;
       if (6<=pt1&&pt1<8) pt1_position=2;
       if (8<=pt1       ) pt1_position=3;
       if (                  fabs(eta1)<1.05) eta1_position=0;
       if (1.05<=fabs(eta1)&&fabs(eta1)<2.00) eta1_position=1;
       if (2.00<=fabs(eta1)                 ) eta1_position=2;
       
       
       //Check if reconstructed muon is matched to truth muon
       //Barcode requirement imposed to ensure only the Starlight truth matched muons are used
       //TODO:: other reconstructed muons from UE are unfortunately counted as fakes!
       const float prob   =muon_truth_prob ->at(imuon1);
       const int   itruth1=muon_truth_index->at(imuon1);
       const int   barcode=(itruth1>=0)? truth_barcode->at(itruth1) :-1;
       const int   id     =(itruth1>=0)? fabs(truth_id->at(itruth1)):-1;

       bool is_matched1 =false;
       h_Prob_Match->Fill(prob);
       if(prob>=0.5 && itruth1>=0 && barcode<=10002 && barcode>10000 && id==13){
          is_matched1=true;
          h_Deta_Match          ->Fill(eta1-truth_eta->at(itruth1));
          h_Dphi_Match          ->Fill(phi1-truth_phi->at(itruth1));
          h_DpT_Match           ->Fill(pt1 -truth_pt ->at(itruth1)/1000.0);
          p_D0_rec              ->Fill(pt1 ,d01);
          p_Z0_rec              ->Fill(pt1 ,z01*sin(sin(2.0*atan(exp(-eta1)))));

          p_Deta_Match->Fill(eta1,eta1-truth_eta->at(itruth1));
          p_DpT_Match ->Fill(pt1 ,pt1 -truth_pt ->at(itruth1)/1000.0);
          h2_DpT_Match->Fill(truth_pt ->at(itruth1)/1000.0,pt1 -truth_pt ->at(itruth1)/1000.0);

          h_deltaP_overP_dist_1D_singles->Fill(deltaP_overP1);
          h_d0_dist_1D_singles          ->Fill(d01);
          h_deltaP_overP_dist_pt_and_eta_1D_singles[pt1_position][eta1_position]->Fill(deltaP_overP1);
       }

       //Fill histograms for fake-rate calculation
                        h_Fake_Den[ieta1]->Fill(pt1);
       if(!is_matched1) h_Fake_Num[ieta1]->Fill(pt1);
     }
     //--------------------------------------------------------------------------------------
   }


   //--------------------------------------------------------------------------------------
   //Rebinned centralities
   for(int icent=0;icent<Bins::NCENT_ADD;icent++){
     int i=icent+Bins::NCENT_ORIGINAL;
     for(int icent_add=Bins::cent_add_lo[icent];icent_add<Bins::cent_add_up[icent];icent_add++){
       h_Eff_Den_2D_cent    [i]->Add(h_Eff_Den_2D_cent    [icent_add]);
       h_Eff_Num_2D_cent    [i]->Add(h_Eff_Num_2D_cent    [icent_add]);
       h_Eff_Num_2D_cent_rec[i]->Add(h_Eff_Num_2D_cent_rec[icent_add]);
       for(int ieta=0;ieta<TruthType::NETA;ieta++){
         h_Eff_Den_cent    [i][ieta]->Add(h_Eff_Den_cent    [icent_add][ieta]);
         h_Eff_Num_cent    [i][ieta]->Add(h_Eff_Num_cent    [icent_add][ieta]);
         h_Eff_Num_cent_rec[i][ieta]->Add(h_Eff_Num_cent_rec[icent_add][ieta]);
       }
     }
   }
   //--------------------------------------------------------------------------------------


   //--------------------------------------------------------------------------------------
   //TGraphs storing efficiency and Fake
   for(int ieta=0;ieta<TruthType::NETA;ieta++){
     gr_Eff [ieta]->BayesDivide(h_Eff_Num [ieta],h_Eff_Den [ieta]);
     gr_Fake[ieta]->BayesDivide(h_Fake_Num[ieta],h_Fake_Den[ieta]);
     for(int icent=0;icent<Bins::NCENT;icent++){
       gr_Eff_cent[icent][ieta]->BayesDivide(h_Eff_Num_cent[icent][ieta],h_Eff_Den_cent[icent][ieta]);
     }
   }
   //TGraphs storing efficiency+Unfolding corrections
   for(int ieta=0;ieta<TruthType::NETA;ieta++){
     TH1D* htemp=(TH1D*)h_Eff_Num_rec [ieta]->Clone(Common::UniqueName().c_str());
     htemp->Divide(h_Eff_Den [ieta]);
     for(int ibin=1;ibin<=htemp->GetNbinsX();ibin++){
       float x  =htemp->GetBinCenter (ibin);
       float val=htemp->GetBinContent(ibin);
       float err=htemp->GetBinError  (ibin);
       gr_Eff_rec [ieta]->SetPoint     (ibin-1,x  ,val);
       gr_Eff_rec [ieta]->SetPointError(ibin-1,err,err,0,0);
     }
     delete htemp;
     for(int icent=0;icent<Bins::NCENT;icent++){
       TH1D* htemp=(TH1D*)h_Eff_Num_cent_rec[icent][ieta]->Clone(Common::UniqueName().c_str());
       htemp->Divide(h_Eff_Den_cent[icent][ieta]);
       for(int ibin=1;ibin<=htemp->GetNbinsX();ibin++){
         float x  =htemp->GetBinCenter (ibin);
         float val=htemp->GetBinContent(ibin);
         float err=htemp->GetBinError  (ibin);
         gr_Eff_cent_rec[icent][ieta]->SetPoint     (ibin-1,x  ,val);
         gr_Eff_cent_rec[icent][ieta]->SetPointError(ibin-1,err,err,0,0);
       }
       delete htemp;
     }
   }
   //--------------------------------------------------------------------------------------
}


void NoOverlayMC::MakeEffFits()
{
   char name[600];
   for(int ieta=0;ieta<TruthType::NETA;ieta++){
      sprintf(name,"tf1_eff_fit_eta%d",ieta);
      tf1_eff_fit[ieta]= new TF1(name,"[0]/(1+exp(-(x-[1])/[2]))",3.5,20);
      tf1_eff_fit[ieta]->SetParameter(0,0.9);
      tf1_eff_fit[ieta]->SetParameter(1,4.0);
      tf1_eff_fit[ieta]->SetParameter(2,0.2);
      gr_Eff[ieta]->Fit(tf1_eff_fit[ieta],"QR");
      gr_Eff[ieta]->GetListOfFunctions()->Clear();
      gr_Eff[ieta]->Fit(tf1_eff_fit[ieta],"QR");
      gr_Eff[ieta]->GetListOfFunctions()->Clear();
   }

   for(int ieta=0;ieta<TruthType::NETA;ieta++){
     for(int icent=0;icent<Bins::NCENT;icent++){
       sprintf(name,"tf1_eff_fit_cent%d_eta%d",icent,ieta);
       tf1_eff_fit_cent[icent][ieta]= new TF1(name,"[0]/(1+exp(-(x-[1])/[2]))",3.5,20);
       tf1_eff_fit_cent[icent][ieta]->SetParameter(0,0.9);
       tf1_eff_fit_cent[icent][ieta]->SetParameter(1,4.0);
       tf1_eff_fit_cent[icent][ieta]->SetParameter(2,0.2);
       gr_Eff_cent[icent][ieta]->Fit(tf1_eff_fit_cent[icent][ieta],"QR");
       gr_Eff_cent[icent][ieta]->GetListOfFunctions()->Clear();
       gr_Eff_cent[icent][ieta]->Fit(tf1_eff_fit_cent[icent][ieta],"QR");
       gr_Eff_cent[icent][ieta]->GetListOfFunctions()->Clear();
     }
   }
}



void NoOverlayMC::CheckRecoEff(){
   Long64_t nentries = fChain->GetEntries();
   for(Long64_t jentry=0; jentry<nentries;jentry++){
     int nb = fChain->GetEntry(jentry);
     const int centrality=Bins::GetCentralityHITight(trk_numqual->at(3));//Multiplicity based centrality for HIJING overlay
     //--------------------------------------------------------------------------------------
     //Reconstructed Muons
     const int NRec=muon_pt->size();
     for(int imuon1=0;imuon1<NRec;imuon1++){
       const float pt1  =fabs(muon_pt ->at(imuon1)/1000.);
       const float eta1 = muon_eta    ->at(imuon1);
       const float phi1 = muon_phi    ->at(imuon1);
       const int qual1  = muon_quality->at(imuon1);
       const int charge1=(muon_pt     ->at(imuon1)>0)? 1:-1;
       if(pt1<3.8) continue;

       if(!Bins::PassQualityCut(qual1,qual1,m_muon_selection_cut)) continue;
       if((qual1&32 )==0) continue;//IDCuts

       int ieta1=GetEtaBin(eta1,charge1);
       if(ieta1<0) continue;

       //Check if reconstructed muon is matched to truth muon
       //Barcode requirement imposed to ensure only the Starlight truth matched muons are used
       //TODO:: other reconstructed muons from UE are unfortunately counted as fakes!
       const float prob   =muon_truth_prob ->at(imuon1);
       const int   itruth1=muon_truth_index->at(imuon1);
       const int   barcode=(itruth1>=0)? truth_barcode->at(itruth1) :-1;
       const int   id     =(itruth1>=0)? fabs(truth_id->at(itruth1)):-1;

       bool is_matched1 =false;
       if(prob>=0.5 && itruth1>=0 && barcode<=10002 && barcode>10000 && id==13){
          is_matched1=true;
          float  eff1=GetReconstructionEfficiency(eta1,pt1,charge1,centrality);
          if(eff1<=0.0001){
            std::cout<<pt1<<" "<<eta1<<" "<<charge1<<" "<<centrality<<std::endl;
            Common::Exception(__LINE__,__FILE__);
          }
          h_Spectra_Reco_EffCorrected->Fill(pt1,1.0/eff1);
       }
     }
     //--------------------------------------------------------------------------------------
   }
}


void NoOverlayMC::ProcessPairs()
{
   Long64_t nentries = fChain->GetEntries();
   for(Long64_t jentry=0; jentry<nentries;jentry++){
     int nb = fChain->GetEntry(jentry);
     if(jentry%10000==0) std::cout<<"Running "<<jentry<<" out of "<<nentries<<std::endl;
     const int centrality=Bins::GetCentralityHITight(trk_numqual->at(3));//Multiplicity based centrality for HIJING overlay

     //--------------------------------------------------------------------------------------
     //--------------------------------------------------------------------------------------
     //Reconstructed Muons
     const int NRec=muon_pt->size();
     for (int imuon1=0; imuon1<(NRec-1); imuon1++){
       for (int imuon2=imuon1+1; imuon2<NRec;imuon2++){
         const float pt1          = fabs(muon_pt     ->at(imuon1)/1000.);
         const float eta1         = muon_eta         ->at(imuon1);
         const float phi1         = muon_phi         ->at(imuon1);
         const float deltaP_overP1= muon_deltaP_overP->at(imuon1);
         const float d01          = muon_d0          ->at(imuon1);
         const float d01_abs      = fabs(d01);
      // const float d01_err      =(muon_d0_err      ->at(imuon1))/(muon_d0_err->at(imuon1));//SOUMYA: BUG 
         const float d01_err      =(muon_d0_err      ->at(imuon1));
         const float d01_sig      = d01/d01_err;
         const int qual1          = muon_quality     ->at(imuon1);
         const int charge1        =(muon_pt          ->at(imuon1)>0)? 1:-1;
         
         const float pt2          = fabs(muon_pt     ->at(imuon2)/1000.);
         const float eta2         = muon_eta         ->at(imuon2);
         const float phi2         = muon_phi         ->at(imuon2);
         const float deltaP_overP2= muon_deltaP_overP->at(imuon2);
         const float d02          = muon_d0          ->at(imuon2);
         const float d02_abs      = fabs(d02);
      // const float d02_err      =(muon_d0_err      ->at(imuon2))/(muon_d0_err->at(imuon2));//SOUMYA: BUG 
         const float d02_err      =(muon_d0_err      ->at(imuon2));
         const float d02_sig      = d02/d02_err;
         const int qual2          = muon_quality     ->at(imuon2);
         const int charge2        =(muon_pt          ->at(imuon2)>0)? 1:-1;

         const float ptdiff       =(pt1-pt2);
         const float dphi         = [phi1,phi2]{double dphi=atan2(sin(phi1-phi2),cos(phi1-phi2));
                                            //Shift dphi
                                            return (dphi <  -Common::PI/2)?dphi+2*Common::PI : dphi;
                                          }();
         const int dphi_position  = (dphi>Common::PI/2)? DPhi::AwaySide:DPhi::NearSide;
         const int ich            = ((charge1*charge2)>0)? Bins::PairSignCut::SAME_SIGN : Bins::PairSignCut::OPPOSITE_SIGN;
         const float dpop_signif1 = CalculateDpopSig(pt1,eta1,deltaP_overP1,centrality);
         const float dpop_signif2 = CalculateDpopSig(pt2,eta2,deltaP_overP2,centrality);
         
         //create pair values
         const float deta         = fabs(eta1-eta2);
         const float d0_pair      = sqrt(d01*d01 + d02*d02);
         const float d0_pair_sig  = sqrt(d01_sig*d01_sig + d02_sig*d02_sig);
         const float dpop_signif_rad = sqrt((pow(dpop_signif1,2.0))+(pow(dpop_signif2,2.0)));


         //pt and eta binning for dp/p
         int pt1_position=0; int eta1_position=0; int pt2_position=0; int eta2_position=0;
         if (4<=pt1&&pt1<5) pt1_position=0;if (4<=pt2&&pt2<5) pt2_position=0;
         if (5<=pt1&&pt1<6) pt1_position=1;if (5<=pt2&&pt2<6) pt2_position=1;
         if (6<=pt1&&pt1<8) pt1_position=2;if (6<=pt2&&pt2<8) pt2_position=2;
         if (8<=pt1       ) pt1_position=3;if (8<=pt2       ) pt2_position=3;
         if (                  fabs(eta1)<1.05) eta1_position=0;if (                  fabs(eta2)<1.05) eta2_position=0;
         if (1.05<=fabs(eta1)&&fabs(eta1)<2.00) eta1_position=1;if (1.05<=fabs(eta2)&&fabs(eta2)<2.00) eta2_position=1;
         if (2.00<=fabs(eta1)                 ) eta1_position=2;if (2.00<=fabs(eta2)                 ) eta2_position=2;
         
         //To randomize filling
         static TDatime *timegen = new TDatime();
         static double time      = timegen->GetTime();
         static TRandom3 *rangen = new TRandom3(time);
         double ran              = rangen->Uniform(0,1);
         
         if(!Bins::PassQualityCut(qual1,qual2,m_muon_selection_cut)) continue;
         if((qual1&qual2&32)==0) continue;//IDCuts
         
         int ieta1=GetEtaBin(eta1,charge1);
         if(ieta1<0) continue;
         int ieta2=GetEtaBin(eta2,charge2);
         if(ieta2<0) continue;

         //Check if reconstructed muon is matched to truth muon
         //Barcode requirement imposed to ensure only the Starlight truth matched muons are used
         //TODO:: other reconstructed muons from UE are unfortunately counted as fakes!
         const float prob1   =muon_truth_prob ->at(imuon1);
         const int   itruth1 =muon_truth_index->at(imuon1);
         const int   barcode1=(itruth1>=0)? truth_barcode->at(itruth1) :-1;
         const int   id1     =(itruth1>=0)? fabs(truth_id->at(itruth1)):-1;

         const float prob2   =muon_truth_prob ->at(imuon2);
         const int   itruth2 =muon_truth_index->at(imuon2);
         const int   barcode2=(itruth2>=0)? truth_barcode->at(itruth2) :-1;
         const int   id2     =(itruth2>=0)? fabs(truth_id->at(itruth2)):-1;
         
         if(prob1>=0.5 && itruth1>=0 && barcode1<=10002 && barcode1>10000 && id1==13 && prob2>=0.5 && itruth2>=0 && barcode2<=10002 && barcode2>10000 && id2==13){
           /*
           h_deltaP_overP_dist_2D->Fill(deltaP_overP1,deltaP_overP2);
           h_deltaP_overP_dist_1D->Fill(deltaP_overP1); h_deltaP_overP_dist_1D->Fill(deltaP_overP2);
           h_d0_dist_2D          ->Fill(d01,d02);
           h_d0_dist_1D          ->Fill(d01); h_d0_dist_1D          ->Fill(d01);
           h_d0_pair_dist_1D     ->Fill(d0_pair);
           */
           if (ran<0.5) {h_deltaP_overP_dist_2D[ich][dphi_position]->Fill(deltaP_overP1,deltaP_overP2);}
           else         {h_deltaP_overP_dist_2D[ich][dphi_position]->Fill(deltaP_overP2,deltaP_overP1);}

           h_deltaP_overP_dist_1D[ich][dphi_position]->Fill(deltaP_overP1); h_deltaP_overP_dist_1D[ich][dphi_position]->Fill(deltaP_overP2);

           //h_d0_dist_2D[ich][dphi_position]->Fill(d01_abs,d02_abs);
           if (ran<0.5) h_d0_dist_2D[ich][dphi_position]->Fill(d01_abs,d02_abs);
           else         h_d0_dist_2D[ich][dphi_position]->Fill(d02_abs,d01_abs);
           
           h_d0_dist_1D[ich][dphi_position]->Fill(d01); h_d0_dist_1D[ich][dphi_position]->Fill(d02);

           h_d0_pair_dist_1D[ich][dphi_position]->Fill(d0_pair);
           h_pt1_minus_pt2_1D->Fill(ptdiff);
           if (ran<0.5) {h_pt1_minus_pt2_rand_1D->Fill( ptdiff);}
           else         {h_pt1_minus_pt2_rand_1D->Fill(-ptdiff);}

           h_deltaP_overPpairsig_dist_1D[ich][dphi_position]->Fill(dpop_signif_rad);

           h_deltaP_overP_dist_pt_and_eta_1D[pt1_position][eta1_position]->Fill(deltaP_overP1);
           h_deltaP_overP_dist_pt_and_eta_1D[pt2_position][eta2_position]->Fill(deltaP_overP2);

           corr_sig_signif3D[Bins::GetCentBin(centrality)][ich]->Fill(dpop_signif_rad,deta,dphi);
           
         }
       }
     }
     //--------------------------------------------------------------------------------------
   }
   //--------------------------------------------------------------------------------------
}






int NoOverlayMC::GetEtaBin(float eta,float charge)
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

float NoOverlayMC::GetReconstructionEfficiency(float eta, float pt, float charge, int centrality_percentile)
{
  int ieta     =GetEtaBin(eta,charge);
  if(ieta<0){
    cout<<__PRETTY_FUNCTION__<<": exception at line "<<__LINE__<<":  eta="<<eta<<std::endl;
    Common::Exception(__LINE__, __FILE__);
  }

  int icent_bin=Bins::GetCentBin(centrality_percentile);
  //TODO for bins that are more peripheral than 80% centrality
  //we use the efficiency for 80-100% centrality to gain statistics
  //this is the proper thing todo, this is still marked as TODO to
  //keep in mind that this is being done
  if(centrality_percentile>=80){
    icent_bin=Bins::GetCentIndex(80,100);
  }
  //TODO for bins that are more central than 10% centrality in the HIJING Overlay
  //we use the efficiency for 0-10% centrality to gain statistics
  //this is the proper thing todo, this is still marked as TODO to
  //keep in mind that this is being done
  else if(centrality_percentile<10){
    icent_bin=Bins::GetCentIndex(0,10);
  }

  if(m_muon_selection_cut==Bins::QualityCut::MEDIUM){
    return tf1_eff_fit_cent[icent_bin][ieta]->Eval(pt);
  }
  //For Tight tracks we have very bad statistics/fits in this case
  //jobs using tight muons should be run without efficiency corrections
  //if efficiencies need to be used, then use the TGraphs directly
  if(m_muon_selection_cut==Bins::QualityCut::TIGHT){
    //return tf1_eff_fit_cent[icent_bin][ieta]->Eval(pt);
		if(pt>20) pt=20;//Fix added for very high-pt muons as graphs donot have statistics here 
    return gr_Eff_cent     [icent_bin][ieta]->Eval(pt);//TODO fix fits
  }
  Common::Exception(__LINE__,__FILE__);
  return 0;
}

void NoOverlayMC::RebinCentralityAndCharge(){
  for(int icent=0;icent<Bins::NCENT;icent++){
    for(int ich:{Bins::PairSignCut::SAME_SIGN,Bins::PairSignCut::OPPOSITE_SIGN}){
      corr_sig_signif3D[icent][Bins::PairSignCut::COMBINED_SIGN]->Add(corr_sig_signif3D[icent][ich]);
    }
  }
  
  for(int icent=0;icent<Bins::NCENT_ADD;icent++){
    int i=icent+Bins::NCENT_ORIGINAL;
    for(int icent_add=Bins::cent_add_lo[icent];icent_add<Bins::cent_add_up[icent];icent_add++){
      for(int ich:{Bins::PairSignCut::SAME_SIGN,Bins::PairSignCut::OPPOSITE_SIGN,Bins::PairSignCut::COMBINED_SIGN}){
        corr_sig_signif3D[i][ich]->Add(corr_sig_signif3D[icent_add][ich]);
      }
    }
  }
  
  for(int ich=0;ich<2;ich++){
    for(int idphi=0;idphi<2;idphi++){
      h_deltaP_overP_dist_2D[ich][2]  ->Add(h_deltaP_overP_dist_2D[ich][idphi]);
      h_deltaP_overP_dist_2D[2][idphi]->Add(h_deltaP_overP_dist_2D[ich][idphi]);
      h_deltaP_overP_dist_2D[2][2]    ->Add(h_deltaP_overP_dist_2D[ich][idphi]);

      h_deltaP_overP_dist_1D[ich][2]  ->Add(h_deltaP_overP_dist_1D[ich][idphi]);
      h_deltaP_overP_dist_1D[2][idphi]->Add(h_deltaP_overP_dist_1D[ich][idphi]);
      h_deltaP_overP_dist_1D[2][2]    ->Add(h_deltaP_overP_dist_1D[ich][idphi]);

      h_d0_dist_2D[ich][2]  ->Add(h_d0_dist_2D[ich][idphi]);
      h_d0_dist_2D[2][idphi]->Add(h_d0_dist_2D[ich][idphi]);
      h_d0_dist_2D[2][2]    ->Add(h_d0_dist_2D[ich][idphi]);

      h_d0_dist_1D[ich][2]  ->Add(h_d0_dist_1D[ich][idphi]);
      h_d0_dist_1D[2][idphi]->Add(h_d0_dist_1D[ich][idphi]);
      h_d0_dist_1D[2][2]    ->Add(h_d0_dist_1D[ich][idphi]);

      h_d0_pair_dist_1D[ich][2]  ->Add(h_d0_pair_dist_1D[ich][idphi]);
      h_d0_pair_dist_1D[2][idphi]->Add(h_d0_pair_dist_1D[ich][idphi]);
      h_d0_pair_dist_1D[2][2]    ->Add(h_d0_pair_dist_1D[ich][idphi]);

      h_deltaP_overPpairsig_dist_1D[ich][2]  ->Add(h_deltaP_overPpairsig_dist_1D[ich][idphi]);
      h_deltaP_overPpairsig_dist_1D[2][idphi]->Add(h_deltaP_overPpairsig_dist_1D[ich][idphi]);
      h_deltaP_overPpairsig_dist_1D[2][2]    ->Add(h_deltaP_overPpairsig_dist_1D[ich][idphi]);
    }
  }
}

TH1D* NoOverlayMC::ProjectDphiSignif(const int ich,const int icent,const double deta_min,const double deta_max,const double dpop_min,const double dpop_max){

  TH3D* hist=corr_sig_signif3D[icent][ich];
  
  Common::Reset3DAxes(hist);

  int bin1=hist->GetYaxis()->FindBin(deta_min+Bins::SmallNumber);
  int bin2=hist->GetYaxis()->FindBin(deta_max-Bins::SmallNumber);

  double diff1=fabs(hist->GetYaxis()->GetBinLowEdge(bin1)-deta_min);
  double diff2=fabs(hist->GetYaxis()->GetBinUpEdge (bin2)-deta_max);
  if(diff1>Bins::SmallNumber || diff2>Bins::SmallNumber) {
    std::cout<<diff1<<" "<<diff2<<std::endl;
    Common::Exception(__LINE__,__FILE__,"Projection range doesnot match histogram bin edges");
  }

  int bin3=hist->GetXaxis()->FindBin(dpop_min+Bins::SmallNumber);
  int bin4=hist->GetXaxis()->FindBin(dpop_max-Bins::SmallNumber);

  double diff3=fabs(hist->GetXaxis()->GetBinLowEdge(bin3)-dpop_min);
  double diff4=fabs(hist->GetXaxis()->GetBinUpEdge (bin4)-dpop_max);
  if(diff3>Bins::SmallNumber || diff4>Bins::SmallNumber) {
    std::cout<<diff3<<" "<<diff4<<std::endl;
    Common::Exception(__LINE__,__FILE__,"Projection range doesnot match histogram bin edges");
  }
  hist->GetYaxis()->SetRange(bin1,bin2);
  hist->GetXaxis()->SetRange(bin3,bin4);
  
  TH1D* ret=(TH1D*) hist->Project3D("z");
  Common::Reset3DAxes(hist);

  return ret;
}


double NoOverlayMC::CalculateDpopSig(double pt_val, double eta_val, double momimb_val, int cent_val){
  double eta_bins   [9]={-2.4,-1.5,-1.05,-0.5,-0.1,0.1,0.5,1.05,1.5};
  double pt_bins    [9]={ 4,4.5 ,5   ,5.5 ,6,8,10,12,15};
  double pt_bins2   [7]={ 4,4.75,5.25,5.75,7,9,11};
  double cent_bins2 [8]={ 0,10,20,30,40,50,60,80};
  int eta_range        =  0;
  int pt_range         =  0;
  int cent_range       =  0;

  //get eta, pt, and cent ranges
  for (int ieta=0;ieta<9;ieta++){
    if (eta_val>eta_bins[ieta]) eta_range=ieta;
  }
  for (int ipt=0;ipt<7;ipt++){
    if (pt_val>pt_bins2[ipt]) pt_range=ipt;
  }
  for (int icent=0;icent<8;icent++){
    if (cent_val>cent_bins2[icent]) cent_range=icent;
  }

  //extract dpop_mean
  double dpop_rise =dpop_mean_byetaandpt[eta_range][pt_range+1]-dpop_mean_byetaandpt[eta_range][pt_range];
  double dpop_run  =(pt_bins[pt_range+1]-pt_bins[pt_range]);
  double dpop_slope=dpop_rise/dpop_run;
  double dpop_int  =0;
  if (pt_range==0) dpop_int=dpop_mean_byetaandpt[eta_range][pt_range]-0.5*dpop_run*dpop_slope;
  else dpop_int=dpop_mean_byetaandpt[eta_range][pt_range];
  double dpop_mean =dpop_int+(pt_val-pt_bins2[pt_range])*dpop_slope;

  //extract dpop_rms
  double dpop_rms_int=dpop_slope_int_byetaandcent[eta_range][cent_range];
  double dpop_rms_slope=dpop_slope_fit[0]*(eta_val*eta_val)+dpop_slope_fit[1];
  //if (eta_range==4) dpop_rms_slope=dpop_slope_fit_n0p1to0p1;
  if (dpop_rms_slope>0) dpop_rms_slope=0;
  double dpop_rms=dpop_rms_int+dpop_rms_slope*(pt_val-4);

  double dpop_signif=(momimb_val-dpop_mean)/dpop_rms;

  return dpop_signif;
}

void NoOverlayMC::MakePlots(int flag){
  if(flag&2){
    PlotEff2D();
    for(int icent_bin=0;icent_bin<Bins::NCENT;icent_bin++) PlotEff(icent_bin);
  }
  if(flag&4){
    PlotMatch (h_Deta_Match,"can_EtaMatch"  );
    PlotMatch (h_Dphi_Match,"can_PhiMatch"  );
    PlotMatch (h_DpT_Match ,"can_pTMatch"   );
    PlotMatch (h_Prob_Match,"can_ProbMatch" );
    PlotMatch (h_count     ,"can_CountMatch");
    //PlotMatch2();
  }

  std::string label=(m_muon_selection_cut==Bins::QualityCut::MEDIUM)?"_medium":"_tight";
  Common::SaveCanvas(m_can_vec,label);
}

void NoOverlayMC::PlotEff2D(){
  //gStyle->SetPalette(kBird);
  char name[600];
  sprintf(name,"Can_recoeff2D");
  TCanvas *C1=Common::StandardCanvas1(name);
  C1->SetRightMargin(0.2);
  m_can_vec.push_back(C1);

  TH1D *heff  =(TH1D*)h_Eff_Num_2D->Clone(Common::UniqueName().c_str());
  Common::FormatHist(heff,Common::StandardFormat());
  heff->Divide(h_Eff_Den_2D);
  heff->GetZaxis()->SetTitle("Efficiency");
  heff->GetYaxis()->SetRangeUser(4,20);
  heff->Draw("COLZ");

  Bins::LabelATLAS2(.2,.90,18,.17,.12,.06);
  Common::myText2  (.2,.78, 1,Bins::QualityCutLabel[m_muon_selection_cut], 18, 43);
  Common::myText2  (.2,.72, 1, "0-100%", 18, 43);
}


void NoOverlayMC::PlotEff(int icent_bin){
  TCanvas *C1=nullptr;
  std::string name="Can_recoeff1D_cent"+std::to_string(icent_bin);
  if(icent_bin==Bins::GetCentIndex(0,100)) name="Can_recoeff1D";
  C1=Common::StandardCanvas9(name);
  m_can_vec.push_back(C1);

  for(int ieta=0;ieta<TruthType::NETA;ieta++){
    C1->cd(ieta+1);

    //Efficiency
    TGraphAsymmErrors*gr     =(TGraphAsymmErrors*)gr_Eff_cent     [icent_bin][ieta]->Clone(Common::UniqueName().c_str());
    TF1*func1                =(TF1*              )tf1_eff_fit_cent[icent_bin][ieta]->Clone(Common::UniqueName().c_str());
    //Efficiency+unfolding
    TGraphAsymmErrors *gr_rec=(TGraphAsymmErrors*)gr_Eff_cent_rec [icent_bin][ieta]->Clone(Common::UniqueName().c_str());
    
    TH1D* hist=new TH1D(Common::UniqueName().c_str(),"",100,0,20);hist->Draw();
    gr->Draw("P");
    gr->SetMarkerStyle(24);
    //Draw TF1 only for Medium muons
    if(m_muon_selection_cut==Bins::QualityCut::MEDIUM) func1->Draw("same");
    //gr_rec->SetMarkerSize (1);
    //gr_rec->SetMarkerStyle(25);
    //gr_rec->SetMarkerColor(2);
    //gr_rec->Draw("P");

    
    Common::FormatHist(hist,Common::StandardFormat());
    hist->GetYaxis()->SetRangeUser(0.01, 1.99);
    hist->GetXaxis()->SetRangeUser(4., 20);
    hist->GetYaxis()->SetTitle("Efficiency");
    hist->GetXaxis()->SetTitle("#it{p}_{T}^{truth}");
    if(ieta%3 !=0 ) hist->GetYaxis()->SetTitle("");
    if(ieta   <=6 ) hist->GetXaxis()->SetTitle("");

    float x=(ieta%3==0)? 0.2:0.10;
    Bins::LabelATLAS2(x, .90, 15, .20, .15, .07);
    Common::myText2  (x, .76, 1, "HIJING Simulation", 15, 43);
    Common::myText2  (x, .69, 1, LabelEta[ieta]+", "+Bins::QualityCutLabel[m_muon_selection_cut], 15, 43);
    Common::myText2  (x, .62, 1, Bins::LabelCent(icent_bin,-1), 15, 43);
  }
}

void NoOverlayMC::PlotMatch(TH1D*hist_,std::string name){
  TCanvas *C1=Common::StandardCanvas1(name);
  m_can_vec.push_back(C1);
  gPad->SetRightMargin(0.02);

  TH1D* hist=(TH1D*)hist_->Clone(Common::UniqueName().c_str());

  Common::FormatHist(hist,Common::StandardFormat());
  hist->Draw();

  if(name=="can_CountMatch") hist->Scale(1.0/hist_->Integral());
  if(name=="can_ProbMatch" ) C1->SetLogy();

  char name1[600],name2[600];
  double mean    =hist->GetMean();
  double mean_err=hist->GetMeanError();
  double rms     =hist->GetRMS();
  double rms_err =hist->GetRMSError();
  sprintf(name1,"mean=%.2g#pm%.2g"  ,mean,mean_err);
  Common::myText2(.2,.58,1,name1,20,43);
  sprintf(name2,"#sigma=%.2g#pm%.2g",rms ,rms_err );
  Common::myText2(.2,.52,1,name2,20,43);

  float x=.60,y=.88,dy=.06;
  if(name=="can_CountMatch") y=.70;
  Bins::LabelATLAS2(x,y     ,20,.14,.10,dy);
  Common::myText2  (x,y-2*dy,1,"HIJING Simulation",20,43);
  Common::myText2  (x,y-3*dy,1,Bins::QualityCutLabel[m_muon_selection_cut]+ "0-100%", 20, 43);
}

void NoOverlayMC::PlotMatch2(){
  std::string name0="can_EtaPtMatch";
  TCanvas *C1=Common::StandardCanvas2(name0);
  m_can_vec.push_back(C1);

  TProfile* prof1=(TProfile*)p_Deta_Match->Clone(Common::UniqueName().c_str());
  TProfile* prof2=(TProfile*)p_DpT_Match ->Clone(Common::UniqueName().c_str());
  Common::FormatHist(prof1,Common::StandardFormat());
  Common::FormatHist(prof2,Common::StandardFormat());
  C1->cd(1);prof1->Draw();
  C1->cd(2);prof2->Draw();
  prof1->GetYaxis()->SetRangeUser(-.005,.005);
  prof2->GetYaxis()->SetRangeUser(-1   ,1 );

  for(int ipad:{1,2}){
    C1->cd(ipad);
    float X=0.25;
    Bins::LabelATLAS2(X,.85,20,.14,.10,.07);
    Common::myText2  (X,.7 ,1,"HIJING Simulation ",20,43);
    Common::myText2  (X,.64,1,Bins::QualityCutLabel[m_muon_selection_cut]+ "0-100%",20,43);
  }
}

#endif
