#include "common.C"

//Compare the Trigger Effs between different datasets
void CompareTrigEff(){


  //gStyle->SetPalette(kBird);
  char name[600];
  sprintf(name,"Can_trigeff");
  TCanvas *C1=Common::StandardCanvas9(name);
  //m_can_vec.push_back(C1);


  //--------------------------------------
  //Medium muon Trigeffs from Qipeng
  TFile *file0  = new TFile("TriggerEfficiency_pp.root","read");
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


  
  for(int idraw:{0,1,2,3}){
    //const int itype = TrigEff::RERUN  ;//rerun=1, ps_corrrected=0
    //const int itrig = TrigEff::HLT_MU4;

    TFile *file0=nullptr;
    TH2D *hNum;
    TH2D *hDen;
    if(idraw==0) {
      file0=new TFile("../01RootFiles/TrigEff_MuonSelectionCut0.root");
      hNum  =(TH2D*)file0->Get("h_Eff_Num_2D_type1_trig0")->Clone(Common::UniqueName().c_str());
      hDen  =(TH2D*)file0->Get("h_Eff_Den_2D_type1_trig0")->Clone(Common::UniqueName().c_str());
    }
    if(idraw==1){
      file0=new TFile("../01RootFiles/TrigEff_MuonSelectionCut1.root");
      hNum  =(TH2D*)file0->Get("h_Eff_Num_2D_type1_trig0")->Clone(Common::UniqueName().c_str());
      hDen  =(TH2D*)file0->Get("h_Eff_Den_2D_type1_trig0")->Clone(Common::UniqueName().c_str());
    } 
    if(idraw==2){
      file0=new TFile("../01RootFiles/TrigEff_MuonSelectionCut0_pp2017.root");
      hNum  =(TH2D*)file0->Get("h_Eff_Num_2D_type1_trig0")->Clone(Common::UniqueName().c_str());
      hDen  =(TH2D*)file0->Get("h_Eff_Den_2D_type1_trig0")->Clone(Common::UniqueName().c_str());
    } 
    if(idraw==3){
      file0=new TFile("../01RootFiles/TrigEff_MuonSelectionCut0_pp2017.root");\
      hNum  =(TH2D*)file0->Get("h_Eff_Num_2D_type0_trig0")->Clone(Common::UniqueName().c_str());
      hDen  =(TH2D*)file0->Get("h_Eff_Den_2D_type0_trig0")->Clone(Common::UniqueName().c_str());
    } 



    for(int ipad=0;ipad<9;ipad++){
      C1->cd(ipad+1);


      const double etabins[]={-2.4,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.4};
      hNum->GetXaxis()->SetRangeUser(etabins[ipad]+0.00001,etabins[ipad+1]-.00001);
      hDen->GetXaxis()->SetRangeUser(etabins[ipad]+0.00001,etabins[ipad+1]-.00001);

      TH1D* h1=hNum->ProjectionY(Common::UniqueName().c_str());
      TH1D* h2=hDen->ProjectionY(Common::UniqueName().c_str());
      h1->Divide(h2);

      int col[]={ 1, 2, 4, 8};
      int sty[]={20,24,25,26};
      Common::format(h1,col[idraw],sty[idraw]);
      
      if(idraw==0){
        h1->Draw();
        h1->GetYaxis()->SetRangeUser(0,1.2);
        h1->GetXaxis()->SetRangeUser(4,10 );
        
        //Qipeng Medium muon Effs
        hist->GetXaxis()->SetRangeUser(etabins[ipad]+0.00001,etabins[ipad+1]-.00001);
        TH1D* htemp=hist->ProjectionY(Common::UniqueName().c_str());
        htemp->Scale( 0.1/(etabins[ipad+1]-etabins[ipad]) );
        Common::SetYError(htemp,0);//Remove Y-errors
        Common::format(htemp,kGreen-2,28);
        htemp->Draw("same");

      }
      else h1->Draw("same");
    }
  }
}