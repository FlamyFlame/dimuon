#ifndef __TrigEffMu4NoL1_h__
#define __TrigEffMu4NoL1_h__

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>
#include "Riostream.h"
#include "vector"
#include "TH2D.h"
#include "common.C"
#include "TrigEff.C"
#include "Bins.h"
#include "set"

class TrigEffMu4NoL1
{

  private:
  enum Cuts{
    NoCuts     =0,
    DpopCut0p1 =1,
    DpopCut0p15=2,
    DpopCut0p2 =4,

    NoIDCuts   =128,

    Data2015   = 8,
    DataBoth   =16,
    DefaultCuts=NoCuts,
  };

  TChain          *fChain;   //!pointer to the analyzed TTree or TChain

  // Declaration of leaf types
  UInt_t          RunNumber;
  UInt_t          lbn;
  Int_t           centrality;

  vector<float>   *muon_pair_muon1_pt          =nullptr;
  vector<float>   *muon_pair_muon1_eta         =nullptr;
  vector<float>   *muon_pair_muon1_phi         =nullptr;
  vector<int>     *muon_pair_muon1_quality     =nullptr;
  vector<float>   *muon_pair_muon1_d0          =nullptr;
  vector<float>   *muon_pair_muon1_z0          =nullptr;
  vector<int>     *muon_pair_muon1_index       =nullptr;

  vector<float>   *muon_pair_muon2_pt          =nullptr;
  vector<float>   *muon_pair_muon2_eta         =nullptr;
  vector<float>   *muon_pair_muon2_phi         =nullptr;
  vector<int>     *muon_pair_muon2_quality     =nullptr;
  vector<float>   *muon_pair_muon2_d0          =nullptr;
  vector<float>   *muon_pair_muon2_z0          =nullptr;
  vector<int>     *muon_pair_muon2_index       =nullptr;

  vector<float>   *muon_deltaP_overP           =nullptr;

  Bool_t          b_HLT_2mu4       =false;
  Bool_t          b_HLT_mu4_mu4noL1=false;
  Bool_t          b_HLT_mu4_mu4noL1_isPrescaled=false;
  Float_t         f_HLT_mu4_mu4noL1_prescale=0;
  Bool_t          b_HLT_mu4=false;
  Bool_t          b_HLT_mu6=false;
  Bool_t          b_HLT_mu8=false;

  vector<bool>    *dimuon_b_HLT_2mu4       =nullptr;
  vector<bool>    *dimuon_b_HLT_mu4_mu4noL1=nullptr;
  vector<bool>    *dimuon_b_HLT_mu4_mu4noL1_swap=nullptr;


  vector<bool>    *muon_b_HLT_mu4=nullptr;
  vector<bool>    *muon_b_HLT_mu6=nullptr;
  vector<bool>    *muon_b_HLT_mu8=nullptr;
 

  public:

  enum{
    NTRIGS=3,
      Mu4       =0,
      Mu4NoL1   =1,
      Both      =2,
    NRefTrigs =4,
      RefMu4=0,
      RefMu6=1,
      RefMu8=2,
      RefAll=3,

  };
  map<int,string> LabelType ={{Mu4       ,"Mu4"},
                              {Mu4NoL1   ,"Mu4NoL1"},
                              {Both      ,"Both"},
                             };

  map<int,string> LabelRef ={{RefMu4   ,"RefMu4"},
                             {RefMu6   ,"RefMu6"},
                             {RefMu8   ,"RefMu8"},
                             {RefAll   ,"RefAll"},
                             };

  static const int NBINSPT=20;
  const double eff_pt_bins[NBINSPT+1]={3.5,3.6,3.7,3.8,3.9,4,4.1,4.25,4.5,5,5.5,6,7,8,9,10,12,14,16,18,20};

  static const int NETA=9;
  const double EtaMin[NETA]={-2.4,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0};
  const double EtaMax[NETA]={-2.0,-1.5,-1.0,-0.5, 0.5,1.0,1.5,2.0,2.4};

  static const int NCENT=8;
  const double CentBins[NCENT+1]={0,5,10,20,30,40,60,80,100};
  static const int NPT=4;
  const double pTbins[NPT+1]={4,5,6,8,20};

  TH2D* h_eff_Den  [NRefTrigs][NTRIGS]      ={{ nullptr }};
  TH2D* h_eff_Num  [NRefTrigs][NTRIGS]      ={{ nullptr }};
  TH2D* h_eff      [NRefTrigs][NTRIGS]      ={{ nullptr }};
  TH1D* h_eff1D    [NRefTrigs][NTRIGS][NETA]={{{nullptr}}};

  //Centrality dependence
  TH2D* h_eff_Den2 [NRefTrigs][NTRIGS]      ={{ nullptr }};
  TH2D* h_eff_Num2 [NRefTrigs][NTRIGS]      ={{ nullptr }};

  //Momentum imbalance
  TH2D* h_Dpop[NCENT]                       ={nullptr};
  TH1D *h_mon_npair,*h_mon_cent,*h_cut_flow; //New

  TFile *m_OutFile=nullptr;
  int m_quality=Bins::QualityCut::MEDIUM;
  int m_cutflags=NoCuts;
  std::vector<TCanvas*> m_can_vec;

  private:
  void Init     ();
  void InitHists(int read_flag=0);
  void Process  ();
  void EvalEffs ();
  int  GetEtaBin(float eta,float charge);
  int  GetCentBin(int centrality);

  public:
  TrigEffMu4NoL1(int read_flag=0,int quality=Bins::QualityCut::MEDIUM,int cutflags=NoCuts);

  float GetTriggerEfficiency    (int itrig, float eta, float pt, float charge, int cent_percentile);
  float CentMultiplier          (int itrig,            float pt,               int cent_percentile);
  float CentMultiplierMedium    (int itrig,            float pt,               int cent_percentile);
  float CentMultiplierSystematic(int itrig,            float pt,               int cent_percentile){
    return CentMultiplier(itrig,pt,cent_percentile) - CentMultiplierMedium(itrig,pt,cent_percentile);
  }

  void Plot();
  void Plot1(int iref);
  void PlotAll(int iref);
  void PlotRelativeEffs(int iref);
  void PlotCentDep(int iref, double ptmin, double ptmax);
  void PlotCentDepV2(int iref, double ptmin, double ptmax);
  void PlotCentDepAll(int iref);
};


void TrigEffMu4NoL1::Plot(){
  //Plot1(RefMu4);
  //Plot1(RefMu6);
  //Plot1(RefMu8);

  
  Plot1(RefAll); //2D and 1D plot for each trigger
  //PlotCentDep(RefAll, 4.0, 4.5);
  //PlotCentDep(RefAll, 4.5, 5.0);
  PlotCentDep(RefAll, 4.0, 5.0);
  PlotCentDep(RefAll, 5.0, 6.0);
  //PlotCentDep(RefAll, 5.5, 6.0);
  PlotCentDep(RefAll, 6.0, 8.0);
  //PlotCentDep(RefAll, 7.0,10.0);
  PlotCentDep(RefAll,8.0,20.0);

  PlotCentDepAll(RefAll);
  PlotAll(RefAll);//Plot all the triggers together
  PlotRelativeEffs(RefAll);//Relative eff of "Both" w.r.t either 
  

  string s="_MuonSelectionCut"+std::to_string(m_quality);
  if(m_cutflags!=DefaultCuts)  s+="_cut"+std::to_string(m_cutflags);
  Common::SaveCanvas(m_can_vec,s);
}

void TrigEffMu4NoL1::Plot1(int iref){
  for(int itrig=0;itrig<NTRIGS;itrig++){//;j<OtherCuts;j++){
    TCanvas *C1=Common::StandardCanvas1("can_eff_"+LabelType[itrig]+"_"+LabelRef[iref]);
    m_can_vec.push_back(C1);
    C1->SetRightMargin(0.2);
    TH2D* hist=(TH2D*) h_eff[iref][itrig]->Clone(Common::UniqueName().c_str());
    Common::FormatHist(hist,Common::StandardFormat());
    hist->Draw("COLZ");
    hist->GetXaxis()->SetRangeUser(-2.4,2.4);
    hist->GetYaxis()->SetRangeUser(4, 20);
    hist->GetZaxis()->SetTitle("Efficiency");

    Bins::LabelATLAS(0,.2,.90,15,.17,.06);
    Common::myText2(.2, .78, 1,LabelType[itrig], 15, 43);
    Common::myText2(.2, .72, 1,Bins::QualityCutLabel[m_quality], 15, 43);
  }


  static TrigEff* l_trigEff=new TrigEff(1,m_quality);

  for(int itrig=0;itrig<NTRIGS;itrig++){//;j<OtherCuts;j++){
    string name1="can_eff1D_"+LabelType[itrig]+"_"+LabelRef[iref];
    string name2="can_eff1D_"+LabelType[itrig]+"_"+LabelRef[iref]+"_ratio";
    TCanvas *C1=Common::StandardCanvas9d(name1);
    TCanvas *C2=Common::StandardCanvas9 (name2);
    m_can_vec.push_back(C1);
    m_can_vec.push_back(C2);
    for(int ieta=0;ieta<NETA;ieta++){
      C1->cd(ieta+1);
      TH1D* h_eff=h_eff1D[iref][itrig][ieta];
      h_eff->Draw();
      h_eff->GetXaxis()->SetRangeUser(4,20);
      h_eff->GetYaxis()->SetRangeUser(0, 1.35);
      TLegend *leg=Common::StandardLegend(.4,.2,.95,.5);
      leg->SetTextSize(0.06);
      leg->AddEntry(h_eff,(LabelType[itrig]+" (From HP)").c_str(),"l");

      Common::FormatHist(h_eff,Common::StandardFormat());
      Bins::LabelATLAS(0, .2, .92, 15, .18, .07);
      if(ieta<6) {
        h_eff->GetXaxis()->SetTitleOffset(100);
        h_eff->GetXaxis()->SetLabelOffset(100);
      }
      if(ieta%3 !=0) {
        h_eff->GetYaxis()->SetTitleOffset(100);
        h_eff->GetYaxis()->SetLabelOffset(100);
      }
      h_eff->GetYaxis()->SetTitle("Efficiency");
      h_eff->GetXaxis()->SetTitle("#it{p}_{T} [GeV]");


      if(itrig==Mu4){
        TH1D* h1=(TH1D*)l_trigEff->h_Eff_Num[TrigEff::TrigEffBins::RERUN][ieta][TrigEff::TrigEffBins::HLT_MU4]->Clone(Common::UniqueName().c_str());
        TH1D* h2=       l_trigEff->h_Eff_Den[TrigEff::TrigEffBins::RERUN][ieta][TrigEff::TrigEffBins::HLT_MU4];
        TGraphAsymmErrors *gr_Eff=new TGraphAsymmErrors(1);
        gr_Eff->BayesDivide(h1,h2);
        gr_Eff->Draw("P");
        leg->AddEntry(gr_Eff,"From Minbias","p");

        C2->cd(ieta+1);
        h1->Divide(h2);
        TH1D* h_ratio=(TH1D*) h1->Clone(Common::UniqueName().c_str());
        h_ratio->Divide(h_eff);
        h_ratio->Draw();
        Common::SetYError(h_ratio,0);
      }

      C1->cd(ieta+1);
      char name1[600];
      sprintf(name1,"%0.1f<#it{q#eta}<%0.1f",EtaMin[ieta],EtaMax[ieta]);
      Common::myText2(.2, .78, 1, name1, 15, 43);

      //C1->GetPad(ieta+1)->SetLogx();
      leg->Draw();
    }
  }
}

void TrigEffMu4NoL1::PlotAll(int iref){
  string name1="can_eff1D_All_"+LabelRef[iref];
  TCanvas *C1=Common::StandardCanvas9d(name1);
  m_can_vec.push_back(C1);

  TLegend *leg=Common::StandardLegend(.4,.2,.95,.5);
  leg->SetTextSize(0.06);
  for(int itrig=0;itrig<NTRIGS;itrig++){//;j<OtherCuts;j++){
    for(int ieta=0;ieta<NETA;ieta++){
      C1->cd(ieta+1);
      TH1D* h_eff=(TH1D*)h_eff1D[iref][itrig][ieta]->Clone(Common::UniqueName().c_str());
      if(itrig==0) h_eff->Draw();
      else         h_eff->Draw("same");
      int col[]={1,2,4};
      h_eff->SetLineColor(col[itrig]);
      h_eff->GetXaxis()->SetRangeUser(4,20);
      h_eff->GetYaxis()->SetRangeUser(0, 1.35);
      if(ieta==0) leg->AddEntry(h_eff,(LabelType[itrig]).c_str(),"l");

      Common::FormatHist(h_eff,Common::StandardFormat());
      
      if(ieta<6) {
        h_eff->GetXaxis()->SetTitleOffset(100);
        h_eff->GetXaxis()->SetLabelOffset(100);
      }
      if(ieta%3 !=0) {
        h_eff->GetYaxis()->SetTitleOffset(100);
        h_eff->GetYaxis()->SetLabelOffset(100);
      }
      h_eff->GetYaxis()->SetTitle("Efficiency");
      h_eff->GetXaxis()->SetTitle("#it{p}_{T} [GeV]");

      if(itrig==0){
        Bins::LabelATLAS(0, .2, .92, 15, .18, .07);
        char name1[600];
        sprintf(name1,"%0.1f<#it{q#eta}<%0.1f",EtaMin[ieta],EtaMax[ieta]);
        Common::myText2(.2, .78, 1, name1, 15, 43);
      }
      //C1->GetPad(ieta+1)->SetLogx();
    }
  }
  leg->Draw();
}

void TrigEffMu4NoL1::PlotRelativeEffs(int iref){
  string name2="can_Relativeff1D_"+LabelRef[iref];
  TCanvas *C2=Common::StandardCanvas9 (name2);
  m_can_vec.push_back(C2);

  TLegend *leg2=Common::StandardLegend(.2,.2,.95,.5);
  leg2->SetTextSize(0.06);
  int idraw=0;
  for(int itrig:{Mu4,Mu4NoL1}){
    int itrig2=(itrig==Mu4)? Mu4NoL1:Mu4;// Other trigger than current loop variable
    for(int ieta=0;ieta<NETA;ieta++){
      C2->cd(ieta+1);

      TH1D* h_eff     =(TH1D*)h_eff1D[iref][itrig][ieta]->Clone(Common::UniqueName().c_str());
      TH1D* h_eff_Both=(TH1D*)h_eff1D[iref][Both ][ieta]->Clone(Common::UniqueName().c_str());

      int col[]={1,2,4};
      h_eff->SetLineColor(col[itrig]);
      h_eff->GetXaxis()->SetRangeUser(4,20);
      h_eff->GetYaxis()->SetRangeUser(0, 1.35);

      if(idraw==0) h_eff->Draw();
      else         h_eff->Draw("same");

      h_eff_Both->Divide(h_eff);
      h_eff_Both->SetLineStyle(2);
      h_eff_Both->SetLineColor(col[itrig2]);//note: itrig2
      h_eff_Both->Draw("same");


      if(ieta==0){
        leg2->AddEntry(h_eff     ,(                             LabelType[itrig]).c_str(),"l");
        leg2->AddEntry(h_eff_Both,(LabelType[itrig2]+" w.r.t. "+LabelType[itrig]).c_str(),"l");
      }



      Common::FormatHist(h_eff,Common::StandardFormat());
      if(ieta<6) {
        h_eff->GetXaxis()->SetTitleOffset(100);
        h_eff->GetXaxis()->SetLabelOffset(100);
      }
      if(ieta%3 !=0) {
        h_eff->GetYaxis()->SetTitleOffset(100);
        h_eff->GetYaxis()->SetLabelOffset(100);
      }
      h_eff->GetYaxis()->SetTitle("(Relative) Efficiency");
      h_eff->GetXaxis()->SetTitle("#it{p}_{T} [GeV]");

      if(idraw==0){
        Bins::LabelATLAS(0, .2, .92, 15, .18, .07);
        char name1[600];
        sprintf(name1,"%0.1f<#it{q#eta}<%0.1f",EtaMin[ieta],EtaMax[ieta]);
        Common::myText2(.2, .78, 1, name1, 15, 43);
      }

      //C1->GetPad(ieta+1)->SetLogx();
      idraw++;
    }
  }
  leg2->Draw();
}


void TrigEffMu4NoL1::PlotCentDep(int iref, double ptmin, double ptmax){
  //for(int itrig=0;itrig<NTRIGS;itrig++)
  for(int itrig:{Mu4,Mu4NoL1,Both})
  {
    string name1="can_eff_centdep_"+LabelType[itrig]+"_"+LabelRef[iref]+std::to_string(int(ptmin*10))+"_"+std::to_string(int(ptmax*10));
    TCanvas *C1=Common::StandardCanvas3(name1);
    m_can_vec.push_back(C1);


    int bin1=h_eff_Num2[iref][itrig]->GetYaxis()->FindFixBin(ptmin+1e-4);
    int bin2=h_eff_Num2[iref][itrig]->GetYaxis()->FindFixBin(ptmax-1e-4);
    TH1D* h_eff   =(TH1D*    ) h_eff_Num2[iref][itrig]->ProjectionX(Common::UniqueName().c_str(),bin1,bin2);
    TH1D* h_den   =(TH1D*    ) h_eff_Den2[iref][itrig]->ProjectionX(Common::UniqueName().c_str(),bin1,bin2);
    TProfile* p_pT=(TProfile*) h_eff_Den2[iref][itrig]->ProfileX   (Common::UniqueName().c_str(),bin1,bin2);

    static const int NCENT2=7;
    const double CentBins2[NCENT2+1]={0,10,20,30,40,60,80,100};
    h_eff=(TH1D*    )h_eff->Rebin(NCENT2,Common::UniqueName().c_str(),CentBins2);
    h_den=(TH1D*    )h_den->Rebin(NCENT2,Common::UniqueName().c_str(),CentBins2);
    p_pT =(TProfile*)p_pT ->Rebin(NCENT2,Common::UniqueName().c_str(),CentBins2);

    C1->cd(1);
    p_pT->Draw();
    p_pT->GetYaxis()->SetTitle("#LT#it{p}_{T}#GT [GeV]");
    char name[600];
    sprintf(name,"#it{p}_{T}#in(%0.1f,%0.1f) GeV",ptmin,ptmax);
    Common::myText2(.5, .50, 1,name, 15, 43);

    //h_eff->Sumw2();
    //h_den->Sumw2();
    Common::FormatHist(h_eff,Common::StandardFormat());
    float scale=h_eff->Integral()/h_den->Integral();// average efficiency across all centralities
    h_eff->Divide(h_den);
    h_eff->GetYaxis()->SetTitle("Efficiency");
    //h_eff->GetYaxis()->SetRangeUser(0, 1.35);

    C1->cd(2);
    h_eff->DrawClone();

    C1->cd(3);
    h_eff->Scale(1.0/scale);//scale by average (centrality-integrated) efficiency
    h_eff->GetYaxis()->SetTitle("Relative change");
    h_eff->Draw();
    h_eff->GetYaxis()->SetRangeUser(.9, 1.2);

    //Draw my parameterizations, for checking    
    TH1D* h_eff_param=(TH1D*)h_eff->Clone(Common::UniqueName().c_str());
    h_eff_param->Reset();
    for(int ibin=1;ibin<=h_eff_param->GetNbinsX();ibin++){
      h_eff_param->SetBinContent(ibin,CentMultiplier( itrig, (ptmin+ptmax)/2.0, h_eff_param->GetBinCenter(ibin)));
    }
    h_eff_param->SetLineStyle(2);
    h_eff_param->SetLineColor(2);
    h_eff_param->DrawClone("same");
    h_eff_param->Scale(1.03  );h_eff_param->SetLineColor(4);h_eff_param->DrawClone("same");
    h_eff_param->Scale(1/1.06);h_eff_param->SetLineColor(4);h_eff_param->DrawClone("same");

    //Draw my parameterizations, for checking    
    TH1D* h_eff_param2=(TH1D*)h_eff->Clone(Common::UniqueName().c_str());
    h_eff_param2->Reset();
    for(int ibin=1;ibin<=h_eff_param->GetNbinsX();ibin++){
      h_eff_param2->SetBinContent(ibin,CentMultiplierMedium( itrig, (ptmin+ptmax)/2.0, h_eff_param->GetBinCenter(ibin)));
    }
    h_eff_param2->SetLineStyle(1);
    h_eff_param2->SetLineColor(6);
    //h_eff_param2->DrawClone("same");


    sprintf(name,"#it{p}_{T}#in(%0.1f,%0.1f) GeV",ptmin,ptmax);
    Bins::LabelATLAS(0, .2, .90, 15, .15, .05);
    Common::myText2 (   .2, .80, 1,LabelType[itrig]+", "+name, 15, 43);
    Common::myText2 (   .2, .75, 1,Bins::QualityCutLabel[m_quality], 15, 43);
    //Common::myText2(.5, .50, 1,name, 15, 43);
  }
}

//Same as PlotCentDep but Renormalize the centraliy dependence scale to have an asymptotic value of 1.0
void TrigEffMu4NoL1::PlotCentDepV2(int iref, double ptmin, double ptmax){
  //for(int itrig=0;itrig<NTRIGS;itrig++)
  for(int itrig:{Both,Mu4NoL1})
  {
    string name1="can_eff_centdep_"+LabelType[itrig]+"_"+LabelRef[iref]+std::to_string(int(ptmin*10))+"_"+std::to_string(int(ptmax*10));
    TCanvas *C1=Common::StandardCanvas3(name1);
    m_can_vec.push_back(C1);


    int bin1=h_eff_Num2[iref][itrig]->GetYaxis()->FindFixBin(ptmin+1e-4);
    int bin2=h_eff_Num2[iref][itrig]->GetYaxis()->FindFixBin(ptmax-1e-4);
    TH1D* h_eff   =(TH1D*    ) h_eff_Num2[iref][itrig]->ProjectionX(Common::UniqueName().c_str(),bin1,bin2);
    TH1D* h_den   =(TH1D*    ) h_eff_Den2[iref][itrig]->ProjectionX(Common::UniqueName().c_str(),bin1,bin2);
    TProfile* p_pT=(TProfile*) h_eff_Den2[iref][itrig]->ProfileX   (Common::UniqueName().c_str(),bin1,bin2);

    static const int NCENT2=7;
    const double CentBins2[NCENT2+1]={0,10,20,30,40,60,80,100};
    h_eff=(TH1D*    )h_eff->Rebin(NCENT2,Common::UniqueName().c_str(),CentBins2);
    h_den=(TH1D*    )h_den->Rebin(NCENT2,Common::UniqueName().c_str(),CentBins2);
    p_pT =(TProfile*)p_pT ->Rebin(NCENT2,Common::UniqueName().c_str(),CentBins2);

    C1->cd(1);
    p_pT->Draw();
    p_pT->GetYaxis()->SetTitle("#LT#it{p}_{T}#GT [GeV]");
    char name[600];
    sprintf(name,"#it{p}_{T}#in(%0.1f,%0.1f) GeV",ptmin,ptmax);
    Common::myText2(.5, .50, 1,name, 15, 43);

    //h_eff->Sumw2();
    //h_den->Sumw2();
    Common::FormatHist(h_eff,Common::StandardFormat());
    float scale=h_eff->Integral()/h_den->Integral();// average efficiency across all centralities
    h_eff->Divide(h_den);
    h_eff->GetYaxis()->SetTitle("Efficiency");
    //h_eff->GetYaxis()->SetRangeUser(0, 1.35);

    C1->cd(2);
    h_eff->DrawClone();

    C1->cd(3);
    h_eff->Scale(1.0/scale);//scale by average (centrality-integrated) efficiency
    h_eff->GetYaxis()->SetTitle("Relative change");
    h_eff->Draw();
    float scale2=(h_eff->GetBinContent(4)+h_eff->GetBinContent(5))/2;//30-60% bin
    h_eff->Scale(1.0/scale2);
    h_eff->GetYaxis()->SetRangeUser(.8, 1.1);

    //Draw my parameterizations, for checking    
    TH1D* h_eff_param=(TH1D*)h_eff->Clone(Common::UniqueName().c_str());
    h_eff_param->Reset();
    for(int ibin=1;ibin<=h_eff_param->GetNbinsX();ibin++){
      h_eff_param->SetBinContent(ibin,CentMultiplier( itrig, (ptmin+ptmax)/2.0, h_eff_param->GetBinCenter(ibin)));
    }
    h_eff_param->SetLineStyle(2);
    h_eff_param->SetLineColor(2);
    h_eff_param->Scale(1.0/h_eff_param->GetBinContent(5));
    h_eff_param->DrawClone("same");
    for(int i:{0,1}){
    	auto *h=(TH1D*)h_eff_param->Clone(Common::UniqueName().c_str()); 
    	//additional uncertainty for 0-10%
      if(i==0) {h->Scale(1.03);h->SetBinContent(1,h->GetBinContent(1)*1.02);}
      if(i==1) {h->Scale(0.97);h->SetBinContent(1,h->GetBinContent(1)*0.98);}
      h->SetLineColor(4);
      h->Draw("same");
    }


    sprintf(name,"#it{p}_{T}#in(%0.1f,%0.1f) GeV",ptmin,ptmax);
    Bins::LabelATLAS(0, .2, .90, 15, .15, .05);
    Common::myText2 (   .2, .80, 1,LabelType[itrig]+", "+name, 15, 43);
    Common::myText2 (   .2, .75, 1,Bins::QualityCutLabel[m_quality], 15, 43);
    //Common::myText2(.5, .50, 1,name, 15, 43);
  }
}

void TrigEffMu4NoL1::PlotCentDepAll(int iref){
  string name1="can_eff_centdepAll_"+LabelRef[iref];
  TCanvas *C1=Common::StandardCanvas1a(name1);
  m_can_vec.push_back(C1);

  TLegend *leg=Common::StandardLegend(.2,.2,.95,.5);
  leg->SetTextSize(0.06);
  int idraw=0;

  for(int itrig=0;itrig<NTRIGS;itrig++){//;j<OtherCuts;j++){
    TH1D* h_eff=(TH1D*)h_eff_Num2[iref][itrig]->ProjectionX(Common::UniqueName().c_str());
    TH1D* h_den=(TH1D*)h_eff_Den2[iref][itrig]->ProjectionX(Common::UniqueName().c_str());
    h_eff->Sumw2();
    h_den->Sumw2();
    float scale=h_eff->Integral()/h_den->Integral();//average eff over all cent
    h_eff->Divide(h_den);
    h_eff->Scale(1.0/scale);
    if(idraw==0) h_eff->Draw();
    else         h_eff->Draw("same");
    h_eff->GetYaxis()->SetRangeUser(0, 1.35);
    h_eff->GetYaxis()->SetTitle("Efficiency");
    Common::FormatHist(h_eff,Common::StandardFormat());
    int sty[]={24,25,28};
    int col[]={ 1, 2, 4};
    Common::format(h_eff,col[idraw],sty[idraw]);
    leg->AddEntry(h_eff,LabelType[itrig].c_str(),"p");
    idraw++;
  }
  leg->Draw();
  Bins::LabelATLAS(0, .2, .92, 15, .18, .07);
}



TrigEffMu4NoL1::TrigEffMu4NoL1(int read_flag,int quality,int cutflags){
  m_quality =quality;
  m_cutflags=cutflags;
  if(read_flag==0){
    Init();
    InitHists(read_flag);
    Process();
    EvalEffs();
    m_OutFile->Write();
  }
  else if(read_flag==1) {
    InitHists(read_flag);
    Plot();
  }
  else if(read_flag==2){
    InitHists(read_flag);
  }
  else if(read_flag==3){
    InitHists(1);
    PlotCentDepV2(RefAll, 4.0, 5.0);
    PlotCentDepV2(RefAll, 5.0, 6.0);
    PlotCentDepV2(RefAll, 6.0, 8.0);
    PlotCentDepV2(RefAll, 4.0,20.0);
  }
  else Common::Exception(__LINE__,__FILE__);
}

void TrigEffMu4NoL1::Init(){
  fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
  fChain->SetMakeClass(1);

  if(m_cutflags&Cuts::Data2015){
    fChain->Add("Data/ForTrigEffs/user.soumya.TrigRates.physics_HP.PbPb2015.15April2022._MYSTREAM/*.root");
  }
  else if(m_cutflags&Cuts::DataBoth){
    fChain->Add("Data/ForTrigEffs/user.soumya.TrigRates.physics_HP.PbPb2018.15April2022._MYSTREAM/*.root");
    fChain->Add("Data/ForTrigEffs/user.soumya.TrigRates.physics_HP.PbPb2015.15April2022._MYSTREAM/*.root");
  }
  else{
    fChain->Add("Data/ForTrigEffs/user.soumya.TrigRates.physics_HP.PbPb2018.15April2022._MYSTREAM/*.root");
  }

  fChain->SetBranchAddress("RunNumber"                   , &RunNumber);
  fChain->SetBranchAddress("lbn"                         , &lbn);
  fChain->SetBranchAddress("centrality"                  , &centrality);

  fChain->SetBranchAddress("muon_pair_muon1_pt"          , &muon_pair_muon1_pt);
  fChain->SetBranchAddress("muon_pair_muon1_eta"         , &muon_pair_muon1_eta);
  fChain->SetBranchAddress("muon_pair_muon1_phi"         , &muon_pair_muon1_phi);
  fChain->SetBranchAddress("muon_pair_muon1_quality"     , &muon_pair_muon1_quality);
  fChain->SetBranchAddress("muon_pair_muon1_d0"          , &muon_pair_muon1_d0);
  fChain->SetBranchAddress("muon_pair_muon1_z0"          , &muon_pair_muon1_z0);
  fChain->SetBranchAddress("muon_pair_muon1_index"       , &muon_pair_muon1_index);

  fChain->SetBranchAddress("muon_pair_muon2_pt"          , &muon_pair_muon2_pt);
  fChain->SetBranchAddress("muon_pair_muon2_eta"         , &muon_pair_muon2_eta);
  fChain->SetBranchAddress("muon_pair_muon2_phi"         , &muon_pair_muon2_phi);
  fChain->SetBranchAddress("muon_pair_muon2_quality"     , &muon_pair_muon2_quality);
  fChain->SetBranchAddress("muon_pair_muon2_d0"          , &muon_pair_muon2_d0);
  fChain->SetBranchAddress("muon_pair_muon2_z0"          , &muon_pair_muon2_z0);
  fChain->SetBranchAddress("muon_pair_muon2_index"       , &muon_pair_muon2_index);

  fChain->SetBranchAddress("muon_deltaP_overP"           , &muon_deltaP_overP);

  fChain->SetBranchAddress("b_HLT_2mu4"                   , &b_HLT_2mu4);
  fChain->SetBranchAddress("b_HLT_mu4_mu4noL1"            , &b_HLT_mu4_mu4noL1);
  fChain->SetBranchAddress("b_HLT_mu4_mu4noL1_isPrescaled", &b_HLT_mu4_mu4noL1_isPrescaled);
  fChain->SetBranchAddress("f_HLT_mu4_mu4noL1_prescale"   , &f_HLT_mu4_mu4noL1_prescale);

  fChain->SetBranchAddress("b_HLT_mu4"                   , &b_HLT_mu4);
  fChain->SetBranchAddress("b_HLT_mu6"                   , &b_HLT_mu6);
  fChain->SetBranchAddress("b_HLT_mu8"                   , &b_HLT_mu8);

  fChain->SetBranchAddress("dimuon_b_HLT_2mu4"           , &dimuon_b_HLT_2mu4);
  fChain->SetBranchAddress("dimuon_b_HLT_mu4_mu4noL1"    , &dimuon_b_HLT_mu4_mu4noL1);
  fChain->SetBranchAddress("dimuon_b_HLT_mu4_mu4noL1_swap", &dimuon_b_HLT_mu4_mu4noL1_swap);

  fChain->SetBranchAddress("muon_b_HLT_mu4"              , &muon_b_HLT_mu4);
  fChain->SetBranchAddress("muon_b_HLT_mu6"              , &muon_b_HLT_mu6);
  fChain->SetBranchAddress("muon_b_HLT_mu8"              , &muon_b_HLT_mu8);


  // SetBranch Status
  fChain->SetBranchStatus("*"                               ,0);
  fChain->SetBranchStatus("RunNumber"                       ,1);
  fChain->SetBranchStatus("lbn"                             ,1);
  fChain->SetBranchStatus("centrality"                      ,1);

  fChain->SetBranchStatus("muon_pair_muon1_pt"              ,1);
  fChain->SetBranchStatus("muon_pair_muon1_eta"             ,1);
  fChain->SetBranchStatus("muon_pair_muon1_phi"             ,1);
  fChain->SetBranchStatus("muon_pair_muon1_quality"         ,1);
  fChain->SetBranchStatus("muon_pair_muon1_d0"              ,1);
  fChain->SetBranchStatus("muon_pair_muon1_z0"              ,1);
  fChain->SetBranchStatus("muon_pair_muon1_index"           ,1);

  fChain->SetBranchStatus("muon_pair_muon2_pt"              ,1);
  fChain->SetBranchStatus("muon_pair_muon2_eta"             ,1);
  fChain->SetBranchStatus("muon_pair_muon2_phi"             ,1);
  fChain->SetBranchStatus("muon_pair_muon2_quality"         ,1);
  fChain->SetBranchStatus("muon_pair_muon2_d0"              ,1);
  fChain->SetBranchStatus("muon_pair_muon2_z0"              ,1);
  fChain->SetBranchStatus("muon_pair_muon2_index"           ,1);

  fChain->SetBranchStatus("muon_deltaP_overP"               ,1);

  fChain->SetBranchStatus("b_HLT_2mu4"                      ,1);
  fChain->SetBranchStatus("b_HLT_mu4_mu4noL1"               ,1);
  fChain->SetBranchStatus("b_HLT_mu4_mu4noL1_isPrescaled"   ,1);
  fChain->SetBranchStatus("f_HLT_mu4_mu4noL1_prescale"      ,1);

  fChain->SetBranchStatus("b_HLT_mu4"                       ,1);
  fChain->SetBranchStatus("b_HLT_mu6"                       ,1);
  fChain->SetBranchStatus("b_HLT_mu8"                       ,1);

  fChain->SetBranchStatus("dimuon_b_HLT_2mu4"               ,1);
  fChain->SetBranchStatus("dimuon_b_HLT_mu4_mu4noL1"        ,1);
  fChain->SetBranchStatus("dimuon_b_HLT_mu4_mu4noL1_swap"   ,1);

  fChain->SetBranchStatus("muon_b_HLT_mu4"                  ,1);
  fChain->SetBranchStatus("muon_b_HLT_mu6"                  ,1);
  fChain->SetBranchStatus("muon_b_HLT_mu8"                  ,1);
}


void TrigEffMu4NoL1::InitHists(int read_flag){
  string s="_cut"+std::to_string(m_cutflags);
  if(m_cutflags==DefaultCuts)  s="";
  string name="01RootFiles/TrigEffMu4NoL1_MuonSelectionCut"+std::to_string(m_quality)+s+".root";
  if(read_flag==0) m_OutFile=new TFile(name.c_str(),"recreate");
  else             m_OutFile=new TFile(name.c_str(),"read");

  for(int iref=0;iref<NRefTrigs;iref++){
    for(int itrig=0;itrig<NTRIGS;itrig++){
      char name0[600],name1[600],name2[600];
      //------------------------------------------------
      if(read_flag==0){
        sprintf(name0,"h_eff_Num_ref%d_trig%d" ,iref,itrig);
        sprintf(name1,"h_eff_Den_ref%d_trig%d" ,iref,itrig);
        h_eff_Num[iref][itrig]=new TH2D(name0,";q*#it{#eta};#it{p}_{T} [GeV];Counts",50,-2.5,2.5,NBINSPT,eff_pt_bins);
        h_eff_Den[iref][itrig]=new TH2D(name1,";q*#it{#eta};#it{p}_{T} [GeV];Counts",50,-2.5,2.5,NBINSPT,eff_pt_bins);

        sprintf(name0,"h_eff_Num2_ref%d_trig%d" ,iref,itrig);
        sprintf(name1,"h_eff_Den2_ref%d_trig%d" ,iref,itrig);
        h_eff_Num2[iref][itrig]=new TH2D(name0,";Centrality [%];#it{p}_{T} [GeV];Counts",NCENT,CentBins,160,4,20);
        h_eff_Den2[iref][itrig]=new TH2D(name1,";Centrality [%];#it{p}_{T} [GeV];Counts",NCENT,CentBins,160,4,20);
        h_eff_Num2[iref][itrig]->Sumw2();
        h_eff_Den2[iref][itrig]->Sumw2();

        if(iref==RefAll && itrig==Mu4NoL1){
          for(int icent=0;icent<NCENT;icent++){
            sprintf(name0,"h_Dpop_cent%d",icent);
            h_Dpop[icent]=new TH2D(name0,";#frac{#Delta p}{p};#it{p}_{T} [GeV];Counts",50,-1,1,NPT,pTbins);
          }
        }
        if(iref==0 && itrig==0){ //New
          h_mon_npair=new TH1D("h_mon_npair",";#pairs;#events",101,-0.5,100.5); //New
          h_mon_cent =new TH1D("h_mon_cent" ,";#cent ;#events",101,-0.5,100.5); //New
          h_cut_flow =new TH1D("h_cut_flow" ,";#cut  ;#pairs" ,16 ,-0.5,15.5 ); //New
        } //New
      }
      //------------------------------------------------
      else if(read_flag==1){
        sprintf(name0,"h_eff_Num_ref%d_trig%d" ,iref,itrig);
        sprintf(name1,"h_eff_Den_ref%d_trig%d" ,iref,itrig);
        sprintf(name2,"h_eff_ref%d_trig%d"     ,iref,itrig);
        h_eff_Num[iref][itrig]=(TH2D*)m_OutFile->Get(name0);
        h_eff_Den[iref][itrig]=(TH2D*)m_OutFile->Get(name1);
        h_eff    [iref][itrig]=(TH2D*)m_OutFile->Get(name2);

        sprintf(name0,"h_eff_Num2_ref%d_trig%d" ,iref,itrig);
        sprintf(name1,"h_eff_Den2_ref%d_trig%d" ,iref,itrig);
        h_eff_Num2[iref][itrig]=(TH2D*)m_OutFile->Get(name0);
        h_eff_Den2[iref][itrig]=(TH2D*)m_OutFile->Get(name1);

        for(int ieta=0;ieta<NETA;ieta++){
          sprintf(name2,"h_eff1D_ref%d_trig%d_eta%d",iref,itrig,ieta);
          h_eff1D[iref][itrig][ieta]=(TH1D*)m_OutFile->Get(name2);
        }
      }
      //------------------------------------------------
      //------------------------------------------------
      else if(read_flag==2){
        if(iref!=RefAll) continue;
        for(int ieta=0;ieta<NETA;ieta++){
          sprintf(name2,"h_eff1D_ref%d_trig%d_eta%d",iref,itrig,ieta);
          h_eff1D[iref][itrig][ieta]=(TH1D*)m_OutFile->Get(name2);
        }
      }
      //------------------------------------------------
    }
  }
}



int TrigEffMu4NoL1::GetEtaBin(float eta,float charge){
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

int TrigEffMu4NoL1::GetCentBin(int centrality){
  for(int ibin=0;ibin<NCENT;ibin++){
    if (centrality<CentBins[ibin+1]) return ibin; 
  }
  return NCENT-1;
}



float TrigEffMu4NoL1::GetTriggerEfficiency(int itrig,float eta, float pt, float charge, int cent_percentile){
  int ieta=GetEtaBin(eta,charge);
  if(ieta<0){
    cout<<__PRETTY_FUNCTION__<<": exception at line "<<__LINE__<<":  eta="<<eta<<std::endl;
    throw std::exception();
  }
  if(pt>=20) pt=19.5;

  float scale=CentMultiplier(itrig,pt,cent_percentile);

  int ibin = h_eff1D[RefAll][itrig][ieta]->GetXaxis()->FindFixBin(pt);
  return     h_eff1D[RefAll][itrig][ieta]->GetBinContent(ibin) * scale;
}


float TrigEffMu4NoL1::CentMultiplier(int itrig, float pt, int cent_percentile){
  float scale=1.0;
  if(itrig==Mu4){
    if(pt<5){  //4-5 
      if     (cent_percentile<10) scale=0.943;
      else if(cent_percentile<20) scale=1.015;
      else if(cent_percentile<30) scale=1.049;
      else if(cent_percentile<40) scale=1.065;
      else if(cent_percentile<60) scale=1.080;
      else                        scale=1.080;
    }
    else if(pt<6){  //5-6
      if     (cent_percentile<10) scale=0.963;
      else if(cent_percentile<20) scale=1.005;
      else if(cent_percentile<30) scale=1.026;
      else if(cent_percentile<40) scale=1.045;
      else if(cent_percentile<60) scale=1.046;
      else                        scale=1.046;
    }
    else if(pt<8){  //6-8
      if     (cent_percentile<10) scale=0.967;
      else if(cent_percentile<20) scale=1.008;
      else if(cent_percentile<30) scale=1.019;
      else if(cent_percentile<40) scale=1.030;
      else if(cent_percentile<60) scale=1.030;
      else                        scale=1.030;
    }
    else{  //>8
      if     (cent_percentile<10) scale=0.985;
      else if(cent_percentile<20) scale=1.000;
      else if(cent_percentile<30) scale=1.005;
      else if(cent_percentile<40) scale=1.008;
      else if(cent_percentile<60) scale=1.010;
      else                        scale=1.010;
    }
  }
  else if(itrig==Mu4NoL1){
    if(pt<5){  //4-5
      if     (cent_percentile<10) scale=0.978;
      else if(cent_percentile<20) scale=1.005;
      else if(cent_percentile<30) scale=1.023;
      else if(cent_percentile<40) scale=1.028;
      else if(cent_percentile<60) scale=1.030;
      else                        scale=1.030;
    }
    else if(pt<6){  //5-6
      if     (cent_percentile<10) scale=0.993;
      else if(cent_percentile<20) scale=1.002;
      else if(cent_percentile<30) scale=1.007;
      else if(cent_percentile<40) scale=1.007;
      else if(cent_percentile<60) scale=1.007;
      else                        scale=1.007;
    }
    else if(pt<8){  //6-8
      if     (cent_percentile<10) scale=0.996;
      else if(cent_percentile<20) scale=1.000;
      else if(cent_percentile<30) scale=1.005;
      else if(cent_percentile<40) scale=1.005;
      else if(cent_percentile<60) scale=1.005;
      else                        scale=1.005;
    }
    else{  //>8
      if     (cent_percentile<10) scale=0.995;
      else if(cent_percentile<20) scale=1.000;
      else if(cent_percentile<30) scale=1.003;
      else if(cent_percentile<40) scale=1.003;
      else if(cent_percentile<60) scale=1.003;
      else                        scale=1.003;
    }
  }
  else if(itrig==Both){
    if(pt<5){  //4-5
      if     (cent_percentile<10) scale=0.938;
      else if(cent_percentile<20) scale=1.018;
      else if(cent_percentile<30) scale=1.057;
      else if(cent_percentile<40) scale=1.075;
      else if(cent_percentile<60) scale=1.090;
      else                        scale=1.090;
    }
    else if(pt<6){  //5-6
      if     (cent_percentile<10) scale=0.960;
      else if(cent_percentile<20) scale=1.005;
      else if(cent_percentile<30) scale=1.029;
      else if(cent_percentile<40) scale=1.048;
      else if(cent_percentile<60) scale=1.048;
      else                        scale=1.048;
    }
    else if(pt<8){  //6-8
      if     (cent_percentile<10) scale=0.967;
      else if(cent_percentile<20) scale=1.009;
      else if(cent_percentile<30) scale=1.020;
      else if(cent_percentile<40) scale=1.033;
      else if(cent_percentile<60) scale=1.033;
      else                        scale=1.033;
    }
    else{  //>8
      if     (cent_percentile<10) scale=0.980;
      else if(cent_percentile<20) scale=1.000;
      else if(cent_percentile<30) scale=1.005;
      else if(cent_percentile<40) scale=1.011;
      else if(cent_percentile<60) scale=1.013;
      else                        scale=1.013;
    }
  }
  return scale;
}


float TrigEffMu4NoL1::CentMultiplierMedium(int itrig, float pt, int cent_percentile){
  float scale=1.0;
  if(itrig==Mu4){
    if(pt<5){  //4-5 
      if     (cent_percentile<10) scale=0.925;
      else if(cent_percentile<20) scale=1.019;
      else if(cent_percentile<30) scale=1.070;
      else if(cent_percentile<40) scale=1.090;
      else if(cent_percentile<60) scale=1.120;
      else                        scale=1.125;
    }
    else if(pt<6){  //5-6
      if     (cent_percentile<10) scale=0.950;
      else if(cent_percentile<20) scale=1.007;
      else if(cent_percentile<30) scale=1.035;
      else if(cent_percentile<40) scale=1.065;
      else if(cent_percentile<60) scale=1.078;
      else                        scale=1.078;
    }
    else if(pt<8){  //6-8
      if     (cent_percentile<10) scale=0.961;
      else if(cent_percentile<20) scale=1.005;
      else if(cent_percentile<30) scale=1.025;
      else if(cent_percentile<40) scale=1.045;
      else if(cent_percentile<60) scale=1.045;
      else                        scale=1.045;
    }
    else{  //>8
      if     (cent_percentile<10) scale=0.974;
      else if(cent_percentile<20) scale=1.000;
      else if(cent_percentile<30) scale=1.010;
      else if(cent_percentile<40) scale=1.015;
      else if(cent_percentile<60) scale=1.017;
      else                        scale=1.017;
    }
  }
  else if(itrig==Mu4NoL1){
    if(pt<5){  //4-5
      if     (cent_percentile<10) scale=0.96 ;
      else if(cent_percentile<20) scale=1.005;
      else if(cent_percentile<30) scale=1.04 ;
      else if(cent_percentile<40) scale=1.06 ;
      else if(cent_percentile<60) scale=1.075;
      else                        scale=1.09 ;
    }
    else if(pt<6){  //5-6
      if     (cent_percentile<10) scale=0.985;
      else if(cent_percentile<20) scale=1.001;
      else if(cent_percentile<30) scale=1.012;
      else if(cent_percentile<40) scale=1.020;
      else if(cent_percentile<60) scale=1.025;
      else                        scale=1.030;
    }
    else if(pt<8){  //6-8
      if     (cent_percentile<10) scale=0.992;
      else if(cent_percentile<20) scale=1.000;
      else if(cent_percentile<30) scale=1.005;
      else if(cent_percentile<40) scale=1.007;
      else if(cent_percentile<60) scale=1.010;
      else                        scale=1.010;
    }
    else{  //>8
      if     (cent_percentile<10) scale=0.995;
      else if(cent_percentile<20) scale=1.000;
      else if(cent_percentile<30) scale=1.003;
      else if(cent_percentile<40) scale=1.003;
      else if(cent_percentile<60) scale=1.003;
      else                        scale=1.003;
    }
  }
  else if(itrig==Both){
    if(pt<5){  //4-5
      if     (cent_percentile<10) scale=0.920;
      else if(cent_percentile<20) scale=1.020;
      else if(cent_percentile<30) scale=1.080;
      else if(cent_percentile<40) scale=1.100;
      else if(cent_percentile<60) scale=1.140;
      else                        scale=1.140;
    }
    else if(pt<6){  //5-6
      if     (cent_percentile<10) scale=0.945;
      else if(cent_percentile<20) scale=1.008;
      else if(cent_percentile<30) scale=1.042;
      else if(cent_percentile<40) scale=1.070;
      else if(cent_percentile<60) scale=1.085;
      else                        scale=1.085;
    }
    else if(pt<8){  //6-8
      if     (cent_percentile<10) scale=0.957;
      else if(cent_percentile<20) scale=1.005;
      else if(cent_percentile<30) scale=1.025;
      else if(cent_percentile<40) scale=1.045;
      else if(cent_percentile<60) scale=1.045;
      else                        scale=1.045;
    }
    else{  //>8
      if     (cent_percentile<10) scale=0.970;
      else if(cent_percentile<20) scale=1.000;
      else if(cent_percentile<30) scale=1.010;
      else if(cent_percentile<40) scale=1.015;
      else if(cent_percentile<60) scale=1.017;
      else                        scale=1.017;
    }
  }
  return scale;
}


void TrigEffMu4NoL1::Process(){
  std::cout<<"Starting"<<std::endl;

  Long64_t nentries = fChain->GetEntries();

  int l_qual=0;
  if      (m_quality==Bins::QualityCut::MEDIUM) l_qual= 8;
  else if (m_quality==Bins::QualityCut::TIGHT ) l_qual=16;
  else    Common::Exception(__LINE__,__FILE__);
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%10000==0) cout<<jentry<<" "<<nentries<<endl;
    fChain->GetEntry(jentry);

    int NPairs=muon_pair_muon1_pt->size();
    std::set<int> l_CheckUnique[NRefTrigs];

    const int icent=GetCentBin(centrality);
    h_mon_npair->Fill(NPairs    ); //New
    h_mon_cent ->Fill(centrality); //New

    if(jentry%10000==0) cout<<jentry<<" "<<nentries<<" "<<NPairs<<endl;
    for(int i=0;i<NPairs;i++){
      const float pt1 = fabs(muon_pair_muon1_pt ->at(i))/1000.0;
      const float pt2 = fabs(muon_pair_muon2_pt ->at(i))/1000.0;
      const float eta1=      muon_pair_muon1_eta->at(i);
      const float eta2=      muon_pair_muon2_eta->at(i);
      const float phi1=      muon_pair_muon1_phi->at(i);
      const float phi2=      muon_pair_muon2_phi->at(i);
      const int q1    =     (muon_pair_muon1_pt ->at(i)>0)? 1:-1;
      const int q2    =     (muon_pair_muon2_pt ->at(i)>0)? 1:-1;

      int icut=0; //New
      h_cut_flow ->Fill(icut++); //New 1
      //exclude nearby muons
      {
        float deta=eta1-eta2;
        float dphi=phi1-phi2;
        dphi=atan2(sin(dphi),cos(dphi));
        float dR2=deta*deta + dphi*dphi;
        if(dR2<0.16) continue; //dR^2<0.14 ==> dr<0.4;
      }
      h_cut_flow ->Fill(icut++); //New 2

      //require medium or tight muons depending on flag
      const int quality1=muon_pair_muon1_quality->at(i);
      const int quality2=muon_pair_muon2_quality->at(i);
      if((quality1&quality2&1     )==0) continue;//combined muon
      h_cut_flow ->Fill(icut++); //New 3
      if((quality1&quality2&l_qual)==0) continue;//Medium/Tight muon
      h_cut_flow ->Fill(icut++); //New 4
      if((m_cutflags&NoIDCuts)==0 )
      {
        if((quality1&quality2&32 )==0) continue;//IDCuts
        if((quality1&quality2&256)==0) continue;//MuonCuts
      }
      h_cut_flow ->Fill(icut++); //New 5

      if( (pt1<3.7) || (pt2<3.7) || (fabs(eta1)>2.4) || (fabs(eta2)>2.4) ) continue;
      h_cut_flow ->Fill(icut++); //New 6

      const float d01   =muon_pair_muon1_d0->at(i);
      const float d02   =muon_pair_muon2_d0->at(i);
      const float z0Sin1=muon_pair_muon1_z0->at(i)*sin(2.0*atan(exp(-eta1)));
      const float z0Sin2=muon_pair_muon2_z0->at(i)*sin(2.0*atan(exp(-eta2)));
      if(fabs(d01)>1 || fabs(d02)>1 || fabs(z0Sin1)>1 || fabs(z0Sin2)>1) continue;
      h_cut_flow ->Fill(icut++); //New 7

      const int index1=muon_pair_muon1_index->at(i);
      const int index2=muon_pair_muon2_index->at(i);
 
      //momentum-imbalance
      const double dpop1= muon_deltaP_overP->at(index1);
      const double dpop2= muon_deltaP_overP->at(index2);
      if     (m_cutflags&DpopCut0p1 ){if(dpop1>0.1     || dpop2>0.1    ) continue;}
      else if(m_cutflags&DpopCut0p15){if(dpop1>0.15    || dpop2>0.15   ) continue;}
      else if(m_cutflags&DpopCut0p2 ){if(dpop1>0.2     || dpop2>0.2    ) continue;}
      if(fabs(dpop1)>1 || fabs(dpop2)>1) continue;
      h_cut_flow ->Fill(icut++); //New 8


      //Trigger match for muon pair
      bool trig_match_2mu4        = dimuon_b_HLT_2mu4       ->at(i);
      bool trig_match_mu4_mu4noL1 = dimuon_b_HLT_mu4_mu4noL1->at(i);
      bool trig_match_both        = trig_match_2mu4 && trig_match_mu4_mu4noL1;
      bool is_prescaled           = b_HLT_mu4_mu4noL1_isPrescaled;
      //For 2015 the b_HLT_mu4_mu4noL1_isPrescaled is incorrect, we instead use the relation below
      if(RunNumber<300000 ) is_prescaled=(f_HLT_mu4_mu4noL1_prescale<0.9999 || f_HLT_mu4_mu4noL1_prescale>1.0001);
      if(trig_match_mu4_mu4noL1!=dimuon_b_HLT_mu4_mu4noL1_swap->at(i)) {Common::Exception(__LINE__,__FILE__);}//Check
      if(is_prescaled) continue;
      h_cut_flow ->Fill(icut++); //New 9
      //cout<<"A1"<<endl;


      for(int iref=0;iref<NRefTrigs;iref++){
        bool pass1=false;
        bool pass2=false;
        if      (iref==RefMu4) {
          if(       !b_HLT_mu4) continue;
          h_cut_flow ->Fill(icut++); //New 10
          pass1=muon_b_HLT_mu4->at(index1);
          pass2=muon_b_HLT_mu4->at(index2);
          if(!(pass1||pass2)) continue;//New
          h_cut_flow ->Fill(icut++); //New 11
          if(!(pt1>4 || pt2>4)) continue;//New
          h_cut_flow ->Fill(icut++); //New 12
        }
        else if (iref==RefMu6) {
          if(       !b_HLT_mu6) continue;
          pass1=muon_b_HLT_mu6->at(index1);
          pass2=muon_b_HLT_mu6->at(index2);
        }
        else if (iref==RefMu8) {
          if(       !b_HLT_mu8) continue;
          pass1=muon_b_HLT_mu8->at(index1);
          pass2=muon_b_HLT_mu8->at(index2);
        }
        else if(iref==RefAll){
          if(!      (b_HLT_mu4             ||      b_HLT_mu6             ||      b_HLT_mu8   )) continue;
          pass1=muon_b_HLT_mu4->at(index1) || muon_b_HLT_mu6->at(index1) || muon_b_HLT_mu8->at(index1);
          pass2=muon_b_HLT_mu4->at(index2) || muon_b_HLT_mu6->at(index2) || muon_b_HLT_mu8->at(index2);
        }
        else Common::Exception(__LINE__,__FILE__);

        //first muon matches single-muon trigger
        if(pass1 && pt1>4){
          if(l_CheckUnique[iref].count(index2)==0)//only use this muon if it was not already used
          {
                                       h_eff_Den[iref][Mu4    ]->Fill(q2*eta2,pt2);h_eff_Den2[iref][Mu4    ]->Fill(centrality,pt2);
                                       h_eff_Den[iref][Mu4NoL1]->Fill(q2*eta2,pt2);h_eff_Den2[iref][Mu4NoL1]->Fill(centrality,pt2);
                                       h_eff_Den[iref][Both   ]->Fill(q2*eta2,pt2);h_eff_Den2[iref][Both   ]->Fill(centrality,pt2);
            if(trig_match_2mu4       ){h_eff_Num[iref][Mu4    ]->Fill(q2*eta2,pt2);h_eff_Num2[iref][Mu4    ]->Fill(centrality,pt2);}
            if(trig_match_mu4_mu4noL1){h_eff_Num[iref][Mu4NoL1]->Fill(q2*eta2,pt2);h_eff_Num2[iref][Mu4NoL1]->Fill(centrality,pt2);}
            if(trig_match_both       ){h_eff_Num[iref][Both   ]->Fill(q2*eta2,pt2);h_eff_Num2[iref][Both   ]->Fill(centrality,pt2);}
            l_CheckUnique[iref].insert(index2);
            if(iref==RefAll){
              h_Dpop[icent]->Fill(dpop2,pt2);
            }
            if(iref==RefMu4) h_cut_flow ->Fill(icut); //New 13
          }
        }
        //Second muon matches single-muon trigger
        if(pass2 && pt2>4){
          if(l_CheckUnique[iref].count(index1)==0)//only use this muon if it was not already used
          {
                                       h_eff_Den[iref][Mu4    ]->Fill(q1*eta1,pt1);h_eff_Den2[iref][Mu4    ]->Fill(centrality,pt1);
                                       h_eff_Den[iref][Mu4NoL1]->Fill(q1*eta1,pt1);h_eff_Den2[iref][Mu4NoL1]->Fill(centrality,pt1);
                                       h_eff_Den[iref][Both   ]->Fill(q1*eta1,pt1);h_eff_Den2[iref][Both   ]->Fill(centrality,pt1);
            if(trig_match_2mu4       ){h_eff_Num[iref][Mu4    ]->Fill(q1*eta1,pt1);h_eff_Num2[iref][Mu4    ]->Fill(centrality,pt1);}
            if(trig_match_mu4_mu4noL1){h_eff_Num[iref][Mu4NoL1]->Fill(q1*eta1,pt1);h_eff_Num2[iref][Mu4NoL1]->Fill(centrality,pt1);}
            if(trig_match_both       ){h_eff_Num[iref][Both   ]->Fill(q1*eta1,pt1);h_eff_Num2[iref][Both   ]->Fill(centrality,pt1);}
            l_CheckUnique[iref].insert(index1);
            if(iref==RefAll){
              h_Dpop[icent]->Fill(dpop1,pt1);
            }
            if(iref==RefMu4) h_cut_flow ->Fill(icut); //New 13
          }
        }
      }
    }
  }
}



void TrigEffMu4NoL1::EvalEffs(){
  char name[600];
  for(int iref=0;iref<NRefTrigs;iref++){
    for(int itrig=0;itrig<NTRIGS;itrig++){

      sprintf(name,"h_eff_ref%d_trig%d",iref,itrig);
      h_eff  [iref][itrig]=(TH2D*) h_eff_Num[iref][itrig]->Clone(name);
      h_eff  [iref][itrig]->Divide(h_eff_Den[iref][itrig]);

      for(int ieta=0;ieta<NETA;ieta++){
        sprintf(name,"h_eff1D_ref%d_trig%d_eta%d",iref,itrig,ieta);
                   h_eff_Num[iref][itrig]->GetXaxis()->SetRangeUser(EtaMin[ieta]+1e-6,EtaMax[ieta]-1e-6);
                   h_eff_Den[iref][itrig]->GetXaxis()->SetRangeUser(EtaMin[ieta]+1e-6,EtaMax[ieta]-1e-6);
        TH1D* hNum=h_eff_Num[iref][itrig]->ProjectionY(name);
        TH1D* hDen=h_eff_Den[iref][itrig]->ProjectionY(Common::UniqueName().c_str());
        hNum->Divide(hDen);
        h_eff1D[iref][itrig][ieta]=hNum;
        delete hDen;
      }
      h_eff_Num[iref][itrig]->GetXaxis()->SetRangeUser(-2.5,2.5);
      h_eff_Den[iref][itrig]->GetXaxis()->SetRangeUser(-2.5,2.5);
    }
  }
}

#endif
