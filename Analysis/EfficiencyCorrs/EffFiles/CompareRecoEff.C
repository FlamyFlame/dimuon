#include "common.C"

//Compare the Reco Effs between the different files
void CompareRecoEff(){

  //MCP pp effs
  //--------------------------------------------------------------------
  TCanvas *C1=Common::StandardCanvas9("Can_Compare_RecoEff");

  TFile *_file0 = TFile::Open("210222_Precision_r21/Reco_Medium_JPsi.root");
  TFile *_file1 = TFile::Open("210222_Precision_r21/Reco_Tight_JPsi.root");
  TH2D  *hist0  =(TH2D*) (_file0->Get("Eff_2017")->Clone(Common::UniqueName().c_str()));
  TH2D  *hist1  =(TH2D*) (_file1->Get("Eff_2017")->Clone(Common::UniqueName().c_str()));

  
  for(int idraw=0;idraw<2;idraw++){
    for(int ipad=0;ipad<9;ipad++){
      C1->cd(ipad+1);

      TH1D* h1=nullptr;
      if(idraw==0) h1=hist0->ProjectionX(Common::UniqueName().c_str(),ipad+1,ipad+1);
      if(idraw==1) h1=hist1->ProjectionX(Common::UniqueName().c_str(),ipad+1,ipad+1);


      int col[]={ 1, 2};
      int sty[]={24,24};
      Common::format(h1,col[idraw],sty[idraw]);
      
      if(idraw==0){
        h1->Draw();
        h1->GetYaxis()->SetRangeUser(0,1.2);
        h1->GetXaxis()->SetRangeUser(4,10 );
      }
      else h1->Draw("same");
    }
  }
  //--------------------------------------------------------------------


  //My PbPb Effs
  //--------------------------------------------------------------------
  int idraw=0;
  for(std::string filename:{"../01RootFiles/MC_MuonSelection0.root","../01RootFiles/MC_MuonSelection1.root"}){
    TFile *file = TFile::Open(filename.c_str());
    TH2D* hist3a=(TH2D*)file->Get("h_Eff_Num_2D");
    TH2D* hist3b=(TH2D*)file->Get("h_Eff_Den_2D");

    for(int ipad=0;ipad<9;ipad++){
      C1->cd(ipad+1);
      //NOTE: Not exact matches with the MCP histograms
      const double etabins[]={-2.5,-2.0,-1.3,-1.0,-0.1,0.1,1.0,1.3,2.0,2.5};

      hist3a->GetXaxis()->SetRangeUser(etabins[ipad]+1e-3,etabins[ipad]-1e-3);
      hist3b->GetXaxis()->SetRangeUser(etabins[ipad]+1e-3,etabins[ipad]-1e-3);

      TH1D* h_num=hist3a->ProjectionY(Common::UniqueName().c_str());
      TH1D* h_den=hist3b->ProjectionY(Common::UniqueName().c_str());
      h_num->Divide(h_den);
      h_num->Draw("same");

      int col[]={ 1, 2};
      int sty[]={25,25};
      Common::format(h_num,col[idraw],sty[idraw]);
    }
    idraw++;
  }
  //--------------------------------------------------------------------

  
}