#ifndef __COMMON_C__
#define __COMMON_C__

#include <TLatex.h>
#include <TMarker.h>
#include "TLine.h"
#include "map"

#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TProfile.h"
#include "THStack.h"
#include "TRandom3.h"
#include "TSystem.h"


//#define LONG_PAPER


#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphErrors.h"

using std::vector;
using std::map;
using std::string;
using std::cout;
using std::endl;
using std::pair;

void MyExit(){
  gSystem->Exec(Form("kill %d", getpid()));//terminate the program to exit faster
}

namespace Common{
double PI=acos(-1.0);
//std::string Internal="<Internal>";
std::string Internal="";
//string Internal="Preliminary";

bool CheckObject(const TObject *obj,std::string name,TFile *file=nullptr){
  std::string file_description="";
  if(file){
    file_description=" in TFile ";
    file_description+=file->GetName();
  }
  if(!obj) {
    std::cout<<"*******************************************"<<std::endl;
    std::cout<<name<<" Object is invalid"<<file_description<<std::endl;
    std::cout<<"*******************************************"<<std::endl;
    return false;
  }
  if(strcmp(obj->GetName(),name.c_str())) {
    std::cout<<"*******************************************"<<std::endl;
    std::cout<<"Expected="<<obj->GetName()<<"  Found="<<name<<file_description<<std::endl;
    std::cout<<"*******************************************"<<std::endl;
    return false;
  }
  return true;
}

bool CheckObject2(const TObject *obj,std::string name,TFile *file=nullptr){
  //if(!CheckObject(obj,name,file)) throw std::exception();
  if(!CheckObject(obj,name,file)) return true;//throw std::exception();
  else return true;
}

void CheckFile(const TFile *file,std::string name){
  if(file->IsZombie()){
    std::cout<<name<<" File is Zombie"<<std::endl;
    throw std::exception();
  }
}

void Exception(int line, std::string file,std::string message="");

std::string UniqueName(){
  static int UniqueNameCounter=0;
  std::string str="UniqueName";
  str+=std::to_string(UniqueNameCounter);
  UniqueNameCounter++;
  return str;
}
std::string UniqueName2(){
  static int UniqueNameCounter=0;
	std::string str=std::to_string(UniqueNameCounter);
  UniqueNameCounter++;
  return str;
}

TGraphErrors *FlipGraphAxes(TGraphErrors*gr1){
  if(!gr1) return nullptr;

  double x[100],y[100],xe[100],ye[100];

  int N=gr1->GetN();
  if(N>100){
    std::cout<<"FlipGraphAxes too many points"<<std::endl;
    throw std::exception();
  }
  Double_t *X1=gr1->GetX();
  Double_t *Y1=gr1->GetY();
  Double_t *XE1=gr1->GetEX();
  Double_t *YE1=gr1->GetEY();

  for(int i=0;i<N;i++){
    x[i]  = Y1[i];
    y[i]  = X1[i];
    ye[i] = XE1[i];
    xe[i] = YE1[i];
  }
  TGraphErrors*gr=  new TGraphErrors(N,x,y,xe,ye);
  gr->SetMarkerColor(gr1->GetMarkerColor());
  gr->SetLineColor  (gr1->GetLineColor  ());
  gr->SetFillColor  (gr1->GetFillColor  ());
  gr->SetMarkerSize (gr1->GetMarkerSize ());
  gr->SetLineWidth  (gr1->GetLineWidth  ());
  gr->SetLineStyle  (gr1->GetLineStyle  ());
  gr->SetMarkerStyle(gr1->GetMarkerStyle());
  return gr;
}




TGraphErrors *h2gr(TH1*h){
  if(!h) return 0;
  int N = h->GetNbinsX();
  if(N>100){
    std::cout<<"h2gr() too many points"<<std::endl;
    throw std::exception();
  }
  double x[100],y[100],ye[100];

  for(int i=0;i<N;i++){
    x[i]  = h->GetBinCenter (i+1);
    y[i]  = h->GetBinContent(i+1);
    ye[i] = h->GetBinError  (i+1);
  }
  TGraphErrors*gr=  new TGraphErrors(N,x,y,0,ye);
  gr->SetMarkerColor(h->GetMarkerColor());
  gr->SetLineColor  (h->GetLineColor  ());
  gr->SetFillColor  (h->GetFillColor  ());
  gr->SetMarkerSize (h->GetMarkerSize ());
  gr->SetLineWidth  (h->GetLineWidth  ());
  gr->SetLineStyle  (h->GetLineStyle  ());
  gr->SetMarkerStyle(h->GetMarkerStyle());
  return gr;
}


TGraphErrors *h2gr_Xerrors(TH1*h){
  if(!h) return 0;
  int N = h->GetNbinsX();
  if(N>100){
    std::cout<<"h2gr_Xerrors() too many points"<<std::endl;
    throw std::exception();
  }
  double x[100],y[100],xe[100],ye[100];

  for(int i=0;i<N;i++){
    x[i]  = h->GetBinCenter (i+1);
    y[i]  = h->GetBinContent(i+1);
    ye[i] = h->GetBinError  (i+1);
    xe[i] = (h->GetBinLowEdge(i+2)-h->GetBinLowEdge(i+1))/2.0;
  }
  TGraphErrors*gr=  new TGraphErrors(N,x,y,xe,ye);
  gr->SetMarkerColor(h->GetMarkerColor());
  gr->SetLineColor  (h->GetLineColor  ());
  gr->SetFillColor  (h->GetFillColor  ());
  gr->SetMarkerSize (h->GetMarkerSize ());
  gr->SetLineWidth  (h->GetLineWidth  ());
  gr->SetLineStyle  (h->GetLineStyle  ());
  gr->SetMarkerStyle(h->GetMarkerStyle());
  return gr;
}


TGraphErrors *ProfToTGraphErrors(TProfile *h){
  if(!h) return nullptr;
  int N = h->GetNbinsX();
  if(N>100){
    std::cout<<"ProfToGraph(TProfile *h) too many points"<<std::endl;
    throw std::exception();
  }

  int ibin=0;
  double x[100],y[100],ye[100];
  for(int i=0;i<N;i++){
    if(h->GetBinEntries(i+1)==0) continue;
    x [ibin] = h->GetBinCenter (i+1);
    y [ibin] = h->GetBinContent(i+1);
    ye[ibin] = h->GetBinError  (i+1);
    ibin++;
  }
  TGraphErrors*gr=  new TGraphErrors(ibin,x,y,0,ye);
  gr->SetMarkerColor(h->GetMarkerColor());
  gr->SetLineColor  (h->GetLineColor  ());
  gr->SetFillColor  (h->GetFillColor  ());
  gr->SetMarkerSize (h->GetMarkerSize ());
  gr->SetLineWidth  (h->GetLineWidth  ());
  gr->SetLineStyle  (h->GetLineStyle  ());
  gr->SetMarkerStyle(h->GetMarkerStyle());
  return gr;
}


TLegend* StandardLegend(float x1,float y1,float x2,float y2){
  TLegend *leg2 = new TLegend(x1,y1,x2,y2);
  leg2->SetTextSize(0.045);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  return leg2;
}


TCanvas* StandardCanvas1(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),600,450);
  Can->SetLeftMargin (0.125);
  Can->SetTopMargin  (0.05);
  Can->SetRightMargin(0.15);
  return Can;
}


TCanvas* StandardCanvas1a(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),600,450);
  Can->SetLeftMargin (0.12);
  Can->SetTopMargin  (0.01);
  Can->SetRightMargin(0.01);
  return Can;
}

TCanvas* StandardCanvas1b(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),600,450);
  Can->SetLeftMargin (0.15);
  Can->SetTopMargin  (0.01);
  Can->SetRightMargin(0.01);
  return Can;
}

TCanvas* StandardCanvas1d(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1000,450);
  Can->SetLeftMargin (0.125);
  Can->SetTopMargin  (0.05);
  Can->SetRightMargin(0.05);
  return Can;
}

TCanvas* StandardCanvas1e(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),800,450);
  Can->SetLeftMargin (0.125);
  Can->SetTopMargin  (0.05);
  Can->SetRightMargin(0.15);
  return Can;
}

TCanvas* StandardCanvas2(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1200,450);
  Can->Divide(2);
  for(int ipad:{0,1}){
    Can->GetPad(ipad+1)->SetLeftMargin (0.15);
    Can->GetPad(ipad+1)->SetTopMargin  (0.05);
    Can->GetPad(ipad+1)->SetRightMargin(0.15);
  }
  return Can;
}

TCanvas* StandardCanvas2a(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1200,800);
  Can->Divide(2);
  for(int ipad:{0,1}){
    Can->GetPad(ipad+1)->SetLeftMargin (0.12);
    Can->GetPad(ipad+1)->SetTopMargin  (0.05);
    Can->GetPad(ipad+1)->SetRightMargin(0.05);
  }
  return Can;
}

TCanvas* StandardCanvas2b(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1200,450);
  Can->Divide(2);
  for(int ipad:{0,1}){
    Can->GetPad(ipad+1)->SetLeftMargin (0.15);
    Can->GetPad(ipad+1)->SetTopMargin  (0.02);
    Can->GetPad(ipad+1)->SetRightMargin(0.05);
  }
  return Can;
}

TCanvas* StandardCanvas2d(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1000,900);
  Can->Divide(1,2);
  for(int ipad:{0,1}){
    Can->GetPad(ipad+1)->SetLeftMargin (0.125);
    Can->GetPad(ipad+1)->SetTopMargin  (0.05);
    Can->GetPad(ipad+1)->SetRightMargin(0.05);
  }
  return Can;
}

TCanvas* StandardCanvas2e(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1200,450);
  Can->Divide(2);
  for(int ipad:{0,1}){
    Can->GetPad(ipad+1)->SetLeftMargin (0.12);
    Can->GetPad(ipad+1)->SetTopMargin  (0.05);
    Can->GetPad(ipad+1)->SetRightMargin(0.02);
  }
  return Can;
}

TCanvas* StandardCanvas3(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1800,450);
  Can->Divide(3);
  for(int ipad:{0,1,2}){
    Can->GetPad(ipad+1)->SetLeftMargin  (0.12);
    Can->GetPad(ipad+1)->SetTopMargin   (0.01);
    Can->GetPad(ipad+1)->SetRightMargin (0.02);
    Can->GetPad(ipad+1)->SetBottomMargin(0.12);
  }
  return Can;
}

TCanvas* StandardCanvas3a(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1000,1400);
  Can->Divide(1,3);
  for(int ipad:{0,1,2}){
    Can->GetPad(ipad+1)->SetLeftMargin  (0.12);
    Can->GetPad(ipad+1)->SetTopMargin   (0.01);
    Can->GetPad(ipad+1)->SetRightMargin (0.02);
  }
  return Can;
}

TCanvas* StandardCanvas3b(std::string can_name){
  TCanvas *C1=new TCanvas(can_name.c_str(),can_name.c_str(),1800,450);
  C1->Divide(3,1,.001,.001);
  for(int ipad:{1,2,3}){
    C1->cd(ipad);
    gPad->SetTopMargin   (0.06);
    gPad->SetRightMargin (0.01);
    gPad->SetLeftMargin  (0.14);
    gPad->SetBottomMargin(0.05);
  }
  return C1;
}

TCanvas* StandardCanvas3c(std::string can_name){
  TCanvas *C1=new TCanvas(can_name.c_str(),can_name.c_str(),1800,450);
  C1->Divide(3,1,.001,.001);
  for(int ipad:{1,2,3}){
    C1->cd(ipad);
    gPad->SetTopMargin   (0.02);
    gPad->SetRightMargin (0.20);
    gPad->SetLeftMargin  (0.11);
    gPad->SetBottomMargin(0.15);
  }
  return C1;
}
TCanvas* StandardCanvas3d(std::string can_name){
  TCanvas *C1=new TCanvas(can_name.c_str(),can_name.c_str(),1800,450);
  C1->Divide(3,1,.001,.001);
  for(int ipad:{1,2,3}){
    C1->cd(ipad);
    gPad->SetTopMargin   (0.01);
    gPad->SetRightMargin (0.01);
    gPad->SetLeftMargin  (0.11);
    gPad->SetBottomMargin(0.15);
  }
  return C1;
}
TCanvas* StandardCanvas3e(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1800,450);
  Can->Divide(3);
  for(int ipad:{0,1,2}){
    Can->GetPad(ipad+1)->SetLeftMargin  (0.12);
    Can->GetPad(ipad+1)->SetTopMargin   (0.02);
    Can->GetPad(ipad+1)->SetRightMargin (0.02);
    Can->GetPad(ipad+1)->SetBottomMargin(0.12);
  }
  return Can;
}

TCanvas* StandardCanvas4(std::string can_name){
  TCanvas *C2=new TCanvas(can_name.c_str(),can_name.c_str(),1200,1000);
  C2->Divide(2,2);
  for(int ican:{1,2,3,4}){
    C2->GetPad(ican)->SetRightMargin(0.02);
    C2->GetPad(ican)->SetLeftMargin (0.15);
    C2->GetPad(ican)->SetTopMargin  (0.05);
    C2->GetPad(ican)->SetBottomMargin(0.15);
  }
  return C2;
}

TCanvas* StandardCanvas4b(std::string can_name){
  TCanvas *C2=new TCanvas(can_name.c_str(),can_name.c_str(),1200,900);
  C2->Divide(2,2);
  for(int ican:{1,2,3,4}){
    C2->GetPad(ican)->SetRightMargin(0.01);
    C2->GetPad(ican)->SetLeftMargin (0.12);
    C2->GetPad(ican)->SetTopMargin  (0.01);
    C2->GetPad(ican)->SetBottomMargin(0.15);
    C2->GetPad(ican)->SetPad(0.5*((ican-1)%2),(ican<=2)?0.5:0,0.5*((ican-1)%2+1),(ican<=2)?1.0:0.5);
  }
  return C2;
}

TCanvas* StandardCanvas4c(std::string can_name){
  TCanvas *C2=new TCanvas(can_name.c_str(),can_name.c_str(),1800,450);
  C2->Divide(4);
  for(int ican:{1,2,3,4}){
    C2->GetPad(ican)->SetRightMargin(0.01);
    C2->GetPad(ican)->SetLeftMargin (0.12);
    C2->GetPad(ican)->SetTopMargin  (0.01);
    C2->GetPad(ican)->SetBottomMargin(0.15);
    //C2->GetPad(ican)->SetPad(0.5*((ican-1)%2),(ican<=2)?0.5:0,0.5*((ican-1)%2+1),(ican<=2)?1.0:0.5);
  }
  return C2;
}
  

TCanvas* StandardCanvas6(std::string can_name){
  TCanvas *C2=new TCanvas(can_name.c_str(),can_name.c_str(),1200,1000);
  C2->Divide(2,3);
  for(int ican:{1,2,3,4,5,6}){
    C2->GetPad(ican)->SetRightMargin(0.02);
    C2->GetPad(ican)->SetLeftMargin (0.15);
    C2->GetPad(ican)->SetTopMargin  (0.02);
    C2->GetPad(ican)->SetBottomMargin(0.15);
  }
  return C2;
}

TCanvas* StandardCanvas6a(std::string can_name){
  TCanvas *C2=new TCanvas(can_name.c_str(),can_name.c_str(),800,900);
  C2->Divide(2,3);
  for(int ican:{1,2,3,4,5,6}){
    C2->GetPad(ican)->SetRightMargin(0.02);
    C2->GetPad(ican)->SetLeftMargin (0.15);
    C2->GetPad(ican)->SetTopMargin  (0.02);
    C2->GetPad(ican)->SetBottomMargin(0.15);
  }
  return C2;
}

TCanvas* StandardCanvas6b(std::string can_name){
  TCanvas *C2=new TCanvas(can_name.c_str(),can_name.c_str(),1200,900);
  C2->Divide(3,2);
  for(int ican:{1,2,3,4,5,6}){
    C2->GetPad(ican)->SetRightMargin(0.06);
    C2->GetPad(ican)->SetLeftMargin (0.18);
    C2->GetPad(ican)->SetTopMargin  (0.02);
    C2->GetPad(ican)->SetBottomMargin(0.18);
    /*
    int irow=(ican-1)/3;
    int icol=(ican-1)%3;

    C2->GetPad(ican)->SetTopMargin  (0.02);
    C2->GetPad(ican)->SetRightMargin(0.06);
    if(icol==0) C2->GetPad(ican)->SetLeftMargin(0.16);
    else        C2->GetPad(ican)->SetLeftMargin(0.06);
    if(irow==1) C2->GetPad(ican)->SetBottomMargin(0.16);
    else        C2->GetPad(ican)->SetBottomMargin(0.055);
    */

  }
  return C2;
}

TCanvas* StandardCanvas6c(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),1800,900);
  Can->Divide(3,2);
  for(int ipad:{0,1,2,3,4,5}){
    Can->GetPad(ipad+1)->SetLeftMargin  (0.12);
    Can->GetPad(ipad+1)->SetTopMargin   (0.01);
    Can->GetPad(ipad+1)->SetRightMargin (0.02);
    Can->GetPad(ipad+1)->SetBottomMargin(0.12);
  }
  return Can;
}


TCanvas*  StandardCanvas1Ratio(std::string can_name){
  TCanvas *Can=new TCanvas(can_name.c_str(),can_name.c_str(),600,450);
  Can->Divide(2);
  for(int i:{1,2}){
    TPad* Pad=(TPad*)Can->GetPad(i);
    Pad->cd();

    float Xmin=0.0,Xmax=1.0;
    float Ymin=0.50,Ymax=1.0;
    if(i==2) Ymin=0.0,Ymax=0.50;
    Pad->SetPad(Xmin,Ymin,Xmax,Ymax);

    Pad->SetLeftMargin  (0.15);
    Pad->SetRightMargin (0.01);
    Pad->SetTopMargin   (0.01);
    Pad->SetBottomMargin(0.0 );
    if(i==1) Pad->SetBottomMargin(0.01);
    else     Pad->SetBottomMargin(0.15);
  }
  return Can;
}


TCanvas* StandardCanvas9(std::string name){
  TCanvas *C1=new TCanvas(name.c_str(),name.c_str(),1000,900);
  C1->Divide(3,3);
  for(int ipad=0;ipad<9;ipad++){
    C1->cd(ipad+1);
    int irow=ipad/3;
    int icol=ipad%3;


    gPad->SetTopMargin (0.015);
    gPad->SetRightMargin(0.125 );
    if(icol==0) gPad->SetLeftMargin(0.16);
    else        gPad->SetLeftMargin(0.06);
    if(irow==2) gPad->SetBottomMargin(0.16);
    else        gPad->SetBottomMargin(0.055);

    double x1=0,x2=.36,x3=.675,x4=.99;
    double y1=0,y2=.36,y3=.675,y4=.99;
    if(ipad==0) gPad->SetPad(x1,y3,x2,y4);
    if(ipad==1) gPad->SetPad(x2,y3,x3,y4);
    if(ipad==2) gPad->SetPad(x3,y3,x4,y4);
    if(ipad==3) gPad->SetPad(x1,y2,x2,y3);
    if(ipad==4) gPad->SetPad(x2,y2,x3,y3);
    if(ipad==5) gPad->SetPad(x3,y2,x4,y3);
    if(ipad==6) gPad->SetPad(x1,y1,x2,y2);
    if(ipad==7) gPad->SetPad(x2,y1,x3,y2);
    if(ipad==8) gPad->SetPad(x3,y1,x4,y2);
  }
  return C1;
}

TCanvas* StandardCanvas9b(std::string name){
  TCanvas *C1=new TCanvas(name.c_str(),name.c_str(),1200,900);
  C1->Divide(3,3);
  for(int ipad=0;ipad<9;ipad++){
    C1->cd(ipad+1);
    int irow=ipad/3;
    int icol=ipad%3;

    gPad->SetTopMargin  (0.015);
    gPad->SetRightMargin(0.05 );
    gPad->SetLeftMargin (0.16);
  }
  return C1;
}

TCanvas* StandardCanvas9c(std::string name){
  TCanvas *C1=new TCanvas(name.c_str(),name.c_str(),1200,900);
  C1->Divide(3,3,0.001,0.001);
  for(int ipad=0;ipad<9;ipad++){
    C1->cd(ipad+1);

    gPad->SetTopMargin   (0.01);
    gPad->SetRightMargin (0.02);
    gPad->SetLeftMargin  (0.17);
    gPad->SetBottomMargin(0.12);

    //std::cout<<ipad <<"  "<< gPad->GetXlowNDC()<<std::endl;
  }
  return C1;
}


TCanvas* StandardCanvas9d(std::string name){
  TCanvas *C1=new TCanvas(name.c_str(),name.c_str(),1000,900);
  C1->Divide(3,3);
  for(int ipad=0;ipad<9;ipad++){
    C1->cd(ipad+1);
    int irow=ipad/3;
    int icol=ipad%3;

    gPad->SetTopMargin  (0.02);
    gPad->SetRightMargin(0.16 );
    if(icol==0) gPad->SetLeftMargin(0.16);
    else        gPad->SetLeftMargin(0.06);
    if(irow==2) gPad->SetBottomMargin(0.16);
    else        gPad->SetBottomMargin(0.055);

    double x1=0,x2=.36,x3=.675,x4=.99;
    double y1=0,y2=.36,y3=.675,y4=.99;
    if(ipad==0) gPad->SetPad(x1,y3,x2,y4);
    if(ipad==1) gPad->SetPad(x2,y3,x3,y4);
    if(ipad==2) gPad->SetPad(x3,y3,x4,y4);
    if(ipad==3) gPad->SetPad(x1,y2,x2,y3);
    if(ipad==4) gPad->SetPad(x2,y2,x3,y3);
    if(ipad==5) gPad->SetPad(x3,y2,x4,y3);
    if(ipad==6) gPad->SetPad(x1,y1,x2,y2);
    if(ipad==7) gPad->SetPad(x2,y1,x3,y2);
    if(ipad==8) gPad->SetPad(x3,y1,x4,y2);
  }
  return C1;
}

TCanvas* StandardCanvas9e(std::string name){
  TCanvas *C1=new TCanvas(name.c_str(),name.c_str(),1000,900);
  C1->Divide(3,3);
  for(int ipad=0;ipad<9;ipad++){
    C1->cd(ipad+1);
    int irow=ipad/3;
    int icol=ipad%3;


    gPad->SetTopMargin (0.015);
    gPad->SetRightMargin(0.125 );
    if(icol==0) gPad->SetLeftMargin(0.16);
    else        gPad->SetLeftMargin(0.10);
    if(irow==2) gPad->SetBottomMargin(0.16);
    else        gPad->SetBottomMargin(0.055);

    double x1=0,x2=.36,x3=.675,x4=.99;
    double y1=0,y2=.36,y3=.675,y4=.99;
    if(ipad==0) gPad->SetPad(x1,y3,x2,y4);
    if(ipad==1) gPad->SetPad(x2,y3,x3,y4);
    if(ipad==2) gPad->SetPad(x3,y3,x4,y4);
    if(ipad==3) gPad->SetPad(x1,y2,x2,y3);
    if(ipad==4) gPad->SetPad(x2,y2,x3,y3);
    if(ipad==5) gPad->SetPad(x3,y2,x4,y3);
    if(ipad==6) gPad->SetPad(x1,y1,x2,y2);
    if(ipad==7) gPad->SetPad(x2,y1,x3,y2);
    if(ipad==8) gPad->SetPad(x3,y1,x4,y2);
  }
  return C1;
}

TCanvas* StandardCanvas9x2(std::string name,std::string title){
  TCanvas *C1=new TCanvas(name.c_str(),title.c_str(),1000,900);
  C1->Divide(6,3);
  for(int ipad=0;ipad<18;ipad++){
    C1->cd(ipad+1);
    int irow=ipad/6;
    int icol=ipad%3;

    gPad->SetTopMargin  (0.015);
    gPad->SetRightMargin(0.125 );
    if(icol==0) gPad->SetLeftMargin(0.16);
    else        gPad->SetLeftMargin(0.06);
    if(irow==2) gPad->SetBottomMargin(0.5);
    else        gPad->SetBottomMargin(0.01);

    double x1 =0.0 ,x2 =.36,x3 =.675,x4 =.99;
    double y1 =0.0 ,y2 =.36,y3 =.675,y4 =.99;
    double y1a=0.16,y2a=.46,y3a=.775;
    double Top=0.05;

    if(ipad== 0) {gPad->SetPad(x1,y3a,x2,y4);gPad->SetBottomMargin(0.01);gPad->SetTopMargin  (Top);}
    if(ipad== 1) {gPad->SetPad(x2,y3a,x3,y4);gPad->SetBottomMargin(0.01);gPad->SetTopMargin  (Top);}
    if(ipad== 2) {gPad->SetPad(x3,y3a,x4,y4);gPad->SetBottomMargin(0.01);gPad->SetTopMargin  (Top);}

    if(ipad== 3) {gPad->SetPad(x1,y3,x2,y3a);gPad->SetTopMargin(0.01);}
    if(ipad== 4) {gPad->SetPad(x2,y3,x3,y3a);gPad->SetTopMargin(0.01);}
    if(ipad== 5) {gPad->SetPad(x3,y3,x4,y3a);gPad->SetTopMargin(0.01);}

    if(ipad== 6) {gPad->SetPad(x1,y2a,x2,y3);gPad->SetBottomMargin(0.01);gPad->SetTopMargin  (Top);}
    if(ipad== 7) {gPad->SetPad(x2,y2a,x3,y3);gPad->SetBottomMargin(0.01);gPad->SetTopMargin  (Top);}
    if(ipad== 8) {gPad->SetPad(x3,y2a,x4,y3);gPad->SetBottomMargin(0.01);gPad->SetTopMargin  (Top);}

    if(ipad== 9) {gPad->SetPad(x1,y2,x2,y2a);gPad->SetTopMargin(0.01);}
    if(ipad==10) {gPad->SetPad(x2,y2,x3,y2a);gPad->SetTopMargin(0.01);}
    if(ipad==11) {gPad->SetPad(x3,y2,x4,y2a);gPad->SetTopMargin(0.01);}

    if(ipad==12) {gPad->SetPad(x1,y1a,x2,y2);gPad->SetBottomMargin(0.01);gPad->SetTopMargin  (Top);}
    if(ipad==13) {gPad->SetPad(x2,y1a,x3,y2);gPad->SetBottomMargin(0.01);gPad->SetTopMargin  (Top);}
    if(ipad==14) {gPad->SetPad(x3,y1a,x4,y2);gPad->SetBottomMargin(0.01);gPad->SetTopMargin  (Top);}

    if(ipad==15) {gPad->SetPad(x1,y1,x2,y1a);gPad->SetTopMargin(0.01);}
    if(ipad==16) {gPad->SetPad(x2,y1,x3,y1a);gPad->SetTopMargin(0.01);}
    if(ipad==17) {gPad->SetPad(x3,y1,x4,y1a);gPad->SetTopMargin(0.01);}
  }
  return C1;
}

TCanvas* StandardCanvas12(std::string name){
  TCanvas *C1=new TCanvas(name.c_str(),name.c_str(),1300,900);
  C1->Divide(4,3);
  for(int ipad=0;ipad<12;ipad++){
    C1->cd(ipad+1);

    gPad->SetTopMargin   (0.01);
    gPad->SetRightMargin (0.02);
    gPad->SetLeftMargin  (0.17);
    gPad->SetBottomMargin(0.12);
  }
  return C1;
}
  
TCanvas* StandardCanvas16(std::string name,std::string title){
  TCanvas *C1=new TCanvas(name.c_str(),title.c_str(),900,900);
  C1->Divide(4,4);
  for(int ipad=0;ipad<16;ipad++){
    C1->cd(ipad+1);
    int irow=ipad/4;
    int icol=ipad%4;

    gPad->SetTopMargin  (0.015);
    gPad->SetRightMargin(0.025 );
    if(icol==0) gPad->SetLeftMargin(0.16);
    else        gPad->SetLeftMargin(0.10);
    if(irow==3){gPad->SetTopMargin  (0.07);gPad->SetBottomMargin(0.16);}
    else        gPad->SetBottomMargin(0.055);



    double x1=.04 +.24*icol;
    double x2=x1  +.24;
    double y1=.99 -.24*(irow+1);
    double y2=y1  +.24;
    if(icol==0) x1=0.01;
    if(irow==3) y1=0.005;
    gPad->SetPad(x1,y1,x2,y2);
  }
  return C1;
}

  
TCanvas* StandardCanvas20(std::string name,std::string title){
  TCanvas *C1=new TCanvas(name.c_str(),title.c_str(),1100,1400);
  C1->Divide(4,5);
  for(int ipad=0;ipad<20;ipad++){
    C1->cd(ipad+1);
    int irow=ipad/4;
    int icol=ipad%4;

    gPad->SetTopMargin  (0.06 );
    gPad->SetRightMargin(0.025);
    if(icol==0) gPad->SetLeftMargin(0.16);
    else        gPad->SetLeftMargin(0.10);
    if(irow==4) gPad->SetBottomMargin(0.16);
    else        gPad->SetBottomMargin(0.055);



    double x1=.04 +.24*icol;
    double x2=x1  +.24;
    double y1=.99 -.193*(irow+1);
    double y2=y1  +.193;
    if(icol==0) x1=0.01;
    if(irow==4) y1=0.005;
    gPad->SetPad(x1,y1,x2,y2);
  }
  return C1;
}


TCanvas* StandardCanvas20b(std::string name,std::string title){
  TCanvas *C1=new TCanvas(name.c_str(),title.c_str(),1100,1400);
  C1->Divide(4,5);
  for(int ipad=0;ipad<20;ipad++){
    C1->cd(ipad+1);
    int irow=ipad/4;
    int icol=ipad%4;

    gPad->SetTopMargin  (0.015);
    gPad->SetRightMargin(0.025 );
    if(icol==0) gPad->SetLeftMargin(0.16);
    else        gPad->SetLeftMargin(0.06);
    if(irow==4) gPad->SetBottomMargin(0.16);
    else        gPad->SetBottomMargin(0.055);



    double x1=.04 +.24*icol;
    double x2=x1  +.24;
    double y1=.99 -.193*(irow+1);
    double y2=y1  +.193;
    if(icol==0) x1=0.01;
    if(irow==4) y1=0.005;
    gPad->SetPad(x1,y1,x2,y2);
  }
  return C1;
}


TCanvas* StandardCanvas24(std::string name){
  TCanvas *C1=new TCanvas(name.c_str(),name.c_str(),1300,900);
  C1->Divide(6,4);
  for(int ipad=0;ipad<24;ipad++){
    C1->cd(ipad+1);
    int irow=ipad/6;
    int icol=ipad%6;

    gPad->SetTopMargin  (0.015);
    gPad->SetRightMargin(0.130);
    if(icol==0) gPad->SetLeftMargin(0.16);
    else        gPad->SetLeftMargin(0.06);
    if(irow==3) gPad->SetBottomMargin(0.16);
    else        gPad->SetBottomMargin(0.055);

    /*
    double x1=.04 +.24*icol;
    double x2=x1  +.24;
    double y1=.99 -.193*(irow+1);
    double y2=y1  +.193;
    if(icol==0) x1=0.01;
    if(irow==4) y1=0.005;
    gPad->SetPad(x1,y1,x2,y2);
    */
  }
  return C1;
}


void SaveCanvas(std::vector<TCanvas*> &can_vec,std::string base,std::string extension=".pdf")
//void SaveCanvas(std::vector<TCanvas*> &can_vec,std::string base,std::string extension=".png")
{
  char name[600];
  for(auto &can:can_vec){
    sprintf(name,"figs/%s%s%s",can->GetName(),base.c_str(),extension.c_str());
    can->SaveAs(name);
  }
  can_vec.clear();
}


void SetYError(TH1* h,float err=0.001){
  int NBins=h->GetNbinsX();
  for(int ibin=1;ibin<=NBins;ibin++){
    h->SetBinError(ibin,err);
  }
}

void ShiftXaxis(TH1* h,float shift, int excluded_lower_bins=0){
  double X[100];
  int NBins=h->GetNbinsX();
  if(NBins>100) {std::cout<<"void ShiftXaxis(TH1* h) Too many bins"<<std::endl;throw std::exception();}

  for(int ibin=0;ibin<excluded_lower_bins;ibin++){
    X[ibin]=h->GetBinLowEdge(ibin+1);
  }
  for(int ibin=excluded_lower_bins;ibin<=NBins;ibin++){
    X[ibin]=h->GetBinLowEdge(ibin+1)+shift;
  }

  h->GetXaxis()->Set(NBins,X);
}

void ScaleXaxis(TH1* h,float scale){
  double X[100];
  int NBins=h->GetNbinsX();
  if(NBins>100) {std::cout<<"void ScaleXaxis(TH1* h) Too many bins"<<std::endl;throw std::exception();}

  for(int ibin=0;ibin<=NBins;ibin++){
    X[ibin]=h->GetBinLowEdge(ibin+1)*scale;
  }

  h->GetXaxis()->Set(NBins,X);
}


//Fold 1D histogram in dphi to 0-PI range
//Only works for a particular binning
TH1D* Fold_1D(TH1 *hist){
  int NBins=hist->GetNbinsX();
  if(NBins!=36) {std::cout<<"Error in Fold()"<<std::endl;throw std::exception();}

  static TH1* hist2=0;
  if(!hist2) hist2=(TH1*)hist->Clone("hist_temp_Fold");

  char name[600];
  sprintf(name,"%s_folded",hist->GetName());
  TH1D* hist3=new TH1D(name,"",18,0.0,PI);
  hist3->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
  hist3->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());

  hist2->Reset();
  for(int ibin=1;ibin<=9;ibin++){
    double val=hist->GetBinContent(ibin);
    double err=hist->GetBinError  (ibin);
    hist2->SetBinContent(18+1-ibin,val);
    hist2->SetBinError  (18+1-ibin,err);

    val=hist->GetBinContent(18+1-ibin);
    err=hist->GetBinError  (18+1-ibin);
    hist2->SetBinContent(ibin,val);
    hist2->SetBinError  (ibin,err);

    val=hist->GetBinContent(18+ibin);
    err=hist->GetBinError  (18+ibin);
    hist2->SetBinContent(36+1-ibin,val);
    hist2->SetBinError  (36+1-ibin,err);

    val=hist->GetBinContent(36+1-ibin);
    err=hist->GetBinError  (36+1-ibin);
    hist2->SetBinContent(18+ibin,val);
    hist2->SetBinError  (18+ibin,err);
  }
  hist2->Add(hist);

  for(int ibin=1;ibin<=18;ibin++){
    double val=hist2->GetBinContent(ibin+9);
    double err=hist2->GetBinError  (ibin+9);
    hist3->SetBinContent(ibin,val);
    hist3->SetBinError  (ibin,err);
  }
  return hist3;
}


bool CheckBinLimits(const TAxis* a1, const TAxis * a2)
{
   const TArrayD * h1Array = a1->GetXbins();
   const TArrayD * h2Array = a2->GetXbins();
   Int_t fN = h1Array->fN;
   if ( fN != 0 ) {
      if ( h2Array->fN != fN ) {
         cout<<"fN="<<fN<<" "<<h2Array->fN<<std::endl;
         Exception(__LINE__,__FILE__,"different fN");
         return false;
      }
      else {
         for ( int i = 0; i < fN; ++i ) {
            // for i==fN (nbin+1) a->GetBinWidth() returns last bin width
            // we do not need to exclude that case
            double binWidth = a1->GetBinWidth(i);
            if ( ! TMath::AreEqualAbs( h1Array->GetAt(i), h2Array->GetAt(i), binWidth*1E-10 ) ) {
               std::cout<<i<<" "<<h1Array->GetAt(i)<<" "<<h2Array->GetAt(i)<<" "<<binWidth*1E-10<<std::endl; 
               Exception(__LINE__,__FILE__,"different vals");
               return false;
            }
         }
      }
   }
 
   return true;
}


//Symmetrize 1D correlation in dphi
//Only works for a particular binning
void Symmetrize_1D(TH1 *hist, bool scale_error_sqrt2=true){
  const int B1=32,B2=16,B3=8;

  int NBins=hist->GetNbinsX();
  if(NBins!=B1) {std::cout<<"Error in Symmetrize()"<<std::endl;throw std::exception();}

  TH1* hist2=(TH1*)hist->Clone("hist_temp_symmetrize");
  //CheckBinLimits(hist->GetXaxis(),hist2->GetXaxis());

  //Make a reflected clone of original histogram
  hist2->Reset();
  for(int ibin=1;ibin<=B3;ibin++){
    double val=hist->GetBinContent(ibin);
    double err=hist->GetBinError  (ibin);
    hist2->SetBinContent(B2+1-ibin,val);
    hist2->SetBinError  (B2+1-ibin,err);

    val=hist->GetBinContent(B2+1-ibin);
    err=hist->GetBinError  (B2+1-ibin);
    hist2->SetBinContent(ibin,val);
    hist2->SetBinError  (ibin,err);

    val=hist->GetBinContent(B2+ibin);
    err=hist->GetBinError  (B2+ibin);
    hist2->SetBinContent(B1+1-ibin,val);
    hist2->SetBinError  (B1+1-ibin,err);

    val=hist->GetBinContent(B1+1-ibin);
    err=hist->GetBinError  (B1+1-ibin);
    hist2->SetBinContent(B2+ibin,val);
    hist2->SetBinError  (B2+ibin,err);
  }

  //Add reflected clone to original histogram
  hist->Add(hist2);
  delete hist2;

  //Scale by 2.0 and errors by sqrt(2.0)
  if(scale_error_sqrt2==true){
    double sqrt2=sqrt(2.0);
    for(int ibin=1;ibin<=B1;ibin++){
      double val=hist->GetBinContent(ibin);
      double err=hist->GetBinError  (ibin);
      hist->SetBinContent(ibin,val/2.0);
      hist->SetBinError  (ibin,err/sqrt2);
    }
  }
  //Scale both errors and values by 2
  else hist->Scale(1/2.0);
}


//Symmetrize 2D correlation (both in deta and dphi)
//Only works for a particular binning
TH2D* Symmetrize_2D(TH2D *hist){

  int NBins_phi =hist->GetNbinsX();
  int NBins_deta=hist->GetNbinsY();
  if(NBins_phi !=36) {std::cout<<"Error in Symmetrize()"<<std::endl;throw std::exception();}
  if(NBins_deta!=50) {std::cout<<"Error in Symmetrize()"<<std::endl;throw std::exception();}

  char histname[600];
  sprintf(histname,"%s_sym",hist->GetName());
  TH2D* hist2= new TH2D(histname,histname,36,-PI/2,1.5*PI,100,-5.0,5.0);

  for(int ibin_deta=1;ibin_deta<=50;ibin_deta++){
    for(int ibin_phi=1;ibin_phi<=9;ibin_phi++){
      double val1=hist->GetBinContent(ibin_phi     ,ibin_deta);
      double err1=hist->GetBinError  (ibin_phi     ,ibin_deta);
      double val2=hist->GetBinContent(18+1-ibin_phi,ibin_deta);
      double err2=hist->GetBinError  (18+1-ibin_phi,ibin_deta);
      double val=val1+val2;
      double err=sqrt(err1*err1 + err2*err2);

      hist2->SetBinContent(ibin_phi     ,50+ibin_deta,val);
      hist2->SetBinError  (ibin_phi     ,50+ibin_deta,err);
      hist2->SetBinContent(18+1-ibin_phi,50+ibin_deta,val);
      hist2->SetBinError  (18+1-ibin_phi,50+ibin_deta,err);
      hist2->SetBinContent(ibin_phi     ,51-ibin_deta,val);
      hist2->SetBinError  (ibin_phi     ,51-ibin_deta,err);
      hist2->SetBinContent(18+1-ibin_phi,51-ibin_deta,val);
      hist2->SetBinError  (18+1-ibin_phi,51-ibin_deta,err);


      val1=hist->GetBinContent(18+  ibin_phi,ibin_deta);
      err1=hist->GetBinError  (18+  ibin_phi,ibin_deta);
      val2=hist->GetBinContent(36+1-ibin_phi,ibin_deta);
      err2=hist->GetBinError  (36+1-ibin_phi,ibin_deta);
      val=val1+val2;
      err=sqrt(err1*err1 + err2*err2);

      hist2->SetBinContent(18+  ibin_phi,50+ibin_deta,val);
      hist2->SetBinError  (18+  ibin_phi,50+ibin_deta,err);
      hist2->SetBinContent(36+1-ibin_phi,50+ibin_deta,val);
      hist2->SetBinError  (36+1-ibin_phi,50+ibin_deta,err);
      hist2->SetBinContent(18+  ibin_phi,51-ibin_deta,val);
      hist2->SetBinError  (18+  ibin_phi,51-ibin_deta,err);
      hist2->SetBinContent(36+1-ibin_phi,51-ibin_deta,val);
      hist2->SetBinError  (36+1-ibin_phi,51-ibin_deta,err);
    }
  }
  return hist2;
}



void myText(float x, float y, Color_t color, std::string text,float tsize=0.06){
  TLatex l; //l.SetTextAlign(12);
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text.c_str());
}

void myText2(float x, float y, Color_t color, std::string text,int size, int font){
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.SetTextSize(size);
  l.SetTextFont(font);
  l.DrawLatex(x,y,text.c_str());
}

void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle,std::string text,Float_t msize,Double_t tsize=0.06)
{
  TMarker *marker = new TMarker(x-(0.4*tsize),y,8);
  marker->SetMarkerColor(color);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(msize);
  marker->Draw();

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.DrawLatex(x,y,text.c_str());
}

void myLineText(Double_t x,Double_t y,Int_t color,Int_t lstyle,std::string text,Float_t lsize,Double_t tsize=0.06)
{
  TLine *line = new TLine(x-(0.65*tsize),y,x-(0.15*tsize),y);
  line->SetLineColor(color);  line->SetNDC();
  line->SetLineStyle(lstyle);
  line->SetLineWidth(lsize);
  line->Draw();

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.DrawLatex(x,y,text.c_str());
}


void format_hist(TH1* hist){
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
}

void format(TH1* hist,int col=1,int sty=20){
  //hist->GetYaxis()->CenterTitle();
  //hist->GetXaxis()->CenterTitle();

  hist->SetMarkerColor(col);
  hist->SetLineColor  (col);
  hist->SetMarkerStyle(sty);
}

void format(TGraph* hist,int col=1,int sty=20){
  //hist->GetYaxis()->CenterTitle();
  //hist->GetXaxis()->CenterTitle();

  hist->SetMarkerColor(col);
  hist->SetLineColor  (col);
  hist->SetMarkerStyle(sty);
}



std::map<std::string,double> StandardFormat(){
  std::map<std::string,double> _format={
      {"XTitleSize"  ,0.06},
      {"YTitleSize"  ,0.06},
      {"XLabelSize"  ,0.05},
      {"YLabelSize"  ,0.05},
      {"XNdivisions" ,505},
      {"YNdivisions" ,505},
      {"YTitleOffset",1.0},
      {"XTitleOffset",1.0},};
  return _format;
}
std::map<std::string,double> StandardFormat2(){
  std::map<std::string,double> _format={
      {"XTitleSize"  ,0.055},
      {"YTitleSize"  ,0.055},
      {"XLabelSize"  ,0.05},
      {"YLabelSize"  ,0.05},
      {"XNdivisions" ,505},
      {"YNdivisions" ,505},
      {"YTitleOffset",1.0},
      {"XTitleOffset",1.0},};
  return _format;
}

void FormatHist(TH1* hist1,std::map<std::string,double> f){
  if(f.find("XLabelSize"  ) !=f.end()) hist1->GetXaxis()->SetLabelSize  (f["XLabelSize"]);
  if(f.find("YLabelSize"  ) !=f.end()) hist1->GetYaxis()->SetLabelSize  (f["YLabelSize"]);
  if(f.find("XLabelOffset") !=f.end()) hist1->GetXaxis()->SetLabelOffset(f["XLabelOffset"]);
  if(f.find("YLabelOffset") !=f.end()) hist1->GetYaxis()->SetLabelOffset(f["YLabelOffset"]);
  if(f.find("XTitleOffset") !=f.end()) hist1->GetXaxis()->SetTitleOffset(f["XTitleOffset"]);
  if(f.find("YTitleOffset") !=f.end()) hist1->GetYaxis()->SetTitleOffset(f["YTitleOffset"]);
  if(f.find("XTitleSize"  ) !=f.end()) hist1->GetXaxis()->SetTitleSize  (f["XTitleSize"]);
  if(f.find("YTitleSize"  ) !=f.end()) hist1->GetYaxis()->SetTitleSize  (f["YTitleSize"]);
  if(f.find("XNdivisions")  !=f.end()) hist1->GetXaxis()->SetNdivisions (int(f["XNdivisions"]));
  if(f.find("YNdivisions")  !=f.end()) hist1->GetYaxis()->SetNdivisions (int(f["YNdivisions"]));

  //hist1->GetXaxis()->CenterTitle();
  //hist1->GetYaxis()->CenterTitle();
  hist1->SetStats(0);
  hist1->SetTitle("");
}

void FormatHist(THStack* hist1,std::map<std::string,double> f){
  if(f.find("XLabelSize"  ) !=f.end()) hist1->GetXaxis()->SetLabelSize  (f["XLabelSize"]);
  if(f.find("YLabelSize"  ) !=f.end()) hist1->GetYaxis()->SetLabelSize  (f["YLabelSize"]);
  if(f.find("XLabelOffset") !=f.end()) hist1->GetXaxis()->SetLabelOffset(f["XLabelOffset"]);
  if(f.find("YLabelOffset") !=f.end()) hist1->GetYaxis()->SetLabelOffset(f["YLabelOffset"]);
  if(f.find("XTitleOffset") !=f.end()) hist1->GetXaxis()->SetTitleOffset(f["XTitleOffset"]);
  if(f.find("YTitleOffset") !=f.end()) hist1->GetYaxis()->SetTitleOffset(f["YTitleOffset"]);
  if(f.find("XTitleSize"  ) !=f.end()) hist1->GetXaxis()->SetTitleSize  (f["XTitleSize"]);
  if(f.find("YTitleSize"  ) !=f.end()) hist1->GetYaxis()->SetTitleSize  (f["YTitleSize"]);
  if(f.find("XNdivisions")  !=f.end()) hist1->GetXaxis()->SetNdivisions (int(f["XNdivisions"]));
  if(f.find("YNdivisions")  !=f.end()) hist1->GetYaxis()->SetNdivisions (int(f["YNdivisions"]));

  //hist1->GetXaxis()->CenterTitle();
  //hist1->GetYaxis()->CenterTitle();
  hist1->SetTitle("");
}


void CopyFormat(TH1*h_target,TH1* h_source){
  h_target->SetMarkerColor(h_source->GetMarkerColor());
  h_target->SetLineColor  (h_source->GetLineColor  ());
  h_target->SetFillColor  (h_source->GetFillColor  ());
  h_target->SetMarkerSize (h_source->GetMarkerSize ());
  h_target->SetLineWidth  (h_source->GetLineWidth  ());
  h_target->SetLineStyle  (h_source->GetLineStyle  ());
  h_target->SetMarkerStyle(h_source->GetMarkerStyle());

  h_target->GetXaxis()->SetLabelSize  (h_source->GetXaxis()->GetLabelSize  ());
  h_target->GetYaxis()->SetLabelSize  (h_source->GetYaxis()->GetLabelSize  ());
  h_target->GetXaxis()->SetLabelOffset(h_source->GetXaxis()->GetLabelOffset());
  h_target->GetYaxis()->SetLabelOffset(h_source->GetYaxis()->GetLabelOffset());
  h_target->GetXaxis()->SetTitleSize  (h_source->GetXaxis()->GetTitleSize  ());
  h_target->GetYaxis()->SetTitleSize  (h_source->GetYaxis()->GetTitleSize  ());
  h_target->GetXaxis()->SetTitleOffset(h_source->GetXaxis()->GetTitleOffset());
  h_target->GetYaxis()->SetTitleOffset(h_source->GetYaxis()->GetTitleOffset());
  h_target->GetXaxis()->SetNdivisions (h_source->GetXaxis()->GetNdivisions ());
  h_target->GetYaxis()->SetNdivisions (h_source->GetYaxis()->GetNdivisions ());
}




void  MyDivide2Panels(TCanvas* Can){
  Can->Divide(2);
  for(int i:{1,2}){
    TPad* Pad=(TPad*)Can->GetPad(i);
    Pad->cd();

    float    Xmin=0.0 ,Xmax=0.55;
    if(i==2) Xmin=0.55,Xmax=1.00;
    float Ymin=0.0,Ymax=1.0;
    Pad->SetPad(Xmin,Ymin,Xmax,Ymax);

    Pad->SetLeftMargin  (0.00);
    Pad->SetRightMargin (0.00);
    Pad->SetTopMargin   (0.01);
    Pad->SetBottomMargin(0.14);
    if(i==1) Pad->SetLeftMargin  (10.0/55 );
    else     Pad->SetRightMargin (0.01);

  }
}



void MyDivide_SixPanels(TCanvas *C1){
    C1->Divide(3,2,0,0);
    float X1[6]={0.00,0.37,0.68,0.00,0.37,0.68};
    float X2[6]={0.37,0.68,0.99,0.37,0.68,0.99};
    float Y1[6]={0.53,0.53,0.53,0.00,0.00,0.00};
    float Y2[6]={0.99,0.99,0.99,0.53,0.53,0.53};

    for(int ipad=0;ipad<6;ipad++){
      TPad *pad=static_cast<TPad*>(C1->GetPad(ipad+1));
      pad->SetPad(X1[ipad],Y1[ipad],X2[ipad],Y2[ipad]);
      pad->SetTopMargin(0.01);
      if(ipad!=0 && ipad!=3) pad->SetLeftMargin  (0.0 );
      else                   pad->SetLeftMargin  (0.16);
      if(ipad==2||ipad==5)   pad->SetRightMargin (0.01);
      else                   pad->SetRightMargin (0.0 );
                             pad->SetTopMargin   (0.0 );
      if(ipad<3)             pad->SetBottomMargin(0.0 );
      else                   pad->SetBottomMargin(0.15);

    }
}
void MyDivide_SixPanels_ver2(TCanvas *C1,float BottomMargin){
    C1->Divide(3,2,0,0);
    float X1[6]={0.00,0.35,0.67,0.00,0.35,0.67};
    float X2[6]={0.35,0.67,0.99,0.35,0.67,0.99};
    float Y1[6]={0.54,0.54,0.54,0.00,0.00,0.00};
    float Y2[6]={0.99,0.99,0.99,0.54,0.54,0.54};

    for(int ipad=0;ipad<6;ipad++){
      TPad *pad=static_cast<TPad*>(C1->GetPad(ipad+1));
      pad->SetPad(X1[ipad],Y1[ipad],X2[ipad],Y2[ipad]);
      pad->SetTopMargin(0.01);
      if(ipad!=0 && ipad!=3) pad->SetLeftMargin  (0.15);
      else                   pad->SetLeftMargin  (0.22);
                             pad->SetRightMargin (0.01);
                             pad->SetTopMargin   (0.0 );
      if(ipad<3)             pad->SetBottomMargin(0.0 );
      else                   pad->SetBottomMargin(BottomMargin);

    }
}


TH1* Take_Sqrt( TH1* MyHist,int flag=0){//Symmetric errors!
  TH1* NewHist;
  if(flag){//create new one
    NewHist =(TH1*)MyHist->Clone();
    NewHist->Reset();
  }else{
    NewHist = MyHist;
  }
  int N=NewHist->GetNbinsX();
  for(int I=1;I<=N;I++){
    float val1     = MyHist->GetBinContent(I);
    float val = (val1>0)? sqrt(val1):-sqrt(-val1);
    {
      float val_err = fabs(MyHist->GetBinError  (I));
      float val_up=val1+val_err;
      float val_dn=val1-val_err;
      val_up  = (val_up>0)? sqrt(val_up):-sqrt(-val_up);
      val_dn  = (val_dn>0)? sqrt(val_dn):-sqrt(-val_dn);

      float val_err0 = fabs(val_up-val_dn)/2.0;
      float val_err1 = fabs(val_up-val);
      float val_err2 = fabs(val_dn-val);
      val_err=val_err0;
      if(val_err1>val_err) val_err=val_err1;
      if(val_err2>val_err) val_err=val_err2;

      NewHist->SetBinContent(I,val    );
      NewHist->SetBinError  (I,val_err);
    }
  }
  return NewHist;
}


TH1* TakePower(const TH1* MyHist,double power){//Symmetric errors!
  TH1* NewHist;
  NewHist =(TH1*)MyHist->Clone(UniqueName().c_str());
  if(power==1.0) return NewHist;
  NewHist->Reset();

  int N=NewHist->GetNbinsX();
  for(int I=1;I<=N;I++){
    float val     = MyHist->GetBinContent(I);
    float err     = MyHist->GetBinError  (I);

    if(val<0) continue;
    float new_val   =pow(fabs(val),power);
    float new_val_up=pow(fabs(val+err),power);
    float new_val_lo=0;
    if((val-err)<0) new_val_lo=0;
    else  new_val_lo=pow(fabs(val-err),power);

    float new_err1=(new_val_up-new_val);
    float new_err2=(new_val   -new_val_lo);
    float new_err =(new_err1>new_err2)? new_err1:new_err2;

    NewHist->SetBinContent(I,new_val);
    NewHist->SetBinError  (I,new_err);
  }
  return NewHist;
}



//Convert a TGraph to a histogram
TH1* GraphToHist(TGraph* gr, TH1*hist, std::string info=""){
  int N=hist->GetNbinsX();
  if(N>gr->GetN()){
    std::cout<<"Exception in GraphToHist(TGraphErrors* gr,TH1*hist,string info): size of TH1>TGraph "<<N<<" "<<gr->GetN()<<" :Info="<<info<<std::endl;
    throw std::exception();
  }
  Double_t *X =gr->GetX();
  Double_t *Y =gr->GetY();

  for(int ibin=1;ibin<=N;ibin++){
    hist->SetBinContent(ibin,Y [ibin-1]);

    double err=0.0;
    std::string ClassName=gr->ClassName();
    if     (ClassName=="TGraph"           ) err=0.0;
    else if(ClassName=="TGraphErrors"     ) err=(gr->GetEY())[ibin-1];
    else if(ClassName=="TGraphAsymmErrors") err=((gr->GetErrorYhigh(ibin-1)) + (gr->GetErrorYlow(ibin-1)))/2.0;
    else {std::cout<<"GraphToHist(TGraph* gr,TH1*hist,string info):: Unknown Class "<<ClassName<<" :Info="<<info<<std::endl; throw std::exception();}
    hist->SetBinError  (ibin,err);


    if(X[ibin-1]<hist->GetBinLowEdge(ibin) || X[ibin-1]>hist->GetBinLowEdge(ibin+1)){
      std::cout<<"Exception in GraphToHist(TGraphErrors* gr,TH1*hist,string info) "<<" :Info="<<info<<std::endl;
      throw std::exception();
    }
  }
  return hist;
}

//if normalize_by_width_only==true, then do not rescale histogram to unit integral
//  but instead scale by bin-widths only
//if min_x<max_x then rescale histogram to be unit integral only over this range
//  caller must be careful and ensure that min_x and max_x correspond to values at edge of bins
void Normalize(TH1*hist, bool normalize_by_width_only=false,float min_x=1,float max_x=-1){
  if(normalize_by_width_only==false){
    double integral=0;
    if(min_x> max_x) integral=hist->Integral();
    else{
      int bin1=hist->GetXaxis()->FindBin(min_x);
      int bin2=hist->GetXaxis()->FindBin(max_x);
      integral=hist->Integral(bin1,bin2);
    }
    hist->Scale(1.0/integral);
  }
  int nbins=hist->GetNbinsX();
  for(int i=1;i<=nbins;i++){
    double width=hist->GetBinWidth(i);
    double val  =hist->GetBinContent(i)/width;
    double err  =hist->GetBinError  (i)/width;
    hist->SetBinContent(i,val);
    hist->SetBinError  (i,err);
  }
}

void ScaleHistogram(TH1*hist,double val,double err){
  int nbins=hist->GetNbinsX();
  TH1* hist2=(TH1*)hist->Clone(Common::UniqueName().c_str());
  for(int ibin=1;ibin<=nbins;ibin++){
    hist2->SetBinContent(ibin,val);
    hist2->SetBinError  (ibin,err);
  }
  hist->Divide(hist2);
  delete hist2;
}

void Resample(TH1D* hist){
  static TRandom3 rnd;
  int Nbins=hist->GetNbinsX();
  for(int ibin=1;ibin<=Nbins;ibin++){
     double val=hist->GetBinContent(ibin);
     double err=hist->GetBinError(ibin);
     double resampled_val=rnd.Gaus(val,err);
     if (resampled_val<0) resampled_val=0;
     hist->SetBinContent(ibin,resampled_val);
  }
}

void Exception(int line, std::string file,std::string message){
  std::cout<<"****************************************************"<<std::endl;
  std::cout<<"Exception at Line "<<line<<" in File "<<file         <<std::endl;
  if(message!="") std::cout<<message<<std::endl;
  std::cout<<"****************************************************"<<std::endl;
  std::cerr<<"****************************************************"<<std::endl;
  std::cerr<<"Exception at Line "<<line<<" in File "<<file         <<std::endl;
  if(message!="") std::cout<<message<<std::endl;
  std::cerr<<"****************************************************"<<std::endl;
  throw std::exception();
}

TObject* GetAndClone(TFile* InFile,std::string name){
  if(!InFile           ) Common::Exception(__LINE__,__FILE__);
  if(InFile->IsZombie()) Common::Exception(__LINE__,__FILE__);
  TObject* obj=InFile->Get(name.c_str());
  if(!obj){
    char temp[600];
    sprintf(temp,"Object named %s not found in TFile %s",name.c_str(),InFile->GetName());
    Common::Exception(__LINE__,__FILE__,temp);
  }
  return obj->Clone(Common::UniqueName().c_str());
}


//Resets the 3D histogram axis ranges
void Reset3DAxes(TH3D* hist){
  hist->GetXaxis()->SetRange(1,hist->GetXaxis()->GetNbins());
  hist->GetYaxis()->SetRange(1,hist->GetYaxis()->GetNbins());
  hist->GetZaxis()->SetRange(1,hist->GetZaxis()->GetNbins());
}

//Scales 2D histograms by bin-width
void ScaleByBinWidth2D(TH2D* hist){
  int NbinsX=hist->GetNbinsX();
  int NbinsY=hist->GetNbinsY();
  for(int ibinx=1;ibinx<=NbinsX;ibinx++){
    double widthx=hist->GetXaxis()->GetBinWidth(ibinx);
    for(int ibiny=1;ibiny<=NbinsY;ibiny++){
      double widthy=hist->GetYaxis()->GetBinWidth(ibiny);
      double val=hist->GetBinContent(ibinx,ibiny)/(widthx*widthy);
      double err=hist->GetBinError  (ibinx,ibiny)/(widthx*widthy);
      hist->SetBinContent(ibinx,ibiny,val);
      hist->SetBinError  (ibinx,ibiny,err);
    }
  }
}


void Divide(TCanvas* can,int x,int y,float marx, float mary){
  double xcoor1[10][10], xcoor2[10][10], ycoor1[10][10], ycoor2[10][10];
  Double_t xlow,ylow,xup,yup;
  double ratx[]={1,1,1,1,1,1,1,1};
  double raty[]={1,1,1,1,1,1,1,1};
  double fracx[10], fracy[10];

  double xsli=0,ysli=0;//for boundary
  //define the slice size;
  ratx[0]   +=marx;
  raty[y-1] +=mary;
  for(int i=0;i<x;i++){
    xsli+=ratx[i];
  }
  for(int i=0;i<y;i++){
    ysli+=raty[i];
  }
  fracx[0]=0;  fracy[0]=1;
  for(int i=1;i<=x;i++){
    fracx[i]=fracx[i-1]+ratx[i-1]/xsli;
  }
  for(int i=1;i<=y;i++){
    fracy[i]=fracy[i-1]-raty[i-1]/ysli;
  }
  //rescale
  double scal=0.995;
  for(int i=0;i<=x;i++){
    fracx[i]= fracx[i]*scal+(1-scal)*(0.5-fracx[i]);
  }
  for(int i=0;i<=y;i++){
    fracy[i]= fracy[i]*scal+(1-scal)*(0.5-fracy[i]);
  }
  can->cd();
  can->Divide(x,y);
  int count=1;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      can->cd(count);
      count++;
      xlow = fracx[j];      xup = fracx[j+1];
      ylow = fracy[i+1];    yup = fracy[i];
      xcoor1[i][j] = xlow;      xcoor2[i][j] = xup;
      ycoor1[i][j] = ylow;      ycoor2[i][j] = yup;
      //std::cout<<xlow<<" "<<ylow<<" "<<xup<<" "<<yup<<std::endl;
      gPad->SetPad(xlow,ylow,xup,yup);
      gPad->SetLeftMargin(0.23);
      gPad->SetRightMargin(0.004);
      gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.15);
      /*      if(j==0){
        gPad->SetLeftMargin(marx/ratx[0]);
      }
      if(i==y-1){
        gPad->SetBottomMargin(mary/raty[y-1]);
        }*/

    }
  }
}

}//namespace common
#endif
