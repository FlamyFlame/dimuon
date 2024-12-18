#ifndef __NoOverlayMC_h__
#define __NoOverlayMC_h__

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "vector"
#include "common.C"
#include "TGraphAsymmErrors.h"
#include "Bins.h"


class NoOverlayMC {
public:

   enum DPhi{
     NDPHI=3,
       NearSide=0,
       AwaySide=1,
       All     =2
   };

   TH1D* ProjectDphiSignif(const int ich,const int icent,const double deta_min,const double deta_max,const double dpop_min,const double dpop_max);
  
   //Used for dp/p significance
   //mean
   double dpop_mean_byetaandpt[9][8]       ={{ 0.0074,0.0075,0.0076,0.0063,0.0067,0.0047,0.0046,0.0059},
					     { 0.0168,0.0184,0.0182,0.0209,0.0184,0.0163,0.0182,0.0105},
					     { 0.0110,0.0129,0.0165,0.0136,0.0141,0.0119,0.0060,0.0074},
					     {-0.0013,0.0058,0.0112,0.0078,0.0095,0.0085,0.0102,0.0042},
					     {-0.0007,0.0036,0.0106,0.0124,0.0125,0.0144,0.0051,0.0043},
					     {-0.0084,0.0028,0.0080,0.0091,0.0120,0.0116,0.0052,0.0125},
					     { 0.0029,0.0111,0.0172,0.0195,0.0165,0.0134,0.0102,0.0077},
					     { 0.0147,0.0214,0.0275,0.0279,0.0250,0.0197,0.0185,0.0123},
					     { 0.0131,0.0120,0.0101,0.0099,0.0080,0.0081,0.0064,0.0061}};
   //rms from MC & data hybrid
   double dpop_slope_fit                [2]= {4.38e-4,-2.166e-3};
   //double dpop_slope_fit_n0p1to0p1         = -1.18e-3;
   double dpop_slope_int_byetaandcent[9][8]={{0.0641,0.0629,0.0605,0.0593,0.0589,0.0578,0.0578,0.0573},
					     {0.0675,0.0633,0.0618,0.0612,0.0614,0.0593,0.0592,0.0595},
					     {0.0618,0.0604,0.0610,0.0600,0.0590,0.0592,0.0592,0.0588},
					     {0.0578,0.0575,0.0579,0.0565,0.0566,0.0558,0.0555,0.0555},
					     {0.0679,0.0674,0.0678,0.0689,0.0673,0.0651,0.0654,0.0655},
					     {0.0579,0.0577,0.0553,0.0534,0.0543,0.0535,0.0546,0.0536},
					     {0.0629,0.0616,0.0575,0.0587,0.0574,0.0584,0.0579,0.0573},
					     {0.0674,0.0686,0.0649,0.0638,0.0629,0.0595,0.0591,0.0610},
					     {0.0652,0.0629,0.0601,0.0591,0.0577,0.0564,0.0578,0.0569}};
  
private :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain

   UInt_t          RunNumber;
   UInt_t          lbn;
   UInt_t          bcid;
   ULong64_t       eventNumber;
   Float_t         ActIntPerXing;
   Float_t         AvgIntPerXing;
   vector<int>     *trk_numqual=nullptr;
   vector<float>   *vtx_z=nullptr;
   vector<float>   *vtx_x=nullptr;
   vector<float>   *vtx_y=nullptr;
   vector<int>     *vtx_ntrk=nullptr;
   vector<float>   *muon_pt=nullptr;
   vector<float>   *muon_eta=nullptr;
   vector<float>   *muon_phi=nullptr;
   vector<int>     *muon_quality=nullptr;
   vector<float>   *muon_deltaP_overP=nullptr;
   vector<float>   *muon_d0=nullptr;
   vector<float>   *muon_d0_err=nullptr;
   vector<float>   *muon_z0=nullptr;
   vector<float>   *muon_mspt=nullptr;
   vector<float>   *muon_msp=nullptr;
   vector<int>     *muon_elosstype=nullptr;
   vector<int>     *muon_truth_index=nullptr;
   vector<float>   *muon_truth_prob=nullptr;
   vector<int>     *muon_truth_barcode=nullptr;
   vector<bool>    *muon_truth_IsPrimary=nullptr;
   vector<float>   *truth_pt=nullptr;
   vector<float>   *truth_eta=nullptr;
   vector<float>   *truth_phi=nullptr;
   vector<float>   *truth_charge=nullptr;
   vector<int>     *truth_id=nullptr;
   vector<int>     *truth_barcode=nullptr;
   vector<int>     *truth_qual=nullptr;
   vector<int>     *truth_muon_index=nullptr;
   //Int_t           centrality=1000;
   Float_t         FCal_Et=0;

   void Init();
   void InitHists(bool read_flag);

   void MakePlots(int flag);
   void PlotMatch           (TH1D*hist_             ,std::string name);
   void PlotMatch2          ();
   //void PlotMeanImpactParameters();
   //void PlotTruthAndReco    (TH1D*hist_[]           ,std::string name);
   //void PlotAcopOrAsym      (std::string type,bool logplot);
   //void PlotD0              (TH1D*hist_,std::string name,int icent_bin=-1);
   //void PlotD0MultipleCents (std::string name);
   //void PlotD0RmsWidths     (std::string name);
   void PlotEff             (int icent_bin=-1);
   //void PlotUnfoldCorrection(int icent_bin=-1);
   //void PlotFake            ();
   void PlotEff2D           ();
   //void PlotAcopAsymKperpResolutions(std::string name="Acop");


   TFile *File=nullptr;
   void ProcessSingles();//Fill histograms for efficiency and fake-rate calculations
   void MakeEffFits   ();//Fit efficiencies with TF1
   void CheckRecoEff  ();//Check effect of reconstruction eff
   void ProcessPairs  ();//Fill truth and reco Acoplanad, d0_pair etc distributions


   void RebinCentralityAndCharge();
   double CalculateDpopSig(double pt_val, double eta_val, double momimb_val, int cent_val);
  
   std::vector<TObject**> m_all_objects;

public:
   Bins::QualityCut m_muon_selection_cut;
   std::string      m_name;
   std::vector<TCanvas*> m_can_vec;

   NoOverlayMC(int read_flag         =1,
               int muon_selection_cut=static_cast<int>(Bins::QualityCut::MEDIUM   ));

   ~NoOverlayMC(){
     for(auto *obj:m_all_objects){
       if(*obj) {
         delete *obj;
         *obj=nullptr;
       }
     }
   }



   enum TruthType
   {
     NTYPE=2,
       TRUTH=0,
       RECO =1,

     NETA=9, //-2.4<q*eta<2, -2<q*eta<-1, -1<q*eta<0, 0<q*eta|<1, 1<q*eta<2, 2<q*eta<2.4
   };
   int   GetEtaBin(float eta,float charge);//Eta bin for reconstruction efficiency calculation
   std::map<int,std::string> Label   ={ {TruthType::TRUTH,"Truth"},{TruthType::RECO,"Reconstructed"}};
public:
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

   float GetReconstructionEfficiency(float eta, float pt, float charge, int centrality_percentile);//returns reconstruction efficiency

   TH1D* h_cent;
   TH1D* h_FCalEt;

   //Spectra for truth and reconstructed muons
   TH1D* h_Spectra_Reco_EffCorrected;
   TH1D *h_Spectra[TruthType::NTYPE]={nullptr};
   TH1D *h_eta    [TruthType::NTYPE]={nullptr};
   TH1D *h_phi    [TruthType::NTYPE]={nullptr};

   //Efficiency as function of truth pT in eta bins
   TH1D              *h_Eff_Num[TruthType::NETA]={nullptr};
   TH1D              *h_Eff_Den[TruthType::NETA]={nullptr};
   TGraphAsymmErrors *gr_Eff   [TruthType::NETA]={nullptr};
   //TF1 storing fit to efficiency
   TF1* tf1_eff_fit[TruthType::NETA]={nullptr};
   //2D eff hists
   TH2D* h_Eff_Num_2D=nullptr;
   TH2D* h_Eff_Den_2D=nullptr;

   //Efficiency as function of truth pT in centrality+eta bins
   //These are filled only when we run with m_mc_type!=Bins::MCType::NOOVERLAY && m_mc_type!=Bins::MCType::NOOVERLAY_2
   TH1D              *h_Eff_Num_cent[Bins::NCENT][TruthType::NETA]={{nullptr}};
   TH1D              *h_Eff_Den_cent[Bins::NCENT][TruthType::NETA]={{nullptr}};
   TGraphAsymmErrors *gr_Eff_cent   [Bins::NCENT][TruthType::NETA]={{nullptr}};
   //TF1 storing fit to efficiency
   TF1* tf1_eff_fit_cent[Bins::NCENT][TruthType::NETA]={{nullptr}};
   //2D eff hists
   TH2D* h_Eff_Num_2D_cent[Bins::NCENT]={nullptr};
   TH2D* h_Eff_Den_2D_cent[Bins::NCENT]={nullptr};

   //same as above histogrems without "_rec" at the end but store Efficiency+Unfolding correction
   TH1D              *h_Eff_Num_rec                     [TruthType::NETA]= {nullptr};
   TGraphAsymmErrors *gr_Eff_rec                        [TruthType::NETA]= {nullptr};
   TH2D              *h_Eff_Num_2D_rec                                   =  nullptr;
   TH1D              *h_Eff_Num_cent_rec   [Bins::NCENT][TruthType::NETA]={{nullptr}};
   TGraphAsymmErrors *gr_Eff_cent_rec      [Bins::NCENT][TruthType::NETA]={{nullptr}};
   TH2D              *h_Eff_Num_2D_cent_rec[Bins::NCENT]                 = {nullptr};

   //Fake-rate as function of reco pT
   TH1D              *h_Fake_Num[TruthType::NETA]={nullptr};
   TH1D              *h_Fake_Den[TruthType::NETA]={nullptr};
   TGraphAsymmErrors *gr_Fake   [TruthType::NETA]={nullptr};

   //Matching between truth and reco
   TH1D *h_Deta_Match=nullptr;
   TH1D *h_Dphi_Match=nullptr;
   TH1D *h_DpT_Match =nullptr;
   TH1D *h_DpT_MatchV2=nullptr;
   TH1D *h_Prob_Match=nullptr;
   TH1D *h_count     =nullptr;
   TH2D* h2_DpT_Match=nullptr;
   TProfile*p_Deta_Match=nullptr,*p_DpT_Match=nullptr;

   //
   TProfile *p_D0_rec=nullptr;
   TProfile *p_Z0_rec=nullptr;


   //Histograms for storing distributions to make template fits
   TH1D *h_deltaP_overP_dist_1D_singles=nullptr;
   TH1D *h_d0_dist_1D_singles          =nullptr;
   TH1D *h_deltaP_overP_dist_pt_and_eta_1D_singles[4][3]={{nullptr}};
  
   TH2D *h_deltaP_overP_dist_2D[3][3]={{nullptr}};
   TH1D *h_deltaP_overP_dist_1D[3][3]={{nullptr}};
   TH2D *h_d0_dist_2D          [3][3]={{nullptr}};
   TH1D *h_d0_dist_1D          [3][3]={{nullptr}};
   TH1D *h_d0_pair_dist_1D     [3][3]={{nullptr}};
   TH1D *h_pt1_minus_pt2_1D     =nullptr;
   TH1D *h_pt1_minus_pt2_rand_1D=nullptr;
   TH1D *h_deltaP_overP_dist_pt_and_eta_1D[4][3]={{nullptr}};

   TH1D *h_deltaP_overPpairsig_dist_1D[3][3]={{nullptr}};

   TH3D* corr_sig_signif3D     [Bins::NCENT][Bins::PairSignCut::NCH]={{nullptr}};
  
};



void NoOverlayMC::Init()
{
   fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
   fChain->SetMakeClass(1);

   fChain->Add("Data/ATLHI313/Overlay_ATLHI313.root");
   fChain->Add("Data/ATLHI304/Overlay_ATLHI304.root");
   //TODO: NOTE the no-overlay sample is added to increase statistics in the UPC bin
   fChain->Add("Data/ATLHI313/SignalOnly_ATLHI313.root");
   fChain->Add("Data/ATLHI304/SignalOnly_ATLHI304.root");

   fChain->SetBranchAddress("trk_numqual", &trk_numqual);
   fChain->SetBranchAddress("RunNumber", &RunNumber);
   fChain->SetBranchAddress("lbn", &lbn);
   fChain->SetBranchAddress("bcid", &bcid);
   fChain->SetBranchAddress("eventNumber", &eventNumber);
   fChain->SetBranchAddress("ActIntPerXing", &ActIntPerXing);
   fChain->SetBranchAddress("AvgIntPerXing", &AvgIntPerXing);
   fChain->SetBranchAddress("vtx_z", &vtx_z);
   fChain->SetBranchAddress("vtx_x", &vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y);
   fChain->SetBranchAddress("vtx_ntrk", &vtx_ntrk);
   fChain->SetBranchAddress("muon_pt", &muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi);
   fChain->SetBranchAddress("muon_quality", &muon_quality);
   fChain->SetBranchAddress("muon_deltaP_overP", &muon_deltaP_overP);
   fChain->SetBranchAddress("muon_d0", &muon_d0);
   fChain->SetBranchAddress("muon_d0_err", &muon_d0_err);
   fChain->SetBranchAddress("muon_z0", &muon_z0);
   fChain->SetBranchAddress("muon_mspt", &muon_mspt);
   fChain->SetBranchAddress("muon_msp", &muon_msp);
   fChain->SetBranchAddress("muon_elosstype", &muon_elosstype);
   fChain->SetBranchAddress("muon_truth_index", &muon_truth_index);
   fChain->SetBranchAddress("muon_truth_prob", &muon_truth_prob);
   fChain->SetBranchAddress("muon_truth_barcode", &muon_truth_barcode);
   fChain->SetBranchAddress("muon_truth_IsPrimary", &muon_truth_IsPrimary);
   fChain->SetBranchAddress("truth_pt", &truth_pt);
   fChain->SetBranchAddress("truth_eta", &truth_eta);
   fChain->SetBranchAddress("truth_phi", &truth_phi);
   fChain->SetBranchAddress("truth_charge", &truth_charge);
   fChain->SetBranchAddress("truth_id", &truth_id);
   fChain->SetBranchAddress("truth_barcode", &truth_barcode);
   fChain->SetBranchAddress("truth_qual", &truth_qual);
   fChain->SetBranchAddress("truth_muon_index", &truth_muon_index);
   //fChain->SetBranchAddress("centrality", &centrality);
   fChain->SetBranchAddress("FCal_Et", &FCal_Et);
}

#endif
