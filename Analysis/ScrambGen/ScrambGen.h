#ifndef ScrambGen_h
#define ScrambGen_h
// Mixed-event combinatoric pair generator (PbPb), modern object model.
// Reads single-muon trees (one TTree `muon_tree`, one branch `MuonObj` of type MuonPbPb,
// one entry = one muon), bins muons by ev_centrality into ParamsSet::nCtrIntvls 5%-wide
// intervals, and forms UNCORRELATED pairs by mixing two muons from DIFFERENT events within
// the SAME interval. Writes muon_pair_tree_sign1 (SS) / muon_pair_tree_sign2 (OS) with a
// `MuonPairObj` (MuonPairPbPb) branch, so the existing RDF reads it as df_ss / df_op and the
// RDF signal_cuts build the mixed-event minv template T_mix.
// NO resonance / photoproduction / dR / minv cuts here (those live in the RDF signal_cuts).
// See docs/tracking/low_mass_dimuon_template_fit.md (Design Decision "ScrambGen rewrite").
#include "../MuonObjectsParamsAndHelpers/MuonPairPbPb.h"   // MuonPbPb, MuonPairPbPb, MuonPairBase/Reco
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"       // ParamsSet::nCtrIntvls, CtrStep
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <vector>
#include <string>
#include <iostream>

class ScrambGen {
public:
   int  oversample = 5;          // N_pairs(interval) = oversample * N_muons(interval); shape-only (N_C floats in fit)
   UInt_t rng_seed = 20260623;   // fixed seed -> reproducible scrambled sample
   ScrambGen() {}
   ~ScrambGen() {}
   void Run(int run_year);       // 23 / 24 / 25

private:
   TRandom3 rng;

   static std::string DataDir(int yr){
      return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + std::to_string(yr) + "/";
   }
   static int NParts(int yr){ return (yr == 23) ? 4 : (yr == 24) ? 2 : 6; }
   static std::string OutputPath(int yr){
      return DataDir(yr) + "muon_pairs_pbpb_20" + std::to_string(yr) + "_single_mu4_scrambled.root";
   }

   // Load all muons from the per-part single-muon trees into per-interval vectors.
   void LoadMuons(int yr, std::vector<std::vector<MuonPbPb>>& by_ctr){
      by_ctr.assign(ParamsSet::nCtrIntvls, {});
      TChain ch("muon_tree");
      for (int p = 1; p <= NParts(yr); ++p){
         ch.Add((DataDir(yr) + "single_muon_trees_pbpb_20" + std::to_string(yr) +
                 "_part" + std::to_string(p) + "_single_mu4_mindR_0_02.root").c_str());
      }
      MuonPbPb* m = nullptr;
      ch.SetBranchAddress("MuonObj", &m);
      const Long64_t nentries = ch.GetEntries();
      for (Long64_t i = 0; i < nentries; ++i){
         ch.GetEntry(i);
         if (!m) continue;
         const int c = m->ev_centrality;
         if (c < 0) continue;                                   // exclude >85% / invalid
         const int idx = c / (int)ParamsSet::CtrStep;           // 5%-wide interval index
         if (idx < 0 || idx >= (int)ParamsSet::nCtrIntvls) continue;
         by_ctr[idx].push_back(*m);                             // copy the muon object
      }
   }

   // Mix muons within each interval and write the SS/OS scrambled pair trees.
   void GeneratePairs(int yr, std::vector<std::vector<MuonPbPb>>& by_ctr){
      TFile* out = new TFile(OutputPath(yr).c_str(), "recreate");
      if (!out || out->IsZombie()){ std::cerr << "[ScrambGen] cannot open output " << OutputPath(yr) << std::endl; return; }
      MuonPairPbPb* pair = new MuonPairPbPb();
      TTree* t_ss = new TTree("muon_pair_tree_sign1", "scrambled same-sign (SS) pairs");
      TTree* t_op = new TTree("muon_pair_tree_sign2", "scrambled opposite-sign (OS) pairs");
      t_ss->Branch("MuonPairObj", &pair);
      t_op->Branch("MuonPairObj", &pair);

      long n_ss = 0, n_op = 0;
      for (int c = 0; c < (int)ParamsSet::nCtrIntvls; ++c){
         std::vector<MuonPbPb>& mu = by_ctr[c];
         const int nm = (int)mu.size();
         if (nm < 2) continue;
         const long target = (long)oversample * nm;
         for (long k = 0; k < target; ++k){
            const int i = (int)rng.Integer(nm);
            int j = i, tries = 0;
            do { j = (int)rng.Integer(nm); } while (mu[j].ev_num == mu[i].ev_num && ++tries < 50);
            if (mu[j].ev_num == mu[i].ev_num) continue;          // could not find a different-event partner
            *pair = MuonPairPbPb();                              // reset
            pair->m1 = mu[i];
            pair->m2 = mu[j];
            pair->year   = yr;
            pair->weight = 1.0;
            pair->Update();   // minv, pair_pt/eta, dr, same_sign
            // Set the mixed-pair centrality label EXPLICITLY from the two binned muons. Both are in the
            // SAME 5% interval so the average is the correct label. This is REQUIRED: for year==25
            // Update()->PairValueCalcHook calls UpdateCentrality() which recomputes centrality from the
            // PAIR-level FCal_Et — which a mixed pair (two different events) does not have (stays -1e6 ->
            // GetCentrality returns -1), so without this override ALL yr25 pairs get avg_centrality=-1 and
            // are dropped by the RDF centrality filter. (yr23/24 take the ev_centrality-average path, so this
            // is a no-op there — same value.)
            pair->avg_centrality = (mu[i].ev_centrality + mu[j].ev_centrality) / 2;
            if (pair->same_sign){ t_ss->Fill(); ++n_ss; }
            else                { t_op->Fill(); ++n_op; }
         }
      }
      out->cd();
      t_ss->Write();
      t_op->Write();
      std::cout << "[ScrambGen yr" << yr << "] SS=" << n_ss << "  OS=" << n_op
                << "  -> " << OutputPath(yr) << std::endl;
      out->Close();
      delete pair;
   }
};

inline void ScrambGen::Run(int run_year){
   rng.SetSeed(rng_seed);
   std::vector<std::vector<MuonPbPb>> by_ctr;
   LoadMuons(run_year, by_ctr);
   long tot = 0;
   for (auto& v : by_ctr) tot += (long)v.size();
   std::cout << "[ScrambGen yr" << run_year << "] loaded " << tot << " muons into "
             << ParamsSet::nCtrIntvls << " centrality intervals (5% each)" << std::endl;
   GeneratePairs(run_year, by_ctr);
}

#endif
