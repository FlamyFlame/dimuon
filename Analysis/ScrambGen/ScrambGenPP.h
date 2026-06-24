#ifndef ScrambGenPP_h
#define ScrambGenPP_h
// Mixed-event combinatoric pair generator (pp), modern object model. Same idea as the PbPb
// ScrambGen but pp has NO centrality binning (single pool) and uses MuonPP / MuonPairPP.
// Reads single-muon trees (TTree `muon_tree`, branch `MuonObj` of type MuonPP, one entry =
// one muon), forms UNCORRELATED pairs by mixing two muons from DIFFERENT events, and writes
// muon_pair_tree_sign1 (SS) / muon_pair_tree_sign2 (OS) with a `MuonPairObj` (MuonPairPP)
// branch. NO physics cuts here (the RDF signal_cuts apply them).
#include "../MuonObjectsParamsAndHelpers/MuonPairReco.h"   // MuonPairPP; pulls Muon.h
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <vector>
#include <string>
#include <iostream>

class ScrambGenPP {
public:
   int  oversample = 5;
   UInt_t rng_seed = 20260623;
   ScrambGenPP() {}
   ~ScrambGenPP() {}
   void Run();

private:
   TRandom3 rng;
   static std::string DataDir(){ return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"; }
   static const int kNParts = 12;
   static std::string OutputPath(){ return DataDir() + "muon_pairs_pp_2024_2mu4_scrambled.root"; }

   void LoadMuons(std::vector<MuonPP>& muons){
      muons.clear();
      TChain ch("muon_tree");
      for (int p = 1; p <= kNParts; ++p){
         ch.Add((DataDir() + "single_muon_trees_pp_2024_part" + std::to_string(p) +
                 "_2mu4_mindR_0_02.root").c_str());
      }
      MuonPP* m = nullptr;
      ch.SetBranchAddress("MuonObj", &m);
      const Long64_t nentries = ch.GetEntries();
      for (Long64_t i = 0; i < nentries; ++i){
         ch.GetEntry(i);
         if (!m) continue;
         muons.push_back(*m);
      }
   }

   void GeneratePairs(std::vector<MuonPP>& muons){
      TFile* out = new TFile(OutputPath().c_str(), "recreate");
      if (!out || out->IsZombie()){ std::cerr << "[ScrambGenPP] cannot open output " << OutputPath() << std::endl; return; }
      MuonPairPP* pair = new MuonPairPP();
      TTree* t_ss = new TTree("muon_pair_tree_sign1", "scrambled same-sign (SS) pairs");
      TTree* t_op = new TTree("muon_pair_tree_sign2", "scrambled opposite-sign (OS) pairs");
      t_ss->Branch("MuonPairObj", &pair);
      t_op->Branch("MuonPairObj", &pair);

      const int nm = (int)muons.size();
      long n_ss = 0, n_op = 0;
      if (nm >= 2){
         const long target = (long)oversample * nm;
         for (long k = 0; k < target; ++k){
            const int i = (int)rng.Integer(nm);
            int j = i, tries = 0;
            do { j = (int)rng.Integer(nm); } while (muons[j].ev_num == muons[i].ev_num && ++tries < 50);
            if (muons[j].ev_num == muons[i].ev_num) continue;
            *pair = MuonPairPP();
            pair->m1 = muons[i];
            pair->m2 = muons[j];
            pair->weight = 1.0;
            pair->Update();   // minv, pair_pt/eta, dr, same_sign
            if (pair->same_sign){ t_ss->Fill(); ++n_ss; }
            else                { t_op->Fill(); ++n_op; }
         }
      }
      out->cd();
      t_ss->Write();
      t_op->Write();
      std::cout << "[ScrambGenPP] SS=" << n_ss << "  OS=" << n_op
                << "  -> " << OutputPath() << std::endl;
      out->Close();
      delete pair;
   }
};

inline void ScrambGenPP::Run(){
   rng.SetSeed(rng_seed);
   std::vector<MuonPP> muons;
   LoadMuons(muons);
   std::cout << "[ScrambGenPP] loaded " << muons.size() << " muons" << std::endl;
   GeneratePairs(muons);
}

#endif
