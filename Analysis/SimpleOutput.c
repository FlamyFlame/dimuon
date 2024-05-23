#ifndef __SimpleOutput_C__
#define __SimpleOutput_C__

#include "SimpleOutput.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"
#include<bits/stdc++.h>
#include <stdio.h>
#include <algorithm>

void SimpleOutput::ProcessData(){

  Long64_t nentries = fChain->GetEntries();//number of events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {//loop over the events
    if (std::find(events_of_interest.begin(),events_of_interest.end(),jentry) == events_of_interest.end())
      continue;
    // std::cout << jentry << std::endl;
    int num_bytes = fChain->GetEntry(jentry);//read in an event
    if(num_bytes==0){
      std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
      throw std::exception();
    }

    int nparticles = truth_pt->size();//number of muon pairs in the event

    for(int i=0;i<nparticles;i++){//loop over all muon-pairs in the event

      int ev_num = jentry;
      int pid = truth_id->at(i);
      int pbarcode = truth_barcode->at(i);
      std::vector<int> pparents = truth_parents->at(i);
      std::vector<int> pchildren = truth_children->at(i);
      float ppt = truth_pt->at(i)/1000.;
      float peta = truth_eta->at(i);
      float pphi = truth_phi->at(i);
      // float weight = EventWeights;

      // outfile << ev_num << "\t" << pbarcode << "\t" << pid << "\t" << ppt << "\t" << peta << "\t" << "{{";
      // for (auto iparent : pparents) outfile << iparent << " ";
      // outfile << "}}" << std::endl;
      fprintf(outfile, "%4d%6d%6d%12.2e%12.2f%12.2f  {{",ev_num,pbarcode,pid,ppt,peta,pphi);
      // fprintf(outfile, "%4d%4d%6d%8.3f%8.2f  {{",ev_num,pbarcode,pid,ppt,peta);
      for (auto iparent : pparents) fprintf(outfile,"%d ",iparent);
      fprintf(outfile, "}}\t{{");
      for (auto ichild : pchildren) fprintf(outfile,"%d ",ichild);
      fprintf(outfile,"}}\n");
    }
  }
}


void SimpleOutput::Run(){
  // if (mc_mode == "mc_truth_bb") filter_effcy = filter_effcy_bb;
  // else if (mc_mode == "mc_truth_cc") filter_effcy = filter_effcy_cc;
  // std::cout << "The MC mode is " << mc_mode << ". Filter efficiency is " << filter_effcy << std::endl;

  // std::cout << "Mode = " << mode << ". Output file is " << m_outfile << std::endl;
  InitInput();
  std::string fname = (mc_mode == "mc_truth_bb")? "truth_print_bb.txt" : "truth_print_cc.txt";
  events_of_interest = (mc_mode == "mc_truth_bb")? bb_events_of_interest : cc_events_of_interest;
  
  // std::cout << events_of_interest.size() << " " << events_of_interest[0] << " " << events_of_interest[-1] << std::endl;
  outfile = fopen ((mcdir + fname).c_str(),"w");
  fprintf(outfile,"#ev\tbarcode\tid\tpt\t\teta\t\tparents\t\tchildren\n");

  ProcessData();
  fclose(outfile);
}

#endif
