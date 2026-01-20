#include "DimuonAlgCoreT.h"
#include <memory>
#include <stdexcept>
#include <string>

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::Run(){
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    Initialize();
    ProcessData();
    HistAdjust();
    Finalize();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << "CPU time used (in seconds) is " << cpu_time_used << std::endl;
}

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::Initialize(){
    
    mpair = std::make_shared<pair_t>();
    mpair_raw_ptr = mpair.get();

    self().InitParamsHook();
    self().InitInputHook();
    InitOutput();
    self().InitTempVariablesHook();
    self().InitializeExtraHook();
}

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::InitOutput(){
    self().InitOutputSettingsHook();
    InitOutputTrees();
    InitOutputHists();
}

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::InitOutputTrees(){
    self().OutputTreePath(); // including output hist file path
    if (!m_outfile) throw std::runtime_error("Output file is nullptr!");

    // Create muon/muon-pair output file
    m_outfile = new TFile(output_file_path.c_str(), "recreate");

    // Create main trees
    if (output_single_muon_tree) {
        muonOutTree = new TTree("muon_tree", "All single muons");
        muonOutTree->Branch("MuonObj", &muon_raw_ptr);
    } else {
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++) {
            muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u", ksign+1), Form("All muon pairs, sign%u", ksign+1));
            muonPairOutTree[ksign]->Branch("MuonPairObj", &mpair_raw_ptr);
        }
    }
    self().InitOutputTreesExtraHook();
}

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::InitOutputHists(){
    self().OutputHistPath(); // including output hist file path
    if (!m_outHistFile) throw std::runtime_error("Output hist file is nullptr!");

    cout << "DimuonAlgCoreT::InitOutputHists: [DEBUG] numCuts = " << numCuts << endl;
    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        h_cutAcceptance[isign] = new TH1D(Form("h_cutAcceptance_sign%d",isign+1),Form("h_cutAcceptance_sign%d",isign+1),numCuts,0,numCuts);
    }
    self().InitOutputHistsExtraHook();
}


template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::FillSingleMuonTree(){
    muonOutTree->Fill();
    self().FillSingleMuonTreeHook();
}

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::FillMuonPairTree(){
    int nsign = (mpair->same_sign)? 0:1;

    try{
      if (!mpair) throw std::runtime_error("FillMuonPairTree: Muon Pair doesn't exist!");
      if (!muonPairOutTree[nsign]) throw std::runtime_error("FillMuonPairTree: Muon Pair output tree doesn't exist!");
      mpair_raw_ptr = mpair.get();
    }catch(const std::runtime_error& e){
      std::cout << "Runtime error caught in function FillMuonPairTree: " << e.what() << std::endl;
      std::cout << "Return without filling the muon pair in the output trees!" << std::endl;
      return;
    }
    muonPairOutTree[nsign]->Fill();
    self().FillMuonPairTreeHook();
}

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::ResonanceTaggingV2(){
    // function that checks if a muon pair comes from the resonance
  
    if (!(mpair->same_sign)){ // opposite sign
        if (mpair->minv > pms.minv_upper){ // if minv is too large - treat as a "resonance" and tag both
            resonance_tagged_muon_index_list_v2.push_back(mpair->m1.ind);
            resonance_tagged_muon_index_list_v2.push_back(mpair->m2.ind);
            return;
        }

        for (array<float,2> ires : pms.minv_cuts_v2){
            if (mpair->minv > ires[0] && mpair->minv < ires[1]){ // if minv fits within a resonance-cut range: tag both as from resonance
                resonance_tagged_muon_index_list_v2.push_back(mpair->m1.ind);
                resonance_tagged_muon_index_list_v2.push_back(mpair->m2.ind);
                return;
            }
        }
    }
    
    return; // not resonance
}

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::ResonanceTagging(){
    // function that checks if a muon pair comes from the resonance
  
    if (!(mpair->same_sign)){ // opposite sign
        if (mpair->minv > pms.minv_upper){ // if minv is too large - treat as a "resonance" and tag both
            resonance_tagged_muon_index_list.push_back(mpair->m1.ind);
            resonance_tagged_muon_index_list.push_back(mpair->m2.ind);
            return;
        }

        for (array<float,2> ires : pms.minv_cuts){
            if (mpair->minv > ires[0] && mpair->minv < ires[1]){ // if minv fits within a resonance-cut range: tag both as from resonance
                resonance_tagged_muon_index_list.push_back(mpair->m1.ind);
                resonance_tagged_muon_index_list.push_back(mpair->m2.ind);
                return;
            }
        }
    }
    
    return; // not resonance
}

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::HistAdjust(){
    // set labels for the cut-acceptance histograms
    try{
        for (int ibin = 0; ibin < numCuts; ibin++){
            h_cutAcceptance[0]->GetXaxis()->SetBinLabel(ibin+1,cutLabels.at(ibin).c_str());
            h_cutAcceptance[1]->GetXaxis()->SetBinLabel(ibin+1,cutLabels.at(ibin).c_str());
        }     
    }catch (const std::out_of_range& e){
        std::cerr << "out-of-range error occured for the cut labels: " << e.what() << endl;
        std::cout << "Return without setting labels for the cut-acceptance histograms!" << std::endl;
        return;
    }catch (...){
        std::cerr << "Some other unknown error occured when setting labels for the cut-acceptance histograms" << endl;
        std::cout << "Return without setting labels for the cut-acceptance histograms!" << std::endl;
        return;
    }
  
    // normalize the cut-acceptance histograms
    try{
        if (h_cutAcceptance[0]->GetBinContent(1) == 0)   
            throw std::runtime_error("The first bin (#events with no cut applied) of the same-sign cut-acceptance histogram has ZERO bin content! Cut acceptance histogram erroneous.");
        if (h_cutAcceptance[1]->GetBinContent(1) == 0)   
            throw std::runtime_error("The first bin (#events with no cut applied) of the opposite-sign cut-acceptance histogram has ZERO bin content! Cut acceptance histogram erroneous.");
        
        h_cutAcceptance[0]->Scale(1./h_cutAcceptance[0]->GetBinContent(1));
        h_cutAcceptance[1]->Scale(1./h_cutAcceptance[1]->GetBinContent(1));
    }catch (const std::runtime_error & e){
        std::cout << "Runtime error caught: " << e.what() << std::endl;
        std::cout << "Return without normalizing the cut-acceptance histograms!" << std::endl;
        return;
    }catch (...){
        std::cerr << "Some other unknown error occured when normalizing the cut-acceptance histograms" << endl;
    }

    self().HistAdjustHook();
}

template <class PairT, class MuonT, class Derived>
void DimuonAlgCoreT<PairT, MuonT, Derived>::Finalize(){
    m_outfile->Write();
    m_outHistFile->Write();
    m_outfile->Close();
    m_outHistFile->Close();
    std::cout << "Output muon-pair trees have been written to: " << output_file_path << std::endl;
    std::cout << "Output histograms have been written to: " << output_hist_file_path << std::endl;

    self().FinalizeHook();
}