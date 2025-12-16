#include "DimuonAnalysisBaseClass.h"
#include <memory>
#include <stdexcept>
#include <string>

void DimuonAnalysisBaseClass::InitOutput(){
    cout << "DimuonAnalysisBaseClass::InitOutput: [DEBUG] numCuts = " << numCuts << endl;
    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        h_cutAcceptance[isign] = new TH1D(Form("h_cutAcceptance_sign%d",isign+1),Form("h_cutAcceptance_sign%d",isign+1),numCuts,0,numCuts);
    }
}

void DimuonAnalysisBaseClass::FillSingleMuonTree(){
  muonOutTree->Fill();
}

void DimuonAnalysisBaseClass::ResonanceTaggingV2(std::shared_ptr<MuonPair> const& mpair){
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

void DimuonAnalysisBaseClass::ResonanceTagging(std::shared_ptr<MuonPair> const& mpair){
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


bool DimuonAnalysisBaseClass::IsPhotoProduction(const std::shared_ptr<MuonPair>& mpair){
  return (!(mpair->same_sign) && mpair->asym < 0.05 && mpair->acop < 0.01);
}


void DimuonAnalysisBaseClass::HistAdjust(){
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

}
