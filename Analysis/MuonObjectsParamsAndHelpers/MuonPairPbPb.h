#pragma once
#include "MuonPairReco.h"
#include <functional>

template <class Derived>
struct PairPbPbExtras {
  	int avg_centrality{};
	  int year{-1};
  	float FCal_Et{-1e6f};
  	float FCal_Et_A{-1e6f};     // FCal_Et_P in raw data (positive z = side A)
  	float FCal_Et_C{-1e6f};     // FCal_Et_N in raw data (negative z = side C)
  	float ZDC_E_tot{-1e6f};     // zdc_ZdcEnergy[0] + [1]
  	float ZDC_t_A{-1e6f};       // zdc_ZdcTime[1]
  	float ZDC_t_C{-1e6f};       // zdc_ZdcTime[0]
  	float ZDC_preamp_A{-1e6f};  // sum zdc_ZdcModulePreSampleAmp[1][0..3]
  	float ZDC_preamp_C{-1e6f};  // sum zdc_ZdcModulePreSampleAmp[0][0..3]
  	int ntrk_HIloose{-1};           // trk_numqual[2]: HI_LOOSE, pT > 400 MeV
  	int ntrk_HItight{-1};           // trk_numqual[3]: HI_TIGHT, pT > 400 MeV
  	int ntrk_HIloose_noPtCut{-1};   // trk_numqual[6]: HI_LOOSE, no pT cut
  	int ntrk_HItight_noPtCut{-1};   // trk_numqual[7]: HI_TIGHT, no pT cut
  	void UpdateCentrality();

protected:
  	static const std::vector<float> FCal_ET_Bins_PbPb2015;
  	static const std::vector<float> FCal_ET_Bins_PbPb2023;
  	static const std::vector<float> FCal_ET_Bins_PbPb2024;
  	  
  	static int GetCentrality(float fcalET, const std::vector<float> & fCal_centr_boundaries);
  	static int GetCentralityPbPb2015(float FCal_Et);
  	static int GetCentralityPbPb2023(float FCal_Et);
  	static int GetCentralityPbPb2024(float FCal_Et);

public:
	void PairValueCalcPbPb() {
      	auto& d = static_cast<Derived&>(*this);
      	d.avg_centrality = (d.m1.ev_centrality + d.m2.ev_centrality)/2;
    }
};

struct MuonPairPbPb
  : MuonPairBaseT<MuonPairPbPb, MuonPbPb>
  , PairRecoExtras<MuonPairPbPb>
  , PairDataExtras<MuonPairPbPb>
  , PairPbPbExtras<MuonPairPbPb>
{
    void PairValueCalcHook() {
        this->PairValueCalcPbPb(); // sets avg_centrality from ev_centrality branch
        // For years where the centrality branch is not filled, recalculate from FCal ET.
        // pbpb2025: centrality branch is all zeros in the skim.
        if (year % 2000 == 25) this->UpdateCentrality();
    }
};



// Glauber v3.2 based: https://cds.cern.ch/record/2672348/files/ATL-COM-PHYS-2019-340.pdf
template <class Derived>
const std::vector<float> PairPbPbExtras<Derived>::FCal_ET_Bins_PbPb2015={
  4.26372, 4.08338, 3.92041, 3.767, 3.6226, 3.48557, 3.35466, 3.22967, 3.11027, 2.99594, // 0-10
  2.88571, 2.77968, 2.6776, 2.57956, 2.48475, 2.39311, 2.30447, 2.21888, 2.13617, 2.05577, // 10-20
  1.97802, 1.90273, 1.82974, 1.75906, 1.69047, 1.62434, 1.56005, 1.49754, 1.43729, 1.37892, // 20-30
  1.32212, 1.26707, 1.21391, 1.16233, 1.1122, 1.06371, 1.01686, 0.971487, 0.927582, 0.885172, // 30-40
  0.844192, 0.804607, 0.766305, 0.729251, 0.693606, 0.659269, 0.626047, 0.59407, 0.563263, 0.533608, // 40-50
  0.505085, 0.477734, 0.451509, 0.426354, 0.402144, 0.378968, 0.356885, 0.335738, 0.315523, 0.29617, // 50-60
  0.277742, 0.260219, 0.243588, 0.227751, 0.212684, 0.198428, 0.184922, 0.172097, 0.160033, 0.148625, // 60-70
  0.137874, 0.127745, 0.118249, 0.109333, 0.100928, 0.093071, 0.085729, 0.078834, 0.072411, 0.066402 // 70-80
};

template <class Derived>
const std::vector<float> PairPbPbExtras<Derived>::FCal_ET_Bins_PbPb2023 = {
  4.51272, 4.32043, 4.15372, 3.99602, 3.84498, 3.69944, 3.55802, 3.42045, 3.28744, 3.15972, // 0-10
  3.03748, 2.92012, 2.80723, 2.69878, 2.59464, 2.49406, 2.39646, 2.3018, 2.21028, 2.12188, // 10-20
  2.03659, 1.95428, 1.87489, 1.79842, 1.72484, 1.65387, 1.58516, 1.51853, 1.45406, 1.39178, // 20-30
  1.33168, 1.27363, 1.21752, 1.16336, 1.11112, 1.06069, 1.0121, 0.965176, 0.919908, 0.876324, // 30-40
  0.834306, 0.793842, 0.754927, 0.71753, 0.681616, 0.647138, 0.61393, 0.581945, 0.55126, 0.52171, // 40-50
  0.493477, 0.466404, 0.440432, 0.415694, 0.392004, 0.369334, 0.347675, 0.327011, 0.30731, 0.288534, // 50-60
  0.270635, 0.253565, 0.237264, 0.221782, 0.207197, 0.19325, 0.179834, 0.167476, 0.155568, 0.144347, // 60-70
  0.13388, 0.12389, 0.114683, 0.105976, 0.0976472, 0.0901451, 0.082643, 0.0758922, 0.0695501, 0.063208, // 70-80
  0.0575959, 0.052731, 0.0478661, 0.0430012, 0.0388482 // 80-85
};

template <class Derived>
const std::vector<float> PairPbPbExtras<Derived>::FCal_ET_Bins_PbPb2024 = { // update to agree with 2023 final
  4.51272, 4.32043, 4.15372, 3.99602, 3.84498, 3.69944, 3.55802, 3.42045, 3.28744, 3.15972, // 0-10
  3.03748, 2.92012, 2.80723, 2.69878, 2.59464, 2.49406, 2.39646, 2.3018, 2.21028, 2.12188, // 10-20
  2.03659, 1.95428, 1.87489, 1.79842, 1.72484, 1.65387, 1.58516, 1.51853, 1.45406, 1.39178, // 20-30
  1.33168, 1.27363, 1.21752, 1.16336, 1.11112, 1.06069, 1.0121, 0.965176, 0.919908, 0.876324, // 30-40
  0.834306, 0.793842, 0.754927, 0.71753, 0.681616, 0.647138, 0.61393, 0.581945, 0.55126, 0.52171, // 40-50
  0.493477, 0.466404, 0.440432, 0.415694, 0.392004, 0.369334, 0.347675, 0.327011, 0.30731, 0.288534, // 50-60
  0.270635, 0.253565, 0.237264, 0.221782, 0.207197, 0.19325, 0.179834, 0.167476, 0.155568, 0.144347, // 60-70
  0.13388, 0.12389, 0.114683, 0.105976, 0.0976472, 0.0901451, 0.082643, 0.0758922, 0.0695501, 0.063208, // 70-80
  0.0575959, 0.052731, 0.0478661, 0.0430012, 0.0388482 // 80-85
};

template <class Derived>
int PairPbPbExtras<Derived>::GetCentrality(float fcalET, const std::vector<float> & fCal_centr_boundaries) {
  // fCal_centr_boundaries is sorted descending
  // We want the first boundary that is <= fcalET.
  // That is exactly what lower_bound with std::greater does:
  auto it = std::lower_bound(fCal_centr_boundaries.begin(),
                             fCal_centr_boundaries.end(),
                             fcalET,
                             std::greater<float>());

  // centrality is the 0-based index
  int centrality = it - fCal_centr_boundaries.begin();

  if (centrality < 0) centrality = 0;
  if (centrality >= (int)fCal_centr_boundaries.size()) { // if centrality >= 85%, return -1
    return -1;
  }

  return centrality;
}

template <class Derived>
int PairPbPbExtras<Derived>::GetCentralityPbPb2015(float FCal_Et){
  return GetCentrality(FCal_Et, FCal_ET_Bins_PbPb2015);
}

template <class Derived>
int PairPbPbExtras<Derived>::GetCentralityPbPb2023(float FCal_Et){
  return GetCentrality(FCal_Et, FCal_ET_Bins_PbPb2023);
}

template <class Derived>
int PairPbPbExtras<Derived>::GetCentralityPbPb2024(float FCal_Et){
  return GetCentrality(FCal_Et, FCal_ET_Bins_PbPb2024);
}

template <class Derived>
void PairPbPbExtras<Derived>::UpdateCentrality(){
  switch (year % 2000){
  case 15:
    avg_centrality = GetCentralityPbPb2015(FCal_Et);
    break;
  case 18: // use 2015 centralities
    avg_centrality = GetCentralityPbPb2015(FCal_Et);
    break;
  case 23:
    avg_centrality = GetCentralityPbPb2023(FCal_Et);
    break;
  case 24:
    avg_centrality = GetCentralityPbPb2024(FCal_Et);
    break;
  case 25:
    avg_centrality = GetCentralityPbPb2023(FCal_Et);  // use pbpb2023 thresholds until pbpb2025 are derived
    break;
  default:
    std::cout << "PairPbPbExtras::UpdateCentrality:    WARNING:: UpdateCentrality called but year is INVALID (must be 2015 / 2023 / 2024 / 2025)" << std::endl;
    std::cout << "Centrality is NOT updated!" << std::endl;
  }
}

