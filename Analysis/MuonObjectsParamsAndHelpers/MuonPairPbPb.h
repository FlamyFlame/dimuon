#include "MuonPairData.h"

class MuonPairPbPb : public MuonPairData{ 
private:
  static const std::vector<float> FCal_ET_Bins_PbPb2015;


  static const std::vector<float> FCal_ET_Bins_PbPb2023;

  static const std::vector<float> FCal_ET_Bins_PbPb2024;
    

  static int GetCentrality(float fcalET, const std::vector<float> & fCal_centr_boundaries);
  static int GetCentralityPbPb2015(float FCal_Et);
  static int GetCentralityPbPb2023(float FCal_Et);
  static int GetCentralityPbPb2024(float FCal_Et);

public:
  int year;
  float FCal_Et;
  void UpdateCentrality();

};


// Glauber v3.2 based: https://cds.cern.ch/record/2672348/files/ATL-COM-PHYS-2019-340.pdf
const std::vector<float> MuonPairPbPb::FCal_ET_Bins_PbPb2015={
  4.26372, 4.08338, 3.92041, 3.767, 3.6226, 3.48557, 3.35466, 3.22967, 3.11027, 2.99594, // 0-10
  2.88571, 2.77968, 2.6776, 2.57956, 2.48475, 2.39311, 2.30447, 2.21888, 2.13617, 2.05577, // 10-20
  1.97802, 1.90273, 1.82974, 1.75906, 1.69047, 1.62434, 1.56005, 1.49754, 1.43729, 1.37892, // 20-30
  1.32212, 1.26707, 1.21391, 1.16233, 1.1122, 1.06371, 1.01686, 0.971487, 0.927582, 0.885172, // 30-40
  0.844192, 0.804607, 0.766305, 0.729251, 0.693606, 0.659269, 0.626047, 0.59407, 0.563263, 0.533608, // 40-50
  0.505085, 0.477734, 0.451509, 0.426354, 0.402144, 0.378968, 0.356885, 0.335738, 0.315523, 0.29617, // 50-60
  0.277742, 0.260219, 0.243588, 0.227751, 0.212684, 0.198428, 0.184922, 0.172097, 0.160033, 0.148625, // 60-70
  0.137874, 0.127745, 0.118249, 0.109333, 0.100928, 0.093071, 0.085729, 0.078834, 0.072411, 0.066402 // 70-80
};

const std::vector<float> MuonPairPbPb::FCal_ET_Bins_PbPb2023 = {
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

const std::vector<float> MuonPairPbPb::FCal_ET_Bins_PbPb2024 = {
  4.52146, 4.31811, 4.13353, 3.96046, 3.79717, 3.64263, 3.49592, 3.35557, 3.22165, 3.09343, // 0-10
  2.97016, 2.8517, 2.73778, 2.62777, 2.52149, 2.41863, 2.3191, 2.22271, 2.12924, 2.03863, // 10-20
  1.95086, 1.86565, 1.78306, 1.70296, 1.62528, 1.55017, 1.47738, 1.40701, 1.33893, 1.27333, // 20-30
  1.21019, 1.14936, 1.09078, 1.0345, 0.980452, 0.928742, 0.879123, 0.831534, 0.786079, 0.74255, // 30-40
  0.700943, 0.661196, 0.623259, 0.587086, 0.552536, 0.519523, 0.488175, 0.45821, 0.429779, 0.402687, // 40-50
  0.377083, 0.352811, 0.329731, 0.307802, 0.287017, 0.26733, 0.248698, 0.231103, 0.214482, 0.198731, // 50-60
  0.184059, 0.170295, 0.157174, 0.145029, 0.133661, 0.122843, 0.113043, 0.103605, 0.0950586, 0.0870635, // 60-70
  0.0792377, 0.0726958, 0.066154, 0.059703, 0.0546939, 0.0496847, 0.0446756, 0.0397848, 0.0365526, 0.0333205, // 70-80
  0.0300883, // 80-81
};

int MuonPairPbPb::GetCentrality(float fcalET, const std::vector<float> & fCal_centr_boundaries) {
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


int MuonPairPbPb::GetCentralityPbPb2015(float FCal_Et){
  return GetCentrality(FCal_Et, FCal_ET_Bins_PbPb2015);
}
int MuonPairPbPb::GetCentralityPbPb2023(float FCal_Et){
  return GetCentrality(FCal_Et, FCal_ET_Bins_PbPb2023);
}
int MuonPairPbPb::GetCentralityPbPb2024(float FCal_Et){
  return GetCentrality(FCal_Et, FCal_ET_Bins_PbPb2024);
}

void MuonPairPbPb::UpdateCentrality(){
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
  default:
    std::cout << "MuonPairPbPb::UpdateCentrality:    WARNING:: UpdateCentrality called but year is INVALID (must be 2015 / 2023 / 2024)" << std::endl;
    std::cout << "Centrality is NOT updated!" << std::endl;
  }
}

