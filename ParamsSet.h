#include <vector>
#ifndef ParamsSet_h
#define ParamsSet_h

class ParamsSet{
public:
   	static const unsigned int ndRselcs=3;
   	static const unsigned int ndphiselcs=3;
   	static const unsigned int nCtrBins=5;
   	static const unsigned int nPtBins=5;
   	static const unsigned int nSigns=2;
   	static const unsigned int nGapCuts=2; //0: no gap cut; 1: having gap cut

	int scaleFactorCtrs[nCtrBins] = {1,1,2,3,3};

   	float deltaP_overP_thrsh;
   	float deltaP_overP_max;
   	float deltaP_overP_step;
   	int deltaP_overP_nbins;
   	
   	float deltaR_thrsh[ndRselcs] = {0.8,1.2,5.75};
   	float 	deltaR_step = 0.01;
   	int deltaR_nbins[ndRselcs];

   	float eta_gap_cut = 0.1;

	// // no log
   	// float minv_min[nSigns] = {0.2,1.06};
   	// float minv_max[ndRselcs] = {10,15,100};
   	// int minv_nbins[ndRselcs] = {100,150,1000};

   	float minv_max[ndRselcs] = {10,15,100};
   	// static const int minv_nbins[ndRselcs] = {25,40,200};
   	static constexpr int minv_nbins[ndRselcs] = {25,40,200};
   		// the C++ standard does not specifiy how floating point should be implemented and is left to the processor. 
   		// To get around this and other limitations constexpr was introduced.
   	float minv_logpow[nSigns][ndRselcs];
   	std::vector<float> minv_bins[nSigns][ndRselcs];

   	std::vector<std::array<float,2>> minv_cuts;

  	double PI;

  	static const int npt_bins = 50;
  	float ptlogpow = 0.02194;
  	float ptmax = 50;
  	double pTBins[npt_bins+1];
	
  	static const int npairPT_bins = 40;
  	float pairPTlogpow[nSigns][ndRselcs];
  	float pairPTmax = 40;
  	double pairPTBins[nSigns][ndRselcs][npairPT_bins+1];

  	ParamsSet();
  	~ParamsSet(){}
};

ParamsSet::ParamsSet(){
	deltaP_overP_thrsh = 0.12;
 	deltaP_overP_max = 0.12 * sqrt(2);
 	deltaP_overP_step = 0.002;
 	deltaP_overP_nbins = static_cast<int>(deltaP_overP_max/deltaP_overP_step);

	// deltaR_thrsh[ndRselcs] = {0.8,1.2,5.75};
	// deltaR_step = 0.01;

	PI=acos(-1.0);
  	for (unsigned int idr = 0; idr < ndRselcs; idr++){
    	deltaR_nbins[idr] = static_cast<int>(deltaR_thrsh[idr]/deltaR_step);
  	}

	for(int i = 0; i <= npt_bins; i++){
    	pTBins[i] = ptmax * pow(10.0, ((float)(i-npt_bins))*ptlogpow);
  	}

	pairPTlogpow[0][0] = 0.0198;
	pairPTlogpow[0][1] = 0.0198;
	pairPTlogpow[0][2] = 0.052;
	pairPTlogpow[1][0] = 0.0198;
	pairPTlogpow[1][1] = 0.0198;
	pairPTlogpow[1][2] = 0.082;


	// ={{0.0198,0.0198,0.052},{0.0198,0.0198,0.082}};
  	for (int isign = 0; isign < nSigns; isign++){
  		for (int idr = 0; idr < ndRselcs; idr++){
  			for(int ipt = 0; ipt <= npairPT_bins; ipt++){
    			pairPTBins[isign][idr][ipt] = pairPTmax * pow(10.0, ((float)(ipt - npairPT_bins)) * pairPTlogpow[isign][idr]);
  			}
  		}
  	}

   	// minv_logpow={{0.0672,0.0465,0.01342},{0.039,0.0288,0.0099}};
	minv_logpow[0][0] = 0.0672;
	minv_logpow[0][1] = 0.0465;
	minv_logpow[0][2] = 0.01342;
	minv_logpow[1][0] = 0.039;
	minv_logpow[1][1] = 0.0288;
	minv_logpow[1][2] = 0.0099;

  	for (int isign = 0; isign < nSigns; isign++){
  		for (int idr = 0; idr < ndRselcs; idr++){
  			for(int ibin = 0; ibin <= minv_nbins[idr]; ibin++){
  				minv_bins[isign][idr].push_back(minv_max[idr] * pow(10.0, ((float)(ibin - minv_nbins[idr])) * minv_logpow[isign][idr]));
			}
		}
	}

  	minv_cuts.push_back({0,1.06});
   	minv_cuts.push_back({2.9,3.3});
   	minv_cuts.push_back({3.59,3.74});
   	minv_cuts.push_back({9,9.8});

}

#endif

