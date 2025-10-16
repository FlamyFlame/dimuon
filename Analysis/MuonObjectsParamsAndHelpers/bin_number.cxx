// helper function that for a given target (edge of region, e.g, for projection)
// and a vector of double for bin edges 
int bin_number(double target, const std::vector<double>& bin_edges, float resolution = 1e-6){

	// bin indices
	int low = 0; // zero-based
	int high = bin_edges.size() - 1; // zero-based, last bin index
	int med = (low + high) / 2; // medium or medium left

	// initial check on whether target is out of range
	if (target < bin_edges.at(low)) return 0;
	if (target > bin_edges.at(high)) return bin_edges.size();

	do {
		if (abs(target - bin_edges.at(med)) < resolution) { // medium value equals target
			return med;
		}
		else{
			if (target > bin_edges.at(med)){ // target > medium
				low = med;
			} else{ // target < medium
				high = med;
			}
			med = (low + high) / 2;
			if (low == med){ // with med = (low + high) / 2, target always >= medium at this stage
				return (target - bin_edges.at(med) <= bin_edges.at(med+1) - target)? med : med + 1; // return the closest 
			}
		}

		if (low >= med || high <= med){ // should not happen: give ERROR and return nonsensible result
			std::cerr << "Either low >= medium or high <= medium; not already returned desired bin number" << std::endl;
			std::cerr << "low: " << low << ", medium " << med << ", high: " << high << std::endl;
			break;
		}
	}while(1);

	return -100;
}

int bin_number(float target, const std::vector<double>& bin_edges){
	return bin_number((double) target, bin_edges);
}

int bin_number(double target, const std::vector<float>& bin_edges){
	std::vector<double> bin_edges_double(bin_edges.begin(), bin_edges.end());
	return bin_number(target, bin_edges_double);
}

int bin_number(float target, const std::vector<float>& bin_edges){
	std::vector<double> bin_edges_double(bin_edges.begin(), bin_edges.end());
	return bin_number((double) target, bin_edges_double);
}

void test(){
	std::vector<double> v = {0,1,2,3,4,5,6,7,8,9,10};
	cout << bin_number(2.3, v) << endl;
	cout << bin_number(2.8, v) << endl;
	cout << bin_number(7.3, v) << endl;
	cout << bin_number(7.7, v) << endl;
	cout << bin_number(8., v) << endl;
	cout << bin_number(8.2, v) << endl;

	std::vector<double> v2 = { -2.4000000, -2.3000000, -2.2000000, -2.1000000, -2.0000000, -1.9000000, -1.8000000, -1.7000000, -1.6000000, -1.5000000, -1.4000000, -1.3000000, -1.2800000, -1.2600000, -1.2400000, -1.2200000, -1.2000000, -1.1800000, -1.1600000, -1.1400000, -1.1200000, -1.1000000, -1.0800000, -1.0600000, -1.0400000, -1.0200000, -1.0000000, -0.98000000, -0.96000000, -0.94000000, -0.92000000, -0.90000000, -0.88000000, -0.86000000, -0.84000000, -0.82000000, -0.80000000, -0.78000000, -0.76000000, -0.74000000, -0.72000000, -0.70000000, -0.68000000, -0.66000000, -0.64000000, -0.62000000, -0.60000000, -0.58000000, -0.56000000, -0.54000000, -0.52000000, -0.50000000, -0.48000000, -0.46000000, -0.44000000, -0.42000000, -0.40000000, -0.38000000, -0.36000000, -0.34000000, -0.32000000, -0.30000000, -0.28000000, -0.26000000, -0.24000000, -0.22000000, -0.20000000, -0.19000000, -0.18000000, -0.17000000, -0.16000000, -0.15000000, -0.14000000, -0.13000000, -0.12000000, -0.11000000, -0.10000000, -0.090000000, -0.080000000, -0.070000000, -0.060000000, -0.050000000, -0.040000000, -0.030000000, -0.020000000, -0.010000000, 0.0000000, 0.010000000, 0.020000000, 0.030000000, 0.040000000, 0.050000000, 0.060000000, 0.070000000, 0.080000000, 0.090000000, 0.10000000, 0.11000000, 0.12000000, 0.13000000, 0.14000000, 0.15000000, 0.16000000, 0.17000000, 0.18000000, 0.19000000, 0.20000000, 0.22000000, 0.24000000, 0.26000000, 0.28000000, 0.30000000, 0.32000000, 0.34000000, 0.36000000, 0.38000000, 0.40000000, 0.42000000, 0.44000000, 0.46000000, 0.48000000, 0.50000000, 0.52000000, 0.54000000, 0.56000000, 0.58000000, 0.60000000, 0.62000000, 0.64000000, 0.66000000, 0.68000000, 0.70000000, 0.72000000, 0.74000000, 0.76000000, 0.78000000, 0.80000000, 0.82000000, 0.84000000, 0.86000000, 0.88000000, 0.90000000, 0.92000000, 0.94000000, 0.96000000, 0.98000000, 1.0000000, 1.0200000, 1.0400000, 1.0600000, 1.0800000, 1.1000000, 1.1200000, 1.1400000, 1.1600000, 1.1800000, 1.2000000, 1.2200000, 1.2400000, 1.2600000, 1.2800000, 1.3000000, 1.3200000, 1.3400000, 1.3600000, 1.3800000, 1.4000000, 1.5000000, 1.6000000, 1.7000000, 1.8000000, 1.9000000, 2.0000000, 2.1000000, 2.2000000, 2.3000000, 2.4000000 };
	cout << bin_number(-0.5, v2) << ", ";
	cout << find(v2.begin(), v2.end(), -0.5) - v2.begin() << endl;
	cout << bin_number(-2.4, v2) << ", ";
	cout << find(v2.begin(), v2.end(), -2.4) - v2.begin() << endl;
	cout << bin_number(-2.3, v2) << ", ";
	cout << find(v2.begin(), v2.end(), -2.3) - v2.begin() << endl;
}
