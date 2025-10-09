// helper function that for a given target (edge of region, e.g, for projection)
// and a vector of double for bin edges 
float bin_number(double target, const std::vector<double>& bin_edges, float resolution = 1e-6){

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

float bin_number(float target, const std::vector<double>& bin_edges){
	return bin_number((double) target, bin_edges);
}

float bin_number(double target, const std::vector<float>& bin_edges){
	std::vector<double> bin_edges_double(bin_edges.begin(), bin_edges.end());
	return bin_number(target, bin_edges_double);
}

float bin_number(float target, const std::vector<float>& bin_edges){
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
}
