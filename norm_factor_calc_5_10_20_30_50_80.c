float single_norm_factor_calc(float taa, int count_5percent){
	// Parameters:
	// taa: TAA in unit of inverse microbarn
	// count_5percent: integer multiple of 5% for the percentage of total data that the current centrality bin covers
	// e.g, count_5percent = 2 means the current centrality bin covers 10% of data

	return (1./taa/(count_5percent * 15. / 20.));
}

float norm_factor_calc_5_10_20_30_50_80(){
	std::cout <<  single_norm_factor_calc(26.0,1) << ", ";
	std::cout <<  single_norm_factor_calc(20.4,1) << ", ";
	std::cout <<  single_norm_factor_calc(14.4,2) << ", ";
	std::cout <<  single_norm_factor_calc(8.77,2) << ", ";
	std::cout <<  single_norm_factor_calc((5.09+2.75)/2,4) << ", ";
	std::cout <<  single_norm_factor_calc((1.35 + 0.601 + 0.239)/3,6) << std::endl;
}