// double* poly_3rd_order_fit (std::vector<int>* nev){
void poly_3rd_order_fit (){
	int nbins = 16;
	// std::vector<int> nev_ss = {216222, 176461, 138893, 105992, 77987, 56988, 40304, 27617, 18393, 11999, 7655, 4810, 2835, 1729, 973, 571};
	std::vector<int> nev = {216222, 176461, 138893, 105992, 77987, 56988, 40304, 27617, 18393, 11999, 7655, 4810, 2835, 1729, 973, 571};
	// std::vector<int> nev_op = {218774, 181482, 146571, 115040, 87534, 65863, 48612, 35240, 24965, 17274, 11861, 8168, 5298, 3421, 2127, 1412};
	// std::vector<int> nev = {218774, 181482, 146571, 115040, 87534, 65863, 48612, 35240, 24965, 17274, 11861, 8168, 5298, 3421, 2127, 1412};

	float ctr_intvl = 5.;
	float ctr_max = ctr_intvl * nbins;

	TH1F *hist = new TH1F("hist", "My Histogram", nbins, 0, ctr_max);

	for (int i = 0; i < nbins; i++){
		hist->SetBinContent(i+1, nev[i]);
	}

	// hist->GetXaxis()->SetTitle("ctr interval (5%)");
	// hist->GetYaxis()->SetTitle("# dimuon events");

	// // Create TGraph object with scatter plot of bin contents
	// TGraph *graph = new TGraph(hist);
	// graph->SetMarkerStyle(20);
	// graph->SetMarkerSize(0.8);
	// graph->Draw("AP");

	// // Perform linear fit and draw it on the same canvas
	TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0, ctr_max);
	hist->Fit(fitFunc, "Q");
	// fitFunc->SetLineColor(kRed);
	// fitFunc->Draw("same");

	// Draw legend with fit equation and chi-square value
	std::cout << fitFunc->GetParameter(0) << std::endl;
	std::cout << fitFunc->GetParameter(1) << std::endl;
	std::cout << fitFunc->GetParameter(2) << std::endl;
	std::cout << fitFunc->GetParameter(3) << std::endl;

	// return fitFunc->GetParameters();

	// double chi2 = fitFunc->GetChisquare();
	// TLegend *legend = new TLegend(0.15, 0.7, 0.5, 0.85);
	// legend->SetFillColor(0);
	// legend->SetFillStyle(0);
	// legend->SetTextSize(0.04);
	// legend->AddEntry(graph, "Data", "p");
	// legend->AddEntry(fitFunc, Form("Fit: y = %.2f + %.2f x + %.2f x * x + %.2f x * x * x, #chi^{2}/NDF = %.2f", a, b, c, d, chi2/6), "l");
	// legend->Draw();

	// // Update canvas
	// gPad->Modified();
	// gPad->Update();

}




