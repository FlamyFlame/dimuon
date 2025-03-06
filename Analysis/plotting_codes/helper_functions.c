void ApplyCutsTo1DHistogram(TH1* h, const std::vector<std::pair<double, double>>& cuts) {
    if (!h) return;

    int nBins = h->GetNbinsX();

    for (int bin = 1; bin <= nBins; ++bin) { // Bins indexed from 1 in ROOT
        double binLowEdge  = h->GetBinLowEdge(bin);
        double binHighEdge = h->GetBinLowEdge(bin + 1); // Correct for non-uniform bins

        for (const auto& cut : cuts) {
            double cutLow = cut.first;
            double cutHigh = cut.second;

            // Check if bin overlaps with the cut range
            if (!(binHighEdge <= cutLow || binLowEdge >= cutHigh)) { 
                h->SetBinContent(bin, 0); // the bins will still be drawn, but it will screw integral calculations
                h->SetBinError(bin, 0);
                break; // No need to check further cuts
            }
        }
    }
}

void ApplyCutsTo1DHistogram(TH1* h, const std::vector<std::array<float,2>>& cuts) {
    if (!h) return;

    int nBins = h->GetNbinsX();

    for (int bin = 1; bin <= nBins; ++bin) { // Bins indexed from 1 in ROOT
        double binLowEdge  = h->GetBinLowEdge(bin);
        double binHighEdge = h->GetBinLowEdge(bin + 1); // Correct for non-uniform bins

        for (const auto& cut : cuts) {
            double cutLow = cut.at(0);
            double cutHigh = cut.at(1);

            // Check if bin overlaps with the cut range
            if (!(binHighEdge <= cutLow || binLowEdge >= cutHigh)) { 
                h->SetBinContent(bin, 0); // the bins will still be drawn, but it will screw integral calculations
                h->SetBinError(bin, 0);
                break; // No need to check further cuts
            }
        }
    }
}
