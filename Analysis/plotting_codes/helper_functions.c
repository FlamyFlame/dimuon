void ApplyCutsTo1DHistogram(TH1* h, const std::vector<std::pair<double, double>>& cuts) {
    if (!h){
        std::cout << "WARNING: Attempting to apply cut to 1D histogram, but histogram pointer is nullptr!" << std::endl;
        std::cout << "Return without applying cuts" << std::endl;
        return;
    }

    if (h->GetDimension() != 1){
        std::cout << "WARNING: Attempting to apply cut to 1D histogram, but histogram dimension is NOT 1!" << std::endl;
        std::cout << "Histogram dimension: " << h->GetDimension() << std::endl;
        std::cout << "Return without applying cuts" << std::endl;
        return;
    }

    int nBins = h->GetNbinsX();

    for (int bin = 1; bin <= nBins; ++bin) { // Bins indexed from 1 in ROOT
        double binLowEdge  = h->GetBinLowEdge(bin);
        double binHighEdge = h->GetBinLowEdge(bin + 1);

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
    if (!h){
        std::cout << "WARNING: Attempting to apply cut to 1D histogram, but histogram pointer is nullptr!" << std::endl;
        std::cout << "Return without applying cuts" << std::endl;
        return;
    }

    if (h->GetDimension() != 1){
        std::cout << "WARNING: Attempting to apply cut to 1D histogram, but histogram dimension is NOT 1!" << std::endl;
        std::cout << "Histogram dimension: " << h->GetDimension() << std::endl;
        std::cout << "Return without applying cuts" << std::endl;
        return;
    }

    int nBins = h->GetNbinsX();

    for (int bin = 1; bin <= nBins; ++bin) { // Bins indexed from 1 in ROOT
        double binLowEdge  = h->GetBinLowEdge(bin);
        double binHighEdge = h->GetBinLowEdge(bin + 1);

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

void ApplyXAxisCutsTo2DHistogram(TH2* h, const std::vector<std::array<float,2>>& cuts) {
    if (!h){
        std::cout << "WARNING: Attempting to apply cut to 2D histogram, but histogram pointer is nullptr!" << std::endl;
        std::cout << "Return without applying cuts" << std::endl;
        return;
    }

    if (h->GetDimension() != 2){
        std::cout << "WARNING: Attempting to apply cut to 2D histogram, but histogram dimension is NOT 2!" << std::endl;
        std::cout << "Histogram dimension: " << h->GetDimension() << std::endl;
        std::cout << "Return without applying cuts" << std::endl;
        return;
    }

    int nBinsX = h->GetNbinsX();
    int nBinsY = h->GetNbinsY();

    std::string hx_name = "hx"+std::to_string(rand()); // randomized hx name to avoid histogram replacement if helper function is applied multiple times in the same macro
    TH1D* hx = h->ProjectionX(hx_name.c_str());
    for (int binx = 1; binx <= nBinsX; ++binx) {
        double binLowEdge  = hx->GetBinLowEdge(binx);
        double binHighEdge = hx->GetBinLowEdge(binx + 1);

        for (const auto& cut : cuts) {
            double cutLow = cut.at(0);
            double cutHigh = cut.at(1);

            // Check if binx overlaps with the cut range
            if (!(binHighEdge <= cutLow || binLowEdge >= cutHigh)) { 
                for (int biny = 1; biny <= nBinsY; ++biny) {
                    h->SetBinContent(binx, biny, 0); // the bins will still be drawn, but it will screw integral calculations
                    h->SetBinError(binx, biny, 0);
                    break; // No need to check further cuts
                }
            }
        }
    }
}

void ApplyYAxisCutsTo2DHistogram(TH2* h, const std::vector<std::array<float,2>>& cuts) {
    if (!h){
        std::cout << "WARNING: Attempting to apply cut to 2D histogram, but histogram pointer is nullptr!" << std::endl;
        std::cout << "Return without applying cuts" << std::endl;
        return;
    }

    if (h->GetDimension() != 2){
        std::cout << "WARNING: Attempting to apply cut to 2D histogram, but histogram dimension is NOT 2!" << std::endl;
        std::cout << "Histogram dimension: " << h->GetDimension() << std::endl;
        std::cout << "Return without applying cuts" << std::endl;
        return;
    }

    int nBinsX = h->GetNbinsX();
    int nBinsY = h->GetNbinsY();

    std::string hy_name = "hy"+std::to_string(rand()); // randomized hy name to avoid histogram replacement if helper function is applied multiple times in the same macro
    TH1D* hy = h->ProjectionY(hy_name.c_str());
    for (int biny = 1; biny <= nBinsY; ++biny) {
        double binLowEdge  = hy->GetBinLowEdge(biny);
        double binHighEdge = hy->GetBinLowEdge(biny + 1);

        for (const auto& cut : cuts) {
            double cutLow = cut.at(0);
            double cutHigh = cut.at(1);

            // Check if biny overlaps with the cut range
            if (!(binHighEdge <= cutLow || binLowEdge >= cutHigh)) { 
                for (int binx = 1; binx <= nBinsX; ++binx) {
                    h->SetBinContent(binx, biny, 0); // the bins will still be drawn, but it will screw integral calculations
                    h->SetBinError(binx, biny, 0);
                    break; // No need to check further cuts
                }
            }
        }
    }
}

