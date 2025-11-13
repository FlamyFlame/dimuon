#include <string>

void hist_helper(TH1* h, float norm, bool norm_unity, std::string title, std::string ytitle=""){

  h->SetStats(0);

    if (norm_unity){
        h->Scale(1.,"width");
        try{
            if (h->Integral("width") == 0) throw std::runtime_error("Histogram interval is ZERO! Cannot normalize to unity!");
            h->Scale(1./h->Integral("width"));
        }catch (const std::runtime_error& e){
            std::cout << "runtime_error caught: " << e.what() << std::endl;
            std::cout << "Proceed without normalizing to unity!" << std::endl; 
        }
        h->GetYaxis()->SetTitle("pdf");
    }else{
        h->Scale(norm,"width");
        if (ytitle.length() != 0){
            h->GetYaxis()->SetTitle(ytitle.c_str());
        }
    }

    // h->SetTitle(title.c_str());
    // h->SetTitleSize(35);
    h->GetYaxis()->SetLabelFont(43);
    h->GetYaxis()->SetLabelSize(28);
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitleSize(28);
    h->GetYaxis()->SetTitleOffset(2.1);
    h->GetXaxis()->SetLabelFont(43);
    h->GetXaxis()->SetLabelSize(28);
    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetXaxis()->SetTitleFont(43);  
    h->GetXaxis()->SetTitleSize(28);
    h->GetXaxis()->SetTitleOffset(1);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.9);
}


std::string CutBeforeBracket(const std::string& s) {
    // find both patterns
    size_t pos1 = s.find(" [");
    size_t pos2 = s.find("[");

    // choose the earliest valid occurrence
    size_t pos = std::string::npos;

    if (pos1 != std::string::npos) pos = pos1;
    if (pos2 != std::string::npos) pos = (pos == std::string::npos ? pos2 : std::min(pos, pos2));

    // if neither found, return original
    if (pos == std::string::npos) return s;

    return s.substr(0, pos);
}

void thstack_helper(THStack* h, std::string kin_name, std::string title){
  h->SetTitle(title.c_str());
  h->GetXaxis()->SetTitle(kin_name.c_str());
  
  // we manually set the y-title titles since unlike histogram drawing
  // where the first histogram drawn sets the properties for the subcanvases (e.g, axis properties)
  // THStack does not "inherit y-axis titles" from the individual histograms in the stack
  std::string ytitle = "";
  const std::string prefix = "FULL: ";
  if (kin_name.rfind(prefix, 0) == 0) {  // starts with "FULL: "
      ytitle = kin_name.substr(prefix.length());
  } else {
      ytitle = "d#sigma/d" + CutBeforeBracket(kin_name) + " [nb]";
  }

  h->GetYaxis()->SetTitle(ytitle.c_str());

  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(23);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitleOffset(2);

  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(23);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(25);
  h->GetXaxis()->SetTitleOffset(1);
}

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
