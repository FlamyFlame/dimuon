const int nbins = 16;
const float ctr_intvl = 5.;
const float ctr_max = ctr_intvl * nbins;

TH1F *hist_op = new TH1F("hist_op", "My Histogram", nbins, 0, ctr_max);

std::vector<int> nev_ss = {216222, 176461, 138893, 105992, 77987, 56988, 40304, 27617, 18393, 11999, 7655, 4810, 2835, 1729, 973, 571};
std::vector<int> nev_op = {218774, 181482, 146571, 115040, 87534, 65863, 48612, 35240, 24965, 17274, 11861, 8168, 5298, 3421, 2127, 1412};

for (int i = 0; i < nbins; i++){
	hist_op->SetBinContent(i+1, nev_op[i]);
}

hist_op->GetXaxis()->SetTitle("ctr interval (5%)");
hist_op->GetYaxis()->SetTitle("# dimuon events");

// Create TGraph object with scatter plot of bin contents
TGraph *graph_op = new TGraph(hist_op);
graph_op->SetMarkerStyle(20);
graph_op->SetMarkerSize(0.8);
graph_op->Draw("AP");

// Perform linear fit and draw it on the same canvas
TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0, ctr_max);
hist_op->Fit(fitFunc, "Q");
fitFunc->SetLineColor(kRed);
fitFunc->Draw("same");

// Draw legend with fit equation and chi-square value
double a = fitFunc->GetParameter(0);
double b = fitFunc->GetParameter(1);
double c = fitFunc->GetParameter(2);
double d = fitFunc->GetParameter(3);
double chi2 = fitFunc->GetChisquare();
TLegend *legend = new TLegend(0.15, 0.7, 0.5, 0.85);
legend->SetFillColor(0);
legend->SetFillStyle(0);
legend->SetTextSize(0.04);
legend->AddEntry(graph_op, "Data", "p");
legend->AddEntry(fitFunc, Form("Fit: y = %.2f + %.2f x + %.2f x * x + %.2f x * x * x, #chi^{2}/NDF = %.2f", a, b, c, d, chi2/6), "l");
legend->Draw();

// Update canvas
gPad->Modified();
gPad->Update();

