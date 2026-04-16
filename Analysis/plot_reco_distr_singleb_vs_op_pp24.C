// plot_reco_distr_singleb_vs_op_pp24.C
// 2×2 canvas: reco pair_pt (top row) + minv (bottom row)
//   left  col: truth single-b pairs passing medium WP  (from_same_b && pair_pass_medium)
//   right col: all op-sign pairs passing medium WP     (pair_pass_medium, reco-only, like data)
// Drawn with error bars; y-axis = estimated cross-section [pb / bin-width].

void plot_reco_distr_singleb_vs_op_pp24(){
    gROOT->SetBatch(kTRUE);

    // Load struct dict (needed for MuonPairObj.xxx sub-branch access in TTree::Draw)
    gROOT->ProcessLine(".L /gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/Analysis/NTupleProcessingCode/PythiaAnalysisClasses.h+");

    const char* fpath  = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/"
                         "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts.root";
    const char* outdir = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/plots/";

    TFile* f = TFile::Open(fpath, "READ");
    if (!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", fpath); return; }
    TTree* t = (TTree*)f->Get("muon_pair_tree_sign2");
    if (!t){ printf("ERROR: muon_pair_tree_sign2 not found\n"); f->Close(); delete f; return; }

    // ---- binnings ----
    // pair_pt: 12 log bins from 8-80 GeV  (matches pT_bins_80)
    const int n_pt = 12;
    double pt_e[n_pt + 1];
    for (int i = 0; i <= n_pt; i++)
        pt_e[i] = 8.0 * pow(80.0 / 8.0, double(i) / n_pt);

    // minv: 50 bins from 1.0-3.0 GeV (data signal window)
    const int n_mv = 50;
    const double mv_lo = 1.0, mv_hi = 3.0;

    // ---- create histograms ----
    TH1D* h_pt_sb = new TH1D("h_pt_sb", "", n_pt, pt_e);  h_pt_sb->Sumw2();
    TH1D* h_pt_op = new TH1D("h_pt_op", "", n_pt, pt_e);  h_pt_op->Sumw2();
    TH1D* h_mv_sb = new TH1D("h_mv_sb", "", n_mv, mv_lo, mv_hi); h_mv_sb->Sumw2();
    TH1D* h_mv_op = new TH1D("h_mv_op", "", n_mv, mv_lo, mv_hi); h_mv_op->Sumw2();

    // ---- fill ----
    // left: truth single-b only (no reco signal cuts)
    const char* sel_sb = "MuonPairObj.weight * MuonPairObj.from_same_b";
    // right: reco signal cuts like data (from RDFBasedHistFillingPP.cxx)
    const char* sel_op = "MuonPairObj.weight * (MuonPairObj.pair_pass_medium"
                         " && MuonPairObj.minv > 1.08 && MuonPairObj.minv < 2.9"
                         " && MuonPairObj.pair_pt > 8"
                         " && MuonPairObj.m1.charge * MuonPairObj.m1.eta < 2.2"
                         " && MuonPairObj.m2.charge * MuonPairObj.m2.eta < 2.2"
                         " && MuonPairObj.dr > 0.05)";

    t->Draw("MuonPairObj.pair_pt >> h_pt_sb", sel_sb, "goff");
    t->Draw("MuonPairObj.pair_pt >> h_pt_op", sel_op, "goff");
    t->Draw("MuonPairObj.minv    >> h_mv_sb", sel_sb, "goff");
    t->Draw("MuonPairObj.minv    >> h_mv_op", sel_op, "goff");

    // ---- divide by bin width -> d sigma / d var ----
    auto divByBinWidth = [](TH1D* h){
        for (int ib = 1; ib <= h->GetNbinsX(); ++ib){
            double w = h->GetBinWidth(ib);
            if (w > 0){
                h->SetBinContent(ib, h->GetBinContent(ib) / w);
                h->SetBinError  (ib, h->GetBinError  (ib) / w);
            }
        }
    };
    divByBinWidth(h_pt_sb); divByBinWidth(h_pt_op);
    divByBinWidth(h_mv_sb); divByBinWidth(h_mv_op);

    // ---- style ----
    auto sty = [](TH1D* h, const char* xtitle){
        h->SetStats(0); h->SetTitle("");
        h->GetXaxis()->SetTitle(xtitle);
        h->GetYaxis()->SetTitle("d#sigma/d(var) [pb/GeV]");
        h->GetYaxis()->SetTitleOffset(1.4);
        h->SetLineColor(kBlack); h->SetMarkerColor(kBlack);
        h->SetLineWidth(2); h->SetMarkerStyle(20); h->SetMarkerSize(0.7);
    };
    sty(h_pt_sb, "p_{T}^{pair,reco} [GeV]");
    sty(h_pt_op, "p_{T}^{pair,reco} [GeV]");
    sty(h_mv_sb, "m_{#mu#mu}^{reco} [GeV]");
    sty(h_mv_op, "m_{#mu#mu}^{reco} [GeV]");

    // ---- canvas: 2 rows × 2 cols ----
    TCanvas c("c_reco_distr", "c_reco_distr", 1300, 900);
    c.Divide(2, 2);

    // column labels (drawn in the top two pads)
    const char* lbl_sb = "single-b (truth) + medium WP";
    const char* lbl_op = "all op-sign, medium WP (reco, like data)";

    // --- pad 1: pair_pt single-b ---
    c.cd(1);
    gPad->SetTicks(1,1); gPad->SetLogx(1); gPad->SetLogy(1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    h_pt_sb->Draw("E1");
    TLatex lab1; lab1.SetNDC(); lab1.SetTextSize(0.052);
    lab1.DrawLatex(0.14, 0.93, lbl_sb);

    // --- pad 2: pair_pt all-op ---
    c.cd(2);
    gPad->SetTicks(1,1); gPad->SetLogx(1); gPad->SetLogy(1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    h_pt_op->Draw("E1");
    TLatex lab2; lab2.SetNDC(); lab2.SetTextSize(0.052);
    lab2.DrawLatex(0.14, 0.93, lbl_op);

    // --- pad 3: minv single-b ---
    c.cd(3);
    gPad->SetTicks(1,1); gPad->SetLogy(1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    h_mv_sb->Draw("E1");

    // --- pad 4: minv all-op ---
    c.cd(4);
    gPad->SetTicks(1,1); gPad->SetLogy(1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    h_mv_op->Draw("E1");

    gSystem->mkdir(outdir, kTRUE);
    c.SaveAs((std::string(outdir) + "reco_distr_singleb_vs_op_mediumWP.png").c_str());
    printf("Saved to %sreco_distr_singleb_vs_op_mediumWP.png\n", outdir);

    // ---- print summary crossx ----
    auto xsec_total = [](TH1D* h){
        double s = 0;
        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
            s += h->GetBinContent(ib) * h->GetBinWidth(ib);
        return s;
    };
    printf("pair_pt xsec (single-b): %.4f pb\n",  xsec_total(h_pt_sb));
    printf("pair_pt xsec (all op):   %.4f pb\n",  xsec_total(h_pt_op));
    printf("minv    xsec (single-b): %.4f pb\n",  xsec_total(h_mv_sb));
    printf("minv    xsec (all op):   %.4f pb\n",  xsec_total(h_mv_op));

    delete h_pt_sb; delete h_pt_op; delete h_mv_sb; delete h_mv_op;
    f->Close(); delete f;
}
