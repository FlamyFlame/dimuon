# Sub-task: MC single-muon mu4 trigger efficiency, 2D map in (q*eta, pT)

Scratch doc for subagent. Append-only. Source of truth for this sub-task.

## Objective
New macro `plotting_codes/trig_effcy/mc_based/plot_mc_singles_2d_effcy.cxx` producing, per WP:
- `step1_eff_2d_pt_vs_q_eta_charge_sepr.png` (2 pads: mu+ left, mu- right)
- `step1_eff_2d_pt_vs_q_eta_charge_comb.png` (1 pad: (num+ + num-)/(den+ + den-))
under `<out_base>/step1_singles_data_mc/`.

## Ownership
Only file I may write: the macro above + this scratch doc. NEVER run git.

## Step 0 — locate the data-side reference (DONE)

Data plot: `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mu4/no_corr/pp24/single_muon_effcy/pt2nd_vs_q_eta2nd_2D_trig_effcy_sepr.png`

Produced by:
- `Analysis/plotting_codes/trig_effcy/TrigEffPlotterPP.cxx:494` `TrigEffPlotterPP::plot2D_SingleMuonEffcy(const std::string& var)`
  - canvas: `TCanvas(...,"",1100,450*nRows)` with nRows=1 -> **1100x450**, `c->Divide(2,1)`  (line 511-512)
  - per pad: `gPad->SetRightMargin(0.15)` (line 530); `gPad->SetLogx(xy.first); gPad->SetLogy(xy.second)` (line 531)
  - title: `Form("%s, %s", "mu4", sign=="_sign1"?"#mu^{+}":"#mu^{-}")` (lines 532-536) -> **sign1 = mu+ (left pad), sign2 = mu- (right pad)**
  - draw: `numH->Draw("COLZ")` (line 537)
  - eff via `numH->Clone(); numH->Divide(denH)` (lines 527-528) -> TH1::Divide leaves bins with 0 denominator at content 0 AND error 0.
  - output: `c->SaveAs(Form("%s/%s_2D_trig_effcy%s.png", outdir, var, fs))` with `outdir = plot_base_dir + "/single_muon_effcy"` (lines 540-543)
  - histogram keys: `hName(var, fs, trg, sign)` = `h_pt2nd_vs_q_eta2nd` + fs(`_sepr`) ... num trg `2mu4`, denom trg `mu4`, sign `_sign1`/`_sign2` (TrigEffPlotterBaseClass.h:498)
- log-axis flags from `trig_effcy_plot_pp.cxx:20`: `{"pt2nd_vs_q_eta2nd",{false,true}}` -> **logx=false (q*eta), logy=true (pT)**.
- z-range: NOT set explicitly in the data code; the palette runs 0->1 because the efficiency ratio populates that range. For the MC plot I set `SetMinimum(0); SetMaximum(1)` explicitly (equivalent look, deterministic).
- Empty cells: the data plot shows white (unfilled) cells at high pT. `TH2::Divide` writes 0 in a 0-denominator cell, which COLZ paints as the bottom palette colour, not white. In the data plot those cells are white because ROOT's COLZ skips cells with content == 0 exactly. To be safe and explicit I will set such cells to `SetBinContent(...,0)` *and* rely on the default `gStyle->SetHistMinimumZero(0)`; ROOT does not draw a COLZ box when content==0 and minimum<=0.  -> verified visually on the produced PNG.

## Step 1 — binning provenance / comparison (DONE)
- MC booking: `Analysis/RDFBasedHistFilling/FillMCTrigEffHists.cxx:958-964`
  `Histo2D({..., ";q#eta;p_{T} [GeV]", bins.q_eta, bins.pt}, "q_eta", "pt", wcol)`
- `Binnings MakeBinnings()` (FillMCTrigEffHists.cxx:229-241):
  - `b.pt = pms.pT_bins_8;  b.pt.insert(end, pms.pT_bins_60...)` -- an EXACT copy of the data
    construction in `RDFBasedHistFilling/RDFBasedHistFillingData.cxx:315-317`
    (`pT_bins_single_muon` = pT_bins_8 ++ pT_bins_60), registered as `hist_binning_map["pT_bins_single_muon"]`.
  - `b.q_eta = ParamsSet::makeEtaTrigEffcyBinning(1)` -- same call as data `RDFBasedHistFillingData.cxx:322`
    (`hist_binning_map["eta_bins_trig_effcy"]`).
  - `ParamsSet.h:544-545`: `fillLogBinningArray(pT_bins_8, 20, 4.0, 8.0)` and
    `fillLogBinningArray(pT_bins_60, 20, 8.0, 60.0)` -> 21 + 21 = 42 edges, 41 bins.

### DUPLICATE 8.0 GeV EDGE — finding
The zero-width bin at pT = 8.0 GeV is **REAL, INTENTIONAL and INHERITED FROM DATA**, not an MC bug.
`pT_bins_8` ends at 8.0 and `pT_bins_60` starts at 8.0; concatenating them without dropping the
shared edge yields `... 7.72749, 8, 8, 8.84796 ...` i.e. one bin of zero width. The MC code does this
DELIBERATELY and says so at FillMCTrigEffHists.cxx:234-237 ("exact data pt2nd construction ...
this duplicates the 8.0 edge -> one zero-width bin, as in data"), so the MC and data 2D maps are
bin-for-bin identical. The same zero-width bin therefore exists in every data single-muon
trigger-efficiency histogram built from `pT_bins_single_muon`.
Consequence: that bin can never be filled (no value satisfies 8 <= pt < 8), so num = denom = 0 there
and it is simply an always-empty cell. It is invisible on a plot (zero width). NOT fixed here
(FillMCTrigEffHists.cxx is not mine); handled by leaving it empty like any other empty cell.

## Step 2 — binning comparison MC vs DATA (DONE, numerical)
Compared `h_mc_pt_vs_q_eta_denom_muplus` (MC, `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/mc_trig_eff_hists_pp24_full.root`)
against `h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr` (data,
`/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_single_mu4_coarse_q_eta_bin_qeta_fid.root`):

    MC   : nx=184  ny=41
    DATA : nx=184  ny=41
    max |edge difference| :  x = 0    y = 0        -> BIT-IDENTICAL binning

X (q*eta): [-2.4, 2.4], 184 variable bins (0.1 coarse steps outside the fine regions).
Y (pT), all 42 edges:
  4 4.14106 4.28709 4.43828 4.59479 4.75683 4.92458 5.09824 5.27803 5.46416 5.65685 5.85634
  6.06287 6.27667 6.49802 6.72717 6.9644 7.21 7.46426 7.72749 **8 8** 8.84796 9.78579 10.823
  11.9702 13.239 14.6423 16.1943 17.9108 19.8092 21.9089 24.2311 26.7995 29.6401 32.7818
  36.2565 40.0995 44.3498 49.0507 54.2498 60
Zero-width Y bin: index 21, [8, 8]  (the duplicate 8.0 edge -- see the finding above).

## Step 3 — macro written and run (DONE)
File (the ONLY source file I created):
`/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/plotting_codes/trig_effcy/mc_based/plot_mc_singles_2d_effcy.cxx`
Signature: `void plot_mc_singles_2d_effcy(const std::string& sample = "pp_full", bool use_tight_wp = true)`.
Sample identity entirely from `GetDrCorrSample` (`dr_correction_sample_cfg.h`) -> works unchanged for
`"overlay"`; no path, no headline and no binning is hardcoded.

Design points:
- eff = num/den per cell, computed by hand (not `TH2::Divide`) so a cell with an EMPTY DENOMINATOR
  can be set to a sentinel (-1) instead of 0. With `SetMinimum(0)` COLZ skips it -> the cell is left
  white/unfilled, matching the data map. A cell with a populated denominator and empty numerator
  keeps a genuine 0 and IS painted (bottom of the palette).
- error: binomial sqrt(eff(1-eff)/D) (display only; nothing downstream reads it).
- combined canvas: sum(num+) + sum(num-) over sum(den+) + sum(den-), then divide -- NOT the average
  of the two efficiencies.
- axes/ranges/titles taken from the input histograms; log-y (log-binned pT), `SetMoreLogLabels()` +
  `SetNoExponent()` because the default log labelling printed only "10" over 4-60 GeV.
- style: `SetAtlasStyle()` (`Analysis/AtlasStyle.C`), OptStat/OptTitle off, 100 contours, default
  (kBird) palette as in the data map. Right margin 0.17 for the palette + its title.
- on-canvas text is ONLY: axis titles with units (`q #upoint #eta`, `p_{T} [GeV]`, `#varepsilon(mu4)`),
  the charge (`#mu^{+}` / `#mu^{-}` / `#mu^{+} and #mu^{-}`) and the sample/WP headline taken from
  `cfg.sample_text` + "Tight/Medium muons". No step index, no file/histogram names, no prose.
  Headline drawn in its own reserved strip (pads stop at y = 0.93) so it cannot land inside a pad.

## Step 4 — outputs + visual check (DONE)
Ran (from `plotting_codes/trig_effcy/mc_based/`):
  root -l -b -q 'plot_mc_singles_2d_effcy.cxx+("pp_full", true)'
  root -l -b -q 'plot_mc_singles_2d_effcy.cxx+("pp_full", false)'
NOT run: `noovl` (out of scope) and `overlay` (PbPb not being rerun this round; its inputs
`mc_trig_eff_hists_hijing_overlay_pbpb23{,_medium_wp}.root` were verified to exist).

Produced (all verified present and non-trivial):
  /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based/step1_singles_data_mc/step1_eff_2d_pt_vs_q_eta_charge_sepr.png   (50.7 kB)
  /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based/step1_singles_data_mc/step1_eff_2d_pt_vs_q_eta_charge_comb.png   (32.1 kB)
  /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based_medium/step1_singles_data_mc/step1_eff_2d_pt_vs_q_eta_charge_sepr.png (50.8 kB)
  /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based_medium/step1_singles_data_mc/step1_eff_2d_pt_vs_q_eta_charge_comb.png (32.3 kB)

Console numbers (plateau = pT > 8 GeV, integrated over all q*eta):
  Tight : mu+ 0.9138 (D=12.997), mu- 0.9155 (D=13.071), comb 0.9147 (D=26.067)
  Medium: mu+ 0.9115 (D=13.171), mu- 0.9135 (D=13.245), comb 0.9125 (D=26.416)
  Cells: 6400 drawn / 1144 undefined, IDENTICAL for both charges and for the combination.
  That identity is structural, not a bug: 24 q*eta columns are empty for EVERY pT because of the
  single-muon fiducial gap cut, plus the one zero-width pT row ->
  40*24 + 184 = 1144 exactly. Measured empty q*eta runs (from the denominator):
      [-1.200, -1.060)   [-0.060, 0.060)   [2.300, 2.400)
  which is exactly `ParamsSet::single_mu_fiducial_gap_cuts` (forward edge now 2.30, commit a6d1d6e).

Visual check (read each PNG back): headline centred in its own strip, no clipping; y axis now
labels 5,6,7,8,10,20,30,40,50,60; palette 0-1 with its title `#varepsilon(mu4)` fully inside the
canvas; empty gap columns and unpopulated high-pT corners render WHITE, not as efficiency 0; low
efficiency at pT < 6 GeV (turn-on) and the |eta| ~ 1.05 / feet structures visible as expected;
no pad empty, no text collision. Iterated once (added the extra log labels) -- clean afterwards.

## Not done / caveats
- Did NOT touch `FillMCTrigEffHists.cxx` (the duplicate 8.0 GeV edge lives there and is intentional).
- Did NOT run the `overlay` or `noovl` samples.
- Did NOT run git (per hard rule); the new macro is uncommitted.
