# Same-sign fraction vs pair pT in data (pp24, Pb+Pb 50–80 %)

## Objective
Decide whether the difficulty of producing high-pT SAME-SIGN MC statistics
(`additional_mc_request_jo_test.md`) matters for the analysis. Measure, in DATA, the raw SS and OS
pair counts in the single-b signal region vs pair pT, and the SS fraction SS/(SS+OS):
pp24 and Pb+Pb 23+24+25+26 combined, 50–80 % centrality. If SS is negligible at high pT, a coarse
SS trigger-efficiency correction with its uncertainty in the systematics is acceptable (user).

Deliverable: one 2×2 PNG — columns pp24 | Pb+Pb 50–80 %; top raw counts (SS, OS), bottom SS/(SS+OS).

## Physics Procedure
1. **Motivation.** SS pairs enter only as background (combinatoric, fake, hadronic, and the SS part
   of g→bb̄) in the OS/SS template fit. The SS fraction in data bounds how much any SS-specific
   correction (trigger efficiency) can move the signal.
2. **Top-level quantity.** Per pair-pT bin: `N_SS`, `N_OS` = UNWEIGHTED pair counts passing the
   exact crossx signal region; `f_SS = N_SS / (N_SS + N_OS)` (user 2026-10-07), binomial error.
3. **Method.** (a) Selection = the nominal crossx `signal_cuts` (mass window from
   `ParamsSet::SignalMinvCutExpr`, pair-pT cut, gap cuts on both muons, pair |η|<2.2) + Tight WP,
   identical for SS and OS — the same RDF nodes that fill the existing SS/OS crossx yields. Mass
   range = SIGNAL WINDOW (user 2026-10-07; the template fit now also uses the signal window).
   The OS-only resonance veto has no hole inside the window, so SS/OS are on equal footing.
   (b) Counts are raw (no trigger/reco/lumi/T_AA weight). (c) Pb+Pb: centrality bin 50–80 from
   the nominal `ctr_bin_edges`, years summed as raw counts. (d) Pair-pT axis = the default crossx
   axis `ParamsSet::pT_bins_150` (same as the existing count histograms), projected over pair η.
4. **Negative constraints.** No new binning; no efficiency correction; no per-year Pb+Pb panels;
   SS and OS must come from the same selection code path.

## Implementation Plan
1. Book SS count histograms beside the existing OS ones (PP: `h2d_ss_counts_pair_pt_pair_eta_binned_w_signal_cuts`;
   PbPb: `h2d_ss_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_<ctr>_counts`). Per §3a,b. `/review-analysis-code`.
2. Back up nominal crossx outputs; rerun crossx stages pp24, pbpb23/24/25/26; verify pre-existing
   histograms unchanged (old vs new OS counts).
3. Plot macro → 2×2 PNG. Per §3c,d. `/review-plot`.
4. Answer: SS % at high pT.

## Progress Log
- **2026-10-07, step 1.** Booked `h2d_ss_counts_pair_pt_pair_eta_binned_w_signal_cuts` (RDFBasedHistFillingPP.cxx,
  on `df_ss_weighted`, unweighted like the OS twin) and `h2d_ss_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_<ctr>_counts`
  (RDFBasedHistFillingPbPb.cxx, `df_ss`.Filter(signal_cuts).Filter(ctr_filter), weight "weight" like the OS twin).
  Nominal outputs backed up to the session scratchpad. Step 2 rerun launched (pp24, pbpb23..26 sequential).
  Plot macro `plotting_codes/single_b_analysis/plot_ss_fraction_vs_pair_pt.cxx` written (asserts integral == entries,
  identical axes across pp/years; CP 68.3 % errors on f_SS).
- **2026-10-07, step 2 DONE.** Reruns pp24 + pbpb23/24/25/26 exit 0 (2.5 min pp, ~0.5–2.5 min each PbPb).
  Old (backup) vs new: every pre-existing histogram present; max relative bin difference ≤ 1.2e-14 (summation
  order of the multithreaded fill), 0 changed entries, 0 count-histogram differences; pbpb24 bit-identical.
  Added keys: 1 (pp), 6 per PbPb year (one per centrality bin). `/review-analysis-code` PASS (whitespace NOTE fixed).
- **2026-10-07, step 3.** Plot + CSV produced: `dimuon_data/plots/single_b_analysis/ss_fraction/
  ss_os_counts_and_ss_fraction_vs_pair_pt.{png,csv}`. `/review-plot` PASS; its W1/W2 applied (bottom row fixed −0.5…20 %, few-pair CP bars clipped, values
  in the CSV; N_SS=0 points lifted off the frame line); label "50–80 % centrality". Plot regenerated.

## Results & Observations
Signal window, Tight, raw counts. Pb+Pb 50–80 % per year OS/SS: 23 2584/213, 24 1817/161, 25 5229/456, 26 5001/456.
f_SS = SS/(SS+OS) in the canonical coarse pair-pT bins (ParamsSet coarse = every 2nd fine edge):

| pair pT [GeV] | pp OS | pp SS | pp f_SS | PbPb 50–80 OS | SS | f_SS |
|---|---|---|---|---|---|---|
| 9–12.8 | 196427 | 18043 | 8.4 % | 5681 | 702 | 11.0 % |
| 12.8–18.2 | 203182 | 14191 | 6.5 % | 5459 | 426 | 7.2 % |
| 18.2–25.9 | 85110 | 4537 | 5.1 % | 2408 | 124 | 4.9 % |
| 25.9–36.7 | 25101 | 1044 | 4.0 % | 764 | 26 | 3.3 % |
| 36.7–52.2 | 5625 | 199 | 3.4 % | 251 | 8 | 3.1 % |
| 52.2–74.2 | 1041 | 18 | 1.7 % | 59 | 0 | 0 |
| 74.2–105.5 | 144 | 2 | 1.4 % | 9 | 0 | 0 |
| 105.5–150 | 16 | 2 | (11 %, 2 pairs) | 0 | 0 | – |

pp above 52 GeV: 22 SS / 1223 pairs = 1.8 %. Pb+Pb 50–80 % above 52 GeV: 0 SS / 68.
f_SS falls monotonically with pair pT (≈9 % → ≈1.5–2 % in pp); pp and peripheral Pb+Pb agree within stats
from 13 GeV up. The data SS includes combinatoric + fake/hadronic + SS g→bb̄, so it is an UPPER bound on the
SS g→bb̄ component.

### Conclusion
Within the signal window the SS fraction in data falls from ≈8–11 % (9–13 GeV) to ≈3 % at 37–52 GeV and
≈1.5–2 % above 52 GeV in pp (22 of 1223 pairs). In Pb+Pb 50–80 % there are 0 SS pairs among 68 above 52 GeV
(Clopper–Pearson 68.3 % interval upper edge ≈2.7 %). The data SS is an upper bound on SS g→bb̄. So a ±X % SS trigger-efficiency
uncertainty moves the OS−SS / template result by ≲0.02·X % above 52 GeV (pp). Statistical caveat: >100 GeV
holds only 2 (pp) / 0 (Pb+Pb) SS pairs — the fraction there is unmeasured, not shown to be negligible.

## Latest Stage
DONE 2026-10-07. Open for the user: adopt the coarse SS trigger-efficiency correction + systematic.
