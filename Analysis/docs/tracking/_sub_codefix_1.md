# Sub-agent scratch: round-9 reviewer fixes C1-C7 (codefix_1)

Append-only. Owner: subagent. Orchestrator merges + deletes.

## Mandate
Fix 5 WARNINGs + 2 INFOs from the round-9 analysis-code review. NO fitted or measured
number may move (no fit function / range / limit / plateau window / selection / binning
change). Editable files ONLY:
- RDFBasedHistFilling/FillMCTrigEffHists.cxx
- plotting_codes/trig_effcy/mc_based/fit_dr_corrections.cxx
- plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff.cxx
- pipelines/run_mc_trigeff_round7.sh + pipelines/run_dr_correction_fits.sh
No git. Samples: pp_full, overlay only (never noovl).

## Step 0 (DONE) -- baseline snapshot
36 fit_report*.txt (pp_trigger_efficiency/mc_based{,_medium} x step{3,4} x
{expo,polyu_fixedRp,interp} x {intgr,same_sign,opposite_sign}) copied to
/tmp/claude-101379/-gpfs-mnt-atlasgpfs01-usatlas-workarea-yuhanguo-dimuon-codes-Analysis/104587a7-5fbd-4e7a-95b1-38298be6d69b/scratchpad/codefix_baseline/

Baseline (BUGGY) "cells marked UNUSABLE" counts, Step 3:
  Tight  (mc_based):        expo intgr 12 / ss 25 / os 15 ; polyu_fixedRp 12 / 31 / 13 ; interp 12 / 11 / 13
  Medium (mc_based_medium): expo intgr 12 / ss 24 / os 13 ; polyu_fixedRp 12 / 30 / 11 ; interp 12 / 11 / 11
Step 4 Tight: expo 2/5/3 ; polyu 2/4/3 ; interp 2/2/3
Step 4 Medium: expo 2/2/2 ; polyu 2/4/2 ; interp 2/2/2

## Step 1 (DONE) -- C1: n_unusable missed both early-continue paths
`fit_dr_corrections.cxx`:
- `++n_unusable` added, guarded by `!inclusive`, in BOTH early exits:
  * the `pnb <= 0 || !DrCorrPlateauUsable(...)` branch ("no plateau -> cell skipped")
  * the `k < M.nfree + 2` branch ("too few points -> no fit")
- declaration comment (was "fit failed / plateau outside tolerance / correction not > 0")
  rewritten to "EVERY cell persisted with fit_ok = 0, by whichever route".
- report line reason-list extended to name the two new routes ("no measurable plateau, too
  few points to fit, ..."). The `cells marked UNUSABLE (h_stepN_fit_ok = 0` prefix is kept;
  verified by grep that NO code anywhere parses this line.

### CORRECTED counts (Step 3, pp_full), old -> new
| WP | series | expo | polyu_fixedRp | interp |
|---|---|---|---|---|
| Tight | sign-integrated | 12 -> **13** | 12 -> **13** | 12 -> **13** |
| Tight | same sign       | 25 -> **32** | 31 -> **38** | 11 -> **15** |
| Tight | opposite sign   | 15 -> **17** | 13 -> **15** | 13 -> **15** |
| Medium | sign-integrated| 12 -> **12** | 12 -> **12** | 12 -> **12** |
| Medium | same sign      | 24 -> **31** | 30 -> **37** | 11 -> **15** |
| Medium | opposite sign  | 13 -> **15** | 11 -> **13** | 11 -> **13** |

Step 4 for completeness, old -> new:
| WP | series | expo | polyu_fixedRp | interp |
|---|---|---|---|---|
| Tight | sign-integrated | 2 -> 2 | 2 -> 2 | 2 -> 2 |
| Tight | same sign       | 5 -> 12 | 4 -> 11 | 2 -> 6 |
| Tight | opposite sign   | 3 -> 3 | 3 -> 3 | 3 -> 3 |
| Medium | sign-integrated| 2 -> 2 | 2 -> 2 | 2 -> 2 |
| Medium | same sign      | 2 -> 9 | 4 -> 11 | 2 -> 6 |
| Medium | opposite sign  | 2 -> 2 | 2 -> 2 | 2 -> 2 |

INDEPENDENT CROSS-CHECK: counted fit_ok == 0 bins directly in the 18 Step-3
`dr_correction_fits_pp24_full*_step3_*{,_ss,_os}.root` `h_step3_fit_ok` TH2Ds -- every one
equals the NEW report number exactly (13/32/17, 13/38/15, 13/15/15 Tight;
12/31/15, 12/37/13, 12/15/13 Medium). The map and its caption now agree.

## Step 2 (DONE) -- C2: charge-blindness rationale
`FillMCTrigEffHists.cxx` ~line 690: comment rewritten. Blindness restricted to STEP 4 (and
R25's f(0) = 1.047 vs 1.148 quoted as the evidence); for STEP 3 it says blindness is REFUTED
by R25 (0.427 at dR = 0.025, surviving in the raw un-inverse-weighted numraw/denom) and that
the sign-integrated Step-3 series stays nominal only PENDING AN OPEN USER DECISION. Comment
only; no executable line touched.

## Step 3 (DONE) -- C3: 2D-singles stage artefact check
`run_mc_trigeff_round7.sh`: the stage no longer trusts the exit code. It extracts the macro's
own two `wrote` lines (`sed -n 's/^ *wrote //p'`), requires exactly 2, `val_file`s each, and
requires the two expected basenames under `step1_singles_data_mc/`. Paths come from the LOG,
not re-derived in shell, because the directory is built in C++ from DrCorrSampleCfg::out_base.
VERIFIED on a real run of plot_mc_singles_2d_effcy("pp_full", true): n=2, both files
non-empty, both basenames matched.

## Step 4 (DONE) -- C4: statistics-table stage
`run_mc_trigeff_round7.sh`: only the benign pre-round-9 case is downgraded to a log line --
detected by grepping the log for `is MISSING from` (the Get2D throw). Everything else
(pair-pT bin-count vs MCTrigEffPairPt::NBins(), axis-edge mismatch among the six TH2Ds,
all-empty binned range) now calls `fail`. VERIFIED: pp_full/tight log has
"CSVs written to" (pass), overlay/tight log has "is MISSING from" (benign note) -- so the
transition case still passes and the guards are live.

## Step 5 (DONE) -- C5: PNG-count check in the fit driver
`run_dr_correction_fits.sh`: main and ratio canvases counted separately
(`! -name "*_ratio.png"` / `-name "*_ratio.png"`); expected NPT+1 main always, NPT+1 ratio in
`sign_sepr`, 0 ratio in `sign_intgr`. VERIFIED on disk after the re-run:
step3/expo/sign_intgr main=9 ratio=0, sign_sepr main=9 ratio=9; driver reported 0 artefact
failures.

## Step 6 (DONE) -- C6: the RAW joint trigger probability figure
`plot_mc_trig_eff.cxx`, new block after the Step-3 per-sign plateau loop. Draws, per sign,
`numraw/denom` = P(both muons fire | dR) with NO eps_MC anywhere, plus the same-sign /
opposite-sign ratio pad. Errors reuse the SHARED SetConditionalRatioErrors (with
A = B = sum_fired w^2 taken from Sumw2, the a_i = w_i specialisation) rather than a second
copy of the formula; the SS/OS ratio uses quadrature because the two pair samples are
disjoint. Canvas: reserved header strip (headline), reserved strip above the frame carrying
the defining equation `P = Sigma w_MC (pairs that fired) / Sigma w_MC (all pairs)` and the
legend ("same sign, P_SS" / "opposite sign, P_OS"); y title
`P(2mu4 | mu pair)` for pp and `P(both mu pass mu4 | mu pair)` for the overlay (the trigger
requirement genuinely differs, FillMCTrigEffHists `trig_cond`); ratio title `P_SS / P_OS`;
PNG only; no code tokens. Graceful skip with a printed note when numraw is absent.

Output (Tight and Medium):
  .../pp_trigger_efficiency/mc_based/step3_dr_correction/step3_raw_joint_trigger_probability_by_sign.png
  .../pp_trigger_efficiency/mc_based_medium/step3_dr_correction/step3_raw_joint_trigger_probability_by_sign.png
Overlay: SKIPPED with the printed note (h_mc_dr_{ss,os}_zoom_numraw not in its Step-3 file --
it has not been refilled since numraw was booked). NO refill run.

Printed numbers reproduce R25's raw table EXACTLY (pp_full Tight):
  dR 0.025 ss 0.2375 os 0.5564 ratio 0.427
  dR 0.075 ss 0.3135 os 0.5489 ratio 0.571
  dR 0.125 ss 0.4714 os 0.5762 ratio 0.818
  dR 0.175 ss 0.5595 os 0.5789 ratio 0.967
  dR 0.225 ss 0.6316 os 0.6147 ratio 1.027
  dR 0.325 ss 0.6420 os 0.6503 ratio 0.987
PNG inspected: series separated, no overlap of legend/equation/data, ratio pad readable.
NEW OBSERVATION visible in the figure (not previously in the doc): the OPPOSITE-sign raw
probability DIPS to ~0.44 around dR = 0.65-0.75 while the same-sign curve stays flat at
~0.62, driving P_SS/P_OS up to ~1.4 there. Plausibly the OS-only resonance veto applied at
the ntuple stage (memory project_os_resonance_veto) removing a kinematic band. NOT
investigated -- flagged for the orchestrator.

## Step 7 (DONE) -- C7 and a sibling instance
`FillMCTrigEffHists.cxx:1296` and `:1373`: the hardcoded `dR in [1,4]` replaced by
`Form("  large-dR average (dR in [%g,%g], weighted): ", MCTrigEffPlateau::kLo,
MCTrigEffPlateau::kHi)` -> prints `[2,3.5]`, the window the sum actually uses.

SAME DEFECT CLASS FOUND AND FIXED in `plot_mc_trig_eff.cxx:1400` and `:2234` (both editable
files): the plateau printf used `%.0f`, so the nominal window [2, 3.5] printed as "[2,4]" and
was indistinguishable from the RETIRED [1,4] systematic window printed right next to it.
Changed to `%g`; now prints "dR in [2,3.5] ... window syst vs dR in [1,4]". stdout label only,
no number touched.

## Step 8 (DONE) -- proof nothing moved
All 36 pp_full fit_report*.txt snapshotted BEFORE the re-run, then diffed after.
TOTAL differing lines that are not the `cells marked UNUSABLE` line: **0** across all 36
files. Every plateau, every fitted parameter, every chi2/ndf, every f(0), every at-limit
record is byte-identical.
Fit-stage re-run: `SAMPLES=pp_full WPS="tight medium" STEPS="3 4"
METHODS="expo polyu_fixedRp interp" SKIP_MEASURE=1 bash pipelines/run_dr_correction_fits.sh`
-> exit 2 (the pre-existing FULL-sample plateau-guard failure, as expected), **0 artefact
failures**.
All four macros recompile clean; `bash -n` passes on both pipeline scripts.
