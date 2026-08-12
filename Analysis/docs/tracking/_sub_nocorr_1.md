# _sub_nocorr_1 -- Step-3 dR fits WITHOUT the plateau correction (subagent scratch, append-only)

## Mandate (from the orchestrator)
Add a `plateau_mode` variant to the Step-3 dR-correction fit chain:
  * `corr`   = today's behaviour (each cell's eps_dR divided by its own [2,3.5] plateau), NOMINAL,
               must be BIT-FOR-BIT unchanged, only relocated one directory deeper.
  * `nocorr` = fit the RAW, un-normalised eps_dR with a FREE additive baseline C. The [2,3.5]
               plateau NEVER enters the fit; the baseline comes from the dR<1 data itself.
Directory layout (plateau mode is the TOP level):
  step3_dr_fit/plateau_corrected/<method>/{sign_intgr,sign_sepr}/
  step3_dr_fit/no_plateau_correction/<method>/{sign_intgr,sign_sepr}/
  step4_dr_fit/plateau_corrected/<method>/{sign_intgr,sign_sepr}/    (step 4: corr ONLY)
Editable files (exactly four): fit_dr_corrections.cxx, plot_dr_correction_fits.cxx,
dr_correction_sample_cfg.h, pipelines/run_dr_correction_fits.sh.  NEVER run git.
Run only `pp_full`.

## Physics motivation (user)
The inversely-weighted dR distribution shows structure even in the full range, worst in the
gap-enclosing pair-eta bins => a plateau determined over a LARGE dR window ([2,3.5]) may not be
the right baseline for the SMALL-dR region of interest. The nocorr variant lets the fit decide
the baseline from the dR<1 data alone.

## Step 0 -- baseline snapshot (DONE, before any edit)
scratchpad/nocorr_baseline/ : 84 .txt + 324 .png copied from
  plots/pp_trigger_efficiency/{mc_based,mc_based_medium}/step{3,4}_dr_fit/**
scratchpad/nocorr_baseline_rootfiles.txt : `ls -la` of all dr_correction_{fits,plateaus}_* in
  /usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/

## Fit-domain audit (asked by the orchestrator, answered from the source)
fit_dr_corrections.cxx: kFitLo=0.0, kFitHi=1.0; the fit graph is built ONLY from the ZOOM 3D
histograms (h_mc_dr_zoom_vs_pt_eta_*), keeping bins with kFitLo <= centre <= kFitHi and a
non-zero error; TF1 constructed on [0,1]; `g->Fit(f,"QRNS")` -- "R" = use the function range.
The full-range histograms (dR>1, incl. the [2,3.5] plateau window) are read ONLY by the PLOT
stage, drawn with open markers, and never enter the fit.  => the nocorr path inherits exactly
this: C is determined only from data in [0,1]; no [2,3.5] bin enters the fit.

## Plan
1. dr_correction_sample_cfg.h : DrCorrPlateauModeDir/Tag helpers + optional `plateau_mode` arg on
   DrCorrFitFile (default "" -> byte-identical existing names).
2. fit_dr_corrections.cxx : plateau_mode param; mode dir in ONE place; nocorr = no division, no
   plateau screen/guard consequence, free-C formulas, C initialised from the upper part of the
   fit domain, fit_ok = valid && f>0 on [0,Rp]; plateau still REPORTED (labelled not used).
3. plot_dr_correction_fits.cxx : plateau_mode param; raw points in nocorr; y-title differs;
   per-panel `C (fitted plateau)`; interp prints its pinned baseline; readback flat reference
   generalised to f(Rp) in nocorr.
4. pipeline : PLATEAU_MODES env (default "corr nocorr"), step 4 forced to corr, artefact paths.
5. Run pp_full, both WPs, steps 3+4, 3 methods; verify corr diff == empty; delete the old
   step{3,4}_dr_fit/<method>/ level.

## Step 1-4 -- implementation (DONE)
* `dr_correction_sample_cfg.h`: `DrCorrPlateauModeDir` ("plateau_corrected/" | "no_plateau_correction/")
  and `DrCorrPlateauModeTag` ("" | "_nocorr"); `DrCorrFitFile` gained a 6th arg `plateau_mode = ""`
  -> every existing call site keeps its filename byte-for-byte (verified: the nominal fit files
  still come out as dr_correction_fits_pp24_full_step3_expo.root etc.).
* `fit_dr_corrections.cxx`: 7th arg `plateau_mode = "corr"`. ONE place puts the mode in the path
  (`sdir = out_base + stepN_dr_fit/ + mode_dir`), used by the guard report and by `mdir`.
  nocorr: no `r->Scale(1/plateau)`, no plateau screen, guard never fatal, verdict "not applicable",
  formulas gain a free C as the LAST parameter (expo [3], polyu [4], powerlaw [3]), C started at
  the mean of the fitted points with dR >= Rp and limited to [0, max(5*C0, 2)], interp pins its
  flat branch to the last measured knot, fit_ok = valid && f>0 on [0,Rp], plateau still reported.
* `plot_dr_correction_fits.cxx`: 6th arg `plateau_mode = "corr"`; `norm_of` (1.0 in nocorr) and
  `cell_drawable` (always true in nocorr) replace the plateau divisions/screens; y-axis title is
  the raw efficiency in nocorr and ".../ plateau" in corr; equation shows C; per-panel
  "C (fitted plateau) = v +- e" (sepr: "C", the header defines it; interp: "C (plateau)");
  the [2,3.5] plateau line is NOT drawn in nocorr (report only) so two numbers called "plateau"
  never share a canvas; read-back flatness/persistence reference generalised from 1 to the
  cell's own C.
* `pipelines/run_dr_correction_fits.sh`: PLATEAU_MODES env (default "corr nocorr"), pmode loop
  between step and method, step 4 forced to corr, all artefact paths/tokens mode-aware, guard
  fatality only in corr, chi2 summary reads the all-cells line in nocorr.

## Step 5 -- full run
SAMPLES=pp_full WPS="tight medium" STEPS="3 4" METHODS="expo polyu_fixedRp interp" SKIP_MEASURE=1

## Step 5 -- RESULTS

### Runs
1. `PLATEAU_MODES` default (corr+nocorr), SAMPLES=pp_full WPS="tight medium" STEPS="3 4"
   METHODS="expo polyu_fixedRp interp" SKIP_MEASURE=1  -> /tmp/drfit_full_run.log
   0 ARTEFACT FAILURES; exit 2 (the pre-existing FULL-sample plateau-guard failure, corr only:
   pp_full/{tight,medium}/step{3,4}).  nocorr never trips the guard (verdict "not applicable").
2. Guard-report banner made nocorr-only so the corr report stays byte-identical -> re-ran
   PLATEAU_MODES=corr (/tmp/drfit_corr_rerun.log): 0 ARTEFACT FAILURES, exit 2.
3. `C (fitted plateau) = ...` was clipping the frame edge at xcol=0.55 -> single-series column
   moved to 0.45 in nocorr ONLY (corr keeps 0.55) and the nocorr plots were re-made
   (SKIP_FIT=1, /tmp/drfit_nocorr_replot.log): 0 ARTEFACT FAILURES.

### corr did NOT move (requirement C)
scratchpad/check_corr_unchanged.sh compares every file of the pre-change snapshot against
<base>/step<N>_dr_fit/plateau_corrected/<same relative path>:
   identical=408  differing=0  missing=0     (324 PNG + 84 TXT)
i.e. every fit_report*.txt, plateau_guard_report*.txt, readback_check*.txt and every PNG
(sign_intgr AND sign_sepr, both WPs, steps 3 and 4) is BYTE-IDENTICAL; only the path changed.
The nominal fit ROOT files kept their exact names (dr_correction_fits_pp24_full[_medium_wp]_
step<N>_<method>[_ss|_os].root); nocorr adds 18 new files with the `_nocorr` token (step 3 only).

### Inclusive cell, step 3, expo, sign-integrated  (fit range dR in [0,1] in BOTH modes)
Tight   corr  : A=-0.1881+-0.0024  lambda=0.2530+-0.0026  p=2.4812+-0.0833
                chi2/ndf=10.248 (chi2=174.22, ndf=17)
Tight   nocorr: A=-0.1848+-0.0030  lambda=0.2513+-0.0031  p=2.5199+-0.0931  C=0.9906+-0.0015
                chi2/ndf=10.824 (chi2=173.18, ndf=16)     [2,3.5] plateau (not used) = 0.9921
Medium  corr  : A=-0.1853+-0.0023  lambda=0.2594+-0.0025  p=2.5364+-0.0825
                chi2/ndf=11.266 (chi2=191.52, ndf=17)
Medium  nocorr: A=-0.1832+-0.0028  lambda=0.2588+-0.0030  p=2.5508+-0.0908  C=0.9920+-0.0014
                chi2/ndf=11.960 (chi2=191.36, ndf=16)     [2,3.5] plateau (not used) = 0.9925
=> nested-model sanity check passes: chi2 DROPS slightly when C is freed (174.22->173.18,
   191.52->191.36) while chi2/ndf rises only because ndf goes 17->16. Inclusively C agrees with
   the [2,3.5] plateau to 0.15%, i.e. the inclusive cell is NOT where the two differ.

### fit_ok, step 3, sign-integrated, 72 (pair pT, pair eta) cells
Tight  expo          : corr 59, nocorr 68 | both 58, only corr 1, only nocorr 10, neither 3
Tight  polyu_fixedRp : corr 59, nocorr 66 | both 59, only corr 0, only nocorr  7, neither 6
Tight  interp        : corr 59, nocorr 70 | both 59, only corr 0, only nocorr 11, neither 2
Medium expo          : corr 60, nocorr 67 | both 58, only corr 2, only nocorr  9, neither 3
Medium polyu_fixedRp : corr 60, nocorr 66 | both 60, only corr 0, only nocorr  6, neither 6
Medium interp        : corr 60, nocorr 70 | both 60, only corr 0, only nocorr 10, neither 2
"only nocorr" = cells corr rejected on |plateau-1| > 0.15 / unmeasurable plateau, which cannot
disqualify anything when nothing is normalized.  "only corr" = 3 cells where the nocorr fit is
rejected by the POSITIVITY screen (f < 0 somewhere on [0,Rp]):
  Tight  pT[50,72.1) x eta[-1.5,-1.0): min f = -0.142 (A=-1.389, lam=0.088, p=0.458, C=1.247)
  Medium pT[50,72.1) x eta[-1.5,-1.0): min f = -0.338 (A=-1.662, lam=0.082, p=0.366, C=1.324)
  Medium pT[50,72.1) x eta[ 1.5, 2.0): min f = -0.671 (A=-1.813, lam=0.021, p=0.367, C=1.142)

### Deleted (old layout, superseded)
408 files listed in scratchpad/deleted_old_layout.txt:
  {mc_based,mc_based_medium}/step{3,4}_dr_fit/{expo,polyu_fixedRp,interp}/  (27 PNG + 6 TXT each)
  {mc_based,mc_based_medium}/step{3,4}_dr_fit/plateau_guard_report{,_same_sign,_opposite_sign}.txt
NOT touched: mc_based_pt4bin*, the pbpb (overlay) trees and the *.bak trees -- they still hold the
old layout because this task did not regenerate them (pp_full nominal binning only).

### Run 4 (FINAL, authoritative) -- /tmp/drfit_final_run.log
The guard report's "# rule: ... the fatal tier is ENFORCED" sentence contradicted the nocorr mode
banner, so the rule paragraph was made mode-aware (nominal branch byte-for-byte unchanged) and the
FULL matrix (corr + nocorr, both WPs, steps 3+4, 3 methods, 3 signs) was re-run from scratch:
  0 ARTEFACT FAILURES; exit 2 from the pre-existing corr plateau-guard FAIL (pp_full tight+medium,
  steps 3 and 4). nocorr verdict everywhere: "not applicable".
  corr-vs-baseline re-check after this run: identical=408, differing=0, missing=0.
Final tree: 486 PNG + 126 TXT under {mc_based,mc_based_medium}/step{3,4}_dr_fit/<mode>/<method>/,
9 PNG per sign_intgr leaf (8 pair-pT + inclusive), 18 per sign_sepr leaf (+ the ratio family).

### Not done / notes
* Step 4 nocorr deliberately NOT produced (user: ignore step 4). The tree shape is the same:
  step4_dr_fit/plateau_corrected/... only.
* overlay / noovl / mc_based_pt4bin* still carry the OLD (mode-less) layout -- they were not
  re-run, so their outputs were left in place rather than deleted with nothing to replace them.
