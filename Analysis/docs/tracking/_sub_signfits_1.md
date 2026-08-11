# _sub_signfits_1 -- sign-separated dR-correction fits + fit plots (subagent scratch doc)

Owner: subagent (signfits_1). Append-only. Orchestrator merges this into the canonical doc.

## Mandate (from the orchestrator)
- Add a trailing `sign` parameter ("" / "ss" / "os") to `fit_dr_corrections(...)` and
  `fit_dr_corrections_all(...)`; sign-tagged input keys, plateau tag, fit-file name, report names.
- Guard: sign-integrated keeps the FULL-sample THROW; sign-separated measures/reports but NEVER throws.
- New plot layout: `step{3,4}_dr_fit/<method>/{sign_intgr,sign_sepr}/`.
- `plot_dr_correction_fits.cxx` gets a mode parameter; `sign_sepr` overlays both signs + both fits.
- Extend `pipelines/run_dr_correction_fits.sh` with `SIGNS` (default "intgr ss os"); derived PNG count.
- Fix the pre-existing missing `MCTrigEffPairPt::FileSuffix()` in the plot macro's `hist_path`.
- Files I may edit: fit_dr_corrections.cxx, plot_dr_correction_fits.cxx, dr_correction_sample_cfg.h,
  pipelines/run_dr_correction_fits.sh. NEVER git. Outputs only under step{3,4}_dr_fit/ and
  dr_correction_fits_*.root.

## Step 0 -- inputs verified (DONE)
- `mc_trig_eff_hists_pp24_full[_medium_wp]_step3.root` has `h_mc_dr_{ss,os}_{zoom,full}[_vs_pt_eta]_{num,denom,errA,errB}`.
- `..._step4.root` has `h_mc_single_dr_{ss,os}_...` incl. covP/covQ.
- `dr_correction_plateaus_pp24_full.root` keys: `h_step{3,4}[_ss|_os]_plateau[_err|_rms|_nbins|_syst|_inclusive|_syst_inclusive]`, `prov_step{3,4}[_ss|_os]`.
- Inclusive plateaus (Tight, pp_full) read back from the file:
  step3 = 0.992076 +- 0.000873, step3_ss = 0.988930 +- 0.001589, step3_os = 0.993315 +- 0.001044;
  step4 = 0.996215 +- 0.000400, step4_ss = 0.994984 +- 0.000718, step4_os = 0.996703 +- 0.000480.
  -> matches the orchestrator's regression anchors.
- prov stamp format identical for all tags: `sample=pp_full; WP=Tight muons; quantity=eps_dR; ...`
  => the existing sample/WP provenance check works verbatim on `prov_step3_ss`.
- Existing plot dirs: `step{3,4}_dr_fit/{expo,interp,polyu_fixedRp}/` with 9 PNGs each
  (8 pair-pT canvases + inclusive) + fit_report.txt + readback_check.txt, and
  `step{3,4}_dr_fit/plateau_guard_report.txt`.
  => the pipeline's `NPNG -lt 5` check is ALREADY wrong for the 8-bin nominal (it happens to pass).

## Design decisions
- D1. ONE place injects the sign into the histogram prefix: `MakeStepCfg(step, sign)` builds
  `h_prefix = "h_mc_dr_" + sign + "_"`.
- D2. TWO tags in the fit stage: `ptag` = `step<N>[_<sign>]` for READING the plateau file, and
  `tag` = `step<N>` for the keys WRITTEN into the fit file. The sign already lives in the fit FILE
  name, so the internal keys stay byte-identical across signs -> plot stage and the pipeline
  artefact check need no per-sign key logic.
- D3. `DrCorrFitFile(..., sign="")` appends `_<sign>` before `.root`; default keeps every existing
  call site byte-for-byte identical.
- D4. Text reports stay where they are today (`<method>/fit_report*.txt`,
  `step<N>_dr_fit/plateau_guard_report*.txt`, `<method>/readback_check*.txt`); the two new
  subdirectories `sign_intgr/` and `sign_sepr/` hold ONLY PNGs.
  Sign-separated report names: `*_same_sign.txt` / `*_opposite_sign.txt` (physical wording).
- D5. Guard policy: the |plateau-1| > 0.15 FAILURE tier throws ONLY for a sign-INTEGRATED series of
  a FULL sample. A sign-separated series measures, reports and lists everything but never throws --
  same sign is ~13% of the pairs, so per-sign failures are a STATISTICS statement about that
  subsample, not a defect of the nominal correction. `DrCorrPlateauUsable` and `fit_ok` are
  unchanged for all three series, so nothing unusable is ever published.

## Step 1 -- BASELINE captured before any edit (DONE)
- Ran the UNMODIFIED `fit_dr_corrections` on the CURRENT (2026-08-11) plateau + hist files, pp_full,
  tight+medium, steps 3+4, methods expo/polyu_fixedRp/interp, `allow_plateau_violation=true`
  (needed: see the guard note below). Reports copied to the session scratchpad
  `.../scratchpad/baseline_newinputs/`.
- All 12 fit reports are **byte-identical** to the Aug-7 reports already on disk
  => the orchestrator's per-sign refill did NOT move the sign-integrated nominal at all.
- **Pre-existing guard state (NOT caused by this work):** pp_full step 3 already FAILS the
  |plateau-1| <= 0.15 tier in 12 cells (50-72 GeV eta[1.0,1.5]; the 72-104 and 104-150 GeV rows),
  step 4 in 2 cells. Identical cells and identical numbers before and after. The pipeline default
  (STRICT_GUARD=0) reports it and re-runs with the explicit override.

## Step 2 -- code changes (DONE)
### dr_correction_sample_cfg.h
- `DrCorrSignText(sign)` -> "" / "same sign" / "opposite sign" (throws on an unknown token).
- `DrCorrSignFileTag(sign)` -> "" / "_same_sign" / "_opposite_sign" (report file names).
- `DrCorrFitFile(s, wp, step, method, sign = "")` appends `_ss` / `_os`; default keeps every
  existing call site byte-for-byte (verified: the sign-integrated fit file name is unchanged).
### fit_dr_corrections.cxx
- `MakeStepCfg(step, sign)` -- the ONE place the sign enters the histogram prefix.
- `fit_dr_corrections(..., allow_plateau_violation, sign = "")`, `fit_dr_corrections_all(..., sign)`.
- `ptag` (plateau keys, sign-tagged, incl. `prov_step3_ss`) vs `tag` (fit-file keys, NOT tagged).
- Availability probe up front: a sign series whose per-sign plateau key or per-sign histogram is
  missing prints `SKIPPED: ...` and RETURNS (non-fatal).
- `guard_is_fatal = cfg.is_full_sample && sign.empty()`; the report carries a paragraph stating the
  policy, and the verdict for a per-sign series is `reported, not enforced (<series> series)`.
- Reports: `plateau_guard_report[_same_sign|_opposite_sign].txt` in `step<N>_dr_fit/`,
  `fit_report[_same_sign|_opposite_sign].txt` in `<method>/`.
### plot_dr_correction_fits.cxx (rewritten around a `Series` vector)
- New trailing `mode` parameter: `sign_intgr` | `sign_sepr`; output PNGs go to
  `step<N>_dr_fit/<method>/<mode>/`. Text artefacts stay in `<method>/`.
- sign_sepr overlays same sign (blue curve kBlue, dark-blue kBlue+3 square markers) and opposite
  sign (red kRed, dark-red kRed+3 round markers); filled = inside the fit domain, open = outside.
- Legend: 4 entries (`same sign`, `opposite sign`, `fit, same sign`, `fit, opposite sign`) in its
  OWN row of the canvas header strip -- inside a panel it covered the measurement.
- Per-panel annotation: shared equation + shared FIXED parameters once, then one colour-coded
  column per sign (name, plateau, free parameters, chi2/ndf). The block is placed in whichever
  half of the frame that panel's markers leave free (off-scale points excluded -- they are a thin
  arrow, not a marker), with a clamp so it can never run into the axis labels.
- **PRE-EXISTING BUG FIXED:** `hist_path` now includes `MCTrigEffPairPt::FileSuffix()`. Without it
  a `MCTRIGEFF_PAIRPT_4BIN=1` run drew the NOMINAL 8-bin histograms against the 4-bin fit file.
- Read-back audit now loops over the series and writes
  `readback_check[_same_sign|_opposite_sign].txt`.
### pipelines/run_dr_correction_fits.sh
- `SIGNS` env var (default `"intgr ss os"`); Stage 2 loops over signs, Stage 3 over MODES
  (`sign_intgr` if the integrated fit exists, `sign_sepr` only if BOTH per-sign fits exist).
- Guard-verdict grep restricted to `sgn == intgr` (a per-sign verdict is never a pipeline failure);
  its `sed` range fixed from the non-existent `^cells with` to `^FAILING cells`.
- A missing per-sign fit whose log says `SKIPPED:` is a printed note, not an ARTEFACT_FAILURE.
- **PNG count DERIVED**, not a literal: `fit_file_npt()` reads `h_step<N>_plateau`'s x-axis nbins
  from the fit file, expected = nbins + 1. (The old literal 5 was the 4-bin count; the 8-bin
  nominal makes 9, so the check was under-counting.)
- **SECOND PRE-EXISTING BUG FOUND AND FIXED:** `root_file_has_objects()` used `root -l -b -q` with
  a heredoc. `root -l -b -q` with no macro argument QUITS BEFORE READING STDIN, so the check never
  ran and the function returned 0 for EVERY input -- every "file missing/incomplete" artefact check
  in this pipeline was a no-op. Verified: `-q` + `gSystem->Exit(7)` on stdin exits 0, without `-q`
  exits 7. Both this helper and the new `fit_file_npt` now use `root -l -b` + explicit
  `gSystem->Exit`. Re-tested: present objects -> pass, missing object -> fail, missing file -> fail.

## Incident (2026-08-11 01:08, no damage)
While testing a shell helper I `source`d `run_dr_correction_fits.sh`, which RAN it with defaults
and therefore executed Stage 1 (`plot_mc_trig_eff`) for pp_full/tight before the 2-minute tool
timeout killed it. It re-wrote step1_*, step2_* and part of step3_dr_correction/ -- all
deterministic regenerations from UNCHANGED inputs with UNCHANGED code. Checked afterwards:
no ROOT file in the sample dir was touched (`dr_correction_plateaus_pp24_full.root` still 00:39),
all 506 PNGs under pp_trigger_efficiency/ are complete (valid PNG header + IEND), no orphan ROOT
process left. Nothing lost; recorded here for the orchestrator. All later runs use SKIP_MEASURE=1.

## Step 3 -- full run for pp_full (DONE)
Command (SKIP_MEASURE=1 on purpose -- Stage 1 / plot_mc_trig_eff is the orchestrator's step and the
plateau file was already current):
`SAMPLES=pp_full WPS="tight medium" STEPS="3 4" METHODS="expo polyu_fixedRp interp" SIGNS="intgr ss os" SKIP_MEASURE=1 bash pipelines/run_dr_correction_fits.sh`
Exit code 2 = plateau-guard failures only (the PRE-EXISTING sign-integrated ones). **No ARTEFACT
failures**: every fit file complete, every PNG count met the DERIVED expectation, every read-back
report shows `0 FAILED persistence`. All 36 series (2 WP x 2 steps x 3 methods x 3 signs) produced.
Log: session scratchpad `pipeline_run.log`.

### 3a. Directory layout (both WP trees, both steps, all 3 methods -- 12 method dirs)
```
step{3,4}_dr_fit/plateau_guard_report[_same_sign|_opposite_sign].txt
step{3,4}_dr_fit/<method>/fit_report[_same_sign|_opposite_sign].txt
step{3,4}_dr_fit/<method>/readback_check[_same_sign|_opposite_sign].txt
step{3,4}_dr_fit/<method>/sign_intgr/  9 PNGs  (8 pair-pT canvases + inclusive)
step{3,4}_dr_fit/<method>/sign_sepr/   9 PNGs
```
= 108 sign_intgr + 108 sign_sepr PNGs. The 9 SUPERSEDED PNGs still sitting directly in each
`<method>/` are for the ORCHESTRATOR to delete (108 files); nothing writes there any more.

### 3b. Sign-integrated results are UNCHANGED (regression proof)
- All 12 `fit_report.txt` files: numeric rows byte-identical to the pre-change baseline
  (`diff` of everything except `#` comment lines is empty). Only the header comment lines changed
  (`series=sign-integrated`, the plateau-key/histogram provenance lines).
- All 108 `sign_intgr/*.png`: **md5-identical** to the pre-existing PNGs in `<method>/`.

### 3c. Same sign vs opposite sign -- INCLUSIVE cell, Tight WP
`f(0)` is the correction at dR = 0 (1 = no correction); the plateau is the large-dR normalization.
| step | method | plateau ss / os | f(0) ss / os | chi2/ndf ss / os |
|---|---|---|---|---|
| 3 | expo          | 0.9889 / 0.9933 | **0.3041 / 0.8155** | 1.31 / 9.80 |
| 3 | polyu_fixedRp | 0.9889 / 0.9933 | **0.1543 / 0.8409** | 2.13 / 8.12 |
| 3 | interp        | 0.9889 / 0.9933 | **0.3570 / 0.8357** | n/a |
| 4 | expo          | 0.9950 / 0.9967 | 1.0888 / 1.1308 | 1.39 / 42.59 |
| 4 | polyu_fixedRp | 0.9950 / 0.9967 | 0.9412 / 1.1187 | 1.22 / 51.93 |
| 4 | interp        | 0.9950 / 0.9967 | 1.0465 / 1.1476 | n/a |
Medium WP behaves the same (step-3 expo inclusive chi2/ndf 1.11 ss vs 10.63 os, etc.).

**PHYSICS FINDING -- the two signs do NOT agree in Step 3.** The pair-level correction at dR -> 0 is
0.30-0.36 for same sign but 0.82-0.84 for opposite sign: the same-sign turn-on is far deeper AND
narrower (expo lambda 0.137 vs 0.258). The sign-INTEGRATED result (f(0) = 0.8119, lambda 0.2530) sits
essentially on top of the opposite-sign one, as it must -- opposite sign is ~87% of the selected
pairs. Step 4 (single leg) agrees between the signs at the few-% level (f(0) 1.09 vs 1.13 for expo),
so the disagreement is specific to the PAIR (two-leg) term. Plateaus themselves agree well
(0.9889 vs 0.9933 step 3; 0.9950 vs 0.9967 step 4), i.e. this is a small-dR SHAPE difference, not a
normalization one. FOR THE ORCHESTRATOR TO JUDGE -- I did not change the nominal in any way.
Caveat on the chi2 comparison: opposite sign has ~7x the statistics, so its error bars are ~2.6x
smaller and the same functional form is judged much more harshly (chi2/ndf 9.8-52 vs 1.2-2.1).

### 3d. Guard outcome per series (all pp_full, FULL production)
| tree | step | sign-integrated | same sign | opposite sign |
|---|---|---|---|---|
| mc_based | 3 | **FAIL (FULL sample)** 12 failing / 5 flagged / 1 unmeasurable | reported, not enforced -- 13 / 4 / 2 | reported, not enforced -- 13 / 4 / 2 |
| mc_based | 4 | **FAIL (FULL sample)** 2 / 6 / 0 | reported, not enforced -- 2 / 5 / 1 | reported, not enforced -- 3 / 4 / 0 |
| mc_based_medium | 3 | **FAIL (FULL sample)** 12 / 6 / 0 | reported, not enforced -- 13 / 5 / 2 | reported, not enforced -- 11 / 6 / 2 |
| mc_based_medium | 4 | **FAIL (FULL sample)** 2 / 8 / 0 | reported, not enforced -- 2 / 7 / 0 | reported, not enforced -- 2 / 6 / 0 |
The sign-integrated FAIL is PRE-EXISTING (identical cells and numbers before any of this work; the
driver re-runs those with the explicit override, as it always has, and still exits 2). The
per-sign series measure and list exactly the same way and never throw -- the designed behaviour.
`fit_ok` and `DrCorrPlateauUsable` behave identically in all three series.

## Step 4 -- plot iteration from LOOKING at the PNGs (DONE)
Four rounds of look-and-fix on the rendered sign_sepr canvases (nothing here touched sign_intgr,
which stayed md5-identical throughout):
1. The 4-entry legend collided with the headline when placed beside it -> it now has its OWN row in
   a taller header strip (160 px on the grid canvases, 80 px on the inclusive one).
2. The parameter columns landed on the markers -> the annotation block now picks, PER PANEL, the
   half of the frame the markers leave free (off-scale points excluded; they are a thin arrow at
   the frame edge), with a clamp that keeps it off the axis labels when no half is free.
3. The `interp` equation quoted R_p without ever giving its value -> R_p is now drawn, READ FROM
   THE FIT FILE's provenance stamp (never retyped), whenever the equation names R_p and no drawn
   parameter carries it. For `polyu_fixedRp` R_p is a FIXED parameter, so it is stated once with
   the equation instead of twice in the two columns.
4. In the sparsest cells the left column ran into the right one -> columns moved to NDC x 0.42 /
   0.71 at text size 0.022 (sized by the widest line these columns ever carry).
   Also added a terse per-column `no fit` when one sign has a usable plateau but too few points to
   fit, so a column is never silently just a plateau.
Also: the mode directory is created only after the inputs load, so a skipped `sign_sepr` leaves no
empty directory behind.

**Canvases inspected with the Read tool (final versions):** step3 expo pair-pT (Tight), step3
interp pair-pT (Tight), step3 expo pair-pT 104-150 GeV (Medium -- the sparsest cell, with two
fully-unmeasurable panels), step4 polyu_fixedRp pair-pT (Tight), step3 expo inclusive (Tight),
step4 expo inclusive (Tight), step4 interp inclusive (Tight). No clipped legends, no text off the
frame, no column-on-column overlap. The only residual overlap is a marker drawn behind text in the
one Medium 104-150 GeV panel whose points fill the entire frame; the numbers stay readable and
there is no free band anywhere in that panel.

## Step 5 -- graceful skip verified (DONE)
`fit_dr_corrections("overlay", true, 3, "expo", false, "ss")` and
`plot_dr_correction_fits("overlay", true, 3, "expo", "sign_sepr")` both print an explicit SKIPPED
note naming the missing key/file, return normally (rc 0) and write nothing. `overlay` therefore
still runs end-to-end in sign-integrated mode.

## Handover to the orchestrator
1. **Delete the 108 superseded PNGs** sitting directly in
   `{pp}_trigger_efficiency/mc_based{,_medium}/step{3,4}_dr_fit/<method>/*.png` -- they are the old
   sign-integrated canvases, now reproduced byte-identically under `<method>/sign_intgr/`. Nothing
   writes to `<method>/*.png` any longer. (Keep the .txt reports there.)
2. **Physics finding to judge (§3c):** same sign and opposite sign disagree strongly in the Step-3
   pair correction at small dR (f(0) ~ 0.30-0.36 vs ~0.82-0.84 inclusive); Step 4 agrees to a few %.
   I changed nothing about the nominal.
3. **Two pre-existing bugs fixed** (both listed in Step 2): the missing `MCTrigEffPairPt::FileSuffix()`
   in the plot macro's `hist_path`, and `root_file_has_objects()` in the driver being a no-op
   because `root -l -b -q` never reads stdin.
4. The sign-integrated plateau guard still FAILS for pp_full (pre-existing, unchanged, 12 cells in
   step 3 and 2 in step 4); the driver exits 2 for that reason alone.
