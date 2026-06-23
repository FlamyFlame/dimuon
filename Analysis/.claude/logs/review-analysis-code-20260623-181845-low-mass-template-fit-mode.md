# Analysis Code Review Log
**Task**: Add low_mass_template_calc data-RDF mode (PP+PbPb): read _no_res_cut, distinct output, fill 0-4 GeV OS+SS minv template spectra (1D + 2D pT/eta, PbPb per-ctr), dsigma-weighted. D_OS/D_SS for the 5b closure + fit.
**Log file**: review-analysis-code-20260623-181845-low-mass-template-fit-mode.md
**Started**: 2026-06-23T22:18:45Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Executor work (pre-review)
Files modified:
- RDFBasedHistFillingData.h: added PUBLIC `bool low_mass_template_calc = false;` (mode flag).
- RDFBasedHistFillingData.cxx TriggerModeSettings: out_file_suffix += "_template_fit" when low_mass_template_calc.
- RDFBasedHistFillingPP.cxx: SetIOPathsHook new top branch (low_mass_template_calc -> _no_res_cut, throw if absent);
  FillHistogramsCrossx early-return template block (after npt setup) filling h1d_crossx_minv_0_4_{op,ss}_dsigma +
  2D {pair_pt_log_150,minv}/{pair_eta,minv} _{op,ss}_dsigma, dsigma weight, selection signal_cuts MINUS minv (no dR), return.
- RDFBasedHistFillingPbPb.cxx: SetIOPathsHook low_mass branch in BOTH 23/24/25 and 15/18 blocks;
  FillHistogramsCrossx early-return per-centrality template block (after dsigma_lumi_factor) filling the same
  histos with _<ctr> suffix.
- New run scripts: run_template_fit_pp24.sh, run_template_fit_pbpb2{3,4,5}.sh (low_mass_template_calc=true,
  output_generic_hists=false [light], same trigger_mode/mu4_nominal as nominal crossx).

INCIDENT + RECOVERY (during execution): I first placed low_mass_template_calc in the PROTECTED section of
RDFBasedHistFillingData.h. The cling macro assignment `pp.low_mass_template_calc = true;` then FAILED SILENTLY
(non-fatal cling error) -> stayed false -> pp.Run() executed in NOMINAL mode (trigger_mode=3, V1) with
output_generic_hists=false and OVERWROTE histograms_real_pairs_pp_2024_2mu4_nominal.root (dropped generic/gapcut
histos, 678KB->658KB). SAME root cause as the 2026-06-22 output_single_muon_tree incident. RECOVERY: (1) moved the
flag to PUBLIC (next to mu4_nominal_pbpb_NO_trig_calc), recompiled PP, verified `pp.low_mass_template_calc=true`
now compiles (FLAG_SETTABLE=1); (2) RESTORED the nominal pp output by re-running run_crossx_hist_filling_pp24.sh
(canonical nominal, output_generic_hists=true) -> histograms_real_pairs_pp_2024_2mu4_nominal.root 678KB, intact.
No lasting damage (nominal output reproduced from V1). LESSON: any member set from a cling run-script macro MUST
be public; protected fails silently and runs in the default (nominal-overwriting) mode.

Compilation: PP + PbPb ACLiC clean (separate sessions; PP .so 18:3x, PbPb .so 18:37).
pp24 template run: SUCCESS. Output histograms_real_pairs_pp_2024_2mu4_nominal_template_fit.root (distinct).
Verification (pp24): 1D OS integral=17130, SS=935.4; OS HAS RESONANCES PRESENT [0,1.06]=3330, [2.9,3.3]=7131
(both EXACTLY 0 under V1 -- the goal); SS smooth (140/102 continuum); all 4 2D histos present (15x50 pT, 24x50 eta),
integrals match 1D. Nominal pp output untouched by the template pass (distinct file).
PbPb 23/24/25 template runs: IN PROGRESS (background b3pe0qk84) -- will verify per-centrality outputs before review.

PbPb 23/24/25 template runs: SUCCESS (exit 0, "template mode completed", distinct outputs
histograms_real_pairs_pbpb_20YY_single_mu4_no_trg_plots_nominal_template_fit.root 324/315/340 KB).
Verification (pbpb23): per-ctr OS has RESONANCES PRESENT (ctr0_5 [0,1.06]=15950 J/psi[2.9,3.3]=21130;
ctr50_80 3779/6949); SS smooth (ctr0_5 SS=27030 large central combinatoric, expected); 2D per-ctr present.
NOMINAL PbPb outputs UNTOUCHED (mtime Jun 22 01:28, distinct template files only). All 4 datasets done.
Numbers reported: pp24 OS=17130 SS=935.4; pbpb23 ctr0_5 OS=78090 SS=27030; resonances present in OS only.

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 0 CRITICAL/WARNING (2 INFO, non-blocking):
 - INFO1: nominal output STILL carries the V1-vetoed h1d_crossx_minv_0_4 (the old end-of-function block runs
   in nominal; my early-return only fires in template mode). Harmless — different file (_template_fit), nominal
   V1 version ~0 in veto windows, unused. (My executor note "gated out of nominal" was imprecise; code correct.)
 - INFO2: PbPb template block Defines q_eta1 unconditionally (PP guards on generic_weight_col); safe because
   run scripts set output_generic_hists=false (no pre-defined columns). Optional symmetry fix.
**Details**: All numbers MATCH (<0.03% on totals; window sub-integrals differ a few % only because 1.06/3.3
aren't bin boundaries). Flag PUBLIC confirmed (incident fixed). _no_res_cut only in template mode; Part-1
nominal/trig-eff branches unchanged. Distinct output. OS=df_op/SS=df_ss correct (OS>>SS, OS has phi/J/psi, SS
doesn't). dsigma weight not T_AA. No minv/dR cut. Binnings real. C1/C2 OK (smooth continuum + resonance peaks).
**Numerical verification**: all MATCH.

**Status**: APPROVED at iteration 1
**Summary**: low_mass_template_calc mode added (PP+PbPb): reads _no_res_cut, distinct _template_fit output,
fills D_OS/D_SS 1D + 2D (pair pT/eta), PbPb per-ctr. OS resonances present (the goal); nominal untouched.
Protected-member incident fixed (flag public) + nominal pp restored. Reviewer PASS iter 1.
