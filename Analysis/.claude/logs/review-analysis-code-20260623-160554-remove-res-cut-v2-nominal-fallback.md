# Analysis Code Review Log
**Task**: Remove the obsolete _res_cut_v2 / _no_res_cut / bare-.root fallback from the NOMINAL/crossx input-file selection in the data RDF (PP + PbPb). Nominal must require ONLY V1 _mindR_0_02.root and throw a descriptive error if absent. Trigger-efficiency branch unchanged.
**Log file**: review-analysis-code-20260623-160554-remove-res-cut-v2-nominal-fallback.md
**Started**: 2026-06-23T20:05:54Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Executor work (pre-review)
Files modified:
- `RDFBasedHistFilling/RDFBasedHistFillingPP.cxx` SetIOPathsHook (L15-45): nominal branch
  (`!trigger_effcy_calc`) now resolves ONLY V1 `pp_base+input_mindR_suffix+".root"`; throws a
  descriptive runtime_error if absent (states V1 _mindR_0_02 required, fallback removed). Trig-eff
  `else` branch candidate list (V2 first) UNCHANGED.
- `RDFBasedHistFilling/RDFBasedHistFillingPbPb.cxx` SetIOPathsHook: BOTH nominal branches
  (`mu4_nominal_pbpb_NO_trig_calc`, run_year 23/24/25 block and 15/18 block) now resolve ONLY V1
  `base+input_mindR_suffix+".root"`; throw descriptive error if absent. Both trig-eff `else` lists
  (V2 first) UNCHANGED.
Compilation (ACLiC, separate sessions per PP/PbPb isolation):
- PP: fresh `.so` 16:09, instantiates (run_year=24), exit 0, no errors (only pre-existing
  ParamsSet.h sign-compare warnings).
- PbPb: "creating shared library", PBPB_OK ctor (run_year=23), fresh `.so` 16:13, exit 0.
No crossx rerun (per task; orchestrator recompiles + sanity-runs after review).
Numbers: none computed (IO-path logic change only).

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 0
**Details**: None found. Verified: trig-eff branch untouched (all 3 places, `_res_cut_v2` first);
nominal branch restricted to V1 `_mindR_0_02` only; `gSystem->AccessPathName` throw sense correct
(throws when file absent); `in_path`/`output_file`/`else`-scoped vector all correct; descriptive
throw messages; grep confirms no external caller relied on the nominal hook returning a
`_res_cut_v2`/`_no_res_cut` file (other hits are independent raw-path readers).
**Numerical verification**: No numbers to verify.

**Status**: APPROVED at iteration 1
**Summary**: Nominal/crossx data-RDF input selection (PP + both PbPb blocks) now requires V1
`_mindR_0_02.root` only and throws a descriptive error if absent; obsolete `_res_cut_v2`/`_no_res_cut`
fallback removed; trigger-efficiency branch unchanged. ACLiC-clean both classes. Reviewer PASS iter 1.
</content>
