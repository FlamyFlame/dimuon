# Analysis Code Review Log
**Task**: Decouple the ad-hoc dR-fallback from `pythia_only_barcode_cache`; add `use_dr_fallback` flag (default false) so the DEFAULT overlay matching is pure prob>0.5 (no dR, no geometric). Preserve pythia_only_barcode_cache's cache-restriction role.
**Log file**: review-analysis-code-20260612-041807-decouple-dr-fallback-default-off.md
**Started**: 2026-06-12T04:18:07Z
**Status**: APPROVED at iteration 1
**Iterations completed**: 1
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 1 INFO (no CRITICAL/WARNING)
**Details**: INFO — helper script `run_pythia_fullsim_overlay_r17662_nodr.sh` had a stale `pythia_only_barcode_cache=false` line and misleading comment (no-op after decoupling, also overridden by InitParamsExtra). Fixed: removed the line, updated comments to state default = pure prob>0.5 (no dR), use_dr_fallback=true for dR.
**Numerical verification**: All MATCH — default reco_match=0.8266/pass_med=0.7426; use_dr_fallback=true reco_match=0.8955/pass_med=0.7430 (6-8 GeV, N=2918). dR recovers reco_match but pass_medium unchanged.

**Key discovery during review**: the original `pythia_only_barcode_cache` dR gate was UN-disableable from run macros because `InitParamsExtra()` (called in `Run()` via `InitParamsHook()`) re-sets it true. So all prior "nodr" r17662 runs (Steps 20-23) actually ran WITH the dR fallback. The new `use_dr_fallback` flag (untouched by InitParamsExtra) correctly defaults off → true pure prob>0.5. The geometric-matching refactor (prior review) was REVERTED — it appeared to "regress" reco_match only because it exposed the true no-dR rate (0.8266) that the old code could never reach; pass_medium (the efficiency) is identical (~0.743) across pure-prob / +dR / geometric, so all physics conclusions hold.

**Status**: APPROVED at iteration 1
**Summary**: Added `use_dr_fallback` (default false) and re-gated the dR fallback on it; default overlay matching is now pure prob>0.5 with no dR, no geometric. Verified on r17662. Geometric refactor reverted (not needed in production; physics result stands).
