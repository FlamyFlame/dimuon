# Note Review Log
**Task**: datasets.tex pp24 only: lumi 400.412 -> 395.836 pb^-1; pp GRL -> physics_ppref2024_all_good.xml (v144-pro42-02); run 488456 LB 84-220 removed.
**Log file**: review-note-20261012-232317-datasets-pp24-grl-lumi.md
**Started**: 2026-10-13T03:23:17Z
**Status**: APPROVED
**Iterations completed**: 1
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 4 (all INFO): long source line (rewrapped); analysis_metadata.md:179 still says prescale ~1.0001 (compatible with 100.00 %); optionally add the v144 tag to the GRL table; no other tex file quotes 400.412.
**Numerical verification**: 395.836 pb^-1 (table, text) == CSV Total == PPBaseClass; prescale fraction 100.00 %; GRL tags; GRL diff only run 488456 (84 vs 221); LB drop matches PPEventSelConfig/PPExtras — all MATCH.

**Status**: APPROVED at iteration 1
**Summary**: datasets.tex pp24 lumi 395.836 pb^-1, new pp GRL and the run 488456 LB 84-220 removal documented; no uncertainty/recalibration claim; no Pb+Pb 2026 added (user scope).
