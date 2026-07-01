# Analysis Code Review Log
**Task**: D1 muon Δp/p templates (real/fake/hadronic) from fullsim+overlay, inclusive + low-mass-pair subset; check within-jet vs inclusive hadronic dp/p shape.
**Log file**: review-analysis-code-20260625-020343-dpop-templates.md
**Started**: 2026-06-25T02:03:43-04:00
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Execution (iteration 1 prep)
Macro dimuon_data/plots/template_fitting/dpop_fake_muon_20260624/code/dpop_templates.C. ACLiC-clean. Ran signal+overlay.
Outputs dpop_templates_{signal,overlay}.root (h_dpop_{real,fake,hadronic} + _lowmasspair_*).
SIGNAL: real mean+0.014 rms0.071 (92.7% in cut); hadronic mean+0.092 rms0.185 (43% beyond). low-mass-pair hadronic
  mean+0.089, KS vs inclusive=1.000 (identical).
OVERLAY: real mean+0.016 (90.8% in cut); hadronic inclusive mean+0.077 rms0.242; low-mass-pair hadronic mean+0.027
  rms0.292, KS vs inclusive=0.000 (DIFFERENT — within-jet hadronic SOFTER/less-separable). fake overlay broad,
  mean-0.033 (low-mass-pair -0.070).
KEY (D4 input): in realistic UE the within-jet (low-mass-pair) hadronic Δp/p differs from inclusive (softer, less
shifted) -> the fit template must be the within-jet shape; Δp/p separation weaker for within-jet -> likely under-removes
-> pair-level residual template likely needed. Code correct; result is physics.

## Iteration 1
**Reviewer verdict**: PASS (0 CRITICAL/WARNING; 1 INFO: pT-hat slices unweighted — inherited design, note for the D2 yield fit where slice weighting can shift the template mixture)
**Numerical verification**: ALL MATCH (signal real mean+0.0136/0.927; hadronic mean+0.0916/43% beyond; lmp hadronic mean+0.089 KS=1.000; overlay lmp hadronic KS=0.000; sentinel handling confirmed).
**Status**: APPROVED at iteration 1
**Summary**: D1 Δp/p templates built (real/fake/hadronic, inclusive + low-mass-pair). Real peaks ρ≈0; hadronic shifted +0.08-0.09. KEY: overlay within-jet hadronic Δp/p SOFTER than inclusive (KS=0, mean +0.027 vs +0.077) -> fit must use the within-jet template; weaker separation -> Δp/p likely under-removes -> pair residual likely needed (D4).
