# Investigation Review Log
**Task**: 5b DATA closure (THE GATE) — develop a clean coupled OS+SS fit (shared N_C, wide clean-SS anchor, OS resonance mask, non-negativity); reach a defensible PASS/FAIL on SS_corr = k*OS_corr. First attempt: k_data=0.296 vs k_MC=0.308 (good) but N_C negative + bad chi2 (fit-setup issues).
**Log file**: review-investigation-20260623-230057-5b-closure-coupled-fit-gate.md
**Started**: 2026-06-24T03:00:57Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Executor findings (pre-review)
Hypotheses tested + evidence (pp24):
1. FIRST attempt (separate fits, [1.08,4]): k_data=0.296+/-0.010 vs k_MC=0.308 (ratio 0.96) BUT N_C NEGATIVE
   (SS -562, OS -1353), chi2/ndf OS=83. -> fit-setup issues, retest with fixes.
2. COUPLED NNLS (shared N_C, SS wide [0.3,4] clean anchor, OS [0.3,4] resonance-masked, yields>=0 via
   brute-force subset NNLS): integrated N_S=5773, **N_C(shared) pinned to EXACTLY 0**, N_G=2723, N_kG=887;
   k_data=0.326; **chi2/ndf=368** (catastrophic). N_C=0 is UNPHYSICAL: SS is combinatoric-DOMINATED, yet the
   fit uses NO combinatoric (G_SS absorbs the SS). The wide-range + coupling + NNLS did NOT rescue it.
3. DECISIVE shape diagnostic (SS, area-norm [0.3,4]): corr(T_mix, G_SS) = **0.978** (near-degenerate
   templates), corr(D_SS, T_mix) = **0.926**. From the per-bin table: T_mix is systematically HARDER than
   D_SS (high-mass excess, low-mass deficit). Physics check: D_SS = C + k*G and G RISES toward high mass, so a
   correct combinatoric C would make D_SS HARDER than C — but the data shows the OPPOSITE (D_SS softer than
   T_mix). => the mixed-event T_mix OVER-POPULATES high mass (the single-muon pool used for mixing is harder
   than the soft muons in real low-mass b-decay pairs). Plot: DIAGNOSTIC_SS_shapes_degeneracy.png.

CONCLUSION / VERDICT: the 5b data closure does NOT cleanly pass. TWO root causes, both real (not just fit
setup): (a) the mixed-event T_mix does NOT reproduce the data combinatoric minv shape (too hard); (b) the
combinatoric (T_mix) and correlated (G_SS) templates are near-degenerate in shape (corr 0.98), so even a
correct T_mix could not robustly separate combinatoric from correlated in the SS fit. The coupled wide-range
NNLS pins N_C to 0 with chi2/ndf=368 — unphysical. The promising integrated k_data~k_MC from attempt 1 was a
coincidence of a poorly-conditioned fit, NOT a clean closure.
RECOMMENDATION (gate): this is a closure FAIL/MARGINAL -> STOP and escalate to the user (the user's gate
condition). Options for the user: (A) refine the mixed-event method (mix within pT/eta classes so T_mix
matches the data combinatoric) then retry the combined fit; (B) switch to MC-only separate SS/OS fits
(code-set B, the pre-agreed fallback); (C) accept the combined fit only with a large systematic. The
near-degeneracy (corr 0.98) is a fundamental concern that even a perfect T_mix may not overcome -> (B) is the
likely robust choice. Quantitative support: numbers_pp.txt, numbers_pp_coupled.txt, the shape correlations.

## Iteration 1
**Reviewer verdict**: FAIL (on the diagnosis WRITEUP; the GATE OUTCOME — closure fails, escalate — was CONFIRMED correct)
**Issues**: (WARNING) my "corr" 0.978/0.926 were uncentered COSINE similarity, not Pearson -> overstated the
degeneracy root-cause. Proper Pearson: corr(T_mix,G_SS)=0.89, corr(D_SS,T_mix)=0.40. (INFO) lead with the
robust shape-INFEASIBILITY argument (data softer than BOTH templates) rather than degeneracy.
**Numerical verification**: all MATCH (cosine values reproduced; Pearson recomputed; N_C=0, chi2=367.74 MATCH).

## Iteration 1 fix (amended diagnosis)
Recomputed PEARSON correlations: corr(T_mix,G_SS)=0.895, corr(D_SS,T_mix)=0.402 (LOW -> T_mix genuinely
mismatches the data combinatoric), corr(D_SS,G_SS)=0.729. **Shape-infeasibility (PRIMARY evidence):** D_SS
peaks at m=1.72 GeV then FALLS; T_mix peaks at 3.96, G_SS at 3.64 -> the DATA IS SOFTER THAN BOTH background
templates. 41/46 bins in [0.3,4] have D_SS outside the [min,max] of the two templates (18 where D_SS<both at
high mass, 23 where D_SS>both at low mass). => NO non-negative combination N_C*T_mix + N_kG*G_SS can fit D_SS
=> the NNLS pins N_C=0 with chi2/ndf=368. This is a genuine shape-feasibility failure, NOT a fit-conditioning
artifact (error/range/NNLS tweaks cannot rescue it). Plot DIAGNOSTIC_SS_shapes_degeneracy.png re-annotated to
lead with this. CORRECTED VERDICT: the combined OS+SS fit does NOT close — the mixed-event T_mix is too hard
(Pearson 0.40 vs data; over-populates high mass) and, more decisively, the data SS combinatoric is SOFTER than
both background templates. Gate = FAIL -> escalate to user. Options: (A) refine mixed-event (pT-class mixing to
soften T_mix) + retry; (B) MC-only separate fits (code-set B, pre-agreed fallback); (C) combined + large syst.

## Iteration 2
**Reviewer verdict**: PASS
**Issues found**: 0
**Details**: All amended claims reproduce exactly (Pearson 0.895/0.402/0.729; peaks D_SS 1.72, T_mix 3.96,
G_SS 3.64; 41/46 bins outside envelope [18 high-mass deficit + 23 low-mass excess]; N_C=0, chi2/ndf=367.74).
Plot re-annotated (leads with shape-infeasibility, shows Pearson not cosine). Gate logic sound: combined fit
does NOT close -> escalate (option B MC-only is the documented fallback). Escalation is the correct call.
**Numerical verification**: all MATCH.

**Status**: APPROVED at iteration 2
**Summary**: 5b DATA CLOSURE FAILS (the GATE). The combined OS+SS template fit does not close for pp24: the
data SS combinatoric is SOFTER than both the mixed-event T_mix (too hard; Pearson 0.40 vs data) and the MC
G_SS (peaks 3.6-4.0 vs data 1.72), so no non-negative N_C*T_mix+N_kG*G_SS describes D_SS (N_C pinned to 0,
chi2/ndf=368). Root causes: (1) mixed-event T_mix over-populates high mass; (2) data softer than both bkg
templates. Gate condition met -> STOP + escalate to user. Options: (A) refine mixed-event (pT-class mixing);
(B) MC-only separate SS/OS fits (code-set B, fallback); (C) combined+large systematic.
