# Analysis Code Review Log
**Task**: Review the D3 change to plot_pp24_fullsim_pair_bin_migration.cxx: stability -> bin-by-bin efficiency (same formula), purity denominator N[R_j] -> sum_{i in T} M_ij, o_j removed from plots (kept in txt), C_j unchanged. Physics Procedure §2 + §4 + D3.
**Log file**: review-analysis-code-20261008-151400-bin-migration-d3-purity.md
**Started**: 2026-10-08T15:14:00-04:00
**Status**: APPROVED
**Iterations completed**: 4
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 4 (2 WARNING, 2 INFO)
**Details**:
1. [DOC/LABEL] tracking doc §3(b) still says ε/s/p (+o), cell purity/stability/ε — WARNING. Fix: update to ε/ε^bin/p, o table-only.
2. [PHYSICS-RESULTS/DOC] Results R1–R3 still use old purity/"stability"; R3 m-vs-ΔR conclusion (ΔR 20× larger low-purity fraction) reverses under D3 (ΔR 0.0000 vs m 0.0002; min purity ΔR 0.846 > m 0.816); ΔR edge feed-in now only in o_j/C_j — WARNING. Fix: rewrite R1–R3 with D3 numbers.
3. [COMMENT] plotter :152 "every filled truth row inside T" — for dr guaranteed by filler throw, not inT — INFO.
4. [STYLE] C_j line default blue vs black markers — INFO (plot review).
**Numerical verification**: all MATCH (purities, bbeff, 3D means/mins, identity sum_{i∈T}M_ij = R_j − Rout_j max rel diff 0, row sums = TR, C_j = p(1−o)/(ε ε^bin)).
**Amend**: Issue 2 is documentation of correct results (not an unphysical result) → amend, no investigation. Fixing 1, 2, 3; 4 deferred to /review-plot.

## Iteration 2
**Reviewer verdict**: FAIL
**Issues found**: 3 (1 WARNING, 2 INFO)
**Details**:
1. WARNING: doc R2/R3 claim m_μμ o_j ≤ 0.05 everywhere contradicts table (48 m-cells > 0.05; max 0.090 for N_eff ≥ 100, top mass bin forward η). Fix: restate.
2. INFO: R1 Medium "within ≈0.002" → ≈0.003.
3. INFO: title/§1 still "stability".
**Numerical verification**: 25 doc numbers vs txt — MATCH except issue 1; C identity on all 25 1D bins within 3e-4.
**Amend**: all three fixed in the doc (Objective line 10 keeps the user's original 2026-10-06 wording).

## Iteration 3
**Reviewer verdict**: FAIL
**Issues found**: 5 (2 WARNING, 3 INFO)
**Details**:
1. WARNING: "38 of 48 cells o>0.05" wrong → Tight 42 / 36 top-mass (executor recount confirms 42/36).
2. WARNING: ΔR max o for N_eff ≥ 100 is 0.179 (pT 12.8–18.2, η −2.2…−2.0, ΔR 0.4–0.6, N_eff 317), 2× m_μμ max 0.090 → "same size" wrong.
3. INFO: "≤ 0.09" → max 0.0903. 4. INFO: "forward η" → mostly forward. 5. INFO: second ΔR high-o group pT 36.7–52.2, ΔR 0.2–0.4, N_eff ≤ 21.
**Numerical verification**: all other R1–R3 numbers MATCH.
**Amend**: R3 and R2 rows rewritten with recounted values (42/36; max o 0.18 ΔR vs 0.090 m; mostly forward; added low-N_eff row incl. second ΔR group).

## Iteration 4
**Reviewer verdict**: PASS
**Issues found**: 2 INFO (R3 "0.3–0.7" → "0.25–0.73"; m max-o cell η label) — both applied.
**Numerical verification**: every R1–R3 number MATCH (Tight + Medium).

**Status**: APPROVED at iteration 4
**Summary**: D3 purity/bin-by-bin-efficiency change correct in code (exact closures); iterations 1–3 fixed tracking-doc text (§3b, R1–R3, feed-in counts).
