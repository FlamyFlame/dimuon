# Plot Review Log
**Task**: Review the D3-updated bin-migration plot set (reco efficiency / bin-by-bin efficiency / purity with purity = M_jj/sum_{i in T} M_ij; o_j lines removed), Tight + Medium.
**Log file**: review-plot-20261008-151427-bin-migration-d3-plots.md
**Started**: 2026-10-08T15:14:27-04:00
**Status**: APPROVED
**Iterations completed**: 2
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 6 (1 WARNING, 5 INFO)
**Details**:
1. [6] compare-plot right-panel y title clipped ("ction of N^truth ... bin-by-bin efficiency < x") — WARNING; fix: shorten title / widen margin.
2. [4] 1D: purity and bin-by-bin efficiency markers coincide — INFO; optional x offset.
3. [R1] 1D 0–1.15 range compresses p/ε^bin departures — INFO; none needed.
4. C_j line default blue — INFO.
5. cell-panel floor from central values; one Wilson bar touches floor — INFO.
6. [19] no on-plot formula for p/ε^bin — INFO; caption must state definitions.
**Numerical verification**: all MATCH (every Tight 1D value, 3D means both WPs, closures M flow rows 0, row sums = TR, Σ_{i∈T}M_ij + Rout = R, C identity 4e-16).
**Amend**: fix 1 (shorter y title + wider margin), 2 (±10% bin-width x offset of series), 4 (C_j line black), 5 (floor from v − e_lo).

## Iteration 2
**Reviewer verdict**: PASS
**Issues found**: 1 (INFO: dashed ΔR curve risers faint on compare plot; optional)
**Numerical verification**: pT bins 1/8, ΔR top bin, a cell point vs txt — MATCH; Tight/Medium distinct.

**Status**: APPROVED at iteration 2
**Summary**: D3-updated bin-migration plots (reco eff / bin-by-bin eff / purity, o_j removed) approved after fixing clipped compare-plot y title, marker overlap, C_j colour, cell floor.
