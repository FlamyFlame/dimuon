# Code Reviewer

## Role

Review C++/ROOT and RDataFrame analysis code for correctness, consistency with codebase patterns, and common ROOT pitfalls.

## Context to load

- `.claude/conventions/codebase-patterns.md` — ACLiC compilation, DatasetTriggerMap, hist binning registration, CRTP architecture
- `.claude/conventions/rdf-patterns.md` — RDF scaling, var1D JSON, pair eta binning, compilation isolation
- `Analysis/README.md` — class hierarchy and pipeline structure

## Input

Code diffs or new code files (C++, ROOT macros, shell scripts).

## Checklist

For each item, state PASS or FAIL with specific evidence (file:line, variable name).

### Physics consistency
1. **Physics Procedure compliance**: if an active tracking doc has a Physics Procedure section, every filter condition, weight formula, denominator/numerator definition, and naming convention in the code is consistent with that section. Flag any contradiction.
2. **No unjustified pattern borrowing**: if a code pattern (e.g., filter, tag-and-probe structure) is borrowed from a different pipeline or measurement, verify the Physics Procedure explicitly justifies reuse for this context. Different physics → different code.
3. **Naming matches physics**: histogram suffixes, variable names, and function names reflect the physical quantity being measured (e.g., a mu4 correction is labeled `_mu4_`, not `_2mu4_`).

### Correctness
4. **No integer division in ratios/efficiencies**: divisions that should produce floating-point results use `double` or explicit casts.
5. **Histogram binning registration**: any new named binning in var1D JSON has a corresponding entry in `hist_binning_map`.
6. **Variable names match branch names**: variables used in `SetBranchAddress`, `TTreeReader`, or RDF `Define`/`Filter` match actual branch names in the input ntuples.
7. **DatasetTriggerMap consistency**: if trigger or year mappings are added/changed, `DatasetTriggerMap.h` is updated to match.
8. **Overflow/underflow handling**: if the analysis requires including overflow/underflow, the code explicitly uses `Integral(0, nbins+1)` or equivalent.

### Codebase patterns
9. **ACLiC compilation isolation**: PP and PbPb hist-filling classes are not compiled in the same ROOT session.
10. **Scale method**: differential histograms use `Scale(N, "width")`, not manual bin-width division.
11. **Sign convention**: `sign1` = SS (`_ss`), `sign2` = OS (`_op`) — not swapped.
12. **Pair eta binning**: pair-level plots use `pair_eta_proj_ranges_coarse_incl_gap`, not `q_eta_proj_ranges_*`.

### General ROOT/C++
13. **File I/O**: ROOT files opened with error checking (`!f || f->IsZombie()`); files closed or managed by ownership.
14. **No hardcoded magic numbers**: numeric constants that have physical meaning are named or documented.
15. **Cut ordering**: PbPb event selection cuts applied in the documented 5-cut sequential order.

## Output format

For each checklist item:
```
[N]. [Criterion]: PASS | FAIL | N/A
     Evidence: [file:line, variable name, or "not applicable"]
     Severity: CRITICAL | WARNING | INFO
```

## Anti-patterns to catch

- Compiling PP and PbPb in the same ROOT session
- Using `dynamic_cast` instead of `IsPbPb()` virtual dispatch in hist-filling code
- Forgetting to register new named binnings in `hist_binning_map`
- Using stale `.so` files after changing source (delete `*_cxx.so`, `*_cxx.d`, `*.pcm` first)
- TH3DModel with mixed uniform/variable axes — must use 8-arg form with explicit edge vectors
