# Analysis Code Review Log
**Task**: Physics+code correctness of the Medium→TIGHT default muon-WP change (pp & PbPb: data crossx wiring, reco-eff Tight, trig-eff Tight turn-on). Nominal tight=unsuffixed; medium=_medium_wp systematic. All files compile clean.
**Log file**: review-analysis-code-20260708-010455-tight-wp-default.md
**Started**: 2026-07-08T05:04:55Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 1 WARNING + 2 INFO
**Details**:
1. [WARNING] Medium-WP trig-eff plumbing clobbers the Tight nominal + name-mismatched: `_medium_wp` suffix applied only when `!isTight && !trigger_effcy_calc` (Data.cxx:106), so a Medium trig-eff run (trigger_effcy_calc=true) writes the UNSUFFIXED graph file → clobbers Tight; and `OpenEffcyPtFitFile` always reads the unsuffixed fit. (Medium systematic only; nominal Tight unaffected.)
2. [INFO] FillHistogramsGeneric control-hists not WP-filtered but use Tight reco-eff keys (diagnostic only, pre-existing).
3. [INFO] pp Tight reco-eff = interim reuse of Medium Fig.31 (labeled, no invented numbers) — accepted gap.
**Reviewer confirmed:** NOMINAL Tight is fully correct + self-consistent (crossx spectrum + reco-eff Tight keys + trig-eff tag-and-probe all Tight; `node` reference valid; only the quality bit changed; C5/C6 unaffected; non-clobber for crossx). Placeholder has 65 tight + 65 medium keys.

**Amendment (iter 1→2):** Fixed the WARNING. `Data.cxx:106` → `if(!isTight) out_file_suffix += "_medium_wp"` (Medium suffix now on ALL Medium runs incl. trig-eff graph output). `OpenEffcyPtFitFile` PP+PbPb → read `single_mu_effcy_pT_fit{wpsuf}.root` + `..._fine_q_eta_bin{wpsuf}.root` with `wpsuf = isTight?"":"_medium_wp"`. Now Medium chain is non-clobbering + name-consistent: RDF medium trig-eff→`..._medium_wp` graph → fitter(`_medium_wp`) → `..._fit_medium_wp.root` → consumed medium. Nominal Tight unchanged (unsuffixed). PP recompiles clean; PbPb compiling. INFO #2/#3 left as-is (diagnostic / accepted gap).

## Iteration 2
**Reviewer verdict**: PASS
**Issues found**: 0 (None found)
**Details**: Focused re-review of the medium-suffix amendment. Nominal Tight byte-identical (isTight=true → wpsuf="" everywhere). Medium round-trip verified non-clobber + name-matched: producer `..._fine_q_eta_bin_medium_wp.root` → fitter(_medium_wp) → `single_mu_effcy_pT_fit_medium_wp.root` → consumer reads same. Tight vs Medium filenames disjoint. No new bug (pure string concat).

**Status**: APPROVED at iteration 2
**Summary**: Medium→Tight default WP change verified physically correct + self-consistent for the Tight nominal (crossx + reco-eff + trig-eff all Tight), and the Medium systematic plumbing is now non-clobbering. 0 CRITICAL/0 WARNING. (INFO: pp-tight reco-eff = labeled interim; FillHistogramsGeneric control-hists Medium-selected — both accepted/known.)
