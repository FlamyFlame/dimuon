# Analysis Code Review Log
**Task**: Stage-2 prompt-content decomposition of bkg pairs (prompt_plus_bkg vs bkg_bkg) in bkg_mc_provenance.C; test near-side charge symmetry of genuinely-combinatoric bkg_bkg.
**Log file**: review-analysis-code-20260624-135708-stage2-prompt-decomp.md
**Started**: 2026-06-24T13:57:08-04:00
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Execution (iteration 1 prep)
Extended bkg_mc_provenance.C: added pair classes prompt_plus_bkg (nprompt==1) and bkg_bkg (nprompt==0),
histos h_bkgmc_minv_{op,ss}_{prompt_plus_bkg,bkg_bkg}, per-class mirror printout + per-file bkg_bkg counts +
consistency check. ACLiC-clean. Ran overlay + signal.
OVERLAY: consistency PASS (prompt+prompt_plus_bkg+bkg_bkg=all=19903; prompt_plus_bkg+bkg_bkg=background=1838).
  background=1838 (prompt_plus_bkg=1805, bkg_bkg=33). MIRROR near-side[0,1.5): prompt_plus_bkg OS=242 SS=48
  SS/OS=0.198 (OS-enhanced); bkg_bkg OS=8 SS=2 SS/OS=0.250 (33 pairs total -> INCONCLUSIVE stats).
SIGNAL: bkg_bkg ~15 pairs total (no UE) -> negligible.
VERDICT(physics): H2 structural claim CONFIRMED (bkg is 98% prompt_plus_bkg carrying the OS-enhancement); the
decisive bkg_bkg charge-symmetry test is STATISTICALLY INCONCLUSIVE (pair_pt>8 suppresses two-soft-bkg pairs)
-> Stage 3 needs higher-stat UE/min-bias.

## Iteration 1
**Reviewer verdict**: PASS (0 CRITICAL/WARNING; 3 INFO: don't overreach on bkg_bkg [won't]; [0,1.5)->[0,1.44) label; counter-vs-integral labeling)
**Numerical verification**: ALL MATCH — class sums exact per OS & SS (prompt+ppb+bb=all; ppb+bb=background) in both overlay & signal; mirror integrals reproduce (ppb near-side OS=242 SS=48 SS/OS=0.198; bkg_bkg OS=8 SS=2).
**Status**: APPROVED at iteration 1
**Summary**: Stage-2 prompt-content split validated. background=98% prompt_plus_bkg (OS-enhanced near-side, SS/OS=0.198); bkg_bkg too rare (overlay 33 / signal ~15) to test charge symmetry. H2 structural claim confirmed; bkg_bkg decisive test inconclusive -> Stage 3 needs higher stats.
