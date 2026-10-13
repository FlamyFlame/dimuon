# Analysis Code Review Log
**Task**: pp24: PPEventSelConfig.h (run 488456 LB 84-220 dropped), PPExtras PassEventSelExtra, PPBaseClass 2mu4 = 1/395.836, TrigRates_CA.py pp24 GRL -> physics_ppref2024_all_good.xml.
**Log file**: review-analysis-code-20261012-231355-pp24-grl-lb-exclusion.md
**Started**: 2026-10-13T03:13:55Z
**Status**: APPROVED
**Iterations completed**: 1
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 6 (all INFO): stale 400.412 in IntNotes/tex/datasets.tex:70,82; comment in plot_mc_data_pair_pt_in_eta.cxx:20; historical docs; min-bias path doesn't bind RunNumber/lbn (not used by pp24); sel_data.py diagnostic not audited; new lumitable has a uniform +0.266 % lumi-calibration shift in the 7 unchanged runs (+0.877 pb^-1) on top of the LB-drop (-5.453 pb^-1).
**Numerical verification**: GRL diff only run 488456 (84->221) MATCH; 974 LBs GRL==lumitable MATCH; Total 395.836 MATCH; 70.592->65.139 MATCH; PPBaseClass 1/395.836 MATCH; PPEventSelConfig logic MATCH.

**Status**: APPROVED at iteration 1
**Summary**: pp24 numerator drops run 488456 LB 84-220 and denominator is the new-GRL lumitable (395.836 pb^-1); hook wiring verified by reading.
