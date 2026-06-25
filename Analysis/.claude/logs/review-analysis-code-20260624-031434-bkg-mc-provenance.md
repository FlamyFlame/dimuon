# Analysis Code Review Log
**Task**: Standalone reco-seeded provenance macro on pythia fullsim → demonstrate low-mass background (fake/hadronic) + SS↔OS mirror (V2 of OS-SS foundation validation).
**Log file**: review-analysis-code-20260624-031434-bkg-mc-provenance.md
**Started**: 2026-06-24T03:14:34-04:00
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Execution (iteration 1 prep)
Macro: dimuon_data/plots/template_fitting/bkg_mc_provenance_20260624/code/bkg_mc_provenance.C (ACLiC-clean after adding <sstream>).
Ran over 24 pythia fullsim NTUP files. RESULTS:
- selected reco muons=364686: prompt 99.32%, fake 0.02% (61), hadronic 0.66% (2408).
- pairs: prompt 103576, fake 23, hadronic 2145, background 2168, all 105744 (prompt+bkg=all, fake+had=bkg: consistent).
- MIRROR (background): OS[0,4]=844 SS=343 SS/OS=0.406 ; OS[0,1.5]=629 SS=183 SS/OS=0.291.
- per-file bkg SS/OS clustered ~0.5-0.67 (high-pTH), consistent.
INTERPRETATION: signal-MC background = correlated HF (prompt+hadronic, OS-enhanced, SS/OS~k~0.3), NOT the
symmetric UE combinatoric of data SS. Check(2) [genuine fake/hadronic bkg surviving quality cuts] CONFIRMED;
check(3) [SS<->OS mirror] NOT demonstrable in low-UE signal MC -> needs HIJING overlay (UE) or data-driven T_mix.
Output: bkg_mc_provenance.root.

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 1 INFO (cosmetic: "[0,1.5]" printout is really [0,1.44) — bin-edge; no physics impact)
**Details**: All code-correctness checks passed (charge=sign(muon_pt); pt MeV→GeV; quality bitmask; z0 sin theta; provenance classify; disjoint pair classes prompt+bkg=all & fake+had=bkg; q*eta<2.2 both; vector<bool> read; size guards; double ratios; zombie checks). barcode>200k correctly avoided.
**Numerical verification**: ALL MATCH (op_background Integral=844, ss=343 SS/OS=0.406; [0,1.44) 629/183=0.291; prompt+bkg=all & fake+had=bkg consistent; op_all>>ss_all OS-dominated).
**Physics note**: non-unity background mirror (SS/OS≈0.3-0.4) confirmed EXPECTED (correlated HF in low-UE signal MC, ≈k), not a defect.

**Status**: APPROVED at iteration 1
**Summary**: Reco-seeded provenance macro on pythia fullsim validated. Background exists+survives quality cuts; mirror NOT symmetric in signal MC (correlated HF, SS/OS≈k≈0.3) → need HIJING overlay (UE) for the symmetric-combinatoric mirror.
