# Plot Review Log
**Task**: r-tag sanity plot set for the HIJING-overlay Pb+Pb-2024 test sample after the contact's
re-run (tracking doc hijing_overlay_pbpb24_test_sample_skim.md, R6-R9). Code:
Analysis/plotting_codes/overlay_rtag_checks/plot_r17864_event_level.cxx
(modes fcal|vtx|ntrk|l1te|centr|bcid|trainpos; centr/bcid/trainpos are NEW).
Outputs: ~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/plots/r17864_rtag_sanity/.
**Log file**: review-plot-20260923-163759-rtag-sanity-bunch-train.md
**Started**: 2026-09-23T16:37:59-04:00
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Execution record (iteration 1 input)
Files modified:
- Analysis/plotting_codes/overlay_rtag_checks/plot_r17864_event_level.cxx
  - header updated (the pbpb24 NTUP is now the RE-RUN production, merge tid52629267; the first
    production survives as ...NTUP.bak_20260923.root)
  - NEW_LABEL "2024 conditions (r17864)" -> "2024 conditions (r17864, re-run)"
  - NEW mode centr  : centrality population table -> centrality_population_*.txt (no plot)
  - NEW mode bcid   : filled-bunch pattern (first 400 crossings, 1 bin per bcid) + <FCal> vs bcid
                      profile (12 bcid per point) -> bcid_structure_*.png
  - NEW mode trainpos: <FCal> and <L1TE> vs position in the bunch train (56 positions)
                      -> train_position_*.png
Plots produced (all 10 000-event samples, both sides the SAME pT-hat 125-300 slice):
  fcal_sum_et_r17864_vs_r17618_pTH125_300.png       (re-run)
  l1te_r17864_vs_r17618_pTH125_300.png              (re-run)
  trk_multiplicity_r17864_vs_r17618_pTH125_300.png  (re-run)
  vertex_z_reco_vs_truth_r17864.png                 (re-run)
  bcid_structure_r17864_vs_r17618_pTH125_300.png    (NEW)
  train_position_r17864_vs_r17618_pTH125_300.png    (NEW)
  centrality_population_r17864_vs_r17618_pTH125_300.txt (NEW, text)
  README.txt  (documents the mixed vintage: the 5 single-muon efficiency PNGs in the same
               directory are from 2026-09-16 and were deliberately NOT regenerated)
KB: no entry covers LAr digitisation / bunch structure (checked .claude/kb/index.md).

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 8 (4 WARNING, 4 INFO)
**Details**:
1. [PHYSICS / DERIVATION ARTEFACT] WARNING — the 2 isolated bunches are assigned tpos = 0 and merged
   into the train-head point, conflating the only two healthy crossings with the (unhealthy) train
   heads. Train heads 5.104+-0.054 FCal / 10.421+-0.188 L1TE; isolated 4.564+-0.127 / 16.671+-0.532.
2. [LABEL SEMANTICS / P1c] WARNING — on-canvas "18 trains of 56 bunches, 2 bcid spacing" misstates the
   pattern: 882 gaps of 2 + 108 of 4, i.e. 18 trains of 7x8 bunches, the 4-bcid gaps after positions
   7/15/23/31/39/47. The 8-bunch sub-burst IS the observed period; "period 8 = 16 bcid" is wrong.
3. [NUMERICAL CLAIM] WARNING — "the 39 good bcids are EXACTLY 18x2 heads + 2 isolated" is false
   (18x2+2 = 38); the 39th is bcid 1993 (position 25), a chance pass in the recovery phase.
4. [BOOKKEEPING] WARNING — README.txt does not list train_position_*.png.
5. [CODE DOC] INFO — macro header claims a progressive degradation; the data are periodic.
6. [PROVENANCE DOC] INFO — header's branch list stale for the new modes (bcid/L1TE/trk_numqual/centrality).
7. [REDUNDANCY] INFO — trk_numqual[0]==[4] and [3]==[7] element-wise; panel titles contradict that.
8. [R3] INFO — trk_multiplicity has no ratio panel (pre-existing, passed 2026-09-17).
**Numerical verification**: every reported number MATCHED except issue 3 (the "exactly 39" claim).
The crux — modulation locked to position-in-train, not absolute bcid — was independently CONFIRMED
(flat in bcid % 16; mod-8 argmax = phase 1 in all 18 trains over 14 distinct start % 16 values).

**Amendments made (all verified independently by the executor before applying):**
- train_positions(): single-bunch groups now flagged kIsolated = -2 instead of position 0; the profile
  covers in-train bunches only and the isolated bunches are drawn as their own dashed reference line
  with the value in the legend.
- Legend header corrected to "18 trains of 7x8 bunches + 2 isolated bunches"; bcid-mode comment now
  states the measured gap spectrum and the 18-bcid sub-burst period.
- Legend moved to a reserved strip above the frame (SetTopMargin 0.24): the previous in-frame legend
  was CLIPPED at the right edge, which the executor caught by looking at the rendered PNG.
- trainpos printf replaced by the mod-8 fold + the isolated-bunch values (both samples).
- Macro header: progressive -> periodic; provenance branch list extended with a bcid clause.
- trk_numqual panel titles now state that [0]==[4] and [3]==[7] (stored collection already has pT>400 MeV).
- README.txt lists both new figures, tagged [R8]/[R9].
- Tracking doc R9 rewritten (isolated-vs-head split, 7x8 pattern, 38-of-39, period = sub-burst = 18 bcid);
  INDEX scope line corrected (it had said "<FCal> vs bcid is FLAT" — wrong).

## Iteration 2
**Reviewer verdict**: FAIL
**Issues found**: 7 (1 WARNING, 6 INFO)
**Details**:
1. [DOC-CONSISTENCY] WARNING — `Remaining Work` bullet 1 still cited "R8's flat <FCal> vs bcid" as report
   evidence, the very claim R9 overturned, and R9 was missing from the evidence lists in Remaining Work,
   Latest Stage and the INDEX scope line. A draft built from that bullet would have told the production
   contact the deficit is flat in bcid and dropped the one reproducible signature.
2. [LAYOUT] INFO — fcal right-panel legend entries ran into the right frame edge / axis ticks.
3. [DOC] INFO — R9's "+15 % / -36 %" had an unstated baseline.
4. [DOC] INFO — "gap of more than 5 empty crossings" vs the code's bcid DIFFERENCE > 5.
5. [CODE-COMMENT] INFO — header doc pointer "(R1, R4, R6-R8)"; "~13k tracks" vs measured max 11 952.
6. [PRESENTATION] INFO — trk_multiplicity still draws two pairs of identical panels (accepted as labelled).
7. [DOC] INFO — pbpb23 <L1TE> span quoted 16.1-16.7, actual 15.86-16.66.
**Numerical verification**: EVERY quantity MATCH — filled pattern, gap spectrum {2:882, 4:108}, train-length
spectrum {56:18, 1:2}, isolated bcids 826/2611, 4-bcid gaps after positions 7/15/23/31/39/47, sub-burst
period 18 bcid (450 ns), head-vs-isolated split, all mod-8 phases with errors, threshold scan, the
"38 of 39" position census {0:18, 1:18, -2:2, 25:1}, <FCal> vs bcid%16 featureless (0.562-0.777 TeV),
mod-8 argmax = phase 1 in all 18 trains over 14 distinct start%16 values, trk_numqual [0]==[4]/[3]==[7].
The reviewer also re-ran `bcid`, `trainpos` and `ntrk` from the working-tree source and got
BYTE-IDENTICAL PNGs.

**Amendments made:**
- Doc: Remaining Work bullet 1 rewritten around R9; Latest Stage -> "report R1/R4/R7/R8/R9"; INDEX scope
  line -> "R1/R4/R7/R8/R9 ... (R9 is the actionable part)". (WARNING 1)
- Doc: "+15 % vs the 2023 reference" / "(-36 % vs it)"; "a bcid DIFFERENCE of more than 5, i.e. more than
  4 empty crossings"; pbpb23 <L1TE> "15.9-16.7 (all-event mean 16.27)". (INFO 3, 4, 7)
- Macro: header pointer "(R1, R4, R6-R9)"; "~12k tracks"; fcal right-panel legend entries reduced to one
  per sample (the header already defines solid = A / dashed = C), with a code comment recording why.
  (INFO 5, 2)
- INFO 6 deliberately NOT changed: accepted in iteration 1 as addressed by honest panel titles.
- Regenerated `fcal` and `ntrk`; inspected the rendered fcal PNG - no clipping, no overlap.

## Iteration 3
**Reviewer verdict**: PASS
**Issues found**: 1 INFO (macro comment still said "gap of more than 5 empty crossings" where the code
splits on a bcid DIFFERENCE > 5) — fixed immediately after the verdict.
**Numerical verification**: pbpb24 median FCal 0.4071 TeV unbinned / 0.4103 binned, 36.17 % negative;
pbpb23 4.3966 TeV, 0.00 % negative — all MATCH the values drawn in the legend. Rest carried over from
iteration 2 (all MATCH). The reviewer also grepped the whole doc + INDEX row for any surviving
flat-in-bcid claim: none, and confirmed `git diff` of the macro contains only the documented changes.

**Status**: APPROVED at iteration 3
**Summary**: The r-tag sanity plot set for the re-run pbpb24 slice is approved. Two review rounds
corrected a real derivation artefact (isolated bunches merged into the train-head bin, which had hidden
the sample's only healthy crossings), a misstated filled-bunch pattern (18 trains of 7x8, not 56 at
uniform 2-bcid spacing; sub-burst period 18 bcid, not 16), an arithmetic error ("exactly 39" vs 38 of 39),
a clipped legend, and a superseded claim ("<FCal> vs bcid is FLAT") that had propagated into the doc's
report-evidence list and the INDEX scope line.
