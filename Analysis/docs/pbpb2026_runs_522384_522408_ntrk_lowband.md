# Pb+Pb 2026 runs 522384 and 522408: reduced N_trk^HItight at fixed FCal ΣE_T — evidence record

**Status:** evidence record for the DQ / detector experts. **No run is excluded and no cut has been
changed by the analysis on this basis** — that decision belongs to the experts
(memory `feedback_run_exclusion_needs_experts`). Produced 2026-09-17 by
`Analysis/plotting_codes/event_selection/plot_pbpb_ntrk_lowband_runs.cxx` from the merged `hi2026`
skim (`~/usatlasdata/dimuon_data/pbpb_2026/data_pbpb26_part{1..7}.root`, 270 087 106 events,
`HLT_mu4_L1MU3V` events only, no event-selection cut applied).

## Summary

| | runs 522384 + 522408 | all other 33 GRL runs |
|---|---|---|
| Events (HLT_mu4, skim) | 6 502 522 + 11 072 507 = 17 575 029 (**6.51 %** of 2026) | 252 512 077 |
| Luminosity (`Prescale Corrected`, HLT_mu4 table) | 63.224 + 107.570 = **170.79 µb⁻¹ = 6.51 %** of 2623.16 µb⁻¹ | 2452.37 µb⁻¹ |
| LB range in the skim (GRL-good LBs) | 133–356 (224 good LBs) and 132–688 (557 good LBs) | — |
| N_trk^HItight / FCal ΣE_T^{A+C} (FCal > 1 TeV), peak | **≈ 470 TeV⁻¹** | ≈ 590 TeV⁻¹ |
| Fraction of events below the cut-5 lower band edge at FCal > 1.5 TeV | 40.1 % / 38.6 % | 0.19 % |
| Fraction of events cut 5 would reject (any FCal) | 28.7 % / 27.6 % | 0.23 % |
| LB dependence | **none** — every one of the 223 / 551 LBs sits at 20–60 % below-edge (mean 40.0 / 37.9 %); no LB block is normal | — |

**What is seen.** In these two runs the number of HItight tracks per unit FCal transverse energy is
~20 % lower than in every other 2026 run (and than in 2023/24/25), for the WHOLE run, in every
lumiblock, and the N_trk–FCal correlation is not a scaled copy of the normal band but bends over
(saturates) at high FCal: N_trk^HItight reaches ~2000 at FCal ΣE_T ≈ 5 TeV instead of ~3000. The
shortfall grows with multiplicity, which points at a tracking-side condition (detector element(s)
missing from the HItight track selection / occupancy-dependent efficiency loss), not at the
calorimeter: FCal ΣE_T itself is normal (the runs populate the same FCal range as the others). The
two runs are consecutive Pb+Pb runs in the same period; the GRL in use is `physics_HI2026_50ns_noIBL.xml`.

**Two neighbouring runs checked and found normal:** 522336 and 522355 have the nominal ratio peak
(≈ 590 TeV⁻¹; below-edge fractions 0.30 % / 0.23 %, same as the bulk). Run 522336 shows a small
shoulder at ratio ≈ 500–530 TeV⁻¹ (visible in the 1D ratio panel), a mild version of the effect in a
subset of events — worth a look by the experts but not comparable to 522384/522408.

## Why it matters for the analysis

Event-selection **cut 5** (N_trk^HItight within µ ± 5σ of the linear N_trk–FCal band, derived per
year from the year's own data) does **not** handle these runs cleanly:
- the band was fitted on the two-population sample, so its lower edge was pulled down for all 2026
  runs (looser than 2025's);
- the lower edge bisects the low band: events of 522384/522408 are **rejected above ≈ 3.5 TeV, partly
  rejected between ≈ 1.5 and 3.5 TeV, and kept below** — a centrality-dependent removal (28 % of the
  two runs' events, 1.56 % of all 2026 events; `centrality_ratio_after_before_cuts_pbpb_2026.png`
  shows an 8–10 % event loss over centrality 0–8 %) while the runs' full luminosity stays in the
  cross-section denominator.
- A direct cut on N_trk^HItight alone cannot isolate them either: their N_trk range overlaps the
  normal band's mid-central range. The variable that separates the population is
  N_trk^HItight / FCal ΣE_T (≈ 470 vs 590 TeV⁻¹, fully resolved at FCal > 1 TeV).

Possible treatments, **for the experts / analysis lead to decide**: (a) exclude both runs from data
AND luminosity (−6.5 % of the 2026 luminosity; `PbPbBadRuns(26)`, `PbPbSampledLumi.h`,
`make_crossx_factors_pbpb_2026`, lumi README in one change, then re-derive the 2026 cuts and rerun);
(b) keep them and add an N_trk/FCal-based selection or a per-run-block cut-5 band, after confirming
the muon reconstruction (combined muons need ID tracks) and the pair yield per luminosity in these
runs are nominal; (c) accept and record the central-class bias.

## Figures and table (all in `~/usatlasdata/dimuon_data/plots/single_b_analysis/event_selection/pbpb_2026/ntrk_lowband_runs/`)

- `ntrk_lowband_ntrk_vs_fcal_pbpb_2026.png` — N_trk^HItight vs FCal ΣE_T^{A+C}, same axes as the
  cut-5 figures ({120,−0.5,5.5} × {100,0,3500}), one panel per run (522336, 522355, 522384, 522408)
  + all other runs, with the cut-5 band (red) overlaid.
- `ntrk_lowband_ratio_1d_pbpb_2026.png` — N_trk^HItight / FCal ΣE_T (FCal > 1 TeV), same panels.
- `ntrk_lowband_ntrk_1d_pbpb_2026.png` — N_trk^HItight distribution, same panels.
- `ntrk_lowband_runs_pbpb_2026.txt` — per-run table for all 35 runs (events, LB range, events at
  FCal > 1.5 TeV, below-edge count and fraction, LB range of below-edge events, cut-5 rejection
  fraction) and the per-LB listing for the four runs above.
- Context: `../event_sel_cut5_nTrk_FCal_band_standalone_pbpb_2026.png`,
  `../event_sel_cut5_nTrk_FCal_band_fail_5panel_pbpb_2026.png`,
  `../centrality_ratio_after_before_cuts_pbpb_2026.png`.

Regenerate: `cd Analysis/plotting_codes/event_selection && root -l -b -q 'plot_pbpb_ntrk_lowband_runs.cxx+(26, {522336, 522355, 522384, 522408})'`.
