# Task 02 — pp24 2mu4 trigger efficiency (Pipelines 2+3)

**Roadmap:** `analysis_roadmap_2026_06.md` §Q3a / chain [3]
**Depends on:** task 01 (new pp24 ntuples with passmu4noL1 fix).
**Reviewer:** `/review-analysis-code`, then `/review-plot`.
**Physics procedure:** `docs/tracking/mu4_trig_effcy_implementation.md`
— read §Top-level equation (PP 2mu4) and §3d before any work. Key
constraints: P(pair passes 2mu4 | dR) = ε₁^nc · ε₂^nc · ε_dR^{2mu4}(dR);
**no tag/probe, no cross-term, numerator uses `pass2mu4`** (hardware
pair decision, not m1.passmu4 && m2.passmu4).

## Objective

Measure the pp24 no-correlation single-muon mu4 efficiency (Pipeline 2)
and the 2mu4 pair-level dR correction (Pipeline 3, §3d), completing
mu4-doc Remaining Work item "PP24 Pipeline 2+3 when re-skim completes".

## Inputs

- pp24 trig-eff ntuples from task 01 (`_mindR_0_02_res_cut_v2` flavor —
  verify the pp NTP run produces these; trig-eff processing needs
  trigger_mode allowing efficiency derivation, mirroring PbPb P2 setup)
- Code already implemented: `RDFBasedHistFillingPP.cxx`
  (FillTrigEffcyHistsInvWeightedbySingleMuonEffcies, OpenEffcyPtFitFile,
  MakeAndWriteDRTrigEffGraphs), `RDFBasedHistFillingData.{h,cxx}`,
  `SingleMuEffcyPtTurnOnFitter.cxx` (PP mode = `erf_plus_log`).

## Steps

1. Read the mu4 tracking doc Physics Procedure fully (INVARIANT).
2. P2: RDF hist fill, cycle=generic, fine q·η (useCoarseQEtaBin=false),
   trigger-efficiency derivation enabled for pp24.
3. pT turn-on fitting with `erf_plus_log` (PP-specific; run interpreted,
   `.L` without `+` — CRTP rootcling issue).
4. Validate fits (TF1 count, TH2D `_divided` fallbacks) as in PbPb
   pipeline stages 6–7.
5. P3: cycle=2 (inv_weight). Note D9 bias caveat: on a 2mu4-selected
   sample the inv-weighted ratio numerator condition equals the
   selection condition — check with the user / Physics Procedure how
   the §3d measurement is interpreted before quoting numbers.
6. Plots via `TrigEffPlotterPP` and `plot_dR_trig_corr.C` (PP paths).
7. Update mu4 tracking doc (mark 2f / Remaining Work), roadmap ledger.

## Verification

- ε^nc fits: clear turn-on, plateau values physically sensible
  (compare PbPb plateau 0.6–0.9).
- mu4_mu4noL1-related fits now constrained (passmu4noL1 fix in skim) if
  that trigger is in scope (see roadmap Q2.4).
