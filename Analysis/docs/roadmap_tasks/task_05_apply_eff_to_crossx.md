# Task 05 — Apply trigger + reco efficiency corrections to crossx

**Roadmap:** `analysis_roadmap_2026_06.md` §Q3c / chain [5],[6]
**Depends on:** task 04 (dummy ε_reco). **Reviewer:** `/review-analysis-code`
(this is the physics-critical step), then `/review-plot`.

> **Status update 2026-06-14:** the **trigger** part of this task is DONE —
> the no-correlation per-pair trigger weight is already applied in the crossx
> pipeline (PbPb `RDFBasedHistFillingPbPb.cxx:926`, pp24
> `RDFBasedHistFillingPP.cxx:373`). This task now covers **only the reco-eff
> weight**.

## Objective

Add the per-pair **reconstruction-efficiency** weight to the crossx
measurement (trigger weight already applied). Target total weight:

```
w⁻¹ = ε_trig^pair · ε_reco(pair pT, pair η, dR)
```

- **ε_reco** is a single **pair** reconstruction efficiency binned in
  **(pair pT, pair η, dR)** — derived from Pythia fullsim (pp) / Pythia
  fullsim HIJING overlay (PbPb, per centrality). It is **NOT** a product of
  two single-muon efficiencies (the two muons are correlated; see
  `docs/tracking/hijing_overlay_reco_effcy_investigation.md` and
  `RDFBasedHistFillingPythiaFullsim.cxx` `reco_effcy_var3Ds`).
- **ε_trig^pair** (already applied): PbPb single-mu4 inclusion–exclusion
  `ε₁^nc + ε₂^nc − ε₁^nc·ε₂^nc`; pp24 2mu4 `ε₁^nc · ε₂^nc`. εᵢ^nc = fitted
  ε^{no-corr}(pTᵢ, q·ηᵢ) per (centrality, μ±) from
  `pbpb_20YY/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root`.
- dR **trigger** correction ε_dR: **dummy ≡ 1** until measured on unbiased
  MC (roadmap Q4); keep a hook for it.

## Design notes

- The lookup machinery already exists for Pipeline 3
  (`OpenEffcyPtFitFile`, `EvaluateSingleMuonEffcyPtFitted`, TH2D
  fallback for q·η gaps) in `RDFBasedHistFillingData.cxx` /
  `RDFBasedHistFillingPbPb.cxx` / `RDFBasedHistFillingPP.cxx`. Reuse it
  in the crossx (cycle=generic, `mu4_nominal_pbpb_NO_trig_calc=true`)
  path — but verify in the Physics Procedure that reuse is justified
  (it is the same ε^nc; different application context).
- Add a config flag (e.g. `apply_pair_eff_weights`) defaulting OFF so
  existing outputs remain reproducible; pipeline scripts turn it ON.
- Weight column: multiply existing per-pair weight (FCal-scaled RAA
  weight path in `FillMuonPairExtra` / `CalculateWeightForRAA`) by w.
  Watch double-counting: confirm no other efficiency-like factor is
  already in the weight.
- Guard 1/ε blow-up near threshold (ε small at pT≈4 GeV): clamp ε at a
  floor (e.g. 0.05) and count clamped pairs.
- Naming: derive from physics terms (eff_pair_trig, invw_pair_eff), per
  tracking-doc naming rules.

## Steps

1. Write plan into mu4 implementation tracking doc (new step) — this is
   trigger-efficiency *application*, governed by its Physics Procedure.
2. Implement in RDF crossx filling for PbPb and PP; ACLiC both.
3. `/review-analysis-code` with the top-level-equation sections quoted.
4. Rerun crossx hist filling for PbPb 23/24/25 + pp24
   (`SKIP_CONDOR=1` pipelines; ntuples unchanged).
5. Rerun combined crossx plotting; produce corrected-vs-uncorrected
   overlay plots for sanity.
6. Update docs (pbpb_pipelines.md, data_analysis.md) + roadmap ledger.

## Verification

- ⟨w⟩ distribution sane (Run 2 analog: smooth, w ≳ 1, no long tail);
  plot ⟨w⟩ vs pT of lower-pT muon as in Run 2 note Fig. check_pair_effs.
- Corrected/uncorrected ratio ≈ 1/⟨w⟩, smooth in pT and centrality.
- With all ε set to 1 (flag off), outputs bit-identical to current.
