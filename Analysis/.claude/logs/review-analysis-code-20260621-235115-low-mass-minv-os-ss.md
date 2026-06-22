# Analysis Code Review Log
**Task**: Add efficiency-corrected, crossx-normalized OS and SS minv histograms over 0–4 GeV with the dimuon mass-window cut REMOVED, to the PP and PbPb crossx RDF hist-filling (data inputs D_OS, D_SS for the low-mass dimuon template fit). Tracking doc: Analysis/docs/tracking/low_mass_dimuon_template_fit.md (Physics Procedure §2,§3a,§4).
**Log file**: review-analysis-code-20260621-235115-low-mass-minv-os-ss.md
**Started**: 2026-06-21 23:51:15
**Status**: APPROVED at iteration 1
**Summary**: OS+SS minv 0–4 GeV no-mass-cut crossx histos added to PP + PbPb RDF (dσ-weighted, identical OS/SS selection/binning/weight). Reviewer PASS iter 1, all numbers MATCH. Surfaced OS-only resonance veto (ntuple-level) as a data characteristic.
**Iterations completed**: 1
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 0 (one INFO: removing minv<2.9 sends high-mass pairs to overflow; harmless, fit uses in-range bins).
**Numerical verification**: ALL MATCH — PP OS=6209, SS=990.6 (50×[0,4]); OS≥SS continuum 0 violations; PbPb ctr0_5 OS=3.672e4/SS=2.655e4, ctr50_80 OS=6736; all 6 PbPb centralities × OS+SS = 12 histos present; OS resonance-veto zeros at [0,1.06],[2.9,3.3],[3.55,3.8] confirmed (explained, not a bug).

## Executor work (pre-review)
**Files modified:**
- `RDFBasedHistFilling/RDFBasedHistFillingPP.cxx` — new block at end of `FillHistogramsCrossx`: `signal_cuts_no_minv` (= signal_cuts minus the two minv conditions), `attach_crossx_weight` lambda (reuses generic w_trig/w_reco if present else defines inline, then crossx_weight·w_reco·w_trig), fills `h1d_crossx_minv_0_4_{op,ss}_dsigma` (50 bins, 0–4 GeV, weight `crossx_weight_trig_corr`).
- `RDFBasedHistFilling/RDFBasedHistFillingPbPb.cxx` — new block at end of `FillHistogramsCrossx`: same `signal_cuts_no_minv`, `attach_dsigma_weight` lambda (OR-logic pair trig ε₁+ε₂−ε₁ε₂, per-ctr reco, `weight_for_dsigma_trig_corr`), per-centrality fills `h1d_crossx_minv_0_4_{op,ss}_dsigma_<ctr>` (50 bins, 0–4 GeV).

**Compilation:** PP `.so` and PbPb `.so` both ACLiC-clean in separate ROOT sessions (PP_COMPILE_OK, PBPB_COMPILE_OK; only pre-existing ParamsSet sign-compare warnings).

**Runtime (verification):** pp24 + pbpb23 crossx reran (read existing ntuples). pp24: "# 1D histograms: 2" (OS+SS); pbpb23: "# 1D histograms: 12" (6 ctr × OS+SS). Both "completed successfully".

**Numbers (sanity):**
- PP OS integral=6209 (1.68M entries), 50 bins 0–4 GeV; PP SS integral=990.6 (552k entries). OS≥SS in ALL continuum bins 1.05–2.0 (0 violations). OS/SS ≈ 6.3.
- PbPb ctr0_5: OS=3.672e4, SS=2.655e4 (SS/OS≈72%, large central combinatoric — expected); ctr50_80 OS=6736.

**KEY PHYSICS FINDING (data characteristic, NOT a code bug):** the PP OS spectrum has content EXACTLY 0 in [0,1.06], [2.9,3.3], [3.55,3.8] GeV while SS is smooth everywhere. Root cause: `DimuonAlgCoreT::ResonanceTaggingImpl` (NTupleProcessingCode/DimuonAlgCoreT.c:128) tags+removes resonance muons **only for op_sign (OS)** — SS is never tagged. Veto windows `pms.minv_cuts` (ParamsSet.h:404-407): [0,1.06] (ρ/ω/φ + below), [2.9,3.3] (J/ψ 3.097), [3.55,3.8] (ψ(2S) 3.686), [9.08,10.5] (Υ). So the OS muon-pair tree is the resonance-vetoed `_no_res` ntuple; SS is unvetoed. The new histos faithfully reflect this (OS and SS use the identical RDF filter — the asymmetry is purely upstream in the trees). Implication for the fit: OS continuum = S+C+G as modeled; but OS−SS is only valid in OS-populated bins (mask the veto windows in BOTH OS and SS), and to SEE the J/ψ in OS one would need a no-veto ntuple pass (trigger_effcy_calc skips resonance cuts).
