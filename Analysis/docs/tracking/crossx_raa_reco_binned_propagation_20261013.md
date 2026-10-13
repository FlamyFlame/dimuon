# Propagate the 3D bin-by-bin reco efficiency (+ muon SFs) to pp/Pb+Pb crossx, MC–data comparison and R_AA

**Created:** 2026-10-13 00:55 · **Mode:** Implementation · **Session:** crossx · **Deadline (user):** plots by 11:00 2026-10-13

## Objective
User request 2026-10-13: propagate the reco efficiency with bin-by-bin corrections in (pair pT, η^pair, m_μμ)
(`pair_reco_eff_bin_by_bin_reco_binned.md`) to all final results: pp crossx, Pb+Pb crossx, pp MC–data comparison,
R_AA. Keep (1) the Run-2 Pb+Pb reco-eff placeholder and (2) the Pb+Pb placeholder trigger correction (no correlation
efficiency). R_AA stays OS−SS (template fit not ready). Use the OLD GRL + lumitable (the NTuple rerun for the new
GRL/ZDC cut, `ntuple_rerun_pp24_pbpb2526_20261013.md`, cannot finish in time). First handle:
1. Pb+Pb placeholder lacks the gap-cut acceptance → wrong η^pair dependence: keep the Run-2 centrality + pair-pT
   dependence, take the η^pair dependence from pp24 (gap acceptance).
2. Pair-level (MC) reco efficiencies must be corrected by the product of single-muon data/MC SFs; the Pb+Pb
   placeholder is MC-derived (colleague's Run-2 MC TF1 fits) → corrected too.
3. Read-only investigation: why R_AA at high pair pT is ≫ 1 (mistakes in the crossx process?). NO code changes from it.

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = (a) SF-weighted reco-binned pp24 ε product + scaled Pb+Pb placeholder product built and verified;
  (b) pp24 crossx RDF refilled with ε^{reco,SF}(pT,η^pair,m) and old-GRL inputs/lumi (400.412 pb⁻¹), crossx plots +
  correction-stage/sanity plots regenerated; (c) Pb+Pb 23/24/25/26 crossx RDF refilled with the scaled placeholder
  (old NTuples, old lumi incl. 2026 = 2.62316 nb⁻¹), combined crossx plots regenerated; (d) pp MC–data comparison
  plots regenerated; (e) R_AA (OS−SS) plots regenerated; (f) investigation findings reported (no code change);
  (g) reviews (/review-analysis-code, /review-plot), docs, INDEX, commits. Plots by 11:00.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

### 1. Motivation
The measured yield is binned in reco kinematics; the result is quoted for the truth population. The pp24 fullsim
3D bin-by-bin correction ε^reco_j (pT, η^pair, m) is the decided reco correction (user 2026-10-09/11). MC efficiencies
must be brought to data with the single-muon data/MC efficiency SFs (pair = product of the two legs). The Pb+Pb
analysis has no overlay-based ε yet; its Run-2 placeholder (single-muon ε(pT, q·η, centrality) product) knows nothing
of the Run-3 detector-gap fiducial cuts (acceptance, folded into ε_reco since 2026-09-17), so its η^pair dependence
is wrong for our signal region.

### 2. Top-level equations
pp (per reco 3D cell j = (pT, η^pair, m) on the canonical axes `ParamsSet::pair_pt_coarse_bins` ×
`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` × `ParamsSet::SignalMinvRecoEffcyBins()`):

    dσ_pp ∝ (1/L_pp) Σ_{data OS pairs} w_trig / ε^{reco,SF}_j
    ε^{reco,SF}_j = Σ_{MC pairs ∈ R_j} w·SF_1·SF_2 / N[T_j]

T_j, R_j exactly as `pair_reco_eff_bin_by_bin_reco_binned.md` §2 (R incl. out-of-T feed-in, no cap at 1).
SF_i = `m_i.eff_sf_<wp>` (skim MuonEfficiencyScaleFactors, WP-matched; Tight nominal). w_trig = the existing pp24
HYBRID 2mu4 weight (unchanged).

SF_i below the SF tool's validity (pT < 5 GeV, no SF provided) = the MC-weighted mean filled SF of WP muons in the
nearest valid slice [5, 6) GeV at the muon's η (D5); SF = 1 there is the systematic variation.

Pb+Pb (per pair, centrality c of the event, pair in reco cell j with pT bin p) — D4, anchored to pp24:

    ε_PbPb = C_Run2(p, c) × ε^{reco,SF}_j
    C_Run2(p, c) = E_Run2(p, c) / E_Run2(p, 60–80 %)
    E_Run2(p, c) = Σ_{T pairs, pT bin p} w·ε_Run2(pT_1, q·η_1; c)·ε_Run2(pT_2, q·η_2; c) / Σ_{T pairs, pT bin p} w

- ε_Run2 = the existing placeholder single-muon fit (`EvaluateSingleMuonRecoEffPlaceholder`, WP-matched key, its
  centrality intervals 0-10,…,50-60,60-80), evaluated at the TRUTH muon (pT, q·η) of the pp24 fullsim T pairs —
  the η^pair-AVERAGED Run-2 pair efficiency.
- Run 2 supplies ONLY the centrality dependence relative to its most peripheral interval (per pT bin); level, pT,
  η^pair (gap acceptance), m and the SFs are the pp24 ε^{reco,SF}_j. 60–80 % Pb+Pb ε = pp24 ε exactly; the
  per-centrality R_AA efficiency ratio is 1/C_Run2.
- Trigger: unchanged placeholder union ε_1 + ε_2 − ε_1ε_2 of the Pb+Pb data single-muon turn-ons, no correlation term.

R_AA = (OS−SS Pb+Pb, per-year ΣN/ΣL combined) / (T_AA × OS−SS pp), unchanged procedure (RAA_plotting.cxx).

### 3. Step-by-step method
(a) **SF cells.** pp24 fullsim `_mc_trig_full` pair tree (the SAME truth-seeded pair sample as the nominal `_full`
tree, plus SF branches; `PythiaFullSimExtras.c:368-375`) → N[T_j], N[R_j], M (unweighted) — REQUIRED to equal the
migration file cell by cell (set identity) — plus SF-weighted N^{SF}[R_j], M^{SF}, TR^{SF}, and N[R^{nogap}] per pT bin.
(b) **pp product** ε^{reco,SF}_j (3D, m-integrated 2D, inclusive); value = N^{SF}[R]/N[T] (exact); error = the
MC-only A/B/D relative error carried over (A/B/D assumes one weight per pair, SF-on-R breaks it; SF ≈ const per cell);
SF statistical/systematic uncertainties not yet propagated (Remaining Work); new definition marker.
(c) **Pb+Pb product** C_Run2(p, c) (TH2D, + absolute E_Run2 for reference) and ε^{reco,SF}_j (TH3D) in one file with its own marker.
(d) RDF: pp `w_reco` = 1/ε^{reco,SF}(reco pT, η^pair, m); Pb+Pb `w_reco` = 1/(E·G) in every place it is defined.
(e) Rerun on frozen OLD NTuples (`~/usatlasdata/dimuon_data/frozen_old_grl_20261013/`) with OLD lumis
(pp 400.412 pb⁻¹; Pb+Pb 23/24/25/26 = 1.17576/0.85112/2.59933/2.62316 nb⁻¹) from a run-only git worktree whose
extra TEMP commit (input dir + lumis) is NEVER merged.

### 4. Negative constraints
1. No binning change; all axes read from the canonical sources. No ΔR anywhere in the reco correction.
2. Do not cap ε^{reco,SF} at 1; do not exclude feed-in.
3. Do not touch the trigger corrections, the NTuples, the GRL/lumi constants on master, or the other session's files.
4. Pb+Pb: do not take the centrality dependence from pp24; do not take the η^pair dependence or the level from Run 2 (D4).
5. R_AA: OS−SS, no template fit.
6. Investigation (item 3) is read-only: no code change from it.

## Design Decisions
- D1 (agent, user-delegated "find a way"): Pb+Pb scaling = §2 factorization E_Run2(p,c)·G_j, G normalised to the
  pp24 NO-GAP η-integrated ε of the same pT bin (a normalised-shape G would drop the gap-acceptance LOSS, ~18 %, from
  the level). G is 3D (incl. m) so pp and Pb+Pb share the same (η^pair, m) structure in R_AA.
- D2 (agent): Pb+Pb SFs = the pp24 Run-3 Tight SFs carried by G (no Pb+Pb SF source exists for the Run-2 placeholder).
- D1 SUPERSEDED by D4 (old: E_Run2 absolute × G = ε^{reco,SF}/ε^{nogap}, keeping the Run-2 level — review it.1 W2:
  even 60–80 % carried ε^{nogap}/E = 1.12–1.24 into R_AA, 0–10 % 1.36–1.59).
- D4 (USER 2026-10-13, "anchor to pp24"): ε_PbPb = [E_Run2(p,c)/E_Run2(p,60–80 %)] × ε^{reco,SF}_j (§2). Physics: a
  Run-2 MC fit vs Run-3 MC level difference is detector/era modelling, not Pb+Pb physics, and must not enter R_AA.
- D5 (USER 2026-10-13, "nearest valid SF"): muons with pT < 5 GeV (SF tool out of validity, skim stores 0 → NTP 1)
  take the [5,6) GeV mean SF at their η (0.1 lookup granularity); SF = 1 = systematic (review it.1 W1: −3.1 % on ε at
  9–12.79 GeV, −1.0 % 12.79–18.18, −0.5 % 18.18–25.85, −1.7 % inclusive; cancels in R_AA).
- D3 (agent): OLD GRL realised by frozen copies of the pre-rerun NTuples + a run-only worktree TEMP commit with the
  old lumis; master keeps the new constants.

## Implementation Plan
| # | Step | § | Review | Status |
|---|---|---|---|---|
| 1 | Freeze old NTuples | §3e | — | DONE 00:45 |
| 2 | SF-cells filler + pp product builder + Pb+Pb product builder | §3a–c | /review-analysis-code (§2–4) | DONE (review running) |
| 3 | RDF pp + Pb+Pb w_reco | §3d | /review-analysis-code | DONE a2121af (review running) |
| 4 | Worktree run: pp RDF + plots, Pb+Pb RDF ×4 + plots, MC–data, R_AA | §3e | /review-plot | |
| 5 | Investigation (subagent, read-only) | item 3 | — (read-only, reported to user) | DONE 01:15 |
| 6 | Docs, INDEX, commits | — | — | |

## Progress Log
- **Step 1 (00:45).** Frozen copies (cp -p, mtimes kept): pp24 2mu4 nominal (10-06) + single_mu4 res_cut_v2 (09-09)
  from `pp_2024/pre_rerun_20261013/`; Pb+Pb 23/24/25/26 `_mindR_0_02` + `_res_cut_v2` (09-09/09-17) →
  `~/usatlasdata/dimuon_data/frozen_old_grl_20261013/`. Old lumis from git: pp 400.412 (0e3942e^), 2026 2.62316 (a9beddc^).
  These NTuples predate the pp LB exclusion, the 2026 bad-run/LB exclusion, D8 and the new cut 2.
  Known: old pp24 NTuple lacks part13 (run 488534 recovery) while 400.412 includes that run (other doc's finding).

- **Step 2 (01:00–01:25).** `plotting_codes/reco_effcy/reco_binned/fill_pp24_fullsim_pair_reco_sf_cells.C` →
  `~/usatlasdata/pythia_fullsim_full_sample/reco_eff/pair_reco_sf_cells_pp24_full.root`: unweighted T/R/TR/M from the
  `_mc_trig_full` tree IDENTICAL to the migration file (max rel diff 0, both WPs) ⇒ same pair sample.
  `build_pp24_fullsim_pair_reco_eff_reco_binned.C(true)` → `pair_reco_eff_reco_binned_sf_pp24_full.root` (marker ValueSF;
  error = MC-only relative error, since A/B/D assumes one weight per pair). `build_pbpb_reco_eff_scaled_placeholder.C`
  (run after `.L RDFBasedHistFillingPbPb.cxx+`, Run-2 lookup re-exported from the data class) →
  `pbpb_reco_eff_scaled_placeholder.root`; evaluator `Utilities/PbPbRecoEffScaledPlaceholderEvaluator.h`.
- **Step 3 (01:30).** RDF: pp `w_reco` = 1/ε^{reco,SF}(reco pT, η^pair, m) (`PairRecoEffRecoBinnedEvaluator`, SF marker);
  Pb+Pb: the 3 per-muon placeholder blocks → one `PbPbRecoEffPair` (E·G); census in `HistPostProcessPbPb`. Commit a2121af.
  Worktree `/usatlas/u/yuhanguo/workarea/dimuon_codes_wt_oldgrl` (branch `crossx_old_grl_run_20261013`) + TEMP 215e921
  (frozen input dirs; pp 400.412; 2026 2.62316 in PbPbSampledLumi.h + PbPbBaseClass.h) — never merged. Previous outputs
  (10-07) backed up to `frozen_old_grl_20261013/prev_outputs_20261007/`. Driver: session scratchpad `run_oldgrl.sh`.

## Results & Observations
### R1 — muon SFs (pp24 fullsim, Tight R pairs)
Nominal (D5) ⟨SF_1 SF_2⟩ = 0.9148 (Tight), 0.9636 (Medium); SF = 1 below 5 GeV variation: 0.9308 / 0.9703. Per-muon
Tight SF ≈ 0.95 (5–6 GeV) → 0.977 (100 GeV); η: ≈0.976 barrel, ≈0.92–0.93 at 1.4<|η|<1.8, 0.99 at η≈0. SF UNAVAILABLE
(tool out of validity; skim `TrigRates.cxx:1322-1341` zeroes it on a non-Ok code, NTP stores 1) for EVERY muon with
pT < 5 GeV and for none above (guard in the filler, both WPs) → D5 nearest-valid [5,6) GeV SF at η. D5 vs SF = 1 on
ε^{reco,SF} (Tight): −3.10 % (9–12.79), −0.94 % (12.79–18.18), −0.53 % (18.18–25.85), −0.33 % … −0.07 % above;
−1.72 % inclusive (review it.2 verified). Cancels in R_AA (both systems carry ε^{reco,SF}).
### R2 — ε products (Tight; Medium in brackets), D4/D5 nominal
pp inclusive ε^{reco,SF} = 0.5415 (0.6202) vs MC-only 0.5919 (0.6436); per pT bin m-integrated ≈ 0.49 → 0.64.
Pb+Pb: C_Run2 = E_Run2/E_Run2(60–80 %) (Tight), per pT bin, centrality 0–10/10–20/20–30/30–40/40–50/50–60/60–80 %:
| pT [GeV] | 0–10 | 10–20 | 20–30 | 30–40 | 40–50 | 50–60 | 60–80 |
|---|---|---|---|---|---|---|---|
| 9–12.79 | 0.829 | 0.932 | 0.972 | 0.995 | 1.008 | 1.027 | 1 |
| 12.79–18.18 | 0.794 | 0.902 | 0.939 | 0.968 | 0.975 | 0.990 | 1 |
| 18.18–25.85 | 0.787 | 0.897 | 0.933 | 0.962 | 0.969 | 0.983 | 1 |
| 25.85–36.74 | 0.784 | 0.895 | 0.930 | 0.960 | 0.966 | 0.980 | 1 |
| 36.74–52.23 | 0.783 | 0.896 | 0.930 | 0.959 | 0.964 | 0.978 | 1 |
| 52.23–74.24 | 0.783 | 0.897 | 0.930 | 0.959 | 0.964 | 0.976 | 1 |
| 74.24–105.53 | 0.783 | 0.898 | 0.930 | 0.959 | 0.963 | 0.975 | 1 |
| 105.53–150 | 0.783 | 0.902 | 0.932 | 0.960 | 0.962 | 0.974 | 1 |
(Medium 0–10 %: 0.822–0.832.) C_Run2 > 1 only at 9–12.79 GeV, 40–60 % (≤ 2.7 %, inherited non-monotonicity of the
Run-2 fits; placeholder systematic). Absolute E_Run2 (reference, h_E_run2_abs): Tight 0–10 % 0.466–0.504, 60–80 %
0.563–0.644, i.e. 0.80–0.89 × pp24 ε^{nogap} even peripheral — the level D4 removes from R_AA.
Superseded D1-era products (E absolute × ε^{reco,SF}/ε^{nogap}, SF=1 below 5 GeV) were used only for the first
(01:00) RDF/plot pass, overwritten by the D4/D5 rerun.
### R3 — R_AA high-pT investigation (subagent, read-only; scratch merged below)
Conclusion (on the 10-07 inputs, before this propagation): NO crossx-process mistake inflates high-pT R_AA.
(B) cells > 1 are statistics: 74–150 GeV integrated R_AA < 1 in every centrality (0.45±0.09 … 0.59±0.21); all 34 cells
> 52 GeV fit flat per centrality with χ² 33.4/28; pp has 112/32/12/4 OS pairs in bins 13–16; the 5–10 % excess is one
year (2025: 6 vs 29 pairs above 74 GeV in 0–5 vs 5–10 %, ~3.9σ local, no artefact found). Symmetric errors understate the
upward tail. (A, pushes R_AA DOWN) pp 2mu4 correction above 74.24 GeV (single-value MC pair eff × SFs) over-corrects pp
by ~1.5× (1.50±0.25 vs single-mu4 shape; 0.95/1.02 below) — consistent with R13 of mu_pt45_gap125_pairpt9_adoption.md.
(A, known) mixed fiducial (pp ε incl. gap acceptance ×0.817, Pb+Pb placeholder without) biased R_AA low ×0.82 —
addressed by THIS doc's Pb+Pb scaling. (C) rising R_AA with pT is expected (KB open_hf_production.md). Ruled out: weight
blow-ups, out-of-range fits/floors, SS subtraction (except pp 105–126 GeV: 2 SS / 12 OS), binning/selection mismatch,
year combination, centrality mis-assignment, photoproduction cut.

#### Appendix — merged subagent scratch doc `_sub_raa_highpt.md` (verbatim)
# Sub-investigation: R_AA high pair-pT values (read-only)
Created 2026-10-13. READ-ONLY. Scratch dir = this dir.

##### Objective
User: R_AA vs pair pT crazy (>1) at high pT. Besides stats, are there mistakes in the crossx process?
Trace each factor of high-pT bins (pp + PbPb), quantify, classify A(bug)/B(stats)/C(physics).

##### Plan
1. Read docs (raa tracking, pp24 crossx rerun, pp trig-eff high-pT jump, mupt45 adoption) + RAA_plotting.cxx.
2. Locate current R_AA PNGs + input ROOT files; dump R_AA numbers bin by bin (macro).
3. Trace factors: OS, SS, counts vs weighted, mean weights per bin pp vs PbPb; normalization; binning.
4. Per-pair weight distributions from frozen pair trees if needed.
5. Rank findings.

##### Findings (append-only)

###### Step 1-2 (docs + numbers) DONE
- Current R_AA PNG: ~/usatlasdata/dimuon_data/plots/single_b_analysis/RAA/raa_pair_pt_pbpb23_24_25_26_combined_pp24_2mu4.png (2026-10-06 11:32), produced by Analysis/RAA_plotting.cxx case 6 mode 1.
  Inputs refilled AFTER it: pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root (2026-10-07 00:58), pbpb_20{23..26}/..._single_mu4_no_trg_plots_nominal.root (2026-10-07 01:01-01:02). PNG mildly stale; recomputed numbers within few % of PNG (e.g. 5-10% bin15 1.87 vs 1.8 on PNG; 30-50% bin16 2.94 vs 2.96).
- Recompute script: raa_dump.C -> raa_dump.out (R_AA per ctr x 16 pT bins, pp OS/SS counts, pp stage factors).
- PbPb factors: pbpb_fac.C -> pbpb_fac.out (per-ctr OS/SS counts summed over years, <w_trig>,<w_reco> from stage hists _corr_raw/_corr_unfolded_trig/_corr_unfolded_reco).
- KEY NUMBERS (pair pT bins of ParamsSet::pT_bins_150, 16 log bins 9-150):
  pp <w_trig> = 1.74 (15-22 GeV) -> 2.20 (44-52) -> 2.69/2.84 (52-74) -> 3.77/4.07/4.59/4.11 (74-150; hybrid MC single-value branch)
  pp <w_reco> flat 1.46-1.6 (above 15 GeV).
  PbPb (0-80 all) <w_trig> flat 1.04-1.06; <w_reco> flat 1.71-1.85 (Run-2 TF1 placeholder).
  => <w_PbPb>/<w_pp> = 0.66 (15-18) -> 0.56 (44-52) -> 0.42 (62-74) -> 0.32/0.33/0.27/0.28 (74-150): the CORRECTIONS pull high-pT R_AA DOWN, not up.
  Raw OS COUNT ratio N_PbPb(0-80, 4 yrs)/N_pp: 0.43 (15-18) 0.73 (37-44) 0.80 (44-52) 0.85 (52-62) 1.14 (62-74) 1.22 (74-88) 1.84 (88-105) 2.17 (105-126) 3.25 (126-150): raw ratio rises x7.5 => the R_AA rise lives in RAW counts, i.e. pp 2mu4 raw suppression at small dR (high pT) only partly compensated by w_trig.
- R_AA >1 cells: 5-10% b14 1.13+-0.38, b15 1.87+-0.94, b16 1.12+-1.0; 10-20% b15 1.20+-0.59; 30-50% b16 2.94+-2.05; 50-80% b12 1.10+-0.25, b14 1.12+-0.69; 20-30% b16 1.12+-0.99. None >1.5 sigma above 1 (5-10% b15: 0.9 sigma above 1).
- pp counts high bins: OS 112/32/12/4, SS 1/1/2/0 (bins 13-16). PbPb OS (0-80) 137/59/26/13.

###### Step 3-4 (factor tracing, data-driven checks) DONE
Scripts/outputs in this dir: trig2mu4.C/.out, selchk.C/.out, raa_int.C/.out, hiyr.C/.out, list25.C + list25_ctr5_10.out, ctrdist.C/.out, neff.C/.out.
- VALIDATION: cut reproduction on frozen pp 2mu4 tree (ParamsSet Signal/Fiducial/PairFiducial exprs + pair_pass_tight) reproduces h2d_counts per bin EXACTLY (67397,129030,...,112,32,12,4). PbPb frozen pT>52 = 1209 vs hist 1198 (centrality acceptance).
- F1 (data-driven 2mu4 closure) pp24 single_mu4 frozen tree (res_cut_v2, apply mu pT>=4.5): P(2mu4|mu4) = 0.583 (15-22) -> 0.448/0.435 (37-52) -> 0.419 (52-62) -> 0.315 (62-74) -> 0.25 (74-88). Barrel-only 0.457 -> 0.29 -> 0.21. Both-muons-mu4-matched pairs: 0.86 -> 0.35 => close-pair (same-RoI) 2mu4 loss is REAL.
  S = N_2mu4*<w_trig>/N_mu4 normalized to [15.25,36.74): 0.947+-0.028 [36.7,52.2), 1.015+-0.066 [52.2,74.2), 1.50+-0.25 [74.2,150) (45 mu4 pairs).
  => pp hybrid trigger correction closes vs single-mu4 shape to ~5% up to 74 GeV; above 74 GeV (single-value MC pair-eff branch) the pp is OVER-corrected by ~1.5x (2 sigma) -> R_AA there biased LOW, not high. Caveat: assumes eps_OR(single mu4) flat in pT (PbPb bare-union weight 1.04-1.06 flat).
- F2 weights: pp N_eff/N = 0.75-0.89 in bins 5-16 (no few-pair domination). pp barrel <w_trig> 2.22 (18-22) -> 3.2 (44-52) -> 4.1-4.5 (52-74) -> 6.6/6.5/6.2/5.8 (74-150); endcap 1.41 -> 1.5 -> 1.65 -> 1.9-2.4. PbPb w flat; reco placeholder pT clamped to [4,19] (RDFBasedHistFillingData.cxx:806-809), TF1 eval fine, floor 0.05 not reached (<w_reco> 1.7-2.3).
- F3 integrated R_AA (raa_int.out): 74.2-150 GeV: 0-5% 0.446+-0.088, 5-10% 0.875+-0.147, 10-20% 0.595+-0.097, 20-30% 0.483+-0.106, 30-50% 0.673+-0.127, 50-80% 0.594+-0.207. 52-150: 0.45-0.77. ALL <1. chi2 of bins 11-16 cells vs per-ctr 52-150 integrated value = 33.4/28 (p~0.2); no cell >1.5 sigma ABOVE; only pull <-1.5 is 5-10% bin11 (-2.9 sigma, low).
- F4 5-10% high-pT excess traced: per year/ctr pT>74 counts (hiyr.out): 2025 0-5%=6 vs 5-10%=29 (15-37 GeV counts 4794 vs 4838) -> local ~3.9 sigma split in ONE year; other years 9/8, 5/7, 17/13. 29 pairs = 29 distinct events across 20 runs (list25_ctr5_10.out); fine-centrality distributions smooth (ctrdist.out) => statistical, no structural artefact found. Now-dropped 2026 runs/LBs contribute 3 pairs >74 GeV.
- F5 selection consistency (selchk.out): both sides min mu pT 4.500, |eta|<2.4, dP/p one-sided <=0.12, |d0|<2, tight, same signal-cut expressions (RDFBasedHistFillingPP.cxx:612-616, PbPb.cxx:992-997), same res-cut mode (mindR_0_02 v1), same axes (16 pT_bins_150 x 48 eta). pp |z0| up to 293 mm vs PbPb <10 mm (pp pileup; not pT specific, not investigated).
- F6 mixed fiducial (known, pair_reco_eff_gap_acceptance.md R3): pp eps_reco includes gap acceptance (x0.817), PbPb placeholder does not -> R_AA biased LOW by ~0.82 globally; fixing it RAISES R_AA ~x1.22. pT-dependence of A_gap alone not measured here; pp <w_reco>/PbPb <w_reco> ratio varies 1.04-1.21 bin to bin (smooth-ish).
- F7 lumi transient: hist files filled 2026-10-07 with L26=2.62316 (crossx factor) but PbPbSampledLumi.h now 2.35549 (D8, 7e46aff 2026-10-13) -> 2026 events weighted 0.898 in the combine; net combined effect <~1%, global; resolves at the in-progress refill.
- F8 presentation: R_AA PNG (10-06) older than inputs (10-07); empty cells (50-80% bins 15-16) drawn as 0+-0; linear x on log bins; Gaussian ratio errors for N_pp=4..12 / N_PbPb=2..9 understate the upward tail (B).
- Physics (KB open_hf_production.md:105-115): D-meson R_AA rises to ~0.7 at 100 GeV; b/c/light degeneracy above ~25 GeV => rising R_AA 0.3->0.5-0.8 with pair pT is EXPECTED (C).

##### CONCLUSION (ranked)
1. (B) R_AA>1 cells are statistical: integrated 74-150 R_AA <=0.88 in all centralities; chi2 33.4/28 vs flat; largest outlier is 2025 0-5/5-10 split (3.9 sigma local in one year, no artefact found).
2. (A, opposite sign) pp 2mu4 hybrid correction above 74.24 GeV over-corrects pp by 1.50+-0.25 vs data single-mu4 closure -> high-pT R_AA biased LOW by ~1/1.5. Up to 74 GeV closes to ~5%.
3. (A, known) mixed fiducial: R_AA globally LOW by ~0.82.
4. (A, minor, transient) 2026 lumi mismatch hist vs header (<1% net); stale PNG; 0+-0 empty cells; Gaussian ratio errors.
5. (C) rising R_AA with pT expected.
Ruled out: weight blow-up (pp N_eff/N>=0.75; PbPb flat), TF1 out-of-range/floors (PbPb reco clamp 19 GeV, trig union flat 1.05), SS over-subtraction (SS/OS<=3% above 37 GeV except pp bin15 2/12), binning mismatch, selection mismatch, year-combination formula, centrality mis-assignment (smooth int-centrality distributions).


## Remaining Work
- SF statistical + systematic uncertainties (incl. D5 SF=1 variation) not propagated.
- Steps 4–6 (rerun with D4/D5, plots, reviews, commit).

## Latest Stage
02:00 — code review it.2 → fixing marker/doc WARNINGs; rebuild products, rerun RDF + plots (D4/D5) from the worktree; then /review-plot on result plots, commit, report.
