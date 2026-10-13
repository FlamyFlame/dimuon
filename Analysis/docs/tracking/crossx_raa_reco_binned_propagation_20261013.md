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

Pb+Pb (per pair, centrality c of the event, pair in reco cell j with pT bin p):

    ε_PbPb = E_Run2(p, c) × G_j
    E_Run2(p, c) = Σ_{T pairs, pT bin p} w·ε_Run2(pT_1, q·η_1; c)·ε_Run2(pT_2, q·η_2; c) / Σ_{T pairs, pT bin p} w
    G_j = ε^{reco,SF}_j / ε^{nogap}_p ,  ε^{nogap}_p = Σ_{MC pairs ∈ R^{nogap}_p} w / N[T_p]

- ε_Run2 = the existing placeholder single-muon fit (`EvaluateSingleMuonRecoEffPlaceholder`, WP-matched key, its
  centrality intervals 0-10,…,50-60,60-80), evaluated at the TRUTH muon (pT, q·η) of the pp24 fullsim T pairs —
  the η^pair-AVERAGED Run-2 pair efficiency, so only its pT and centrality dependence survive.
- R^{nogap}_p: like R but WITHOUT the single-muon q·η gap windows (reco mass window, reco pair pT, reco |η^pair| < 2.2
  kept = inside the axis), no SF, reco pT in bin p: the efficiency the Run-2 placeholder describes (no gap cuts, MC).
- G_j therefore carries, relative to the Run-2 level: the gap-cut acceptance + the η^pair (and m) shape + the SF
  product. E·G reduces to ε^{reco,SF}_j if E_Run2 = ε^{nogap}_p (pp24 ↔ Run-2 level difference is the only thing kept
  from Run 2, together with its centrality dependence).
- Trigger: unchanged placeholder union ε_1 + ε_2 − ε_1ε_2 of the Pb+Pb data single-muon turn-ons, no correlation term.

R_AA = (OS−SS Pb+Pb, per-year ΣN/ΣL combined) / (T_AA × OS−SS pp), unchanged procedure (RAA_plotting.cxx).

### 3. Step-by-step method
(a) **SF cells.** pp24 fullsim `_mc_trig_full` pair tree (the SAME truth-seeded pair sample as the nominal `_full`
tree, plus SF branches; `PythiaFullSimExtras.c:368-375`) → N[T_j], N[R_j], M (unweighted) — REQUIRED to equal the
migration file cell by cell (set identity) — plus SF-weighted N^{SF}[R_j], M^{SF}, TR^{SF}, and N[R^{nogap}] per pT bin.
(b) **pp product** ε^{reco,SF}_j (3D, m-integrated 2D, inclusive) with the existing A/B/D error formula on the
SF-weighted pieces; new definition marker.
(c) **Pb+Pb product** E_Run2(p, c) (TH2D) and G_j (TH3D) in one file with its own marker.
(d) RDF: pp `w_reco` = 1/ε^{reco,SF}(reco pT, η^pair, m); Pb+Pb `w_reco` = 1/(E·G) in every place it is defined.
(e) Rerun on frozen OLD NTuples (`~/usatlasdata/dimuon_data/frozen_old_grl_20261013/`) with OLD lumis
(pp 400.412 pb⁻¹; Pb+Pb 23/24/25/26 = 1.17576/0.85112/2.59933/2.62316 nb⁻¹) from a run-only git worktree whose
extra TEMP commit (input dir + lumis) is NEVER merged.

### 4. Negative constraints
1. No binning change; all axes read from the canonical sources. No ΔR anywhere in the reco correction.
2. Do not cap ε^{reco,SF} at 1; do not exclude feed-in.
3. Do not touch the trigger corrections, the NTuples, the GRL/lumi constants on master, or the other session's files.
4. Pb+Pb: do not take the pT or centrality dependence from pp24; do not take the η^pair dependence from Run 2.
5. R_AA: OS−SS, no template fit.
6. Investigation (item 3) is read-only: no code change from it.

## Design Decisions
- D1 (agent, user-delegated "find a way"): Pb+Pb scaling = §2 factorization E_Run2(p,c)·G_j, G normalised to the
  pp24 NO-GAP η-integrated ε of the same pT bin (a normalised-shape G would drop the gap-acceptance LOSS, ~18 %, from
  the level). G is 3D (incl. m) so pp and Pb+Pb share the same (η^pair, m) structure in R_AA.
- D2 (agent): Pb+Pb SFs = the pp24 Run-3 Tight SFs carried by G (no Pb+Pb SF source exists for the Run-2 placeholder).
- D3 (agent): OLD GRL realised by frozen copies of the pre-rerun NTuples + a run-only worktree TEMP commit with the
  old lumis; master keeps the new constants.

## Implementation Plan
| # | Step | § | Review | Status |
|---|---|---|---|---|
| 1 | Freeze old NTuples | §3e | — | DONE 00:45 |
| 2 | SF-cells filler + pp product builder + Pb+Pb product builder | §3a–c | /review-analysis-code (§2–4) | |
| 3 | RDF pp + Pb+Pb w_reco | §3d | /review-analysis-code | |
| 4 | Worktree run: pp RDF + plots, Pb+Pb RDF ×4 + plots, MC–data, R_AA | §3e | /review-plot | |
| 5 | Investigation (subagent, read-only) | item 3 | /review-investigation | running |
| 6 | Docs, INDEX, commits | — | — | |

## Progress Log
- **Step 1 (00:45).** Frozen copies (cp -p, mtimes kept): pp24 2mu4 nominal (10-06) + single_mu4 res_cut_v2 (09-09)
  from `pp_2024/pre_rerun_20261013/`; Pb+Pb 23/24/25/26 `_mindR_0_02` + `_res_cut_v2` (09-09/09-17) →
  `~/usatlasdata/dimuon_data/frozen_old_grl_20261013/`. Old lumis from git: pp 400.412 (0e3942e^), 2026 2.62316 (a9beddc^).
  These NTuples predate the pp LB exclusion, the 2026 bad-run/LB exclusion, D8 and the new cut 2.
  Known: old pp24 NTuple lacks part13 (run 488534 recovery) while 400.412 includes that run (other doc's finding).

## Results & Observations

## Remaining Work
Steps 2–6.

## Latest Stage
00:55 — Step 2: write SF-cells filler (reco_binned/), extend builder (SF variant), Pb+Pb product builder.
