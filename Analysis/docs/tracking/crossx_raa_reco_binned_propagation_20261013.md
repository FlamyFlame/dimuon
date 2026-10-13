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
(d) RDF: pp `w_reco` = 1/ε^{reco,SF}(reco pT, η^pair, m); Pb+Pb `w_reco` = 1/(C_Run2·ε^{reco,SF}_j) in every place it is defined.
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
- D6 (USER 2026-10-13, after /review-plot it.1 C1): the MC–data SIGNAL family compares the trigger+reco(SF)-CORRECTED
  data (= truth signal region WITHOUT gap cuts, |η^pair| < 2.2) with MC truth in the SAME population: new filters
  `_<cat>_pass_signal_truth_nogap` (Pythia fullsim single_b/op/ss; POWHEG fullsim single_b) = truth mass window +
  truth pair pT + |truth η^pair| < 2.2, no single-muon gap windows; plotter keys switched. `_gapcut` kept (uncorrected
  comparison). Corrects the 09-17 premise in `pair_reco_eff_gap_acceptance.md` (d).
- D7 (USER 2026-10-13): the R_AA-vs-η^pair W-shape (×1.6, centrality-independent, from the trigger corrections, not
  the reco eff) → read-only investigation, no trigger change (subagent; results → R4).
- D7 resolved (USER 2026-10-13): R4 W-shape = document only (placeholder limitation; fix with the Run-3 overlay ε and
  a pp-hybrid per-|η^pair| calibration/systematic later).
- D8 (USER 2026-10-13, binning): R_AA uses the crossx binning — pair pT on `ParamsSet::pT_bins_150` (16 bins) inside the
  9 coarse η^pair cells `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` (RAA_plotting mode 4, one PNG per
  centrality: `raa_pair_pt_in_eta_subplots_ctr<lo>_<hi>_…png`); R_AA vs η^pair (mode 2) on the coarse cells (was 48 fine).
- D9 (USER 2026-10-13): MC–data signal η^pair plot on the coarse cells (the bin-by-bin ε is valid per coarse cell only).
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

- **Step 4 (01:50–02:10).** Rerun with D4/D5 (`run_oldgrl_all2.log`): RDF all 5 OK; pp ∫ = 4483.59 pb (×1.0758 vs
  10-07's 4167.62; entries 516649 identical); Pb+Pb ×1.093–1.097 per year (entries identical). Code review APPROVED it.3
  (`.claude/logs/review-analysis-code-20261013-005603-reco-binned-sf-propagation.md`).
- **Step 4b (02:10–).** /review-plot it.1 FAIL: (1) R_AA-vs-η W-shape [trigger, → D7 investigation]; (2) MC–data
  fiducial mismatch [→ D6, MC refill: Pythia fullsim FULL + POWHEG fullsim pp17, backups
  `pythia_fullsim_full_sample/backup/…_pre_nogap`, `powheg_full_sample/histograms_powheg_fullsim_pp17.root.bak_20261013_pre_nogap`];
  (3) Pb+Pb pT-in-η header covers points → `SingleBCrossxPlotterBase::subplot_ceil_pad` (Pb+Pb 30, pp default).
  (A first refill attempt was killed: the no-gap filter lacked |truth η^pair| < 2.2; fixed before any output was written.)

- **Step 4c (12:50–13:25).** MC refills (Pythia fullsim FULL: 467 old keys bit-identical + 20 no-gap keys; op no-gap ∫ =
  6.73888 = N[T]; gap/no-gap 0.817. POWHEG fullsim pp17: fixed a latent duplicate `powheg_n_gen_sample` Define (base class
  per-sample N_gen since 022b9bd 09-11); refill also refreshed 402 stale pre-09-11-normalised keys (×~2, MC–data key
  unchanged 2204.65 pb; no-gap 2699.76 pb)). Final data RDF with the D5 SF marker 13:00–13:03 (pp ∫ 4483.6 pb). Plots
  13:05–13:20 incl. D8/D9. /review-plot it.2 running. NOTE: wall clock — the 11:00 deadline passed while a user question
  was pending (≈02:15 → ≈12:35); the D4/D5 crossx/R_AA plots had been at the canonical paths since ≈02:00.

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
| 9–12.79 | 0.829 | 0.932 | 0.972 | 0.995 | 1.007 | 1.027 | 1 |
| 12.79–18.18 | 0.794 | 0.902 | 0.939 | 0.967 | 0.975 | 0.990 | 1 |
| 18.18–25.85 | 0.787 | 0.897 | 0.933 | 0.962 | 0.969 | 0.983 | 1 |
| 25.85–36.74 | 0.784 | 0.895 | 0.930 | 0.960 | 0.966 | 0.980 | 1 |
| 36.74–52.23 | 0.783 | 0.896 | 0.929 | 0.959 | 0.964 | 0.978 | 1 |
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


### R4 — R_AA vs η^pair W-shape (D7 investigation, read-only; scratch merged below)
Normalised factorisation (|η^pair| 0–0.5/0.5–1/1–1.5/1.5–2/2–2.2, pT 9–150): reco eff cancels (Pb+Pb/pp mean reco weight
flat in η: 1.238–1.245 at 0–5 %, 0.991–1.000 at 50–80 %). SAME-trigger test (Pb+Pb mu4 vs pp24 single_mu4, both union-
corrected): 0.947/1.150/1.054/0.866/0.907, identical at 0–5 % and 50–80 %, every year, every pT slice → (D, DOMINANT,
≈×1.33 of the ×1.45 contrast) a Pb+Pb-vs-pp η-dependent reco/selection efficiency difference that D4 cannot represent
(D4 takes the η shape from pp24 by construction; the Run-2 placeholder had a W-like η shape, 60–80 %:
0.75/0.85/0.81/0.78/0.80, that partly compensated before). Corrected pp data match Pythia reco η shape (1.042/1.000/
0.976/0.987/0.975); corrected Pb+Pb do not (0.987/1.150/1.028/0.855/0.884). Ruled out: WP, momentum scale/Δp/p, Pb+Pb
trigger turn-ons, physics (centrality-independent). (B/D, ≈7 %) pp 2mu4 hybrid under-corrects barrel vs endcap:
P(2mu4|mu4) obs/pred 0.936/0.920 (|η^pair|<1) vs 1.015/1.006 (1–2), 0.974 (2–2.2). (D, minor, opposite sign) Pb+Pb union
lacks correlation: both-legs obs/pred 1.04–1.08 low pT, ≈1.4 barrel >26 GeV. Fix options (USER): η-dependent Pb+Pb ε from
the Run-3 overlay, or a Pb+Pb/pp η-shape stopgap factor + systematic; pp hybrid per-|η^pair| calibration or systematic.

#### Appendix — merged subagent scratch doc `_sub_raa_eta.md` (verbatim)
# _sub_raa_eta — R_AA vs eta^pair W-shape root cause (READ-ONLY investigation)
Mandate: read-only; outputs only in this scratch dir. Done = ranked root-cause findings w/ evidence + classification (A/B/C/D) + fix proposal.

##### Plan
S1 read tracking docs + code (PairTrigEffCrossxEvaluator.h, PP AddPairEfficiencyWeightColumns, PbPb union weight, RAA_plotting case 6)
S2 inspect RDF hist files (*_corr_* histos) per eta^pair: raw, trig-corr, reco-corr pp vs PbPb
S3 pp24 data closure: single_mu4 res_cut_v2 ntuple -> fraction of sig-region OS pairs firing 2mu4 vs eta^pair vs hybrid prediction
S4 PbPb union weight eta shape; peripheral vs pp shapes
S5 history: old plots / 09-10 version

##### Findings (append-only)

###### S1/S2 (done) — docs + hist decomposition
- Read crossx_raa_reco_binned_propagation (D4: eps_PbPb = C_Run2(p,c) x eps_pp^{reco,SF}; trig PbPb = bare union), pp24_trig_eff_hybrid_application (§2 hybrid), mc_trigger_efficiency §1-2 (eps^nc = data T&P at dR>0.8; pp product eps1 eps2 eps_dR^2mu4).
- eta.C (copy of plotreview1/eta.C) -> eta.out. CAUTION: its "trig" column = crossx/no_trig_corr, and no_trig_corr = 1/L-normalised RAW counts (recopp=1/400.4) -> that column is the TOTAL (trig x reco) correction, mislabelled.
- R_AA(norm. to integral) new vs old (10-07 inputs) per eta^pair bin: OLD shape is the opposite in 1-1.5 vs 1.5-2 (Rold 0.65-0.85 at 1-1.5, 1.05-1.27 at 1.5-2; Rnew 0.9-1.24 / 0.72-0.9). The change = reco correction swap (D4, pp eps in both) -> old PbPb Run-2 placeholder lacked gap acceptance.
- Raw OS count ratio PbPb/pp (normalised): 0-5%: 0.62-0.70 at |eta|1.5-2, 1.34-1.51 at 0.5-1, 0.75-1.1 at |eta|<0.5 (dip 0.75 at -0.1..0), ~0.8-1.2 at 1-1.5, 0.83-1.05 at 2-2.2. 50-80% same shape (0.62-0.79 / 1.30-1.70).
- stages.C -> stages.out (stage hists _corr_raw/_corr_unfolded_trig/_corr_unfolded_reco_trig, OS only, pT-integrated, |eta^pair| groups 0/.5/1/1.5/2/2.2):
  0-5%: raw AA/pp 1.129/1.451/0.949/0.676/0.891; <w_trig>pp 2.192/2.340/1.729/1.495/1.815; <w_trig>AA 1.138/1.178/1.153/1.081/1.074; <w_reco> ratio AA/pp flat 1.238-1.245; R(norm) 0.966/1.205/1.041/0.802/0.865.
  50-80%: raw 1.131/1.422/0.907/0.701/0.977; <w_trig>AA 1.082/1.093/1.057/1.019/1.023; reco ratio 0.991-1.000; R(norm) 0.985/1.171/0.986/0.846/0.974.
  => reco cancels (flat). Raw barrel(0.5-1)/endcap(1.5-2) contrast = 2.07 (0-5) / 2.03 (50-80); trigger-weight ratio (w_pp/w_AA) contrast predicted = (2.340/1.178)/(1.495/1.081) = 1.44 (0-5), (2.340/1.093)/(1.495/1.019)=1.46 (50-80). Residual ~1.4-1.5 = the W.

###### S3 plan — pp data-driven closure (closure.C)
On frozen pp24 single_mu4 res_cut_v2 tree (sign2, signal cuts as RDFBasedHistFillingPP.cxx:612-616 + pair_pass_tight):
 per |eta^pair| group: N_mu4, N(pass2mu4), N(both legs passmu4), predicted <eps_hyb/eps_union>, <e1e2/union>.
 eps_hyb from Utilities/PairTrigEffCrossxEvaluator.h (pp_full Tight) with pp data T&P TF1s (reimplemented lookup, cap1 floor .01).
 Also Y_mu4 = Sum 1/union (pp mu4 sample, PbPb-style correction) vs Y_2mu4 = Sum w_hyb (2mu4 nominal tree) shape ratio.

###### S3 result — closure.C -> closure.out (pp24 frozen single_mu4 res_cut_v2 sample, sig cuts + tight, mu pT>=4.5)
pair pT [9,150), |eta^pair| groups 0/.5/1/1.5/2/2.2:
 f_obs=P(2mu4|mu4):          0.465 0.435 0.614 0.692 0.557
 f_pred=<eps_hyb/eps_union>: 0.496 0.473 0.605 0.689 0.571 -> obs/pred 0.936 0.920 1.015 1.006 0.974
 P(both legs passmu4) obs/pred(e1e2/union): 1.065 1.034 1.046 1.023 0.989 (25.85-74 GeV barrel 1.37-1.39: close-by positive correlation, as eps_dR models)
 <eps_dR> 0.885 0.886 0.958 0.958 0.898
 KEY: Y_2mu4(hybrid, 2mu4 tree)/Y_mu4(union w/ pp T&P, mu4 tree), normalised: 0.961 0.966 1.052 1.038 1.000 (same in pT slices 9-12.8, 12.8-25.9, 25.9-74: within 0.93-1.07)
 => pp hybrid vs an independent pp mu4-union correction agree to <=+-5% in eta^pair; barrel/endcap(1.5-2) 0.93. The hybrid under-corrects pp barrel by ~7% relative to 1.5-2 at most (f_obs/f_pred 0.92 vs 1.006) -> pushes R_AA barrel UP by ~7-8%: a minor piece of the W (W is ~x1.45). The pp hybrid is NOT the main cause.

###### S4 result — closure2.C -> closure2.out: SAME trigger (HLT_mu4) in both systems
PbPb (frozen _mindR_0_02, 4 yrs, tight, sig cuts) vs pp24 single_mu4 (frozen res_cut_v2), |eta^pair| groups 0/.5/1/1.5/2/2.2, pT 9-150:
 ctr 0-80: rawA/rawP(norm) 0.954 1.146 1.040 0.874 0.940 ; union-corrected A/P(norm) 0.947 1.150 1.054 0.866 0.907
 <w_union> PbPb 1.109 1.136 1.097 1.044 1.040 ; pp 1.086 1.100 1.053 1.024 1.047  (both ~flat, differ <=4%)
 Same pattern in every pT slice (9-12.8 / 12.8-25.9 / 25.9-150) and in 0-5% AND 50-80% (50-80% 12.8-25.9: 0.915 1.131 1.047 0.892 0.996).
 P(both legs mu4)|mu4 obs/pred (e1e2/union): PbPb 1.08/1.06/1.04/1.00/0.94 ~ pp 1.07/1.03/1.05/1.02/0.99 -> per-leg T&P self-consistent in both systems.
 => ~x1.32 of the barrel(0.5-1)/endcap(1.5-2) contrast is in the RAW single-mu4 yields of PbPb vs pp and is NOT a trigger-correction effect (union weights only 1.02-1.14, cannot move 30%). The remaining ~x1.1 of the x1.45 W is the 2mu4-vs-mu4 trigger-correction difference in pp (S3: 0.96/1.04).
 Next: test selection pieces whose PbPb-vs-pp eta efficiency is NOT modelled (D4 copies pp reco eta shape): Tight WP fraction, Medium yields.
- wp.C -> wp.out: Tight fraction of sig-region OS pairs per |eta^pair| group: pp 0.824/0.866/0.921/0.983/0.991, PbPb50-80 0.832/0.875/0.927/0.974/0.989, PbPb0-10 0.831/0.868/0.914/0.968/0.978 -> WP not the cause (<1.5% diff). Raw ratio with NO WP: PbPb50-80/pp 0.927/1.127/1.023/0.907/0.999 (same W). Ruled out: Tight WP.
###### S5 plan — locate the raw deficit in single-muon eta / phi / year (mueta.C)
- mueta.C -> mueta.out: per-MUON eta (0.1 bins, both legs of sig-region tight OS pairs, mu4 samples) PbPb/pp_mu4 normalised ratio, per year 23/24/25/26:
  barrel |eta| 0.6-1.1: 1.11-1.29 (all years); 1.5-2.0: 0.74-0.95 (all years; 1.7-2.0 A side 0.74-0.87); eta~0-0.3: 0.87-1.0. Smooth broad modulation, SAME in every PbPb year (2023..2026) -> not a single-year detector condition. No sharp step at 1.3 (NSW edge).
  endcap phi (both mu |eta|>1.3): some sectors 0.71-0.78 (A: phi -2.55,-2.16; C: phi 1.77, 2.16) vs others 1.1-1.18 -> also sector structure.
  pp 2mu4/pp mu4 per-muon ratio: barrel 0.69-0.94, endcap 1.17-1.32 (expected: 2mu4 requires both legs -> barrel L1 ineff squared).
###### S6 plan — momentum-scale / dp/p differences (mom.C): <pt/trk_pt>, dP_overP distribution per |eta| pp vs PbPb
- mom.C -> mom.out: <pt/trk_pt> per |eta| identical pp vs PbPb23/25 (0.994-1.000, diffs <0.4%); dP/p distributions similar -> momentum scale / dp/p cut NOT the cause.
- Run-2 PbPb placeholder single-muon reco eff (run2_reco_eff_placeholder.root, tight, pT 6/10): 60-80%: |qeta|<0.5 0.75, 0.5-1 0.83-0.86, 1-1.5 0.81, 1.5-2 0.77-0.80, 2-2.4 0.80-0.82; 0-10%: 0.65 / 0.75-0.78 / 0.75 / 0.64-0.67 / 0.69-0.72.
  -> Run-2 PbPb MC eps has the SAME W pattern (pair ~eps^2: 0.5-1 vs 1.5-2 ~1.18 (60-80), ~1.3 (0-10); vs |eta|<0.5 ~1.25). This is why the 09-10 (Run-2-placeholder-reco) R_AA vs eta was flat apart from eta~0, and D4 (pp24-anchored eta shape, 2026-10-13) removed it.
- Exact factorisation of R(norm): R = [PbPb_mu4 union-corr / pp_mu4 union-corr] / [Y_2mu4/Y_mu4]
  = (0.947/0.961, 1.150/0.966, 1.054/1.052, 0.866/1.038, 0.907/1.000) = 0.985 1.190 1.002 0.834 0.907, reproducing stages R(norm) (0-5: 0.966 1.205 1.041 0.802 0.865; 50-80: 0.985 1.171 0.986 0.846 0.974).
- pp MC-data eta^pair plot (plots/mc_data_compr/signal/pair_eta_mc_data_compr.png, 10-13): MC/data sawtooth with steps at the coarse eta^pair reco-cell edges 0.5/1/1.5 (0.5-1: 1.3-1.6; 1-1.5: 0.8-0.9; 1.5-2: 1.15-1.35) -> within-cell eta shape of the reco-binned correction (cancels in R_AA, both sides use same eps).

###### S7 — which system deviates? data vs Pythia pp24 fullsim RECO-level eta^pair shape
MC: histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root h_truth_pair_eta_op_pass_tight_and_signal_truth_and_reco (reco-passing, tight, sig+gap), group fractions 0.2183 0.2846 0.1735 0.2552 0.0684.
pp24 data mu4 union-corrected (closure2, 9-150, N*<wu>): fractions 0.2276 0.2845 0.1694 0.2519 0.0667 -> pp/MC = 1.042 1.000 0.976 0.987 0.975  (pp data == Pythia Run-3 reco shape within +-4%)
PbPb data 0-80 union-corrected: fractions 0.2155 0.3274 0.1784 0.2182 0.0605 -> PbPb/MC = 0.987 1.150 1.028 0.855 0.884
=> the Pb+Pb data deviate from the Run-3 pp detector-level eta shape; pp does not. Centrality-independent (0-5 == 50-80), year-independent (23..26), WP-, trigger-, momentum-scale-independent -> a Pb+Pb-vs-pp reconstruction/selection efficiency eta dependence (HI reco / Pb+Pb-run detector configuration) that D4 cannot carry (D4 copies the pp24 eta shape into Pb+Pb by construction; C_Run2 depends on pT & centrality only).
Run-2 Pb+Pb MC placeholder carries a qualitatively similar W (single-mu 60-80%: 0.75/0.85/0.81/0.78/0.80 by |q eta| group) vs Run-2 pp (HF note Fig.31) barrel 0.95-0.97 / endcap 0.92-0.95 -> Run-2 PbPb/pp per-muon barrel/endcap contrast ~1.07 (pair ~1.14): explains roughly half of the observed x1.33 and is why the 09-10 R_AA(eta) was ~flat.
Overlay MC check (pythia_fullsim_hijing_overlay_test_sample[_pbpb23]) INCONCLUSIVE: truth denominators not defined like the pp file (pbpb23 truth shape incl. gap-like dips), stats 681-3017 pairs.

##### CONCLUSION (ranked)
1. [D, dominant, ~x1.33 of the x1.45 barrel(0.5-1)/endcap(1.5-2) contrast] Pb+Pb reco/selection-efficiency eta dependence NOT modelled: D4 ("anchor to pp24") gives Pb+Pb exactly the pp24 eta^pair shape. Data: same-trigger (HLT_mu4) Pb+Pb/pp yields after the union correction = 0.947/1.150/1.054/0.866/0.907; pp data agree with Pythia Run-3 reco shape to +-4%, Pb+Pb deviate +15%/-15%. Same in every year and centrality -> detector/reco, not physics. Fix: eta-dependent Pb+Pb efficiency from the Run-3 HIJING overlay (proper), or as interim a Pb+Pb/pp eta-shape factor (e.g. data-driven Pb+Pb J/psi or Z T&P vs pp; or Run-2 PbPb/Run-2 pp ratio as a stopgap with a systematic). User decision (it reverses part of D4).
2. [B/D, secondary, ~x1.07-1.08] pp 2mu4 hybrid under-corrects the barrel relative to the endcap: data closure P(2mu4|mu4) obs/pred = 0.936/0.920 (|eta^pair|<1) vs 1.015/1.006 (1-2), 0.974 (2-2.2); all pT slices (9-12.8: 0.95/0.93 vs 1.03/1.02; 25.9-74: 0.92/0.90 vs 0.99/0.97). Y_2mu4(hybrid)/Y_mu4(union) = 0.961/0.966/1.052/1.038/1.000. The MC eps_dR (|eta^pair|<1 group mean 0.885) / data eps^nc product is ~7% too optimistic in the barrel (barrel L1 close-by loss larger in data than MC). Fix: data-driven calibration of the hybrid with the per-|eta^pair|(,pT) ratio f_obs/f_pred from the single_mu4 sample, or a systematic of that size.
3. [D, minor, opposite sign] Pb+Pb union has no correlation term: P(both legs)|mu4 obs/pred 1.04-1.08 at low pT, ~1.4 above 26 GeV in the barrel (RoI sharing) -> union overestimates P(>=1) by ~2% (low pT) to ~14% (barrel, >26 GeV) -> pushes barrel R_AA DOWN; does not cause the W.
4. [C, expected] reco-eff cancels by construction (AA/pp <w_reco> ratio 1.238-1.245 at 0-5 %, 0.991-1.000 at 50-80 %, flat).
Ruled out: Tight WP (tight fractions equal to <1.5 %), momentum scale / dp/p (pt/trk_pt equal <0.4 %), Pb+Pb trigger turn-ons (union weights 1.02-1.14, self-consistent P(both|>=1)), single PbPb year, centrality-dependent effect.
History: 09-10 (3 yrs) flat because the Run-2 Pb+Pb placeholder (absolute, own eta shape ~W) partially compensated (and pre-hybrid pp trigger); 10-07 inputs: mixed fiducial (pp eps incl gap acceptance, PbPb Run-2 without) -> R low in gap cells (1-1.5, eta~0), high at 1.5-2 (Rold 0.65-0.85 vs 1.05-1.27); 10-13: D4 removes the Run-2 eta shape -> raw Pb+Pb/pp detector W exposed.


## Remaining Work
- SF statistical + systematic uncertainties (incl. D5 SF=1 variation) not propagated.
- Steps 4–6 (rerun with D4/D5, plots, reviews, commit).

## Latest Stage
13:25 — /review-plot it.2 running (all result plots incl. D8 R_AA pT-in-η and D9). Then: close review, commit docs, INDEX row, report.
