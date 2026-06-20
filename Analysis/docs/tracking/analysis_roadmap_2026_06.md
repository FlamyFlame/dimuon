# Analysis Roadmap — Run 3 Dimuon Internal Note (ANA-HION-2023-07-INT1)

**Status date: 2026-06-10.** Living document; update statuses as steps
complete. Companion checklist tracking doc:
`analysis_roadmap_buildout.md`. Per-task agent instructions:
`Analysis/docs/roadmap_tasks/`.

**Physics ground truth (stable, no status):** `Analysis/docs/analysis_overview.md`
— objective, observables, scientific methodology, efficiency strategy, sample
roles. This roadmap owns all *status*; the overview owns the *what/why*.

**Reference:** Run 2 note summary in
`.claude/kb/analysis/run2_dimuon_note.md` (full note:
`IntNotesRun2DimuonReference/`). The Run 3 note skeleton is at
`dimuon_codes/IntNotes/` (ANA-HION-2023-07-INT1) — currently bare bones
(Intro / Datasets / Analysis / Results chapters; `tex/datasets.tex` lists
MC JIRAs ATLHI-658, ATLHI-596, ATLHI-666).

**The measurement (as implemented in this repo):** single-b dimuon
analysis at √s_NN = 5.36 TeV — opposite-sign muon pairs from the same
b-decay chain (B → μμX via b→μ + b→c→μ; signal region truth-side:
minv ∈ [1.08, 2.9] GeV, pair pT > 8 GeV, |pair η| < 2.2, dR > 0.05) —
differential cross-sections in pp24 and PbPb 23/24/25 vs (pair pT,
pair η, centrality), and R_AA. *(Whether a Δφ-correlation/width
measurement à la Run 2 is also in scope is an open question for the
user — see §Q2.9.)*

---

## Q1. IntNote sections that can already be written

Proposed section layout mirrors the Run 2 note (see kb summary,
"Note structure"). Status legend: **READY** = results exist and are
current; **DRAFT+PH** = methodology settled, write now with placeholder
numbers/figures; **BLOCKED** = see Q4.

| # | Section (proposed) | Status | Basis |
|---|---|---|---|
| 1 | Introduction (physics motivation, HF energy loss, single-b method) | DRAFT+PH | Run 2 note intro + Pythia/Powheg truth studies |
| 2 | Datasets: PbPb 23/24/25 data, triggers, skim | DRAFT+PH | May 2026 skim record (`data-merging-record.txt`); event counts 124.5M/92.6M/260.4M; exact HLT names/GRLs need user confirmation (Q2) |
| 3 | Datasets: pp24 data | DRAFT+PH | May 2026 skim (12 parts); lumi needs confirmation |
| 4 | Datasets: MC samples | DRAFT+PH | ATLHI-596 (Pythia truth), ATLHI-666 (Pythia fullsim pp24 — reco-eff & det-resp), overlay r-tags (ATLHI-576/r17662 — reco-eff overlay), Powheg **truth** bb/cc (NLO template fit only). Powheg **fullsim** (ATLHI-658) obsolete — not requested for Run 3. Final production tags pending |
| 5 | PbPb event selection (5-cut sequential: banana, ZDC time, preamp, nTrk cuts) | **READY** | Cuts + plots rerun 2026-06-08 with May skim, all 3 years (`plots/single_b_analysis/event_selection*`) |
| 6 | Centrality determination (FCal ET, 2023 Glauber thresholds, FCal cross-year scaling) | **READY** (method) | `fcal_scale_pbpb_20YY.root`, centrality recalc documented; official T_AA refs needed (Q2) |
| 7 | Muon & pair selections, cutflow | **READY** (method+plots) | NTuple processing cuts; `hists_cut_acceptance_*`, signal acceptance cutflow plotter |
| 8 | Trigger efficiency: PbPb single-mu4 ε^nc(pT, q·η, ctr) + Fermi+log fits | **READY** | P2 rerun 2026-06-09 with May skim, all 3 years; 132 plots + fit canvases |
| 9 | Trigger efficiency: dR correlation corrections | DRAFT+PH | Method settled (D9: cross-term as MC reference); final numbers need unbiased MC (Q4) |
| 10 | Trigger efficiency: pp24 2mu4 | **READY** | pp24 2mu4 ε^nc + erf+log fits done 2026-06-10 (P2; PP has no Pipeline 3 by design) |
| 11 | MC truth studies (muon sources analog): Pythia flavor/origin, Powheg ancestor groups | **READY** | Pythia truth 5.36/5.02 + Powheg truth bb/cc processed, plots exist |
| 12 | Reco efficiency (PbPb overlay + pp fullsim) | BLOCKED (Q4) | Method settled; pair ε_reco(pair pT, pair η, dR). Index-mapping/HIJING bugs fixed; residual low-pT overlay deficit shown PHYSICAL (see investigation Step 23–24). Still test samples only → placeholder figures |
| 13 | Detector response / momentum imbalance | BLOCKED (Q4) | Now from **Pythia fullsim** (pp24 + overlay); pp17 Powheg det-resp kept as method demo only (Powheg fullsim obsolete). Full Pythia fullsim samples needed |
| 14 | Cross-sections (pp24, PbPb combined) | DRAFT+PH | Plots exist on May skim (PbPb 2026-06-08, pp24 2026-06-10) **with no-correlation trigger-efficiency correction applied**; still **missing reco-eff correction** + confirmed lumi/T_AA (esp. PbPb25) — preliminary |
| 15 | R_AA | DRAFT+PH | `RAA_plotting.cxx` exists; needs corrected crossx + 2025 lumi/T_AA |
| 16 | Background / fake-muon purity (Δp/p templates) | BLOCKED (Q4) | No code yet; method from Run 2 note §DpopAna |
| 17 | Systematics | DRAFT+PH (list only) | Source list adaptable from Run 2; no evaluations yet |
| 18 | Results / Summary | BLOCKED | Needs corrected results |

**Recommendation:** populate the IntNotes skeleton now with sections
1–8, 11 (and method text for 9, 14, 15, 17 with `\todo` placeholders).
That is roughly half the note.

---

## Q2. Information needed from the user (do NOT assume)

Numbers found in the repo are listed with their location, but their
provenance is unverified — **treat none of these as known facts until
confirmed.**

> **2026-06-15 — Q2.1, Q2.3, Q2.4, Q2.5 ANSWERED by the user.** Authoritative
> values now documented in **`IntNotes/analysis_metadata.md`** (private CERN
> GitLab IntNote repo), backed by `IntNotes/data/luminosity/` (per-run lumi
> CSVs) and `IntNotes/data/centrality/TaaValues2023.txt`. Summary:
> - **Q2.1 lumi** (`Prescale Corrected`, analysis triggers): pp24 2mu4 = 400.412 pb⁻¹;
>   PbPb23 mu4 = 1.17576 nb⁻¹ (corrected 2026-06-19: correct-GRL v120 total
>   1183.650457 µb⁻¹ = 1.18365 nb⁻¹ − the TWO b-hadron runs 461674 (2.623332) +
>   462964 (5.262407) µb⁻¹, both excluded at event level via PbPbBadRuns; old
>   1.02426 used a wrong/old GRL total 1029.52 and excluded only 462964 because
>   461674 was below that table's range);
>   PbPb24 mu4 = 0.85112 nb⁻¹ (851.118 µb⁻¹, GRL ≥489703 — old 1.59663 used a
>   stale GRL incl. bad runs <489703; corrected 2026-06-19); PbPb25 mu4 = 2.59933 nb⁻¹; PbPb26 = placeholder.
>   **Status:** `PbPbBaseClass.h` + `Utilities/PbPbSampledLumi.h` now carry
>   1.17576 / 0.85112 / 2.59933 nb⁻¹ (DONE via /review-analysis-code). `PPBaseClass.h`
>   pp24 (410.815 → 400.412) still to reconcile. Per-year official lumi
>   **uncertainty** still needed.
> - **Q2.3 T_AA**: 2023 values in `TaaValues2023.txt`; **2024 & 2025 reuse 2023
>   T_AA as a placeholder** (official centralities unavailable) — see reminder below.
> - **Q2.4 HLT chains** and **Q2.5 GRL names**: per-year tables in
>   `analysis_metadata.md`, from `SkimCode/scripts/TrigRates_CA.py` (consistent
>   with the `SkimCode/run_*/` XMLs). Analysis triggers: PbPb `HLT_mu4_L1MU3V`,
>   pp24 `HLT_2mu4_L12MU3V`.
>
> **Q2.2 (total PbPb hadronic σ = 7.8 b) is STILL OPEN** — needs a citable
> 5.36 TeV reference.
>
> **INTERNAL-NOTE REMINDER (Q2.3):** the 2024 and 2025 cross-section / R_AA /
> T_AA-dependent results use **2023 Glauber T_AA as a placeholder**. This MUST be
> stated explicitly in the internal note (and revisited once official 2024/2025
> centralities are released). The code carries the matching remark in
> `PbPbBaseClass.h::make_crossx_factors_pbpb_2024/2025`.

1. **Integrated luminosities (per year, per trigger, GRL-matched to our
   skim):**
   - In code (`MuonObjectsParamsAndHelpers/PPBaseClass.h`):
     pp24 2mu4 = 410.815 pb⁻¹; pp24 mu4_mu4noL1 = 113.999 pb⁻¹
     (prescale-corrected?); pp17 2mu4 = 256.8 pb⁻¹.
   - In code (`PbPbBaseClass.h` crossx factors): PbPb23 = 1.3896 nb⁻¹,
     PbPb24 = 1.5411 nb⁻¹. **PbPb25 = explicit TODO — currently a
     placeholder reusing 2023 values.** ⟵ most urgent
   - Needed: confirmed values + official lumi uncertainty per year.
2. **Total PbPb hadronic cross-section at 5.36 TeV:** code uses 7.8 b.
   Confirm value + reference (Run 2 note used 7.66 b at 5.02 TeV).
3. **Glauber ⟨T_AA⟩ / N_part for 5.36 TeV** per centrality bin +
   uncertainties + citable reference. Code hardcodes T_AA = {26.1428,
   20.3241, 14.0502, 8.5074, 3.7733, 0.6716} mb⁻¹ (source?); 2025
   values explicitly missing.
4. **Exact HLT chain names + prescale status per year:** PbPb23/24/25
   mu4 chain (and mu6/mu8 support chains), mu4_mu4noL1, 2mu4; pp24
   2mu4. Especially: why PbPb25 passmu4mu4noL1 ≈ 2.4% vs ~70–80% in
   23/24 (menu/prescale change?) — affects whether mu4_mu4noL1 is
   usable at all. *(2026-06-14: mu4_mu4noL1 is currently NOT used for pp
   or PbPb; only single-mu4 (PbPb) and 2mu4 (pp24). Revisit only if that
   decision changes.)*
5. **GRL names per year** (PbPb23/24/25, pp24). Skim investigation doc
   mentions v113 vs v120 for PbPb24 — which is final for the note?
6. **MC production details for the datasets section:** final container
   names/AMI tags, generator versions/tunes, filter efficiencies, event
   counts for ATLHI-596 (Pythia truth) and ATLHI-666 (Pythia fullsim pp24);
   HIJING overlay production status (ATLHI-576 — r17662 StandardSignalOnlyTruth
   recommended by our investigation; has full production been requested?
   expected stats & timeline?). (Powheg **fullsim** — ATLHI-658 — is obsolete
   for Run 3 and NOT part of the final MC sets; only Powheg **truth** bb/cc is
   used, for the NLO template fit. No Powheg fullsim overlay is planned.)
7. **MCP recommendations for Run 3 HI muons** (scale factors, reco-eff
   uncertainties): which release/recommendation applies, or do we go
   fully MC-driven?
8. **Muon working point for the note's nominal** (Run 2 used Tight with
   Medium as systematic; current RDF code fills medium & tight — which
   is nominal here?).
9. **Scope of the note:** crossx + R_AA only, or also Δφ away-side
   correlation widths (Run 2 style)? This decides whether the Δφ-fit
   machinery (flow-modulated Lorentzian fits) must be ported/written.
10. **Note metadata:** title, author list, abstract, reference code
    confirmation (skeleton says ANA-HION-2023-07-INT1).

*Guesses recorded (NOT facts):* pp24 2mu4 lumi ~411 pb⁻¹ and PbPb23/24
lumis above are probably what was used for the existing preliminary
crossx plots; σ_PbPb = 7.8 b is plausibly the 5.36 TeV Glauber value.
They must be confirmed before any number enters the note.

---

## Q3. Capable now, but not performed / not finished

### (a) Not yet performed (in sequence; all inputs on disk)

*(none currently)*

> **Removed 2026-06-14 (finished — see status ledger):** pp24 ntuple
> reprocessing, pp24 crossx, pp24 trig eff P2 (`task_01`, `task_02`) — all done
> on the May skim 2026-06-10. Overlay reco-eff Bug #1/#2 fixes + pipeline rerun
> (`task_03`) — done & committed (investigation Steps 11–13). PbPb23 part1
> hadd-inclusion check — verified. The cross-term dR correction is **not** an
> active reco/trig step (PbPb P3 cross-term plot dropped as biased, D9; PP has
> no Pipeline 3); dR trigger correction stays a dummy ≡ 1 (Q4).
>
> **Dropped 2026-06-14:** PbPb passmu4noL1 ntuple reprocessing — **not needed**.
> mu4_mu4noL1 is **not used** for either pp or PbPb (current decision); only
> single-mu4 (PbPb) and 2mu4 (pp24) drive the analysis.

### (b) Performed but preliminary / not refined

- **PbPb & pp24 crossx plots** — produced on the May skim (PbPb
  2026-06-08, pp24 2026-06-10) **with the no-correlation trigger-efficiency
  correction applied** (`w_trig = 1/ε_trig^pair`); still **missing the
  reco-efficiency correction**; normalization (lumi/T_AA) unconfirmed (Q2).
  Treat current crossx/R_AA numbers as preliminary until reco-eff + lumi/T_AA.
- **Overlay reco efficiency** — three bugs fixed; pair ε now ~49–55%, and the
  residual low-pT deficit vs pp is shown PHYSICAL (investigation Steps 23–24).
  Final numbers still await the full signal-only-truth (r17662) sample (Q4).
- **Powheg fullsim pp17 chain** — Run 2 conditions; **obsolete for Run 3**
  (reco-eff & det response now come from Pythia fullsim — see `docs/powheg.md`
  role note). Retained for reference only.

### (c) Precise but not propagated to final results

- **PbPb single-mu4 ε^nc + pp24 2mu4 ε^nc per-pair trigger weights** —
  **now wired into and applied in the crossx pipeline** (done 2026-06-10):
  PbPb uses inclusion–exclusion `1/(ε₁+ε₂−ε₁ε₂)` (`RDFBasedHistFillingPbPb.cxx:926`,
  `weight_for_RAA_trig_corr`), pp24 uses `1/(ε₁·ε₂)`
  (`RDFBasedHistFillingPP.cxx:373`, `crossx_weight_trig_corr`). The only
  remaining unpropagated correction is the **pair reco efficiency**
  ε_reco(pair pT, pair η, dR) → `task_05` (reco-eff part only).
- **Event-selection cuts & FCal scaling** — derived and applied
  (propagated; no action).

---

## Q4. Cannot perform yet — missing prerequisites

All are MC-sample-availability limited. For each, a **dummy strategy**
keeps the chain runnable end-to-end.

> **Consolidated placeholder list:** `docs/placeholder.md` (symlinked as
> `IntNotes/placeholder.md`) indexes every standing placeholder — 2024/2025
> centrality & ⟨T_AA⟩, reco efficiency, σ_PbPb, etc. — with code locations and
> note-disclosure requirements.

| Step | Missing prerequisite | Dummy strategy until available |
|---|---|---|
| PbPb pair reco efficiency ε_reco(pair pT, pair η, dR) | Full Pythia fullsim HIJING-overlay production (r17662 recommended; only 60k-evt r17618 + 10k-evt r17662 test samples exist). Bug #1/#2 fixes DONE (task 03) | Use bug-fixed r17618 60k-evt pair efficiencies (large stat errors) as dummy; or flat ε per Run 2-like values |
| pp pair reco efficiency ε_reco(pair pT, pair η, dR) | Full Pythia fullsim pp24 sample (`pythia_fullsim_full_sample/` empty; ATLHI-666 test only) | Test-sample efficiency or flat dummy |
| Det response / unfolding inputs | **Pythia fullsim** pp24 (+ HIJING overlay) — only test samples exist. (Powheg fullsim is obsolete: pp17-only, no pTHat → poor high-pT stats; see `docs/powheg.md`) | Use Pythia fullsim test-sample det-resp shapes as placeholder |
| ~~PbPb NLO overlay~~ | Powheg fullsim overlay — **dropped** (Powheg fullsim obsolete for Run 3) | n/a |
| Δp/p significance + template fit (purity) | Signal templates need fullsim/overlay true muons; bkg templates need π/K-enriched MC (Run 2 used jetjet JZ slices — Run 3 equivalent not identified) **and code does not exist yet** | Build the framework now (task 07) on test fullsim samples; quote purity as placeholder |
| dR trigger-correlation correction (final) | Unbiased trigger decision in MC → needs fullsim overlay with trigger sim | Cross-term on test sample; or assume ε_dR ≡ 1 with a systematic |
| Final R_AA normalization | 2025 lumi + T_AA (Q2.1/2.3) | Keep placeholder = 2023 values, flagged |
| MCP scale factors | Run 3 HI recommendations (Q2.7) | Omit (pure MC-driven) with note |

**Note on sequencing:** Bug #1/#2 fixes (task 03) and the Δp/p framework
(task 07) are *code* work we can do now; only the *final numbers* await
samples. They are deliberately placed in the "perform now with dummies"
column of the roadmap.

---

## Q5. Full analysis chain roadmap (dummies where needed)

```
                     ┌─────────────────────────────────────────────┐
 DATA (done)         │ MC (partly dummy)                            │
 ───────────         │ ──────────────                               │
 May 2026 skim       │ Pythia truth 5.36/5.02 ........ DONE         │
   PbPb23/24/25 ✓    │ Powheg truth bb/cc ............ DONE         │
   pp24 ✓            │ Pythia fullsim pp24 ........... TEST ONLY    │
        │            │ Pythia fullsim HIJING overlay . TEST ONLY    │
        ▼            │ (Powheg fullsim ............... OBSOLETE)    │
 [1] PbPb event sel ✓└─────────────────────────────────────────────┘
        │
        ▼
 [2] NTuple processing
       PbPb 23/24/25 ✓ (2026-06-08)        pp24 ✓ (2026-06-10)
        │
        ▼
 [3] Trigger efficiency
       PbPb P2 ε^nc + fits ✓ (06-09)       pp24 2mu4 P2 ✓ (06-10)
       dR corr: dropped (biased, D9 / no PP P3) → dummy ε_dR=1
        │
        ▼
 [4] Reco efficiency                        ← critical dummy
       Bug#1/#2 fixed ✓ (task_03); dummy ε_reco files → task_04
        │
        ▼
 [5] Apply per-pair weights w⁻¹ = ε_trig^pair · ε_reco(pair pT, pair η, dR)
       trigger part ✓ applied (06-10); reco-eff part → task_05
        │
        ▼
 [6] Cross-sections: pp24 + PbPb (combined years) → rerun in task_05
        │
        ▼
 [7] R_AA (needs 2025 lumi/T_AA from user) → task_06
        │
        ▼
 [8] Purity: Δp/p significance template fits (framework now, dummy
       templates) → task_07
        │
        ▼
 [9] Systematics framework (variations: WP, eff ±1σ, binning, fit
       model, ...) → task_08
        │
        ▼
 [10] Internal note writing (sections per Q1) → task_09  [can start NOW]
```

### Task files (instructions for Claude Code agents)

In `Analysis/docs/roadmap_tasks/` — one .md per step, each with
objective, inputs, exact commands/files, verification, and which
reviewer skill to invoke:

| File | Step | Can start | Blocking input |
|---|---|---|---|
| `task_04_dummy_reco_eff.md` | Dummy ε_reco(pair pT, pair η, dR) correction files | NOW | — |
| `task_05_apply_eff_to_crossx.md` | Wire **ε_reco** pair weight into crossx; rerun (trigger weight already applied) | after 04 | — |
| `task_06_raa.md` | R_AA from corrected crossx | after 05 | 2025 lumi/T_AA (user) |
| `task_07_dpop_template_framework.md` | Δp/p significance + template-fit code | NOW | bkg MC choice (user) |
| `task_08_systematics.md` | Systematics variation framework | after 05 | WP choice (user) |
| `task_09_intnote_sections.md` | Write IntNote sections (Q1 READY/DRAFT+PH) | NOW | Q2 numbers (placeholders OK) |

> **Done & removed 2026-06-14:** `task_01` (pp24 NTP + crossx),
> `task_02` (pp24 2mu4 trig eff), `task_03` (overlay reco-eff Bug #1/#2 +
> rerun). The trigger-weight half of `task_05` is also done; `task_05` now
> covers the reco-eff weight only. Task instruction files kept on disk.

Immediately startable in parallel: **04, 07** (and **09** for note text).

---

## Status ledger (update as work proceeds)

| Date | Step | Status change |
|---|---|---|
| 2026-06-10 | Roadmap created | — |
| 2026-06-14 | task_01 (pp24 NTP+crossx), task_02 (pp24 2mu4 trig eff) | DONE (May skim, 2026-06-10) — removed from Q3a/Q5 |
| 2026-06-14 | task_03 (overlay reco-eff Bug #1/#2 + rerun) | DONE (investigation Steps 11–13) — removed |
| 2026-06-14 | Trigger no-correlation correction in crossx | DONE & APPLIED (PbPb + pp24); task_05 reduced to reco-eff only |
| 2026-06-14 | Reco-eff formula corrected | pair ε_reco(pair pT, pair η, dR), not ε₁·ε₂ |
| 2026-06-14 | Powheg fullsim | marked obsolete; reco-eff & det-resp → Pythia fullsim; Powheg truth → template fit |
| 2026-06-15 | Reco-eff PLACEHOLDER applied (chain [4]/[5], task_05 reco part) | Run 2 single-muon ε_reco proxy (ε₁·ε₂) wired into crossx as a correction STAGE (`CorrectionStages.h`); PbPb from dimuon note F.2, pp from HF R_AA Fig.31; crossx reran pp+PbPb; before/after 3-line plots; both `/review-*` PASS. Proper 3D pair ε_reco still pending MC (Q4). Docs: `reco_eff_placeholder_run2.md`, `placeholder.md` item 3 |
| 2026-06-16 | Reco-eff PLACEHOLDER promoted to NOMINAL (user request; Run 3 MC ≥2-3 mo away) | `w_reco` folded into nominal corrected weight (`*_trig_corr` = base·w_reco·w_trig) → all crossx histos + R_AA 3D input now reco+trig corrected (== validated reco_trig stage; verified pp 5739, PbPb 3D ×1.67). `/review-analysis-code` PASS. Crossx reran + nominal crossx plots reran (pp24, pbpb_23_24_25_combined). Pre-reco backup `crossx_hist_backup_20260616_pre_reco_nominal/`. Crossx/R_AA remain placeholder/preliminary |
| 2026-06-16 | **task_06 R_AA DONE** (reco-corrected R_AA runs) | Added SS signal-region histos to RDF crossx (PbPb 3D `h3d_ss_..._vs_centr`, pp 2D `h2d_ss_...`) for OS−SS combinatorial subtraction (`/review-analysis-code` PASS). Modernized `RAA_plotting.cxx` (cluster paths, RDF inputs case 6, combined PbPb 23+24+25 vs pp24, OS−SS, 15-bin pT, index-wise ratio, segfault+mode-3 off-by-one fixed; `/review-analysis-code` PASS iter2). R_AA plots vs pair pT/η/centrality (`/review-plot` PASS iter4). 2025 lumi verified 2.59933 nb⁻¹; 2023 T_AA placeholder. Tracking: `raa_from_rdf_crossx.md`. |
| 2026-06-18 | **PbPb reco placeholder → colleague's exact Run 2 Medium fits** | Replaced eyeball F.2 arrays in `run2_reco_eff_placeholder.root` with dense samples of `MuonRecoEffcyRun2MC_medium.root::tf1_eff_fit_cent{C}_eta{E}` (ctr map {12,13,4,5,6,7,8}, q·η slice i↔eta{i}); same 63 graph names ⇒ no lookup change. Builder `/review-analysis-code` PASS; PbPb crossx RDF (23/24/25)+crossx/R_AA/stage plots reran, `/review-plot` PASS. dσ shifts ctr0_5 +4.2%/ctr30_50 −5.6%/ctr50_80 −4.3%; reco/raw 1.43–1.81. pp unchanged. Tracking `reco_eff_placeholder_run2.md` Step 11 |
| 2026-06-16 | **Year-combination normalization FIXED** (crossx + R_AA) | Switched to luminosity-weighted average `Σ(L_y·h_y)/ΣL_y` (HF R_AA note HION-2019-58 §4.1 Eq.3) via single-source `Utilities/PbPbSampledLumi.h` (2023=1.02426, 2024=1.59663, 2025=2.59933 nb⁻¹). Applied in R_AA `HistRetrieve`, crossx combined plotter `GetHistObject` (`_counts`=simple sum), and the before/after stage plotter. `/review-analysis-code` + `/review-plot` PASS. **R_AA scale now physical (~0.1–0.9; was ~1.5–2.5).** Remaining R_AA-scale caveats: 2023 T_AA placeholder, σ_PbPb 5.36 TeV guess, Run 2 reco placeholder |
| 2026-06-19 | **PbPb 2024 lumi GRL correction** | 1.59663→0.85112 nb⁻¹ (GRL ≥489703; old GRL over-counted, events already correct). PbPbBaseClass.h + PbPbSampledLumi.h + docs; 2024 crossx refilled; combined R_AA ×1.17. `/review-analysis-code` + `/review-plot` PASS. Tracking `raa_from_rdf_crossx.md`. |
| 2026-06-19 | **PbPb 2023 lumi GRL correction + b-hadron-run subtraction** | 1.02426→1.17576 nb⁻¹ (correct-GRL v120 total 1.18365 nb⁻¹ minus the two b-hadron runs 461674+462964, both already excluded at event level via PbPbBadRuns). PbPbBaseClass.h (6 factors) + PbPbSampledLumi.h case 23 + docs (metadata, lumi README, datasets.tex, roadmap, status); 2023 crossx refilled. `/review-analysis-code` PASS. Combined R_AA × 0.967 (ΣL 4.475→4.626 nb⁻¹). Tracking `raa_from_rdf_crossx.md`. |
