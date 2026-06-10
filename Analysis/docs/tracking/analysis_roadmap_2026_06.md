# Analysis Roadmap — Run 3 Dimuon Internal Note (ANA-HION-2023-07-INT1)

**Status date: 2026-06-10.** Living document; update statuses as steps
complete. Companion checklist tracking doc:
`analysis_roadmap_buildout.md`. Per-task agent instructions:
`Analysis/docs/roadmap_tasks/`.

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
| 4 | Datasets: MC samples | DRAFT+PH | ATLHI-596 (Pythia truth), ATLHI-666 (Pythia fullsim pp24), ATLHI-658 (Powheg fullsim pp24), overlay r-tags (ATLHI-576 test); final production tags pending |
| 5 | PbPb event selection (5-cut sequential: banana, ZDC time, preamp, nTrk cuts) | **READY** | Cuts + plots rerun 2026-06-08 with May skim, all 3 years (`plots/single_b_analysis/event_selection*`) |
| 6 | Centrality determination (FCal ET, 2023 Glauber thresholds, FCal cross-year scaling) | **READY** (method) | `fcal_scale_pbpb_20YY.root`, centrality recalc documented; official T_AA refs needed (Q2) |
| 7 | Muon & pair selections, cutflow | **READY** (method+plots) | NTuple processing cuts; `hists_cut_acceptance_*`, signal acceptance cutflow plotter |
| 8 | Trigger efficiency: PbPb single-mu4 ε^nc(pT, q·η, ctr) + Fermi+log fits | **READY** | P2 rerun 2026-06-09 with May skim, all 3 years; 132 plots + fit canvases |
| 9 | Trigger efficiency: dR correlation corrections | DRAFT+PH | Method settled (D9: cross-term as MC reference); final numbers need unbiased MC (Q4) |
| 10 | Trigger efficiency: pp24 2mu4 | BLOCKED→Q3 | Code ready (§3d); needs pp24 ntuple rerun (task 02) |
| 11 | MC truth studies (muon sources analog): Pythia flavor/origin, Powheg ancestor groups | **READY** | Pythia truth 5.36/5.02 + Powheg truth bb/cc processed, plots exist |
| 12 | Reco efficiency (PbPb overlay + pp fullsim) | BLOCKED (Q4) | Test samples only + 2 open bugs; write method text with placeholder figures |
| 13 | Detector response / momentum imbalance | BLOCKED (Q4) | Powheg pp17 det-resp exists as method demo; pp24/overlay samples needed |
| 14 | Cross-sections (pp24, PbPb combined) | DRAFT+PH | Plots exist (2026-06-08 PbPb; pp24 from old skim) but **uncorrected** for trig/reco eff — preliminary |
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
   usable at all.
5. **GRL names per year** (PbPb23/24/25, pp24). Skim investigation doc
   mentions v113 vs v120 for PbPb24 — which is final for the note?
6. **MC production details for the datasets section:** final container
   names/AMI tags, generator versions/tunes, filter efficiencies, event
   counts for ATLHI-596/666/658; HIJING overlay production status
   (ATLHI-576 — r17662 StandardSignalOnlyTruth recommended by our
   investigation; has full production been requested? expected stats &
   timeline?). Also Powheg fullsim **overlay** (PbPb conditions) plans.
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

| Step | What | Why pending | Task file |
|---|---|---|---|
| pp24 ntuple reprocessing | NTP on May 2026 skim (12 parts); current `muon_pairs_pp_2024_*` are Feb/Mar 2026 (old skim) | Just not run since skim landed | `task_01` |
| pp24 crossx | RDF hist fill + `plot_single_b_crossx_pp` on new ntuples | Depends on previous | `task_01` |
| pp24 trig eff P2+P3 | 2mu4 ε^nc + ε_dR^{2mu4} (code implemented, §3d) | Depends on pp24 ntuples | `task_02` |
| PbPb passmu4noL1 ntuple reprocessing | Skimmer-side fix applied in code; ntuples not reprocessed | Only needed if mu4_mu4noL1 is used (see Q2.4) | — (decide first) |
| Overlay reco-eff Bug #1 + #2 fixes | Index-mapping fix + dR-fallback matching (specs written in `hijing_overlay_reco_effcy_investigation.md` Step 10) | Awaiting implementation decision | `task_03` |
| Cross-term dR correction on MC | Run P3-style cross-term on fullsim (unbiased) sample | Test-sample stats only for now → preliminary | `task_03` (validation step) |
| PbPb23 part 1 raw re-download | `data_pbpb23_part1.root` was missing at one point (mu4 doc Step 16); May skim record lists part1 (24.9M, v2) — verify current hadd includes it | Verification, not blocked | `task_01` (pre-check) |

### (b) Performed but preliminary / not refined

- **PbPb & pp24 crossx plots** — produced (PbPb 2026-06-08 with May
  skim; pp24 with old skim) but **without trigger- or reco-efficiency
  corrections**; normalization (lumi/T_AA) unconfirmed (Q2). Treat all
  current crossx/R_AA numbers as shape-level preliminary.
- **Overlay reco efficiency** — produced but unphysical (~20%); three
  bugs identified, one fixed; numbers invalid until task 03.
- **Powheg fullsim pp17 chain** (single-muon det resp + mixed pairs) —
  complete, but Run 2 (2017) conditions: serves as method demonstration
  only; final needs pp24 sample (Q4).

### (c) Precise but not propagated to final results

- **PbPb single-mu4 ε^nc(pT, q·η, ctr) fits** (P2, all 3 years, May
  skim, fine q·η) — exist in
  `pbpb_20YY/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root`
  but are **not applied as per-pair weights** in the crossx pipeline
  (P1 deliberately applies no trigger-efficiency correction). The
  inclusion–exclusion pair weight 1/(ε₁+ε₂−ε₁ε₂) needs to be wired
  into the crossx hist filling. → `task_05`.
- **Event-selection cuts & FCal scaling** — derived and applied
  (propagated; no action).

---

## Q4. Cannot perform yet — missing prerequisites

All are MC-sample-availability limited. For each, a **dummy strategy**
keeps the chain runnable end-to-end.

| Step | Missing prerequisite | Dummy strategy until available |
|---|---|---|
| PbPb pair reco efficiency ε_reco(pT, q·η, ctr) | Full Pythia fullsim HIJING-overlay production (r17662 recommended; only 60k-evt r17618 + 1k-evt r17662 test samples exist) **and** Bug #1/#2 fixes (task 03 — doable now) | After task 03: use bug-fixed r17618 60k-evt efficiencies (large stat errors) as dummy; or flat ε per Run 2-like values |
| pp pair reco efficiency | Full Pythia fullsim pp24 sample (`pythia_fullsim_full_sample/` empty; ATLHI-666 test only) | Test-sample efficiency or flat dummy |
| NLO det response / unfolding inputs, NLO eff cross-check | Powheg fullsim pp24 (ATLHI-658) — not produced; only pp17 exists | Use pp17 Powheg det-resp shapes as placeholder |
| PbPb NLO overlay | Powheg fullsim overlay samples (classes exist for yr 18, 23–26; no samples) | Skip; not on critical path |
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
        ▼            │ Powheg fullsim ................ pp17 ONLY    │
 [1] PbPb event sel ✓└─────────────────────────────────────────────┘
        │
        ▼
 [2] NTuple processing
       PbPb 23/24/25 ✓ (2026-06-08)        pp24 ✗ → task_01
        │
        ▼
 [3] Trigger efficiency
       PbPb P2 ε^nc + fits ✓ (06-09)       pp24 P2/P3 ✗ → task_02
       dR corr: cross-term ref ✓ (D9);  final from MC → dummy ε_dR=1
        │
        ▼
 [4] Reco efficiency                        ← critical dummy
       fix Bug#1/#2 → task_03; dummy ε_reco files → task_04
        │
        ▼
 [5] Apply per-pair weights w⁻¹ = ε_trig_pair · ε_reco,1 · ε_reco,2
       → task_05 (new: currently NO corrections applied in crossx)
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
| `task_01_pp24_reprocess.md` | pp24 NTP + crossx on May skim | NOW | — |
| `task_02_pp24_trig_eff.md` | pp24 2mu4 P2+P3 | after 01 | — |
| `task_03_overlay_recoeff_bugfix.md` | Bug #1/#2 fixes + rerun overlay pipeline | NOW | — |
| `task_04_dummy_reco_eff.md` | Dummy ε_reco correction files | after 03 | — |
| `task_05_apply_eff_to_crossx.md` | Wire ε_trig (+ ε_reco) weights into crossx; rerun | after 02, 04 | — |
| `task_06_raa.md` | R_AA from corrected crossx | after 05 | 2025 lumi/T_AA (user) |
| `task_07_dpop_template_framework.md` | Δp/p significance + template-fit code | NOW | bkg MC choice (user) |
| `task_08_systematics.md` | Systematics variation framework | after 05 | WP choice (user) |
| `task_09_intnote_sections.md` | Write IntNote sections (Q1 READY/DRAFT+PH) | NOW | Q2 numbers (placeholders OK) |

Three tasks are immediately startable in parallel: **01, 03, 07** (and
**09** for note text).

---

## Status ledger (update as work proceeds)

| Date | Step | Status change |
|---|---|---|
| 2026-06-10 | Roadmap created | — |
