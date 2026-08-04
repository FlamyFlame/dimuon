# Systematic Uncertainties — authoritative registry

**Created:** 2026-08-03. **Status:** living document.

This is the single place that records **what** each systematic uncertainty is, **where its
inputs live**, and **how it is evaluated**. It is a registry, not an evaluation: sizes appear
here only once they have actually been measured, and every number carries a pointer to the
file or tracking doc it came from. An entry with no measured size says so explicitly — do not
infer one.

**Companion documents**
- `docs/roadmap_tasks/task_08_systematics.md` — the *variation framework* plan (how variations
  are wired as config flags, the `_syst_<source>_<updown>` output convention, the ratio +
  smoothing plotter). That file is the implementation plan; this file is the registry.
- `docs/muon_wp_registry.md` — every site where the muon working point is configured
  (the Tight↔Medium variation touches all of them).
- `docs/signal_selection_change_impact.md` — **MUST-READ before varying any signal-selection
  cut**: the full recompile → refill → replot blast radius of a selection change.
- `docs/tracking/INDEX.md` — the tracking docs that own each measurement.

**Method (inherited from the Run 2 dimuon note, kb `run2_dimuon_note.md` §Systematics):**
vary → take the ratio to nominal → smooth/fit the ratio → sum sources in quadrature.

---

## 1. Trigger efficiency

### 1a. ΔR-correlation correction — large-ΔR plateau normalization  ⚠ ACTIVE, size to be re-derived

**What it is.** The MC-derived ΔR corrections are ratios that must, by construction, tend to 1
for well-separated muons (`mc_trigger_efficiency.md` §3.3/§3.4 diagnostic 2): a plateau ≠ 1
measures the quality of the ε_MC turn-on parameterization used in the inverse weighting, not
physics. The corrections are therefore **plateau-normalized** before use, and the size of the
offset that was normalized away is taken as a systematic on the correction.

**Where the numbers live — READ THESE, DO NOT COPY THEM INTO PROSE:**
- **Machine-readable (authoritative): the plateau ROOT file** written by the Step-3/Step-4
  stage of `plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff.cxx`, one per sample and
  working point, next to that sample's histogram files:
  `<sample dir>/dr_correction_plateaus_<label>[_medium_wp].root`
  (`<sample dir>` = `~/usatlasdata/pythia_fullsim_full_sample/` for pp24,
  `~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/` for the HIJING overlay). It holds
  the large-ΔR plateau **per (pair pT, pair η) cell** — value, statistical error, RMS scatter
  and the number of contributing ΔR bins — for both Step 3 and Step 4, plus the inclusive
  plateau. The ΔR-correction fit stage reads it directly; nothing downstream may hardcode a
  plateau value or parse it out of a .md/.txt.
- **Human-readable mirrors** (same content, for eyeballing only):
  `step3_plateau_pair_eta_pt.txt`, `step3_plateau_fluctuation_pair_eta_pt.txt` and the Step-4
  equivalents, under
  `~/usatlasdata/dimuon_data/plots/{pp,pbpb}_trigger_efficiency/mc_based/step{3,4}_*/`.
- **Narrative + current values:** `docs/tracking/mc_trigger_efficiency.md` (R9, R12, R13 and
  Remaining Work 6).

**Guard (TWO TIERS, user 2026-08-04).** For any **FULL** sample the fit stage
(`fit_dr_corrections.cxx`) classifies every (pair pT, pair η) cell:
`|plateau − 1| > 0.10` → **FLAGGED** (reported loudly in `plateau_guard_report.txt` and on
stdout, allowed, and the cell carries `|plateau − 1|` as a systematic);
`|plateau − 1| > 0.15` → **FAILURE**, the stage throws and exits. The 0.10 tier is already a
deliberately loose closure requirement; 0.15 is the "this cell is not usable" line.
A cell is additionally marked **unusable** (`h_stepN_fit_ok = 0`) if its fit did not converge or
if the fitted correction is not positive over the whole fitted range — **consumers must require
`fit_ok == 1`**. The HIJING overlay and the r17663 no-overlay sample
are 10 000-event **TEST** samples whose high-pair-pT cells are noise-dominated, so they are
exempt from the guard (their offending cells are reported, not fatal). The guard lifts to the
overlay once the full-statistics PbPb production lands.

**⚠ TWO pp24 CELLS ARE FLAGGED (0.10 < |plateau−1| ≤ 0.15) — accepted and normalized (user
decision 2026-08-04); NONE fails.** Numbers below are on the CANONICAL pair-pT binning
{8, 13.75, 23.63, 40.62, 120} (see `.claude/CLAUDE.md` §Binnings; the pre-2026-08-04 numbers were
measured on the retired {8,15,27,50,150} edges and are superseded). 34 of the 36 pp24
(pair pT × pair η) cells have plateaus within 7% of 1; two are flagged, both **Step 3** in the
top pair-pT bin, and **Step 4 has none**:

| step | cell | plateau | deviation from 1 |
|---|---|---|---|
| 3 | pT_pair[41,120) × η_pair[−2.4,−2.0) | 1.1490 ± 0.0105 | 14.2σ, \|Δ\| = 0.149 |
| 3 | pT_pair[41,120) × η_pair[1.0,1.5) | 0.8648 ± 0.0277 | 4.9σ, \|Δ\| = 0.135 |

Neighbouring cells in the same pair-pT row are 0.95–1.00, so these are isolated outliers rather
than a trend, and neither is a statistical artefact. **Decision: normalize each cell by its own
measured plateau as usual and carry `|plateau − 1|` in these two cells as a systematic on the
ΔR correction.** The guard stays in place and reports loudly.

**Unusable cells (`h_stepN_fit_ok = 0`).** pp24: **0** for both steps — every cell is usable.
HIJING overlay: **31 of 36** (Step 3) and **13 of 36** (Step 4) — with 10 000 events the per-cell
plateaus are not measurable, so **only the inclusive overlay correction is usable today**
(inclusive plateau 0.8681 ± 0.0189 Step 3, 0.9476 ± 0.0088 Step 4, Tight). Per-cell PbPb
corrections need the full-statistics overlay.

**Open lead on the cause (NOT yet demonstrated):** the MC single-muon turn-on TF1s are fitted
only over `pT ∈ [4, 60] GeV` and the evaluator clamps above 60 GeV. In a 50–150 GeV *pair* pT bin
a large fraction of legs exceed 60 GeV, so their ε is a clamped extrapolation rather than a
measurement — the same family as the out-of-range TF1 trap in
`docs/tracking/pp_trig_eff_highpt_jump.md`. **Caveat: the sign does not obviously match** — a
clamp that under-estimates ε at high pT biases the inverse-weighted ratio *high*, while these
cells sit *low*. Extending the turn-on fits above 60 GeV and re-measuring would settle it.

**⚠ Size must be RE-DERIVED (2026-08-03).** Round 7 corrected the statistical errors on these
corrections: they had been over-stated by 1.5–3.2× by an independent-error `TH1::Divide` on a
ratio whose numerator is a subset of its denominator (`mc_trigger_efficiency.md` R12). With the
correct errors the plateau offset is no longer consistent with 1 — e.g. pp24 Step 4 inclusive
went from "0.993 ± 0.0008" to **0.9911 ± 0.0004, i.e. 22σ from 1**. Any systematic previously
sized by comparing |plateau − 1| with its error is therefore **too small** and must be redone.

### 1b. MC over-efficiency of the second muon leg

MC over-estimates the second-leg (L1 RoI × HLT) efficiency by ~20%: MC/data
`f(2mu4|mu4) = 1.23` against `f(mu4_mu4noL1|mu4) = 1.02` (`mc_trigger_info_skim.md` §8d, CLOSED;
per-leg L1 ratio ≈ 1.13, turn-on dominated). This is why the analysis keeps the **data-derived**
ε^nc and uses MC only for the ΔR **ratios**, in which a smooth per-leg offset largely cancels.
The residual is carried as a systematic on the ΔR correction. **Size not yet evaluated.**

### 1c. Forward-endcap trigger configuration (r16578)

The forward-endcap L1/HLT anomaly at `|q·η| > 2`, `pT ≲ 6 GeV` is a property of the r16578
production **configuration**, not of HIJING occupancy — established by the r17663 no-overlay
sample (`mc_trigger_efficiency.md` R8, R10). The failing step is deliberately NOT localized
(our observables cannot separate an L1 threshold/coincidence effect from the HLT hypothesis or
the chain matching), and a question to the trigger group is pending (Remaining Work 3).
**Open:** whether this is handled by a selection requirement (rejecting muons with
`pT < 7 GeV && q·η < −2`) or by a systematic — see `mc_trigger_efficiency.md` round 7.

### 1d. ε^nc turn-on parameterization

Fit-parameter ±1σ on the data tag-and-probe turn-on fits, and fit-vs-interpolation below
8 GeV. Run 2 precedent: small (`task_08_systematics.md`). **Not yet evaluated.**

**Related trap (not a systematic — a bug to avoid):** compiled TF1 turn-on fits return 0 above
their range on read-back, which floors ε at 0.01 and produced a spurious pp cross-section jump
at pair pT 50–60 GeV (`docs/tracking/pp_trig_eff_highpt_jump.md`, CLOSED — root cause found,
fix awaiting a user decision). Persist fits as TFormula-based TF1s and clamp at the point of use.

## 2. Muon working point

Tight (nominal, since 2026-07-08) ↔ Medium. Every plot set, efficiency and fit must expose a WP
config variable defaulting to Tight; the complete list of sites is
**`docs/muon_wp_registry.md`**. Implementation history: `docs/tracking/tight_wp_default_change.md`.
Run 2 size: 0.5–2% (`task_08_systematics.md`). **Run 3 size not yet evaluated.**

## 3. Reconstruction efficiency

Currently a **placeholder** in both systems, so its uncertainty is not yet meaningful:
- PbPb: a colleague's exact Run 2 Medium TF1 fits (`MuonRecoEffcyRun2MC_medium.root`), stored as
  TF1s and evaluated at the exact pT.
- pp: HF-muon note HION-2019-58 Fig. 31, with peripheral PbPb as the fallback.
See `docs/tracking/reco_eff_placeholder_run2.md` (CLOSED) and
`docs/tracking/pp_reco_eff_placeholder.md`. The nominal target is a genuine 3D **pair**
efficiency ε_reco(pair pT, pair η, ΔR) — not ε₁·ε₂ — from the Run 3 MC.
**Uncertainty cannot be assigned until the Run 3 measurement replaces the placeholder.**

## 4. Background subtraction / signal extraction

Owned by the `low_mass_dimuon_template_fit.md` umbrella and its sub-docs A–E. Sources include
the OS→SS factor k and the mixed-event model (sub-doc B), the abandoned combined OS+SS fit
retained *as* a systematic (sub-doc B), the background composition and normalization (sub-doc C),
and the Δp/p fake-muon yield program (sub-doc D). Run 2 purity systematic: <2%.
**Run 3 sizes not yet evaluated.**

## 5. Signal-selection variations

Vary the m_μμ window, the pair-pT threshold and the ΔR requirement. **Before changing any of
these, read `docs/signal_selection_change_impact.md`** — a selection change invalidates the
hist-filling outputs and every downstream cross-section/R_AA plot, and the doc enumerates
exactly what goes stale.

## 6. Event selection (PbPb) and centrality

- Event selection: nominal vs the `_alt` banana cuts, already produced as
  `event_sel_cuts_pbpb_20YY_alt.root` (5-cut sequential selection;
  `docs/tracking/PbPb_JEDI...`/`pbpb_pipelines.md`).
- Centrality calibration: vary the FCal scale factors (`fcal_scale_pbpb_20YY.root`, applied in
  `FillMuonPairExtra`, with centrality recalculated via `GetCentralityPbPb2023` for 2024/2025).
**Neither evaluated yet.**

## 7. Normalization

Luminosity (per-run CSVs in `IntNotes/data/luminosity`, "Prescale Corrected" column; PbPb in
µb⁻¹ → normalize to nb⁻¹, pp mixed pb⁻¹/nb⁻¹) and T_AA for R_AA. Run 2 luminosity uncertainty:
1.5–1.6%. Affects the normalization only, not the shapes.

**MC normalization trap (not a systematic — a hard blocker):** AMI cross-section files are keyed
by beam + pT-hat slice only, so pointing the code at a new MC production with a stale `ami_info/`
silently produces wrong, slice-dependent weights that cancel **nowhere**. See
`docs/ami_weights.md` (BLOCKING) before adopting any new MC dataset.

---

## Maintenance rule

When a systematic is evaluated, add its size **and the path of the file the number came from**
to the relevant section above, and bump the date. When a measurement that feeds a systematic is
re-run, check whether the entry's "size" line is now stale and say so explicitly rather than
leaving a number whose provenance no longer holds.
