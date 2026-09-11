# PbPb 2026 Analysis-Code Support (all stages) + 2026 Luminosity

**Created:** 2026-09-10 · **Mode:** Implementation

## Objective

Make the whole Pb+Pb analysis chain year-complete for a **fourth Run-3 heavy-ion
period, Pb+Pb 2026**, so that **combined 2023 + 2024 + 2025 + 2026** cross-section and
$R_{AA}$ results can be produced exactly as the 2023+24+25 combination is produced today.

Three deliverables:
1. Register the 2026 `HLT_mu4` luminosity table (`IntNotes/data/luminosity/pbpb_2026/`)
   and update the luminosity README.
2. Update the analysis status / metadata documentation for the new period.
3. Add 2026 to **every stage** of the analysis code (NTuple processing, event selection,
   centrality/FCal, RDF hist filling, efficiency lookups, pipelines, plotting, R_AA
   combination), so every Pb+Pb pipeline runs for 2026 the moment the skim NTUPs land.

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. `IntNotes/data/luminosity/README.md` documents `pbpb_2026` (file, unit, total,
     per-year table, GRL cross-check, bad-run caveat).
  2. `IntNotes/analysis_metadata.md` carries the real 2026 lumi / GRL / trigger-chain /
     prescale / T_AA-placeholder rows (no "placeholder (data not yet in skim)" left for
     items that are now known), and `Analysis/docs/tracking/analysis_status_summary.md`
     has a 2026 row/section stating what exists and what is pending.
  3. Every code site that is keyed by Pb+Pb year — enumerated exhaustively and
     cross-checked by four independent sweeps — either handles 2026 or is explicitly
     recorded here as deliberately unchanged, with the reason.
  4. Everything that must compile, compiles (ACLiC / the pipelines' own build steps);
     every Pb+Pb pipeline script accepts `26` and its dry-run/pre-flight passes.
  5. Every input that cannot be known until the 2026 skim finishes (number of part
     files, entry counts, bad-run list, per-run ZDC preamp cuts, FCal scale factors)
     is implemented as an explicit, clearly-labelled **placeholder guess** listed in
     one table in this doc, so it can be confirmed and corrected in one pass.
  6. No 2026 code path silently falls back to another year's constants: every year
     switch either has a real 2026 branch or throws.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

### 1. Motivation

The measurement's statistical reach is set by the integrated Pb+Pb luminosity. Adding
2026 (`Prescale Corrected` = 2.62316 nb⁻¹) to 2023+24+25 (4.62621 nb⁻¹) raises the total
to **7.24937 nb⁻¹**, ≈ 3.7× the Run-2 dimuon analysis (1.94 nb⁻¹) and ≈ 1.57× the
current Run-3 combination. The years are combined, never plotted separately
([[feedback_pbpb_crossx_combined]]), so a year that is silently omitted from a combined
loop is invisible in the output.

### 2. Top-level equation (what the year enters)

For each centrality class $c$ and each analysis bin $i$, the combined Pb+Pb
differential yield/cross-section is a **luminosity-weighted sum over years**:

$$\frac{d\sigma}{dx}\bigg|_{i,c} \;=\; \frac{\sum_{y}\, n^{y}_{i,c}\;/\;
\big(\varepsilon^{y}_{\rm reco}\,\varepsilon^{y}_{\rm trig}\big)}
{\sum_{y} L^{y}_{\rm int}\;\cdot\;\Delta x_i \cdot f^{y}_{c}}\,,
\qquad y \in \{23,24,25,\mathbf{26}\}$$

and $R_{AA}$ divides that by the pp reference scaled by $\langle T_{AA}\rangle_c$.
Symbols: $n^{y}_{i,c}$ = raw pair count of year $y$; $\varepsilon^{y}_{\rm reco/trig}$ =
that year's reconstruction / trigger efficiency; $L^{y}_{\rm int}$ = that year's
`Prescale Corrected` luminosity **minus the luminosity of any run excluded at event
level**; $f^{y}_{c}$ = the centrality-class event fraction; $\Delta x_i$ = bin width.

**Consequence for the code:** a year is not a cosmetic label. It enters (a) the yield
numerator, (b) the luminosity denominator, (c) the efficiency lookups, (d) the
centrality definition. Adding 2026 to some of these and not others produces a wrong
number with no error message.

### 3. Step-by-step method — what 2026 must inherit, and from where

a. **Skim → NTuple processing.** 2026 raw NTUPs come from the `hi2026` skim run mode,
   which is *procedurally identical* to 2023/24/25 (same `TrigRates` algorithm, muon
   tools, trigger chain lists, stored branches; only the GRL and run list differ) — see
   `pbpb2026_skim_and_lgd_storage.md` §Physics Procedure. The NTuple processing must
   therefore treat 2026 exactly as it treats 2025: same cuts, same trigger mode
   (`trigger_mode=1`, single `mu4`), same resonance-cut mode, same derived branches.
b. **Event selection.** The 5-cut sequential Pb+Pb event selection
   ([[project_pbpb_event_selection]]) is re-derived per year from that year's data
   (ZDC/FCal distributions are year-dependent). 2026 needs its own cuts ROOT file.
   Cut 3 (ZDC preamp) is per-run $\mu+7\sigma$ for 2025 and a hard scalar for 23/24
   ([[project_pbpb_preamp_cut]]); which form 2026 takes is decided from the 2026 data.
c. **Centrality.** All Run-3 Pb+Pb years use the **PbPb2023 FCal-$E_{\rm T}$ → centrality
   thresholds** (registered user decision **D6**, 2026-09-10). `MuonPairPbPb::PairValueCalcHook`
   recomputes the percentile from FCal $E_{\rm T}$ for 2025 and 2026, which **enforces** that
   calibration rather than trusting whatever the skim wrote.

   **CORRECTED 2026-09-11 — the previous justification was false.** This section, and a
   long-standing comment in `MuonPairPbPb.h`, claimed the Pb+Pb 2025 skim writes an *all-zero*
   `centrality` branch, so the recompute was a rescue. Measured directly on
   `data_pbpb25_part1.root` (2 M events): **mean 20.454, max 84, only 4.08 % zeros** —
   populated, and indistinguishable in character from 2023 (19.947 / 84 / 4.20 %). It is
   moreover **already on the 2023 calibration**: it agrees with the 2023-table recompute in
   **299 992 of 300 000** events, the 8 exceptions differing by exactly one unit and
   attributable to the rounded mirror table (see the parked Glauber item below). So the
   recompute is a **no-op for 2025**, not a repair — and the 2026 branch was added on a premise
   that does not hold, even though the resulting behaviour is exactly what D6 asks for.
   No number changes; the recorded reason does.

   **No cross-year FCal rescaling exists in the code** — `fcal_scale` / `fcal_corr_weight`
   appear in documentation only, never in a source file (removed in `40e67c9`; confirmed by
   four independent sweeps). 2026 needs no `fcal_scale_pbpb_2026.root`. Note also that
   `FCal_ET_Bins_PbPb2024` is byte-identical to the 2023 vector, so `case 24:` calling
   `GetCentralityPbPb2024` puts 2024 on the same calibration as the rest.
d. **Luminosity.** 2026 `Prescale Corrected` total = **2623.16 µb⁻¹ = 2.62316 nb⁻¹**
   from `lumitable_pbpb_26_HLT_mu4.csv`, whose run list is byte-identical to the skim
   GRL `physics_HI2026_50ns_noIBL.xml` (35 runs, 522041–523437). Enters
   `PbPbBaseClass.h` crossx factors and `Utilities/PbPbSampledLumi.h`.
e. **$T_{AA}$.** No official 2026 Glauber values; 2026 reuses the **2023 $T_{AA}$
   placeholder**, exactly as 2024 and 2025 already do (`analysis_metadata.md` §3).
   This is an existing, already-flagged placeholder — not a new one.
f. **Efficiencies.** Trigger efficiency is measured per year from that year's own data
   (tag-and-probe $\varepsilon_{\mu4}(p_{\rm T}, q\!\cdot\!\eta)$, pipelines P2+P3), so
   2026 needs its own P2/P3 run. Reconstruction efficiency is currently the Run-2
   placeholder keyed by centrality × $q\!\cdot\!\eta$ only (not by year)
   ([[project_pbpb_reco_eff_placeholder]]), so it needs no 2026 entry.

### 4. Negative constraints

- Do **NOT** change any binning to accommodate 2026. The pair-p_T, pair-$\eta$ and
  single-muon-p_T axes come from `ParamsSet` / `CommonEffcyConfig` and are the same for
  every year — that is what makes the years addable (project CLAUDE.md BLOCKING rule).
- Do **NOT** let 2026 inherit another year's *measured* constants silently. Luminosity,
  FCal scale factor, event-selection cuts and trigger-efficiency fits are all measured
  per year; a `default:` branch that quietly hands 2026 the 2025 values is a silent
  physics error. Where a 2026 measurement does not exist yet, the code must either throw
  or carry a loudly-labelled placeholder that is registered in this doc.
- Do **NOT** apply the 2023 b-hadron-run exclusion logic to 2026 by analogy. 2026 has no
  bad-run list yet; if one is later defined, its luminosity must leave
  `PbPbSampledLumi` / `PbPbBaseClass` in the same commit (numerator and denominator must
  stay consistent).
- Do **NOT** produce per-year 2026 crossx plots as a new output family; Pb+Pb crossx is
  always the combined result ([[feedback_pbpb_crossx_combined]]).
- Do **NOT** change the `hi2026` skim (owned by `pbpb2026_skim_and_lgd_storage.md`).

## Context

- **Sibling ACTIVE doc:** `pbpb2026_skim_and_lgd_storage.md` owns the skim, the grid
  submission/download, and the LOCALGROUPDISK migration. It explicitly declares
  analysis-code support OUT of its scope and hands off these values:

  | Item | Value |
  |---|---|
  | Download dir | `~/usatlasdata/dimuon_data/pbpb_2026/` |
  | Merged NTUP names | `data_pbpb26_part<N>.root` |
  | Tree name | `HeavyIonD3PD` |
  | Skim run mode | `hi2026` |
  | GRL | `SkimCode/xmls/physics_HI2026_50ns_noIBL.xml` |
  | grid_monitor mapping | `PbPb2026data...partN._EXT0` → `pbpb_2026/data_pbpb26_partN.root` |

  That doc also records that `scripts/grid_monitor.sh::get_code_update_info()`
  deliberately returns **empty** for `pbpb_2026` — i.e. the automatic
  `file_batch_max` / `.sub` `queue N` update is **skipped** for 2026 and must be set by
  hand here once the part count is known.
- Existing period reference points: PbPb 2023 = 4 parts / 124,473,950 events;
  2024 = 2 parts / 92,650,031; 2025 = 6 parts / 260,392,022.
- Rerun blast radius for any selection change: `Analysis/docs/signal_selection_change_impact.md`.

## Scope

**IN:** the 2026 lumi table registration + README; `analysis_metadata.md` and
`analysis_status_summary.md` updates; 2026 support in NTuple processing, event
selection, centrality/FCal, RDF hist filling, efficiency lookups, Condor run scripts,
pipeline drivers, plotting and R_AA combination; compile checks; the placeholder
registry below.

**OUT:** the 2026 skim itself and the LGD migration (sibling doc); producing 2026
physics results (the data does not exist yet); internal-note prose (`IntNotes/tex/`),
which is gated by the academic-writing chain and is listed under Remaining Work.

## Placeholder Registry (guesses to CONFIRM after the skim finishes)

Every value here is a **guess** made because the 2026 skim is still running. Each must
be revisited once `~/usatlasdata/dimuon_data/pbpb_2026/` is populated.

| # | Item | Placeholder value used | How to confirm |
|---|------|------------------------|----------------|
| P1 | Number of 2026 NTUP part files (`file_batch_max{26}`, Condor `queue N`, `QUEUE_COUNTS[26]`, `ScrambGen::NParts(26)`, every plotting part list) | **5** — the skim is submitted as 5 grid tasks (`SkimCode/run_26hi/InDstxt_PbPb2026_5p36TeV_part1..5.txt`, 45 datasets over the 35 GRL runs). NOT a copy of 2025's 6. Can end up LARGER: `grid_monitor`'s chunked-hadd fallback splits an over-large task into extra parts. | `pipelines/preflight_pbpb_year.sh 26` checks **all five declaration classes** against each other and against `ls` (extended 2026-09-11; it previously covered only three), treating an UNDER-count as an error and an OVER-count as a warning. **`grid_monitor` auto-updates only `PbPbExtras.c` and `run_pbpb_26.sub`** — the other five `.sub`, `ScrambGen::NParts` and the six plotting lists are manual, so run the preflight after the skim lands. |
| P2 | Total recorded 2026 events | *(unknown)* | entry count of the hadded NTUPs / `data-merging-record.txt` |
| P3 | 2026 bad-run list | **empty** (no run excluded) — `PbPbBadRuns` has no 26 entry, so the R_AA luminosity is the full GRL total 2.62316 nb⁻¹ | 2026 DQ review. If it becomes non-empty, subtract those runs' `Prescale Corrected` in `PbPbSampledLumi.h`, `make_crossx_factors_pbpb_2026()` AND the lumi README **in the same change** — numerator and denominator must cover the same runs |
| P4 | What calibration is the 2026 skim's `centrality` branch on? | **Irrelevant to the result** — under D6 the 2023 calibration is enforced by `UpdateCentrality` regardless. (The old entry asked whether the branch is "zero-filled"; that premise was false even for 2025 — see §3c.) | Open `data_pbpb26_part1.root` and check whether the branch is on a **non-2023** calibration. If it is, the override is intentional under D6 and **must be disclosed in the note**, not silently applied. |
| P5 | Does 2026 need per-run µ+7σ ZDC preamp cuts (2025-style) or a hard scalar (23/24-style)? | **assumed per-run**, i.e. `UsesPerRunPreampCuts(26) = true` in the cuts macro. In the *NTuple processing* this is no longer a guess at all: the loader is now **presence-driven** — it uses per-run cuts iff the year's cuts file contains `t_preamp_per_run` | inspect the 2026 preamp distributions as was done for 2025 |
| P6 | 2026 hard scalar preamp cut (A, C) [ADC] | **{850, 700}** — copied from 2025 (23: {300,300}, 24: {420,420}) | re-derive from the 2026 1D preamp Gaussian-to-tail turnover |
| P7 | 2026 ZDC background-band search window | **{170, 320}** — copied from 2025 | read off the 2026 ZDC-vs-FCal 2D |
| P8 | 2026 background-band µ(x) quadratic coefficients | copied from 2025 | 3 calibration points on the 2026 2D |
| P9 | 2026 alt-banana calibration points | **{3.4, 200, 4.8, 152}** — copied from 2025 | re-derive on 2026 |
| P10 | 2026 event-selection cuts file `event_sel_cuts_pbpb_2026.root` | *(does not exist)* — **hard prerequisite**: `PbPbExtras::InitEventSel` throws without it, so no 2026 NTuple job can run | produce it with `plot_pbpb_event_sel_event_level.cxx(26)` then `plot_pbpb_event_sel_cuts.cxx(26)` |
| P11 | A 2026 minimum-bias sample for the data-driven trigger efficiency | *(unknown)* — `TrigEffPlotterPbPb` needs `histograms_real_pairs_pbpb_2026_MB.root` | as for the other years |
| P12 | 2026 ⟨T_AA⟩ | **2023 Glauber values** — the existing convention, identical to what 2024 and 2025 already do; flagged in `placeholder.md` and required to be disclosed in the note | official 2026 Glauber calibration |
| P13 | 2026 FCal→centrality thresholds | **PbPb2023 thresholds — REGISTERED USER DECISION 2026-09-10 (D6), no longer a guess.** Same as 2024 and 2025 (whose vectors are byte-identical). Interim, not final. | an official 2026 Glauber centrality calibration; until then nothing to confirm |

**Every placeholder is labelled as such in the code**, with a comment pointing back at this
doc, so `grep -rn "PLACEHOLDER" ` over the 2026 sites enumerates them.

## Design Decisions

**D1 — A year switch either has a real 2026 branch or it throws; no silent fallback.**
Rationale: four independent sweeps found that every dangerous 2026 site failed *silently*
(centrality 0, luminosity weight 0, crossx factors −1, invented cut constants), producing
complete and plausible plots. A loud failure costs one rerun; a silent one costs a wrong
result that a plot review cannot catch.

**D2 — Where a year gate was standing in for a property of the data, test the property
instead.** The ZDC per-run preamp cut is the case: `PbPbExtras::InitEventSel` now loads
`t_preamp_per_run` iff the year's cuts file contains it, rather than asking `year == 25`.
Verified no-op for 2023/2024 (their cuts files have no such key) and identical for 2025.
This removes the guess for 2026 entirely — the code follows whatever the cuts file says.

**D3 — 2026 inherits 2025's treatment exactly, and every inherited constant is labelled a
PLACEHOLDER.** Per the sibling doc's Physics Procedure the `hi2026` skim is procedurally
identical to 2023/24/25, so the analysis must treat 2026 as it treats 2025. Constants that
are *measured* per year (preamp cuts, background band, alt-banana points) are seeded from
2025 and marked for re-derivation; constants that are *conventions* (⟨T_AA⟩ = 2023 Glauber,
σ_PbPb = 7.8 b) are reused unchanged because deviating would itself be the anomaly.

**D6 — Pb+Pb 2026 uses the PbPb2023 FCal→centrality calibration (USER DECISION 2026-09-10).**
Registered on the user's instruction: "for now, we will use the pbpb2023 fcal for pbpb26, just
like for pbpb24 & 25". This does not change any number — `UpdateCentrality` `case 26:` already
resolved to `GetCentralityPbPb2023` — it changes the choice's **status**, from an agent
placeholder awaiting confirmation to a registered decision. Verified while registering it that
the premise holds exactly: `FCal_ET_Bins_PbPb2024` is **byte-identical** to
`FCal_ET_Bins_PbPb2023` (all 85 entries), so although `case 24:` calls
`GetCentralityPbPb2024`, all four Run-3 years sit on the same 2023 calibration today. Flagged
at the 2024 vector that editing it alone would silently desynchronise 2024 from the others.
Still **interim**: revisit when an official 2026 Glauber calibration exists. ⟨T_AA⟩ remains on
the same 2023 placeholder, consistent with 2024 and 2025.

**D5 — For year INPUTS the rule is probe-skip-relabel; only year CONSTANTS throw.**
*Caveat recorded 2026-09-11:* where the probe targets an **earlier stage's output**
(`run_all_crossx.sh`, `run_scrambgen.sh`, `RAA_plotting` mode 6), a skip covers both "not
produced yet" and "that stage failed or its output was deleted". It is mitigated rather than
silent — the combined directory name, the R_AA filename suffix and every legend are rebuilt
from the surviving years, so a shrunk combination changes the output path and is visible. Where
the distinction is decidable it IS made: see `plot_forward_qeta_edge_scan`, which throws when a
period has data but is missing the requested working-point variant.
Added after the round-1 review. A missing luminosity, crossx factor or centrality mapping is
a configuration error and must throw (D1). A missing *data file* is a normal transient state
while a year is being produced, and must not take down the years that are ready — nor be
silently claimed in a caption. Probe, skip with an `[INFO]` line, and rebuild every label and
output path from the years that actually contributed.

**D4 — Part count comes from the grid partition, not from 2025.** `SkimCode/run_26hi/
InDstxt_PbPb2026_5p36TeV_part1..5.txt` gives 5, not the 6 that copying 2025 would have given.
`grid_monitor`'s auto-update for `pbpb_2026` was enabled so the value self-corrects from
disk (it had been deliberately disabled while the 2026 analysis code did not exist).

## Implementation Plan

| # | Step | Status |
|---|------|--------|
| 1 | Tracking doc + INDEX registration | DONE |
| 2 | Lumi README: register `pbpb_2026`, per-year table, GRL cross-check | DONE |
| 3 | `analysis_metadata.md` + `analysis_status_summary.md` 2026 update | DONE |
| 4 | Exhaustive enumeration of year-keyed code sites (4 independent sweeps, cross-checked) | DONE |
| 5 | Implement 2026 in the NTuple-processing stage + shared headers | DONE |
| 6 | Implement 2026 in the RDF hist-filling stage + efficiency lookups | DONE |
| 7 | Implement 2026 in pipelines + Condor run scripts | DONE |
| 8 | Implement 2026 in plotting + R_AA combination | DONE |
| 9 | Compile + pre-flight every Pb+Pb pipeline for year 26 | DONE (`pipelines/preflight_pbpb_year.sh`) |
| 10 | `/review-analysis-code` on the C++/RDF changes | round 1 FAILED -> all findings fixed + executor-verified; **independent round-2 review still OWED** (reviewer hit an API session rate limit) |
| 11 | Commit; update INDEX; final summary | pending |

## Progress Log

### 2026-09-10 — Step 2 DONE (2026 luminosity registered)

`IntNotes/data/luminosity/pbpb_2026/lumitable_pbpb_26_HLT_mu4.csv` (user-supplied)
inspected and registered.

- **35 runs, 522041 … 523437.** Verified **byte-identical** to the run list in the skim
  GRL `SkimCode/xmls/physics_HI2026_50ns_noIBL.xml` (`diff` of the two sorted run lists
  → identical). Lumi denominator and skim numerator therefore cover the same runs.
- **Unit = µb⁻¹**, consistent with every other Pb+Pb table (2722.21 "delivered" is
  2.72 nb⁻¹, the right order for a Pb+Pb year; nb⁻¹ or pb⁻¹ would be unphysical).
- Totals: `LDelivered` 2722.21, `LRecorded` 2627.51, `LAr Corrected` 2627.51,
  **`Prescale Corrected` 2623.16 µb⁻¹ = 2.62316 nb⁻¹**; live fraction 96.52 %,
  LAr fraction 100.00 %, prescale fraction 100.17 % → average prescale 1.0017.
- Only run **522658** is prescaled (191 bad LBs, prescale fraction 104.80 %). Every
  other run: 0 bad LBs, prescale 1.000.
- New Run-3 Pb+Pb totals: GRL sum 1.18365 + 0.85112 + 2.59933 + 2.62316 =
  **7.25726 nb⁻¹**; R_AA sum (2023 minus its two b-hadron runs) = **7.24937 nb⁻¹**,
  ≈ 3.7× Run 2, ≈ 1.57× the current 23+24+25 combination.

`README.md` updated: `pbpb_2026` added to the dataset/trigger table and the units
table; a new **per-year totals table** (GRL total, event-level exclusions, R_AA
luminosity, Σ 2023–26); a new **pbpb_2026 section** with the run range, GRL
cross-check, column totals and the prescaled-run note; and an explicit warning that the
2026 R_AA luminosity assumes no 2026 run is excluded at event level (registry item P3).

### 2026-09-10 — Step 3 DONE (status + metadata docs)

- `IntNotes/analysis_metadata.md`: the two "placeholder (data not yet in skim)" 2026 rows
  replaced with real values — §1 lumi table (2623.16 µb⁻¹, prescale 1.0017, 2.62316 nb⁻¹),
  §4 GRL (`physics_HI2026_50ns_noIBL.xml`), §5 HLT chains (a 2026 row, byte-identical to
  2025 by design, `TrigRates_CA.py:190-196`) and the prescale line. New §"Pb+Pb 2026"
  subsection; §3 T_AA placeholder note extended to 2026; code-sync status records that the
  value was not yet wired into the analysis code at the time of writing (it is now).
- `analysis_status_summary.md`: 2026 added to the skim reference, plus a **PbPb 2026 status
  table** (what is DONE / IN PROGRESS / BLOCKED and on what). Also added an explicit staleness
  banner: the rest of that doc has not been revised since 2026-06-22 and predates the
  muon-pT 4→4.5 GeV / gap / pair-pT 9 GeV adoption, so its "DONE" marks describe the OLD
  selection.

### 2026-09-10 — Step 4 DONE (enumeration, 4 independent sweeps)

Four **read-only** subagents enumerated the year-keyed sites: three partitioned by directory
(NTuple stage + shared headers; RDF + pipelines + efficiency; plotting + event selection +
R_AA) and one independent whole-repo sweep partitioned by *pattern class* as a cross-check.
Their scratch docs (`_sub_pbpb26_code_1..4.md`) were merged here and deleted.

**Cross-check outcome:** the three partitioned sweeps agreed site-for-site on the shared
headers, and the pattern sweep found nothing the partitioned ones had missed. The pattern
sweep DID contribute one thing none of the others had: the real 2026 grid partition.

**Corrections the cross-check produced:**
1. **Part count is 5, not 6.** The first three sweeps proposed "guess 6, like 2025". The
   pattern sweep found `SkimCode/run_26hi/InDstxt_PbPb2026_5p36TeV_part1..5.txt` — the skim
   is submitted as **5 tasks** (45 datasets, 35 GRL runs). Verified directly. All sites were
   corrected from 6 to 5, and both writer subagents were messaged mid-task.
2. **`fcal_scale` does not exist in the code.** All four sweeps independently confirmed that
   `fcal_scale` / `fcal_corr_weight` appear in documentation only, never in a source file
   (the FCal reweighting was removed in `40e67c9`). The Physics Procedure §3c of this doc was
   written on the opposite assumption and has been corrected; so have
   `analysis_overview.md` §(e) and `systematic_uncertainties.md` §6, both of which claimed a
   cross-year FCal scaling that is not implemented. The stale auto-memory index line was fixed.

**The four SILENT failure modes found (each would have produced complete, plausible, wrong
plots with no crash and no warning):**

| Site | Unknown-year behaviour |
|---|---|
| `MuonPairPbPb.h` `PairValueCalcHook` | 2026 keeps the skim's all-zero `centrality` branch → **every pair classified 0-1 %, most central** |
| `PbPbSampledLumi.h` `default: return 0.` | 2026 enters the luminosity-weighted year combination at **weight zero** (silently excluded) while the legend still claims it; and `1.0/L` in the RDF becomes **+inf** |
| `PbPbBaseClass.h` `SetCrossxFactorsPbPbCtrBinned` | unregistered year → **crossx factors all −1** (negative cross sections); `SanityCheckPbPb` only checks the vector's SIZE, so it passed |
| `RAA_plotting.cxx` `pbpb_lumis[iy] : 1.0` | a 4th input file with a 3-entry lumi vector → that year silently weighted **1.0 nb⁻¹** |

Plus, in the event-selection derivation: `GetPreampCuts` → an invented **(385, 385) ADC**,
`GetBgSearchRange` / `GetBgMuGuess` / `GetAltCutPts` → **2024's constants**, `FilesForYear`
→ an **empty TChain**, and `ScrambGen::NParts` → **2025's part count**.

**Every one of these now throws instead.**

### 2026-09-10 — Steps 5-8 DONE (implementation)

Design rule applied throughout: *a year switch either has a real 2026 branch or it throws.*
Where the year gate was standing in for a property of the data, it was replaced by a test of
that property, which removes the guess entirely:
- **ZDC per-run preamp cuts (NTuple stage) are now presence-driven, not year-gated.**
  `PbPbExtras::InitEventSel` loads `t_preamp_per_run` whenever the year's cuts file contains
  it. Verified that only the 2025 cuts file has that key (2023 and 2024 do not), so 23/24
  keep the scalar cut byte-for-byte and 2026 picks up per-run cuts iff its cuts file was
  derived with them.
- **`plot_npairs_vs_centrality.cxx` was de-magicked**: it used `TH1D* h_os[4]` where index 3
  was *simultaneously* the array bound and the "combined" slot. Now `kNY = years.size()`,
  `kComb = kNY`, `std::vector<TH1D*>(kNY+1)`; no index literal survives.
- **The three `plot_zdc_preamp_*` macros** now size every array, loop, `TGraph`, axis limit,
  `SetNdivisions`, `TLine`, label and `min/max_element` from `kYearsPP[]`/`kNYears`, so a
  3-year figure under a 4-year legend is structurally impossible.

Files changed: see the commit. New files: six `run_pbpb_26*.{sh,sub}`, three
`run_*_pbpb26.sh` RDF runners, and `pipelines/preflight_pbpb_year.sh`.

**A pre-existing bug fixed in passing (changes an existing 2025 figure):**
`plot_pbpb_event_sel_event_level.cxx::MakeZDCTimeCentralityPlot` read the raw `centrality`
branch with **no year gate at all**, so for 2025 (whose branch is zero-filled) every event
landed in the 0-1 % panel. It now recomputes from FCal E_T for 25/26 and captions the panel
accordingly. `ZDC_time_AC_corr_top_5_centr_pbpb_2025.png` will therefore change when
regenerated. **No cut is derived or saved by that function**, so nothing downstream moves.

### 2026-09-10 — Step 9 DONE (compile + pre-flight)

- ACLiC/`g++` compile, 0 errors: `DataAnalysisClasses.h` (NTuple stage),
  `RDFBasedHistFillingPbPb.cxx`, `RAA_plotting.cxx`, `plot_single_b_crossx_pbpb.cxx`,
  `plot_crossx_reco_eff_stages.C`, `plot_crossx_trig_corr_sanity.C`, and all six
  event-selection macros. Only pre-existing `-Wsign-compare` warnings from `ParamsSet.h`.
  (`plot_dR_trig_corr.C` does not build under ACLiC and did not before — it is an
  interpreted macro missing `#include <iostream>`; verified by interpreted `.L`.)
- **Behavioural check on year 26 through the real NTuple-processing entry point:**
  `PbPbAnalysis(26, 1)` is accepted, resolves its paths to `pbpb_2026/`, and fails **loudly**
  at the expected point (missing `data_pbpb26_part1.root`). `PbPbAnalysis(26, 9)` throws on
  the file-batch range; `PbPbAnalysis(27, 1)` throws on the year whitelist.
- **New `pipelines/preflight_pbpb_year.sh <yr>`**: checks that every repo artifact the Pb+Pb
  pipelines reference for that year exists, and cross-checks the three declarations of the
  part count (`file_batch_max`, `queue N`, `QUEUE_COUNTS`) against each other and against the
  files on disk — the one mismatch class that silently processes only part of the data.
  Result: **26 passes** (15 repo artifacts present, 9 data artifacts pending as expected,
  0 count mismatches), and **23 / 24 / 25 also pass against their real on-disk data**, so the
  check is not vacuously green.

### 2026-09-10 — Step 10 round 1: `/review-analysis-code` returned **FAIL**, findings fixed

**Root cause (accepted):** design decision D1 ("a year switch either has a real 2026 branch
or it throws") was correct for year **constants** but was applied indiscriminately to year
**inputs**. Because the 2026 data does not exist yet, that made seven existing 2023/24/25
workflows unproducible, and left three more silently overwriting correct figures with a
4-year caption over 3-year data (`TChain::Add` only warns on a missing file).

**D5 (new) — for year INPUTS the rule is probe-skip-relabel, not throw.** The pattern already
existed in `plot_single_b_crossx_pbpb.cxx`: probe each candidate with
`gSystem->AccessPathName`, skip an absent year with an `[INFO]` line, and rebuild **both** the
labels and the output path from the years that actually contributed — so a figure can never
claim a year it did not use, and a year still in production cannot take down the years that
are ready. Applied uniformly. A year *constant* that is missing is still a configuration
error and still throws (D1 stands for those).

Fixed, each verified:
- **`RAA_plotting` mode 6** now discovers its year set. Verified end-to-end with `pbpb_2026/`
  absent: logs the skip, reports "combining 3 Pb+Pb year(s)", regenerates the same three
  filenames. The suffix is built with no separator before the first year so the 23/24/25 set
  reproduces `_pbpb23_24_25_combined_pp24_2mu4` exactly — otherwise every existing R_AA PNG
  would have been orphaned.
- **Run-2 regression (the subtlest finding).** `PbPbMu4SampledLumiNb`'s new throwing default
  broke years 15/18, which are *not* unknown: `RDFBasedHistFillingPbPb` still carries a full
  Run-2 branch and `PbPbBaseClass` still registers `{18,"default"}`. Added
  `case 15/18 = 1.817182 nb⁻¹`, verified at runtime to equal exactly the luminosity already
  baked into `make_crossx_factors_pbpb_run2()` (0.436882 + 1.3803). The throw message, which
  the first pass had rewritten to "23/24/25/26 (Run 3 Pb+Pb)", was restored to include 15/18.
- **Four pipeline drivers + `run_scrambgen.sh` + `run_all_crossx.sh`** filter the year list on
  input presence and announce every skip. Two bugs in the first attempt at that filter, both
  caught by testing rather than assumption: `exit` inside `<(...)` only ends the subshell, and
  `printf '%s\n'` with **zero** arguments still prints one blank line, which `mapfile` turns
  into a one-element array containing `""` so the caller's emptiness test passed with nothing
  to do. Verified: default → `23 24 25` with a loud skip; `YEARS="26"` → FATAL, exit 1.
- **ScrambGen** would have written an *empty* scrambled file, which passes the downstream
  existence check and feeds an empty `T_mix` into the template fit. `LoadMuons` now throws on
  zero entries and `GeneratePairs` refuses to create the output when there is nothing to mix.
- **`run_all_crossx.sh`** had been missed by all four enumeration sweeps: no 2026 filling step,
  and a hard-coded `pbpb_23_24_25_combined` validation path while the plotter derives that name
  from the years it finds. Both now come from one discovered list (verified: still resolves to
  `pbpb_23_24_25_combined` today).
- **Six multi-year plotting macros** and **the six event-selection macros** converted to
  probe-skip-relabel, with the silent-empty traps closed (empty chain, non-null-but-empty
  histogram, unchecked `TFile::Get`, single-part "availability").
- **The Glauber mirror was single-sourced with ZERO numerical change.** The two hand-typed
  copies now include one `plotting_codes/event_selection/PbPbCentralityFCalMirror.h` whose
  values are byte-identical to the pre-edit macro text (machine-verified against `81ad7ef`).
  Two copies became one; no number moved; the reconcile-with-canonical decision stays parked.

**A pre-existing plotting bug found by actually running a macro** (not by reading it):
`plot_npairs_vs_centrality` pad 4 ("OS / SS combined") rendered EMPTY on a 0-50000 axis, and
pad 3's linear panel was squashed into the bottom sixth of its frame. Both pads clone
`h_os[kComb]` *after* `draw_log()` has called `SetMaximum(5x peak)` on it. The obvious repair
fails, and instructively: **`TH1::GetMaximum()` returns the STORED `fMaximum` once
`SetMaximum` has been called**, so `GetMaximum()*1.2` just returns the same 5e4. Fixed by
taking the largest bin content. Pad 3 now peaks at ~10700 pairs near 7 % centrality and OS/SS
rises from ~1.0 central to ~2.5 peripheral — the expected shape (combinatorial-dominated
centrally, correlated fraction growing peripherally). Pre-existing at `81ad7ef`, identical
code with the literal index `[3]`.

### 2026-09-11 — Step 10 round 2: independent review returned FAIL; findings fixed

Round 2 ran after the API session limit reset. **3 WARNING + 8 INFO, zero CRITICAL** — all
eight round-1 CRITICALs verified closed. Every reported number re-derived independently and
MATCHED (2026 luminosity, GRL run-list identity, crossx factor ratio 0.990915536986 in all six
bins, Run-3 sums, Run-2 1.817182, the 5-task partition, the Glauber mirror byte-identity, and
`FCal_ET_Bins_PbPb2024` ≡ `_PbPb2023`).

**W1 — the stated premise for the centrality recompute was FALSE (the important one).** The
reviewer read the real NTUPs; I reproduced it independently. See the correction now in §3c:
the 2025 `centrality` branch is populated (mean 20.454, max 84, 4.08 % zeros) and already on
the 2023 calibration (299 992/300 000 agreement with the 2023-table recompute, the 8 exceptions
off by one unit and traceable to the rounded mirror table). The recompute is a **no-op** for
2025, and the 2026 branch rested on a premise that does not hold. **Behaviour unchanged and
still correct under D6** — it is an *override* enforcing the 2023 calibration, not a repair.
Corrected in `MuonPairPbPb.h`, both event-selection macros (including the runtime log line and
the in-loop comments), §3c and registry P4. P4's confirm-step was the wrong test and is
rewritten: ask whether the 2026 branch is on a *non-2023* calibration, and if so disclose the
deliberate override in the note.

**W2 — `plot_forward_qeta_edge_scan` Medium WP silently shrank to one year.** `_medium_wp`
histograms exist only for 2023 (verified on disk); the new probe-skip turned a loud failure
into a 2023-only Medium figure sitting beside a three-period Tight figure — a mismatched pair
for the working-point systematic, presented as finished. `AvailablePbPbYears` now distinguishes
**(a)** a period with no data at all (skip, D5) from **(b)** a period that *has* data but is
missing this WP variant (**throw**, naming the file to produce). This is the "not produced yet
vs genuinely failed" distinction D5 did not make.

**W3 — the preflight overstated its own guard.** Registry P1 claimed it cross-checks all
declarations of the part count; it checked three of five. Extended to assert
`ScrambGen::NParts(<yr>)` and the per-part file lists in the six event-selection/preamp/FCal
macros, and all `.sub` variants are now checked (the deliberately-partial one-offs are named
explicitly, so a NEW variant is checked by default rather than silently exempt). It now also
distinguishes the two directions, which are not equally dangerous: an **under**-count silently
processes a subset (error), an **over**-count throws loudly on the missing part (warning).

**The extended preflight immediately found real pre-existing drift** (all pre-dating this work,
none introduced by it, none fixed unilaterally):
- The three `plot_zdc_preamp_*` macros list **3 of 2023's 4** part files, so those cross-year
  diagnostic figures are built from a subset of 2023. (The *cut derivation* macros correctly
  list 4 — verified — so no derived cut is affected.)
- `run_pbpb_23_no_res_mu4_mu4noL1.sub` (queue 6 vs 4), `run_pbpb_24_no_res_mu4_mu4noL1.sub` and
  `run_pbpb_24_no_res_mu4.sub` (queue 9 vs 2) are stale. Over-counts, so they fail loudly; and
  `mu4_mu4noL1` is not maintained.

INFO items fixed: the legacy `case 26:` carries 25's TODO; `plot_dR_trig_corr`'s pT-slice header
is derived from `kPtSliceYear` instead of a literal; the two crossx sanity macros probe with
`AccessPathName` before `TFile::Open` so an unproduced year no longer prints a raw ROOT error;
the wrong `// [0]=A, [1]=C` comment in `PbPbExtras.c:99` is flipped to the project convention.

**Declared exceptions to "nothing changed for existing years" — now three, not two:**
(a) `MakeZDCTimeCentralityPlot` (2025 centrality panel; derives no cut, nothing consumes it);
(b) `plot_npairs_vs_centrality` pads 3/4 (axis ranges only);
(c) **cosmetic year relabelling `"PbPb 23+24+25"` → `"PbPb 2023+2024+2025"`** in five macros,
required by the no-typed-year-string rule; filenames unchanged.
Also noted: regenerating the 2025 cuts file changes a `TTree` *title* string (`(PbPb25)` →
`(PbPb2025)`) and transient histogram names; the key name and every cut value are unchanged, so
a byte-diff of a regenerated file is not a physics change.

## Results & Observations

### Open questions for the user (none blocked the work; all recorded)

1. **`plot_dR_trig_corr.C:196`** reads `files[25]` by hard-coded map key for the pT-slice
   plot set, with a "PbPb 2025, " header. This is a deliberate single-year choice, so it was
   **not** silently repointed at 2026. Given the standing "Pb+Pb always combined" preference,
   it arguably should sum over years instead of being pinned to one.
2. **`plot_pbpb_fcal_comparison.cxx`**: now that `LoadCutsFC` throws on missing inputs, an
   unconditional 2026 column would make the existing 23/24/25 figure unregenerable until the
   2026 skim lands. It currently warns and prints "PbPb 2026: not available" in that pad
   instead of throwing.
3. **`TrigEffPlotterPbPb.cxx:560`** draws a hard-coded `"Pb+Pb 2023, "` pad label — a
   pre-existing bug that already mislabels the 2024 and 2025 output. Untouched.

### Pre-existing defect found, NOT fixed (needs a decision — it touches existing years)

**The Glauber FCal→centrality table is duplicated, and the copies are rounded.**
`PairPbPbExtras::FCal_ET_Bins_PbPb2023` (`MuonPairPbPb.h`) is the canonical 85-entry table.
`plot_pbpb_event_sel_cuts.cxx` carries a hand-typed mirror `kFCalBinsPbPb2023[85]` (this
predates the 2026 work — it is in `HEAD`), and the 2026 work propagated that same mirror into
`plot_pbpb_event_sel_event_level.cxx`. Checked numerically: **40 of the 85 entries differ
from the canonical values**, all by rounding to 5 significant figures — e.g. index 37
`0.965176` → `0.96518`, index 79 (the 80 % boundary) `0.063208` → `0.06321`, worst relative
difference 8.2e-5 at index 82.

The two mirrors agree with each other exactly, so the event-selection derivation is
self-consistent across years; the divergence is against the table the *analysis* uses. The
practical effect is confined to events sitting within ~1e-5 TeV of a threshold, i.e.
statistically negligible — but it is exactly the silent binning divergence the project's
BLOCKING binning rule exists to prevent, and index 79 is the 0-80 % acceptance boundary used
by `is_ctr80` in the cut derivation.

**Not fixed here**, because the correct fix (make both macros read the canonical table
instead of mirroring it) perturbs the event-selection derivation for 2023/2024/2025, which is
outside the requested 2026 scope and is a user decision.

## Remaining Work

- **Step 10:** `/review-analysis-code` on the C++/RDF changes and `/review-pipeline` on the
  pipeline changes (per the project's per-step protocol).
- **After the skim lands:** confirm every Placeholder Registry item P1–P12; run
  `pipelines/preflight_pbpb_year.sh 26` first — it will catch a part-count drift immediately.
- **Then, in order:** derive `event_sel_cuts_pbpb_2026.root` (hard prerequisite) → NTuple
  processing → hadd → trigger efficiency (needs a 2026 MB sample) → RDF crossx →
  the combined 23+24+25+26 crossx and R_AA plots.
- **User decisions** listed under Results & Observations: the `plot_dR_trig_corr.C` single-year
  pT-slice set, the FCal-comparison missing-year behaviour, and the duplicated/rounded Glauber
  threshold table.
- `IntNotes/tex/datasets.tex` needs 2026 (prose "2023, 2024, and 2025", the recorded-event
  table, the luminosity table, the GRL table). Deferred: it is internal-note prose governed by
  the academic-writing gate chain (`/review-note`), and its event-count row cannot be filled
  until the skim finishes.

## Latest Stage

**Code is 2026-ready and every Pb+Pb workflow runs today with the 2026 data still absent.**
Steps 1-9 DONE. Step 10 round 1 returned FAIL; all findings are fixed, committed and verified
by execution (year filter run from the real script under four configs; 13 macros compile clean;
`RAA_plotting` and `plot_npairs_vs_centrality` run end-to-end and their figures inspected;
per-year constants diffed against `81ad7ef` and found byte-identical for 23/24/25).

**OWED: an independent round-2 `/review-analysis-code` pass.** The round-2 reviewer subagent
terminated early on an API session rate limit, so round 2 was verified by the executor, not by
an independent reviewer. The loop is NOT closed with a PASS — rerun it before treating this
work as reviewed.

**Then:** wait for the 2026 skim, run `pipelines/preflight_pbpb_year.sh 26` first (it will catch
a part-count drift immediately), and work the Placeholder Registry P1-P12.

**Open user decisions** (all recorded above, none blocking): the Glauber mirror reconcile, the
naive year sum in `plot_crossx_trig_corr_sanity.C`, the scalar-vs-per-run preamp cut drawn in
the preamp/FCal figures, the single-year `files[25]` pT-slice set in `plot_dR_trig_corr.C`, and
the FCal-comparison canvas widening at unchanged filenames.
