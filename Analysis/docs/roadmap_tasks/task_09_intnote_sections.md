# Task 09 — Write IntNote sections (ANA-HION-2023-07-INT1)

**Roadmap:** `analysis_roadmap_2026_06.md` §Q1 / chain [10]
**Can start:** immediately. **Reviewer:** `/review-note` (every
section). **Template:** Run 2 note — structure and per-section content
summarized in `.claude/kb/analysis/run2_dimuon_note.md`; exact wording
in `IntNotesRun2DimuonReference/tex/*.tex`.

## Objective

Populate `dimuon_codes/IntNotes/` (currently bare skeleton:
Intro/Datasets/Analysis/Results chapters + `tex/datasets.tex` with 3
JIRA refs) with all sections rated READY or DRAFT+PH in roadmap Q1.

## Section plan (write in this order)

1. **Skeleton restructure:** convert the bare chapters into the Run 2
   note's section flow: Introduction / Datasets (PbPb data, pp data,
   triggers+GRL, MC samples, centrality) / Event selection / Muon &
   pair selections / Efficiency corrections (trigger, reco) / Analysis
   & observables / Background estimate / Systematics / Results /
   Appendices. Create one `tex/<Section>.tex` per section, `\input`
   from the master file (mirror Run 2 conventions; document class is
   already atlasdoc).
2. **Introduction** — adapt Run 2 intro to single-b crossx/R_AA at
   5.36 TeV; cite Run 2 note/paper as the predecessor measurement.
3. **Datasets** — May 2026 skim: PbPb23 (4 parts, 124,473,950 evts),
   PbPb24 (2 parts, 92,650,031), PbPb25 (6 parts, 260,392,022), pp24
   (12 parts). `\todo` placeholders for: dataset container names, GRLs,
   HLT chain names/prescales, luminosities (roadmap Q2 — NEVER fill
   these from guesses). MC: ATLHI-596/666/658 + overlay r-tags with
   production status.
4. **PbPb event selection** — 5-cut sequence with derivation method
   per cut (kb overview §PbPb event selection); figures from
   `plots/single_b_analysis/event_selection/`. Fully writable.
5. **Centrality** — FCal ET, 2023 Glauber thresholds for all years,
   PbPb25 recalculation, FCal cross-year scale factors. Writable;
   T_AA table = `\todo`.
6. **Muon & pair selections** — selections from NTupleProcessingCode
   (quality WPs, pT/η cuts, resonance vetoes V1/V2, photoproduction
   Aco/Asym cuts, signal-region definition); cutflow figures.
7. **Trigger efficiency** — PbPb: Physics Procedure from
   `mu4_trig_effcy_implementation.md` (tag-and-probe dR>0.8, role-swap
   combination, Fermi+log fits, per-centrality); figures from
   `plots/pbpb_trigger_efficiency/`. State the D9 finding (dR
   corrections biased on triggered data; MC-based ε_dR planned). pp24
   2mu4 = placeholder until task 02.
8. **MC truth studies appendix** — Pythia flavor/origin and Powheg
   ancestor-group results (plots exist under `plots/pythia`,
   `plots/powheg`, `plots/flavor_categoried`).
9. **Reco efficiency** — method text (overlay-based, match-prob > 0.5,
   per-centrality) with placeholder figures labeled test-sample/DUMMY.
10. **Cross-sections / R_AA** — method + preliminary figures with
    "uncorrected/preliminary normalization" caveats until task 05/06.
11. **Systematics** — source list table (task 08), all values `\todo`.

## Rules

- Every number lacking confirmed provenance → `\todo{}` with what is
  needed; never copy Run 2 values silently (different energy/datasets).
- Per CLAUDE.md auto-dispatch: each writing step goes through
  `/review-note`.
- Keep a build check (`make` in IntNotes/) green after every section.
- Update IntNotes git (separate repo) with one commit per section.
