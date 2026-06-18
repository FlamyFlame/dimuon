---
description: >
  Read-only sweep checking the internal note is in sync with the current analysis:
  numbers in prose vs ROOT/table outputs, figure freshness, citation hygiene, and
  placeholder honesty. Use when the user says "is the note up to date", "check note
  sync", "does the note match the analysis", or before committing/freezing the note.
  Aggregates Gates G2(existence)+G3(spot)+G4(staleness)+G7 of
  Analysis/docs/academic_writing_workflow.md.
---

A single **read-only** consistency sweep of the internal note against the current
analysis. Modify nothing.

Scope: $ARGUMENTS  (a section, a file, or empty = whole note)

## Setup
Read `Analysis/docs/academic_writing_workflow.md` and `Analysis/docs/placeholder.md`.

## Checks (report each as MATCH / MISMATCH / STALE / MISSING / WARNING)

1. **Numbers ↔ analysis (G3 spot-check).** Extract numeric claims in the in-scope
   prose (luminosities, kinematic ranges, bin edges, cut values, yields,
   efficiencies, cross-sections, R_AA, T_AA). For each, identify the authoritative
   source: `IntNotes/analysis_metadata.md`, the code constants
   (`PbPbBaseClass.h`/`PPBaseClass.h`), or a ROOT output. Spawn `numerical-verifier`
   for values backed by ROOT files; compare metadata-backed values directly. Report
   MATCH/MISMATCH with file:line + source.
2. **Figures fresh (G4).** Run `/sync-note-figures --check`; fold in its
   STALE / MISSING-SOURCE / UNTRACKED findings.
3. **Citations resolve (G2 existence).** Quick pass: every `\cite` key is defined in
   a bib resource; flag UNDEFINED. (For full claim-support, point the user to
   `/verify-citations`.)
4. **Placeholder honesty (G7).** Cross-check every placeholder/preliminary statement
   against `placeholder.md`: each standing placeholder (2024/2025 ⟨T_AA⟩, reco-eff,
   σ_PbPb, etc.) that appears in the note must be explicitly flagged as such; any
   `\todo`/TODO/FIXME left in the source is listed.
5. **Stale-result guard.** Flag any figure whose provenance `analysis_git_sha`
   predates the current `Analysis` HEAD, or any number whose source ROOT file is
   older than the most recent crossx/RAA run referenced in the status summary.

## Output
A grouped report (one block per check) and a final:
- **IN SYNC** — all MATCH/UP-TO-DATE, no MISMATCH/STALE/MISSING/UNDEFINED.
- **OUT OF SYNC** — enumerated discrepancies, each with the exact fix and which
  command to run (`/sync-note-figures`, `/verify-citations`, or a re-run of the
  relevant analysis step).

This is a diagnostic. It never edits the note; it tells the user (or the
`/review-note` loop) what to fix.
