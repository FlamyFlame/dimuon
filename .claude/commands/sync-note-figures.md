---
description: >
  Sync the internal note's figures with the latest analysis plot outputs, via a
  copy-based provenance manifest (CI-safe), and report stale/missing/orphaned
  figures. Use when the user says "sync the figures", "update the note figures",
  "are the note figures current", or before compiling/freezing the note. Implements
  Gate G4 (sync half) of Analysis/docs/academic_writing_workflow.md.
---

Sync `IntNotes/figures/` with current analysis plots and report freshness.

Arguments: $ARGUMENTS
- (empty) → copy-sync all manifest entries + write provenance + staleness report
- `--dry-run` → report only; copy nothing
- `--symlink` → **local fast-path**: symlink instead of copy (live preview;
  NEVER commit — symlinks dangle on the CERN GitLab CI runner and on standalone
  clones; re-run a plain copy-sync before committing/freezing)
- `--check` → staleness/orphan report only (same as `--dry-run` but no provenance write)

## Why copy, not symlink (do not change without re-reading this)
`IntNotes/` is a standalone CERN GitLab submodule compiled **remotely** by the ATLAS
PO CI, where neither the `Analysis/` repo nor the data area is checked out. `\graphicspath`
needs figures **physically present** inside the submodule. A committed symlink to
`../Analysis/plots/...` (or the data area) dangles there and on any collaborator clone —
the same trap already documented for `IntNotes/placeholder.md`. So the committed figures
are **real copies**, with a provenance record proving they are the latest output.

## The manifest
`IntNotes/figures/figure_manifest.yaml` — one entry per note figure:
```yaml
# note_name: filename used in \includegraphics (relative to figures/)
# source:    absolute path to the current analysis plot output
figures:
  - note_name: reco_eff_placeholder_pp.png
    source: /gpfs/.../Analysis/plots/reco_effcy/placeholder/reco_eff_placeholder_pp.png
    note: "PP reco-eff placeholder (G7: placeholder, see placeholder.md)"
```
If the manifest does not exist, create a stub from the figures currently referenced
by `\includegraphics` in the note and ask the user to fill in each `source`.

## Implementation
The mechanism is implemented in `IntNotes/scripts/sync_note_figures.py` (stdlib
only — no pip, no PyYAML). **Run it** rather than re-deriving the logic:
```
python3 IntNotes/scripts/sync_note_figures.py [--check|--dry-run|--symlink]
```
It copies current sources into `IntNotes/figures/`, writes
`figure_provenance.json` (with the `Analysis` repo git SHA), and prints the status
table + orphan/coverage findings, exiting non-zero on STALE/MISSING-SOURCE/UNTRACKED.
Relay its output and the final `SYNC: CLEAN` / `SYNC: ACTION NEEDED` verdict. The
steps below document what the script does (and what to do if the manifest is absent).

## Procedure (what the script does)
1. `mkdir -p IntNotes/figures`. Read `figure_manifest.yaml`; parse the list of
   `(note_name, source, note)` entries.
2. Collect the set of figures actually used: grep `\includegraphics` across
   `IntNotes/*.tex` and `IntNotes/tex/*.tex` for basenames.
3. For each manifest entry:
   - If `source` missing on disk → **MISSING-SOURCE** (CRITICAL); skip copy.
   - Compute `source` sha256 + mtime; compare to the committed
     `IntNotes/figures/<note_name>` (if present) and to the last provenance record.
   - Status: **UP-TO-DATE** (identical sha), **STALE** (source newer / sha differs),
     **NEW** (no committed copy yet).
   - Unless `--dry-run`/`--check`: copy `source` → `IntNotes/figures/<note_name>`
     (`cp -p`), or symlink it under `--symlink`.
4. **Orphan/coverage checks:**
   - figure used by `\includegraphics` but absent from the manifest → **UNTRACKED** (WARNING).
   - manifest entry whose `note_name` is never `\includegraphics`'d → **UNUSED** (INFO).
   - committed file in `figures/` neither in manifest nor a logo → **ORPHAN** (INFO).
5. Unless `--dry-run`/`--check`: write `IntNotes/figures/figure_provenance.json`:
   ```json
   {"synced_at":"<ISO8601>","analysis_git_sha":"<git -C Analysis rev-parse --short HEAD>",
    "figures":[{"note_name":"...","source":"...","source_mtime":"...","sha256":"...","status":"..."}]}
   ```
   Use a stdlib-only `python3 -c` (hashlib/json/os) snippet — no pip, no pyyaml needed
   (you parse the YAML yourself when reading the manifest).
6. The committed artifacts are the **real copied figures**, `figure_manifest.yaml`,
   and `figure_provenance.json` (do not gitignore these). `--symlink` replaces a
   figure with a symlink **at the same path** — there is no separate extension to
   gitignore, so a symlinked figure must NEVER be committed: under `--symlink` the
   script stamps `"mode":"symlink"` in the provenance and you must re-run a plain
   copy-sync (no flag) before any `git add`/commit/freeze. Do not commit anything yourself.

## Output
A table `| note_name | status | source mtime | committed sha (short) |` plus the
orphan/coverage findings, then:
- **SYNC: CLEAN** — every manifest figure UP-TO-DATE, no MISSING-SOURCE, no UNTRACKED.
- **SYNC: ACTION NEEDED** — list STALE (now refreshed unless dry-run), MISSING-SOURCE,
  UNTRACKED, with the exact next step.

Never edit `.tex`/`.bib`. Only write under `IntNotes/figures/` (+ the `.gitignore` line).
