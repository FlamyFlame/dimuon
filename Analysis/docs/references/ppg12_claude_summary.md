# `ppg12/.claude` — HEP/HI Reference Summary and Cherry-Pick

**Source:** https://github.com/Shuonli/ppg12 (`.claude/`) — Shuon Li, sPHENIX PPG12
isolated-photon cross-section analysis. A real, HEP/HI-specific Claude Code setup
on a working collaboration analysis. Studied from a local clone.

Companion to `academic_research_skills_summary.md`; both feed
`../academic_writing_workflow.md`. ARS supplies generic *rigor protocols*; ppg12
supplies *HEP-shaped concrete mechanisms* that map cleanly onto our codebase.

---

## 1. What ppg12 ships (`.claude/`)

- **Commands:** `sync-check` (constants/filenames/features/bins in sync across
  pipeline codes), `check-results` (systematic-variation outputs complete/healthy),
  `compile-note` (pdflatex+bibtex build + log grep), `make-report`, `session-report`,
  `make-plots`, `run-systematics`, `run-efficiency`, `inspect-root`, `validate-config`,
  `deploy-pages` (publish reports+figures to GitHub Pages), `ask-shuonli` (advisor).
- **Agents:** `note-critic` (read-only LaTeX note critique vs current configs),
  `numerical-rederivation` (independently re-extract a cited number from ROOT, never
  trusting the producing code), `physics-reviewer`, `report-writer`, `code-writer`,
  `data-explorer`, `fix-validator`, `plot-cosmetics-reviewer`, `prompt-refiner`.
- **Rules:** `root-macros.md`, `python-analysis.md`, `yaml-config.md`, `wiki-usage.md`,
  `double-interaction.md`.
- **Pattern:** several plotting macros **write figures directly into the note's
  `Figures/` directory**, so the note's figures are always the latest plot output.

---

## 2. Relevance: high and direct

ppg12 is the same kind of object we are building — a HEP collaboration analysis with
a LaTeX note (sPHENIX `sphenix` class / our ATLAS `atlasdoc`), ROOT outputs,
systematic variations, and a need to keep prose synchronised with code. Its agents
already encode HEP-correct verification discipline (ROOT binning conventions, "never
trust the producing code, re-extract from the output file", config-vs-prose
accuracy). It needs **no adaptation of philosophy**, only renaming to our pipeline.

---

## 3. Cherry-pick map (ppg12 → our implementation)

| ppg12 element | Ported into | Notes |
|---|---|---|
| `sync-check` (values in sync across codes) | `commands/check-note-sync.md` | reframed: note prose ↔ current analysis ROOT/table outputs + figure staleness + placeholder inventory |
| `numerical-rederivation` agent | enhanced `agents/numerical-verifier.md` | re-extract from ROOT with own minimal script; per-quantity tolerances; `CANNOT LOCATE` not guess |
| `note-critic` checklist | enhanced `agents/note-reviewer.md` | config/number accuracy vs source; `defs.sty` macro consistency; sentence-length; vague/hedge language; figure-data fidelity |
| `compile-note` | `commands/compile-note.md` | wraps `IntNotes/Makefile` (latexmk/biber, not bibtex); greps undefined refs / missing figures / overfull boxes |
| figures-written-into-note-dir | `commands/sync-note-figures.md` | we use copy+provenance manifest (CI-safe) instead of macros writing in-place, because plots live outside the submodule and the note compiles remotely |
| `deploy-pages` | noted future (thesis/results website) | out of scope now |
| `session-report` / `make-report` | partially covered by our tracking docs + `wrap-up` | not re-implemented |

What we **do not** adopt: ppg12's analysis-execution commands (`run-efficiency`,
`train-bdt`, `apply-bdt`, `run-systematics`) — they are PPG12-pipeline-specific; our
equivalents are the existing `/steer-pipeline` and `/review-pipeline`.

---

## 4. Key divergence from ppg12: figure sync

ppg12 can let plotting macros write straight into `PPG12-analysis-note/Figures/`
because the note and the analysis live in **one repo**. Our `IntNotes/` is a
**separate CERN GitLab submodule compiled remotely**, and our plots live outside it
(in `Analysis/plots/` and the data area, untracked). So a write-in-place or symlink
approach would dangle on the CI runner. We therefore use a **copy + provenance
manifest + staleness check** (`/sync-note-figures`) — same goal (note figures always
current), CI-safe mechanism. See `../academic_writing_workflow.md` §G4 and the
tracking doc Design Decision D2.
