# Academic Writing Production Chain — Implementation Tracking

> **Mode:** Implementation. **Note:** This task is *tooling/workflow*, not physics.
> The "Physics Procedure" slot is replaced by the **Workflow Specification**
> below, which is the authoritative reference for all design/implementation
> choices here. Flagged to the user per CLAUDE.md (no physics procedure exists
> for a meta-workflow task).

## Objective

Establish a single, coherent, academic-grade writing production chain for the
Run 3 dimuon analysis (internal note → journal paper → thesis) that guarantees
(1) **rigor / hallucination prevention** — every statement backed by evidence,
every citation real *and* supporting, captions following from data, prose at
ATLAS/publisher standard — and (2) **auto-sync** — the note/paper always reflect
the latest analysis outputs (numbers + figures). The quality procedures must
**actually fire** on real writing tasks, even a single section, with nothing
missing.

## Workflow Specification (AUTHORITATIVE — replaces Physics Procedure)

### Motivation
LLM-assisted scientific writing fails in characteristic ways: fabricated or
non-supporting citations, numbers/figures drifting from the current analysis,
plausible-but-wrong methodology prose, and sycophantic self-review. The chain
below closes each gap with a named gate, ported from two references
(`academic-research-skills`, `ppg12`) and grafted onto the user's existing
executor↔reviewer infrastructure.

### Top-level chain (every writing task passes the gates that apply)
```
[research/KB] → [draft] → G1 evidence-support → G2 citation existence+support
   → G3 numerical re-derivation → G4 figure-data fidelity + figure-sync
   → G5 reviewer loop (anti-sycophancy) → G6 compile → G7 disclosure/epistemic
```
The chain is enforced through the existing `/review-note` and `/review-paper`
commands (which embed the gates) plus three new read-only commands
(`/verify-citations`, `/check-note-sync`, `/compile-note`) and one sync command
(`/sync-note-figures`). `Analysis/docs/academic_writing_workflow.md` is the
human-readable ground-truth spec; the commands reference it.

### The gates (what each checks; what is and is NOT required)
- **G1 Evidence-support (anti-hallucination).** Every factual/quantitative claim
  in prose traces to either (a) a real citation, or (b) a verifiable
  analysis-output source (ROOT file + histogram/bin, table, or plot). Claims with
  neither are flagged `[UNSUPPORTED]`. Anti-leakage: prefer session/analysis
  materials over model memory; flag `[MATERIAL GAP]` rather than filling from
  memory. *(ARS anti-leakage + claim-verification.)*
- **G2 Citation existence + support.** Every `\cite{key}` resolves in a
  `\addbibresource` `.bib` (biblatex/biber). Each cited reference must actually
  support the adjacent claim (verdict taxonomy: VERIFIED / MINOR / MAJOR
  DISTORTION / UNVERIFIABLE). HEP adaptation: resolve via INSPIRE-HEP / arXiv /
  DOI where egress allows, else flag `UNVERIFIED — needs manual`. NOT required:
  re-verifying the standard ATLAS/CMS curated bib entries' existence (they are
  INSPIRE exports); DO verify the *support* relationship and any
  user-added entries in `ANA-HION-2023-07-INT1.bib`.
- **G3 Numerical re-derivation.** Every number in prose is independently
  re-extracted from the source ROOT file by `numerical-verifier` — never trust
  the prose or the producing code. Tolerances per quantity type; `CANNOT LOCATE`
  rather than guess. *(ppg12 numerical-rederivation.)*
- **G4 Figure-data fidelity + sync.** Each `\includegraphics` figure (a) is
  referenced in text, (b) has a caption that follows from the data and states
  dataset+selection, (c) is the *current* analysis output via the figure
  provenance manifest (staleness check). *(ARS figure-table fidelity + ppg12
  figures-in-note-dir.)*
- **G5 Reviewer loop (anti-sycophancy).** The existing executor↔reviewer loop,
  with the reviewer doing a paper-blind pre-commitment of expected
  criteria before reading the draft (generator-evaluator discipline), so it does
  not rubber-stamp. PASS only with zero CRITICAL and zero WARNING. *(ARS
  generator-evaluator.)*
- **G6 Compile.** `latexmk`/`pdflatex`+`biber` build is clean; undefined refs,
  missing figures, overfull boxes are reported. *(ppg12 compile-note.)*
- **G7 Disclosure / epistemic status.** Placeholders are honestly labelled and
  cross-checked against `docs/placeholder.md`; preliminary vs verified vs
  placeholder is explicit; AI-assistance disclosure available for the paper.

### Negative constraints (must NOT do)
- Must NOT enable the ARS plugin as an auto-routing skill (would collide with and
  displace `/review-note` `/review-paper`; wrong format/domain; deps/egress
  blocked). Port protocols only. See Design Decision D1.
- Must NOT implement figure sync as raw committed symlinks into `IntNotes/`
  (dangle on the CERN GitLab CI runner and standalone clones — the note compiles
  remotely with neither the Analysis repo nor the data area present). Use
  copy+provenance. See Design Decision D2.
- Must NOT require pip installs or live API egress beyond the allowlist.

## Context
- `IntNotes/` is a CERN GitLab submodule (ANA-HION-2023-07-INT1), biblatex+biber,
  master `ANA-HION-2023-07-INT1.tex`, `\graphicspath{{logos/}{figures/}}`,
  sections in `tex/`, bib in `bib/*.bib` + own `ANA-HION-2023-07-INT1.bib`.
  Compiles remotely via `atlas-physics-office/gitlab-integration` CI.
- Plots live in `Analysis/plots/` (and the data area), NOT git-tracked.
- Existing infra: `commands/{review-note,review-paper}.md`, agents
  `{executor,note-reviewer,numerical-verifier,physics-reviewer,plot-reviewer}`,
  CLAUDE.md Auto-Dispatch Rules, KB, `placeholder.md` ledger.
- References cloned to `/tmp/{acad-research-skills,ppg12}` (study copies).

## Scope
In: ARS summary doc; cherry-pick protocols; new commands `verify-citations`,
`sync-note-figures`, `check-note-sync`, `compile-note`; new agent
`citation-verifier`; enhance `numerical-verifier`, `note-reviewer`; upgrade
`review-note`/`review-paper`; CLAUDE.md + settings wiring; vendor ARS reference.
Out: the future generic agentic-workflow task; porting ARS Python scripts to run
live; thesis/website deploy (noted future).

## Design Decisions
- **D1 — Do not enable ARS plugin; vendor + cherry-pick.** Rationale in plan +
  summary doc: routing collision with the HEP note/paper commands; generic
  markdown/IMRAD/APA7 vs ATLAS LaTeX/biblatex/INSPIRE; no-pip + restricted egress
  block its automated gates; CSL-JSON-vs-bib citation model mismatch. Port the
  *protocols* into the existing workflow instead.
- **D2 — Copy-based figure sync + provenance manifest, not symlinks.** The note
  compiles remotely in isolation; symlinks to `../Analysis`/data dangle there and
  on clones (user already documented this for `placeholder.md`). Copy current
  plots into committed `IntNotes/figures/` + `figure_provenance.json`; optional
  local-only `--symlink` fast path (gitignored).
- **D3 — Fix stale path drift.** Existing commands/agents reference
  `IntNote/mydocument.tex`; the real path is `IntNotes/ANA-HION-2023-07-INT1.tex`.
  Correct during the upgrade.

## Implementation Plan (status)
1. [done] Tracking doc (this file) + register in CLAUDE.md Active Tracking Docs.
2. [done] ARS summary doc + ppg12 companion + vendored ARS protocol subset.
3. [done] `academic_writing_workflow.md` ground-truth spec + Quality Gate Checklist.
4. [done] New agent `citation-verifier.md`; enhanced `numerical-verifier.md`,
   `note-reviewer.md`.
5. [done] New commands `verify-citations`, `sync-note-figures`, `check-note-sync`,
   `compile-note` (+ tested helper `IntNotes/scripts/sync_note_figures.py`).
6. [done] Upgraded `review-note`/`review-paper` (embedded G1/G3/G4/G5/G7 +
   pre-flight + generator-evaluator pre-commitment; fixed D3 path drift).
7. [done] Wired CLAUDE.md Auto-Dispatch + Documentation-References pointer + Active
   Tracking Docs; added WebFetch inspirehep.net/doi.org to settings.local.json.
8. [done] Verified: figure-sync end-to-end (8 NEW→copied→CLEAN; re-run UP-TO-DATE;
   modified source→STALE→exit 1; restored intact); settings JSON valid; helper
   py_compile OK; no stale `IntNote/mydocument` refs remain.

## Progress Log
- 2026-06-16 — Step 1: created tracking doc with Workflow Specification (gates
  G1–G7), Design Decisions D1–D3. Plan approved (plan file
  `~/.claude/plans/abundant-humming-ocean.md`). References studied: ARS v3.12.1
  (4 skills, claim/citation/anti-leakage/figure-fidelity protocols), ppg12
  (sync-check, numerical-rederivation, note-critic, compile-note). Verified:
  IntNotes is biblatex+biber CERN submodule; plots untracked in `Analysis/plots/`;
  symlink approach unsafe for CI (D2).
- 2026-06-16 — Steps 2–8 (all done in one session):
  - Docs: `docs/references/academic_research_skills_summary.md` (full ARS summary +
    install decision + cherry-pick map), `ppg12_claude_summary.md`,
    `academic-research-skills-vendored/` (6 curated protocol files + LICENSE + NOTICE),
    `docs/academic_writing_workflow.md` (ground-truth G1–G7 + Quality Gate Checklist).
  - Agents: new `citation-verifier.md`; enhanced `numerical-verifier.md` (re-extract
    from ROOT, per-quantity tolerances, M1/M3 flags, CANNOT LOCATE) and
    `note-reviewer.md` (G1/G4/G7 + macro/writing-quality checklist; path fix).
  - Commands: new `verify-citations`, `sync-note-figures`, `check-note-sync`,
    `compile-note`; upgraded `review-note`/`review-paper` (Step-0 pre-flight,
    G1 anti-leakage in executor, reviewer pre-commitment + gate criteria, D3 fix).
  - Wiring: CLAUDE.md Auto-Dispatch (4 new rules + note/paper now reference the gate
    chain), Documentation-References pointer, Active Tracking Docs entry;
    settings.local.json += WebFetch inspirehep.net/doi.org.
  - Helper: `IntNotes/scripts/sync_note_figures.py` (stdlib-only, tested).
  - Verification: figure-sync 3-run test PASS (copy/clean/stale-detect, source
    restored intact); settings JSON valid; py_compile OK.

## Results & Observations
- **Environment:** cluster has `pdflatex`+`bibtex` but NOT `latexmk`/`biber` (even
  after `setup.sh`). Note is biblatex+biber → canonical full build is the CERN
  GitLab CI; local `/compile-note` does a structural pdflatex pass + reports bib
  UNRESOLVED. Encoded in `compile-note.md`.
- **Figure sync verified working** (copy-based, CI-safe): seeded manifest from the
  8 reco-eff placeholder plots; provenance stamps `Analysis` git SHA; staleness via
  sha256; STALE/MISSING-SOURCE/UNTRACKED/UNUSED/ORPHAN classes all exercised.
- **Note is a skeleton** (Intro/Datasets/Analysis/Results placeholders; no `\cite`
  yet; `tex/datasets.tex` uses `\href` JIRA links, not citations) → `/verify-citations`
  and `/check-note-sync` will be exercised on real content as sections get written.
- **Created artifacts in `IntNotes/figures/`** (working tree, NOT committed):
  `figure_manifest.yaml`, `figure_provenance.json`, 8 copied PNGs — ready for the
  reco-eff section to `\includegraphics`.

## Remaining Work (build complete; natural follow-ups)
- Fill `IntNotes/figures/figure_manifest.yaml` with the rest of the analysis figures
  as sections get written (currently seeded with reco-eff placeholders only).
- First real exercise of the full chain when an actual section is drafted; tune gate
  wording from that experience.
- Optional future: thesis/results website deploy (ppg12 `deploy-pages` analog).
- Decide when to `git add` the `IntNotes/figures/` artifacts into the CERN submodule
  (user's call — not committed by this task).

## Live test (2026-06-16) — first real sections written through the chain

**Wrote & gated the Datasets chapter + Introduction.** Per-section gate matrix:

| Section | G1 evid | G2 cite | G3 numbers | G4 fig | G5 review | G6 compile | G7 placeholder |
|---|---|---|---|---|---|---|---|
| `tex/datasets.tex` (PbPb/pp/MC) | PASS | PASS (4/4 verified, citation-verifier agent) | PASS (event-count sums re-derived MATCH; lumi MATCH) | N/A | FAIL→fixed (1 WARNING: NLO; +2 INFO addressed) | macros validated standalone; full build CI-only | PASS (`\PHtext` lumi-unc, `\PHbox` MC test-only/reco proxy) |
| `tex/introduction.tex` | PASS | PASS (HION-2019-58 verified; gaps `\PHtext`-flagged, not fabricated) | N/A | N/A | self+macro check | compiles standalone | PASS (2 citation gaps flagged red) |

**Placeholder highlighting (user req):** added `\PH`/`\PHtext`/`\PHbox` to
`ANA-HION-2023-07-INT1-defs.sty` (red + bold + large). 6 `\PH*` markers now in `tex/`
(greppable). Missing-but-needed citations flagged with `\PHtext` rather than fabricated.

**Workflow frictions found & fixed during the live test (steering):**
1. `\PH` used `\bfseries` → crashed in math mode (`$...\PH{}$`). **Fixed** math-safe
   (`\ifmmode`) in defs.sty. (Caught by a standalone compile.)
2. Local full build impossible: not only `biber`/`latexmk` missing but the ATLAS class
   font dep `newtxtext.sty` is absent → bare `pdflatex` of the master fails before
   content. **Codified** in `compile-note.md`: full build is CERN-CI-only; added a
   **standalone-wrapper section compile** as the local G6 substitute (used here to
   validate macros + both sections).
3. Reviewer (anti-sycophancy pre-commitment) caught a real WARNING (NLO undefined) +
   INFO (2024 GRL "unofficial"; pp event-count rationale) — all addressed. Confirms the
   G5 gate is not rubber-stamping.

**Artifacts:** `IntNotes/tex/{datasets,introduction}.tex`, master wired
(`\input{tex/introduction}`), `IntNotes/figures/` (8 synced placeholder plots +
manifest + provenance from the earlier build). Nothing committed (user's call).

**Deferred:** PbPb event-selection section — needs careful code-reading of the exact
5-cut values to meet G1/G3; a focused session, not rushed here.

## Latest Stage
**Build complete + first live test passed (2026-06-16).** Doc kept active pending the
user's review; natural next writing targets: event selection, centrality, trigger
efficiency. Close (remove from CLAUDE.md Active Tracking Docs) once those land or the
user signs off.
