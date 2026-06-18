---
description: >
  Compile the ATLAS internal note (latexmk/pdflatex+biber) via its Makefile and
  report build health (undefined refs/citations, missing figures, overfull boxes).
  Use when the user says "compile the note", "build the note", "does the note
  compile". Implements Gate G6 of Analysis/docs/academic_writing_workflow.md.
---

Compile the internal note and report build health.

Arguments: $ARGUMENTS
- (empty) → full build via the Makefile default (`latexmk -pdf`)
- `clean` → `make clean` first, then full build
- `quick` → single `pdflatex` pass (fast, no bibliography resolution)

Working directory: `IntNotes/` (submodule; master `ANA-HION-2023-07-INT1.tex`,
biblatex+biber).

## Environment reality (verified 2026-06-16, live test)
This cluster has `pdflatex` and `bibtex` but **NOT `latexmk` or `biber`** (even after
`source /usatlas/u/yuhanguo/setup.sh`), **and the local TeX is missing ATLAS class
dependencies** (e.g. `newtxtext.sty`), so even a bare `pdflatex` of the full master
`ANA-HION-2023-07-INT1.tex` **fails before reaching any content**. The note is
biblatex+biber, so the **canonical full build runs on the CERN GitLab CI**
(`.gitlab-ci.yml` → `atlas-physics-office/gitlab-integration`), not locally.

**Local G6 substitute (standalone section/macro check)** — when the full build can't run
locally, validate the LaTeX of a *section* (and the `\PH*` placeholder macros) in a
minimal wrapper that does not load the ATLAS class:
```latex
\documentclass{article}\usepackage{xcolor}
\providecommand{\TeV}{\ensuremath{\,\mathrm{TeV}}}\providecommand{\cite}[1]{[#1]}
\input{IntNotes/ANA-HION-2023-07-INT1-defs.sty}   % the real placeholder macros
\begin{document}\input{IntNotes/tex/<section>.tex}\end{document}
```
`pdflatex` this in `/tmp`. It catches undefined macros, math-mode errors, and broken
tables/`\ref` syntax (it will not resolve `\cite`/cross-document `\ref` or ATLAS-class
styling — report those as "CI-only"). This is the recommended local G6 check here.

## Steps
1. Probe the toolchain: `which latexmk pdflatex biber`. If only `pdflatex` is
   present, proceed in **structural mode** and say so in the report. If even
   `pdflatex` is missing, `source /usatlas/u/yuhanguo/setup.sh` and re-probe; if
   still missing, report that the note cannot be compiled here and stop (do not fake a build).
2. Build:
   - `clean`: `make clean` then proceed.
   - `quick` / structural (no biber): `pdflatex -interaction=nonstopmode ANA-HION-2023-07-INT1.tex`
     (run twice to settle `\ref`s). Report bibliography as **UNRESOLVED (needs biber/CI)**.
   - full (biber present): `make` (`latexmk -pdf`); fallback
     `latexmk -pdf ANA-HION-2023-07-INT1` or `pdflatex`+`biber`+`pdflatex`×2.
3. Report:
   - Whether `ANA-HION-2023-07-INT1.pdf` was produced (+ size, page count if available).
   - **Errors**: lines beginning `!` in the `.log`.
   - **Undefined references/citations**: grep the `.log` for
     `LaTeX Warning: Reference .* undefined`, `Citation .* undefined`,
     `There were undefined references`.
   - **Missing figures**: grep for `File .* not found` / `Unable to load picture` —
     cross-reference with the figure manifest (these usually mean `/sync-note-figures`
     was not run).
   - **Overfull boxes**: count `Overfull \hbox` > 10pt and list the worst few.
4. End with **BUILD: CLEAN** (PDF produced, no errors, no undefined refs/citations,
   no missing figures) or **BUILD: ISSUES** with the specific list.

This command may modify only LaTeX build artifacts (aux/log/pdf) — never `.tex`/`.bib`.
