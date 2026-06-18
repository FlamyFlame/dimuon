# Academic Writing Production Chain — Ground-Truth Specification

**Status:** authoritative spec for all internal-note, journal-paper, and thesis
writing in this analysis. The writing commands (`/review-note`, `/review-paper`)
and the support commands (`/verify-citations`, `/check-note-sync`,
`/sync-note-figures`, `/compile-note`) **must honor this document on every
invocation, including a single sentence or a single section.**

**Provenance:** ported from `academic-research-skills` (rigor protocols) and
`ppg12` (HEP mechanisms); see `docs/references/academic_research_skills_summary.md`
and `docs/references/ppg12_claude_summary.md`. We did **not** install the ARS
plugin (reasons there). Tracking: `docs/tracking/academic_writing_workflow.md`.

---

## 0. Principles

1. **Evidence or it doesn't ship.** Every factual or quantitative statement is
   backed by (a) a real, supporting citation, or (b) a verifiable analysis output
   (ROOT histogram/bin, table, or plot). No third category.
2. **Materials over memory (anti-leakage).** Write from the analysis outputs, KB,
   tracking docs, and `analysis_metadata.md` — not from the model's parametric
   memory. If the materials don't cover something the section needs, write
   `[MATERIAL GAP]`, never plausible-sounding filler.
3. **Honesty about status.** Preliminary / placeholder / verified is always
   explicit and consistent with `docs/placeholder.md`.
4. **Latest version always.** Numbers and figures in the note reflect the *current*
   analysis outputs; staleness is detected, not assumed away.
5. **No rubber-stamping.** The reviewer commits to its criteria before seeing the
   draft and fails on any CRITICAL/WARNING.

These restate, for HEP, the ARS anti-leakage + 7-failure-mode philosophy.

---

## 1. The chain

```
[research / KB / metadata]
   → DRAFT (materials-grounded)
   → G1  evidence-support audit
   → G2  citation existence + support      (/verify-citations)
   → G3  numerical re-derivation            (numerical-verifier agent)
   → G4  figure-data fidelity + sync        (/sync-note-figures, note-reviewer)
   → G5  reviewer loop (anti-sycophancy)    (/review-note · /review-paper)
   → G6  compile clean                      (/compile-note)
   → G7  disclosure / epistemic status
```

`/review-note` and `/review-paper` are the **entry points**: they run the draft +
embed G1, G3, G4, G5, G7, and call the support commands for G2, G4-sync, G6 as
pre-flight/post-flight steps. A standalone support command can also be run alone.

---

## 2. Quality Gate Checklist (apply every gate that is in scope for the task)

> Scope rule: even a one-sentence edit runs G1, G2, G7 on the touched text, and
> G3/G4 on any number/figure it introduces or relies on. G5 (reviewer loop) always
> runs via the review command. G6 runs whenever the build could be affected.

### G1 — Evidence-support (anti-hallucination)
- [ ] Every claim in the touched text traces to a citation **or** a named analysis
      output (file + histogram/bin/table/plot). Anchor the trace explicitly.
- [ ] Methodology prose describes only what the code/pipeline actually does
      (cross-check against `Analysis/README.md`, the pipeline docs, the tracking
      docs). No invented procedures. *(ARS failure mode 6.)*
- [ ] Uncovered required content is marked `[MATERIAL GAP]`, not filled from memory.
- [ ] No overclaiming: hedged where the evidence is preliminary/placeholder. *(M5.)*

### G2 — Citation existence + support  (`/verify-citations`)
- [ ] Every `\cite{key}` resolves in a `\addbibresource` `.bib`
      (`bib/*.bib` + `ANA-HION-2023-07-INT1.bib`); biblatex/biber.
- [ ] No undefined keys (`[?]` in build); flag bib entries never cited.
- [ ] For each cited reference, the adjacent claim is actually supported by it —
      verdict VERIFIED / MINOR_DISTORTION / MAJOR_DISTORTION / UNVERIFIABLE.
      MAJOR or UNVERIFIABLE ⇒ FAIL. Resolve/support-check via INSPIRE-HEP / arXiv /
      DOI where egress allows; otherwise mark `UNVERIFIED — needs manual` (do not
      silently pass). *(ARS citation existence + claim verification.)*
- [ ] No fabricated/placeholder references presented as real.

### G3 — Numerical re-derivation  (`numerical-verifier` agent)
- [ ] Every number stated in prose is **independently re-extracted from the source
      ROOT file** by the verifier — never trust the prose or the producing code.
- [ ] ROOT binning honored (bin 1 = first physical bin; overflow/underflow handled);
      `Scale(N,"width")` vs `Scale(N)`; OS/SS (`_op`/`_ss`) histogram picked correctly.
- [ ] Suspiciously round numbers / identical error bars across conditions are
      challenged (could be a constant leaking through a broken pipeline). *(M1/M3.)*
- [ ] Unlocatable source ⇒ `CANNOT LOCATE`, never a guessed value.

### G4 — Figure-data fidelity + sync  (`/sync-note-figures`, `note-reviewer`)
- [ ] Every `\includegraphics` figure is referenced by `\ref` in the body.
- [ ] Caption states what is shown, dataset, and selection, and its interpretation
      **follows from the data** (no claim the figure doesn't support). *(ARS fig-table fidelity.)*
- [ ] The committed figure in `IntNotes/figures/` is the **current** analysis output:
      `/sync-note-figures` provenance manifest shows source not newer than the copy,
      source present, and the figure is listed in the manifest.
- [ ] Axis labels readable; units present; ATLAS style; no misleading scale. *(M7.)*

### G5 — Reviewer loop (anti-sycophancy)  (`/review-note` · `/review-paper`)
- [ ] Reviewer subagent pre-commits its expected criteria/red-flags **before**
      reading the draft, then evaluates against them (generator-evaluator discipline).
- [ ] Verdict PASS only with **zero CRITICAL and zero WARNING**; loop ≤ 5 iterations,
      else ESCALATE to the user.
- [ ] Prose-quality deltas applied: vague/hedge language flagged, AI-typical phrasing
      (e.g. "delve/pivotal/robust"-as-filler, em-dash overuse, throat-clearing
      openers) flagged, `defs.sty`/`atlasphysics` macros used consistently,
      acronyms defined on first use. *(ARS writing quality check + ppg12 note-critic.)*
- [ ] Closure/ablation sanity for any results claim (does a control rule out the
      obvious artifact?). *(M4 shortcut reliance.)*

### G6 — Compile clean  (`/compile-note`)
- [ ] `latexmk -pdf` (or `pdflatex`+`biber`) builds with no errors.
- [ ] No undefined references / citations; no missing-figure errors; overfull
      `\hbox` > 10pt flagged.

### G7 — Disclosure / epistemic status
- [ ] Placeholders labelled and consistent with `docs/placeholder.md`
      (2024/2025 T_AA, reco-eff, σ_PbPb, etc. — must be stated in the note).
- [ ] Preliminary vs final results explicitly distinguished.
- [ ] (Paper/thesis) AI-assistance disclosure prepared per venue policy when required.

---

## 3. Tooling map (what enforces each gate)

| Gate | Enforcer |
|---|---|
| G1 | `/review-note` `/review-paper` executor preamble + reviewer checklist |
| G2 | `/verify-citations` → `agents/citation-verifier.md` |
| G3 | `agents/numerical-verifier.md` (invoked by the review loop) |
| G4 | `/sync-note-figures` (provenance/staleness) + `agents/note-reviewer.md` (caption fidelity) |
| G5 | `/review-note` `/review-paper` reviewer subagent |
| G6 | `/compile-note` |
| G7 | review commands + `docs/placeholder.md` |
| all | `/check-note-sync` runs G2(existence)+G3(spot)+G4(staleness)+G7 inventory as one read-only sweep |

## 4. Standard order of operations for a writing task
1. Identify the section + gather materials (KB, metadata, tracking docs, ROOT
   outputs, plots). 
2. `/sync-note-figures` for any figures the section uses (ensures latest).
3. Draft via `/review-note` (note) or `/review-paper` (publication-grade), which
   internally runs G1/G3/G4/G5/G7 and calls `/verify-citations` (G2).
4. `/compile-note` (G6).
5. `/check-note-sync` before any commit/freeze of the note (full read-only sweep).

## 5. Out of scope (here)
Generic agentic-workflow improvements (future task); running ARS Python gates live;
thesis/results website deploy. Figure sync is copy-based, not symlink (Design
Decision D2 in the tracking doc).
