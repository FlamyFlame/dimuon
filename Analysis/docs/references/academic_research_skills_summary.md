# `academic-research-skills` (ARS) — Summary, Relevance, and Cherry-Pick Decision

**Source:** https://github.com/imbad0202/academic-research-skills — "Academic
Research Skills", author Cheng-I Wu, **v3.12.1** (2026-06-15), license
**CC-BY-NC-4.0**. Studied from a local clone; a curated protocol subset is
vendored read-only under `academic-research-skills-vendored/` (with attribution).

This document (1) summarises ARS's full procedure, (2) assesses its relevance to
our ATLAS Run 3 dimuon analysis, (3) records the decision **not to enable it as a
plugin**, and (4) maps exactly which of its protocols we port. Companion:
`ppg12_claude_summary.md` (the HEP-specific reference). Our integrated workflow:
`../academic_writing_workflow.md`.

---

## 1. What ARS is

A **field-agnostic Claude Code plugin** for the full scholarly lifecycle —
*research → write → review → revise → finalize* — built as **4 skills + 27 modes +
~39 agents + 16 `/ars-*` slash commands + hooks**, backed by a large Python
`scripts/` library (≈100 checkers/tests) and a CI suite. Default reference domain
is higher-education / social-science / clinical; output is markdown/LaTeX/DOCX with
APA7 (+ Chicago/MLA/IEEE/Vancouver) citations.

### 1.1 The four skills

| Skill | Role | Modes (selected) |
|---|---|---|
| **deep-research** | upstream research engine (literature search, fact-check, systematic review) | full, quick, socratic, review, lit-review, three-way-scan, fact-check, systematic-review (PRISMA) |
| **academic-paper** | downstream writing engine (12 agents) | full, plan, outline-only, revision, revision-coach, abstract-only, lit-review, format-convert, citation-check, disclosure, rebuttal-audit |
| **academic-paper-reviewer** | multi-perspective peer review (5 reviewers: EIC + 3 peers + Devil's Advocate) | full, re-review, quick, methodology-focus, guided, calibration |
| **academic-pipeline** | orchestrator coordinating the above into a 10-stage flow | (dispatch only) |

### 1.2 The full pipeline (10 stages with integrity gates)
```
deep-research (socratic/full)
 → academic-paper (plan/full)
   → INTEGRITY CHECK (Stage 2.5, blocking)
     → academic-paper-reviewer (full/guided)
       → academic-paper (revision)
         → academic-paper-reviewer (re-review, ≤2 loops)
           → FINAL INTEGRITY CHECK (Stage 4.5, blocking)
             → academic-paper (format-convert → output)
               → Process Summary + AI Self-Reflection Report
```
Each stage requires user confirmation. Two **mandatory blocking integrity gates**
(2.5 pre-review, 4.5 final) must pass before proceeding.

### 1.3 The rigor machinery (the valuable core)

ARS's distinctive contribution is a stack of **anti-hallucination / integrity**
protocols, motivated by recent literature (Lu et al. *Nature* 2026 on AI-research
failure modes; Zhao et al. 2026 documenting 146,932 hallucinated citations across
preprints; Song et al. PaperOrchestra 2026):

1. **Anti-Leakage Protocol** — when the user supplies research materials, the
   writer builds the paper **primarily from those materials, not parametric
   memory**. Every claim must trace to a provided source; uncovered topics are
   flagged `[MATERIAL GAP]` instead of filled from memory. Directly attacks
   methodology fabrication and unintentional plagiarism.
2. **7-mode AI-research failure taxonomy** (blocking at 2.5/4.5): (1) implementation
   bug passing self-review, (2) hallucinated citation, (3) hallucinated experimental
   result, (4) shortcut reliance, (5) overclaiming/generalization, (6) methodology
   fabrication, (7) misleading visualization. The integrity reviewer must explicitly
   rule out or flag each.
3. **Citation existence verification** — deterministic, LLM-independent gate
   (v3.11): four-index resolution (Semantic Scholar + OpenAlex + Crossref + arXiv),
   cross-index *triangulation* to reduce false positives, a 90-day SQLite cache, and
   an opt-in `strict` terminal policy that hard-blocks ID-keyed unmatched citations.
4. **Claim verification (claim → source support)** — Phase E: existence is not
   enough; each numerical/categorical/trend/causal claim is traced to the specific
   passage of its cited source, with a verdict taxonomy (VERIFIED / MINOR_DISTORTION
   / MAJOR_DISTORTION / UNVERIFIABLE / UNVERIFIABLE_ACCESS). MAJOR or UNVERIFIABLE ⇒ FAIL.
5. **Three-layer citation emission** — every cited claim carries
   `claim → <!--ref:slug--> → <!--anchor:quote|page|section-->`, so the *locator*
   (exact quote/page) backing each citation is explicit and checkable.
6. **Generator-evaluator contract (anti-sycophancy)** — the reviewer commits to its
   expected criteria/scores **paper-blind** (Phase 1) before seeing the draft
   (Phase 2), so it cannot rubber-stamp; Devil's Advocate scores rebuttals 1–5 and
   does not concede below 4/5; cross-model divergence is flagged, not averaged.
7. **VLM figure verification** — a vision-LLM checks the *rendered* figure against
   the source data + standards (truncated labels, wrong scale, data mismatch,
   missing series), in a ≤2-iteration refinement loop. v3.12 adds a *figure/table
   fidelity* prose contract: does the caption's interpretation follow from the data,
   and does the manuscript cite the artifact for a claim it actually supports?
8. **Epistemic status / contamination signals + AI disclosure** — preprint-post-LLM
   and index-unmatched advisory flags; venue-specific AI-usage disclosure
   statements; a "Material Passport" ledger that records verification provenance and
   supports cross-session resume.

### 1.4 Packaging
Installs via `/plugin marketplace add imbad0202/academic-research-skills` +
`/plugin install`, or by symlinking the four skill dirs into `.claude/skills/`. Ships
a `SessionStart` announce hook and a `PreToolUse` write-scope guard. The Python
gates need `pyyaml, ruamel.yaml, jsonschema, pypdf, defusedxml` and live HTTP to the
four bibliographic indexes.

---

## 2. Relevance to our analysis

**Conceptually: high.** ARS targets exactly the failure classes that threaten an
LLM-assisted internal note / paper / thesis — fabricated or non-supporting
citations, numbers drifting from the real result, methodology prose invented from
memory, and self-review that rubber-stamps. Its protocol designs are well-motivated
and directly informative.

**Mechanically: low fit as-is**, for four reasons (these drive the decision below):
- **Format/domain mismatch.** ARS produces generic markdown/IMRAD/APA7 papers for
  social-science/clinical/higher-ed venues. We write **ATLAS LaTeX** notes/papers
  with **biblatex+biber**, ATLAS-curated `bib/*.bib` (INSPIRE-HEP exports), `atlasdoc`
  class, collaboration boilerplate, and `\todo`/placeholder discipline.
- **Citation-model mismatch.** ARS verifies CSL-JSON/markdown citations against DOI
  indexes; our citations are biblatex keys resolved against `.bib` resources and the
  HEP literature lives on **INSPIRE-HEP / arXiv**, not primarily Crossref.
- **Dependency + egress constraints.** Cluster policy is **no pip install**; egress
  is restricted (only arxiv.org / gitlab.cern.ch / a few allowlisted). ARS's
  load-bearing *deterministic* gates therefore cannot run here as shipped.
- **Routing collision.** ARS skills auto-trigger on "write a paper / review my paper
  / draft an abstract / write a section / write the methodology" — the **exact**
  intents that must fire our HEP-specific `/review-note` and `/review-paper`. Enabling
  ARS globally would create competing skills on the same triggers and could
  **displace** the correct ATLAS-LaTeX workflow — the opposite of the user's goal
  that the right procedure *actually* runs.

---

## 3. Decision: do NOT enable ARS as a plugin; vendor + port protocols

**We do not install/enable ARS as an auto-routing Claude Code plugin.** The four
reasons in §2 mean a wholesale install would, at best, sit dormant behind our
commands and, at worst, hijack note/paper routing with the wrong output format and
unrunnable gates. This is the user's explicitly-offered "give a reason + cherry-pick"
branch.

Instead we **port ARS's protocols** (prompt-level, no pip, no extra egress) into our
**single existing executor↔reviewer workflow**, and keep a curated, attributed,
read-only copy of the source protocols under `academic-research-skills-vendored/`.
CC-BY-NC-4.0 permits this non-commercial academic reuse with attribution (recorded
in the vendored `NOTICE.md`).

---

## 4. Cherry-pick map (ARS protocol → our implementation)

| ARS protocol | Ported into | As gate |
|---|---|---|
| Anti-Leakage (materials > memory; `[MATERIAL GAP]`) | `review-note`/`review-paper` executor preamble; workflow spec | **G1** evidence-support |
| Claim verification (claim→source support verdicts) | `agents/citation-verifier.md`; `commands/verify-citations.md` | **G2** citation support |
| Citation existence (4-index, adapted to INSPIRE/arXiv/DOI) | `commands/verify-citations.md` (biblatex `\cite`↔`.bib`) | **G2** citation existence |
| 7-mode failure taxonomy (esp. M1 impl-bug, M3 hallucinated result, M4 shortcut) | `agents/numerical-verifier.md` (M1/M3: number must trace to saved ROOT output; suspiciously-round / identical-error-bar checks); `review-*` integrity step (M4: closure/ablation) | **G3** + reviewer |
| Three-layer emission (claim→ref→anchor) | evidence-support gate requires each claim cite a *specific* figure/table/ROOT-number anchor | **G1/G2** |
| Generator-evaluator (paper-blind pre-commitment) | reviewer-subagent prompt in `review-note`/`review-paper` | **G5** anti-sycophancy |
| VLM / figure-table fidelity | `note-reviewer` figure checks + `commands/sync-note-figures.md` provenance | **G4** figure-data fidelity |
| Epistemic status / disclosure | G7 + existing `docs/placeholder.md` ledger | **G7** disclosure |
| Writing Quality Check (AI-typical phrasing) | `note-reviewer` style checklist deltas | **G5** (style) |

What we deliberately **do not** port: the markdown/APA7 formatter, the deep-research
literature-search engine (our `/kb-build` covers KB ingestion), the Material Passport
schema machinery, and all Python deterministic gates (deps/egress).

---

## 5. Bottom line
ARS is a strong *design reference* and a weak *drop-in tool* for this repo. We adopt
its rigor philosophy and concrete protocols, expressed as prompt-level gates inside
the ATLAS-LaTeX workflow we already have, so the quality procedures run on every real
writing task without importing an incompatible, unrunnable, routing-colliding plugin.
