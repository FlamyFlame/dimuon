---
description: >
  Review/validate the knowledge base under .claude/kb/. Use after /kb-build, or
  when the user asks to review, audit, or validate the KB or a KB entry. Two
  modes: SINGLE (one/few new sources) and FULL (whole-KB structure + use-case
  testing). Auto-dispatched per CLAUDE.md.
---

You are the **reviewer-orchestrator** for the KB. You run an executor-reviewer
loop: the builder (`/kb-build`) is the executor; you evaluate and send it back to
fix issues until it passes.

Task / target: $ARGUMENTS

## Setup (read first)
1. `.claude/kb/KB_BUILDING_GUIDE.md` — the criteria you enforce.
2. `.claude/kb/gotchas/reading_pdfs.md` — you also read PDFs with `gs` to verify.
3. `Analysis/docs/analysis_overview.md` + `docs/tracking/analysis_roadmap_2026_06.md`
   — relevance ground truth.

## Pick the mode
- **SINGLE** — one or a few sources were just added/changed.
- **FULL** — whole-KB review (after a bulk build, a reorg/step-3, or on request).

### Logging
- Slug; `LOG_FILE = .claude/logs/kb-review-YYYYMMDD-HHMMSS-<slug>.md`; write path to
  `.claude/logs/active-session.txt`.
- `tracking.jsonl` start event: `{"event":"start","command":"kb-review","mode":"SINGLE|FULL",...}`.
- MAX_ITERATIONS = 4.

---

## SINGLE mode — per-source review

Spawn a **fresh reviewer subagent** (`general-purpose`) per the prompt below. It
independently reads the new doc AND the source (via `gs` on the committed PDF, or
`curl`+`gs` for a supportive URL) and the source's reference list.

Reviewer prompt (fill brackets):
```
You review one knowledge-base source summary for an ATLAS Run-3 single-b dimuon
analysis. Criteria live in .claude/kb/KB_BUILDING_GUIDE.md — read it first.
Read PDFs with `gs -sDEVICE=txtwrite` (see gotchas/reading_pdfs.md); never WebFetch.

New/changed doc: [path]
Source PDF/URL: [path or url]
Analysis relevance ground truth: Analysis/docs/analysis_overview.md + roadmap.

### (I) Standalone criteria
1. Content correct & objective vs the source (spot-check ≥5 specific claims/numbers
   against the PDF via gs; report MATCH/MISMATCH with values + page).
2. Relevance to THIS analysis captured specifically (names concrete step / IntNote
   or thesis section / observable + use-type), not vague. (GUIDE §1,§4)
3. No hallucination / no physics assumptions; condition differences acknowledged
   but not quantified. (GUIDE §5)
4. Concision/detail balance proper — self-sufficient for common questions, no bulky
   irrelevant generic recap. (GUIDE §3,§7)
5. References fully scanned: every reference offering new crucial relevant info is
   listed (≤3) OR an explicit "none" is stated. Independently skim the source's
   reference list and flag any important ref the summarizer missed, and any listed
   future-read that is NOT actually worth including. (GUIDE §6)

### (II) Relationship to the full KB
1. Correct area directory? (GUIDE §8)
2. Relationships to other docs handled via `[[...]]`/links (bidirectional)?
3. Findable via index.md (entry exists, accurate)?
4. Deduplicated — shared concepts not re-explained but linked to a concept doc?

Output EXACTLY:
## Issues  (N. [CATEGORY] desc / Evidence / Severity CRITICAL|WARNING|INFO / Fix)
## Numerical Verification  (Quantity | summary value | verified value | page | MATCH/MISMATCH)
## What Passed
VERDICT: PASS   (only if zero CRITICAL and zero WARNING)  | VERDICT: FAIL
```

---

## FULL mode — structure + use-case validation

### Part A — structural audit (reviewer subagent, `general-purpose`)

> **Primary risk to hunt:** the KB was built one-subagent-per-paper, and those
> builders were blind to each other, so cross-paper duplication and missing
> cross-paper links are the *expected* failure mode. Actively look for two source
> docs on the same topic that (a) re-explain the same concept instead of sharing a
> concept doc, or (b) fail to link to each other. Do not assume integration was
> done just because each doc looks fine on its own.

Check, reporting Issues with severities + a verdict:
1. **Directory structure** sound and matches the step-3 design (one folder/area;
   per-paper + concept-doc balance). (GUIDE §8)
2. **Index** complete & accurate: every doc listed; no dead links; classifications
   shown.
3. **Knowledge graph**: related docs cross-linked bidirectionally; **no orphans**.
4. **Dedup**: no concept re-explained across docs without a canonical concept doc.
5. **FUTURE_READ.md** sane (≤3/source, deduped, not bloated).
6. **Self-sync (GUIDE §9):** every analysis step in the roadmap and every source
   area has a covering use-case scenario in Part B; flag gaps.

### Part B — use-case "unit tests" (the core of FULL review)

Treat the KB like software: run blind, scenario-based tests, then grade.

**B1. Maintain the scenario list (living — GUIDE §9).** Derive scenarios from the
current roadmap. Baseline set (update/extend as the analysis and KB grow):
1. **Implement an analysis procedure in code** — e.g. "apply the per-pair reco-eff
   weight ε_reco(pair pT, pair η, dR) in the crossx pipeline."
2. **Perform an investigation** — e.g. "why is the HIJING-overlay reco efficiency
   low at low pair pT?"
3. **Write a paper/IntNote methodology section** — e.g. "write the muon
   reconstruction & efficiency methodology paragraph for the IntNote." (Should find
   the detector docs + `concepts/muon_source_template_fits`.)
4. **Review code / plots / physics results** — e.g. "review a trigger-efficiency
   turn-on plot and its fit for physics correctness."
5. **Write the Introduction / physics motivation** — e.g. "draft the IntNote/thesis
   Introduction motivating b-quark in-medium energy loss in the QGP." (Should find
   the `physics/heavy_ion/` reviews — hub `open_hf_production` + `hi_big_picture`.)
6. **Reason about backgrounds / templates** — e.g. "explain the dominant correlated
   background to the single-b dimuon signal and how the template is built." (Should
   find `physics/background/gluon_splitting_flavour_excitation` + analysis_overview §6.)
7. **Centrality / R_AA normalization** — e.g. "justify the ⟨T_AA⟩ values used for
   the PbPb R_AA normalization and the Glauber method behind them." (Should find
   `physics/centrality/` — hub `glauber_modeling` + `atlas_centrality_2023`, which
   is the source of our hardcoded T_AA.)
(Add scenarios when new steps appear — unfolding, template fitting, systematics —
or when new sources are added; a scenario must exercise them. Scenarios 5–6 added
2026-06-15 (HI-field + background); scenario 7 added 2026-06-15 (centrality area).)

**B2. Run each scenario with a FRESH, BLIND test agent.** Spawn a
`general-purpose` agent given ONLY the realistic task — **no hint it is a test**,
no instruction to use the KB. After it works the task, ask it to report:
(1) all info/knowledge it needed, (2) every source/doc it consulted, (3) its
thought process and the path it followed — **including what it checked and found
NOT useful.** Save each transcript to the log.

**B3. Grade with a separate evaluator agent** (fresh `general-purpose`; the
orchestrator built/owns the KB, so grading is delegated to avoid bias). Give it
the transcripts + the KB index + the GUIDE. For each scenario judge:
- **Discovery:** did it consult the KB as an important source? *If it never checked
  the KB → CLAUDE.md / discoverability problem.*
- **Efficiency:** did it find the relevant info quickly? *If it checked the KB but
  missed relevant info that is present → indexing/graph problem.*
- **Comprehensiveness:** did it find ALL the relevant info (not just the first hit)?
Output per scenario: PASS/FAIL + the specific diagnosis → concrete fix (update
CLAUDE.md pointer / add index entry / add graph link / rewrite a Relevance
section / move a doc).

### FULL verdict
PASS only if structural audit passes AND every scenario passes. Otherwise FAIL
with the concrete fix list.

---

## Loop & exit (both modes)
1. Run the review (spawn the subagent(s) above). Append the verdict + Issues +
   (FULL) per-scenario diagnoses to the log.
2. **PASS** → Exit: Approved. **FAIL** & iter < MAX → send ONLY the listed fixes
   back to `/kb-build` (or fix index/graph directly if orchestrator-owned), then
   re-review with a fresh subagent. **FAIL** & iter ≥ MAX → Exit: Escalate.
3. **Exit:** update log status; `tracking.jsonl` end event
   `{"event":"end","command":"kb-review","mode":...,"status":"approved|escalated","iterations":<n>,...}`;
   delete `.claude/logs/active-session.txt`; report to the user.
</content>
