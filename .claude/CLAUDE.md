- Remote cluster (BNL SDCC); run `/usatlas/u/yuhanguo/setup.sh` for ROOT environment
- No pip install — use what's in the release
- For C++ class testing, use ROOT ACLiC: `.L MyClass.cxx+`
- **Reading PDFs:** the Read tool and WebFetch CANNOT read PDF content here (no poppler; WebFetch returns binary/abstract-only). Use `gs -sDEVICE=txtwrite` to extract text. See `.claude/kb/gotchas/reading_pdfs.md`.
- Read `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/README.md` and `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/` for analysis context (class hierarchy, pipelines, sample types)
- For any analysis change: always update and maintain the relevant documentation in those files

## Tracking Documents

**INVARIANT:** For every active tracking doc, every step MUST: (1) write plan to Latest Stage BEFORE work, (2) verify against Physics Procedure, (3) append results to Progress Log AFTER work, with physics motivation where applicable. No exceptions, including after compaction. **First action in any conversation or after compaction MUST be reading all active tracking docs — no code changes, no tool calls (other than Read), no planning until this is done.**

Create a tracking doc when: (a) investigating an unknown root cause with
multiple hypotheses, or (b) the user requests documentation for a new
analysis or major rework spanning many editing cycles.

### Document structure

**Start:** Create `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/<name>.md`. Register in Active Tracking
Docs below. Structure depends on mode:
- *Investigation:* Objective, Context, Sub-steps, Accumulated Findings
  (append-only), Ruled Out (append-only), Latest Stage.
- *Implementation:* Objective, **Physics Procedure** (REQUIRED — see below),
  Context, Scope, Design Decisions (with rationale), Implementation Plan
  (numbered, with status), Progress Log (append-only), Results &
  Observations (organized, mutable), Remaining Work, Latest Stage.

**Physics Procedure (implementation docs only):** Every implementation doc
MUST have this section immediately after Objective. It is the
**authoritative reference** — all design decisions, implementation choices,
variable names, and filter conditions must follow from it. If anything
contradicts it, flag to the user before proceeding. Contents:
1. **Motivation:** Why this measurement/correction is needed.
2. **Top-level equation:** The final result applied to the analysis.
   Define every symbol.
3. **Step-by-step method:** Each pipeline/step with (a) what physical
   quantity it measures, (b) the mathematical procedure, (c) explicit
   statements of what is and is NOT required (e.g., "no trigger
   requirement on the other muon").
4. **Negative constraints:** What the code must NOT do where confusion
   with similar-but-different procedures is likely (e.g., "this differs
   from mu4_mu4noL1 where the other muon must pass mu4").

**Implementation-specific sections:**
- **Results & Observations** (replaces Accumulated Findings): organize by
  topic, combine related items, delete obsolete entries when fixed. Not
  append-only.
- **No Ruled Out section** in implementation docs (keep for investigation
  only).

### Per-step protocol

**Before work:**
- Update Latest Stage with current step, plan, and files involved.
- Verify the planned step is consistent with the Physics Procedure section.
- Implementation Plan steps must reference the Physics Procedure section
  they implement (e.g., "per §3a") and specify which reviewer command
  applies (e.g., /review-analysis-code, /review-plot) and which Physics
  Procedure sections to include in the task prompt.

**During work:**
- Naming conventions in code must derive from physics terminology in the
  Physics Procedure, not from code convenience or other pipeline patterns.
- If reusing a code pattern from another pipeline, verify in the Physics
  Procedure that the physics justifies reuse. Different physics →
  different code, even if the structure looks similar.
- Steps that write or modify C++/ROOT/RDF code → /review-analysis-code.
  Include the relevant Physics Procedure section(s) in the task prompt.
  The reviewer must check physics correctness against the procedure.
- Steps that create or modify plots or fitting output → /review-plot.
- These reviewer rules apply to both user requests and agent-initiated
  work within an implementation plan, including after compaction.

**After work:**
- Append results to Progress Log (implementation) or Accumulated Findings
  (investigation) with step number. Keep exact values: paths, line numbers,
  names, numbers.
- For investigations: add ruled-out approaches to Ruled Out with reason.
- For implementations: mark step done and update Remaining Work.

**Design changes (implementation):** Record old approach, new approach, and
reason in Design Decisions before proceeding. When a physics motivation
applies, reference the Physics Procedure section and include a concise
physics reason (not just a code rationale). Verify the new approach is
consistent with the Physics Procedure. If it contradicts, update the
procedure first (with user approval).

### Lifecycle

**Continuity (CRITICAL):** At the start of every new conversation, and
before resuming after any context compression, check Active Tracking Docs.
If any exist, Read each fully before doing anything else — then re-read
the Per-step protocol and INVARIANT above. For implementation docs,
re-read the Physics Procedure section first. The doc is ground truth — if
conversation history or compaction summaries conflict, trust the doc.

**How to detect compaction:** If you cannot recall reading the tracking
doc's full text in this conversation (i.e., there is no Read tool call
for it in your visible history), treat it as a compaction event and
re-read before proceeding. When in doubt, re-read — reading the doc is
cheap, skipping it risks contradicting the Physics Procedure.

**Completion:** Write final summary, clear Latest Stage, remove from Active
Tracking Docs.

**Never:** Keep findings/progress only in conversation (write to doc before
next action); start a step without writing plan to doc first; declare
complete without re-reading Objective and Physics Procedure; implement
code that contradicts the Physics Procedure without user approval.

## Documentation References

Before working on any task, check these existing docs:
- **High-level analysis overview (objective, observables, physics methodology, sample roles): `Analysis/docs/analysis_overview.md`** — the stable conceptual ground truth for implementation and academic writing (no status; status lives in the roadmap).
- Knowledge base: `.claude/kb/index.md` (analysis overview, decisions, samples, variables, gotchas)
- Class hierarchy, code architecture, pipeline stages: `Analysis/README.md`
- Per-pipeline docs: `Analysis/docs/` (pythia_truth, pythia_fullsim_pp, pythia_fullsim_overlay, powheg, data_analysis)
- Skimming code, branches, grid workflow: `SkimCode/README.md`
- Data paths and directory layout: root `README.md`
- Internal note structure: `IntNote/tex/` (section files), `IntNote/mydocument.tex` (master)

**Cross-reference other tracking docs.** Before making a factual statement on a
topic that is NOT the subject of the current tracking doc, search for and
cross-reference other existing tracking docs in `Analysis/docs/tracking/` —
including CLOSED ones (those removed from "Active Tracking Docs"). They are not
deleted and often already settled the question. Cite the doc you relied on, and
do not assert from inference what a sibling doc has already established or refuted.

## Auto-Dispatch Rules

When the user asks to create, fix, or modify a plot → invoke `/review-plot`
When the user asks to write, modify, or fix C++/ROOT/RDF analysis code → invoke `/review-analysis-code`
When the user asks to investigate a discrepancy, debug, or understand an unexpected result → invoke `/review-investigation`
When the user asks to write or edit an internal note section → invoke `/review-note`
When the user asks to write or polish paper text for publication → invoke `/review-paper`
When the user asks to review, audit, or validate a Claude Code plugin or skill → invoke `/review-plugin`
When the user asks to write, modify, review, or audit a multi-step analysis pipeline script → invoke `/review-pipeline`
When the user asks to run, execute, or steer a pipeline, or wants autonomous end-to-end pipeline execution → invoke `/steer-pipeline`
When the user asks to summarize a paper/source into the knowledge base, add references to the KB, or build/reorganize the KB → invoke `/kb-build` (criteria: `.claude/kb/KB_BUILDING_GUIDE.md`)
When the user asks to review, audit, or validate the knowledge base or a KB entry → invoke `/kb-review`

## Active Tracking Docs

- `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/analysis_status_summary.md` — Current analysis status: which steps updated with May 2026 skim
- `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/analysis_roadmap_2026_06.md` — Analysis roadmap (2026-06-10): IntNote readiness, missing inputs, full chain with dummies; task files in Analysis/docs/roadmap_tasks/
- `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/kb_building.md` — KB-building system (4 steps): step 1 design DONE (criteria GUIDE + /kb-build + /kb-review); awaiting paper list for step 2 bulk build

<!-- PARKED (2026-06-12), reopen when inputs ready — do NOT auto-load, but DO consult per the "Cross-reference other tracking docs" rule:
- Analysis/docs/tracking/hijing_overlay_reco_effcy_investigation.md — HIJING overlay reco efficiency: RESOLVED (deficit is PHYSICAL, method/sample-independent; dR fallback default-off); awaiting full ~10M-event signal-truth-only sample (≥2 months). -->

