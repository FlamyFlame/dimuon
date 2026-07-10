<!-- Setup basics (remote cluster, no pip, ACLiC) inherited from parent .claude/CLAUDE.md — not repeated here. -->
- **Reading PDFs:** the Read tool and WebFetch CANNOT read PDF content here (no poppler; WebFetch returns binary/abstract-only). Use `gs -sDEVICE=txtwrite` to extract text. See `.claude/kb/gotchas/reading_pdfs.md`.
- Read `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/README.md` and `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/` for analysis context (class hierarchy, pipelines, sample types)
- For any analysis change: always update and maintain the relevant documentation in those files

## NTuple-Processing Provenance (BLOCKING — for any task that reads analysis data/MC)

Recurring, high-impact failure mode: an agent writes standalone code that reads the **raw
NTUPs** and silently diverges from the ntuple-processing procedure (a missing cut, a
truth-navigation step, a weight, a mode flag) → wrong physics that a plot/number review can
miss. To prevent it, for ANY task that reads analysis data/MC (plots, numbers, templates,
studies, investigations):

1. **Prefer ntuple-processing OUTPUT.** Use the output of the ntuple processing
   (`Analysis/NTupleProcessingCode/*` → muon_pairs / single-muon / pythia-truth trees and
   their derived branches: flavor/origin categories, `from_same_b`, traced HF ancestry, etc.)
   whenever it can answer the question. If the needed quantity is not produced yet, **add a
   configuration mode/flag to the ntuple processing and rerun it** (distinct output suffix,
   never clobber nominal) — this is PREFERRED over re-deriving from raw NTUPs.
2. **If standalone-from-raw is genuinely unavoidable**, you MUST first READ the relevant
   ntuple-processing code and reproduce its procedure EXACTLY — every selection cut, truth
   step (e.g. the generator-block barcode cutoff for HIJING overlay,
   `PythiaTruthExtras.c` `pythia_only_barcode_cache` = first `truth_barcode>200000`), weight,
   and convention, identical to nominal. The ONLY settings that may differ are those the
   request explicitly changes for a stated physics purpose (e.g. loosen Δp/p for a Δp/p
   distribution). Even then, if only ONE setting changes, prefer (1) — add a flag + rerun.
3. **STOP AND ASK (hard blocking point).** If the request does not specify whether a given
   cut/procedure should be mirrored, and you believe a standalone deviation is justified for a
   physics reason, DO NOT proceed on your own judgment — stop and ask the user first.

Binds the main agent AND every subagent — quote this rule in delegated data/MC-analysis task
prompts. Reviewers MUST enforce it: see `.claude/conventions/ntuple-provenance.md`.

## Tracking Documents

**INVARIANT:** For every in-scope tracking doc, every step MUST: (1) write plan to Latest Stage BEFORE work, (2) verify against Physics Procedure, (3) append results to Progress Log AFTER work, with physics motivation where applicable. No exceptions, including after compaction.

**Doc triage — first action in any new conversation and after every compaction, before any non-Read tool call, planning, or code change:**

1. Read `Analysis/docs/tracking/INDEX.md` — one line per doc: status + scope.
2. Read fully every ACTIVE doc whose scope overlaps the request; usually 1–2. Read a CLOSED doc only if its scope clearly bears on the request — closed docs settle past questions, they do not govern current work.
3. When ambiguous, or the request is broad ("what's the analysis status?"), read it. Reading is cheap; contradicting a Physics Procedure is not.
4. Read nothing if the request is not analysis work (CLAUDE.md/meta questions, plugin or skill work, environment setup).
5. If the task later moves **outside** the scope of the docs you read — a new topic, a different pipeline, a different measurement — re-read INDEX.md and triage again before working on it.
6. State the triage in one line in your first message: which docs you read, and which **ACTIVE** docs you skipped.

**Scope accuracy:** INDEX.md scopes are the *only* thing triage sees. If a doc's content grows past what its scope line describes, rewrite that line (and bump `Updated`) in the same step that grows the doc. A stale scope silently hides the doc from every future conversation.

For the `low_mass_dimuon_template_fit.md` umbrella (62 kB): read its Physics Procedure + index sections, then only the sub-doc (A–E) that owns your thread.

Create a tracking doc when: (a) investigating an unknown root cause with
multiple hypotheses, or (b) the user requests documentation for a new
analysis or major rework spanning many editing cycles.

### Autonomy Contract (write early; survives compaction; prevents false early-stop)

**Trigger:** the user asks to run a task autonomously to completion ("until it's
fully done", "until the bug is truly gone", "produce all the final plots/numbers").
Such a request is itself a reason to have a tracking doc (case (b) above) — if none
exists, create one.

**As one of the FIRST steps** — after Doc triage / reading the relevant tracking
doc(s), *before* detailed planning — write a pinned **Autonomy Contract** block at
the top of the doc (right after Objective), with exactly these three fields:

```markdown
## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = <concrete final outputs / acceptance checks, derived from the request>
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).
```

Rules governing the block:
- **Agent fills `Done` itself** from the request — the requested deliverables are
  usually stated clearly (e.g. "3 steps each producing a plot set" → Done = those 3
  plot sets regenerated at their paths). This is how a post-compaction agent knows it
  is *not yet* done.
- **No pre-enumerated ambiguity list — copy the `Stop-and-ask` line VERBATIM.** Only
  `Done` is task-specific (angle-bracket placeholder); the `Stop-and-ask` line is fixed
  literal text — do NOT specialize it into a concrete list of anticipated ambiguities.
  Two reasons: (1) at this early write-time your knowledge of the task is limited and
  most real ambiguities surface only mid-task, so any list you write now is guessing;
  (2) a fixed list would both miss the task-specific points (usually the very ones not
  foreseen in the original request) and license reckless push-through on anything
  off-list. The stop-condition stays a runtime judgment call, biased toward stopping.
- **If you cannot pin `Done`** from the request, or you already see an ambiguity of
  unclear blocking-status at this initial stage → **STOP and ask before doing any
  work.** Ambiguities that instead surface later (during planning/implementation,
  after exploration) → stop and `AskUserQuestion` exactly as agents already do well;
  record the resolution and, if it moves the target, update `Done`.
- **Investigation `Done` includes the fix and its blast radius** — but *conditionally*:
  apply the fix only if confident (tests confirm it removes the root cause rather than
  masking it, and changes nothing else), then regenerate **every** result the fix
  affects (rerun affected pipelines; see `signal_selection_change_impact.md`). If you
  cannot resolve the issue, or are unsure the proposed fix is correct, that is itself a
  blocking ambiguity → STOP and ask; proceed-to-done does NOT apply.
- **On compaction:** the block is re-read as part of Doc triage (see Lifecycle);
  resume from the first unmet `Done` item — do not restart with a fresh confirmation.
- **Completion clears it:** when Done is met, mark the block `DONE` (or delete it)
  alongside the normal Completion step.
- **Delegated autonomous work:** put the same Mandate/Done/Stop-and-ask in the
  subagent's scratch doc *and* its task prompt (per §Delegated subagent memory).

### Document structure

**Start:** Create `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/<name>.md`. Register it in the
**Active** section of `Analysis/docs/tracking/INDEX.md` (top of the section) with a one-line
scope. Structure depends on mode:
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
- Bump the doc's `Updated` date in `Analysis/docs/tracking/INDEX.md`. If the work
  pushed the doc beyond what its INDEX scope line describes, rewrite that scope
  line now — a stale scope hides the doc from every future triage.

**Design changes (implementation):** Record old approach, new approach, and
reason in Design Decisions before proceeding. When a physics motivation
applies, reference the Physics Procedure section and include a concise
physics reason (not just a code rationale). Verify the new approach is
consistent with the Physics Procedure. If it contradicts, update the
procedure first (with user approval).

### Delegated subagent memory

The tracking-doc INVARIANT above binds the **main (orchestrator) agent**.
Reviewer subagents (/review-*) are stateless: they get everything in the task
prompt and **return a structured verdict** — no doc writes, low compaction
risk, nothing to change here. The rule below covers the *other* case:
delegating a **semi-complex investigation or implementation** (not a review)
to a subagent that may run long enough to hit its own context limit and
compact mid-task, silently losing in-flight findings before they ever reach
the orchestrator.

When you (the orchestrator) delegate such work via the Agent tool:

- **Scratch doc per subagent.** Instruct the subagent to checkpoint to its
  OWN scratch tracking doc, never the canonical one. Use a distinct path per
  subagent to avoid write races:
  `Analysis/docs/tracking/_sub_<task>_<n>.md` (or under the session scratchpad
  dir). One file per subagent — concurrent subagents MUST NOT share a file.
- **Append-only, as it works.** The subagent writes its plan before each
  step and appends findings/results (exact paths, line numbers, names,
  numbers) after each step — same discipline as the main protocol, so the
  doc survives the subagent's own compaction. Tell the subagent to re-read
  its scratch doc first if it detects compaction.
- **Subagents NEVER run git** and never edit the canonical tracking doc or
  other shared files (mirrors the kb-build lesson: concurrent `git`/shared-file
  writes collide on `.git/index.lock` and clobber siblings). The orchestrator
  owns all shared state.
- **Merge then clean up.** When the subagent returns, the orchestrator merges
  its scratch doc into the canonical tracking doc (Progress Log / Accumulated
  Findings, per the Per-step protocol) and THEN deletes the scratch file. If
  the subagent died, recover its findings by reading its scratch doc before
  retrying.

The subagent's returned summary is a convenience, not the source of truth —
the scratch doc is. Treat anything only in the return text (not in the doc)
as at risk.

### Lifecycle

**Continuity (CRITICAL):** At the start of every new conversation, and before
resuming after any context compression, run **Doc triage** (above). Scope it
from the current request; after compaction, from the task described in the
summary. Then re-read the Per-step protocol and INVARIANT above. For
implementation docs, re-read the Physics Procedure section first. **If the doc
has an ACTIVE Autonomy Contract, re-read it and resume from the first unmet
`Done` item — an autonomous task is not finished until Done is met; a compaction
is not a reason to stop and re-confirm.** The doc is ground truth — if
conversation history or compaction summaries conflict, trust the doc.

**How to detect compaction:** If you cannot recall reading the tracking
doc's full text in this conversation (i.e., there is no Read tool call
for it in your visible history), treat it as a compaction event and
re-read before proceeding. When in doubt, re-read — reading the doc is
cheap, skipping it risks contradicting the Physics Procedure.

**Completion:** Write final summary, clear Latest Stage, move the doc from the
**Active** to the **Closed** section of `Analysis/docs/tracking/INDEX.md`, and
make its scope line describe what the doc *concluded* (not what it set out to
do) — that line is all a future conversation will see. Never delete the file.
A doc that is finished-for-now but awaiting external inputs is **PARKED**: keep
it in Closed, say so in the scope line, and note what unblocks it.

**Never:** Keep findings/progress only in conversation (write to doc before
next action); start a step without writing plan to doc first; declare
complete without re-reading Objective and Physics Procedure; implement
code that contradicts the Physics Procedure without user approval.

## Documentation References

Before working on any task, check these existing docs:
- **High-level analysis overview (objective, observables, physics methodology, sample roles): `Analysis/docs/analysis_overview.md`** — the stable conceptual ground truth for implementation and academic writing (no status; status lives in the roadmap).
- **Signal-selection change impact / rerun map: `Analysis/docs/signal_selection_change_impact.md`** — **MUST-READ before adding, removing, or changing the value of ANY single-b signal-selection cut (minv, pair pT, q·η, ΔR, …), including selection systematics.** It enumerates the full recompile→rerun-hist-filling→replot blast radius (which code, which outputs go stale, what stays unchanged). The signal region itself is defined in `analysis_overview.md` §2.
- **Academic writing production chain (rigor + auto-sync gates G1–G7): `Analysis/docs/academic_writing_workflow.md`** — ground-truth spec that `/review-note`, `/review-paper`, `/verify-citations`, `/sync-note-figures`, `/check-note-sync`, `/compile-note` enforce on EVERY writing task (even one section). Reference material: `Analysis/docs/references/academic_research_skills_summary.md` (why the ARS plugin is NOT installed) + `ppg12_claude_summary.md`.
- **Knowledge base — index of physics references: `.claude/kb/index.md`** — the curated literature/physics reference library, NOT just analysis bookkeeping. It holds: the two highest-priority Run 2 reference analyses ours derives from (HF-muon R_AA/v_n note+paper; back-to-back dimuon note+Letter), heavy-ion physics (especially heavy-flavor background), ATLAS muon detector (reco + trigger), centrality (ATLAS 2023 + Glauber), plus analysis bookkeeping (decisions, samples, variables, gotchas). **Consult the index for EVERY physics question/task/investigation/decision** — see the required-use rule below.
- Class hierarchy, code architecture, pipeline stages: `Analysis/README.md`
- Per-pipeline docs: `Analysis/docs/` (pythia_truth, pythia_fullsim_pp, pythia_fullsim_overlay, powheg, data_analysis)
- Skimming code, branches, grid workflow: `SkimCode/README.md`
- Data paths and directory layout: root `README.md`
- Internal note structure: `IntNotes/tex/` (section files), `IntNotes/ANA-HION-2023-07-INT1.tex` (master; biblatex+biber, CERN GitLab submodule)

**Use the KB for all physics work (REQUIRED).** The knowledge base index
(`.claude/kb/index.md`) is the mandatory entry point for every physics
question, task, investigation, and analysis decision. Workflow: consult the
index → see what references are available → pull the specific entries relevant
to the task (index-driven lookup — do NOT read every paper; that is the whole
point of the index). The two Run 2 reference analyses ours derives from are the
highest priority: HF-muon R_AA/v_n (note + paper) and back-to-back dimuon
(note + Letter). Beyond those, any question touching centrality (ATLAS 2023 +
Glauber), heavy-ion-specific physics (especially heavy-flavor background), or
ATLAS muon-detector specifics (reco/trigger) MUST be grounded in the relevant
KB entries before answering. Do not answer a physics question from memory when
a KB entry covers it — give physically grounded, reference-backed answers.

**Cross-reference other tracking docs.** Before making a factual statement on a
topic that is NOT the subject of the current tracking doc, consult
`Analysis/docs/tracking/INDEX.md` and cross-reference the sibling docs whose scope
covers it — **including CLOSED ones**. They are never deleted and often already
settled the question. Cite the doc you relied on, and do not assert from inference
what a sibling doc has already established or refuted.

## Auto-Dispatch Rules

Before writing/modifying ANY code that reads analysis data/MC, apply the **NTuple-Processing Provenance (BLOCKING)** rule above (prefer ntuple-processing output; if standalone-from-raw, mirror the procedure exactly; if unclear whether to mirror, STOP AND ASK). The `/review-*` reviewers verify provenance per `.claude/conventions/ntuple-provenance.md`.
When the user asks to create, fix, or modify a plot → invoke `/review-plot`
When the user asks to write, modify, or fix C++/ROOT/RDF analysis code → invoke `/review-analysis-code`
When the user asks to investigate a discrepancy, debug, or understand an unexpected result → invoke `/review-investigation`
When the user asks to write or edit an internal note section → invoke `/review-note` (enforces the Academic Writing gate chain G1–G7; see `Analysis/docs/academic_writing_workflow.md`)
When the user asks to write or polish paper text for publication → invoke `/review-paper` (publication-grade gate chain)
When the user asks to check/verify citations or references (real & supporting the claim) → invoke `/verify-citations`
When the user asks to sync, update, or check the note's figures against the latest analysis → invoke `/sync-note-figures`
When the user asks whether the note is up to date / matches the analysis → invoke `/check-note-sync`
When the user asks to compile or build the internal note → invoke `/compile-note`
When the user asks to review, audit, or validate a Claude Code plugin or skill → invoke `/review-plugin`
When the user asks to write, modify, review, or audit a multi-step analysis pipeline script → invoke `/review-pipeline`
When the user asks to run, execute, or steer a pipeline, or wants autonomous end-to-end pipeline execution → invoke `/steer-pipeline`
When the user asks to summarize a paper/source into the knowledge base, add references to the KB, or build/reorganize the KB → invoke `/kb-build` (criteria: `.claude/kb/KB_BUILDING_GUIDE.md`)
When the user asks to review, audit, or validate the knowledge base or a KB entry → invoke `/kb-review`

## Tracking Doc Index

All tracking docs — ACTIVE and CLOSED — are registered with a one-line scope in
**`Analysis/docs/tracking/INDEX.md`**. Doc triage (see §Tracking Documents) reads
that file first and pulls only the docs whose scope overlaps the request. Do not
maintain a duplicate list here; INDEX.md is the single source of truth for both
status and scope.

