<!-- Setup basics (remote cluster, no pip, ACLiC) inherited from parent .claude/CLAUDE.md — not repeated here. -->
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
topic that is NOT the subject of the current tracking doc, search for and
cross-reference other existing tracking docs in `Analysis/docs/tracking/` —
including CLOSED ones (those removed from "Active Tracking Docs"). They are not
deleted and often already settled the question. Cite the doc you relied on, and
do not assert from inference what a sibling doc has already established or refuted.

## Auto-Dispatch Rules

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

## Active Tracking Docs

- `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/analysis_status_summary.md` — Current analysis status: which steps updated with May 2026 skim
- `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/analysis_roadmap_2026_06.md` — Analysis roadmap (2026-06-10): IntNote readiness, missing inputs, full chain with dummies; task files in Analysis/docs/roadmap_tasks/
- `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/academic_writing_workflow.md` — Academic writing production chain (rigor + auto-sync): building the G1–G7 gates, new commands/agents, ARS-not-installed decision, copy-based figure sync
- `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/low_mass_dimuon_template_fit.md` — **Umbrella: full background-subtraction program that REPLACES provisional OS−SS in crossx/R_AA.** Low-mass (0–4 GeV) OS/SS minv template fit (single-b vs g→QQ̄/FE vs mixed-event combinatoric; SS-anchored G norm via k=G_SS/G_OS) + origin-blind eff/unfolding BEFORE fit + signal acceptance AFTER. Step 2 (data histos) in progress; acceptance + R_AA-integration steps (5–9) await user approval (scope-extended 2026-06-21).
- `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/raa_from_rdf_crossx.md` — task_06 R_AA from RDF crossx. REOPENED 2026-06-19: PbPb 2024 lumi corrected 1.59663→0.85112 nb⁻¹ (old GRL didn't exclude runs <489703; events already correct, lumi-only fix; combined R_AA ×1.17) + R_AA equation rewritten to common notation (n_AA raw yield, explicit 1/N_evt) + y-title relabel.
<!-- COMPLETED (2026-06-16), do NOT auto-load:
- Analysis/docs/tracking/reco_eff_placeholder_run2.md — Reco-eff placeholder (F.2 PbPb + HF R_AA Fig.31 pp); ε₁·ε₂ proxy in nominal crossx + R_AA. Follow-ups Q1+Q2 DONE 2026-06-16: Q1 PbPb genuine differential cross-section dσ/dp_T=1/L·dN (nb/GeV, differential_crossx dir; T_AA-weighted kept as R_AA input); Q2 reco folded into pp generic weight (generic_weight_col=w_reco_trig) so MC-data comparison reflects reco + INVARIANT: rerun MC-data comparison after any pp eff/det-resp/unfolding change. CLOSED. Follow-up: pt_150 differential crossx; proper 3D pair ε_reco when MC lands. -->

<!-- COMPLETED (2026-06-15), do NOT auto-load:
- Analysis/docs/tracking/kb_building.md — KB-building system (4 steps ALL DONE): /kb-build + /kb-review + GUIDE built; 12-source bulk build done & review-validated. Further sources via /kb-build ADD mode. -->


<!-- PARKED (2026-06-12), reopen when inputs ready — do NOT auto-load, but DO consult per the "Cross-reference other tracking docs" rule:
- Analysis/docs/tracking/hijing_overlay_reco_effcy_investigation.md — HIJING overlay reco efficiency: RESOLVED (deficit is PHYSICAL, method/sample-independent; dR fallback default-off); awaiting full ~10M-event signal-truth-only sample (≥2 months). -->

