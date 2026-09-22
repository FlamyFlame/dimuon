# Dimuon project guidance for Codex

## Scope and compatibility

This repository contains a Run 3 ATLAS heavy-ion dimuon physics analysis.
Preserve Claude Code's workflow: never modify `.claude/` or write Codex events
to `.claude/logs/`. Read Claude-owned context directly when useful, but apply
Codex-native tools, skills, hooks, and path semantics.

The Git root is the canonical Codex project root. Local plots outside the
repository are at `../plots/` from here and `../../plots/` from `Analysis/`.
BNL paths in `.claude/CLAUDE.md` are cluster equivalents, not valid local paths.

## Authority and context routing

For scientific claims, use this precedence:

1. The relevant implementation tracking document's `Physics Procedure` or
   academic workflow's `Workflow Specification`.
2. `Analysis/docs/analysis_overview.md` for stable physics methodology.
3. Current status and unresolved inputs in
   `Analysis/docs/tracking/analysis_roadmap_2026_06.md` and
   `Analysis/docs/tracking/analysis_status_summary.md`.
4. Repository code and pipeline documentation.
5. `.claude/kb/index.md` and the KB article it routes to.
6. Conversation history or Codex memory.

Never use memory as evidence for a number, dataset property, correction,
selection, uncertainty, or publication claim. Re-read its canonical source.
Treat anything labeled placeholder, guess, preliminary, open, or unverified as
such in code, plots, prose, and summaries.

Before domain work, read the relevant portions of `.claude/CLAUDE.md`, including
Tracking Documents, Documentation References, Auto-Dispatch Rules, and Active
Tracking Docs. Read every active tracking document fully before changing code
or planning an analysis step. The session-start hook supplies the current list,
but the files remain ground truth.

Use `rg` to route context instead of loading the whole KB. Start from:

- `.claude/kb/index.md` for domain knowledge.
- `Analysis/README.md` for architecture and pipelines.
- `Analysis/docs/analysis_overview.md` for scientific methodology.
- `Analysis/docs/tracking/` for current and historical decisions.
- `IntNotes/analysis_metadata.md` for confirmed note metadata and inputs.

## Long-task externalization

Use the `track-long-task` skill for investigations, multi-step implementations,
pipeline runs, publication work, and Codex infrastructure changes.

- Write the plan and `Latest stage` before implementation.
- Verify every analysis step against the authoritative Physics Procedure.
- Append exact results, paths, commands, values, and validation evidence after
  each step.
- Record rejected hypotheses and design changes with reasons.
- Re-read the objective and authoritative procedure before declaring complete.
- After compaction, trust the tracking document and checkpoint over a compressed
  conversational recollection.

For analysis work, use or create the canonical document under
`Analysis/docs/tracking/` when the task authorizes documentation changes. For
Codex configuration and other agent-infrastructure work, use
`.codex/tracking/`. Generated checkpoints belong in `.codex/state/` and are not
version controlled.

## Review and safety

Use the matching Codex skill when changing analysis code, plots, investigations,
pipelines, or publication text. Claude slash commands named in existing docs are
workflow references, not callable Codex commands. Do not claim an independent
review occurred unless it actually did. Spawn specialized subagents only when
the user explicitly requests delegation or parallel agents.

- No `pip install`; use the available release/runtime.
- Ask before grid submissions, Rucio rule changes, `pbook` actions, commits, or
  pushes unless the user explicitly requested that state change.
- Never initialize a VOMS proxy for the user.
- Preserve unrelated changes and submodule state.
- Prefer minimal, reversible changes and validate in proportion to scientific
  risk.
- For PDFs, use the available PDF tooling or Ghostscript extraction; do not
  infer content from filenames or abstracts.

## Completion standard

A task is complete only when requested outputs exist, relevant checks pass,
scientific claims are source-backed, tracking state is current, placeholders and
limitations are explicit, and unrelated repository state is untouched.
