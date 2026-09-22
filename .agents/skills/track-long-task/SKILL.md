---
name: track-long-task
description: Externalize plans, progress, decisions, evidence, and recovery state for long or multi-step work. Use for investigations with multiple hypotheses, analysis implementations spanning several edits, pipeline runs, publication workflows, agent-infrastructure changes, or any task likely to cross turns or context compaction.
---

# Track a long task

## Select the record

- Reuse the relevant active tracking document when one exists.
- For analysis work, create `Analysis/docs/tracking/<slug>.md` only when the task
  authorizes repository documentation changes.
- For Codex infrastructure, create `.codex/tracking/<slug>.md`.
- Never create a duplicate record merely because a different agent is working.

## Initialize before implementation

Write these sections before changing task artifacts:

1. `Objective`
2. `Physics Procedure` for analysis implementation, or `Authoritative
   specification` for other fragile workflows
3. `Context and scope`
4. `Safety invariants`
5. `Design decisions`
6. Numbered implementation plan with checkboxes and validation requirements
7. Append-only `Progress log`
8. `Remaining work`
9. `Latest stage`

For investigations, replace the implementation plan with hypotheses/substeps and
include append-only `Accumulated findings` and `Ruled out` sections.

## Per-step protocol

Before each step:

- Update `Latest stage` with the step, plan, files, and expected validation.
- Verify the plan against the authoritative procedure.
- Record design changes before implementing them; include physics rationale.

After each step:

- Mark checklist status.
- Append exact commands, paths, line locations, values, observations, and test
  results to the progress log.
- Update remaining work and ruled-out hypotheses.
- Do not leave a material finding only in conversation.

After compaction or resumption, read the record fully and resume from `Latest
stage`. On completion, re-read the objective and procedure, record a completion
summary, resolve every checklist item, and set `Latest stage` to `cleared`.
