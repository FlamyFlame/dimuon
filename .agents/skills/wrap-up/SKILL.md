---
name: wrap-up
description: Close a substantial dimuon task with tracking, documentation, verification, repository-state, and agent-workflow checks. Use when the user says wrap up, finalize, clean up, prepare to commit, or after completing a multi-step implementation or investigation.
---

# Wrap up a task

Read `.claude/commands/wrap-up.md` as the latest read-only housekeeping
reference, then apply this Codex-safe checklist:

1. Re-read the task objective, authoritative procedure, implementation plan,
   remaining work, and latest stage.
2. Confirm requested artifacts exist and run the relevant verification again.
3. Update only documentation made inaccurate by the authorized task. Preserve
   unrelated files and do not modify `.claude/`.
4. Finish the active tracking document: exact evidence, unresolved limitations,
   downstream reruns, completion summary, resolved checklist, and `Latest stage:
   cleared`.
5. Inspect `git status`, focused diffs, submodules, untracked artifacts, and
   accidental generated files. Do not absorb unrelated user changes.
6. Run `claude-to-codex-setting-sync` in dry-run mode if the task exposed agent
   workflow drift. Do not silently rewrite adapters.
7. Commit or push only if the user explicitly requested it. Stage explicit
   logical groups; never use `git add .` or `git add -A`.
8. Report artifacts, checks, documentation/tracking changes, remaining risks,
   and any restart or user action required.
