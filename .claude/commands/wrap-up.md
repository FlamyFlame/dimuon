---
description: >
  End-of-task housekeeping: update docs, commit to git, and review .claude
  setup for improvements. Use after completing a substantial task — when
  the user says "wrap up", "finalize", "clean up and commit", or at the
  end of a multi-step implementation. NOT for mid-task use.
---

Task context: $ARGUMENTS
(If blank, infer from the conversation history and `git diff` what was done.)

Run the following checklist in order. Skip any section marked with a gate
condition that isn't met. Be concise — this is housekeeping, not a new task.

---

## 1. Documentation updates

### 1a. Tracking docs

Read the Active Tracking Docs list in CLAUDE.md. For each active doc:
- If this session's work touched that doc's topic: update Progress Log,
  mark completed steps, update Remaining Work and Latest Stage.
- If the task is fully complete: write final summary, clear Latest Stage,
  remove from Active Tracking Docs in CLAUDE.md.

### 1b. Repo docs

Check whether any of these need updates based on this session's changes:
- `Analysis/README.md` (class hierarchy, pipeline list, variable list)
- `Analysis/docs/*.md` (pipeline-specific docs)
- `SkimCode/README.md` (if skimming code was touched)

Only update docs where the content is now factually wrong or missing
something a future reader would need. Don't rewrite for style.

---

## 2. Git commit

Run `git status` and `git diff --stat` to see all changes.

- Group changes into logical commits (one per logical change — see
  feedback_git_commits memory).
- Stage and commit with clear messages. Include all modified files that
  belong to the completed work.
- Push to remote.
- If there are untracked files that look like temporary artifacts (test
  scripts, scratch files), list them for the user but don't commit.

---

## 3. CLAUDE.md review

**Gate: skip if this session only did a small, routine task.**

Read the full CLAUDE.md. Check for:

1. **Behavior consistency**: Did any behavior drift or failure occur this
   session that a new CLAUDE.md rule would prevent? If so, draft the
   minimal rule addition.
2. **Redundancy**: Are there rules that say the same thing in different
   words? If so, merge them.
3. **Staleness**: Are there references to things that no longer exist
   (removed files, old conventions, completed tracking docs still listed)?
4. **Conciseness**: Can any verbose section be shortened without losing
   meaning? Aim for the shortest phrasing that a model will still follow.

Present proposed changes to the user before applying. Don't rewrite
sections that are working fine.

---

## 4. Skill & command review

**Gate: skip if no executor-reviewer loop or skill was used this session.**

For each skill/command used this session:

### 4a. Triggering description
- Was the skill auto-triggered correctly, or did the user have to invoke
  it manually? If manual invocation was needed, update the `description:`
  field to cover the missed trigger pattern.
- Is the description concise? Remove trigger phrases that overlap with
  other skills.

### 4b. Reviewer criteria (for executor-reviewer commands only)
- Did the reviewer miss an issue that was later caught manually? If so,
  add a specific criterion for it.
- Did the reviewer flag false positives repeatedly? If so, tighten the
  criterion wording or add an exception.
- Are there criteria that are never actionable (always N/A)? Consider
  removing or consolidating them.

### 4c. Instruction clarity
- Were any steps skipped or misinterpreted during execution? If so,
  rewrite the instruction to be unambiguous.
- Are steps ordered correctly for the typical workflow?

Present proposed changes to the user before applying.

---

## 5. Summary

Output a short summary:
- Docs updated: [list or "none"]
- Commits made: [list with short messages]
- CLAUDE.md changes: [list or "none needed"]
- Skill/command changes: [list or "none needed"]
