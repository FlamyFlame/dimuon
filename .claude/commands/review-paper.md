---
description: Executor-reviewer loop for paper writing (publication-level quality)
---

You are running a structured executor-reviewer loop.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (general executor behavior)
- `.claude/agents/note-reviewer.md` (document consistency — ALL checklist items 1-12, including paper-specific 9-12)
- `.claude/agents/physics-reviewer.md` (physics content checks)
- `.claude/agents/plot-reviewer.md` (for any embedded figures)
- `.claude/kb/index.md` (check for relevant KB articles)
- `IntNote/mydocument.tex` (master document structure)

Configuration:
- MAX_ITERATIONS = 5
- LOG_FILE = `.claude/logs/review-paper-$(date +%Y%m%d-%H%M%S).md`

Initialize the log file with a header:
```
# Paper Review Log
**Task**: $ARGUMENTS
**Started**: [timestamp]
**Status**: IN PROGRESS
```

## Loop

### Iteration [N]

#### Execute
Write or edit the paper section described in $ARGUMENTS. Follow the existing document structure. Apply publication-level language quality: precise, formal, no colloquialisms, consistent tense.

#### Review
Switch to reviewer mode. Apply ALL three review checklists:

**Note review — full paper checklist** (from `note-reviewer.md`, items 1-12): all internal note criteria PLUS publication-specific items (language quality, systematic uncertainty discussion, comparison with prior measurements, collaboration conventions).

**Physics review** (from `physics-reviewer.md`): are physics statements correct and complete?

**Plot review** (from `plot-reviewer.md`): for any figures referenced or included — apply the 6-item plot checklist.

For each criterion:
- State the criterion
- State PASS or FAIL
- If FAIL: cite specific evidence (file:line, figure/table number)
- Classify: CRITICAL / WARNING / INFO

#### Decision
- All PASS (no CRITICAL or WARNING) → go to [Exit: Approved]
- CRITICAL or WARNING items AND iteration < MAX_ITERATIONS →
  log findings, state specific amendments needed, go to [Amend]
- iteration >= MAX_ITERATIONS → go to [Exit: Escalate]

#### Amend
Address ONLY the specific issues identified by the reviewer.
Do not refactor or change anything not flagged.
Log what was changed and why.
Go back to [Review] (not [Execute] — only re-review, don't redo from scratch).

## Exit: Approved

Update log:
```
**Status**: APPROVED at iteration [N]
**Summary**: [1-2 sentences]
```

Present final output to user.

## Exit: Escalate

Update log:
```
**Status**: ESCALATED at iteration [N] (max iterations reached)
**Unresolved issues**:
- [Issue 1]: [what it is, what was tried, why unresolved]
- [Issue 2]: ...
**Recommended human action**: [specific suggestions]
```

Present unresolved issues to user with:
"These issues need your judgment: [list with context]"
