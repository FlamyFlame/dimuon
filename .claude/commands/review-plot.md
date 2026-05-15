---
description: Executor-reviewer loop for creating or fixing plots
---

You are running a structured executor-reviewer loop.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (general executor behavior)
- `.claude/agents/plot-reviewer.md` (plot review criteria)
- `.claude/agents/physics-reviewer.md` (physics content checks)
- `.claude/conventions/atlas-plotting.md` (plotting conventions)
- `.claude/kb/index.md` (check for relevant KB articles)

Configuration:
- MAX_ITERATIONS = 5
- LOG_FILE = `.claude/logs/review-plot-$(date +%Y%m%d-%H%M%S).md`

Initialize the log file with a header:
```
# Plot Review Log
**Task**: $ARGUMENTS
**Started**: [timestamp]
**Status**: IN PROGRESS
```

## Loop

### Iteration [N]

#### Execute
Create or modify the plot described in $ARGUMENTS. Follow the plotting conventions in `.claude/conventions/atlas-plotting.md`. Save output to the appropriate figures/plots directory. Report all files created and any computed values.

#### Review
Switch to reviewer mode. Apply BOTH review checklists:

**Plot review** (from `plot-reviewer.md`): check all 6 criteria (completeness, non-empty, legends, no obscured data, ATLAS style if requested, axis readability).

**Physics review** (from `physics-reviewer.md`): check applicable items — is the distribution physically reasonable? Are axis ranges sensible for this observable? Are sign conventions correct?

For each criterion:
- State the criterion
- State PASS or FAIL
- If FAIL: cite specific evidence (file path, plot element, code line)
- Classify: CRITICAL / WARNING

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
