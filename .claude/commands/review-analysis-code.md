---
description: Executor-reviewer loop for writing or modifying analysis code
---

You are running a structured executor-reviewer loop.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (general executor behavior)
- `.claude/agents/code-reviewer.md` (code review criteria)
- `.claude/agents/physics-reviewer.md` (physics correctness)
- `.claude/agents/numerical-verifier.md` (number verification)
- `.claude/conventions/codebase-patterns.md` (codebase conventions)
- `.claude/conventions/rdf-patterns.md` (RDataFrame conventions)
- `.claude/kb/index.md` (check for relevant KB articles)
- `Analysis/README.md` (class hierarchy and pipeline structure)

Configuration:
- MAX_ITERATIONS = 5
- LOG_FILE = `.claude/logs/review-analysis-code-$(date +%Y%m%d-%H%M%S).md`

Initialize the log file with a header:
```
# Analysis Code Review Log
**Task**: $ARGUMENTS
**Started**: [timestamp]
**Status**: IN PROGRESS
```

## Loop

### Iteration [N]

#### Execute
Write or modify the analysis code described in $ARGUMENTS. Follow existing code patterns in the repository. Compile with ACLiC and run to verify. Report all files modified, compilation status, and any computed values.

#### Review
Switch to reviewer mode. Apply ALL three review checklists:

**Code review** (from `code-reviewer.md`): check all 12 criteria (integer division, binning registration, branch names, DatasetTriggerMap, overflow, compilation isolation, Scale method, sign convention, pair eta, file I/O, magic numbers, cut ordering).

**Physics review** (from `physics-reviewer.md`): check applicable items — are physics choices correct for this code?

**Numerical verification** (from `numerical-verifier.md`): for any numbers reported by the executor (yields, efficiencies, counts), independently verify by reading the source ROOT file or re-deriving the computation.

For each criterion:
- State the criterion
- State PASS or FAIL
- If FAIL: cite specific evidence (file:line, variable name, value)
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
