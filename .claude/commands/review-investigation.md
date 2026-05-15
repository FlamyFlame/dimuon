---
description: Executor-reviewer loop for investigating discrepancies or complex physics questions
---

You are running a structured executor-reviewer loop.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (general executor behavior)
- `.claude/agents/physics-reviewer.md` (physics review criteria)
- `.claude/agents/numerical-verifier.md` (number verification)
- `.claude/kb/index.md` (check for relevant KB articles)
- `Analysis/README.md` (analysis overview)

Configuration:
- MAX_ITERATIONS = 5
- LOG_FILE = `.claude/logs/review-investigation-$(date +%Y%m%d-%H%M%S).md`

Initialize the log file with a header:
```
# Investigation Review Log
**Task**: $ARGUMENTS
**Started**: [timestamp]
**Status**: IN PROGRESS
```

## Loop

### Iteration [N]

#### Execute
Investigate the question posed in $ARGUMENTS. Structure the investigation:
1. State the hypothesis or question clearly
2. Identify what evidence is needed
3. Examine the evidence (read files, run queries, produce diagnostic plots)
4. Test alternative explanations
5. Report conclusions with quantitative support

Follow the tracking document protocol in CLAUDE.md if the investigation is complex (multiple hypotheses or spanning many steps).

#### Review
Switch to reviewer mode. Apply BOTH review checklists:

**Physics review** (from `physics-reviewer.md`): are conclusions physically sensible? Are the right variables and conventions used? Is the event selection correct for the context?

**Numerical verification** (from `numerical-verifier.md`): independently verify every number cited in the investigation. Re-derive computed quantities step by step.

Additionally check:
- Are conclusions supported by the evidence presented?
- Are alternative explanations considered and either tested or explicitly noted as untested?
- Are all cited numbers traceable to a specific source (file, histogram, tree)?

For each criterion:
- State the criterion
- State PASS or FAIL
- If FAIL: cite specific evidence
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
