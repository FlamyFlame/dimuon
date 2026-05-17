---
description: >
  Executor-reviewer loop for investigating discrepancies or complex physics
  questions. Use when the user asks to investigate, debug, diagnose, or
  understand an unexpected result, discrepancy, anomaly, or says "why is this
  happening", "what's causing this", "look into this", "check if X explains Y",
  or any physics investigation/root-cause analysis.
---

You are the **executor** in a structured executor-reviewer loop. A separate
reviewer subagent will independently evaluate your work each iteration.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (your behavior as executor)
- `.claude/kb/index.md` (check for relevant KB articles)
- `Analysis/README.md` (analysis overview)

Configuration:
- MAX_ITERATIONS = 5

### Initialize log file

1. Derive a short slug from the task (3-5 words, kebab-case).
2. Set `LOG_FILE = .claude/logs/review-investigation-YYYYMMDD-HHMMSS-<slug>.md`
3. Write the log file path to `.claude/logs/active-session.txt` (overwrite).
4. Initialize the log file:

```
# Investigation Review Log
**Task**: $ARGUMENTS
**Log file**: <LOG_FILE basename>
**Started**: [timestamp]
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5
```

### Emit tracking event

Append to `.claude/logs/tracking.jsonl`:
```json
{"event":"start","command":"review-investigation","task":"$ARGUMENTS","slug":"<slug>","timestamp":"<ISO 8601>"}
```

Follow the tracking document protocol in CLAUDE.md if the investigation is
complex (multiple hypotheses or spanning many steps).

## Compaction recovery

If you have lost context (e.g., after compaction), recover by:
1. Read `.claude/logs/active-session.txt` to get the log file path.
2. If that file is missing, run `ls -t .claude/logs/review-investigation-*.md | head -1`.
3. Read the log file fully — it is ground truth. Resume from where it left off.

## Loop

### Step 1: Execute

Investigate the question posed in the task:
1. State the hypothesis or question clearly
2. Identify what evidence is needed
3. Examine the evidence (read files, run queries, produce diagnostic plots)
4. Test alternative explanations
5. Report conclusions with quantitative support

Record:
- Hypotheses tested and evidence for/against each
- Files examined (with paths)
- Numbers extracted (values, sources)
- Diagnostic plots produced (with paths)

### Step 2: Spawn reviewer subagent

Use the **Agent tool** to spawn a reviewer subagent with `subagent_type: "general-purpose"`.
The reviewer runs in its own context window — give it everything it needs.

Compose the reviewer prompt as follows (fill in the bracketed parts):

```
You are a physics reviewer for an ATLAS heavy-ion dimuon analysis investigation.

## Your task
Review the following investigation and return a structured verdict.

## Investigation question
[paste the original $ARGUMENTS]

## Executor's findings
[paste a complete summary of findings: hypotheses, evidence, conclusions]

## Numbers cited
[list every quantitative value cited, with its source file/histogram/tree]

## Files examined
[list all files read or produced, with full paths — the reviewer will Read them to verify]

## Review criteria

### Physics checklist (apply only items relevant to this investigation)
1. Sign convention: OS=sign2/_op, SS=sign1/_ss.
2. Trigger mode: PbPb23/24/25=single_mu4/mode1, pp24=2mu4/mode3.
3. Event selection: PbPb 5-cut sequential order.
4. ZDC/FCal indexing: 0=C-side, 1=A-side. FCal_Et_P=A, FCal_Et_N=C.
5. Centrality: all PbPb years use PbPb2023 FCal-ET thresholds; PbPb25 must recalculate.
6. Preamp cut: PbPb23/24=hard scalar, PbPb25=per-run mu+7sigma from TTree.
7. Pair eta binning: pair_eta_proj_ranges_coarse_incl_gap, NOT q_eta_proj_ranges_*.
8. Differential scaling: Scale(N, "width"), not manual bin-width division.

### Investigation quality checklist
1. **Conclusions supported by evidence**: every conclusion cites specific data.
2. **Alternative explanations considered**: at least one alternative tested or noted as untested.
3. **Numbers traceable**: every cited number traces to a specific source (file, histogram, tree).
4. **Correct conventions used**: sign, trigger, eta binning, centrality conventions followed.

### Numerical verification
For every number cited:
- Identify the source ROOT file and histogram/tree.
- Read the source file independently and extract the same quantity.
- Report MATCH (<0.1% relative difference) or MISMATCH with both values.
- For computed quantities, re-derive step by step.

## Response format

You MUST structure your response EXACTLY as follows:

## Issues
[numbered list of issues, or "None found." if all pass]

For each issue:
N. [CATEGORY] Description of the problem
   Evidence: [what you found that contradicts the executor]
   Severity: CRITICAL | WARNING | INFO
   Fix: [what the executor should do differently]

## Numerical Verification
[for each number: Quantity, Executor's value, Verified value, Source, MATCH/MISMATCH]
[or "No numbers to verify." if none reported]

## What Passed
[brief list of criteria that passed]

VERDICT: PASS

or

VERDICT: FAIL
```

### Step 3: Parse verdict and update log

Read the reviewer's response. Extract the verdict by searching for the line
starting with `VERDICT:`.

Append to the log file:

```
## Iteration N
**Reviewer verdict**: PASS | FAIL
**Issues found**: [count]
**Details**:
[paste the reviewer's Issues section verbatim]
**Numerical verification**:
[paste the reviewer's Numerical Verification section verbatim]
```

Update the `**Iterations completed**` count in the log header.

### Step 4: Decide

Count iterations by reading the log file (count `## Iteration` headers).

- **VERDICT: PASS** → go to **Exit: Approved**
- **VERDICT: FAIL** AND iteration count < MAX_ITERATIONS → go to **Step 5: Amend**
- **VERDICT: FAIL** AND iteration count >= MAX_ITERATIONS → go to **Exit: Escalate**

### Step 5: Amend

Address ONLY the specific issues listed in the reviewer's response.
Do not change conclusions or add investigations not flagged.
Log what you changed in the log file under the current iteration.
Go back to **Step 2** (spawn a fresh reviewer subagent for the amended work).

## Exit: Approved

Update log:
```
**Status**: APPROVED at iteration [N]
**Summary**: [1-2 sentences]
```

Append to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"review-investigation","status":"approved","iterations":<N>,"timestamp":"<ISO 8601>","criticals":0,"warnings":0,"slug":"<slug>"}
```

Delete `.claude/logs/active-session.txt`.
Present final conclusions to user with quantitative support.

## Exit: Escalate

Update log:
```
**Status**: ESCALATED at iteration [N] (max iterations reached)
**Unresolved issues**:
- [Issue 1]: [what it is, what was tried, why unresolved]
**Recommended action**: [specific suggestions]
```

Append to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"review-investigation","status":"escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<count>,"warnings":<count>,"unresolved":["<short issue description>",...],"slug":"<slug>"}
```

Delete `.claude/logs/active-session.txt`.
Present unresolved issues to user with:
"These issues need your judgment: [list with context]"
