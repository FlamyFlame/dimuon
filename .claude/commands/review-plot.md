---
description: >
  Executor-reviewer loop for creating or fixing plots. Use when the user asks
  to make a plot, fix a plot, update plot styling, add a new histogram panel,
  change axis labels/ranges, or says "plot X vs Y", "make a figure for",
  "the plot looks wrong", "update the figure", or any histogram/figure/canvas work.
---

You are the **executor** in a structured executor-reviewer loop. A separate
reviewer subagent will independently evaluate your work each iteration.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (your behavior as executor)
- `.claude/conventions/atlas-plotting.md` (plotting conventions)
- `.claude/kb/index.md` (check for relevant KB articles)

Configuration:
- MAX_ITERATIONS = 5

### Initialize log file

1. Derive a short slug from the task (3-5 words, kebab-case, e.g. `fcal-centrality-distribution`).
2. Set `LOG_FILE = .claude/logs/review-plot-YYYYMMDD-HHMMSS-<slug>.md`
3. Write the log file path to `.claude/logs/active-session.txt` (overwrite).
4. Initialize the log file:

```
# Plot Review Log
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
{"event":"start","command":"review-plot","task":"$ARGUMENTS","slug":"<slug>","timestamp":"<ISO 8601>"}
```

## Compaction recovery

If you have lost context (e.g., after compaction), recover by:
1. Read `.claude/logs/active-session.txt` to get the log file path.
2. If that file is missing, run `ls -t .claude/logs/review-plot-*.md | head -1`.
3. Read the log file fully — it is ground truth. Resume from where it left off.

## Loop

### Step 1: Execute

Create or modify the plot described in the task. Follow `.claude/conventions/atlas-plotting.md`.
Save output to the appropriate figures/plots directory.

After completing the work, record what you did:
- Files created or modified (with paths)
- Any computed values with their source
- Plots produced (with paths)

### Step 2: Spawn reviewer subagent

Use the **Agent tool** to spawn a reviewer subagent with `subagent_type: "general-purpose"`.
The reviewer runs in its own context window — it cannot see your conversation
history, so you must give it everything it needs in the prompt.

Compose the reviewer prompt as follows (fill in the bracketed parts):

```
You are a plot reviewer for an ATLAS heavy-ion dimuon analysis.

## Your task
Review the following plot work and return a structured verdict.

## Task description
[paste the original $ARGUMENTS]

## Files to review
[list every plot file produced, with full paths — the reviewer will Read them]
[list the code file(s) that produced them, with full paths — the reviewer will Read them]

## Numbers reported by executor
[list any computed values — integrals, yields, bin counts — or "None." if none reported]

## Review criteria

### Plot checklist
For each item, state PASS or FAIL with specific evidence.

1. **All required plots are made**: every plot specified in the task exists.
2. **All plots are non-empty**: visible data in every plot.
3. **Legend for multi-dataset canvases**: legend present, labels comprehensible.
4. **Legends/textboxes do not obscure data**: no overlap with data points.
5. **ATLAS style (if required)**: SetAtlasStyle() called, "ATLAS Internal" present. Mark N/A if not requested.
6. **Axis labels readable**: not cropped, not overlapping, appropriate size.

### Physics checklist (apply only items relevant to this plot)
1. Sign convention: OS=sign2/_op, SS=sign1/_ss.
2. Trigger mode: PbPb23/24/25=single_mu4/mode1, pp24=2mu4/mode3.
3. ZDC/FCal indexing: 0=C-side, 1=A-side.
4. Centrality: all years use PbPb2023 FCal-ET thresholds.
5. Pair eta binning: pair_eta_proj_ranges_coarse_incl_gap, NOT q_eta_proj_ranges_*.
6. Differential scaling: Scale(N, "width"), not manual bin-width division.
7. PbPb cross-section plots always combined (all years), never per year.

### Numerical verification (if numbers reported)
For any numbers reported by the executor:
- Identify the source ROOT file and histogram/tree.
- Read the source independently and verify. Report MATCH or MISMATCH.

### Anti-patterns to catch
- Saving PDF when only PNG was requested
- Default ROOT color palette instead of distinct readable colors
- Raw branch names as axis titles instead of physics labels
- Missing log scale when dynamic range > ~10x
- Ratio panels with Y-axis range hiding data points

## Response format

You MUST structure your response EXACTLY as follows:

## Issues
[numbered list of issues, or "None found." if all pass]

For each issue:
N. [CATEGORY] Description of the problem
   File: [path:line or plot file name]
   Severity: CRITICAL | WARNING | INFO
   Fix: [specific instruction for what to change]

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
Do not refactor or change anything not flagged.
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
{"event":"end","command":"review-plot","status":"approved","iterations":<N>,"timestamp":"<ISO 8601>","criticals":0,"warnings":0,"slug":"<slug>"}
```

Delete `.claude/logs/active-session.txt`.
Present final output to user.

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
{"event":"end","command":"review-plot","status":"escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<count>,"warnings":<count>,"unresolved":["<short issue description>",...],"slug":"<slug>"}
```

Delete `.claude/logs/active-session.txt`.
Present unresolved issues to user with:
"These issues need your judgment: [list with context]"
