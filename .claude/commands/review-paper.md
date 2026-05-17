---
description: >
  Executor-reviewer loop for paper writing at publication-level quality. Use
  when the user asks to write, edit, or polish a paper section, prepare text
  for publication, or says "write the abstract", "polish this section for
  publication", "prepare the systematics discussion", or any publication-level
  LaTeX writing that requires formal language, complete uncertainty discussion,
  and comparison with prior measurements.
---

You are the **executor** in a structured executor-reviewer loop. A separate
reviewer subagent will independently evaluate your work each iteration.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (your behavior as executor)
- `.claude/kb/index.md` (check for relevant KB articles)
- `IntNote/mydocument.tex` (master document structure)

Configuration:
- MAX_ITERATIONS = 5

### Initialize log file

1. Derive a short slug from the task (3-5 words, kebab-case).
2. Set `LOG_FILE = .claude/logs/review-paper-YYYYMMDD-HHMMSS-<slug>.md`
3. Write the log file path to `.claude/logs/active-session.txt` (overwrite).
4. Initialize the log file:

```
# Paper Review Log
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
{"event":"start","command":"review-paper","task":"$ARGUMENTS","slug":"<slug>","timestamp":"<ISO 8601>"}
```

## Compaction recovery

If you have lost context (e.g., after compaction), recover by:
1. Read `.claude/logs/active-session.txt` to get the log file path.
2. If that file is missing, run `ls -t .claude/logs/review-paper-*.md | head -1`.
3. Read the log file fully — it is ground truth. Resume from where it left off.

## Loop

### Step 1: Execute

Write or edit the paper section described in the task. Follow the existing
document structure. Apply publication-level language quality: precise, formal,
no colloquialisms, consistent tense.

Record:
- Files created or modified (with paths)
- Figures/tables added or referenced
- Any numbers stated in prose

### Step 2: Spawn reviewer subagent

Use the **Agent tool** to spawn a reviewer subagent with `subagent_type: "general-purpose"`.
The reviewer runs in its own context window — give it everything it needs.

Compose the reviewer prompt as follows (fill in the bracketed parts):

```
You are a publication-level reviewer for an ATLAS paper on dimuon continuum
measurement in pp and Pb+Pb collisions.

## Your task
Review the following paper section and return a structured verdict.
Apply the FULL checklist (items 1-12), including publication-specific criteria.

## Task description
[paste the original $ARGUMENTS]

## Files to review
[list every LaTeX file modified/created, with full paths — the reviewer will Read them]
[list the master document: IntNote/mydocument.tex]

## Numbers stated in prose
[list any numeric values written in the text, with their source]

## Review criteria

### Document checklist (items 1-8, same as internal note)
1. **Every figure is referenced in text**: each figure environment has at least one \ref in the body.
2. **Numbers in text match tables/figures**: numeric values in prose match corresponding data.
3. **Acronyms defined on first use**: first occurrence is written in full with abbreviation in parentheses.
4. **References complete**: all \cite commands resolve; no [?] in compiled output.
5. **Table captions self-contained**: each table caption is understandable without surrounding text.
6. **Figure captions describe content**: each caption states what is shown, which dataset, which selection.
7. **No placeholders remain**: no TODO, FIXME, XXX, or TBD markers in the source.
8. **LaTeX compiles cleanly**: no errors; overfull hbox > 10pt flagged.

### Publication-specific checklist (items 9-12)
9. **Publication-level language**: no colloquialisms, consistent tense, passive voice where appropriate, precise word choices.
10. **Systematic uncertainties discussed**: each source of systematic uncertainty is identified, quantified, and its impact on the result stated.
11. **Comparison with prior measurements**: results compared to relevant published measurements with proper citations.
12. **ATLAS collaboration conventions**: author list, acknowledgements, and boilerplate sections present (if applicable to this section).

### Physics checklist (apply only items relevant to this section)
1. Sign convention: OS=sign2/_op, SS=sign1/_ss.
2. Trigger mode: PbPb23/24/25=single_mu4/mode1, pp24=2mu4/mode3.
3. Event selection: PbPb 5-cut sequential order.
4. ZDC/FCal indexing: 0=C-side, 1=A-side.
5. Centrality: all PbPb years use PbPb2023 FCal-ET thresholds.
6. Differential scaling: Scale(N, "width"), not manual bin-width division.
7. PbPb cross-section always combined (all years), never per year.

### Plot checklist (for any embedded figures)
1. All required plots exist and are non-empty.
2. Legends present for multi-dataset canvases, labels comprehensible.
3. Legends/textboxes do not obscure data.
4. Axis labels readable (not cropped, not overlapping).
5. ATLAS style applied (SetAtlasStyle, "ATLAS Preliminary" or similar).

### Numerical verification (if numbers stated in prose)
For any numeric values stated in the text:
- Identify the source (ROOT file, table, figure, or computation).
- Verify the value matches the source. Report MATCH or MISMATCH.

### Anti-patterns to catch
- Figures included but never referenced in text
- Inconsistent notation (mixing p_T and pT, or GeV and GeV/c)
- Missing units on axis labels in embedded figures
- Vague language ("large", "significant") without quantification
- Passive voice used inconsistently
- Systematic uncertainty sources mentioned but not quantified

## Response format

You MUST structure your response EXACTLY as follows:

## Issues
[numbered list of issues, or "None found." if all pass]

For each issue:
N. [CATEGORY] Description of the problem
   File: [path:line or figure/table number]
   Severity: CRITICAL | WARNING | INFO
   Fix: [specific instruction for what to change]

## Numerical Verification
[for each number: Quantity, Stated value, Verified value, Source, MATCH/MISMATCH]
[or "No numbers to verify." if none stated in prose]

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
{"event":"end","command":"review-paper","status":"approved","iterations":<N>,"timestamp":"<ISO 8601>","criticals":0,"warnings":0,"slug":"<slug>"}
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
{"event":"end","command":"review-paper","status":"escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<count>,"warnings":<count>,"unresolved":["<short issue description>",...],"slug":"<slug>"}
```

Delete `.claude/logs/active-session.txt`.
Present unresolved issues to user with:
"These issues need your judgment: [list with context]"
