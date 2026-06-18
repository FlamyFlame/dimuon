---
description: >
  Executor-reviewer loop for writing or editing internal note sections. Use
  when the user asks to write, edit, or update a section of the internal note,
  add a figure/table to the note, update captions, or says "write the
  introduction", "update the event selection section", "add this figure to
  the note", or any IntNote/ LaTeX editing work.
---

You are the **executor** in a structured executor-reviewer loop. A separate
reviewer subagent will independently evaluate your work each iteration.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `Analysis/docs/academic_writing_workflow.md` (**ground-truth gate spec G1–G7 — this loop enforces it**)
- `.claude/agents/executor.md` (your behavior as executor)
- `.claude/kb/index.md` (check for relevant KB articles)
- `IntNotes/ANA-HION-2023-07-INT1.tex` (master document structure; biblatex+biber, `\graphicspath{{logos/}{figures/}}`, sections in `IntNotes/tex/`)
- `Analysis/docs/placeholder.md` (≡ `IntNotes/placeholder.md`) — standing placeholders (G7)

Configuration:
- MAX_ITERATIONS = 5

## Step 0: Pre-flight (run before drafting — keeps the note on the latest version)

1. For any figure the section uses or will add: run `/sync-note-figures --check`
   (or full sync if you are adding/refreshing a figure) and resolve STALE/MISSING-SOURCE first.
2. Note that `/verify-citations` (G2) and `numerical-verifier` (G3) run inside the
   loop below; if the section is citation-heavy you may run `/verify-citations` now.

### Initialize log file

1. Derive a short slug from the task (3-5 words, kebab-case).
2. Set `LOG_FILE = .claude/logs/review-note-YYYYMMDD-HHMMSS-<slug>.md`
3. Write the log file path to `.claude/logs/active-session.txt` (overwrite).
4. Initialize the log file:

```
# Note Review Log
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
{"event":"start","command":"review-note","task":"$ARGUMENTS","slug":"<slug>","timestamp":"<ISO 8601>"}
```

## Compaction recovery

If you have lost context (e.g., after compaction), recover by:
1. Read `.claude/logs/active-session.txt` to get the log file path.
2. If that file is missing, run `ls -t .claude/logs/review-note-*.md | head -1`.
3. Read the log file fully — it is ground truth. Resume from where it left off.

## Loop

### Step 1: Execute

Write or edit the internal note section described in the task. Follow the
existing note structure in `IntNotes/`. Match the LaTeX style and conventions
already used in `IntNotes/tex/` and the ATLAS `atlasphysics`/`defs.sty` macros.

**Gate G1 — evidence-support + anti-leakage (apply while drafting):**
- Write from the analysis materials (ROOT outputs, `analysis_metadata.md`, KB,
  tracking docs), NOT from parametric memory. Every factual/quantitative claim must
  trace to a citation OR a named analysis output (file + histogram/bin/table/plot).
- For each number/figure, record the exact source anchor (you will hand numbers to
  `numerical-verifier` and figures to the figure-fidelity check).
- If the section needs content the materials do not cover, write `[MATERIAL GAP]` —
  do not invent procedures or numbers. Describe only what the pipeline actually does.
- Label any placeholder/preliminary result honestly, consistent with `placeholder.md` (G7).

Record:
- Files created or modified (with paths)
- Figures/tables added or referenced (+ their manifest source paths)
- Any numbers stated in prose (+ the ROOT file/histogram each came from)

### Step 2: Spawn reviewer subagent

Use the **Agent tool** to spawn a reviewer subagent with `subagent_type: "general-purpose"`.
The reviewer runs in its own context window — give it everything it needs.

Compose the reviewer prompt as follows (fill in the bracketed parts):

```
You are a document reviewer for an ATLAS internal note on dimuon continuum
measurement in pp and Pb+Pb collisions.

## Your task
Review the following note section and return a structured verdict.

## Pre-commitment (do this FIRST, before reading the draft — anti-sycophancy)
Based ONLY on the task description below, write down the 3–6 specific things a
correct version of this section MUST get right (expected claims, the numbers/figures
it should cite, the placeholders it must disclose, the likely failure modes). THEN
read the draft and check it against your pre-committed list. Do not lower the bar to
match what the draft happens to contain.

## Task description
[paste the original $ARGUMENTS]

## Files to review
[list every LaTeX file modified/created, with full paths — the reviewer will Read them]
[list the master document: IntNotes/ANA-HION-2023-07-INT1.tex]

## Numbers stated in prose
[list any numeric values written in the text, with their source]

## Review criteria

### Note checklist (items 1-8)
1. **Every figure is referenced in text**: each figure environment has at least one \ref in the body.
2. **Numbers in text match tables/figures**: numeric values in prose match corresponding data.
3. **Acronyms defined on first use**: first occurrence is written in full with abbreviation in parentheses.
4. **References complete**: all \cite commands resolve; no [?] in compiled output.
5. **Table captions self-contained**: each table caption is understandable without surrounding text.
6. **Figure captions describe content**: each caption states what is shown, which dataset, which selection.
7. **No placeholders remain**: no TODO, FIXME, XXX, or TBD markers in the source.
8. **LaTeX compiles cleanly**: no errors; overfull hbox > 10pt flagged.

### Evidence, citation & sync checklist (gates G1, G2, G4, G7 — apply to touched text)
9. **Evidence-support (G1)**: every factual/quantitative claim traces to a citation OR a named analysis output (ROOT file+histogram/bin, table, plot). A claim with neither → CRITICAL `[UNSUPPORTED]`.
10. **No memory-leakage (G1)**: methodology prose matches what the pipeline actually does (`Analysis/README.md` + pipeline docs); invented procedures / un-sourced facts → CRITICAL. `[MATERIAL GAP]` markers are acceptable, silent filler is not.
11. **Citations (G2)**: every `\cite` key resolves in a bib resource; each cited reference actually supports its claim. If not already run, instruct the executor to run `/verify-citations`; treat its MAJOR_DISTORTION / UNVERIFIABLE / UNDEFINED as CRITICAL.
12. **Figure-data fidelity + freshness (G4)**: each caption's interpretation follows from the data and states dataset+selection; each figure is current per `/sync-note-figures --check` (STALE/MISSING-SOURCE → CRITICAL).
13. **Placeholder honesty (G7)**: every placeholder/preliminary result is labelled and consistent with `placeholder.md`; an unlabelled placeholder presented as final → CRITICAL.
14. **Writing quality**: flag vague/hedging ("roughly", "seems"), AI-typical filler ("delve", "pivotal", "it is important to note that", em-dash overuse, throat-clearing openers), sentences > 40 words, inconsistent macro use.

### Physics checklist (apply only items relevant to this section)
1. Sign convention: OS=sign2/_op, SS=sign1/_ss.
2. Trigger mode: PbPb23/24/25=single_mu4/mode1, pp24=2mu4/mode3.
3. Event selection: PbPb 5-cut sequential order.
4. ZDC/FCal indexing: 0=C-side, 1=A-side.
5. Centrality: all PbPb years use PbPb2023 FCal-ET thresholds.
6. Differential scaling: Scale(N, "width"), not manual bin-width division.

### Plot checklist (for any embedded figures)
1. All required plots exist and are non-empty.
2. Legends present for multi-dataset canvases, labels comprehensible.
3. Legends/textboxes do not obscure data.
4. Axis labels readable (not cropped, not overlapping).

### Numerical verification (if numbers stated in prose)
For any numeric values stated in the text:
- Identify the source (ROOT file, table, figure, or computation).
- Verify the value matches the source. Report MATCH or MISMATCH.

### Anti-patterns to catch
- Figures included but never referenced in text
- Inconsistent notation (mixing p_T and pT, or GeV and GeV/c)
- Missing units on axis labels in embedded figures
- Stale figure files that don't match current analysis results

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

## Verdict rules (MANDATORY — follow exactly)

- **VERDICT: PASS** ONLY if zero CRITICAL AND zero WARNING issues.
- **VERDICT: FAIL** if ANY CRITICAL or WARNING issue exists.
- Do NOT return PASS with outstanding WARNINGs.

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

**MANDATORY — do ALL three in order:**

1. Update log:
```
**Status**: APPROVED at iteration [N]
**Summary**: [1-2 sentences]
```

2. **Emit end event** — append to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"review-note","status":"approved","iterations":<N>,"timestamp":"<ISO 8601>","criticals":0,"warnings":0,"slug":"<slug>"}
```

3. Delete `.claude/logs/active-session.txt`, then present final output to user.

## Exit: Escalate

**MANDATORY — do ALL three in order:**

1. Update log:
```
**Status**: ESCALATED at iteration [N] (max iterations reached)
**Unresolved issues**:
- [Issue 1]: [what it is, what was tried, why unresolved]
**Recommended action**: [specific suggestions]
```

2. **Emit end event** — append to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"review-note","status":"escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<count>,"warnings":<count>,"unresolved":["<short issue description>",...],"slug":"<slug>"}
```

3. Delete `.claude/logs/active-session.txt`, then present unresolved issues to user with:
"These issues need your judgment: [list with context]"
