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
- `Analysis/docs/academic_writing_workflow.md` (**ground-truth gate spec G1–G7 — this loop enforces it at publication grade**)
- `.claude/agents/executor.md` (your behavior as executor)
- `.claude/conventions/physics-results-review.md` (physics-results criteria C1–C4 — apply to any quoted result/figure)
- `.claude/kb/index.md` (check for relevant KB articles)
- `IntNotes/ANA-HION-2023-07-INT1.tex` (master document structure; biblatex+biber, sections in `IntNotes/tex/`)
- `Analysis/docs/placeholder.md` — standing placeholders (G7); for a publication these must be resolved or explicitly scoped, not silently carried

Configuration:
- MAX_ITERATIONS = 5

## Step 0: Pre-flight
Run `/sync-note-figures --check` for any figure in scope (resolve STALE/MISSING first)
and `/verify-citations` for the section (publication grade ⇒ 100% citation
support-check, not a sample). Numerical re-derivation (G3) runs inside the loop.

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

**Gate G1 — evidence-support + anti-leakage (publication grade):** write strictly
from the analysis materials, not parametric memory; every claim traces to a citation
or a named analysis output; describe only what the pipeline actually does; `[MATERIAL
GAP]` over invention. Quantify every comparative/uncertainty statement (no vague
"large"/"significant"). Hedge nothing that the evidence supports, and overclaim nothing.

Record:
- Files created or modified (with paths)
- Figures/tables added or referenced (+ manifest source paths)
- Any numbers stated in prose (+ the ROOT file/histogram each came from)

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

## Pre-commitment (do this FIRST, before reading the draft — anti-sycophancy)
From the task description alone, write the 3–6 things a publication-grade version
MUST get right (expected claims, required citations + prior-measurement comparison,
the numbers/figures it must cite, systematic-uncertainty completeness, likely failure
modes). THEN check the draft against that list. Do not lower the bar to fit the draft.

## Task description
[paste the original $ARGUMENTS]

## Files to review
[list every LaTeX file modified/created, with full paths — the reviewer will Read them]
[list the master document: IntNotes/ANA-HION-2023-07-INT1.tex]

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

### Evidence, citation & sync checklist (gates G1, G2, G4, G7 — publication grade)
E1. **Evidence-support (G1)**: every factual/quantitative claim traces to a citation OR a named analysis output. Un-sourced claim → CRITICAL.
E2. **No memory-leakage (G1)**: methodology matches the actual pipeline; no invented procedures/numbers.
E3. **Citations (G2)**: 100% of cited references support their claims (VERIFIED), every `\cite` resolves; MAJOR_DISTORTION / UNVERIFIABLE / UNDEFINED → CRITICAL. (Run `/verify-citations`.)
E4. **Figure-data fidelity + freshness (G4)**: captions follow from the data; figures current per `/sync-note-figures --check`.
E5. **Epistemic status (G7)**: any preliminary/placeholder input is resolved or explicitly scoped for publication — an unlabelled placeholder in a paper draft → CRITICAL.

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

### Physics-results review (MANDATORY for any quoted physics result — publication grade)
For every physics result stated or shown (efficiency, cross-section, yield, R_AA,
ratio, purity, correction factor), **read `.claude/conventions/physics-results-review.md`
and apply C1–C4** — at publication grade an internally-consistent-but-unphysical
result is a hard stop:
- **C2 Shape & magnitude (rubric-first):** form the physics expectation (scale and
  shape) BEFORE checking, then verify (R_AA O(0.1–1.5), efficiency ∈[0,1] at a
  sensible plateau, falling cross-section spectrum). Scale/shape miss → **CRITICAL**.
- **C1 Discontinuity:** unexplained jumps in any quoted series / shown distribution → **CRITICAL**.
- **C3 Run 2 cross-check:** compare magnitude & shape against the Run 2 references
  via `.claude/kb/index.md` (HF-muon R_AA, back-to-back dimuon; trigger/reco
  magnitudes), accounting for analysis and Run 2→Run 3 differences (efficiency
  **much lower** than Run 2 for the same trigger → problem). For a publication a
  comparison plot/number vs prior measurements is expected where one exists; if not
  comparable in-session → `RUN2-CROSSCHECK UNVERIFIED` (never invent a value).
- Label such CRITICALs `PHYSICS-RESULTS` so the loop routes them to C4.

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
First **classify the FAIL**: does any CRITICAL issue carry the `PHYSICS-RESULTS`
label or arise from criteria C1/C2/C3 (discontinuity, shape/magnitude violation,
Run 2 inconsistency)?

- **VERDICT: PASS** → go to **Exit: Approved**
- **VERDICT: FAIL with a PHYSICS-RESULTS (C1/C2/C3) issue** → go to **Step 5b: Investigate**.
- **VERDICT: FAIL, prose/citation/figure issues only** AND iteration count < MAX_ITERATIONS → go to **Step 5: Amend**
- **VERDICT: FAIL** AND iteration count >= MAX_ITERATIONS → go to **Exit: Escalate**

### Step 5: Amend

Address ONLY the specific issues listed in the reviewer's response.
Do not refactor or change anything not flagged.
Log what you changed in the log file under the current iteration.
Go back to **Step 2** (spawn a fresh reviewer subagent for the amended work).

### Step 5b: Investigate (physics-results failure — criterion C4)

A C1/C2/C3 failure means a result destined for publication is physically wrong —
do NOT polish the prose around it.
1. Log the physics finding (quantity, expected vs observed, failing C-item).
2. **Invoke `/review-investigation`** on the underlying result (tracking doc + its
   own loop); pass the reviewer's evidence. Fix the root cause via
   `/review-analysis-code` / `/review-plot` and refresh the quoted numbers/figures.
3. Re-run this review (back to **Step 2**) once corrected.
   - **If the investigation concludes the result is EXPECTED physics** (not a bug),
     record the justification and include it in the next reviewer prompt so it is not
     re-flagged.
4. **If the investigation cannot resolve it**, go to **Exit: Escalate** with the
   **full investigation report** — never publish a result known to be off.

## Exit: Approved

**MANDATORY — do ALL three in order:**

1. Update log:
```
**Status**: APPROVED at iteration [N]
**Summary**: [1-2 sentences]
```

2. **Emit end event** — append to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"review-paper","status":"approved","iterations":<N>,"timestamp":"<ISO 8601>","criticals":0,"warnings":0,"slug":"<slug>"}
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
{"event":"end","command":"review-paper","status":"escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<count>,"warnings":<count>,"unresolved":["<short issue description>",...],"slug":"<slug>"}
```

3. Delete `.claude/logs/active-session.txt`, then present unresolved issues to user with:
"These issues need your judgment: [list with context]"
