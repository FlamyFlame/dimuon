---
description: >
  Executor-reviewer loop for writing or modifying analysis code. Use when the
  user asks to write, modify, fix, or refactor C++/ROOT analysis code,
  RDataFrame pipelines, hist-filling classes, event selection code, or says
  "add a new variable", "fix the compilation", "update the cuts", or any
  analysis code changes in the Analysis/ directory.
---

You are the **executor** in a structured executor-reviewer loop. A separate
reviewer subagent will independently evaluate your work each iteration.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (your behavior as executor)
- `.claude/conventions/codebase-patterns.md` (codebase conventions)
- `.claude/conventions/rdf-patterns.md` (RDataFrame conventions)
- `.claude/conventions/physics-results-review.md` (physics-results criteria C1–C4 — apply when the code produces a physics result)
- `.claude/kb/index.md` (check for relevant KB articles)
- `Analysis/README.md` (class hierarchy and pipeline structure)

Configuration:
- MAX_ITERATIONS = 5

### Initialize log file

1. Derive a short slug from the task (3-5 words, kebab-case).
2. Set `LOG_FILE = .claude/logs/review-analysis-code-YYYYMMDD-HHMMSS-<slug>.md`
3. Write the log file path to `.claude/logs/active-session.txt` (overwrite).
4. Initialize the log file:

```
# Analysis Code Review Log
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
{"event":"start","command":"review-analysis-code","task":"$ARGUMENTS","slug":"<slug>","timestamp":"<ISO 8601>"}
```

## Compaction recovery

If you have lost context (e.g., after compaction), recover by:
1. Read `.claude/logs/active-session.txt` to get the log file path.
2. If that file is missing, run `ls -t .claude/logs/review-analysis-code-*.md | head -1`.
3. Read the log file fully — it is ground truth. Resume from where it left off.

## Loop

### Step 1: Execute

Write or modify the analysis code described in the task. Follow existing code
patterns in the repository. Compile with ACLiC and run to verify.

After completing the work, record what you did:
- Files created or modified (with paths)
- Compilation and runtime status
- Any numbers computed (yields, efficiencies, etc.) with their source

### Step 2: Spawn reviewer subagent

Use the **Agent tool** to spawn a reviewer subagent with `subagent_type: "general-purpose"`.
The reviewer runs in its own context window — give it everything it needs.

Compose the reviewer prompt as follows (fill in the bracketed parts):

```
You are a code reviewer for an ATLAS heavy-ion dimuon analysis codebase (C++/ROOT/RDataFrame).

## Your task
Review the following code changes and return a structured verdict.

## Task description
[paste the original $ARGUMENTS]

## Files to review
[list every file modified/created, with full paths — the reviewer will Read them]
[for each file, briefly state what was changed]

## Compilation/runtime result
[paste compilation output and any runtime output/errors]

## Numbers reported
[list any yields, efficiencies, or computed values — the reviewer will verify these]

## Review criteria

### Code correctness checklist
1. **No integer division** in ratios/efficiencies — uses double or explicit casts.
2. **Histogram binning registration**: new named binnings in var1D JSON have entries in hist_binning_map.
3. **Variable names match branch names** in SetBranchAddress/TTreeReader/RDF Define/Filter.
4. **DatasetTriggerMap consistency**: trigger/year mappings match DatasetTriggerMap.h.
5. **Overflow/underflow handling**: Integral(0, nbins+1) used where overflow must be included.

### Codebase pattern checklist
6. **ACLiC compilation isolation**: PP and PbPb hist-filling classes NOT compiled in same ROOT session.
7. **Scale method**: differential histograms use Scale(N, "width"), not manual bin-width division.
8. **Sign convention**: sign1=SS (_ss), sign2=OS (_op) — not swapped.
9. **Pair eta binning**: pair-level plots use pair_eta_proj_ranges_coarse_incl_gap, NOT q_eta_proj_ranges_*.

### General ROOT/C++ checklist
10. **File I/O**: ROOT files opened with error checking (!f || f->IsZombie()); files closed properly.
11. **No hardcoded magic numbers**: numeric constants with physical meaning are named or documented.
12. **Cut ordering**: PbPb event selection cuts in documented 5-cut sequential order.

### Physics checklist (apply only items relevant to this code)
1. Sign convention: OS=sign2/_op, SS=sign1/_ss.
2. Trigger mode: PbPb23/24/25=single_mu4/mode1, pp24=2mu4/mode3.
3. Event selection: PbPb 5-cut sequential order.
4. ZDC/FCal indexing: 0=C-side, 1=A-side. FCal_Et_P=A, FCal_Et_N=C.
5. Centrality: all PbPb years use PbPb2023 FCal-ET thresholds; PbPb25 must recalculate.
6. Preamp cut: PbPb23/24=hard scalar, PbPb25=per-run mu+7sigma from TTree.
7. Differential scaling: Scale(N, "width"), not manual bin-width division.

### Numerical verification
For any numbers reported above (yields, efficiencies, counts):
- Identify the source ROOT file and histogram/tree.
- Read the source file independently and extract the same quantity.
- Report MATCH (<0.1% relative difference) or MISMATCH with both values.
- For computed quantities, re-derive step by step.

### Physics-results review (MANDATORY when the code produces or changes a physics result)
If this code produces or alters a physics result (efficiency, cross-section,
yield, R_AA, ratio, correction factor, spectrum), **read
`.claude/conventions/physics-results-review.md` and apply C1–C4**, not just
"does the number match the code's own output" (the code may carry the bug):
- **C2 Shape & magnitude (rubric-first):** write the physics expectation for the
  produced quantity (shape AND scale) BEFORE checking, then verify the output
  against it — e.g. a 1/ε correction factor must be ≥1 and smoothly varying; an
  efficiency ∈[0,1] plateauing sensibly; a cross-section a smoothly falling
  spectrum. Order-of-magnitude scale miss or shape violation → **CRITICAL**.
- **C1 Discontinuity:** if the output is/feeds a distribution, scan for unexplained
  jumps across bins (or between corrected/uncorrected) → **CRITICAL**.
- **C3 Run 2 cross-check:** for final/near-final results, sanity-check magnitude
  vs the Run 2 references via `.claude/kb/index.md` (HF-muon R_AA, back-to-back
  dimuon; trigger/reco magnitudes), accounting for analysis and Run 2→Run 3
  differences (efficiency **much lower** than Run 2 for the same trigger → problem).
  Not comparable in-session → `RUN2-CROSSCHECK UNVERIFIED` (never invent a value).
- Label any such CRITICAL `PHYSICS-RESULTS` so the loop routes it to C4.

### Anti-patterns to catch
- Compiling PP and PbPb in the same ROOT session
- Using dynamic_cast instead of IsPbPb() virtual dispatch
- Forgetting to register new named binnings in hist_binning_map
- Using stale .so files after changing source
- TH3DModel with mixed uniform/variable axes — must use 8-arg form

## Response format

You MUST structure your response EXACTLY as follows:

## Issues
[numbered list of issues, or "None found." if all pass]

For each issue:
N. [CATEGORY] Description of the problem
   File: [path:line]
   Severity: CRITICAL | WARNING | INFO
   Fix: [specific instruction for what to change]

## Numerical Verification
[for each number: Quantity, Executor's value, Verified value, Source, MATCH/MISMATCH]
[or "No numbers to verify." if none reported]

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
- **VERDICT: FAIL, code/logic issues only** AND iteration count < MAX_ITERATIONS → go to **Step 5: Amend**
- **VERDICT: FAIL** AND iteration count >= MAX_ITERATIONS → go to **Exit: Escalate**

### Step 5: Amend

Address ONLY the specific issues listed in the reviewer's response.
Do not refactor or change anything not flagged.
Recompile and rerun to verify the fix.
Log what you changed in the log file under the current iteration.
Go back to **Step 2** (spawn a fresh reviewer subagent for the amended work).

### Step 5b: Investigate (physics-results failure — criterion C4)

A C1/C2/C3 failure means the produced result is unphysical (possibly a
methodology/procedure bug), not just a code-style fix. Do NOT patch the number to
look right.
1. Log the physics finding (which quantity/bins, expected vs observed, failing C-item).
2. **Invoke `/review-investigation`** on the issue (creates/uses a tracking doc, runs
   its own loop) to find the root cause; pass the reviewer's evidence.
3. Apply the root-cause fix, recompile/rerun, then go back to **Step 2** for a fresh review.
   - **If the investigation concludes the result is EXPECTED physics** (not a bug),
     record the justification and include it in the next reviewer prompt so it is not
     re-flagged.
4. **If the investigation cannot resolve it** (escalates / hits max iterations), go to
   **Exit: Escalate** with the **full investigation report** — never force a PASS.

## Exit: Approved

**MANDATORY — do ALL three in order:**

1. Update log:
```
**Status**: APPROVED at iteration [N]
**Summary**: [1-2 sentences]
```

2. **Emit end event** — append to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"review-analysis-code","status":"approved","iterations":<N>,"timestamp":"<ISO 8601>","criticals":0,"warnings":0,"slug":"<slug>"}
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
{"event":"end","command":"review-analysis-code","status":"escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<count>,"warnings":<count>,"unresolved":["<short issue description>",...],"slug":"<slug>"}
```

3. Delete `.claude/logs/active-session.txt`, then present unresolved issues to user with:
"These issues need your judgment: [list with context]"
