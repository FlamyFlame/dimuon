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
- `.claude/conventions/physics-results-review.md` (MANDATORY physics-results criteria C1–C4 the reviewer enforces)
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
- If fitting was performed: report fit parameters, χ²/ndf, and fit range for each fit. These will be passed to the reviewer for numerical verification.

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

## Directory structure context
[If the user requested a specific directory structure or plot organization, paste that request here verbatim. If plots for a related physics procedure already exist elsewhere, list the existing paths and whether they should be reorganized. Otherwise write "No specific directory structure requested."]

## Review criteria

### Physics-results review (MANDATORY — apply FIRST, before the plot checklist)
**Read `.claude/conventions/physics-results-review.md` now and apply EVERY
applicable criterion (C1–C4).** These catch *physics* problems that a
readability/binning review misses (e.g. an efficiency-corrected spectrum that
jumps by ×10–100 at high p_T). The plot checklist below does NOT replace these.

For each plotted physics result (efficiency, cross-section, yield, R_AA, ratio,
correction factor, spectrum, mass distribution):
- **C1 Discontinuity / smoothness** — scan every panel/curve for any jump, kink,
  spike, or drop-to-zero (up or down). Physical distributions are smooth. Any
  **unexplained** discontinuity (one not from a named kinematic threshold or stated
  selection boundary, and not covered by the error bars) → **CRITICAL FAIL**.
- **C2 Shape & magnitude vs expectation (rubric-first)** — BEFORE judging the plot,
  write down what is plotted and an explicit list of physics expectations for its
  shape and magnitude (e.g. a cross-section vs p_T must be a smoothly, steeply
  *falling* spectrum; an efficiency ∈[0,1] plateauing at a sensible value; R_AA
  O(0.1–1.5)). Check each. Any shape violation or order-of-magnitude scale miss →
  **CRITICAL FAIL**.
- **C3 Run 2 reference cross-check** — for final/near-final results, compare
  magnitude & shape against the Run 2 references via the KB index
  (`.claude/kb/index.md`: `analysis/run2_hf_muon_raa.md`, `run2_dimuon_note.md`,
  `run2_dimuon_backtoback_paper.md`; trigger/reco magnitudes in
  `physics/detector/atlas_run2_muon_trigger.md`). Account for the analysis
  differences (nearby muon **pair**, not single HF-muon / not back-to-back) and
  Run 2→Run 3 differences (perf should improve *moderately*; an efficiency **much
  lower** than Run 2 for the same trigger signals a problem). Clear unexplained
  inconsistency → **CRITICAL FAIL**; if not comparable in-session, report
  `RUN2-CROSSCHECK UNVERIFIED` (never invent a Run 2 value).
- **C4** — any C1/C2/C3 CRITICAL is a *physics* failure: the calling command will
  auto-trigger an investigation (see the command's Decide step). In your report,
  label such issues clearly as `PHYSICS-RESULTS` so the loop routes them correctly.

### Plot checklist
For each item, state PASS or FAIL with specific evidence.

1. **All required plots are made**: every plot specified in the task exists.
2. **All plots are non-empty**: visible data in every plot.
3. **Legend for multi-dataset canvases**: legend present, labels comprehensible.
4. **Legends/textboxes do not obscure data**: no overlap with data points.
5. **ATLAS style (if required)**: SetAtlasStyle() called, "ATLAS Internal" present. Mark N/A if not requested.
6. **Axis labels readable**: not cropped, not overlapping, appropriate size.
7. **Plot output directory**: plots saved under the designated area for the sample type (data → `.../dimuon_data/plots/`, Pythia truth → `.../pythia_truth_full_sample/plots/`, Pythia fullsim PP → `.../pythia_fullsim_full_sample/plots/`, Pythia fullsim HIJING overlay → `.../pythia_fullsim_hijing_overlay_test_sample/plots/`, Powheg → `.../powheg_full_sample/plots/`). Check SaveAs paths in code.
8. **Directory structure clear**: subdirectories organized by physics topic; no flat dump of many PNGs. If "Directory structure context" above specifies a user-requested layout, verify it is implemented exactly.
9. **Related plots co-located**: if existing plots for the same physics procedure live elsewhere, verify both old and new plotting code are updated to use a common parent directory.

### Fitting checks (apply when plots include fits)
10. **Fit overlay visible**: fit curve drawn on data, visually distinguishable (different color/style).
11. **Fit follows data**: no wild divergence or large systematic residual pattern. Do NOT fail on moderate scatter — goodness-of-fit depends on statistics and model choice.
12. **Fit not obviously unconstrained**: flag fits driven by parameter bounds rather than data (featureless shape despite turn-on model, identical across many bins, data consistent with zero).
13. **Fit range and asymptote reasonable**: not extrapolated far beyond data; efficiency plateau ≤ 1; dR corrections approach expected asymptote at large separation.

### Efficiency-specific checks (apply when plots show efficiencies)
14. **Efficiency values in expected range**: standard efficiencies in [0, 1]. Inverse-weighted ratios (dR corrections from TH1::Divide with 1/ε weights) can fluctuate above 1 — flag only if far above 1 beyond what statistical error bars explain.
15. **Turn-on shape qualitatively correct**: pT turn-on should trend increasing (sigmoid-like). Do not require strict monotonic increase — raw unfitted data has statistical fluctuations, especially at high pT. Flag only if overall trend is non-physical.

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
- Plots saved to ad-hoc or flat directories instead of designated plot areas
- New plots for a procedure placed in a different directory from existing plots for the same procedure

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

First, **classify the FAIL** (if any): does any CRITICAL issue carry the
`PHYSICS-RESULTS` label or arise from criteria C1/C2/C3 (unexplained
discontinuity, shape/magnitude violation, or Run 2 inconsistency)?

- **VERDICT: PASS** → go to **Exit: Approved**
- **VERDICT: FAIL with a PHYSICS-RESULTS (C1/C2/C3) issue** → go to
  **Step 5b: Investigate** (do NOT cosmetically amend a physics problem).
- **VERDICT: FAIL, cosmetic/plotting issues only** AND iteration count < MAX_ITERATIONS
  → go to **Step 5: Amend**
- **VERDICT: FAIL** AND iteration count >= MAX_ITERATIONS → go to **Exit: Escalate**

### Step 5: Amend

Address ONLY the specific issues listed in the reviewer's response.
Do not refactor or change anything not flagged.
Log what you changed in the log file under the current iteration.
Go back to **Step 2** (spawn a fresh reviewer subagent for the amended work).

### Step 5b: Investigate (physics-results failure — criterion C4)

A C1/C2/C3 failure is a *physics* problem (possibly a methodology/procedure bug),
not a cosmetic plot fix. Do NOT just restyle or re-bin to hide it.
1. Log the physics finding in the plot-review log (which bins/panels, expected vs
   observed, the failing C-item).
2. **Invoke `/review-investigation`** on the flagged issue, passing the reviewer's
   evidence. It creates/uses a tracking doc and runs its own loop to find the root
   cause.
3. When the investigation returns a root cause, **fix it** (via `/review-analysis-code`
   if code changes are needed), regenerate the plot, then go back to **Step 2** for a
   fresh review.
   - **If the investigation concludes the feature is EXPECTED physics** (a named
     kinematic threshold or stated selection boundary, not a bug), there is nothing
     to fix: record that justification in the log AND include it in the next
     reviewer prompt (under "Numbers reported" / context) so C1's named-exception
     applies and the reviewer does not re-flag the same feature.
4. **If the investigation cannot resolve it** (escalates / hits its max iterations),
   go to **Exit: Escalate** and present the user the **full investigation report**
   (root-cause status, hypotheses tried, what remains) — never force a PASS.

## Exit: Approved

**MANDATORY — do ALL three in order:**

1. Update log:
```
**Status**: APPROVED at iteration [N]
**Summary**: [1-2 sentences]
```

2. **Emit end event** — append to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"review-plot","status":"approved","iterations":<N>,"timestamp":"<ISO 8601>","criticals":0,"warnings":0,"slug":"<slug>"}
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
{"event":"end","command":"review-plot","status":"escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<count>,"warnings":<count>,"unresolved":["<short issue description>",...],"slug":"<slug>"}
```

3. Delete `.claude/logs/active-session.txt`, then present unresolved issues to user with:
"These issues need your judgment: [list with context]"
