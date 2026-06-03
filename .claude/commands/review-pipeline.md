---
description: >
  Executor-reviewer loop for writing or reviewing analysis pipeline scripts.
  Use when the user asks to write, modify, fix, or review a multi-step bash
  pipeline that orchestrates NTuple processing, RDF hist filling, fitting,
  plotting, or any end-to-end analysis chain. Also use when the user says
  "review this pipeline", "check the pipeline", "write a pipeline for", or
  points to an existing pipeline script for audit.
---

You are the **executor** in a structured executor-reviewer loop. A separate
reviewer subagent will independently evaluate your work each iteration.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (your behavior as executor)
- `.claude/conventions/codebase-patterns.md` (codebase conventions)
- `.claude/kb/index.md` (check for relevant KB articles)
- `Analysis/README.md` (class hierarchy and pipeline structure)
- `Analysis/docs/` — read the pipeline doc(s) relevant to this task

Configuration:
- MAX_ITERATIONS = 5

**MANDATORY SETUP**: You MUST complete ALL steps below (log file + tracking
event) BEFORE entering the Loop. Do not skip any step. If you are a subagent
spawned by a parent agent, these obligations still apply to you.

### Initialize log file

1. Derive a short slug from the task (3-5 words, kebab-case).
2. Set `LOG_FILE = .claude/logs/review-pipeline-YYYYMMDD-HHMMSS-<slug>.md`
3. Write the log file path to `.claude/logs/active-session.txt` (overwrite).
4. Initialize the log file:

```
# Pipeline Review Log
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
{"event":"start","command":"review-pipeline","task":"$ARGUMENTS","slug":"<slug>","timestamp":"<ISO 8601>"}
```

## Compaction recovery

If you have lost context (e.g., after compaction), recover by:
1. Read `.claude/logs/active-session.txt` to get the log file path.
2. If that file is missing, run `ls -t .claude/logs/review-pipeline-*.md | head -1`.
3. Read the log file fully — it is ground truth. Resume from where it left off.

## Loop

### Step 1: Execute

If the task targets an **existing pipeline**: read the script fully. If writing
a **new pipeline**: create it following patterns in `Analysis/pipelines/`.

After writing or reading, produce this analysis:
- List all stages the pipeline implements (numbered), mapping each to its
  analysis function (NTuple processing, hadd, RDF hist filling, fitting, etc.)
- List env vars, skip flags, and configuration knobs
- Identify expected input files and final output files
- Read the relevant pipeline documentation in `Analysis/docs/` and check
  whether the pipeline's stages match the documented procedure
- List files created or modified (with paths)

### Step 2: Spawn reviewer subagent

**CRITICAL**: The reviewer MUST be an independent agent — use the **Agent tool**
with `subagent_type: "general-purpose"` to spawn it. Do NOT review the code
yourself; do NOT skip the Agent spawn; do NOT combine executor and reviewer
roles. The reviewer runs in its own context window with no access to your
conversation history. You must give it everything it needs in the prompt.

Compose the reviewer prompt as follows. Use this template VERBATIM — fill in
the `[bracketed parts]` but do NOT omit, reorder, or summarize any section,
especially the review criteria (all 6 dimensions, all 29 items) and the
response format. The reviewer must evaluate every numbered item individually.

```
You are a pipeline reviewer for an ATLAS heavy-ion dimuon analysis. You review
multi-step bash pipeline scripts that orchestrate end-to-end physics analyses
on a Condor batch system with ROOT/RDF processing.

## Your task
Review the following pipeline and return a structured verdict.

## Task description
[paste the original $ARGUMENTS]

## Files to review
[list the pipeline script(s) with full paths — the reviewer will Read them]
[list any auxiliary files: .sub files, shell scripts called by the pipeline]

## Stage map
[paste the numbered stage list from Step 1, mapping each stage to its analysis function]

## Pipeline documentation
[paste the relevant section from Analysis/docs/ or Analysis/README.md that
describes the expected stages for this analysis pipeline. If no documentation
exists, state "No pipeline documentation found."]

## Review criteria

Evaluate EVERY item below. For each, state PASS or FAIL with specific evidence
(line numbers, function names, file paths). If an item is not applicable, state N/A.

### Dimension 1: Stage completeness

1. **All documented stages present**: every stage in the pipeline documentation
   (or Analysis/README.md) has a corresponding implementation in the script.
   List any missing stages.
2. **Stage ordering correct**: stages run in the correct dependency order.
   No stage reads an output that hasn't been produced yet.
3. **Input-output chaining**: each stage's expected input files match the
   previous stage's actual output files (filenames, paths, naming conventions).
4. **Year/sample coverage**: if the pipeline loops over years or samples, verify
   all expected values are covered and the loop variables are consistent
   (QUEUE_COUNTS, filename patterns, RDF configurations).
5. **Final outputs sufficient**: the pipeline's final outputs are enough to
   produce the physics result (plots, fitted parameters, tables for the note).
   Identify any gaps.

### Dimension 2: Safety & error handling

6. **set -Eeuo pipefail**: script uses strict bash mode.
7. **ERR trap**: errors are caught with context (line number, command, exit code).
8. **Per-stage validation**: every stage that produces ROOT files validates them
   before the next stage uses them (file exists, non-empty, non-zombie,
   has expected keys/trees).
9. **Fail-fast on bad output**: if validation fails, the script exits
   immediately with a clear error message naming the stage and the file.
   It must NOT continue to the next stage.
10. **No silent overwrite of good outputs**: check whether `rm -f` + `hadd -f`
    or similar patterns could destroy valid previous outputs when a re-run
    partially fails. Specifically: if hadd's input files are bad but the
    combined file from a previous run exists, does the pipeline delete the
    good combined file before validating the inputs?
11. **Condor failure handling**: held jobs are detected and reported (not
    silently ignored). Cluster IDs are parsed with error checking.
12. **Timeout support**: long-running waits have configurable timeouts.
13. **Clear error messages**: every `fail` call includes both WHAT failed
    (stage name, file path) and WHY (exit code, validation detail).

### Dimension 3: Idempotency & resumability

14. **Skip mechanisms**: the pipeline supports skipping completed expensive
    stages (e.g., SKIP_CONDOR for NTuple processing). Document what skip
    flags exist and whether additional ones are needed.
15. **Safe re-run**: running the pipeline twice with the same inputs produces
    the same outputs without corruption. Check for append-vs-overwrite issues.
16. **Partial re-run**: can individual stages be re-run without re-running
    everything? If not, flag as WARNING (not CRITICAL — some pipelines
    are designed as all-or-nothing).

### Dimension 4: Environment & resource handling

17. **Required commands checked**: all external tools (root, hadd, condor_submit,
    condor_q) are checked with `require_cmd` or equivalent before first use.
18. **Environment setup**: the script sources the ROOT environment reliably
    (handles set -eu conflicts with setup scripts).
19. **Configurable parameters**: hardcoded values that should be configurable
    are exposed as env vars with defaults (thread counts, poll intervals,
    timeout, year list).
20. **Working directory safety**: the script uses absolute paths or
    pushd/popd correctly. No reliance on implicit cwd.

### Dimension 5: Documentation

21. **Header documents all stages**: the script header lists every stage the
    pipeline performs, numbered, with a brief description.
22. **Usage section**: the header shows how to run the script and lists all
    env vars with their defaults and meanings.
23. **Pipeline doc exists**: a corresponding document in `Analysis/docs/`
    describes this pipeline's physics purpose, stage-by-stage procedure,
    expected inputs/outputs, and how to interpret results.
24. **Doc matches implementation**: the documentation accurately reflects the
    current script. No stale references to removed stages or renamed files.
25. **Agent-readable**: an agent encountering this pipeline for the first time
    can determine from the script + docs: (a) what physics result it produces,
    (b) what each stage does and what it outputs, (c) how to re-run a
    specific stage, (d) how to verify the outputs are correct.

### Dimension 6: Physics integration (pipeline-level)

These checks verify physics correctness at the ORCHESTRATION level — not
the internal C++ logic (that's /review-analysis-code's job).

26. **Correct mode/trigger settings**: RDF and NTuple stages use the right
    trigger_mode, mindR, hist_filling_cycle for this pipeline's physics goal.
    Cross-reference with DatasetTriggerMap conventions.
27. **Validation thresholds physically motivated**: magic numbers in validation
    (e.g., "at least 10 TF1s", "at least 4 TH2D") correspond to the actual
    binning scheme. State what the expected counts should be and whether they
    match.
28. **Stage-to-stage physics consistency**: the same physics assumptions
    (trigger mode, sign convention, year-specific cuts) are used consistently
    across all stages. No stage uses different settings than what the previous
    stage produced.
29. **Final result chain complete**: trace the full chain from raw input to
    final physics plot/number. Identify any step where a physics quantity is
    computed but never validated or used.

### Anti-patterns to catch
- `rm -f combined_output` BEFORE validating per-batch inputs
- Validation that only checks file existence (not ROOT integrity)
- `set +e` without restoring `set -e` afterward
- Hardcoded year lists that miss new data-taking periods
- ROOT macros embedded in heredocs without proper escape handling
- `pushd`/`popd` mismatches (especially in error paths)
- Condor wait loops that silently continue when `condor_q` fails
- Pipeline stages that share output filenames (stage N overwrites stage M output)

## Response format

You MUST structure your response EXACTLY as follows:

## Findings by Dimension

### Dimension 1: Stage completeness
[numbered findings for items 1-5, each PASS/FAIL/N/A with evidence]

### Dimension 2: Safety & error handling
[numbered findings for items 6-13]

### Dimension 3: Idempotency & resumability
[numbered findings for items 14-16]

### Dimension 4: Environment & resource handling
[numbered findings for items 17-20]

### Dimension 5: Documentation
[numbered findings for items 21-25]

### Dimension 6: Physics integration
[numbered findings for items 26-29]

## Issues
[consolidated numbered list of all FAIL findings, or "None found."]

For each issue:
N. [DIMENSION] Description of the problem
   File: [path:line]
   Severity: CRITICAL | WARNING | INFO
   Fix: [specific instruction for what to change]

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

Read the reviewer's response. First verify it follows the required format:
- Has a `## Findings by Dimension` section with per-item PASS/FAIL/N/A
- Has a `## Issues` section
- Has a `VERDICT:` line

If the reviewer did not follow the format (e.g., gave a summary instead of
per-item findings), note this in the log and re-spawn the reviewer with the
same prompt. Do NOT accept a verdict without per-item evidence.

Extract the verdict by searching for the line starting with `VERDICT:`.

**MANDATORY**: Append the following to the log file. Paste the reviewer's
full output — do NOT summarize or paraphrase. The log is the audit trail.

```
## Iteration N
**Reviewer verdict**: PASS | FAIL
**Issues found**: [count]
**Reviewer full output**:
[paste the reviewer's ENTIRE response verbatim — all 6 dimension sections,
Issues, What Passed, and Verdict. Do not truncate or summarize.]
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

**MANDATORY — do ALL four steps in order. Do NOT skip any step, even if you
are a subagent. The log file and tracking event are required for auditability.**

1. Update log:
```
**Status**: APPROVED at iteration [N]
**Summary**: [1-2 sentences]
**Dimension scorecard**:
[for each dimension: "Dim N (name): PASS — X/Y criteria met"]
```

2. **Emit end event** — append to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"review-pipeline","status":"approved","iterations":<N>,"timestamp":"<ISO 8601>","criticals":0,"warnings":0,"slug":"<slug>"}
```

3. Verify the log file contains the full reviewer output for every iteration
   (not just summaries). If any iteration is missing the verbatim reviewer
   output, append it now.

4. Delete `.claude/logs/active-session.txt`, then present final output to user.

## Exit: Escalate

**MANDATORY — do ALL four steps in order. Do NOT skip any step, even if you
are a subagent.**

1. Update log:
```
**Status**: ESCALATED at iteration [N] (max iterations reached)
**Unresolved issues**:
- [Issue 1]: [what it is, what was tried, why unresolved]
**Recommended action**: [specific suggestions]
```

2. **Emit end event** — append to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"review-pipeline","status":"escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<count>,"warnings":<count>,"unresolved":["<short issue description>",...],"slug":"<slug>"}
```

3. Verify the log file contains the full reviewer output for every iteration.

4. Delete `.claude/logs/active-session.txt`, then present unresolved issues to user with:
"These issues need your judgment: [list with context]"
