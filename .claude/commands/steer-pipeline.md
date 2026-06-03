---
description: >
  Autonomous pipeline runner: pre-flight source-change detection, execution
  monitoring, error recovery, post-run output review, and bookkeeping
  (docs + git). Use when the user says "run this pipeline", "steer",
  "execute pipeline", or wants autonomous end-to-end pipeline execution
  with physics-aware guardrails.
---

You are the **steering agent** for an analysis pipeline. You autonomously
run a pipeline script with pre-flight safety checks, monitor execution,
recover from errors, review outputs, and handle bookkeeping.

Task: $ARGUMENTS

## Argument parsing

Parse `$ARGUMENTS` for:
- `<script-path>` (required): path to the pipeline script, relative to
  `Analysis/pipelines/` or absolute. Resolve to full path.
- `[pipeline-args...]`: any positional args forwarded to the pipeline
  (e.g., `hijing` for overlay, `--skip-mixing` for mixed pairs)
- `--from-stage N` (optional): skip stages before N using native skip flags
- `--full-review` (optional): force full plot review loop on outputs
- `--dry-run` (optional): run pre-flight only, do not execute
- `KEY=VALUE` pairs: passed as env vars to the pipeline

## Setup

Read these files before starting:
- The pipeline script (fully)
- `Analysis/README.md` (class hierarchy and pipeline structure)
- `Analysis/docs/` — the pipeline doc relevant to this script
- `.claude/kb/index.md` (check for relevant KB articles)

**MANDATORY SETUP**: Complete ALL steps below BEFORE entering Phase 1.
If you are a subagent spawned by a parent agent, these still apply to you.

### Initialize run log

1. Derive a slug from the pipeline name (e.g., `pbpb-trig-eff`).
   If pipeline args exist, append them (e.g., `pythia-overlay-hijing`).
2. Set `RUN_LOG = .claude/logs/pipeline-runs/<slug>-YYYYMMDD-HHMMSS.log`
3. Set `OUTPUT_LOG = .claude/logs/pipeline-runs/<slug>-YYYYMMDD-HHMMSS-output.log`
4. Initialize the run log:

```
# Pipeline Run Log
**Pipeline**: <script basename>
**Args**: <pipeline args, or "none">
**Task**: $ARGUMENTS
**Started**: [timestamp]
**Status**: IN PROGRESS
**Git SHA**: [current HEAD sha — run `git rev-parse HEAD`]
**Phases completed**: []
**Recovery attempts**: 0
**Dependencies**: [populated in Phase 1]
```

### Emit tracking event

Append to `.claude/logs/tracking.jsonl`:
```json
{"event":"start","command":"steer-pipeline","pipeline":"<script>","task":"$ARGUMENTS","slug":"<slug>","timestamp":"<ISO 8601>","git_sha":"<HEAD>"}
```

## Compaction recovery

If you have lost context (after compaction or conversation compression):
1. `ls -t .claude/logs/pipeline-runs/*.log | head -1` → latest run log
2. Read it fully — it is ground truth. Resume from the last completed phase.
3. If Phase 2 was in progress, check for the output log and whether the
   pipeline process is still running (`ps aux | grep <script-name>`).

## Reference tables

### Cross-pipeline dependency map

Used for downstream warnings in Phase 5:
- `pipeline_pbpb_trig_eff.sh` → `pipeline_pbpb_crossx.sh`
  (trig eff corrections feed crossx nominal processing)
- `pipeline_powheg_fullsim_single_muon.sh` → `pipeline_powheg_fullsim_mixed_pairs.sh`
  (single muon combined trees used by mixer)
- `pipeline_pythia_truth.sh` ↔ `pipeline_pythia_fullsim_overlay.sh`
  (truth is reference for validating overlay results)

### Native skip flag map

Map `--from-stage N` to native pipeline skip mechanisms:
| Pipeline | Flag | Stages skipped |
|---|---|---|
| `pipeline_pbpb_trig_eff.sh` | `SKIP_CONDOR=1` | 1-4 (submit, wait, validate batch, hadd) |
| `pipeline_pbpb_crossx.sh` | `SKIP_CONDOR=1` | 1-4 (submit, wait, validate batch, hadd) |
| `pipeline_pythia_truth.sh` | `SKIP_CONDOR=1` | 1-3 (submit, wait, validate batch) |
| `pipeline_pythia_fullsim_overlay.sh` | `--dry-run` | 1-2 (submit, wait) |
| `pipeline_powheg_fullsim_single_muon.sh` | `--skip-condor` | 1-2 (submit, wait) |
| `pipeline_powheg_fullsim_mixed_pairs.sh` | `--skip-mixing` | mixing stage only |

If `--from-stage N` targets a stage that the native flags can skip, set the
appropriate flag. If the target stage has no native skip mechanism, warn the
user: "Cannot skip to stage N directly; the nearest skip point is stage M.
Proceed from stage M?" If user agrees, set the flag and continue.

---

## Phase 1: Pre-flight

### 1.1 Extract dependency list

Read the pipeline script and extract ALL referenced source files:

| Pattern | File type |
|---|---|
| `.L <file>` or `.L <file>+` | C++ source (ACLiC) |
| `root -l -b -q "<file>(...)"` | ROOT macro |
| `condor_submit "<file>.sub"` | Condor submit file |
| `source <file>` or `. <file>` | Sourced shell script |
| `bash <file>` or `./<file>` | Called shell script |
| Variables like `*_DIR`/`*_BASE` + filename | Computed paths |

For each `.cxx` or `.c` file, check if a corresponding `.h` file exists
and add it. Also add any helper files referenced by Condor `.sub` files
(e.g., the `executable` line in `.sub` files).

Resolve all paths using the script's directory variables (`SCRIPT_DIR`,
`NTP_DIR`, `RDF_DIR`, `PLOT_DIR`, etc.).

Log the full dependency list to the run log under `**Dependencies**:`.

### 1.2 Check source changes

1. Read run history: `.claude/logs/pipeline-runs/<slug>.jsonl`
2. Find last entry with `"outcome":"success"` → extract `git_sha`
3. If no previous run → log "First run, skipping diff" → skip to 1.5
4. Run: `git diff --stat <last-sha>..HEAD -- <dep1> <dep2> ...`
5. Also check: `git diff --stat <last-sha>..HEAD -- <pipeline-script>`
6. If no changes → log "No source changes since last run (SHA <last-sha>)" → skip to 1.5

### 1.3 Classify source changes

For each changed file, read the full diff with `git diff <last-sha>..HEAD -- <file>`.

Classify every hunk:

**SAFE** (log and continue automatically):
- Whitespace, formatting, indentation changes
- Comment-only changes (not physics-documenting comments)
- Logging, print, echo statement changes
- Variable renaming with identical semantics
- Error message text changes

**PHYSICS-AFFECTING** (block — present to user and ask):
- Cut values, selection criteria, filter conditions, veto logic
- Binning definitions (bin edges, nbins, range)
- Trigger logic: trigger_mode, mindR, trigger matching, trigger names
- Physics formulas: efficiency calculations, corrections, weights, scale factors
- Sign conventions, charge selection, OS/SS logic
- Sample paths, input file paths, dataset definitions
- Any change inside a function that computes a physics quantity
- Changes to Condor `.sub` files that alter job count or input file lists

**AMBIGUOUS** (block — present to user and ask):
- Structural refactors that might change execution order
- New or removed function calls in the main pipeline flow
- Changes to validation thresholds
- Anything you're not confident classifying

For SAFE changes: log file name and one-line summary, continue.

For PHYSICS-AFFECTING or AMBIGUOUS:
1. Present a summary table to the user:

```
| File | Classification | What changed | Affected stages |
|---|---|---|---|
```

2. For each physics-affecting change, explain: what the old value was,
   what the new value is, and which final physics outputs might change.

3. Ask: **"Proceed with pipeline execution? (yes / abort / review first)"**
   - yes → continue to Phase 1.4
   - abort → exit with status "aborted-preflight"
   - review first → suggest the appropriate `/review-*` command and exit

### 1.4 Check pipeline-docs consistency

If the pipeline script itself changed since the last run:
1. Read the pipeline's doc file in `Analysis/docs/`
2. Quick check: does the script's header stage list match the doc?
3. If mismatch → log "WARNING: pipeline header and docs out of sync"
   (Phase 5 will reconcile)

### 1.5 Environment validation

Run these checks — fail fast before touching data:

1. **ROOT**: `which root` succeeds. If not, run `source ~/setup.sh` and retry.
2. **hadd**: `which hadd` succeeds.
3. **Condor** (if pipeline uses `condor_submit`): `condor_q > /dev/null 2>&1`
   succeeds. If Condor is down, and `--from-stage` can skip Condor stages,
   suggest that. Otherwise fail.
4. **Disk space**: `df -BG <output-directory> | awk 'NR==2{print $4}'` → warn
   if < 10 GB free.
5. **Syntax check**: `bash -n <pipeline-script>` → if fails, report the syntax
   error and exit.

Log: `"Pre-flight complete. [N] deps checked, [M] changed ([K] safe, [J] physics). Environment OK."`

If `--dry-run` → update run log status to `"DRY RUN COMPLETE"`, emit end
tracking event, exit with the pre-flight report. Do not proceed to Phase 2.

---

## Phase 2: Execute + Monitor

### 2.1 Configure execution

1. If `--from-stage N` was set: apply the native skip flag from the map.
2. Collect all env var overrides from arguments.
3. Build the full command line:
   ```
   <ENV_VARS> bash <script-path> [pipeline-args] 2>&1 | tee <OUTPUT_LOG>
   ```

### 2.2 Run pipeline

**Determine execution mode:**
- If the pipeline has Condor stages AND those stages are NOT skipped →
  use `run_in_background: true` (pipelines with Condor can take hours).
- If Condor is skipped (SKIP_CONDOR, --skip-condor, --dry-run) →
  run in foreground with `timeout: 600000` (10 min).

Execute the command. Log: `"Phase 2: pipeline started at [timestamp]"`

### 2.3 On completion

1. Read exit code from Bash tool result.
2. Read the output log (`<OUTPUT_LOG>`) fully.
3. Parse stage transitions — search for lines matching patterns like:
   `Stage N:`, `Step N:`, `=== Stage`, `--- Stage`, or the pipeline's
   specific stage markers.
4. Build a stage completion list.
5. Append to run log:
   ```
   **Phase 2 result**: exit_code=[N]
   **Stages completed**: [list]
   **Output log**: <OUTPUT_LOG path>
   ```

- Exit 0 → log "Pipeline completed successfully" → go to **Phase 4**
- Exit non-zero → log "Pipeline FAILED at stage [N]" → go to **Phase 3**

---

## Phase 3: Error Recovery

**Entry**: pipeline exited non-zero.
**Recovery budget**: max 2 auto-recovery attempts per invocation.
**Physics guardrail**: NEVER auto-fix anything that changes physics outputs.
When in doubt, ask the user.

### 3.1 Data corruption check — DO THIS FIRST

Before ANY fix attempt, verify the integrity of existing outputs:

1. From the pipeline script, identify all output files that are produced
   stage-by-stage. Focus on "combined" or "final" outputs (hadd results,
   RDF histogram files, fitted parameter files, plots).

2. For outputs from stages that COMPLETED SUCCESSFULLY in this run:
   - Verify files exist and are valid ROOT files (non-zombie)
   - These are fine — the pipeline validated them before proceeding

3. For outputs from PREVIOUS SUCCESSFUL RUNS that this pipeline MIGHT
   have touched before failing:
   - Check the pipeline for `rm -f` or `hadd -f` patterns that delete
     old outputs before creating new ones
   - If the pipeline deleted a previous good output but failed to create
     the replacement → **DATA LOSS**

4. **If data loss detected**:
   - Log: `"DATA CORRUPTION: [file] was deleted but replacement failed"`
   - Report to user immediately with full details
   - Do NOT attempt auto-recovery — the user must decide how to restore
   - Exit with status "data-loss"

5. If no data loss → proceed to 3.2

### 3.2 Root cause analysis

Read the output log tail (last 50 lines) plus any lines containing
"ERROR", "FAIL", "Segmentation", "error:", or "fatal".

Classify the error:

#### Auto-fixable (TRANSIENT / TECHNICAL — no physics impact)

| Error | Auto-fix |
|---|---|
| Condor job held | `condor_release <cluster_id>`, wait for completion |
| Condor job evicted/removed | Resubmit the same `.sub` file |
| Condor timeout | Double `CONDOR_TIMEOUT_SECONDS`, resubmit |
| ROOT segfault on one batch | Log the bad file, skip it if pipeline supports `--skip-batch` |
| ACLiC stale `.pcm`/`.d` | `rm *.pcm *.d *.so` in the source dir, retry |
| `source ~/setup.sh` failure | Retry once; if still fails, report env issue |
| Disk quota warning (not full) | Report usage, suggest cleanup, wait for user |
| Network/CVMFS transient | Wait 60s, retry |
| `hadd` failure on one file | Check if the input ROOT file is zombie, skip it |

#### Needs user input (LOGIC / PHYSICS)

| Error | Why user input needed |
|---|---|
| Compilation error in physics code | Fix might change physics logic |
| Missing input file from another pipeline | Cross-pipeline dependency — user decides which to re-run |
| Validation failure (wrong # of histograms) | Might indicate binning or config change |
| Wrong ROOT macro arguments | Macro interface might have changed intentionally |
| Any error you cannot confidently classify | Safety default |

#### Unknown

Any error not matching the above → report full context to user.

### 3.3 Apply fix and re-run

**For auto-fixable errors:**

1. Apply the fix (release Condor, clean stale files, etc.)
2. Determine which stages completed successfully (from 2.3 stage list)
3. Set native skip flags to skip completed stages
4. Append to run log:
   ```
   **Recovery attempt [N]**: [error description]
   **Auto-fix**: [what was done]
   **Re-running from**: stage [M] (skip flags: [flags])
   ```
5. Increment recovery counter
6. Go back to **Phase 2.2** (re-run with skip flags)

**If recovery counter > 2:**
- Log: "Auto-recovery exhausted after 2 attempts"
- Report to user: "Pipeline failed 3 times. Errors: [list all 3].
  Auto-recovery could not resolve this. Please investigate."
- Exit with status "recovery-exhausted"

**For errors needing user input:**
- Present: the error, your root cause analysis, the suggested fix, and
  explicitly state whether the fix has physics implications
- Wait for user decision
- If user provides a fix → apply it, re-run
- If user says abort → exit

---

## Phase 4: Post-run Review

**Entry**: pipeline completed successfully (exit 0).

### 4.1 Output file validation (always)

Read the pipeline script to identify ALL expected output files. Look for:
- Paths in validation stages (e.g., `validate_root_file_quick "$FILE"`)
- Final output assignments (e.g., `COMBINED=...`, `HIST_FILE=...`)
- Plot output paths (e.g., plot macros that `SaveAs` or `Print`)

For each expected ROOT output file, run a validation ROOT macro:

```cpp
{
  auto f = TFile::Open("<path>");
  if (!f || f->IsZombie()) { cout << "FAIL:zombie" << endl; exit(1); }
  auto keys = f->GetListOfKeys();
  cout << "OK: " << keys->GetEntries() << " keys" << endl;
  for (int i = 0; i < TMath::Min(keys->GetEntries(), 5); i++)
    cout << "  key: " << keys->At(i)->GetName() << endl;
}
```

For each expected plot file (`.png`):
- Check exists and file size > 0

Log all validation results. If any file is missing or invalid → report as
WARNING to user (pipeline passed its own validation, so this is unexpected).

### 4.2 Output regression check (if previous successful run exists)

1. Read run history → find last successful run's `output_files` list
2. Compare against current run's output files:

| Check | Threshold | Severity |
|---|---|---|
| File present last time, missing now | any | WARNING |
| File > 50% smaller than last time | 50% | WARNING |
| New file not in previous run | any | INFO (log only) |

3. For ROOT files containing histograms, open both current and previous
   versions and compare a sample of histograms:

```cpp
{
  auto f1 = TFile::Open("<current>");
  auto f2 = TFile::Open("<previous>");
  // For each TH1/TH2 key in f1:
  //   compare Integral, GetMean, GetRMS
  //   flag if relative difference > 10%
}
```

Only do this comparison for the pipeline's final output files (not
intermediate files that may legitimately change).

4. If regressions found:
   - If source changes were detected in Phase 1 → report regressions as
     INFO: "Output changed, likely due to source modifications: [list]"
   - If NO source changes → report as WARNING: "Output changed with no
     source modifications — possible non-determinism or environment change: [list]"

### 4.3 Plot review (conditional)

**Trigger full plot review** if ANY of:
- `--full-review` flag was set
- This is the first run (no previous successful run in history)
- PHYSICS-AFFECTING source changes were detected in Phase 1

**If triggered:**
1. Identify the pipeline's output plot directory from the script (look for
   the plotting stage's output path — `PLOT_DIR`, `SaveAs` calls, etc.)
2. Spawn a plot reviewer subagent using the **Agent tool** with
   `subagent_type: "general-purpose"`. Give it:
   - The plot directory path
   - The pipeline's physics goal (from the doc)
   - What the plots should show (from the pipeline's plotting stage)
   - Instructions to check: axis labels, units, legend entries, scale
     (log vs linear), empty bins, unexpected features, and whether the
     plots are consistent with each other
3. Parse the reviewer's assessment
4. Log the full review result

**If NOT triggered:**
- Log: "Plot review skipped (no triggering condition). Previous run at
  SHA [sha] was reviewed."

---

## Phase 5: Bookkeeping

### 5.1 Append to run history

Append one JSON line to `.claude/logs/pipeline-runs/<slug>.jsonl`:

```json
{
  "timestamp": "<ISO 8601>",
  "git_sha": "<HEAD at start>",
  "outcome": "success | failed | recovered | aborted-preflight | data-loss | recovery-exhausted",
  "stages_run": [1, 2, 3],
  "skip_flags": ["SKIP_CONDOR"],
  "args": "<pipeline args>",
  "errors": [],
  "fixes": [],
  "recovery_attempts": 0,
  "review": "lightweight | full | skipped",
  "regressions": [],
  "output_files": ["<paths>"],
  "duration_s": 1234,
  "dependencies": ["<source file paths>"],
  "run_log": "<RUN_LOG path>"
}
```

### 5.2 Doc update

If the pipeline script was modified during this session (error recovery,
skip flag changes, any edits):

1. Read the pipeline's doc in `Analysis/docs/`
2. Check: does the doc accurately describe the current script?
   - Stage list matches?
   - Env vars / flags documented?
   - Input/output paths correct?
3. If updates needed → apply them
4. Log what was updated

### 5.3 Git commit + push

If ANY files were modified during this session (pipeline script, docs,
helper scripts):

1. `git diff --name-only` to list changes
2. Stage ONLY: pipeline scripts, docs, helper scripts, `.sub` files.
   Do NOT stage: output ROOT files, plots, log files, `.jsonl` files.
3. Commit:
   ```
   steer: <slug> — <1-line summary>

   Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
   ```
4. `git push`
5. Log the commit SHA

If no files modified → skip.

### 5.4 Downstream dependency warning

Check the cross-pipeline dependency map in the Reference Tables section.

If this pipeline has downstream dependents AND this run's outcome was
`"success"` or `"recovered"`:

1. Check if output files changed compared to the previous successful run
   (from Phase 4.2 results)
2. If outputs changed (or this is the first run):
   - Log and report to user:
     ```
     Downstream pipeline may need re-running:
       <downstream-pipeline> depends on <this-pipeline>'s outputs.
       Run: /steer-pipeline <downstream-script>
     ```
3. If outputs are unchanged → log "Downstream pipelines unaffected."

---

## Exit

Update the run log:
```
**Status**: <COMPLETE | FAILED | ESCALATED | ABORTED>
**Outcome**: <success | failed | recovered | data-loss | recovery-exhausted | aborted-preflight>
**Phases completed**: [1, 2, 4, 5]
**Recovery attempts**: [N]
**Duration**: [total elapsed time]
**Summary**: [1-2 sentences: what ran, what was found, what needs attention]
```

Emit end event to `.claude/logs/tracking.jsonl`:
```json
{"event":"end","command":"steer-pipeline","pipeline":"<script>","status":"<outcome>","phases":[1,2,4,5],"recovery_attempts":0,"regressions":0,"timestamp":"<ISO 8601>","slug":"<slug>"}
```

**Final report to user** (format depends on outcome):

**Success:**
```
Pipeline <name> completed successfully.
- Stages: [list]
- Outputs: [key output files]
- Review: [lightweight/full — summary]
- [Downstream warning if any]
```

**Recovered:**
```
Pipeline <name> completed after [N] recovery attempt(s).
- Error: [what failed]
- Fix: [what was auto-fixed]
- Outputs: [validated — summary]
```

**Failed / Escalated:**
```
Pipeline <name> needs your input:
- Error: [what failed and why]
- Attempted: [what was tried]
- Recommendation: [specific next step]
```
