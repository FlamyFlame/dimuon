---
description: >
  Executor-reviewer loop for reviewing Claude Code plugins and skills.
  Use when the user asks to review, audit, validate, or check a plugin,
  or says "review the plugin", "check the plugin for portability",
  "is this plugin ready to share", "audit the skill", or any plugin/skill
  quality check. Also auto-triggered after plugin creation or modification.
---

You are the **executor** in a structured executor-reviewer loop.

**CRITICAL CONSTRAINT:** You must NOT evaluate the plugin checklist yourself.
Your job is to prepare context, then delegate ALL quality evaluation to a
separate reviewer subagent via the **Agent tool**. You fix issues; the
reviewer finds them. If you catch yourself applying the checklist criteria
to the plugin files, STOP — that is the reviewer's job. Proceed to Step 2
and spawn the subagent.

Task: $ARGUMENTS

## Setup

Read these files before starting:
- `.claude/agents/executor.md` (your behavior as executor)
- `.claude/docs/command-design.md` (loop architecture)

Then identify:
1. The plugin directory (from $ARGUMENTS or ask user)
2. Any reference tracking docs or procedure docs the plugin encodes
   - Check `## Active Tracking Docs` in CLAUDE.md for relevant docs
   - If the user provided reference docs, note them

List all files in the plugin directory (for the reviewer prompt).
Read the plugin's `plugin.json` to confirm the plugin name and description.

Configuration:
- MAX_ITERATIONS = 5

### Initialize log file

1. Derive a short slug from the task (3-5 words, kebab-case).
2. Set `LOG_FILE = .claude/logs/review-plugin-YYYYMMDD-HHMMSS-<slug>.md`
3. Write the log file path to `.claude/logs/active-session.txt` (overwrite).
4. Initialize the log file:

```
# Plugin Review Log
**Task**: $ARGUMENTS
**Plugin**: <plugin directory path>
**Log file**: <LOG_FILE basename>
**Started**: [timestamp]
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5
```

### Emit tracking event

Append to `.claude/logs/tracking.jsonl`:
```json
{"event":"start","command":"review-plugin","task":"$ARGUMENTS","slug":"<slug>","timestamp":"<ISO 8601>"}
```

## Compaction recovery

If you have lost context (e.g., after compaction), recover by:
1. Read `.claude/logs/active-session.txt` to get the log file path.
2. If that file is missing, run `ls -t .claude/logs/review-plugin-*.md | head -1`.
3. Read the log file fully — it is ground truth. Resume from where it left off.

## Loop

### Step 1: Prepare or amend

**First iteration:** Do NOT read or evaluate the plugin files yourself.
Just gather the file list and reference doc paths — the reviewer will read
and evaluate them. Proceed immediately to Step 2.

**Subsequent iterations (after reviewer FAIL):** Address ONLY the specific
issues the reviewer flagged. Edit the plugin files to fix them. Do not
refactor or change anything not flagged. Record what you changed, then
proceed to Step 2.

### Step 2: Spawn reviewer subagent (MANDATORY — use the Agent tool)

You MUST use the **Agent tool** to spawn a reviewer subagent with
`subagent_type: "general-purpose"`. Do NOT skip this step. Do NOT
evaluate the checklist yourself instead of spawning. The reviewer runs in
its own context window with fresh eyes — that is the entire point.

Compose the reviewer prompt as follows (fill in the bracketed parts):

```
You are a plugin reviewer for Claude Code plugins. Your job is to
independently evaluate a plugin for portability, procedure fidelity,
self-sufficiency, and operational robustness.

## Task description
[paste the original $ARGUMENTS]

## Plugin directory
[full path to plugin root]

## Files to review
[list every file in the plugin, with full paths — you will Read each one]

## Reference documents
[list any tracking docs, procedure docs, or external references that the
 plugin encodes. Include full paths. If none provided, note that and check
 whether the plugin is self-contained enough to not need them.]

## Walkthrough scenario

Use this concrete scenario to trace through each skill step by step. For
every step, check: does the skill tell you exactly what to run, what output
to expect, what to do if it fails, and when to proceed? If you get stuck
at any step — that is a SELF-SUFFICIENCY or ROBUSTNESS issue.

[Paste the worked example from the plugin's README if one exists. If not,
 paste a scenario from the tracking doc's pilot test, with personal
 identifiers replaced by generic placeholders. If neither exists, flag
 the absence as a SELF-SUFFICIENCY WARNING: "no walkthrough scenario
 available — standalone usability cannot be concretely verified."]

When applying the walkthrough:
- Mentally execute each skill using the scenario's values as arguments.
- At each step, ask: "would I know what to do next with ONLY the skill
  text and this scenario?" If not, that is a gap.
- Check that error paths mentioned in the scenario (e.g., quota failures,
  DID conflicts) are handled by the skill.

## Review criteria

### Portability (no hardcoded environment)

1. **No hardcoded usernames**: no literal usernames in paths, scope prefixes,
   Rucio accounts, or examples that would break for another user. Must derive
   from runtime (e.g., `rucio whoami`, `$USER`, argument).
   - Exception: BNL-specific plugins may assume generic BNL SDCC infrastructure
     paths (`/cvmfs/atlas.cern.ch/`, `/pnfs/usatlas.bnl.gov/`) but NOT personal
     paths (`/usatlas/u/yuhanguo/`, `~/dcachearea`).
2. **No hardcoded file paths**: no absolute paths to specific data files, output
   directories, or working areas. Paths must come from arguments or be
   constructed at runtime.
3. **No hardcoded dataset names or DIDs**: Rucio dataset names, rule IDs,
   container names must be arguments or derived, not baked in.
4. **No hardcoded credentials or tokens**: no proxy paths, certificate paths,
   or auth tokens.
5. **Environment assumptions documented**: if the plugin requires a specific
   environment (e.g., BNL SDCC, CVMFS), the README states this clearly.

### Procedure fidelity

6. **Faithful to source procedure**: every step in the plugin's skills matches
   the procedure in the reference documents. No steps omitted, added, or
   reordered in ways that change semantics. If no reference docs provided,
   check internal consistency only and flag that fidelity could not be fully
   verified.
7. **Error handling matches source**: failure behavior at each step matches
   what the source procedure specifies.
8. **Key facts preserved**: constants, RSE names, path prefixes, timing
   expectations, quota limits from the source are present and correct.
9. **Rollback procedure included**: if the procedure has destructive steps,
   rollback is documented.

### Self-sufficiency

10. **Standalone usability**: a user installing this plugin with zero additional
    context can follow every skill end-to-end. All commands, flags, expected
    outputs, and decision criteria are in the skill files.
11. **No implicit knowledge**: the plugin does not assume the user knows
    project-specific conventions unless documented in the plugin itself.
12. **Arguments fully specified**: every skill with arguments documents what
    each argument is, its format, and provides an example. Missing-argument
    behavior is defined.
13. **README covers full workflow**: README explains purpose, skill ordering,
    prerequisites, and troubleshooting.
14. **README entry point is obvious**: a new user knows immediately which
    skill to run first, whether one skill handles everything or they must
    orchestrate multiple skills, and when to use the other skills. Read the
    README as a new user: would you know what to type after installation?
15. **Worked example in README**: concrete end-to-end walkthrough with
    generic placeholders, expected output at each step, common failure
    modes. Must cover all major paths (e.g., simple mode and advanced mode).
    Mentally trace through every skill step — if you get stuck, that's a gap.

### Skill quality (for skill-based plugins)

16. **Description triggers correctly**: each skill's `description:` field
    contains enough specific terms to trigger correctly and avoid false
    positives.
17. **allowed-tools minimal**: each skill lists only the tools it actually needs.
18. **Argument hints usable**: `argument-hint` clearly shows expected format.
19. **disable-model-invocation appropriate**: set to `true` for purely
    procedural skills; `false` or absent for skills requiring judgment.
20. **Inter-skill references correct**: cross-references to other skills use
    correct names and exist.

### Operational robustness

21. **Idempotency or guard checks**: state-modifying skills check for existing
    state or document re-run behavior.
22. **Destructive actions flagged**: irreversible steps and side-effects on the
    user's working state (git stash, git checkout) are clearly marked. The
    skill must ask the user before touching uncommitted work and offer
    alternatives (e.g., commit first vs stash).
23. **Long-running operations handled**: steps that can take hours document
    expected timing, monitoring, and stuck-vs-normal criteria.
24. **Output verification**: automatic verification (entry counts, file counts,
    readability checks), not just "ask the user to check."
25. **No hardcoded git branch names**: rollback/checkout commands must detect
    the base branch at runtime, not assume `main` or `master`.
26. **Crash-safe state transitions**: multi-step mutations (rename + rebuild)
    must not leave a broken state if interrupted. Check for staging/atomic
    swap patterns.
27. **Cleanup on all exit paths**: side-effects (git stash, temp branches,
    temp dirs) cleaned up on both success and failure. If cleanup can fail
    (e.g., stash pop merge conflict), the skill handles or escalates.
28. **Thorough codebase search**: path/pattern searches must not rely on
    literal string grep alone. If no matches found, skill asks user or
    exits with clear manual instructions — never silently skips.
29. **Agent-appropriate instruction level**: skills should describe the goal
    and edge cases, not prescribe exact shell commands for tasks agents can
    figure out (searching, code reading, testing). Hardcoded commands are
    appropriate for non-obvious operations (Rucio flags, PFN stripping).
    Over-specified = agent can't adapt to the user's codebase.
    Under-specified = agent misses edge cases. Check the balance.
30. **Exhaustive failure-mode coverage**: using the worked example, walk
    through every step and enumerate what can go wrong. For each failure,
    check if the skill handles it (error message, fallback, rollback, or
    escalation). Flag unhandled failures. Include: external service failures,
    user-state conflicts, data integrity issues, mid-step interruptions.
31. **Maximal automation at minimal risk**: automate aggressively when safe
    (file checks, compilation, read-only scans, temp-dir output). Ask the
    user when risky (code edits, git stash, directory renames, file deletion).
    Under-automated (asks for something the agent could safely do) and
    under-protected (silently performs risky action) are both failures.

## Severity definitions

| Level | Meaning | Blocks PASS? |
|-------|---------|:------------:|
| CRITICAL | Plugin would fail, corrupt data, lose user work, crash mid-procedure, silently produce wrong results, or leave user in a broken state due to an unhandled failure. | Yes |
| WARNING | Confusing, suboptimal but functional, missing docs that don't prevent operation. | Yes |
| INFO | Minor suggestion, style preference. | No |

**CRITICAL triggers** (must be CRITICAL, not WARNING):
- Unhandled failure in multi-step procedure (crash leaves broken state)
- Risky action without user awareness (git stash, code edit w/o confirmation)
- Missing rollback for destructive step
- Hardcoded username/path that breaks for other users
- Silent skip on error instead of stopping
- Non-atomic multi-step state mutation
- Side-effect not cleaned up on a failure path

## How to review

1. Read every file in the plugin directory.
2. If reference documents are listed, read them and compare against the plugin's
   procedure step by step.
3. For each checklist item, evaluate with specific evidence.
4. Pay special attention to:
   - Literal usernames (grep for common patterns: `yuhang`, `yuhanguo`, specific usernames)
   - Absolute paths to personal directories
   - Examples that use specific values without marking them as examples
   - Steps in reference docs that are missing from the plugin
   - Skills that reference external docs the plugin user won't have

## Response format

You MUST structure your response EXACTLY as follows:

## Issues
[numbered list of issues, or "None found." if all pass]

For each issue:
N. [CATEGORY] Description of the problem
   File: [path:line or element name]
   Severity: CRITICAL | WARNING | INFO
   Fix: [specific instruction for what to change]

Categories: PORTABILITY, FIDELITY, SELF-SUFFICIENCY, SKILL-QUALITY, ROBUSTNESS

## What Passed
[brief list of criteria that passed, grouped by category]

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
{"event":"end","command":"review-plugin","status":"approved","iterations":<N>,"timestamp":"<ISO 8601>","criticals":0,"warnings":0,"slug":"<slug>"}
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
{"event":"end","command":"review-plugin","status":"escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<count>,"warnings":<count>,"unresolved":["<short issue description>",...],"slug":"<slug>"}
```

3. Delete `.claude/logs/active-session.txt`, then present unresolved issues to user with:
"These issues need your judgment: [list with context]"
