# Executor-Reviewer Command Design

Design criteria for the five `/review-*` commands that implement multi-agent
executor-reviewer loops in Claude Code.

## Architecture

```
Main agent (executor)
  - Holds full task context, code state, iteration history
  - Does the work: writes code, produces plots, writes LaTeX, investigates
  - Composes reviewer prompt with all context the subagent needs
  |
  |--- spawns ---> Reviewer subagent (via Agent tool)
  |                  - Fresh context window each iteration (no anchoring bias)
  |                  - Receives: task, files to review, criteria, output format
  |                  - Returns: structured verdict (PASS/FAIL) + issues
  |<-- result ------|
  |
  - Parses verdict, iterates or exits
```

Communication is parent-child only. Subagents cannot talk to each other.
The executor retains full history; the reviewer gets fresh eyes each iteration.

## Design Criteria Checklist

Every `/review-*` command must satisfy all of the following:

### 0. Trigger Description

- Multi-line YAML `description:` in frontmatter
- Lists concrete trigger phrases users might say
- Covers both direct requests ("make a plot") and indirect signals ("the plot
  looks wrong", "fix the axis")
- Enables auto-dispatch without the user typing `/command-name`

### 1. Role Assignment

- Main agent = executor (stated explicitly: "You are the **executor**")
- Reviewer = subagent spawned via `Agent` tool with `subagent_type: "general-purpose"`
- Reviewer prompt is self-contained: includes task, file paths, criteria, and
  output format (the subagent has no conversation history)

### 2. Review Criteria

- Full criteria embedded directly in the reviewer subagent prompt
- Criteria sourced from `.claude/agents/*.md` files (plot-reviewer, code-reviewer,
  physics-reviewer, numerical-verifier, note-reviewer)
- Reviewer does not need to read external files for its checklist — it's in its prompt
- Reviewer IS instructed to Read source files to verify content

### 3. Iteration Control & Breaking

- `MAX_ITERATIONS = 5` (configurable per command)
- Iteration count tracked on disk in the log file (`## Iteration N` headers)
- Breaking logic:
  - `VERDICT: PASS` → exit approved
  - `VERDICT: FAIL` AND count < MAX → amend and re-review
  - `VERDICT: FAIL` AND count >= MAX → exit escalated
- Iteration count determined by re-reading the log file (not context memory)

### 4. Structured Tracking (JSONL)

- Start event emitted at setup (before any work)
- End event emitted at exit (both approved and escalated paths)
- All events appended to `.claude/logs/tracking.jsonl`
- Event schema:
  ```json
  {"event":"start","command":"<name>","task":"<task>","slug":"<slug>","timestamp":"<ISO 8601>"}
  {"event":"end","command":"<name>","status":"approved|escalated","iterations":<N>,"timestamp":"<ISO 8601>","criticals":<n>,"warnings":<n>,"slug":"<slug>"}
  ```
- Escalated events additionally include `"unresolved":["<short descriptions>"]`
- Analysis: `python3 .claude/scripts/analyze-tracking.py`

### 5. Log File Naming & Recovery

- Pattern: `review-<type>-YYYYMMDD-HHMMSS-<slug>.md`
- Slug: 3-5 words from task description, kebab-case
- Active session pointer: `.claude/logs/active-session.txt` (overwritten at start,
  deleted at exit)
- Compaction recovery procedure:
  1. Read `active-session.txt`
  2. Fallback: `ls -t .claude/logs/review-<type>-*.md | head -1`
  3. Read log file fully — it is ground truth

## Reviewer Output Format

The reviewer subagent must structure its response with these sections in order:

```
## Issues
[numbered list, or "None found."]

N. [CATEGORY] Description
   File: [path:line or element name]
   Severity: CRITICAL | WARNING | INFO
   Fix: [specific instruction]

## Numerical Verification
[structured per-number verification, or "No numbers to verify."]

## What Passed
[brief list]

VERDICT: PASS | FAIL
```

### Section Roles

| Section | Purpose | Drives verdict? |
|---------|---------|:---------------:|
| Issues | All problems: style, logic, conventions, content, correctness | Yes — CRITICAL/WARNING block PASS |
| Numerical Verification | Structured evidence for number accuracy specifically | Indirectly — MISMATCHes appear as CRITICAL in Issues |
| What Passed | Confirmation of correct aspects | No |
| VERDICT | Final binary decision | Controls loop flow |

### Severity Levels

| Level | Meaning | Blocks PASS? | Examples |
|-------|---------|:------------:|---------|
| CRITICAL | Wrong result, data corruption, logic error, physics mistake | Yes | Sign convention swapped, integer division in efficiency, missing cut |
| WARNING | Readability issue, style violation, potential confusion | Yes | Legend obscures data, axis label cropped, inconsistent notation |
| INFO | Observation, minor suggestion, note for awareness | No | Could use a comment here, consider reordering for clarity |

### Numerical Verification Applicability

Not all tasks produce numbers. The section adapts:

| Command | Typical numbers | When "No numbers to verify" |
|---------|----------------|----------------------------|
| review-plot | Integrals, bin counts, scale factors | Purely cosmetic plot fixes |
| review-analysis-code | Yields, efficiencies, event counts | Refactoring with no output change |
| review-investigation | All cited quantities | Never (investigations always cite numbers) |
| review-note | Numbers in prose (e.g., "85% efficiency") | Qualitative sections (motivation, detector description) |
| review-paper | Same as note + systematic uncertainties | Same as note |

When numbers ARE present, the reviewer:
1. Identifies the source (ROOT file, histogram, tree, computation)
2. Independently extracts the same quantity (via Read/Bash)
3. Reports MATCH (<0.1% relative difference) or MISMATCH
4. Any MISMATCH also appears as a CRITICAL item in ## Issues

## CLAUDE.md Integration

The commands are auto-dispatched via explicit rules in CLAUDE.md:

```markdown
## Auto-Dispatch Rules
When the user asks to create/fix/modify a plot → invoke /review-plot
When the user asks to write/modify/fix C++/ROOT/RDF code → invoke /review-analysis-code
...
```

This provides a two-layer trigger:
1. **Description matching** — Claude matches user request against command descriptions
2. **Explicit CLAUDE.md rules** — direct instruction to invoke specific commands

## Adding a New Command

To create a new `/review-*` command:

1. Create `.claude/commands/review-<type>.md`
2. Write multi-line `description:` with trigger phrases
3. Follow the template structure:
   - Setup (read refs, init log, emit tracking start)
   - Compaction recovery
   - Loop (Execute → Spawn reviewer → Parse → Decide → Amend)
   - Exit: Approved (update log, emit tracking end, delete active-session)
   - Exit: Escalate (same + unresolved list)
4. Compose reviewer prompt with:
   - Role statement
   - Task description
   - File list with "the reviewer will Read them"
   - Full review criteria inline
   - Exact response format specification
5. Add auto-dispatch rule to CLAUDE.md
6. Test with a real task and verify tracking output
