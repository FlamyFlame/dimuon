---
name: hep-skill-creator
description: >
  Create, test, and optimize Claude Code skills tailored to heavy-ion physics analysis
  workflows. Use this skill whenever the user wants to create a new analysis skill (e.g.
  publication-quality plotting, template fitting, RooFit workspace setup, event selection
  cuts, FCal/ZDC diagnostics, RDataFrame pipelines, Athena job submission, grid job
  management, or any ATLAS/HEP-specific workflow), or wants to improve or optimize an
  existing analysis skill's triggering description. Also trigger when the user says
  "make a skill for", "turn this into a skill", "create a skill that", or refers to
  automating a recurring analysis pattern. Even if the user just says "I keep doing X
  manually, can we automate it?" — that's a skill request. This skill wraps the upstream
  anthropics/skills skill-creator with physics-aware context gathering and conventions
  from the user's repository, CLAUDE.md, and analysis environment.
---

# HEP Analysis Skill Creator

A meta-skill that creates Claude Code skills grounded in your specific heavy-ion physics
analysis environment. It bridges two things: your domain knowledge (repo conventions,
ROOT/RDataFrame patterns, ATLAS style, data paths, cut definitions) and the upstream
skill-creator's structured eval/optimization pipeline.

## Prerequisites

Before first use, ensure the upstream skill-creator is available locally:

```bash
# Check if already cloned
ls ~/tools/skill-creator/SKILL.md 2>/dev/null || {
  echo "Cloning skill-creator from anthropics/skills..."
  mkdir -p ~/tools
  cd ~/tools
  git clone --depth 1 --filter=blob:none --sparse \
    https://github.com/anthropics/skills.git skills-upstream
  cd skills-upstream
  git sparse-checkout set skills/skill-creator
  ln -sf ~/tools/skills-upstream/skills/skill-creator ~/tools/skill-creator
}
```

The scripts require Python with `pyyaml` installed. The description optimization
scripts use `claude -p` (Claude Code CLI) as a subprocess — no separate API key needed.

## Workflow Overview

The process has four phases. Always confirm with the user before moving to the next phase.

```
Phase 1: Context Gathering    → Read repo, CLAUDE.md, memories for domain grounding
Phase 2: Skill Drafting       → Write SKILL.md + bundled references using physics context
Phase 3: Test & Iterate       → Run test cases, review outputs, improve the skill
Phase 4: Optimize & Deliver   → Run description optimizer, git commit, summarize
```

---

## Phase 1: Context Gathering

Before writing anything, gather the context that makes this skill physics-aware rather
than generic. Read these sources in order — stop and discuss with the user if anything
is missing or ambiguous.

### 1.1 Read the project CLAUDE.md

```bash
cat CLAUDE.md
```

Extract: coding conventions, preferred patterns (first-principles, minimal code),
safety rules (never automate passwords, ask before irreversible actions), any
analysis-specific instructions.

### 1.2 Read the repository structure

```bash
find . -maxdepth 3 -type f \( -name "*.C" -o -name "*.cxx" -o -name "*.py" \
  -o -name "*.h" -o -name "*.sh" -o -name "Makefile" -o -name "CMakeLists.txt" \) \
  | head -80
```

Identify: existing macros, plotting code, fitting code, event selection scripts,
job submission scripts. These are the patterns the new skill should encode.

### 1.3 Read relevant existing code

Based on what the user wants the skill to do, read the actual source files that
implement that workflow today. For example:
- Plotting skill → read the plotting macros, style headers, axis label conventions
- Template fitting → read existing RooFit workspace code, PDF definitions
- Event selection → read cut implementation code, diagnostic plot generators

Extract concrete patterns: function signatures, histogram naming conventions,
file path conventions, canvas layout patterns, style settings.

### 1.4 Check for existing skills

```bash
ls -la .claude/skills/ 2>/dev/null
```

Avoid duplicating existing skills. If the new skill overlaps, discuss with the
user whether to extend the existing one or create a separate skill.

### 1.5 Summarize context to user

Present a brief summary of what you found and confirm:
- "Here's what I understand about your [plotting/fitting/etc.] conventions: ..."
- "The new skill should encode these patterns: ..."
- "Does this look right? Anything I'm missing?"

---

## Phase 2: Skill Drafting

### 2.1 Choose a skill name and location

Skills live in `.claude/skills/<skill-name>/`. Use descriptive kebab-case names:
- `publication-plots` — ATLAS-style publication-quality plots
- `template-fitting` — RooFit template fits with systematic variations
- `event-selection` — Event cleaning cuts and diagnostic grids
- `rdataframe-pipeline` — RDataFrame analysis chain setup
- `grid-jobs` — PANDA/rucio grid job submission and monitoring

### 2.2 Write the SKILL.md

Follow the upstream skill-creator anatomy:

```
<skill-name>/
├── SKILL.md              (required — under 500 lines)
├── references/           (analysis-specific reference docs)
│   ├── style-guide.md    (e.g., ATLAS style settings, axis labels)
│   ├── conventions.md    (e.g., histogram naming, file paths)
│   └── examples.md       (e.g., annotated code examples from the repo)
└── scripts/              (deterministic helper scripts if needed)
```

When writing the SKILL.md body:

**Ground instructions in the actual codebase.** Don't write generic "use ROOT to
make a histogram" instructions. Write "use the TStyle settings from
`Analysis/plotStyle.h`, apply the ATLAS label with `ATLASLabel()`, use the canvas
division pattern from the 3×2 diagnostic grids." Reference actual file paths in the
repo.

**Explain the physics why.** Instead of "always use log scale on the y-axis for
FCal distributions", write "FCal ΣE_T distributions span several orders of magnitude
in Pb+Pb due to the wide range of centralities, so log scale on y is needed to see
both peripheral and central events."

**Encode the user's CLAUDE.md preferences.** The skill should inherit: first-principles
thinking, precise minimal code, stop to discuss unclear goals, ask before long-running
or irreversible actions.

**Make the description pushy.** The description field drives triggering. Include
specific physics terms, ROOT class names, and casual phrasings a physicist would use.
See the examples in the description optimization section below.

### 2.3 Write bundled reference files

For each reference file, extract and distill the relevant patterns from the repo.
Keep individual reference files under 300 lines. If a reference would be longer,
split by topic (e.g., `style-guide.md` vs `canvas-layouts.md`).

Include a table of contents at the top of any reference file over 100 lines.

### 2.4 Review draft with user

Present the full SKILL.md and reference files. Ask:
- "Does this capture how you actually do [X] in practice?"
- "Are there edge cases or variations I missed?"
- "Any conventions that should be stricter or more flexible?"

---

## Phase 3: Test & Iterate

Follow the upstream skill-creator's test loop. Read the full instructions at:
```bash
cat ~/tools/skill-creator/SKILL.md
```

The key steps adapted for this environment:

### 3.1 Write test cases

Create 2-3 realistic prompts — the kind of thing you'd actually type at the terminal.
Save to `<skill-name>/evals/evals.json`.

For a publication-plots skill, examples might be:
```json
{
  "skill_name": "publication-plots",
  "evals": [
    {
      "id": 1,
      "prompt": "make a plot of the dimuon invariant mass spectrum for Pb+Pb 2024 data, ATLAS style, log y, with the J/psi and Upsilon peaks labeled",
      "expected_output": "A ROOT macro or script that produces an ATLAS-style plot with correct axis labels, log scale, peak annotations, and uses the project's style conventions",
      "files": []
    },
    {
      "id": 2,
      "prompt": "I need the 3x2 diagnostic grid for the ZDC pre-sample cut showing before and after distributions",
      "expected_output": "A 6-panel canvas following the Mason-style diagnostic grid pattern with before/after overlays for the ZDC pre-sample amplitude cut",
      "files": []
    }
  ]
}
```

Confirm test cases with the user before running.

### 3.2 Run test cases

Since Claude Code on BNL doesn't have subagents, run test cases sequentially:

```bash
# For each test case, invoke Claude with the skill:
claude -p "<test prompt>" --skill-path .claude/skills/<skill-name>
```

Save outputs to `<skill-name>-workspace/iteration-1/eval-<ID>/`.

### 3.3 Review with user

Present each output. For physics skills, review criteria include:
- Does the code compile and run with the project's ROOT version?
- Are ATLAS style conventions followed correctly?
- Are file paths and histogram names consistent with repo conventions?
- Is the physics correct (axis labels, units, scale choices)?

### 3.4 Iterate

Based on feedback, improve the SKILL.md. Common issues for physics skills:
- Instructions too generic (not grounded in the actual repo patterns)
- Missing edge cases (e.g., pp vs Pb+Pb require different style choices)
- Over-specified (too rigid for variations the user might want)

Rerun test cases after each revision.

---

## Phase 4: Optimize & Deliver

### 4.1 Description optimization

Generate trigger eval queries — 8-10 should-trigger, 8-10 should-not-trigger.

For physics skills, **should-trigger** examples should include:
- Casual physics shorthand: "make me an FCal plot", "set up the template fit"
- Specific technical requests: "create a RooFit workspace with crystal ball + exponential"
- Workflow requests: "I need to redo the ZDC banana cut plots for the new dataset"

**Should-not-trigger** examples (the tricky near-misses):
- Generic ROOT questions: "how do I create a TH1F?" (too simple, no skill needed)
- ML work on the Mac: "train the PyTorch model on the Olist dataset"
- Different analysis domain: "make a plot of the neural network loss curve"
- Athena/grid questions that overlap but need a different skill

Save as `<skill-name>/evals/trigger-eval.json`, confirm with user, then run:

```bash
cd ~/tools/skill-creator
python -m scripts.run_loop \
  --eval-set <path-to-trigger-eval.json> \
  --skill-path <path-to-skill> \
  --model claude-sonnet-4-6 \
  --max-iterations 5 \
  --verbose
```

Apply the `best_description` from the output to the SKILL.md frontmatter.

### 4.2 Git commit and push

```bash
cd <repo-root>
git add .claude/skills/<skill-name>/
git commit -m "feat: add <skill-name> skill for [brief description]"
```

Ask the user before pushing: "Ready to push to remote?"

### 4.3 Summarize

Present a summary:
- Skill name and location
- What it does (1-2 sentences)
- Trigger phrases (from the optimized description)
- Test results (what worked, what was tricky)
- Bundled references and scripts
- Suggestions for future improvements or companion skills

---

## Important Notes

**Cost awareness:** The description optimization loop in `run_loop.py` invokes
`claude -p` many times (3 runs per query × ~20 queries × up to 5 iterations). On a
Pro subscription this uses your included usage, but be aware it can take 15-30 minutes
and significant context. Warn the user before starting Phase 4.1.

**BNL environment specifics:**
- ROOT and RDataFrame are available via ATLAS setup scripts, not pip
- Data files live on `/usatlas/u/yuhanguo/usatlasdata/` — skills should reference
  paths via conventions, not hardcode absolute paths
- Grid certificates require `voms-proxy-init` which must NEVER be automated
  (consistent with CLAUDE.md rules)

**Skill count:** If you accumulate many skills (40+), triggering reliability can
degrade. Keep the skill set focused on genuinely multi-step workflows. Simple
one-liner tasks don't need skills.
