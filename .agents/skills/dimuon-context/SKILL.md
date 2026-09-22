---
name: dimuon-context
description: Load and route authoritative context for the Run 3 ATLAS dimuon analysis. Use before analysis code, plots, pipelines, numerical claims, datasets, selections, efficiencies, R_AA, internal-note text, paper text, or any task where project-specific physics knowledge and current tracking state matter.
---

# Load dimuon context

1. Resolve the Git root and read the root `AGENTS.md` plus the closest nested
   `AGENTS.md` for the working directory.
2. Read `.claude/CLAUDE.md` sections `Tracking Documents`, `Documentation
   References`, `Auto-Dispatch Rules`, and `Active Tracking Docs`. Treat Claude
   slash commands as workflow references, not Codex commands.
3. Read every active tracking document fully. For implementation work, re-read
   `Objective`, `Physics Procedure`, `Remaining Work`, and `Latest Stage` first.
4. Route additional context narrowly:
   - methodology: `Analysis/docs/analysis_overview.md`
   - architecture: `Analysis/README.md`
   - pipeline behavior: the matching file in `Analysis/docs/`
   - current status: roadmap and status tracking documents
   - domain knowledge: `.claude/kb/index.md`, then only relevant linked entries
   - confirmed note inputs: `IntNotes/analysis_metadata.md`
5. Search historical tracking documents before inferring a previously settled
   fact. Use `rg -n` with physics terms, variable names, sample tags, and paths.

## Evidence discipline

- Prefer the authoritative Physics Procedure over code patterns or memory.
- Re-derive quoted numbers from source files or logs when practical.
- Label placeholders, guesses, preliminary results, missing inputs, and stale
  outputs explicitly.
- Never use Codex memory as the sole support for a scientific claim.
- Do not modify `.claude/`; share its context by direct read-only reference.
