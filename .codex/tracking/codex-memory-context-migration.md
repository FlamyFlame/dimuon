# Codex Memory and Long-Task Context Migration

## Objective

Add persistent, project-aware Codex context management for the dimuon physics
analysis while preserving Claude Code's existing workflow and files. Enable
native Codex memories, share the repository's canonical scientific context by
reference, externalize long-task state, and add Codex-native skills and guards
for publication-level work.

## Safety invariants

- Do not modify anything under `.claude/`.
- Do not modify any pre-existing file under `Analysis/docs/`.
- Do not modify the dirty `IntNotes` submodule or unrelated user changes.
- Treat `Analysis/docs/analysis_overview.md`, relevant Physics Procedures, and
  canonical tracking documents as scientific ground truth; memories are never
  authoritative for numerical or current factual claims.
- Keep generated Codex state and logs unversioned.
- Make synchronization from Claude to Codex read-only on the Claude side and
  dry-run by default on the Codex side.

## Design decisions

1. The Git repository root (`dimuon/`) is the Codex project root. Analysis work
   may start in `Analysis/`; local slide work reads plots from `../plots/` at
   repository root or `../../plots/` from `Analysis/`.
2. Project behavior is versioned in this repository. User-local memories and
   machine-specific profiles remain under `~/.codex`.
3. Existing Claude knowledge, conventions, and tracking documents remain
   canonical and are read directly. They are not copied or symlinked.
4. Claude commands and agents require Codex-native adapters because their
   formats and execution semantics differ.
5. Hooks provide deterministic start/compaction/stop context guards; skills
   provide reusable task workflows; `AGENTS.md` provides always-on policy.

## Implementation plan and checklist

### 1. Native memory

- [x] Back up `~/.codex/config.toml`.
- [x] Enable `[features].memories` without disturbing existing settings.
- [x] Configure conservative memory use/generation settings.
- [x] Verify configuration parses and `codex features list` reports memories.

### 2. Versioned project instruction spine

- [x] Add root `AGENTS.md` with context hierarchy, safety, tracking, review,
      path, and completion rules.
- [x] Add `Analysis/AGENTS.md` with physics-procedure and publication-level
      analysis requirements.
- [x] Add `.codex/config.toml` with project-safe settings only.
- [x] Narrow `.gitignore` so declarative `.codex` files are tracked while
      generated state/logs stay ignored.

### 3. Long-task context hooks

- [x] Add a deterministic context helper that parses active tracking documents
      from `.claude/CLAUDE.md` without modifying it.
- [x] Add checkpoint validation/loading helpers for startup and compaction.
- [x] Configure `SessionStart`, `PreCompact`, `PostCompact`, and `Stop` hooks.
- [x] Ensure hooks degrade safely when cluster-only paths are unavailable.

### 4. Codex skills and agents

- [x] Create and validate `dimuon-context`.
- [x] Create and validate `track-long-task`.
- [x] Create and validate the priority workflow adapters:
      `review-analysis-code`, `review-plot`, `review-investigation`,
      `steer-pipeline`, `wrap-up`, and `academic-writing-gates`.
- [x] Create and validate `claude-to-codex-setting-sync` with a deterministic
      drift-report script and dry-run default.
- [x] Add Codex custom-agent adapters for code, physics, numerical, and plot
      review, while preserving explicit-user-request requirements for
      subagent use.

### 5. Permission and environment compatibility

- [x] Add project command rules for high-value safe/blocked command patterns.
- [x] Add a pre-tool policy hook for protected data and prohibited installs.
- [x] Add local and BNL profile templates without embedding machine-local
      settings in active project config.

### 6. Verification

- [x] Validate TOML, JSON, Python, skills, and shell/rule files.
- [x] Run sync drift reporting and hook smoke tests.
- [x] Verify instruction discovery from repository root and `Analysis/`.
- [x] Verify `.claude/`, existing `Analysis/docs/`, and `IntNotes` were not
      changed by this migration.
- [x] Review the final diff for accidental secrets or absolute local paths in
      versioned configuration.

## Progress log

- 2026-06-23: Plan written before implementation. Read-only repository and
  Codex capability audit completed in the preceding planning turn. Current
  Codex CLI is 0.139.0; memories are installed but disabled. Existing Claude
  context, commands, agents, KB, logs, and active tracking structure were
  inventoried without modification.
- 2026-06-23: Backed up `~/.codex/config.toml` to
  `~/.codex/config.toml.pre-dimuon-memory-20260623`; enabled memory use and
  generation with external-context generation disabled. TOML parses and
  `codex features list` reports `memories ... true`. Marked the Git root trusted.
- 2026-06-23: Added the root and Analysis instruction chain, project config,
  and narrow `.gitignore` exceptions. Existing Claude and analysis documents
  remain untouched.
- 2026-06-23: Added lifecycle and pre-tool hooks. Smoke tests parsed the four
  active analysis tracking documents, wrote an ignored checkpoint, restored
  routing context, blocked false completion with an open checklist, blocked
  package installation, and allowed read-only access to `.claude`.
- 2026-06-23: Created nine Codex skills and four read-only custom reviewer
  agents. All skills pass `quick_validate.py`. Added tested command rules,
  local/BNL profile templates, and an accepted SHA-256 sync baseline. The sync
  script reports zero current drift while explicitly listing unsupported
  Claude-only workflows for later manual translation.
- 2026-06-23: A fresh read-only Codex 0.142.0 process launched from `Analysis/`
  loaded both root and nested `AGENTS.md`, discovered all nine repository skills,
  and received the exact `DIMUON PROJECT CONTEXT GUARD` developer context with
  the four canonical active tracking paths. Updated the standalone CLI from
  0.139.0 to 0.142.0 to match the desktop app and eliminate hook behavior drift.
- 2026-06-23: Final acceptance passed: user/project TOML, hook JSON/Python,
  profile templates, all skill validators, rule unit checks, sync `--check`,
  `git diff --check`, generated-state ignores, secret-pattern scan, and local
  absolute-path scan. `git status -- .claude Analysis/docs` is empty. The
  pre-existing dirty `IntNotes` submodule remains untouched.

## Results and compatibility limits

- Native memories are enabled globally with use/generation on and memory
  generation disabled for threads using external context. The user will restart
  Codex after this task.
- Project instructions, skills, custom agents, hooks, rules, profiles, and the
  sync baseline are versionable in the dimuon repository; generated state and
  logs are ignored.
- Active Claude tracking context and conventions synchronize by direct read-only
  reference. No symlink is useful for these files because Codex needs a semantic
  adapter rather than another path to Claude-specific syntax.
- Dedicated Codex adapters now cover the highest-value analysis, plotting,
  investigation, steering, wrap-up, and academic-writing workflows.
- `kb-build`, `kb-review`, `review-pipeline`, and `review-plugin`, plus Claude's
  executor/KB/plugin agents and `hep-skill-creator`, remain explicitly flagged
  by the sync report as unsupported dedicated adapters. Their sources remain
  readable, but they are not falsely presented as native Codex workflows. The
  HEP skill creator in particular invokes `claude -p` and writes `.claude`, so a
  symlink would interfere with Claude ownership.
- Hooks require one-time review/trust through `/hooks` after restart. Changed
  hook definitions will require re-trust by design.

## Completion summary

The Codex memory and long-context system is implemented and validated without
modifying Claude-owned configuration, existing analysis tracking/docs, or the
dirty IntNotes submodule. The only remaining action is the user's planned Codex
restart followed by one-time hook trust review.

## Latest stage

cleared — implementation and validation complete; restart and `/hooks` trust are
user-side activation steps.
