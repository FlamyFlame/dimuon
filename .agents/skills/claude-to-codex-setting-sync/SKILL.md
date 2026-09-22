---
name: claude-to-codex-setting-sync
description: Detect and reconcile drift between the dimuon repository's Claude Code instructions, commands, agents, conventions, settings, and Codex-native adapters. Use when Claude configuration changes, Codex behavior seems stale, agent workflows are reviewed, a migration audit is requested, or before finalizing updates to AGENTS.md, .codex, or .agents/skills.
---

# Synchronize Claude behavior into Codex adapters

Preserve `.claude/` as read-only. Synchronization is directional:
Claude source behavior → reviewed Codex adapters.

1. Run the deterministic dry-run report:

   ```bash
   python3 .agents/skills/claude-to-codex-setting-sync/scripts/sync_settings.py
   ```

2. For every changed source, read the source and its mapped Codex targets. Classify
   the change as:
   - direct reference: no copy needed
   - mechanical metadata/path update
   - semantic translation required
   - unsupported/manual design needed
3. Translate semantics, not filenames. Replace Claude slash commands,
   `$ARGUMENTS`, Task/subagent assumptions, tool names, absolute BNL paths, and
   `.claude/logs` writes with Codex-native equivalents.
4. Do not auto-apply behavior changes from a hook. Present material changes for
   review, then edit only Codex-owned targets.
5. Validate affected skills, TOML, JSON, hooks, and rule tests.
6. After reviewing all drift, update the accepted baseline:

   ```bash
   python3 .agents/skills/claude-to-codex-setting-sync/scripts/sync_settings.py --accept
   ```

7. Run `--check`; it must exit successfully before declaring settings synchronized.

The script never edits `.claude/` and never rewrites Codex adapters. `--accept`
only records source hashes after a human-reviewed translation.
