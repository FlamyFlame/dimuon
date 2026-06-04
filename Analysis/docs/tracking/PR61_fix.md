# PR #61 Fix Requirements

Changes needed to pass CodeRabbit review and `check-skills` CI.

## All 4 SKILL.md files

### 1. Frontmatter — replace entirely

Remove all non-compliant keys. Only `name` and `description` are allowed.
`description` must start with "Use when" and describe triggering conditions only
(not a summary of the skill's content or workflow). Max ~1024 characters total.

**Required format:**
```yaml
---
name: <skill-name>          # letters, numbers, hyphens only
description: >-
  Use when <triggering condition 1>, <condition 2>, ...
---
```

**Keys to remove from all files:**
- `disable-model-invocation`
- `arguments`
- `argument-hint`
- `allowed-tools`

**Per-file values:**

| File | `name` | `description` |
|------|--------|---------------|
| `skills/migrate/SKILL.md` | `bnl-localgroupdisk-migrate` | `Use when you need to migrate a directory of local ROOT files to BNL-OSG2_LOCALGROUPDISK via Rucio on BNL SDCC nodes.` |
| `skills/preflight/SKILL.md` | `bnl-localgroupdisk-preflight` | `Use when running pre-flight checks before a LOCALGROUPDISK migration to verify Rucio account, grid proxy, quotas, RSE names, and pnfs mount availability on BNL SDCC.` |
| `skills/check-rule/SKILL.md` | `bnl-localgroupdisk-check-rule` | `Use when monitoring a Rucio replication rule until it completes (state OK) or gets stuck after adding a rule to BNL-OSG2_LOCALGROUPDISK.` |
| `skills/build-symlinks/SKILL.md` | `bnl-localgroupdisk-build-symlinks` | `Use when you need to build a symlink farm from an already-replicated Rucio dataset on BNL-OSG2_LOCALGROUPDISK for transparent proxy-free access on BNL SDCC nodes.` |

---

### 2. Section order — reorder body content

Canonical order (enforced by linter):

Overview
When to Use
Key Concepts
Canonical Patterns
Gotchas
Interop
Docs

Deep skills (`migrate`) may also include `## Worked Example` and
`## Troubleshooting` between Canonical Patterns and Gotchas.
Target length: ~150–250 lines for medium skills, ~300–400 for deep skills.

**Per-file mapping** (where to move existing content):

#### `migrate/SKILL.md`
| Existing content | → Section |
|-----------------|-----------|
| What the skill does (intro paragraph) | `## Overview` |
| Decision Point descriptions | `## When to Use` |
| Phase/Decision terminology | `## Key Concepts` |
| Phase 1–3 step sequences + bash snippets | `## Canonical Patterns` |
| Worked example (if included) | `## Worked Example` |
| Troubleshooting table | `## Troubleshooting` |
| `_orig` guard, DID conflict, FTS queue delay warnings | `## Gotchas` |
| Proxy/SDCC/symlink interop notes | `## Interop` |
| Rollback + Key facts | fold into `## Gotchas` / `## Interop` |

#### `preflight/SKILL.md`
| Existing content | → Section |
|-----------------|-----------|
| Intro paragraph | `## Overview` |
| When to run this vs migrate | `## When to Use` |
| Rucio scope, VOMS group explanation | `## Key Concepts` |
| The 5 bash check blocks + summary table | `## Canonical Patterns` |
| Common failure modes (NO LGD QUOTA, missing VOMS) | `## Gotchas` |
| Notes on pnfs mount path variant | `## Interop` |

#### `check-rule/SKILL.md`
| Existing content | → Section |
|-----------------|-----------|
| Intro | `## Overview` |
| When to run standalone vs inside migrate | `## When to Use` |
| State/Locks field meanings | `## Key Concepts` |
| Quick check bash block | `## Canonical Patterns` |
| Long-running monitor loop | `## Canonical Patterns` (subsection) |
| STUCK / Error interpretation | `## Gotchas` |
| FTS queue timing notes | `## Interop` |

#### `build-symlinks/SKILL.md`
| Existing content | → Section |
|-----------------|-----------|
| "How it works" paragraph | `## Overview` |
| When to use standalone vs inside migrate | `## When to Use` |
| LOCALGROUPDISK hash-path explanation | `## Key Concepts` |
| Step 1–3 bash blocks | `## Canonical Patterns` |
| Non-empty `$farm_dir` warning | `## Gotchas` |
| No proxy needed on SDCC note | `## Interop` |

---

### 3. Add `## Docs` section (all 4 files)

Every SKILL.md must end with a `## Docs` section as its **last** section,
linking canonical upstream documentation.

```markdown
## Docs

- [bnl-localgroupdisk plugin](https://github.com/FlamyFlame/claude-bnl-localgroupdisk)
- [BNL SDCC storage documentation](https://usatlas.github.io/af-docs/bnl/storage/)
```

---

## Summary checklist

- [ ] `migrate/SKILL.md` — fix frontmatter, reorder sections, add `## Docs`
- [ ] `preflight/SKILL.md` — fix frontmatter, reorder sections, add `## Docs`
- [ ] `check-rule/SKILL.md` — fix frontmatter, reorder sections, add `## Docs`
- [ ] `build-symlinks/SKILL.md` — fix frontmatter, reorder sections, add `## Docs`
- [ ] Push to `add-bnl-localgroupdisk-plugin` branch — CodeRabbit will re-review automatically
- [ ] Confirm `check-skills` CI passes (pixi linter)
