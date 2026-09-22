#!/usr/bin/env python3
"""Report Claude-to-Codex adapter drift without changing Claude-owned files."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Mapping:
    classification: str
    targets: tuple[str, ...]
    note: str


COMMAND_TARGETS = {
    "review-analysis-code.md": (".agents/skills/review-analysis-code/SKILL.md",),
    "review-plot.md": (".agents/skills/review-plot/SKILL.md",),
    "review-investigation.md": (".agents/skills/review-investigation/SKILL.md",),
    "steer-pipeline.md": (".agents/skills/steer-pipeline/SKILL.md",),
    "wrap-up.md": (".agents/skills/wrap-up/SKILL.md",),
    "review-note.md": (".agents/skills/academic-writing-gates/SKILL.md",),
    "review-paper.md": (".agents/skills/academic-writing-gates/SKILL.md",),
    "verify-citations.md": (".agents/skills/academic-writing-gates/SKILL.md",),
    "sync-note-figures.md": (".agents/skills/academic-writing-gates/SKILL.md",),
    "check-note-sync.md": (".agents/skills/academic-writing-gates/SKILL.md",),
    "compile-note.md": (".agents/skills/academic-writing-gates/SKILL.md",),
}

AGENT_TARGETS = {
    "code-reviewer.md": (".codex/agents/code-reviewer.toml",),
    "physics-reviewer.md": (".codex/agents/physics-reviewer.toml",),
    "numerical-verifier.md": (".codex/agents/numerical-verifier.toml",),
    "plot-reviewer.md": (".codex/agents/plot-reviewer.toml",),
    "note-reviewer.md": (".agents/skills/academic-writing-gates/SKILL.md",),
    "citation-verifier.md": (".agents/skills/academic-writing-gates/SKILL.md",),
}


def repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in (current, *current.parents):
        if (candidate / ".git").exists():
            return candidate
    raise SystemExit("Not inside a Git repository")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources(root: Path) -> list[Path]:
    claude = root / ".claude"
    found = [claude / "CLAUDE.md", claude / "settings.local"]
    found.extend(sorted((claude / "commands").glob("*.md")))
    found.extend(sorted((claude / "agents").glob("*.md")))
    found.extend(sorted((claude / "conventions").glob("*.md")))
    found.extend(sorted((claude / "skills").glob("*/SKILL.md")))
    return [path for path in found if path.is_file()]


def mapping_for(root: Path, path: Path) -> Mapping:
    rel = path.relative_to(root).as_posix()
    if rel == ".claude/CLAUDE.md":
        return Mapping("translated", ("AGENTS.md", "Analysis/AGENTS.md"), "Instructions need Codex semantics; active docs are read directly by hook.")
    if rel == ".claude/settings.local":
        return Mapping("translated", (".codex/hooks/pre_tool_policy.py", ".codex/rules/dimuon.rules", ".codex/profiles/"), "Claude permissions map to hooks, rules, and user profile templates.")
    if rel.startswith(".claude/commands/"):
        targets = COMMAND_TARGETS.get(path.name)
        if targets:
            return Mapping("translated", targets, "Claude command is a read-only workflow reference for this Codex skill.")
        return Mapping("unsupported", (), "No dedicated Codex adapter yet; keep visible for manual triage.")
    if rel.startswith(".claude/agents/"):
        targets = AGENT_TARGETS.get(path.name)
        if targets:
            return Mapping("translated", targets, "Codex adapter uses TOML or a consolidated skill.")
        return Mapping("unsupported", (), "No dedicated Codex custom-agent adapter yet.")
    if rel.startswith(".claude/conventions/"):
        return Mapping("direct-reference", ("AGENTS.md",), "Codex reads this convention directly when routed by a skill.")
    if rel.startswith(".claude/skills/"):
        return Mapping("unsupported", (), "Claude skill may invoke Claude-specific CLI or write .claude; review manually before symlinking.")
    return Mapping("unsupported", (), "Unclassified source.")


def load_baseline(path: Path) -> dict[str, str]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    hashes = value.get("source_hashes", {}) if isinstance(value, dict) else {}
    return hashes if isinstance(hashes, dict) else {}


def write_baseline(path: Path, hashes: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"schema_version": 1, "source_hashes": dict(sorted(hashes.items()))}
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Repository root or a path inside it")
    parser.add_argument("--accept", action="store_true", help="Record current Claude source hashes after reviewed translation")
    parser.add_argument("--check", action="store_true", help="Exit nonzero when source drift or missing mapped targets exist")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    root = repo_root(args.root)
    baseline_path = root / ".agents/skills/claude-to-codex-setting-sync/references/sync-baseline.json"
    baseline = load_baseline(baseline_path)
    current: dict[str, str] = {}
    rows = []
    problems = 0

    for path in sources(root):
        rel = path.relative_to(root).as_posix()
        sha = digest(path)
        current[rel] = sha
        old = baseline.get(rel)
        status = "unchanged" if old == sha else ("new" if old is None else "changed")
        mapping = mapping_for(root, path)
        missing_targets = [target for target in mapping.targets if not (root / target).exists()]
        if status != "unchanged" or missing_targets:
            problems += 1
        rows.append(
            {
                "source": rel,
                "status": status,
                "classification": mapping.classification,
                "targets": list(mapping.targets),
                "missing_targets": missing_targets,
                "note": mapping.note,
            }
        )

    removed = sorted(set(baseline) - set(current))
    problems += len(removed)

    if args.accept:
        missing_required = [row for row in rows if row["classification"] == "translated" and row["missing_targets"]]
        if missing_required:
            print("Refusing --accept because translated targets are missing:", file=sys.stderr)
            for row in missing_required:
                print(f"- {row['source']}: {', '.join(row['missing_targets'])}", file=sys.stderr)
            return 2
        write_baseline(baseline_path, current)
        baseline = current
        for row in rows:
            row["status"] = "unchanged"
        removed = []
        problems = 0

    report = {
        "root": str(root),
        "baseline": str(baseline_path.relative_to(root)),
        "rows": rows,
        "removed_sources": removed,
        "drift_count": problems,
    }

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("# Claude → Codex setting sync report")
        print()
        print(f"Baseline: `{report['baseline']}`")
        print(f"Drift or missing-target count: **{problems}**")
        print()
        print("| Source | Status | Class | Targets |")
        print("|---|---|---|---|")
        for row in rows:
            targets = ", ".join(f"`{target}`" for target in row["targets"]) or "manual"
            if row["missing_targets"]:
                targets += " (missing: " + ", ".join(row["missing_targets"]) + ")"
            print(f"| `{row['source']}` | {row['status']} | {row['classification']} | {targets} |")
        if removed:
            print()
            print("Removed Claude sources recorded in baseline:")
            for rel in removed:
                print(f"- `{rel}`")

    return 1 if args.check and problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
