#!/usr/bin/env python3
"""Codex lifecycle context and checkpoint guard for the dimuon repository."""

from __future__ import annotations

import datetime as dt
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


def read_event() -> dict[str, Any]:
    try:
        value = json.load(sys.stdin)
    except (json.JSONDecodeError, OSError):
        return {}
    return value if isinstance(value, dict) else {}


def git_root(cwd: str | None) -> Path | None:
    start = Path(cwd or ".").resolve()
    try:
        result = subprocess.run(
            ["git", "-C", str(start), "rev-parse", "--show-toplevel"],
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    root = Path(result.stdout.strip())
    return root if root.is_dir() else None


def git_sha(root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() or None


def active_claude_tracking_docs(root: Path) -> list[str]:
    claude = root / ".claude" / "CLAUDE.md"
    try:
        text = claude.read_text(encoding="utf-8")
    except OSError:
        return []
    marker = "## Active Tracking Docs"
    if marker not in text:
        return []
    section = text.split(marker, 1)[1]
    section = re.sub(r"<!--.*?-->", "", section, flags=re.DOTALL)
    paths: list[str] = []
    for line in section.splitlines():
        if line.startswith("## "):
            break
        if not line.lstrip().startswith("-"):
            continue
        match = re.search(r"`([^`]+\.md)`", line)
        if not match:
            continue
        raw = match.group(1)
        if "/Analysis/" in raw:
            rel = "Analysis/" + raw.split("/Analysis/", 1)[1]
        else:
            candidate = Path(raw)
            rel = str(candidate if not candidate.is_absolute() else candidate.name)
        if rel not in paths:
            paths.append(rel)
    return paths


def latest_stage(text: str) -> str:
    match = re.search(r"^## Latest stage\s*$([\s\S]*?)(?=^## |\Z)", text, re.MULTILINE | re.IGNORECASE)
    return match.group(1).strip() if match else ""


def active_codex_tracking_docs(root: Path) -> list[str]:
    tracking_dir = root / ".codex" / "tracking"
    active: list[str] = []
    for path in sorted(tracking_dir.glob("*.md")):
        try:
            stage = latest_stage(path.read_text(encoding="utf-8"))
        except OSError:
            continue
        if stage and not re.search(r"\b(cleared|complete|closed)\b", stage, re.IGNORECASE):
            active.append(str(path.relative_to(root)))
    return active


def checkpoint_path(root: Path, session_id: str) -> Path:
    safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", session_id or "unknown")
    return root / ".codex" / "state" / "checkpoints" / f"{safe_id}.json"


def write_checkpoint(root: Path, event: dict[str, Any]) -> Path:
    session_id = str(event.get("session_id") or "unknown")
    path = checkpoint_path(root, session_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "session_id": session_id,
        "turn_id": event.get("turn_id"),
        "trigger": event.get("trigger"),
        "cwd": event.get("cwd"),
        "git_sha": git_sha(root),
        "transcript_path": event.get("transcript_path"),
        "active_claude_tracking_docs": active_claude_tracking_docs(root),
        "active_codex_tracking_docs": active_codex_tracking_docs(root),
        "recovery_protocol": [
            "Read AGENTS.md and the closest nested AGENTS.md.",
            "Read every active tracking document fully.",
            "Resume from the tracking document's Latest stage, not memory.",
            "Verify the next analysis step against its Physics Procedure.",
        ],
    }
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path


def load_checkpoint(root: Path, session_id: str) -> dict[str, Any] | None:
    path = checkpoint_path(root, session_id)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def startup_context(root: Path, event: dict[str, Any]) -> str:
    claude_docs = active_claude_tracking_docs(root)
    codex_docs = active_codex_tracking_docs(root)
    lines = [
        "DIMUON PROJECT CONTEXT GUARD:",
        "Before domain work, read AGENTS.md and every active tracking document below fully.",
        "For analysis implementations, re-read the authoritative Physics Procedure before planning or editing.",
        "Never treat Codex memory as evidence for scientific facts or numbers.",
    ]
    if claude_docs:
        lines.append("Canonical active analysis tracking documents:")
        lines.extend(f"- {path}" for path in claude_docs)
    else:
        lines.append("No active analysis documents could be parsed; inspect .claude/CLAUDE.md manually.")
    if codex_docs:
        lines.append("Active Codex infrastructure tracking documents:")
        lines.extend(f"- {path}" for path in codex_docs)
    checkpoint = load_checkpoint(root, str(event.get("session_id") or "unknown"))
    if checkpoint:
        lines.append(
            "A pre-compaction checkpoint exists for this session at "
            f".codex/state/checkpoints/{checkpoint_path(root, str(event.get('session_id') or 'unknown')).name}; "
            "use it only as a routing aid and trust tracking documents over it."
        )
    return "\n".join(lines)


def has_incomplete_codex_tracking(root: Path) -> list[str]:
    incomplete: list[str] = []
    for rel in active_codex_tracking_docs(root):
        path = root / rel
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        if re.search(r"^- \[ \]", text, re.MULTILINE):
            incomplete.append(rel)
    return incomplete


def emit(value: dict[str, Any]) -> None:
    print(json.dumps(value, ensure_ascii=False))


def main() -> int:
    event = read_event()
    name = str(event.get("hook_event_name") or "")
    root = git_root(event.get("cwd"))
    if root is None:
        emit({"continue": True})
        return 0

    if name in {"SessionStart", "SubagentStart"}:
        context = startup_context(root, event)
        emit(
            {
                "continue": True,
                "hookSpecificOutput": {
                    "hookEventName": name,
                    "additionalContext": context,
                },
            }
        )
        return 0

    if name == "PreCompact":
        write_checkpoint(root, event)
        emit({"continue": True})
        return 0

    if name == "PostCompact":
        emit(
            {
                "continue": True,
                "systemMessage": (
                    "Dimuon context was compacted. Re-read active tracking documents and "
                    "resume from Latest stage before further analysis work."
                ),
            }
        )
        return 0

    if name == "Stop":
        if event.get("stop_hook_active"):
            emit({"continue": True})
            return 0
        message = str(event.get("last_assistant_message") or "")
        claims_full_completion = bool(
            re.search(
                r"\b(task|implementation|migration|work) (?:is |was )?(?:fully )?(?:complete|completed|done)\b|\ball done\b",
                message,
                re.IGNORECASE,
            )
        )
        incomplete = has_incomplete_codex_tracking(root) if claims_full_completion else []
        if incomplete:
            emit(
                {
                    "decision": "block",
                    "reason": (
                        "Before claiming full completion, update the active Codex tracking "
                        "checklist and clear Latest stage: " + ", ".join(incomplete)
                    ),
                }
            )
        else:
            emit({"continue": True})
        return 0

    emit({"continue": True})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
