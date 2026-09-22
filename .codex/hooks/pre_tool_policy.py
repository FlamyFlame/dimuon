#!/usr/bin/env python3
"""Narrow Codex pre-tool guard for Claude-owned config and protected data."""

from __future__ import annotations

import json
import re
import sys
from typing import Any


def read_event() -> dict[str, Any]:
    try:
        value = json.load(sys.stdin)
    except (json.JSONDecodeError, OSError):
        return {}
    return value if isinstance(value, dict) else {}


def deny(reason: str) -> int:
    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "deny",
                    "permissionDecisionReason": reason,
                }
            }
        )
    )
    return 0


def main() -> int:
    event = read_event()
    tool_name = str(event.get("tool_name") or "")
    tool_input = event.get("tool_input")
    if not isinstance(tool_input, dict):
        return 0
    command = str(tool_input.get("command") or "")
    lowered = command.lower()

    if re.search(r"(?:^|[;&|\s])(?:pip|pip3)\s+install\b", lowered) or re.search(
        r"(?:python|python3)\s+-m\s+pip\s+install\b", lowered
    ):
        return deny("Package installation is prohibited for this analysis; use the configured release/runtime.")

    if re.search(r"(?:^|[;&|\s])voms-proxy-init\b", lowered):
        return deny("VOMS proxy initialization requires the user to act interactively.")

    if re.search(r"\bgit\s+add\s+(?:\.|-a|--all)(?:\s|$)", lowered):
        return deny("Do not stage the entire worktree; stage explicit logical file groups.")

    if tool_name == "apply_patch" and re.search(
        r"^\*\*\* (?:Add|Update|Delete) File: .*[/\\]\.claude[/\\]",
        command,
        re.MULTILINE,
    ):
        return deny("The .claude directory is Claude-owned and read-only for Codex.")

    write_like = re.search(
        r"(?:^|[;&|\s])(?:rm|mv|cp|touch|mkdir|rmdir|truncate|dd|chmod|chown|tee)\b|"
        r"\bsed\s+-i\b|\bperl\s+-pi\b|(?:^|\s)(?:>|>>)",
        lowered,
    )
    if ".claude/" in lowered and write_like:
        return deny("Commands that mutate .claude are blocked; read it directly and write Codex adapters elsewhere.")

    protected_data = (
        "/atlasgpfs01/usatlas/data/yuhanguo/pythia_private_sample/0",
        "/atlasgpfs01/usatlas/data/yuhanguo/pythia_truth_full_sample/",
        "/atlasgpfs01/usatlas/data/yuhanguo/powheg_full_sample/",
    )
    destructive = re.search(r"(?:^|[;&|\s])(?:rm|mv|truncate|dd)\b|(?:^|\s)(?:>|>>)", lowered)
    if destructive and any(path in lowered for path in protected_data):
        return deny("The command would mutate protected BNL analysis data.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
