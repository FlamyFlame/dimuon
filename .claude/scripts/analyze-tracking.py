#!/usr/bin/env python3
"""Analyze executor-reviewer loop usage from tracking.jsonl.

Run from repo root:
    python .claude/scripts/analyze-tracking.py
    python .claude/scripts/analyze-tracking.py --since 2026-05-01
    python .claude/scripts/analyze-tracking.py --command review-plot
"""
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

def load_events(log_path, since=None, command_filter=None):
    if not log_path.exists():
        return []
    events = []
    for line in log_path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            e = json.loads(line)
        except json.JSONDecodeError:
            continue
        if since and e.get("timestamp", "") < since:
            continue
        if command_filter and e.get("command") != command_filter:
            continue
        events.append(e)
    return events


def print_summary(events):
    starts = [e for e in events if e["event"] == "start"]
    ends = [e for e in events if e["event"] == "end"]

    print(f"Total invocations: {len(starts)}")
    print(f"Completed: {len(ends)}")
    orphans = len(starts) - len(ends)
    if orphans > 0:
        print(f"Incomplete (started but no end event): {orphans}")
    print()

    if not starts:
        return

    print("=" * 60)
    print("BY COMMAND")
    print("=" * 60)
    for cmd, count in Counter(e["command"] for e in starts).most_common():
        cmd_ends = [e for e in ends if e["command"] == cmd]
        approved = sum(1 for e in cmd_ends if e["status"] == "approved")
        escalated = sum(1 for e in cmd_ends if e["status"] == "escalated")
        iters = [e["iterations"] for e in cmd_ends if "iterations" in e]
        avg_iter = sum(iters) / len(iters) if iters else 0
        approval_rate = (approved / len(cmd_ends) * 100) if cmd_ends else 0

        print(f"\n  {cmd}:")
        print(f"    Runs: {count}")
        print(f"    Approved: {approved} | Escalated: {escalated}")
        print(f"    Approval rate: {approval_rate:.0f}%")
        print(f"    Avg iterations to converge: {avg_iter:.1f}")
        if iters:
            print(f"    Iteration distribution: {sorted(iters)}")

    print()
    print("=" * 60)
    print("ESCALATION ANALYSIS")
    print("=" * 60)
    escalated = [e for e in ends if e["status"] == "escalated"]
    if not escalated:
        print("  No escalations.")
    else:
        all_issues = []
        for e in escalated:
            all_issues.extend(e.get("unresolved", []))
        if all_issues:
            print("  Most common unresolved issues:")
            for issue, count in Counter(all_issues).most_common(10):
                print(f"    {count}x: {issue}")
        print(f"\n  Escalated tasks:")
        for e in escalated:
            print(f"    - [{e['command']}] {e.get('slug', 'unknown')} "
                  f"({e['iterations']} iters, "
                  f"{e.get('criticals', '?')}C/{e.get('warnings', '?')}W)")

    print()
    print("=" * 60)
    print("TIMELINE (last 10)")
    print("=" * 60)
    recent = sorted(events, key=lambda e: e.get("timestamp", ""))[-10:]
    for e in recent:
        ts = e.get("timestamp", "?")[:16]
        if e["event"] == "start":
            print(f"  {ts}  START  {e['command']}  {e.get('slug', '')}")
        else:
            status = "OK" if e["status"] == "approved" else "ESCALATED"
            print(f"  {ts}  {status:9s}  {e['command']}  "
                  f"({e.get('iterations', '?')} iters)")


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--since", help="Filter events after this date (YYYY-MM-DD)")
    parser.add_argument("--command", help="Filter to a specific command")
    parser.add_argument("--log", default=".claude/logs/tracking.jsonl",
                        help="Path to tracking.jsonl")
    args = parser.parse_args()

    log_path = Path(args.log)
    events = load_events(log_path, since=args.since, command_filter=args.command)

    if not events:
        print("No tracking data found.")
        if not log_path.exists():
            print(f"  (File does not exist: {log_path})")
        elif args.since or args.command:
            print(f"  (Try without --since/--command filters)")
        sys.exit(0)

    print_summary(events)


if __name__ == "__main__":
    main()
