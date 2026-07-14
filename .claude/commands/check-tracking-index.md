---
description: >
  Check and maintain Analysis/docs/tracking/INDEX.md — the tracking-doc index that
  Doc triage reads first. Verifies every doc is listed exactly once, statuses and
  Updated dates are right, sections are newest-first, and each scope line still
  describes its doc. Use when the user says "check the index", "check tracking
  index", "are the tracking scopes current", "maintain INDEX.md", after closing or
  creating a tracking doc, or when triage seems to be missing a relevant doc.
---

Arguments: $ARGUMENTS
- (empty) → check, apply fixes, report
- `check` → report only; make **no** edits
- `all` → deep-check every doc (expensive; only when scopes are broadly suspect)
- `<doc.md> [...]` → deep-check only the named docs

`check` and `all` are keywords only as bare words. An argument ending in `.md` is
always a doc name, never a keyword.

Index: `Analysis/docs/tracking/INDEX.md`. Docs: `Analysis/docs/tracking/*.md`.
Never commit — leave that to the user or `/wrap-up`. Never delete a tracking doc.

**Why this exists:** the scope column is the *only* thing Doc triage sees
(CLAUDE.md §Tracking Documents). A doc whose scope line is stale is invisible to
every future conversation — worse than one that was never indexed, because
nothing looks wrong. Scope accuracy, not row count, is the deliverable.

---

## 1. Structural checks (cheap — all docs, always)

Run these first; they need no doc reading.

1. **Coverage, 1:1.** The coverage universe is every `*.md` in
   `Analysis/docs/tracking/` **except `INDEX.md` and `_sub_*.md`** (scratch docs —
   see §1.4). Every doc in that universe appears in exactly one row. Report rows
   with no file (orphans) and files with no row (missing). A missing file is the
   serious one — it is unreachable by triage. Never add a row for a `_sub_*.md`.

2. **Section integrity.** Exactly two sections, `## Active` and `## Closed`.
   Rows sorted `Updated` descending within each. No doc in both.

3. **Updated dates.** For each doc compute its **effective last-changed date**:

   ```
   git status --porcelain -- <file>
   git log -1 --format=%ad --date=short -- <file>
   ```

   - If `git status --porcelain` prints **any** entry at all — whatever the code
     (`??` untracked, ` M` dirty, `M ` staged, `MM`, `A ` added, `R ` renamed, …) —
     the doc is **not clean**: its effective date is **today**, and it **always
     enters the deep-check set** (§2). Do not enumerate codes and do not treat
     staged-but-working-tree-clean as unchanged; any porcelain output means the
     content moved after it was last indexed, which is exactly when a scope goes
     stale.
   - Only a doc that is **tracked and clean** (no porcelain output) uses the
     committed `git log` date.
   - Never write an empty value into `Updated`, and never treat an empty `git log`
     as "unchanged" — it means untracked.

   Flag any row whose `Updated` is older than the doc's effective date.

   This matters most on the `/wrap-up` path: wrap-up delegates here in §1a but
   commits only in §2, so the session's own edits are still uncommitted when this
   command runs. Using the committed date there would bump `Updated` backwards and
   skip re-scoping the very docs the session just grew.

4. **Leftover scratch docs.** Any `_sub_*.md` is a subagent scratch doc that
   should have been merged and deleted (CLAUDE.md §Delegated subagent memory).
   Do not index them; report them for the user to merge or remove.

## 2. Deep check — scope accuracy (prioritized, never all-by-default)

Reading every doc would cost what triage was built to avoid. Select:

- **All ACTIVE docs** — they govern current work.
- **Any doc on disk with no row** (§1.1 "missing") — it has no scope at all, so one
  must be derived before it can be indexed.
- **Any doc flagged stale by §1.3** — including CLOSED ones, and including docs
  merely dirty in the working tree.
- **The 3 most recently updated CLOSED docs** — a rolling audit of the rest.

Skip the remaining CLOSED docs unless `all` or explicit doc names were passed.

For each selected doc, read **only** what is needed to judge its scope — not the
whole file:
- line 1 (title) and the `## Objective` block;
- `Latest Stage` (current thread, or empty if complete);
- the last ~30 lines of `Progress Log` / `Accumulated Findings` (what actually
  happened most recently).

**Not every doc has those headings** (e.g. `analysis_roadmap_2026_06.md` is
organized as Q1–Q5 + a status ledger; `analysis_status_summary.md` has only an
Objective; `event_selection_banana_cut_comparison.md` is slide-structured). When a
heading is absent, read the section that plays its role — the doc's opening
statement of purpose, and its most recent results/conclusions. If the doc is small
(< ~15 kB), just read it whole.

**Context cap (hard).** Within one run, read any single doc at most once in full;
never re-read a doc to re-derive a scope.

The ACTIVE docs are the always-read baseline — they are never the reason to abort.
Neither are dirty or unindexed docs: those are targeted, individually cheap, and
are precisely the ones known to have moved, so never skip them for cost.

The cap governs only the **rolling-audit CLOSED sample** (the 3 above). Before
reading anything, `wc -c` both sets. If that sample's full-file byte total exceeds
**half the ACTIVE set's full-file byte total**, shrink the sample (drop to the
single most recent CLOSED doc, or none) and say so in the report. Both sides are
full-file bytes from `wc -c`, computed up front — never compare a partial-excerpt
estimate against a full-file total.

If the dirty set is itself very large — many CLOSED docs touched at once — read
them anyway, but note the total in the report so the cost is visible rather than
silent. Correctness beats cost: a doc known to have changed must never go
unre-scoped.

Nothing is written before §3, so stopping or shrinking here leaves INDEX.md
untouched.

Then ask, in order:

1. **Does the scope line still describe the doc?** A scope goes stale when the
   doc grows a sub-thread, changes method, or reaches a conclusion that reverses
   its original objective. The scope must say what the doc *concluded*, not only
   what it set out to do.
2. **Is the status right?** An ACTIVE doc with a final summary and empty Latest
   Stage should be CLOSED. A CLOSED doc with an open Latest Stage or unfinished
   Remaining Work should be ACTIVE, or PARKED if it is blocked on external
   inputs — say what unblocks it.
3. **Would triage find this doc?** Read the scope line cold, as an agent who has
   read nothing else. Does it name the physics terms, sample names, pipelines, or
   observables someone would ask about? A scope that is accurate but unsearchable
   ("various follow-ups") fails this check. Rewrite it with the nouns a future
   request would use.
4. **Does it carry what only it knows?** If the doc established a standing
   invariant, a superseded default, or a live open item (e.g. "fix not applied,
   awaiting user decision"), the scope line must say so — closed docs are read by
   scope alone unless something pulls a reader in.

## 3. Apply

Unless `check` was passed:

- Rewrite stale scope lines; fix statuses; bump `Updated` to each doc's
  **effective date** from §1.3; re-sort each section newest-first.
- **Add a row for every doc found on disk without one** (§1.1 "missing"), in the
  right section, with a freshly derived scope from §2 and its effective date. This
  is the whole point of the coverage check — detecting an unindexed doc and
  leaving it unindexed fixes nothing.
- **Remove orphan rows** whose file no longer exists, and name them in the report.
- Preserve the umbrella / sub-doc annotations (`low_mass_dimuon_template_fit.md` A–E).

If no row changed, do not rewrite the file and take no backup — a second
consecutive run is then a true no-op.

Moving rows between sections and re-sorting is a multi-row mutation of a
load-bearing file. Immediately before an actual write, **copy INDEX.md to
`INDEX.md.bak`** — INDEX.md is often uncommitted (it changes in the same sessions
that change the docs), so `git checkout` is not a reliable restore path. Then
**write INDEX.md once, as a whole file**, from the fully reconciled row set — never
as a sequence of partial edits, so an interruption cannot leave the index
half-sorted or a row dropped. (`INDEX.md.bak` does not end in `.md`, so it never
enters the §1.1 coverage universe. A stray `.bak` from an interrupted run is
overwritten by the next apply.)

Do not invent content. If a doc's real scope cannot be determined from the
sections read, read the rest of that doc once (within the §2 cap) — do not guess,
and do not leave the old line in place hoping it is still true.

## 4. Report

**Before reporting, re-run the §1.1 coverage check against the file you just
wrote.** Assert N rows == N docs in the coverage universe, and that every doc that
had a row before still has exactly one. Never emit a verdict from the pre-apply
check — a row can be lost in the rewrite, and that is precisely the failure this
command exists to catch.

**If the post-apply check fails, do not emit a verdict.** Distinguish the two
causes — they need opposite responses:

- **A row was lost or duplicated by the rewrite** (it existed in the reconciled
  set but not in the file). This is corruption. Rewrite INDEX.md once more from
  the reconciled set and re-verify. If it still fails, restore `INDEX.md.bak` over
  INDEX.md and report **INDEX: ISSUES** naming the exact rows. A reported failure
  with the old index intact is recoverable; a silently dropped row is not.
- **A doc legitimately has no row and none could be derived** (e.g. §2 could not
  determine its scope), or a doc's file **disappeared mid-run** so its row is now
  an orphan. Do **not** roll back — that would discard every correct scope and
  date fix in the same file. Keep the applied fixes and report **INDEX: ISSUES**
  naming the doc or row and why.

On success, delete `INDEX.md.bak` if it exists. On failure, leave it and say so.

In `check` mode nothing is written and no backup is taken, so §4 collapses to the
report bullets: never delete a stray `INDEX.md.bak` there — `check` makes no edits.

- **Coverage:** N rows / N docs, 1:1 or the specific orphans + missing.
- **Deep-checked:** which docs, and why each was selected (active / stale / rolling).
- **Scopes rewritten:** doc → one-line before/after, with the reason.
- **Status changes:** doc → Active↔Closed↔Parked, with the evidence.
- **Dates bumped:** doc → old → new, and whether the date came from git or from
  an uncommitted working-tree change.
- **Scratch docs found:** list, or "none".
- **Not checked:** the CLOSED docs skipped — say how many. Never imply full
  coverage when §2 sampled.

End with **INDEX: CLEAN** (1:1 verified after apply, all checked scopes accurate,
statuses right) or **INDEX: FIXED** / **INDEX: ISSUES** with the specific list.
