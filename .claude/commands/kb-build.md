---
description: >
  Build or extend the physics/literature knowledge base under .claude/kb/.
  Use when the user asks to summarize a paper/source into the KB, add one or a
  few references to the KB, or bulk-build/reorganize the KB from a list of
  sources. Says "summarize this paper", "add to the knowledge base", "build the
  KB", "add these references". Auto-dispatched per CLAUDE.md.
---

You are the **orchestrator/executor** for KB building. A separate reviewer is run
via `/kb-review` afterward.

Task: $ARGUMENTS

## Setup (read first, in order)
1. `.claude/kb/KB_BUILDING_GUIDE.md` — the authoritative criteria. Obey it.
2. `.claude/kb/gotchas/reading_pdfs.md` — read PDFs with `gs`, never WebFetch/Read.
3. `.claude/agents/kb-builder.md` — per-source builder behavior.
4. `.claude/kb/index.md`, `.claude/kb/FUTURE_READ.md` — current state.
5. `Analysis/docs/analysis_overview.md` + `docs/tracking/analysis_roadmap_2026_06.md`
   — physics ground truth + which analysis steps/sections define "relevance."

## Pick the mode

- **ADD mode** — one or a few sources into an existing KB.
- **FULL mode** — bulk build/reorg from a provided list (roadmap step 2/3).

If a list of sources is given with intent to build the whole KB → FULL; otherwise
→ ADD. If ambiguous, ask the user.

### Logging (both modes)
- Slug from task (kebab-case). `LOG_FILE = .claude/logs/kb-build-YYYYMMDD-HHMMSS-<slug>.md`.
- Write the path to `.claude/logs/active-session.txt`.
- Append to `.claude/logs/tracking.jsonl`:
  `{"event":"start","command":"kb-build","mode":"ADD|FULL","task":"$ARGUMENTS","slug":"<slug>","timestamp":"<ISO 8601>"}`
- For FULL mode (many editing cycles) follow the tracking-doc protocol in CLAUDE.md.

---

## ADD mode

> **INVARIANT — no source is ever added in isolation.** Adding even ONE paper is
> not "write its doc and stop." Every ADD **must** also (a) integrate it into the
> KB — dedup against existing docs, add **bidirectional** graph links, designate/
> defer to hubs, update `index.md` + `FUTURE_READ.md`, reorganize if it creates a
> new area; and (b) **be reviewed** via `/kb-review` SINGLE (which checks exactly
> this integration). A doc that is written but not linked/deduped/reviewed is an
> **incomplete add** — the isolated-doc failure mode the whole system exists to
> prevent. Steps 2–4 below are mandatory, not optional polish.

For each source (do directly if 1–2; spawn one `kb-builder` subagent per source
if ≥3 to parallelize):
1. Per `kb-builder.md`: get PDF (PRIMARY → commit; SUPPORTIVE → URL), read with
   `gs`, classify, write the **one** source doc per GUIDE §4, scan refs (≤3).
2. **Orchestrator-only writes (never the subagent):** apply the proposed
   `index.md` line, add `[[...]]` graph links both ways, append future-read
   entries to `FUTURE_READ.md` (dedupe; enforce ≤3/source), and if the source
   duplicates a concept already in the KB, factor the shared part into a concept
   doc and link both (GUIDE §8).
   - **Bidirectional checklist (do not skip — this is the recurring miss):** for
     **every** doc the new source links TO, open that doc and add the **reciprocal
     back-link** — especially hubs, whose "analyses/sources that use this" lists
     must gain the new entry. A new doc that only links outward (one-way) is an
     isolated add. Use bare `[[basename]]` (GUIDE §8). When a source is now
     summarized, also **delist it** from any doc's future-read section and move it
     to FUTURE_READ "Done."
3. **Self-sync (GUIDE §9):** if this source covers an analysis step/topic not yet
   in the `/kb-review` full-KB test scenarios, add/extend the scenario.
4. Hand off to **`/kb-review` (single/few-source mode)**.

## FULL mode (bulk build / reorg)

Main agent (opus 4.8 / fable 5) is the **orchestrator**; subagents do isolated
per-source reads. This split exists to avoid concurrency conflicts on shared files.

**Phase 1 — parallel per-source summaries.**
- Partition the source list. Spawn **one `kb-builder` subagent per source**
  (cap concurrency ~4–6). Give each: the source citation/URL, its target **area
  directory**, the GUIDE, and the relevance targets (overview + roadmap).
- **Concurrency contract (state it in every subagent prompt):** the subagent
  writes **only its own source doc** (+ commits its PDF if PRIMARY) and
  **returns** proposed index line, graph links, future-read entries, and noticed
  duplications **as text**. It must **not** touch `index.md`, `FUTURE_READ.md`,
  concept docs, or any other source's doc, and must **not** reorganize anything
  beyond its own file.
- **Subagents are BLIND to each other — this is a DRAWBACK, not a feature.** The
  per-paper split buys parallelism and avoids concurrency conflicts on shared
  files, but the **price** is that each subagent reads only its one paper with no
  knowledge of siblings working papers on the same topic. So a subagent **cannot**
  detect cross-paper duplication, cannot know which other docs to link to, and
  cannot build the graph. Its "noticed duplications" can only flag overlap with
  docs that *already existed* before this build — not with sibling papers built in
  the same pass. **This blind spot is the central cost of the design and the
  orchestrator must consciously pay it down in Phase 2 / step 3** — it does not
  resolve itself; integration is the orchestrator's job and only the orchestrator
  has the global view to do it.

**Phase 2 — integration (orchestrator only, sequential). This is the heart of a
FULL build, not a cleanup pass.** Because the subagents were blind to each other
(Phase 1), **no cross-paper duplication, linking, or graph edge exists until you
create it here** — nobody else can or will. Be **hands-on and exhaustive**: read
every new doc yourself, build a topic map across all of them, and do not assume a
pair of papers is unrelated just because neither subagent mentioned the other
(they could not have). Concretely:
- **Read all new docs and group them by topic** (detector subsystem, observable,
  method, background, system/energy). This global view is the prerequisite for
  the rest — you cannot dedup or link papers you have not cross-read.
- **Dedup (exhaustive) — prefer DESIGNATING a hub over creating new docs.** For
  any concept recurring across ≥2 sources: **first check whether one existing doc
  is already the comprehensive treatment** and *designate it the canonical hub*
  (add a "this is the canonical hub for X" note + back-links to its spokes). Only
  **create a new `concepts/` doc when the concept spans ≥2 docs with no natural
  owner.** (In the 2026-06 build, designating `open_hf_production` and
  `ATLAS_Run2_muon_reconstruction` as hubs beat creating many small concept docs;
  exactly one new concept doc — `muon_source_template_fits` — was genuinely needed.)
  Then make each source doc **defer the shared mechanism to the hub via a pointer**
  while **retaining its own specifics** (provenance, numbers, systematics, per-paper
  usage) — this resolves the GUIDE §3-self-sufficiency vs §8-dedup tension. Check
  every topic group, not just overlaps a subagent happened to flag.
- **Hub back-links (critical):** a doc designated as a hub was often built *before*
  its spokes and will be a **one-way sink** — it links to nothing. The orchestrator
  MUST add the hub→spoke back-links (the blind builders could not).
  **Do not dedup from the summaries alone.** Per-paper summaries were written
  blind and may describe the *same* concept in different words or, worse, describe
  *subtly different* things that only look identical. When merging a concept across
  papers — or whenever the summaries are ambiguous or seem to conflict — **go back
  to the raw PDFs** (`gs` on the relevant pages) and confirm the papers really mean
  the same thing before collapsing them into one concept doc. Preserve genuine
  per-paper differences rather than flattening them.
- **Graph (exhaustive):** add **bidirectional** `[[...]]` links between every pair
  of related docs (same topic / method-and-its-application / measurement-and-its-
  comparison). **No orphans** — every doc must be reachable. Sibling papers on one
  topic must end up mutually linked even though their builders never saw each other.
- **FUTURE_READ.md:** merge all entries, dedupe, sanity-check global size.
- **Self-sync (GUIDE §9):** ensure `/kb-review` full-KB scenarios cover every
  analysis step and every source area now present.
- **Git (orchestrator-only, sequential):** subagents never commit. Stage `.md`
  normally; **`git add -f` the PDFs** (`*.pdf` is gitignored). Commit in a few
  logical groups (sources; then structure/index/concepts), not one mega-commit.
- **Hand-off to step 3:** this Phase-2 integration is the *first* pass; the
  dedicated directory/indexing/graph optimization (roadmap step 3) deepens it and
  inherits the same rule — **the orchestrator owns all cross-paper structure
  because the per-paper builders are blind to each other.** In step 3, especially
  when deduplicating, the orchestrator must again **return to the raw papers**
  (`gs`) when the blind-written summaries are insufficient to decide whether two
  treatments of a topic are truly the same.

> Directory structure + indexing design is roadmap **step 3**. In FULL mode,
> place files in the agreed area folders; if structure is not yet decided, stop
> and resolve step 3 first (do not invent a structure ad hoc).

**Phase 3 — review.** Hand off to **`/kb-review` (full-KB mode)**.

---

## Exit
- Update the log with files created/modified, classifications, future-read count.
- Append to `.claude/logs/tracking.jsonl`:
  `{"event":"end","command":"kb-build","mode":"ADD|FULL","docs":<n>,"timestamp":"<ISO 8601>","slug":"<slug>"}`
- Delete `.claude/logs/active-session.txt`.
- Tell the user which `/kb-review` mode to run next.
</content>
