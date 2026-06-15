# KB-Builder — per-source behavior

Read this file before summarizing a source into the knowledge base in a
build/reviewer loop. You summarize **one source** into one KB doc. You do **not**
reorganize the rest of the KB.

## Context to load
1. `.claude/kb/KB_BUILDING_GUIDE.md` — the authoritative criteria (follow exactly).
2. `.claude/kb/gotchas/reading_pdfs.md` — how to read PDFs here (`gs`, never WebFetch/Read).
3. `Analysis/docs/analysis_overview.md` — physics ground truth (to judge relevance).
4. `Analysis/docs/tracking/analysis_roadmap_2026_06.md` — which analysis steps/sections exist (relevance targets).
5. `.claude/kb/index.md` — existing docs (to link to, and to avoid duplicating).

## Procedure (one source)
1. **Get the PDF.** PRIMARY → `curl -sSL <url> -o <area-dir>/<file>.pdf` (just
   **write the file — do NOT git-add/commit it**; the orchestrator commits, with
   `git add -f` since `*.pdf` is gitignored). SUPPORTIVE → record the URL; `curl`
   to `/tmp` only if you need to read it.
   - **arXiv works** (`https://arxiv.org/pdf/<id>`). **Non-arXiv hosts may block
     `curl`** (e.g. indico/CDS sit behind an Anubis anti-bot proof-of-work →
     `curl` returns an HTML challenge, not the PDF). If blocked: record the URL,
     summarize only what you can legitimately establish, and **never fabricate**
     content you could not read (GUIDE §5). Flag it for the orchestrator.
2. **Read it with ghostscript** (`gs -sDEVICE=txtwrite`, page ranges). Never use
   WebFetch or the Read tool on the PDF. `gs` may segfault on an odd page — extract
   per-page-range and cover gaps from adjacent pages; note any gap.
3. **Classify** PRIMARY vs SUPPORTIVE per GUIDE §2.
4. **Write the summary** to its doc using the GUIDE §4 template: relevance-first
   (§1), self-sufficient (§3), no physics assumptions / condition-difference
   warnings (§5), key tables/numbers inline (§3.4), concision (§7).
5. **Scan the reference list** (GUIDE §6): list ≤3 future-read refs with
   class/new-info/serves, or state "None worth adding."

## Hard limits (concurrency safety)
- Write **only your one source doc** (and, if PRIMARY, write — not commit — its PDF).
- **NEVER run any git command** (no add/commit/reset). Concurrent subagent commits
  collide on `.git/index.lock` and sweep up siblings' files — the orchestrator does
  **all** git, sequentially, after the parallel phase. (Learned the hard way in the
  2026-06 bulk build.)
- **Do NOT** edit `index.md`, `FUTURE_READ.md`, concept docs, or any other
  source's doc. **Return** proposed index line, graph links, and future-read
  entries **as text** to the orchestrator, which owns those shared files.
- **Link convention:** write graph links as `[[basename]]` (no path, no `.md`) —
  e.g. `[[open_hf_production]]`, not `[[physics/heavy_ion/open_hf_production.md]]`.

## No-hallucination rule (GUIDE §5)
Attribute every claim/number to a page/table/eq. Mark uncertain claims
`⟨unverified — check PDF p.N⟩`. Never speculate on the physics implication of a
methodology difference or the size/direction of a condition difference.

## Output (return to orchestrator)
- Source doc path; classification; whether PDF committed (path).
- Proposed `index.md` line; proposed `[[...]]` graph links to related docs.
- Future-read entries (≤3, or "None") for `FUTURE_READ.md`.
- Any duplication you noticed (a concept already covered elsewhere → candidate
  for a shared concept doc).
</content>
