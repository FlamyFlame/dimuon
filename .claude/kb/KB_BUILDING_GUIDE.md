# Knowledge-Base Building Guide (authoritative criteria)

This is the **constitution** for the literature/physics knowledge base under
`.claude/kb/`. The `/kb-build` (executor) and `/kb-review` (reviewer) commands
operationalize it. If anything in those commands or in a KB doc contradicts this
guide, this guide wins — flag the conflict.

> **Scope.** This guide governs *how sources are summarized and organized*.
> The concrete directory layout and indexing/graph design are finalized in
> roadmap **step 3**; §8 here states the principles those must satisfy.

---

## 0. Why this KB exists (the three jobs)

Every entry must serve at least one of these, and say which:

1. **The analysis itself, at publication rigor** — methods/numbers we use or
   must reproduce (efficiencies, fit forms, selection definitions, normalization
   inputs).
2. **Journal-paper / Internal-Note writing** (e.g. PLB, ANA-HION-2023-07-INT1) —
   motivation, prior measurements to cite/compare, methodology phrasing.
3. **PhD thesis writing** — broader field background and pedagogical context.

The KB is a **working reference for THIS analysis**, not a generic literature
archive. Analysis context lives in `docs/analysis_overview.md` (physics ground
truth) and `docs/tracking/analysis_roadmap_2026_06.md` (status); read those to
judge relevance.

---

## 1. Relevance-first principle (generic vs analysis-relevant)

Default to **analysis-relevant depth**, not generic completeness. For any source:

- Summarize a generic topic **only to the depth it serves one of the three jobs**.
- Explicitly tag each block as either **`[method-we-use]`** (something the
  analysis or its writing directly relies on — summarize thoroughly) or
  **`[background-for-writing]`** (intro/motivation/thesis context — summarize
  briefly, enough to cite and paraphrase).
- A bulky generic recap that serves none of the three jobs does not belong in
  the KB. The longer the KB, the harder it is for a task agent to find the one
  thing it needs — **brevity is a correctness property, not a style preference.**

---

## 2. Source classification: PRIMARY vs SUPPORTIVE

Classify every source. Classification drives storage (§3) and search weight.

- **PRIMARY** — high relevance/referenceability to this analysis: the analyses
  we directly derive upon; performance papers whose methods or numbers we use;
  key field reviews we will cite heavily; the measurements we compare against.
  Prioritized in KB search/recall.
  → **Requires a committed local PDF (§3) + a thorough, self-sufficient summary
  + the URL.**
- **SUPPORTIVE / complementary** — genuinely useful (else excluded) and
  complements a primary on the same topic, but either less relevant to this
  analysis or less comprehensive (e.g. narrower scope than a review).
  → **URL to the PDF + a summary of ONLY the relevant slice.** Local PDF
  optional (store only if small and likely to be re-consulted).

Notes:
- A reference **discovered through** a supportive source can itself be PRIMARY if
  highly relevant. Classification follows relevance, not how it was found.
- A source relevant for **one narrow reason only** (e.g. just a muon-pair
  reco-efficiency method): summarize **only that slice**, keep the rest generic
  or omit it, and **link it from the relevant topic** in the index/graph so it is
  found when (and only when) that topic comes up.

---

## 3. The URL-vs-PDF rule (decided empirically — do not relitigate)

**Finding (2026-06-14, this cluster):**
- `WebFetch` on an arXiv **abstract** page returns only the abstract — no tables,
  cuts, or fit forms.
- `WebFetch` on an arXiv **PDF URL** returns undecoded binary — **content is not
  retrievable**.
- The **Read** tool on a local PDF **fails**: rendering needs `pdftoppm`
  (poppler), which is **not installed and cannot be installed** (no pip/apt).
- `gs -sDEVICE=txtwrite` (ghostscript 9.54, present) **does** extract usable text
  (minor ligature/spacing artifacts; tables/figures extract poorly).
- Outbound HTTPS works (so a PDF can be `curl`-ed on demand), but **WebFetch will
  still not parse it** — only `gs` will.

**Consequences (mandatory):**
1. **A URL is NOT a substitute for content.** No agent can pull a specific
   number/method from a URL via WebFetch. Therefore **every summary must be
   self-sufficient for the analysis-relevant content** — an agent should be able
   to answer the common analysis/writing question *without opening the PDF*.
2. **PRIMARY refs commit the PDF in-repo** (next to the summary). The PDF is the
   **fallback for rare deep details** the summary omits.
3. **To read any PDF, use ghostscript, never WebFetch/Read:**
   ```
   gs -q -dNOPAUSE -dBATCH -sDEVICE=txtwrite \
      -dFirstPage=N -dLastPage=M -sOutputFile=out.txt input.pdf
   ```
   For a supportive ref with no committed PDF, `curl` it to `/tmp` first, then
   `gs`. (See `gotchas/reading_pdfs.md`.)
4. Capture key **tables/numbers in the summary text**, because `gs` extracts
   tables badly — the PDF fallback is weakest exactly where numbers live.
5. **arXiv `curl` works; many non-arXiv hosts do not.** indico/CDS sit behind an
   anti-bot proof-of-work, so `curl` returns an HTML challenge, not the PDF. When a
   source cannot be fetched, record the URL, summarize only what you can
   legitimately establish, and **never fabricate** unread content — add a citeable
   arXiv anchor as a future-read instead.

---

## 4. Mandatory summary structure (one file per source)

```markdown
# <Short title>

**Source:** <full citation> — <journal/report no.>
**arXiv / DOI:** <id>
**PDF:** `./<file>.pdf` (PRIMARY)   |   **URL:** <url> (SUPPORTIVE)
**Classification:** PRIMARY | SUPPORTIVE
**Added:** YYYY-MM-DD

## Relevance to this analysis   (REQUIRED — specific, not vague)
| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| <e.g. tag-and-probe reco-eff method> | <e.g. pp24 reco-eff, roadmap step 4 / IntNote §12> | [method-we-use] |
| <e.g. b-quark energy-loss motivation> | <Intro / thesis ch. 1> | [background-for-writing] |

## Scope & condition-difference warnings
<system, √s, dataset, run period>. Then explicit warnings, e.g.:
"Run 2 pp 13 TeV; our analysis is Run 3 5.36 TeV with HIJING overlay —
detector conditions and performance differ. ACKNOWLEDGE the difference; DO NOT
assume its size or direction (§5)."

## Content summary
<structured, relevance-first per §1; key tables/numbers inline per §3.4>

## References worth future reading   (§6; ≤3, or "None")

## Related KB docs   (knowledge graph; §8)
- [[<other-doc>]] — relation
```

The **Relevance** section is the single most important part. "Relevant to muon
efficiency" is **not** acceptable; name the concrete analysis step, the
IntNote/thesis section, or the observable, and the use type.

---

## 5. No hallucination, no physics assumptions (hard rule)

- **Never invent physics implications.** If another analysis uses a different
  method/binning/sample than ours, **do not** speculate on what that difference
  implies for our measurement.
- **Condition/performance differences** (Run 2↔Run 3, pp↔PbPb, √s, pile-up,
  occupancy): **acknowledge and warn** future agents, but **do not estimate the
  size or direction** of the difference unless the source itself states it.
- **Attribute everything.** Distinguish "stated in the source" from "our analysis
  context." Numbers require a traceable anchor (page / table / equation / figure).
- When unsure whether a claim is in the source, mark it `⟨unverified — check
  PDF p.N⟩` rather than asserting it.

---

## 6. Reference scanning & the Future-Read list (the ≤3 rule)

For every source, the summarizer **must scan its reference list** and identify
references offering **new, crucial, analysis-relevant info not already covered**
by the current source.

- **Cap ≤ 3 per paper** added to `FUTURE_READ.md`. Keep the KB lean: a bulkier
  KB is more duplicative and harder to search.
- If **none** are worth adding, **state so explicitly** ("No new references worth
  adding") — the reviewer checks this.
- For each future-read entry record: citation, tentative PRIMARY/SUPPORTIVE,
  one line on *what new info* it adds and *which part of our work* it serves.
- Record each entry **both** in the source's own KB doc (§4 "References worth
  future reading") **and** appended to `FUTURE_READ.md`.

The reviewer independently checks the source's references to (a) confirm listed
future-reads are genuinely worth it and (b) catch important refs the summarizer
missed.

---

## 7. Concision vs detail (what to keep, what to compress)

- **Keep:** numbers/tables/fit-forms/selections we use or cite; definitions;
  conventions; methodology we adopt; results we will compare against.
- **Compress:** generic narrative, unrelated channels/results, history.
- **Test:** a summary passes if a task agent can answer the *common* analysis or
  writing question from the summary alone; the PDF is only for rare deep dives.

---

## 8. Placement, dedup & knowledge graph (principles; structure finalized in step 3)

- **One folder per area** (e.g. analysis papers / detector-performance /
  heavy-ion-field reviews / methods). Default **one file per source**.
- **Dedup shared concepts — prefer a designated HUB over a new doc.** When a
  concept recurs across sources (a detector subsystem, b-quark energy loss, a
  method), first check whether one existing doc already treats it comprehensively
  and **designate that doc the canonical hub** (mark it as such; other docs defer
  to it). Only **create a new `concepts/` doc when the concept spans ≥2 docs with
  no natural owner.** Either way: explain once, and have each source doc **defer the
  shared mechanism via a pointer while keeping its own specifics** (provenance,
  numbers, systematics, per-paper usage) — this is how §3 self-sufficiency and §8
  dedup coexist.
- **Every doc** is registered in `index.md` and **cross-linked** to related docs so
  an agent finds **all** relevant info, not just the first hit. **No orphan docs;
  links are bidirectional** — in particular, a doc designated as a hub is often
  built before its spokes and starts as a one-way sink, so back-links into it must
  be added during integration. **Link syntax: `[[basename]]`** (no path, no `.md`).
- The structure must make information **fast to find** (clear, indexed) and
  **hard to miss** (graph-linked, deduped). These are the acceptance criteria the
  step-3 structure and the `/kb-review` full-KB mode test.

---

## 9. Self-sync (the KB and its tests are living)

When (a) the analysis gains steps/sources (e.g. unfolding, template fitting,
systematics) or (b) new sources are added to the KB, then:
- the **Relevance** mappings of affected docs, and
- the **use-case test scenarios** in `/kb-review` (full-KB mode)

**must be updated to include them.** A source added without updating the relevant
test scenario, or an analysis step added without a covering scenario, is an
incomplete change and should be flagged by the reviewer.
</content>
</invoke>
