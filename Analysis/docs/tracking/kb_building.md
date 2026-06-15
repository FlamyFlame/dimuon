# Knowledge-Base Building System

**Type:** Implementation (infrastructure/tooling, **not** a physics measurement).
There is no Physics Procedure section because this builds the KB *machinery*, not
a measurement/correction. **In its place, the authoritative reference is
`.claude/kb/KB_BUILDING_GUIDE.md`** — all design decisions, command behavior, and
review criteria must follow from it; if anything contradicts it, flag and fix the
guide first (with user approval).

## Objective

Build a physics/literature knowledge base under `.claude/kb/` that serves, at
publication rigor: (1) the analysis itself, (2) journal-paper + Internal-Note
writing, and (3) PhD-thesis writing. Plus the tooling to build, review, and
maintain it. Four roadmap steps (below).

## Authoritative Reference (in lieu of Physics Procedure)

`.claude/kb/KB_BUILDING_GUIDE.md` (the "constitution"): relevance-first; PRIMARY
vs SUPPORTIVE classification; URL/PDF rule; mandatory summary template; no
hallucination / no physics assumptions; ≤3 future-reads per source; dedup +
knowledge graph; self-sync. The `/kb-build` and `/kb-review` commands operationalize
it.

## Context

- KB lives at `.claude/kb/` (existing: index, analysis/{overview,decisions,
  run2_dimuon_note}, data/{samples,variables}, physics/detector/ATLAS muon reco
  + its PDF, procedures, gotchas).
- Relevance ground truth: `Analysis/docs/analysis_overview.md` (physics) +
  `docs/tracking/analysis_roadmap_2026_06.md` (status / which steps exist).
- External reference for structure: Shuonli ppg12 wiki
  (https://github.com/Shuonli/ppg12/tree/main/wiki) — inspiration for step 3.

## Scope (the four steps)

1. **Design phase** — criteria + `/kb-build` + `/kb-review` skills. *(this doc's
   current step)*
2. **Bulk build** — from a user-provided paper list (analysis-derived / ATLAS
   detector & muon perf / HI-field reviews / specific methods).
3. **Indexing & structural optimization** — directory structure, indexing system,
   dedup, knowledge graph.
4. **Review & self-improvement** — of the KB-building procedure itself.

## Design Decisions (with rationale)

- **URL is NOT sufficient → PDF for PRIMARY refs.** Empirically (2026-06-14, this
  cluster) WebFetch returns only the abstract (abs page) or undecoded binary (PDF
  URL); the Read tool fails (no poppler); only `gs -sDEVICE=txtwrite` extracts
  text. ⇒ summaries must be self-sufficient; PRIMARY refs commit the PDF in-repo
  (deep-dive fallback via `gs`); SUPPORTIVE refs carry the URL + relevant-slice
  summary. (GUIDE §3; gotcha `reading_pdfs.md`.)
- **PRIMARY vs SUPPORTIVE classification** drives storage + search weight; a ref
  found via a supportive one can itself be PRIMARY. (GUIDE §2.)
- **Relevance-first, tagged** `[method-we-use]` / `[background-for-writing]`;
  Relevance section must name a concrete step/section/observable, never vague.
- **No physics assumptions:** acknowledge but never quantify Run2↔Run3 / pp↔PbPb
  condition differences. (GUIDE §5.)
- **≤3 future-reads per source** (or explicit "none") to keep the KB lean. (§6.)
- **Bulk orchestration = one subagent per paper**, orchestrator owns shared files
  (index/FUTURE_READ/concept docs). Rationale: parallelism + no concurrency
  conflicts. **Known DRAWBACK (the price):** subagents are blind to each other →
  cannot do cross-paper dedup/linking/graph; the orchestrator must pay this down
  exhaustively in Phase 2 / step 3, returning to **raw PDFs** when blind-written
  summaries are insufficient to decide if two treatments are truly the same.
- **kb-review FULL mode = structural audit + blind use-case "unit tests"**: fresh
  test agents (no hint it's a test) run scenarios (implement code / investigate /
  write section / review results); a *separate* evaluator grades Discovery /
  Efficiency / Comprehensiveness (avoids orchestrator self-grading bias).
- **Auto-trigger** added to CLAUDE.md, scoped to explicit "summarize/add/build KB".
- **Self-sync:** scenarios + relevance mappings must update as analysis steps and
  sources grow. (GUIDE §9.)

## Implementation Plan

1. [DONE] Design: GUIDE + `/kb-build` (+ `kb-builder` agent) + `/kb-review`;
   CLAUDE.md auto-dispatch + PDF gotcha; index/FUTURE_READ stubs; memories.
2. [TODO] Bulk build from user's paper list → `/kb-build` FULL → `/kb-review` FULL.
   Blocking input: the paper list (grouped).
3. [TODO] Directory structure + indexing + knowledge graph (Shuonli wiki as ref);
   re-review FULL mode incl. structure + use-case sims.
4. [TODO] Review & self-improve the KB-building procedure (skills, guidelines).

## Progress Log (append-only)

- **2026-06-14 — Step 1 DONE.** Tested PDF retrieval empirically (decisive: URL
  not runtime-retrievable; `gs` works). Wrote `.claude/kb/KB_BUILDING_GUIDE.md`,
  `.claude/commands/kb-build.md`, `.claude/commands/kb-review.md`,
  `.claude/agents/kb-builder.md`, `.claude/kb/FUTURE_READ.md`,
  `.claude/kb/gotchas/reading_pdfs.md`. Updated `.claude/kb/index.md` and
  `.claude/CLAUDE.md` (auto-dispatch for `/kb-build` + `/kb-review`; PDF gotcha).
  Saved memories `project_kb_building`, `reference_reading_pdfs`.
- **2026-06-14 — Step 2 batch 1 DONE** (5 subagents, parallel). Created:
  `analysis/run2_dimuon_backtoback_paper.md` (+2308.16652.pdf; updated
  `run2_dimuon_note.md` with publication-level framing + cross-link),
  `analysis/run2_hf_muon_raa.md` (+2109.00411.pdf),
  `analysis/run2_hf_muon_vn.md` (+2003.03565.pdf),
  `physics/heavy_ion/hi_big_picture.md` (+1802.04801.pdf),
  `physics/heavy_ion/open_hf_production.md` (+1903.07709.pdf).
  - **Dedup signals (for step 3 concept docs):** ρ=Δp/p momentum-imbalance
    fake-muon (π/K) template fit (note+raa+vn); d0 c/b template fit (raa+vn);
    tag-and-probe muon eff (↔ detector/muon_reco); FCal/Glauber centrality +
    event-plane (heavy_ion docs); R_AA & v_n definitions (overview+both HI
    reviews); per-pair dimuon trig+reco eff (backtoback+note).
  - **Future-read candidates (raw, pre-dedup):** 2206.12594 (γγ→μμ non-UPC PbPb,
    PRIMARY — per-pair trig+reco eff method); 1805.05220 (HF μ 2.76 TeV R_AA+v2,
    ×3, SUPPORTIVE); 1909.01650 (c/b μ v2 pp13, origin of d0 method); 1705.01974
    (jets in HI RMP, SUPP); Andronic EPJC76:107 (HF+quarkonium review, SUPP);
    Beraudo NPA979 (HF E-loss models, SUPP); CMS bottom R_AA Sirunyan (PRIMARY
    comparison data). → dedup in step 3.
  - **CONCURRENCY BUG (for step 4):** subagents were told to commit their PDFs;
    concurrent `git add/commit` collided on `.git/index.lock` and one swept up
    siblings' staged files (recovered via `git reset`). FIX: subagents must NOT
    run git at all — only write files; the **orchestrator does all git**. Also
    `*.pdf` is gitignored → orchestrator uses `git add -f`. Update
    `kb-builder.md` + `kb-build.md` accordingly in step 4.

- **2026-06-14 — Step 1 refinement.** Per user feedback, made the per-paper
  "subagents blind to each other" explicit in `kb-build.md` as a **drawback/price**
  (not a feature) requiring exhaustive orchestrator integration in Phase 2/step 3;
  required dedup to **return to raw PDFs** when summaries are ambiguous; added a
  matching "primary risk to hunt" note to `kb-review.md` FULL Part A.

## Results & Observations

- **PDF handling on this cluster (verified 2026-06-14):** WebFetch(abs)=abstract
  only; WebFetch(pdf URL)=binary, not retrievable; Read tool=fails (no
  pdftoppm/poppler, can't install); `gs -sDEVICE=txtwrite` (gs 9.54)=usable text
  (tables/figures extract poorly); outbound HTTPS works (curl 200, no proxy).
- Existing convention already pairs a local PDF with its `.md` summary
  (`physics/detector/ATLAS_Run2_muon_reco_eff_EPJC81_578.pdf`) — consistent with
  the decided PRIMARY-ref policy.

## Remaining Work

- Steps 2–4 (above). Step 2 needs the user's grouped paper list. Step 3 should
  consult the Shuonli ppg12 wiki and may reshuffle directories/indexing.

## Latest Stage

**Steps 2–4 in progress (2026-06-14).** Approach chosen by user: **dogfood** —
orchestrate the build myself (one subagent per paper), then improve `/kb-build` +
`/kb-review` from the experience (feeds step 4).

### Step-2 source plan (classification + target dir)
Analysis (PRIMARY) → `kb/analysis/`:
- A. 2308.16652 Run2 dimuon back-to-back az. corr (journal/PRL). Also re-read
  internal note `IntNotesRun2DimuonReference/ATL-COM-PHYS-2021-1094.pdf` and
  update `run2_dimuon_note.md` re journal relationship.
- B. 2109.00411 Run2 HF-muon R_AA + local note `Run2 HF muons nuclear
  modification factor internal notes.pdf`.
- C. 2003.03565 Run2 HF-muon v_n + local note `Run2 HF muons azimuthal
  anisotropies internal notes.pdf`.

HI field → `kb/physics/heavy_ion/`:
- D. 1802.04801 HI big picture (PRIMARY; weight open-HF refs).
- E. 1903.07709 Open HF production in HI (PRIMARY).
- F. local `Heavy-flavor production ... hot QCD matter.pdf` (PRIMARY).
- G. 2105.11656 RHIC open-HF review (SUPPORTIVE).
- H. 1901.01606 HF theory overview (SUPPORTIVE).

Detector → `kb/physics/detector/`:
- I. 2004.13447 Run2 muon trigger perf (PRIMARY).
- J. 2012.00578 Run2 muon reco — already summarized (`ATLAS_Run2_muon_reconstruction.md`); integrate only.
- K. WEB SEARCH Run3 ATLAS muon perf (reco+trigger) (PRIMARY if official source found).
- L. WEB SEARCH ATLAS inner detector/tracking (run2+run3) relevant to muon/HI; skip if nothing.

Background → `kb/physics/background/`:
- M. Gluon splitting / flavour excitation (indico 2308 link is anti-bot blocked;
  capture concept + future-read citeable source). Orchestrator-handled.

### Empirical learnings so far (for step 4)
- arXiv `curl https://arxiv.org/pdf/<id>` works (PDF). indico is behind an
  Anubis anti-bot PoW → curl returns HTML; cannot auto-download. → GUIDE/command
  should note: non-arXiv hosts may block curl; fall back to URL-only + concept
  capture, never fabricate.
</content>
