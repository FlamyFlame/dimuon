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
2. [DONE] Bulk build (12 sources, 2 batches of subagents). → batches 1+2.
3. [DONE] Structure + index + knowledge graph + `/kb-review` FULL (audit + blind
   use-case tests). Re-audit iter 2 = PASS.
4. [DONE] Self-improved the procedure from real build experience (see step-4 log).

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

- **2026-06-15 — Step 2 batch 2 DONE** (6 subagents, no-git). Created:
  `physics/heavy_ion/hf_hot_qcd_matter.md` (Averbeck PPNP70 2013, no arXiv),
  `physics/heavy_ion/rhic_open_hf_review.md` (2105.11656, SUPP),
  `physics/heavy_ion/hf_theory_overview.md` (Gossiaux 1901.01606, SUPP),
  `physics/detector/atlas_run2_muon_trigger.md` (2004.13447, +PDF),
  `physics/detector/atlas_run3_muon_performance.md` (2401.06630 trigger paper +
  Run3 reco slides, +2 PDFs; **no full Run-3 reco paper / HI muon perf exists**),
  `physics/detector/atlas_inner_detector_tracking.md` (1704.07983 dense-env
  tracking — justifies our pair ε_reco ΔR axis, SUPP). No-git instruction worked
  (no concurrency collisions this batch).
- **2026-06-15 — Step 3 DONE (integration; orchestrator).** Read all 12 new docs.
  Decision: subagents already deduped by deferring shared concepts → **designated
  canonical hubs** (`open_hf_production` = HF-physics; `ATLAS_Run2_muon_reconstruction`
  = muon types/WP/tag-and-probe) rather than splitting them. Created the one
  genuinely cross-paper concept doc `concepts/muon_source_template_fits.md`
  (Δp/p + d0 fits, preserving per-paper use differences). Wrote
  `physics/background/gluon_splitting_flavour_excitation.md` (indico talk anti-bot
  blocked → recorded URL + concept + citeable anchor arXiv:1812.09283 as top
  future-read). Rebuilt `index.md` (Meta/Analysis/Concepts/HI/Detector/Background/
  Data/Procedures/Gotchas). Deduped+prioritized `FUTURE_READ.md` (5 HIGH + 6 MED).
  Added bidirectional hub/back-links. Self-synced kb-review scenarios (+Intro
  motivation, +background/template). Committed (30a907d sources, 8485bac structure).
- **2026-06-15 — Step 3 review DONE (the loop).** Ran `/kb-review` FULL: (A)
  structural-audit subagent → FAIL with 4 WARNINGs (all blind-builder artifacts:
  one-way hub sink, uncrosslinked duplicate-topic pair + stale claim, one dead
  `[[...]]`, un-trimmed concept duplication) + INFOs. (B) **3 blind use-case test
  agents** (write Intro motivation / explain background+template / implement Δp/p
  purity) — **all PASS**: each discovered the KB via CLAUDE.md→index, found the
  right docs efficiently (process notes praised the hub designation as "highest
  leverage"), produced grounded output; the Δp/p agent even found a skim-branch
  gating issue via KB pointers. Fixed all 4 WARNINGs + INFOs; **re-audit iter 2 =
  PASS** (0 CRITICAL/0 WARNING; lone INFO = cosmetic `[[...]]` syntax mix).
- **2026-06-15 — Step 4 DONE (self-improvement from real experience).** Codified
  build learnings into the tooling: (1) **subagents NEVER run git** — orchestrator
  does all git, `git add -f` for gitignored PDFs (fixes the batch-1 index.lock
  collision); (2) **non-arXiv hosts (indico/CDS) block curl** (Anubis anti-bot) →
  URL-only fallback, never fabricate; (3) **prefer DESIGNATING an existing
  comprehensive doc as a hub** over creating many concept docs (+ add hub→spoke
  back-links, since hubs start as one-way sinks); (4) the **deferral-pointer
  pattern** resolves GUIDE §3-self-sufficiency vs §8-dedup; (5) standardize
  `[[basename]]` links; (6) `gs` may segfault per-page. Edits: `kb-builder.md`,
  `kb-build.md` (Phase 2), `KB_BUILDING_GUIDE.md` (§3, §8), `gotchas/reading_pdfs.md`.

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

**All four steps complete.** Ongoing (not part of this effort — via `/kb-build`
ADD mode as needed):
- Summarize queued `FUTURE_READ.md` items — top priority arXiv:2206.12594 (the
  actual per-pair trig+reco eff method) and arXiv:1812.09283 (g→bb̄ at small
  opening angle, our regime).
- Add docs as new analysis steps (unfolding, template fitting, systematics) appear;
  `/kb-review` FULL scenarios are stubbed for self-sync.

## Final KB inventory (2026-06-15)
24 docs. analysis (5: overview, decisions, dimuon note + backtoback paper, HF-muon
R_AA, HF-muon v_n) · concepts (1: muon_source_template_fits) · physics/heavy_ion
(5: open_hf_production [HUB], hi_big_picture, hf_hot_qcd_matter, rhic_open_hf_review,
hf_theory_overview) · physics/detector (4: ATLAS_Run2_muon_reconstruction [HUB],
run2 trigger, run3 perf, inner_detector_tracking) · physics/background (1:
gluon_splitting_flavour_excitation) · data (2) · procedures (1) · gotchas (1) ·
Meta (GUIDE, FUTURE_READ). 12 PDFs committed (`git add -f`). Commits: 30a907d
(sources), 8485bac (structure/tooling), f44111a (review fixes), + step-4 tooling.

## Latest Stage

— (complete; removed from Active Tracking Docs 2026-06-15.) The KB-building system
(`/kb-build`, `/kb-review`, `KB_BUILDING_GUIDE`, `kb-builder` agent) and the initial
bulk build (12 sources) are done and review-validated (re-audit PASS; blind
use-case tests PASS). Further sources are added incrementally via `/kb-build` +
`/kb-review`.
</content>
