# Analysis Roadmap Buildout (Tracking)

## Objective

Strategic task (2026-06-10):
1. Read the Run 2 dimuon internal note (`IntNotesRun2DimuonReference/`,
   ATL-COM-PHYS-2021-1094, 51 pp) and write a comprehensive summary into
   the knowledge base (`.claude/kb/`), pointing back to the note for
   reference.
2. Read the entire repo and produce a dedicated **analysis roadmap doc**
   (dated 2026-06-10), one section per question:
   1) Which IntNote (ANA-HION-2023-07-INT1) sections can already be
      written (results ready, or methodology clear → placeholders)?
   2) What info is missing and must come from the user (lumi, prescales,
      MC info, ...)? Never assume — list guesses as guesses.
   3) Which analysis steps are we capable of but haven't performed /
      finished (preliminary-only OR precise-but-not-propagated)?
   4) Which steps cannot yet be performed + lacking prerequisites
      (mostly fullsim / fullsim-overlay sample availability); dummies
      will be used to keep the chain moving.
   5) Full roadmap of the analysis chain (with dummies), plus per-step
      .md instruction files for future Claude Code agents.

## Plan / Checklist

- [x] Read all 4 active tracking docs (protocol)
- [x] Survey IntNotes skeleton + Run 2 note file structure
- [x] Read Run 2 note (LaTeX sources — no poppler for PDF on this node)
- [x] Write kb summary file `.claude/kb/analysis/run2_dimuon_note.md`, registered in kb index
- [x] Read Analysis/README.md + Analysis/docs/*.md (pipelines, samples)
- [x] Read IntNotes skeleton main tex (bare bones: Intro/Datasets/Analysis/Results; tex/datasets.tex has 3 JIRA refs ATLHI-658/596/666)
- [x] Survey datasets on disk (data + MC) and existing outputs/plots
- [x] Write roadmap doc `Analysis/docs/tracking/analysis_roadmap_2026_06.md`
      (sections Q1-Q5, dated 2026-06-10)
- [x] Write per-step instruction .md files (`Analysis/docs/roadmap_tasks/`
      task_01..task_09)
- [x] Register roadmap in Active Tracking Docs (CLAUDE.md)

## Key survey findings (2026-06-10)

- **pp24 stale:** raw May 2026 skim (12 parts) on disk, but
  `muon_pairs_pp_2024_*` ntuples date Feb 22–Mar 31, 2026 (old skim).
  pp24 NTP + crossx + P2/P3 all need rerun.
- **Luminosity values in code (provenance UNVERIFIED — ask user):**
  `PPBaseClass.h`: pp17 2mu4 256.8 pb⁻¹; pp24 mu4_mu4noL1 113.999 pb⁻¹;
  pp24 2mu4 410.815 pb⁻¹. `PbPbBaseClass.h`: PbPb23 1.3896 nb⁻¹,
  PbPb24 1.5411 nb⁻¹, **PbPb25 = TODO placeholder (uses 2023 values)**;
  σ_PbPb = 7.8 b; T_AA per ctr bin hardcoded (26.1428, 20.3241, 14.0502,
  8.5074, 3.7733, 0.6716).
- **Fullsim full samples ABSENT:** `pythia_fullsim_full_sample/` is empty;
  overlay has only test samples (r17618 60k evts + r17662/r17663 ~1k).
  Powheg fullsim exists for **pp17 only** (Run 2 conditions); pp24
  (ATLHI-658) not produced.
- **No template-fitting code** (Δp/p or otherwise) exists in the repo.
- **Crossx pipeline applies NO trig/reco efficiency corrections** (P1
  uses trigger only for selection). P2 ε^nc fits exist for PbPb 23/24/25
  but are not propagated into pair weights.

## Progress Log (cont.)

### 2026-06-10: kb summary written
- `.claude/kb/analysis/run2_dimuon_note.md`: full summary of
  ATL-COM-PHYS-2021-1094 (physics goal, datasets/triggers/GRL, MC,
  selections, corrections, Δp/p templates, Δφ fits, systematics tables,
  note structure, Run2→Run3 differences). Registered in kb index.

## Progress Log

### 2026-06-10: Setup
- Read all 4 active tracking docs. Key statuses:
  - PbPb P1/P2/P3 all rerun with May 2026 skim (all 3 years) as of 06-09.
  - dR corrections: single-muon & pair-level removed (D9, biased);
    cross-term kept as MC reference. PP24 P2/P3 blocked on re-skim.
  - Overlay reco efficiency: 3 bugs found; Bug #3 (HIJING truth in denom)
    fixed; Bug #1 (index mapping) + Bug #2 (AOD truth-match dR fallback)
    have written specs, NOT yet implemented.
  - r-tag comparison: r17662 (StandardSignalOnlyTruth) recommended for
    production overlay.
- IntNotes skeleton = `dimuon_codes/IntNotes/` (ANA-HION-2023-07-INT1),
  essentially empty (only MC request JIRA list found in main tex).
- Run 2 note tex sections: Introduction, Data_and_selections, MuonAnalysis,
  MediumMuons, DpopAna, MuonSources, Templates, Analysis, Corrections,
  DrellYan, McClosure, Results, Systematics, SystematicsPP, StatErrs,
  DphiFits, DphiFitsDeta0p9, Conclusion + 5 appendices.

### 2026-06-10: Roadmap + task files delivered — TASK COMPLETE

Deliverables:
1. kb: `.claude/kb/analysis/run2_dimuon_note.md` (+ index entry)
2. Roadmap: `Analysis/docs/tracking/analysis_roadmap_2026_06.md`
   (Q1 IntNote readiness table; Q2 missing-input list with unverified
   in-code values flagged; Q3 capable-but-pending; Q4 blocked + dummy
   strategies; Q5 chain diagram + task table; status ledger)
3. Task instructions: `Analysis/docs/roadmap_tasks/task_01..09_*.md`
4. Roadmap registered in CLAUDE.md Active Tracking Docs.

Immediately startable in parallel: task_01 (pp24 reprocess), task_03
(overlay reco-eff bugfix), task_07 (Δp/p framework), task_09 (note
sections). This buildout doc is complete; ongoing status lives in the
roadmap doc's ledger.

## Latest Stage

(complete — see final summary above)
