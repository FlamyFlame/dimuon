# Task 01 — pp24 NTuple reprocessing + crossx on May 2026 skim

**Roadmap:** `docs/tracking/analysis_roadmap_2026_06.md` §Q3a / chain [2],[6]
**Can start:** immediately. **Reviewer:** `/review-analysis-code` for any
code change; `/review-plot` for the crossx plots.

## Objective

Bring pp24 to the same skim generation as PbPb: reprocess the May 2026
pp24 skim through NTuple processing, hadd, RDF crossx hist filling, and
crossx plotting.

## Why

`muon_pairs_pp_2024_*` files date Feb 22–Mar 31 2026 (old skim), while
raw `data_pp24_part1..12.root` are the May 2026 skim (see
`~/usatlasdata/dimuon_data/data-merging-record.txt`; parts are chunks of
2 grid tasks). The new skim also carries the passmu4noL1 fix and
mu6/mu8-disable (`DimuonDataAlgCoreT.c`), required by task 02.

## Inputs

- Raw: `~/usatlasdata/dimuon_data/pp_2024/data_pp24_part{1..12}.root`
- Code: `NTupleProcessingCode/DataAnalysisClasses.h` (`PPAnalysis`),
  `DimuonDataAlgCoreT.{h,c}`, `PPExtras.{h,c}`
- Condor: `NTupleProcessingCode/run_pp_24.sub` — **check queue count
  matches 12 parts** before submitting (old skim may have had fewer).

## Steps

1. Pre-checks: verify all 12 raw parts open and have entries matching
   the merging record; confirm `use_mu6_for_trg_eff=false` and the
   passmu4noL1 legacy-derivation fix are present in
   `DimuonDataAlgCoreT.c`. Also verify (for the PbPb side) that the
   current `muon_pairs_pbpb_2023_*` hadd includes part1 (24.9M entries).
2. Back up existing `muon_pairs_pp_2024_*single_mu4*` and `*2mu4*`
   combined files with a dated suffix (do not delete).
3. Submit Condor NTP jobs (per user rules: submit without asking,
   monitor until done). pp24 is trigger_mode=3 (2mu4) for nominal
   crossx — confirm against `DatasetTriggerMap.h` (`{24,"pp"} → 2mu4`).
4. Validate per-part outputs (non-zombie, both sign trees non-empty),
   hadd, validate combined.
5. RDF crossx hist filling (see `RDFBasedHistFilling/test_crossx_pp24.sh`
   — pp24 has no production pipeline; consider promoting this script to
   `pipelines/pipeline_pp24_crossx.sh` with the same validation pattern
   as the PbPb pipelines).
6. Run `plotting_codes/single_b_analysis/plot_single_b_crossx_pp.cxx`.
7. Update `data_analysis.md` doc + record in roadmap status ledger.

## Verification

- Combined sign1/sign2 entry counts comparable to (likely larger than)
  the old-skim files (old combined `_res_cut_v2` was 479 MB).
- Crossx plots regenerate; normalization uses
  `PPBaseClass::GetCrossxFactor(24,"2mu4")` (=1/410.815 pb⁻¹ —
  provenance unconfirmed, flag on plots as preliminary).
- m1/m2.passmu4noL1 branch non-trivial (fix validation, needed by task 02).
