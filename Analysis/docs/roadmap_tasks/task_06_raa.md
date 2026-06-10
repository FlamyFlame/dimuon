# Task 06 — R_AA from corrected cross-sections

**Roadmap:** `analysis_roadmap_2026_06.md` chain [7]
**Depends on:** task 05. **Blocked input from user:** confirmed
luminosities (esp. PbPb25 — currently a placeholder copying 2023) and
5.36 TeV T_AA values (roadmap Q2.1–2.3). Proceed with placeholders,
clearly flagged. **Reviewer:** `/review-plot`.

## Objective

Produce R_AA(pair pT, pair η, centrality) from efficiency-corrected
PbPb (23+24+25 combined) and pp24 cross-sections.

## Inputs

- Corrected `histograms_real_pairs_pbpb_20YY_*` (task 05) and pp24
  equivalents; PbPb crossx always combined across years
  (memory: feedback_pbpb_crossx_combined).
- `Analysis/RAA_plotting.cxx` — exists but **check `base_dir`: it is
  hardcoded to a macOS path** (`/Users/yuhanguo/...`); port to
  `DATA_BASE=/usatlas/u/yuhanguo/usatlasdata/dimuon_data` or make it
  configurable before use.
- Normalization: `PbPbBaseClass.h` crossx factors (1/(Δcent · σ_PbPb ·
  T_AA · L)); `PPBaseClass.h` for pp. Update PbPb25 entries when the
  user provides values; until then the 2025 contribution carries the
  2023 placeholder and **every output plot must carry a
  "PRELIMINARY NORMALIZATION" label**.

## Steps

1. Fix/port `RAA_plotting.cxx` paths; ACLiC compile.
2. Verify the 3D hist names it reads
   (`h3d_{ss,op}_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt`)
   still exist in the corrected output files.
3. Produce R_AA vs pair pT (per centrality), vs centrality (per pT
   bin), vs pair η; SS and OS separately; log axes per plotting rules
   for log-binned variables.
4. `/review-plot`.
5. Update roadmap ledger; record placeholder-normalization caveats.

## Verification

- R_AA peripheral (50–80%) closer to unity than central; statistical
  errors propagate from both pp and PbPb.
- Cross-check one bin by hand: N_PbPb/(N_evt·T_AA·σ...) against the
  plotted value.
