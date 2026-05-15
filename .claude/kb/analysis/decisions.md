# Key Analysis Decisions

Log of technique choices and their rationale. Append new decisions as they are made.

## Differential histogram scaling

**Decision**: Use `h->Scale(NORM_FACTOR, "width")` for differential plots.
**Rationale**: ROOT divides each bin by its own width automatically. Works correctly for non-uniform binning. Never manually compute bin_width and divide.

## Pair eta binning

**Decision**: Pair eta binning uses `pair_eta_proj_ranges_coarse_incl_gap` from `CommonEffcyConfig.h`.
**Rationale**: Signal cuts changed from |pair_eta| < 2.2 to |q*eta| < 2.2 for each muon, so pair eta extends to 2.4. The `q_eta_proj_ranges_*` binnings are for single-muon trigger efficiency fitting only. Mixing them silently cuts off the [2.0, 2.4] pair eta bin.

## Plot output format

**Decision**: Save all plots as PNG only. No PDF unless explicitly requested.

## PbPb25 per-run preamp cut

**Decision**: Use per-run mu+7sigma Gaussian-fit ZDC preamp cut for PbPb25, instead of a single hard-coded scalar.
**Rationale**: Beam conditions and ZDC calibration shift significantly run-to-run in PbPb25, making a single global cut ineffective.

## Overlay barcode handling

**Decision**: Use first-writer-wins (`emplace`) for barcode cache in overlay, not last-writer-wins.
**Rationale**: Overlay has ~570 duplicate barcodes/event in bc 1-510 (Pythia and HIJING generators both start from bc=1). Last-writer-wins returned the HIJING particle instead of the Pythia particle, causing ~68% s_light misclassification. Fixed Apr 2026.
