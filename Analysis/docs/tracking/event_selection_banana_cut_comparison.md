# PbPb ZDC–FCal Banana Cut: Two Procedures Compared

> **How to use:** Upload each image file listed under its slide, then paste this document into Claude.  
> All image paths are on the BNL SDCC cluster at  
> `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/`

---

## Slide 1 — Nominal procedure: two-band Gaussian fitting (PbPb 2024)

**Image:** `event_selection/pbpb_2024/cut1_ZDC_FCal_2graph_support_pbpb_2024.png`

**Caption:**  
The nominal ZDC–FCal banana cut is derived through a two-band Gaussian fitting procedure applied independently to each 0.1 TeV FCal E_T slice.

**Main band (red):** In each slice the ZDC total-energy projection is fitted with a two-pass Gaussian. The first pass uses a ±3 RMS window around the distribution peak; the second refines to ±3σ₁ around the first-pass mean. The banana cut is placed at μ + 5σ (upper edge of the main hadronic interaction band).

**Pileup / out-of-time band (blue/yellow):** For FCal E_T > 2 TeV, a second Gaussian is fitted inside a per-year search window ([220, 350] TeV for 2023; [180, 300] TeV for 2024; [170, 320] TeV for 2025). An initial quadratic guess constrains the pass-1 centre. The final veto threshold is μ_bg − 3σ_bg (lower edge of the pileup band). Where the pileup fit does not converge, the cut is extended via wider-bin retry (0.2 or 0.4 TeV slices) and TGraph extrapolation from converged slices.

**Final cut:** max(main μ + 5σ, pileup μ − 3σ), so events in the out-of-time pileup band are rejected. Yellow markers show the pileup μ ± σ; orange dashes show the quadratic initial guess.

---

## Slide 2 — Nominal two-band cut applied: PbPb 2023

**Image:** `event_selection/pbpb_2023/ZDC_E_tot_vs_FCal_Et_AC_with_cut_pbpb_2023.png`

**Caption:**  
ZDC total energy vs FCal E_T for PbPb 2023 data (HLT single-muon trigger). The red curve is the final nominal banana cut: max(main μ + 5σ, pileup μ − 3σ). Events above the red curve are rejected. The two visible bands — the main hadronic band curving downward with centrality and the flatter pileup band at higher ZDC — are both incorporated into the cut.

---

## Slide 3 — Nominal two-band cut applied: PbPb 2024

**Image:** `event_selection/pbpb_2024/ZDC_E_tot_vs_FCal_Et_AC_with_cut_pbpb_2024.png`

**Caption:**  
ZDC total energy vs FCal E_T for PbPb 2024 data. Same two-band procedure as 2023, with the pileup search window [180, 300] TeV. The pileup band is clearly separated from the main band at high FCal E_T, and the cut follows its lower edge where it rises above the main band's upper edge.

---

## Slide 4 — Nominal two-band cut applied: PbPb 2025

**Image:** `event_selection/pbpb_2025/ZDC_E_tot_vs_FCal_Et_AC_with_cut_pbpb_2025.png`

**Caption:**  
ZDC total energy vs FCal E_T for PbPb 2025 data. Pileup search window [170, 320] TeV. With the largest dataset (~33 M trigger-passing events), the pileup band is well-populated and the cut closely tracks its lower edge across the full FCal range.

---

## Slide 5 — FCal E_T shape comparison after nominal two-band cuts

**Image:** `event_selection/fcal_comparison_pbpb_nominal.png`

**Caption:**  
FCal E_T distributions after applying all five nominal event-selection cuts (ZDC banana, ZDC time, ZDC preamp, nTrk fraction, nTrk–FCal band), normalized to unit area. **Left:** PbPb 2024 (red, open) overlaid on 2023 (black), with the 24/23 ratio below. **Right:** PbPb 2025 (blue, open) overlaid on 2023, with the 25/23 ratio. Each year uses its own per-year cut file derived from the two-band procedure.

---

## Slide 6 — Alternative procedure: quadratic banana cut applied: PbPb 2023

**Image:** `event_selection_alternative_banana/pbpb_2023/ZDC_E_tot_vs_FCal_Et_AC_with_cut_pbpb_2023.png`

**Caption:**  
**Alternative banana cut procedure:** The pileup-band Gaussian fit is replaced by a simple second-order polynomial. Only the main hadronic band is fitted (same two-pass Gaussian, cut at μ + 5σ). For FCal E_T < 2 TeV the cut follows the main-band TGraph directly. For FCal E_T ≥ 2 TeV, the cut is a quadratic a·x² + b·x + c whose three coefficients are fully determined by three points: (1) the main-band cut evaluated at x = 2 TeV (derived from data), and (2, 3) two hard-coded per-year calibration points chosen in the gap between the main band and pileup band:

| Year | Point 2 | Point 3 |
|------|---------|---------|
| 2023 | (3.4 TeV, 225 TeV) | (4.8 TeV, 165 TeV) |
| 2024 | (3.4 TeV, 202 TeV) | (4.8 TeV, 153 TeV) |
| 2025 | (3.4 TeV, 200 TeV) | (4.8 TeV, 152 TeV) |

This avoids the fit instability at FCal E_T < 3 TeV where the main-band tail contaminates the pileup search window. Shown here for PbPb 2023.

---

## Slide 7 — Alternative quadratic cut applied: PbPb 2024

**Image:** `event_selection_alternative_banana/pbpb_2024/ZDC_E_tot_vs_FCal_Et_AC_with_cut_pbpb_2024.png`

**Caption:**  
ZDC total energy vs FCal E_T for PbPb 2024 with the alternative quadratic banana cut. The red curve smoothly transitions from the fitted main-band TGraph below 2 TeV to the quadratic above, passing through (2.0, 251.3), (3.4, 202), and (4.8, 153) TeV.

---

## Slide 8 — Alternative quadratic cut applied: PbPb 2025

**Image:** `event_selection_alternative_banana/pbpb_2025/ZDC_E_tot_vs_FCal_Et_AC_with_cut_pbpb_2025.png`

**Caption:**  
ZDC total energy vs FCal E_T for PbPb 2025 with the alternative quadratic banana cut. The quadratic passes through (2.0, 242.0), (3.4, 200), and (4.8, 152) TeV.

---

## Slide 9 — FCal E_T shape comparison after alternative quadratic cuts

**Image:** `event_selection_alternative_banana/fcal_comparison_pbpb_alt_banana.png`

**Caption:**  
FCal E_T distributions after applying all five event-selection cuts with the alternative quadratic banana cut, normalized to unit area. Layout is identical to Slide 5. Comparing the 24/23 and 25/23 ratios between the two procedures (Slide 5 vs this slide) quantifies the systematic effect of the banana-cut choice on the centrality-weighted event sample.
