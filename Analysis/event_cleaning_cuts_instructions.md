# Instructions: Successive Event Cleaning Cuts for Pb+Pb Event Selection

## Context

We are performing event selection (cleaning) for 2024 Pb+Pb collision data at √s_NN = 5.36 TeV (ATLAS). The goal is to remove pathological events (pileup, out-of-time backgrounds, non-hadronic interactions) using a sequence of four cuts applied one after another. At each stage we produce a set of **six 2D diagnostic histograms** to visualize the effect of each cut.

The code is ROOT C++ using TH2 histograms and TCanvas for plotting.

---

## Overview of the Four Successive Cuts

The cuts are applied **cumulatively in this order**:

1. **Track quality fraction cut** — ratio of HITight tracks to all tracks (`n_trk_HITight / n_trk_all`). Events with an anomalously low fraction are rejected.
2. **ZDC presample amplitude sum cut** — the sum of ZDC C-side presample amplitudes. Events with large presample signals (out-of-time pileup) are rejected.
3. **ZDC–FCal "banana" cut** — a 2D cut on the correlation between ZDC energy (`E_ZDC`) and FCal ΣE_T. Events falling outside the expected correlation band are rejected.
4. **ZDC timing cut** — a 2D cut on ZDC A-side time vs C-side time. Events with out-of-time ZDC signals are rejected.

---

## The Six Diagnostic 2D Histograms

For every stage (all events, fail cut N, pass cut N, pass cuts 1..N fail cut N+1, etc.), we produce the **same six plots**:

| # | X-axis | Y-axis | Description |
|---|--------|--------|-------------|
| 1 | `n_trk_all` | `n_trk_HITight / n_trk_all` | Track quality fraction vs total track multiplicity |
| 2 | `FCal Σ_{η<0} E_T` [TeV] | `FCal Σ_{η>0} E_T` [TeV] | FCal E_T A-side vs C-side |
| 3 | `FCal ΣE_T` [TeV] | `E_ZDC` [TeV] | ZDC energy vs FCal total E_T |
| 4 | `C Time` [ns] | `A Time` [ns] | ZDC timing: A-side vs C-side |
| 5 | `FCal ΣE_T` [TeV] | `n_trk_HITight` | HITight track multiplicity vs FCal E_T |
| 6 | `Primary vtx N_trk` | `σ_z²` [mm²] | Primary vertex z-width squared vs vertex track count |

All histograms use a **log z-axis color scale** (prescale-weighted events). Use the `COLZ` draw option. Each plot should include an `ATLAS Internal` label, `Pb+Pb √s_NN = 5.36 TeV`, and `2024 PC+CC`.

---

## Cut Definitions

Below are the cut definitions. You will need to adapt variable names to match whatever branch/variable names exist in the existing codebase. Search the codebase for these variables to find the correct names.

### Cut 1: Track quality fraction

```
pass_trackfrac = (n_trk_HITight / n_trk_all) > THRESHOLD
```

- The threshold is a function of `n_trk_all`. From the plots, it appears to be a constant or simple function. Look in the existing code for a track fraction cut — it likely already exists. If not, a typical value is around 0.2–0.3 depending on multiplicity, or a centrality-dependent parameterization.
- **Action**: Search the codebase for any existing track quality fraction cut or variable like `trackFraction`, `trkQualFrac`, `HITight`, `ntrk_ratio`, or similar.

### Cut 2: ZDC presample amplitude sum

```
pass_zdcpresample = (ZDC_C_presample_amp_sum < THRESHOLD)
```

- This cuts on the ZDC C-side presample amplitude sum. Events with large presample values have out-of-time pileup.
- **Action**: Search for variables like `ZDC_presample`, `presampleAmp`, `ZdcModulePresampleAmp`, or similar. The threshold should be tuned from the distribution — typically a value that removes the tail above the main peak.

### Cut 3: ZDC–FCal "banana" cut

```
pass_banana = E_ZDC < f(FCal_SumET)
```

- This is a 2D cut on the ZDC energy vs FCal ΣE_T plane. The boundary follows the upper edge of the "banana"-shaped correlation. Events above this boundary (high ZDC energy for their FCal E_T) are electromagnetic dissociation or UPC events.
- **Action**: Search for any existing banana cut, ZDC–FCal correlation cut, or parameterization like `zdcFcalCut`, `banana`, `E_ZDC_cut`.

### Cut 4: ZDC timing

```
pass_zdctime = (|ZDC_A_time| < T_CUT) && (|ZDC_C_time| < T_CUT)
```

- Or possibly a 2D circular/elliptical cut around (0, 0) in the ZDC A-time vs C-time plane.
- **Action**: Search for variables like `ZDC_time`, `zdcTimeA`, `zdcTimeC`, `timeCut`.

---

## Plot Production Procedure

### Step 0: Setup

Find and read the existing histogram filling code that was modified in the relevant commit. Identify:
- The input TTree/TChain and branch names for all six diagnostic observables
- The existing event loop structure
- How prescale weights are applied
- How ATLAS labels are drawn (look for `ATLASLabel`, `myText`, or similar helper functions)
- The output directory for plots

### Step 1: Define the cut booleans

In the event loop, compute four boolean flags per event:

```cpp
bool pass_trackfrac    = /* track quality fraction cut */;
bool pass_zdcpresample = /* ZDC presample cut */;
bool pass_banana       = /* ZDC-FCal banana cut */;
bool pass_zdctime      = /* ZDC timing cut */;
```

### Step 2: Create histogram sets

Create **nine sets** of the six 2D histograms (each set = 6 TH2D/TH2F). Name them with a suffix indicating the selection stage:

| Set | Title / suffix | Selection |
|-----|---------------|-----------|
| 0 | `_allEvents` | No cuts (all events) |
| 1 | `_fail_trackfrac` | Fail cut 1 |
| 2 | `_pass_trackfrac` | Pass cut 1 |
| 3 | `_pass_trackfrac_fail_zdcpresample` | Pass cut 1, fail cut 2 |
| 4 | `_pass_trackfrac_zdcpresample` | Pass cuts 1+2 |
| 5 | `_pass_trackfrac_zdcpresample_fail_banana` | Pass cuts 1+2, fail cut 3 |
| 6 | `_pass_trackfrac_zdcpresample_banana` | Pass cuts 1+2+3 |
| 7 | `_pass_trackfrac_zdcpresample_banana_fail_zdctime` | Pass cuts 1+2+3, fail cut 4 |
| 8 | `_pass_trackfrac_zdcpresample_banana_zdctime` | Pass all 4 cuts |

That gives 9 × 6 = 54 histograms total.

### Step 3: Fill histograms in the event loop

For each event, compute the four cut booleans, then fill the appropriate histogram sets:

```cpp
// Always fill set 0
fillHistSet(hset_allEvents, ...);

// Cut 1 split
if (!pass_trackfrac) {
    fillHistSet(hset_fail_trackfrac, ...);
} else {
    fillHistSet(hset_pass_trackfrac, ...);

    // Cut 2 split (only among events passing cut 1)
    if (!pass_zdcpresample) {
        fillHistSet(hset_pass_trackfrac_fail_zdcpresample, ...);
    } else {
        fillHistSet(hset_pass_trackfrac_zdcpresample, ...);

        // Cut 3 split (only among events passing cuts 1+2)
        if (!pass_banana) {
            fillHistSet(hset_pass_trackfrac_zdcpresample_fail_banana, ...);
        } else {
            fillHistSet(hset_pass_trackfrac_zdcpresample_banana, ...);

            // Cut 4 split (only among events passing cuts 1+2+3)
            if (!pass_zdctime) {
                fillHistSet(hset_pass_trackfrac_zdcpresample_banana_fail_zdctime, ...);
            } else {
                fillHistSet(hset_pass_trackfrac_zdcpresample_banana_zdctime, ...);
            }
        }
    }
}
```

Use prescale weights in the `Fill()` calls.

### Step 4: Draw and save canvases

For each of the 9 histogram sets, create a **single TCanvas divided into a 3×2 grid** (or six separate canvases — match the existing code style). Each sub-pad shows one of the six diagnostic plots.

Drawing settings for each 2D histogram:
- `SetStats(0)` — no stats box
- Draw with `"COLZ"` option
- `gPad->SetLogz(1)` — log color scale
- For the vertex σ_z² plot: also `gPad->SetLogx(1)` and `gPad->SetLogy(1)`
- Add ATLAS label and dataset info text on each pad

Title each canvas with the selection stage, e.g.:
- "All events"
- "Fail: track fraction"
- "Pass: track fraction"
- "Pass: track fraction / Fail: ZDC presample"
- "Pass: track fraction, ZDC presample"
- "Pass: track fraction, ZDC presample / Fail: banana"
- "Pass: track fraction, ZDC presample, banana"
- "Pass: track fraction, ZDC presample, banana / Fail: ZDC timing"
- "Pass: track fraction, ZDC presample, banana, ZDC timing"

Save each canvas as both `.pdf` and `.png`.

---

## Histogram Binning (approximate, adjust as needed)

| Observable | Typical range | Suggested bins |
|-----------|---------------|----------------|
| `n_trk_all` | 0 – 12000 | 120 bins |
| `n_trk_HITight / n_trk_all` | 0 – 1 | 100 bins |
| `FCal Σ_{η<0} E_T` | 0 – 3.5 TeV | 70 bins |
| `FCal Σ_{η>0} E_T` | 0 – 3.5 TeV | 70 bins |
| `FCal ΣE_T` | 0 – 7 TeV | 70 bins |
| `E_ZDC` | 0 – 400 TeV | 80 bins |
| `C Time` | −12 – 12 ns | 120 bins |
| `A Time` | −12 – 12 ns | 120 bins |
| `n_trk_HITight` | 0 – 4000 | 80 bins |
| `Primary vtx N_trk` | 1 – 2×10⁴ | 100 log-spaced bins |
| `σ_z²` | 10⁻⁵ – 1 mm² | 100 log-spaced bins |

---

## Implementation Notes

- **Do not hardcode cut values** — define them as constants or configurable parameters at the top of the file so they can be tuned.
- **Search the existing codebase first** for any existing implementations of these cuts, variable names, helper functions, and plotting utilities before writing new code.
- **Match the existing code style** — follow the same conventions for histogram creation, canvas styling, output file naming, and directory structure as the rest of the codebase.
- **Prescale weights** — make sure `Fill(x, y, weight)` uses the prescale weight for data histograms.
- **Canvas size** — use a wide canvas for the 3×2 grid, e.g., `TCanvas("c", "c", 2400, 1600)` divided into `Divide(3,2)`.
- The plots in the reference presentation use the **viridis-like** ROOT color palette. Set it with `gStyle->SetPalette(kBird)` or `gStyle->SetPalette(55)`.

---

## Summary of Deliverables

1. Modified C++ source file(s) with histogram booking, filling logic (nested if/else for successive cuts), and plotting.
2. Nine output canvases (PDF + PNG), each containing the six diagnostic 2D histograms for one selection stage.
3. Optionally, a summary printout of the number of events (prescale-weighted) passing and failing each successive cut.
