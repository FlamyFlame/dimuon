# Slide Deck Instructions: PbPb Event Selection
**Topic:** Event selection for the PbPb dimuon analysis (Pb+Pb 2024 data, √s_NN = 5.36 TeV)

---

## Background & Context

These slides document the event-level quality cuts applied to raw PbPb collision data before the dimuon physics analysis. The data comes from the ATLAS detector at CERN. Events are first required to fire the HLT_mu4_L1MU3V single-muon trigger. Then five sequential cuts are applied, each requiring all prior cuts to already pass.

**Key detector components referenced in plots:**
- **FCal ET**: Forward Calorimeter transverse energy, measured separately on side A (η > 0) and side C (η < 0). The sum FCal_ET^{A+C} is a proxy for collision centrality — more central (head-on) collisions produce more FCal energy (up to ~5 TeV total).
- **ZDC**: Zero Degree Calorimeter, positioned at ±140 m from the interaction point. Measures forward neutrons from nuclear breakup. Used to tag genuine hadronic PbPb collisions and remove backgrounds.
- **nTrk HItight**: Number of tracks passing tight HI (heavy-ion) quality criteria. Should correlate strongly with FCal ET for genuine collisions.

**Centrality:** "Top 80% centrality" (centrality 0–79, where 0 = most central) is the analysis acceptance. Survival rates quoted below are always for this centrality range.

**The 5 cuts in order:**

| Cut | Variable | Method | Survival (top 80%) |
|-----|----------|--------|-------------------|
| 1 | ZDC total energy vs FCal ET | Upper boundary μ+5σ per FCal slice (TGraph) | 99.64% |
| 2 | ZDC time, both sides | Hard box: \|t_A\| < 1.8 ns AND \|t_C\| < 1.8 ns | 99.99% |
| 3 | ZDC pre-sample amplitude (A and C) | EMG fit, upper cut at μ+6σ per side | 94.76% |
| 4 | nTrk HItight fraction vs total nTrk | Lower bound μ−5σ per nTrk slice (TGraph) | 99.96% |
| 5 | nTrk HItight vs FCal ET | Band μ±5σ per FCal slice (TGraph pair) | 99.92% |

**What the 5-panel diagnostic plots show:** Each 5-panel canvas contains these subplots:
1. (top-left) ZDC total energy vs FCal ET sum — shows the ZDC banana structure
2. (top-center) ZDC time side A vs side C — should cluster at (0,0) for good events
3. (top-right) FCal ET side A vs side C — should be symmetric for genuine collisions
4. (bottom-left) nTrk HItight fraction vs total nTrk — should be near constant ~0.8–0.9
5. (bottom-center) nTrk HItight vs FCal ET — tight linear correlation for good events

These are shown at each stage (no cuts, after each cut for pass and fail populations) so you can see what the passing and failing events look like across all variables simultaneously.

**Local plot directory:** `/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/plots/single_b_analysis/event_selection/`
**Alternative plot directory:** `/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/plots/single_b_analysis/event_selection_alternative/`

All plots are for Pb+Pb 2024 data (suffix `_pbpb_2024.png`).

---

## Slide-by-Slide Instructions

---

### SLIDE 1 — Overview: Sequential Event Selection

**Title:** Sequential Event Selection for Pb+Pb 2024 Data

**No figure.** Create a clean text slide with the following content:

- Opening sentence: "Five sequential quality cuts are applied to events passing the HLT_mu4_L1MU3V single-muon trigger. Each cut requires all prior cuts to pass."
- Present the 5 cuts as a numbered vertical flow (arrow connecting each step), with a one-line description per cut:
  1. **ZDC–FCal banana** — remove events with anomalously high ZDC energy for their FCal ET (non-hadronic background)
  2. **ZDC time box** — require both ZDC sides to fire within ±1.8 ns of collision time (remove out-of-time pile-up)
  3. **ZDC pre-sample amplitude** — remove events with saturated/noisy ZDC pre-amplifiers (both sides < threshold)
  4. **nTrk HItight fraction** — require the fraction of tight-quality tracks to be consistent with genuine collisions
  5. **nTrk HItight vs FCal ET band** — require track multiplicity to be consistent with the measured FCal energy
- Footer note: "Survival rates shown are for top 80% centrality events (centrality 0–79), the analysis acceptance."

---

### SLIDE 2 — Distributions Before Any Cuts

**Title:** Event Distributions Before Any Quality Cuts

**Figure:** `event_sel_nocuts_5panel_pbpb_2024.png`

**Caption:** "Two-dimensional distributions of key event-level variables for all events passing the HLT_mu4_L1MU3V trigger, before any quality cuts. The five panels show: ZDC energy vs FCal ET (banana structure from nuclear breakup neutrons), ZDC time correlation (should cluster near origin), FCal ET side A vs C (symmetric for genuine collisions), nTrk HItight fraction vs total nTrk, and nTrk HItight vs FCal ET. All genuine Pb+Pb collision events are expected to populate well-defined correlation bands in each panel."

---

### SLIDE 3 — Cut 1: ZDC–FCal Banana

**Title:** Cut 1 — ZDC–FCal Banana

**Figure:** `event_sel_cut1_ZDC_FCal_banana_standalone_pbpb_2024.png`

**Caption:** "ZDC total energy vs FCal ET for events before cut 1. The red curve is the derived upper boundary: for each FCal ET slice, the ZDC distribution is fit with a Gaussian and the cut is placed at μ+5σ. Events above this boundary (high ZDC energy for their FCal ET) are inconsistent with genuine hadronic collisions and are removed. Survival rate for top 80% centrality: **99.64%**."

---

### SLIDE 4 — Cut 1: Events Failing the ZDC–FCal Banana

**Title:** Cut 1 Rejected Events — ZDC–FCal Banana

**Figure:** `event_sel_cut1_ZDC_FCal_banana_fail_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for events rejected by the ZDC–FCal banana cut. The rejected population sits in the high-ZDC-energy tail, visible in the top-left panel. Note these events show no significant anomaly in other variables, confirming the cut selects purely on ZDC–FCal correlation."

---

### SLIDE 5 — Cut 1: Events Passing the ZDC–FCal Banana

**Title:** Cut 1 Accepted Events — ZDC–FCal Banana

**Figure:** `event_sel_cut1_ZDC_FCal_banana_pass_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for events passing the ZDC–FCal banana cut. The ZDC–FCal correlation is now cleanly bounded, and all other distributions look nominal."

---

### SLIDE 6 — Cut 2: ZDC Time Box

**Title:** Cut 2 — ZDC Time Box

**Figure:** `event_sel_cut2_ZDC_time_standalone_pbpb_2024.png`

**Caption:** "ZDC time on side A vs side C for events passing cut 1. The red dashed lines mark the ±1.8 ns acceptance box. Genuine collision events produce simultaneous ZDC signals on both sides and cluster tightly at the origin. Events outside the box are out-of-time backgrounds (beam–gas interactions, satellite bunches). Survival rate for top 80% centrality: **99.99%** — nearly all central events have good ZDC timing."

---

### SLIDE 7 — Cut 2: Events Failing the ZDC Time Box

**Title:** Cut 2 Rejected Events — ZDC Time Box

**Figure:** `event_sel_cut2_ZDC_time_fail_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for events rejected by the ZDC time cut. The ZDC time panel shows the rejected events lie outside the ±1.8 ns box. Other variables (FCal, nTrk) can still look collision-like, motivating this explicit timing requirement."

---

### SLIDE 8 — Cut 2: Events Passing the ZDC Time Box

**Title:** Cut 2 Accepted Events — ZDC Time Box

**Figure:** `event_sel_cut2_ZDC_time_pass_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for events passing the ZDC time cut. The ZDC time correlation is now tightly clustered at the origin."

---

### SLIDE 9 — Cut 3: ZDC Pre-sample Amplitude (EMG Fit)

**Title:** Cut 3 — ZDC Pre-sample Amplitude (Fit-derived Cut)

**Figure:** `event_sel_cut3_ZDC_preamp_standalone_pbpb_2024.png`

**Caption:** "Distribution of the ZDC pre-amplifier ADC sum for side A (left) and side C (right), for events entering cut 3 (passing cuts 1 and 2). The pre-sample amplitude measures the signal before the main ZDC amplifier; anomalously large values indicate saturated or noisy channels that should be excluded.

**Cut derivation:** The distribution is fit with an Exponentially Modified Gaussian (EMG) — a Gaussian core with an exponential right tail — described by parameters μ (peak position), σ (width), and λ (tail rate). The upper cut is placed at **μ + 6σ**, shown as the blue dashed vertical line. The fitted function is shown in red over the fit range.

**Survival rate for top 80% centrality: 94.76%** — this is the largest rejection of the five cuts.

**Note on fit quality:** The EMG fit appears to slightly underestimate the width of the Gaussian core: the fitted curve sits below the data points on the left (low-ADC) side of the peak. This means the fitted σ is smaller than the true Gaussian width, causing μ + 6σ to fall lower than it should and resulting in a tighter-than-intended cut. This likely contributes to the relatively low survival rate of 94.76%."

---

### SLIDE 10 — Cut 3: Events Failing the Pre-sample Cut

**Title:** Cut 3 Rejected Events — ZDC Pre-sample Amplitude

**Figure:** `event_sel_cut3_ZDC_preamp_fail_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for events rejected by the ZDC pre-sample cut. These events have large pre-amplifier ADC sums on at least one ZDC side, indicating noisy or saturated ZDC modules."

---

### SLIDE 11 — Cut 3: Events Passing the Pre-sample Cut

**Title:** Cut 3 Accepted Events — ZDC Pre-sample Amplitude

**Figure:** `event_sel_cut3_ZDC_preamp_pass_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for events passing the ZDC pre-sample cut. The ZDC–FCal banana and time structures are preserved, and the nTrk–FCal correlation remains tight."

---

### SLIDE 12 — Cut 3 Alternative: Hard-coded Pre-sample Cut at 385 ADC

**Title:** Cut 3 Alternative — ZDC Pre-sample Hard Cut at 385 ADC

**Figure (from alternative directory):** `event_selection_alternative/event_sel_cut3_ZDC_preamp_standalone_pbpb_2024.png`

**Caption:** "Same ZDC pre-sample distributions as the previous slide, now with a hard-coded cut at **385 ADC for both sides** (blue dashed line), shown without fitting. This value is a direct estimate of the physical 'turning point' of the distribution — the ADC value above which the pre-sample sum is inconsistent with a clean ZDC signal. It is also consistent with independent μ + 6σ results from Mason's analysis.

**Survival rate for top 80% centrality: 95.90%** — higher than the fit-derived cut (94.76%), because 385 ADC is a looser (higher) threshold that does not suffer from the EMG underestimating the Gaussian width. This alternative is therefore a more robust and physically motivated choice."

---

### SLIDE 13 — Cut 4: nTrk HItight Fraction

**Title:** Cut 4 — nTrk HItight Fraction vs Total nTrk

**Figure:** `event_sel_cut4_nTrk_frac_standalone_pbpb_2024.png`

**Caption:** "Ratio of HItight-quality tracks to total tracks, as a function of total track count, for events passing cuts 1–3. For genuine collision events, this ratio is stable (~0.85) across all multiplicities. The red curve is the derived lower boundary: for each slice in total nTrk, the fraction distribution is fit with a Gaussian and the cut is placed at **μ − 5σ**. Events below this boundary have an anomalously low fraction of tight-quality tracks, indicating possible beam-induced backgrounds or reconstruction failures. Survival rate for top 80% centrality: **99.96%**."

---

### SLIDE 14 — Cut 4: Events Failing the nTrk Fraction Cut

**Title:** Cut 4 Rejected Events — nTrk HItight Fraction

**Figure:** `event_sel_cut4_nTrk_frac_fail_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for events rejected by the nTrk fraction cut. The bottom-left panel shows these events have abnormally low HItight fraction."

---

### SLIDE 15 — Cut 4: Events Passing the nTrk Fraction Cut

**Title:** Cut 4 Accepted Events — nTrk HItight Fraction

**Figure:** `event_sel_cut4_nTrk_frac_pass_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for events passing the nTrk fraction cut. The track quality distribution is now cleanly within the expected band."

---

### SLIDE 16 — Cut 5: nTrk HItight vs FCal ET Band

**Title:** Cut 5 — nTrk HItight vs FCal ET Correlation Band

**Figure:** `event_sel_cut5_nTrk_FCal_band_standalone_pbpb_2024.png`

**Caption:** "Number of HItight tracks vs FCal ET for events passing cuts 1–4. For genuine Pb+Pb collisions, nTrk HItight follows a tight linear correlation with FCal ET. The two red curves are the derived upper and lower boundaries of a **μ ± 5σ band**: for each FCal ET slice, the nTrk distribution is fit with a Gaussian and the cut keeps events within 5σ of the mean. Events outside this band have anomalous track multiplicity for their centrality and are removed. Survival rate for top 80% centrality: **99.92%**."

---

### SLIDE 17 — Cut 5: Events Failing the nTrk–FCal Band

**Title:** Cut 5 Rejected Events — nTrk HItight vs FCal ET

**Figure:** `event_sel_cut5_nTrk_FCal_band_fail_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for events rejected by the nTrk–FCal band cut. The bottom-center panel shows these events lie above or below the expected nTrk–FCal correlation."

---

### SLIDE 18 — Cut 5: Events Passing the nTrk–FCal Band

**Title:** Cut 5 Accepted Events — nTrk HItight vs FCal ET

**Figure:** `event_sel_cut5_nTrk_FCal_band_pass_5panel_pbpb_2024.png`

**Caption:** "Five-panel distributions for the final selected event sample, passing all five cuts. All correlation structures are clean and consistent with genuine Pb+Pb collisions."

---

### SLIDE 19 — Survival Rate Summary

**Title:** Event Selection Survival Rates — Top 80% Centrality

**Figure:** `event_sel_survival_pbpb_2024.png`

**Caption:** "Per-cut survival fraction for events with centrality 0–79 (top 80%, the analysis acceptance). Each point shows the fraction of events entering that cut which pass it, with binomial statistical error bars. Cuts 1, 2, 4, and 5 each pass >99.9% of events. Cut 3 (ZDC pre-sample, fit-derived) passes 94.76%, making it the dominant rejection cut. The alternative hard-coded cut at 385 ADC raises this to 95.90% (shown separately in slides 12)."

---

## Formatting Notes for the Assistant

- Use a clean, minimal slide design. White or light background.
- Figures should fill most of the slide area; keep titles and captions readable but compact.
- The 5-panel figures are wide (1500×900 px); display them landscape filling the slide width.
- The standalone cut figures are 1200×600 px (two-panel); display similarly landscape.
- The survival rate figure (800×600 px) can be centered with some white space on the sides.
- Dataset label visible on all plots: "Pb+Pb 2024 data".
- Total slides: 19.
