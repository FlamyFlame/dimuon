# HIJING Overlay: Combinatorial Background and Single-B Statistics

Test sample: HIJING 0–5 fm impact parameter (ultra-central Pb+Pb), Pythia bb→μμ overlaid.

---

## 1. Combinatorial Background

### Pair structure in the overlay NTUP

The overlay `truth_mupair_*` branches store pairs of **one Pythia signal truth muon** paired with **one reconstructed muon from the combined (Pythia+HIJING) event**. This is not the same as PP fullsim, where both muons in every pair are Pythia-generated truth muons.

Because most reconstructed "second muons" come from HIJING soft muon production (π/K decays), the pair pool is dominated by combinatorial (signal×background) pairs.

### Counts (test sample: 60k events total, 6 kn bins × 10k events)

| Centrality | OS pairs (`_op`) | SS pairs (`_ss`) | SS/OS |
|-----------|-----------------|-----------------|-------|
| ctr0–5 | 67,362 | 20,790 | 31% |
| ctr5–10 | 14,888 | 4,842 | 31% |
| ctr10–20 | 554 | — | — |
| ctr20+ | 0 | — | — |

Only the two most central bins are populated because the HIJING sample uses b < 5 fm impact parameter selection.

### Effect on reconstruction efficiency

The `_op` reco efficiency uses all OS pairs as denominator. Since the denominator is dominated by HIJING combinatorial background (Pythia×HIJING pairs), the apparent OS efficiency in ctr0–5 is:

| Category | Denom | Pass medium | Apparent efficiency |
|----------|-------|-------------|---------------------|
| `_op_ctr0_5` | 67,362 | 612 | ~0.9% |

This does **not** represent the signal reconstruction efficiency; it measures the probability that a random HIJING combinatorial OS pair passes medium selection. The `_single_b` efficiency is the physically meaningful quantity (see §2).

### SS/OS symmetry

The SS/OS ratio is ~31% (not ~100%) because the genuine Pythia signal pair (μ+×μ−) contributes only to OS. With N_HIJING soft muons per event, all with equal probability of + or −:
- OS ≈ N_HIJING + 1 (HIJING×Pythia cross-pairs + signal pair)
- SS ≈ N_HIJING (HIJING×Pythia same-sign cross-pairs)

The imbalance reflects that Pythia contributes exactly one OS pair per event with no SS counterpart.

### Strategies to deal with combinatorial background

| Method | Signal purity | Status |
|--------|--------------|--------|
| `_single_b` (same-b truth flag) | 100% | Implemented; see §2 for structural limitation |
| OS − SS subtraction | High (leading order) | SS per-centrality histograms now filled |
| `_op` (all OS) | ~0.01% | Misleading; documents combinatorial level only |

---

## 2. Single-B Statistics and Structural Limitation

### Comparison: overlay vs PP fullsim

| Sample | Events (total) | `single_b` (truth, all) | `single_b` pass medium | Rate per event |
|--------|----------------|------------------------|----------------------|---------------|
| PP fullsim | 240k (40k/kn × 6) | 93,400 | 67,516 (72%) | 0.39 |
| Overlay (HIJING) | 60k (10k/kn × 6) | **6** | 0 | ~0.0001 |
| Expected if scaled | 60k | ~23,350 | ~17,000 | 0.39 |
| Actual / expected | — | **×4,000 lower** | — | — |

### What `from_same_b` means

`from_same_b = true` requires both muons in a pair to trace to the **same b-hadron** (same `m1_eldest_bhadron_barcode == m2_eldest_bhadron_barcode`). The dominant topology is the b→c→μ cascade: a single B meson produces two muons — one direct (B→Dμν) and one secondary (D→μν). This is the open-bottom dimuon signal of interest (~39% of in-acceptance events in PP fullsim).

### Root cause of the 4,000× suppression

The suppression is **not from combinatorial background**. Combinatorial pairs add more entries to the OS denominator but do not remove signal pairs. The root cause is structural: the overlay NTUP `truth_mupair_*` branches store **(Pythia signal μ) × (HIJING reco μ)** pairs, not **(Pythia μ₁) × (Pythia μ₂)** pairs.

For these pairs:
- `m1_eldest_bhadron_barcode`: valid (Pythia signal muon traces to a b-hadron)
- `m2_eldest_bhadron_barcode`: **−10** (HIJING soft muon traces to π/K, not a b-hadron; ancestry trace terminates without finding a b-flavored hadron)
- Condition `m1_bcode != -10 && m1_bcode == m2_bcode` → **always false**

The 6 observed `from_same_b` events are rare HIJING events where the second muon happens to come from a b-quark produced in the HIJING heavy-ion collision itself.

### Per-centrality single-B truth counts

| Bin | `single_b` denom | pass_medium | pass_tight | det response |
|-----|-----------------|-------------|------------|--------------|
| ctr0–5 | 4 | 0 | 0 | 0 |
| ctr5–10 | 2 | 0 | 0 | 0 |
| ctr10–20 | 0 | 0 | 0 | 0 |

- The 4+2 pairs in the top two bins are the rare HIJING-b cases, not Pythia signal.
- Zero pass medium/tight: HIJING-origin muons embedded in ultra-central events fail quality selection.
- Detector response histograms require `_single_b_pass_medium` (non-zero reco) → identically empty.

### Path forward

The current overlay pair structure does not permit `from_same_b` efficiency measurements. Two approaches exist:

1. **Store Pythia-only truth pairs in the overlay NTUP** (skimmer change): add a separate branch set `truth_mupair_signal_*` containing only (Pythia truth μ₁) × (Pythia truth μ₂) pairs. This is the clean solution.

2. **Process Pythia truth pairs separately in the NTP** (analysis-side): after the standard overlay pair loop, iterate over the Pythia truth muon pair (two signal muons per event), match each to a reco muon, and record the efficiency for that signal pair separately. Requires access to `truth_mupair_bar1/bar2` that specifically references both Pythia signal barcodes.

Both require that the full sample (1.6M–3.2M events/kn) is available, at which point the signal yield at full statistics would give ~{630k–1.25M} true single-b pairs in ctr0–5 — ample for efficiency maps.
