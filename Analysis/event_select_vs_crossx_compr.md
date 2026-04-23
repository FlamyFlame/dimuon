# Event Selection vs Cross-Section: File Chain Comparison (PbPb 2024)

## 1. Pair ntuple production (`NTupleProcessingCode/`)

Two condor jobs write to **different output files**:

| Job | Extra flag | `resonance_cut_mode` | Output pair ntuple |
|-----|-----------|----------------------|--------------------|
| `run_pbpb_24.sub` | *(none)* | 2 (`_res_cut_v2`) forced by `trigger_effcy_calc=true` | `muon_pairs_pbpb_2024_single_mu4_mindR_0_02_res_cut_v2.root` |
| `run_pbpb_24_nominal.sub` | `pbpb24_mu4_NO_trig_calc = true` | 1 (old cuts, empty suffix) | `muon_pairs_pbpb_2024_single_mu4_mindR_0_02.root` |

Both use the same trigger (`single_mu4`, trigger_mode=1) and the same photoproduction/resonance filtering
(`filter_out_photo_resn_for_trig_effcy = true` by default), so pairs are filtered identically.
The only difference is the resonance cut version applied (v1 vs v2).

**The nominal sub is the cross-section path** (res cut v1, no `_res_cut_v2` suffix).

---

## 2. Hist-filling input (pair ntuple picked up)

`RDFBasedHistFillingPbPb::SetIOPathsHook` tries these paths in order (with default `mindR_trig=0.02`):

```
1. muon_pairs_pbpb_2024_single_mu4_mindR_0_02_res_cut_v2.root   ← picked if exists
2. muon_pairs_pbpb_2024_single_mu4_mindR_0_02.root
3. muon_pairs_pbpb_2024_single_mu4_no_res_cut.root
4. muon_pairs_pbpb_2024_single_mu4.root
```

Currently `_res_cut_v2.root` exists (from `run_pbpb_24.sub`), so **hist-filling picks that up first**,
regardless of whether it is a crossx or event-selection run.

**Both crossx and event-selection hist-filling read the same input file** — they run inside the same
`Run()` → `FillHistograms()` call, which calls both `FillHistogramsCrossx()` and
`FillHistogramsEventSelection()` sequentially.

---

## 3. Hist-filling runs and output files

`FillHistogramsCrossx` and `FillHistogramsEventSelection` always run together in one invocation.
What differs between the two invocations (crossx vs event-selection) is the `doTrigEffcy` flag,
which changes the **output filename** via `trig_suffix`:

| Invocation | `doTrigEffcy` | `trig_suffix` | Output histogram file |
|------------|---------------|---------------|-----------------------|
| Cross-section | `true` (default) | `_single_mu4` | `histograms_real_pairs_pbpb_2024_single_mu4_coarse_q_eta_bin.root` |
| Event selection | `false` | `_single_mu4_no_trg_plots` | `histograms_real_pairs_pbpb_2024_single_mu4_no_trg_plots_coarse_q_eta_bin.root` |

Both files are written to `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024/`.

Each output file contains **both** crossx and event-selection histograms (since both fill functions
run in the same call); the distinction is only in the output path.

---

## 4. Plotting inputs

| Plotter | Input histogram file |
|---------|---------------------|
| `plot_single_b_crossx_pbpb.cxx` | `histograms_real_pairs_pbpb_2024_single_mu4_coarse_q_eta_bin.root` (falls back to `_single_mu4.root`) |
| `plot_pbpb_event_selection.cxx` | `histograms_real_pairs_pbpb_2024_single_mu4_no_trg_plots_coarse_q_eta_bin.root` (hardcoded) |

The cross-section plotter resolves the trigger string via `DatasetTriggerMap::GetTrigger(24, "PbPb") = "single_mu4"`.

---

## 5. Summary of the full chain

```
run_pbpb_24_nominal.sub
  └─> muon_pairs_pbpb_2024_single_mu4_mindR_0_02.root          (res cut v1, nominal)

run_pbpb_24.sub
  └─> muon_pairs_pbpb_2024_single_mu4_mindR_0_02_res_cut_v2.root  (res cut v2, trig effcy)
         │
         └─ [picked first by hist-filling if it exists]

Hist-filling (doTrigEffcy=true, crossx run):
  input:  _res_cut_v2.root  (priority 1)
  output: histograms_real_pairs_pbpb_2024_single_mu4_coarse_q_eta_bin.root
            └─> plot_single_b_crossx_pbpb.cxx

Hist-filling (doTrigEffcy=false, event-selection run):
  input:  _res_cut_v2.root  (same priority 1)
  output: histograms_real_pairs_pbpb_2024_single_mu4_no_trg_plots_coarse_q_eta_bin.root
            └─> plot_pbpb_event_selection.cxx
```

**Note:** As long as `_res_cut_v2.root` exists on disk, both hist-filling runs pick it up regardless
of which sub was intended. To force hist-filling onto the nominal (`_mindR_0_02.root`) file,
the `_res_cut_v2.root` file must either not exist or `mindR_trig` must be set to a value that
produces no `_res_cut_v2` match.
