# Branch Names, Conventions, and Units

## Sign convention (NTuple processing output trees)

| Tree name | Meaning | Output suffix |
|-----------|---------|---------------|
| `muon_pair_tree_sign1` | Same-sign pairs (SS) | `_ss` |
| `muon_pair_tree_sign2` | Opposite-sign pairs (OS) | `_op` |

Verified by from_same_b fractions: sign1 ~0% (SS), sign2 ~25% (OS — B-meson sequential decay).

## ZDC array indexing

All ZDC arrays (`zdc_ZdcEnergy[2]`, `zdc_ZdcTime[2]`, `zdc_ZdcModulePreSampleAmp[2][4]`, etc.):
- **Index 0 = Side C** (eta < 0)
- **Index 1 = Side A** (eta > 0)

Source: `TrigRates.cxx::ProcessZdc`: `iside = (zdcside > 0) ? 1 : 0` where zdcSide() > 0 = A-side.

## FCal variables

- `FCal_Et_P` (eta > 0) = Side A
- `FCal_Et_N` (eta < 0) = Side C
- FCal_Et_P/N are stored in **MeV** in the skim; divide by 1e6 for TeV

## ZDC preamp sum computation

Declare `Float_t zdc_ZdcModulePreSampleAmp[2][4]` (with SetMakeClass=1).
- preamp_A = sum of `[1][0..3]` (side A)
- preamp_C = sum of `[0][0..3]` (side C)

## Event selection cut storage

All derived cuts saved to `event_sel_cuts_pbpb_20YY.root`:

| Key | Type | Description |
|-----|------|-------------|
| `g_ZDC_FCal_cut` | TGraph | ZDC vs FCal banana cut |
| `g_ntrk_frac_cut_lo` | TGraph | nTrk fraction lower bound |
| `g_ntrk_fcal_cut_lo` | TGraph | nTrk vs FCal lower band |
| `g_ntrk_fcal_cut_hi` | TGraph | nTrk vs FCal upper band |
| `ZDC_time_cut_ns` | TParameter | Time cut threshold |
| `ZDC_preamp_A_cut_ADC` | TParameter | Preamp A scalar cut |
| `ZDC_preamp_C_cut_ADC` | TParameter | Preamp C scalar cut |
| `t_preamp_per_run` | TTree | Per-run preamp cuts (PbPb25 only) |

## Skim output branches

See `SkimCode/README.md` sections "Output branch reference" for `StoreTracks`, `StoreZdc`, and `StoreEventInfo` bitmaps and their branch layouts.
