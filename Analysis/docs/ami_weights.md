# AMI weights — the MC cross-section registry (**BLOCKING**)

**Authoritative registry of the AMI cross-section weights for every MC dataset this analysis
uses, and the rule that governs them.** Verified against `pyami` on **2026-07-14**; evgen-keyed
layout and the "weight = the Pythia evgen's" rule from the user ruling of **2026-09-17**.

---

## THE RULE (hard blocker — read before touching any MC dataset)

> **AMI weights are a HARD BLOCK on every MC dataset.**
> **When you move to a NEW MC dataset — test → full sample, a new overlay, a different
> generator, a new production tag, ANY new DSID — you MUST fetch that dataset's OWN AMI
> weights from `pyami` BEFORE running any analysis on it.**
> **NEVER reuse the old dataset's weights. That is a SILENT failure that propagates all the
> way to the final results.**

### The weight is a property of the PYTHIA EVGEN, not of the AOD chain (user ruling 2026-09-17)

The per-event weight of every Pythia sample — truth-only, pp24 fullsim, HIJING overlay, a future
data overlay — is `σ · ε_filt` of the **Pythia evgen dataset the event was generated in**. A fullsim
or overlay production does not re-generate events; it takes the already-generated EVNT as input to
simulation. That is why the Pythia e-tag does not even appear in an overlay's AMI tag chain: in
`e8613_e8586_s4684_r17864_r17855` (the pbpb24 overlay test sample) **`e8613` is the HIJING evgen**
(`860250.Hijing_PbPb_UCC_Flow_JJFV6_ip0_5.e8613_s4684_s4688`), not a Pythia tag. Querying AMI on the
AOD dataset and calling the result "this sample's own AMI weights" (done 2026-09-16, registry B2,
reverted 2026-09-17) is therefore wrong in principle, even where the numbers happen to agree.

The same principle is why only **Pythia** truth muons enter the overlay's reconstruction- and
trigger-efficiency numerators and denominators, never HIJING truth muons: each event carries the
Pythia evgen weight, and HIJING is background / environment / underlying event for the Pythia truth
particles. A HIJING muon counted with the Pythia weight would carry a wrong weight.

**Two Pythia evgen productions exist at 5.36 TeV** (`FullSimSampleType.h`, `enum PythiaEvgen`):

| evgen | e-tag | DSIDs | isospin | PDF | suitable for | used by |
|---|---|---|---|---|---|---|
| **nPDF** | e8599 | 802758–802781 | pp/pn/np/nn, Pb ratio 4:6:6:9 | `LHAPDF6:nNNPDF30_nlo_as_0118_A208_Z82/0001` (nuclear) | **Pb+Pb** conditions | truth-only analysis; pp24 **TEST** fullsim (produced on it by mistake); every HIJING overlay (pbpb23 r17618/r17662, pbpb24 r17864, the full production); future data overlay |
| **PDF** | e8599 (`_pdf`) | 803015–803020 | pp only | `NNPDF23_lo_as_0119_qed` (proton) | **pp** conditions | pp24 fullsim **FULL** sample; a truth-only skim of it, should one be made |

Their `σ·ε_filt` differ **slice-dependently** (table B, FULL/TEST 0.879–1.545), so choosing the
wrong evgen cancels in no ratio — the same failure mode as reusing another production's file.

**The AMI files live with the evgen, in ONE place each** (fetch scripts next to them):

```
~/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/ami_info_nPDF/  ami_info_mc23_5p36TeV_Py8EG_A14_<beam>_hQCD_DiMu_pTH<lo>_<hi>.txt      (fetch_ami_info_nPDF.sh)
~/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/ami_info_PDF/   ami_info_mc23_5p36TeV_Py8EG_A14_pp_hQCD_DiMu_pTH<lo>_<hi>_pdf.txt    (fetch_ami_info_PDF.sh)
```

`PythiaAlgCoreT` derives `ami_evgen` from the same switch as the input directory and the isospin
treatment (`FullSimSampleEvgen(type, isTestSample)`: pp24 FULL → PDF, everything else → nPDF) and
reads `PythiaEvgenAmiDir(ami_evgen)`. **A copy of AMI files inside a fullsim sample directory is
never read by code** (the 2026-09-16 "prefer the sample's own `ami_info/`" logic is gone).

### Why it is silent, and why it is worse than it looks

The AMI info files are named by **beam + pT-hat slice ONLY**:

```
ami_info_mc23_5p36TeV_Py8EG_A14_<beam>_hQCD_DiMu_pTH<lo>_<hi>.txt
```

The filename is **byte-identical between productions**. Nothing in the path, the file name, or
the NTUP names distinguishes the test sample from the full sample. Point the code at a new
dataset while leaving the old `ami_info/` in place and **every file opens successfully, every
number comes out, and every one of them is wrong.** There is no crash, no warning, nothing to
notice in a plot review.

**And it does NOT cancel.** This is the crucial point:

| error type | slice-dependent? | cancels in ratios? |
|---|---|---|
| isospin weight (the old 4/25 bug) | no — same constant for every slice | **YES** — drops out of reco-eff, det-response, MC trig-eff |
| **AMI cross-section from the wrong production** | **YES** | **NO — cancels NOWHERE** |

A wrong σ·ε_filt **reweights the pT-hat mixture**. The MC trigger efficiency is a
σ-weighted average over slices *within each (pT, η) bin*, so a per-slice error survives the
numerator/denominator ratio. Every cross-section is wrong outright. Measured for the pp24
sample: the FULL/TEST ratio of σ·ε_filt spans **0.879 – 1.545** across the six slices — a
**76 % spread**, not a harmless overall scale.

### The guard (enforced in code)

`PythiaAlgCoreT::InitInputFullsim` parses `datasetNumber` out of each AMI file and **THROWS**
if it is not in `expected_ami_dsids`:

```
InitInputFullsim: AMI PROVENANCE MISMATCH for pTH8_14 beam pp: <path>
  has datasetNumber=802781, which is NOT in the expected DSID list for this sample {803015,...}
```

- `ami_evgen` — which Pythia evgen's AMI directory to read; derived from
  `(fullsim_sample_type, isTestSample)`, never set by hand (`ami_info_dir_override` is a
  diagnostic escape hatch only).
- `expected_ami_dsids` — the DSIDs this run is allowed to see; defaults to the evgen's list,
  run scripts may declare/narrow it.
- Both follow from the single `isTestSample` switch (default `false` = the FULL production), so
  the input files and their cross-sections can never come from different productions.
- A **missing** AMI file is now **fatal** too (it used to leave `ami_weight = 0`, silently giving
  that pT-hat slice **zero weight**).

### Checklist when adopting a NEW MC dataset

1. Get the DSIDs (`rucio list-dids`).
2. `lsetup pyami`; `ami show dataset info <dataset>` for **every** DSID.
3. Identify the **Pythia evgen** the dataset was simulated from (the AMI `ldn` of the EVNT, not of
   the AOD). If it is one of the two known evgens (nPDF / PDF above), its files already exist. If it
   is a NEW Pythia evgen, add a `PythiaEvgen` value + `ami_info_<name>/` next to the two existing
   ones, and fetch the EVNT dataset's info for **every** DSID into it — never overwrite another
   evgen's files, never copy them into the fullsim sample directory.
4. Make `FullSimSampleEvgen()` map the new sample to its evgen, and declare `expected_ami_dsids` in
   the run script (defaults to the evgen's DSID list).
5. Add the numbers to the table below, with the date they were verified.
6. Re-run **everything** that consumes that sample's weights (see §Blast radius).

---

## Units

**AMI `crossSection` is in nb, NOT pb** (`project_ami_crosssection_nb_units`). The comment in
`PythiaAlgCoreT` that once said pb was wrong. The per-event MC weight is

```
w_s = σ_s · ε_filt,s · r_isospin / N_s      [nb]
```

with `r_isospin` = 1 for a single-beam (pp-conditions) sample and the Pb ratio 4:6:6:9 for the
4-beam PbPb overlay (`ami_weights` ↔ `Analysis/docs/tracking/pythia_fullsim_pp24_full_sample_skim.md`
D2). Comparing an absolute MC dσ to **pp data** (dσ = N/L, L in pb⁻¹) requires **nb → pb, ×1000**.

---

## Registry (verified against `pyami` 2026-07-14 — both sets: **no drift**)

#### A. **nPDF evgen** `e8599`, DSIDs 802758–802781 — truth-only analysis, pp24 fullsim **TEST** sample, EVERY HIJING overlay (pbpb23 r17618/r17662, pbpb24 r17864, the full production)

The HIJING overlays are built on the **pp-beam** evgen DSIDs (802776–802781) — the AMI weight comes
from the **evgen** dataset, so neither the overlay's r-tag (r17618 / r17662 / r17864) nor its HIJING
e-tag (e8613) changes it.

| slice | beam | DSID | σ [nb] | genFiltEff | **σ·ε_filt [nb]** |
|---|---|---|---|---|---|
| pTH8_14 | pp | 802781 | 4.816e+06 | 4.938681e-06 | **23.7847** |
| pTH8_14 | pn | 802775 | 4.816e+06 | 4.937969e-06 | **23.7813** |
| pTH8_14 | np | 802769 | 4.816e+06 | 4.959556e-06 | **23.8852** |
| pTH8_14 | nn | 802763 | 4.816e+06 | 4.958479e-06 | **23.8800** |
| pTH14_24 | pp | 802777 | 726090 | 4.392067e-05 | **31.8904** |
| pTH14_24 | pn | 802771 | 726100 | 4.388425e-05 | **31.8644** |
| pTH14_24 | np | 802765 | 726090 | 4.385502e-05 | **31.8427** |
| pTH14_24 | nn | 802758 | 726090 | 4.390211e-05 | **31.8769** |
| pTH24_40 | pp | 802778 | 100200 | 0.0001929893 | **19.3375** |
| pTH24_40 | pn | 802772 | 100200 | 0.0001929752 | **19.3361** |
| pTH24_40 | np | 802766 | 100200 | 0.0001927405 | **19.3126** |
| pTH24_40 | nn | 802760 | 100200 | 0.0001930126 | **19.3399** |
| pTH40_70 | pp | 802779 | 13860 | 0.0005569396 | **7.7192** |
| pTH40_70 | pn | 802773 | 13860 | 0.0005616423 | **7.7844** |
| pTH40_70 | np | 802767 | 13860 | 0.0005613345 | **7.7801** |
| pTH40_70 | nn | 802761 | 13860 | 0.0005581628 | **7.7361** |
| pTH70_125 | pp | 802780 | 1296.8 | 0.001280053 | **1.6600** |
| pTH70_125 | pn | 802774 | 1296.7 | 0.001280214 | **1.6601** |
| pTH70_125 | np | 802768 | 1296.8 | 0.001283622 | **1.6646** |
| pTH70_125 | nn | 802762 | 1296.7 | 0.001282946 | **1.6636** |
| pTH125_300 | pp | 802776 | 89.529 | 0.002342579 | **0.2097** |
| pTH125_300 | pn | 802770 | 89.534 | 0.002341105 | **0.2096** |
| pTH125_300 | np | 802764 | 89.529 | 0.00234476 | **0.2099** |
| pTH125_300 | nn | 802759 | 89.521 | 0.002344699 | **0.2099** |

AMI files: `~/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/ami_info_nPDF/` (this IS the
Pythia truth-only production's evgen).

#### B. **PDF evgen** `_pdf` (e8599), DSIDs 803015–803020 — pp24 fullsim **FULL** sample

**pp beam ONLY**, proton PDF `NNPDF23_lo_as_0119_qed` (this evgen was re-generated for pp
conditions; isospin weight 1 — see D2 of the tracking doc).

| slice | beam | DSID | σ [nb] | genFiltEff | **σ·ε_filt [nb]** | FULL/TEST(pp) |
|---|---|---|---|---|---|---|
| pTH8_14 | pp | 803020 | 6.1995e+06 | 5.926805e-06 | **36.7432** | **1.545** |
| pTH14_24 | pp | 803016 | 816620 | 4.363898e-05 | **35.6365** | 1.117 |
| pTH24_40 | pp | 803017 | 101810 | 0.0001821213 | **18.5418** | 0.959 |
| pTH40_70 | pp | 803018 | 13114 | 0.0005262908 | **6.9018** | 0.894 |
| pTH70_125 | pp | 803019 | 1177 | 0.00123943 | **1.4588** | **0.879** |
| pTH125_300 | pp | 803015 | 82.548 | 0.002395971 | **0.1978** | 0.943 |

AMI files: `~/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/ami_info_PDF/` (`_pdf.txt`
suffix). The copy that lived in `pythia_fullsim_full_sample/ami_info/` (fetched 2026-07-13 on the
AOD dataset; identical σ and ε_filt) was renamed `ami_info_AOD_record_20260713/` on 2026-09-17 and
is not read by anything.

> **The FULL/TEST ratio spans 0.879 – 1.545.** It is **slice-dependent**, so it cancels in
> **nothing** — not in the MC trigger efficiency, not in any cross-section. This is exactly the
> silent failure the rule above exists to prevent.

#### B2 (RETIRED 2026-09-17). "pbpb24 overlay test sample, its own AMI"

Registered 2026-09-16 from `ami show dataset info` on the **AOD** dataset
`802776...merge.AOD.e8613_e8586_s4684_r17864_r17855` (σ 89.541 nb, ε_filt 2.342314e-3, 1e-4 from
table A). Wrong in principle: e8613 is the HIJING evgen tag; the sample's Pythia evgen is 802776
e8599 = table A, which is what it is weighted with now. The fetched file was moved to
`pythia_fullsim_hijing_overlay_test_sample/ami_info_AOD_record_20260915/` (record only).

#### C. Not yet fetched — **fetch before use**

| sample | status |
|---|---|
| HIJING-overlay **FULL** sample (PbPb conditions, 4 isospin beams) | **in production** on the nPDF evgen 802758–802781 → table A applies (4 beams × 6 slices). **Verify the EVNT in the production's input chain, do not assume**; a new Pythia evgen would need its own `ami_info_<name>/`. |
| Truth-only skim of the **PDF** evgen (pp-suitable generator-level sample) | not skimmed; if ever made, table B applies unchanged (same EVNT datasets). |
| POWHEG (truth, NLO template) | uses its own weights (`weight_norm`); not in this registry yet. |

---

## Blast radius — what to re-run when an AMI weight changes

An AMI change is **slice-dependent**, so it invalidates **everything weighted**, not just absolute
cross-sections:

- **NTuple processing** of that sample (the weight is baked into the pair/single-muon trees at
  `PythiaFullSimExtras.c` `event_weight = fullsim_weight_factor`), then
- **RDF hist filling**, then
- **every plot**: absolute dσ/dp_T, the MC trigger efficiency (a σ-weighted average over slices),
  MC/data comparisons, and any template built from that sample.

Contrast with a **slice-independent** factor (e.g. the isospin weight): that one cancels in
reco-eff, detector response, MC trig-eff and area-normalized templates, and only moves absolute
normalizations.

---

## References

- Code guard: `Analysis/MuonObjectsParamsAndHelpers/FullSimSampleType.h` (`PythiaEvgen`,
  `FullSimSampleEvgen`, `PythiaEvgenAmiDir/FileName/Dsids`),
  `Analysis/NTupleProcessingCode/PythiaAlgCoreT.{h,c}` (`ami_evgen`, `expected_ami_dsids`,
  `isTestSample`; truth-only reader `InitInputCentrProd` → nPDF, missing file fatal).
- Ruling record: `Analysis/docs/tracking/hijing_overlay_pbpb24_test_sample_skim.md` §"AMI-weight
  correction (2026-09-17)".
- Tracking: `Analysis/docs/tracking/pythia_fullsim_pp24_full_sample_skim.md` (D2 isospin, D3 AMI guard).
- Units: `project_ami_crosssection_nb_units` — AMI σ is in **nb**.
- Provenance convention: `.claude/conventions/ntuple-provenance.md`.
