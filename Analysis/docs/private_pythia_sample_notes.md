# Private Pythia samples — conventions, weights, and pitfalls (reference)

Written 2026-10-07 for the no-muon-p_T-filter acceptance sample (`docs/tracking/signal_acceptance_mass_window.md`
PD4/PD9). Read before generating, processing, or reading ANY private Pythia sample. Paths: `PSC` =
`dimuon_codes/PythiaSimulationCode`, `NTP` = `Analysis/NTupleProcessingCode`, `JO` =
`dimuon_codes/PythiaCentralProdJO/evgen_high_pT_additional`.

There are TWO kinds of private Pythia sample in this project. They follow DIFFERENT conventions:

| | (A) standalone-Pythia "private sample" | (B) private Athena Gen_tf sample | ATLAS central |
|---|---|---|---|
| generator | `PSC/RunHardQCD_dimuon.cc` (own C++ main) | Gen_tf + our JO (`JO/run/`) | Gen_tf + JO, AMI-registered |
| format | `PyTree` + `meta_tree` + TEnv `pythia_params` | EVNT → TRUTH0 → SkimCode `TruthTrigRates` NTUP (`HeavyIonD3PD`) | same as (B) |
| NTuple mode | `is_private=true` | EXTERNAL mode (`truth_input_override` + `truth_raw_counts`, or production-weighted `truth_external_sigma_eff_nb`; `docs/pythia_truth.md`) | nonprivate |
| b/c-hadron decays | Pythia internal tables | **EvtGen** (`2014Inclusive` dec) | EvtGen |
| normalisation | per-file `eventWeight` in µb | σ, ε_filt from the sample's OWN `log.generate` (nb), applied downstream | AMI σ·ε_filt (nb) |

## 1. The previous standalone private sample (A) — what it really is

- **It is 5.02 TeV, nCTEQ15 nuclear PDF, 4 isospin beams** (`pythia_params`: `Beams:eCM = 5020`,
  `nCTEQ15npFullNuc_208_82`), Pythia 8.307, 5 pTHat slices {5,10,25,60,120,3200}, muon filter ≥2 muons
  p_T > **3.5** GeV, |η| < 2.5. `docs/pythia_truth.md` used to say 5.36 TeV — wrong (fixed 2026-10-07).
- **It was produced by the 2023 code in the OLD repo** (`/usatlas/u/yuhanguo/workarea/pythia/SimulationCode/RunHardQCD_dimuon.cc`),
  NOT by the current `PSC/RunHardQCD_dimuon.cc` (since edited: 7 slices 5–8…125–300, p_T > 3.7, nNNPDF30 nPDF,
  eCM still 5020, test event counts; the compiled binary 2025-05-21 is older than the .cc 2026-03-09).
- Batches used (`NTP/PythiaAlgCoreT.c:139-169`): batch 1 `_allto0318` = {0317_all_k, 0318_all_k, 0318_k0}
  (meta branch `totalSigma`), batch 2 `_after0322` = {0322_k0_k1, 0323_k0_k1, 0325_all_k, 0401_all_k, 0429_all_k}
  (meta branch `eventWeight`; same meaning). 2025 batches (052025, 052125_nNNPDF30, 051425_single_muon) are unused.
  8 069 910 events total (k0 119 910 from 1.73e11 trials, filter eff 7e-7). Raw files have LGD replicas
  (`*_lgd_staging` symlink farms) but the local originals were never deleted.
- `PSC/submission_scripts/submit_hard_qcd_no_muon_filter.sh` is STALE: its `RunHardQCD_noMuonFilter` source no
  longer exists anywhere. **No no-muon-filter Pythia sample has ever been produced.**

## 2. Conventions that differ between (A) and (B)/central — must get right

1. **Barcode / index.** (A): `truth_barcode` = Pythia event-record index, **0-based, entry 0 = system (id 90,
   status −11)**, 1–2 = beams; `truth_mother1/2`, `truth_daughter1/2` index the same array (verified: barcode ==
   vector index in every event). (B)/central: HepMC barcodes, **1-based, == index+1**, no system entry, beams
   status 4; parents as `truth_parents` = list of parent BARCODES. `PythiaTruthExtras::GetParticleIndex`
   (`NTP/PythiaTruthExtras.c:334-371`) switches on `isPrivate` — never mix.
2. **Status.** (A): native Pythia status (final muons 91 ≈ 89 %, 92 ≈ 8 %, 62 = γ→μμ in the shower); "final" =
   `status > 0`. (B)/central: HepMC status — 1 final, 2 decayed, 4 beam, **otherwise |Pythia status|** kept for
   intermediates. The ancestry code relies on native 21/23/31/33 (hard-scatter/ISR/FSR) in BOTH modes; a record
   with only 1/2/4 would break that classification silently. (A) files on disk have NO `truth_statusHepMC`.
3. **Parent navigation.** (A): `mother1..mother2` with Pythia's range convention for |status| 81–86/101–106
   (all indices in [mother1, mother2]); `mother1 == 0` is an error (`PythiaTruthExtras.c:711-741`). (B): `truth_parents`.
4. **Units & types.** (A): doubles in GeV, keeps pTHat/mHat (branch `QHard` actually = `info.QFac()`).
   (B)/central: floats in MeV (÷1000 in `FillMuonPair`), pTHat = mHat = −1000.
5. **The HIJING barcode cutoff** (`pythia_only_barcode_cache`, first barcode > 200000) is overlay-only; it is
   OFF on every truth/private path and must stay off.
6. `from_same_b` (= same eldest B-hadron barcode for both muons) is convention-independent once
   `GetParticleIndex` is consistent.

## 3. Event weights — what a private weight MEANS

- (A) generator: `meta_tree.eventWeight = info.sigmaGen()[mb]·1e3 / info.weightSum()` = σ_slice / N_tried in
  **µb**, per file (`PSC/RunHardQCD_dimuon.cc:576`). Since ε_filt = N_acc/N_tried, this = σ·ε_filt / N_acc:
  each ACCEPTED event carries the cross section it represents. Verified: eventWeight·nTried = Σ processSigma·1e3
  to 1e-6 in every file ⇒ `weightSum == getCounter(4)` (number of successful `next()`), `event_weight` = 1.
- (A) NTuple: per pair `weight = eventWeight_file / N_files_total(slice)` (`NTP/PythiaAlgCoreT.c:1201`) =
  average over files. Correct ONLY if every file has exactly the hard-coded `nevents_per_file`
  (`kjob` lookup `:1168`), no file is missing (the skip at `:225` does NOT reduce the hard-coded
  `njobs_all_files_combined`), and both batches are hadded. `nevents_max` truncation does NOT renormalise.
  **Isospin 4:6:6:9 exists only as file counts per beam** — there is no isospin weight.
- **Units: (A) µb, central/AMI nb** (×1000). RDF uses the `weight` column as is.
- (B): either the external raw-counts mode (weight = 1, `nproc` per slice in `meta_tree_out`, slice weight
  applied downstream) or — since 2026-10-07 — the production-weighted external mode
  (`truth_external_sigma_eff_nb`, weight = σ·ε_filt/N_processed in nb at the NTuple stage; one run per slice
  over all its files; `docs/pythia_truth.md`). In both cases the slice numbers come from: w_slice = σ_slice · ε_filt / N_acc with σ = mean over jobs of
  `MetaData: cross-section (nb)` and ε_filt = Σacc/Σgen from `EvgenFilterSeq ... Filter Efficiency` in THIS
  sample's own `log.generate` (template: `JO/analysis/plot_test_sample_sanity.C:64-83`). **Never** use the central
  AMI files for a private sample: the AMI-rule failure mode (slice-dependent, cancels nowhere,
  `docs/ami_weights.md`) applies identically to slice weights from the wrong production.
- Slices combine by σ: Σ_slices w_slice·(counts). For a RATIO such as A_sig the global scale cancels but the
  RELATIVE slice weights do not — a wrong ε_filt in one slice biases A_sig.

## 4. Generator settings vs ATLAS central (for any new standalone run)

The (A) code's tune block = ATLAS `Pythia8_A14_NNPDF23LO_Common.py` verbatim (Tune:ee 7, Tune:pp 14 + A14
parameters) EXCEPT: `ColourReconnection:range = 1.71` missing (Pythia default 1.8); ATLAS base fragment
(`ParticleDecays:limitTau0 = on`, `tau0Max = 10` — K0S/Λ left to Geant4; SM masses) missing; **no EvtGen**. PDF:
A14-native is `NNPDF23_lo_as_0130_qed`, but **the pp24 central PDF evgen (803015–803020) overrides to
`NNPDF23_lo_as_0119_qed`** (`docs/ami_weights.md` L39). Seeds: `OFFSET + 10000·k + 1000·beam + job`; the 2023
batches used offsets 100000…900000 (900000 is already taken by 0429 and the 2025 tests); job ≥ 1000 collides
with the next beam. The submit script `rm -rf`s its output path — always a new path.

## 5. Generator filters truncate acceptance denominators

A generator muon filter removes every event whose muons fail it; an acceptance whose denominator extends
below the filter threshold is then biased (the event appears only if OTHER muons pass). Central: ≥2 muons
p_T > 3.7, |η| < 2.5; old private: p_T > 3.5. A filter is harmless iff EVERY denominator object passes it.
Note the filter acts on ANY muon pair in the event (any origin), so it is a superset of a single-b condition.

Athena `xAODDiLeptonMassFilter` (`JO/source/athena/Generators/GeneratorFilters/src/xAODDiLeptonMassFilter.cxx`):
stable leptons from the container, `pt >= MinPt`, `|eta| <= MaxEta`, `AllowSameCharge=False` → OS same-flavour,
`MinMass < m < MaxMass`, `pT(ll) > MinDilepPt`. `TruthMuons` (`xAODTruthParticleSlimmerMuon`) = all stable muons,
no p_T cut; it is created by `createxAODSlimmedContainer("TruthMuons", prefiltSeq)` (normally pulled in by
`xAODMultiMuonFilter_Common.py`).

## 6. Low pTHat: HardQCD is not regularised below a few GeV

`HardQCD:all` with `pTHatMin → 0` is cut only at `PhaseSpace:pTHatMinDiverge` = 1 GeV: σ(0–5) = **2.4 b**
at 5.36 TeV (measured 2026-10-07), ≫ σ_inel ≈ 70 mb — an interaction count, not an event cross section; the
"hardest-interaction" slicing is invalid there. The physical low-pT complement of HardQCD slices is
`SoftQCD:nonDiffractive` (σ ≈ 48 mb) with the hardest subcollision `info.pTHat() < pTHatMin(lowest slice)`.
The pp24 central PDF evgen starts at pTHat 8 (no 5–8 slice).

## 7. Processing pitfalls (code hook points)

- `PassCuts_PythiaCore` (`NTP/PythiaAlgCoreT.c:902-915`, truth muon |η| < 2.4, p_T > 4.5) is shared by truth,
  private AND the fullsim truth-pair loop (`PythiaFullSimExtras.c:510`, reco-eff denominators): any relaxed
  fiducial must be opt-in, default nominal, and refuse fullsim.
- The central-format skim stores ALL status-1 muon pairs, no p_T/η cut (`SkimCode TrigRates.cxx:2440-2500`;
  `truth_mupair_pass` is informational) → a relaxed truth fiducial is possible at the pair-list level; (A)'s
  `truth_mupair_*` holds only filter-passing muons.
- External mode: one kn per run into the pp beam slot; throws with private/fullsim; RDF has no reader for its
  outputs by default.
- `pipeline_pythia_truth.sh private` hadds the GLOB `pythia_private_sample/muon_pairs_pythia_*_no_data_resonance_cuts.root`
  — never write a new-mode output there. `plot_signal_acceptance_pythia.cxx` has a HARDCODED output dir (the
  nominal acceptance plots). Nonprivate diagnostic files (`ancestor_tracing_output_diagnostic_files/`) have an
  empty `batch_suffix` → shared across jobs.
- Missing central NTUP beam file = warning + skip (`PythiaAlgCoreT.c:329-333`) — silently drops an isospin
  component; `pp_only` still applies the 4/25 ratio to pp.

## 8. The no-muon-p_T-filter acceptance sample (2026-10-07)

See the tracking doc PD9 for decisions; the production record (path, slices, statistics, per-slice σ and ε_filt)
is in §8a below once produced.

### 8a. Production record
(filled when produced)
