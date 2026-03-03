# Pythia NTuple Processing Transformation Plan

**Overview:** Transform PythiaNTupleFirstPass from a private-only, fixed-structure code to support both private and ATLAS centrally produced samples, using CRTP + extras (mirroring Powheg/Data), dynamic kinematic ranges, event weighting for non-private, truth_parents branch handling, and per-kinematic-range batch processing.

---

## Current State Summary

- **PythiaNTupleFirstPass** inherits from `DimuonAnalysisBaseClass` (non-template), uses `MuonPairPythiaOld` (file not found; likely `MuonPairPythia`/`MuonPairPythiaTruth`), fixed `nKinRanges=5`, `PyTree` + `meta_tree` for private sample.
- **InitInputCentrProd** is declared but **not implemented** (would cause link error).
- **Private sample**: `/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/` with `job_dirs/kin_dirs/beam_dirs/pytree_N.root`, tree `PyTree`.
- **Non-private sample**: `/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/` with files `pythia_5p36TeV_{pp,pn,np,nn}_hQCD_DiMu_pTH{lo}_{hi}.TRUTH0.NTUP.root`, tree `HeavyIonD3PD`. `ami_info/` contains `ami_info_mc23_5p36TeV_Py8EG_A14_{beam}_hQCD_DiMu_pTH{lo}_{hi}.txt` with `crossSection` and `genFiltEff`.
- **Powheg/Data** use CRTP: `PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>`, `DimuonAlgCoreT`, `enable_and_bind`, `truth_parents` (vector of vectors) vs Pythia private `truth_mother1`/`truth_mother2`.

---

## Phase 1: Backups and Prerequisites

1. **Backup current code**
   - Copy `PythiaNTupleFirstPass.c` → `PythiaNTupleFirstPassOld.c`
   - Copy `PythiaNTupleFirstPass.h` → `PythiaNTupleFirstPassOld.h`
   - Backup current muon pair output files before testing (for cross-checks)

2. **Fix MuonPairPythiaOld reference**
   - Current header includes `MuonPairPythiaOld.h` (file not found). Replace with `MuonPairPythia.h` and use `MuonPairPythiaTruth` / `MuonPythiaTruth`. Ensure `FillMuonPair` and `PassCuts` use `truth_pt`, `truth_eta`, `truth_phi`, `truth_charge` (and `truth_same_sign`) to match `MuonMCTruthKinExtra` / `PairMCTruthExtras`.

---

## Phase 2: Dynamic nKinRanges and Data Structures

1. **Replace fixed nKinRanges with vectors**
   - Remove `static const int nKinRanges = 5`; use `std::vector`-based configuration.
   - Convert all fixed arrays to vectors:
     - `kin_dirs`, `beam_dirs`, `nevents`, `njobs`, `nevents_per_file`, `kinRanges`, `njobs_all_files_combined`
     - `muonPairOutTreeKinRange[nKinRanges][ParamsSet::nSigns]` → `std::vector<std::vector<TTree*>>` or `std::map`
     - `nentries_k0`..`nentries_k4` → `std::vector<long>` or `std::map`
   - For private: `nKinRanges = 5`, `kinRanges = {5, 10, 25, 60, 120, 3200}` (GeV).
   - For non-private: `nKinRanges = 6`, `kinRanges = {8, 14, 24, 40, 70, 125, 300}` (GeV).

---

## Phase 3: isPrivate Flow and Non-Private I/O

1. **isPrivate branching**
   - `isPrivate == true`: existing private logic (PyTree, meta_tree, job_dirs, etc.).
   - `isPrivate == false`: central-production logic with 6 kinematic ranges, 4 beam types.

2. **Non-private I/O structure**
   - Base path: `/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV` (5.36) or `pythia_5TeV` (5.02). Test only 5.36 for now.
   - Input files: `pythia_5p36TeV_{pp,pn,np,nn}_hQCD_DiMu_pTH{lo}_{hi}.TRUTH0.NTUP.root`.
   - Use 2D structure for trees and entry counts:
     - `std::vector<std::vector<TChain*>> evChains` or `std::map<std::pair<int,int>, TChain*>` for (kinematic range, beam type).
     - `std::vector<std::vector<Long64_t>> nentries_kn_beam` for #entries per (kn, beam).
     - `std::vector<Long64_t> nentries_kn_sum` for sum over beam types per kn.
   - Implement `InitInputCentrProd()`: build TChains from file pattern, fill `nentries_kn_beam` and `nentries_kn_sum`, set branch addresses for `HeavyIonD3PD`.

---

## Phase 4: Event Weighting (!isPrivate)

1. **Three-step event weight for non-private**
   - **(i) AMI weight**: For each (kn, beam), read `ami_info/ami_info_mc23_5p36TeV_Py8EG_A14_{beam}_hQCD_DiMu_pTH{lo}_{hi}.txt`, grep `crossSection` and `genFiltEff`. Event weight factor = `crossSection * genFiltEff`. Error if file missing.
   - **(ii) Beam-type correction**: Nominal ratios `pp:4/25, np:6/25, pn:6/25, nn:9/25` (map beam type → ratio). For each kn, compute actual ratio = `N_{beam}^{kn} / N_{sum}^{kn}`. Correction factor = `(actual_ratio / nominal_ratio)^{-1}`.
   - **(iii) Normalization**: Divide by `N_{sum}^{kn}`.
   - Final weight: `(crossSection * genFiltEff) * correction_factor / N_sum`.
   - Write `mpair->weight` for each muon pair in the event.

---

## Phase 5: Branch Differences (truth_mother1/2 vs truth_parents)

*(Truth particle branches and all access to them live in PythiaTruthExtras; see Phase 7.)*

1. **Branch address setup (in PythiaTruthExtras InitInputExtra)**
   - `isPrivate`: bind `truth_mother1`, `truth_mother2` (existing).
   - `!isPrivate`: bind `truth_parents` (std::vector<std::vector<int>>*).

2. **Truth origin tracing**
   - Introduce helper or inline logic at each use of `truth_mother1->at(IND)` / `truth_mother2->at(IND)`:
     - `isPrivate`: `mother1 = truth_mother1->at(IND)`, `mother2 = truth_mother2->at(IND)`.
     - `!isPrivate`: `mother1 = (truth_parents->at(IND).empty()) ? 0 : truth_parents->at(IND).at(0)` (with out_of_range check); `mother2` = last entry of `truth_parents->at(IND)`.
   - Affected locations (when moved into PythiaTruthExtras): UpdateCurParents, GluonHistoryTracking, termination check. Also check GetPtEtaPhiMFromBarcode and any `truth_id->at(bar)` usage: private uses barcode as index; non-private may use `truth_barcode` to map barcode → index.

3. **Branch compatibility check**
   - Inspect one private and one non-private ROOT file (e.g. `root -l file.root`, `tree->Print()`) to list branches. Document differences. Handle or flag any mismatches; clarify with user if needed.

---

## Phase 6: Per-Kinematic-Range Batching (!isPrivate)

1. **kn_batch for non-private**
   - Add `kn_batch` (or `kin_batch`) constructor parameter: process only kinematic range `kn_batch` (0..5 for 6 ranges).
   - For `!isPrivate`, load only trees for that kn, process events from that kn only.
   - Create `run_pythia_truth_5p36.sh` and `run_pythia_truth_5p36.sub` (queue 6 jobs, one per kn).
   - Defer isPrivate batching (batch_num × nKinRanges) to a later step if it adds risk.

---

## Phase 7: CRTP + Extras Refactor

**Core vs TruthExtras split (mirror Powheg):** History/origin tracing and generic truth-particle data live in **PythiaTruthExtras**, not in PythiaAlgCoreT.

1. **PythiaAlgCoreT (lean, like PowhegAlgCoreT)**
   - Create `PythiaAlgCoreT.h` (mirror `PowhegAlgCoreT.h`):
     - `template <class PairT, class MuonT, class Derived, class... Extras> class PythiaAlgCoreT : public DimuonAlgCoreT<PairT, MuonT, Derived>`
   - **In Core only:** isPrivate, E_CoM (5.36/5.02), kn_batch, py_dir, InitInput (chain building), muon-pair branches (truth_mupair_* or muon_pair_muon*), event-level branches used for weighting/selection (e.g. QHard, pTHat, mHat only if needed for core flow), meta_tree (private), EventWeights / weight-related members (non-private). No truth particle arrays (truth_id, truth_barcode, truth_parents/mothers). No ancestor-tracing or origin logic.

2. **PythiaTruthExtras (all history/origin + generic truth)**
   - Create `PythiaTruthExtras.h/.c` (mirror `PowhegTruthExtras`):
   - **In TruthExtras:** All generic truth particle branch pointers: `truth_id`, `truth_barcode`, `truth_status`, `truth_m`, `truth_pt`, `truth_eta`, `truth_phi`; and either `truth_mother1`/`truth_mother2` (isPrivate) or `truth_parents` (non-private). All ancestor/history logic and related state:
     - **Functions:** MuonPairAncestorTracing, SingleMuonAncestorTracing, UpdateCurParents, GluonHistoryTracking, FindHeavyQuarks, GetPtEtaPhiMFromBarcode, PrintHistory, ParentGrouping, HardAnalysisCategr, FillCategoryHistograms, HFMuonPairAnalysis (and any helpers they need).
     - **State:** m1_history, m2_history, m1_history_particle, m2_history_particle, cur_prt_bars/ids, m1_ancestor_is_incoming, m2_ancestor_is_incoming, m1_from_tau, m2_from_tau, resonance_ids, m1_resonance_barcode, m2_resonance_barcode, parent group / origin category labels, and Pythia-specific histograms (parent_groups, origin categories, QQ plots, etc.).
   - Truth origin tracing and mother/parent access (including isPrivate vs truth_parents handling) are implemented only inside TruthExtras.

3. **Placeholder extras**
   - Create `PythiaFullSimExtras.h/.c` and `PythiaFullSimOverlayExtras.h/.c` as placeholders (input paths blank).

4. **Analysis classes (template-based, not hard-coded)**
   - Pythia analysis must be a **template class** mirroring Powheg: `PairT` and `MuonT` are template parameters, not hard-coded `MuonPairPythiaTruth` / `MuonPythiaTruth`.
   - Example: `PythiaTruthAnalysis` : `PythiaAlgCoreT<MuonPairPythiaTruth, MuonPythiaTruth, PythiaTruthAnalysis, PythiaTruthExtras<MuonPairPythiaTruth, PythiaTruthAnalysis>>` — the concrete analysis class *instantiates* the template with those types; the core `PythiaAlgCoreT` is `template <class PairT, class MuonT, class Derived, class... Extras>`.
   - `PythiaFullSimAnalysisNoTruth`, `PythiaFullSimAnalysisWTruth`, etc. as placeholders with their respective PairT/MuonT.

5. **PythiaAnalysisClasses.h**
   - Include `PythiaAlgCoreT.c`, `PythiaTruthExtras.c`, `PythiaFullSimExtras.c`, `PythiaFullSimOverlayExtras.c`. Define all Pythia analysis classes (mirror `PowhegAnalysisClasses.h`).

---

## Phase 8: Scripts and Testing

1. **Scripts**
   - `run_pythia_truth_bb.sh` (analog of `run_powheg_truth_bb.sh`): `source ~/setup.sh`, `root -b -l`, `.L PythiaAnalysisClasses.h`, `PythiaTruthAnalysis pw(...); pw.Run();`
   - `run_pythia_truth_bb.sub`: queue jobs (for private: batch_num; for non-private: kn_batch).
   - `run_pythia_truth_5p36.sh` and `run_pythia_truth_5p36.sub` for non-private 5.36 TeV.

2. **Testing**
   - **Private**: `source ~/setup.sh`, `root -b -l`, `.L PythiaAnalysisClasses.h`, `PythiaTruthAnalysis pw(1, true); pw.nevents_max = 1000; pw.Run();` (or equivalent for batch_num=1, isPrivate=true). Do NOT compile; use ROOT `.L` only.
   - **Non-private**: Same pattern with `isPrivate=false`, `kn_batch=0`, `nevents_max=1000`.
   - Compare output with backed-up muon pair files for private sample.

---

## Clarifications and Risks

| Item                         | Note                                                                                                                                                                                                                                                         |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **MuonPairPythiaOld**        | File missing; use `MuonPairPythiaTruth` and ensure `pt`/`eta`/`phi`/`charge` → `truth_pt`/`truth_eta`/`truth_phi`/`truth_charge` in FillMuonPair/PassCuts.                                                                                                   |
| **Private index vs barcode** | Private code uses `truth_pt->at(bar)` (barcode as index). Non-private typically uses `truth_barcode` to map. Verify indexing in both samples.                                                                                                                |
| **ProcessData loop**         | Private loops over kn, then events; uses `evChain->GetEntry(jevent)`, `metaChain->GetEntry(kjob)`. Non-private: no meta_tree; loop over (kn, beam) or use a single chain per kn built from 4 beam files. Need clear event-ordering strategy for non-private. |
| **ami_info grep**            | Parse `crossSection` and `genFiltEff` from AMI txt (key: value format). Implement robust parsing.                                                                                                                                                          |
| **Branch check**             | Run `tree->Print()` on one file per sample to list branches; document and handle differences before full testing.                                                                                                                                            |

---

## Suggested Implementation Order

1. Backups (Phase 1)
2. Fix MuonPairPythiaOld and vectorize nKinRanges (Phases 1–2)
3. Implement InitInputCentrProd and non-private I/O (Phase 3)
4. truth_mother1/2 vs truth_parents handling (Phase 5)
5. Event weighting (Phase 4)
6. CRTP + extras refactor (Phase 7)
7. Per-kn batching and scripts (Phases 6, 8)
8. Testing (Phase 8)
