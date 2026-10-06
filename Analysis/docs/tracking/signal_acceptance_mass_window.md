# Signal acceptance A_sig in (pair pT, η^pair): truth mass window + single-muon p_T/η cuts

(File name kept for continuity; since 2026-10-06 A_sig ALSO contains the single-muon cuts — PD4.)

Mode: **IMPLEMENTATION**. Opened 2026-10-05.

## Objective

Redefine and rebuild the truth signal acceptance A_sig (Pythia + POWHEG truth) per the user's
measurement definition (2026-10-05/06), removing the gap cuts it double counts with ε_reco
(`pair_reco_eff_gap_acceptance.md`), and plot its pair-pT AND pair-η dependence. **2026-10-06
(user decision, PD4):** A_sig must ALSO correct the single-muon truth cuts p_T > 4.5 GeV and
|η| < 2.4, which requires a NEW private Pythia truth sample WITHOUT a muon-p_T generator filter.

## Physics Procedure

### 1. Motivation
The final result is an UNFOLDED cross section of truth single-b OS dimuon pairs,
dσ/d(p_T^pair, η^pair), measured in truth pair p_T > 9 GeV and |truth η^pair| < 2.2 (user,
2026-10-05). The detector chain (ε_trig, ε_reco, unfolding) turns the reconstructed yield into the
number of truth pairs in the region T_reco-def = ε_reco's truth denominator. A_sig then converts
that to the number of truth single-b pairs in the TARGET phase space of each (p_T^pair, η^pair)
cell. A_sig must therefore contain exactly the truth cuts that (i) a pair inside a measured cell
can fail and (ii) ε_reco's denominator does not already correct — no more (double count), no less
(gap). The decisions below (§5) fix that list.

### 2. Top-level equation
    dσ/dX (truth X, after unfolding) = N_unfolded(X) / (L · A_sig(X)),   X = (p_T^pair, η^pair)

    TARGET (decided 2026-10-06, PD4):
    A_sig(cell) = N[ T_target AND both μ truth p_T > 4.5 GeV & |η| < 2.4 AND m_μμ^truth ∈ window ]
                / N[ T_target ]
    T_target = truth single-b OS pair (`from_same_b`), truth p_T^pair > 9 GeV, |truth η^pair| < 2.2,
               both muons |truth η| < 2.5 (ID coverage, PD4), NO muon-p_T cut, ALL masses.
    The numerator region is exactly ε_reco's truth denominator (+ the measured range), so
    A_sig · ε_reco = N_reco-selected / N[T_target]: no double count, no gap.

    CURRENTLY IMPLEMENTED (interim, 2026-10-06): only the mass-window factor —
    A_sig^mass = N[T_meas AND window] / N[T_meas], T_meas = T_target with the muon cuts
    p_T > 4.5 & |η| < 2.4 already applied (NTuple-stage truth fiducial). Plots/numbers R1–R4 are
    this interim quantity. The muon-cut factor needs the new private sample (§3e).
Binned in TRUTH (p_T^pair, η^pair): `pT_bins_150` (16 log bins 9–150 GeV) × 44 η^pair bins over
[−2.2, 2.2] (PD7), MC-event weighted.

### 3. Step-by-step method
(a) RDF truth fill (`RDFBasedHistFillingPythiaTruth.cxx`, `RDFBasedHistFillingPowhegTruth.cxx`,
`FillHistogramsSignalAcceptance`): den = `SignalAcceptanceDefinition::DenomExpr()` =
`from_same_b && pair p_T > 9 && |η^pair| < 2.2`; num = `NumExpr()` = den && truth mass window.
TNamed marker; the plotter refuses files without the current marker.
(b) Plots: 2D map; A_sig vs p_T^pair in coarse η^pair panels; A_sig vs η^pair in coarse p_T^pair
panels (`pair_pt_coarse_bins`); cutflow STARTING from the measured range with only the mass cuts
(its last cumulative bin = integrated A_sig).
(c) POWHEG: the truth NTuple must be the `_no_data_resonance_cuts` one (PD8).
(d) A_sig is applied LAST, to the unfolded yield (not yet wired: unfolding is a placeholder).
(e) **Muon-cut factor (TO IMPLEMENT, PD4/PD9):** new PRIVATE Pythia truth sample with NO muon-p_T
generator filter → truth NTuple processing with a NEW opt-in flag that replaces the PassCuts truth
fiducial (p_T > 4.5, |η| < 2.4) by |η| < 2.5 only (distinct output suffix; nominal untouched;
everything else identical — provenance rule) → A_sig producer with DenomExpr = T_target and
NumExpr = + muon cuts + mass window; the cutflow then lists (after "no cut" = T_target) the muon
cuts and the mass cuts. NOT implemented yet (user, 2026-10-06: "no code implement yet").

### 4. Negative constraints
1. Do NOT include any detector-gap cut (single-muon q·η windows, the reco |η^pair| gap) in A_sig.
2. Do NOT acceptance-correct p_T^pair > 9 or |η^pair| < 2.2 — apply them to num AND den, and never
   show them as an acceptance loss (cutflow starts from the measured range).
3. Do NOT compute the muon-cut factor from the central (3.7 GeV filter) or old private (3.5 GeV
   filter) Pythia samples — their denominators are truncated by the generator filter (PD9). Do NOT
   relax the NOMINAL NTuple truth cuts; the no-muon-p_T-cut truth processing is a NEW opt-in mode.
4. Do NOT compute A_sig from a resonance-vetoed NTuple (POWHEG nominal truth NTuple).
5. Do NOT change the crossx binning; the 44-bin η axis is an opt-in suffixed variant for A_sig only.

### 5. Physics Decisions (what A_sig includes / excludes, with reasoning)

**PD1 — EXCLUDED: truth p_T^pair > 9 GeV and |truth η^pair| < 2.2 (they define the measured range).**
(user, 2026-10-05.) Reasoning: these are cuts on the very variables the cross section is
differential in. Every reported cell lies entirely inside them, so no signal pair of any reported
cell is removed by them — there is nothing to correct. "Correcting" for them would mean quoting the
cross section in cells we do not measure (p_T^pair < 9, |η^pair| > 2.2): a pure model extrapolation
into unprobed phase space, which would also land in the axis under/overflow. This is universal HF
practice — the measured rapidity/p_T range is the result's definition, never acceptance-corrected:
CMS non-prompt D⁰ arXiv:2212.01636 §3 p.5 (correction = reco / generated "in a certain kinematic
range", i.e. the |y| < 1 analysis bin), ATLAS quarkonium arXiv:1709.03089 §5.2 pp.5–6 (A as a
function of the parent p_T, y — the p_T/y bins themselves are not corrected). Edge effect: truth
outside / reco inside (e.g. truth |η^pair| ≥ 2.2 reconstructed inside, 0.02 % of the numerator,
`pair_reco_eff_gap_acceptance.md` R1) is bin migration → the response matrix (unfolding), not A_sig.

**PD2 — EXCLUDED: detector-gap cuts (single-muon q·η windows ×2, reco |η^pair| gap).**
(user, 2026-09-17, `pair_reco_eff_gap_acceptance.md`.) Reasoning: they are imposed because the
detector gaps make the muon reconstruction/trigger efficiency low and rapidly varying there — a
DETECTOR effect, applied only to reconstructed muons. A truth signal pair lost to a gap is lost
because of the detector, i.e. an inefficiency, carried by ε_reco (gap cuts on the reco numerator
only; factor 0.817 pp24). Having them also in A_sig double counts: the pre-2026-10-05 A_sig did, and
showed dips to ≈0.45–0.53 in the gap η^pair cells that now vanish.

**PD3 — INCLUDED: the truth mass window 1.08 < m_μμ < 2.9 GeV ("extrapolate all mass").**
(user, 2026-10-05.) Reasoning: the window is a signal-SELECTION cut (removes the φ region below and
the J/ψ region above, where resonances dominate the data), not an axis of the observable and not a
detector effect. Inside every measured (p_T^pair, η^pair) cell there are truth single-b pairs outside
the window; ε_reco's truth denominator contains the window, so nothing else corrects for them. A_sig
answers "in each 2D cell, what fraction of truth single-b pairs does the mass window keep". It is
model-dependent through the b-hadron decay kinematics and branching fractions (Pythia: 15.9 % below
1.08 GeV incl. φ 2.4 %; 9.2 % above 2.9 GeV incl. B → J/ψ X → μμ 3.1 %, which therefore counts as
single-b signal) → a generator/decay-table systematic (Pythia vs POWHEG+Pythia8, PD8).

**PD4 — INCLUDED (DECIDED, user 2026-10-06): single-muon truth p_T > 4.5 GeV and |η| < 2.4.**
(Supersedes the 2026-10-05 "TO BE CLARIFIED / keep fiducial" entry; the user corrected that the
filter limits do NOT bound the extrapolation.) Reasoning: the muon cuts are not axes of the
observable. A truth single-b pair can sit inside a measured (p_T^pair, η^pair) cell while one muon has
p_T < 4.5 GeV (the b → c → μ cascade muon is typically soft) or |η| > 2.4; such a pair is NOT in
ε_reco's truth denominator, so — by the same logic as PD3 — only A_sig can correct for it.
- **Target muon phase space: NO muon-p_T cut, |truth η_μ| < 2.5** (user Remark 2: the inner detector
  covers only |η| < 2.5 and an ID track is required for a combined muon, so |η_μ| < 2.5 stays part
  of the target definition; A_sig extrapolates |η| 2.4 → 2.5 and p_T 4.5 GeV → 0).
- **Model dependence** (to be quantified as a systematic): the soft-muon extrapolation rests on the
  MC b-hadron decay kinematics (b → c → μ cascade spectrum); expected to be a small effect overall
  (user), but it must be measured, not assumed.
- Literature analogue: ATLAS quarkonium arXiv:1709.03089 §5.2 pp.5–6 — A = probability that both
  decay muons fall in p_T > 4 GeV, |η| < 2.4, as a function of the parent p_T, y (single-muon cuts
  inside A; parent kinematic range = the measured range).

**PD9 — Sample for the muon-cut factor: a NEW PRIVATE Pythia truth sample with NO muon-p_T filter
(user, 2026-10-06).** Reasoning: a generator filter truncates the denominator — a single-b pair
whose soft muon is below the filter threshold appears only if two OTHER muons pass the filter, so
the sub-threshold part of T_target is missing/biased. The ATLAS central samples (DSIDs 802758–802781,
e8599, and the pp24 PDF slices) filter at muon p_T > 3.7 GeV (`xAODMultiMuonFilter`, |η| < 2.5);
the previous private Pythia sample (`docs/pythia_truth.md` mode `private`) filters at p_T > 3.5 GeV.
Neither can give the 4.5 → 0 GeV extrapolation. User Remark 1: ATLAS central MC is the default for
nearly every task; a private sample is acceptable ONLY for this purpose (extending below the
central filter for the acceptance), and the effect may be small. PROPOSED requirements for the
private sample (Claude's proposal, NOT yet confirmed by the user — confirm before generating): same
physics configuration as the central
production (Pythia 8 A14 NNPDF2.3LO, HardQCD:all, pT-hat slices, √s = 5.36 TeV, pp beam) so that
A_sig's mass-window part can be cross-checked against the central sample above 4.5 GeV; generator
filter at most |η_μ| < 2.5 with no p_T cut (or none); per-slice σ·ε weights taken from THIS
sample's own generator output — never another sample's (AMI rule spirit, `docs/ami_weights.md`:
slice-dependent weights do not cancel in A_sig); enough statistics at high p_T^pair given the
unfiltered rate.

**PD5 — Applied to the UNFOLDED yield, in TRUTH variables.** Reasoning: A_sig is a pure
generator-level ratio (both legs truth) defined in truth (p_T^pair, η^pair); detector migration
belongs to the response matrix. Literature: A at generator level and separate from the detector
correction C — ATLAS W/Z arXiv:1603.09222 §2 Eq.(1) p.2. Counter-example (allowed when migration is
small): ATLAS quarkonium arXiv:1709.03089 §5.2 applies generator-level A per candidate at reco p_T
plus an FSR correction.

**PD6 — OS only; ΔR not cut.** A_sig uses the OS truth tree (`df_op`, `from_same_b`); no ΔR cut exists
in the signal selection (removed 2026-06-22), so none appears in A_sig.

**PD7 — η^pair axis = 44 bins over [−2.2, 2.2] (user, 2026-10-05).** Reasoning: the measured range
ends at 2.2, so bins beyond it are structurally empty; width 0.1 keeps every edge identical to the
crossx fine-axis edges inside the range (A_sig cell k ↔ crossx `pair_eta_crossx` cell k+2) and every
coarse-panel boundary an edge. Opt-in suffixed ParamsSet variant `pair_eta_measured_bins`
(`N_PAIR_ETA_MEASURED_BINS`, `PAIR_ETA_MEASURED_MAX`, checked == `pair_eta_fiducial_max`); the crossx
axis is unchanged.

**PD8 — POWHEG from the bb NO-resonance-cut truth NTuple.** Reasoning: A_sig's denominator is ALL
masses, but the nominal POWHEG truth NTuple drops every pair with a resonance-tagged muon (R2),
removing precisely the regions the mass window cuts → biased-high A_sig (0.965). New flag
`PowhegAlgCoreT::turn_data_resonance_cuts_on` (default true = nominal; false keeps the tagged pairs,
mirroring Pythia's `turn_data_resonance_cuts_on = false`; suffix `_no_data_resonance_cuts`). bb only:
the cc sample has exactly ZERO single-b OS pairs in the measured range (checked 2026-10-06: N = 0),
so cc cannot enter A_sig (checked on the nominal cc NTuple: `from_same_b` = 0 at ALL masses, 2 441 289 OS pairs; a no-res cc file could add no single-b pairs).

## Design Decisions
- D1 — single source of the definition: `Utilities/SignalAcceptanceDefinition.h` (DenomExpr/NumExpr
  from ParamsSet helpers + marker), used by both producers, the cutflow and the plotter guard.
- D2 — η-dependence view projects the SAME 2D histograms in `pair_pt_coarse_bins` panels (CLAUDE.md
  Binnings rule 2; coarse edges checked to be fine-axis edges).
- D3 — cutflow (round 2, user): starts from the measured range (first bin labelled "no cut", user
  2026-10-06: the measured range is not an acceptance cut) and lists ONLY the mass cuts → last cumulative bin = A_sig (PD1).
- D4 — POWHEG no-res-cut path is bb-only end to end: NTuple `muon_pairs_powheg_bb_truth_no_data_resonance_cuts.root`
  → RDF `RunRDFBasedHistFillingPowhegTruth(false)` (`run_rdf_powheg_truth.sh nores`) →
  `histograms_powheg_truth_bb_no_data_resonance_cuts.root` (PD8). Nominal POWHEG files untouched.

## Implementation Plan
| # | Step | Physics § | Review | Status |
|---|---|---|---|---|
| 1 | Producers: new num/den + marker (Pythia + POWHEG truth) | §3a, §4.1–2 | /review-analysis-code | DONE (APPROVED it.2) |
| 2 | Plotter: marker guard, η-dependence view; cutflow kCuts | §3b | /review-analysis-code (plots checked C1/C2) | DONE |
| 3 | Compile; rerun Pythia 5.36 truth RDF + POWHEG truth RDF; replot | §3a–b | /review-plot | Pythia DONE; POWHEG blocked (R2) |
| 4 | Docs (analysis_overview, template-fit §3f pointer, signal_selection_change_impact, powheg.md), INDEX, commit | — | — | DONE |
| 5 | Round 2: cutflow from measured range, mass cuts only | PD1, D3 | /review-analysis-code | DONE |
| 6 | Round 2: 44-bin η^pair axis (ParamsSet variant) in both producers + marker bump; refill Pythia | PD7 | /review-analysis-code | DONE (Pythia 0.7493 unchanged, axis 44 bins) |
| 7 | Round 2: POWHEG NTuple no-res flag + Condor bb rerun (clusters 9, 10, 11) + hadd + RDF nores + plots | PD8 | /review-analysis-code | DONE (R3), review running |
| 8 | Round 2: Physics Decisions section (this doc) | §5 | — | DONE |
| 9 | Private no-muon-p_T-filter Pythia truth sample: spec (confirm with user) + production | PD9 | /review-pipeline | TODO |
| 10 | Opt-in truth-NTuple flag: truth fiducial → |η_μ|<2.5 only (distinct suffix) | PD4, §3e | /review-analysis-code | TODO |
| 11 | A_sig with T_target denominator; muon cuts + mass window in numerator; marker bump; plots | §2, PD4 | /review-analysis-code + /review-plot | TODO |

## Progress Log

- **Steps 1–2 code (2026-10-05).** NEW `Utilities/SignalAcceptanceDefinition.h` (Key/Value marker +
  `DenomExpr()` = `from_same_b && SignalPairPtCutExpr(truth_pair_pt) && PairFiducialEtaCutExpr(truth_pair_eta)`,
  `NumExpr()` = Denom && `SignalMinvCutExpr(truth_minv)`). `RDFBasedHistFillingPythiaTruth.cxx` /
  `RDFBasedHistFillingPowhegTruth.cxx` `FillHistogramsSignalAcceptance` use them (gap cuts removed);
  marker written in `WriteOutputExtra` (new override for POWHEG truth, `RDFBasedHistFillingPowheg.h`).
  `SignalAcceptancePlotter.cxx`: marker guard in `Init`, 2D map and η view restricted to |η^pair|<2.2,
  binomial ("B") ratio errors, NEW `DrawAcceptanceEtaByPt` (A_sig vs η^pair in the 8
  `pair_pt_coarse_bins` panels, edge-alignment checked → `*_sig_accept_pair_eta_in_pt_subplots.png`,
  default pT_bins_150 axis only). `plot_sig_accept_cutflow_above_60GeV.cxx`: kCuts = pair pT>9,
  |η^pair|<2.2, m>1.08, m<2.9 (no gap); prints integrated A_sig. Verified both truth NTuples apply
  truth muon pT>4.5 & |η|<2.4 (`PythiaAlgCoreT.c:906-911`, `PowhegAlgCoreT.c:252-257`; NTuples dated
  2026-09-09/10/11, post 4.5-GeV adoption). Backups:
  `.../pythia_5p36TeV/backup/histograms_pythia_5p36TeV_no_data_resonance_cuts.bak_20261005_sigacc.root`,
  `.../powheg_full_sample/backup/histograms_powheg_truth.bak_20261005_sigacc.root`.

- **Review (2026-10-05):** `/review-analysis-code` APPROVED at iteration 2
  (`.claude/logs/review-analysis-code-20261005-225359-signal-acceptance-mass-window.md`). It.1 WARNING:
  `SetRangeUser(±2.2f)` displayed ±2.3 (float rounding) → fixed with `SetMeasuredEtaRange` (bin
  index); cutflow index → `kNMeasuredRangeCuts`. All numbers independently reproduced. Per coarse
  pair-pT bin A_sig = 0.758, 0.750, 0.738, 0.723, 0.710, 0.711, 0.714, 0.721. C3: RUN2-CROSSCHECK
  UNVERIFIED (no comparable Run 2 quantity). INFO left: zero binomial error at A=0/1 edge cells;
  stored `h2d_sig_accept_pt_eta` uses plain Divide (not plotted with errors).

- **Round 2 review (2026-10-06):** `/review-analysis-code` PASS at iteration 1
  (`.claude/logs/review-analysis-code-*-signal-acceptance-round2.md`). Provenance verified pair-by-pair on
  part6: nores = nominal + exactly the vetoed pairs (identical weights/minv, no duplicates, cut-acceptance
  hists identical up to the resonance step). All numbers MATCH; A_sig edges = crossx edges to 4e-16. INFO:
  POWHEG has negative weights (0.09 % of pairs) → a few high-pT cells with num > den (one plotted cell
  A = 1.26 at 126–150 GeV) and N_eff/N = 0.17 → POWHEG A_sig meaningful only for p_T^pair ≲ 50 GeV; the
  nores POWHEG NTuple stores no veto tag (cannot rebuild nominal from it); `pair_eta_measured_bins` vector
  unused (axis built from the constants).
- **Round 2 code (2026-10-06).** Cutflow rewritten (D3). `ParamsSet.h`: `N_PAIR_ETA_MEASURED_BINS = 44`,
  `PAIR_ETA_MEASURED_MAX = 2.2` (double, checked == `pair_eta_fiducial_max`), `pair_eta_measured_bins`;
  both producers book `h2d_sig_accept_*` on it; marker Value now names the axis (old 48-bin files refused).
  `PowhegAlgCoreT.{h,c}`: `turn_data_resonance_cuts_on` (+ public setter), suffix
  `_no_data_resonance_cuts`, throws on fullsim; new `run_powheg_truth_bb_no_res_cuts.{sh,sub}`.
  `RDFBasedHistFillingPowheg.{h,cxx}` + Truth: `with_data_resonance_cuts` (bb-only nores input + distinct
  output); `run_rdf_powheg_truth.sh nores`. `plot_signal_acceptance_powheg.cxx` default input → the nores
  output. Local 20k-event test (batch 1): nores ⊇ nominal (OS 7844 vs 6174 pairs; OS J/ψ [3.0,3.2] 814
  vs 0; m<1 385 vs 0; SS 2672 vs 2589), meta `nentries_before_cuts` identical (20000).
  ⚠ INCIDENT: a first test call hit a protected setter and ran in NOMINAL mode, overwriting the
  intermediate `bb_evgen_truth_full_sample/muon_pairs_powheg_bb_truth_part1.root` (+ its hists file)
  with a 20k-event version (no backup existed). The hadded `muon_pairs_powheg_bb_truth.root` read by
  all consumers is untouched. Regenerated by a nominal Condor job (cluster 10, batch 1).
  Condor: cluster 9 = 6 nores bb batches.
- **Condor (2026-10-06).** Cluster 10 (nominal part1 regen) done: entries part1 = combined − parts2–6
  EXACTLY (SS 127551, OS 304683, meta 1; diff 0) → incident fully repaired AND the nominal code path
  is verified unchanged. Cluster 9: batches 1,3,5,6 done in ~10–15 min; 9.1 and 9.3 hung on the /pnfs
  input read (4 s CPU in 2 h, no output after 00:01) → condor_hold/release at 02:00.

## Results & Observations

### R1 — Pythia 8.3 (5.36 TeV truth, 2026-10-05 refill)
Integrated A_sig (measured range, all pT) = **0.7493** (num 6.11509 / den 8.16069 in
`h2d_sig_accept_{num,denom}_pt_eta`; cutflow 0.749335), pT>60: 0.7064. Single-b truth mass in the
measured range: m<1.08 → 15.9 % (φ peak [0.98,1.06] 2.4 %), m>2.9 → 9.2 % (J/ψ [3.0,3.2] 3.1 %:
B→J/ψX→μμ counts as single-b, user: "extrapolate all mass"). Shape: ≈0.76 at 9 GeV falling smoothly
to ≈0.71 at ~100 GeV; FLAT in η^pair (old gap dips 0.45/0.53 gone — they were the double-counted
gap acceptance). Plots: `plots/single_b_analysis/pythia/pythia_sig_accept_{pair_pt_vs_pair_eta,
pair_pt_in_eta_subplots,pair_eta_in_pt_subplots,above_60GeV,above_60GeV_accum}.png` (+ `pythia_pt_120/`).

### R2 — ★ POWHEG A_sig INVALID: the POWHEG truth NTuple is resonance-vetoed
POWHEG refill integrated A_sig = 0.9654 (num 2691.44 / den 2787.97) vs Pythia 0.749. Root cause:
`PowhegAlgCoreT.c:327-351` applies `ResonanceTagging` UNCONDITIONALLY (no switch), so
`muon_pairs_powheg_bb_truth.root` OS single-b pairs have ZERO weight for m<1.0 and in [2.9,3.3],
[3.5,3.8] (measured: m<1.08 0.8 %, m>2.9 2.7 %, [3.0,3.2] 0.0 %) — exactly the regions the mass
window removes, so the all-mass denominator is missing them. Pythia reads the `_no_data_resonance_cuts`
NTuple. Fix (provenance rule 1): add a no-resonance-cut flag to POWHEG NTuple processing (distinct
suffix) + rerun bb truth on Condor. POWHEG A_sig plots NOT regenerated (old ones are stale: old gap
definition). AWAITING USER.

### R3 — POWHEG A_sig from the bb no-resonance-cut NTuple (2026-10-06)
NTuple: 6 Condor batches (cluster 9; batches 2/4 hung on /pnfs NFS reads of e.g. `mc_truth_bb_06.root`
— LGD pool unreachable over NFS while xrootd works — rerun as cluster 11 with the existing
`POWHEG_USE_XROOTD=1` switch; physics identical), hadded →
`bb_evgen_truth_full_sample/muon_pairs_powheg_bb_truth_no_data_resonance_cuts.root`. Checks: N_gen
4 999 000 = nominal; OS 1 958 925 vs nominal 1 525 248, SS 661 775 vs 640 565 (superset); single-b in
1.1<m<2.85 weight +0.28 % (nominal also drops pairs sharing a muon with a resonance-mass pair).
RDF (`run_rdf_powheg_truth.sh nores`) → `histograms_powheg_truth_bb_no_data_resonance_cuts.root`.
**POWHEG integrated A_sig = 0.513** (num 2699.0 / den 5260.65; pT>60: 0.541), flat in η^pair and
p_T^pair within stats. Plots `plots/single_b_analysis/powheg/powheg_sig_accept_*.png`.

### R4 — ★ A_sig is dominated by the B → J/ψ X → μμ fraction, which differs 10× between generators
Single-b OS mass fractions in the measured range: J/ψ [3.0,3.2) **31.4 % POWHEG vs 3.1 % Pythia**;
[1.08,2.9) 51.3 % vs 74.9 %; m<1.08 11.3 % vs 15.9 %. 94 % of POWHEG J/ψ-window pairs have both
muons `parent_group == b_to_c` (J/ψ = cc̄ from the b) → genuine B → J/ψ decays in the sample, not a
tag artefact. Consequence: with "extrapolate all mass" (PD3) A_sig is ~32 % generator-dependent,
driven almost entirely by whether / how often B → J/ψ X → μμ is counted as single-b signal.
**Physics question for the user:** (a) which generator's B → J/ψ rate is right (check against PDG
BR(B → J/ψ X) ≈ 1.1 % × BR(J/ψ → μμ) ≈ 6 %, and the kinematic enhancement of two hard J/ψ muons vs
the soft cascade muon); (b) whether B → J/ψ X → μμ should be part of the signal at all (PD3 counts
it, since it is "single-b"). AWAITING USER.

## Remaining Work
- **PD4/PD9 implementation (NEXT, not started — user said no code yet):** (1) define + produce the
  private no-muon-p_T-filter Pythia truth sample (open items in PD9); (2) opt-in truth-NTuple flag
  (truth fiducial → |η_μ| < 2.5 only, distinct suffix); (3) A_sig producer/def with T_target
  denominator and muon cuts + mass window in the numerator (marker bump); (4) cross-check vs the
  central sample's mass-window factor; (5) plots + review.
- R4 B → J/ψ question (user).
- POWHEG: its truth sample also has a generator muon filter (filter_effcy 0.003 bb) → POWHEG cannot
  give the muon-cut factor either; POWHEG A_sig^mass stays a cross-check of the mass-window part only.
- Wire A_sig into the crossx after unfolding exists.

## Latest Stage
**2026-10-06 — DECISION RECORDED, NO CODE YET.** PD4 decided: A_sig includes the single-muon truth
cuts p_T > 4.5 GeV and |η| < 2.4; target muon phase space = no p_T cut, |η_μ| < 2.5. PD9: requires a
NEW private Pythia truth sample with no muon-p_T filter (central filters 3.7 GeV, old private 3.5 GeV).
Current outputs (R1–R4: Pythia 0.749, POWHEG 0.513) are the INTERIM mass-window-only factor.
Next step when the user says go: Remaining Work "PD4/PD9 implementation" (1)–(5), starting with the
private-sample specification (PD9 open items) — confirm the spec with the user before generating.
Also OPEN: R4 (B → J/ψ fraction 31 % POWHEG vs 3 % Pythia).
