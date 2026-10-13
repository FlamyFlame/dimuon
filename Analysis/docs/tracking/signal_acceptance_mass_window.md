# Signal acceptance A_sig in (pair pT, η^pair): truth mass window + single-muon p_T/η cuts

(File name kept for continuity; since 2026-10-06 A_sig ALSO contains the single-muon cuts — PD4.)

Mode: **IMPLEMENTATION**. Opened 2026-10-05.

## Objective

Redefine and rebuild the truth signal acceptance A_sig (Pythia + POWHEG truth) per the user's
measurement definition (2026-10-05/06), removing the gap cuts it double counts with ε_reco
(`pair_reco_eff_gap_acceptance.md`), and plot its pair-pT AND pair-η dependence. **2026-10-06
(user decision, PD4):** A_sig must ALSO correct the single-muon truth cuts p_T > 4.5 GeV and
|η| < 2.4, which requires a NEW private Pythia truth sample WITHOUT a muon-p_T generator filter.

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = (1) reference doc on private-sample generation/processing subtleties written
  (`Analysis/docs/private_pythia_sample_notes.md`), linked here, and summarised to the user;
  (2) NEW private pp-beam NNPDF2.3LO Pythia truth sample with NO muon-pT filter (filter: one OS
  muon pair, both |eta_mu|<2.5, pair pT>8 GeV, no mass cut), all pTHat slices incl. a lowest
  slice (0-5 unless filter eff < 1e-7, then 5-8 or suitable), statistics chosen and justified,
  generated (grid or Condor), truth-NTuple-processed with an opt-in flag (truth muon cut |eta|<2.5
  only, distinct suffix; no main/previous-private output overwritten); (3) sanity plots of the new
  sample (pair pT all slices summed, pair eta, origin-categorised mass 0-4 GeV, per-slice checks)
  under a distinct plot dir; (4) A_sig with T_target denominator incl. muon pT>4.5 & |eta|<2.4 +
  mass window; plots of EACH single-muon cut factor and the total vs pair pT (in pair-eta panels)
  and vs pair eta (in pair-pT panels); /review-analysis-code + /review-plot PASS; docs + commit;
  then remind the user of R4.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

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
NumExpr = + muon cuts + mass window. In place of the cutflow, the target plots show EACH factor
separately (muon p_T, muon η, both muon cuts, mass window on T_target) and the total, vs p_T^pair in
η^pair panels and vs η^pair in p_T^pair panels (`plot_signal_acceptance_pythia_target.cxx`); a target
cutflow can be added later if wanted. NOT implemented yet (user, 2026-10-06: "no code implement yet").

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
central filter for the acceptance), and the effect may be small. CONFIRMED by the user 2026-10-07: pp beam ONLY, NNPDF2.3LO (match the pp24 central PDF production);
filter = one OS muon pair with both |η_μ| < 2.5, NO muon-p_T cut, pair p_T > 8 GeV (1 GeV margin
below the 9 GeV measured-range edge), no mass cut — every T_target pair passes it, so the
denominator is unbiased; all pTHat slices plus a lowest slice pTHat 0–5 GeV (user: if its filter
efficiency < 1e-7 use 5–8 GeV or a suitable range as the lowest). Original proposal text: same
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

- D5 (2026-10-07) — GENERATOR = Athena Gen_tf (private grid production), NOT the standalone
  `PythiaSimulationCode` program. Reason (physics): the user fixed "match the pp24 central PDF
  production"; that production decays b/c hadrons with EvtGen (`Pythia8_A14_NNPDF23LO_EvtGen_Common.py`),
  which sets the semileptonic B/D branching fractions and muon spectra — precisely what the soft-muon
  extrapolation of A_sig rests on (PD4 model dependence). The standalone program has no EvtGen and
  also lacks `ColourReconnection:range=1.71` and the ATLAS base fragment. Gen_tf with the central A14 +
  `NNPDF23_lo_as_0119_qed` + EvtGen fragments reproduces the central configuration by construction,
  and its output (EVNT → TRUTH0 → SkimCode NTUP) is in the CENTRAL format, so the nominal
  (nonprivate) truth navigation (HepMC barcodes, `truth_parents`, status) applies unchanged — none
  of the standalone-format conventions (0-based index barcode, Pythia status, µb weights) enter.
  Chain = the tested `PythiaCentralProdJO/evgen_high_pT_additional/grid/` chain. JO:
  `PythiaCentralProdJO/evgen_nomupt_acceptance/run/Py8Control_hQCD_OSdimu_nomupt.py`. Local test
  2026-10-07 (pTH70_125, 20 ev): exit 0, σ = 1184 nb, ε_filt = 2.0e-3 (standalone estimate 2.4e-3), ~40 ms/trial.
  Reference: `docs/private_pythia_sample_notes.md`.
- D6 (2026-10-07) — slice weights: w_slice = σ_slice·ε_filt/N_processed [nb], σ = generated-event-weighted mean of the
  jobs' `MetaData: cross-section (nb)`, ε_filt = Σacc/Σgen over ALL jobs (one NTuple run per slice) of `EvgenFilterSeq Filter Efficiency`, both from
  THIS production's own `log.generate` (AMI-rule spirit). Applied at the NTuple stage through a new
  external-mode option (alternative to raw counts) so the per-pair `weight` is in nb exactly as in
  the nonprivate path, and the hadded file is RDF-ready.
- D7 (2026-10-07) — target truth fiducial: new opt-in NTuple flag `truth_muon_fiducial_target`
  (|η_μ| < 2.5 only, no p_T cut), allowed ONLY with external input (refuses nominal/private/fullsim);
  output name carries `_muon_fid_target`. The muon cuts (p_T > 4.5, |η| < 2.4) move into the A_sig
  numerator. Single-muon cut constants move to ParamsSet (no retyping).

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

- **Step 9a — investigation + reference doc (2026-10-07).** Two read-only subagents (scratch docs merged
  here, then deleted): previous standalone private sample = 5.02 TeV, nCTEQ15 nPDF, 4 beams, Pythia 8.307,
  5 slices, muon filter p_T > 3.5 (made by the 2023 code in the OLD repo; `docs/pythia_truth.md` said
  5.36 TeV — fixed); barcode = 0-based Pythia index (entry 0 = system), Pythia status, `eventWeight` =
  σ/N_tried in µb per file, isospin only via file counts; central/Athena = 1-based HepMC barcodes,
  HepMC status (native |status| for intermediates), nb AMI weights. Full record:
  `Analysis/docs/private_pythia_sample_notes.md` (NEW reference doc).
- **Step 9b — filter-efficiency scan (2026-10-07, standalone Pythia 8.310 with the ATLAS A14 list +
  `NNPDF23_lo_as_0119_qed`, pp 5.36 TeV, no EvtGen; Condor clusters 12, 13;
  `~/usatlasdata/pythia_private_nomupt/effest_standalone/`).** Filter = OS pair, both |η_μ| < 2.5,
  pair p_T > 8 GeV:
  | slice | σ_gen | ε_filt |
  |---|---|---|
  | HardQCD 0–5 | 2422 mb (unphysical: only the 1 GeV `pTHatMinDiverge` cut) | 0 / 1e7 |
  | SoftQCD:nonDiffractive, hardest pTHat < 5 (physical complement) | 48.3 mb × 0.703 | 1 / 4.0e7 = 2.5e-8; 0 single-b target pairs |
  | 5–8 | 28.3 mb | 9.0e-7 |
  | 8–14 | 6.20 mb | 1.4e-5 |
  | 14–24 | 0.816 mb | 1.2e-4 |
  | 24–40 | 0.102 mb | 4.1e-4 |
  | 40–70 | 13.1 µb | 1.2e-3 |
  | 70–125 | 1.18 µb | 2.4e-3 |
  | 125–300 | 82.5 nb | 3.9e-3 |
  **Lowest slice = pTHat 5–8** by the user's rule (ε of the physical pTHat < 5 complement 2.5e-8 < 1e-7;
  its σ·ε = 8.5e-7 mb ≈ 0.3 % of the filtered σ, no single-b target pair seen). HardQCD below ~5 GeV is not
  a valid event cross section (σ ≫ σ_inel) — not generated.
- **Step 10 — code (2026-10-07).** `PythiaAlgCoreT.{h,c}`: external-mode options
  `truth_external_sigma_eff_nb` + `truth_external_weight_source` (weight = σ·ε_filt/N_processed nb; exactly
  one of raw-counts / weighted), `truth_external_kin_edges`, `truth_muon_fiducial_target` (PassCuts |η|≤2.5,
  no p_T; output tag `_muon_fid_target`; all refuse non-external runs); meta_tree_out gets
  `external_sigma_eff_nb` + UserInfo provenance. `ParamsSet.h`: `SINGLE_MU_PT_MIN/ETA_MAX`,
  `TARGET_MU_ETA_MAX` + expression helpers (PassCuts now reads them; nominal values identical).
  `SignalAcceptanceDefinition.h`: `TargetValue()` + `Target{Denom,NumMuPt,NumMuEta,NumMu,NumMass,Num}Expr`.
  RDF: `input/output_file_override` + `muon_fiducial_target` (name must carry `_muon_fid_target`),
  `RunRDFBasedHistFillingPythiaTruthTargetFiducial(in,out)`; target mode books
  `h2d_sig_accept_{denom,num,num_mupt,num_mueta,num_mu,num_mass}_pt_eta` and stamps `TargetValue()`.
  Local end-to-end test (pTH70_125, 20 events: Gen_tf → TRUTH0 (Athena 25.0.72) → TruthTrigRates skim
  (AthAnalysis 25.2.89) → NTuple target+weighted → RDF): 20/20 events have a filter-passing OS pair, muons
  down to 0.12 GeV kept; weight = 2.394 nb / 20; guards throw as designed.
- **Step 9c — target fraction + muon-cut factor per slice (standalone, 2026-10-07, cluster 13 final).**
  f_tgt = events with a single-b T_target pair / filter-passing events; A_μ = fraction of T_target pairs
  with both muons p_T > 4.5 & |η| < 2.4 (single-b = same eldest b-hadron via mother1):
  | slice | ε_filt | f_tgt | σ·ε·f_tgt [mb] | A_μ (N_tgt) |
  |---|---|---|---|---|
  | 5–8 | 8.3e-7 | 0.034 | 8.1e-7 (≈1 %) | — (1) |
  | 8–14 | 1.31e-5 | 0.267 | 2.17e-5 | 0.17 (35) |
  | 14–24 | 1.17e-4 | 0.337 | 3.22e-5 | 0.08 (195) |
  | 24–40 | 4.66e-4 | 0.378 | 1.80e-5 | 0.17 (335) |
  | 40–70 | 1.19e-3 | 0.319 | 5.0e-6 | 0.34 (667) |
  | 70–125 | 2.38e-3 | 0.291 | 8.2e-7 | 0.41 (1178) |
  | 125–300 | 3.94e-3 | 0.276 | 9.0e-8 | 0.54 (1524) |
  ★ σ-weighted A_μ ≈ 0.15: most single-b pairs with p_T^pair > 9 GeV have one muon below 4.5 GeV (the
  soft b→c→μ cascade muon). The PD4 extrapolation is a factor ≈ 2–10, NOT small — Athena/EvtGen check
  running (Condor cluster 14, slices 24–40 / 40–70) before raising it with the user.
- **Step 9d — Athena/EvtGen confirmation (2026-10-07).** Condor cluster 14: 20 Gen_tf jobs × 20 events,
  pTH24_40, the production JO → 400 accepted / 911 747 generated (ε = 4.39e-4; standalone 4.66e-4),
  σ = 101 805 nb; TRUTH0 + skim + NTuple (`truth_muon_fiducial_target`, weighted) →
  `PythiaCentralProdJO/evgen_nomupt_acceptance/condor_check/chain_24_40/`. With the nominal ancestry code
  (`from_same_b`): 120 T_target pairs, **A_μ = 0.19 ± 0.04** (standalone 0.17) — the small muon-cut factor
  is CONFIRMED with EvtGen. No upstream |η| cut: 13/460 OS pairs have a muon at 2.4 < |η| ≤ 2.5.
  pTH40_70 batch failed (Gen_tf requires nEventsPerJob ∈ {1,2,5,10,20,25,50,100,200,500,1000}); resubmitted
  (cluster 15). Grid CPU per trial from the 2026-09 request logs: 21–26 ms (41 ms at 125–300).
- **PROPOSED production plan (NOT submitted — awaiting the user on the A_μ finding):**
  | slice | acc. events | evts/job × jobs | ≈ T_target pairs | ≈ CPU-h |
  |---|---|---|---|---|
  | 5–8 | 200 | 1 × 200 | ~7 (≈1 % of σ_tgt) | 1 340 |
  | 8–14 | 10 000 | 20 × 500 | 2 700 | 4 240 |
  | 14–24 | 100 000 | 200 × 500 | 34 000 | 4 750 |
  | 24–40 | 100 000 | 500 × 200 | 38 000 | 1 200 |
  | 40–70 | 100 000 | 1000 × 100 | 32 000 | 470 |
  | 70–125 | 100 000 | 1000 × 100 | 29 000 | 230 |
  | 125–300 | 100 000 | 1000 × 100 | 28 000 | 290 |
  Total ≈ 12 500 CPU-h. Rationale: ≥ ~3e4 T_target pairs per slice where cheap (A_μ to ≈ 1–2 % per coarse
  cell); 8–14 (28 % of σ_tgt, 0.4 CPU-h per accepted event) to ±0.007 absolute; 5–8 (≈1 % of σ_tgt,
  6.7 CPU-h per event) only enough to bound its share.
- **Step 9e — grid production SUBMITTED (2026-10-07 17:2x, user approved the plan as proposed).**
  `PythiaCentralProdJO/evgen_nomupt_acceptance/grid/{slices.txt,submit_evgen.log}`; outDS
  `user.yuhang.nomuptAcc.Py8EG_A14_pp_hQCD_OSDiMu8_nomupt_pTH<lo>_<hi>.v1`; jediTaskIDs 5–8 52929567,
  8–14 52929573, 14–24 52929578, 24–40 52929583, 40–70 52929587, 70–125 52929591, 125–300 52929592.
  Chain driver `grid/chain.py` (nohup, PID in `grid/chain.pid`, log `grid/chain.log`) → TRUTH0 → skim →
  download to `~/usatlasdata/pythia_private_nomupt/grid_download/<cfg>/`. PD4 question still open with the user.

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
- R4 B → J/ψ question — DEFERRED by the user 2026-10-07: remind the user after the muon cuts are done.
- POWHEG: its truth sample also has a generator muon filter (filter_effcy 0.003 bb) → POWHEG cannot
  give the muon-cut factor either; POWHEG A_sig^mass stays a cross-check of the mass-window part only.
- Wire A_sig into the crossx after unfolding exists.

## Latest Stage
**2026-10-07 — STOPPED FOR USER (physics): the muon-cut factor is A_μ ≈ 0.15 (σ-weighted; 0.08–0.54 by
slice; EvtGen-confirmed 0.19 ± 0.04 at pTHat 24–40)** — PD4's "expected small" premise does not hold: A_sig
would extrapolate the measured yield by ×≈5–10 into soft-muon phase space (b→c→μ cascade muon < 4.5 GeV),
with the model dependence that implies. Asked the user whether PD4 stands; the user approved the production plan (submitted 2026-10-07,
Step 9e) and asked what PD4 is (explained). Code (NTuple/RDF/plotter target modes) done + tested locally; not yet reviewed/committed.
Lowest slice = 5–8 (user rule). **RESUME (new session):** production SUBMITTED + code committed (Step 9e,
commits 7113829/df171f3/b0065cf/7ed2769). `PythiaCentralProdJO/evgen_nomupt_acceptance/grid/chain.py` runs
under nohup (PID in `grid/chain.pid`; state `grid/chain_state.json`; log `grid/chain.log`). On resume: check
the log (`STOP:` = a task ended non-clean → inspect on BigPanDA, retry failed jobs or accept, then restart
`nohup python3 grid/chain.py`), check `voms-proxy-info -timeleft` (user renews), start a log waiter. Remaining
after download: NTuple per slice (weighted, σ·ε from that slice's log.generate: ε=Σacc/Σgen, σ gen-weighted
mean; `truth_external_kin_edges={5,8,14,24,40,70,125,300}`, `truth_muon_fiducial_target`) → hadd →
`RunRDFBasedHistFillingPythiaTruthTargetFiducial` → sanity plots + `plot_signal_acceptance_pythia_target.cxx`
→ /review-plot. PD4 question open with the user; R4 deferred. Superseded plan text: `grid/slices.txt` from the plan → `grid/submit_evgen.sh` →
`grid/chain.py` → NTuple (weighted, target) per slice → hadd → RDF target → sanity + A_sig plots → reviews.
