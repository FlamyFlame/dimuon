# NTuple-Processing Provenance (MANDATORY)

**Why this exists.** A recurring, high-impact failure: an agent writes standalone code that
reads the **raw NTUPs** directly and silently diverges from the ntuple-processing procedure —
a missing selection cut, a skipped truth-navigation step, a wrong/absent weight, or an
unreplicated mode flag — producing wrong physics that a readability/label review passes. Two
real instances: (a) a standalone |d0|/Δp/p provenance macro that traced truth parents WITHOUT
the generator-block barcode cutoff → in the HIJING overlay the duplicated barcodes
(Geant4/HIJING at `truth_barcode>200000`) mis-mapped signal-HF muons' parents → real-HF
collapsed from ~92% to ~5%; (b) the AMI `crossSection` nb-vs-pb unit slip in a standalone
weighting macro. The processing code already handles these; re-deriving from raw repeats
every mistake it already solved.

## Rule (the ladder — try in order)

1. **Use ntuple-processing OUTPUT.** The processing (`Analysis/NTupleProcessingCode/*`) emits
   trees with all cuts/weights/truth-tracing already applied: muon_pairs (data + MC),
   single-muon trees, pythia-truth trees, and derived branches (flavor/origin categories,
   `from_same_b`, `_last_b_hadron`/HF ancestry, per-pair efficiency/weight columns). If the
   output can answer the question, USE IT.
2. **Extend the processing.** If the needed quantity is not produced, ADD a configuration
   mode/flag to the ntuple processing and rerun it (distinct output suffix; NEVER clobber
   nominal — see the clobber-incident history). Preferred over standalone-from-raw, especially
   when only one setting changes (e.g. a loosened-Δp/p single-muon output for a Δp/p plot).
3. **Standalone-from-raw (last resort).** Only if (1) and (2) are genuinely impossible. Then
   you MUST first READ the relevant processing code and reproduce its procedure EXACTLY:
   - every single-muon / pair selection cut (quality WP bits, pt, |η|, Δp/p, |d0|, |z0 sinθ|, …),
   - every truth-navigation step (esp. the **generator-block barcode cutoff** for HIJING
     overlay: `PythiaTruthExtras.c` `pythia_only_barcode_cache`, map+trace only over `[0,lim)`
     with `lim` = first `truth_barcode>200000`),
   - every weight (AMI nb→pb, luminosity, efficiency, T_AA), binning, and sign/convention,
   identical to nominal. The ONLY settings that may differ are those the request EXPLICITLY
   changes for a stated physics purpose.
4. **Unclear → STOP AND ASK (hard blocking point).** If the request does not say whether a
   given cut/procedure should be mirrored, and you think a standalone deviation is justified
   for a physics reason, do NOT decide unilaterally — ask the user before writing the code.

## Reviewer checklist (apply in `/review-analysis-code`, `/review-plot`, `/review-investigation`)

When the work under review is standalone code reading raw NTUPs (tree `HeavyIonD3PD` or raw
`*.root` skim files, rather than processed muon_pairs/single-muon/pythia-truth outputs):

- **P1 — Output-vs-raw justification.** Is there a reason the ntuple-processing output (or a new
  processing mode) could not be used? If processed output would have sufficed → **FAIL** (should
  not have gone standalone).
- **P2 — Procedure parity.** Read the corresponding processing code and diff the procedure. Every
  cut, truth step, weight, binning, and convention must match nominal EXCEPT settings the request
  explicitly changed. Any unstated divergence (missing cut, skipped truth step like the barcode
  cutoff, wrong/absent weight) → **CRITICAL FAIL** (this is a physics-results failure — route to
  investigation like C1–C4).
- **P3 — Sanity vs the processed sample.** Where feasible, cross-check a yield/fraction/shape
  against the processed output or the other sample (e.g. signal-vs-overlay HF fraction should be
  comparable; a class fraction that swings ×10+ between samples is a red flag for a provenance bug).
- **P4 — Deviation authorization.** Any intentional deviation must be explicitly requested by the
  user (or approved) AND documented in the tracking doc with its physics motivation. An
  agent-chosen deviation without user sign-off → **FAIL** (the rule requires STOP-AND-ASK).

A P2 CRITICAL is a *physics* failure — do not close it with a cosmetic amend; fix the code to
mirror the procedure (or add the processing mode) and re-review.
