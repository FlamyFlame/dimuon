# Physics-Results Review Criteria (MANDATORY for every physics result)

> **Why this exists.** A trigger-efficiency bug (compiled TF1s evaluating to 0
> above the fit range → ε floored to 0.01 → pp 2mu4 product weight blow-up)
> produced a ×13–88 discontinuity in the pp cross-section and R_AA above
> pair p_T ~58 GeV, and the plot review PASSED it. The reviewer checked
> readability/binning/conventions but never asked *"is this distribution
> physically sensible?"*. These criteria close that gap. They apply to **any
> review of a physics result** — a plotted distribution (`/review-plot`) or a
> quoted number (numerical verification inside `/review-analysis-code`,
> `/review-note`, `/review-paper`). The reviewer MUST apply every applicable
> item below and treat failures at the stated severity.

A "physics result" = any efficiency, cross-section, yield, R_AA, ratio,
correction factor, spectrum, or other measured/derived quantity (distribution or
single number). Apply the items flagged for your reviewer type (PLOT and/or
NUMERICAL). When an item applies and fails, it is **CRITICAL** unless stated
otherwise — a PASS verdict is forbidden with any CRITICAL open.

---

## C1 — Discontinuity / smoothness check  *(PLOT; also NUMERICAL when a number is one point of a distribution)*

Scan every distribution (each panel, each curve) for any **discontinuity, sudden
jump, kink, spike, or drop-to-zero** — up OR down — between adjacent bins or
between a corrected and uncorrected version of the same quantity.

For each one found, ask: **is it expected from the physics?** Almost always it is
NOT. Physical spectra, efficiencies, and ratios are smooth functions of kinematics
(a kinematic threshold, a known bin-edge effect, or an explicitly stated selection
boundary are the only legitimate sources — and those must be named).

- An **unexplained** discontinuity → **CRITICAL FAIL**. Do not rationalize it as
  "low statistics" unless the error bars genuinely span the jump AND the neighbors;
  a localized step far larger than its error bars is never just statistics.
- A discontinuity may signal a **physics-procedure / methodology problem**, not
  merely a cosmetic bug. It MUST trigger the investigation protocol (C4).

## C2 — Shape & magnitude vs physics expectation  *(PLOT and NUMERICAL)*

Do this **rubric-first**, before looking at the result:
1. From the task description alone, write down **what is plotted/quoted** (observable,
   axes, dataset, selection, units).
2. Form an **explicit list of physics expectations** for it — shape AND magnitude.
   Examples of expectation rubrics:
   - A cross-section / yield vs p_T is a **steeply, smoothly falling** spectrum
     (roughly power-law / exponential); it must not rise at high p_T.
   - A trigger or reconstruction **efficiency** is in [0,1], rises through a
     turn-on, then **plateaus near a sensible value** (single-muon HLT_mu4/2mu4-leg
     plateau is typically high, ~0.8–0.98 in the central acceptance; far lower is
     suspicious — see C3).
   - **R_AA** is O(0.1–1.5) for this measurement; a smooth, slowly varying function
     of p_T/centrality; central more suppressed than peripheral.
   - A **correction factor** (1/ε) is ≥ 1 and smoothly varying; values reaching
     tens–hundreds in isolated bins are a red flag.
   - A mass spectrum has the **expected resonance/continuum structure** in the
     stated window.
3. Check the actual result against **each** rubric item. Note every miss.

Any result that **violates its expected shape or sits orders of magnitude off the
expected scale** → **CRITICAL FAIL** (it is a physics-results check failure, even
if the code "runs" and the binning/labels are fine).

## C3 — Run 2 reference cross-check  *(PLOT and NUMERICAL; for final/near-final physics results)*

Sanity-check magnitude and shape against the **Run 2 reference analyses**, sourced
through the **KB index** (`.claude/kb/index.md` — the reference channel):
- **Run 2 HF-muon R_AA** — `.claude/kb/analysis/run2_hf_muon_raa.md`
  (arXiv:2109.00411 + internal note): R_AA scale/shape, T_AA normalization,
  single-muon efficiencies, b/c hierarchy.
- **Run 2 back-to-back dimuon** — `.claude/kb/analysis/run2_dimuon_note.md` +
  `run2_dimuon_backtoback_paper.md` (arXiv:2308.16652): dimuon yields, per-pair
  efficiency method, signal region.
- **Detector performance** — `.claude/kb/physics/detector/atlas_run2_muon_trigger.md`
  and `.../atlas_run3_muon_performance.md`: trigger/reco efficiency magnitudes
  (mu4/2mu4) for the trigger/magnitude comparison.

Procedure: identify whether a **corresponding Run 2 plot or number exists**; if so,
compare against it. **Account for the known differences** (do not expect equality —
expect *consistency after accounting for them*):
- **Analysis differences:** this is a **nearby muon-pair** measurement, not single
  HF-muon nor back-to-back dimuon → different signal region, different background,
  **non-negligible pair correlations** (vs. back-to-back), a muon **pair** rather
  than one HF-decay muon, and a different b-hadron kinematic region (vs. HF-muon).
- **Run 2 → Run 3 differences:** detector conditions, performance, and luminosity
  differ; detector performance should, in principle, **improve moderately** from
  Run 2. So for the **same trigger**, efficiencies should be comparable-to-somewhat-
  better than Run 2. If a trigger or reconstruction efficiency differs **greatly**
  from Run 2 — **especially much lower** — that clearly signals a problem.

A clear, unexplained magnitude/shape inconsistency with Run 2 (after accounting for
the above) → **CRITICAL FAIL**. If the Run 2 reference cannot be located/compared
in-session, or no Run 2 analog exists for this quantity, report it explicitly as
`RUN2-CROSSCHECK UNVERIFIED` (do not silently skip, do not invent a Run 2 value).
**`RUN2-CROSSCHECK UNVERIFIED` is a neutral status annotation, NOT a WARNING or
CRITICAL — it must not block a PASS.** It only blocks PASS if there genuinely IS a
comparable Run 2 result and it is inconsistent (then it is a C3 CRITICAL above).

## C4 — Auto-trigger investigation on a physics issue  *(control flow — the command, not the reviewer subagent, executes this)*

If the reviewer flags a **physics-results** issue under C1/C2/C3 (an unphysical
discontinuity, a shape/magnitude violation, or a Run 2 inconsistency it cannot
explain), the failure must NOT be closed with a cosmetic amend. Instead the
calling command MUST:
1. **Launch a physics investigation with tracking** — invoke `/review-investigation`
   on the flagged issue (it creates/uses a tracking doc per CLAUDE.md and runs its
   own executor↔reviewer loop). Pass the reviewer's evidence (which bins, which
   panels, expected vs observed, the C-item that failed).
2. **Fix per the investigation's root cause**, then re-run the originating review.
3. **Escalate if unresolved** — if the investigation cannot resolve it even after a
   thorough loop (escalates / hits max iterations), stop and **escalate to the user**
   with the full investigation report (root-cause status, what was tried, what
   remains), rather than forcing a PASS or an unexplained amend.

A purely cosmetic finding (legend overlap, axis label, directory) does NOT trigger
C4 — fix it with a normal amend. C4 is for *physics* failures only.
