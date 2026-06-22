# Gluon Splitting (g→QQ̄) and Flavour Excitation — correlated dimuon backgrounds

**Sources:**
- PI-provided talk: M. Mangano / "Gluon splitting to bottom quarks at the LHC",
  indico https://indico.cern.ch/event/557400/contributions/2372773/attachments/1375225/2087952/FCCeegbb.pdf
  — **content not auto-retrievable** (indico is behind an Anubis anti-bot
  proof-of-work; `curl` returns an HTML challenge, not the slides). URL recorded
  as a human pointer; **no slide content is summarized below** (GUIDE §5 — do not
  fabricate). The mechanism below is established QCD + our own
  `analysis_overview` §6.
- Citeable anchor (see Future-read): ATLAS, *Properties of g→bb̄ at small opening
  angles in pp at √s=13 TeV*, arXiv:1812.09283 — **directly our nearby-pair
  regime.**
**Classification:** SUPPORTIVE / background-concept (URL-only sources).
**Added:** 2026-06-15

> **Why this doc exists.** The two correlated backgrounds to our single-b → dimuon
> signal both produce **two muons from two *different* heavy quarks** with low pair
> mass / small opening angle, mimicking the single-b topology. This doc names and
> separates them so they are not conflated with the signal (one b chain → b→μ and
> b→c→μ). It is the background-concept hub; specific measurements are linked.

## Relevance to this analysis   (specific)

| What | Where it serves OUR work | Use type |
|---|---|---|
| **g→QQ̄ (final-state radiation, FSR):** a hard gluon splits to a QQ̄ pair; both Q decay semileptonically → two muons, **small opening angle / low mass** | The dominant *correlated* background to the single-b dimuon signal — `analysis_overview` §6; background/template discussion (IntNote §16, roadmap §16/task_07) | [background-for-writing] |
| **Flavour excitation (initial-state radiation, ISR):** a heavy quark from the proton sea (from a g→QQ̄ ISR splitting) is put on shell by the hard scatter; its partner Q is produced too → two heavy quarks, also low-mass/small-angle dimuons. Often grouped loosely (and imprecisely) with "gluon splitting" | Second correlated background; must be modeled/separated; Intro/background framing | [background-for-writing] |
| At **NLO** the classic flavour-creation / flavour-excitation / gluon-splitting "topologies" merge into the 2→3 matrix element — the low-mass/low-kT g→QQ̄ tail is the piece NLO (Powheg-HQ) under-populates | Justifies the **stitch** of Powheg-HQ NLO + a **Pythia g→QQ̄ template** for the background (and signal-mix) templates — `analysis_overview` §6 (`docs/powheg.md`) | [method-we-use] (template construction) |
| **Charm handle:** requiring **both muons pT > 4 GeV** kinematically suppresses the cc̄ (charm) dimuon contribution by ~an order of magnitude vs bottom (Powheg+Pythia8) | Why our signal region is bottom-dominated; quantified in [[run2_dimuon_backtoback_paper]] | [background-for-writing] |
| **OS/SS charge structure** of signal (OS-only) vs g→cc̄ (OS-only) vs g→bb̄ (both signs) → the **k = G_SS/G_OS** like-sign handle | **SS-anchored normalization** of the correlated background in the low-mass dimuon template fit — `docs/tracking/low_mass_dimuon_template_fit.md`, IntNote §16, roadmap step 16 | [method-we-use] |

## Scope & condition-difference warnings
- The mechanism is process-level QCD, energy-independent in concept; specific
  *rates/shapes* of g→bb̄ are measured at √s = 13 TeV pp (arXiv:1812.09283) and
  carry **large theoretical uncertainties** (the reason that paper exists). Do NOT
  assume the size/direction of any rate when transferring to our 5.36 TeV pp or
  Pb+Pb (GUIDE §5).
- "Flavour excitation" vs "gluon splitting" terminology is LO-language; at NLO the
  separation is not well defined (see above). Use carefully in writing.

## Content (established mechanism — attributed, not from the blocked talk)
- **Signal vs these backgrounds (the key distinction):** our **signal** = one
  b-hadron's decay chain gives both muons (b→μ and b→c→μ) → genuinely a *single*
  parent. These **backgrounds** = two muons from **two different** heavy quarks
  (a QQ̄ pair from one gluon, or two excited heavy quarks). Both can give low pair
  mass and small ΔR, so kinematics alone do not fully separate them — hence the
  template/purity program (`analysis_overview` §6).
- **Connection to the IMR DD̄ dimuon continuum:** the same two-heavy-quark dimuon
  source appears in the heavy-ion literature as the intermediate-mass-region (IMR)
  continuum from correlated DD̄ semimuonic decays — see [[hf_hot_qcd_matter]]
  (Averbeck §2.1/§5). In-medium effects modify both the energy loss and the
  opening angle of the QQ̄ pair.
- **Template strategy (our analysis):** Powheg-HQ supplies the NLO-hard spectrum
  but is strongly suppressed in the low-mass/low-kT g→QQ̄ tail (kinematic stacking:
  the splitting gluon needs pT ≳ 2m_Q, each Q hadronizes and decays to a muon
  carrying a fraction of the momentum, and each muon must pass pT > 4 GeV). The
  anticipated approach stitches Powheg-HQ with a Pythia g→QQ̄ template to supply
  that tail (`analysis_overview` §6; `docs/powheg.md`).

## OS/SS charge structure of signal vs g→QQ̄ background (the like-sign handle)  [method-we-use]

**Why this matters.** Our combinatorial subtraction and the *normalization* of the
correlated g→QQ̄ background both rest on how the muon-pair **charge** (opposite-sign
OS vs same-sign SS) distributes across sources. This is established decay physics
(semileptonic branchings + B⁰ mixing), and it is the basis for using the SS
spectrum to **data-anchor** the g→QQ̄ background in the low-mass dimuon template fit
(`docs/tracking/low_mass_dimuon_template_fit.md`; template-fit context in
[[muon_source_template_fits]]).

- **Signal (single-b → μμ): OS-only.** Both muons come from one b-hadron decay
  chain — direct b→μ⁻ and cascade b→c→μ⁺ — opposite-sign by construction. B⁰
  oscillation does not break this (both muons originate in the one b-hadron decay).
  SS signal ≈ 0 (confirmed in our MC — memory `project_overlay_pair_structure`).
- **g→cc̄: OS-only.** c→μ⁺ and c̄→μ⁻; charm has no sign-flipping cascade and D⁰–D̄⁰
  mixing is negligible (PDG x_D ≈ 0.4%), so the two muons are opposite-sign.
- **g→bb̄: BOTH OS and SS.** Two separate b-hadrons, each muon direct (b→μ) or
  cascade (b→c→μ): direct+direct and cascade+cascade → OS; direct+cascade → SS.
  Additionally **B⁰–B̄⁰ oscillation** flips a direct muon's sign → moves OS↔SS
  (time-integrated mixing probability χ_d ≈ 0.187, PDG).

**Consequence — the k ratio.** Writing the OS open-HF yield G_OS = (bb̄ OS)+(cc̄ OS)
and the SS yield G_SS = (bb̄ SS) only:
```
k ≡ G_SS / G_OS = (bb̄ SS) / [(bb̄ OS) + (cc̄ OS)] ,   0 < k < 1
```
The uncorrelated combinatorial is charge-symmetric (comb_OS ≈ comb_SS), so the
low-mass dimuon model is `OS = S + comb + G_OS`, `SS = comb + k·G_OS`.

**Why k is theory-robust (the payoff).** k is a *ratio* drawn from the **same**
g→QQ̄ population, so the absolute g→QQ̄ cross-section — which carries large theory
uncertainty (this doc's body; [[atlas_g_to_bb_small_angle]]) — largely cancels. k
is set mainly by *measured* inputs (semileptonic BRs for b→μ vs b→c→μ; χ_d). Its
one genuinely theory-sensitive ingredient is the **cc̄:bb̄ ratio in OS** (cc̄ enters
OS but not SS). This is why anchoring the correlated-background normalization via
**SS + k from data** beats taking the g→QQ̄ normalization from Pythia directly. The
cc̄:bb̄ sensitivity is mitigated: cc̄ and bb̄ have different m_inv shapes (the OS fit
partly constrains the mix), and the residual is taken as a systematic.

**Negative constraints (for the code).** k is bb̄-only — do NOT model g→cc̄ as
contributing to SS beyond the lumped approximation. The **leading** bb̄ SS
mechanisms are **direct+cascade pairing** (one b→μ direct, the partner b̄→c̄→μ →
same sign) and **single B⁰ mixing**; these are *captured* by a measured k, not
neglected. (Note: cascade+cascade is OS, not SS.) k must be **measured in MC**
(Pythia g→bb̄/g→cc̄ + flavour excitation), including its (pair pT, η, m_inv)
dependence. **k is EMPIRICAL — it absorbs everything with the QQ̄ charge
structure:** flavour excitation (ISR g→QQ̄, identical charge structure to FSR
gluon-splitting) AND all bb̄ charge-flip mechanisms (direct+cascade, single AND
double B⁰ mixing, b-baryons) are *in* the measured k, not separately modeled or
"neglected." `SS = C + k·G_OS` can fail only from: (a) **SS sources outside
{correlated-HF, combinatoric}** — residual jet/fake muons and **>2 HF μ in a
jet** (the latter expected rare → if the data closure shows the OS/SS factor
roughly holds, carry it as a systematic rather than modeling it); (b) **MC
mismodeling of the ratio** (cc̄:bb̄ and
FE:GS fractions, mixing) → k_MC ≠ k_data, the systematic, checked by a data
closure test (predict SS = C_mixed-event + k·G_OS vs actual SS); (c)
**combinatoric charge asymmetry** C_OS ≠ C_SS; (d) uncaptured k(pair pT, η). Do
NOT leave `OS − SS` as the final background subtraction — it cancels
combinatorial but leaves `(1−k)·G_OS`.

## References worth future reading   (≤3)
1. **ATLAS, *Properties of g→bb̄ at small opening angles in pp at √s=13 TeV*,
   arXiv:1812.09283.** — ✅ **Now summarized in the KB: [[atlas_g_to_bb_small_angle]]**
   (the citeable anchor for this background, in our nearby-pair regime).
2. **ATLAS, *Measurement of b-hadron pair production in pp at √s=8 TeV*,
   arXiv:1705.03374.** — SUPPORTIVE. New info: bb̄ pair kinematics incl.
   flavour-creation vs flavour-excitation vs gluon-splitting decomposition;
   small-ΔR / collinear bb̄ region. Serves: background composition framing.

## Related KB docs
- [[muon_source_template_fits]] — template-fit concept hub; the OS/SS charge
  structure here underpins the **SS-anchored normalization** of the correlated
  background, used alongside the Δp/p fake-muon purity fit in the low-mass dimuon
  template fit.
- [[atlas_g_to_bb_small_angle]] — the ATLAS g→bb̄ small-opening-angle measurement (citeable anchor / our regime).
- `Analysis/docs/analysis_overview.md` §6 — our background & template program (ground truth; external to KB).
- [[run2_dimuon_backtoback_paper]] — quantifies the charm-vs-bottom dimuon handle (both muons pT>4 GeV).
- [[hf_hot_qcd_matter]] — the DD̄ IMR dimuon continuum (same two-heavy-quark source in HI).
- [[open_hf_production]] — HF energy-loss/production field hub.
- [[decisions]] — analysis decisions; `Analysis/docs/powheg.md` has the Powheg-truth NLO template role & the stitch with a Pythia g→QQ̄ template.
</content>
