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

## References worth future reading   (≤3)
1. **ATLAS, *Properties of g→bb̄ at small opening angles in pp at √s=13 TeV*,
   arXiv:1812.09283 (Phys. Rev. D 99, 052004).** — **PRIMARY candidate (strong).**
   New info: a direct measurement of g→bb̄ at **small opening angle** — exactly our
   nearby-pair regime — using R=0.2 track-jets / b-tagging; the region with large
   theory uncertainty our background template must get right. Serves: background
   template validation (roadmap §16), Intro background discussion. *Worth a full KB
   doc later.*
2. **ATLAS, *Measurement of b-hadron pair production in pp at √s=8 TeV*,
   arXiv:1705.03374.** — SUPPORTIVE. New info: bb̄ pair kinematics incl.
   flavour-creation vs flavour-excitation vs gluon-splitting decomposition;
   small-ΔR / collinear bb̄ region. Serves: background composition framing.

## Related KB docs
- `Analysis/docs/analysis_overview.md` §6 — our background & template program (ground truth; external to KB).
- [[run2_dimuon_backtoback_paper]] — quantifies the charm-vs-bottom dimuon handle (both muons pT>4 GeV).
- [[hf_hot_qcd_matter]] — the DD̄ IMR dimuon continuum (same two-heavy-quark source in HI).
- [[open_hf_production]] — HF energy-loss/production field hub.
- [[decisions]] — analysis decisions; `Analysis/docs/powheg.md` has the Powheg-truth NLO template role & the stitch with a Pythia g→QQ̄ template.
</content>
