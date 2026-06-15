# Open Heavy Flavor & Quarkonium at RHIC (experimental review)

**Source:** Z. Tang, W. Zha, Y. Zhang (USTC), "An experimental review of open heavy flavor and quarkonium production at RHIC"
**arXiv / DOI:** arXiv:2105.11656 [nucl-ex] (v1, 25 May 2021)
**PDF:** *(not committed — SUPPORTIVE)* — `curl -sSL https://arxiv.org/pdf/2105.11656 -o /tmp/2105.11656.pdf` then read with `gs`
**URL:** https://arxiv.org/abs/2105.11656
**Classification:** SUPPORTIVE / complementary (RHIC-energy counterpart to the LHC-focused field reviews already in the KB)
**Added:** 2026-06-14

RHIC (STAR/PHENIX) experimental survey of open heavy flavor (charm + bottom)
and quarkonium in Au+Au and p+p at √s_NN = 200 GeV. Its value to us is twofold:
(1) it provides the **RHIC-energy** open-HF picture to contrast against the
LHC-energy reviews ([[open_hf_production]], [[hi_big_picture]]); (2) its
**references** point to concrete b/c-separation and bottom-RAA methodologies. We
keep the **open-bottom** slice thorough, open charm moderate, quarkonium brief.

---

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| **Open-bottom measured indirectly via non-prompt decay products** (B→J/ψ, B→D0, b→e) — RHIC has no direct open-b hadron reco (low σ_b, tiny hadronic BR) | Motivation/context for our **indirect single-b → dimuon** proxy; IntNote §1 Intro; thesis ch.1 | [background-for-writing] |
| **Data-driven b/c isolation from HF decay electrons** (template fit; ref [69]) and DCA/impact-parameter template fit (refs [65–67]) | Analogy for our signal/background separation & template-fit framework (analysis_overview §6; roadmap §16); methodology lead | [background-for-writing] / method-context |
| **b→e less suppressed than c→e at high p_T** (RHIC) → mass-ordered energy loss ΔE_b<ΔE_c | Comparison/contrast point for our b-quark R_AA(pair p_T); IntNote Results / Intro energy-loss framing | [background-for-writing] / interpretation |
| **Bottom v₂ small / beauty not thermalized at RHIC** (vb→e deviates from B-meson NCQ scaling, 99% CL) | Context for our *potential* v₂ observable (analysis_overview §3.3, roadmap §Q2.9); RHIC↔LHC contrast | [background-for-writing] |
| R_AA(p_T) = y_AA / (N_coll · y_pp) definition (Eq. 4) | Same nuclear-modification ratio as our primary observable (analysis_overview §3) | [background-for-writing] |
| D0 R_AA at 200 GeV **comparable to LHC 2.76 TeV despite large √s difference** | Concrete RHIC↔LHC energy-comparison data point for open HF | [background-for-writing] / interpretation |

---

## Scope & condition-difference warnings

- System/energy: **RHIC Au+Au √s_NN = 200 GeV** (p+p baseline at 200/500/510 GeV).
  Our analysis is **LHC Pb+Pb Run 3, √s_NN = 5.36 TeV**. ACKNOWLEDGE the large
  energy gap (200 GeV → 5.36 TeV); the review *notes* charm R_AA is similar at
  RHIC and LHC despite the gap, but **DO NOT assume the size/direction** of any
  difference for spectra, σ_b, T_AA, or suppression in our case.
- Channel: heavy flavor here is accessed via **D-mesons, Λ_c, and decay
  electrons** with silicon-vertex (STAR HFT) tagging — **not** our single-b →
  dimuon (two muons from one b chain). Treat all results as comparison/context,
  **not** a method we reproduce one-to-one.
- Beam-energy regime differs qualitatively: at RHIC thermal HQ production is
  "negligible," and bottom statistics are very limited (mostly STAR Preliminary
  / proceedings). Review date 2021 — pre-Run 3.

## Content summary   (relevance-first)

### Open bottom production (§II.B) — most relevant
- **Why indirect:** bottom mass ~3× charm → smaller energy loss expected, but the
  **low b production cross section and tiny hadronic-decay BR prevent direct
  open-b hadron reconstruction at RHIC.** STAR instead uses the **HFT** silicon
  vertex detector to separate the long-lived b-decay products from prompt charm
  via **impact parameter (DCA)** distributions — i.e. a **template fit** on the
  differing DCA shapes of signal vs background (refs [65–67]).
- **Measured channels:** non-prompt **B→J/ψ**, **B→D0**, and **b→e** at
  mid-rapidity in √s_NN = 200 GeV Au+Au (Fig. 7). Centralities differ per channel:
  non-prompt J/ψ 0–80%, non-prompt D0 0–10%, electrons 0–80%.
- **R_AA result (Fig. 7):** B→J/ψ suppressed across the whole 2–8 GeV/c range;
  B→D0 and b→e suppressed at high p_T → bottom-quark energy loss in the medium.
  Crucially, **b→e is systematically *less* suppressed than c→e**, and non-prompt
  D0 at ~4 GeV/c shows ~no suppression → **bottom loses less energy than charm**
  (mass-dependent parton energy loss). The **DUKE transport model** [68]
  reproduces the data within uncertainties.
- **Bottom v₂ (Fig. 8, ref [69]):** data-driven extraction of beauty v₂ by
  subtracting the (precisely measured) open-charm electron contribution from the
  inclusive HF electron spectrum. Electron v₂ from beauty hadrons at p_T > 3 GeV/c
  is ~**4σ** above zero (χ²/ndf = 29.7/6). At p_T < 4 GeV/c **vb→e < vc→e**
  (larger beauty mass), and vb→e **deviates from the B-meson NCQ-scaling
  hypothesis at 99% CL** (χ²/ndf = 14.3/4 over 2.5–4.5 GeV/c) → **beauty is
  unlikely thermalized and too heavy to follow the collective flow** of lighter
  partons at RHIC energy. (Contrast: charm v₂ *does* follow light-flavor NCQ
  scaling — see open charm.)

### Open charm production (§II.A) — moderate
- **D0 R_AA** in 0–10% central Au+Au at 200 GeV is **comparable to LHC Pb+Pb
  2.76 TeV despite the large √s difference**; significant suppression at
  p_T > 5 GeV/c (similar level to light-flavor mesons → strong charm–medium
  coupling), with a characteristic **low-p_T "bump"** at p_T < 5 GeV/c. Duke and
  LBT transport models (sizable charm collective motion) describe the STAR data;
  p+p reference uncertainty dominates the R_AA systematics.
- **D0 v₂** is large and follows **NCQ scaling** like light/strange mesons →
  charm flows and **may reach thermalization**; transport diffusion coefficient
  extracted but with large uncertainties.
- **Hadronization:** **Λ_c/D0 and D_s/D0 enhanced** in Au+Au vs p+p → charm
  **coalescence/recombination** with thermal partons plays an important role.
  Integrated D cross section per binary collision shown vs centrality (N_part).

### Quarkonium (§III–V) — brief
- p+p baseline: inclusive J/ψ production cross section, p_T dependence, and
  **polarization** (HX / Collins-Soper frames) measured; compared to CEM, NRQCD,
  CGC+NRQCD. **Feeddown** to inclusive J/ψ characterized (χ_c ~14% at 2 GeV/c;
  non-prompt J/ψ←B fraction rising with p_T/energy).
- In-medium: two competing hot effects — **color-screening dissociation**
  (suppression; T_d ordering, excited states melt first) vs **(re)combination**
  of uncorrelated heavy quarks (enhancement); plus CNM effects (nPDF, breakup).
  **Υ sequential melting** measured by STAR: Υ(1S) and Υ(2S+3S) R_AA vs N_part,
  compared to TAMU (Rapp) transport and Rothkopf calculations. J/ψ v₂ also shown.
- This is the bound-state (hidden-HF) program — distinct from our **open**-b
  signal; kept here only as Intro contrast.

---

## References worth future reading   (≤3)

1. **F. Si, X.-L. Chen, …, X. Dong, N. Xu, "Charm and beauty isolation from heavy
   flavor decay electrons in Au+Au collisions at √s_NN = 200 GeV," Phys. Lett. B
   805 (2020) 135465, arXiv:1906.08974** [ref 69] — *SUPPORTIVE*. New info: a
   concrete **data-driven b/c separation** of semileptonic HF decay electrons via
   template fit. Serves our **signal/background template separation**
   (analysis_overview §6, roadmap §16) and bottom-vs-charm interpretation — the
   closest published analog to isolating a bottom contribution from inclusive HF.
2. **S. Cao, G.-Y. Qin, S. A. Bass, "Energy loss, hadronization and hadronic
   interactions of heavy flavors in relativistic heavy-ion collisions," Phys. Rev.
   C92 (2015) 024907, arXiv:1505.01413** [ref 68] — *SUPPORTIVE*. New info: the
   **DUKE** heavy-flavor transport model that reproduces the RHIC open-bottom R_AA
   and v₂ here. Serves interpretation/systematics of our b-quark R_AA (a specific
   energy-loss model; complements the EMMI/Beraudo model comparison already queued
   under [[open_hf_production]]).
3. **STAR open-bottom R_AA via displaced J/ψ, D0, electrons** (X. Chen, PoS Hard
   Probes 2018 158; S. Zhang, Int. J. Mod. Phys. Conf. Ser. 46 (2018) 1860014)
   [refs 65–67] — *SUPPORTIVE (proceedings)*. New info: the **DCA / impact-
   parameter template-fit** method to tag non-prompt (bottom) products, and the
   measured RHIC bottom R_AA itself. Serves vertex-based b-tagging methodology
   context and the RHIC bottom-R_AA comparison point (look for the journal
   version when citing).

## Related KB docs   (knowledge graph)
- [[open_hf_production]] — Dong/Lee/Rapp open-HF review (RHIC **and** LHC). This
  doc is the RHIC-weighted experimental complement; that one is the broader
  theory+data review (diffusion, 2πT D_s, energy-loss mechanisms, dead cone).
- [[hi_big_picture]] — QGP/HIC big-picture motivation (R_AA, v_n, Glauber);
  supplies the field framing this RHIC survey sits inside.
- *(physics/ also holds the "Heavy-flavor production … hot QCD matter" review PDF,
  not yet summarized — another LHC-focused open-HF complement.)*
- [[analysis/overview]] — our single-b → dimuon goal, R_AA observable, and the
  background/template separation these RHIC b/c-isolation methods inform.
