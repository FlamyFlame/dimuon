# Future-Read List

References discovered while summarizing KB sources that offer **new, crucial,
analysis-relevant** info not covered by the source they were found in. Governed by
`KB_BUILDING_GUIDE.md` §6 (≤3 per source; keep lean). Deduped across the 2026-06
bulk build (≥1 paper was flagged by several builders).

When an item is summarized into the KB, move it to "Done" with a link.

## To read — HIGH priority (method/comparison-critical, mostly PRIMARY)

| Citation | arXiv/DOI | Class | New info it adds | Serves | Found in |
|---|---|---|---|---|---|
| Iterative Bayesian (D'Agostini) unfolding + RooUnfold | NIM A 362 (1995) 487; Adye arXiv:1105.1160 | SUPPORTIVE | Standard ATLAS-HI **unfolding** method — **no KB doc covers unfolding yet** | Detector response / unfolding (roadmap step 13, overview §4c) | atlas_gammagamma_mumu_pbpb |
| ATLAS c/b muon v2 in **pp** 13 TeV | arXiv:1909.01650 (PRL 124 (2020) 082301) | PRIMARY | **Origin of the d0 charm/bottom template-separation** method; pp HF-muon flow baseline | c/b separation concept; pp baseline | run2_hf_muon_raa, run2_hf_muon_vn |
| CMS bottom R_AA (B±, non-prompt J/ψ, non-prompt D0) | Sirunyan et al. (refs behind 1903.07709 Fig.5) | PRIMARY (comparison) | Measured **bottom** R_AA(pT) at LHC — closest existing comparison to our single-b R_AA | R_AA results comparison (step 15); Intro | open_hf_production |
| ALICE HF decay-**muon** R_AA, fwd-y, Pb+Pb 2.76 TeV | Averbeck ref [222] (©2012 APS) | PRIMARY (comparison) | The measured **muon-channel** HF R_AA (same decay channel as ours), explicit high-pT bottom fraction | R_AA results comparison; thesis prior measurements | hf_hot_qcd_matter |

## To read — MEDIUM priority (SUPPORTIVE background/method context)

| Citation | arXiv/DOI | Class | New info | Serves | Found in |
|---|---|---|---|---|---|
| ATLAS HF-muon suppression+v2, Pb+Pb 2.76 TeV | arXiv:1805.05220 (PRC 98 (2018) 044905) | SUPPORTIVE | Predecessor inclusive HF-muon R_AA (no c/b split) + original event-plane method + the π/K **Δp/p purity** procedure | History/Intro; Δp/p purity (step 16) | backtoback, raa, vn |
| ATLAS b-hadron pair production, pp 8 TeV | arXiv:1705.03374 | SUPPORTIVE | bb̄-pair kinematics with flavour-creation/excitation/gluon-splitting decomposition; small-ΔR region | Background composition framing | gluon_splitting_flavour_excitation |
| Andronic et al., HF + quarkonium in the LHC era | arXiv:1506.03981 (EPJC 76:107 (2016)) | SUPPORTIVE | More detailed dedicated HF+quarkonium pp→AA review (baseline, mechanisms, models) | IntNote §1 Intro; thesis ch.1 | open_hf_production, hf_theory_overview |
| Beraudo et al. (EMMI task force), HF transport coefficients | arXiv:1803.03824 / Nucl.Phys.A979:21 (2018) | SUPPORTIVE | Quantitative model-to-model HF energy-loss/transport comparison incl. bottom; H_AA uncertainty | b R_AA interpretation & systematics framing | open_hf_production, hf_theory_overview |
| Connors/Nattrass/Reed/Salur, jet measurements in HI | arXiv:1705.01974 (RMP 90 (2018) 025005) | SUPPORTIVE | Comprehensive **experimental** hard-probe/HF R_AA review | Intro; R_AA discussion | hi_big_picture |
| Si et al., charm/beauty isolation from HF decay electrons | arXiv:1906.08974 (PLB 805 (2020) 135465) | SUPPORTIVE | Concrete **data-driven b/c separation** of semileptonic HF leptons via template fit | Signal/bkg template separation (§6/§16) | rhic_open_hf_review |
| Loizides, Kamin, d'Enterria, Glauber improved params | arXiv:1710.07098 (PRC 97 (2018) 054910) | SUPPORTIVE | Origin of the Woods-Saxon / σ_NN parameter set + Glauber systematics used by the ATLAS centrality notes | Centrality systematics (overview §4e, roadmap Q2.3) | atlas_centrality_2023 |
| Ilten, Rodd, Thaler, Williams, disentangling heavy flavor | arXiv:1702.02947 (PRD 96 (2017) 054019) | SUPPORTIVE | Phenomenological statistical separation of g→QQ̄ from other HF production | Background composition / template separation (§6) | atlas_g_to_bb_small_angle |

> Deferred/merged to keep the list lean: DUKE (Cao 1505.01413) and Prino&Rapp
> (1603.00529) HF-transport models — folded into the Beraudo/EMMI comparison;
> Majumder&van Leeuwen (jet-quenching theory) — covered by 1705.01974; STAR
> open-b proceedings — await journal versions.

## Done
| Citation | KB doc | Date |
|---|---|---|
| ATLAS γγ→μμ non-UPC Pb+Pb (arXiv:2206.12594) | [[atlas_gammagamma_mumu_pbpb]] | 2026-06-15 |
| ATLAS g→bb̄ small opening angle (arXiv:1812.09283) | [[atlas_g_to_bb_small_angle]] | 2026-06-15 |
</content>
