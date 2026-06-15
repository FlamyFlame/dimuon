# Knowledge Base Index

A working reference for the ATLAS Run-3 single-b → dimuon analysis. Built/reviewed
via `/kb-build` + `/kb-review` per `KB_BUILDING_GUIDE.md`. **PRIMARY** = high
relevance, PDF committed; **SUPPORTIVE** = URL-only, relevant-slice summary.
Physics ground truth lives in `Analysis/docs/analysis_overview.md`.

## Meta (how this KB is built)
- [KB_BUILDING_GUIDE.md](KB_BUILDING_GUIDE.md) — authoritative criteria (PRIMARY/SUPPORTIVE, relevance-first, no-hallucination, ≤3 future-reads, PDF rule, dedup + knowledge graph, self-sync)
- [FUTURE_READ.md](FUTURE_READ.md) — queued references worth summarizing later (≤3/source)

## Analysis (papers ours derives from + internal notes)
- [overview.md](analysis/overview.md) — Physics goal, signal, backgrounds, datasets, analysis strategy (older brief; see `Analysis/docs/analysis_overview.md` for the maintained version)
- [decisions.md](analysis/decisions.md) — Key analysis decisions and rationale
- [run2_dimuon_note.md](analysis/run2_dimuon_note.md) — **PRIMARY.** Run 2 dimuon internal note (ATL-COM-PHYS-2021-1094): datasets, selections, per-pair corrections, Δp/p purity, Δφ fits, systematics; **the template for the Run 3 note** + publication-level framing
- [run2_dimuon_backtoback_paper.md](analysis/run2_dimuon_backtoback_paper.md) — **PRIMARY.** Published Letter of that analysis (PRL 132 (2024) 202301, arXiv:2308.16652): publication framing, away-side Δφ widths, per-pair eff method we adapt (+ΔR), σ_int≲0.2 limit (PDF in dir)
- [run2_hf_muon_raa.md](analysis/run2_hf_muon_raa.md) — **PRIMARY.** ATLAS Run 2 HF-muon charm/bottom **R_AA** (PLB 829 (2022) 137077, arXiv:2109.00411) — the single-HF-muon analysis ours derives from: R_AA/T_AA norm, FCal-Glauber centrality, Δp/p + d0 template fits, efficiencies, systematics (PDF + internal note in dir)
- [run2_hf_muon_vn.md](analysis/run2_hf_muon_vn.md) — **PRIMARY.** ATLAS Run 2 HF-muon azimuthal anisotropy **v_n** (PLB 807 (2020) 135595, arXiv:2003.03565): event-plane flow method, ρ/d0 template fits; charm v_n > bottom v_n (PDF + internal note in dir)
- [atlas_gammagamma_mumu_pbpb.md](analysis/atlas_gammagamma_mumu_pbpb.md) — **PRIMARY.** ATLAS γγ→μμ in non-UPC Pb+Pb 5.02 TeV (PRC 107 (2023) 054907, arXiv:2206.12594) — Ref [60] of the back-to-back Letter: the **actual per-pair trigger+reco efficiency method** (single-μ ε(pT,qη) product, valid for back-to-back; CONTRAST our pair ε_reco(pair pT,η,ΔR)) + γγ α/A/k⊥ removal (also our OS background). **+ companion internal note ANA-HION-2020-10** (full per-pair eff implementation) (PDFs in dir)
- [yun_tian_thesis.md](analysis/yun_tian_thesis.md) — **PRIMARY.** Yun Tian Columbia PhD thesis (CERN-THESIS-2018-078, 2018): in-depth companion to the **2.76 TeV HF-muon R_AA+v_n** analysis (arXiv:1805.05220 — Δp/p "Eloss" template fit Run-1 origin, MC-truth efficiencies, EP+scalar-product v_n); also **Ch7 HF dimuon correlations 5.02 TeV (Run-1 predecessor of OUR analysis)** and Ch8 γγ→μμ (arXiv:1806.08708, first observation, superseded by 2206.12594); **a writing template for the PI's own thesis** (PDF in dir)

## Concepts (cross-paper method hubs)
- [concepts/muon_source_template_fits.md](concepts/muon_source_template_fits.md) — Δp/p momentum-imbalance (real vs π/K fake) and d0 (charm vs bottom) template fits; how they differ across the 3 analysis sources; **our Δp/p purity framework**

## Physics — Heavy Ion (field background)
> Canonical hub for HF physics (energy loss, R_AA, dead-cone, transport): **open_hf_production**. Others defer shared concepts to it.
- [physics/heavy_ion/open_hf_production.md](physics/heavy_ion/open_hf_production.md) — **PRIMARY (hub).** Open HF in HIC review (Dong/Lee/Rapp, arXiv:1903.07709): R_AA & v2, collisional/radiative energy loss, dead-cone mass hierarchy (b<c<q<g), bottom-charm R_AA hierarchy, HQ diffusion 2πT D_s≈2–4, hadronization, CNM (PDF in dir)
- [physics/heavy_ion/hi_big_picture.md](physics/heavy_ion/hi_big_picture.md) — **PRIMARY.** QGP/HIC big-picture review (Busza/Rajagopal/van der Schee, arXiv:1802.04801): QGP-as-liquid (η/s≈1/4π), jet quenching, R_AA def, Glauber centrality, flow; Intro/thesis motivation (PDF in dir)
- [physics/heavy_ion/hf_hot_qcd_matter.md](physics/heavy_ion/hf_hot_qcd_matter.md) — **PRIMARY.** Earlier open-HF review (Averbeck, PPNP 70 (2013)): RHIC-era discovery, decay-lepton (e/μ) channel, DD̄ IMR dimuon continuum (PDF at physics/ root)
- [physics/heavy_ion/rhic_open_hf_review.md](physics/heavy_ion/rhic_open_hf_review.md) — SUPPORTIVE. RHIC (STAR/PHENIX) open HF + quarkonium (arXiv:2105.11656): open-b via DCA template fit, b→e less suppressed than c→e, bottom v2 small; RHIC↔LHC contrast
- [physics/heavy_ion/hf_theory_overview.md](physics/heavy_ion/hf_theory_overview.md) — SUPPORTIVE. HF-transport theory overview (Gossiaux, QM2018, arXiv:1901.01606): Langevin vs Boltzmann models, transport coefficients, R_AA–v2 puzzle, hadronization

## Physics — Detector (ATLAS muon & tracking performance)
> Canonical hub for muon types / WPs / tag-and-probe: **ATLAS_Run2_muon_reconstruction**. Trigger & Run-3 docs defer to it.
- [physics/detector/ATLAS_Run2_muon_reconstruction.md](physics/detector/ATLAS_Run2_muon_reconstruction.md) — **PRIMARY (hub).** Run 2 pp muon reco/ID/isolation/vertex efficiency (EPJC 81 (2021) 578, arXiv:2012.00578): 5 muon types, WPs, tag-and-probe & scale factors, results (PDF in dir)
- [physics/detector/atlas_run2_muon_trigger.md](physics/detector/atlas_run2_muon_trigger.md) — **PRIMARY.** Run 2 muon **trigger** performance (JINST 15 (2020) P09015, arXiv:2004.13447): L1+HLT chains (mu4/mu8/2mu4), tag-and-probe trig-eff, dimuon factorization ε_leg1·ε_leg2·ρ_μμ + ρ_ΔR close-by correction (PDF in dir)
- [physics/detector/atlas_run3_muon_performance.md](physics/detector/atlas_run3_muon_performance.md) — **PRIMARY.** Run 3 muon reco + trigger (trigger paper JINST 19 P06029/arXiv:2401.06630; reco = ATLAS-Preliminary): NSW added to WP, Run-3 thresholds, Pb+Pb trigger menu; **no full Run-3 reco paper or HI muon perf yet** (PDFs in dir)
- [physics/detector/atlas_inner_detector_tracking.md](physics/detector/atlas_inner_detector_tracking.md) — SUPPORTIVE. ID dense-environment tracking (EPJC 77 (2017) 673): merged/shared-cluster + ambiguity solving — the mechanism behind close-track correlated reco (**justifies our pair ε_reco ΔR axis**)

## Physics — Background (signal-mimicking processes)
- [physics/background/gluon_splitting_flavour_excitation.md](physics/background/gluon_splitting_flavour_excitation.md) — **concept hub.** g→QQ̄ (FSR) and flavour excitation (ISR): the two correlated low-mass/small-angle dimuon backgrounds (two muons from two heavy quarks vs our single-b); template-stitch context
- [physics/background/atlas_g_to_bb_small_angle.md](physics/background/atlas_g_to_bb_small_angle.md) — **PRIMARY.** ATLAS g→bb̄ at small opening angle, pp 13 TeV (PRD 99 (2019) 052004, arXiv:1812.09283): differential ΔR/z/mass via double-b-tagged R=0.2 track-jets; Pythia mismodels collinear g→bb̄ (large theory uncert.) → validates/cautions our g→QQ̄ background template (PDF in dir)

## Physics — Centrality (Pb+Pb geometry / Glauber)
> Canonical hub for Glauber methodology (N_part/N_coll/T_AA, percentile centrality, systematics): **glauber_modeling**. ATLAS centrality papers defer the mechanism to it.
- [physics/centrality/glauber_modeling.md](physics/centrality/glauber_modeling.md) — **PRIMARY (hub).** Glauber methodology review (Miller/Reygers/Sanders/Steinberg, ARNPS 57 (2007) 205, arXiv:nucl-ex/0701025): N_part/N_coll/T_AB definitions, Woods-Saxon + σ_NN inputs, optical vs MC Glauber, centrality from percentiles, ⟨T_AB⟩=⟨N_coll⟩/σ_NN, systematics (PDF in dir)
- [physics/centrality/atlas_centrality_2023.md](physics/centrality/atlas_centrality_2023.md) — **PRIMARY.** ATLAS centrality for **5.36 TeV 2023** Pb+Pb (internal note GROUP-2024-033): FCal-ΣE_T + PHOBOS Glauber; **Table 7 ⟨T_AA⟩ = the source of our hardcoded 6-bin T_AA (confirmed)**; centile cuts, N_part/N_coll, systematics; roadmap Q2.3 (σ_PbPb=7.8b is NOT here → Q2.2 open) (PDF in dir)
- [physics/centrality/atlas_centrality_2015.md](physics/centrality/atlas_centrality_2015.md) — **PRIMARY.** Earliest ATLAS FCal-ΣE_T centrality methodology note (5.02 TeV / 2015, ATL-COM-PHYS-2016-XXX): Glauber two-component fit to ΣE_T (x=0.09), 40 GeV gap WP, centile cuts, ⟨N_part/N_coll/T_AA⟩ + systematics; **method-depth** companion (numbers are 5.02 TeV, don't transfer) (PDF in dir)

## Data
- [samples.md](data/samples.md) — Dataset container names, runs, luminosity, data paths
- [variables.md](data/variables.md) — Branch names, sign conventions, ZDC/FCal indexing, units

## Procedures
- [procedures/setup_new_machine.md](procedures/setup_new_machine.md) — environment setup

## Gotchas
- [gotchas/reading_pdfs.md](gotchas/reading_pdfs.md) — read PDFs with `gs -sDEVICE=txtwrite`; WebFetch/Read can't
</content>
