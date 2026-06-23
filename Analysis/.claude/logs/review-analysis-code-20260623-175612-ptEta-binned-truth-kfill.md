# Analysis Code Review Log
**Task**: Add (pair pT, pair eta)-binned 2D truth-template fills to FillHistogramsTemplateMinvSignalRegion (RDFBasedHistFillingPythiaTruth.cxx) for k(m,pT,eta). 2D: {pair_pt_log_150,minv_zoomin}+{pair_eta,minv_zoomin}, _sigsel, signal selection unchanged.
**Log file**: review-analysis-code-20260623-175612-ptEta-binned-truth-kfill.md
**Started**: 2026-06-23T21:56:12Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Executor work (pre-review)
File modified: RDFBasedHistFilling/RDFBasedHistFillingPythiaTruth.cxx, FillHistogramsTemplateMinvSignalRegion.
Added local `minv_template_2D = {{"pair_pt_log_150","minv_zoomin"},{"pair_eta","minv_zoomin"}}` and passed it
as the 5th (vars2D) arg of FillHistogramsSingleDataFrame in BOTH the flavor and origin loops. kin_cuts,
sign->df map, 1D fill, _sigsel suffix all UNCHANGED. Axes resolve via Var1DSearch from var1D_pythia_truth.json
(pair_pt_log_150->pT_bins_150 [base-class hist_binning_map:139]; pair_eta->24 bins; minv_zoomin->50 bins 0-4).
Compile: ACLiC clean, fresh .so 17:57, TRUTH_OK ctor.
Rerun: run_rdf_pythia_truth.sh nonprivate 5.36 nocuts, exit 0 (63s); 1D=562, 2D=475, 3D=475 histos.
Verification:
- 1D _sigsel UNCHANGED: G_OS=2.2138, G_SS=0.6814 (matches pre-change 2.214/0.681).
- New 2D present + non-empty: h_minv_zoomin_vs_pair_pt_log_150_sign{1,2}_<cat>_sigsel (15x50) and
  h_minv_zoomin_vs_pair_eta_sign{1,2}_<cat>_sigsel (24x50), cats {bb,cc,one_b_one_c,single_b}.
- 2D integrals match 1D per category (bb OS 1.3159=1.3159; G as a whole consistent) -> no under/overflow loss.
- Sign mapping correct: k_bb=G_SS_bb/G_OS_bb=0.6734/1.3159=0.512; k_cc=0.0081/0.8979=0.009; one_b_one_c empty.
Backup: histograms_pythia_5p36TeV_no_data_resonance_cuts_pre_ptEta_kfill_backup.root.
Numbers reported: G_OS=2.2138, G_SS=0.6814 (1D, unchanged); 2D per-cat integrals as above.

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 0 (1 INFO non-blocking: executor's hist-count print 562/475/475 = RResultPtr-map sizes incl. not-written; on-disk TH keys 498/423/0 — no TH3 here. Does not affect correctness.)
**Details**: None blocking. All numbers independently re-extracted MATCH: G_OS=2.2138, G_SS=0.6814 (1D unchanged); 2D present 15x50(pT)/24x50(eta); 2D in-range integrals == 1D per cat; k_bb=0.5117, k_cc=0.0090, one_b_one_c empty. Sign mapping confirmed NOT swapped (k_bb~0.51 not ~1.95). Binnings registered. kin_cuts unchanged. Per-bin k spot-checks sane (bb ~0.48-0.50 well-populated, mild high-pT/edge drift = low-stat, not a violation).
**Numerical verification**: all MATCH (hist-count INFO artifact only).

**Status**: APPROVED at iteration 1
**Summary**: (pT,eta)-binned 2D truth template fills added (h_minv_zoomin_vs_{pair_pt_log_150,pair_eta}_sign{1,2}_<cat>_sigsel); 1D templates unchanged; sign mapping correct; ACLiC-clean; truth file regenerated. Raw material for k(m,pT,eta) ready.
