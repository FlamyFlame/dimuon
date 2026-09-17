#ifndef ParamsSet_h
#define ParamsSet_h

#include <cstdio>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <utility> 
#include <functional>
#include <cmath>

struct AxisInfo{ // a simplified version of axis information - for constructing TH1DModel
    int        	nbins		= 1;         ///<  Number of bins
    double     	min			= 0.;          ///<  Low edge of first bin
    double     	max			= 1.;          ///<  Upper edge of last bin
    const double*     bin_edges		= nullptr;         ///<  Bin edges array in X
};

class ParamsSet{
public:
	static const unsigned int nCtrBins=6; // number of coarse bins; each can be studied by itself
   	static const unsigned int nPtBins=10;
   	static const unsigned int nPythiaKinRanges=5;
   	static const unsigned int nSigns=2;

   	static const unsigned int ndRselcs=3;
   	static const unsigned int ndphiselcs=3;
   	static const unsigned int ndphiRegions=3;
   	static const unsigned int ndetaRegions=2;

   	static const unsigned int nEffCorr=2; // 0: no efficiency correction; 1: having effiency correction
   	static const unsigned int nPhotoProdCuts=2; //0: no gap cut; 1: having gap cut
   	static const unsigned int nGapCuts=2; //0: no gap cut; 1: having gap cut
  	static const unsigned int CtrStep = 5;
  	static const unsigned int nCtrIntvls = 20; // number of small intervals; the intervals not to be studied by themselves, but to allow easy combinations

	static std::vector<std::function<bool(float)>> dphi_cut_funcs;
	static std::vector<std::function<bool(float)>> deta_cut_funcs;

  	// float pTbins[nPtBins] = {4.,5.,6.,7.,8.,9.,10.,12.,15.,20.};
	// int scaleFactorCtrs[nCtrBins] = {1,1,2,3,3};
	// std::vector<float> ctrbins = {0, 5, 10, 20, 30, 50, 80};
    std::vector<double> pT_bins_40;
    std::vector<double> pT_bins_8;
    std::vector<double> pT_bins_60;
    std::vector<double> pT_bins_80;
    std::vector<double> pT_bins_120;  // OPT-IN `_pt_120` alternative view (16 bins, 9-120 GeV)
    std::vector<double> pT_bins_150;  // DEFAULT FINE pair-pT crossx axis (16 bins, 9-150 GeV);
                                      // nests 2:1 inside pair_pt_coarse_bins -- see constructor

    // ---- FINE pair-eta axis of the cross-section 2D/3D views (SINGLE SOURCE OF TRUTH) ----
    // 48 uniform bins over [-2.4, +2.4] (was 44 until 2026-08-25; see the note below). This is
    // the axis the pp/PbPb crossx histograms
    // h2d_crossx_*_pair_eta_binned_* and h3d_crossx_* are binned on, and therefore the axis
    // every 1D pair-eta view of the SAME quantity must use (.claude/CLAUDE.md Binnings rule 2:
    // 1D/2D/3D views of one quantity share one binning -- project the SAME histogram).
    // It is a DIFFERENT object from the 9 coarse analysis panels
    // CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, which are the physics cells.
    // Introduced 2026-08-25 to retire the retyped copies of "44, -2.4, 2.4"
    // (RDFBasedHistFillingPP.cxx x5+3).
    //
    // 2026-08-25: **44 -> 48** (user decision). Bin width becomes exactly 0.1, which makes every
    // boundary of the 9 coarse analysis panels
    // (CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap: +-0.5, +-1.0, +-1.5, +-2.0) an
    // EXACT bin edge. With 44 bins the width was 0.109090..., none of the 8 internal boundaries
    // was an edge, and every panel projection `FindBin(lo+eps)..FindBin(hi-eps)` therefore shared
    // bins 4, 9, 13, 18, 27, 32, 36, 41 with its neighbour: the 9 panels summed to 6942.93 pb
    // against a true total of 5784.06 pb (+20.0 %), and a panel labelled [-1.0,-0.5] actually
    // drew [-1.0909,-0.4364] -- a label/binning disagreement, i.e. .claude/CLAUDE.md Binnings
    // rule 5. Affected the pp24 crossx plot set (SingleBCrossxPlotterBase), the signal-acceptance
    // plots and the MC-data comparison alike. The TOTAL is unaffected (sigma_fid = 5784.06 pb);
    // only the per-panel split moves.
    // ADOPTED BY Pb+Pb ON 2026-09-08 (mu_pt45_gap125_pairpt9_adoption.md D1), together with the
    // fiducial + pair-level gap cuts, in the SAME step this comment used to demand.
    // RDFBasedHistFillingPbPb.cxx no longer retypes "44, -2.4, 2.4" in either spelling, and nor
    // do the Pythia/POWHEG truth signal-acceptance producers.
    // ** The CODE is migrated; the on-disk Pb+Pb histograms are NOT yet refilled. ** Until that
    // rerun completes, Pb+Pb panel plots and R_AA must not be quoted -- the guards in
    // Utilities/PairEtaPanelBins.h will throw on the stale files rather than mislabel them.
    static const int N_PAIR_ETA_CROSSX_BINS = 48;
    static constexpr double PAIR_ETA_CROSSX_MIN = -2.4;
    static constexpr double PAIR_ETA_CROSSX_MAX =  2.4;
    std::vector<double> pair_eta_crossx_bins;   // N_PAIR_ETA_CROSSX_BINS + 1 = 49 edges, filled in the constructor

    // ---- NOMINAL COARSE pair-pT binning (SINGLE SOURCE OF TRUTH) ----
    // The coarse pair-pT bins used for the low-mass template fit / R_AA coarse binning and
    // any coarsely-binned pair-pT plot set. THIS is the only place these values live: read
    // `pair_pt_coarse_bins` (and `N_COARSE_PAIR_PT_BINS`) from here, never re-invent per plot
    // so every plot set stays mutually consistent. Both the coarse and the fine-log
    // (pT_bins_150) pair-pT binnings may be re-optimised for the final analysis — change them
    // HERE ONLY.
    // The CANONICAL binning is **8 LOGARITHMIC bins from 9 to 150 GeV**, generated by
    // fillLogBinningArray in the constructor (never retyped).
    // PHYSICS: pair pT is the analysis's key observable. With only 4 coarse bins the pair-pT
    // dependence of the efficiency corrections is smeared out, and that smearing propagates
    // straight into the corrected pair-pT spectrum and into R_AA. 8 bins resolve it.
    // HISTORY (every previous definition is DELETED, by user instruction):
    //   2026-07-01  introduced as {8,15,27,50,150}
    //   2026-08-04  briefly the pT_bins_120 group edges {8, 13.75, 23.63, 40.62, 120}, with
    //               {8,15,27,50,150} retained as an opt-in `pair_pt_coarse_bins_pt150`
    //   2026-08-04  -> 8 log bins 8-150 GeV; `pair_pt_coarse_bins_pt150` REMOVED (it had no
    //               readers), and the coarse axis was no longer derived from the fine one, so
    //               NO interior edge coincided with it.
    //   2026-09-08  -> 8 log bins **9**-150 GeV, the low edge following the signal-region
    //               pair-pT cut 8 -> 9 GeV (user). At the same time the DEFAULT fine crossx
    //               axis became pT_bins_150 = 16 log bins over the SAME 9-150 range, so the
    //               coarse edges ARE fine edges again (coarse k = fine 2k, two fine bins per
    //               coarse cell). The 2026-08-04 warning above is therefore no longer true of
    //               pT_bins_150 -- but it REMAINS true of the `_pt_120` alternative axis.
    // If N changes, change N_COARSE_PAIR_PT_BINS only -- the edges follow automatically.
    static const int N_COARSE_PAIR_PT_BINS = 8;   // = pair_pt_coarse_bins.size() - 1
    std::vector<double> pair_pt_coarse_bins;      // CANONICAL edges, set in the constructor

    // OPT-IN COMPARISON VARIANT, never a second default: 4 logarithmic bins over the SAME
    // 9 -> 150 GeV range, for a like-for-like 4-vs-8 bin comparison of the MC trigger
    // efficiency Steps 2/3/4 with everything else held fixed (user, 2026-08-04). Selected
    // ONLY through Utilities/MCTrigEffPairPtBinning.h, which also supplies the matching
    // `_pt4bin` file token and `pt4bin/` plot subdirectory so a 4-bin run can never
    // overwrite or be confused with the nominal. Everything else reads
    // `pair_pt_coarse_bins`.
    static const int N_COARSE_PAIR_PT_BINS_4BIN = 4;
    std::vector<double> pair_pt_coarse_bins_4bin;

    // ---- NOMINAL COARSE single-muon pT binning (SINGLE SOURCE OF TRUTH) ----
    // Coarse binning in SINGLE-MUON pT (GeV), ~half the pair_pt_coarse_bins scale (each muon
    // carries ~half the pair pT); starts at the pt>4.5 GeV single-muon cut; last bin open-ended
    // (spectrum falls steeply). Used for single-muon |d0| / Delta-p/p distributions by
    // provenance (upfront hadronic/fake background reduction study, docs/tracking/
    // tf_upfront_bkg_reduction.md). Read from HERE, never re-invent per plot. Change HERE ONLY;
    // if N changes update BOTH the edges (constructor) and N_COARSE_SINGLE_MU_PT_BINS.
    static const int N_COARSE_SINGLE_MU_PT_BINS = 4;  // = single_mu_pt_coarse_bins.size() - 1
    std::vector<double> single_mu_pt_coarse_bins;     // edges, initialised in the constructor
	
	static std::vector<float> pTbins;
	static std::vector<double> pTHatbins_pythia;
    static std::vector<int> ctrbins;

    std::vector<std::string> pt_bin_labels;
    std::vector<std::string> pt_bin_exprs;
    std::vector<std::string> pt_titles = {"#bar{p}_T 4.5-5 GeV", "#bar{p}_T 5-6 GeV", "#bar{p}_T 6-7 GeV", "#bar{p}_T 7-8 GeV", "#bar{p}_T 8-9 GeV", "#bar{p}_T 9-10 GeV", "#bar{p}_T 10-12 GeV", "#bar{p}_T 12-15 GeV", "#bar{p}_T 15-20 GeV", "#bar{p}_T > 20 GeV"};
    
    std::vector<std::string> ctr_bin_labels;
    std::vector<std::string> ctr_bin_exprs;
	std::vector<std::string> ctr_titles = {"Centrality 0-5", "Centrality 5-10", "Centrality 10-20", "Centrality 20-30", "Centrality 30-50", "Centrality 50-80"};
	std::vector<std::string> ctrNpp_titles = {"Centrality 0-5", "Centrality 5-10", "Centrality 10-20", "Centrality 20-30", "Centrality 30-50", "Centrality 50-80", "pp"};


	std::string signs[nSigns] = {"ss","op"};
	std::vector<std::string> sign_labels = {"_ss","_op"};
	std::vector<std::string> sign_titles = {"same sign", "opposite sign"};

	std::vector<std::string> deta_cut_labels = {"_deta_lt0_8", "_deta_gt0_8"};
	std::vector<std::string> deta_titles = {"#Delta #eta < 0.8","#Delta #eta #geq 0.8"};

	std::vector<std::string> dphi_cut_labels = {"_near", "_away", "_dphi_lt0_6"};
	std::vector<std::string> dphi_titles = {"#Delta #phi < #pi/2","#Delta #phi #geq #pi/2", "#Delta #phi < 0.6"};

	std::vector<std::string> gapcut_labels = {"_nogapcut", "_wgapcut"};
	std::vector<std::string> gapcut_titles = {"no gap cut", "with gap cut"};

	std::vector<std::string> photocut_labels = {"_nophotocut", "_wphotocut"};
	std::vector<std::string> photocut_titles = {"no photoproduction cut", "with photoproduction cut"};


   	float deltaP_overP_thrsh; // single-muon dP/P cut
   	float deltaP_overP_max;
   	float deltaP_overP_step;
   	int deltaP_overP_nbins;
   	
   	float d0cut; // single-muon |d0| cut (unit: mm)
   	float z0cut; // single-muon |z0 sin(theta)| cut (unit: mm)

   	static float deltaR_thrsh[ndRselcs];
   	static float deltaR_thrsh_zoomin;
   	float 	deltaR_step = 0.01;
   	int deltaR_nbins[ndRselcs];

   	// float eta_gap_cut1 = 0.1;
   	static constexpr float eta_gap_cut1 = 0.135;
   	// float eta_gap_cut2[2] = {1.05,1.29};

	// // no log
   	// float minv_min[nSigns] = {0.2,1.06};
   	// float minv_max[ndRselcs] = {10,15,100};
   	// int minv_nbins[ndRselcs] = {100,150,1000};

   	float minv_max[ndRselcs] = {10,15,60};
   	// static constexpr int minv_nbins[ndRselcs] = {25,40,120};
   	static constexpr int minv_nbins[ndRselcs] = {20,30,120};
   		// the C++ standard does not specifiy how floating point should be implemented and is left to the processor. 
   		// To get around this and other limitations constexpr was introduced.
   	// float minv_logpow[nSigns][ndRselcs];
   	// std::vector<float> minv_bins[nSigns][ndRselcs];

   	std::vector<std::array<float,2>> minv_cuts;
   	std::vector<std::array<float,2>> minv_cuts_v2;
   	static std::vector<std::array<float,2>> charge_eta_gap_cuts;

   	// ---- single-muon detector-gap FIDUCIAL cut (q*eta) ----
   	// Reject a muon whose q*eta falls inside any of these windows. Compensated by an acceptance
   	// efficiency eps_acc in the correction chain (eps_acc itself is NOT built yet).
   	// STATUS: LIVE in the trigger-efficiency chain AND in the pp24 SIGNAL SELECTION. The two
   	// STATUS blocks at the end of this comment, not this header, are the current site list.
   	// Changing any window here changes the MC trig-eff, the data tag-and-probe efficiencies, the
   	// pp24 cross-section and everything downstream, and requires rerunning all of it
   	// (docs/signal_selection_change_impact.md). NOT a free-to-edit constant. See
   	// Analysis/docs/tracking/muon_gap_cuts_acceptance.md (F6/F7/F11/F17).
   	//
   	// Derived from the measured single-muon q*eta spectra (pp24 + PbPb 23/24/25):
   	//   {-0.10, +0.06} barrel crack at eta~0, ASYMMETRIC (user, 2026-08-12). The crack
   	//                  profile is NOT symmetric about 0: the minimum sits at q*eta ~ -0.02
   	//                  and the depletion is one-sided, in the SAME sense for both charges
   	//                  (+/- yield ratio at |q*eta| = 0.05 is 2.4 in pp / 2.0 in PbPb).
   	//                  Measured in 0.1-wide slices, N(slice)/N(mirror slice) is 0.51 (pp)
   	//                  / 0.58 (PbPb) for (-0.10, 0.00) -- the only depleted slice near
   	//                  zero -- while every slice from -1.0 to -0.2 is at or ABOVE 1.0,
   	//                  i.e. enriched. So the window is extended on the negative side to
   	//                  -0.10 to cover the whole depletion, and left at +0.06 on the
   	//                  positive side, where there is none. The earlier symmetric
   	//                  {-0.06, +0.06} cut only part of the depleted region.
   	//                  It is still narrower than the legacy eta_gap_cut1 = 0.135, which
   	//                  is symmetric in |eta| and so cuts a healthy +0.06..+0.135 slice.
   	//   {-1.25, -1.05} barrel/endcap transition. This is a q*eta (toroid bending
   	//                  direction) effect, NOT a fixed-|eta| geometric gap: each charge
   	//                  dips on ONE side only (mu+ at eta=-1.15, mu- at eta=+1.15, both
   	//                  q*eta=-1.15) with no mirror dip, a factor ~7-9 asymmetry. So it
   	//                  MUST be cut one-sided in q*eta; an |eta| window would discard
   	//                  twice the phase space for no gain. The window brackets the
   	//                  measured half-depth region (-1.17,-1.08). EDGE HISTORY: -1.20
   	//                  (2026-08-04) -> -1.30 (user, 2026-09-07, to cover the full
   	//                  negative-side extent) -> -1.25 (user, 2026-09-08). -1.30 reached
   	//                  ~0.13 beyond the measured half-depth into a region only 7-33%
   	//                  depleted (F6: -1.30..-1.22 runs 0.93 -> 0.67 of plateau), i.e. it
   	//                  cut away largely good acceptance; -1.25 still covers the dip. In pp
   	//                  the dip bottoms at 32% of plateau, in PbPb only 67% and centred
   	//                  nearer -1.16 -- a single common window is a pp/PbPb compromise.
   	//                  NOTE: -1.25 is NOT an edge of the FINE q*eta axis
   	//                  (makeEtaTrigEffcyBinning: ..., -1.26, -1.24, ... in 0.02 steps), so
   	//                  it splits bin [-1.26,-1.24). Accepted (user, 2026-09-08): the cut is
   	//                  applied to the muon's own q*eta value, and the CORRECTION uses the
   	//                  coarse contiguous bins, so only the shaded band on the fine-axis
   	//                  diagnostic figures and the 2D (q*eta,pT) map cell edge are affected.
   	//   {+2.20, +2.40} forward acceptance edge, ONE-SIDED: the yield collapses to 8.8%
   	//                  (pp) / 52% (PbPb) of its q*eta=2.2 value by 2.4, while the
   	//                  NEGATIVE side still holds ~71% out to -2.4. So there is
   	//                  deliberately no mirror window at -2.4.
   	//                  EDGE CHOICE (2.20, user 2026-09-07; was 2.30 from 2026-08-04):
   	//                  2.20 coincides with the historical standalone per-muon signal cut
   	//                  `m*.charge*m*.eta < 2.2` that this vector replaces, so the swap is
   	//                  exactly yield-neutral at the forward edge. The coarse turn-on
   	//                  binning's top edge MUST follow it down to 2.20
   	//                  (CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap), which is
   	//                  enforced by the startup throw in RDFBasedHistFillingData.
   	// Deliberately NOT included:
   	//   - the legacy "feet" window (0.56,0.67) and positive-side (1.064,1.29): the
   	//     measured structure there is too shallow to be worth the acceptance.
   	//   - any pT dependence: unlike charge_eta_gap_cuts, which only applies for
   	//     pT < 6 GeV, this is pT-INDEPENDENT so that the acceptance factors as a pure
   	//     function of q*eta.
   	// SCOPE: COMPLETE as of 2026-09-08. This vector is the SINGLE source of the gap/fiducial
   	// definition and has now REPLACED the standalone `q*eta < 2.2` in every live selection.
   	// The last five holdouts -- RDFBasedHistFillingPbPb.cxx (crossx + both template blocks),
   	// RDFBasedHistFillingPythiaTruth.cxx, the per-centrality nodes of
   	// RDFBasedHistFillingPythiaFullsimOverlay.cxx, and the Pythia/POWHEG truth signal
   	// ACCEPTANCE -- were migrated by decisions D1/D7 of
   	// docs/tracking/mu_pt45_gap125_pairpt9_adoption.md. A grep for a retyped `q*eta < 2.2` in
   	// RDFBasedHistFilling/, NTupleProcessingCode/ and Utilities/ now returns nothing live.
   	// ** The CODE is migrated; the affected on-disk histograms are NOT yet refilled. ** Any
   	// change to the windows themselves remains a signal-selection change with the full rerun
   	// blast radius in Analysis/docs/signal_selection_change_impact.md -- get the user's
   	// go-ahead first.
   	// COST, per muon. ** THE TABLES BELOW ARE FOR THE SUPERSEDED -1.30 SET AND ARE STALE FOR
   	// THE CURRENT -1.25 SET. ** They are kept because they are the last MEASURED values and
   	// because they bracket the current one from above (a narrower window can only cost less).
   	// Re-measured by plot_muon_q_eta_spectrum.cxx / plot_muon_truth_q_eta_spectrum.cxx as part
   	// of the 2026-09-08 rerun (docs/tracking/mu_pt45_gap125_pairpt9_adoption.md); update here
   	// when those numbers land. NOTE they are also stale in a SECOND way: they were measured on
   	// single-muon trees cut at pT > 4 GeV, and the threshold moved to 4.5 GeV in the same change.
   	// Measured 2026-09-07 for {-1.30,-1.05}, {-0.10,+0.06}, {2.20,2.40} on the then-current
   	// (pT > 4 GeV) single-muon trees (fraction of RECONSTRUCTED muons rejected):
   	//                 Tight     Medium
   	//   pp24          6.93%     7.03%
   	//   PbPb 0-80%    9.03%     9.77%
   	// (centrality-flat in PbPb: 9.00 / 9.07 / 9.04 / 9.06 / 9.15% over 0-10/10-20/20-30/30-50/
   	//  50-80%.) The p_T dependence is mild and monotonic in pp (5.30% in 4.0-4.5 GeV -> 8.00%
   	//  above 6 GeV) and flat in PbPb -- plot_muon_q_eta_pt_dependence.cxx.
   	// The earlier set {-1.20,-1.05},{-0.10,+0.06},{2.30,2.40} cost 3.79% (pp, per-window
   	// 1.62/1.80/0.37) / 5.60% (PbPb, 1.81/2.88/0.91), i.e. 7.44% / 10.89% of PAIRS. (The
   	// pair-level cost is not measured in either table; the PAIR-LEVEL |eta^pair| < 2.2 cut
   	// declared further down adds to it.)
   	// TRUTH-level cost per window -- what the acceptance eps_acc must carry -- from
   	// plot_muon_truth_q_eta_spectrum.cxx on the pp24 Pythia fullsim / PbPb overlay truth muons,
   	// again for the SUPERSEDED -1.30 set at truth pT > 4 GeV:
   	//                 [-1.30,-1.05]  [-0.10,+0.06]  [+2.20,+2.40]   eps_acc
   	//   pp24 fullsim      5.52%          4.16%          2.44%        0.8789
   	//   PbPb overlay      5.58%          4.15%          2.62%        0.8765
   	// (was eps_acc = 0.9133 / 0.9121 for the {-1.20,...,2.30} set.)
   	// STATUS 2026-08-04: WIRED IN, for the TRIGGER EFFICIENCY (user instruction).
   	//   - MC trig-eff sample: FillMCTrigEffHists.cxx, all of Steps 1-4, both legs of a pair.
   	//   - DATA tag-and-probe: applied to the PROBE only (user decision) -- eps^nc is a per-muon
   	//     efficiency and is only ever evaluated for muons outside the gaps, so eps(probe | probe
   	//     outside gap) is exactly the object the analysis applies. The tag is left uncut.
   	// STATUS 2026-08-17: WIRED IN as the pp24 SIGNAL SELECTION, on BOTH muons of a pair,
   	// REPLACING the standalone per-muon `q*eta < 2.2` (user instruction; tracking doc
   	// docs/tracking/pp24_crossx_rerun_2026_08.md; blast radius per
   	// docs/signal_selection_change_impact.md):
   	//   - pp24 data crossx: RDFBasedHistFillingPP.cxx `signal_cuts` + both `signal_cuts_no_minv`.
   	//   - pp24 fullsim reco efficiency: RDFBasedHistFillingPythiaFullsim.cxx
   	//     `pass_signal_truth` (TRUTH q*eta) and `pass_signal_reco` (RECO q*eta).
   	//   - the data-like mirror Utilities/MCTrigEffPairSelection.h SingleBSignalCutsReco().
   	// Because the ntuple stage already requires |eta| < 2.4, the forward window {2.20, 2.40}
   	// makes the effective forward edge 2.20 == the top edge of
   	// CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap, so every surviving muon has a
   	// fitted turn-on (no sentinel, no silent pair drop).
   	// STILL NOT APPLIED (future to-do, each with its own rerun):
   	// the PbPb data crossx and the PbPb overlay reco efficiency, the Pythia/Powheg TRUTH
   	// signal acceptance, the data ntuple processing, and the template fits.
   	// NOTE the resulting pp24 cross-section is a FIDUCIAL (gap-cut) cross-section: the
   	// truth-level gap acceptance eps_acc = 0.8789 (pp24 fullsim) / 0.8765 (PbPb overlay)
   	// (muon_gap_cuts_acceptance.md F12/F17) is a SEPARATE factor, not applied anywhere yet --
   	// and it does NOT yet include the cost of the PAIR-LEVEL |eta^pair| < 2.2 window declared
   	// below.
   	static std::vector<std::pair<float,float>> single_mu_fiducial_gap_cuts;

   	// --- The fiducial gap cut, in the TWO forms the analysis needs -----------------------
   	// A muon is REJECTED when q*eta = charge*eta falls inside any window above.
   	// Deliberately a NEW helper, not an extension of PassSingleMuonGapCut/MuPairPassGapCut
   	// (RDFBasedHistFillingData.cxx): those implement a DIFFERENT object (a symmetric
   	// |eta| < eta_gap_cut1 veto at all pT plus charge_eta_gap_cuts only below 6 GeV) and are the
   	// definition behind every existing `_wgapcut` diagnostic histogram -- redefining them would
   	// silently change the meaning of already-produced outputs whose names would not change.
   	// This cut is pT-INDEPENDENT by requirement, so the acceptance factorises in q*eta alone.
   	// The windows are rejected CLOSED, `[lo, hi]`, not open. That is not a detail:
   	// the forward window's lower edge (2.20) is ALSO the top edge of the contiguous coarse q*eta
   	// turn-on binning (CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap, whose bins are
   	// half-open `[lo, hi)`), and the ntuple keeps |eta| <= 2.4, so with OPEN windows a muon at
   	// exactly q*eta = 2.20 or 2.40 survived the cut and then had NO fitted turn-on -- which
   	// EvaluateSingleMuonEffcyPtFitted throws on, by design, rather than silently dropping the
   	// pair. Closing the interval makes the surviving region exactly [-2.4, 2.20) minus the two
   	// interior windows, i.e. precisely the region the turn-on fits cover. (Found 2026-08-17 when
   	// a real pp24 muon landed on q*eta = 2.300; the change is measure-zero everywhere else,
   	// since it only moves exact boundary values.)
   	static bool PassSingleMuFiducialGap(float eta, int charge) {
   		const float q_eta = charge * eta;
   		for (const auto& w : single_mu_fiducial_gap_cuts)
   			if (q_eta >= w.first && q_eta <= w.second) return false;
   		return true;
   	}
   	// RDF/JIT string form, so the window numbers are never retyped into a Filter expression.
   	// `q_eta_expr` is any expression evaluating to q*eta (e.g. "charge*eta", "lg_charge*lg_eta").
   	//
   	// EVERYTHING HERE IS IN FLOAT, deliberately: the expression is cast to float and the edges
   	// carry an `f` suffix. Without that, the JIT compares a float q*eta PROMOTED TO DOUBLE
   	// against a DOUBLE decimal literal -- and 2.30f = 2.2999999523... is strictly less than the
   	// double 2.3 = 2.2999999999..., so a muon sitting exactly on the forward edge was NOT
   	// rejected, while `FindBinReturnStr` (which compares in float, against the float bin edge
   	// 2.30f) then found no q*eta bin for it and EvaluateSingleMuonEffcyPtFitted threw. Observed
   	// on real pp24 data, 2026-08-18. The float cast makes this Filter and the float-valued
   	// binning agree bit for bit -- and makes it agree with PassSingleMuFiducialGap, which has
   	// always compared in float.
   	static std::string FiducialGapCutExpr(const std::string& q_eta_expr) {
   		std::string s = "!(";
   		for (size_t i = 0; i < single_mu_fiducial_gap_cuts.size(); ++i) {
   			const auto& w = single_mu_fiducial_gap_cuts[i];
   			if (i) s += " || ";
   			s += "((float)(" + q_eta_expr + ") >= " + std::to_string(w.first) + "f"
   			   + " && (float)(" + q_eta_expr + ") <= " + std::to_string(w.second) + "f)";
   		}
   		return s + ")";
   	}

   	// ---- PAIR-LEVEL detector-gap fiducial cut: |eta^pair| < 2.2 --------------------------
   	// The single-muon fiducial cut above removes q*eta in [2.20, 2.40], ONE-SIDED. A PAIR can
   	// still reach |eta^pair| > 2.2 with both muons passing that cut -- e.g. two muons at
   	// eta = +2.3 whose charges both give q*eta = -2.3. But in that region the pair acceptance is
   	// carved out by the SINGLE-MUON q*eta window in a charge- and configuration-dependent way, so
   	// the pair efficiency there is a complicated function of the single-muon cut rather than a
   	// smooth detector response. Rather than model that, the region is removed (user instruction,
   	// 2026-09-07): |eta^pair| < 2.2, SYMMETRIC -- unlike the single-muon windows, which are
   	// one-sided in q*eta because they track the toroid bending direction; eta^pair carries no
   	// charge, so there is nothing for a one-sided window to track.
   	// This is a GAP CUT AT PAIR LEVEL and belongs everywhere the single-muon gap cut is applied to
   	// both legs of a pair. Like the single-muon windows, its cost must eventually be carried by
   	// the acceptance factor eps_acc (muon_gap_cuts_acceptance.md F12).
   	// Comparison is STRICT (`<`), so the surviving region is the OPEN interval (-2.2, +2.2), and
   	// it is done in FLOAT for the same JIT float/double reason documented above. Strict is the
   	// correct direction here: the coarse pair-eta bins are half-open [lo,hi) with outer edges at
   	// exactly +-2.2, so a pair sitting on 2.2 must be rejected or it would land in overflow.
   	// NOTE PassPairFiducialEta is the non-RDF (plain C++) form and currently has NO call sites --
   	// every live consumer is an RDF Filter and uses PairFiducialEtaCutExpr. It is kept as the
   	// predicate a non-RDF consumer should call rather than re-implementing the test.
   	static constexpr float pair_eta_fiducial_max = 2.2f;
   	static bool PassPairFiducialEta(float pair_eta) {
   		return std::fabs(pair_eta) < pair_eta_fiducial_max;
   	}
   	// RDF/JIT string form; `pair_eta_expr` is any expression evaluating to the pair eta.
   	static std::string PairFiducialEtaCutExpr(const std::string& pair_eta_expr) {
   		return "(fabs((float)(" + pair_eta_expr + ")) < "
   		     + std::to_string(pair_eta_fiducial_max) + "f)";
   	}

   	// ---- SIGNAL-REGION pair-pT threshold (SINGLE SOURCE OF TRUTH) -------------------------
   	// 8.0 -> 9.0 GeV (user, 2026-09-08), together with the single-muon cut 4 -> 4.5 GeV.
   	// PHYSICS: with both muons above 4.5 GeV the region just above 8 GeV is a sculpted corner
   	// of phase space whose acceptance is set by the single-muon cut rather than by a smooth
   	// detector response. It is removed rather than corrected.
   	// This value was RETYPED as a bare `8` at ~10 sites (data crossx signal region, the
   	// template blocks, the generic family, and every truth/reco MC analog) until 2026-09-08.
   	// Read it from HERE -- never retype it -- so data and the MC analogs cannot drift apart.
   	// COUPLED: the canonical pair-pT axes start at this value (pair_pt_coarse_bins,
   	// pT_bins_150, pT_bins_120, pT_bins_80 and the 4-bin variant), so changing it without
   	// moving them leaves a partly-empty first bin. Full blast radius:
   	// docs/signal_selection_change_impact.md.
   	static constexpr float signal_pair_pt_min = 9.0f;
   	static bool PassSignalPairPt(float pair_pt) { return pair_pt > signal_pair_pt_min; }
   	// RDF/JIT string form; `pair_pt_expr` is the pair-pT variable (e.g. "pair_pt" for reco,
   	// "truth_pair_pt" for the truth analogs).
   	static std::string SignalPairPtCutExpr(const std::string& pair_pt_expr = "pair_pt") {
   		return "((float)(" + pair_pt_expr + ") > "
   		     + std::to_string(signal_pair_pt_min) + "f)";
   	}

   	// ---- SIGNAL-REGION dimuon MASS WINDOW (SINGLE SOURCE OF TRUTH) --------------------------
   	// 1.08 < m_mumu < 2.9 GeV (docs/analysis_overview.md §2): above the phi and its radiative
   	// tail, below the J/psi. Until 2026-09-17 it was RETYPED as a bare `minv > 1.08 && minv < 2.9`
   	// at ~12 sites (data crossx signal region, every truth/reco MC analog, the MC trig-eff
   	// signal selection, the single-value pair trigger efficiency's `sig` window). Read it from
   	// HERE -- never retype it.
   	// COUPLED CONSUMERS (docs/tracking/pp24_trig_eff_hybrid_application.md D3/D5): the
   	// single-value pair trigger efficiency (Utilities/PairTrigEffEvaluator.h) is MEASURED inside
   	// this window and applied to pairs above ParamsSet::pair_pt_coarse_bins[N-2]; the low-mass
   	// template fit is to be performed inside this SAME window (future work). Changing the window
   	// therefore requires re-measuring the pair efficiency (run_mc_trigeff_pair_eff.sh) and
   	// refilling every crossx / MC analog -- docs/signal_selection_change_impact.md.
   	// DOUBLES, formatted with %g in the cut expression: that reproduces the literal strings the
   	// sites carried before ("1.08", "2.9"), so the JIT selection is byte-identical and every
   	// output bit-identical to the retyped era (the float form "1.080000f" would not be).
   	static constexpr double signal_minv_min = 1.08;
   	static constexpr double signal_minv_max = 2.9;
   	static bool PassSignalMinv(double minv) { return minv > signal_minv_min && minv < signal_minv_max; }
   	// RDF/JIT string form; `minv_expr` is the mass variable ("minv" for reco, "truth_minv" for
   	// the truth analogs). Open interval, exactly as every site wrote it.
   	static std::string SignalMinvCutExpr(const std::string& minv_expr = "minv") {
   		char buf[128];
   		snprintf(buf, sizeof(buf), "%s > %g && %s < %g",
   		         minv_expr.c_str(), signal_minv_min, minv_expr.c_str(), signal_minv_max);
   		return std::string(buf);
   	}

   	float minv_upper = 60;

  	static double PI;

  	static AxisInfo ptlead_axis;
	static AxisInfo dr_axis;
	static AxisInfo dphi_axis;


  	static const int npt_bins = 50;
  	float ptlogpow = 0.02194;
  	float ptmax = 50;
  	double pTBins[npt_bins+1];
	
  	static const int npairPT_bins = 40;
  	float pairPTlogpow[nSigns][ndRselcs];
  	float pairPTmax = 80;
  	double pairPTBins[nSigns][ndRselcs][npairPT_bins+1];

  	static const int n_hq_pt_bins = 40;
  	float hq_ptlogpow = 0.034;
  	float hq_ptmax = 110;
  	double hq_pTBins[n_hq_pt_bins+1];

  	static const int n_hq_minv_bins = 40;
  	float hq_minvlogpow = 0.049;
  	float hq_minvmax = 220;
  	double hq_minvBins[n_hq_minv_bins+1];
  	
  	ParamsSet();
  	~ParamsSet(){}
    static std::vector<double> makeEtaTrigEffcyBinning(int rebin_factor = 1);
  	void fillLogBinningArray(std::vector<double>& bins, int nBins, double low, double high);
  	template <typename T>
  	std::string write_single_bin_expr (std::string kinvar, T a, T b);
  	template <typename T>
	std::string write_single_bin_label (std::string kinvar, T a, T b);
	template <typename T>
	void write_cut_label_strs (std::string kinvar_expr, std::string kinvar_label, const std::vector<T> & var_bin_bdrys, std::vector<std::string> & var_bin_exprs, std::vector<std::string> & var_bin_labels, bool open_ended);
};

void ParamsSet::fillLogBinningArray(std::vector<double>& bins, int nBins, double low, double high) {
    bins.clear();

    double logLow = std::log10(low);
    double logHigh = std::log10(high);
    double logStep = (logHigh - logLow) / nBins;

    for (int i = 0; i <= nBins; ++i) {
        bins.push_back(std::pow(10, logLow + i * logStep));
    }
}

#include <vector>
#include <algorithm>
#include <cmath>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

std::vector<double> ParamsSet::makeEtaTrigEffcyBinning(int rebin_factor)
{
    const double minEdge = -2.4;
    const double maxEdge =  2.4;

    // Base steps (before rebinning)
    const double ultraFineBase = 0.01;  // adjust to 0.005 if you really want that
    const double fineBase      = 0.02;
    const double coarseStep    = 0.10;

    const double ultraFineStep = ultraFineBase * rebin_factor;
    const double fineStep      = fineBase      * rebin_factor;

    // Fine ranges (to be merged)
    std::vector<std::pair<double,double>> fineRanges = {
        {-1.3,  1.0},
        {-0.8,  1.4},
        { 2.2,  2.4}
    };

    // Ultra-fine range
    const std::pair<double,double> ultraFineRange = {-0.2, 0.2};

    // --- merge overlapping fine ranges ---
    std::sort(fineRanges.begin(), fineRanges.end());
    std::vector<std::pair<double,double>> merged;
    for (auto &r : fineRanges) {
        if (merged.empty() || r.first > merged.back().second) {
            merged.push_back(r);
        } else {
            merged.back().second = std::max(merged.back().second, r.second);
        }
    }

    std::vector<double> edges;
    edges.reserve(200); // just to avoid reallocations
    const double eps = 1e-10;

    auto push_edge = [&](double x) {
        // round to kill FP noise
        x = std::round(x * 1e12) / 1e12;
        if (edges.empty() || std::fabs(edges.back() - x) > eps) {
            edges.push_back(x);
        }
    };

    // Add segment [start, end] with given step
    auto add_segment = [&](double start, double end, double step) {
        if (end <= start + eps) {
            push_edge(start);
            return;
        }
        push_edge(start);
        double x = start;
        while (x + step < end - eps) {
            x += step;
            push_edge(x);
        }
        push_edge(end);
    };

    double x = minEdge;

    // Walk over all merged fine ranges in order
    for (const auto &fr : merged) {
        const double a = fr.first;
        const double b = fr.second;

        // 1) Coarse region up to the start of this fine range
        if (x < a - eps) {
            add_segment(x, a, coarseStep);
            x = a;
        }

        // 2) Inside the fine range [a, b]
        double uf_lo = std::max(a, ultraFineRange.first);
        double uf_hi = std::min(b, ultraFineRange.second);

        // 2a) Fine region before ultra-fine part
        if (x < uf_lo - eps) {
            add_segment(x, uf_lo, fineStep);
            x = uf_lo;
        }

        // 2b) Ultra-fine region (overlap)
        if (uf_lo < uf_hi - eps) {
            add_segment(x, uf_hi, ultraFineStep);
            x = uf_hi;
        }

        // 2c) Fine region after ultra-fine part within this fine range
        if (x < b - eps) {
            add_segment(x, b, fineStep);
            x = b;
        }
    }

    // 3) After the last fine range: coarse step to maxEdge
    if (x < maxEdge - eps) {
        add_segment(x, maxEdge, coarseStep);
    }

    return edges;
}


double ParamsSet::PI = acos(-1.0);

float ParamsSet::deltaR_thrsh[ndRselcs] = {0.8,1.2,5.75};
float ParamsSet::deltaR_thrsh_zoomin = 0.8;

std::vector<std::array<float,2>> ParamsSet::charge_eta_gap_cuts = {{0.56,0.67}, {1.064,1.29}, {-1.29,-1.12}};
// Fiducial gap windows in q*eta -- see the declaration for the derivation.
// HISTORY of the two windows that have moved:
//   barrel/endcap: {-1.20,-1.05} -> {-1.30,-1.05} (user, 2026-09-07) -> {-1.25,-1.05}
//   (user, 2026-09-08; -1.30 reached ~0.13 beyond the measured half-depth into a region only
//   7-33 % depleted, i.e. it cut away largely good acceptance). Widened on the
//     negative-q*eta side to cover the full extent of the one-sided transition dip.
//   forward edge: 2.20 -> 2.30 (user, 2026-08-04, on the edge scan of
//     mc_trigger_efficiency.md R22) -> back to 2.20 (user, 2026-09-07). At 2.20 the
//     forward window again coincides with the historical standalone per-muon
//     `q*eta < 2.2` signal cut, so the fiducial swap is yield-neutral there.
// COUPLED: CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap's top bin must END at this
// same 2.20, or muons between the two edges would survive the cut with no fitted turn-on. The
// coupling is not left to this comment: RDFBasedHistFillingData::SetIOPaths checks it at startup
// and THROWS if the two numbers disagree.
std::vector<std::pair<float,float>> ParamsSet::single_mu_fiducial_gap_cuts = {{-1.25f,-1.05f}, {-0.10f,0.06f}, {2.20f,2.40f}};
// Pair-average-pT bin edges. Low edge follows the single-muon pT cut: 4 -> 4.5 GeV
// (user, 2026-09-08) so the first bin is not left partly empty.
std::vector<float> ParamsSet::pTbins = {4.5,5.,6.,7.,8.,9.,10.,12.,15.,20.};
std::vector<int> ParamsSet::ctrbins = {0, 5, 10, 20, 30, 50, 80};
// std::vector<double> ParamsSet::pTHatbins_pythia = {5, 10, 25, 60, 120, 3200};
std::vector<double> ParamsSet::pTHatbins_pythia = {4.5, 10, 25, 60, 120, 3200}; // lower to 4.5 since a few events have pTHat < 5GeV (leakage)

AxisInfo ParamsSet::ptlead_axis = {100, 0, 30., nullptr};
AxisInfo ParamsSet::dr_axis = {100, 0, ParamsSet::deltaR_thrsh[2], nullptr};
AxisInfo ParamsSet::dphi_axis = {32, -ParamsSet::PI/2., ParamsSet::PI * 3. / 2., nullptr};

// Named helper functions instead of lambdas: cling cannot JIT-link lambdas
// defined at file scope because closure types lack external linkage.
inline bool ParamsSet_dphi_near  (float x){ return std::abs(x) < 1.5708f; }
inline bool ParamsSet_dphi_away  (float x){ return std::abs(x) >= 1.5708f; }
inline bool ParamsSet_dphi_lt0_6 (float x){ return std::abs(x) < 0.6f; }
inline bool ParamsSet_deta_lt0_8 (float x){ return std::abs(x) < 0.8f; }
inline bool ParamsSet_deta_ge0_8 (float x){ return std::abs(x) >= 0.8f; }

std::vector<std::function<bool(float)>> ParamsSet::dphi_cut_funcs {ParamsSet_dphi_near, ParamsSet_dphi_away, ParamsSet_dphi_lt0_6};
std::vector<std::function<bool(float)>> ParamsSet::deta_cut_funcs {ParamsSet_deta_lt0_8, ParamsSet_deta_ge0_8};
const unsigned int ParamsSet::nCtrBins; // number of coarse bins; each can be studied by itself
const unsigned int ParamsSet::nPtBins;
const unsigned int ParamsSet::nPythiaKinRanges;
const unsigned int ParamsSet::nSigns;
const unsigned int ParamsSet::ndRselcs;
const unsigned int ParamsSet::ndphiselcs;
const unsigned int ParamsSet::ndphiRegions;
const unsigned int ParamsSet::ndetaRegions;
const unsigned int ParamsSet::nPhotoProdCuts; //0: no gap cut; 1: having gap cut
const unsigned int ParamsSet::nGapCuts; //0: no gap cut; 1: having gap cut
const unsigned int ParamsSet::CtrStep;
const unsigned int ParamsSet::nCtrIntvls; // number of small intervals; the intervals not to be studied by themselves, but to allow easy combinations

template <typename T>
std::string ParamsSet::write_single_bin_expr (std::string kinvar, T a, T b){
    return kinvar + " >= " + std::to_string(a) + " && " + kinvar + " < " + std::to_string(b);
}

template <typename T>
std::string ParamsSet::write_single_bin_label (std::string kinvar, T a, T b){
    return "_" + kinvar + std::to_string(static_cast<int>(a)) + "to" + std::to_string(static_cast<int>(b));
}

template <typename T>
void ParamsSet::write_cut_label_strs (std::string kinvar_expr, std::string kinvar_label, const std::vector<T> & var_bin_bdrys, std::vector<std::string> & var_bin_exprs, std::vector<std::string> & var_bin_labels, bool open_ended){
    var_bin_exprs.clear();
    var_bin_labels.clear();
    for (typename std::vector<T>::const_iterator it = var_bin_bdrys.begin() + 1; it < var_bin_bdrys.end(); it++){
        var_bin_exprs.push_back(write_single_bin_expr(kinvar_expr, *(it-1), *it));
        var_bin_labels.push_back(write_single_bin_label(kinvar_label, *(it-1), *it));
    }

    if (open_ended){
        var_bin_exprs.push_back(kinvar_expr + " >= " + std::to_string(var_bin_bdrys[var_bin_bdrys.size()-1]));
        var_bin_labels.push_back("_" + kinvar_label + "_above" + std::to_string(static_cast<int>(var_bin_bdrys[var_bin_bdrys.size()-1])));
    }
}

ParamsSet::ParamsSet(){
	deltaP_overP_thrsh = 0.12;
 	deltaP_overP_max = 0.12 * sqrt(2);
 	deltaP_overP_step = 0.002;
 	deltaP_overP_nbins = static_cast<int>(deltaP_overP_max/deltaP_overP_step);

	d0cut = 2; // 2mm
	z0cut = 2; // 2mm

// ------------------------------------ pt & ctr labels & expressions ------------------------------------
	    
    // write_cut_label_strs("(m1.pt + m2.pt) / 2.", "ptavg", pTbins, pt_bin_exprs, pt_bin_labels, true);
    write_cut_label_strs("pt_avg", "ptavg", pTbins, pt_bin_exprs, pt_bin_labels, true);
    write_cut_label_strs("avg_centrality", "ctr", ctrbins, ctr_bin_exprs, ctr_bin_labels, false);

  	for (unsigned int idr = 0; idr < ndRselcs; idr++){
    	deltaR_nbins[idr] = static_cast<int>(deltaR_thrsh[idr]/deltaR_step);
  	}

	for(int i = 0; i <= npt_bins; i++){
    	pTBins[i] = ptmax * pow(10.0, (static_cast<float>(i-npt_bins))*ptlogpow);
  	}

	for(int i = 0; i <= n_hq_pt_bins; i++){
    	hq_pTBins[i] = hq_ptmax * pow(10.0, (static_cast<float>(i-n_hq_pt_bins)) * hq_ptlogpow);
  	}

	for(int i = 0; i <= n_hq_minv_bins; i++){
    	hq_minvBins[i] = hq_minvmax * pow(10.0, (static_cast<float>(i-n_hq_minv_bins)) * hq_minvlogpow);
  	}

	pairPTlogpow[0][0] = 0.0198;
	pairPTlogpow[0][1] = 0.0198;
	// pairPTlogpow[0][2] = 0.052;
	pairPTlogpow[0][2] = 0.078;
	pairPTlogpow[1][0] = 0.0198;
	pairPTlogpow[1][1] = 0.0198;
	pairPTlogpow[1][2] = 0.068;


	// ={{0.0198,0.0198,0.052},{0.0198,0.0198,0.082}};
  	for (int isign = 0; isign < nSigns; isign++){
  		for (int idr = 0; idr < ndRselcs; idr++){
  			for(int ipt = 0; ipt <= npairPT_bins; ipt++){
    			pairPTBins[isign][idr][ipt] = pairPTmax * pow(10.0, (static_cast<float>(ipt - npairPT_bins)) * pairPTlogpow[isign][idr]);
  			}
  		}
  	}

    // ---- pT axes. Every LOW EDGE tracks the cut that defines it (user, 2026-09-08) ----
    // SINGLE-MUON axes start at the single-muon cut, 4 -> 4.5 GeV.
    // PAIR-pT axes start at the signal-region pair-pT cut,  8 -> 9 GeV.
    // Bin COUNTS and the log-spacing rule are unchanged except where noted.
    fillLogBinningArray(pT_bins_40,  18, 4.5,  40.0);  // 18 log bins from 4.5 to 40  GeV (single mu)

    // pT_bins_8 -- the low half of pT_bins_single_muon, i.e. the trigger-efficiency `pt2nd`
    // MEASUREMENT axis. It is deliberately NOT an analysis axis: `single_mu_pt_coarse_bins` and
    // `pTbins` carry the 4.5 GeV analysis threshold, this one does not.
    // TWO SEGMENTS, with 4.5 FORCED as an interior bin edge (user, 2026-09-08, D10):
    //   * the DATA tag-and-probe keeps probes down to 4.0 GeV (the trigger-efficiency NTuple mode
    //     is cut at 4.0 while the nominal analysis mode is cut at 4.5 -- D8), because eps^nc is a
    //     PER-MUON efficiency: its measurement population need not equal the population it is
    //     applied to, and cutting probes at 4.5 would delete ~40 % of the 4->6 GeV mu4 rise and
    //     leave the turn-on midpoint outside the fitted range in 63-80 % of q*eta cells,
    //     degenerate in (mean, sigma) with NO existing warning able to detect it.
    //   * MC has nothing below 4.5 at all (D9), so it must START on a bin edge. Were 4.5 to fall
    //     inside a bin, MC would get one ~29 %-filled bin whose graph point sits at the bin CENTRE
    //     while its survivors' mean is higher -- an x-shift exactly where the turn-on is steepest.
    // Total bin count is unchanged (20), and above 4.5 the data and MC maps stay bit-identical,
    // so the Step-1 MC/data ratio is unaffected.
    {
        std::vector<double> seg_lo, seg_hi;
        fillLogBinningArray(seg_lo,  2, 4.0, 4.5);   //  2 log bins 4.0 -> 4.5 (DATA probes only)
        fillLogBinningArray(seg_hi, 18, 4.5, 8.0);   // 18 log bins 4.5 -> 8.0 (data AND MC)
        pT_bins_8 = seg_lo;
        pT_bins_8.pop_back();                        // drop the duplicated 4.5 edge
        pT_bins_8.insert(pT_bins_8.end(), seg_hi.begin(), seg_hi.end());
    }
    fillLogBinningArray(pT_bins_60,  20, 8.0,  60.0);  // 20 log bins from 8  to 60  GeV (single mu,
                                                       //   high half; 8 is not a cut -> unchanged)
    fillLogBinningArray(pT_bins_80,  12, 9.0,  80.0);  // 12 log bins from 9  to 80  GeV (pair pT)
    // ---- FINE pair-pT axes of the cross-section views ----
    // pT_bins_150 is the NOMINAL/default crossx axis (user, 2026-09-08): 16 log bins 9 -> 150,
    // i.e. EXACTLY TWO fine bins per coarse bin of pair_pt_coarse_bins (8 log bins 9 -> 150),
    // so every coarse edge IS a fine edge (coarse edge k = fine edge 2k). That nesting is the
    // whole reason for 16 rather than 15, and it is what stops a panel/cell projection from
    // straddling bins the way the 44-vs-48 pair-eta axis did (see N_PAIR_ETA_CROSSX_BINS).
    // pT_bins_120 is the OPT-IN `_pt_120` alternative view, same 16 bins over the shorter
    // 9 -> 120 reach. It deliberately does NOT nest with the coarse cells -- it is a display
    // variant only and must never be used to bin a correction.
    fillLogBinningArray(pT_bins_120, 16, 9.0, 120.0);  // 16 log bins from 9  to 120 GeV (alternative)
    fillLogBinningArray(pT_bins_150, 16, 9.0, 150.0);  // 16 log bins from 9  to 150 GeV (DEFAULT)

    // FINE pair-eta axis of the crossx 2D/3D views -- see the declaration. Generated, so the
    // edges follow N_PAIR_ETA_CROSSX_BINS automatically and are never retyped.
    pair_eta_crossx_bins.clear();
    pair_eta_crossx_bins.reserve(N_PAIR_ETA_CROSSX_BINS + 1);
    for (int i = 0; i <= N_PAIR_ETA_CROSSX_BINS; ++i)
        pair_eta_crossx_bins.push_back(PAIR_ETA_CROSSX_MIN
            + (PAIR_ETA_CROSSX_MAX - PAIR_ETA_CROSSX_MIN)
              * static_cast<double>(i) / N_PAIR_ETA_CROSSX_BINS);

    // NOMINAL COARSE pair-pT binning (single source of truth; see the declaration above).
    // 8 LOGARITHMIC bins, 9 -> 150 GeV. Generated, never retyped, so N_COARSE_PAIR_PT_BINS is
    // the only number to change. Low edge 8 -> 9 with the signal-region pair-pT cut
    // (user, 2026-09-08).
    // NESTING (restored 2026-09-08, and now deliberate): the default fine crossx axis
    // pT_bins_150 is 16 log bins over the SAME 9 -> 150 range, so coarse edge k IS fine edge 2k
    // -- exactly two fine bins per coarse cell. This REVERSES the 2026-08-04 note that "nothing
    // may assume the coarse edges are a subset of the fine ones": for pT_bins_150 they now are.
    // It is still FALSE for the `_pt_120` alternative axis (9 -> 120), which does not nest.
    fillLogBinningArray(pair_pt_coarse_bins, N_COARSE_PAIR_PT_BINS, 9.0, 150.0);
    // Opt-in 4-bin comparison variant: SAME range and SAME log spacing rule, so the only
    // difference from the nominal is the number of cells (its edges are the nominal ones
    // at indices 0, 2, 4, 6, 8 -- every 4-bin edge is also an 8-bin edge, which is what
    // makes the comparison a clean merge rather than a re-binning).
    fillLogBinningArray(pair_pt_coarse_bins_4bin, N_COARSE_PAIR_PT_BINS_4BIN, 9.0, 150.0);

    // NOMINAL COARSE single-muon pT binning (single source of truth; see the declaration above).
    // 4 bins: [4.5,8), [8,14), [14,25), [25,100]. ~half the pair_pt_coarse_bins scale; starts at
    // the single-muon pT cut (4 -> 4.5 GeV, user 2026-09-08); last bin open-ended. Keep the edge
    // count consistent with N_COARSE_SINGLE_MU_PT_BINS.
    single_mu_pt_coarse_bins = {4.5, 8.0, 14.0, 25.0, 100.0};


  	// minv cut V1 - cut off all below 1.06GeV
  	minv_cuts.push_back({0,1.06});
   	minv_cuts.push_back({2.9,3.3});
   	minv_cuts.push_back({3.55,3.8});
   	minv_cuts.push_back({9.08,10.5}); // previously 9 - 9.8

   	// minv cut V2 - cut off narrower windows for individual sub-GeV light resonance peaks
  	minv_cuts_v2.push_back({0.,0.6});
  	minv_cuts_v2.push_back({0.72,0.85});
  	minv_cuts_v2.push_back({0.94,1.06});
   	minv_cuts_v2.push_back({2.9,3.3});
   	minv_cuts_v2.push_back({3.55,3.8});
   	minv_cuts_v2.push_back({9.08,10.5}); // previously 9 - 9.8
}

#endif
