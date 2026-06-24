# Analysis Code Review Log
**Task**: Rewrite ScrambGen (mixed-event combinatoric pair generator) to the object model — read MuonObj single-muon trees, mix across events within 5% ctr intervals, write muon_pair_tree_sign1/sign2 (MuonPairObj). Produces T_mix for the 5b closure + fit.
**Log file**: review-analysis-code-20260623-202339-scrambgen-objectmodel-rewrite.md
**Started**: 2026-06-24T00:23:39Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Executor work (pre-review)
Files written (replacing old non-compiling content):
- ScrambGen/ScrambGen.h (header-only PbPb class `ScrambGen`): LoadMuons (TChain over per-part single-muon
  trees, bin by ev_centrality/CtrStep into nCtrIntvls=20, exclude ev_centrality<0); GeneratePairs (per
  interval, target=oversample(5)*N_muons; pick i, pick j with ev_num[j]!=ev_num[i], cap 50 retries; build
  MuonPairPbPb m1/m2/year/weight=1/Update(); route same_sign->sign1(SS), else sign2(OS)); Run(year). Fixed
  seed TRandom3. Includes MuonPairPbPb.h + ParamsSet.h.
- ScrambGen/ScrambGen.c: #include "ScrambGen.h" (ACLiC compile unit + dict).
- ScrambGen/ScrambGenPP.h (PP class `ScrambGenPP`): single pool (no centrality), MuonPP/MuonPairPP, 12 parts,
  same mixing + sign routing, no year/avg_centrality. ScrambGen/ScrambGenPP.c: include.
- ScrambGen/run_scrambgen.sh: builds + runs PbPb 23/24/25 then pp (separate ROOT sessions).
Tree/branch names: muon_pair_tree_sign1 (SS) / muon_pair_tree_sign2 (OS), branch MuonPairObj — EXACTLY what
the RDF reads (RDFBasedHistFillingBaseClass.h:82-83). NO physics cuts in ScrambGen (RDF applies signal_cuts).
Compilation: ACLiC clean both (separate sessions): PbPb SCRAMB_PBPB_OK, PP SCRAMB_PP_OK (after adding the
explicit ParamsSet.h include — MuonPairPbPb.h did not pull it transitively).
Run: IN PROGRESS (background byjitwg6l; ~100M scrambled pairs with object I/O). Verify outputs + the key
physics check (OS minv = smooth combinatoric continuum, NO resonance peaks) before review.

Run RESULTS (all exit 0): yr23 2.82M muons -> SS=7.046M/OS=7.049M; yr24 1.99M -> SS=4.972M/OS=4.971M;
yr25 5.63M -> SS=14.075M/OS=14.073M; pp 9.26M muons -> output written. SS≈OS every year (combinatoric is
charge-symmetric). Output files: pbpb 1.9/1.4/3.7 GB, pp 1.7 GB (object trees). Both trees present + non-empty.
PHYSICS VERIFICATION (pbpb23 OS minv 0-4): mean=2.656; **NO resonance peaks** — J/psi[3.0,3.2]/sidebands
ratio=0.997, phi[1.0,1.06]/sidebands=0.968 (both ~1 = smooth combinatoric continuum). This is the defining
correctness check for a mixed-event sample (uncorrelated events -> no real J/psi/phi). Mixing correct.
Numbers reported: per-year SS/OS counts above; pbpb23 OS J/psi-ratio=0.997, phi-ratio=0.968 (no peaks).

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 1 CRITICAL (+1 INFO)
**Details**: CRITICAL [PHYSICS/LOGIC] yr25 scrambled pairs get avg_centrality=-1 for ALL pairs.
Root cause (reviewer-diagnosed): for year==25, MuonPairPbPb::Update()->PairValueCalcHook calls
UpdateCentrality() which recomputes centrality from the PAIR-level FCal_Et; ScrambGen never sets pair->FCal_Et
(mixed pair has no single event FCal) -> stays -1e6 -> GetCentralityPbPb2023(-1e6)=-1 -> RDF centrality filter
drops ALL yr25 pairs. yr23/24 correct (ev_centrality-average path). INFO: MuonPairPbPb.h:51 yr25
"centrality all zeros -> recalc from FCal" premise is stale (yr25 single-muon ev_centrality now valid, mean
15.48) — flag to user, outside ScrambGen scope.
**Numerical verification**: OS/SS entries, OS minv mean 2.656, J/psi-ratio 0.997, phi-ratio 0.968 all MATCH
(no resonance peaks — mixing correct); yr25 avg_centrality MISMATCH (all -1) → the bug.

## Iteration 1 fix
Added explicit `pair->avg_centrality = (mu[i].ev_centrality + mu[j].ev_centrality)/2;` AFTER Update()
(ScrambGen.h GeneratePairs). This is the original Design-Decision resolution #4 (I had wrongly simplified to
rely on Update()). No-op for yr23/24 (same value); fixes yr25. Recompile + re-run PbPb 23/24/25 (fixed seed
-> yr23/24 byte-identical; yr25 fixed). pp unaffected (no centrality). C4 note: clear root-caused code bug,
not a methodology issue -> direct fix + re-review (no separate /review-investigation needed).

Fix verification (re-run, all exit 0; counts identical to iter-0, fixed seed): yr25 OS avg_centrality
mean=15.22, frac at -1 = 0.000 (was 100% at -1 -> FIXED). yr23 mean=14.59, yr24 mean=14.96, both 0% at -1
(unchanged). All three years now have valid mixed-pair centrality labels. pp unaffected. Spawning iter-2 reviewer.

## Iteration 2
**Reviewer verdict**: PASS
**Issues found**: 0
**Details**: None. Fix verified: avg_centrality=(m1.ev_centrality+m2.ev_centrality)/2 (int) after Update(),
not overwritten; both muons same 5% interval -> valid label. yr25 mean=15.22 (0% at -1, was 100%); yr23
mean=14.59 unchanged; OS minv still smooth (no J/psi peak); per-centrality population central-weighted
([0,10)=5.81M ... [50,84]=0.36M yr25, monotonic). No other behavior changed.
**Numerical verification**: all MATCH.

**Status**: APPROVED at iteration 2
**Summary**: ScrambGen rewritten to object model (read MuonObj single-muon trees, mix across events within
5% ctr intervals, write muon_pair_tree_sign1/sign2 MuonPairObj). Combinatoric T_mix correct: SS≈OS, NO
resonance peaks (uncorrelated mixing), valid per-pair centrality (yr25 bug fixed). Scrambled muon_pairs
produced for pbpb 23/24/25 + pp24. Ready for the RDF mixed-event T_mix fill.
