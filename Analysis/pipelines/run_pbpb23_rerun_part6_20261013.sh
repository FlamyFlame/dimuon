#!/usr/bin/env bash
# One-off (2026-10-13): Pb+Pb 2023 rerun after adding data_pbpb23_part6.root (events of runs 463315/462759 that the
# May-2026 grid merge lost; docs/tracking/data_per_run_signal_yield_qa.md Round 5). PbPb23-ONLY: nothing combined.
#   1 event-selection cuts (event_level -> cuts -> cuts_alt)   2 Condor NTuples (run_ntuple_rerun_20261013.sh pbpb23)
#   3 trig-eff Tight (pipeline_pbpb_trig_eff.sh stages 5-8) + Medium (run_data_trigeff_medium_wp.sh)
#   4 crossx RDF fill (run_crossx_hist_filling_pbpb23.sh) + OS-SS twin   5 data QA
# Same commands as those scripts, except every ACLiC compile goes to a PRIVATE build dir (other sessions compile the
# same classes concurrently). START_AT=<n> resumes at step n. Last log line: "ALL DONE pbpb23" or "FATAL pbpb23: ...".
source /usatlas/u/yuhanguo/setup.sh >/dev/null 2>&1
set -uo pipefail
export X509_USER_PROXY=/tmp/x509up_u101379
A=/gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/Analysis
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data; Y=$D/pbpb_2023
L=/gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/IntNotes/data/luminosity
B=$A/pipelines/logs/build_pbpb23_rerun_$$; LG=$A/pipelines/logs/pbpb23_rerun_20261013; mkdir -p $B $LG
START_AT=${START_AT:-1}
log(){ echo "[$(date '+%F %T')] $*"; }
fail(){ log "FATAL pbpb23: $*"; exit 1; }
fresh(){ [[ -s $1 && $(stat -c %Y $1) -gt $2 && $(stat -c %s $1) -gt 100 ]] || fail "not fresh/empty: $1"; }
BD="gSystem->SetBuildDir(\"$B\", true);"

if (( START_AT <= 1 )); then
  T=$(date +%s); mkdir -p $Y/pre_rerun_20261013
  for f in event_sel_cuts_pbpb_2023.root event_sel_cuts_pbpb_2023_alt.root; do
    [[ -f $Y/pre_rerun_20261013/$f ]] || cp -p $Y/$f $Y/pre_rerun_20261013/
  done
  cd $A/plotting_codes/event_selection
  log "1a event_level(23)"
  root -l -b -q "plot_pbpb_event_sel_event_level.cxx(23)" > $LG/event_level.log 2>&1 < /dev/null || fail "event_level rc"
  grep -q "Cut TGraphs saved to" $LG/event_level.log && grep -q "Alt cut TGraphs saved to" $LG/event_level.log || fail "event_level: cut TGraphs not saved"
  log "1b cuts(23)"
  root -l -b -q "plot_pbpb_event_sel_cuts.cxx(23)" > $LG/cuts.log 2>&1 < /dev/null || fail "cuts rc"
  grep -q "All cuts saved to" $LG/cuts.log || fail "cuts: not saved"
  log "1c cuts_alt(23)"
  root -l -b > $LG/cuts_alt.log 2>&1 <<EOF || fail "cuts_alt rc"
.L plot_pbpb_event_sel_cuts.cxx
plot_pbpb_event_sel_cuts_alt(23);
gSystem->Exit(0);
EOF
  grep -q "All cuts saved to" $LG/cuts_alt.log || fail "cuts_alt: not saved"
  for f in event_sel_cuts_pbpb_2023.root event_sel_cuts_pbpb_2023_alt.root; do fresh $Y/$f $T; done
  grep -h "Excluding" $LG/event_level.log $LG/cuts.log $LG/cuts_alt.log | sort | uniq -c
  log "1 DONE evsel"
fi

if (( START_AT <= 2 )); then
  T=$(date +%s)
  log "2 NTuples (log $LG/ntuple.log)"
  CUTS_AFTER=$(date -d '2026-10-13 16:30' +%s) \
    bash $A/pipelines/run_ntuple_rerun_20261013.sh pbpb23 > $LG/ntuple.log 2>&1
  tail -1 $LG/ntuple.log | grep -q "ALL DONE pbpb23" || fail "NTuple: $(tail -1 $LG/ntuple.log)"
  grep "hadd .* OK" $LG/ntuple.log
  log "2 DONE NTuples"
fi

if (( START_AT <= 3 )); then
  cd $A/RDFBasedHistFilling
  T=$(date +%s); log "3a trig-eff RDF fill Tight"
  root -l -b > $LG/trigeff_rdf.log 2>&1 <<EOF || fail "trig-eff RDF rc"
$BD
ROOT::EnableImplicitMT(2);
.L RDFBasedHistFillingPbPb.cxx+
RDFBasedHistFillingPbPb pbpb(23);
pbpb.trigger_mode = 1;
pbpb.mindR_trig = 0.02;
pbpb.Run();
gSystem->Exit(0);
EOF
  RDF=$Y/histograms_real_pairs_pbpb_2023_single_mu4_coarse_q_eta_bin_qeta_fid.root; fresh $RDF $T
  log "3b turn-on fit Tight"; cd $A
  root -l -b > $LG/trigeff_fit.log 2>&1 <<EOF || fail "fit rc"
.L SingleMuEffcyPtTurnOnFitter.cxx
single_muon_trig_effcy_pT_fitting_PbPb(23);
gSystem->Exit(0);
EOF
  FIT=$Y/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root; fresh $FIT $T
  log "3c validate TF1 + TH2D fallback"
  root -l -b > $LG/trigeff_check.log 2>&1 <<EOF || fail "TF1/TH2D validation: $(tail -2 $LG/trigeff_check.log)"
TFile *ff = TFile::Open("$FIT"); int ntf1 = 0; TIter n1(ff->GetListOfKeys()); TKey *k;
while ((k = (TKey*)n1())) if (std::string(k->GetClassName()) == "TF1") ntf1++;
TFile *fh = TFile::Open("$RDF"); int n2 = 0; TIter n2i(fh->GetListOfKeys());
while ((k = (TKey*)n2i())) { std::string s = k->GetName(); if (std::string(k->GetClassName()) == "TH2D" && s.find("_divided") != std::string::npos && s.find("pt2nd_vs_q_eta2nd") != std::string::npos) n2++; }
std::cout << "TF1=" << ntf1 << " TH2D_fallback=" << n2 << std::endl;
gSystem->Exit((ntf1 >= 10 && n2 >= 4) ? 0 : 1);
EOF
  tail -1 $LG/trigeff_check.log
  log "3d trig-eff plots"; cd $A/plotting_codes/trig_effcy
  root -l -b -q "trig_effcy_plot_PbPb.cxx(23)" > $LG/trigeff_plot.log 2>&1 < /dev/null || fail "trig-eff plot rc"
  T=$(date +%s); log "3e Medium RDF fill"; cd $A/RDFBasedHistFilling
  root -l -b > $LG/trigeff_medium_rdf.log 2>&1 <<EOF || fail "Medium RDF rc"
$BD
ROOT::EnableImplicitMT(2);
.L RDFBasedHistFillingPbPb.cxx+
RDFBasedHistFillingPbPb pbpb(23);
pbpb.trigger_mode = 1;
pbpb.mindR_trig = 0.02;
pbpb.isTight = false;
pbpb.Run();
gSystem->Exit(0);
EOF
  fresh $Y/histograms_real_pairs_pbpb_2023_single_mu4_coarse_q_eta_bin_qeta_fid_medium_wp.root $T
  log "3f Medium fit"; cd $A
  root -l -b > $LG/trigeff_medium_fit.log 2>&1 <<EOF || fail "Medium fit rc"
.L SingleMuEffcyPtTurnOnFitter.cxx
single_muon_trig_effcy_pT_fitting_PbPb(23, "_medium_wp");
gSystem->Exit(0);
EOF
  fresh $Y/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit_medium_wp.root $T
  log "3 DONE trig-eff"
fi

if (( START_AT <= 4 )); then
  cd $A/RDFBasedHistFilling
  T=$(date +%s); log "4a crossx RDF fill"
  root -l -b > $LG/crossx_rdf.log 2>&1 <<EOF || fail "crossx RDF rc"
$BD
.L RDFBasedHistFillingPbPb.cxx+
RDFBasedHistFillingPbPb pbpb(23);
pbpb.trigger_mode = 1;
pbpb.mu4_nominal_pbpb_NO_trig_calc = true;
pbpb.mindR_trig   = 0.02;
pbpb.Run();
gSystem->Exit(0);
EOF
  for f in $Y/histograms_real_pairs_pbpb_2023_single_mu4_no_trg_plots_nominal.root; do
    fresh $f $T
    root -l -b >/dev/null 2>&1 <<EOF || fail "$f: h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr0_5 missing/empty (RDF threw)"
TFile* f = TFile::Open("$f"); if (!f || f->IsZombie()) gSystem->Exit(2);
TH1* h = (TH1*)f->Get("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr0_5");
gSystem->Exit((h && h->GetEntries() > 0) ? 0 : 3);
EOF
  done
  T=$(date +%s); log "4b OS-SS crossx RDF fill"
  BUILD_DIR=$B/ss ./run_crossx_hist_filling_pbpb_ss_subtr.sh 23 > $LG/crossx_ss_subtr.log 2>&1 < /dev/null || fail "ss_subtr rc"
  f=$Y/histograms_real_pairs_pbpb_2023_single_mu4_no_trg_plots_nominal_ss_subtr_crossx.root; fresh $f $T
  root -l -b >/dev/null 2>&1 <<EOF || fail "$f: OS signal histogram missing/empty"
TFile* f = TFile::Open("$f"); if (!f || f->IsZombie()) gSystem->Exit(2);
TH1* h = (TH1*)f->Get("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr0_5");
gSystem->Exit((h && h->GetEntries() > 0) ? 0 : 3);
EOF
  log "4 DONE crossx"
fi

if (( START_AT <= 5 )); then
  log "5 data QA (log $LG/qa.log)"
  bash $A/pipelines/qa_pbpb_year.sh 23 6 false $B/qa > $LG/qa.log 2>&1 < /dev/null
  tail -1 $LG/qa.log | grep -q "ALL DONE pbpb23 QA" || fail "QA: $(tail -1 $LG/qa.log)"
  log "5 DONE QA"
fi
log "ALL DONE pbpb23"
