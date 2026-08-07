#!/usr/bin/env bash
# Round 8: regenerate the MC trigger-efficiency chain on the NOMINAL 8 log pair-pT bins and then
# again on the OPT-IN 4-bin comparison variant, with EVERYTHING else identical (same samples, same
# q*eta gap cut, same coarse q*eta binning, same plateau window, same fit methods). The only
# difference between the two passes is the number of pair-pT cells, which is the whole point --
# it isolates whether the finer binning buys resolution or noise.
#
# The 4-bin pass is selected ONLY by MCTRIGEFF_PAIRPT_4BIN, which drives the binning AND the
# output naming from one place (Utilities/MCTrigEffPairPtBinning.h), so the two passes cannot
# collide: hists get a _pt4bin token, plots a pt4bin/ subdirectory.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
mkdir -p pipelines/logs_round8

for mode in 8bin 4bin; do
  if [[ $mode == 4bin ]]; then export MCTRIGEFF_PAIRPT_4BIN=1; else unset MCTRIGEFF_PAIRPT_4BIN; fi
  echo "================ MC chain, ${mode} ================"
  SAMPLES="pp_full overlay noovl" WPS="tight medium" \
      bash pipelines/run_mc_trigeff_round7.sh   > "pipelines/logs_round8/mc_chain_${mode}.log" 2>&1
  echo "  chain rc=$?  -> pipelines/logs_round8/mc_chain_${mode}.log"
  SAMPLES="pp_full overlay noovl" WPS="tight medium" \
      bash pipelines/run_dr_correction_fits.sh  > "pipelines/logs_round8/dr_fits_${mode}.log" 2>&1
  echo "  fits  rc=$?  -> pipelines/logs_round8/dr_fits_${mode}.log"
done
echo "BOTH_BINNINGS_DONE"
