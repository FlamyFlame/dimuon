#!/usr/bin/env bash
# Build + run the mixed-event combinatoric pair generator (ScrambGen) for the low-mass
# dimuon template fit. PbPb per-year (23/24/25/26) and pp24, reading the single-muon trees and
# writing muon_pairs_*_scrambled.root (muon_pair_tree_sign1/sign2, MuonPairObj). In-memory,
# fast, local. PbPb and pp run in SEPARATE ROOT sessions (different object models).
set -Eeuo pipefail
set +e; set +u
source ~/setup.sh
setup_status=$?
set -e; set -u
if [[ ${setup_status} -ne 0 ]]; then echo "Environment setup failed (source ~/setup.sh)."; exit ${setup_status}; fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

echo "[ScrambGen] PbPb (years with data) ..."
# Only run years whose single-muon trees exist: ScrambGen now THROWS on an empty input
# (rather than writing an empty scrambled file that downstream existence checks accept),
# so a year still being produced would otherwise abort the whole run.
DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"
SG_YEARS=()
for yr in 23 24 25 26; do
  if compgen -G "${DATA_BASE}/pbpb_20${yr}/single_muon_trees_pbpb_20${yr}_part*_single_mu4_mindR_0_02.root" > /dev/null; then
    SG_YEARS+=("${yr}")
  else
    echo "[SKIP] Pb+Pb 20${yr}: no single-muon trees on disk -- ScrambGen skipped." >&2
  fi
done
if [[ ${#SG_YEARS[@]} -eq 0 ]]; then
  echo "[FATAL] no Pb+Pb year has single-muon trees on disk." >&2
  exit 1
fi
echo "[INFO] ScrambGen Pb+Pb years: ${SG_YEARS[*]}"
{
  echo '.L ScrambGen.c+'
  echo 'ScrambGen g;'
  for yr in "${SG_YEARS[@]}"; do echo "g.Run(${yr});"; done
  echo 'gSystem->Exit(0);'
} | root -l -b

echo "[ScrambGen] pp24 ..."
root -l -b <<EOF
.L ScrambGenPP.c+
ScrambGenPP g;
g.Run();
gSystem->Exit(0);
EOF

echo "[ScrambGen] done. Output files:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20{23,24,25,26}/muon_pairs_pbpb_*_scrambled.root \
       /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024_2mu4_scrambled.root 2>/dev/null
