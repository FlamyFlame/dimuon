#!/usr/bin/env bash
# Build + run the mixed-event combinatoric pair generator (ScrambGen) for the low-mass
# dimuon template fit. PbPb per-year (23/24/25) and pp24, reading the single-muon trees and
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

echo "[ScrambGen] PbPb 23/24/25 ..."
root -l -b <<EOF
.L ScrambGen.c+
ScrambGen g;
g.Run(23);
g.Run(24);
g.Run(25);
gSystem->Exit(0);
EOF

echo "[ScrambGen] pp24 ..."
root -l -b <<EOF
.L ScrambGenPP.c+
ScrambGenPP g;
g.Run();
gSystem->Exit(0);
EOF

echo "[ScrambGen] done. Output files:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20{23,24,25}/muon_pairs_pbpb_*_scrambled.root \
       /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024_2mu4_scrambled.root 2>/dev/null
