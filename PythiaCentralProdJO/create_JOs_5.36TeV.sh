#!/usr/bin/env bash
set -euo pipefail

template="mc.template.5.36TeV.py"
outdir="JOs_5.36TeV";  mkdir -p "$outdir"

# ---------------- isospin settings ----------------
isospins=(pp pn np nn)                      # output labels
beam1=(PROTON PROTON NEUTRON NEUTRON)       # Beam1 values
beam2=(PROTON NEUTRON PROTON NEUTRON)       # Beam2 values

# ---------------- pTHat windows -------------------
labels=(5_8 8_14 14_24 24_40 40_70 70_125 125_300)
pmin=(5. 8. 14. 24. 40. 70. 125.)
pmax=(8. 14. 24. 40. 70. 125. 300.)
nevt=(1 10 50 200 500 2000 2000)

# ---------------- main loop -----------------------
for i in "${!isospins[@]}"; do
  iso=${isospins[$i]}
  b1=${beam1[$i]}
  b2=${beam2[$i]}

  for j in "${!labels[@]}"; do
    tag=${labels[$j]}  min=${pmin[$j]}  max=${pmax[$j]}  nev=${nevt[$j]}
    out="${outdir}/mc.Py8EG_A14_5.36_${iso}_hQCD_DiMu_pTH${tag}.py"

    sed  -e "s/\(genSeq\.Pythia8\.Beam1[[:space:]]*=[[:space:]]*\"\)[^\"]*\(\"\)/\1${b1}\2/" \
         -e "s/\(genSeq\.Pythia8\.Beam2[[:space:]]*=[[:space:]]*\"\)[^\"]*\(\"\)/\1${b2}\2/" \
         -e "s/\(pTHatMin[[:space:]]*=[[:space:]]*\)[0-9.+-eE]\+/\1${min}/g" \
         -e "s/\(pTHatMax[[:space:]]*=[[:space:]]*\)[0-9.+-eE]\+/\1${max}/g" \
         -e "s/\(evgenConfig\.nEventsPerJob[[:space:]]*=[[:space:]]*\)[0-9]\+/\1${nev}/" \
         "$template" > "$out"

    printf 'Created %-60s  (Beam1=%s, Beam2=%s, pTHat %s-%s GeV, nEvt/job=%s)\n' \
           "$out" "$b1" "$b2" "$min" "$max" "$nev"
  done
done
