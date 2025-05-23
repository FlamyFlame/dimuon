#!/usr/bin/env bash

# make_JOs.sh  ─ generate 24 Pythia8 job‑option files from mc.template.5.02TeV.py
template="mc.template.5.02TeV.py"              # source JO (pp, 60‑120 slice)

# -------- isospin settings --------------------------------------------------
isospins=(pp  pn  np  nn)              # output labels
beam1=(PROTON PROTON NEUTRON NEUTRON)  # Beam1 values to write
beam2=(PROTON NEUTRON PROTON NEUTRON)  # Beam2 values to write

# -------- pTHat windows & per‑job stat --------------------------------------
labels=(5_8 8_14 14_24 24_40 40_70 70_125 125_300)
pmin=(5. 8. 14. 24. 40. 70. 125.)
pmax=(8. 14. 24. 40. 70. 125. 300.)
nevt=(1 10 50 200 500 2000 2000)

# -------- main loop ---------------------------------------------------------
for i in "${!isospins[@]}"; do
  iso=${isospins[$i]}
  b1=${beam1[$i]}
  b2=${beam2[$i]}

  for j in "${!labels[@]}"; do
    tag=${labels[$j]}
    min=${pmin[$j]}
    max=${pmax[$j]}
    nev=${nevt[$j]}

    out="JOs_5.02TeV/mc.Py8EG_A14_5.02_${iso}_hQCD_DiMu_pTH${tag}.py"

    sed -e "s/genSeq.Pythia8.Beam1 = \"[A-Z]*\"/genSeq.Pythia8.Beam1 = \"${b1}\"/" \
        -e "s/genSeq.Pythia8.Beam2 = \"[A-Z]*\"/genSeq.Pythia8.Beam2 = \"${b2}\"/" \
        -e "s/pTHatMin = [0-9.]\+/pTHatMin = ${min}/" \
        -e "s/pTHatMax = [0-9.]\+/pTHatMax = ${max}/" \
        -e "s/evgenConfig.nEventsPerJob = [0-9]\+/evgenConfig.nEventsPerJob = ${nev}/" \
        "$template" > "$out"

    echo "Created  $out   (Beam1=${b1}, Beam2=${b2}, pT^hat ${min}-${max} GeV, nEvt/job=${nev})"
  done
done
