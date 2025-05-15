#!/usr/bin/env bash

# make_JOs.sh  ─ generate 24 Pythia8 job‑option files from mc.template.py
template="mc.template.py"              # source JO (pp, 60‑120 slice)

# -------- isospin settings --------------------------------------------------
isospins=(pp  pn  np  nn)              # output labels
beam1=(PROTON PROTON NEUTRON NEUTRON)  # Beam1 values to write
beam2=(PROTON NEUTRON PROTON NEUTRON)  # Beam2 values to write

# -------- pTHat windows & per‑job stat --------------------------------------
labels=(5_7 7_10 10_25 25_60 60_120 120_3200)
pmin=(5. 7. 10. 25. 60. 120.)
pmax=(7. 10. 25. 60. 120. 3200.)
# nevt=(20 500 500 5000 10000 10000)
nevt=(10 2 2 20 100 200)

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

    out="JOs/mc.Py8EG_A14_${iso}_HardQCD_DiMu_pTH${tag}.py"

    if [[ "$j" -eq 0 ]]; then
      sed -e "s/genSeq.Pythia8.Beam1 = \"[A-Z]*\"/genSeq.Pythia8.Beam1 = \"${b1}\"/" \
          -e "s/genSeq.Pythia8.Beam2 = \"[A-Z]*\"/genSeq.Pythia8.Beam2 = \"${b2}\"/" \
          -e "s/pTHatMin = [0-9.]\+/pTHatMin = ${min}/" \
          -e "s/pTHatMax = [0-9.]\+/pTHatMax = ${max}/" \
          -e "s/evgenConfig.nEventsPerJob = [0-9]\+/evgenConfig.nEventsPerJob = ${nev}/" \
          -e "s/filtSeq.xAODMultiMuonFilter.NMuons = [0-9]\+/filtSeq.xAODMultiMuonFilter.NMuons = 1/" \
          "$template" > "$out"
    else
      sed -e "s/genSeq.Pythia8.Beam1 = \"[A-Z]*\"/genSeq.Pythia8.Beam1 = \"${b1}\"/" \
          -e "s/genSeq.Pythia8.Beam2 = \"[A-Z]*\"/genSeq.Pythia8.Beam2 = \"${b2}\"/" \
          -e "s/pTHatMin = [0-9.]\+/pTHatMin = ${min}/" \
          -e "s/pTHatMax = [0-9.]\+/pTHatMax = ${max}/" \
          -e "s/evgenConfig.nEventsPerJob = [0-9]\+/evgenConfig.nEventsPerJob = ${nev}/" \
          "$template" > "$out"
    fi

    echo "Created  $out   (Beam1=${b1}, Beam2=${b2}, pT^hat ${min}-${max} GeV, nEvt/job=${nev})"
  done
done
