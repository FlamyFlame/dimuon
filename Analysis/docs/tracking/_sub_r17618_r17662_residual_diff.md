# SUB: r17618 vs r17662 residual reco-eff difference (scratch, append-only)

## Objective
Explain the small residual reco-eff disagreement between r17618 (full HIJING truth,
barcode collision Pythia vs HIJING) and r17662 (StandardSignalOnlyTruth) on the SAME
10000 events / SAME reco. Confirm or rule out the "Pythia truth muon barcode collides
with a HIJING truth MUON barcode -> spurious reco match" hypothesis.

## Constraints
- READ-ONLY on raw NTUPs. No git. No edits to tracked analysis files.
- Mirror `PythiaFullSimExtras.c::ProcessEventFullsim` matching EXACTLY:
  - boundary = index of FIRST truth_barcode > 200000
  - GetNPythiaTruthMuons = count(status==1 && |truth_id|==13) in [0, boundary)
  - take first that many entries of truth_muon_* vectors
  - real list = reco muons with muon_truth_prob > 0.5
  - for each Pythia truth muon in order: first UNCLAIMED reco muon in list with
    muon_truth_barcode == truth barcode -> claim. No dR fallback.
- Key events by eventNumber (entry order DIFFERS between files).

## Plan
- Step 0: inspect branches/types, build eventNumber->entry maps.
- Step 1: verify reco muon collections identical; verify Pythia truth muon content identical;
          explain the ~4-muon fiducial count difference.
- Step 2: enumerate disagreements (reco_match, matched reco index, pass_medium/tight).
- Step 3: test HIJING-muon-barcode-collision hypothesis on each disagreement.
- Step 4: if not explained, check muon_truth_prob/muon_truth_barcode differing between r-tags.
- Step 5: verdict + bias direction + fix proposal.

## Log
(append below)
