# ---------------------------------------------------------------
# Hard-process definition & phase-space slice
# ---------------------------------------------------------------
pTHatMin = 24.
pTHatMax = 40.

# ---------------------------------------------------------------
# Muon-pair filter (MeV)
# ---------------------------------------------------------------
PairMassMax  = 10000.     # m(mumu)  < 10 GeV
PairPtMin    = -1.        # NO pair-pT cut in this slice (user decision 2026-09-22:
                          # with pT(mumu) > 20 GeV the filter efficiency is ~2e-7)
PairSameSign = True

include('Py8Control_hQCD_dimu_pdf_dilepton_filters.py')

# re-tune after the test run, once the filter efficiency is known
evgenConfig.nEventsPerJob = 500

genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"