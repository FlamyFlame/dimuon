# ---------------------------------------------------------------
# Hard-process definition & phase-space slice
# ---------------------------------------------------------------
pTHatMin = 24.
pTHatMax = 40.

# ---------------------------------------------------------------
# Muon-pair filter (MeV)
# ---------------------------------------------------------------
PairMassMax  = 10000.     # m(mumu)  < 10 GeV
PairPtMin    = 20000.     # pT(mumu) > 20 GeV
PairSameSign = True

include('Py8Control_hQCD_dimu_pdf_dilepton_filters.py')

# re-tune after the test run, once the filter efficiency is known
evgenConfig.nEventsPerJob = 500

genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"