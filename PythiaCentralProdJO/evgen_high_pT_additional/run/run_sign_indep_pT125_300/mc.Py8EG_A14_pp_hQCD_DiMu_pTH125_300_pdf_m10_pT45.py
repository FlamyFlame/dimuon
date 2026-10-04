# ---------------------------------------------------------------
# Hard-process definition & phase-space slice
# ---------------------------------------------------------------
pTHatMin = 125.
pTHatMax = 300.

# ---------------------------------------------------------------
# Muon-pair filter (MeV)
# ---------------------------------------------------------------
PairMassMax  = 10000.     # m(mumu)  < 10 GeV
PairPtMin    = 45000.     # pT(mumu) > 45 GeV
PairSameSign = False

include('Py8Control_hQCD_dimu_pdf_dilepton_filters.py')

# filter efficiency ~2.7e-4 -> ~0.04 CPU-h per accepted event; 1000 events/job ~ 40 CPU-h
evgenConfig.nEventsPerJob = 1000

genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"