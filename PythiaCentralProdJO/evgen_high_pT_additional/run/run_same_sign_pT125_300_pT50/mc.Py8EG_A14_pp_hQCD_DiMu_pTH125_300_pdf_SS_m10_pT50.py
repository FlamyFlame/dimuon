# ---------------------------------------------------------------
# Hard-process definition & phase-space slice
# ---------------------------------------------------------------
pTHatMin = 125.
pTHatMax = 300.

# ---------------------------------------------------------------
# Muon-pair filter (MeV)
# ---------------------------------------------------------------
PairMassMax  = 10000.     # m(mumu)  < 10 GeV
PairPtMin    = 50000.     # pT(mumu) > 50 GeV
PairSameSign = True

include('Py8Control_hQCD_dimu_pdf_dilepton_filters.py')

# filter efficiency ~1e-5 -> ~1 CPU-h per accepted event; 50 events/job ~ 50 CPU-h
evgenConfig.nEventsPerJob = 50

genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"
