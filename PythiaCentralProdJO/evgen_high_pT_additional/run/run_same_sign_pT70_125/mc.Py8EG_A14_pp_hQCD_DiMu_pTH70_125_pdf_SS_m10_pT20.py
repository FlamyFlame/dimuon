# ---------------------------------------------------------------
# Hard-process definition & phase-space slice
# ---------------------------------------------------------------
pTHatMin = 70.
pTHatMax = 125.

# ---------------------------------------------------------------
# Muon-pair filter (MeV)
# ---------------------------------------------------------------
PairMassMax  = 10000.     # m(mumu)  < 10 GeV
PairPtMin    = 20000.     # pT(mumu) > 20 GeV
PairSameSign = True

include('Py8Control_hQCD_dimu_pdf_dilepton_filters.py')

# filter efficiency ~3.2e-5 -> ~0.24 CPU-h per accepted event; 200 events/job ~ 48 CPU-h
evgenConfig.nEventsPerJob = 200

genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"