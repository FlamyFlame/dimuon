# ---------------------------------------------------------------
# Hard-process definition & phase-space slice
# ---------------------------------------------------------------
pTHatMin = 125.
pTHatMax = 300.

# ---------------------------------------------------------------
# Muon-pair filter (MeV)
# ---------------------------------------------------------------
PairMassMax  = 10000.     # m(mumu)  < 10 GeV
PairPtMin    = 20000.     # pT(mumu) > 20 GeV
PairSameSign = True

include('Py8Control_hQCD_dimu_pdf_dilepton_filters.py')

# filter efficiency ~1.2e-4 -> ~0.092 CPU-h per accepted event; 500 events/job ~ 46 CPU-h
evgenConfig.nEventsPerJob = 500

genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"