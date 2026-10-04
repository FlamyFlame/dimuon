# ---------------------------------------------------------------
# Hard-process definition & phase-space slice
# ---------------------------------------------------------------
pTHatMin = 70.
pTHatMax = 125.

# ---------------------------------------------------------------
# Muon-pair filter (MeV)
# ---------------------------------------------------------------
PairMassMax  = 10000.     # m(mumu)  < 10 GeV
PairPtMin    = 45000.     # pT(mumu) > 45 GeV
PairSameSign = False

include('Py8Control_hQCD_dimu_pdf_dilepton_filters.py')

# filter efficiency ~4.3e-5 -> ~0.15 CPU-h per accepted event; 200 events/job ~ 31 CPU-h
evgenConfig.nEventsPerJob = 200

genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"