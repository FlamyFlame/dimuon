# ---------------------------------------------------------------
# Hard-process definition & phase-space slice
# ---------------------------------------------------------------
pTHatMin = 40.
pTHatMax = 70.

# ---------------------------------------------------------------
# Muon-pair filter (MeV)
# ---------------------------------------------------------------
PairMassMax  = 10000.     # m(mumu)  < 10 GeV
PairPtMin    = 20000.     # pT(mumu) > 20 GeV
PairSameSign = True

include('Py8Control_hQCD_dimu_pdf_dilepton_filters.py')

# filter efficiency ~4.0e-6 -> ~2.3 CPU-h per accepted event; 20 events/job ~ 47 CPU-h
evgenConfig.nEventsPerJob = 20

genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"