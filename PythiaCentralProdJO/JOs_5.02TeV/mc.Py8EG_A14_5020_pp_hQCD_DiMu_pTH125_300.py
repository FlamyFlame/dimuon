## ------------------------------------------------------------------
##  Pythia8 Hard‑QCD → all, √s = 5.02TeV, pTHatMin<p̂T<pTHatMax,
##  require ≥2 muons with pT >3.7GeV (no charge requirement).
## ------------------------------------------------------------------

evgenConfig.description  = 'Pythia HardQCD:All with dimuon filter; Py8 A14 NNPDF23LO'
evgenConfig.keywords     += ['QCD','SM']
evgenConfig.generators    = ['Pythia8']
evgenConfig.contact       = ['Yuhan Guo']
evgenConfig.process       = 'HardQCD -> all'

evgenConfig.nEventsPerJob = 2000

# ---------------------------------------------------------------
#   Base fragments (A14 tune + EvtGen hooks kept minimal)
# ---------------------------------------------------------------
include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")

genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"
runArgs.ecmEnergy = 5020. # direct setting of genSeq.Pythia8.CollisionEnergy gets overwritten by runArgs.ecmEnergy in EvgenJobTransforms/Generate_ecmenergies.py

genSeq.Pythia8.Commands += [
    'PDF:pSet = LHAPDF6:nNNPDF30_nlo_as_0118_A208_Z82/0001',
]

# ---------------------------------------------------------------
#   Hard‑process definition & phase‑space slice
# ---------------------------------------------------------------
pTHatMin = 125.
pTHatMax = 300.

genSeq.Pythia8.Commands += [
  'HardQCD:all = on',
  f'PhaseSpace:pTHatMin = {pTHatMin:.1f}',
  f'PhaseSpace:pTHatMax = {pTHatMax:.1f}',
  # miscellaneous hygiene
  'Init:showChangedSettings = on',
  'Next:numberCount = 0',
]

# ---------------------------------------------------------------
# Configure TestHepMC to expect the correct E_{CM} energy
# ---------------------------------------------------------------

if hasattr(testSeq, 'TestHepMC'):
    testSeq.remove(testSeq.TestHepMC)


# ---------------------------------------------------------------
#   Dimuon truth filter (≥2 muons with pT > 3.7GeV)
# ---------------------------------------------------------------

# muon filter
include('GeneratorFilters/xAODMultiMuonFilter_Common.py')
filtSeq.xAODMultiMuonFilter.Ptcut = 3700.
filtSeq.xAODMultiMuonFilter.Etacut = 2.5
filtSeq.xAODMultiMuonFilter.NMuons = 2


