##  mc.Py8EG_A14NNPDF23LO_HardQCD_DiMuon_5p02TeV_pTHat15_45.py
## ------------------------------------------------------------------
##  Pythia8 Hard‑QCD → all, √s = 5.02TeV, pTHatMin<p̂T<pTHatMax,
##  require ≥2 muons with pT >3.7GeV (no charge requirement).
## ------------------------------------------------------------------

# ---------------------------------------------------------------
#   Metadata (update run number, contact, XS & ε once validated)
# ---------------------------------------------------------------
evgenConfig.description  = 'Pythia HardQCD:All with dimuon filter; Py8 A14 NNPDF23LO'
evgenConfig.keywords     += ['QCD','SM']
evgenConfig.contact       = ['']
evgenConfig.generators    = ['Pythia8']
# evgenConfig.nEventsPerJob = 2000  # typical value


# ---------------------------------------------------------------
#   Base fragments (A14 tune + EvtGen hooks kept minimal)
# ---------------------------------------------------------------
include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")

# ---------------------------------------------------------------
#   Retrieve & set parameters from run arguments
# ---------------------------------------------------------------

beam1 = "PROTON"
if hasattr(runArgs, 'beam1'):
    beam1 = runArgs.beam1

beam2 = "PROTON"
if hasattr(runArgs, 'beam2'):
    beam2 = runArgs.beam2

ecmEnergy = 5020.
if hasattr(runArgs, 'ecmEnergy'):
    ecmEnergy = runArgs.ecmEnergy

pTHatMin = 10.
if hasattr(runArgs, 'pTHatMin'):
    pTHatMin = runArgs.pTHatMin

pTHatMax = 25.
if hasattr(runArgs, 'pTHatMax'):
    pTHatMax = runArgs.pTHatMax

# ---------------------------------------------------------------
#   Beam energy and PDF for 5.02TeV Pb–Pb (nCTEQ15 optional)
# ---------------------------------------------------------------
genSeq.Pythia8.Beam1 = beam1
genSeq.Pythia8.Beam2 = beam2

genSeq.Pythia8.Commands += [
  f'Beams:eCM = {ecmEnergy:.1f}',           # √s NN
]

genSeq.Pythia8.Commands += [
    # Nominal PDF (already in your config)
    'PDF:pSet = LHAPDF6:nNNPDF30_nlo_as_0118_A208_Z82/0001',
    # Enable LHAPDF and reweighting
    'Variations:doVariations = on',  # Enable weight variations
]

# ---------------------------------------------------------------
#   Hard‑process definition & phase‑space slice
# ---------------------------------------------------------------
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
    testSeq.TestHepMC.CmEnergy = ecmEnergy * 1000.  # Beam COM energy in MeV


# ---------------------------------------------------------------
#   Dimuon truth filter (≥2 muons with pT > 3.7GeV)
# ---------------------------------------------------------------

# muon filter
include("GeneratorFilters/MultiMuonFilter.py")
MultiMuonFilter = filtSeq.MultiMuonFilter
MultiMuonFilter.Ptcut = 3700.
MultiMuonFilter.Etacut = 2.5
MultiMuonFilter.NMuons = 2
