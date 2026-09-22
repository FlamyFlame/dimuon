## ------------------------------------------------------------------
## Pythia8 Hard-QCD -> all, pTHatMin < pTHat < pTHatMax,
## require >=2 muons with pT > 3.7 GeV forming a pair with
## m(mumu) < PairMassMax and pT(mumu) > PairPtMin
## (optionally same-sign only).
## ------------------------------------------------------------------
## Set in the top-level mc.*.py, before the include:
##   pTHatMin, pTHatMax          [GeV]
##   PairMassMax, PairPtMin      [MeV]
##   PairSameSign                [bool]
## ------------------------------------------------------------------

evgenConfig.description = 'Pythia HardQCD:All with dimuon pair mass/pT filter; Py8 A14 NNPDF23LO'
evgenConfig.keywords += ['QCD','SM']
evgenConfig.contact = ['yuhan.guo@cern.ch']
evgenConfig.process = 'HardQCD -> all'

# ---------------------------------------------------------------
# Base fragments (A14 tune + EvtGen hooks kept minimal)
# ---------------------------------------------------------------
include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")

genSeq.Pythia8.Commands += [
    'PDF:pSet = LHAPDF6:NNPDF23_lo_as_0119_qed',
]

genSeq.Pythia8.Commands += [
    'HardQCD:all = on',
    f'PhaseSpace:pTHatMin = {pTHatMin:.1f}',
    f'PhaseSpace:pTHatMax = {pTHatMax:.1f}',
]

# ---------------------------------------------------------------
# Dimuon truth filter (>=2 muons with pT > 3.7GeV)
# ---------------------------------------------------------------

# muon filter
include('GeneratorFilters/xAODMultiMuonFilter_Common.py')
filtSeq.xAODMultiMuonFilter.Ptcut = 3700.
filtSeq.xAODMultiMuonFilter.Etacut = 2.5
filtSeq.xAODMultiMuonFilter.NMuons = 2

# ---------------------------------------------------------------
# muon-pair filter -- pair invariant mass and pair pT
#
# MinMass/MaxMass and MinDilepPt act on the mu-mu 4-vector sum
# (xAODDiLeptonMassFilter.cxx: vec.m() and vec.perp()), i.e. the pair
# quantities; MinPt/MaxEta are the per-muon acceptance cuts and are
# kept identical to the xAODMultiMuonFilter ones above.
#
# The container is switched from the default "TruthLightLeptons"
# (electrons + muons) to the slimmed "TruthMuons" collection that the
# xAODMultiMuonFilter fragment above already creates, so that only muon
# pairs can satisfy the filter -- otherwise an e+e- pair could accept
# the event, since with AllowSameCharge=True the pair condition is
# |pdgId1| == |pdgId2|.
#
# RequireSameCharge is the new property (see athena_patch/); it needs a
# release that contains the patched filter.
# ---------------------------------------------------------------
include('GeneratorFilters/xAODDiLeptonMassFilter_Common.py')
filtSeq.xAODDiLeptonMassFilter.TruthParticleContainerKey = "TruthMuons"
filtSeq.xAODDiLeptonMassFilter.MinPt             = 3700.
filtSeq.xAODDiLeptonMassFilter.MaxEta            = 2.5
filtSeq.xAODDiLeptonMassFilter.MinMass           = 0.
filtSeq.xAODDiLeptonMassFilter.MaxMass           = PairMassMax
filtSeq.xAODDiLeptonMassFilter.MinDilepPt        = PairPtMin
filtSeq.xAODDiLeptonMassFilter.AllowElecMu       = False
filtSeq.xAODDiLeptonMassFilter.AllowSameCharge   = True

if PairSameSign:
    # only set when needed, so the sign-independent samples still run in a
    # release that does not yet have the patched filter
    filtSeq.xAODDiLeptonMassFilter.RequireSameCharge = True
