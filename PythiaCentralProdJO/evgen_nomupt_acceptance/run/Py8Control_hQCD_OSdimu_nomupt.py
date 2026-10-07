## ------------------------------------------------------------------
## Signal-acceptance sample (tracking doc signal_acceptance_mass_window.md PD9, 2026-10-07):
## Pythia8 HardQCD:all, pTHatMin < pTHat < pTHatMax, pp beams, A14 + NNPDF23_lo_as_0119_qed + EvtGen
## (= the pp24 central PDF evgen 803015-803020 configuration), with a filter that keeps EVERY
## target pair: >= 1 opposite-sign muon pair, both |eta_mu| <= 2.5, NO muon-pT cut, pT(mumu) > 8 GeV,
## no mass cut. (Target T_target = single-b OS pair, pT^pair > 9, |eta^pair| < 2.2, |eta_mu| < 2.5.)
## Set pTHatMin, pTHatMax [GeV] in the top-level mc.*.py before the include.
## ------------------------------------------------------------------
evgenConfig.description = 'Pythia HardQCD:All, OS dimuon pT>8 GeV filter, no muon pT cut; Py8 A14 NNPDF23LO EvtGen'
evgenConfig.keywords += ['QCD','SM']
evgenConfig.contact = ['yuhan.guo@cern.ch']
evgenConfig.process = 'HardQCD -> all'

include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
genSeq.Pythia8.Commands += ['PDF:pSet = LHAPDF6:NNPDF23_lo_as_0119_qed']
genSeq.Pythia8.Beam1 = "PROTON"
genSeq.Pythia8.Beam2 = "PROTON"
genSeq.Pythia8.Commands += [
    'HardQCD:all = on',
    f'PhaseSpace:pTHatMin = {pTHatMin:.1f}',
    f'PhaseSpace:pTHatMax = {pTHatMax:.1f}',
]

# all stable muons (no pT cut) into the slimmed TruthMuons container
include("GeneratorFilters/CreatexAODSlimContainers.py")
createxAODSlimmedContainer("TruthMuons", prefiltSeq)
prefiltSeq.xAODCnv.AODContainerName = 'GEN_EVENT'

include('GeneratorFilters/xAODDiLeptonMassFilter_Common.py')
filtSeq.xAODDiLeptonMassFilter.TruthParticleContainerKey = "TruthMuons"   # muon pairs only
filtSeq.xAODDiLeptonMassFilter.MinPt           = 0.        # NO muon-pT cut
filtSeq.xAODDiLeptonMassFilter.MaxEta          = 2.5
filtSeq.xAODDiLeptonMassFilter.MinMass         = -1.       # no mass cut (filter tests MinMass < m)
filtSeq.xAODDiLeptonMassFilter.MaxMass         = 1.e9      # [MeV]
filtSeq.xAODDiLeptonMassFilter.MinDilepPt      = 8000.     # pT(mumu) > 8 GeV
filtSeq.xAODDiLeptonMassFilter.AllowElecMu     = False
filtSeq.xAODDiLeptonMassFilter.AllowSameCharge = False     # opposite sign only
