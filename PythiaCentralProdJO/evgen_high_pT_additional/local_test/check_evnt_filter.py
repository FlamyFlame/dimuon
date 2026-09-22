# Event-by-event check of the generator filter on an EVNT file (truth, HepMC GEN_EVENT).
# Muon/pair definitions = the FILTER's own (xAODMultiMuonFilter + xAODDiLeptonMassFilter in
# run/Py8Control_hQCD_dimu_pdf_dilepton_filters.py): stable (status 1) muons, pT > 3.7 GeV,
# |eta| < 2.5; a passing pair has 0 < m < MAXMASS and pT(mumu) > MINPAIRPT (same sign if SAMESIGN).
# Usage: athena.py -c "INPUT='x.EVNT.pool.root'; MINPAIRPT=20000.; SAMESIGN=True" check_evnt_filter.py
import math
import AthenaPoolCnvSvc.ReadAthenaPool
from AthenaCommon.AppMgr import ServiceMgr as svcMgr, theApp
from AthenaCommon.AlgSequence import AlgSequence
from AthenaPython import PyAthena
svcMgr.EventSelector.InputCollections = [INPUT]
theApp.EvtMax = -1
MAXMASS = globals().get('MAXMASS', 10000.)

class CheckFilter(PyAthena.Alg):
    def initialize(self):
        self.n = 0; self.nfail = 0
        return PyAthena.StatusCode.Success
    def execute(self):
        coll = self.evtStore['GEN_EVENT']
        evt = coll[0]
        mus = []
        for p in evt.particles():
            if abs(p.pdg_id()) != 13 or p.status() != 1: continue
            m = p.momentum()
            pt = math.hypot(m.px(), m.py())
            if pt <= 3700.: continue
            eta = m.eta()
            if abs(eta) >= 2.5: continue
            mus.append((p.pdg_id(), m.px(), m.py(), m.pz(), m.e(), pt, eta))
        best = []
        for i in range(len(mus)):
            for j in range(i+1, len(mus)):
                a, b = mus[i], mus[j]
                ss = a[0]*b[0] > 0
                px, py, pz, e = a[1]+b[1], a[2]+b[2], a[3]+b[3], a[4]+b[4]
                m2 = e*e-px*px-py*py-pz*pz
                mass = math.sqrt(max(m2, 0.)); ppt = math.hypot(px, py)
                ok = (0. < mass < MAXMASS) and ppt > MINPAIRPT and (ss or not SAMESIGN)
                best.append((ok, ss, mass/1e3, ppt/1e3))
        passed = any(b[0] for b in best)
        self.n += 1
        if not passed or len(mus) < 2: self.nfail += 1
        okpairs = [b for b in best if b[0]]
        print("CHECKEVT evt=%d nmu=%d npairs=%d npass=%d PASS=%s passing(sign,m,pT)=%s" % (
            self.n, len(mus), len(best), len(okpairs), passed and len(mus) >= 2,
            ["%s/%.2f/%.1f" % ("SS" if b[1] else "OS", b[2], b[3]) for b in okpairs]))
        return PyAthena.StatusCode.Success
    def finalize(self):
        print("CHECKSUMMARY events=%d failing=%d" % (self.n, self.nfail))
        return PyAthena.StatusCode.Success

topSeq = AlgSequence()
topSeq += CheckFilter("CheckFilter")
