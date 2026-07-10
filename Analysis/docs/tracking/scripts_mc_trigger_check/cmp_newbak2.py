import ROOT, os
ROOT.gROOT.SetBatch(True)
D="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/"
BEAMS=["pp","pn","np","nn"]
SL=["pTH8_14","pTH14_24","pTH24_40","pTH40_70","pTH70_125","pTH125_300"]
VEC=["muon_pt","muon_eta","muon_quality","muon_trk_pt","muon_d0","muon_truth_pt","truth_muon_pt"]
def snap(t):
    t.SetBranchStatus("*",0)
    for v in VEC+["eventNumber"]: t.SetBranchStatus(v,1)
    d={}
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        d[int(t.eventNumber)] = tuple(tuple(getattr(t,v)) for v in VEC)
    return d
hdr=f"{'sample':<26}{'N_new':>7}{'N_bak':>7}{'evt sets':>10}   " + "".join(f"{v.replace('muon_','mu_').replace('truth_','tr_'):>11}" for v in VEC)
print(hdr); print("-"*len(hdr))
allok=True
for b in BEAMS:
    for s in SL:
        base=f"Pythia_5p36TeV_{b}_hQCD_DiMu_{s}.FullSimPP24.NTUP"
        pnew, pbak = D+base+".root", D+base+".bak_20260709.root"
        if not os.path.exists(pbak): continue
        fn=ROOT.TFile.Open(pnew); tn=fn.Get("HeavyIonD3PD")
        fb=ROOT.TFile.Open(pbak); tb=fb.Get("HeavyIonD3PD")
        Nn,Nb=tn.GetEntries(),tb.GetEntries()
        dn, db = snap(tn), snap(tb)
        same_set = set(dn)==set(db)
        diffs=[0]*len(VEC)
        if same_set:
            for e,vn in dn.items():
                vb=db[e]
                for k in range(len(VEC)):
                    if vn[k]!=vb[k]: diffs[k]+=1
        ok = (Nn==Nb) and same_set and all(d==0 for d in diffs)
        allok &= ok
        print(f"{b+' '+s:<26}{Nn:>7}{Nb:>7}{'SAME' if same_set else 'DIFF':>10}   " + "".join(f"{d:>11}" for d in diffs) + ("" if ok else "   <-- DIFFERS"))
        fn.Close(); fb.Close()
print("\nRESULT: " + ("ALL IDENTICAL — re-skim changed nothing outside the trigger branches" if allok else "MISMATCHES FOUND"))
