import ROOT, json, os
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare(open("/tmp/claude-101379/sel.h").read())
D="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/"
ami=json.load(open("/tmp/claude-101379/ami.json"))
SL=["pTH8_14","pTH14_24","pTH24_40","pTH40_70","pTH70_125","pTH125_300"]
MU4="b_HLT_mu4_L1MU3V"; DI="b_HLT_2mu4_L12MU3V"; NOL1="b_HLT_mu4_mu4noL1_L1MU3V"
edges=[4,5,6,8,10,15,25,50]; bins=ROOT.std.vector('double')(edges)
BEAMS=[b for b in ["pp","pn","np","nn"]
       if all(os.path.exists(D+f"Pythia_5p36TeV_{b}_hQCD_DiMu_{s}.FullSimPP24.NTUP.root") for s in SL)]
# only beams whose files are ALL re-skimmed (have trigger branches)
def has_trig(p):
    f=ROOT.TFile.Open(p); t=f.Get("HeavyIonD3PD")
    r=any(b.GetName()=="b_HLT_2mu4_L12MU3V" for b in t.GetListOfBranches()); f.Close(); return r
BEAMS=[b for b in BEAMS if all(has_trig(D+f"Pythia_5p36TeV_{b}_hQCD_DiMu_{s}.FullSimPP24.NTUP.root") for s in SL)]
print("beams with full trigger-enabled set:",BEAMS)
ISO={"pp":4/25.,"pn":6/25.,"np":6/25.,"nn":9/25.}
out={}
for mode in [0,1]:
    per={}
    comb_den=[0.0]*7; comb_di=[0.0]*7; comb_nl=[0.0]*7
    for b in BEAMS:
        den=[0.0]*7; di=[0.0]*7; nl=[0.0]*7
        for s in SL:
            xs,eff=ami[f"{b}_{s}"]; w=xs*eff/10000.0
            df=ROOT.RDataFrame("HeavyIonD3PD", D+f"Pythia_5p36TeV_{b}_hQCD_DiMu_{s}.FullSimPP24.NTUP.root") \
                .Define("ok",f"passmu(muon_pt,muon_eta,muon_quality,muon_deltaP_overP,muon_d0,muon_z0,{mode})") \
                .Define("nsel","Sum(ok)").Filter("nsel>=2") \
                .Define("subpt","subleadpt_sel(muon_pt,ok)").Filter(MU4)
            hd=df.Histo1D(ROOT.RDF.TH1DModel("d","",7,bins.data()),"subpt")
            hi=df.Filter(DI).Histo1D(ROOT.RDF.TH1DModel("i","",7,bins.data()),"subpt")
            hn=df.Filter(NOL1).Histo1D(ROOT.RDF.TH1DModel("n","",7,bins.data()),"subpt")
            hd=hd.GetValue(); hi=hi.GetValue(); hn=hn.GetValue()
            for k in range(7):
                den[k]+=w*hd.GetBinContent(k+1); di[k]+=w*hi.GetBinContent(k+1); nl[k]+=w*hn.GetBinContent(k+1)
        per[b]={"den":den,"di":di,"nl":nl}
        for k in range(7):
            comb_den[k]+=ISO[b]*den[k]; comb_di[k]+=ISO[b]*di[k]; comb_nl[k]+=ISO[b]*nl[k]
    out[mode]={"edges":edges,"per_beam":per,"beams":BEAMS,
               "comb":{"den":comb_den,"di":comb_di,"nl":comb_nl},
               "iso_complete": sorted(BEAMS)==sorted(["pp","pn","np","nn"])}
json.dump(out, open("/tmp/claude-101379/sel_mc.json","w"))
print("SEL_MC_DONE")
