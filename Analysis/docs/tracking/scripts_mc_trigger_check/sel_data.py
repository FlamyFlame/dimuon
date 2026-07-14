import ROOT, glob, json
ROOT.EnableImplicitMT(8)
ROOT.gInterpreter.Declare(open("/tmp/claude-101379/sel.h").read())
fs=sorted(glob.glob('/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part*.root'))
MU4="b_HLT_mu4_L1MU3V"; DI="b_HLT_2mu4_L12MU3V"; NOL1="b_HLT_mu4_mu4noL1_L1MU3V"
edges=[4,5,6,8,10,15,25,50]; bins=ROOT.std.vector('double')(edges)
out={}
for mode in [0,1]:
    df=ROOT.RDataFrame("HeavyIonD3PD", fs) \
        .Define("ok",f"passmu(muon_pt,muon_eta,muon_quality,muon_deltaP_overP,muon_d0,muon_z0,{mode})") \
        .Define("nsel","Sum(ok)").Filter("nsel>=2") \
        .Define("subpt","subleadpt_sel(muon_pt,ok)")
    d=df.Filter(MU4)
    hd=d.Histo1D(ROOT.RDF.TH1DModel("d","",7,bins.data()),"subpt")
    hi=d.Filter(DI).Histo1D(ROOT.RDF.TH1DModel("i","",7,bins.data()),"subpt")
    hn=d.Filter(NOL1).Histo1D(ROOT.RDF.TH1DModel("n","",7,bins.data()),"subpt")
    tot=df.Count(); nm=d.Count()
    hd=hd.GetValue(); hi=hi.GetValue(); hn=hn.GetValue()
    out[mode]={"edges":edges,"nevt_sel":tot.GetValue(),"nevt_mu4":nm.GetValue(),
               "den":[hd.GetBinContent(k) for k in range(1,8)],
               "di":[hi.GetBinContent(k) for k in range(1,8)],
               "nl":[hn.GetBinContent(k) for k in range(1,8)],
               "den_tot":hd.Integral(0,9),"di_tot":hi.Integral(0,9),"nl_tot":hn.Integral(0,9)}
    print(f"mode={mode} events>=2 sel muons: {tot.GetValue()}  of which mu4-fired: {nm.GetValue()}")
json.dump(out, open("/tmp/claude-101379/sel_data.json","w"))
print("SEL_DATA_DONE")
