import ROOT
from glob import glob
fs=glob("./v4_run_job_[0-9]/analysis_strong.root")
d=ROOT.RDataFrame("t",fs)
print(d.Count().GetValue())
d=d.Filter("lep_cut>0")
print(d.Count().GetValue())

d=d.Filter("fjet_cut>0")
print(d.Count().GetValue())

d=d.Filter("ujet_cut>0")
print(d.Count().GetValue())

d=d.Filter("MET_cut>0")
print(d.Count().GetValue())

d=d.Filter("reco_cut>0")
print(d.Count().GetValue())

d=d.Filter("REGCUT>0")
print(d.Count().GetValue())

print("good")

fs=glob("./v4_run_job_[0-9]/analysis_weak.root")
d=ROOT.RDataFrame("t",fs)
print(d.Count().GetValue())
d=d.Filter("lep_cut>0")
print(d.Count().GetValue())

d=d.Filter("fjet_cut>0")
print(d.Count().GetValue())

d=d.Filter("ujet_cut>0")
print(d.Count().GetValue())

d=d.Filter("MET_cut>0")
print(d.Count().GetValue())

d=d.Filter("reco_cut>0")
print(d.Count().GetValue())

d=d.Filter("REGCUT>0")
print(d.Count().GetValue())