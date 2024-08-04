#!/usr/bin/env python
# coding: utf-8

import ROOT
ROOT.gInterpreter.ProcessLine('.L classes/DelphesClasses.cc')
from math import cos as COS
from math import pi as PI
from uncertainties import ufloat
import glob
import os
import numpy as np
import pickle
DEBUG=False
DEBUG_NUM=100000
 
if DEBUG:
    ROOT.ROOT.DisableImplicitMT()
else:
    # ROOT.ROOT.DisableImplicitMT()
    ROOT.ROOT.EnableImplicitMT(16)

# define the kernel functions ########################################
if not os.path.exists("kernel_common_h.so"):
    ROOT.gInterpreter.ProcessLine('.L kernel_common.h+') # do only once
else:
    ROOT.gSystem.Load("kernel_common_h.so") # for all and no portable way.
    
with open("./kernel_common.C") as f:
    ROOT.gInterpreter.Declare(f.read())
with open("./kernel_truth_ana.C") as f:
    ROOT.gInterpreter.Declare(f.read())
with open("./kernel_reco_ana.C") as f:
    ROOT.gInterpreter.Declare(f.read())

def Rdefine(df,name,content):
    if df.HasColumn(name):
        print(f"Warning redefine {name}")
        return df.Redefine(name,content)
    else:
        return df.Define(name,content)

def P_Rdefine(df,name,content): # truth variables
    name=f"truth_{name}"
    if df.HasColumn(name):
        print(f"Warning redefine {name}")
        return df.Redefine(name,content)
    else:
        return df.Define(name,content)

# note since we need both truth and reco, should not fileter out, but just set cut flag
# def cut_comb(cts): # only combine the last one should be enough
#     return "&&".join([f"({c})" for c in cts])
# def cut(d,c,n,f_def,previous_cuts):
#     cts=previous_cuts+[c]
#     n=f"c{len(cts)}_{n}"
#     return f_def(d,n,cut_comb(cts)),n,cts

# http://cp3.irmp.ucl.ac.be/downloads/RootTreeDescription.html
# https://pythia.org/latest-manual/ParticleProperties.html
   
def doRegionSel(d,sum_evt,sum_weight,is_parton,lumi_ref=139,plot=False):
    suffix="" if not is_parton else "truth_" 
    DEF=Rdefine if not is_parton else P_Rdefine
    s=""
    if not plot:
        s=f"{suffix}REGCUT"
        if not is_parton:
            s=f"({s})&&(reco_cut)"
            before_evt=d.Count().GetValue()
            before_weight=d.Sum('w').GetValue()
    else:
        s="reco_cut"
    hs=[]
    if plot:
        suffix="truth_" 
        DEF=P_Rdefine
        d=DEF(d,"theta",f"acos({suffix}cos_theta)")
        print("before",d.Sum('w').GetValue())
        hs.append(d.Histo2D(
                ROOT.RDF.TH2DModel("mtt_theta_before", f"mtt-theta(before);theta(before);mtt", 20, 0, +PI, 17, 300, 2000,),
                f"{suffix}theta",f"{suffix}m_tt","w"))
    d=d.Filter(s)
    if plot:
        print("after",d.Sum('w').GetValue())
        hs.append(d.Histo2D(
                ROOT.RDF.TH2DModel("mtt_theta_after", f"mtt-theta(after);theta(after);mtt", 20, 0, +PI, 17, 300, 2000,),
                f"{suffix}theta",f"{suffix}m_tt","w"))
        return d,hs
    
    if is_parton:
        print(f"### Phase space (truth) cut applied: {s}")
        print(f"Ratio {d.Count().GetValue()}/{sum_evt}={d.Count().GetValue()/sum_evt}")
    else:
        print(f"### Reconstruction cut applied: {s}")
        print(f"Efficiency= {d.Count().GetValue()/before_evt}")
        print(f"Efficiency(weighted)= {d.Sum('w').GetValue()/before_weight}")
    xs=d.Sum('w').GetValue()/sum_evt
    print(f"Eqiv XS {xs} pb, {xs*1000*lumi_ref:.2e}@{lumi_ref}/fb")
    return d,hs

def define_obs(df,op,is_parton):
    DEF=Rdefine if not is_parton else P_Rdefine
    for n,v in op.items():
        df=DEF(df, n, v)
    return df

################################ ok then define the calculation
def getAsymetry(n,d,sum_evt, sum_weight, lumi_fb, eff=1.0, k_factor=1.8,use_weight=True,debug=False): # note the input is pb
    if use_weight:
        # the weight sum to cross-section
        Np=d.Filter(f"{n}>0").Sum("w").GetValue()/sum_evt*1000*lumi_fb*k_factor*eff
        Nn=d.Filter(f"{n}<0").Sum("w").GetValue()/sum_evt*1000*lumi_fb*k_factor*eff
    else:
        Np=d.Filter(f"{n}>0").Count().GetValue()/sum_evt*sum_weight/sum_evt*1000*lumi_fb*k_factor*eff
        Nn=d.Filter(f"{n}<0").Count().GetValue()/sum_evt*sum_weight/sum_evt*1000*lumi_fb*k_factor*eff
    Np=ufloat(Np,Np**0.5)
    Nn=ufloat(Nn,Nn**0.5)
    if debug:
        print(f"Asy analyzing... Np={Np}, Nn={Nn}")
    denom=Np+Nn if Np+Nn > 0 else 1
    return 1.*(Np-Nn)/denom

def getMean(n,d,factor=1.):
    return d.Mean(n).GetValue()*factor

def getCij(Asy,KA,KB): # note the input is pb # always top antitop so -1*ka*kb
    return -1.*Asy*4/KA/KB

def doCij(df,lumi,sum_evt,sum_weight,use_weight, typ,is_parton,eff=1.0):
    suffix="" if not is_parton else "truth_" 
    Cij={}
    KB={
        "opt":0.64,
        "jet" :1.0,
        "soft":0.5,   
        }
    for i in ["n","r","k"]:
        for j in ["n","r","k"]:
            Cij[f'{i}{j}']=getCij(getAsymetry(
                f"{suffix}cosAcosB_{i}{j}_{typ}",df,
                sum_evt=sum_evt, sum_weight=sum_weight,
                use_weight=use_weight,
                eff=eff,lumi_fb=lumi),KA=1.0,KB=KB[typ])
    return Cij

def printCij(C):
    print("###### Cij #########")
    print(f"    n       \t        r\t        k")
    print(f"n: {C['nn']}\t{C['nr']}\t{C['nk']}")
    print(f"r: {C['rn']}\t{C['rr']}\t{C['rk']}")
    print(f"k: {C['kn']}\t{C['kr']}\t{C['kk']}")
    print("####################")

# ## Individual QE variables
# Thresold
# 2C=2*(-(C11+C22+C33)-1)
# boosted
# 2C=2*(-(-C11-C22+C33)-1)
# 1 2 3 = r k n
def calQE_individual(Cij,boosted): # now we first implment the individual method, from Cij. Remain later the direct way
    if boosted:
        return -1*(-Cij["rr"]-Cij["kk"]+Cij["nn"])-1
    else:
        return -1*(Cij["rr"]+Cij["kk"]+Cij["nn"])-1

# B=sqrt(2)
def calBel_individual(Cij): # now we first implment the individual method, from Cij. Remain later the direct way
    return -Cij["nn"]+Cij["rr"]-2**0.5
    
def doQE_individual(df,sum_evt,sum_weight,boosted,is_parton,lumi=139,typ="opt",use_weight=True): # now we first implment the individual method, from Cij. Remain later the direct way
    Cij=doCij(df,lumi=lumi,sum_evt=sum_evt,sum_weight=sum_weight,use_weight=use_weight,typ=typ,is_parton=is_parton)
    C=calQE_individual(Cij,boosted=boosted)
    B=calBel_individual(Cij)
    return Cij,C,B

####################################### ok let' go
def read_files(FS,debug=DEBUG):
    _df = ROOT.RDataFrame("p",FS)
    # ROOT.RDF.Experimental.AddProgressBar(_df)
    if debug:
        _df=_df.Range(0,DEBUG_NUM)
    EVTN=_df.Count().GetValue()
    WEIGHT=_df.Sum("w").GetValue()
    print(f"Loaded {EVTN} enents, sumwt={WEIGHT}, XS={WEIGHT/EVTN}pb")
    return _df,EVTN,WEIGHT

################################################################# resolved
def do_QE_samples(FS, prefix="threshold", save=False): # note you need to do the parton and reco at the same time to derive the efficiency
    df,EVT,SWT=read_files(FS)
    df,_=doRegionSel(df,sum_evt=EVT,sum_weight=SWT,is_parton=True)
    Cij,C,_ = doQE_individual(df,boosted=("boost" in prefix),sum_evt=EVT,sum_weight=SWT,is_parton=True) # no need to apply eff again
    printCij(Cij)
    print(f"QE observable ({prefix}) truth: {C}, reference value ")
    #now do the reco level
    df,_=doRegionSel(df,sum_evt=EVT,sum_weight=SWT,is_parton=False)
    # print(f"Reference eff= evt@139/fb=")
    Cij,C,_ = doQE_individual(df,boosted=("boost" in prefix),sum_evt=EVT,sum_weight=SWT,is_parton=False) # no need to apply eff again
    printCij(Cij)

    print(f"QE observable ({prefix}): {C}, reference value ")
    
    # some plotting?
    
################################################################# boosted
def do_BEL_sample(FS, prefix="weak", save=False): # note you need to do the parton and reco at the same time to derive the efficiency
    df,EVT,SWT=read_files(FS)
    df,_=doRegionSel(df,sum_evt=EVT,sum_weight=SWT,is_parton=True)
    Cij,_,B = doQE_individual(df,boosted=True,sum_evt=EVT,sum_weight=SWT,is_parton=True) # no need to apply eff again
    printCij(Cij)
    print(f"BELL observable ({prefix}) truth: {B}, reference value ")
    #now do the reco level
    df,_=doRegionSel(df,sum_evt=EVT,sum_weight=SWT,is_parton=False)
    # print(f"Reference eff= evt@139/fb=")
    Cij,_,B = doQE_individual(df,boosted=True,sum_evt=EVT,sum_weight=SWT,is_parton=False) # no need to apply eff again
    printCij(Cij)

    print(f"BELL observable ({prefix}): {B}, reference value ")
    
    # some plotting?
###################################################################
# plotting of selection efficiency
def do_QE_eff_plotting(FS,n=""): # note you need to do the parton and reco at the same time to derive the efficiency
    # denom: all
    # nom: reco_cut only
    df,EVT,SWT=read_files(FS)
    _,hs=doRegionSel(df,sum_evt=EVT,sum_weight=SWT,is_parton=False,plot=True)
    f=ROOT.TFile(f"{n}resolve_analysis_eff.root","recreate")
    # should not runGraphs!! since the tow plotting could not be run in parallel(sequential!)
    for i,h in enumerate(hs):
        h.GetValue().Write()
    f.Close()
        
def do_BELL_eff_plotting(FS,n=""): # note you need to do the parton and reco at the same time to derive the efficiency
    # denom: all
    # nom: reco_cut only
    df,EVT,SWT=read_files(FS)
    _,hs=doRegionSel(df,sum_evt=EVT,sum_weight=SWT,is_parton=False,plot=True)
    f=ROOT.TFile(f"{n}boosted_analysis_eff.root","recreate")
    for i,h in enumerate(hs):
        h.GetValue().Write()
    f.Close()

FS_inclusive=list(glob.glob(f"<PATH_TO_ANALYSIS_OUTPUT>/run_job_[1-2][0-9][0-9]/analysis_threshold.root"))
FS_boosted=list(glob.glob(f"<PATH_TO_ANALYSIS_OUTPUT>/run_job_[1-2][0-9][0-9]/analysis_boosted.root"))
FS_weak=list(glob.glob(f"<PATH_TO_ANALYSIS_OUTPUT>/run_job_[1-2][0-9][0-9]/analysis_weak.root"))
FS_strong=list(glob.glob(f"<PATH_TO_ANALYSIS_OUTPUT>/run_job_[1-2][0-9][0-9]/analysis_strong.root"))
# print("############# resolved(threshold) ############")
# do_QE_samples(FS_inclusive,"thresold")
# print("############# boosted ############")
# do_QE_samples(FS_boosted,"boosted")
# print("############# weak ############")
# do_BEL_sample(FS_weak,"weak")
# print("############# strong ############")
# do_BEL_sample(FS_strong,"strong")

# do_QE_eff_plotting(FS_inclusive)
# do_BELL_eff_plotting(FS_weak)
# cross-check
# do_QE_eff_plotting(FS_inclusive,"check_")
# do_BELL_eff_plotting(FS_weak,"check_")

print("ALL DONE")
