#!/usr/bin/env python
# coding: utf-8
# determine nu first

import ROOT
import os
DEBUG=False

# ROOT.gSystem.Load("libDelphes.so")
ROOT.gInterpreter.ProcessLine('.L classes/DelphesClasses.cc')
from math import cos as COS
from math import pi as PI
from uncertainties import ufloat
import glob
import numpy as np
import pickle
DEBUG_NUM=100
SAVEONLY=True # do save only and not do selection and Cij. All the calcualted variable will be saved.
DUMP_TRUTH=False # dump all raw objects
DUMP_RECO=False
DUMP_WTS=False # dump all weights variation and sumwt
DUMP_META=False # dump all weights variation and sumwt
import sys
assert len(sys.argv)==2
seed=int(sys.argv[1])
print(f"Run on {seed}")
 
# if DEBUG:
#     ROOT.ROOT.DisableImplicitMT()
# else:
#     ROOT.ROOT.EnableImplicitMT(8)

# define the kernel functions ########################################
# https://root-forum.cern.ch/t/how-to-make-snapshot-of-rdataframe-that-contains-stl-collections/38838
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

def Rdefine(df,name,content,debug=DEBUG):
    if df.HasColumn(name):
        print(f"Warning redefine {name}")
        return df.Redefine(name,content)
    else:
        return df.Define(name,content)

def P_Rdefine(df,name,content,debug=DEBUG): # truth variables
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
def process_resolved(_df,debug=DEBUG):
    assert False # not consider resolved for now.
    DEF=Rdefine
    # MET cut
    _df = DEF(_df,"MET_cut",
                f"MissingET.MET.size()>0 && MissingET.MET.at(0)>30")
    # lepton
    _df = DEF(_df,"charged_lep",
                f"get_lepton(Electron.PT,Electron.Eta,Electron.Phi,Electron.Charge,Muon.PT,Muon.Eta,Muon.Phi,Muon.Charge)")
    _df = DEF(_df,"charge",
                "charged_lep.first")
    _df = DEF(_df,"lep",
                "charged_lep.second")
    _df = DEF(_df,"lep_cut",
                f"lep.Pt()>0")
    _df = DEF(_df,"lep_pt",
                "lep.Pt()")
    # print(_df.AsNumpy(["lep_pt"]))
    # exit()
    # solve nu
    _df = DEF(_df,"nu",
                "solve_nu(lep,MissingET.MET,MissingET.Phi)")
    _df = DEF(_df,"nu_cut",
                f"nu.Pt()>0")
    # _df = DEF(_df,"n",
    #             "nu.Pz()")
    # print(_df.AsNumpy(["n"]))
    # exit()
    # jet
    _df = DEF(_df,"blep_bhad_j1_j2",
                f"get_blep_bhad_j1_j2(lep,nu,Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass,Jet.BTag)")
    _df = DEF(_df,"b_lep",
                "get_jet(blep_bhad_j1_j2.at(0),Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass)")
    _df = DEF(_df,"b_had",
                "get_jet(blep_bhad_j1_j2.at(1),Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass)")
    _df = DEF(_df,"jet1",
                "get_jet(blep_bhad_j1_j2.at(2),Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass)")
    _df = DEF(_df,"jet2",
                "get_jet(blep_bhad_j1_j2.at(3),Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass)")
    _df = DEF(_df,"jet_cut",
                "blep_bhad_j1_j2.at(0)>=0")
    # _df = DEF(_df,"debug_1",
                # "blep_bhad_j1_j2.size()")
    # _df = DEF(_df,"i0",
    #             "blep_bhad_j1_j2.at(0)")
    # _df = DEF(_df,"i1",
    #             "blep_bhad_j1_j2.at(1)")
    # _df = DEF(_df,"i2",
    #             "blep_bhad_j1_j2.at(2)")
    # _df = DEF(_df,"i3",
    #             "blep_bhad_j1_j2.at(3)")
    # _df = DEF(_df,"e",
    #             "Event.Number.at(0)")
    # print(_df.AsNumpy(["e"]))
    # print(_df.AsNumpy(["i0"]))
    # print(_df.AsNumpy(["i1"]))
    # print(_df.AsNumpy(["i2"]))
    # print(_df.AsNumpy(["i3"]))
    # print(_df.AsNumpy(["debug_1"]))
    # exit()
    
    # ok then others should be fine
    _df = DEF(_df,"lep_top",
                "b_lep+nu+lep")
    _df = DEF(_df,"had_top",
                "b_had+jet1+jet2")
    _df = DEF(_df,"top",
                "charge ? lep_top : had_top")
    _df = DEF(_df,"anti_top",
                "charge ? had_top : lep_top")
    _df = DEF(_df,"had_W",
                "jet1+jet2")
    
    _df = DEF(_df,"reco_cut",
                f"jet_cut")
    
    # print(_df.AsNumpy(["reco_cut"]))
    # exit()
    return _df

def process_boosted(_df,debug=DEBUG):
    DEF=Rdefine
    # lepton
    _df = DEF(_df,"charged_lep",
                f"get_lepton(Electron.PT,Electron.Eta,Electron.Phi,Electron.Charge,Muon.PT,Muon.Eta,Muon.Phi,Muon.Charge)")
    _df = DEF(_df,"charge",
                "charged_lep.first")
    _df = DEF(_df,"lep",
                "charged_lep.second")
    _df = DEF(_df,"lep_cut",
                f"lep.Pt()>0")
    # dump variables
    _df = DEF(_df,"lep_pt",
                f"lep.Pt()")
    _df = DEF(_df,"lep_eta",
                f"lep.Eta()")
    _df = DEF(_df,"lep_phi",
                f"lep.Phi()")
    _df = DEF(_df,"lep_energy",
                f"lep.E()")
    _df = DEF(_df,"lep_charge","charge")
    
    # jet
    _df = DEF(_df,"fjet",
                "get_Fjet(FatJet.PT,FatJet.Eta,FatJet.Phi,FatJet.Mass,Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass,Jet.BTag)")
    _df = DEF(_df,"fjet_cut",
                f"fjet.Pt()>0")
    _df = DEF(_df,"ujet_cut",
                f"ujet_cut(Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass,Jet.BTag)")
    # define the fjet to export (pt, eta, phi, j)
    _df = DEF(_df,"fjet_pt",
                f"fjet.Pt()")
    _df = DEF(_df,"fjet_eta",
                f"fjet.Eta()")
    _df = DEF(_df,"fjet_phi",
                f"fjet.Phi()")
    _df = DEF(_df,"fjet_energy",
                f"fjet.E()")
        
    # MET cut
    _df = DEF(_df,"MET_cut",
                f"MissingET.MET.size()>0 && MissingET.MET.at(0)>30")
    _df = DEF(_df,"met_met",
                f"MissingET.MET.size()>0?MissingET.MET.at(0):0")
    _df = DEF(_df,"met_phi",
                f"MissingET.Phi.size()>0?MissingET.Phi.at(0):0")
    _df = DEF(_df,"met_eta",
                f"MissingET.Eta.size()>0?MissingET.Eta.at(0):0")
    
    # solve nu
    _df = DEF(_df,"nu",
                "solve_nu(lep,MissingET.MET,MissingET.Phi)")
    _df = DEF(_df,"nu_cut",
                f"nu.Pt()>0")
    #  jets    
    _df = DEF(_df,"blep_bhad_j1_j2",
                "get_blep_bhad_j1_j2_boosted(lep,nu,fjet,Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass,Jet.BTag)")
    _df = DEF(_df,"b_lep",
                "get_jet(blep_bhad_j1_j2.at(0),Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass)")
    _df = DEF(_df,"b_had",
                "get_jet(blep_bhad_j1_j2.at(1),Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass)")
    _df = DEF(_df,"jet1",
                "get_jet(blep_bhad_j1_j2.at(2),Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass)")
    _df = DEF(_df,"jet2",
                "get_jet(blep_bhad_j1_j2.at(3),Jet.PT,Jet.Eta,Jet.Phi,Jet.Mass)")
    _df = DEF(_df,"jet_cut",
                "blep_bhad_j1_j2.at(0)>=0")
    
    _df = DEF(_df,"lep_top",
                "b_lep+nu+lep")
    _df = DEF(_df,"had_top",
                "b_had+jet1+jet2")
    # _df = DEF(_df,"had_top",
    #             "fjet")
    _df = DEF(_df,"top",
                "charge ? lep_top : had_top")
    _df = DEF(_df,"anti_top",
                "charge ? had_top : lep_top")
    
    _df = DEF(_df,"had_W",
                "jet1+jet2")
    
    _df = DEF(_df,"reco_cut",
                f"jet_cut")
    
    # export the jets to use
    _df = DEF(_df,"jet_pt",
                f"Jet.PT")
    _df = DEF(_df,"jet_eta",
                f"Jet.Eta")
    _df = DEF(_df,"jet_phi",
                f"Jet.Phi")
    _df = DEF(_df,"jet_mass",
                f"Jet.Mass")
    
    _df = DEF(_df,"jet_n",
                f"Jet.PT.size()")
    
    return _df

def process_parton(df,debug=DEBUG):
    DEF=P_Rdefine
    suffix="truth_"
    
    df = DEF(df,"top_index",
                "get_top(Particle_size, 0, Particle.PID, Particle.Status)")
    df = DEF(df,"anti_top_index",
            "get_top(Particle_size, 1, Particle.PID, Particle.Status)")
    
    df = DEF(df,"top",
                f"get_particle_idx({suffix}top_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    df = DEF(df,"anti_top",
                f"get_particle_idx({suffix}anti_top_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    
    # now there are a lot of partons!! 
    df = DEF(df,"lep_index", "find_lep(Particle_size, Particle.PID,Particle.Status,Particle.M1,Particle.M2)")
    # t --> b W+ --> b j j, t~ --> b~ W- --> b~ l nu, charge<0 PID>0 top_is_had
    # t~ --> b~ W- --> b~ j j, t --> b W+ --> b l~ nu
    # print(df.AsNumpy(["truth_top_index","truth_anti_top_index","truth_lep_index"]))
    df = DEF(df,"top_is_had", f"Particle.PID.at({suffix}lep_index)>0")
    df = DEF(df,"had_W_index", f"get_W(Particle_size, {suffix}top_is_had, Particle.PID, Particle.Status)")
    df = DEF(df,"lep_W_index", f"get_W(Particle_size, !{suffix}top_is_had, Particle.PID, Particle.Status)")
    df = DEF(df,"lep_top_index", f"{suffix}top_is_had?{suffix}anti_top_index:{suffix}top_index")
    df = DEF(df,"had_top_index", f"{suffix}top_is_had?{suffix}top_index:{suffix}anti_top_index")
    df = DEF(df,"lep", f"get_particle_idx({suffix}lep_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    df = DEF(df,"had_W", f"get_particle_idx({suffix}had_W_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    df = DEF(df,"had_top", f"get_particle_idx({suffix}had_top_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    df = DEF(df,"lep_top", f"get_particle_idx({suffix}lep_top_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    
    df = DEF(df,"ujet_index", "find_u(Particle_size,Particle.PID,Particle.Status,Particle.M1,Particle.M2)")
    df = DEF(df,"djet_index", "find_d(Particle_size,Particle.PID,Particle.Status,Particle.M1,Particle.M2)")
    
    df = DEF(df,"j1_index", f"Particle.PT.at({suffix}ujet_index)>Particle.PT.at({suffix}djet_index)?{suffix}ujet_index:{suffix}djet_index")
    df = DEF(df,"j2_index", f"Particle.PT.at({suffix}ujet_index)<=Particle.PT.at({suffix}djet_index)?{suffix}ujet_index:{suffix}djet_index")
    
    df = DEF(df,"jet1", f"get_particle_idx({suffix}j1_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    df = DEF(df,"jet2", f"get_particle_idx({suffix}j2_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    
    df = DEF(df,"ujet", f"get_particle_idx({suffix}ujet_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    df = DEF(df,"djet", f"get_particle_idx({suffix}djet_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    # for dump
    df = DEF(df,"nu_index",
                f"find_nu(Particle_size,Particle.PID,Particle.Status,Particle.M1,Particle.M2)")
    df = DEF(df,"nu", f"get_particle_idx({suffix}nu_index, Particle.E, Particle.PT, Particle.Eta, Particle.Phi)")
    
    df = DEF(df,"nu_pt", f"{suffix}nu.Pt()")
    df = DEF(df,"nu_eta", f"{suffix}nu.Eta()")
    df = DEF(df,"nu_phi", f"{suffix}nu.Phi()")
    df = DEF(df,"nu_e", f"{suffix}nu.E()")
    
    return df
    
def define_presel(_df,is_parton,debug=DEBUG):
    suffix="" if not is_parton else "truth_" 
    DEF=Rdefine if not is_parton else P_Rdefine
    _df = DEF(_df,"ttbar", f"{suffix}top + {suffix}anti_top")
    _df = DEF(_df,"cos_theta", f"get_cos_theta({suffix}top,{suffix}ttbar)")
    _df = DEF(_df,"m_tt", f"({suffix}top + {suffix}anti_top).M()")
    _df = DEF(_df,"p_tt", f"({suffix}top + {suffix}anti_top).P()")
    _df = DEF(_df,"abs_beta", f"abs({suffix}p_tt/{suffix}m_tt)")
    return _df

Rthres=[
    [(COS(PI/20*1),COS(PI/20*0)),(300,600)],
    [(COS(PI/20*2),COS(PI/20*1)),(300,500)],
    [(COS(PI/20*3),COS(PI/20*2)),(300,500)],
    [(COS(PI/20*4),COS(PI/20*3)),(300,400)],
    [(COS(PI/20*5),COS(PI/20*4)),(300,400)],
    [(COS(PI/20*6),COS(PI/20*5)),(300,400)],
    [(COS(PI/20*7),COS(PI/20*6)),(300,400)],
    [(COS(PI/20*8),COS(PI/20*7)),(300,400)],
    [(COS(PI/20*9),COS(PI/20*8)),(300,400)],
    [(COS(PI/20*10),COS(PI/20*9)),(300,400)],
]

Rboost=[
    [(COS(PI/20*6),COS(PI/20*5)),(1100,1600)],
    [(COS(PI/20*7),COS(PI/20*6)),(800,1600)],
    [(COS(PI/20*8),COS(PI/20*7)),(700,1600)],
    [(COS(PI/20*9),COS(PI/20*8)),(700,1600)],
    [(COS(PI/20*10),COS(PI/20*9)),(700,1600)],
]

Rweak=[
    [(COS(PI/20*9),COS(PI/20*8.5)),(1200,1600)],
    [(COS(PI/20*10),COS(PI/20*9)),(1000,1600)],
]

Rstrong=[
    [(COS(PI/20*9),COS(PI/20*8.5)),(1500,1600)],
    [(COS(PI/20*10),COS(PI/20*9)),(1050,1600)],
]
    
def doRegionSel(d,sel,sum_evt,sum_weight,is_parton,beta_cut, lumi_ref=139, add_only=False,debug=DEBUG):
    suffix="" if not is_parton else "truth_" 
    DEF=Rdefine if not is_parton else P_Rdefine
    d=define_presel(d,is_parton=is_parton)
    s=[f"(({suffix}cos_theta>={k[0][0]})*({suffix}cos_theta<{k[0][1]})*({suffix}m_tt>={k[1][0]})*({suffix}m_tt<{k[1][1]}))" for k in sel]
    s="+".join(s)
    s=f"({s})>0"
    if beta_cut:
        s=f"({s})&&({suffix}abs_beta<=0.9)"
    if not add_only and not is_parton:
        # apply the reco region cut
        s=f"({s})&&(reco_cut)"
        before_evt=d.Count().GetValue()
        before_weight=d.Sum('w').GetValue()
    
    d=DEF(d,"REGCUT", s)
    if add_only:
        return d
    # use lambda
    
    d=d.Filter(f"{suffix}REGCUT")
    
    if is_parton:
        print(f"### Phase space (truth) cut applied: {s}")
        print(f"Ratio {d.Count().GetValue()}/{sum_evt}={d.Count().GetValue()/sum_evt}")
    else:
        print(f"### Reconstruction cut applied: {s}")
        print(f"Efficiency= {d.Count().GetValue()/before_evt}")
        print(f"Efficiency(weighted)= {d.Sum('w').GetValue()/before_weight}")
    xs=d.Sum('w').GetValue()/sum_evt
    print(f"Eqiv XS {xs} pb, {xs*1000*lumi_ref:.2e}@{lumi_ref}/fb")
    return d

def define_var(df,is_parton,debug=DEBUG):
    suffix="" if not is_parton else "truth_" 
    DEF=Rdefine if not is_parton else P_Rdefine
    df = DEF(df,"helix_frame", f"get_top_frame({suffix}top,{suffix}ttbar)")
    df = DEF(df,"cos_theta_W", f"CalculatePoptThetaW({suffix}had_W, {suffix}had_top, {suffix}ttbar, {suffix}jet1, {suffix}jet2)")

    for axis in ["k", "r", "n"]:
        df = DEF(df,f"cos_theta_lep_{axis}",
                    f'get_theta_lep("{axis}", {suffix}lep_top, {suffix}helix_frame, {suffix}ttbar, {suffix}lep)') 
        df = DEF(df,f"cos_theta_opt_{axis}",
                    f'get_theta_jet_opt("{axis}", {suffix}had_top, {suffix}helix_frame, {suffix}ttbar, {suffix}cos_theta_W, {suffix}jet1, {suffix}jet2)') 
        # ok let's save all
        if is_parton:
            df = DEF(df,f"cos_theta_jet_{axis}",
                    f'get_theta_jet("{axis}", {suffix}had_top, {suffix}helix_frame, {suffix}ttbar, {suffix}djet)') 
        df = DEF(df,f"cos_theta_soft_{axis}",
            f'get_theta_jet("{axis}", {suffix}had_top, {suffix}helix_frame, {suffix}ttbar, {suffix}jet2)') 
    return df

def define_obs(df,op,is_parton,debug=DEBUG):
    DEF=Rdefine if not is_parton else P_Rdefine
    for n,v in op.items():
        df=DEF(df, n, v)
    return df

################################ ok then define the calculation
def define_cosAcosB(df,is_parton,typ=["opt","jet","soft"],debug=DEBUG):
    suffix="" if not is_parton else "truth_" 
    op={
    }
    for i in ["n","r","k"]:
        for j in ["n","r","k"]:
            for t in typ:
                if t=="jet" and not is_parton: continue
                op.update({
                    f"cosAcosB_{i}{j}_{t}":f"{suffix}cos_theta_lep_{i}*{suffix}cos_theta_{t}_{j}",
                })
    df=define_obs(df,op,is_parton=is_parton)
    return df

def getAsymetry(n,d,sum_evt, sum_weight, lumi_fb, eff=1.0, k_factor=1.8,use_weight=True,debug=DEBUG): # note the input is pb
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

def getMean(n,d,factor=1.,debug=DEBUG):
    return d.Mean(n).GetValue()*factor

def getCij(Asy,KA,KB,debug=DEBUG): # note the input is pb # always top antitop so -1*ka*kb
    return -1.*Asy*4/KA/KB

def doCij(df,lumi,sum_evt,sum_weight,use_weight, typ,is_parton,eff=1.0,debug=DEBUG):
    suffix="" if not is_parton else "truth_" 
    Cij={}
    KB={
        "opt":0.64,
        "jet" :1.0,
        "soft":0.5,   
        }
    for i in ["n","r","k"]:
        for j in ["n","r","k"]:
            for t in typ:
                if t=="jet" and not is_parton: continue
                Cij[f'{i}{j}']=getCij(getAsymetry(
                    f"{suffix}cosAcosB_{i}{j}_{t}",df,
                    sum_evt=sum_evt, sum_weight=sum_weight,
                    use_weight=use_weight,
                    eff=eff,lumi_fb=lumi),KA=1.0,KB=KB[t])
    return Cij

def printCij(C):
    print(f"    n       \t        r\t        k")
    print(f"n: {C['nn']}\t{C['nr']}\t{C['nk']}")
    print(f"r: {C['rn']}\t{C['rr']}\t{C['rk']}")
    print(f"k: {C['kn']}\t{C['kr']}\t{C['kk']}")

# ## Individual QE variables
# Thresold
# 2C=2*(-(C11+C22+C33)-1)
# boosted
# 2C=2*(-(-C11-C22+C33)-1)
# 1 2 3 = r k n
def calQE_individual(Cij,boosted,debug=DEBUG): # now we first implment the individual method, from Cij. Remain later the direct way
    if boosted:
        return -1*(-Cij["rr"]-Cij["kk"]+Cij["nn"])-1
    else:
        return -1*(Cij["rr"]+Cij["kk"]+Cij["nn"])-1

# B=sqrt(2)
def calBel_individual(Cij,debug=DEBUG): # now we first implment the individual method, from Cij. Remain later the direct way
    return -Cij["nn"]+Cij["rr"]-2**0.5
    
def doQE_individual(df,sum_evt,sum_weight,boosted,is_parton,lumi=139,typ=["opt","jet","soft"],use_weight=True,debug=DEBUG): # now we first implment the individual method, from Cij. Remain later the direct way
    Cij=doCij(df,lumi=lumi,sum_evt=sum_evt,sum_weight=sum_weight,use_weight=use_weight,typ=typ,is_parton=is_parton)
    C=calQE_individual(Cij,boosted=boosted)
    B=calBel_individual(Cij)
    return Cij,C,B

####################################### ok let' go
def read_files(FS,debug=DEBUG):
    _df = ROOT.RDataFrame("Delphes",FS)
    # ROOT.RDF.Experimental.AddProgressBar(_df)
    if debug:
        _df=_df.Range(0,DEBUG_NUM)
    EVTN=_df.Count().GetValue()
    _df = Rdefine(_df,"w", "Event.Weight.at(0)")
    WEIGHT=_df.Sum("w").GetValue()
    print(f"Loaded {EVTN} enents, sumwt={WEIGHT}, XS={WEIGHT/EVTN}pb")
    return _df,EVTN,WEIGHT

####################################### storage
def build_save_list(typ=["opt","jet","soft"],has_fatjet=False,dump_truth=DUMP_TRUTH,dump_reco=DUMP_RECO,dump_meta=DUMP_META,dump_weights=DUMP_WTS, debug=DEBUG):
    branchList = ROOT.std.vector('std::string')()
    to_save_parton=[
        # regiosn varuables
        "cos_theta",
        "m_tt",
        "p_tt",
        "abs_beta",
        "REGCUT",
        # charge
        "top_is_had",
        # indices
        "lep_index",
        "had_W_index",
        "lep_W_index",
        "lep_top_index",
        "had_top_index",
        "ujet_index",
        "djet_index",
        "j1_index",
        "j2_index",
        # selected object
        # "top",
        # "anti_top",
        # "lep",
        # "had_W",
        # "had_top",
        # "lep_top",
        # "jet1",
        # "jet2",
        # "ujet",
        # "djet",
        # #
        # "helix_frame",
        # "cos_theta_W",
        # dump
        "nu_pt",
        "nu_eta",
        "nu_phi",
        "nu_e",
    ]
    for i in ["n","r","k"]:
        to_save_parton.append(f"cos_theta_lep_{i}")
        for t in typ:
            to_save_parton.append(f"cos_theta_{t}_{i}") # others could be soft or jet(truth d jet)
    for i in ["n","r","k"]:
        for j in ["n","r","k"]:
            for t in typ:
                to_save_parton.append(f"cosAcosB_{i}{j}_{t}")
            
    suffix="truth_"
    to_save_parton=[f"{suffix}{i}" for i in to_save_parton]
    if dump_truth:
        to_save_parton+=["Particle"]
        
    to_save_reco=[
        # processed variables: region
        "cos_theta",
        "abs_beta", # no need for this but for debug
        "m_tt",
        "p_tt",
        "REGCUT",
        # charge
        "charge",
        # # processed objects
        # "lep",
        # "nu",
        # "lep_top",
        # "had_top",
        # "top",
        # "anti_top",
        # "had_W",
        # "blep_bhad_j1_j2",
        # # processed variables: cut
        "MET_cut",
        "lep_cut",
        # "jet_cut",
        # "nu_cut",
        "reco_cut",
        # #
        # "helix_frame",
        # "cos_theta_W",
    ]
    if has_fatjet:
        to_save_reco+=[
            "fjet_cut",
            "ujet_cut",
            # "fjet",
            # the dump variables
            "lep_pt",
            "lep_eta",
            "lep_phi",
            "lep_energy",
            "lep_charge",
            "fjet_pt",
            "fjet_eta",
            "fjet_phi",
            "fjet_energy",
            "met_met",
            "met_phi",
            "met_eta",
            "jet_pt",
            "jet_eta",
            "jet_phi",
            "jet_mass",
            "jet_n",
        ]
    
    for i in ["n","r","k"]:
        to_save_reco.append(f"cos_theta_lep_{i}")
        for t in typ:
            if t=="jet": continue # no truth d in reco level
            to_save_reco.append(f"cos_theta_{t}_{i}")
    for i in ["n","r","k"]:
        for j in ["n","r","k"]:
            for t in typ:
                if t=="jet": continue # no truth d in reco level
                to_save_reco.append(f"cosAcosB_{i}{j}_{t}")
            
    if dump_reco:
        to_save_reco+=[
            # raw reco objects
            "MissingET",
            "Jet",
            "Electron",
            "Muon",
        ]
        if has_fatjet:
            to_save_reco+=["FatJet"]
    
    to_save=to_save_parton + to_save_reco
    if dump_meta:
        to_save+=[
            "Event",
        ]
    
    if dump_weights:
        to_save+=[
            "Weight",
        ]
    
    to_save.append("w")
    
    for i in to_save:
        branchList.push_back(i)
    return branchList

############### 
def do_resolved_boosted(FS,is_boosted, save=False, debug=DEBUG): # note you need to do the parton and reco at the same time to derive the efficiency
    N="boosted" if is_boosted else "threshold"
    R=Rboost if is_boosted else Rthres
    
    df,EVT,SWT=read_files(FS)
    df=process_parton(df)
    #now you have all the parton level variable
    df=doRegionSel(df,R,sum_evt=EVT,sum_weight=SWT,is_parton=True,beta_cut=not is_boosted,add_only=SAVEONLY)
    
    df=define_var(df,is_parton=True)
    df=define_cosAcosB(df,is_parton=True)
    
    if not SAVEONLY:
        print(f"{N} weighted parton Cij: truth")
        Cij,C,_ = doQE_individual(df,boosted=is_boosted,sum_evt=EVT,sum_weight=SWT,is_parton=True) # no need to apply eff again
        printCij(Cij)
        print(f"QE observable ({N}) truth: {C}, reference value ")
    #now do the reco level
    df=process_resolved(df)
    df=doRegionSel(df,R,sum_evt=EVT,sum_weight=SWT,is_parton=False,beta_cut=not is_boosted,add_only=SAVEONLY)
    # print(f"Reference eff= evt@139/fb=")
    
    df=define_var(df,is_parton=False)
    df=define_cosAcosB(df,is_parton=False)

    if not SAVEONLY:
        print(f"{N} weighted parton Cij")
        Cij,C,_ = doQE_individual(df,boosted=is_boosted,sum_evt=EVT,sum_weight=SWT,is_parton=False) # no need to apply eff again
        printCij(Cij)
        
        # reference_Cij='''
        # parton-level  n              r              k

        # '''
        # print("reference Cij")
        # print(reference_Cij)

        print(f"QE observable ({N}): {C}, reference value ")
    if SAVEONLY or save:
        df.Snapshot("p",f"analysis_{N}.root",build_save_list())


def do_weak_strong(FS,is_strong, save=False,debug=DEBUG): # note you need to do the parton and reco at the same time to derive the efficiency
    N="strong" if is_strong else "weak"
    R=Rstrong if is_strong else Rweak
    
    df,EVT,SWT=read_files(FS)
    df=process_parton(df)
    #now you have all the parton level variable
    df=doRegionSel(df,R,sum_evt=EVT,sum_weight=SWT,is_parton=True,beta_cut=False,lumi_ref=300,add_only=SAVEONLY)
    
    df=define_var(df,is_parton=True)
    df=define_cosAcosB(df,is_parton=True)
    if not SAVEONLY:
        print(f"{N} weighted parton Cij: truth")
        Cij,_,B = doQE_individual(df,boosted=True,sum_evt=EVT,sum_weight=SWT,is_parton=True) # no need to apply eff again
        printCij(Cij)
        print(f"BELL observable ({N}) truth: {B}, reference value ")
    
    #now do the reco level
    df=process_boosted(df)
    df=doRegionSel(df,R,sum_evt=EVT,sum_weight=SWT,is_parton=False,beta_cut=False,lumi_ref=300,add_only=SAVEONLY)
    # print(f"Reference eff= evt@139/fb=")
    
    df=define_var(df,is_parton=False)
    df=define_cosAcosB(df,is_parton=False)
    if not SAVEONLY:
        print(f"{N} weighted parton Cij")
        Cij,_,B = doQE_individual(df,boosted=True,sum_evt=EVT,sum_weight=SWT,is_parton=False) # no need to apply eff again
        printCij(Cij)
        
        # reference_Cij='''
        # parton-level  n              r              k

        # '''
        # print("reference Cij")
        # print(reference_Cij)

        print(f"QE observable ({N}): {B}, reference value ")
    if SAVEONLY or save:
        df.Snapshot("t",f"analysis_{N}.root",build_save_list(has_fatjet=True))
        
# FS_inclusive=list(glob.glob(f"/lustre/collider/liuqibin/MC/QE_bulk_resolved/run_job_{seed}/*.reco.root"))
FS_boosted=list(glob.glob(f"<PATH_TO_DELPHES_OUTPUT>/run_job_{seed}/*.reco.root"))
print("############# resolved(threshold) ############")
# do_resolved_boosted(FS_inclusive,is_boosted=False)
print("############# boosted ############")
# do_resolved_boosted(FS_inclusive,is_boosted=True)
print("############# weak ############")
do_weak_strong(FS_boosted,is_strong=False)
print("############# strong ############")
do_weak_strong(FS_boosted,is_strong=True)

print("ALL DONE")
