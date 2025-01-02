# draw delayedE for both reactor IBD and accidental
import ROOT
import sys
sys.path.append('/home/huangp/BiPo/Analysis')
from PyROOT_Style import set_style
from PyROOT_Style import Line_Draw_Setting
from PyROOT_Style import TLegend_Setting
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
matplotlib.rcParams.update({'font.size': 22})
set_style()

    




if __name__ == "__main__":

    
    # new ttree
    IBDT = ROOT.TChain("IBDT")
    
    AccT = ROOT.TChain("AccT")
    IBDT.Add("/home/huangp/AntiNu/Analysis/AccIBDClassifierold.root")
    AccT.Add("/home/huangp/AntiNu/Analysis/AccIBDClassifierold.root")
    Expected_Events = {"ReactorIBD": 27.9,"Accidental":0.9325}

    '''
    IBDh1 = ROOT.TH1F("IBDh1","",170,-7.,10.)
    Acch1  = ROOT.TH1F("Acch1","",270,-20.,7)
    var = "post_prob"
    IBDT.Project("IBDh1",var); AccT.Project("Acch1",var)
    IBDh1.Scale(Expected_Events["ReactorIBD"]/float(IBDT.GetEntries()), "nosw2")
    Acch1.Scale(Expected_Events["Accidental"]/float(AccT.GetEntries()), "nosw2")

    '''
    post_cut = np.arange(-8,4.1,0.1)
    FOM = np.zeros(len(post_cut))
    
    # calculate FOM = S/(S+B)**0.5
    scaleibd = Expected_Events["ReactorIBD"]/float(IBDT.GetEntries())
    scaleacc = Expected_Events["Accidental"]/float(AccT.GetEntries())
    maxFOMindex = 0; maxFOM = 0
    for icut in range(len(post_cut)):
        cut = f"post_prob > {post_cut[icut]}" 
        print(cut)
        IBDEntries = scaleibd*float(IBDT.GetEntries(cut))
        AccEntries = scaleacc*float(AccT.GetEntries(cut))
        FOM_temp = IBDEntries/(AccEntries+IBDEntries)**0.5
        if maxFOM < FOM_temp:
            maxFOM      = FOM_temp
            maxFOMindex = icut
            print(FOM_temp,maxFOMindex)
        FOM[icut] = FOM_temp

    print(f" Best FOM occurs at {post_cut[maxFOMindex]}, with value {FOM[maxFOMindex]}")

    EffIBD = float(IBDT.GetEntries(f"post_prob > {post_cut[maxFOMindex]}"))/float(IBDT.GetEntries())
    EffAcc = float(AccT.GetEntries(f"post_prob > {post_cut[maxFOMindex]}"))/float(AccT.GetEntries())
    print(f" Eff of IBD at {post_cut[maxFOMindex]} is {EffIBD} ,whereas that of Acc is  {EffAcc}")
    fig, ax = plt.subplots()
    ax.plot(post_cut,FOM,'bo',lw=2)
    ax.set_xlabel("Posterior Probability")
    ax.set_ylabel(r"FOM($S/\sqrt{S+B}$)")
    fig.show()

    input()