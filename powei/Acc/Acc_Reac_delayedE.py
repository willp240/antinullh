# draw delayedE for both reactor IBD and accidental
import ROOT
import sys
sys.path.append('/home/huangp/BiPo/Analysis')
from PyROOT_Style import set_style
from PyROOT_Style import Line_Draw_Setting
from PyROOT_Style import TLegend_Setting
import glob

set_style()
antinueventtype = ["Geoibd_U","Geoibd_Th","Tl210","Alphan_Lab_13c","Alphan_Lab_Avout_Av_18o","Alphan_Lab_Avin_Av_13c","Alphan_Lab_Avout_Av_13c","Alphan_Lab_Avin_Av_18o","Reactoribd","Bipo214","Bipo212" ]
filepath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_sim"
datafilepath = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_Acc/accidental/"

def loadfile(treename,fileloc,evttype):
    for f in glob.glob(f"{fileloc}scaled_oscillated*.root"):
        #print(f)
        if int(f[-24:-18]) > 300600:
            treename.Add(f)
    return
def loaddatafile(treename,fileloc):
    for f in glob.glob(f"{fileloc}*.root"):
        #print(f)
        treename.Add(f)
    return


    




if __name__ == "__main__":

    
    # new ttree
    IBDDelayT = ROOT.TChain("DelayT")
    for itype in antinueventtype:
        fileloc = f"{filepath}/{itype}/" 
        print(fileloc)
    
        if itype == "Reactoribd":#"scaled_oscillated_Reactoribd":
            loadfile(IBDDelayT,fileloc,itype)
    # load data
    DataDelayT = ROOT.TChain("AccT")
    loaddatafile(DataDelayT,datafilepath)
        
    # new histogram
    simbins = 42; databins = 42
    IBDDelay_energyh1 = ROOT.TH1D("IBDDelay_energyh1","IBD_delayedEcorr",simbins,1.85,2.5)
    Acc_Delayenergyh1  = ROOT.TH1D("Acc_Delayenergyh1","Acc_delayedEcorr",databins,1.85,2.5)
    # Project Tree to histo
    var = "delayedEcorr"
    IBDDelayT.Project("IBDDelay_energyh1","delayedEcorr","1.85<delayedEcorr && delayedEcorr<2.5")
    DataDelayT.Project("Acc_Delayenergyh1","delayedEcorr","1.85<delayedEcorr && delayedEcorr<2.5")


    # Print entries information
    Num_ReactorIBD = IBDDelay_energyh1.GetEntries()
    Num_Data       = Acc_Delayenergyh1.GetEntries()
    print("*************** Info about each backgrounds entries *********")
    print(f"*************** ReactorIBD: {Num_ReactorIBD} *********")
    print(f"*************** Data:       {Num_Data} *********")

    # Scale the distribution to the expected yields from Tiny's thesis(After Cuts)
    #Expected_Events = {"ReactorIBD": 27.9}
    IBDDelay_energyh1.Scale(1./IBDDelay_energyh1.Integral(), "nosw2")
    Acc_Delayenergyh1.Scale(1./Acc_Delayenergyh1.Integral(), "nosw2")

    # Canvas 1: Plot every backgrounds/signal distribution
    c1 = ROOT.TCanvas("c1","",800,500)
    # Set Fill Color
    IBDDelay_energyh1.SetLineWidth(2)
    Acc_Delayenergyh1.SetLineWidth(2)
    IBDDelay_energyh1.SetLineColor(1)
    Acc_Delayenergyh1.SetLineColor(2);  
    #IBDDelay_energyh1.SetMaximum(18.5)
    
    IBDDelay_energyh1.GetXaxis().SetTitle("Delayed Reconstructed Energy [MeV]")
    IBDDelay_energyh1.GetYaxis().SetTitle("Asimov Events/0.01 [MeV]")
    IBDDelay_energyh1.Draw()
    Acc_Delayenergyh1.Draw("same")
 
    legend0 = ROOT.TLegend(0.65,0.65,0.85,0.85)

    TLegend_Setting(legend0,[IBDDelay_energyh1,Acc_Delayenergyh1],["Reactor antinu","Data"],"L")
    legend0.Draw()
    

    c1.SaveAs("/home/huangp/AntiNu/Plots/Acc_Reac_delayE.png")
    c1.SaveAs("/home/huangp/AntiNu/Plots/Acc_Reac_delayE.pdf")

    output_root_address = "/home/huangp/AntiNu/Analysis/delayedPDFs.root"
    f_out = ROOT.TFile( output_root_address,"RECREATE")
    IBDDelay_energyh1.Write()
    
    Acc_Delayenergyh1.Write()
    f_out.Close()

    input()