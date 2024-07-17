import ROOT
#from ROOT import gROOT
import numpy as np
import os
import argparse

ROOT.gROOT.SetStyle("thesisStyle")
ROOT.gROOT.ForceStyle()

parser = argparse.ArgumentParser()
parser.add_argument("maindir")
args = parser.parse_args()

PDF_TYPES = ["2d", "3d_composite_full", "3d_full", "composite_full"]#, "2d"]#, "composite_full", "3d_full", "3d_composite_full"]

DATASET_DIRECTORIES = ["/data/snoplus2/kroupova/bb_march20/150mev/independant_datasets", "/data/snoplus2/kroupova/bb_march20/150mev/datasets_same_repeating_externals_1", "/data/snoplus2/kroupova/bb_march20/150mev/datasets_same_repeating_externals_2", "/data/snoplus2/kroupova/bb_march20/150mev/datasets_same_repeating_externals_3"] 

SIG_EFF = [0.742, 0.739, 0.739, 0.739]

directory = args.maindir
os.chdir(MAIN_DIR)

outfile = ROOT.TFile(os.path.join(directory, "bias.root"), "recreate")
outfile.cd()

#loop over pdf types
for j in range(0, len(PDF_TYPES)):

    itype = PDF_TYPES[j]
    #Graph
    count_histo = ROOT.TH1D("", "", 1000, 0., 1000)
    true_count_histo = ROOT.TH1D("", "", 1000, 0., 1000)
    bias_histo = ROOT.TH1D("", "", 1000, -1., 1)
    pull_histo = ROOT.TH1D("", "", 1000, -10., 10)
    
   

    #add the histos from the chains
    for i in range(0,len(DATASET_DIRECTORIES)):
        
        for k in range(0,4):
            idir = os.path.join(DATASET_DIRECTORIES[i],"dataset_{0}".format(k) )
            kdir = os.path.join(idir,itype )

            txtfile = open(os.path.join(DATASET_DIRECTORIES[i], "fake_data_lt_3__{0}.txt".format(k)))
            line = txtfile.readline()
            print line
            true_count = float(line.split()[-1])
            print true_count
            txtfile.close()
            true_count_histo.Fill(float(true_count))
            
            try: 
                projection_file = ROOT.TFile.Open(os.path.join(kdir, "summed_1dlh_proj/0v.root"))
                projection_file.cd()
                ihisto=projection_file.Get("")
                ihisto.Scale(1.0/(ihisto.Integral()))
                #ihisto.SetDirectory(0)

                mean_histo = ihisto.GetMean()
                rms_histo = ihisto.GetRMS()
                ihisto.Fit("gaus", "+", "", 0, 1000)
                ifit = ihisto.GetFunction("gaus")
                imean = ifit.GetParameter(1)
                ibin_mean = ihisto.FindBin(imean)
                print "mean: ", imean
                print "rms: ", rms_histo
                print "man histo: ", mean_histo
                print "bin: ", ibin_mean
                print "i: ", i, ", pdf: ", itype
                


                count_histo.Fill(imean)

                ibias = (imean-(SIG_EFF[j]*true_count))/(SIG_EFF[j]*true_count)
                ipull = (imean-(SIG_EFF[j]*true_count))/(rms_histo)
                print "bias: ", ibias
                bias_histo.Fill(ibias)
                pull_histo.Fill(ipull)
                
                projection_file.Close()
            
            except:
                print "Couldn't add point for pdf type: ", itype

    outfile.cd()
    count_histo.Write("fitted_counts_"+itype)
    
    true_count_histo.Write("true_counts_"+itype)
    bias_histo.Write("bias_"+itype)
    pull_histo.Write("pull_"+itype)
    count_histo.Draw()
    

outfile.Close()


