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

PDF_TYPES = ["2d", "3d_full", "3d_composite_full", "composite_full"]

MAIN_DIR = args.maindir
#MAIN_DIR = "/home/kroupova/bb_analysis/azimov_3yrs_nosig/composite_full"

NO_SIG_DIR = "/home/kroupova/bb_analysis/azimov_3yrs_nosig"
MASSES = [25.0, 40.0, 50.0, 60.0, 65.0, 70.0, 75.0, 85.0, 90.0, 100.0, 110.0, 125.0, 150.0, 200.0]
DIRECTORIES = ["/home/kroupova/bb_analysis/azimov_3yrs_25mev", "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_40mev", "/home/kroupova/bb_analysis/azimov_3yrs_50mev", "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_60mev", "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_65mev", "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_70mev", "/home/kroupova/bb_analysis/azimov_3yrs_75mev", "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_85mev", "/home/kroupova/bb_analysis/azimov_3yrs_90mev", "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_100mev", "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_110mev", "/home/kroupova/bb_analysis/azimov_3yrs_125mev", "/home/kroupova/bb_analysis/azimov_3yrs_150mev", "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_200mev/"]



directory = MAIN_DIR
os.chdir(MAIN_DIR)
discovery_mgraph = ROOT.TMultiGraph()
leg = ROOT.TLegend(0.59, 0.21, 0.86, 0.44)

outfile = ROOT.TFile(os.path.join(directory, "discovery_potential_bayes.root"), "recreate")
outfile.cd()

#loop over pdf types
for j in range(0, len(PDF_TYPES)):


    discovery_graph = ROOT.TGraph()
    itype = PDF_TYPES[j]
    #add the histos from the chains
    for i in range(0,len(DIRECTORIES)):
        idir = os.path.join(DIRECTORIES[i], itype)

        try:
            projection_file = ROOT.TFile.Open(os.path.join(idir, "summed_0v_1dlh.root"))
            print projection_file
            projection_file.cd()
            ihisto=projection_file.Get("0v_1dlh")
            ihisto.Scale(1.0/(ihisto.Integral()))
            #ihisto.SetDirectory(0)

            ihisto.Fit("gaus", "+", "", 0, 1000)
            ifit = ihisto.GetFunction("gaus")
            imean = ifit.GetParameter(1)
            print "mean: ", imean
            print "itype: ", itype, ", mass: ", MASSES[i]

            area = ifit.Integral(0,2*imean)
            print "area fit: ", area
            
            bin1 = ihisto.FindBin(0)
            bin2 = ihisto.FindBin(imean*2)
            
            area_histo = ihisto.Integral(bin1,bin2)
            print "area histo: ", area_histo
            
            ipoint_y = area_histo
            ipoint_x = MASSES[i]

            if (imean<0):
                ipoint_y = 0
            
            discovery_graph.SetPoint(discovery_graph.GetN(), ipoint_x, ipoint_y)
            projection_file.Close()

        except:
            ipoint_x = MASSES[i]
            #discovery_graph.SetPoint(discovery_graph.GetN(), ipoint_x, 0.0)
            print "Couldn't add point for ", itype
    #discovery_graph.SetLineColor(ROOT.GetColor(COLOURS[j]))
    outfile.cd()
    discovery_graph.Write(itype)
    discovery_mgraph.Add(discovery_graph, "PL")
    leg.AddEntry(discovery_graph, itype, "pel")

outfile.cd()
c = ROOT.TCanvas()
c.cd()

discovery_mgraph.Draw("A")
discovery_mgraph.GetXaxis().SetTitle("mev")
discovery_mgraph.GetYaxis().SetTitle("Confidence Level")

leg.SetTextSize(0.045);
leg.Draw()
        
#save to

discovery_mgraph.Write("multi")
c.Write("canvas_multi")
outfile.Close()


