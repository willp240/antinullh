''' Gets plot of discovery potential for various values in terms of CI
by integrating the posterior up to a point where the probability reaches the same value as counts = 0'''

import ROOT
import os
import argparse

ROOT.gROOT.SetStyle("thesisStyle")
ROOT.gROOT.ForceStyle()

parser = argparse.ArgumentParser()
parser.add_argument("maindir")
args = parser.parse_args()

PDF_TYPES = ["2d", "3d_full", "3d_composite_full", "composite_full"]

MASS_DIR_DICT = {"/home/kroupova/bb_analysis/azimov_3yrs_25mev" : 25.0, "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_40mev" : 40.0, "/home/kroupova/bb_analysis/azimov_3yrs_50mev" : 50.0, "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_60mev" : 60.0, "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_65mev" : 65.0, "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_70mev" : 70.0, "/home/kroupova/bb_analysis/azimov_3yrs_75mev" : 75.0, "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_85mev" : 85.0, "/home/kroupova/bb_analysis/azimov_3yrs_90mev" : 90.0, "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_100mev" : 100.0, "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_110mev" : 110.0, "/home/kroupova/bb_analysis/azimov_3yrs_125mev" : 125.0, "/home/kroupova/bb_analysis/azimov_3yrs_150mev" : 150.0, "/data/snoplus2/kroupova/bb_march20/azimov_3yrs_200mev/" : 200.0}


directory = args.maindir
os.chdir(directory)
discovery_mgraph = ROOT.TMultiGraph()
leg = ROOT.TLegend(0.59, 0.21, 0.86, 0.44)

outfile = ROOT.TFile(os.path.join(directory, "discovery_potential_bayes_assym.root"), "recreate")
outfile.cd()

#loop over pdf types
for j in range(0, len(PDF_TYPES)):

    discovery_graph = ROOT.TGraph()
    itype = PDF_TYPES[j]

    #loop over masses
    for idirectory, imass in MASS_DIR_DICT.items():
        idir = os.path.join(idirectory, itype)

        try:
            projection_file = ROOT.TFile.Open(os.path.join(idir, "summed_1dlh_proj/0v.root"))
            print projection_file
            projection_file.cd()
            ihisto=projection_file.Get("")
            ihisto.Scale(1.0/(ihisto.Integral()))

            #fit gaussian to have some idea of the mean - actually not gonna use it really
            ihisto.Fit("gaus", "+", "", 0, 1000)
            ifit = ihisto.GetFunction("gaus")
            imean = ifit.GetParameter(1)
            print "mean: ", imean
            print "itype: ", itype, ", mass: ", imass

            
            #get the height at intercept
            bin1 = ihisto.FindBin(0)
            zero_value = ihisto.GetBinContent(bin1)
            print "zero value: ", zero_value
            

            # now find the last bin after the mean which is still higher then the bin at the intercept
            # and take that as the point where you wanna integrate
            bin_mean = ihisto.FindBin(imean)
            bin2 = bin_mean
            for k in range(bin_mean, ihisto.GetNbinsX()):
                if (ihisto.GetBinContent(k)<zero_value):
                    break
                bin2 = k

            # get the integral between the two bins
            area_histo = ihisto.Integral(bin1,bin2)
            
            ipoint_y = area_histo
            ipoint_x = imass

            if (imean<0):
                ipoint_y = 0
            
            discovery_graph.SetPoint(discovery_graph.GetN(), ipoint_x, ipoint_y)
            projection_file.Close()

        except:
            print "Couldn't add point for ", itype

    outfile.cd()
    discovery_graph.Write(itype)
    discovery_mgraph.Add(discovery_graph, "PL")
    leg.AddEntry(discovery_graph, itype, "pel")

outfile.cd()
c = ROOT.TCanvas()
c.cd()

discovery_mgraph.Draw("A")
discovery_mgraph.GetXaxis().SetTitle("meV")
discovery_mgraph.GetYaxis().SetTitle("Confidence Level")

leg.SetTextSize(0.045);
leg.Draw()
        
#save to
discovery_mgraph.Write("multi")
c.Write("canvas_multi")
outfile.Close()


