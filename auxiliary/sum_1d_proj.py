#################################################################
# Sums the 1dlhproj histograms across all the chains            #
# Assumes the chains were run in subdirectories called fit_*    #
# as output by submit_fits.py                                   #
#                                                               #
# input: main directory with all the fit chains                 #
#################################################################

import ROOT
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("maindir")
args = parser.parse_args()

directory = args.maindir

# specific for your batch system, assumes higher numbers after this identifier are the most recent ones
OUTPUT_FILE_IDENTIFIER = ".sh.o" 
# ideal acceptance is ~0.6-0.7, but depends on your tuning
ACCEPTANCE_MIN = 0.3
ACCEPTANCE_MAX = 0.8


def is_acceptance_ok(subdirectory):
    subdirlist = os.listdir(subdirectory)
    outfiles = [x for x in subdirlist if OUTPUT_FILE_IDENTIFIER in x]
    latest_file = outfiles[0]
    for ifile in outfiles:
        if (float(ifile.split(OUTPUT_FILE_IDENTIFIER)[1]) > float(latest_file.split(OUTPUT_FILE_IDENTIFIER)[1])):
            latest_file = ifile

    open_file = open(os.path.join(subdirectory, latest_file), "r")
    lines = open_file.readlines()

    for line in lines:
        if line.startswith("MCMC:: acceptance"):
            acceptance = float(line.split("=")[1])
            if (acceptance>ACCEPTANCE_MIN) and (acceptance<ACCEPTANCE_MAX):
                return True
            else:
                return False


os.chdir(directory)

#get the list of all of the fit subdirectories
fitdir_list = [s for s in os.listdir(directory) if "fit_" in s]

#take the projection names and histogram binnings from the first subdirectory
example_proj_dir = os.path.join(fitdir_list[0], "1dlhproj")
proj_list =os.listdir(example_proj_dir)
bckg_list = []
histo_list = []
for i in range(0, len(proj_list)):
    bckg_list.append(((proj_list[i]).split("/"))[-1]) # just the name
    projection_file = ROOT.TFile.Open(os.path.join(example_proj_dir, proj_list[i]))
    projection_file.cd()
    histo = projection_file.Get("").Clone()
    histo.SetDirectory(0)
    histo.Reset() #just keep the axises etc, don't wanna double count the first chain
    histo_list.append(histo)
    del histo
    projection_file.Close()

print "Summing projections for the following distributions: ", bckg_list


#loop over the various background projections and sum
for j in range (0, len(bckg_list)):
    ibck = bckg_list[j]
    histo_summed = histo_list[j]

    #add the histos from the chains
    for subdirectory in fitdir_list: 
        
        #only include projection from chains with acceptable acceptance (haha)
        if (is_acceptance_ok(subdirectory) == False):
            print "Ignoring chain with acceptance out of bounds!"
            continue
        
        projection_file = ROOT.TFile.Open(os.path.join(subdirectory, "1dlhproj/{0}".format(ibck)))
        projection_file.cd()
        histo = projection_file.Get("")
        histo_summed.Add(histo)
        projection_file.Close()
        
    histo_summed.SetName("")

    #save the summed projections to its own directory
    summed_dir = os.path.join(directory, "summed_1dlh_proj")
    if not os.path.isdir(summed_dir):                                                      
        os.mkdir(summed_dir)
    outfile = ROOT.TFile(os.path.join(summed_dir, ibck), "recreate")
    outfile.cd()
    histo_summed.Write()
    outfile.Close()

