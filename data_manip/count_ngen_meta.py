'''
Gives counts based on meta information. Important especially for modules with scintEDepCut. Returned is the sum of meta counts across of all given root files. This should be used as an input in event config as n_generated.
'''


from rat import ROOT

def count_ev_indices(*filenames):

    total_ngen = 0

    
    for fn in filenames:

        f = ROOT.TFile(fn)
        meta = f.Get("meta")
        total_ngen += (meta.GetEventsGeneratedCounts().at(meta.GetCurrentPass()))
        
    print "Len: ", len(filenames)    
    return total_ngen
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile_names", type=str, nargs="+")
    args = parser.parse_args()

    entries = count_ev_indices(*args.infile_names)
    
print "Number of generated events from metadata = ", entries

