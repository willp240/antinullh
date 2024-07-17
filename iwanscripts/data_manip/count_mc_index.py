import ROOT

def count_ev_indices(t_name, *filenames):

    total_mc = 0
    
    for fn in filenames:
        ch = ROOT.TChain(t_name)
        ch.Add(fn)

        mc_max = 0
        for x in ch:
            #print "file: ", fn
            #if x.mcIndex <=0:
                #print "x.mcIndex: ", x.mcIndex
                #print "file: ", fn
                #print "x.mcIndex: ", x.mcIndex

            if (x.mcIndex) > mc_max:
                mc_max = x.mcIndex

        total_mc += mc_max+1
    return total_mc
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile_names", type=str, nargs="+")
    parser.add_argument("--tree_name", type=str, default = "output")
    args = parser.parse_args()

    entries = count_ev_indices(args.tree_name, *args.infile_names)
    
print "Number of mc entries = ", entries

