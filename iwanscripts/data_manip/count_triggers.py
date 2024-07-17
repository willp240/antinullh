import ROOT

def count_ev_indices(t_name, *filenames):
    ch = ROOT.TChain(t_name)
    for fn in filenames:
        ch.Add(fn)
    
    triggers = [0] * 10
    for x in ch:
        if (x.evIndex) < 9:
            triggers[x.evIndex + 1] += 1
        else:
            triggers[9] += 1
    return triggers, ch.GetEntries()
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile_names", type=str, nargs="+")
    parser.add_argument("--tree_name", type=str, default = "output")
    args = parser.parse_args()

    event_indices, entries = count_ev_indices(args.tree_name, *args.infile_names)
    for i, count in enumerate(event_indices):
        print i-1, count

print "Number of non-retriggers = ", event_indices[0] + event_indices[1]
print "Number of entries = ", entries

if entries:
    print "Correction = ", (event_indices[0] + event_indices[1])*1./entries

