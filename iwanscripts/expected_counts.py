import re 
import ConfigParser
import ROOT 
import os
import errno

# this is the general form of an entry line, match these to get efficiencies
rex = ".*\s+(\d*\.?\d*)\s+(\d*\.?\d*e?-?\d+)\s+(\d*\.?\d*)\s+(\d*\.?\d*e?-?\d+)\s*"

def count_events(filenames, tree_name):
    ch = ROOT.TChain(tree_name)
    ch.Add(filenames)
    return ch.GetEntries()

def get_n_gen_correction(event_config, name, tree_name):
    cparser = ConfigParser.ConfigParser()
    cparser.read(event_config)

    ngen = int(cparser.get(name, "n_generated"))

    try:
        filenames = os.path.join(cparser.get("summary", "split_ntup_dir_fake")  , name + ".root")
    except ConfigParser.NoOptionError:
        print "Failed to read directory from config file!"
        return 1
    
    
    if ngen == 0:
        return 1

    count = count_events(filenames, tree_name)
    corr = float(count)/ngen
    return corr


def expected_rate(event_config, name):
    cparser = ConfigParser.ConfigParser()
    cparser.read(event_config)
    return float(cparser.get(name, "rate"))

def get_pdf(pdf_dir, name, counts):
    f = ROOT.TFile(os.path.join(pdf_dir,  name + ".root"))
    h = f.Get(f.GetListOfKeys()[0].GetName())
    h.SetDirectory(0)
    if h.Integral():
        h.Scale(counts/h.Integral())
    return h

def create_dir(dirname):
    try:
        os.mkdir(dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

def read_file(event_config, livetime, pdf_dir, tree_name, data_fraction, outdir):
    # make the output directories
    create_dir(outdir)
    create_dir(os.path.join(outdir, "scaled_dists"))

    counts = {}
    efficiencies = {}
    rates = {}
    pdfs = {}

    # load up the config file
    cparser = ConfigParser.ConfigParser()                                                                            
    cparser.read(event_config)


    # mask in the correct contributers
    active = [x.strip() for x in cparser.get("summary", "active").split(",")]
    inactive = [x.rstrip() for x in cparser.get("summary", "inactive").split(",")]

    for name in cparser.sections():
        if name == "summary":
            continue
 
        if not "all" in active and not name in active:
            print "skipping " + name
            continue

        if name in inactive:
            print "skipping " + name
            continue

        rate = expected_rate(event_config, name)

        try:
            with open(os.path.join(pdf_dir, name + ".txt")) as f:
                # the final matching line is the final efficiency
                last_matching_line = [line for line in f.read().splitlines() if re.match(rex, line) is not None][-1]
                
                match = re.match(rex, last_matching_line)
                corr = get_n_gen_correction(event_config, name, tree_name) / data_fraction
                eff = float(match.group(4))/100 * corr 
                counts[name] = eff * rate * livetime
                efficiencies[name] = eff
                rates[name] = rate

                try:
                    pdfs[name] = get_pdf(pdf_dir, name, counts[name])
                except:
                    pass
                    
        except Exception as e:
            print "Exception encountered for event type {0} : {1}".format(name, e)
            pass

    for x in pdfs.keys():
        pdfs[x].SaveAs(os.path.join(outdir, "scaled_dists",  x + ".root"))


    with open(os.path.join(outdir, "expected_counts.dat"), "w") as logfile:
        row_format = "{0:>20}   {1:>20}   {2:>20}   {3:>20}    {4:>20}\n"
        logfile.write(row_format.format("Name", "Total", "Eff", "In box/yr", "In box total"))
        logfile.write("\n")
        for x, y in sorted(counts.items(), key = lambda x: x[1], reverse = True):
            logfile.write( row_format.format(x, rates[x], efficiencies[x], (counts[x]/livetime), counts[x]))
            
        logfile.write("\n\nExpected counts in fit region ={0}".format( sum(counts.values())))
        logfile.write("\n\n number of types with expected rate > 1/year  = {0}".format(len([x for x in counts.values() if x > (1*livetime)])))


if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser()
    parser.add_argument("event_config_file", type=str)
    parser.add_argument("livetime", type=float)
    parser.add_argument("pdf_dir", type=str)
    parser.add_argument("result_dir", type=str)
    parser.add_argument("--data_fraction", type=float, default=1.)
    parser.add_argument("--tree_name", type=str, default="pruned")

    args = parser.parse_args()
    read_file(args.event_config_file,
              args.livetime,
              args.pdf_dir,
              args.tree_name,
              args.data_fraction,
              args.result_dir)
