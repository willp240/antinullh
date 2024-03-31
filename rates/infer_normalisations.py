from glob import glob
import os
import re

def needs_rescale(filename, cut_names):
    ''' Check if this macro includes a cut that should prompt us to rescale
    :param str filename: path to macro
    :param list(str) cut_names:
    :returns bool: True if one of the cut names is in the file
    '''
    with open(filename) as f:
        return any(x in f.read() for x in cut_names)

def grab_rate(filename):
    ''' Parse the macro in filename and extract the rate
    :param str filename:
    :returns float:
    '''
    with open(filename) as f:
        for x in f.readlines():
            if "/generator/rate/set" in x:
                return float(x.split()[1])

def print_dict(rates):
    ''' 
    '''
    for k,v in rates.iteritems():
        print k, v
    
def save_dict(rates, outname):
    ''' Save a dictionary to a text file in the form 
    val0a--val0b
    etc.
    :param dict(string, float) rates: 
    :param str outname: to save the data to
    '''
    if not len(rates):
        return

    with open(outname, "w") as f:
        for pair in rates.iteritems():
            f.write("{0}\t{1}\n".format(*pair))

def count_inside(wildcard_str):
    ''' Count the number of files corresponding in this path (supports *)
    '''
    return len(glob(os.path.join(wildcard_str, "*.root")))
    

def list_top_dir(dirname):
    ''' Get all the of MC file counts
    :param str dirname: this is the top directory of the MC
    :returns dict(str, int): keyed by the subdirectory name
    
    dirname --- directory for bg1
            --- directory for bg2
            etc.
    
    '''
    file_counts = {}
    for x in os.listdir(dirname):
        x = os.path.join(dirname, x)
        if os.path.isdir(x):
            file_counts[x] = count_inside(x)
    return file_counts

def find_dir(top_dir, mac_name):
    ''' Find which MC sub dir corresponds to this name 
    :param str top_dir: this is the top directory of the MC
    :param mac_name: 
    :returns: the path to the subdir relative to topdir
    '''
    for x in os.listdir(top_dir):
        x = os.path.join(top_dir, x)
        base_name = os.path.splitext(os.path.basename(mac_name))[0]

        try:            
            if os.path.isdir(x) and re.search("TeLoaded{0}_r\d+_s\d+_p\d+.ntuple.root".format(base_name).lower(), glob(os.path.join(x, "*.root"))[0].lower()) is not None:
                return x
        except IndexError as e:
#            print "Warning: no ROOT files in {0}".format(x)
            pass

def map_macros_to_dirs(top_dir, mac_names):
    ''' Find the mc subdirectory that corresponds to a set of macros
    :param str top_dir: this is the top directory of the MC
    :param list(str) mac_names: this is the list of macros we want to match up
    :returns dict(str, str): the paths to the subdir relative to topdir keyed by macro name
    ''' 
    mapping = {}
    for x in mac_names:
        mapping[x] = find_dir(top_dir, x)
    return mapping

def dump_to_json(jsfile, result):
    ''' You guessed it
    '''
    with open(jsfile, "w") as f:
        f.write(json.dumps(dict( (x, (y , z)) for (x, y, z) in result )))
            
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("mc_top_dir", type=str)
    parser.add_argument("macro_files", type=str, help="if you're using a wildcard be sure to wrap this argument in double quotes *needs to be abspath*")
    parser.add_argument("run_duration", type=float, help="How long were the simulations run for? Historically this has been one year")
    args = parser.parse_args()

    cuts_requiring_rescale = ["scintEdepCut", ""]

    mac_mapping = map_macros_to_dirs(args.mc_top_dir, glob(os.path.abspath(args.macro_files)))

    file_count_mapping = {}
    for k,v in mac_mapping.iteritems():
        if v is None:
            continue
        file_count_mapping[k] = count_inside(v)

    print "Assuming a simulation time of {0}yr ..\n\n".format(args.run_duration)
    row_format ="{0:>30}   {1:>30}   {2:>30}   {3:>30}"
    print row_format.format("Macro", "Macro rate", "File Count", "Estimated simulated event count")
    warnings = []
    for mac, dir_name in mac_mapping.iteritems():
        if dir_name is None:
            continue
        if needs_rescale(mac, cuts_requiring_rescale):
            file_count = file_count_mapping[mac]
            macro_rate = grab_rate(mac)
            if macro_rate is None:
                warnings.append("Can't interpret a rate for {0}! skipping".format(os.path.basename(mac)))
                continue
            print row_format.format(os.path.basename(mac), macro_rate, file_count, int(file_count * macro_rate * (365 * 24 * 3600 * args.run_duration)))


    print "\n\n\n Warnings:\n\t" + "\n\t".join(warnings)
