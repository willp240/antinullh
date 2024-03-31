''' Count the files and events for the production data
'''
import glob
import sys
import ROOT
import os
import json

def count_files(filestr):
    '''Count the files matching the string given (wildcards included)
    :param str filestr: pattern to match
    :returns int: the count
    '''
    return len(glob.glob(filestr))

def count_events(filestr, tree):
    ''' Count the events inside all the root files that match this pattern
    in this tree
    :param str filestr: pattern to match
    :param str tree: the name of the tree in the ROOT directory
    :returns int: the count
    '''
    ch = ROOT.TChain(tree)
    ch.Add(filestr)
    return ch.GetEntries()

def examine_directory(dirname, tree_name):
    ''' Walk down the directory 
    :param dirname:
    :param str tree: key of the tree in ROOT directory
    :returns list( ("name", int , int) ): dirname, nfiles, nevents
    '''
    returnlist = []
    for root, dirs, files in os.walk(dirname):
        root_files = os.path.join(root, "*.root")
        print "Counting ", root_files
        sys.stdout.flush()
        files = count_files(root_files)
        evs = 0
        if files:
            evs = count_events(root_files, tree_name)
        returnlist.append((os.path.basename(root), files, evs))
    return returnlist

def examine_files(path, tree):
    ''' Look at files matching this path
    :param str path:
    :param str tree: key of the tree in ROOT directory
    :returns list( ("name", int , int) ): dirname, nfiles, nevents
    '''
    print "Counting", path, "..."
    sys.stdout.flush()
    return [(os.path.basename(path), count_files(path), count_events(path, tree))]

def format_result(result):
    '''Take the result and make it a nice string
    :param list( ("name", int, int) )
    :returns str:
    '''
    return "\n".join("{0: <40} | {1: <20} | {2: <20}".format(*x) for x in result)

def examine_path(path, tree):
    ''' If this is a directory look inside, if its a file string just match
    that string
    :param str path: to the files or directory
    :returns list( ("name", int, int) ):
    '''
    if os.path.isdir(path):
        return examine_directory(path, tree)
    else:
        return examine_files(path, tree)
    
def dump_to_json(jsfile, result):
    ''' You guessed it
    '''
    with open(jsfile, "w") as f:
        f.write(json.dumps(dict( (x, (y , z)) for (x, y, z) in result )))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", type = str, nargs="+")
    parser.add_argument("-tree", type = str, default = "output")
    parser.add_argument("-to_json", type = str, dest = "jsfile", default = None)
    args = parser.parse_args()

    result = []
    for path in args.paths:
        print "Top level path = ", path, "\n\n"
        result += examine_path(path, args.tree)

    if args.jsfile is not None:
        dump_to_json(args.jsfile, result)

    print format_result(result)
