#!/usr/bin/env python3
import sys
import os
import re
import ROOT

def main():
    if len(sys.argv) != 2:
        print("Usage: python get_integrals.py <directory>")
        sys.exit(1)

    indir = sys.argv[1]
    if not os.path.isdir(indir):
        print(f"Error: {indir} is not a directory")
        sys.exit(1)

    # Loop over ROOT files in the directory
    for fname in sorted(os.listdir(indir)):
        if not fname.endswith(".root"):
            continue

        fullpath = os.path.join(indir, fname)
        f = ROOT.TFile.Open(fullpath)
        if not f or f.IsZombie():
            print(f"{fname}: could not open file")
            continue

        # Determine histogram name by removing _2p2ppo.root or _bismsb.root
        hist_name = re.sub(r"(_2p2ppo|_bismsb)\.root$", "", fname)

        hist = f.Get(hist_name)
        if not hist or not isinstance(hist, ROOT.TH1):
            print(f"{fname}: histogram '{hist_name}' not found")
            f.Close()
            continue

        # Bin range 8–150 (inclusive)
        bin_min = 8
        bin_max = 150
        integral = hist.Integral(bin_min, bin_max)
        print(hist.Integral())
        print(f"{fname:35s}  {hist_name:30s}  integral(8–150) = {integral:.6g}")

        f.Close()

if __name__ == "__main__":
    main()
