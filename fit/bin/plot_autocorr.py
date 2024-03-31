import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str)
parser.add_argument("outfile", type=str)
args = parser.parse_args()

d = np.loadtxt(args.infile)
plt.xlabel("lag")
plt.ylabel("Auto-correlation")
plt.plot(d[:, 0], d[:, 1])
plt.savefig(args.outfile)
