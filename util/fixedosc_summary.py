"""Summarise the best fit from a fixed-oscillation grid"""

import os
import re
import json
import sys
import subprocess
import math
import csv

# --------------------------------------------------------------
# Run get_integrals.py and plotFixedOscParams
# --------------------------------------------------------------
def run_get_integrals(postfit_dir):
    cmd = ["python", "util/get_integrals.py", postfit_dir]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("❌ Error running get_integrals.py:\n", result.stderr)
        sys.exit(1)
    lines = result.stdout.strip().splitlines()
    pattern = re.compile(r"^(\S+)\s+integral\([^)]*\)\s*=\s*([0-9.+-eE]+)")
    integrals = {}
    for line in lines:
        m = pattern.search(line)
        if m:
            fname = m.group(1)
            val = float(m.group(2))
            integrals[os.path.splitext(fname)[0]] = val
    return integrals


def run_plot_fixed_osc_params(fit_result_tree, corr, alphan):
    cmd = ["root", "-l", "-b", "-q", f'plotting/plotFixedOscParams.C("{fit_result_tree}",{int(corr)},{int(alphan)})']
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("❌ Error running plotFixedOscParams:\n", result.stderr)
        sys.exit(1)
    return result.stdout.strip().splitlines()


# --------------------------------------------------------------
# Parse plotFixedOscParams output
# --------------------------------------------------------------
def parse_plot_params(lines):
    """
    Returns:
      params: {param: (val, err)}
      corr: float or dict {"ppo": x, "bismsb": y}
    """
    params = {}
    geo_corr = {}
    current = None
    for line in lines:
        # Geo correlations
        if "Geo correlation" in line:
            parts = line.split(":")
            if len(parts) == 2:
                label = parts[0].strip().split()[0]  # PPO or bisMSB or nothing
                val = float(parts[1].strip())
                if label in ["PPO", "bisMSB"]:
                    geo_corr[label.lower()] = val
                else:
                    geo_corr["shared"] = val

        # Parameters
        if line.startswith("Par:"):
            current = line.split("Par:")[1].strip()
        elif current and line.startswith("Fit:"):
            parts = line.split()
            if len(parts) >= 3:
                value = float(parts[1])
                uncert = float(parts[2])
                params[current] = (value, uncert)
            current = None

    # Determine style
    if any("_norm" in k for k in params):
        mode = "norm"
    else:
        mode = "direct"

    # Default correlations if missing
    if not geo_corr:
        geo_corr = {"shared": 0.0}

    return params, geo_corr, mode


# --------------------------------------------------------------
# Combine results (norm-style)
# --------------------------------------------------------------
from minuit_results import (
    normalisation_rate_error, combine_results_normstyle, combine_results_directstyle,
    add_geo_ratios, write_csv,
)


# --------------------------------------------------------------
# Main
# --------------------------------------------------------------
def main():
    if len(sys.argv) != 5:
        print("Usage: python util/fixedosc_summary.py /path/to/postfit_dists /path/to/fit_result_tree.root correlated-fit-bool (alpha,n)-classifier-bool")
        sys.exit(1)

    postfit_dir, fit_result_tree, corr, alphan = sys.argv[1:5]
    corr = corr.lower() in ("1", "true", "yes")
    alphan = alphan.lower() in ("1", "true", "yes")
    print("Running get_integrals.py...")
    integrals = run_get_integrals(postfit_dir)
    print("Running plotFixedOscParams...")
    lines = run_plot_fixed_osc_params(fit_result_tree, corr, alphan)
    fit_params, geo_corr, mode = parse_plot_params(lines)
    if not fit_params:
        raise RuntimeError("No fitted parameters were recovered from plotFixedOscParams; check the saved covariance errors")
    print(f"Detected mode: {mode}-style")
    print(f"Geo correlation(s): {geo_corr}")

    print("Combining results...")
    if mode == "norm":
        central, uncert = combine_results_normstyle(integrals, fit_params, geo_corr, alphan)
    else:
        central, uncert = combine_results_directstyle(integrals, fit_params, geo_corr, alphan)

    outdir = os.path.dirname(postfit_dir.rstrip("/"))
    with open(os.path.join(outdir, "postfit_central.json"), "w") as f:
        json.dump(central, f, indent=2)
    with open(os.path.join(outdir, "postfit_uncert.json"), "w") as f:
        json.dump(uncert, f, indent=2)

    write_csv(outdir, central, uncert, mode, alphan)
    print(f"\n✅ All results written under {outdir}")

if __name__ == "__main__":
    main()
