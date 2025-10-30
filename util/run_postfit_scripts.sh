#!/usr/bin/env bash
# Usage: ./util/run_postfit_scripts.sh <dirname> <corr> <alphan> <data>
# Example: ./util/run_postfit_scripts.sh /data/snoplus/parkerw/antinu/Oct22_corr_fit_noprior/ true false false

# Normalise boolean-like values (accept true/false or 1/0)
normalise_bool() {
  case "$1" in
    1|true|True|TRUE|yes|Yes|YES) echo 1 ;;
    0|false|False|FALSE|no|No|NO) echo 0 ;;
    *) echo "Error: invalid bool '$1' (use 1/0 or true/false)" >&2; exit 1 ;;
  esac
}

set -euo pipefail

# Parse arguments
if [ $# -ne 4 ]; then
  echo "Usage: $0 <dirname> <corr> <alphan> <data>"
  exit 1
fi

dirname="$1"
corr="$2"
alphan="$3"
data="$4"

# Output log file
logfile="${dirname}/run_postfit_scripts.log"
mkdir -p "$dirname"
exec > >(tee -a "$logfile") 2>&1

echo "---------------------------------------"
echo "Starting run_fixed_osc.sh"
echo "Directory: $dirname"
echo "corr=$corr, alphan=$alphan, data=$data"
echo "Log file: $logfile"
echo "---------------------------------------"

corr=$(normalise_bool "$corr")
alphan=$(normalise_bool "$alphan")
data=$(normalise_bool "$data")

# Choose which config to run
if [ "$corr" -eq 1 ] && [ "$alphan" -eq 0 ]; then
    cmd="./bin/makeFixedOscTree cfg/bismsb_ppo_corr/fit_config.ini cfg/bismsb_ppo_corr/oscgrid_config.ini"
elif [ "$corr" -eq 1 ] && [ "$alphan" -eq 1 ]; then
    cmd="./bin/makeFixedOscTree cfg/bismsb_ppo_corr_alphan/fit_config.ini cfg/bismsb_ppo_corr_alphan/oscgrid_config.ini"
elif [ "$corr" -eq 0 ] && [ "$alphan" -eq 0 ]; then
    cmd="./bin/makeFixedOscTree cfg/bismsb_ppo_uncorr/fit_config.ini cfg/bismsb_ppo_uncorr/oscgrid_config.ini"
else
    echo "Unsupported combination: corr=$corr, alphan=$alphan"
    exit 1
fi

echo "Running command: $cmd"
$cmd | tee tmp_output.log

# Extract TTree path
ttree_path=$(grep -oE 'TTree saved to [^ ]+' tmp_output.log | awk '{print $4}')
if [ -z "$ttree_path" ]; then
  echo "ERROR: Could not find TTree path in output."
  exit 1
fi
echo "Found TTree path: $ttree_path"

# Extract fit result directory
outputfitdir=$(grep -oE 'Saved fit result to [^ ]+' tmp_output.log | awk '{print $5}' | sed 's|/fit_result.txt||' | head -n1)
if [ -z "$outputfitdir" ]; then
  echo "ERROR: Could not find output fit directory."
  exit 1
fi
echo "Found output fit directory: $outputfitdir"

echo "---------------------------------------"
echo "Running ROOT plotting macros"
echo "---------------------------------------"

root -l -b -q "plotting/plotFixedOscLLH.C(\"${ttree_path}\")"
root -l -b -q "plotting/plotFixedOscParams.C(\"${ttree_path}\", ${corr}, ${alphan})"

if [ "$data" = "false" ]; then
    for i in 0 1 2; do
        root -l -b -q "plotting/plotFixedOscDist.C(\"${ttree_path}\", ${i})"
    done
elif [ "$data" = "true" ] && [ "$alphan" = "false" ]; then
    pdfpath="/data/snoplus/weiiiiiii/antinuFit/all/prunepdfs"
    for i in 0 1 2; do
        root -l -b -q "plotting/plotFixedOscDist.C(\"${ttree_path}\", ${i}, \"${pdfpath}\")"
    done
elif [ "$data" = "true" ] && [ "$alphan" = "true" ]; then
    pdfpath="/data/snoplus/weiiiiiii/antinuFit/all_withANClassifer/prunepdfs"
    for i in 0 1 2; do
        root -l -b -q "plotting/plotFixedOscDist.C(\"${ttree_path}\", ${i}, \"${pdfpath}\")"
    done
fi

echo "---------------------------------------"
echo "Running Python post-processing"
echo "---------------------------------------"

python util/combine_postfit_results.py "${outputfitdir}/postfit_dists" "${ttree_path}" "${corr}" "${alphan}"
python plotting/plotFixedOscFullLLH.py -i "${ttree_path}"

echo "---------------------------------------"
echo "All steps complete!"
echo "Results summary:"
echo "  TTree path:       ${ttree_path}"
echo "  Output fit dir:   ${outputfitdir}"
echo "  Log file:         ${logfile}"
echo "---------------------------------------"
