#!/usr/bin/env bash
# Usage: ./util/run_postfit_scripts.sh <cfgdirname> <corr> <alphan> <data>
# Example: ./util/run_postfit_scripts.sh /home/parkerw/Software/antinu_llh2/cfg/reacAnalysis2025/asmv_correlated_noprior_noalphan true false false

# Normalise boolean-like values (accept true/false or 1/0)
normalise_bool() {
  case "$1" in
    1|true|True|TRUE|yes|Yes|YES) echo 1 ;;
    0|false|False|FALSE|no|No|NO) echo 0 ;;
    *) echo "Error: invalid bool '$1' (use 1/0 or true/false)" >&2; exit 1 ;;
  esac
}

# Get output directory from a fit config file
get_output_directory() {
    local config_file="$1"

    if [[ ! -f "$config_file" ]]; then
        echo "Error: file '$config_file' not found." >&2
        return 1
    fi

    # Extract the value
    local outdir
    outdir=$(grep -E '^output_directory[[:space:]]*=' "$config_file" \
             | head -n 1 \
             | sed -E 's/^output_directory[[:space:]]*=[[:space:]]*//')

    echo "$outdir"
}

set -euo pipefail

# Parse arguments
if [ $# -ne 4 ]; then
  echo "Usage: $0 <config_dir> <corr> <alphan> <data>"
  exit 1
fi

cfg_dir="$1"
corr="$2"
alphan="$3"
data="$4"

fit_config="$cfg_dir/fit_config.ini"
osc_config="$cfg_dir/oscgrid_config.ini"
dirname=$(get_output_directory "$fit_config")

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

cmd="./bin/makeFixedOscTree ${fit_config} ${osc_config}"

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

if [ "$data" -eq 0 ]; then
    for i in 0 1 2; do
        root -l -b -q "plotting/plotFixedOscDist.C(\"${ttree_path}\", ${i})"
    done
elif [ "$data" -eq 1 ] && [ "$alphan" -eq 0 ]; then
    pdfpath="/data/snoplus/weiiiiiii/antinuFit/all/prunepdfs"
    for i in 0 1 2; do
        root -l -b -q "plotting/plotFixedOscDist.C(\"${ttree_path}\", ${i}, \"${pdfpath}\")"
    done
elif [ "$data" -eq 1 ] && [ "$alphan" -eq 1 ]; then
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
