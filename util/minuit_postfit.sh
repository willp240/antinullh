#!/usr/bin/env bash
# Process a completed floated Minuit fit without changing the fitted inputs
set -euo pipefail
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 /path/to/completed_fit" >&2
    exit 1
fi
fit_dir=$(realpath "$1")
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_dir"
python_bin=${PYTHON:-python3}
exec > >(tee -a "$fit_dir/minuit_postfit.log") 2>&1
echo "Post-processing Minuit fit: $fit_dir"
# Quote paths as C++ string literals, then run the existing ROOT drawing tools
"$python_bin" - "$fit_dir" <<'PY'
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, 'util')
from minuit_results import read_minuit_result, write_minuit_metadata, write_minuit_summary
result = read_minuit_result(sys.argv[1])
metadata = str(write_minuit_metadata(result))
write_minuit_summary(result)

def root_call(expression):
    subprocess.run(['root', '-l', '-b', '-q', expression], check=True)

root_call('plotting/plotFitParams.C(' + json.dumps(metadata) + ')')
choices = ([0] if len(result['datasets']) == 2 else [])
choices += [1 if dataset == '2p2ppo' else 2 for dataset in result['datasets']]
for choice in choices:
    root_call('plotting/plotFitSpectra.C(' + json.dumps(metadata) + ',' + str(choice) +
              ',"",' + json.dumps(str(result['directory'])) + ',' +
              ('true' if result['observed_data'] else 'false') + ')')
plots = result['directory'] / 'plots'
expected = ['params.pdf', 'params.root', 'correlations.pdf']
for choice in choices:
    suffix = {0: 'all', 1: '2p2ppo', 2: 'bismsb'}[choice]
    expected += [prefix + '_' + suffix + extension for prefix in ('dist', 'distgroup')
                 for extension in ('.pdf', '.root')]
for filename in expected:
    path = plots / filename
    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError('Missing postfit output: ' + str(path))
print('Postfit plots: ' + str(plots))
PY
echo "Postfit processing complete: $fit_dir"
