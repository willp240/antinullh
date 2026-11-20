#!/usr/bin/env python3
"""Submit an MCMC parameter/count summaries or diagnostics calculation to Condor"""

import argparse
from pathlib import Path
import shlex
import subprocess
import tempfile


REPO = Path(__file__).resolve().parents[1]


def condor_path(path):
    value = str(path)
    if any(character.isspace() or character == '"' for character in value):
        raise ValueError('Condor submission paths must not contain whitespace or double quotes')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='Directory containing completed chains')
    parser.add_argument('--stage', choices=['parameters', 'diagnostics'], required=True,
                        help='Choose the independent calculation to submit')
    parser.add_argument('--burn-in', type=int, required=True, help='Keep Step >= this value in each chain')
    parser.add_argument('--output-dir', type=Path, help='Defaults to output/ inside the input directory')
    parser.add_argument('--max-lag', type=int, default=10000)
    parser.add_argument('-e', '--env-file', type=Path, required=True, help='Complete batch environment setup script')
    parser.add_argument('--python', default='python3', help='Python executable in the batch environment')
    parser.add_argument('--memory', type=int, default=16384, help='Requested memory in MB (default: 16384)')
    parser.add_argument('--wall-time', type=int, default=14400, help='Maximum execution time in seconds (default: 14400)')
    parser.add_argument('--dry-run', action='store_true', help='Write submission files without submitting')
    args = parser.parse_args()
    source = args.input.expanduser().resolve()
    environment = args.env_file.expanduser().resolve()
    output = args.output_dir.expanduser().resolve() if args.output_dir else source / 'output'
    if not source.is_dir() or not next(source.glob('*/cfg/fit_config.ini'), None):
        parser.error('Input must contain completed job directories with cfg/fit_config.ini')
    if not environment.is_file():
        parser.error('Environment setup script does not exist')
    if (args.burn_in is not None and args.burn_in < 0) or args.max_lag < 0:
        parser.error('Burn-in and maximum lag must be non-negative')
    if args.memory <= 0 or args.wall_time <= 0:
        parser.error('Memory and wall time must be positive')
    script_name = {'parameters': 'mcmc_summary.py', 'diagnostics': 'mcmc_diagnostics.py'}[args.stage]
    command = [args.python, '-u', str(REPO / 'util' / script_name)]
    command += [str(source), '--burn-in', str(args.burn_in), '--output-dir', str(output)]
    if args.stage == 'diagnostics':
        command += ['--max-lag', str(args.max_lag)]
    job_root = REPO / 'output/condor_mcmc_postfit'
    job_root.mkdir(parents=True, exist_ok=True)
    job = Path(tempfile.mkdtemp(prefix=source.name + '_' + args.stage + '_', dir=job_root))
    script = job / 'run.sh'
    script.write_text('\n'.join([
        '#!/bin/bash', 'set -e',
        'source ' + shlex.quote(str(environment)),
        'cd ' + shlex.quote(str(REPO)),
        'if [[ -f env_my.sh ]]; then source env_my.sh; fi',
        'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
        'echo "Starting MCMC postfit on $(hostname) at $(date -Is)"',
        'exec ' + shlex.join(command), '',
    ]))
    script.chmod(0o755)
    submission = job / 'job.submit'
    submission.write_text('\n'.join([
        'universe = vanilla',
        'executable = ' + condor_path(script),
        'initialdir = ' + condor_path(REPO),
        'output = ' + condor_path(job / 'stdout.log'),
        'error = ' + condor_path(job / 'stderr.log'),
        'log = ' + condor_path(job / 'condor.log'),
        'getenv = True', 'should_transfer_files = NO',
        'request_cpus = 1', 'request_memory = ' + str(args.memory),
        'allowed_execute_duration = ' + str(args.wall_time),
        'notification = never', 'queue 1', '',
    ]))
    print('Command: ' + shlex.join(command), flush=True)
    print('Condor files and logs: ' + str(job), flush=True)
    print('Postfit results: ' + str(output), flush=True)
    submit = ['condor_submit', '-batch-name', 'postfit_' + source.name + '_' + args.stage, str(submission)]
    if args.dry_run:
        print('Dry run; submit with: ' + shlex.join(submit))
    else:
        subprocess.run(submit, check=True)


if __name__ == '__main__':
    main()
