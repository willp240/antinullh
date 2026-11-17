#!/usr/bin/env python3
"""Analyse independent MCMC chains and make standard postfit products"""

import argparse
import json
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import shutil
import tempfile

# Keep numerical libraries from creating large thread pools on a login node
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/antinu-matplotlib')

import numpy as np
import uproot

from mcmc_results import (read_draws, describe_chains, discover,
                          load_config, posterior_summary, root_open, summarise, write_root_summary)


REPO = Path(__file__).resolve().parents[1]


def reconstruct_fit(repo, output, arrays, manifest, result, count=None):
    warnings = []
    total = len(arrays['samples'])
    count = 0 if count is None else total if count == 0 else min(count, total)
    if not count:
        (output / 'plots/yield_posteriors.pdf').unlink(missing_ok=True)
    rng = np.random.default_rng(1729)
    # One random draw per equal-sized stratum, preserving posterior dwell-time weights
    indices = (np.floor((np.arange(count) + rng.random(count)) * total / count).astype(np.int64)
               if count else np.empty(0, dtype=np.int64))
    selected = np.asarray(arrays['samples'][indices])
    unique, inverse = np.unique(selected, axis=0, return_inverse=True)
    best = np.asarray(arrays['samples'][result['best_sample']['sample_index']])
    requests = np.vstack([best, unique])
    cfg = Path(manifest['chains'][0]['file']).parent / 'cfg'
    print(f'Reconstructing best sample and {len(unique)} unique posterior yield samples', flush=True)
    # Exchange files are temporary; retain the reconstructed spectra and pass numerical results in memory
    with tempfile.TemporaryDirectory(prefix='antinu-yields-') as temporary:
        work = Path(temporary)
        with uproot.recreate(work / 'requests.root') as handle:
            fields = {'sample_id': 'int64', **{name: 'float64' for name in manifest['names']}}
            tree = handle.mktree('requests', fields)
            tree.extend({'sample_id': np.r_[-1, np.arange(len(unique))].astype('int64'),
                         **{name: requests[:, j] for j, name in enumerate(manifest['names'])}})
        fit_path = cfg / 'fit_config.ini'
        recorded_seed = (manifest['chains'][0].get('metadata') or {}).get('random_seed')
        if recorded_seed:
            # Reuse the actual run seed when the saved config requested automatic selection
            fit_config = load_config(fit_path)
            fit_config.set('summary', 'random_seed', str(recorded_seed))
            fit_path = work / 'fit_config.ini'
            with fit_path.open('w') as handle:
                fit_config.write(handle)
        command = [str(repo / 'bin/reconstruct_mcmc_fit'), str(fit_path)]
        command += [str(cfg / name) for name in ('event_config.ini', 'pdf_config.ini', 'syst_config.ini')]
        command += [str(work / 'requests.root'), str(work)]
        with (output / 'reconstruct.log').open('w') as log:
            subprocess.run(command, check=True, cwd=repo, stdout=log, stderr=subprocess.STDOUT)
        metadata = json.loads((work / 'rate_metadata.json').read_text())
        with root_open(work / 'rates.root') as source:
            raw = source['rates'].arrays(library='np')
        if not np.array_equal(raw['sample_id'], np.r_[-1, np.arange(len(unique))]):
            raise ValueError('Rate reconstruction returned unexpected sample IDs')
        shutil.copytree(work / 'bestfit', output / 'bestfit', dirs_exist_ok=True)
    derived = {}
    nominal = {}
    best_yields = {}
    bases = {}
    all_components = {column['component'] for column in metadata['columns'].values()}
    for branch, description in metadata['columns'].items():
        values = raw[branch][1:][inverse]
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError(f'Invalid derived yield in {branch}')
        component = description['component']
        key = component + '_' + description['dataset']
        derived[key] = values
        nominal[key] = description['nominal']
        best_yields[key] = float(raw[branch][0])
        if not np.isfinite(best_yields[key]) or best_yields[key] < 0:
            raise ValueError(f'Invalid best-sample yield in {branch}')
        base = component[:-1] if component.endswith('2') and component[:-1] in all_components else component
        bases.setdefault(base, []).append(key)
    for base, parts in bases.items():
        derived[base + '_total'] = sum(derived[key] for key in parts)
        nominal[base + '_total'] = sum(nominal[key] for key in parts)
        best_yields[base + '_total'] = sum(best_yields[key] for key in parts)
    for dataset, suffix in [('2p2ppo', 'ppo'), ('bismsb', 'bismsb'), ('total', 'total')]:
        if dataset == 'total':
            uranium, thorium = 'geonu_U_total', 'geonu_Th_total'
        else:
            number = '2' if dataset == 'bismsb' else ''
            uranium, thorium = 'geonu_U' + number + '_' + dataset, 'geonu_Th' + number + '_' + dataset
        if uranium in derived and thorium in derived:
            if best_yields[thorium] > 0:
                best_yields['geo_ratio_' + suffix] = best_yields[uranium] / best_yields[thorium]
            if np.any(derived[thorium] <= 0):
                warnings.append(f'U/Th ratio {suffix} omitted because a sampled Th yield is zero')
            else:
                derived['geo_ratio_' + suffix] = derived[uranium] / derived[thorium]
                nominal['geo_ratio_' + suffix] = nominal[uranium] / nominal[thorium]
    summaries = {name: posterior_summary(values) for name, values in derived.items()} if count else {}
    rate_result = {'method': 'Full model evaluation, including shapes and oscillation',
                   'selection': 'Stratified random subsample with original multiplicities' if count else 'Best retained sample only',
                   'selected_draws': count, 'retained_draws': total, 'all_draws_evaluated': count == total,
                   'selection_seed': 1729, 'benchmark': metadata,
                   'parameters': summaries, 'best_fit': best_yields, 'nominal': nominal, 'warnings': warnings,
                   'interval_note': 'Check interval stability as --rate-samples increases; zero evaluates every draw' if count else 'No physical-yield uncertainty calculation requested'}
    if count:
        chain_index = np.searchsorted(np.cumsum([r['retained'] for r in manifest['chains']]), indices, side='right')
        derived.update(sample_index=indices, chain_index=chain_index.astype('int32'),
                       Step=np.asarray(arrays['steps'][indices]))
    else:
        derived = {}
    return rate_result, derived


def write_summary_table(output, result, rates=None):
    parameters = {name: values for name, values in result['parameters'].items() if name != 'LogL'}
    if rates:
        parameters.update(rates['parameters'])
    # Use the established table row names and labels for consistency with Minuit tables
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tables'))
    from make_old_new_fit_tables import ROWS, fmt_value
    rows = []
    grouped = {'osc': [], 'syst': []}
    for key, label, category in ROWS:
        if key in parameters:
            values = parameters[key]
            cell = fmt_value(key, values['median'], (values['error_low'], values['error_high']), category)
            rows.append(label + ' & ' + cell + r' \\')
            grouped[category].append(rows[-1])
    table = '\n'.join([r'\begin{tabular}{lr}', r'\hline', r'Parameter & MCMC posterior \\', r'\hline',
                       *rows, r'\hline', r'\end{tabular}',
                       '% Posterior medians and equal-tailed 68% credible intervals'])
    (output / 'parameter_table.tex').write_text(table + '\n')
    for category, category_rows in grouped.items():
        (output / ('parameter_table_' + category + '.tex')).write_text('\n'.join([
            r'\begin{tabular}{lr}', r'\hline', r'Parameter & MCMC posterior \\', r'\hline',
            *category_rows, r'\hline', r'\end{tabular}']) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', nargs='+', help='Batch directories or original ROOT chain files')
    parser.add_argument('--burn-in', type=int, required=True, help='Keep Step >= this value in each chain')
    parser.add_argument('--output-dir', type=Path,
                        help='Defaults to output/ inside the first input directory, or beside the first input file')
    parser.add_argument('--max-lag', type=int, default=10000)
    parser.add_argument('--rate-samples', type=int, default=None,
                        help='Also calculate yield uncertainties from N retained draws; zero uses all draws')
    parser.add_argument('--skip-rates', action='store_true', help='Parameter analysis and diagnostics only')
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--pairs', nargs='*', default=[], help='Parameter pairs a:b, or all')
    args = parser.parse_args()
    if args.burn_in < 0 or args.max_lag < 0 or (args.rate_samples is not None and args.rate_samples < 0):
        parser.error('Burn-in, maximum lag and rate sample count must be non-negative')
    if args.skip_rates and args.rate_samples is not None:
        parser.error('--skip-rates cannot be combined with --rate-samples')
    inputs = [str(Path(path).resolve()) for path in args.inputs]
    first_input = Path(inputs[0])
    input_directory = first_input.parent if first_input.is_file() else first_input
    output = (args.output_dir if args.output_dir is not None else input_directory / 'output').resolve()
    # Relative input config paths retain the same meaning as during fitting
    os.chdir(REPO)
    files = discover(inputs)
    names, records, warnings = describe_chains(files, args.burn_in)
    output.mkdir(parents=True, exist_ok=True)
    for warning in warnings:
        print('Warning: ' + warning, file=sys.stderr, flush=True)
    (output / 'cfg').mkdir(exist_ok=True)
    for path in (Path(records[0]['file']).parent / 'cfg').glob('*.ini'):
        shutil.copy2(path, output / 'cfg' / path.name)
    if len(records) < 2:
        warnings.append('Only one original chain supplied; between-chain agreement cannot be assessed')
    print(f'Found {len(records)} chains and {len(names)} free parameters', flush=True)
    arrays, manifest = read_draws(names, records, args.burn_in)
    result, slices = summarise(arrays, manifest)
    result['warnings'] = warnings + result['warnings']
    if len({r['retained'] for r in records}) != 1:
        result['warnings'].append('Diagnostics use equal-length prefixes; posterior summaries use all retained samples')
    rates = None
    derived = {}
    if not args.skip_rates:
        if not (REPO / 'bin/reconstruct_mcmc_fit').is_file():
            raise RuntimeError('Build bin/reconstruct_mcmc_fit before reconstructing spectra and yields')
        rates, derived = reconstruct_fit(REPO, output, arrays, manifest, result,
                                          args.rate_samples)
        result['warnings'] += rates['warnings']
    write_summary_table(output, result, rates)
    source_versions = {}
    for name, directory in [('antinu_llh', REPO), ('oxo', Path(os.environ.get('OXO_DIR', '/home/parkerw/Software/oxsx3')))]:
        revision = subprocess.run(['git', '-C', str(directory), 'rev-parse', 'HEAD'], capture_output=True, text=True)
        source_versions[name] = revision.stdout.strip() if revision.returncode == 0 else 'unknown'
    settings = {'schema': 1, 'command': sys.argv,
        'source_versions': source_versions, 'source_hashes': {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in list((REPO / 'util').glob('mcmc_*.py')) + [Path(__file__).resolve()]},
        'inputs': inputs, 'burn_in': args.burn_in, 'rate_samples': args.rate_samples,
        'score': 'LogL = -NLL including constraints; larger is better',
        'parameter_intervals': 'Median and equal-tailed 68% interval',
        'rates': 'Reconstructed' if rates else 'Not reconstructed'}
    write_root_summary(output, manifest, result, rates, derived, settings)
    if not args.no_plots:
        from mcmc_plots import render
        render(REPO, output, arrays, manifest, result, slices, rates, derived,
               args.max_lag, args.pairs)
    report = [f"Analysed {len(records)} chains; {len(arrays['samples'])} retained draws; Step >= {args.burn_in}",
              f"Highest-score retained sample: {result['best_sample']['file']}, entry {result['best_sample']['entry']}",
              '', 'Diagnostics:', *['- ' + warning for warning in result['warnings']]]
    if rates:
        benchmark = rates['benchmark']
        report += ['', 'Physical yields and spectra evaluated at the highest-score retained sample',
                   f"Model setup: {benchmark['setup_seconds']:.1f} s; evaluation: {benchmark['seconds_per_evaluation']:.3f} s/point"]
        if rates['selected_draws']:
            report += [f"Yield uncertainties: {rates['selected_draws']} draws from {rates['retained_draws']}",
                       'Check derived intervals for stability as --rate-samples increases']
    (output / 'report.txt').write_text('\n'.join(report) + '\n')
    print('\n'.join(report), flush=True)
    print(f'Postfit outputs: {output}', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as error:
        print(f'MCMC postfit failed: {error}', file=sys.stderr)
        sys.exit(1)
