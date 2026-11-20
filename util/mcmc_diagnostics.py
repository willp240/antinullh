#!/usr/bin/env python3
"""Calculate MCMC convergence diagnostics independently of parameter summaries"""

import argparse
import json
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import tempfile

import numpy as np
from scipy.stats import norm, rankdata
import uproot

from mcmc_read_chains import prepare_inputs, validate_sources
from mcmc_file_utils import (REPO, check_output_identity, finite, publish_products,
                            root_open, run_cli, settings, stage_directories, stage_lock)
from mcmc_plot_utils import render_diagnostics


def autocorrelation(values, max_lag=None):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.all(np.isfinite(values)):
        raise ValueError('Autocorrelation needs at least two finite scalar samples')
    if max_lag is not None and max_lag < 0:
        raise ValueError('max-lag must be non-negative')
    if np.ptp(values) == 0:
        return None
    # Subtract a reference before averaging to retain precision for nearly constant values
    offsets = values - values[0]
    centered = offsets - offsets.mean()
    length = len(values)
    # Zero padding prevents correlation between the end and the start of the chain
    size = 1 << (2 * length - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=size)
    covariance = np.fft.irfft(spectrum * spectrum.conjugate(), n=size)[:length]
    if not np.isfinite(covariance[0]) or covariance[0] <= 0:
        raise ValueError('Invalid variance for autocorrelation')
    return (covariance / covariance[0])[:None if max_lag is None else max_lag + 1]


def write_autocorrelations(cache, manifest, max_lag):
    """Calculate and save each chain separately using the full retained sample sequence"""
    for i, record in enumerate(manifest['chains']):
        print(f"Autocorrelations for chain {i + 1}/{len(manifest['chains'])}", flush=True)
        with root_open(Path(record['file'])) as source:
            for name in manifest['names']:
                values = source['posteriors'][name].array(entry_start=manifest['burn_in'], library='np')
                try:
                    correlation = autocorrelation(values, max_lag)
                except ValueError as error:
                    raise ValueError(f"{record['file']}: {name}: {error}") from error
                if correlation is None:
                    print(f'Skipping constant parameter {name}', flush=True)
                    continue
                cache[f'autocorrelations/chain_{i}/{name}'] = (
                    correlation, np.arange(len(correlation) + 1, dtype=float) - .5)



def split_chains(chains):
    length = min(len(chain) for chain in chains)
    half = length // 2
    return np.array([part for chain in chains for part in
                     (np.asarray(chain)[:half], np.asarray(chain)[length - half:length])])



def rank_normalise(values):
    # Blom normal scores with average ranks for ties
    ranks = rankdata(values.ravel()).reshape(values.shape)
    return norm.ppf((ranks - 0.375) / (values.size + 0.25))



def rhat(values):
    n = values.shape[1]
    within = np.var(values, axis=1, ddof=1).mean()
    if within == 0:
        return np.nan
    between = n * np.var(values.mean(axis=1), ddof=1)
    return np.sqrt(((n - 1) / n * within + between / n) / within)



def effective_size(values):
    m, n = values.shape
    centered = values - values.mean(axis=1, keepdims=True)
    length = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=length, axis=1)
    acov = np.fft.irfft(spectrum * spectrum.conjugate(), n=length, axis=1)[:, :n] / n
    within = acov[:, 0].mean() * n / (n - 1)
    variance = (n - 1) / n * within + np.var(values.mean(axis=1), ddof=1)
    if variance <= 0:
        return np.nan
    rho = 1 - (within - acov.mean(axis=0)) / variance
    rho[0] = 1.0
    # Geyer initial-positive, initial-monotone sequence of paired autocorrelations
    pairs = rho[:n - n % 2].reshape(-1, 2).sum(axis=1)
    nonpositive = np.flatnonzero(pairs <= 0)
    stop = int(nonpositive[0]) if len(nonpositive) else len(pairs)
    positive = np.minimum.accumulate(pairs[:stop])
    tau = max(1.0, -1 + 2 * positive.sum())
    return min(m * n, m * n / tau)



def diagnostics(chains):
    values = split_chains(chains)
    if values.shape[1] < 4 or np.ptp(values) == 0:
        return {'rhat': None, 'ess_bulk': None, 'ess_tail': None}
    ranked = rank_normalise(values)
    folded = rank_normalise(np.abs(values - np.median(values)))
    low, high = np.quantile(values, [0.05, 0.95])
    return {'rhat': finite(max(rhat(ranked), rhat(folded))),
            'ess_bulk': finite(effective_size(ranked)),
            'ess_tail': finite(min(effective_size((values <= low).astype(float)),
                                    effective_size((values <= high).astype(float))))}



def diagnostic_values(manifest):
    values = {}
    warnings = []
    names = manifest['names'] + ['LogL']
    for j, name in enumerate(names):
        print(f'Convergence diagnostics for {name} ({j + 1}/{len(names)})', flush=True)
        chains = []
        for record in manifest['chains']:
            with root_open(Path(record['file'])) as source:
                chains.append(source['posteriors'][name].array(entry_start=manifest['burn_in'], library='np'))
        values[name] = diagnostics(chains)
        value = values[name]
        if value['rhat'] is None or value['rhat'] > 1.01:
            warnings.append(f"{name}: R-hat is {value['rhat']}")
        if value['ess_bulk'] is None or value['ess_bulk'] < 100 * len(chains):
            warnings.append(f"{name}: low bulk effective sample size {value['ess_bulk']}")
    return values, warnings



def run_diagnostics(args):
    directory, manifest, input_warnings = prepare_inputs(args)
    with stage_lock(directory, 'diagnostics'):
        check_output_identity(directory, manifest)
        values, warnings = diagnostic_values(manifest)
        warnings = input_warnings + warnings
        if len({record['retained'] for record in manifest['chains']}) != 1:
            warnings.append('Diagnostics use equal-length prefixes; parameter summaries use all retained draws')
        if len(manifest['chains']) < 2:
            warnings.append('Only one chain; between-chain agreement cannot be assessed')
        with tempfile.TemporaryDirectory(prefix='.diagnostics-', dir=directory) as temporary:
            work = Path(temporary) / 'products'
            work.mkdir()
            stage_directories(work, Path(manifest['chains'][0]['file']).parent / 'cfg')
            with uproot.recreate(work / 'diagnostics.root') as cache:
                table = cache.mktree('convergence', {'name': 'string', 'rhat': 'float64', 'ess_bulk': 'float64', 'ess_tail': 'float64'})
                table.extend({'name': list(values), **{key: np.array([
                    value[key] if value[key] is not None else np.nan for value in values.values()])
                    for key in ('rhat', 'ess_bulk', 'ess_tail')}})
                stage_settings = settings('diagnostics', manifest)
                stage_settings['max_lag'] = args.max_lag
                stage_settings['autocorrelation_method'] = 'NumPy zero-padded FFT; normalised by full zero-lag sum of squared deviations'
                cache['input_manifest'] = json.dumps(manifest, allow_nan=False)
                cache['analysis'] = json.dumps({'settings': stage_settings, 'warnings': warnings}, allow_nan=False)
                write_autocorrelations(cache, manifest, args.max_lag)
            report = '\n'.join(['Convergence diagnostics', f"Step >= {manifest['burn_in']}",
                                 *warnings, '', *[f'{name}: {value}' for name, value in values.items()]])
            (work / 'logs/diagnostics.txt').write_text(report + '\n')
            if not args.no_plots:
                render_diagnostics(work, manifest, work / 'cfg')
            validate_sources(manifest)
            check_output_identity(directory, manifest)
            publish_products(work, directory, manifest)
        print('Diagnostic products: ' + str(directory), flush=True)



def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', nargs='+', type=Path, help='Chain directories or original ROOT chain files')
    parser.add_argument('--burn-in', type=int, required=True, help='Keep Step >= this value in each chain')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--max-lag', type=int, default=10000)
    args = parser.parse_args(argv)
    if args.burn_in < 0 or args.max_lag < 0:
        parser.error('burn-in and max-lag must be non-negative')
    os.chdir(REPO)
    run_diagnostics(args)


if __name__ == '__main__':
    run_cli(main, 'MCMC diagnostics')
