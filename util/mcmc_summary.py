#!/usr/bin/env python3
"""Summarise MCMC parameters and event counts from saved best-point spectra"""

import argparse
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import shutil
import tempfile
import time
import warnings

import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning

from mcmc_read_chains import BOOKKEEPING, prepare_inputs, read_draws, validate_sources
from mcmc_file_utils import (REPO, check_output_identity, comma_list, load_config,
                            publish_products, root_open, run_cli, settings,
                            stage_directories, stage_lock, write_root_summary, write_summary_table)
from mcmc_plot_utils import render_parameters, render_rates


def posterior_summary(values):
    low, median, high = np.quantile(values, [0.16, 0.5, 0.84])
    return {'median': float(median), 'interval_low': float(low), 'interval_high': float(high),
            'error_low': float(median - low), 'error_high': float(high - median),
            'mean': float(np.mean(values)), 'sd': float(np.std(values, ddof=1))}



def highest_density_regions(counts, edges, probability=.68):
    """Select the densest equal-width histogram bins covering at least the requested mass"""
    order = np.argsort(-counts, kind='stable')
    stop = np.searchsorted(np.cumsum(counts[order]), probability * counts.sum()) + 1
    selected = np.zeros(len(counts), dtype=bool)
    selected[order[:stop]] = True
    transitions = np.diff(np.r_[False, selected, False].astype(int))
    starts, stops = np.flatnonzero(transitions == 1), np.flatnonzero(transitions == -1)
    regions = [[float(edges[start]), float(edges[end])] for start, end in zip(starts, stops)]
    return regions, float(counts[selected].sum() / counts.sum())



def posterior_estimates(values):
    """Additional histogram-based estimates for comparison with sample summaries"""
    counts, edges = np.histogram(values, bins=80)
    peak = int(np.argmax(counts))
    mode = float((edges[peak] + edges[peak + 1]) / 2)
    mean, sd = float(np.mean(values)), float(np.std(values, ddof=1))
    regions, mass = highest_density_regions(counts, edges)
    result = {'mode': mode, 'hpd_regions': regions, 'hpd_mass': mass,
              'gaussian_mean': None, 'gaussian_sigma': None, 'gaussian_amplitude': None}
    if not np.isfinite(sd) or sd == 0:
        result['mode'] = mean
        result['hpd_regions'] = [[mean, mean]]
        result['hpd_mass'] = 1.0
        return result
    # Standardise the fit coordinates so very small oscillation parameters are well conditioned
    centres = ((edges[:-1] + edges[1:]) / 2 - mean) / sd
    def gaussian(x, amplitude, centre, width):
        return amplitude * np.exp(-0.5 * ((x - centre) / width) ** 2)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error', OptimizeWarning)
            fitted, unused = curve_fit(gaussian, centres, counts,
                p0=[float(counts.max()), 0, 1], sigma=np.sqrt(np.maximum(counts, 1)),
                bounds=([0, centres[0], 0.02], [np.inf, centres[-1], 2 * (centres[-1] - centres[0])]),
                absolute_sigma=True, maxfev=10000)
        if np.all(np.isfinite(fitted)):
            result['gaussian_mean'] = float(mean + fitted[1] * sd)
            result['gaussian_sigma'] = float(fitted[2] * sd)
            result['gaussian_amplitude'] = float(fitted[0])
    except (RuntimeError, ValueError, OptimizeWarning):
        pass
    return result



def summarise(arrays, manifest):
    names = manifest['names']
    samples = arrays['samples']
    records = manifest['chains']
    slices = []
    offset = 0
    chain_info = []
    for record in records:
        selected = slice(offset, offset + record['retained'])
        slices.append(selected)
        chain_info.append({'chain_id': record['chain_id'], 'file': record['file'],
                           'retained': record['retained'],
                           'acceptance_after_cut': float(arrays['accepted'][selected].mean())})
        offset += record['retained']
    summaries = {}
    warnings = []
    for j, name in enumerate(names):
        print(f'Summarising {name} ({j + 1}/{len(names) + 1})', flush=True)
        values = samples[:, j]
        summary = posterior_summary(values)
        summary.update(posterior_estimates(values))
        summary.update({'rhat': None, 'ess_bulk': None, 'ess_tail': None})
        summary['per_chain'] = [posterior_summary(samples[s, j]) for s in slices]
        summary['best_fit'] = None
        summaries[name] = summary
    print(f'Summarising LogL ({len(names) + 1}/{len(names) + 1})', flush=True)
    summaries['LogL'] = posterior_summary(arrays['log_score'])
    summaries['LogL'].update(posterior_estimates(arrays['log_score']))
    summaries['LogL'].update({'rhat': None, 'ess_bulk': None, 'ess_tail': None})
    best_index = int(np.argmax(arrays['log_score']))
    chain_index = next(i for i, s in enumerate(slices) if s.start <= best_index < s.stop)
    record = records[chain_index]
    step = int(arrays['steps'][best_index])
    with root_open(Path(record['file'])) as source:
        best_arrays = source['posteriors'].arrays(entry_start=step, entry_stop=step + 1, library='np')
    best = {'sample_index': best_index, 'chain_index': chain_index, 'file': record['file'],
            'entry': step, 'step': step, 'log_score': float(arrays['log_score'][best_index]),
            'parameters': {name: float(values[0]) for name, values in best_arrays.items() if name not in BOOKKEEPING}}
    # Accumulate covariance in chunks to avoid copying the full pooled sample matrix
    mean = np.array([summaries[name]['mean'] for name in names])
    covariance = np.zeros((len(names), len(names)))
    for start in range(0, len(samples), 100000):
        centered = samples[start:start + 100000] - mean
        covariance += centered.T @ centered
    covariance /= len(samples) - 1
    scale = np.sqrt(np.diag(covariance))
    correlation = np.divide(covariance, np.outer(scale, scale),
                            out=np.zeros_like(covariance), where=np.outer(scale, scale) > 0)
    for name in names:
        summaries[name]['best_fit'] = best['parameters'][name]
    summaries['LogL']['best_fit'] = best['log_score']
    result = {'method': 'Highest-LogL retained sample and central 68% interval (16%, 84%); interval need not contain the best point',
              'burn_in': manifest['burn_in'], 'chains': chain_info, 'parameters': summaries,
              'diagnostics_draws_per_chain': None,
              'diagnostics_method': 'Stored separately in diagnostics.root',
              'best_sample': best, 'warnings': warnings}
    result['covariance'] = covariance
    result['correlation'] = correlation
    return result, slices



def check_saved_best(record, best):
    metadata = record.get('metadata') or {}
    saved = metadata.get('best_fit')
    if saved is not None:
        if any(name not in saved or saved[name] != value for name, value in best['parameters'].items()):
            raise ValueError('Saved distributions do not match the highest-LogL retained sample: ' + record['file'])
        return
    # Legacy runs only recorded the distribution parameter set as rounded text
    path = Path(record['file']).parent / 'fit_result.txt'
    if not path.is_file():
        raise ValueError('Cannot verify saved distributions without best_fit metadata or ' + str(path))
    saved = {}
    for line in path.read_text().split('Best Fit Values:', 1)[-1].split('LLH:', 1)[0].splitlines():
        fields = line.split()
        if len(fields) == 2:
            try:
                saved[fields[0]] = float(fields[1])
            except ValueError:
                pass
    if any(name not in saved or not np.isclose(saved[name], value, rtol=5.1e-5, atol=0.)
           for name, value in best['parameters'].items()):
        raise ValueError('Saved distributions do not match the highest-LogL retained sample: ' + record['file'])



def integrate_nonbuffer(path, pdf, metadata):
    axes = comma_list(pdf, 'summary', 'data_axes')
    with root_open(path) as source:
        histograms = [obj for _, obj in source.items() if hasattr(obj, 'to_numpy') and obj.classname.startswith('TH')]
        if len(histograms) != 1:
            raise ValueError('Expected exactly one component histogram in ' + str(path))
        values, *edges = histograms[0].to_numpy(flow=False)
    if values.ndim != len(axes):
        raise ValueError('Saved histogram dimensions do not match configured data axes: ' + str(path))
    selected = []
    for axis, bins in zip(axes, edges):
        n = pdf.getint(axis, 'n_bins')
        expected = np.linspace(pdf.getfloat(axis, 'min'), pdf.getfloat(axis, 'max'), n + 1)
        if len(bins) != len(expected) or not np.allclose(bins, expected, rtol=1e-10, atol=1e-10):
            raise ValueError('Saved histogram binning does not match configuration: ' + str(path))
        buffers = metadata.get('llh_buffer_bins', {}).get(axis)
        if buffers is None:
            low = pdf.getint(axis, 'llh_buffer_bin_low', fallback=None)
            high = pdf.getint(axis, 'llh_buffer_bin_high', fallback=None)
            if axis == 'energy' and (low is None or high is None):
                raise ValueError('Legacy energy buffer settings are missing from saved pdf_config.ini; '
                                 'set llh_buffer_bin_low and llh_buffer_bin_high to the values used by that run')
            buffers = [low or 0, high or 0]
        low, high = buffers
        # Match the fitter's ROOT bin convention: low is the first included bin
        start, stop = max(1, low) - 1, min(n, n - high)
        if low < 0 or high < 0 or stop <= start:
            raise ValueError('Invalid likelihood buffer settings for ' + axis)
        selected.append(slice(start, stop))
    value = float(values[tuple(selected)].sum())
    if not np.isfinite(value) or value < 0:
        raise ValueError('Invalid non-buffer event count in ' + str(path))
    return value



def saved_distribution_rates(output, arrays, manifest, result):
    """Integrate saved spectra and propagate normalisations at fixed shapes and oscillation"""
    best = result['best_sample']
    record = manifest['chains'][best['chain_index']]
    check_saved_best(record, best)
    source = Path(record['file']).parent
    config = load_config(source / 'cfg/fit_config.ini')
    events = load_config(source / 'cfg/event_config.ini')
    pdf = load_config(source / 'cfg/pdf_config.ini')
    syst = load_config(source / 'cfg/syst_config.ini')
    metadata = record.get('metadata') or {}
    if 'best_fit' not in metadata:
        result['warnings'].append('Saved spectrum parameters checked against rounded legacy fit_result.txt values')
    active_systs = comma_list(syst, 'summary', 'active')
    if active_systs == ['all']:
        active_systs = [name for name in syst.sections() if name != 'summary']
    derived, best_yields, nominal, bases = {}, {}, {}, {}
    inputs = []
    for dataset in comma_list(events, 'summary', 'datasets'):
        for component in comma_list(events, dataset, 'active'):
            path = source / 'postfit_dists' / (component + '_' + dataset + '.root')
            count = integrate_nonbuffer(path, pdf, metadata)
            inputs.append(path)
            groups = comma_list(events, component, 'groups')
            factors = [component]
            for name in active_systs:
                kind = syst.get(name, 'type')
                if kind.startswith('Shape:') and kind.split(':', 1)[1].endswith('_norm'):
                    datasets = comma_list(syst, name, 'dataset')
                    group = syst.get(name, 'group', fallback='')
                    if (not datasets or dataset in datasets) and (not group or group in groups):
                        factors += comma_list(syst, name, 'param_names') or [name]
            ratios = np.ones(len(arrays['samples']))
            nominal_ratio = 1.
            for name in dict.fromkeys(factors):
                if name not in manifest['names']:
                    continue
                value = best['parameters'][name]
                if value <= 0:
                    raise ValueError('Cannot scale count uncertainty from a non-positive best normalisation: ' + name)
                ratios *= arrays['samples'][:, manifest['names'].index(name)] / value
                nominal_ratio *= config.getfloat(name, 'nom') / value
            key = component + '_' + dataset
            derived[key] = count * ratios
            best_yields[key] = count
            nominal[key] = count * nominal_ratio
            base = component[:-1] if component.endswith('2') and events.has_section(component[:-1]) else component
            bases.setdefault(base, []).append(key)
    for base, parts in bases.items():
        key = base + '_total'
        derived[key] = sum(derived[part] for part in parts)
        best_yields[key] = sum(best_yields[part] for part in parts)
        nominal[key] = sum(nominal[part] for part in parts)
    for dataset, suffix in [('2p2ppo', 'ppo'), ('bismsb', 'bismsb'), ('total', 'total')]:
        number = '2' if dataset == 'bismsb' else ''
        uranium, thorium = (('geonu_U_total', 'geonu_Th_total') if dataset == 'total' else
                            ('geonu_U' + number + '_' + dataset, 'geonu_Th' + number + '_' + dataset))
        if uranium in derived and thorium in derived and np.all(derived[thorium] > 0) and best_yields[thorium] > 0:
            key = 'geo_ratio_' + suffix
            derived[key] = derived[uranium] / derived[thorium]
            best_yields[key] = best_yields[uranium] / best_yields[thorium]
            nominal[key] = nominal[uranium] / nominal[thorium] if nominal[thorium] > 0 else 0.
    summaries = {}
    for name, values in derived.items():
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError('Invalid scaled event counts: ' + name)
        summaries[name] = {**posterior_summary(values), **posterior_estimates(values), 'best_fit': best_yields[name]}
    # Copy the selected chain's saved inputs and spectra for standalone redraws
    destination = output / 'dists'
    destination.mkdir()
    for directory in ('postfit_dists', 'unscaled_dists'):
        if (source / directory).is_dir():
            shutil.copytree(source / directory, destination / directory)
    for path in source.glob('data_*.root'):
        shutil.copy2(path, destination / path.name)
    return {'method': 'Saved non-buffer spectrum integrals; normalisation uncertainties at fixed shapes and oscillation',
            'parameters': summaries, 'best_fit': best_yields, 'nominal': nominal,
            'source_chain': record['file'], 'source_step': best['step'],
            'distribution_files': [str(path) for path in inputs],
            'legacy_best_fit_check': 'best_fit' not in metadata}, derived



def run_summary(args):
    directory, manifest, warnings = prepare_inputs(args)
    names, records = manifest['names'], manifest['chains']
    with stage_lock(directory, 'parameters'):
        started = time.monotonic()
        check_output_identity(directory, manifest)
        print(f'Found {len(records)} chains and {len(names)} free parameters', flush=True)
        arrays, manifest_draws = read_draws(names, records, args.burn_in)
        print(f"Summarising all {len(arrays['samples'])} retained draws; convergence diagnostics run separately", flush=True)
        result, slices = summarise(arrays, manifest_draws)
        result['analysis_id'] = manifest['analysis_id']
        result['warnings'] = warnings
        if len(records) < 2:
            result['warnings'].append('Only one chain; between-chain agreement cannot be assessed')
        with tempfile.TemporaryDirectory(prefix='.parameters-', dir=directory) as temporary:
            work = Path(temporary)
            stage_directories(work, Path(records[0]['file']).parent / 'cfg')
            rates_result, scaled_counts = saved_distribution_rates(work, arrays, manifest, result)
            write_root_summary(work, manifest, result, rates_result, scaled_counts, settings('parameters', manifest),
                               plot_data=(arrays, slices, args.pairs, args.oscillation_bins))
            write_summary_table(work, result)
            write_summary_table(work, result, rates_result)
            if not args.no_plots:
                print('Rendering parameter plots from saved numerical products', flush=True)
                render_parameters(REPO, work, manifest, result)
                render_rates(REPO, work, rates_result)
            validate_sources(manifest)
            publish_products(work, directory, manifest)
        (directory / 'report.txt').unlink(missing_ok=True)
        print(f'Parameter products: {directory} ({time.monotonic() - started:.1f} s)', flush=True)
        for warning in result['warnings']:
            print('Warning: ' + warning, flush=True)



def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', nargs='+', type=Path, help='Chain directories or original ROOT chain files')
    parser.add_argument('--burn-in', type=int, required=True, help='Keep Step >= this value in each chain')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--pairs', nargs='*', default=[])
    parser.add_argument('--oscillation-bins', type=int, default=160, help='Bins per axis for the oscillation posterior (default: 160)')
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args(argv)
    if args.burn_in < 0:
        parser.error('burn-in must be non-negative')
    if args.oscillation_bins < 2:
        parser.error('oscillation-bins must be at least two')
    os.chdir(REPO)
    run_summary(args)


if __name__ == '__main__':
    run_cli(main, 'MCMC summary')
