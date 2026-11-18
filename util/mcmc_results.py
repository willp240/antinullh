"""Read independent ROOT chains and summarise retained draws"""

import configparser
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm, rankdata
from scipy.optimize import curve_fit, OptimizeWarning
import warnings
import uproot


CONFIG_FILES = ('fit_config.ini', 'event_config.ini', 'pdf_config.ini', 'syst_config.ini')
BOOKKEEPING = {'Step', 'StepTime', 'Accepted', 'LogL'}
# Ignore legacy burn_in and save_projections settings in older chain configs
SAMPLER_OPTIONS = {'iterations', 'burn_in', 'output_directory', 'random_seed',
                   'save_outputs', 'save_projections', 'sigma_scale'}
ORDER = ['deltam21', 'sinsqtheta12', 'sintheta12', 'theta12', 'reactor_nubar_norm',
         'geonu_U_norm', 'geonu_Th_norm', 'alphan_CScatter_norm', 'alphan_OExcited_norm',
         'alphan_PRecoil_norm', 'bipolike_norm', 'atmospheric_norm', 'energy_scale',
         'energy_conv', 'birks_constant', 'p_recoil_energy_scale', 'energy_scale2',
         'energy_conv2', 'birks_constant2', 'p_recoil_energy_scale2', 'class_a_ppo',
         'class_a_bismsb', 'class_s_ppo', 'class_s_bismsb']


def root_open(path):
    # File objects and synchronous executors also work under restricted batch environments
    return uproot.open(path.open('rb'), decompression_executor=uproot.TrivialExecutor(),
                       interpretation_executor=uproot.TrivialExecutor())


def finite(value):
    return float(value) if np.isfinite(value) else None


def load_config(path):
    config = configparser.ConfigParser(interpolation=None)
    if not config.read(path):
        raise ValueError(f'Missing config {path}')
    return config


def target_config(config, fit=False):
    return {section: {key: value for key, value in config[section].items()
                      if not (fit and ((section == 'summary' and key in SAMPLER_OPTIONS)
                                      or (section != 'summary' and key == 'sig')))}
            for section in config.sections()}


def discover(inputs):
    files = []
    for item in inputs:
        path = Path(item).resolve()
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            for config in sorted(path.rglob('cfg/fit_config.ini')):
                for candidate in sorted(config.parent.parent.glob('*.root')):
                    with root_open(candidate) as source:
                        if 'posteriors' in source and hasattr(source['posteriors'], 'num_entries'):
                            files.append(candidate)
        else:
            raise ValueError(f'Input does not exist: {path}')
    if not files:
        raise ValueError('No posteriors trees found beside saved cfg directories')
    if len(set(files)) != len(files):
        raise ValueError('A chain was supplied more than once')
    return sorted(files)


def describe_chains(files, burn_in):
    records = []
    reference = None
    names = None
    warnings = []
    for path in files:
        directory = path.parent
        configs = {name: load_config(directory / 'cfg' / name) for name in CONFIG_FILES}
        signature = {name: target_config(config, name == 'fit_config.ini')
                     for name, config in configs.items()}
        if reference is not None and signature != reference:
            raise ValueError(f'Incompatible likelihood configs: {path}')
        reference = signature
        fit = configs['fit_config.ini']
        with root_open(path) as source:
            if 'posteriors' not in source or not hasattr(source['posteriors'], 'num_entries'):
                raise ValueError(f'No posteriors tree in {path}')
            tree = source['posteriors']
            metadata = json.loads(str(source['mcmc_metadata'])) if 'mcmc_metadata' in source else None
            if not BOOKKEEPING <= set(tree.keys()):
                raise ValueError(f'Missing bookkeeping branches in {path}')
            parameters = set(tree.keys()) - BOOKKEEPING
            unknown = parameters - set(fit.sections())
            if unknown:
                raise ValueError(f'Parameters missing from fit config: {unknown}')
            free = [name for name in fit.sections() if name != 'summary'
                    and not fit.getboolean(name, 'fix', fallback=False)]
            if metadata and 'free_parameters' in metadata:
                free = metadata['free_parameters']
            if not set(free) <= parameters:
                raise ValueError(f'Missing free parameters in {path}: {set(free) - parameters}')
            free.sort(key=lambda name: (ORDER.index(name) if name in ORDER else len(ORDER), name))
            if names is not None and names != free:
                raise ValueError(f'Inconsistent free parameters in {path}')
            names = free
            steps = tree['Step'].array(library='np')
            if len(steps) != fit.getint('summary', 'iterations'):
                raise ValueError(f'Incomplete chain: {path}: {len(steps)} entries')
            if len(steps) == 0 or steps[0] != 0 or not np.all(np.diff(steps) == 1):
                raise ValueError(f'Step must run consecutively from zero in each original chain: {path}')
            retained = int(np.count_nonzero(steps >= burn_in))
            if retained < 4:
                raise ValueError(f'Fewer than four retained samples in {path}')
        logs = list((directory / 'log').glob('*.log'))
        if logs and any('Normal termination (return value 0)' not in p.read_text() for p in logs):
            raise ValueError(f'Job log does not show successful completion: {directory}')
        if not metadata:
            warnings.append(f'{directory.name}: legacy chain has no seed/version metadata')
        records.append({'file': str(path), 'chain_id': directory.name, 'entries': len(steps),
                        'retained': retained, 'size': path.stat().st_size,
                        'mtime_ns': path.stat().st_mtime_ns, 'metadata': metadata,
                        'configs': {name: (directory / 'cfg' / name).read_text() for name in CONFIG_FILES}})
    seeds = [r['metadata']['random_seed'] for r in records
             if r['metadata'] and r['metadata'].get('random_seed') is not None]
    if len(seeds) != len(set(seeds)):
        warnings.append('Multiple chains report the same random seed; check independence')
    for key in ('antinu_revision', 'oxo_revision', 'antinu_diff_hash', 'root_version'):
        versions = {r['metadata'].get(key) for r in records if r['metadata'] and r['metadata'].get(key)}
        if len(versions) > 1:
            warnings.append(f'Mixed {key} values in the supplied chains')
    return names, records, warnings


def read_draws(names, records, burn_in):
    manifest = {'burn_in': burn_in, 'names': names, 'chains': records}
    total = sum(r['retained'] for r in records)
    arrays = {name: np.empty(shape, dtype=dtype)
              for name, dtype, shape in zip(('samples', 'steps', 'log_score', 'accepted'),
                  ('float64', 'int64', 'float64', 'bool'),
                  ((total, len(names)), (total,), (total,), (total,)))}
    arrays['traces'] = []
    offset = 0
    for record in records:
        preview = {name: [] for name in names + ['Step', 'LogL']}
        stride = max(1, record['entries'] // 3000)
        with root_open(Path(record['file'])) as source:
            expressions = names + ['Step', 'LogL', 'Accepted']
            for chunk in source['posteriors'].iterate(expressions, library='np', step_size='32 MB'):
                draw = chunk['Step'] % stride == 0
                for name in preview:
                    preview[name].append(chunk[name][draw])
                keep = chunk['Step'] >= burn_in
                count = int(keep.sum())
                selected = slice(offset, offset + count)
                matrix = np.column_stack([chunk[name][keep] for name in names])
                if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(chunk['LogL'][keep])):
                    raise ValueError(f'Non-finite posterior samples in {record["file"]}')
                arrays['samples'][selected] = matrix
                for key, branch in [('steps', 'Step'), ('log_score', 'LogL'), ('accepted', 'Accepted')]:
                    arrays[key][selected] = chunk[branch][keep]
                offset += count
        arrays['traces'].append({name: np.concatenate(parts) for name, parts in preview.items()})
    return arrays, manifest


def autocorrelation(values, max_lag=None):
    values = np.asarray(values, dtype=float)
    centered = values - values.mean()
    if np.ptp(values) == 0:
        return None
    length = len(values)
    size = 1 << (2 * length - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=size)
    covariance = np.fft.irfft(spectrum * spectrum.conjugate(), n=size)[:length]
    return (covariance / covariance[0])[:None if max_lag is None else max_lag + 1]


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
        values = samples[:, j]
        summary = posterior_summary(values)
        summary.update(posterior_estimates(values))
        summary.update(diagnostics([samples[s, j] for s in slices]))
        summary['per_chain'] = [posterior_summary(samples[s, j]) for s in slices]
        summaries[name] = summary
        if summary['rhat'] is None or summary['rhat'] > 1.01:
            warnings.append(f'{name}: R-hat is {summary["rhat"]}')
        if summary['ess_bulk'] is None or summary['ess_bulk'] < 100 * len(records):
            warnings.append(f'{name}: low bulk effective sample size {summary["ess_bulk"]}')
    summaries['LogL'] = posterior_summary(arrays['log_score'])
    summaries['LogL'].update(posterior_estimates(arrays['log_score']))
    summaries['LogL'].update(diagnostics([arrays['log_score'][s] for s in slices]))
    if summaries['LogL']['rhat'] is None or summaries['LogL']['rhat'] > 1.01:
        warnings.append(f'LogL: R-hat is {summaries["LogL"]["rhat"]}')
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
    result = {'method': 'MCMC posterior; median and equal-tailed 68% interval (16%, 84%)',
              'burn_in': manifest['burn_in'], 'chains': chain_info, 'parameters': summaries,
              'diagnostics_draws_per_chain': min(r['retained'] for r in records),
              'diagnostics_method': 'Rank-normalised folded split R-hat; bulk and 5%/95% tail ESS; Geyer positive monotone sequence; ESS capped at draw count',
              'best_sample': best, 'warnings': warnings}
    result['covariance'] = covariance
    result['correlation'] = correlation
    return result, slices


def write_root_summary(output, manifest, result, rates, derived, settings):
    # Import PyROOT only for the existing ROOT metadata contract and renderers
    import ROOT
    ROOT.gROOT.SetBatch(True)
    names = manifest['names']
    config = load_config(Path(manifest['chains'][0]['file']).parent / 'cfg/fit_config.ini')
    path = output / 'postfit.root'
    source = ROOT.TFile(str(path), 'RECREATE')
    for key, values, type_name in [
        ('param_names', names, 'string'),
        ('all_param_names', [name for name in config.sections() if name != 'summary'], 'string'),
        ('all_tex_labels', [config.get(name, 'tex_label', fallback=name).strip('"')
                            for name in config.sections() if name != 'summary'], 'string'),
        ('tex_labels', [config.get(name, 'tex_label', fallback=name).strip('"') for name in names], 'string'),
        ('param_fit_values', [result['parameters'][name]['median'] for name in names], 'double'),
        ('param_fit_err', [result['parameters'][name]['sd'] for name in names], 'double'),
        ('posterior_mean', [result['parameters'][name]['mean'] for name in names], 'double'),
        ('posterior_sd', [result['parameters'][name]['sd'] for name in names], 'double'),
        ('posterior_gaussian_mean', [result['parameters'][name].get('gaussian_mean') if result['parameters'][name].get('gaussian_mean') is not None else np.nan for name in names], 'double'),
        ('posterior_gaussian_sigma', [result['parameters'][name].get('gaussian_sigma') if result['parameters'][name].get('gaussian_sigma') is not None else np.nan for name in names], 'double'),
        ('posterior_mode', [result['parameters'][name].get('mode', np.nan) for name in names], 'double'),
        ('posterior_best', [result['best_sample']['parameters'][name] for name in names], 'double'),
        ('posterior_error_low', [result['parameters'][name]['error_low'] for name in names], 'double'),
        ('posterior_error_high', [result['parameters'][name]['error_high'] for name in names], 'double'),
        ('param_nom_values', [config.getfloat(name, 'nom') for name in names], 'double'),
        ('constr_mean_values', [config.getfloat(name, 'constraint_mean', fallback=config.getfloat(name, 'nom')) for name in names], 'double'),
        ('constr_sigma_values', [config.getfloat(name, 'constraint_sigma', fallback=0) for name in names], 'double'),
    ]:
        vector = ROOT.std.vector(type_name)()
        for value in values:
            vector.push_back(value)
        source.WriteObject(vector, key)
    region_names, region_low, region_high = [], [], []
    for name in names:
        for low, high in result['parameters'][name].get('hpd_regions', []):
            region_names.append(name)
            region_low.append(low)
            region_high.append(high)
    for key, values, type_name in [('posterior_mode_hpd_names', region_names, 'string'),
                                    ('posterior_mode_hpd_low', region_low, 'double'),
                                    ('posterior_mode_hpd_high', region_high, 'double')]:
        vector = ROOT.std.vector(type_name)()
        for value in values:
            vector.push_back(value)
        source.WriteObject(vector, key)
    for key, values in [('covMatrix', result['covariance']), ('corrMatrix', result['correlation'])]:
        matrix = ROOT.TMatrixD(len(names), len(names))
        for i in range(len(names)):
            for j in range(len(names)):
                matrix[i][j] = float(values[i, j])
        source.WriteObject(matrix, key)
    ROOT.TNamed('uncertainty_method', result['method']).Write()
    source.Close()
    # Numerical tables are directly readable by uproot and the comparison table tool
    parameters = dict(result['parameters'])
    if rates:
        parameters.update(rates['parameters'])
    with uproot.update(path) as handle:
        columns = ['median', 'interval_low', 'interval_high', 'error_low', 'error_high',
                   'mean', 'sd', 'rhat', 'ess_bulk', 'ess_tail', 'mode', 'gaussian_mean', 'gaussian_sigma', 'gaussian_amplitude', 'hpd_mass']
        table = handle.mktree('parameter_summary', {'name': 'string', **{key: 'float64' for key in columns}})
        table.extend({'name': list(parameters), **{
            key: np.array([values.get(key) if values.get(key) is not None else np.nan
                           for values in parameters.values()]) for key in columns}})
        region_rows = [(name, low, high) for name, values in parameters.items()
                       for low, high in values.get('hpd_regions', [])]
        if region_rows:
            regions = handle.mktree('highest_density_regions', {'name': 'string', 'low': 'float64', 'high': 'float64'})
            regions.extend({'name': [row[0] for row in region_rows],
                            'low': np.array([row[1] for row in region_rows]),
                            'high': np.array([row[2] for row in region_rows])})
        nominal = {name: config.getfloat(name, 'nom') for name in config.sections()
                   if name != 'summary' and config.has_option(name, 'nom')}
        if rates:
            nominal.update(rates['nominal'])
        table = handle.mktree('nominal', {'name': 'string', 'value': 'float64'})
        table.extend({'name': list(nominal), 'value': np.array(list(nominal.values()))})
        best = dict(result['best_sample']['parameters'])
        if rates:
            best.update(rates['best_fit'])
        table = handle.mktree('best_fit', {'name': 'string', 'value': 'float64'})
        table.extend({'name': list(best), 'value': np.array(list(best.values()))})
        # Store descriptive settings and chain details together inside the ROOT file
        details = {key: value for key, value in result.items()
                   if key not in {'parameters', 'covariance', 'correlation'}}
        details['per_chain'] = {name: values['per_chain'] for name, values in result['parameters'].items()
                                if 'per_chain' in values}
        details['settings'] = settings
        if rates:
            details['reconstruction'] = {key: value for key, value in rates.items()
                                         if key not in {'parameters', 'nominal', 'best_fit'}}
        handle['analysis'] = json.dumps(details, allow_nan=False)
        if derived:
            fields = {name: 'float64' for name in derived}
            fields.update(sample_index='int64', chain_index='int32', Step='int64')
            handle.mktree('derived_samples', fields).extend(derived)
    return path


