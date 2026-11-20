"""Read and write saved MCMC results and manage their output directories"""

from contextlib import contextmanager
import configparser
import csv
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import uproot

REPO = Path(__file__).resolve().parents[1]


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



def comma_list(config, section, key):
    return [value.strip() for value in config.get(section, key, fallback='').split(',') if value.strip()]



def write_root_summary(output, manifest, result, rates, derived, settings, plot_data=None):
    # Import PyROOT only for the existing ROOT metadata contract and renderers
    import ROOT
    ROOT.gROOT.SetBatch(True)
    names = manifest['names']
    config = configparser.ConfigParser(interpolation=None)
    config.read_string(manifest['chains'][0]['configs']['fit_config.ini'])
    path = output / 'postfit.root'
    source = ROOT.TFile(str(path), 'RECREATE')
    for key, values, type_name in [
        ('param_names', names, 'string'),
        ('all_param_names', [name for name in config.sections() if name != 'summary'], 'string'),
        ('all_tex_labels', [config.get(name, 'tex_label', fallback=name).strip('"')
                            for name in config.sections() if name != 'summary'], 'string'),
        ('tex_labels', [config.get(name, 'tex_label', fallback=name).strip('"') for name in names], 'string'),
        ('param_fit_values', [result['best_sample']['parameters'][name] for name in names], 'double'),
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
        ('constr_sigma_values', [config.getfloat(name, 'constraint_sigma', fallback=-999) for name in names], 'double'),
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
        columns = ['best_fit', 'median', 'interval_low', 'interval_high', 'error_low', 'error_high',
                   'mean', 'sd', 'rhat', 'ess_bulk', 'ess_tail', 'mode', 'gaussian_mean', 'gaussian_sigma', 'gaussian_amplitude', 'hpd_mass']
        table = handle.mktree('parameter_summary', {'name': 'string', **{key: 'float64' for key in columns}})
        if parameters:
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
            details['rates'] = {key: value for key, value in rates.items()
                                         if key not in {'parameters', 'nominal', 'best_fit'}}
        handle['analysis'] = json.dumps(details, allow_nan=False)
        handle['input_manifest'] = json.dumps(manifest, allow_nan=False)
        if rates and derived:
            for name in rates['parameters']:
                handle['posteriors/' + name] = np.histogram(derived[name], bins=80)
        if plot_data is not None:
            save_parameter_data(handle, manifest, *plot_data)
    return path



def read_root_result(path):
    """Read saved parameter and count estimates without a separate JSON cache"""
    with root_open(path) as source:
        result = json.loads(str(source['analysis']))
        manifest = json.loads(str(source['input_manifest']))
        rows = source['parameter_summary'].arrays(library='np')
        parameters = {str(name): {key: finite(rows[key][i]) for key in rows if key != 'name'}
                      for i, name in enumerate(rows['name'])}
        for name, value in parameters.items():
            value['hpd_regions'] = []
            if name in result.get('per_chain', {}):
                value['per_chain'] = result['per_chain'][name]
        if 'highest_density_regions' in source:
            regions = source['highest_density_regions'].arrays(library='np')
            for name, low, high in zip(regions['name'], regions['low'], regions['high']):
                parameters[str(name)]['hpd_regions'].append([float(low), float(high)])
        result['parameters'] = parameters
        result['analysis_id'] = result['settings'].get('analysis_id')
        for name, key in [('covariance', 'covMatrix'), ('correlation', 'corrMatrix')]:
            if key in source:
                matrix = source[key]
                result[name] = matrix.member('fElements').reshape(matrix.member('fNrows'), matrix.member('fNcols'))
        rates = None
        if 'rates' in result or 'reconstruction' in result:
            rates = dict(result.get('rates', result.get('reconstruction')))
            parameter_names = set(manifest['names']) | {'LogL'}
            rates['parameters'] = {name: value for name, value in parameters.items() if name not in parameter_names}
            for key, tree in [('nominal', 'nominal'), ('best_fit', 'best_fit')]:
                data = source[tree].arrays(library='np')
                rates[key] = {str(name): float(value) for name, value in zip(data['name'], data['value'])}
    return manifest, result, rates



def save_parameter_data(cache, manifest, arrays, slices, pairs, oscillation_bins):
    names = manifest['names']
    chosen = list(itertools.combinations(names, 2)) if pairs == ['all'] else [tuple(pair.split(':')) for pair in pairs]
    densities = {}
    if {'sinsqtheta12', 'deltam21'} <= set(names):
        densities['oscillation_posterior'] = ('sinsqtheta12', 'deltam21', oscillation_bins)
    for pair in chosen:
        if len(pair) != 2 or not all(name in names for name in pair):
            raise ValueError(f'Unknown parameter pair: {pair}')
        densities['pair_' + pair[0] + '_' + pair[1]] = (*pair, 80)
    cache['analysis_id'] = manifest['analysis_id']
    posterior_cache = cache.mkdir('posteriors')
    density_cache = cache.mkdir('densities')
    trace_cache = cache.mkdir('traces')
    acceptance_cache = cache.mkdir('acceptance')
    for j, name in enumerate(names + ['LogL']):
        values = arrays['samples'][:, j] if name != 'LogL' else arrays['log_score']
        posterior_cache[name] = np.histogram(values, bins=80)
    for key, (x, y, bins) in densities.items():
        density_cache[key] = np.histogram2d(
            arrays['samples'][:, names.index(x)], arrays['samples'][:, names.index(y)], bins=bins)
    cache['density_axes'] = json.dumps({key: [x, y] for key, (x, y, bins) in densities.items()})
    for i, (trace, selected) in enumerate(zip(arrays['traces'], slices)):
        trace_cache.mktree(f'chain_{i}', {name: 'float64' for name in trace}).extend(trace)
        accepted, steps = arrays['accepted'][selected], arrays['steps'][selected]
        blocks = min(40, len(accepted))
        data = {'Step': np.array([part.mean() for part in np.array_split(steps, blocks)]),
                'fraction': np.array([part.mean() for part in np.array_split(accepted, blocks)])}
        acceptance_cache.mktree(f'chain_{i}', {'Step': 'float64', 'fraction': 'float64'}).extend(data)



def write_summary_table(output, result, rates=None):
    parameters = rates['parameters'] if rates is not None else result['parameters']
    fields = ['parameter', 'best_fit', 'interval_low', 'interval_high']
    filename = 'rate_table.csv' if rates is not None else 'parameter_table.csv'
    best = rates['best_fit'] if rates is not None else result['best_sample']['parameters']
    with (output / filename).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for name, values in parameters.items():
            if name != 'LogL':
                writer.writerow({'parameter': name, 'best_fit': best[name],
                                 'interval_low': values['interval_low'],
                                 'interval_high': values['interval_high']})



def fingerprint(manifest):
    identity = {key: manifest[key] for key in ('names', 'chains', 'burn_in')}
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()



def analysis_directory(path):
    path = path.expanduser().resolve()
    return path if any((path / name).is_file() for name in
                       ('postfit.root', 'rates.root', 'diagnostics.root')) else path / 'output'



def load_analysis(directory):
    manifest, result, _ = read_root_result(directory / 'postfit.root')
    if manifest['analysis_id'] != fingerprint(manifest):
        raise ValueError('Invalid MCMC input metadata')
    return manifest, result



def check_output_identity(directory, manifest):
    for filename in ('postfit.root', 'rates.root', 'diagnostics.root'):
        path = directory / filename
        if path.is_file():
            with root_open(path) as source:
                if 'analysis' not in source:
                    continue
                previous = json.loads(str(source['analysis']))['settings'].get('analysis_id')
            if previous != manifest['analysis_id']:
                raise ValueError('Output contains another analysis; use a new --output-dir for changed chains or burn-in')



@contextmanager
def stage_lock(directory, stage, wait=False):
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / ('.' + stage + '.lock')).open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | (0 if wait else fcntl.LOCK_NB))
        except BlockingIOError as error:
            raise RuntimeError('This stage is already running: ' + stage) from error
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)



@contextmanager
def read_stage_lock(directory, stage):
    path = directory / ('.' + stage + '.lock')
    if not path.exists():
        yield
        return
    with path.open('rb') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError('Saved products are being updated: ' + stage) from error
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)



def settings(stage, manifest):
    versions = {}
    for name, path in [('antinu_llh', REPO), ('oxo', Path(os.environ.get('OXO_DIR', '/home/parkerw/Software/oxsx3')))]:
        revision = subprocess.run(['git', '-C', str(path), 'rev-parse', 'HEAD'], capture_output=True, text=True)
        versions[name] = revision.stdout.strip() if revision.returncode == 0 else 'unknown'
    return {'stage': stage, 'analysis_id': manifest['analysis_id'], 'burn_in': manifest['burn_in'],
            'source_versions': versions, 'source_hashes': {
                path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (REPO / 'util').glob('mcmc_*.py')},
            'parameter_intervals': 'Highest-LogL retained point and central 68% interval (16%, 84%)',
            'parameter_plot': 'Mode and all 68% highest-density regions'}



def stage_directories(work, config):
    (work / 'logs').mkdir()
    shutil.copytree(config, work / 'cfg')



def publish_products(work, directory, manifest):
    """Publish one completed calculation, restoring previous files if publication fails"""
    paths = []
    for path in work.iterdir():
        if path.name in {'cfg', 'logs', 'plots'}:
            paths.extend(item for item in path.rglob('*') if item.is_file())
        else:
            paths.append(path)
    with stage_lock(directory, 'publish', wait=True):
        check_output_identity(directory, manifest)
        with tempfile.TemporaryDirectory(prefix='.publish-', dir=directory) as temporary:
            backup = Path(temporary)
            changed = []
            try:
                for path in paths:
                    relative = path.relative_to(work)
                    target, previous = directory / relative, backup / relative
                    target.parent.mkdir(parents=True, exist_ok=True)
                    if target.exists():
                        previous.parent.mkdir(parents=True, exist_ok=True)
                        target.rename(previous)
                    changed.append((target, previous))
                    path.rename(target)
            except BaseException:
                for target, previous in reversed(changed):
                    if target.is_dir():
                        shutil.rmtree(target)
                    else:
                        target.unlink(missing_ok=True)
                    if previous.exists():
                        previous.rename(target)
                raise



def run_cli(main, description):
    """Report command failures consistently without hiding the original error"""
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as error:
        print(f'{description} failed: {error}', file=sys.stderr)
        sys.exit(1)
