"""Find, validate and read MCMC chains using an explicit burn-in cut"""

import json
from pathlib import Path

import numpy as np

from mcmc_file_utils import fingerprint, load_config, root_open

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
    for chain_number, record in enumerate(records, start=1):
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
        if chain_number % 10 == 0 or chain_number == len(records):
            print(f'Read {chain_number}/{len(records)} chains', flush=True)
    return arrays, manifest



def prepare_inputs(args):
    inputs = [str(path.expanduser().resolve()) for path in args.inputs]
    first = Path(inputs[0])
    directory = args.output_dir.expanduser().resolve() if args.output_dir else (first.parent if first.is_file() else first) / 'output'
    names, records, warnings = describe_chains(discover(inputs), args.burn_in)
    manifest = {'schema': 1, 'names': names, 'chains': records, 'burn_in': args.burn_in, 'inputs': inputs}
    manifest['analysis_id'] = fingerprint(manifest)
    return directory, manifest, warnings



def validate_sources(manifest):
    for record in manifest['chains']:
        path = Path(record['file'])
        stat = path.stat()
        if stat.st_size != record['size'] or stat.st_mtime_ns != record['mtime_ns']:
            raise ValueError('Original chain changed: ' + str(path))
        for name, contents in record['configs'].items():
            if (path.parent / 'cfg' / name).read_text() != contents:
                raise ValueError('Original chain config changed: ' + str(path.parent / 'cfg' / name))

