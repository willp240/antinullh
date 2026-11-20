#!/usr/bin/env python3
"""Draw MCMC plots from saved results without reading chains or rerunning calculations"""

import argparse
from contextlib import ExitStack
import json
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path

import numpy as np

from mcmc_file_utils import (REPO, analysis_directory, load_analysis, read_root_result,
                            read_stage_lock, root_open, run_cli, stage_lock)
from mcmc_plot_utils import (PLOT_GROUPS, PlotOptions, render_parameters,
                             render_rates, render_diagnostics)


def plot(args):
    directory = analysis_directory(args.analysis)
    selected = tuple(args.only or ())
    stages = [stage for stage, groups in PLOT_GROUPS.items() if not selected or set(groups) & set(selected)]
    options = PlotOptions(selected, tuple(args.parameters or ()), args.x_range, args.y_range,
                          args.title, args.x_title, args.y_title)
    destination = args.output_dir.expanduser().resolve() if args.output_dir else None
    with ExitStack() as locks:
        missing = []
        available_names = set()
        products = {}
        for stage in stages:
            path = directory
            filename = {'parameters': 'postfit.root', 'diagnostics': 'diagnostics.root'}[stage]
            if not (path / filename).is_file():
                if selected:
                    missing.append(str(path / filename))
                continue
            locks.enter_context(stage_lock(destination or directory, stage))
            if destination and destination != directory:
                locks.enter_context(read_stage_lock(directory, stage))
            if stage == 'parameters':
                manifest, result = load_analysis(directory)
                products[stage] = (path, manifest, result)
                _, _, counts = read_root_result(path / filename)
                available_names.update(manifest['names'] + ['LogL'] + list(counts['parameters'] if counts else []))
            else:
                with root_open(path / filename) as source:
                    manifest = json.loads(str(source['input_manifest']))
                products[stage] = (path, manifest, None)
                available_names.update(manifest['names'])
        identities = {manifest.get('analysis_id') for _, manifest, _ in products.values() if manifest.get('analysis_id')}
        if len(identities) > 1:
            raise ValueError('Saved stages use different chains or burn-in cuts')
        if missing:
            raise ValueError('Missing saved products: ' + ', '.join(missing))
        if not products:
            raise ValueError('No saved MCMC products found in ' + str(directory))
        unknown = set(options.parameters) - available_names
        if unknown:
            raise ValueError('Unknown saved parameters: ' + ', '.join(sorted(unknown)))
        for stage, (path, manifest, result) in products.items():
            target = destination or path
            if stage == 'parameters':
                render_parameters(REPO, path, manifest, result, options, target)
                _, _, counts = read_root_result(path / 'postfit.root')
                if counts:
                    render_rates(REPO, path, counts, options, target)
            else:
                render_diagnostics(path, manifest, options=options, destination=target)



def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis', type=Path)
    parser.add_argument('--only', nargs='+', choices=[name for groups in PLOT_GROUPS.values() for name in groups], help='Draw only these figures')
    parser.add_argument('--parameters', nargs='+', help='Select parameter names for summaries, marginals, traces and autocorrelations')
    parser.add_argument('--x-range', nargs=2, type=float, metavar=('MIN', 'MAX'))
    parser.add_argument('--y-range', nargs=2, type=float, metavar=('MIN', 'MAX'))
    parser.add_argument('--title')
    parser.add_argument('--x-title')
    parser.add_argument('--y-title')
    parser.add_argument('--output-dir', type=Path, help='Write plots elsewhere without changing the saved analysis')
    args = parser.parse_args(argv)
    for name in ('x_range', 'y_range'):
        bounds = getattr(args, name)
        if bounds and (not np.all(np.isfinite(bounds)) or bounds[0] >= bounds[1]):
            parser.error(name.replace('_', '-') + ' requires finite MIN < MAX')
    os.chdir(REPO)
    plot(args)


if __name__ == '__main__':
    run_cli(main, 'MCMC plotting')
