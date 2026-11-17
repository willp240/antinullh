"""Draw MCMC diagnostics and reuse the established ROOT postfit renderers"""

import itertools
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot

from mcmc_results import load_config, root_open


def root_call(repo, expression, log):
    with log.open('a') as handle:
        subprocess.run(['root', '-l', '-b', '-q', expression], cwd=repo, check=True,
                       stdout=handle, stderr=subprocess.STDOUT)


def plot_diagnostics(repo, output, arrays, manifest, result, slices, max_lag):
    plots = output / 'plots'
    plots.mkdir(exist_ok=True)
    acfs = []
    with tempfile.TemporaryDirectory(prefix='antinu-acf-') as temporary:
        with (output / 'autocorrelations.log').open('w') as log:
            for i, record in enumerate(manifest['chains']):
                path = Path(temporary) / f'chain_{i}.root'
                subprocess.run([str(repo / 'bin/auto_corrs'), record['file'], '--burn-in', str(manifest['burn_in']),
                                '--max-lag', str(max_lag), '--output', str(path)],
                               check=True, stdout=log, stderr=subprocess.STDOUT)
                with root_open(path) as source:
                    acfs.append({name.split(';')[0]: histogram.values() for name, histogram in source['Auto_Corrs'].items()})
    with uproot.update(output / 'postfit.root') as handle:
        for i, chain in enumerate(acfs):
            for name, values in chain.items():
                handle[f'autocorrelations/chain_{i}/{name}'] = values, np.arange(len(values) + 1) - .5
    config = load_config(Path(manifest['chains'][0]['file']).parent / 'cfg/fit_config.ini')
    with uproot.update(output / 'postfit.root') as histograms, PdfPages(plots / 'posteriors.pdf') as pdf, PdfPages(plots / 'traces.pdf') as traces, PdfPages(plots / 'autocorrelations.pdf') as correlations:
        for j, name in enumerate(manifest['names'] + ['LogL']):
            summary = result['parameters'][name]
            values = arrays['samples'][:, j] if name != 'LogL' else arrays['log_score']
            fig, ax = plt.subplots(figsize=(9, 5), layout='constrained')
            edges = np.histogram_bin_edges(values, bins=80)
            histograms['posteriors/' + name] = np.histogram(values, bins=edges)
            ax.hist(values, bins=edges, density=True, color='steelblue', alpha=.35, label='Pooled retained draws')
            for i, selected in enumerate(slices):
                ax.hist(values[selected], bins=edges, density=True, histtype='step', lw=.8, label=f'Chain {i}')
            ax.axvspan(summary['interval_low'], summary['interval_high'], color='steelblue', alpha=.15, label='68% equal-tailed interval')
            ax.axvline(summary['median'], color='navy', ls='--', label='Posterior median')
            if name != 'LogL':
                ax.axvline(config.getfloat(name, 'nom'), color='black', ls=':', label='Nominal')
                if config.has_option(name, 'constraint_sigma'):
                    sigma = config.getfloat(name, 'constraint_sigma')
                    if sigma > 0:
                        from scipy.stats import norm
                        x = np.linspace(edges[0], edges[-1], 300)
                        ax.plot(x, norm.pdf(x, config.getfloat(name, 'constraint_mean'), sigma), c='grey', ls='--', label='Gaussian constraint')
            ax.set(xlabel=name, ylabel='Posterior density', title=name)
            ax.legend(fontsize=8)
            pdf.savefig(fig); plt.close(fig)
            fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
            for i, record in enumerate(manifest['chains']):
                preview = arrays['traces'][i]
                ax.plot(preview['Step'], preview[name], lw=.6, label=f'Chain {i}')
            ax.axvline(manifest['burn_in'], color='black', ls='--', label='Burn-in cut')
            ax.set(xlabel='Step within chain', ylabel=name, title=name)
            ax.legend(fontsize=8)
            traces.savefig(fig); plt.close(fig)
            if name == 'LogL':
                continue
            fig, ax = plt.subplots(figsize=(9, 4), layout='constrained')
            for i, chain in enumerate(acfs):
                if name in chain:
                    ax.plot(np.arange(len(chain[name])), chain[name], label=f'Chain {i}', lw=1)
            ax.axhline(0, c='black', lw=.5)
            ax.set(xlabel='Lag (steps)', ylabel='Autocorrelation', title=name, ylim=(-.3, 1.02))
            ax.legend(fontsize=8)
            correlations.savefig(fig); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4), layout='constrained')
    for i, selected in enumerate(slices):
        accepted = arrays['accepted'][selected]
        for_plot = np.array([chunk.mean() for chunk in np.array_split(accepted, min(40, len(accepted)))])
        steps = np.array([chunk.mean() for chunk in np.array_split(arrays['steps'][selected], len(for_plot))])
        ax.plot(steps, for_plot, label=f'Chain {i}')
    ax.set(xlabel='Step', ylabel='Acceptance fraction in blocks', title='Acceptance after burn-in')
    ax.legend(); fig.savefig(plots / 'acceptance.pdf'); plt.close(fig)


def density_plot(output, x, y, xname, yname, filename, effective_draws=None):
    # Short correlated chains do not support the same binning detail as long chains
    bins = min(70, max(12, int(4 * effective_draws ** .25))) if effective_draws else 40
    counts, xedges, yedges = np.histogram2d(x, y, bins=bins)
    sorted_counts = np.sort(counts.ravel())[::-1]
    enclosed = np.cumsum(sorted_counts) / sorted_counts.sum()
    thresholds = [float(sorted_counts[min(np.searchsorted(enclosed, mass), len(sorted_counts) - 1)]) for mass in (.68, .95)]
    fig, ax = plt.subplots(figsize=(8, 6), layout='constrained')
    mesh = ax.pcolormesh(xedges, yedges, counts.T, shading='auto', cmap='Blues')
    levels = sorted(set(value for value in thresholds if value > 0 and value < counts.max()))
    if levels:
        contour = ax.contour((xedges[:-1] + xedges[1:]) / 2, (yedges[:-1] + yedges[1:]) / 2,
                             counts.T, levels=levels, colors='black', linewidths=1)
        labels = {value: ('68%' if value == thresholds[0] else '95%') for value in levels}
        ax.clabel(contour, fmt=labels)
    ax.set(xlabel=xname, ylabel=yname, title='Binned posterior density; 68% and 95% credible regions')
    fig.colorbar(mesh, ax=ax, label='Retained draws')
    fig.savefig(output / 'plots' / filename); plt.close(fig)
    with uproot.update(output / 'postfit.root') as handle:
        handle['densities/' + filename.removesuffix('.pdf')] = counts, xedges, yedges


def render(repo, output, arrays, manifest, result, slices, rates, derived, max_lag, pairs):
    (output / 'plots').mkdir(exist_ok=True)
    plot_diagnostics(repo, output, arrays, manifest, result, slices, max_lag)
    metadata = output / 'postfit.root'
    log = output / 'root_plots.log'
    root_call(repo, 'plotting/plotFitParams.C(' + json.dumps(str(metadata)) + ',true)', log)
    names = manifest['names']
    if 'deltam21' in names and 'sinsqtheta12' in names:
        effective = min(result['parameters'][name]['ess_bulk'] or 1 for name in ('deltam21', 'sinsqtheta12'))
        density_plot(output, arrays['samples'][:, names.index('deltam21')],
                     arrays['samples'][:, names.index('sinsqtheta12')],
                     'deltam21 (eV squared)', 'sinsqtheta12', 'oscillation_posterior.pdf', effective)
    chosen = list(itertools.combinations(names, 2)) if pairs == ['all'] else [tuple(pair.split(':')) for pair in pairs]
    for pair in chosen:
        if len(pair) != 2 or not all(name in names for name in pair):
            raise ValueError(f'Unknown parameter pair: {pair}')
        density_plot(output, arrays['samples'][:, names.index(pair[0])], arrays['samples'][:, names.index(pair[1])],
                     pair[0], pair[1], 'pair_' + pair[0] + '_' + pair[1] + '.pdf',
                     min(result['parameters'][name]['ess_bulk'] or 1 for name in pair))
    if rates and rates['selected_draws']:
        with PdfPages(output / 'plots' / 'yield_posteriors.pdf') as pdf, uproot.update(output / 'postfit.root') as histograms:
            for name, summary in rates['parameters'].items():
                values = derived[name]
                fig, ax = plt.subplots(figsize=(8, 5), layout='constrained')
                edges = np.histogram_bin_edges(values, bins=50)
                histograms['yield_posteriors/' + name] = np.histogram(values, bins=edges)
                ax.hist(values, bins=edges, density=True, histtype='step', lw=1.5)
                ax.axvspan(summary['interval_low'], summary['interval_high'], color='steelblue', alpha=.2, label='68% interval')
                ax.axvline(summary['median'], c='navy', ls='--', label='Median')
                ax.axvline(rates['nominal'][name], c='black', ls=':', label='Nominal')
                ax.set(xlabel=name, ylabel='Posterior density', title=f'{name}: {rates["selected_draws"]} selected draws')
                ax.legend(); pdf.savefig(fig); plt.close(fig)
    if rates:
        fit = load_config(Path(manifest['chains'][0]['file']).parent / 'cfg/fit_config.ini')
        events = load_config(Path(manifest['chains'][0]['file']).parent / 'cfg/event_config.ini')
        datasets = [name.strip() for name in events.get('summary', 'datasets').split(',')]
        data = not (fit.getboolean('summary', 'asimov') or fit.getboolean('summary', 'fake_data', fallback=False))
        supported = set(datasets) <= {'2p2ppo', 'bismsb'}
        if not supported:
            raise ValueError('Spectrum renderer currently supports 2p2ppo and bismsb datasets')
        choices = ([0] if len(datasets) == 2 else []) + [1 if name == '2p2ppo' else 2 for name in datasets]
        bestfit = output / 'bestfit'
        for choice in choices:
            root_call(repo, 'plotting/plotFitSpectra.C(' + json.dumps(str(metadata)) + ',' + str(choice)
                      + ',"",' + json.dumps(str(bestfit)) + ',' + ('true' if data else 'false') + ')', log)
        # The existing renderer puts spectra beside its reconstructed fit directory
        for path in (bestfit / 'plots').glob('*'):
            shutil.copy2(path, output / 'plots' / path.name)
    required = ['params.pdf', 'params.root', 'correlations.pdf', 'posteriors.pdf',
                'traces.pdf', 'autocorrelations.pdf', 'acceptance.pdf']
    if rates:
        required += ['dist_' + suffix + '.pdf' for suffix in
                     (['all'] if len(datasets) == 2 else []) + datasets]
        if rates['selected_draws']:
            required += ['yield_posteriors.pdf']
    for filename in required:
        path = output / 'plots' / filename
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError('Missing expected plot ' + str(path))
    if shutil.which('pdflatex'):
        with tempfile.TemporaryDirectory(prefix='antinu-tables-') as temporary:
            document = Path(temporary) / 'tables_document.tex'
            document.write_text(r'''\documentclass{article}
\usepackage[margin=1.5cm]{geometry}
\usepackage{amsmath}
\usepackage[version=4]{mhchem}
\begin{document}
\small
\input{parameter_table_osc.tex}
\newpage
\input{parameter_table_syst.tex}
\end{document}
''')
            for table in output.glob('parameter_table*.tex'):
                shutil.copy2(table, Path(temporary) / table.name)
            with (output / 'tables_compile.log').open('w') as handle:
                subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', document.name],
                               cwd=temporary, stdout=handle, stderr=subprocess.STDOUT, check=True)
            shutil.copy2(Path(temporary) / 'tables_document.pdf', output / 'tables.pdf')
