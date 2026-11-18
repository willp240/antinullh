"""Draw MCMC diagnostics and reuse the established ROOT postfit renderers"""

import itertools
import json
import re
from pathlib import Path
import shutil
import subprocess
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.legend_handler import HandlerTuple
import numpy as np
import uproot

from mcmc_results import load_config, root_open, posterior_estimates


def root_call(repo, expression, log):
    with log.open('a') as handle:
        subprocess.run(['root', '-l', '-b', '-q', expression], cwd=repo, check=True,
                       stdout=handle, stderr=subprocess.STDOUT)



def parameter_label(config, name):
    """Use config labels for parameters and the corresponding physical event counts"""
    if name == 'LogL':
        label = 'Log Likelihood'
    elif config.has_section(name):
        label = config.get(name, 'tex_label', fallback=name).strip('" ')
    elif name.startswith('geo_ratio_'):
        dataset = {'ppo': 'PPO', 'bismsb': 'bis-MSB', 'total': 'Total'}[name.removeprefix('geo_ratio_')]
        label = 'U/Th Ratio (' + dataset + ')'
    elif name.endswith('_total'):
        base = name.removesuffix('_total')
        if config.has_section(base + '_norm'):
            label = config.get(base + '_norm', 'tex_label', fallback=base + ' Total').strip('" ')
        else:
            label = config.get(base, 'tex_label', fallback=base).strip('" ')
            label = re.sub(r'\s+(?:PPO|bis-MSB)$', '', label) + ' Total'
    else:
        component = name.rsplit('_', 1)[0]
        label = config.get(component, 'tex_label', fallback=name.replace('_', ' ')).strip('" ')
    # Protect mathematical expressions while capitalising the surrounding words
    pattern = r"#[A-Za-z]+(?:\{(?:[^{}]|\{[^{}]*\})*\})*(?:[_^]\{[^{}]*\})*|[A-Za-z]+(?:[_^]\{[^{}]*\})+"
    pieces = []
    start = 0
    for match in re.finditer(pattern, label):
        pieces.append(re.sub(r"\b[a-z][A-Za-z]+", lambda word: word[0][0].upper() + word[0][1:], label[start:match.start()]))
        math = match[0].replace('#', chr(92))
        if math.startswith('sin^'):
            math = chr(92) + math
        pieces.append('$' + math + '$')
        start = match.end()
    pieces.append(re.sub(r"\b[a-z][A-Za-z]+", lambda word: word[0][0].upper() + word[0][1:], label[start:]))
    return ''.join(pieces)


def posterior_plot(values, summary, label, nominal=None, best=None, constraint=None):
    """Draw one pooled distribution with the estimation methods and numerical values"""
    if 'hpd_regions' not in summary:
        summary = {**summary, **posterior_estimates(values)}
    def number(value):
        return format(value, '.2e' if 0 < abs(value) < .01 else '.2f')
    fig, ax = plt.subplots(figsize=(12, 5.5), layout='constrained')
    histogram, edges = np.histogram(values, bins=80)
    ax.stairs(histogram, edges, fill=True, color='steelblue', alpha=.35,
              label='Posterior')
    median_band = ax.axvspan(summary['interval_low'], summary['interval_high'], color='firebrick', alpha=.10)
    median_line = ax.axvline(summary['median'], color='firebrick', ls='--', lw=1.5,
               label=f"Median = {number(summary['median'])} (-{number(summary['error_low'])}, +{number(summary['error_high'])})\n68% Equal-Tailed")
    ax.axvline(summary['mean'], color='darkorange', ls='-.', lw=1.5,
               label=f"Mean = {number(summary['mean'])}; SD = {number(summary['sd'])}")
    regions = summary['hpd_regions']
    mode_region = next(((low, high) for low, high in regions if low <= summary['mode'] <= high), None)
    if mode_region is not None:
        low, high = mode_region
        mode_label = f"Mode = {number(summary['mode'])} (-{number(summary['mode'] - low)}, +{number(high - summary['mode'])})"
    else:
        mode_label = f"Mode = {number(summary['mode'])}"
    mode_label += '\n68% Highest Density'
    if len(regions) > 1:
        mode_label += f' ({len(regions)} Regions)'
    ax.axvline(summary['mode'], color='teal', ls=':', lw=1.5, label=mode_label)
    x = np.linspace(edges[0], edges[-1], 400)
    if summary.get('gaussian_mean') is not None:
        mean, sigma = summary['gaussian_mean'], summary['gaussian_sigma']
        ax.plot(x, summary['gaussian_amplitude'] * np.exp(-.5 * ((x - mean) / sigma) ** 2),
                color='forestgreen', lw=1.5,
                label=f'Gaussian Fit: Mean = {number(mean)}; Sigma = {number(sigma)}')
    else:
        ax.plot([], [], color='forestgreen', label='Gaussian Fit: Not Available')
    if best is not None:
        ax.axvline(best, color='darkmagenta', ls='--', lw=1.5, label=f'Maximum LLH = {number(best)}')
    if nominal is not None:
        ax.axvline(nominal, color='black', ls=':', lw=1.5, label=f'Nominal = {number(nominal)}')
    if constraint is not None and constraint[1] > 0:
        from scipy.stats import norm
        ax.plot(x, norm.pdf(x, *constraint) * len(values) * (edges[1] - edges[0]),
                color='grey', ls='--', label='Gaussian Constraint')
    ax.set(xlabel=label, ylabel='Steps', title=label)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend([(median_band, handle) if handle is median_line else handle for handle in handles], labels,
              handler_map={tuple: HandlerTuple(ndivide=1)},
              loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=9)
    return fig, (histogram, edges)

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
            label = parameter_label(config, name)
            nominal = config.getfloat(name, 'nom') if name != 'LogL' else None
            best = result['best_sample']['parameters'][name] if name != 'LogL' else result['best_sample']['log_score']
            constraint = None
            if name != 'LogL' and config.has_option(name, 'constraint_sigma'):
                constraint = (config.getfloat(name, 'constraint_mean'), config.getfloat(name, 'constraint_sigma'))
            fig, histogram = posterior_plot(values, summary, label, nominal, best, constraint=constraint)
            histograms['posteriors/' + name] = histogram
            pdf.savefig(fig); plt.close(fig)
            fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
            for i, record in enumerate(manifest['chains']):
                preview = arrays['traces'][i]
                ax.plot(preview['Step'], preview[name], lw=.6, label=f'Chain {i}')
            ax.axvline(manifest['burn_in'], color='black', ls='--', label='Burn-In Cut')
            ax.set(xlabel='Step Within Chain', ylabel=label, title=label)
            ax.legend(fontsize=8)
            traces.savefig(fig); plt.close(fig)
            if name == 'LogL':
                continue
            fig, ax = plt.subplots(figsize=(9, 4), layout='constrained')
            for i, chain in enumerate(acfs):
                if name in chain:
                    ax.plot(np.arange(len(chain[name])), chain[name], label=f'Chain {i}', lw=1)
            ax.axhline(0, c='black', lw=.5)
            ax.set(xlabel='Lag (Steps)', ylabel='Autocorrelation', title=label, ylim=(-.3, 1.02))
            ax.legend(fontsize=8)
            correlations.savefig(fig); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4), layout='constrained')
    for i, selected in enumerate(slices):
        accepted = arrays['accepted'][selected]
        for_plot = np.array([chunk.mean() for chunk in np.array_split(accepted, min(40, len(accepted)))])
        steps = np.array([chunk.mean() for chunk in np.array_split(arrays['steps'][selected], len(for_plot))])
        ax.plot(steps, for_plot, label=f'Chain {i}')
    ax.set(xlabel='Step', ylabel='Acceptance Fraction In Blocks', title='Acceptance After Burn-In')
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
        config = load_config(Path(manifest['chains'][0]['file']).parent / 'cfg/fit_config.ini')
        with PdfPages(output / 'plots' / 'rate_posteriors.pdf') as pdf, uproot.update(output / 'postfit.root') as histograms:
            for name, summary in rates['parameters'].items():
                fig, histogram = posterior_plot(derived[name], summary, parameter_label(config, name),
                    rates['nominal'].get(name), rates['best_fit'].get(name))
                histograms['rate_posteriors/' + name] = histogram
                pdf.savefig(fig); plt.close(fig)
        for previous in ['yield_posteriors.pdf', 'event_count_posteriors.pdf']:
            (output / 'plots' / previous).unlink(missing_ok=True)
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
            required += ['rate_posteriors.pdf']
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
