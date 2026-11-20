"""Draw cached MCMC products as editable ROOT canvases and matching PDFs

All display defaults and drawing routines live here; this module never reads chains
"""

from array import array
from dataclasses import dataclass
from functools import lru_cache
import json
import os
import shutil
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from mcmc_file_utils import load_config, root_open


PLOT_GROUPS = {
    'parameters': ('params', 'violins', 'posteriors', 'oscillation', 'pairs', 'correlations', 'traces', 'acceptance', 'spectra', 'rate-posteriors'),
    'diagnostics': ('autocorrelations',),
}


@dataclass
class PlotOptions:
    only: tuple = ()
    parameters: tuple = ()
    x_range: tuple = None
    y_range: tuple = None
    title: str = None
    x_title: str = None
    y_title: str = None

    def wants(self, name):
        return not self.only or name in self.only

    def names(self, names):
        return [name for name in names if not self.parameters or name in self.parameters]


@lru_cache(maxsize=1)
def root():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetErrorX(0)
    ROOT.gStyle.SetPaintTextFormat('.2f')
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.TH1.AddDirectory(False)
    return ROOT


def label(config, name):
    if name == 'LogL':
        return 'Log Likelihood'
    if config.has_section(name):
        key = name
    elif name.startswith('geo_ratio_'):
        dataset = {'ppo': 'PPO', 'bismsb': 'bis-MSB', 'total': 'Total'}[name.removeprefix('geo_ratio_')]
        return 'U/Th Ratio (' + dataset + ')'
    elif name.endswith('_total'):
        base = name.removesuffix('_total')
        key = base + '_norm' if config.has_section(base + '_norm') else base
    else:
        key = name.rsplit('_', 1)[0]
    return config.get(key, 'tex_label', fallback=name.replace('_', ' ')).strip('" ')


def canvas(name, title='', categorical=False):
    R = root()
    c = R.TCanvas(name, title, 1500 if categorical else 1100, 750 if categorical else 650)
    c.SetLeftMargin(.09)
    c.SetRightMargin(.04)
    c.SetBottomMargin(.29 if categorical else .13)
    c.SetTopMargin(.16 if categorical else .10)
    c.SetGridy(True)
    c._keep = []
    return c


def keep(c, obj):
    c._keep.append(obj)
    return obj


def graph(c, name, x, y, color, width=1):
    R = root()
    g = keep(c, R.TGraph(len(x), array('d', x), array('d', y)) if len(x) else R.TGraph())
    g.SetName(name)
    g.SetLineColor(color)
    g.SetMarkerColor(color)
    g.SetLineWidth(width)
    return g


def line(c, name, x1, y1, x2, y2, color, style=1, width=1):
    g = graph(c, name, [x1, x2], [y1, y2], color, width)
    g.SetLineStyle(style)
    g.Draw('L SAME')
    return g


def histogram(c, name, counts, edges):
    R = root()
    h = keep(c, R.TH1D(name, '', len(counts), array('d', edges)))
    for i, count in enumerate(counts):
        h.SetBinContent(i + 1, float(count))
    return h


def axes(h, options, x_title, y_title, default_y=None, title=''):
    h.SetTitle(options.title if options.title is not None else title)
    h.GetXaxis().SetTitle(options.x_title if options.x_title is not None else x_title)
    h.GetYaxis().SetTitle(options.y_title if options.y_title is not None else y_title)
    h.GetYaxis().SetTitleOffset(1.15)
    if options.x_range:
        h.GetXaxis().SetRangeUser(*options.x_range)
    ylim = options.y_range or default_y
    if ylim:
        h.SetMinimum(ylim[0])
        h.SetMaximum(ylim[1])


def legend(c, entries, columns=1, x2=.96):
    R = root()
    leg = keep(c, R.TLegend(.10, .86 if columns > 1 else .67, x2, .98 if columns > 1 else .89))
    leg.SetNColumns(columns)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(.025)
    for obj, text, style in entries:
        leg.AddEntry(obj, text, style)
    leg.Draw()


def write_figures(plots, name, figures, root_name=None):
    """Write all pages to a ROOT file and PDF from the very same canvases"""
    R = root()
    plots.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.draw-', dir=plots) as temporary:
        pdf = Path(temporary) / (name + '.pdf')
        root_path = Path(temporary) / ((root_name or name) + '.root')
        existing = plots / root_path.name
        if root_name and existing.exists():
            shutil.copy2(existing, root_path)
        target = R.TFile(str(root_path), 'UPDATE' if root_path.exists() else 'RECREATE')
        opened = False
        try:
            for c in figures:
                c.Modified()
                c.Update()
                target.cd()
                c.Write(c.GetName(), R.TObject.kOverwrite)
                for obj in getattr(c, '_export', []):
                    obj.Write(obj.GetName(), R.TObject.kOverwrite)
                for folder, objects in getattr(c, '_directories', {}).items():
                    target.Delete(folder + ';*')
                    directory = target.mkdir(folder)
                    directory.cd()
                    for key, obj in objects.items():
                        obj.Write(key)
                    target.cd()
                if not opened:
                    c.Print(str(pdf) + '[')
                    opened = True
                c.Print(str(pdf))
                c.Close()
            if opened:
                closer = R.TCanvas('pdf_close', '', 10, 10)
                closer.Print(str(pdf) + ']')
                closer.Close()
        finally:
            target.Close()
        if opened:
            os.replace(root_path, plots / root_path.name)
            os.replace(pdf, plots / pdf.name)
            print('Saved ' + str(plots / root_path.name), flush=True)


def summary_plot(cache, config, manifest, result, options, violins=False):
    R = root()
    names = options.names(manifest['names'])
    if not names:
        return
    c = canvas('violins' if violins else 'params', categorical=True)
    frame = histogram(c, 'axis', np.zeros(len(names)), np.arange(len(names) + 1) - .5)
    for i, name in enumerate(names):
        frame.GetXaxis().SetBinLabel(i + 1, label(config, name))
    frame.GetXaxis().LabelsOption('v')
    frame.GetXaxis().SetLabelSize(.028)
    default_y = (0., 3.) if violins else (0., 2.)
    axes(frame, options, '', 'Relative to Nominal', default_y)
    frame.Draw('AXIS')
    ylim = options.y_range or default_y
    prefit = keep(c, R.TH1D('constraints', '', len(names), -.5, len(names) - .5))
    nominal = keep(c, prefit.Clone('nominal'))
    postfit = keep(c, prefit.Clone('postfit'))
    modes = keep(c, R.TGraphAsymmErrors())
    modes.SetName('posterior_intervals')
    regions = keep(c, R.TGraphAsymmErrors())
    regions.SetName('posterior_mode_hpd')
    alternatives = {}
    for key, color, marker in [('median', R.kRed, 24), ('mean', R.kOrange + 7, 21),
                                ('gaussian', R.kGreen + 2, 22), ('best', R.kMagenta + 1, 33)]:
        g = keep(c, R.TGraphAsymmErrors())
        g.SetName('posterior_' + key)
        g.SetLineColor(color)
        g.SetMarkerColor(color)
        g.SetMarkerStyle(marker)
        alternatives[key] = g
    for i, name in enumerate(names):
        nom = config.getfloat(name, 'nom')
        scale, shift = (nom, 0.) if nom else (1., 1.)
        transform = lambda value: (value + shift) / scale
        s = result['parameters'][name]
        for h in (prefit, nominal, postfit):
            h.GetXaxis().SetBinLabel(i + 1, label(config, name))
        nominal.SetBinContent(i + 1, 1.)
        postfit.SetBinContent(i + 1, transform(s['mode']))
        sigma = config.getfloat(name, 'constraint_sigma', fallback=-999)
        centre = transform(config.getfloat(name, 'constraint_mean', fallback=nom))
        if sigma <= 0:
            centre = 1.
            error = max(abs(ylim[0] - centre), abs(ylim[1] - centre))
        else:
            error = sigma / abs(scale)
        prefit.SetBinContent(i + 1, centre)
        prefit.SetBinError(i + 1, error)
        if violins:
            counts, edges = cache['posteriors/' + name].to_numpy()
            density = counts / np.diff(edges)
            widths = .32 * np.repeat(density / density.max(), 2)
            y = transform(np.repeat(edges, 2)[1:-1])
            g = graph(c, 'violin_' + name, np.r_[i - widths, (i + widths)[::-1]],
                      np.r_[y, y[::-1]], R.kAzure + 1)
            # Standard opaque colours and an outline remain visible in native ROOT windows
            g.SetTitle(label(config, name))
            g.SetFillColor(R.kAzure - 9)
            g.SetFillStyle(1001)
            g.SetLineWidth(1)
            g.Draw('FL SAME')
        mode_y = transform(s['mode'])
        modes.SetPoint(i, i, mode_y)
        low, high = next((lo, hi) for lo, hi in s['hpd_regions'] if lo <= s['mode'] <= hi)
        modes.SetPointError(i, 0, 0, (s['mode'] - low) / abs(scale), (high - s['mode']) / abs(scale))
        for low, high in s['hpd_regions']:
            lo, hi = sorted((transform(low), transform(high)))
            j = regions.GetN()
            regions.SetPoint(j, i, (lo + hi) / 2)
            regions.SetPointError(j, 0. if violins else .12, 0. if violins else .12, (hi - lo) / 2, (hi - lo) / 2)
        values = {'median': (s['median'], s['error_low'], s['error_high']),
                  'mean': (s['mean'], s['sd'], s['sd']),
                  'gaussian': (s.get('gaussian_mean'), s.get('gaussian_sigma'), s.get('gaussian_sigma')),
                  'best': (result['best_sample']['parameters'][name], 0., 0.)}
        for key, (value, low, high) in values.items():
            if value is not None and low is not None:
                g = alternatives[key]
                j = g.GetN()
                g.SetPoint(j, i, transform(value))
                g.SetPointError(j, 0, 0, low / abs(scale), high / abs(scale))
    prefit.SetLineColor(R.kBlue)
    prefit.SetLineWidth(1)
    prefit.Draw('E1 SAME')
    line(c, 'nominal_line', -.5, 1., len(names) - .5, 1., R.kGray + 1, 2)
    regions.SetFillColorAlpha(R.kCyan + 2, .22)
    regions.SetLineColor(R.kCyan + 2)
    regions.SetLineWidth(3)
    regions.Draw('Z SAME' if violins else '2 SAME')
    if not violins:
        for g in alternatives.values():
            g.Draw('PZ SAME')
    modes.SetLineColor(R.kCyan + 2)
    modes.SetMarkerColor(R.kCyan + 2)
    modes.SetMarkerStyle(20)
    modes.SetMarkerSize(.8)
    modes.Draw('PZ SAME')
    mode_key = keep(c, modes.Clone('mode_legend'))
    mode_key.SetFillColorAlpha(R.kCyan + 2, .22)
    mode_key.SetFillStyle(1001)
    entries = [(prefit, 'Prefit', 'l'), (mode_key, 'Mode, 68% Highest Density', 'flp' if not violins else 'lp')]
    if violins:
        entries.append((next(obj for obj in c._keep if obj.GetName().startswith('violin_')), 'Posterior Density', 'f'))
    else:
        entries += [(alternatives[key], text, 'lp') for key, text in
                    [('median', 'Median, 68% Equal-Tailed'), ('mean', 'Mean, Standard Deviation'),
                     ('gaussian', 'Gaussian Mean, Fitted Width'), ('best', 'Maximum LLH')]]
    if violins:
        legend(c, entries[:2], columns=2, x2=.72)
        density_legend = keep(c, R.TLegend(.48, .86, .63, .98))
        density_legend.SetBorderSize(0)
        density_legend.SetFillStyle(0)
        density_legend.SetTextSize(.025)
        density_legend.AddEntry(*entries[2])
        density_legend.Draw()
    else:
        legend(c, entries, columns=3)
    postfit.SetLineColor(R.kRed)
    method = keep(c, R.TNamed('uncertainty_method', 'posterior_intervals: mode markers with asymmetric 68% HD errors; posterior_mode_hpd: all disconnected regions'))
    saved_mode = keep(c, modes.Clone('posterior_mode'))
    objects = [prefit, nominal, postfit, modes, saved_mode, regions, method, *alternatives.values()]
    if violins:
        shapes = [obj for obj in c._keep if obj.GetName().startswith('violin_')]
        c._directories = {'violin_shapes': {obj.GetName().removeprefix('violin_'): obj for obj in shapes},
                           'violin_data': {obj.GetName(): obj for obj in [frame, *objects]}}
    else:
        c._export = objects
    yield c


def posterior_canvas(name, counts, edges, summary, title, options, nominal=None, best=None, prior=None):
    R = root()
    c = canvas(name, title)
    c.SetCanvasSize(1500, 690)
    c.SetRightMargin(.29)
    h = histogram(c, 'posterior', counts, edges)
    h.SetFillColorAlpha(R.TColor.GetColor('#4682b4'), .35)
    h.SetLineColor(R.TColor.GetColor('#4682b4'))
    ymax = float(max(counts)) * 1.15
    axes(h, options, title, 'Steps', (0., ymax), title)
    h.Draw('HIST')
    mode_color = R.TColor.GetColor('#008080')
    for j, (low, high) in enumerate(summary['hpd_regions']):
        band = graph(c, 'hpd_region_' + str(j), [low, low, high, high, low], [0, ymax, ymax, 0, 0], mode_color)
        band.SetFillColorAlpha(mode_color, .14)
        band.Draw('F SAME')
    h.Draw('HIST SAME')
    def number(value):
        return format(value, '.2e' if 0 < abs(value) < .01 else '.2f')
    mode_text = f"Mode = {number(summary['mode'])}"
    mode_region = next(((low, high) for low, high in summary['hpd_regions']
                        if low <= summary['mode'] <= high), None)
    if mode_region:
        low, high = mode_region
        mode_text += f" (-{number(summary['mode'] - low)}, +{number(high - summary['mode'])})"
    region_text = '68% Highest Density'
    if len(summary['hpd_regions']) > 1:
        region_text += f" ({len(summary['hpd_regions'])} Regions)"
    entries = [(h, 'Posterior', 'f')]
    for key, value, color, style, text in [
        ('median', summary['median'], R.TColor.GetColor('#b22222'), 2,
         f"#splitline{{Median = {number(summary['median'])} (-{number(summary['error_low'])}, +{number(summary['error_high'])})}}{{68% Equal-Tailed}}"),
        ('mean', summary['mean'], R.TColor.GetColor('#ff8c00'), 4,
         f"Mean = {number(summary['mean'])}; SD = {number(summary['sd'])}"),
        ('mode', summary['mode'], mode_color, 3, f'#splitline{{{mode_text}}}{{{region_text}}}')]:
        g = line(c, key, value, 0, value, ymax, color, style, 2)
        if key == 'mode':
            g.SetFillColorAlpha(mode_color, .14)
        entries.append((g, text, 'fl' if key == 'mode' else 'l'))
    gaussian_color = R.TColor.GetColor('#228b22')
    if summary.get('gaussian_mean') is not None and summary.get('gaussian_sigma', 0) > 0:
        mean, sigma = summary['gaussian_mean'], summary['gaussian_sigma']
        x = np.linspace(edges[0], edges[-1], 400)
        y = summary['gaussian_amplitude'] * np.exp(-.5 * ((x - mean) / sigma) ** 2)
        g = graph(c, 'gaussian_fit', x, y, gaussian_color, 2)
        g.Draw('L SAME')
        entries.append((g, f'Gaus. Fit: #mu = {number(mean)}; #sigma = {number(sigma)}', 'l'))
    else:
        g = graph(c, 'gaussian_fit_unavailable', [], [], gaussian_color, 2)
        entries.append((g, 'Gaus. Fit: Not Available', 'l'))
    for key, value, color, style, text in [
        ('best', best, R.TColor.GetColor('#8b008b'), 2, 'Maximum LLH'),
        ('nominal', nominal, R.kBlack, 3, 'Nominal')]:
        if value is not None:
            entries.append((line(c, key, value, 0, value, ymax, color, style, 2),
                            f'{text} = {number(value)}', 'l'))
    if prior and prior[1] > 0:
        mean, sigma = prior
        x = np.linspace(edges[0], edges[-1], 400)
        y = counts.sum() * np.diff(edges).mean() / (sigma * np.sqrt(2 * np.pi)) * np.exp(-.5 * ((x - mean) / sigma) ** 2)
        g = graph(c, 'gaussian_prior', x, y, R.kGray + 1)
        g.SetLineStyle(2)
        g.Draw('L SAME')
        entries.append((g, 'Gaussian Constraint', 'l'))
    leg = keep(c, R.TLegend(.73, .90 - .055 * (len(entries) + 2), .99, .90))
    leg.SetBorderSize(1)
    leg.SetFillColor(R.kWhite)
    leg.SetTextSize(.023)
    leg.SetMargin(.15)
    for obj, text, style in entries:
        leg.AddEntry(obj, text, style)
    leg.Draw()
    return c


def marginal_figures(cache, config, names, summaries, options, nominal, best):
    for name in options.names(names):
        prior = (config.getfloat(name, 'constraint_mean'), config.getfloat(name, 'constraint_sigma')) if config.has_option(name, 'constraint_mean') else None
        yield posterior_canvas(name, *cache['posteriors/' + name].to_numpy(), summaries[name],
                               label(config, name), options, nominal.get(name), best.get(name), prior)


def line_canvas(name, series, x_title, y_title, options, default_y=None, cut=None, title=''):
    R = root()
    R.gStyle.SetPalette(R.kBird)
    c = canvas(name)
    multi = keep(c, R.TMultiGraph())
    multi.SetName('curves')
    for i, (x, y) in enumerate(series):
        g = graph(c, 'chain_' + str(i), x, y, R.TColor.GetColorPalette((i * 17) % R.TColor.GetNumberOfColors()))
        multi.Add(g, 'L')
    multi.Draw('A')
    axes(multi, options, x_title, y_title, default_y, title)
    c.Update()
    if cut is not None:
        line(c, 'burn_in', cut, c.GetUymin(), cut, c.GetUymax(), R.kBlack, 2, 2)
    return c


def density_canvas(name, counts, xedges, yedges, x_title, y_title, options, contours=True):
    R = root()
    c = canvas(name)
    c.SetRightMargin(.14)
    h = keep(c, R.TH2D('density', '', len(xedges) - 1, array('d', xedges), len(yedges) - 1, array('d', yedges)))
    for i, row in enumerate(counts):
        for j, value in enumerate(row):
            h.SetBinContent(i + 1, j + 1, float(value))
    h.SetTitle(options.title or '')
    h.GetXaxis().SetTitle(options.x_title if options.x_title is not None else x_title)
    h.GetYaxis().SetTitle(options.y_title if options.y_title is not None else y_title)
    h.GetZaxis().SetTitle('Steps' if contours else 'Correlation')
    if options.x_range:
        h.GetXaxis().SetRangeUser(*options.x_range)
    if options.y_range:
        h.GetYaxis().SetRangeUser(*options.y_range)
    if not contours:
        h.SetMinimum(-1)
        h.SetMaximum(1)
    # Save the palette action with the canvas so reopened plots retain their colours
    if contours:
        palette_command = 'gStyle->SetPalette(kBird);'
    else:
        palette_command = ('{ Double_t stops[5] = {0., .25, .5, .75, 1.}; '
                           'Double_t red[5] = {0., .25, 1., 1., .5}; '
                           'Double_t green[5] = {0., .25, 1., .25, 0.}; '
                           'Double_t blue[5] = {.5, 1., 1., .25, 0.}; '
                           'TColor::CreateGradientColorTable(5, stops, red, green, blue, 255); '
                           'gStyle->SetNumberContours(255); }')
    palette = keep(c, R.TExec('set_palette', palette_command))
    h.Draw('AXIS')
    palette.Draw()
    h.Draw('COLZ SAME')
    if contours:
        contour_entries = []
        sorted_counts = np.sort(counts.ravel())[::-1]
        cumulative = np.cumsum(sorted_counts) / sorted_counts.sum()
        for mass in (.68, .95):
            level = sorted_counts[min(np.searchsorted(cumulative, mass), len(sorted_counts) - 1)]
            if 0 < level < counts.max():
                contour = keep(c, h.Clone('contour_' + str(int(mass * 100))))
                contour.SetContour(1, array('d', [float(level)]))
                contour.SetLineColor(R.kBlack if mass == .68 else R.kMagenta + 1)
                contour.SetLineStyle(1 if mass == .68 else 2)
                contour.SetLineWidth(2)
                contour.Draw('CONT3 SAME')
                contour_entries.append((contour, f'{int(mass * 100)}% Posterior', 'l'))
        if contour_entries:
            legend(c, contour_entries, columns=2)
    return c


def render_parameters(repo, output, manifest, result, options=None, destination=None):
    options = options or PlotOptions()
    plots = (destination or output) / 'plots'
    config = load_config(output / 'cfg/fit_config.ini')
    with root_open(output / 'postfit.root') as cache:
        if str(cache['analysis_id']) != manifest['analysis_id']:
            raise ValueError('Saved plot data belongs to another analysis')
        for group, filename in [('params', 'params'), ('violins', 'params_violins')]:
            if options.wants(group):
                write_figures(plots, filename, summary_plot(cache, config, manifest, result, options, group == 'violins'), root_name='params')
                if group == 'violins':
                    (plots / 'params_violins.root').unlink(missing_ok=True)
        if options.wants('posteriors'):
            nominal = {name: config.getfloat(name, 'nom') for name in manifest['names']}
            best = {**result['best_sample']['parameters'], 'LogL': result['best_sample']['log_score']}
            write_figures(plots, 'posteriors', marginal_figures(cache, config, manifest['names'] + ['LogL'], result['parameters'], options, nominal, best))
        if options.wants('traces') or options.wants('acceptance'):
            def series(group, name, y):
                for i in range(len(manifest['chains'])):
                    data = cache[f'{group}/chain_{i}'].arrays(['Step', y or name], library='np')
                    yield data['Step'], data[y or name]
            if options.wants('traces'):
                write_figures(plots, 'traces', (line_canvas(name, series('traces', name, None), 'Step', label(config, name), options,
                                                          cut=manifest['burn_in']) for name in options.names(manifest['names'] + ['LogL'])))
            if options.wants('acceptance'):
                write_figures(plots, 'acceptance', [line_canvas('acceptance', series('acceptance', '', 'fraction'), 'Step', 'Acceptance Fraction', options, (0., 1.))])
        for key, (x, y) in json.loads(str(cache['density_axes'])).items():
            if options.wants('oscillation' if key == 'oscillation_posterior' else 'pairs'):
                write_figures(plots, key, [density_canvas(key, *cache['densities/' + key].to_numpy(), label(config, x), label(config, y), options)])
    if options.wants('correlations'):
        names = options.names(manifest['names'])
        if names:
            indexes = [manifest['names'].index(name) for name in names]
            counts = result['correlation'][np.ix_(indexes, indexes)]
            c = density_canvas('correlations', counts, np.arange(len(names) + 1), np.arange(len(names) + 1), '', '', options, False)
            h = c._keep[0]
            for i, name in enumerate(names):
                h.GetXaxis().SetBinLabel(i + 1, label(config, name))
                h.GetYaxis().SetBinLabel(i + 1, label(config, name))
            c.SetCanvasSize(1500, 1000)
            c.SetBottomMargin(.18)
            c.SetLeftMargin(.18)
            c.SetRightMargin(.16)
            h.GetXaxis().SetLabelSize(.03)
            h.GetYaxis().SetLabelSize(.03)
            h.GetXaxis().SetLabelOffset(.007)
            c.SetGrid(False, False)
            for edge in range(1, len(names)):
                line(c, f'grid_x_{edge}', edge, 0., edge, len(names), root().kGray + 1, 3)
                line(c, f'grid_y_{edge}', 0., edge, len(names), edge, root().kGray + 1, 3)
            write_figures(plots, 'correlations', [c])


def render_diagnostics(output, manifest, config_directory=None, options=None, destination=None):
    options = options or PlotOptions()
    if not options.wants('autocorrelations'):
        return
    config = load_config(output / 'cfg/fit_config.ini')
    with root_open(output / 'diagnostics.root') as cache:
        def figures():
            for name in options.names(manifest['names']):
                series = []
                for i in range(len(manifest['chains'])):
                    key = f'autocorrelations/chain_{i}/{name}'
                    if key in cache:
                        values, edges = cache[key].to_numpy()
                        series.append(((edges[:-1] + edges[1:]) / 2, values))
                if series:
                    yield line_canvas(name, series, 'Lag (Steps)', 'Autocorrelation', options, (-.3, 1.02),
                                      title=label(config, name))
        write_figures((destination or output) / 'plots', 'autocorrelations', figures())


def render_rates(repo, output, rates, options=None, destination=None):
    options = options or PlotOptions()
    plots = (destination or output) / 'plots'
    config = load_config(output / 'cfg/fit_config.ini')
    if options.wants('rate-posteriors') and rates['parameters']:
        with root_open(output / 'postfit.root') as cache:
            write_figures(plots, 'rate_posteriors', marginal_figures(cache, config, list(rates['parameters']), rates['parameters'], options, rates['nominal'], rates['best_fit']))
    if options.wants('spectra'):
        render_spectra(repo, output, options, plots)


def render_spectra(repo, output, options, plots):
    """Reuse the shared spectrum renderer on the selected chain's saved distributions"""
    R = root()
    config = load_config(output / 'cfg/fit_config.ini')
    events = load_config(output / 'cfg/event_config.ini')
    datasets = [name.strip() for name in events.get('summary', 'datasets').split(',')]
    if not set(datasets) <= {'2p2ppo', 'bismsb'}:
        raise ValueError('Spectrum renderer currently supports 2p2ppo and bismsb datasets')
    observed = not (config.getboolean('summary', 'asimov') or config.getboolean('summary', 'fake_data', fallback=False))
    choices = ([0] if len(datasets) == 2 else []) + [1 if name == '2p2ppo' else 2 for name in datasets]
    # The shared macro writes beside its metadata file; use a temporary path for redraws
    with tempfile.TemporaryDirectory(prefix='antinu-spectra-') as temporary:
        work = Path(temporary)
        (work / 'postfit.root').symlink_to(output / 'postfit.root')
        for choice in choices:
            expression = ('plotting/plotFitSpectra.C(' + json.dumps(str(work / 'postfit.root')) + ',' + str(choice)
                          + ',"",' + json.dumps(str(output / 'dists')) + ',' + ('true' if observed else 'false') + ')')
            completed = subprocess.run(['root', '-l', '-b', '-q', expression], cwd=repo, capture_output=True, text=True)
            if completed.returncode:
                raise RuntimeError('Spectrum drawing failed: ' + completed.stdout[-1500:] + completed.stderr[-1500:])
        for path in sorted((work / 'plots').glob('*.root')):
            source = R.TFile.Open(str(path))
            try:
                def figures():
                    for key in source.GetListOfKeys():
                        c = source.Get(key.GetName())
                        if c.InheritsFrom('TCanvas'):
                            c._keep = []
                            def adjust(pad):
                                for obj in pad.GetListOfPrimitives():
                                    if obj.InheritsFrom('TPad'):
                                        adjust(obj)
                                    elif obj.InheritsFrom('TH1'):
                                        title_offset = obj.GetYaxis().GetTitleOffset()
                                        axes(obj, options, obj.GetXaxis().GetTitle(), obj.GetYaxis().GetTitle())
                                        obj.GetYaxis().SetTitleOffset(title_offset)
                                pad.Modified()
                            adjust(c)
                            yield c
                write_figures(plots, path.stem, figures())
            finally:
                source.Close()
