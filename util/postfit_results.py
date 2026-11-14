"""Read single-fit products and share postfit summary calculations"""

import configparser
import csv
import json
import math
from pathlib import Path
import os
import re

def normalisation_rate_error(rate, normalisation, normalisation_error):
    if not math.isfinite(normalisation) or normalisation == 0:
        raise ValueError("Cannot infer a component rate error from a zero or nonfinite normalisation")
    return abs(rate * normalisation_error / normalisation)


def combine_results_normstyle(integrals, fit_params, geo_corr, alphan):
    postfit_central = {}
    postfit_uncert = {}

    # Oscillation and systematics (direct from fit)
    oscillation_params = {"deltam21", "sinsqtheta12"}
    systematics = {}
    if not alphan:
        systematics = {
        "energy_scale", "energy_conv", "birks_constant", "p_recoil_energy_scale",
        "energy_scale2", "energy_conv2", "birks_constant2", "p_recoil_energy_scale2"
        }
    else:
        systematics = {
        "energy_scale", "energy_conv", "birks_constant", "p_recoil_energy_scale",
        "energy_scale2", "energy_conv2", "birks_constant2", "p_recoil_energy_scale2",
        "class_a_ppo", "class_a_bismsb", "class_s_ppo", "class_s_bismsb"
        }
    for p, (val, err) in fit_params.items():
        if p in oscillation_params or p in systematics:
            postfit_central[p] = val
            postfit_uncert[p] = err

    # Normalisations
    base_to_parts = {}
    for stem, integral in integrals.items():
        base = re.sub(r"(_2p2ppo|_bismsb)$", "", stem)
        base = re.sub(r"2$", "", base)
        norm_param = base + "_norm"
        if norm_param not in fit_params:
            continue
        norm_val, norm_err = fit_params[norm_param]
        part = "ppo" if stem.endswith("_2p2ppo") else "bismsb"
        postfit_central[stem] = integral
        postfit_uncert[stem] = normalisation_rate_error(integral, norm_val, norm_err)
        base_to_parts.setdefault(base, {})
        base_to_parts[base][part] = (integral, norm_val, norm_err)

    # Totals
    for base, parts in base_to_parts.items():
        if "ppo" in parts and "bismsb" in parts:
            (ppo_val, norm_val, norm_err) = parts["ppo"]
            (bismsb_val, _, _) = parts["bismsb"]
            total_val = ppo_val + bismsb_val
            total_unc = normalisation_rate_error(total_val, norm_val, norm_err)
            postfit_central[base + "_total"] = total_val
            postfit_uncert[base + "_total"] = total_unc

    # Geo ratios (shared corr)
    corr = geo_corr.get("shared", 0.0)
    add_geo_ratios(postfit_central, postfit_uncert, corr)
    return postfit_central, postfit_uncert


# --------------------------------------------------------------
# Combine results (direct-style)
# --------------------------------------------------------------
def combine_results_directstyle(integrals, fit_params, geo_corr, alphan):
    postfit_central = {}
    postfit_uncert = {}

    oscillation_params = {"deltam21", "sinsqtheta12"}
    systematics = {}
    if not alphan:
        systematics = {
        "energy_scale", "energy_conv", "birks_constant", "p_recoil_energy_scale",
        "energy_scale2", "energy_conv2", "birks_constant2", "p_recoil_energy_scale2"
        }
    else:
        systematics = {
        "energy_scale", "energy_conv", "birks_constant", "p_recoil_energy_scale",
        "energy_scale2", "energy_conv2", "birks_constant2", "p_recoil_energy_scale2",
        "class_a_ppo", "class_a_bismsb", "class_s_ppo", "class_s_bismsb"
        }

    for p, (val, err) in fit_params.items():
        if p in oscillation_params or p in systematics:
            postfit_central[p] = val
            postfit_uncert[p] = err

    # Normalisation-like parameters (use directly)
    base_to_parts = {}
    for stem, integral in integrals.items():
        base = re.sub(r"(_2p2ppo|_bismsb)$", "", stem)
        if base not in fit_params:
            continue
        fit_val, fit_err = fit_params[base]
        part = "ppo" if stem.endswith("_2p2ppo") else "bismsb"
        postfit_central[stem] = integral
        postfit_uncert[stem] = integral * fit_err / fit_val
        base = re.sub(r"2$", "", base)
        base_to_parts.setdefault(base, {})
        base_to_parts[base][part] = (integral, fit_err)

    # Compute totals
    for base, parts in base_to_parts.items():
        if "ppo" in parts and "bismsb" in parts:
            val_ppo, err_ppo = parts["ppo"]
            val_bi, err_bi = parts["bismsb"]
            total_val = val_ppo + val_bi
            total_unc = err_ppo+err_bi
            postfit_central[base + "_total"] = total_val
            postfit_uncert[base + "_total"] = total_unc

    # Geo ratios (separate corr for ppo/bismsb)
    ppo_corr = geo_corr.get("ppo", 0.0)
    bi_corr = geo_corr.get("bismsb", 0.0)
    add_geo_ratios(postfit_central, postfit_uncert, ppo_corr, bi_corr)

    return postfit_central, postfit_uncert


# --------------------------------------------------------------
# Geo ratio helper
# --------------------------------------------------------------
def add_geo_ratios(central, uncert, corr_ppo, corr_bi=None):
    def safe_get(name):
        return central[name], uncert[name]

    try:
        # PPO
        geo_u_ppo, geo_u_ppo_unc = safe_get("geonu_U_2p2ppo") if "geonu_U_2p2ppo" in central else safe_get("geonu_U")
        geo_th_ppo, geo_th_ppo_unc = safe_get("geonu_Th_2p2ppo") if "geonu_Th_2p2ppo" in central else safe_get("geonu_Th")
        # bisMSB
        geo_u_bi, geo_u_bi_unc = safe_get("geonu_U2_bismsb") if "geonu_U2_bismsb" in central else safe_get("geonu_U2")
        geo_th_bi, geo_th_bi_unc = safe_get("geonu_Th2_bismsb") if "geonu_Th2_bismsb" in central else safe_get("geonu_Th2")

        def ratio_and_unc(n, n_unc, d, d_unc, corr):
            r = n / d
            rel = math.sqrt((n_unc / n)**2 + (d_unc / d)**2 - 2*(n_unc / n)*(d_unc / d)*corr)
            return r, r * rel

        r_ppo, r_ppo_unc = ratio_and_unc(geo_u_ppo, geo_u_ppo_unc, geo_th_ppo, geo_th_ppo_unc, corr_ppo)
        r_bi, r_bi_unc = ratio_and_unc(geo_u_bi, geo_u_bi_unc, geo_th_bi, geo_th_bi_unc, corr_bi if corr_bi else corr_ppo)

        # Totals
        geo_u_total = geo_u_ppo + geo_u_bi
        geo_th_total = geo_th_ppo + geo_th_bi
        if corr_bi is None:
            geo_u_total_unc = uncert["geonu_U_total"]
            geo_th_total_unc = uncert["geonu_Th_total"]
        else:
            geo_u_total_unc = math.hypot(geo_u_ppo_unc, geo_u_bi_unc)
            geo_th_total_unc = math.hypot(geo_th_ppo_unc, geo_th_bi_unc)
        r_tot, r_tot_unc = ratio_and_unc(geo_u_total, geo_u_total_unc, geo_th_total, geo_th_total_unc, corr_ppo)

        central["geo_ratio_ppo"], uncert["geo_ratio_ppo"] = r_ppo, r_ppo_unc
        central["geo_ratio_bismsb"], uncert["geo_ratio_bismsb"] = r_bi, r_bi_unc
        central["geo_ratio_total"], uncert["geo_ratio_total"] = r_tot, r_tot_unc
    except KeyError:
        pass


# --------------------------------------------------------------
# CSV writer (same as before)
# --------------------------------------------------------------
def write_csv(outdir, central, uncert, mode, alphan):

    ordered_labels = []

    if not alphan:
        ordered_labels = [
            ("deltam21", "Deltam"),
            ("sinsqtheta12", "Theta"),
            ("reactor_nubar_2p2ppo", "Reac osc PPO"),
            ("reactor_nubar2_bismsb", "Reac osc bisMSB"),
            ("reactor_nubar_total", "Reac osc total"),
            ("geonu_U_2p2ppo", "Geo U PPO"), ("geonu_U2_bismsb", "Geo U bisMSB"), ("geonu_U_total", "Geo U total"),
            ("geonu_Th_2p2ppo", "Geo Th PPO"), ("geonu_Th2_bismsb", "Geo Th bisMSB"), ("geonu_Th_total", "Geo Th total"),
            ("geo_ratio_ppo", "Geo Ratio PPO"), ("geo_ratio_bismsb", "Geo Ratio bisMSB"), ("geo_ratio_total", "Geo Ratio total"),
            ("alphan_CScatter_2p2ppo", "AN CS PPO"), ("alphan_CScatter2_bismsb", "AN CS bisMSB"), ("alphan_CScatter_total", "AN CS total"),
            ("alphan_OExcited_2p2ppo", "AN OE PPO"), ("alphan_OExcited2_bismsb", "AN OE bisMSB"), ("alphan_OExcited_total", "AN OE total"),
            ("alphan_PRecoil_2p2ppo", "AN PR PPO"), ("alphan_PRecoil2_bismsb", "AN PR bisMSB"), ("alphan_PRecoil_total", "AN PR total"),
            ("bipolike_2p2ppo", "BiPo like PPO"), ("bipolike2_bismsb", "BiPo like bisMSB"), ("bipolike_total", "BiPo like total"),
            ("atmospheric_2p2ppo", "Atmos PPO"), ("atmospheric2_bismsb", "Atmos bisMSB"), ("atmospheric_total", "Atmos total"),
            ("energy_scale", "E Scale PPO"), ("energy_scale2", "E Scale bisMSB"),
            ("birks_constant", "Birk's PPO"), ("birks_constant2", "Birk's bisMSB"),
            ("energy_conv", "E Conv PPO"), ("energy_conv2", "E Conv bisMSB"),
            ("p_recoil_energy_scale", "E Scale PR PPO"), ("p_recoil_energy_scale2", "E Scale PR bisMSB")
        ]
    else:
        ordered_labels = [
            ("deltam21", "Deltam"),
            ("sinsqtheta12", "Theta"),
            ("reactor_nubar_2p2ppo", "Reac osc PPO"),
            ("reactor_nubar2_bismsb", "Reac osc bisMSB"),
            ("reactor_nubar_total", "Reac osc total"),
            ("geonu_U_2p2ppo", "Geo U PPO"), ("geonu_U2_bismsb", "Geo U bisMSB"), ("geonu_U_total", "Geo U total"),
            ("geonu_Th_2p2ppo", "Geo Th PPO"), ("geonu_Th2_bismsb", "Geo Th bisMSB"), ("geonu_Th_total", "Geo Th total"),
            ("geo_ratio_ppo", "Geo Ratio PPO"), ("geo_ratio_bismsb", "Geo Ratio bisMSB"), ("geo_ratio_total", "Geo Ratio total"),
            ("alphan_CScatter_2p2ppo", "AN CS PPO"), ("alphan_CScatter2_bismsb", "AN CS bisMSB"), ("alphan_CScatter_total", "AN CS total"),
            ("alphan_OExcited_2p2ppo", "AN OE PPO"), ("alphan_OExcited2_bismsb", "AN OE bisMSB"), ("alphan_OExcited_total", "AN OE total"),
            ("alphan_PRecoil_2p2ppo", "AN PR PPO"), ("alphan_PRecoil2_bismsb", "AN PR bisMSB"), ("alphan_PRecoil_total", "AN PR total"),
            ("bipolike_2p2ppo", "BiPo like PPO"), ("bipolike2_bismsb", "BiPo like bisMSB"), ("bipolike_total", "BiPo like total"),
            ("atmospheric_2p2ppo", "Atmos PPO"), ("atmospheric2_bismsb", "Atmos bisMSB"), ("atmospheric_total", "Atmos total"),
            ("energy_scale", "E Scale PPO"), ("energy_scale2", "E Scale bisMSB"),
            ("birks_constant", "Birk's PPO"), ("birks_constant2", "Birk's bisMSB"),
            ("energy_conv", "E Conv PPO"), ("energy_conv2", "E Conv bisMSB"),
            ("p_recoil_energy_scale", "E Scale PR PPO"), ("p_recoil_energy_scale2", "E Scale PR bisMSB"),
            ("class_a_ppo", "Classifier A PPO"), ("class_a_bismsb", "Classifier A bisMSB"),
            ("class_s_ppo", "Classifier S PPO"), ("class_s_bismsb", "Classifier S bisMSB")
        ]

    csv_path = os.path.join(outdir, "postfit_summary.csv")
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Label", "Parameter", "Central", "Uncertainty", "Rel. Uncertainty (%)"])
        for key, label in ordered_labels:
            if key in central:
                c = central[key]
                u = uncert.get(key, float("nan"))
                rel = (abs(u / c) * 100) if (isinstance(c, (int, float)) and c != 0) else float("nan")
                writer.writerow([label, key, f"{c:.6g}", f"{u:.6g}", f"{rel:.2f}"])
    print(f"✅ CSV written to {csv_path}")



def load_config(path):
    config = configparser.ConfigParser(interpolation=None)
    if not config.read(path):
        raise FileNotFoundError(path)
    return config


def split_names(value):
    return [name.strip() for name in value.split(',') if name.strip()]


def read_minuit_result(fit_dir):
    """Read the original result without rerunning the fitter or changing its configs"""
    import ROOT
    import numpy as np

    directory = Path(fit_dir).resolve()
    fit = load_config(directory / 'cfg/fit_config.ini')
    events = load_config(directory / 'cfg/event_config.ini')
    pdf = load_config(directory / 'cfg/pdf_config.ini')
    text = (directory / 'fit_result.txt').read_text()
    valid = re.findall(r'^Fit Valid:\s*([01])\s*$', text, re.M)
    llh = re.findall(r'^LLH:\s*(\S+)', text, re.M)
    if not valid or valid[-1] != '1' or not llh or not math.isfinite(float(llh[-1])):
        raise ValueError(f'{directory}: missing or invalid fit status/NLLH')
    try:
        source = ROOT.TFile.Open(str(directory / 'fit_result.root'))
    except OSError as error:
        raise ValueError(f'{directory}: cannot open fit_result.root') from error
    if not source or source.IsZombie():
        raise ValueError(f'{directory}: cannot open fit_result.root')
    try:
        for key in ('paramNames', 'paramVals', 'paramErr', 'covMatrix', 'allParamNames'):
            if not source.Get(key):
                raise ValueError(f'{directory}: missing {key}')
        names = list(map(str, source.Get('paramNames')))
        values = list(source.Get('paramVals'))
        errors = list(source.Get('paramErr'))
        all_names = list(map(str, source.Get('allParamNames')))
        matrix = source.Get('covMatrix')
        n = len(names)
        if not n or len(set(names)) != n or len(values) != n or len(errors) != n or matrix.GetNrows() != n or matrix.GetNcols() != n:
            raise ValueError(f'{directory}: inconsistent parameter/covariance dimensions')
        covariance = np.array([[matrix[i][j] for j in range(n)] for i in range(n)])
        if (not np.all(np.isfinite(values)) or not np.all(np.isfinite(errors)) or
                not np.all(np.array(errors) > 0) or not np.all(np.isfinite(covariance)) or
                not np.allclose(covariance, covariance.T)):
            raise ValueError(f'{directory}: invalid values, errors or covariance')
        if not np.allclose(np.diag(covariance), np.square(errors), rtol=1e-5, atol=0):
            raise ValueError(f'{directory}: errors disagree with covariance diagonal')
        correlation = covariance / np.outer(errors, errors)
        if np.linalg.eigvalsh(correlation).min() <= 0:
            raise ValueError(f'{directory}: covariance is not positive definite')
    finally:
        source.Close()
    theta = [name for name in ('theta12', 'sintheta12', 'sinsqtheta12') if name in names]
    if len(theta) != 1 or 'deltam21' not in names:
        raise ValueError('Single-fit Minuit post-processing requires both oscillation parameters floated')
    for name in all_names:
        if not fit.has_section(name):
            raise ValueError(f'{directory}: parameter {name} missing from saved fit config')
    datasets = split_names(events.get('summary', 'datasets'))
    if not datasets or set(datasets) - {'2p2ppo', 'bismsb'}:
        raise ValueError('Distribution plotting currently supports the 2p2ppo and bismsb datasets')
    components = {}
    for dataset in datasets:
        active = split_names(events.get(dataset, 'active'))
        if active == ['all']:
            inactive = split_names(events.get(dataset, 'inactive', fallback=''))
            active = [section for section in events.sections()
                      if events.has_option(section, 'ntup_files') and section not in inactive]
        if not active:
            raise ValueError(f'No active components for {dataset}')
        components[dataset] = active
        for name in active:
            if events.getboolean(name, 'oscillated', fallback=False):
                raise ValueError('Floated post-processing expects the reactor oscillation shape, not pre-oscillated PDFs')
    if pdf.get('summary', 'data_axes').strip() != 'energy':
        raise ValueError('Distribution plotting currently requires the energy data axis')
    low = pdf.getint('energy', 'llh_buffer_bin_low', fallback=8)
    high = pdf.getint('energy', 'llh_buffer_bin_high', fallback=20)
    result = dict(directory=directory, fit=fit, events=events, pdf=pdf, names=names,
                  values=values, errors=errors, all_names=all_names, covariance=covariance,
                  theta=theta[0], datasets=datasets, components=components,
                  observed_data=not (fit.getboolean('summary', 'asimov') or fit.getboolean('summary', 'fake_data', fallback=False)),
                  llh=float(llh[-1]), buffer_low=low, buffer_high=high)
    # Check all required spectra before writing any derived products
    for dataset in datasets:
        for name in components[dataset] + ['postfitdist']:
            read_histogram(directory / 'postfit_dists' / f'{name}_{dataset}.root')
        read_histogram(directory / f'data_{dataset}.root')
    return result


def read_histogram(path):
    import ROOT
    try:
        source = ROOT.TFile.Open(str(path))
    except OSError as error:
        raise ValueError(f'Cannot open saved spectrum {path}') from error
    if not source or source.IsZombie():
        raise ValueError(f'Cannot open saved spectrum {path}')
    try:
        hist = source.Get('oxsx_saved')
        if not hist or not hist.InheritsFrom('TH1'):
            raise ValueError(f'Missing oxsx_saved histogram in {path}')
        if hist.GetDimension() == 2:
            hist = hist.ProjectionX()
        if hist.GetDimension() != 1:
            raise ValueError(f'Expected an energy spectrum in {path}')
        hist.SetDirectory(0)
        return hist
    finally:
        source.Close()


def write_minuit_metadata(result, output_dir=None):
    import ROOT
    directory = Path(output_dir).resolve() if output_dir else result['directory']
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / 'postfit.root'
    output = ROOT.TFile(str(path), 'RECREATE')
    if output.IsZombie():
        raise ValueError(f'Cannot write {path}')
    fit = result['fit']
    names = result['names']
    vectors = {
        'param_names': ('string', names),
        'param_fit_values': ('double', result['values']),
        'param_fit_err': ('double', result['errors']),
        'param_nom_values': ('double', [fit.getfloat(name, 'nom') for name in names]),
        'constr_mean_values': ('double', [fit.getfloat(name, 'constraint_mean', fallback=fit.getfloat(name, 'nom')) for name in names]),
        'constr_sigma_values': ('double', [fit.getfloat(name, 'constraint_sigma', fallback=-999) for name in names]),
        'tex_labels': ('string', [fit.get(name, 'tex_label').strip('"') for name in names]),
        'all_param_names': ('string', result['all_names']),
        'all_tex_labels': ('string', [fit.get(name, 'tex_label').strip('"') for name in result['all_names']]),
        'datasets': ('string', result['datasets']),
    }
    for dataset, components in result['components'].items():
        vectors['components_' + dataset] = ('string', components)
    for key, (kind, entries) in vectors.items():
        vector = ROOT.std.vector(kind)()
        for entry in entries:
            vector.push_back(entry)
        output.WriteObject(vector, key)
    matrix = ROOT.TMatrixD(len(names), len(names))
    for i in range(len(names)):
        for j in range(len(names)):
            matrix[i][j] = result['covariance'][i, j]
    output.WriteObject(matrix, 'covMatrix')
    ROOT.TNamed('fit_directory', str(result['directory'])).Write()
    ROOT.TNamed('uncertainty_method', 'Minuit Hessian; event yields propagate normalisations only').Write()
    ROOT.TParameter(int)('observed_data', int(result['observed_data'])).Write()
    ROOT.TParameter(float)('nllh', result['llh']).Write()
    output.Close()
    return path


def write_minuit_summary(result):
    """Use covariance gradients for normalisation-only yield and U/Th errors"""
    import numpy as np
    directory = result['directory']
    parameters = dict(zip(result['names'], result['values']))
    central = dict(parameters)
    errors = dict(zip(result['names'], result['errors']))
    gradients = {}
    totals = {}
    for dataset, components in result['components'].items():
        for component in components:
            stem = component + '_' + dataset
            hist = read_histogram(directory / 'postfit_dists' / (stem + '.root'))
            low = max(1, result['buffer_low'])
            high = min(hist.GetNbinsX(), hist.GetNbinsX() - result['buffer_high'])
            if high < low:
                raise ValueError(f'Invalid summary bin range for {stem}: {low}..{high}')
            count = float(hist.Integral(low, high))
            base = component[:-1] if component.endswith('2') else component
            norm = base + '_norm' if base + '_norm' in parameters else component
            gradient = np.zeros(len(parameters))
            if norm in parameters:
                if parameters[norm] == 0:
                    raise ValueError(f'Cannot infer yield derivative at zero normalisation: {norm}')
                gradient[result['names'].index(norm)] = count / parameters[norm]
            gradients[stem] = gradient
            central[stem] = count
            totals.setdefault(base, []).append(stem)
    for base, parts in totals.items():
        key = base + '_total'
        central[key] = sum(central[part] for part in parts)
        gradients[key] = sum((gradients[part] for part in parts), np.zeros(len(parameters)))
    for suffix in ('_2p2ppo', '2_bismsb', '_total'):
        u, th = 'geonu_U' + suffix, 'geonu_Th' + suffix
        if u in central and th in central:
            if central[th] <= 0:
                raise ValueError(f'Cannot calculate U/Th ratio with nonpositive Th yield: {th}')
            key = {'_2p2ppo': 'geo_ratio_ppo', '2_bismsb': 'geo_ratio_bismsb', '_total': 'geo_ratio_total'}[suffix]
            central[key] = central[u] / central[th]
            gradients[key] = gradients[u] / central[th] - central[u] * gradients[th] / central[th]**2
    for key, gradient in gradients.items():
        variance = float(gradient @ result['covariance'] @ gradient)
        errors[key] = math.sqrt(max(0, variance))
    for name, values in (('central', central), ('uncert', errors)):
        (directory / ('postfit_' + name + '.json')).write_text(json.dumps(values, indent=2, allow_nan=False) + '\n')
    write_csv(str(directory), central, errors, 'norm' if 'reactor_nubar_norm' in parameters else 'direct',
              any(name.startswith('class_') for name in parameters))
    provenance = dict(fit_result=str(directory / 'fit_result.root'), fit_valid=True, nllh=result['llh'],
                      parameter_errors='Joint Minuit Hessian', event_errors='Normalisation covariance only; explicit shape/oscillation yield dependence omitted',
                      summary_bins=dict(low=result['buffer_low'], upper_buffer=result['buffer_high']),
                      configs={path.name: str(path) for path in (directory / 'cfg').glob('*.ini')})
    (directory / 'postfit_summary_metadata.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return central, errors
