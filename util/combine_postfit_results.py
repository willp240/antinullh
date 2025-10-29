import os
import re
import json
import sys
import subprocess
import math
import csv

# --------------------------------------------------------------
# Run get_integrals.py and plotFixedOscParams
# --------------------------------------------------------------
def run_get_integrals(postfit_dir):
    cmd = ["python", "util/get_integrals.py", postfit_dir]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("❌ Error running get_integrals.py:\n", result.stderr)
        sys.exit(1)
    lines = result.stdout.strip().splitlines()
    pattern = re.compile(r"^(\S+)\s+integral\([^)]*\)\s*=\s*([0-9.+-eE]+)")
    integrals = {}
    for line in lines:
        m = pattern.search(line)
        if m:
            fname = m.group(1)
            val = float(m.group(2))
            integrals[os.path.splitext(fname)[0]] = val
    return integrals


def run_plot_fixed_osc_params(fit_result_tree, corr, alphan):
    cmd = ["root", "-l", "-b", "-q", f'plotting/plotFixedOscParams.C("{fit_result_tree}",{corr},{alphan})']
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("❌ Error running plotFixedOscParams:\n", result.stderr)
        sys.exit(1)
    return result.stdout.strip().splitlines()


# --------------------------------------------------------------
# Parse plotFixedOscParams output
# --------------------------------------------------------------
def parse_plot_params(lines):
    """
    Returns:
      params: {param: (val, err)}
      corr: float or dict {"ppo": x, "bismsb": y}
    """
    params = {}
    geo_corr = {}
    current = None
    for line in lines:
        # Geo correlations
        if "Geo correlation" in line:
            parts = line.split(":")
            if len(parts) == 2:
                label = parts[0].strip().split()[0]  # PPO or bisMSB or nothing
                val = float(parts[1].strip())
                if label in ["PPO", "bisMSB"]:
                    geo_corr[label.lower()] = val
                else:
                    geo_corr["shared"] = val

        # Parameters
        if line.startswith("Par:"):
            current = line.split("Par:")[1].strip()
        elif current and line.startswith("Fit:"):
            parts = line.split()
            if len(parts) >= 3:
                value = float(parts[1])
                uncert = float(parts[2])
                params[current] = (value, uncert)
            current = None

    # Determine style
    if any("_norm" in k for k in params):
        mode = "norm"
    else:
        mode = "direct"

    # Default correlations if missing
    if not geo_corr:
        geo_corr = {"shared": 0.0}

    return params, geo_corr, mode


# --------------------------------------------------------------
# Combine results (norm-style)
# --------------------------------------------------------------
def combine_results_normstyle(integrals, fit_params, geo_corr):
    postfit_central = {}
    postfit_uncert = {}

    # Oscillation and systematics (direct from fit)
    oscillation_params = {"deltam21", "sinsqtheta12"}
    systematics = {
        "energy_scale", "energy_conv", "birks_constant", "p_recoil_energy_scale",
        "energy_scale2", "energy_conv2", "birks_constant2", "p_recoil_energy_scale2"
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
        _, norm_err = fit_params[norm_param]
        part = "ppo" if stem.endswith("_2p2ppo") else "bismsb"
        postfit_central[stem] = integral
        postfit_uncert[stem] = integral * norm_err
        base_to_parts.setdefault(base, {})
        base_to_parts[base][part] = (integral, norm_err)

    # Totals
    for base, parts in base_to_parts.items():
        if "ppo" in parts and "bismsb" in parts:
            (ppo_val, norm_err) = parts["ppo"]
            (bismsb_val, _) = parts["bismsb"]
            total_val = ppo_val + bismsb_val
            total_unc = total_val * norm_err
            postfit_central[base + "_total"] = total_val
            postfit_uncert[base + "_total"] = total_unc

    # Geo ratios (shared corr)
    corr = geo_corr.get("shared", 0.0)
    add_geo_ratios(postfit_central, postfit_uncert, corr)
    return postfit_central, postfit_uncert


# --------------------------------------------------------------
# Combine results (direct-style)
# --------------------------------------------------------------
def combine_results_directstyle(integrals, fit_params, geo_corr):
    postfit_central = {}
    postfit_uncert = {}

    oscillation_params = {"deltam21", "sinsqtheta12"}
    systematics = {
        "energy_scale", "energy_conv", "birks_constant", "p_recoil_energy_scale",
        "energy_scale2", "energy_conv2", "birks_constant2", "p_recoil_energy_scale2"
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
            rel = math.sqrt((n_unc / n)**2 + (d_unc / d)**2 - (n_unc / n)*(d_unc / d)*corr)
            return r, r * rel

        r_ppo, r_ppo_unc = ratio_and_unc(geo_u_ppo, geo_u_ppo_unc, geo_th_ppo, geo_th_ppo_unc, corr_ppo)
        r_bi, r_bi_unc = ratio_and_unc(geo_u_bi, geo_u_bi_unc, geo_th_bi, geo_th_bi_unc, corr_bi if corr_bi else corr_ppo)

        # Totals
        geo_u_total = geo_u_ppo + geo_u_bi
        geo_th_total = geo_th_ppo + geo_th_bi
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
def write_csv(outdir, central, uncert, mode):

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
            ("class_a_ppo", "Classifier A PPO"), ("class_a_bismsb", "Classifier A bisMSB"),
            ("class_s_ppo", "Classifier S PPO"), ("class_s_bismsb", "Classifier S bisMSB"),
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


# --------------------------------------------------------------
# Main
# --------------------------------------------------------------
def main():
    if len(sys.argv) != 3 and len(sys.argv) != 5:
        print("Usage: python combine_postfit_results.py /path/to/postfit_dists /path/to/fit_result_tree.root correlated-fit-bool (alpha,n)-classifier-bool")
        sys.exit(1)

    postfit_dir, fit_result_tree, corr, alphan = sys.argv[1:5]
    print("Running get_integrals.py...")
    integrals = run_get_integrals(postfit_dir)
    print("Running plotFixedOscParams...")
    lines = run_plot_fixed_osc_params(fit_result_tree, corr, alphan)
    fit_params, geo_corr, mode = parse_plot_params(lines)
    print(f"Detected mode: {mode}-style")
    print(f"Geo correlation(s): {geo_corr}")

    print("Combining results...")
    if mode == "norm":
        central, uncert = combine_results_normstyle(integrals, fit_params, geo_corr)
    else:
        central, uncert = combine_results_directstyle(integrals, fit_params, geo_corr)

    outdir = os.path.dirname(postfit_dir.rstrip("/"))
    with open(os.path.join(outdir, "postfit_central.json"), "w") as f:
        json.dump(central, f, indent=2)
    with open(os.path.join(outdir, "postfit_uncert.json"), "w") as f:
        json.dump(uncert, f, indent=2)

    write_csv(outdir, central, uncert, mode)
    print(f"\n✅ All results written under {outdir}")

if __name__ == "__main__":
    main()
