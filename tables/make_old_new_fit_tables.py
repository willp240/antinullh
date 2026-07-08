#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


ROWS = [
    ("deltam21", r"$\Delta m_{21}^{2}$ (eV$^2$)", "osc"),
    ("sinsqtheta12", r"sin$^2\theta_{12}$", "osc"),
    ("reactor_nubar_2p2ppo", r"Reactor-$\bar{\nu_{e}}$ (PPO)", "osc"),
    ("reactor_nubar2_bismsb", r"Reactor-$\bar{\nu_{e}}$ (bis-MSB)", "osc"),
    ("reactor_nubar_total", r"Reactor-$\bar{\nu_{e}}$ (Total)", "osc"),
    ("geonu_U_2p2ppo", r"Geo. $\bar{\nu_{e}}$ \ce{U} (PPO)", "osc"),
    ("geonu_U2_bismsb", r"Geo. $\bar{\nu_{e}}$ \ce{U} (bis-MSB)", "osc"),
    ("geonu_U_total", r"Geo. $\bar{\nu_{e}}$ \ce{U} (Total)", "osc"),
    ("geonu_Th_2p2ppo", r"Geo. $\bar{\nu_{e}}$ \ce{Th} (PPO)", "osc"),
    ("geonu_Th2_bismsb", r"Geo. $\bar{\nu_{e}}$ \ce{Th} (bis-MSB)", "osc"),
    ("geonu_Th_total", r"Geo. $\bar{\nu_{e}}$ \ce{Th} (Total)", "osc"),
    ("geo_ratio_ppo", r"\ce{U}/\ce{Th} Ratio (PPO)", "osc"),
    ("geo_ratio_bismsb", r"\ce{U}/\ce{Th} Ratio (bis-MSB)", "osc"),
    ("geo_ratio_total", r"\ce{U}/\ce{Th} Ratio (Total)", "osc"),
    ("alphan_CScatter_2p2ppo", r"($\alpha$, n) \ce{C} S. (PPO)", "osc"),
    ("alphan_CScatter2_bismsb", r"($\alpha$, n) \ce{C} S. (bis-MSB)", "osc"),
    ("alphan_CScatter_total", r"($\alpha$, n) \ce{C} S. (Total)", "osc"),
    ("alphan_OExcited_2p2ppo", r"($\alpha$, n) \ce{O} E. (PPO)", "osc"),
    ("alphan_OExcited2_bismsb", r"($\alpha$, n) \ce{O} E. (bis-MSB)", "osc"),
    ("alphan_OExcited_total", r"($\alpha$, n) \ce{O} E. (Total)", "osc"),
    ("alphan_PRecoil_2p2ppo", r"($\alpha$, n) P. R. (PPO)", "osc"),
    ("alphan_PRecoil2_bismsb", r"($\alpha$, n) P. R. (bis-MSB)", "osc"),
    ("alphan_PRecoil_total", r"($\alpha$, n) P. R. (Total)", "osc"),
    ("bipolike_2p2ppo", r"($\alpha$, p) (PPO)", "osc"),
    ("bipolike2_bismsb", r"($\alpha$, p) (bis-MSB)", "osc"),
    ("bipolike_total", r"($\alpha$, p) (Total)", "osc"),
    ("atmospheric_2p2ppo", r"Atmospheric (PPO)", "osc"),
    ("atmospheric2_bismsb", r"Atmospheric (bis-MSB)", "osc"),
    ("atmospheric_total", r"Atmospheric (Total)", "osc"),
    ("energy_scale", r"Energy Scale (PPO)", "syst"),
    ("energy_scale2", r"Energy Scale (bis-MSB)", "syst"),
    ("birks_constant", r"Birk's Const. (PPO)", "syst"),
    ("birks_constant2", r"Birk's Const. (bis-MSB)", "syst"),
    ("energy_conv", r"Energy Conv. (PPO)", "syst"),
    ("energy_conv2", r"Energy Conv. (bis-MSB)", "syst"),
    ("p_recoil_energy_scale", r"E. Scale P. R. (PPO)", "syst"),
    ("p_recoil_energy_scale2", r"E. Scale P. R. (bis-MSB)", "syst"),
    ("class_a_ppo", r"($\alpha$, n) Class. A (PPO)", "syst"),
    ("class_a_bismsb", r"($\alpha$, n) Class. A (bis-MSB)", "syst"),
    ("class_s_ppo", r"($\alpha$, n) Class. S (PPO)", "syst"),
    ("class_s_bismsb", r"($\alpha$, n) Class. S (bis-MSB)", "syst"),
]

NOMINAL_KEYS = {
    "deltam21": "deltam",
    "sinsqtheta12": "theta",
    "reactor_nubar_2p2ppo": "reacoscppo",
    "reactor_nubar2_bismsb": "reacoscbismsb",
    "reactor_nubar_total": "reacosctotal",
    "geonu_U_2p2ppo": "geouppo",
    "geonu_U2_bismsb": "geoubismsb",
    "geonu_U_total": "geoutotal",
    "geonu_Th_2p2ppo": "geothppo",
    "geonu_Th2_bismsb": "geothbismsb",
    "geonu_Th_total": "geothtotal",
    "geo_ratio_ppo": "georatioppo",
    "geo_ratio_bismsb": "georatiobismsb",
    "geo_ratio_total": "georatiototal",
    "alphan_CScatter_2p2ppo": "ancsppo",
    "alphan_CScatter2_bismsb": "ancsbismsb",
    "alphan_CScatter_total": "ancstotal",
    "alphan_OExcited_2p2ppo": "anoeppo",
    "alphan_OExcited2_bismsb": "anoebismsb",
    "alphan_OExcited_total": "anoetotal",
    "alphan_PRecoil_2p2ppo": "anprppo",
    "alphan_PRecoil2_bismsb": "anprbismsb",
    "alphan_PRecoil_total": "anprtotal",
    "bipolike_2p2ppo": "bipolikeppo",
    "bipolike2_bismsb": "bipolikebismsb",
    "bipolike_total": "bipoliketotal",
    "atmospheric_2p2ppo": "atmosppo",
    "atmospheric2_bismsb": "atmosbismsb",
    "atmospheric_total": "atmostotal",
    "energy_scale": "escaleppo",
    "energy_scale2": "escalebismsb",
    "birks_constant": "birksppo",
    "birks_constant2": "birksbismsb",
    "energy_conv": "econvppo",
    "energy_conv2": "econvbismsb",
    "p_recoil_energy_scale": "escaleprppo",
    "p_recoil_energy_scale2": "escaleprbismsb",
    "class_a_ppo": "classifierappo",
    "class_a_bismsb": "classifierabismsb",
    "class_s_ppo": "classifiersppo",
    "class_s_bismsb": "classifiersbismsb",
}


def clean_label(label):
    return re.sub(r"[^a-z0-9]", "", label.lower())


def find_summary(path):
    path = Path(path)
    if path.is_file():
        return path

    direct = path / "postfit_summary.csv"
    if direct.exists():
        return direct

    log_path = path / "run_postfit_scripts.log"
    if log_path.exists():
        for line in log_path.read_text(errors="ignore").splitlines():
            marker = "CSV written to "
            if marker not in line:
                continue
            candidate = Path(line.split(marker, 1)[1].strip())
            if candidate.exists():
                return candidate

    matches = sorted(path.glob("th*/th*_dm*/postfit_summary.csv"))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(f"No postfit_summary.csv found under {path}")
    raise RuntimeError(
        "Multiple postfit_summary.csv files found under "
        f"{path}. Pass the exact CSV path instead."
    )


def read_summary(path):
    data = {}
    with open(find_summary(path), newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row["Parameter"].strip()
            data[key] = (
                float(row["Central"]),
                float(row["Uncertainty"]),
            )
    return data


def read_nominal_csv(path):
    data = {}
    keyed_rows = {}
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 3 or not row[0].strip():
                continue
            keyed_rows[clean_label(row[0])] = row

    for param_key, nominal_key in NOMINAL_KEYS.items():
        match = next((row for key, row in keyed_rows.items() if nominal_key in key), None)
        if not match:
            continue
        value = float(match[1]) if match[1] else None
        err = float(match[2]) if match[2] else None
        data[param_key] = (value, err)

    return data


def default_nominal_csv(name):
    clean_name = name.lower()
    if "alphan" in clean_name and "noalphan" not in clean_name:
        return Path("tables/asimov_alphan_priors.csv")
    return Path("tables/asimov_priors.csv")


def decimals(key, table_type):
    if table_type == "syst":
        return 3
    if key == "sinsqtheta12":
        return 3
    return 2


def fmt_value(key, value, err, table_type, include_error=True):
    if value is None:
        return ""

    ndp = decimals(key, table_type)
    if key == "deltam21":
        value *= 1e5
        if err is not None:
            err *= 1e5
        if include_error and err is not None:
            return f"$({value:.2f}\\pm{err:.2f})\\times10^{{-5}}$"
        return f"${value:.2f}\\times10^{{-5}}$"
    if include_error and err is not None:
        return f"${value:.{ndp}f}\\pm{err:.{ndp}f}$"
    return f"${value:.{ndp}f}$"


def lookup(data, key):
    return data.get(key, (None, None))


def split_label(label):
    for suffix in ("PPO", "bis-MSB", "Total"):
        tag = f" ({suffix})"
        if label.endswith(tag):
            return label[: -len(tag)], suffix
    return label, None


def first_column(label, key):
    main, sublabel = split_label(label)
    if sublabel:
        return (
            f"\\textbf{{{main}}}",
            f"\\textbf{{{sublabel}}}",
        )
    if key in ("deltam21", "sinsqtheta12"):
        return f"\\multirow{{2}}{{*}}{{{label}}}", ""
    return f"\\multirow{{2}}{{*}}{{\\textbf{{{label}}}}}", ""


def row_lines(key, label, table_type, nominal, original, reprocessed):
    nv, ne = lookup(nominal, key)
    ov, oe = lookup(original, key)
    rv, re = lookup(reprocessed, key)

    nom = fmt_value(key, nv, ne, table_type)
    old = fmt_value(key, ov, oe, table_type)
    new = fmt_value(key, rv, re, table_type)

    col1, col1b = first_column(label, key)
    lines = [
        f"\n        {col1} & \\multirow{{2}}{{*}}{{{nom}}} & "
        f"\\multirow{{2}}{{*}}{{{old}}} & \\multirow{{2}}{{*}}{{{new}}} \\\\"
    ]

    if col1b:
        lines.append(f"        {col1b} & & & \\\\")
    else:
        lines.append("        & & & \\\\")

    lines.append("        \\hline")
    return lines


def build_table(nominal, original, reprocessed, table_type, caption, label):
    header = (
        "\\begin{table}\n"
        "    \\small\n"
        "    \\centering\n"
        f"    \\caption{{{caption}}}\n"
        f"    \\label{{{label}}}\n"
        "    \\begin{tabular}{lccc}\n"
        "        \\hline\\hline\n"
        "        \\textbf{Parameter Name} & \\textbf{Nominal Value} & "
        "\\textbf{Original Fit} & \\textbf{Reprocessed Fit} \\\\\n"
        "        \\hline"
    )
    lines = [header]

    for key, latex_label, kind in ROWS:
        if kind != table_type:
            continue
        if key not in nominal and key not in original and key not in reprocessed:
            continue
        lines.extend(row_lines(key, latex_label, table_type, nominal, original, reprocessed))

    lines.append("\n        \\hline\\hline\n    \\end{tabular}\n\\end{table}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Build old-vs-reprocessed fit comparison LaTeX tables."
    )
    parser.add_argument("original", help="Original fit dir or postfit_summary.csv")
    parser.add_argument("reprocessed", help="Reprocessed fit dir or postfit_summary.csv")
    parser.add_argument(
        "--nominal-dir",
        default=None,
        help="Nominal/asimov dir or postfit_summary.csv. Kept for fallback/debug use.",
    )
    parser.add_argument(
        "--nominal-csv",
        default=None,
        help=(
            "Existing nominal TSV/CSV in the tables/*.csv style. Defaults to "
            "tables/asimov_priors.csv, or tables/asimov_alphan_priors.csv for alphan fits."
        ),
    )
    parser.add_argument("--name", default="fit", help="Short name for captions/labels")
    parser.add_argument("-o", "--output", default=None, help="Write LaTeX to this file")
    args = parser.parse_args()

    if args.nominal_csv:
        nominal = read_nominal_csv(args.nominal_csv)
    elif args.nominal_dir:
        nominal = read_summary(args.nominal_dir)
    else:
        nominal = read_nominal_csv(default_nominal_csv(args.name))
    original = read_summary(args.original)
    reprocessed = read_summary(args.reprocessed)

    name = args.name.replace("_", " ")
    label_stub = args.name.replace(" ", "_").lower()
    tex = "\n\n".join(
        [
            build_table(
                nominal,
                original,
                reprocessed,
                "osc",
                f"Comparison of oscillation and normalisation parameter values for the {name} fit.",
                f"tab:{label_stub}_osc_norm_old_new",
            ),
            build_table(
                nominal,
                original,
                reprocessed,
                "syst",
                f"Comparison of systematic parameter values for the {name} fit.",
                f"tab:{label_stub}_syst_old_new",
            ),
        ]
    )

    if args.output:
        Path(args.output).write_text(tex + "\n")
    else:
        print(tex)


if __name__ == "__main__":
    main()
