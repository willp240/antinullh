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


def row_lines(key, label, table_type, nominal, fit1, fit2, fit3=None):
    nv, ne = lookup(nominal, key)
    value1, error1 = lookup(fit1, key)
    value2, error2 = lookup(fit2, key)

    nom = fmt_value(key, nv, ne, table_type)
    formatted1 = fmt_value(key, value1, error1, table_type)
    formatted2 = fmt_value(key, value2, error2, table_type)

    col1, col1b = first_column(label, key)
    extra = ""
    if fit3 is not None:
        av, ae = lookup(fit3, key)
        extra = f" & \\multirow{{2}}{{*}}{{{fmt_value(key, av, ae, table_type)}}}"
    lines = [
        f"\n        {col1} & \\multirow{{2}}{{*}}{{{nom}}} & "
        f"\\multirow{{2}}{{*}}{{{formatted1}}} & \\multirow{{2}}{{*}}{{{formatted2}}}{extra} \\\\"
    ]

    empty_cells = " &" * (4 if fit3 is not None else 3)
    if col1b:
        lines.append(f"        {col1b}{empty_cells} \\\\")
    else:
        lines.append(f"       {empty_cells} \\\\")

    lines.append("        \\hline")
    return lines


def build_table(nominal, fit1, fit2, table_type, caption, label,
                fit3=None, label1="Fit 1", label2="Fit 2",
                label3="Fit 3", totals_only=False):
    columns = "lcccc" if fit3 is not None else "lccc"
    extra_header = f" & \\textbf{{{label3}}}" if fit3 is not None else ""
    header = (
        "\\begin{table}\n"
        "    \\small\n"
        "    \\centering\n"
        f"    \\caption{{{caption}}}\n"
        f"    \\label{{{label}}}\n"
        f"    \\begin{{tabular}}{{{columns}}}\n"
        "        \\hline\\hline\n"
        "        \\textbf{Parameter Name} & \\textbf{Nominal Value} & "
        f"\\textbf{{{label1}}} & \\textbf{{{label2}}}{extra_header} \\\\\n"
        "        \\hline"
    )
    lines = [header]

    for key, latex_label, kind in ROWS:
        if kind != table_type:
            continue
        if totals_only and key not in ("deltam21", "sinsqtheta12") and not key.endswith("_total"):
            continue
        if key not in nominal and key not in fit1 and key not in fit2 and (fit3 is None or key not in fit3):
            continue
        lines.extend(row_lines(key, latex_label, table_type, nominal, fit1, fit2, fit3))

    lines.append("\n        \\hline\\hline\n    \\end{tabular}\n\\end{table}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Build comparison LaTeX tables for two or three fits."
    )
    parser.add_argument("fit1", help="First fit dir or postfit_summary.csv")
    parser.add_argument("fit2", help="Second fit dir or postfit_summary.csv")
    parser.add_argument("--fit3", help="Optional third fit dir or postfit_summary.csv")
    parser.add_argument("--label1", default="Fit 1")
    parser.add_argument("--label2", default="Fit 2")
    parser.add_argument("--label3", default="Fit 3")
    parser.add_argument("--totals-only", action="store_true", help="Show oscillation parameters and summed event rates")
    parser.add_argument("--table-type", choices=("all", "osc", "syst"), default="all")
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
    fit1 = read_summary(args.fit1)
    fit2 = read_summary(args.fit2)
    fit3 = read_summary(args.fit3) if args.fit3 else None

    name = args.name.replace("_", " ")
    label_stub = args.name.replace(" ", "_").lower()
    tables = []
    for kind, description in (("osc", "oscillation and normalisation"), ("syst", "systematic")):
        if args.table_type not in ("all", kind):
            continue
        caption = f"Comparison of {description} parameter values for the {name} fit."
        label_kind = "osc_norm" if kind == "osc" else "syst"
        tables.append(build_table(
            nominal, fit1, fit2, kind, caption, f"tab:{label_stub}_{label_kind}_comparison",
            fit3, args.label1, args.label2, args.label3, args.totals_only))
    tex = "\n\n".join(tables)

    if args.output:
        Path(args.output).write_text(tex + "\n")
    else:
        print(tex)


if __name__ == "__main__":
    main()
