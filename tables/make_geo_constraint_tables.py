#!/usr/bin/env python3
import argparse
import csv
import math
import re
from pathlib import Path

from make_old_new_fit_tables import read_summary


PAPER_ROWS = [
    ("deltam21", r"$\Delta m_{21}^{2}$ (eV$^2$)"),
    ("sinsqtheta12", r"sin$^2\theta_{12}$"),
    ("reactor_nubar_total", r"Reactor-$\bar{\nu}_{e}$"),
    ("geonu_U_total", r"Geo. $\bar{\nu}_{e}$ \ce{U}"),
    ("geonu_Th_total", r"Geo. $\bar{\nu}_{e}$ \ce{Th}"),
    ("geo_ratio_total", r"\ce{U}/\ce{Th} Ratio"),
    ("alphan_PRecoil_total", r"($\alpha$, n) p-scattering"),
    ("alphan_other_total", r"Other ($\alpha$, n)"),
    ("bipolike_total", r"($\alpha$, p)"),
    ("atmospheric_total", r"Atmospheric"),
]


JAMES_PAIR = {
    "noprior": 1,
    "thprior": 2,
    "thdmprior": 3,
}

LEGACY_KEYS = {
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
}


def clean_label(label):
    return re.sub(r"[^a-z0-9]", "", label.lower())


def read_legacy(path, james_pair):
    keyed_rows = {}
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if row and row[0].strip():
                keyed_rows[clean_label(row[0])] = row

    predicted = {}
    james = {}
    value_index = 1 + 2 * james_pair
    for parameter, legacy_key in LEGACY_KEYS.items():
        row = next(
            (values for key, values in keyed_rows.items() if legacy_key in key),
            None,
        )
        if row is None:
            continue
        predicted[parameter] = (
            float(row[1]) if len(row) > 1 and row[1] else None,
            float(row[2]) if len(row) > 2 and row[2] else None,
        )
        james[parameter] = (
            float(row[value_index]) if len(row) > value_index and row[value_index] else None,
            float(row[value_index + 1])
            if len(row) > value_index + 1 and row[value_index + 1]
            else None,
        )
    return predicted, james


def lookup(values, key):
    return values.get(key, (None, None))


def format_value(key, value, error):
    if value is None:
        return ""

    if key == "deltam21":
        value *= 1e5
        if error is None:
            return f"${value:.2f}\\times10^{{-5}}$"
        return f"$({value:.2f}\\pm{error * 1e5:.2f})\\times10^{{-5}}$"

    if key == "sinsqtheta12":
        decimals = 3
    elif key == "geo_ratio_total":
        decimals = 2
    else:
        decimals = 0

    if error is None:
        return f"${value:.{decimals}f}$"
    return f"${value:.{decimals}f}\\pm{error:.{decimals}f}$"


def add_paper_composites(values):
    cs_value, cs_error = lookup(values, "alphan_CScatter_total")
    oe_value, oe_error = lookup(values, "alphan_OExcited_total")
    if cs_value is not None and oe_value is not None:
        error = None
        if cs_error is not None and oe_error is not None:
            error = math.hypot(cs_error, oe_error)
        values["alphan_other_total"] = (cs_value + oe_value, error)
    return values


def make_row(key, label, predicted, james, ratio, fractional):
    columns = []
    for values in (predicted, james, ratio, fractional):
        value, error = lookup(values, key)
        columns.append(format_value(key, value, error))

    return [
        f"\n        {label} & " + " & ".join(columns) + " \\\\",
    ]


def build_table(predicted, james, ratio, fractional, name):
    display_name = name.replace("_", " ")
    label_name = name.lower().replace(" ", "_")
    lines = [
        "\\begin{table}\n"
        "    \\small\n"
        "    \\centering\n"
        f"    \\caption{{Comparison of the main post-fit parameter values for the {display_name} fit.}}\n"
        f"    \\label{{tab:{label_name}_geo_constraints}}\n"
        "    \\resizebox{\\textwidth}{!}{%\n"
        "    \\begin{tabular}{lcccc}\n"
        "        \\hline\\hline\n"
        "        \\textbf{Parameter Name} & \\textbf{Predicted} & \\textbf{James} & "
        "\\textbf{OXO Ratio Geo-Constraint} & \\textbf{OXO Fractional Geo-Constraint} \\\\\n"
        "        \\hline"
    ]

    for key, label in PAPER_ROWS:
        if not any(key in values for values in (predicted, james, ratio, fractional)):
            continue
        lines.extend(make_row(key, label, predicted, james, ratio, fractional))

    lines.append(
        "\n        \\hline\\hline\n"
        "    \\end{tabular}%\n"
        "    }\n"
        "\\end{table}"
    )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Build predicted/James/OXO geo-constraint comparison table."
    )
    parser.add_argument("legacy_csv", help="Historical James TSV file")
    parser.add_argument("ratio", help="OXO ratio fit directory or postfit_summary.csv")
    parser.add_argument("fractional", help="OXO fractional fit directory or postfit_summary.csv")
    parser.add_argument("--james-fit", choices=JAMES_PAIR, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    predicted, james = read_legacy(args.legacy_csv, JAMES_PAIR[args.james_fit])
    predicted = add_paper_composites(predicted)
    james = add_paper_composites(james)
    ratio = add_paper_composites(read_summary(args.ratio))
    fractional = add_paper_composites(read_summary(args.fractional))
    table = build_table(predicted, james, ratio, fractional, args.name)

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(table + "\n")
    else:
        print(table)


if __name__ == "__main__":
    main()
