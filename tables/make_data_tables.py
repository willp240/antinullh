#!/usr/bin/env python3
import csv
import sys
import re

def clean_label(label):
    """Normalize label names for flexible lookup."""
    return re.sub(r'[^a-z0-9]', '', label.lower())

def fmt(val, err, ndp=2, sci=False):
    """Return LaTeX string of value ± error with given decimal places."""
    if val == "":
        return ""
    try:
        v = float(val)
        e = float(err) if err not in ("", None, "") else None
    except ValueError:
        return val

    # special scientific format for deltam
    if sci:
        v_scaled = v * 1e5
        e_scaled = e * 1e5 if e else 0
        return f"${v_scaled:.2f}\\pm{e_scaled:.2f}$"

    if e is None:
        return f"${v:.{ndp}f}$"
    return f"${v:.{ndp}f}\\pm{e:.{ndp}f}$"

def read_csv(path):
    """Read TSV file and return dict of label→values."""
    data = {}
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or not row[0].strip():
                continue
            key = clean_label(row[0])
            data[key] = row[1:] + [""] * (16 - len(row))
    return data

rows = [
    ("deltam",  r"$\Delta m_{21}^{2}$", True),
    ("theta",   r"sin$^2\theta_{12}$", False),
    ("reacoscppo", r"Reactor-$\bar{\nu_{e}}$ (PPO)", False),
    ("reacoscbismsb", r"Reactor-$\bar{\nu_{e}}$ (bis-MSB)", False),
    ("reacosctotal", r"Reactor-$\bar{\nu_{e}}$ (Total)", False),
    ("geouppo", r"Geo. $\bar{\nu_{e}}$ \ce{U} (PPO)", False),
    ("geoubismsb", r"Geo. $\bar{\nu_{e}}$ \ce{U} (bis-MSB)", False),
    ("geoutotal", r"Geo. $\bar{\nu_{e}}$ \ce{U} (Total)", False),
    ("geothppo", r"Geo. $\bar{\nu_{e}}$ \ce{Th} (PPO)", False),
    ("geothbismsb", r"Geo. $\bar{\nu_{e}}$ \ce{Th} (bis-MSB)", False),
    ("geothtotal", r"Geo. $\bar{\nu_{e}}$ \ce{Th} (Total)", False),
    ("georatioppo", r"\ce{U}/\ce{Th} Ratio (PPO)", False),
    ("georatiobismsb", r"\ce{U}/\ce{Th} Ratio (bis-MSB)", False),
    ("georatiototal", r"\ce{U}/\ce{Th} Ratio (Total)", False),
]

def build_table(data):
    header = (
        "\\begin{table}\n"
        "    \\small\n"
        "    \\centering\n"
        "    \\caption{Fitted oscillation parameters and signal normalisations for all data fits with correlated normalisations.}\n"
        "    \\resizebox{\\textwidth}{!}{%\n"
        "    \\begin{tabular}{lcccccccc}\n"
        "        \\hline\\hline\n"
        "        \\multirow{3}{*}{\\textbf{Parameter}} & \n"
        "        \\multicolumn{4}{c}{\\textbf{Without ($\\alpha$, n) Classifier}} & \\multicolumn{4}{c}{\\textbf{With ($\\alpha$, n) Classifier}} \\\\\n"
        "        & \\multirow{2}{*}{\\textbf{Asimov}} & \\textbf{No} & \\textbf{sin$^2\\theta_{12}$} & \\textbf{$\\Delta m_{21}^{2} \\&$ sin$^2\\theta_{12}$} & \\multirow{2}{*}{\\textbf{Asimov}} & \\textbf{No} & \\textbf{sin$^2\\theta_{12}$} & \\textbf{$\\Delta m_{21}^{2} \\&$ sin$^2\\theta_{12}$} \\\\\n"
        "        & & \\textbf{Prior} & \\textbf{Prior} & \\textbf{Prior} & & \\textbf{Prior} & \\textbf{Prior} & \\textbf{Prior} \\\\\n"
        "        \\hline"
    )
    lines = [header]

    for key_sub, label, is_deltam in rows:
        k = clean_label(key_sub)
        match = next((v for n, v in data.items() if k in n), None)
        if not match:
            continue

        vals = match[:16]
        if len(vals) < 16:
            vals += [""] * (16 - len(vals))

        # extract 8 pairs
        left = [fmt(vals[i], vals[i+1], ndp=(3 if "theta" in k else 2), sci=is_deltam) for i in range(0, 8, 2)]
        right = [fmt(vals[i], vals[i+1], ndp=(3 if "theta" in k else 2), sci=is_deltam) for i in range(8, 16, 2)]

        # first row
        if is_deltam:
            lines.append(
                f"\n        {label} & \\multirow{{2}}{{*}}{{{left[0]}}} & \\multirow{{2}}{{*}}{{{left[1]}}} & \\multirow{{2}}{{*}}{{{left[2]}}} & \\multirow{{2}}{{*}}{{{left[3]}}} & \\multirow{{2}}{{*}}{{{right[0]}}} & \\multirow{{2}}{{*}}{{{right[1]}}} & \\multirow{{2}}{{*}}{{{right[2]}}} & \\multirow{{2}}{{*}}{{{right[3]}}} \\\\"
            )
            lines.append("        ($\\times10^{-5}$eV$^2$) & & & & & & & & \\\\")
        elif "theta" in k:
            lines.append(
                f"\n        \\multirow{{2}}{{*}}{{{label}}} & \\multirow{{2}}{{*}}{{{left[0]}}} & \\multirow{{2}}{{*}}{{{left[1]}}} & \\multirow{{2}}{{*}}{{{left[2]}}} & \\multirow{{2}}{{*}}{{{left[3]}}} & \\multirow{{2}}{{*}}{{{right[0]}}} & \\multirow{{2}}{{*}}{{{right[1]}}} & \\multirow{{2}}{{*}}{{{right[2]}}} & \\multirow{{2}}{{*}}{{{right[3]}}} \\\\"
            )
            lines.append("        & & & & & & & & \\\\")
        else:
            lines.append(
                f"\n        \\textbf{{{label.split('(')[0].strip()}}} & \\multirow{{2}}{{*}}{{{left[0]}}} & \\multirow{{2}}{{*}}{{{left[1]}}} & \\multirow{{2}}{{*}}{{{left[2]}}} & \\multirow{{2}}{{*}}{{{left[3]}}} & \\multirow{{2}}{{*}}{{{right[0]}}} & \\multirow{{2}}{{*}}{{{right[1]}}} & \\multirow{{2}}{{*}}{{{right[2]}}} & \\multirow{{2}}{{*}}{{{right[3]}}} \\\\"
            )
            sublabel = label.split('(')[1].split(')')[0] if '(' in label else ""
            lines.append(f"        \\textbf{{{sublabel}}} & & & & & & & & \\\\")

        lines.append("        \\hline")

    lines.append(
        "\n        \\hline\\hline\n"
        "    \\end{tabular}%\n"
        "    } % end resizebox\n"
        "    \\label{tab:alldatfitres}\n"
        "\\end{table}"
    )
    return "\n".join(lines)

def main():
    if len(sys.argv) < 2:
        print("Usage: python make_corr_table.py <csvfile>")
        sys.exit(1)

    data = read_csv(sys.argv[1])
    print(build_table(data))

if __name__ == "__main__":
    main()
