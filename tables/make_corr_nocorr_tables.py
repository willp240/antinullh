#!/usr/bin/env python3
import csv
import sys
import re

def clean_label(label):
    """Normalize label names for flexible lookup."""
    return re.sub(r'[^a-z0-9]', '', label.lower())

def fmt(val, err, ndp=2):
    """Return LaTeX string of value ± error with given decimal places."""
    if val == "":
        return ""
    try:
        v = float(val)
        e = float(err) if err not in ("", None, "") else None
    except ValueError:
        return val

    if e is None:
        return f"${v:.{ndp}f}$"
    if ((abs(v) < 1e-4 and v > 0) or abs(v) >= 1e4):
        exp = int(f"{v:e}".split('e')[1])
        mant = v / (10**exp)
        mant_e = e / (10**exp)
        return f"$({mant:.{ndp}f}\\pm{mant_e:.{ndp}f})\\times10^{{{exp}}}$"
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
            data[key] = row[1:] + [""] * (6 - len(row))
    return data

rows = [
    ("deltam", r"$\Delta m_{21}^{2}$ (eV$^2$)", "osc"),
    ("theta", r"sin$^2\theta_{12}$", "osc"),
    ("reacoscppo", r"Reactor-$\bar{\nu_{e}}$ (PPO)", "osc"),
    ("reacoscbismsb", r"Reactor-$\bar{\nu_{e}}$ (bis-MSB)", "osc"),
    ("reacosctotal", r"Reactor-$\bar{\nu_{e}}$ (Total)", "osc"),
    ("geouppo", r"Geo. $\bar{\nu_{e}}$ \ce{U} (PPO)", "osc"),
    ("geoubismsb", r"Geo. $\bar{\nu_{e}}$ \ce{U} (bis-MSB)", "osc"),
    ("geoutotal", r"Geo. $\bar{\nu_{e}}$ \ce{U} (Total)", "osc"),
    ("geothppo", r"Geo. $\bar{\nu_{e}}$ \ce{Th} (PPO)", "osc"),
    ("geothbismsb", r"Geo. $\bar{\nu_{e}}$ \ce{Th} (bis-MSB)", "osc"),
    ("geothtotal", r"Geo. $\bar{\nu_{e}}$ \ce{Th} (Total)", "osc"),
    ("georatioppo", r"\ce{U}/\ce{Th} Ratio (PPO)", "osc"),
    ("georatiobismsb", r"\ce{U}/\ce{Th} Ratio (bis-MSB)", "osc"),
    ("georatiototal", r"\ce{U}/\ce{Th} Ratio (Total)", "osc"),
    ("ancsppo", r"($\alpha$, n) \ce{C} S. (PPO)", "osc"),
    ("ancsbismsb", r"($\alpha$, n) \ce{C} S. (bis-MSB)", "osc"),
    ("ancstotal", r"($\alpha$, n) \ce{C} S. (Total)", "osc"),
    ("anoeppo", r"($\alpha$, n) \ce{O} E. (PPO)", "osc"),
    ("anoebismsb", r"($\alpha$, n) \ce{O} E. (bis-MSB)", "osc"),
    ("anoetotal", r"($\alpha$, n) \ce{O} E. (Total)", "osc"),
    ("anprppo", r"($\alpha$, n) P. R. (PPO)", "osc"),
    ("anprbismsb", r"($\alpha$, n) P. R. (bis-MSB)", "osc"),
    ("anprtotal", r"($\alpha$, n) P. R. (Total)", "osc"),
    ("bipolikeppo", r"($\alpha$, p) (PPO)", "osc"),
    ("bipolikebismsb", r"($\alpha$, p) (bis-MSB)", "osc"),
    ("bipoliketotal", r"($\alpha$, p) (Total)", "osc"),
    ("atmosppo", r"Atmospheric (PPO)", "osc"),
    ("atmosbismsb", r"Atmospheric (bis-MSB)", "osc"),
    ("atmostotal", r"Atmospheric (Total)", "osc"),

    # Systematics
    ("escaleppo", r"Energy Scale (PPO)", "syst"),
    ("escalebismsb", r"Energy Scale (bis-MSB)", "syst"),
    ("birksppo", r"Birk's Const. (PPO)", "syst"),
    ("birksbismsb", r"Birk's Const. (bis-MSB)", "syst"),
    ("econvppo", r"Energy Conv. (PPO)", "syst"),
    ("econvbismsb", r"Energy Conv. (bis-MSB)", "syst"),
    ("escaleprppo", r"E. Scale P. R. (PPO)", "syst"),
    ("escaleprbismsb", r"E. Scale P. R. (bis-MSB)", "syst"),
]

def build_table(data, table_type):
    if table_type == "osc":
        caption = "The fitted oscillation and normalisation parameter values from the Asimov fits with correlated and uncorrelated normalisations."
    else:
        caption = "The fitted systematic parameter values from the Asimov fits with correlated and uncorrelated normalisations."

    header = (
        "\\begin{table}\n"
        "    \\small\n"
        "    \\centering\n"
        f"    \\caption{{{caption}}}\n"
        "    \\begin{tabular}{lccc}\n"
        "        \\hline\\hline\n"
        "        \\textbf{Parameter} & \\textbf{Asimov Value} & \\textbf{Correlated Normalisations} & \\textbf{Uncorrelated Normalisations} \\\\\n"
        "        \\hline"
    )
    lines = [header]

    for key_sub, label, ttype in rows:
        if ttype != table_type:
            continue
        k = clean_label(key_sub)
        match = next((v for n, v in data.items() if k in n), None)
        if not match:
            continue

        asv, ase, corrv, corre, uncorrv, uncorre = match[:6]

        if table_type == "syst":
            ndp = 3
        elif "theta" in k:
            ndp = 3
        else:
            ndp = 2

        av = fmt(asv, ase, ndp)
        c = fmt(corrv, corre, ndp)
        u = fmt(uncorrv, uncorre, ndp)

        # Multirow formatting
        if "(PPO)" in label:
            main = label.replace(" (PPO)", "")
            lines.append(f"\n        \\textbf{{{main}}} & \\multirow{{2}}{{*}}{{{av}}} & \\multirow{{2}}{{*}}{{{c}}} & \\multirow{{2}}{{*}}{{{u}}} \\\\")
            lines.append("        \\textbf{PPO} & & & \\\\")
        elif "(bis-MSB)" in label:
            main = label.replace(" (bis-MSB)", "")
            lines.append(f"\n        \\textbf{{{main}}} & \\multirow{{2}}{{*}}{{{av}}} & \\multirow{{2}}{{*}}{{{c}}} & \\multirow{{2}}{{*}}{{{u}}} \\\\")
            lines.append("        \\textbf{bis-MSB} & & & \\\\")
        elif "(Total)" in label:
            main = label.replace(" (Total)", "")
            lines.append(f"\n        \\textbf{{{main}}} & \\multirow{{2}}{{*}}{{{av}}} & \\multirow{{2}}{{*}}{{{c}}} & \\multirow{{2}}{{*}}{{{u}}} \\\\")
            lines.append("        \\textbf{Total} & & & \\\\")
        else:
            firstcol = (
                f"\\multirow{{2}}{{*}}{{{label}}}"
                if ("theta" in k or "deltam" in k)
                else f"\\multirow{{2}}{{*}}{{\\textbf{{{label}}}}}"
            )
            lines.append(f"\n        {firstcol} & \\multirow{{2}}{{*}}{{{av}}} & \\multirow{{2}}{{*}}{{{c}}} & \\multirow{{2}}{{*}}{{{u}}} \\\\")
            lines.append("        & & & \\\\")

        lines.append("        \\hline")

    lines.append("\n        \\hline\\hline\n    \\end{tabular}\n\\end{table}")
    return "\n".join(lines)

def main():
    if len(sys.argv) < 2:
        print("Usage: python make_corr_nocorr_tables.py <csvfile>")
        sys.exit(1)

    data = read_csv(sys.argv[1])
    print(build_table(data, "osc"))
    print("\n\n")
    print(build_table(data, "syst"))

if __name__ == "__main__":
    main()
