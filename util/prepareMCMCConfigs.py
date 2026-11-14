"""Prepare matching MCMC configs with proposal widths from a data-fit table"""

import argparse
import configparser
import csv
import math
from pathlib import Path
import re
import shutil


def set_value(text, section, key, value):
    pattern = rf"(\[{re.escape(section)}\]\n)(.*?)(?=\n\[|\Z)"

    def update(match):
        body, count = re.subn(
            rf"^{re.escape(key)}\s*=.*$", f"{key} = {value}", match[2], flags=re.M
        )
        if not count:
            body = body.rstrip() + f"\n{key} = {value}\n"
        return match[1] + body

    result, count = re.subn(pattern, update, text, flags=re.S)
    if count != 1:
        raise ValueError(f"Expected one section [{section}]")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--burn-in", type=int, default=0)
    parser.add_argument("--scale", type=float, default=0.5)
    parser.add_argument("--save-outputs", type=int, choices=(0, 1), default=1)
    parser.add_argument("--config-dir", type=Path, default=Path("cfg/reacAnalysis2025_mcmc_validation"))
    parser.add_argument("--output-dir", type=Path, default=Path("output/mcmc_validation_pilot"))
    args = parser.parse_args()
    if args.steps <= 0 or not 0 <= args.burn_in < args.steps:
        parser.error("Require steps > 0 and 0 <= burn-in < steps")
    if not math.isfinite(args.scale) or args.scale <= 0:
        parser.error("Scale must be finite and positive")

    reference = Path("cfg/reacAnalysis2025_mcmc_reference/data_proposal_widths.tsv")
    with reference.open() as source:
        widths = {
            row["parameter"]: (float(row["error_low"]) + float(row["error_high"])) / 2
            for row in csv.DictReader(source, delimiter="\t")
        }
    for name in ("asmv_correlated_noprior_noalphan", "data_correlated_thdmprior_alphan"):
        source = Path("cfg/reacAnalysis2025_floated_validation") / name
        destination = args.config_dir / name
        destination.mkdir(parents=True, exist_ok=True)
        output = (args.output_dir / name).resolve()
        output.mkdir(parents=True, exist_ok=True)
        for filename in ("event_config.ini", "pdf_config.ini", "syst_config.ini"):
            shutil.copy2(source / filename, destination / filename)
        text = (source / "fit_config.ini").read_text()
        config = configparser.ConfigParser()
        config.read_string(text)
        for parameter in config.sections():
            if parameter == "summary" or config.getboolean(parameter, "fix", fallback=False):
                continue
            width = widths[parameter]
            if not math.isfinite(width) or width <= 0:
                raise ValueError(f"Invalid width for {parameter}: {width}")
            text = set_value(text, parameter, "sig", format(width, ".12g"))
        for key, value in {
            "iterations": args.steps,
            "burn_in": args.burn_in,
            "sigma_scale": args.scale,
            "save_outputs": args.save_outputs,
            "output_directory": output,
        }.items():
            text = set_value(text, "summary", key, value)
        (destination / "fit_config.ini").write_text(text)
        print(destination)


if __name__ == "__main__":
    main()
