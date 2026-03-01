#!/usr/bin/env python3
"""
Parse OXSX fit output logs (many files, many fits per file) and plot:

- x = BestfitIntegral(geonu_U)  + BestfitIntegral(geonu_U2)
- y = BestfitIntegral(geonu_Th) + BestfitIntegral(geonu_Th2)
- colour = ΔLLH = (LLH - LLHmin)  [global minimum after filtering FitValid if requested]

Plot style:
- main panel: triangulated filled contours (works for scattered points; no fixed binning)
- overlaid contour lines at chosen ΔLLH levels
- top and right 1D "profiles": min ΔLLH in x/y bins
- colourbar in the top-right panel
"""

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as mticker


# --- Regex helpers ---
RE_FIT_SPLIT = re.compile(r"^[-]*\s*OXSX Fit Result:\s*(.+?)\s*$", re.MULTILINE)
RE_LLH = re.compile(r"^\s*LLH:\s*([0-9]*\.?[0-9]+)\s*$", re.MULTILINE)
RE_FITVALID = re.compile(r"^\s*FitValid:\s*([01])\s*$", re.MULTILINE)

RE_INTEGRAL_ROW = re.compile(
    r"^\s*(?P<pdf>\S+)\s+\[(?P<dataset>[^\]]+)\]\s+"
    r"(?P<asimov>[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+"
    r"(?P<bestfit>[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*$",
    re.MULTILINE,
)


def iter_files_in_dir(d: Path) -> List[Path]:
    files = []
    for p in d.rglob("*"):
        if p.is_file():
            if p.suffix.lower() in {".root", ".png", ".pdf", ".jpg", ".jpeg", ".gif", ".zip", ".gz"}:
                continue
            files.append(p)
    return sorted(files)


def split_into_fit_blocks(text: str) -> List[str]:
    matches = list(RE_FIT_SPLIT.finditer(text))
    if not matches:
        return []
    blocks = []
    for i, m in enumerate(matches):
        start = m.start()
        end = matches[i + 1].start() if (i + 1) < len(matches) else len(text)
        blocks.append(text[start:end])
    return blocks


def extract_bestfit_integrals(block: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for m in RE_INTEGRAL_ROW.finditer(block):
        pdf = m.group("pdf")
        if pdf in {"geonu_U", "geonu_U2", "geonu_Th", "geonu_Th2"}:
            out[pdf] = float(m.group("bestfit"))
    return out


def extract_llh(block: str) -> Optional[float]:
    m = RE_LLH.search(block)
    return float(m.group(1)) if m else None


def extract_fitvalid(block: str) -> Optional[int]:
    m = RE_FITVALID.search(block)
    return int(m.group(1)) if m else None


def parse_dir(directory: Path, include_invalid: bool) -> List[Tuple[float, float, float, int, str, str]]:
    """
    Returns list of tuples:
      (U_sum, Th_sum, LLH, FitValid, source_file, fit_timestamp_str)
    """
    rows: List[Tuple[float, float, float, int, str, str]] = []

    files = iter_files_in_dir(directory)
    if not files:
        print(f"No files found under: {directory}", file=sys.stderr)
        return rows

    for fp in files:
        try:
            text = fp.read_text(errors="replace")
        except Exception as e:
            print(f"Skipping unreadable file {fp}: {e}", file=sys.stderr)
            continue

        blocks = split_into_fit_blocks(text)
        if not blocks:
            continue

        for b in blocks:
            llh = extract_llh(b)
            fitvalid = extract_fitvalid(b)

            ts_match = RE_FIT_SPLIT.search(b)
            ts = ts_match.group(1).strip() if ts_match else ""

            if llh is None or fitvalid is None:
                continue
            if (not include_invalid) and fitvalid == 0:
                continue

            ints = extract_bestfit_integrals(b)
            required = {"geonu_U", "geonu_U2", "geonu_Th", "geonu_Th2"}
            if not required.issubset(ints.keys()):
                continue

            u_sum = ints["geonu_U"] + ints["geonu_U2"]
            th_sum = ints["geonu_Th"] + ints["geonu_Th2"]

            rows.append((u_sum, th_sum, llh, fitvalid, str(fp), ts))

    return rows


def write_csv(rows, out_csv: Path, llh0_raw: float) -> None:
    """
    rows are assumed to contain ΔLLH in column 3 (index 2).
    We also store LLHmin (raw) as a header comment for provenance.
    """
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8") as f:
        f.write(f"# LLHmin_raw,{llh0_raw}\n")
        f.write("U_sum,Th_sum,deltaLLH,FitValid,source_file,fit_timestamp\n")
        for u, th, dllh, fv, src, ts in rows:
            src2 = '"' + src.replace('"', '""') + '"'
            ts2 = '"' + ts.replace('"', '""') + '"'
            f.write(f"{u},{th},{dllh},{fv},{src2},{ts2}\n")


def profile_min_binned(x: np.ndarray, z: np.ndarray, edges: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return bin centers and min(z) in each x-bin (NaN if empty)."""
    centers = 0.5 * (edges[:-1] + edges[1:])
    prof = np.full(len(centers), np.nan, dtype=float)

    idx = np.digitize(x, edges) - 1
    for i in range(len(centers)):
        mask = idx == i
        if np.any(mask):
            prof[i] = np.nanmin(z[mask])
    return centers, prof


def make_plot(
    rows,
    out_img: Path,
    max_z: Optional[float],
    title: str = "",
    xmin: float = 0.0,
    xmax: float = 180.0,
    ymin: float = 0.0,
    ymax: float = 100.0,
    profile_bins: int = 60,
    contour_levels: Optional[List[float]] = None,
) -> None:
    # unpack (now deltaLLH is stored in r[2])
    xs = np.array([r[0] for r in rows], dtype=float)  # U_sum
    ys = np.array([r[1] for r in rows], dtype=float)  # Th_sum
    zs = np.array([r[2] for r in rows], dtype=float)  # ΔLLH

    # Restrict to plotting window for triangulation/profiles/colour scaling
    mask = np.isfinite(xs) & np.isfinite(ys) & np.isfinite(zs)
    mask &= (xs >= xmin) & (xs <= xmax) & (ys >= ymin) & (ys <= ymax)

    xs, ys, zs = xs[mask], ys[mask], zs[mask]
    if xs.size < 3:
        raise RuntimeError("Not enough points inside x/y window to triangulate (need >= 3).")
    
    best_idx = int(np.nanargmin(zs))
    best_x = float(xs[best_idx])
    best_y = float(ys[best_idx])

    # ΔLLH should start at ~0
    vmin = float(np.nanmin(zs))
    if vmin < -1e-6:
        # shouldn't happen, but don't crash; just proceed
        pass
    if contour_levels is None:
        contour_levels = [0.5, 2.0, 4.5]

    # Clamp/saturate for plotting if requested
    if max_z is not None:
        vmax = float(max_z)
        zs_plot = np.minimum(zs, vmax)  # values > max_z drawn at max_z
    else:
        vmax = float(np.nanmax(zs))
        zs_plot = zs

    # Layout: top profile + main + right profile + colourbar (top-right)
    fig = plt.figure(figsize=(9.0, 7.5))
    gs = GridSpec(
        nrows=2,
        ncols=2,
        width_ratios=[4.8, 1.6],
        height_ratios=[1.6, 4.8],
        wspace=0.05,
        hspace=0.05,
    )

    ax_top = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[1, 0], sharex=ax_top)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
    ax_cbar = fig.add_subplot(gs[0, 1])
    # Shrink colourbar width
    pos = ax_cbar.get_position()
    ax_cbar.set_position([pos.x0 + 0.2 * pos.width, pos.y0, 0.35 * pos.width, pos.height])

    # Main panel: triangulated filled contours + contour lines
    tri = mtri.Triangulation(xs, ys)
    # Smooth triangle shading for the filled field (no contour-band seams)
    m = ax_main.tripcolor(
    tri,
    zs_plot,
    shading="gouraud",
    vmin=vmin,
    vmax=vmax,
    rasterized=True,   # rasterize only the colour field in PDF
    )

    # Colourbar should use this mappable now
    cbar = fig.colorbar(m, cax=ax_cbar)
    cbar.set_label("ΔLLH (LLH - LLHmin)", labelpad=10)

    # Contour colours: 1σ=white, 2σ=red, 3σ=yellow
    default_colour_map = {
        0.5: "white",
        2.0: "red",
        4.5: "yellow",
    }

    # Build colour list in same order as contour_levels
    contour_colours = [default_colour_map.get(lv, "white") for lv in contour_levels]

    cs = ax_main.tricontour(
        tri,
        zs,
        levels=contour_levels,
        linewidths=1.4,
        colors=contour_colours,
    )

    ax_main.plot(best_x, best_y, marker="x", markersize=7, markeredgewidth=1.6, color="red")
    # optional label:
    # ax_main.annotate("best", (best_x, best_y), xytext=(5, 5), textcoords="offset points", color="red", fontsize=10)

    sigma_map = {
        0.5: r"$1\sigma$",
        2.0: r"$2\sigma$",
        4.5: r"$3\sigma$",
    }
    fmt = {lv: sigma_map.get(lv, f"{lv:g}") for lv in cs.levels}

    ax_main.clabel(cs, fmt=fmt, fontsize=11)

    m = ax_main.tripcolor(
        tri,
        zs_plot,
        shading="gouraud",
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )

    cbar = fig.colorbar(m, cax=ax_cbar)
    cbar.set_label("ΔLLH (LLH - LLHmin)", labelpad=10)

    # Nice tick spacing
    cbar.locator = mticker.MultipleLocator(5)
    cbar.update_ticks()

    ax_cbar.yaxis.set_ticks_position("right")
    ax_cbar.yaxis.set_label_position("right")
    ax_cbar.set_xticks([])  

    # Profiles: min ΔLLH in bins along x and y
    x_edges = np.linspace(xmin, xmax, profile_bins + 1)
    y_edges = np.linspace(ymin, ymax, profile_bins + 1)

    x_cent, prof_x = profile_min_binned(xs, zs, x_edges)  # min over y in each x bin
    y_cent, prof_y = profile_min_binned(ys, zs, y_edges)  # min over x in each y bin

    ax_top.plot(x_cent, prof_x)
    ax_top.set_ylabel("Profile min ΔLLH", fontsize=12)
    ax_top.tick_params(labelbottom=False)
    ax_top.grid(True, alpha=0.3)
    ax_top.set_ylim(bottom=0.0)

    ax_right.plot(prof_y, y_cent)
    ax_right.set_xlabel("Profile min ΔLLH", fontsize=12)
    ax_right.tick_params(labelleft=False)
    ax_right.xaxis.set_major_locator(mticker.MaxNLocator(nbins=4, prune="lower"))
    ax_right.grid(True, alpha=0.3)
    ax_right.set_xlim(left=0.0)

    # Labels/limits (note: x is U_sum, y is Th_sum)
    ax_main.set_xlabel("Total U Geo Rate", fontsize=16)
    ax_main.set_ylabel("Total Th Geo Rate", fontsize=16)
    ax_main.set_xlim(xmin, xmax)
    ax_main.set_ylim(ymin, ymax)
    ax_main.grid(True, alpha=0.25)

    if title:
        fig.suptitle(title, fontsize=16, y=0.95)

    fig.tight_layout()
    fig.savefig(out_img, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Parse OXSX fit logs and plot ΔLLH vs summed bestfit U/Th integrals (triangulated + profiles)."
    )
    ap.add_argument("directory", help="Directory containing output text files")
    ap.add_argument("output", help="Output image filename (e.g. plot.png or plot.pdf)")

    ap.add_argument(
        "--include-invalid",
        action="store_true",
        help="Include fits with FitValid: 0 (default: only FitValid: 1)",
    )
    ap.add_argument(
        "--max-z",
        type=float,
        default=None,
        help="Clamp colour scale maximum (values above are saturated). Interpreted in ΔLLH units.",
    )

    ap.add_argument("--xmin", type=float, default=0.0, help="Minimum x-axis (default 0)")
    ap.add_argument("--xmax", type=float, default=100.0, help="Maximum x-axis (default 180)")
    ap.add_argument("--ymin", type=float, default=0.0, help="Minimum y-axis (default 0)")
    ap.add_argument("--ymax", type=float, default=100.0, help="Maximum y-axis (default 100)")

    ap.add_argument(
        "--profile-bins",
        type=int,
        default=60,
        help="Number of bins used for 1D profile minima (default 60)",
    )
    ap.add_argument(
        "--contours",
        type=float,
        nargs="*",
        default=[0.5, 2.0, 4.5],
        help="Contour ΔLLH levels (default: 0.5 2.0 4.5)",
    )
    ap.add_argument("--title", default="", help="Optional plot title")

    args = ap.parse_args()

    d = Path(args.directory).expanduser().resolve()
    out_img = Path(args.output).expanduser().resolve()

    if not d.exists() or not d.is_dir():
        print(f"Not a directory: {d}", file=sys.stderr)
        sys.exit(2)

    rows = parse_dir(d, include_invalid=args.include_invalid)
    if not rows:
        print("No fit rows parsed. Check logs contain 'OXSX Fit Result:' blocks and the integral table.", file=sys.stderr)
        sys.exit(1)

    # Shift LLH so global best (minimum) is 0
    llh0_raw = min(r[2] for r in rows)
    rows = [(u, th, llh - llh0_raw, fv, src, ts) for (u, th, llh, fv, src, ts) in rows]

    # Write CSV next to plot
    out_csv = out_img.with_suffix(".csv")
    write_csv(rows, out_csv, llh0_raw=llh0_raw)

    title = args.title

    make_plot(
        rows,
        out_img,
        max_z=args.max_z,
        title=title,
        xmin=args.xmin,
        xmax=args.xmax,
        ymin=args.ymin,
        ymax=args.ymax,
        profile_bins=args.profile_bins,
        contour_levels=args.contours,
    )

    print(f"Wrote plot: {out_img}")
    print(f"Wrote data: {out_csv}")
    print(f"LLHmin_raw subtracted: {llh0_raw:.6g}")
    print(f"Points parsed (after FitValid filter): {len(rows)}")


if __name__ == "__main__":
    main()
