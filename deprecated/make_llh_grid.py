#!/usr/bin/env python3
"""
Parse OXSX fit output logs and build a 2D ΔLLH grid plot with 1D profiles.
"""
from __future__ import annotations

import argparse
import glob
import math
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams.update({
    "axes.labelsize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
})

FIT_START_RE = re.compile(r"^.*OXSX Fit Result:\s*(.+)\s*$")
BEST_FIT_HEADER_RE = re.compile(r"^\s*Best Fit Values:\s*$")
BEST_FIT_LINE_RE = re.compile(r"^\s*([A-Za-z0-9_]+)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$")
LLH_RE = re.compile(r"^\s*LLH:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$")
FITVALID_RE = re.compile(r"^\s*FitValid:\s*([01])\s*$")

REQUIRED_PARAMS = (
    "geonu_Th",
    "geonu_Th2",
    "geonu_Th_norm",
    "geonu_U",
    "geonu_U2",
    "geonu_U_norm",
)

@dataclass
class FitRecord:
    source_file: str
    fit_timestamp: str = ""
    bestfit: Dict[str, float] = field(default_factory=dict)
    llh: Optional[float] = None
    fitvalid: Optional[int] = None

    def has_required(self) -> bool:
        if self.llh is None or self.fitvalid is None:
            return False
        return all(k in self.bestfit for k in REQUIRED_PARAMS)

    def derived_axes(self) -> Tuple[float, float]:
        th = (self.bestfit["geonu_Th"] + self.bestfit["geonu_Th2"]) * self.bestfit["geonu_Th_norm"]
        uu = (self.bestfit["geonu_U"]  + self.bestfit["geonu_U2"])  * self.bestfit["geonu_U_norm"]
        return uu, th


def parse_one_file(path: str) -> List[FitRecord]:
    records: List[FitRecord] = []

    current: Optional[FitRecord] = None
    in_bestfit_block = False
    bestfit_started = False

    def flush_current():
        nonlocal current, in_bestfit_block
        if current is not None:
            records.append(current)
        current = None
        in_bestfit_block = False
        
    with open(path, "r", errors="replace") as f:
        for line in f:
            # Start of a fit
            m = FIT_START_RE.match(line)
            if m:
                # flush previous fit (even if incomplete)
                flush_current()
                current = FitRecord(source_file=os.path.basename(path), fit_timestamp=m.group(1).strip())
                continue

            if current is None:
                continue

            # Best fit values section header
            if BEST_FIT_HEADER_RE.match(line):
                in_bestfit_block = True
                bestfit_started = False
                continue

            # Parse best-fit lines until we hit a blank line or a non-matching line after we've started reading values
            if in_bestfit_block:
                if line.strip() == "":
                    # Many logs have a blank line immediately after "Best Fit Values:"
                    # Only treat blank as end-of-block once we've actually read some values.
                    if bestfit_started:
                        in_bestfit_block = False
                    continue
                bm = BEST_FIT_LINE_RE.match(line)
                if bm:
                    bestfit_started = True
                    name = bm.group(1)
                    val = float(bm.group(2))
                    current.bestfit[name] = val
                # else: ignore lines like separators
                continue

            # LLH line
            lm = LLH_RE.match(line)
            if lm:
                current.llh = float(lm.group(1))
                continue

            # FitValid line; often this is effectively the end-of-fit marker
            vm = FITVALID_RE.match(line)
            if vm:
                current.fitvalid = int(vm.group(1))
                continue

    # flush last
    flush_current()
    return records


def build_grid(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return (u_vals, th_vals, grid) where grid[j,i] corresponds to th_vals[j], u_vals[i].
    """
    u_vals = np.array(sorted(df["U_total"].unique()))
    th_vals = np.array(sorted(df["Th_total"].unique()))

    u_index = {v: i for i, v in enumerate(u_vals)}
    th_index = {v: j for j, v in enumerate(th_vals)}

    grid = np.full((len(th_vals), len(u_vals)), np.nan, dtype=float)
    for _, row in df.iterrows():
        i = u_index[row["U_total"]]
        j = th_index[row["Th_total"]]
        grid[j, i] = row["deltaLLH"]

    return u_vals, th_vals, grid


def profile_min_over_axis(grid: np.ndarray, axis: int) -> np.ndarray:
    """Min ignoring NaN along axis."""
    return np.nanmin(grid, axis=axis)


def plot_with_profiles(u_vals, th_vals, grid, out_pdf: str, title: str = "", max_delta: Optional[float] = None):
    # Profiles:
    prof_u = profile_min_over_axis(grid, axis=0)   # min over Th (rows) -> per U (cols)
    prof_th = profile_min_over_axis(grid, axis=1)  # min over U (cols) -> per Th (rows)

    # Set colour scale upper bound
    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        raise RuntimeError("No finite ΔLLH values to plot.")
    if max_delta is None:
        max_delta = float(np.nanpercentile(finite, 99.0))
        if max_delta <= 0:
            max_delta = float(np.nanmax(finite))

    # Layout: main + top profile + right profile + palette
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(9.0, 7.5))
    gs = GridSpec(
        nrows=2, ncols=2,
        width_ratios=[4.8, 1.6],
        height_ratios=[1.6, 4.8],
        wspace=0.05, hspace=0.05
    )

    ax_top = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[1, 0], sharex=ax_top)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
    ax_cbar  = fig.add_subplot(gs[0, 1])

    # Make bin edges for pcolormesh
    def edges(vals: np.ndarray) -> np.ndarray:
        if len(vals) == 1:
            # arbitrary single-bin width
            d = 1.0
            return np.array([vals[0] - d/2, vals[0] + d/2])
        mids = 0.5 * (vals[:-1] + vals[1:])
        first = vals[0] - (mids[0] - vals[0])
        last = vals[-1] + (vals[-1] - mids[-1])
        return np.concatenate([[first], mids, [last]])

    u_edges = edges(u_vals)
    th_edges = edges(th_vals)

    # Main heatmap
    m = ax_main.pcolormesh(u_edges, th_edges, grid, shading="auto", vmin=0.0, vmax=max_delta)
    cbar = fig.colorbar(m, cax=ax_cbar)
    cbar.set_label(r"$\Delta \mathrm{LLH}$", labelpad=10)

    ax_cbar.yaxis.set_ticks_position("right")
    ax_cbar.yaxis.set_label_position("right")
    ax_cbar.set_xticks([])
    
    pos = ax_cbar.get_position()
    ax_cbar.set_position([pos.x0, pos.y0, 0.3 * pos.width, pos.height])

    ax_main.set_xlabel(r"Total U Geo Rate")
    ax_main.set_ylabel(r"Total Th Geo Rate")

    # Draw contours
    try_levels = [0.5, 2.0, 4.5]
    sigma_map = {
    0.5: r"$1\sigma$",
    2.0: r"$2\sigma$",
    4.5: r"$3\sigma$",
    }
    # Contours need grid centers
    Uc, THc = np.meshgrid(u_vals, th_vals)
    finite_mask = np.isfinite(grid)
    if np.any(finite_mask):
        # Only contour if there are enough points
        if len(u_vals) >= 2 and len(th_vals) >= 2:
            cs = ax_main.contour(Uc, THc, grid, levels=try_levels, linewidths=1.0, colors=["white", "red", "yellow"])
            
    ax_main.clabel(cs,fmt=sigma_map,fontsize=10)

    # Top profile
    ax_top.plot(u_vals, prof_u)
    ax_top.set_ylabel(r"Profile U min $\Delta\mathrm{LLH}$", fontsize = 12)
    ax_top.tick_params(labelbottom=False)
    ax_top.grid(True, alpha=0.3)

    # Right profile
    ax_right.plot(prof_th, th_vals)
    ax_right.set_xlabel(r"Profile Th min $\Delta\mathrm{LLH}$", fontsize = 12)
    ax_right.tick_params(labelleft=False)
    ax_right.grid(True, alpha=0.3)

    # Nice limits
    ax_main.set_xlim(u_edges[0], u_edges[-1])
    ax_main.set_ylim(th_edges[0], th_edges[-1])

    if title:
        fig.suptitle(title, fontsize=18, y = 0.95)

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Make 2D ΔLLH grid + 1D profiles from OXSX output logs.")
    ap.add_argument("--pattern", default="*.output", help="Glob pattern for input files (default: *.output)")
    ap.add_argument("--out-pdf", default="llh_grid.pdf", help="Output PDF filename")
    ap.add_argument("--out-csv", default="llh_points.csv", help="Output CSV filename of all parsed points")
    ap.add_argument("--include-invalid", action="store_true", help="Include FitValid==0 points (default: drop them)")
    ap.add_argument("--max-delta", type=float, default=None, help="Cap colour scale at this ΔLLH (default: auto)")
    ap.add_argument("--title", default="", help="Optional plot title")
    args = ap.parse_args()

    paths = sorted(glob.glob(args.pattern))
    if not paths:
        raise SystemExit(f"No files matched pattern: {args.pattern}")

    all_recs: List[FitRecord] = []
    for p in paths:
        all_recs.extend(parse_one_file(p))

    # Build rows
    rows = []
    n_incomplete = 0
    for r in all_recs:
        if not r.has_required():
            n_incomplete += 1
            continue
        u_total, th_total = r.derived_axes()
        rows.append(
            dict(
                source_file=r.source_file,
                fit_timestamp=r.fit_timestamp,
                LLH=r.llh,
                FitValid=r.fitvalid,
                U_total=u_total,
                Th_total=th_total,
                geonu_U=r.bestfit["geonu_U"],
                geonu_U2=r.bestfit["geonu_U2"],
                geonu_U_norm=r.bestfit["geonu_U_norm"],
                geonu_Th=r.bestfit["geonu_Th"],
                geonu_Th2=r.bestfit["geonu_Th2"],
                geonu_Th_norm=r.bestfit["geonu_Th_norm"],
            )
        )

    if not rows:
        raise SystemExit("Parsed zero complete fit records (missing required params/LLH/FitValid).")

    df = pd.DataFrame(rows)

    if not args.include_invalid:
        df = df[df["FitValid"] == 1].copy()

    if df.empty:
        raise SystemExit("After filtering, no rows left to plot.")

    # ΔLLH relative to global best
    best_llh = df["LLH"].min()
    df["deltaLLH"] = df["LLH"] - best_llh

    # Save CSV of points
    df.sort_values(["Th_total", "U_total"], inplace=True)
    df.to_csv(args.out_csv, index=False)

    # Build grid
    u_vals, th_vals, grid = build_grid(df)

    # Plot
    title = args.title
    if not title:
        title = f"ΔLLH for Fixed Geo Rates (min LLH = {best_llh:.6g})"
    plot_with_profiles(u_vals, th_vals, grid, out_pdf=args.out_pdf, title=title, max_delta=args.max_delta)

    # Small stdout summary
    print(f"Read {len(paths)} files, parsed {len(all_recs)} fits.")
    if n_incomplete:
        print(f"Skipped {n_incomplete} incomplete fits (missing required fields).")
    print(f"Keeping {len(df)} fits for plot. Best LLH = {best_llh:.6g}")
    print(f"Wrote: {args.out_csv}")
    print(f"Wrote: {args.out_pdf}")


if __name__ == "__main__":
    main()
