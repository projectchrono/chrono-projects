# =============================================================================
# PROJECT CHRONO - http://projectchrono.org
#
# Copyright (c) 2026 projectchrono.org
# All rights reserved.
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file at the top level of the distribution and at
# http://projectchrono.org/license-chrono.txt.
#
# =============================================================================
# Author: Radu Serban
# =============================================================================
#
# Gravity: stress ratio against shear displacement for direct shear runs at the
# same normal stress under different gravities (one panel per normal stress, one
# line per gravity), and a table of the plateau and residual ratios, the initial
# porosity, the lid rise and how the normal load splits between the lid and the
# beads above the shear plane.
#
# Usage:
#   python ds_gravity.py --out <file without extension> <run folder> [<run folder> ...]
#
# Optional:
#   --tau KEY             time series column of the shear stress (default tau, the Jenike cell)
#   --size W H            figure size in inches (default 4.47 2.98)
#   --dilation FILE       also plot the lid rise against shear displacement (FILE without extension)
#   --smooth MM           running mean over this much shear displacement (default 0: none)
#
# =============================================================================

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dsplot  # noqa: E402
from dsplot import plt  # noqa: E402

# One style per gravity, from the highest: DEM orange, then dark brown dashed
STYLES = [(dsplot.DEM_ORANGE, "-"), (dsplot.DEM_RAMP[3], "--"), (dsplot.DEM_RAMP[1], ":")]


def running_mean(x, y, width):
    """Mean of y over a window of the given width in x (x evenly spaced); y itself for a zero width."""
    if width <= 0 or len(x) < 2:
        return y
    n = max(1, int(round(width / (x[1] - x[0]))))
    kernel = np.ones(n) / n
    return np.convolve(np.pad(y, (n // 2, n - 1 - n // 2), mode="edge"), kernel, mode="valid")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tau", default="tau")
    ap.add_argument("--size", nargs=2, type=float, default=[4.47, 2.98])
    ap.add_argument("--dilation", default=None)
    ap.add_argument("--smooth", type=float, default=0.0)
    args = ap.parse_args()

    dsplot.set_style()

    # Per normal stress, the runs by gravity (highest first); each with its curves on an even 0.01 mm grid
    groups = {}
    for r in args.runs:
        meta = dsplot.read_meta(r)
        ts = dsplot.read_timeseries(r)
        idx = dsplot.shear_part(ts)
        d = np.array([1e3 * ts["shear_disp"][i] for i in idx])
        grid = np.linspace(0, d[-1], int(round(d[-1] / 0.01)) + 1)
        ratio = np.interp(grid, d, [ts[args.tau][i] / meta["sigma"] for i in idx])
        rise = np.interp(grid, d, [1e3 * ts["dilation"][i] for i in idx])
        run = {"name": os.path.basename(os.path.normpath(r)), "meta": meta, "grid": grid,
               "ratio": running_mean(grid, ratio, args.smooth), "rise": rise,
               "summary": dsplot.summarize(ts, meta, args.tau)}
        groups.setdefault(meta["sigma"], []).append(run)
    sigmas = sorted(groups)
    gravities = sorted({run["meta"]["gravity"] for runs in groups.values() for run in runs}, reverse=True)

    def figure(key, ylabel, ylim, path, legend_panel, legend_loc, legend_anchor):
        fig, axes = plt.subplots(1, len(sigmas), figsize=args.size, sharex=True, sharey=True, squeeze=False)
        axes = axes[0]
        for ax, sigma in zip(axes, sigmas):
            for run in sorted(groups[sigma], key=lambda r: -r["meta"]["gravity"]):
                g = run["meta"]["gravity"]
                color, ls = STYLES[gravities.index(g) % len(STYLES)]
                ax.plot(run["grid"], run[key], color=color, ls=ls, lw=1.2, label="gravity %.2f m/s²" % g)
            ax.text(0.04, 0.97, "%.1f kPa" % (1e-3 * sigma), transform=ax.transAxes, ha="left", va="top", fontsize=9)
            ax.set_xlabel("shear displacement (mm)")
            ax.grid(color=dsplot.LIGHT_GRAY, lw=0.5)
        axes[0].set_ylabel(ylabel)
        axes[0].set_xlim(0, 12)
        axes[0].xaxis.set_major_locator(plt.MultipleLocator(3))
        axes[0].set_ylim(*ylim)
        axes[legend_panel].legend(loc=legend_loc, bbox_to_anchor=legend_anchor, handlelength=1.8, frameon=False)
        fig.tight_layout()
        dsplot.save(fig, path)

    # Legends where the curves leave room: low in the first panel for the stress ratio, high in the last for the rise
    figure("ratio", "shear stress / normal stress", (0, 0.7), args.out, 0, "lower right", (1, 0))
    if args.dilation:
        figure("rise", "lid rise (mm)", (-0.2, 3.5), args.dilation, -1, "upper left", (0, 0.88))

    lines = ["| normal stress (kPa) | gravity (m/s²) | run | lid's pressure on the sample (Pa) | "
             "beads above the shear plane (share of sigma_n) | initial porosity | plateau tau/sigma_n (4 to 7.5 mm) | "
             "residual tau/sigma_n (last 2 mm) | lid rise at end (mm) |", "|---|---|---|---|---|---|---|---|---|"]
    for sigma in sigmas:
        for run in sorted(groups[sigma], key=lambda r: -r["meta"]["gravity"]):
            m, s = run["meta"], run["summary"]
            area = m["nominal_area"]
            lines.append("| %.1f | %.2f | %s | %.0f | %.2f | %.3f | %.3f | %.3f | %.2f |" % (
                1e-3 * sigma, m["gravity"], run["name"], m["lid_mass"] * m["gravity"] / area,
                m["top_ring_bead_mass"] * m["gravity"] / (sigma * area), s["porosity0"], s["plateau_ratio"],
                s["residual_ratio"], 1e3 * s["dilation_end"]))
    table = "\n".join(lines)
    with open(args.out + "_table.md", "w", encoding="utf-8") as f:
        f.write(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
