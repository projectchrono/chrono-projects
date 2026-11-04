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
# Bulk friction against initial porosity: the plateau stress ratio of direct
# shear runs prepared at different porosities (loose, rainfall, dense), next to
# Hartl and Ooi's (2008) single-sphere tests with three filling methods (their
# Tables 2 and 3: secant bulk friction per normal stress, porosity before
# shearing per filling method). One panel per normal stress, with the model's
# runs from the loosest to the densest and the experiment's three fillings.
#
# Usage:
#   python ds_porosity.py --out <file without extension> <run folder> [<run folder> ...]
#
# Optional:
#   --tau KEY       time series column of the shear stress (default tau, the Jenike cell)
#   --size W H      figure size in inches (default 4.47 2.98)
#   --no-lab        leave out the experiment
#   --names A B C   names of the model's preparations, loosest first, written in the first panel
#                   (default loose rainfall dense; --names alone: none)
#
# =============================================================================

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dsplot  # noqa: E402
from dsplot import plt  # noqa: E402

# Hartl and Ooi (2008), single spheres: porosity before shearing (Table 3) and secant bulk friction at 3.1, 6.3, 12.5
# and 24.2 kPa (Table 2), for the central, rainfall and compacted fillings
LAB_STRESSES = [3100, 6300, 12500, 24200]
TEXT_GRAY = "#404040"
LAB = {
    "central": (0.411, [0.54, 0.51, 0.52, 0.49]),
    "rainfall": (0.400, [0.61, 0.58, 0.56, 0.56]),
    "compacted": (0.382, [0.60, 0.55, 0.56, 0.56]),
}


def lab_index(sigma):
    """Index of the experiment's normal stress closest to sigma (6.4 kPa runs match the paper's 6.3)."""
    return int(np.argmin([abs(s - sigma) for s in LAB_STRESSES]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tau", default="tau")
    ap.add_argument("--size", nargs=2, type=float, default=[4.47, 2.98])
    ap.add_argument("--no-lab", action="store_true")
    ap.add_argument("--names", nargs="*", default=["loose", "rainfall", "dense"],
                    help="names of the model's preparations, from the loosest to the densest (none: no names)")
    args = ap.parse_args()

    dsplot.set_style()

    # Per normal stress, the runs sorted by initial porosity
    groups = {}
    for r in args.runs:
        meta = dsplot.read_meta(r)
        ts = dsplot.read_timeseries(r)
        s = dsplot.summarize(ts, meta, args.tau)
        groups.setdefault(meta["sigma"], []).append((s["porosity0"], s, os.path.basename(os.path.normpath(r))))
    sigmas = sorted(groups)

    fig, axes = plt.subplots(1, len(sigmas), figsize=args.size, sharex=True, sharey=True, squeeze=False)
    axes = axes[0]
    for j, (ax, sigma) in enumerate(zip(axes, sigmas)):
        pts = sorted(groups[sigma], key=lambda p: p[0])
        ax.plot([p[0] for p in pts], [p[1]["plateau_ratio"] for p in pts], "-o", color=dsplot.DEM_ORANGE, lw=1.2, ms=4.5,
                label="model")
        if j == 0 and args.names:  # name the model's preparations once, from the loosest to the densest
            # the loosest below its point, the densest above, the ones between below and to the right (clear of the line)
            ordered = list(reversed(pts))
            for k, ((n, s, _), name) in enumerate(zip(ordered, args.names)):
                if k == 0:
                    offset, ha, va = (0, -6), "center", "top"
                elif k == len(ordered) - 1:
                    offset, ha, va = (0, 5), "center", "bottom"
                else:
                    offset, ha, va = (5, -5), "left", "top"
                ax.annotate(name, (n, s["plateau_ratio"]), xytext=offset, textcoords="offset points", ha=ha, va=va,
                            fontsize=8, color=TEXT_GRAY)
        if not args.no_lab:
            k = lab_index(sigma)
            lab = sorted((n, f[k], name) for name, (n, f) in LAB.items())
            ax.plot([p[0] for p in lab], [p[1] for p in lab], "--s", color=dsplot.EXP_BLACK, lw=1.0, ms=4.5, mfc="white",
                    label="Härtl and Ooi")
            if j == 0:  # name the experiment's fillings once (central, the loosest, beside its point: it is at the edge)
                for n, f, name in lab:
                    right = name == "central"
                    ax.annotate(name, (n, f), xytext=(6, 0) if right else (0, 5), textcoords="offset points",
                                ha="left" if right else "center", va="center" if right else "bottom", fontsize=8,
                                color=TEXT_GRAY)
        ax.text(0.04, 0.97, "%.1f kPa" % (1e-3 * sigma), transform=ax.transAxes, ha="left", va="top", fontsize=9)
        ax.set_xlabel("initial porosity")
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.01))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
        ax.grid(color=dsplot.LIGHT_GRAY, lw=0.5)
    axes[0].set_ylabel("plateau stress ratio")
    axes[0].set_xlim(0.415, 0.365)  # shared: from the loosest sample on the left to the densest on the right
    axes[0].set_ylim(0.40, 0.67)
    axes[0].legend(loc="lower right", handlelength=1.4, frameon=False)
    fig.tight_layout()
    dsplot.save(fig, args.out)

    lines = ["| normal stress (kPa) | run | initial porosity | plateau tau/sigma_n (4 to 7.5 mm) | "
             "residual tau/sigma_n (last 2 mm) | lid rise at end (mm) |", "|---|---|---|---|---|---|"]
    for sigma in sigmas:
        for n, s, name in sorted(groups[sigma], key=lambda p: p[0]):
            lines.append("| %.1f | %s | %.3f | %.3f | %.3f | %.2f |" % (1e-3 * sigma, name, n, s["plateau_ratio"],
                                                                     s["residual_ratio"], 1e3 * s["dilation_end"]))
    table = "\n".join(lines)
    with open(args.out + "_table.md", "w", encoding="utf-8") as f:
        f.write(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
