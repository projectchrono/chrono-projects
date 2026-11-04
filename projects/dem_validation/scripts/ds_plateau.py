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
# Plateau stress ratio against normal stress for groups of direct shear runs
# (for example the same cell with and without rolling friction), next to Hartl
# and Ooi's (2008) single-sphere tests (their Table 2: secant bulk friction of
# the rainfall and central fillings). Runs of a group at the same normal stress
# (several seeds) are drawn as their mean, with a band over their range.
#
# Usage:
#   python ds_plateau.py --out <file without extension> --group <label> <run folder> [<run folder> ...]
#                        [--group <label> <run folder> ...]
#
# Optional:
#   --tau KEY       time series column of the shear stress (default tau, the Jenike cell)
#   --size W H      figure size in inches (default 4.47 2.98)
#   --no-lab        leave out the experiment
#
# =============================================================================

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dsplot  # noqa: E402
from dsplot import plt  # noqa: E402

# Hartl and Ooi (2008), Table 2, single spheres: secant bulk friction at 3.1, 6.3, 12.5 and 24.2 kPa
LAB_STRESSES = [3100, 6300, 12500, 24200]
LAB = [("Härtl and Ooi, rainfall", [0.61, 0.58, 0.56, 0.56], dsplot.EXP_BLACK, "s"),
       ("Härtl and Ooi, central filling", [0.54, 0.51, 0.52, 0.49], dsplot.EXP_GRAY, "D")]

# One style per group: color, line style, marker
STYLES = [(dsplot.DEM_ORANGE, "-", "o"), (dsplot.DEM_RAMP[3], "--", "^"), (dsplot.DEM_RAMP[1], ":", "v")]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group", nargs="+", action="append", required=True, metavar=("LABEL", "RUN"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--tau", default="tau")
    ap.add_argument("--size", nargs=2, type=float, default=[4.47, 2.98])
    ap.add_argument("--no-lab", action="store_true")
    args = ap.parse_args()

    dsplot.set_style()

    # Per group, per normal stress, the plateau ratios of its runs
    groups = []
    for label, *runs in args.group:
        by_sigma = {}
        for r in runs:
            meta = dsplot.read_meta(r)
            s = dsplot.summarize(dsplot.read_timeseries(r), meta, args.tau)
            by_sigma.setdefault(meta["sigma"], []).append(s["plateau_ratio"])
        groups.append((label, by_sigma))
    sigmas = sorted({sigma for _, by_sigma in groups for sigma in by_sigma})

    fig, ax = plt.subplots(figsize=args.size)
    for k, (label, by_sigma) in enumerate(groups):
        color, ls, marker = STYLES[k % len(STYLES)]
        x = sorted(by_sigma)
        v = [np.array(by_sigma[sigma]) for sigma in x]
        xk = [1e-3 * sigma for sigma in x]
        if any(len(a) > 1 for a in v):
            ax.fill_between(xk, [a.min() for a in v], [a.max() for a in v], color=color, alpha=0.25, lw=0)
        ax.plot(xk, [a.mean() for a in v], color=color, ls=ls, marker=marker, lw=1.2, ms=4.5, label=label)
    if not args.no_lab:
        for label, values, color, marker in LAB:
            ax.plot([1e-3 * s for s in LAB_STRESSES], values, color=color, marker=marker, lw=1.0, ms=4, mfc="white",
                    label=label)
    ax.set_xscale("log")
    ax.set_xticks([1e-3 * sigma for sigma in sigmas])
    ax.set_xticklabels(["%.1f" % (1e-3 * sigma) for sigma in sigmas])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlim(2.6, 29)
    ax.set_ylim(0.40, 0.68)  # room above the experiment for the legend
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.05))
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
    ax.set_xlabel("normal stress (kPa)")
    ax.set_ylabel("plateau stress ratio")
    ax.grid(color=dsplot.LIGHT_GRAY, lw=0.5)
    ax.legend(loc="upper right", handlelength=1.8, frameon=False)
    fig.tight_layout()
    dsplot.save(fig, args.out)

    # Table: per normal stress, each group's mean (and standard deviation over several runs), and the experiment
    head = "| normal stress (kPa) | " + " | ".join("%s: plateau tau/sigma_n" % label for label, _ in groups)
    lines = [head + (" | " + " | ".join(label for label, *_ in LAB) if not args.no_lab else "") + " |"]
    lines.append("|---" * (1 + len(groups) + (0 if args.no_lab else len(LAB))) + "|")
    for sigma in sigmas:
        cells = []
        for _, by_sigma in groups:
            a = np.array(by_sigma.get(sigma, []))
            if len(a) == 0:
                cells.append("")
            else:
                cells.append("%.3f" % a.mean() + (" ± %.3f (%d runs)" % (a.std(ddof=1), len(a)) if len(a) > 1 else ""))
        if not args.no_lab:
            k = int(np.argmin([abs(s - sigma) for s in LAB_STRESSES]))
            cells += ["%.2f" % values[k] for _, values, _, _ in LAB]
        lines.append("| %.1f | %s |" % (1e-3 * sigma, " | ".join(cells)))
    table = "\n".join(lines)
    with open(args.out + "_table.md", "w", encoding="utf-8") as f:
        f.write(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
