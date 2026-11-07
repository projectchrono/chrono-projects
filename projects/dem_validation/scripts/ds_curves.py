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
# Shear stress against shear displacement for direct shear runs at several
# normal stresses, on the axes of Hartl and Ooi (2008), Fig. 2 (0 to 9 mm,
# 0 to 16 kPa), and a table of peak and residual stress ratios and initial
# void ratios. Runs at the same normal stress (several seeds) are drawn as
# their mean curve with a band over their range, and the table gives the mean
# and standard deviation over them.
#
# Usage:
#   python ds_curves.py --out <file without extension> <run folder> [<run folder> ...]
#
# Optional:
#   --tau KEY                 time series column of the shear stress (default tau_bin: square box, force on the
#                             lower box; tau_upper: on the upper box and lid; tau: Jenike cell, on top ring and lid)
#   --size W H                figure size in inches (default 4.47 2.98)
#   --dilation FILE           also plot the lid rise against shear displacement (FILE without extension)
#
# =============================================================================

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dsplot  # noqa: E402
from dsplot import plt  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tau", default="tau_bin")
    ap.add_argument("--size", nargs=2, type=float, default=[4.47, 2.98])
    ap.add_argument("--dilation", default=None)
    args = ap.parse_args()

    dsplot.set_style()

    # Runs at the same normal stress (several seeds) are one group: the mean curve, with the range of the runs as a band
    groups = {}
    for r in args.runs:
        meta = dsplot.read_meta(r)
        ts = dsplot.read_timeseries(r)
        groups.setdefault(meta["sigma"], []).append((meta, ts, dsplot.summarize(ts, meta, args.tau)))
    sigmas = sorted(groups)
    colors = dsplot.DEM_RAMP if len(sigmas) <= len(dsplot.DEM_RAMP) else [dsplot.DEM_ORANGE] * len(sigmas)

    def curves(key, scale):
        """Common shear displacement grid (mm) and, per group, the curves of 'key' on it."""
        out = []
        for sigma in sigmas:
            runs = []
            for meta, ts, s in groups[sigma]:
                idx = dsplot.shear_part(ts)
                runs.append((np.array([1e3 * ts["shear_disp"][i] for i in idx]), np.array([scale * ts[key][i] for i in idx])))
            d_end = min(d[-1] for d, _ in runs)
            grid = np.linspace(0, d_end, int(round(d_end / 0.01)) + 1)
            out.append((grid, np.array([np.interp(grid, d, v) for d, v in runs])))
        return out

    def plot(ax, data):
        for k, (sigma, (grid, values)) in enumerate(zip(sigmas, data)):
            label = "normal stress = %.1f kPa" % (1e-3 * sigma)
            if len(values) > 1:
                ax.fill_between(grid, values.min(axis=0), values.max(axis=0), color=colors[k], alpha=0.25, lw=0)
            ax.plot(grid, values.mean(axis=0), color=colors[k], lw=1.2, label=label)

    fig, ax = plt.subplots(figsize=args.size)
    plot(ax, curves(args.tau, 1e-3))
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 16)
    ax.set_xlabel("shear displacement (mm)")
    ax.set_ylabel("shear stress (kPa)")
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    dsplot.save(fig, args.out)

    if args.dilation:
        fig, ax = plt.subplots(figsize=args.size)
        plot(ax, curves("dilation", 1e3))
        ax.set_xlim(0, 9)
        ax.set_xlabel("shear displacement (mm)")
        ax.set_ylabel("lid rise (mm)")
        ax.axhline(0, color=dsplot.LIGHT_GRAY, lw=0.6)
        ax.legend(loc="upper left", frameon=False)
        fig.tight_layout()
        dsplot.save(fig, args.dilation)

    # Table: per normal stress, the mean over the runs (and the standard deviation, for several runs)
    def stat(sigma, key, scale=1.0, digits=3):
        v = np.array([scale * s[key] for _, _, s in groups[sigma]])
        text = "%.*f" % (digits, v.mean())
        return text + (" ± %.*f" % (digits, v.std(ddof=1)) if len(v) > 1 else "")

    lines = ["| normal stress (kPa) | runs | peak tau/sigma_n | at (mm) | plateau tau/sigma_n (4 to 7.5 mm) | "
             "plateau tau (kPa) | residual tau/sigma_n (last 2 mm) | residual tau (kPa) | initial void ratio | "
             "initial porosity | lid rise at end (mm) | history overflows |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for sigma in sigmas:
        lines.append("| %.1f | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %d |" % (
            1e-3 * sigma, len(groups[sigma]), stat(sigma, "peak_ratio"), stat(sigma, "peak_disp", 1e3, 2),
            stat(sigma, "plateau_ratio"), stat(sigma, "plateau_tau", 1e-3, 2), stat(sigma, "residual_ratio"),
            stat(sigma, "residual_tau", 1e-3, 2), stat(sigma, "e0"), stat(sigma, "porosity0"),
            stat(sigma, "dilation_end", 1e3, 2), sum(s["overflows"] for _, _, s in groups[sigma])))
    table = "\n".join(lines)
    with open(args.out + "_table.md", "w", encoding="utf-8") as f:
        f.write(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
