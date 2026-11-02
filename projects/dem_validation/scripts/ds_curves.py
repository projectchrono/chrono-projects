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
# void ratios.
#
# Usage:
#   python ds_curves.py --out <file without extension> <run folder> [<run folder> ...]
#
# Optional:
#   --tau tau_bin|tau_upper   shear force on the lower box (2015 definition) or on the upper box and lid
#   --size W H                figure size in inches (default 4.47 2.98, the slot on slide 13)
#   --dilation FILE           also plot the lid rise against shear displacement (FILE without extension)
#
# =============================================================================

import argparse
import os
import sys

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
    runs = []
    for r in args.runs:
        meta = dsplot.read_meta(r)
        ts = dsplot.read_timeseries(r)
        runs.append((meta, ts, dsplot.summarize(ts, meta, args.tau)))
    runs.sort(key=lambda x: x[0]["sigma"])
    colors = dsplot.DEM_RAMP if len(runs) <= len(dsplot.DEM_RAMP) else None

    fig, ax = plt.subplots(figsize=args.size)
    for k, (meta, ts, s) in enumerate(runs):
        idx = dsplot.shear_part(ts)
        d = [1e3 * ts["shear_disp"][i] for i in idx]
        tau = [1e-3 * ts[args.tau][i] for i in idx]
        c = colors[k] if colors else dsplot.DEM_ORANGE
        ax.plot(d, tau, color=c, lw=1.2, label="normal stress = %.1f kPa" % (1e-3 * meta["sigma"]))
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 16)
    ax.set_xlabel("shear displacement (mm)")
    ax.set_ylabel("shear stress (kPa)")
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    dsplot.save(fig, args.out)

    if args.dilation:
        fig, ax = plt.subplots(figsize=args.size)
        for k, (meta, ts, s) in enumerate(runs):
            idx = dsplot.shear_part(ts)
            d = [1e3 * ts["shear_disp"][i] for i in idx]
            dy = [1e3 * ts["dilation"][i] for i in idx]
            c = colors[k] if colors else dsplot.DEM_ORANGE
            ax.plot(d, dy, color=c, lw=1.2, label="normal stress = %.1f kPa" % (1e-3 * meta["sigma"]))
        ax.set_xlim(0, 9)
        ax.set_xlabel("shear displacement (mm)")
        ax.set_ylabel("lid rise (mm)")
        ax.axhline(0, color=dsplot.LIGHT_GRAY, lw=0.6)
        ax.legend(loc="upper left", frameon=False)
        fig.tight_layout()
        dsplot.save(fig, args.dilation)

    # Table
    lines = ["| normal stress (kPa) | peak tau/sigma_n | at (mm) | plateau tau/sigma_n (4 to 7.5 mm) | plateau tau (kPa) | "
             "residual tau/sigma_n (last 2 mm) | residual tau (kPa) | initial void ratio | initial porosity | "
             "lid rise at end (mm) | history overflows |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for meta, ts, s in runs:
        lines.append("| %.1f | %.3f | %.2f | %.3f | %.2f | %.3f | %.2f | %.3f | %.3f | %.2f | %d |" % (
            1e-3 * s["sigma"], s["peak_ratio"], 1e3 * s["peak_disp"], s["plateau_ratio"], 1e-3 * s["plateau_tau"],
            s["residual_ratio"], 1e-3 * s["residual_tau"], s["e0"], s["porosity0"], 1e3 * s["dilation_end"],
            s["overflows"]))
    table = "\n".join(lines)
    with open(args.out + "_table.md", "w") as f:
        f.write(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
