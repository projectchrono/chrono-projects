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
# Contact history comparison: stress ratio against shear displacement (in sphere diameters) for runs sheared from
# the same prepared sample with different tangential displacement models (MultiStep, OneStep, None), and a table of
# peak and residual ratios, dilation, sliding fraction and coordination at a given shear displacement.
#
# Usage:
#   python ds_history.py --out <file without extension> <run folder> [<run folder> ...] [--at 0.5]
#                        [--dilation <file without extension>]
#
# =============================================================================

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dsplot  # noqa: E402
from dsplot import plt  # noqa: E402

LABELS = {"multistep": "MultiStep (stored tangential displacement)",
          "onestep": "OneStep (displacement over one step)",
          "none": "None (no tangential spring)"}
STYLES = {"multistep": (dsplot.DEM_ORANGE, "-"), "onestep": ("#7F3A0B", "--"), "none": ("#E58B4C", ":")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tau", default="tau_bin")
    ap.add_argument("--at", type=float, default=0.5, help="shear displacement (diameters) for the table")
    ap.add_argument("--size", nargs=2, type=float, default=[6.19, 3.40])
    ap.add_argument("--dilation", default=None, help="also plot the lid rise against shear displacement")
    args = ap.parse_args()

    dsplot.set_style()
    fig, ax = plt.subplots(figsize=args.size)
    rows = []
    for r in args.runs:
        meta = dsplot.read_meta(r)
        ts = dsplot.read_timeseries(r)
        mode = meta["history_shear"]
        d = 2 * meta["radius"]
        idx = dsplot.shear_part(ts)
        x = [ts["shear_disp"][i] / d for i in idx]
        y = [ts[args.tau][i] / meta["sigma"] for i in idx]
        c, ls = STYLES.get(mode, (dsplot.DEM_ORANGE, "-"))
        ax.plot(x, y, color=c, ls=ls, lw=1.4, label=LABELS.get(mode, mode))
        s = dsplot.summarize(ts, meta, args.tau)
        k = min(idx, key=lambda i: abs(ts["shear_disp"][i] / d - args.at))
        rows.append((mode, s, d, ts["frac_sliding"][k], ts["coord_mech"][k], ts["dilation"][k] / d,
                     ts["dilation"][idx[-1]] / d))
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2.5)
    ax.set_xlabel("shear displacement / sphere diameter")
    ax.set_ylabel("shear stress / normal stress")
    ax.grid(color=dsplot.LIGHT_GRAY, lw=0.5)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    dsplot.save(fig, args.out)

    if args.dilation:
        fig, ax = plt.subplots(figsize=args.size)
        for r in args.runs:
            meta = dsplot.read_meta(r)
            ts = dsplot.read_timeseries(r)
            mode = meta["history_shear"]
            d = 2 * meta["radius"]
            idx = dsplot.shear_part(ts)
            c, ls = STYLES.get(mode, (dsplot.DEM_ORANGE, "-"))
            ax.plot([ts["shear_disp"][i] / d for i in idx], [ts["dilation"][i] / d for i in idx], color=c, ls=ls, lw=1.4,
                    label=LABELS.get(mode, mode))
        ax.axhline(0, color=dsplot.LIGHT_GRAY, lw=0.6)
        ax.set_xlim(0, 2)
        ax.set_xlabel("shear displacement / sphere diameter")
        ax.set_ylabel("lid rise / sphere diameter")
        ax.grid(color=dsplot.LIGHT_GRAY, lw=0.5)
        ax.legend(loc="upper left", frameon=False)
        fig.tight_layout()
        dsplot.save(fig, args.dilation)

    lines = ["| model | peak tau/sigma_n | at (d) | residual tau/sigma_n | sliding fraction at %.2g d | "
             "mech. coordination at %.2g d | lid rise at %.2g d (d) | lid rise at end (d) |" % ((args.at,) * 3),
             "|---|---|---|---|---|---|---|---|"]
    for mode, s, d, fs, zc, dil, dil_end in rows:
        lines.append("| %s | %.2f | %.2f | %.2f | %.3f | %.2f | %.3f | %.3f |" % (
            mode, s["peak_ratio"], s["peak_disp"] / d, s["residual_ratio"], fs, zc, dil, dil_end))
    table = "\n".join(lines)
    with open(args.out + "_table.md", "w", encoding="utf-8") as f:
        f.write(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
