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
# mu(I) from the annular shear runs (test_VDEM_annularShear, coarse-grained by
# cg_annular.py): the bulk cells of each run, averaged over the height in each
# radial bin, give points (I, mu = tau/p, solid fraction). A law is fitted to all
# points by least squares, and to the points of each pressure for the table:
#   --law jop:    mu(I) = mu_s + (mu_2 - mu_s) / (1 + I_0 / I)
#                 (Jop, Forterre and Pouliquen 2006);
#   --law linear: mu(I) = mu_s + b I, the limit of the Jop law for I << I_0
#                 (b = (mu_2 - mu_s) / I_0), drawn over the range of the points.
# Optionally adds the quasi-static point of the direct shear band (cg_fields.py
# table) and the glass-bead law of Jop et al. (2006) for reference.
#
# Colors carry roles: DEM results in orange, the fitted law in purple,
# literature in gray.
#
# Usage:
#   python mu_of_I.py --out <file prefix> [--law jop|linear] [--width w1d]
#                     [--qs_table <cg_fields_table.md>] [--literature] <run dir> [<run dir> ...]
#
# =============================================================================

import argparse
import math
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cg  # noqa: E402
import dsplot  # noqa: E402

FIT_PURPLE = "#6B3FA0"  # laws and fits derived from the DEM results
MARKERS = ["o", "s", "^", "D", "v"]

# Jop, Forterre and Pouliquen (2006), glass beads on a rough inclined plane
JOP2006 = (math.tan(math.radians(20.9)), math.tan(math.radians(32.76)), 0.279)


def mu_law(I, mu_s, mu_2, I0):
    return mu_s + (mu_2 - mu_s) / (1 + I0 / I)


def linear_law(I, mu_s, b):
    return mu_s + b * I


LAWS = {
    "jop": (mu_law, [0.3, 0.6, 0.3], ([0, 0, 1e-4], [2, 3, 10]), ["mu_s", "mu_2", "I_0"]),
    "linear": (linear_law, [0.3, 0.5], ([0, -10], [2, 10]), ["mu_s", "b"]),
}


def fit(law, I, mu):
    """Least-squares fit of a law to the points: parameters, their standard errors and the rms residual."""
    f, p0, bounds, _ = LAWS[law]
    popt, pcov = curve_fit(f, I, mu, p0=p0, bounds=bounds, maxfev=20000)
    return popt, np.sqrt(np.diag(pcov)), math.sqrt(np.mean((mu - f(I, *popt)) ** 2))


def read_cells(path):
    with open(path) as f:
        lines = [line for line in f if not line.startswith("#")]
    names = lines[0].strip().split(",")
    arr = np.loadtxt(lines[1:], delimiter=",", ndmin=2)
    return {n: arr[:, k] for k, n in enumerate(names)}


def run_points(run, width):
    """Points of one run: the bulk cells averaged over the height in each radial bin."""
    c = read_cells(os.path.join(run, "cg", width + "_cells.csv"))
    meta = cg.read_meta(run)
    d = meta["mean_diameter"]
    rho_s = meta["density"]
    pts = []
    for r in np.unique(c["r"]):
        sel = c["r"] == r
        p = c["p"][sel].mean()
        tau = c["tau"][sel].mean()
        rate = c["rate"][sel].mean()
        pts.append(
            {
                "r": r,
                "p": p,
                "tau": tau,
                "mu": tau / p,
                "mu_spread": (c["mu"][sel].min(), c["mu"][sel].max()),
                "I": rate * d / math.sqrt(p / rho_s),
                "nu": c["nu"][sel].mean(),
                "s_ratio": c["s_phiy"][sel].mean() / c["s_yy"][sel].mean(),
                "kin": c["pk_over_pc"][sel].mean(),
            }
        )
    return meta, pts


def read_qs_table(path):
    """Quasi-static points from the cg_fields.py table: (I, tau/p) of every snapshot and width."""
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            cells = [s.strip() for s in line.strip().strip("|").split("|")]
            if len(cells) == 13 and cells[0].startswith("snap_"):
                rows.append((cells[2], float(cells[11]), float(cells[9])))
    return rows


def main():
    ap = argparse.ArgumentParser(description="mu(I) from the annular shear runs")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--law", choices=sorted(LAWS), default="jop")
    ap.add_argument("--width", default="w1d")
    ap.add_argument("--qs_table", default="")
    ap.add_argument("--literature", action="store_true")
    ap.add_argument("--size", type=float, nargs=2, default=[7.40, 2.47])
    args = ap.parse_args()

    runs = []
    for run in args.runs:
        meta, pts = run_points(run, args.width)
        runs.append((run, meta, pts))
    pressures = sorted({m["pressure"] for _, m, _ in runs})

    I_all = np.array([pt["I"] for _, _, pts in runs for pt in pts])
    mu_all = np.array([pt["mu"] for _, _, pts in runs for pt in pts])
    p_all = np.array([m["pressure"] for _, m, pts in runs for pt in pts])
    law = LAWS[args.law][0]
    popt, perr, rms = fit(args.law, I_all, mu_all)
    fits = [("all", len(I_all), popt, perr, rms)]
    for pr in pressures:
        sel = p_all == pr
        fits.append(("%g kPa" % (1e-3 * pr), int(sel.sum())) + fit(args.law, I_all[sel], mu_all[sel]))

    qs = read_qs_table(args.qs_table) if args.qs_table else []

    dsplot.set_style()
    I_line = np.logspace(math.log10(min(I_all.min(), 1e-4) / 1.5), math.log10(max(I_all.max(), 0.5) * 1.5), 200)
    # The linear law is the small-I limit of the Jop law: drawn only up to the largest I of the points
    I_fit = I_line if args.law == "jop" else I_line[I_line <= I_all.max()]
    fit_label = {"jop": r"fit: $\mu_s$ = %.3f, $\mu_2$ = %.3f, $I_0$ = %.3f",
                 "linear": r"fit: $\mu$ = %.3f + %.3f $I$"}[args.law] % tuple(popt)

    def mu_panel(ax, legend=True):
        for k, pr in enumerate(pressures):
            pts = [pt for _, m, pp in runs if m["pressure"] == pr for pt in pp]
            ax.plot([pt["I"] for pt in pts], [pt["mu"] for pt in pts], MARKERS[k % len(MARKERS)], color=dsplot.DEM_ORANGE, mfc=dsplot.DEM_ORANGE if k == 0 else "white",
                    ms=4, mew=1.0, lw=0, label="DEM, annular shear, %g kPa" % (1e-3 * pr))
        ax.plot(I_fit, law(I_fit, *popt), color=FIT_PURPLE, lw=1.5, label=fit_label)
        if qs:
            I_q = [r[1] for r in qs]
            mu_q = [r[2] for r in qs]
            mid = [r for r in qs if r[0] == "w = d"] or qs
            Im = np.mean([r[1] for r in mid])
            mum = np.mean([r[2] for r in mid])
            ax.errorbar([Im], [mum], xerr=[[Im - min(I_q)], [max(I_q) - Im]], yerr=[[mum - min(mu_q)], [max(mu_q) - mum]], fmt="D", color=dsplot.DEM_ORANGE, mfc="white", ms=5,
                        mew=1.0, elinewidth=0.8, capsize=2, label="DEM, direct shear band, 3.1 kPa")
        if args.literature:
            ax.plot(I_line, mu_law(I_line, *JOP2006), color=dsplot.EXP_GRAY, lw=1.0, ls="--", label="Jop et al. (2006), glass beads")
        ax.set_xscale("log")
        ax.set_xlabel("inertial number $I$")
        ax.set_ylabel(r"$\mu = \tau / p$")
        if legend:
            ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=2)

    fig, ax = plt.subplots(figsize=tuple(args.size))
    mu_panel(ax)
    fig.tight_layout()
    dsplot.save(fig, args.out)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.47, 2.98))
    for k, pr in enumerate(pressures):
        pts = [pt for _, m, pp in runs if m["pressure"] == pr for pt in pp]
        ax.plot([pt["I"] for pt in pts], [pt["nu"] for pt in pts], MARKERS[k % len(MARKERS)], color=dsplot.DEM_ORANGE, mfc=dsplot.DEM_ORANGE if k == 0 else "white", ms=4,
                mew=1.0, lw=0, label="%g kPa" % (1e-3 * pr))
    ax.set_xscale("log")
    ax.set_xlabel("inertial number $I$")
    ax.set_ylabel("solid fraction")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    dsplot.save(fig, args.out + "_phi")
    plt.close(fig)

    names = LAWS[args.law][3]
    law_text = {"jop": "mu(I) = mu_s + (mu_2 - mu_s) / (1 + I_0 / I)", "linear": "mu(I) = mu_s + b I"}[args.law]
    lines = [
        "# mu(I) from the annular shear runs",
        "",
        "Fit of %s to the points (bulk cells averaged over the height in each radial bin; coarse-graining %s), all "
        "together and for each pressure; ± one standard error:" % (law_text, args.width),
        "",
        "| points | n | " + " | ".join(names) + " | rms residual |",
        "|---|---|" + "---|" * len(names) + "---|",
    ]
    for label, n, pf, pe, r in fits:
        lines.append("| %s | %d | " % (label, n) + " | ".join("%.3f ± %.3f" % (a, e) for a, e in zip(pf, pe)) + " | %.4f |" % r)
    n_head = len(lines)
    lines += [
        "",
        "| run | pressure (Pa) | I nominal | r (m) | I | mu = tau/p (cells: min..max) | s_phiy/s_yy | solid fraction | kinetic/contact pressure |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for run, m, pts in runs:
        for pt in pts:
            lines.append(
                "| %s | %g | %g | %.4f | %.3g | %.3f (%.3f..%.3f) | %.3f | %.3f | %.1e |"
                % (os.path.basename(os.path.normpath(run)), m["pressure"], m["I_nominal"], pt["r"], pt["I"], pt["mu"], pt["mu_spread"][0], pt["mu_spread"][1], abs(pt["s_ratio"]),
                   pt["nu"], pt["kin"])
            )
    if qs:
        lines += ["", "Quasi-static points of the direct shear band (cg_fields.py): " + "; ".join("%s: I %.2g, tau/p %.3f" % r for r in qs)]
    with open(args.out + "_table.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines[:n_head]))


if __name__ == "__main__":
    main()
