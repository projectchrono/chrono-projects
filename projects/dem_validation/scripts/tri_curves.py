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
# Curves of triaxial runs against axial strain: the stress ratio
# (sigma1 - sigma3)/(sigma1 + sigma3), the dilation, q/p or the void ratio, next
# to the laboratory tests of Cui, O'Sullivan and O'Neill (2007), Fig. 4. Runs
# that differ only in their seed are drawn as their mean curve with a band over
# their range. Also writes a table of void ratios, peak and end values.
#
# Usage:
#   python tri_curves.py --out <file without extension> <run folder> [<run folder> ...]
#
# Optional:
#   --panels LIST       quantities, one panel each, side by side: ratio (stress ratio), dilation (-eps_v),
#                       eta (q/p) or e (void ratio); default ratio
#   --by_specimen       one panel per specimen (mono, poly) for each quantity; the quantities then go in rows
#   --area0             axial stress on the area at the start of the test (no area correction, as in the
#                       laboratory); default: on the current area (the membrane volume over the height)
#   --lab FILE          laboratory curves, CSV with columns specimen, e0, eps_a_pc, stress_ratio (drawn with --area0);
#                       data/cui2007_fig4_lab.csv holds those of Cui et al. (2007), Fig. 4
#   --labels TEXT...    legend labels of the groups of runs, in the order of the table (default: from the metadata)
#   --xmax X            axial strain range in percent (default: the longest run)
#   --end_window W      the end values are medians over the last W percent of axial strain (default 2)
#   --size W H          figure size in inches (default 4.47 2.98)
#   --legend LOC        legend location (default: lower right)
#
# The void ratios in the legends are those of the specimen as deposited, on its
# nominal volume (mold radius times the height just after the top platen is
# placed): the laboratory computed its void ratios this way. Compression is
# positive (eps_v > 0 contracts).
#
# =============================================================================

import argparse
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dsplot  # noqa: E402
from dsplot import plt  # noqa: E402

SPECIMEN_NAMES = {"mono": "Monodisperse", "poly": "Polydisperse"}
HISTORY_NAMES = {"multistep": "MultiStep", "onestep": "OneStep", "none": "None"}
HISTORY_ORDER = {"multistep": 0, "onestep": 1, "none": 2}


def specimen_of(meta):
    """'mono' for a single sphere size, 'poly' otherwise."""
    return "mono" if len(meta["radii"]) == 1 else "poly"


def group_of(meta):
    """Runs that differ only in their seed form one group."""
    return (specimen_of(meta), meta["mu_fill"], meta["sigma3"], meta["strain_rate"], meta["history_test"],
            meta["membrane_model"])


def nominal_void_ratio(ts, meta):
    """Void ratio on the nominal volume as deposited: mold radius and the platen height just after placement.

    The laboratory computed its void ratios this way, from the specimen's nominal dimensions and the sphere count.
    """
    R = 0.5 * meta["diameter"]
    for i, ph in enumerate(ts["phase"]):
        if ph >= 1 and not math.isnan(ts["height"][i]):
            # first row with the platen resting on the sample: skip the row at placement itself
            j = min(i + 1, len(ts["phase"]) - 1)
            return math.pi * R * R * ts["height"][j] / meta["solid_volume"] - 1
    return math.nan


def test_part(ts):
    return [i for i, ph in enumerate(ts["phase"]) if ph == 2]


def summarize(ts, meta, xmax, end_window):
    idx = test_part(ts)
    ea = np.array([ts["eps_a"][i] for i in idx])
    cur = np.array([ts["sin_phi"][i] for i in idx])
    a0 = np.array([ts["sin_phi_a0"][i] for i in idx])
    k = int(np.nanargmax(cur))
    in_range = ea <= 0.01 * xmax + 1e-9
    k0 = int(np.nanargmax(np.where(in_range, a0, -np.inf)))
    end = ea >= ea[-1] - 0.01 * end_window

    def end_value(key):
        return float(np.nanmedian([ts[key][i] for i, e in zip(idx, end) if e]))

    return {
        "e_nominal": nominal_void_ratio(ts, meta),
        "e_start": meta.get("void_ratio_test_start", math.nan),
        "peak": cur[k],
        "eps_peak": ea[k],
        "peak_a0": a0[k0],
        "eps_peak_a0": ea[k0],
        "end": float(np.nanmedian(cur[end])),
        "eta_end": end_value("eta"),
        "p_end": end_value("p"),
        "q_end": end_value("q"),
        "e_end": end_value("e"),
        "eps_end": ea[-1],
        "eps_v_end": ts["eps_v"][idx[-1]],
        "overflows": meta.get("history_overflows_total", 0),
    }


def read_lab(path, specimen):
    curves = {}
    with open(path) as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            if r["specimen"] == specimen:
                curves.setdefault(float(r["e0"]), []).append((float(r["eps_a_pc"]), float(r["stress_ratio"])))
    return curves


def ramp(n):
    """n colors along the deck's orange-to-brown ramp of DEM series."""
    if n == 1:
        return [dsplot.DEM_ORANGE]
    rgb = np.array([[int(c[i:i + 2], 16) for i in (1, 3, 5)] for c in dsplot.DEM_RAMP], dtype=float)
    s = np.linspace(0, len(rgb) - 1, n)
    out = []
    for x in s:
        i = min(int(x), len(rgb) - 2)
        c = rgb[i] + (x - i) * (rgb[i + 1] - rgb[i])
        out.append("#%02X%02X%02X" % tuple(int(round(v)) for v in c))
    return out


def main():
    ap = argparse.ArgumentParser(description="Curves of triaxial runs against axial strain")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--panels", default="ratio")
    ap.add_argument("--by_specimen", action="store_true")
    ap.add_argument("--area0", action="store_true")
    ap.add_argument("--lab", default=None)
    ap.add_argument("--labels", nargs="+", default=None)
    ap.add_argument("--xmax", type=float, default=None)
    ap.add_argument("--end_window", type=float, default=2.0)
    ap.add_argument("--size", nargs=2, type=float, default=[4.47, 2.98])
    ap.add_argument("--legend", default="lower right")
    args = ap.parse_args()

    dsplot.set_style()
    quantities = args.panels.split(",")

    # Groups of runs, per specimen and pressure, densest deposition first, then by rate and contact history
    groups = {}
    for r in args.runs:
        meta = dsplot.read_meta(r)
        ts = dsplot.read_timeseries(r)
        groups.setdefault(group_of(meta), []).append([r, meta, ts])
    order = sorted(groups, key=lambda g: (g[0], g[2], g[1], g[3], HISTORY_ORDER.get(g[4], 9), g[5]))
    xmax = args.xmax
    if xmax is None:
        xmax = max(100 * max(ts["eps_a"]) for g in order for _, _, ts in groups[g])
    for g in order:
        for run in groups[g]:
            run.append(summarize(run[2], run[1], xmax, args.end_window))
    specimens = sorted({g[0] for g in order}) if args.by_specimen else [None]
    # colors along the ramp among the groups drawn in the same panel
    colors = {}
    for specimen in specimens:
        panel_groups = [g for g in order if specimen is None or g[0] == specimen]
        colors.update(zip(panel_groups, ramp(len(panel_groups))))
    if args.labels and len(args.labels) != len(order):
        raise SystemExit("--labels needs %d labels (groups: %s)" % (len(order), order))

    def varies(k):
        return len({g[k] for g in order}) > 1

    def label(g):
        """The group's legend: what tells it apart from the other groups, its void ratio as deposited, its seeds."""
        if args.labels:
            return args.labels[order.index(g)]
        s = [x[3] for x in groups[g]]
        parts = []
        if varies(0) and not args.by_specimen:
            parts.append(SPECIMEN_NAMES[g[0]].lower())
        if varies(2):
            parts.append(r"$\sigma_3$ = %g kPa" % (1e-3 * g[2]))
        if varies(1):
            parts.append("deposition friction %g" % g[1])
        if varies(3):
            parts.append("%g %%/s" % (100 * g[3]))
        if varies(4):
            parts.append(HISTORY_NAMES.get(g[4], g[4]))
        if varies(5):
            parts.append("%s membrane" % g[5])
        e_n = [x["e_nominal"] for x in s]
        all_e = [x[3]["e_nominal"] for h in order for x in groups[h]]
        if args.lab or varies(1) or max(all_e) - min(all_e) > 0.005:
            parts.append("e = %.3f" % np.mean(e_n))
        text = "DEM" + (", " + ", ".join(parts) if parts else "")
        return text + ("" if len(s) == 1 else ", %d seeds" % len(s))

    def series(g, column, scale):
        """Common axial strain grid (percent) and the group's curves of 'column' on it."""
        runs = []
        for _, _, ts, _ in groups[g]:
            idx = test_part(ts)
            runs.append((np.array([100 * ts["eps_a"][i] for i in idx]), np.array([scale * ts[column][i] for i in idx])))
        e_end = min(min(ea[-1] for ea, _ in runs), xmax)
        grid = np.linspace(0, e_end, int(round(e_end / 0.01)) + 1)
        return grid, np.array([np.interp(grid, ea, v) for ea, v in runs])

    ratio_key = "sin_phi_a0" if args.area0 else "sin_phi"
    spec = {
        "ratio": (ratio_key, 1.0, r"$(\sigma_1 - \sigma_3)/(\sigma_1 + \sigma_3)$"),
        "dilation": ("eps_v", -100.0, r"Dilation, $-\varepsilon_v$ (%)"),
        "eta": ("eta", 1.0, r"$q/p$"),
        "e": ("e", 1.0, "Void ratio"),
    }

    nrows = len(quantities) if args.by_specimen else 1
    ncols = len(specimens) if args.by_specimen else len(quantities)
    fig, axes = plt.subplots(nrows, ncols, figsize=args.size, layout="constrained", squeeze=False)
    for qi, quantity in enumerate(quantities):
        column, scale, ylabel = spec[quantity]
        for si, specimen in enumerate(specimens):
            ax = axes[qi, si] if args.by_specimen else axes[0, qi]
            for g in order:
                if specimen is not None and g[0] != specimen:
                    continue
                grid, values = series(g, column, scale)
                if len(values) > 1:
                    ax.fill_between(grid, values.min(axis=0), values.max(axis=0), color=colors[g], alpha=0.25, linewidth=0)
                ax.plot(grid, values.mean(axis=0), color=colors[g], linewidth=1.4, label=label(g))
            if quantity == "ratio" and args.lab and args.area0:
                lab_specimens = [specimen] if specimen else sorted({g[0] for g in order})
                for ls in lab_specimens:
                    lab = read_lab(args.lab, ls)
                    for c, e0 in zip([dsplot.EXP_BLACK, dsplot.EXP_GRAY], sorted(lab)):
                        x, y = zip(*sorted(lab[e0]))
                        ax.plot(x, y, color=c, linewidth=1.0, linestyle=(0, (4, 2)),
                                label="Lab (Cui et al. 2007), e = %.3f" % e0)
            if quantity == "ratio":
                ax.set_ylim(0, 0.4)
            if quantity == "dilation":
                ax.axhline(0, color=dsplot.LIGHT_GRAY, linewidth=0.8, zorder=0)
            ax.set_xlim(0, xmax)
            ax.set_xlabel("Axial strain (%)")
            if si == 0 or not args.by_specimen:
                ax.set_ylabel(ylabel)
            if specimen is not None and qi == 0:
                ax.text(0.03, 0.95, SPECIMEN_NAMES[specimen], transform=ax.transAxes, ha="left", va="top")
            if (args.by_specimen and qi == 0) or (not args.by_specimen and qi == 0):
                ax.legend(loc=args.legend, frameon=False)
    dsplot.save(fig, args.out)

    # Table
    with open(args.out + "_table.md", "w", encoding="utf-8") as f:
        f.write("| run | σ3 (kPa) | rate (%%/s) | contact history | e as deposited (nominal) | e at σ3 (membrane) | "
                "peak stress ratio | ε_a at peak (%%) | φ peak (°) | max. without area correction, to %g %% | ε_a there (%%) | "
                "stress ratio at end | φ end (°) | q/p at end | e at end | ε_v at end (%%) | history overflows |\n" % xmax)
        f.write("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
        for g in order:
            for r, meta, ts, s in groups[g]:
                f.write("| %s | %g | %g | %s | %.3f | %.3f | %.3f | %.2f | %.1f | %.3f | %.2f | %.3f | %.1f | %.3f | %.3f | "
                        "%.2f | %d |\n" % (
                            os.path.basename(os.path.normpath(r)), 1e-3 * g[2], 100 * g[3], HISTORY_NAMES.get(g[4], g[4]),
                            s["e_nominal"], s["e_start"], s["peak"], 100 * s["eps_peak"], math.degrees(math.asin(s["peak"])),
                            s["peak_a0"], 100 * s["eps_peak_a0"], s["end"], math.degrees(math.asin(s["end"])), s["eta_end"],
                            s["e_end"], 100 * s["eps_v_end"], s["overflows"]))
        if args.lab:
            for sp in sorted({g[0] for g in order}):
                lab = read_lab(args.lab, sp)
                f.write("\nLab (Cui, O'Sullivan and O'Neill 2007, %s): peak stress ratio %s at e = %s (nominal dimensions, "
                        "no area correction).\n" % (sp, ", ".join("%.3f" % max(y for _, y in v) for _, v in sorted(lab.items())),
                                                    ", ".join("%.3f" % e for e in sorted(lab))))
        f.write("\nStress ratio (σ1 − σ3)/(σ1 + σ3) = sin φ on the current area (membrane volume over height), except the "
                "column without area correction (on the area at the start of the test, as in the laboratory). End values: "
                "medians over the last %g %% of axial strain. Compression positive (ε_v > 0 contracts). e as deposited: "
                "nominal volume (mold radius, height after the top platen is placed), the laboratory's definition.\n"
                % args.end_window)


if __name__ == "__main__":
    main()
