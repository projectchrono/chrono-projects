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
# Shared helpers for the direct shear post-processing scripts: reading run
# folders (meta.json, timeseries.csv) and the figure style of the talk figures.
#
# =============================================================================

import csv
import json
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Colors carry roles in the talk deck: DEM orange (with an orange-to-brown ramp for several DEM series),
# experiments in black and grays. Never a saturated red.
DEM_ORANGE = "#C55A11"
DEM_RAMP = ["#F4B183", "#E58B4C", "#C55A11", "#7F3A0B"]
EXP_BLACK = "#000000"
EXP_GRAY = "#737373"
LIGHT_GRAY = "#BFBFBF"


def set_style():
    """Plot style of the talk figures: white background, Calibri, no title."""
    plt.rcParams.update({
        "font.family": "Calibri",
        "font.size": 10,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 0.8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
    })


def save(fig, path_no_ext):
    """Save a figure as PNG (300 dpi) and SVG."""
    fig.savefig(path_no_ext + ".png", dpi=300)
    fig.savefig(path_no_ext + ".svg")


def read_meta(run_dir):
    with open(os.path.join(run_dir, "meta.json")) as f:
        return json.load(f)


def read_timeseries(run_dir):
    """Return the time series as a dict of float lists (NaN for missing values)."""
    with open(os.path.join(run_dir, "timeseries.csv")) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    data = {}
    for key in rows[0].keys():
        col = []
        for r in rows:
            try:
                col.append(float(r[key]))
            except ValueError:
                col.append(math.nan)
        data[key] = col
    return data


def shear_part(ts):
    """Indices of the rows in the shear phase."""
    return [i for i, ph in enumerate(ts["phase"]) if ph == 2]


def median(values):
    v = sorted(x for x in values if not math.isnan(x))
    if not v:
        return math.nan
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def summarize(ts, meta, tau_key="tau_bin", residual_window=0.002, plateau=(0.004, 0.0075)):
    """Peak and residual stress ratios and the initial void ratio of a direct shear run.

    The residual is the median over the last 'residual_window' of shear (Suhr and Six's convention). The plateau is
    the mean over the 'plateau' range of shear displacement (by default 4 to 7.5 mm, where the measured curves of
    Hartl and Ooi (2008) are level). The initial void ratio is the lid value at the start of shearing.
    """
    idx = shear_part(ts)
    sigma = meta["sigma"]
    disp = [ts["shear_disp"][i] for i in idx]
    tau = [ts[tau_key][i] for i in idx]
    ratio = [t / sigma for t in tau]
    i_peak = max(range(len(ratio)), key=lambda k: ratio[k])
    d_end = disp[-1]
    resid = median([r for d, r in zip(disp, ratio) if d >= d_end - residual_window])
    plat = [r for d, r in zip(disp, ratio) if plateau[0] <= d <= plateau[1]]
    plat = sum(plat) / len(plat) if plat else math.nan
    e0 = ts["e_lid"][idx[0]]
    e_end = ts["e_lid"][idx[-1]]
    return {
        "sigma": sigma,
        "peak_ratio": ratio[i_peak],
        "peak_disp": disp[i_peak],
        "plateau_ratio": plat,
        "plateau_tau": plat * sigma,
        "residual_ratio": resid,
        "residual_tau": resid * sigma,
        "e0": e0,
        "porosity0": e0 / (1 + e0),
        "e_end": e_end,
        "travel": d_end,
        "dilation_end": ts["dilation"][idx[-1]],
        "overflows": meta.get("history_overflows_total", 0),
    }
