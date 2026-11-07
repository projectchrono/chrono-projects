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
# Critical state line and isotropic compression from triaxial runs, in the two
# planes of Modified Cam Clay:
# - (p, q): the stress paths of drained triaxial tests at several cell
#   pressures, their states at the end of the test, and the critical state line
#   q = M p fitted through those states;
# - (ln p, v), v = 1 + e: the isotropic compression with unload-reload loops,
#   with the slopes of its loading (lambda) and unloading (kappa) branches, and
#   the end states of the triaxial tests with the critical state line in that
#   plane (slope lambda_cs).
# Also writes a table of the end states and the fitted parameters.
#
# Usage:
#   python tri_csl.py --pq <file without extension> --lnp <file without extension>
#                     --iso <isotropic run folder> <triaxial run folder> [...]
#
# Optional:
#   --end_window W   end states: medians over the last W percent of axial strain (default 2)
#   --range P0 P1    pressure range of the fits of lambda and kappa, kPa (default 100 400)
#   --size_pq W H    size of the (p, q) figure in inches (default 6.60 2.18)
#   --size_lnp W H   size of the (ln p, v) figure in inches (default 5.77 2.18)
#
# Stresses on the current area (membrane volume over height); p and q from the
# time series: p = (sigma1 + 2 sigma3)/3, q = sigma1 - sigma3, compression
# positive. The void ratio is that of the volume enclosed by the membrane.
#
# =============================================================================

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dsplot  # noqa: E402
import tri_curves  # noqa: E402
from dsplot import plt  # noqa: E402

FIT_PURPLE = "#6B3FA0"  # laws and fits derived from the DEM results


def end_state(ts, window):
    """Medians of p, q, e over the last 'window' percent of axial strain of the test, and their drift over the last 5 %."""
    idx = tri_curves.test_part(ts)
    ea = np.array([ts["eps_a"][i] for i in idx])
    sel = ea >= ea[-1] - 0.01 * window

    def med(key, mask):
        return float(np.nanmedian(np.array([ts[key][i] for i in idx])[mask]))

    before = (ea >= ea[-1] - 0.05 - 0.01 * window) & (ea < ea[-1] - 0.05)
    return {
        "p": med("p", sel), "q": med("q", sel), "e": med("e", sel), "eta": med("eta", sel),
        "eta_drift": med("eta", sel) - med("eta", before), "e_drift": med("e", sel) - med("e", before),
        "eps_end": ea[-1],
    }


def branches(ts):
    """Ramps of the isotropic path (phase 2) between holds: (rows, +1 loading or -1 unloading)."""
    idx = [i for i, ph in enumerate(ts["phase"]) if ph == 2]
    s3 = np.array([ts["sigma3"][i] for i in idx])
    out = []
    start = 0
    direction = 0
    for k in range(1, len(idx)):
        d = int(np.sign(s3[k] - s3[k - 1]))
        if d != direction:
            if direction != 0:
                out.append((idx[start:k], direction))  # a ramp ends at a hold or a reversal
            start = k - 1
            direction = d
    if direction != 0:
        out.append((idx[start:], direction))
    return out


def slope(p, v, p0, p1):
    """-dv/dln p fitted over p0 <= p <= p1 (least squares), with the number of points."""
    p = np.asarray(p)
    v = np.asarray(v)
    sel = (p >= p0) & (p <= p1)
    if sel.sum() < 3:
        return math.nan, int(sel.sum())
    a, _ = np.polyfit(np.log(p[sel]), v[sel], 1)
    return -a, int(sel.sum())


def main():
    ap = argparse.ArgumentParser(description="Critical state line and isotropic compression")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--iso", required=True)
    ap.add_argument("--pq", required=True)
    ap.add_argument("--lnp", required=True)
    ap.add_argument("--end_window", type=float, default=2.0)
    ap.add_argument("--range", nargs=2, type=float, default=[100.0, 400.0])
    ap.add_argument("--size_pq", nargs=2, type=float, default=[6.60, 2.18])
    ap.add_argument("--size_lnp", nargs=2, type=float, default=[5.77, 2.18])
    args = ap.parse_args()

    dsplot.set_style()

    # Triaxial runs, by cell pressure
    runs = []
    for r in args.runs:
        meta = dsplot.read_meta(r)
        ts = dsplot.read_timeseries(r)
        runs.append((meta["sigma3"], r, meta, ts, end_state(ts, args.end_window)))
    runs.sort(key=lambda x: x[0])
    colors = tri_curves.ramp(len(runs))
    p_end = np.array([s["p"] for *_, s in runs])
    q_end = np.array([s["q"] for *_, s in runs])
    v_end = 1 + np.array([s["e"] for *_, s in runs])
    M = float(np.sum(p_end * q_end) / np.sum(p_end * p_end))
    sin_cs = 3 * M / (6 + M)
    phi_cs = math.degrees(math.asin(sin_cs))
    lam_cs, gamma = np.polyfit(np.log(1e-3 * p_end), v_end, 1)
    lam_cs = -lam_cs

    # (p, q): stress paths, end states, critical state line
    fig, ax = plt.subplots(figsize=args.size_pq, layout="constrained")
    for c, (s3, r, meta, ts, s) in zip(colors, runs):
        idx = tri_curves.test_part(ts)
        ax.plot([1e-3 * ts["p"][i] for i in idx], [1e-3 * ts["q"][i] for i in idx], color=c, linewidth=0.6, alpha=0.8)
    ax.plot([], [], color=colors[len(colors) // 2], linewidth=0.8,
            label=r"DEM stress paths, $\sigma_3$ = %s kPa" % ", ".join("%g" % (1e-3 * x[0]) for x in runs))
    ax.scatter(1e-3 * p_end, 1e-3 * q_end, s=22, color=colors, edgecolor="black", linewidth=0.5, zorder=3,
               label="DEM, at %g %% axial strain" % round(100 * runs[0][4]["eps_end"]))
    pmax = 1.15 * max(1e-3 * max(ts["p"][i] for i in tri_curves.test_part(ts)) for *_, ts, _ in runs)
    ax.plot([0, pmax], [0, M * pmax], color=FIT_PURPLE, linewidth=1.4,
            label=r"Critical state line, $M$ = %.2f ($\varphi_{cs}$ = %.1f°)" % (M, phi_cs))
    ax.set_xlim(0, pmax)
    ax.set_ylim(0, None)
    ax.set_xlabel(r"Mean stress $p$ (kPa)")
    ax.set_ylabel(r"Deviator stress $q$ (kPa)")
    ax.legend(loc="upper left", frameon=False)
    dsplot.save(fig, args.pq)

    # Isotropic compression: branches and slopes
    meta_iso = dsplot.read_meta(args.iso)
    ts_iso = dsplot.read_timeseries(args.iso)
    p_iso = 1e-3 * np.array(ts_iso["p"])
    v_iso = 1 + np.array(ts_iso["e"])
    p0, p1 = args.range
    fits = []  # (kind, p from, p to, slope over the branch, slope within [p0, p1])
    p_max = 0.0
    for rows, direction in branches(ts_iso):
        p = p_iso[rows]
        v = v_iso[rows]
        lo, hi = float(p.min()), float(p.max())
        if direction > 0 and hi > 1.05 * p_max:
            # loading beyond the largest pressure reached so far (the first loading included): only that part
            kind = "first loading" if p_max == 0 else "loading beyond the previous maximum"
            virgin = p >= p_max
            p, v = p[virgin], v[virgin]
            lo = max(lo, p_max)
        else:
            kind = "reloading" if direction > 0 else "unloading"
        k_all = slope(p, v, lo, hi)[0]
        a, b = max(lo, p0), min(hi, p1)
        k_rng = slope(p, v, a, b)[0] if b >= 1.5 * a else math.nan
        fits.append((kind, lo, hi, k_all, k_rng, direction))
        p_max = max(p_max, hi)
    lam = [f[4] for f in fits if f[0] == "loading beyond the previous maximum" and not math.isnan(f[4])]
    kap = [f[4] for f in fits if f[0] == "unloading" and not math.isnan(f[4])]
    lam_v = lam[0] if lam else math.nan
    kap_v = float(np.mean(kap)) if kap else math.nan

    # (ln p, v)
    fig, ax = plt.subplots(figsize=args.size_lnp, layout="constrained")
    iso = [i for i, ph in enumerate(ts_iso["phase"]) if ph == 2]
    ax.plot(p_iso[iso], v_iso[iso], color=dsplot.DEM_ORANGE, linewidth=1.0,
            label="DEM, isotropic compression with unload-reload loops")
    ax.scatter(1e-3 * p_end, v_end, s=22, color=colors, edgecolor="black", linewidth=0.5, zorder=3,
               label="DEM triaxial tests at %g %% axial strain" % round(100 * runs[0][4]["eps_end"]))
    pp = np.geomspace(0.8 * min(1e-3 * p_end), 1.25 * max(1e-3 * p_end), 50)
    ax.plot(pp, gamma - lam_cs * np.log(pp), color=FIT_PURPLE, linewidth=1.2,
            label=r"Critical state line, $\lambda_{cs}$ = %.3f" % lam_cs)
    ax.set_xscale("log")
    ticks = [t for t in (5, 10, 20, 50, 100, 200, 500) if 0.7 * min(p_iso[iso]) <= t <= 1.3 * max(p_iso[iso])]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["%g" % t for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel(r"Mean stress $p$ (kPa, log scale)")
    ax.set_ylabel(r"Specific volume $v = 1 + e$")
    ax.legend(loc="lower left", frameon=False)
    ax.text(0.98, 0.95, r"$-dv/d\ln p$, %g to %g kPa: loading %.3f, unloading %.3f" % (p0, p1, lam_v, kap_v),
            transform=ax.transAxes, ha="right", va="top", color=dsplot.DEM_ORANGE)
    dsplot.save(fig, args.lnp)

    # Table
    with open(args.pq + "_table.md", "w", encoding="utf-8") as f:
        f.write("| run | σ3 (kPa) | p at end (kPa) | q at end (kPa) | q/p at end | change of q/p over the last 5 % | e at end | "
                "change of e over the last 5 % |\n|---|---|---|---|---|---|---|---|\n")
        for s3, r, meta, ts, s in runs:
            f.write("| %s | %g | %.1f | %.1f | %.3f | %+.3f | %.4f | %+.4f |\n" % (
                os.path.basename(os.path.normpath(r)), 1e-3 * s3, 1e-3 * s["p"], 1e-3 * s["q"], s["eta"], s["eta_drift"],
                s["e"], s["e_drift"]))
        f.write("\nCritical state line through the end states (medians over the last %g %% of axial strain, at %g %%): "
                "q = M p with M = %.3f (least squares through the origin), sin φ_cs = 3M/(6 + M) = %.3f, φ_cs = %.1f°; "
                "CRM's μ = τ̄/p (τ̄ = sqrt(½ τ:τ) = q/√3 in triaxial compression) = M/√3 = %.3f. In (ln p, v): v = Γ − "
                "λ_cs ln p with λ_cs = %.4f, Γ = %.4f (p in kPa).\n\n" % (
                    args.end_window, round(100 * runs[0][4]["eps_end"]), M, sin_cs, phi_cs, M / math.sqrt(3), lam_cs, gamma))
        f.write("Isotropic compression (%s): path %s kPa.\n\n| branch | p from (kPa) | p to (kPa) | −dv/dln p over the "
                "branch | −dv/dln p between %g and %g kPa |\n|---|---|---|---|---|\n" % (
                    os.path.basename(os.path.normpath(args.iso)), ", ".join("%g" % (1e-3 * x) for x in meta_iso["path"]),
                    p0, p1))
        for kind, lo, hi, k_all, k_rng, direction in fits:
            a, b = (lo, hi) if direction > 0 else (hi, lo)
            f.write("| %s | %.1f | %.1f | %.4f | %s |\n" % (kind, a, b, k_all, "" if math.isnan(k_rng) else "%.4f" % k_rng))
        before = [i for i, ph in enumerate(ts_iso["phase"]) if ph == 1][-1]
        f.write("\nλ = %.4f from the loading beyond the previous maximum, κ = %.4f from the unloading branches (mean), both "
                "between %g and %g kPa. Void ratio at %g kPa: %.4f before the loops, %.4f after them.\n" % (
                    lam_v, kap_v, p0, p1, 1e-3 * meta_iso["path"][0], v_iso[before] - 1, v_iso[iso[-1]] - 1))
    print("M = %.3f (phi_cs %.1f deg, mu = M/sqrt3 = %.3f), lambda_cs = %.4f, lambda = %.4f, kappa = %.4f" % (
        M, phi_cs, M / math.sqrt(3), lam_cs, lam_v, kap_v))


if __name__ == "__main__":
    main()
