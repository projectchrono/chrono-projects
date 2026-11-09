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
# Coarse-grained fields of an annular shear run (test_VDEM_annularShear), averaged
# over the azimuth and over the snapshots of the steady shear (nominal strain from
# --strain_min on), on an (r, y) grid, in cylindrical components (r, phi, y) about
# the axis (cg.py; Gaussian kernel of standard deviation w, renormalized by its
# mass between the side walls):
#   solid fraction, velocity, shear rate (the deviatoric strain-rate invariant
#   sqrt(2 D':D'), which is |dV_phi/dy| in simple shear), the pressure p and
#   tau = sqrt(s:s/2) of the total (contact plus kinetic) stress, mu = tau/p and
#   I = rate d / sqrt(p / rho_s);
#   their values in cells of the bulk, away from the plates (margin_plate
#   diameters from the centers of the glued layers) and the side walls
#   (margin_wall diameters).
# Writes <run>/cg/w<width>_cells.csv (one row per bulk cell) and
# <run>/cg/w<width>_fields.npz (the fields on the grid). Needs numpy only.
#
# Usage:
#   python cg_annular.py --run <dir> [--width 1] [--strain_min 2] [--grid 1e-3]
#                        [--cell 1] [--margin_plate 2] [--margin_wall 2]
#   width, cell and margins in mean bead diameters; grid in m.
#
# =============================================================================

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cg  # noqa: E402


def strain_rate(V, grid_r, grid_y):
    """Deviatoric strain-rate invariant sqrt(2 D':D') of an axisymmetric velocity field V = (V_r, V_phi, V_y)."""
    Vr, Vp, Vy = V
    R = grid_r[:, None]
    D = np.zeros((3, 3) + Vr.shape)
    D[0, 0] = np.gradient(Vr, grid_r, axis=0)
    D[1, 1] = Vr / R
    D[2, 2] = np.gradient(Vy, grid_y, axis=1)
    D[0, 2] = D[2, 0] = 0.5 * (np.gradient(Vr, grid_y, axis=1) + np.gradient(Vy, grid_r, axis=0))
    D[0, 1] = D[1, 0] = 0.5 * R * np.gradient(Vp / R, grid_r, axis=0)
    D[1, 2] = D[2, 1] = 0.5 * np.gradient(Vp, grid_y, axis=1)
    tr = (D[0, 0] + D[1, 1] + D[2, 2]) / 3
    for a in range(3):
        D[a, a] -= tr
    return np.sqrt(2 * np.einsum("ab...,ab...->...", D, D))


def snapshot_value(comment, key):
    """Value of 'key = <number>' in the comment line of a snapshot."""
    text = " ".join(comment)
    return float(text.split(key + " =")[1].split()[0].rstrip(","))


def main():
    ap = argparse.ArgumentParser(description="Coarse-grained fields of an annular shear run")
    ap.add_argument("--run", required=True)
    ap.add_argument("--width", type=float, default=1.0)
    ap.add_argument("--strain_min", type=float, default=2.0)
    ap.add_argument("--grid", type=float, default=1e-3)
    ap.add_argument("--cell", type=float, default=1.0)
    ap.add_argument("--margin_plate", type=float, default=2.0)
    ap.add_argument("--margin_wall", type=float, default=2.0)
    args = ap.parse_args()

    run = args.run
    meta = cg.read_meta(run)
    density = meta["density"]
    first = meta["first_particle_index"]
    d = meta["mean_diameter"]
    r_in, r_out = meta["r_in"], meta["r_out"]
    r_g = meta["glued_radius"]
    w = args.width * d

    names = sorted(f[:-4] for f in os.listdir(os.path.join(run, "particles")) if f.startswith("snap_g"))
    use = [n for n in names if float(n[6:]) >= args.strain_min - 1e-9]
    if not use:
        sys.exit("No snapshots at strain >= %g in %s" % (args.strain_min, run))

    grid_r = np.arange(r_in, r_out + 0.5 * args.grid, args.grid)
    sum_part = None
    sum_cont = None
    faces = []
    for name in use:
        s = cg.Snapshot(run, name, density, first)
        y_face = snapshot_value(s.comment, "y_face")
        faces.append(y_face)
        if sum_part is None:
            grid_y = np.arange(0.0, y_face + 2 * d, args.grid)
        ppos, pw = s.particle_sources_cyl()
        cpos, cw = s.contact_sources_cyl()
        part = cg.splat_axisym(ppos, pw, w, grid_r, grid_y)
        cont = cg.splat_axisym(cpos, cw, w, grid_r, grid_y)
        sum_part = part if sum_part is None else sum_part + part
        sum_cont = cont if sum_cont is None else sum_cont + cont
        print("  %s: y_face %.5f m, %d beads, %d contacts" % (name, y_face, len(s.r), len(s.c_fn)), flush=True)
    n = len(use)
    norm = np.repeat(cg.interval_mass(grid_r, r_in, r_out, w)[:, None], len(grid_y), axis=1)
    F = cg.fields(sum_part / n, sum_cont / n, density, np.maximum(norm, 1e-3))

    rate = strain_rate(F["V"], grid_r, grid_y)
    stot = F["sc"] + F["sk"]
    p, tau = cg.invariants(stot)
    pc, _ = cg.invariants(F["sc"])
    pk, _ = cg.invariants(F["sk"])
    s_py = 0.5 * (stot[1, 2] + stot[2, 1])
    s_yy = stot[2, 2]
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = tau / p
        I = rate * d / np.sqrt(np.where(p > 0, p, np.nan) / density)

    # Bulk cells
    y_lo = r_g + args.margin_plate * d
    y_hi = min(faces) - r_g - args.margin_plate * d
    rr_lo = r_in + args.margin_wall * d
    rr_hi = r_out - args.margin_wall * d
    cell = args.cell * d
    r_edges = np.arange(rr_lo, rr_hi + 1e-9, cell)
    y_edges = np.arange(y_lo, y_hi + 1e-9, cell)
    rows = []
    for i in range(len(r_edges) - 1):
        ir = (grid_r >= r_edges[i]) & (grid_r < r_edges[i + 1])
        for j in range(len(y_edges) - 1):
            jy = (grid_y >= y_edges[j]) & (grid_y < y_edges[j + 1])
            sel = np.ix_(ir, jy)
            pm = p[sel].mean()
            taum = tau[sel].mean()
            ratem = rate[sel].mean()
            rows.append(
                [
                    0.5 * (r_edges[i] + r_edges[i + 1]),
                    0.5 * (y_edges[j] + y_edges[j + 1]),
                    F["nu"][sel].mean(),
                    pm,
                    taum,
                    taum / pm,
                    ratem,
                    ratem * d / math.sqrt(pm / density),
                    s_py[sel].mean(),
                    s_yy[sel].mean(),
                    s_py[sel].mean() / s_yy[sel].mean(),
                    pk[sel].mean() / pc[sel].mean(),
                    F["V"][1][sel].mean(),
                ]
            )
    rows = np.array(rows)

    tag = ("w%gd" % args.width).replace(".", "p")
    out_dir = os.path.join(run, "cg")
    os.makedirs(out_dir, exist_ok=True)
    header = "r,y,nu,p,tau,mu,rate,I,s_phiy,s_yy,s_phiy_over_s_yy,pk_over_pc,V_phi"
    comment = (
        "# %s: %d snapshots from strain %g (%s to %s); width %g d (d = %.5g m); bulk cells %g d, %g d from the glued layers "
        "and %g d from the side walls\n" % (os.path.basename(os.path.normpath(run)), n, args.strain_min, use[0], use[-1], args.width, d, args.cell, args.margin_plate,
                                            args.margin_wall)
    )
    with open(os.path.join(out_dir, tag + "_cells.csv"), "w") as f:
        f.write(comment)
        f.write(header + "\n")
        for row in rows:
            f.write(",".join("%.8g" % v for v in row) + "\n")
    np.savez_compressed(os.path.join(out_dir, tag + "_fields.npz"), grid_r=grid_r, grid_y=grid_y, nu=F["nu"], V=F["V"], rate=rate, p=p, tau=tau, mu=mu, I=I, s_phiy=s_py,
                        s_yy=s_yy, pk=pk, pc=pc, faces=np.array(faces), snapshots=np.array(use), width=w, d=d, y_lo=y_lo, y_hi=y_hi, r_lo=rr_lo, r_hi=rr_hi)

    # Summary: bulk averages per radial bin
    print("%s, width %g d, %d snapshots; bulk y %.4f to %.4f m" % (run, args.width, n, y_lo, y_hi))
    print("     r        I       mu      nu      p    tau/p over y-cells (min..max)")
    for rc in np.unique(rows[:, 0]):
        sel = rows[:, 0] == rc
        print("  %.4f  %.2e  %.3f  %.3f  %6.0f  %.3f..%.3f" % (rc, np.mean(rows[sel, 7]), np.mean(rows[sel, 4]) / np.mean(rows[sel, 3]), np.mean(rows[sel, 2]), np.mean(rows[sel, 3]),
                                                           rows[sel, 5].min(), rows[sel, 5].max()))


if __name__ == "__main__":
    main()
