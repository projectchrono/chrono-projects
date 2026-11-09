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
# Coarse-grained fields of one snapshot of the square-box direct shear
# (test_VDEM_directShear_squareBox: the lower box moves along +z, the upper box is
# fixed, Y up), at several coarse-graining widths w (cg.py; Gaussian kernel of
# standard deviation w, renormalized by its mass inside the box near the walls):
#   maps on the vertical plane x = x0, seen from -x as the clips: solid fraction,
#   contact pressure, shear rate |dVz/dy + dVy/dz|;
#   profiles: averages over horizontal planes of the box against height: solid
#   fraction, velocity, shear rate, tau/p, inertial number;
#   the contact and kinetic parts of the stress (profiles, log scale);
#   a table of the values at the center of the shear band, for one or more
#   snapshots and widths.
# The stresses and the kinetic stress use the velocities of the snapshot. The
# velocity field and shear rate use them too, unless --states names a directory
# of saved states (snap_dXX.XXmm.bin, as written by --snap_states or
# --frame_states): then the velocities are the displacements between the states
# --state_offset before and after the snapshot, over that time (the mean motion,
# without the rattling of single steps).
#
# Usage:
#   python cg_fields.py --run <run dir> --snap <name> --out <file prefix>
#                       [--states <dir>] [--state_offset 0.04]
#                       [--map_widths 0.5,2] [--table_widths 0.5,1,2]
#                       [--table_snaps <name>,<name>,...]
#   widths in mean particle diameters; the run directory holds meta.json,
#   particles/<name>.csv and contacts/<name>.csv.
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
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cg  # noqa: E402
import dsplot  # noqa: E402

# Sequential ramp of the DEM orange family, from white
ORANGE_MAP = LinearSegmentedColormap.from_list("dem_orange", ["#FFFFFF"] + dsplot.DEM_RAMP)
ORANGE_MAP.set_bad("#FFFFFF")


def width_tag(wd):
    return "d2" if abs(wd - 0.5) < 1e-9 else ("%gd" % wd).replace(".", "p")


def width_label(wd):
    if abs(wd - 0.5) < 1e-9:
        return "w = d/2"
    if abs(wd - 1) < 1e-9:
        return "w = d"
    return "w = %gd" % wd


class Box:
    """Material domain of the square box: the upper box (fixed) above the middle of the gap, the lower box (moved by
    the shear displacement along +z) below it."""

    def __init__(self, Lx, Lz, y_floor, y_lid, disp, y_split):
        self.Lx, self.Lz, self.y_floor, self.y_lid, self.disp, self.y_split = Lx, Lz, y_floor, y_lid, disp, y_split

    def z_range(self, y):
        shift = np.where(y < self.y_split, self.disp, 0.0)
        return -self.Lz / 2 + shift, self.Lz / 2 + shift

    def inside(self, y, z):
        z_lo, z_hi = self.z_range(y)
        return (y >= self.y_floor) & (y <= self.y_lid) & (z >= z_lo) & (z <= z_hi)

    def map_norm(self, ys, zs, x0, w):
        Y, Z = np.meshgrid(ys, zs, indexing="ij")
        z_lo, z_hi = self.z_range(Y)
        fx = cg.interval_mass(x0, -self.Lx / 2, self.Lx / 2, w)
        fy = cg.interval_mass(Y, self.y_floor, self.y_lid, w)
        fz = cg.interval_mass(Z, z_lo, z_hi, w)
        return np.maximum(fx * fy * fz, 1e-3)

    def profile_norm(self, ys, w):
        return np.maximum(cg.interval_mass(ys, self.y_floor, self.y_lid, w), 1e-3)


def compute_map(snap, density, w, ys, zs, x0, box, nq=8):
    ppos, pw = snap.particle_sources()
    cpos, cw = snap.contact_sources(nq)
    part = cg.splat_plane(ppos, pw, w, ys, zs, 0, x0, (1, 2))
    cont = cg.splat_plane(cpos, cw, w, ys, zs, 0, x0, (1, 2))
    return cg.fields(part, cont, density, box.map_norm(ys, zs, x0, w))


def compute_profile(snap, density, w, ys, box, nq=8):
    ppos, pw = snap.particle_sources()
    cpos, cw = snap.contact_sources(nq)
    area = box.Lx * box.Lz
    part = cg.splat_profile(ppos, pw, w, ys, 1, area)
    cont = cg.splat_profile(cpos, cw, w, ys, 1, area)
    return cg.fields(part, cont, density, box.profile_norm(ys, w))


def band_values(prof, ys, d, density, y_lo, y_hi):
    """Values at the center of the shear band (largest |dVz/dy| between y_lo and y_hi) of a profile."""
    rate = np.abs(np.gradient(prof["Vm"][2], ys))
    inside = (ys > y_lo) & (ys < y_hi)
    k = int(np.argmax(np.where(inside, rate, -1)))
    half = rate[k] / 2
    lo = k
    while lo > 0 and rate[lo] > half and inside[lo]:
        lo -= 1
    hi = k
    while hi < len(ys) - 1 and rate[hi] > half and inside[hi]:
        hi += 1
    s = prof["sc"] + prof["sk"]
    p, tau = cg.invariants(s)
    pc, _ = cg.invariants(prof["sc"])
    pk, _ = cg.invariants(prof["sk"])
    return {
        "y": ys[k],
        "width": ys[hi] - ys[lo],
        "nu": prof["nu"][k],
        "rate": rate[k],
        "p": p[k],
        "tau": tau[k],
        "mu": tau[k] / p[k],
        "syz_syy": 0.5 * (s[1, 2][k] + s[2, 1][k]) / s[1, 1][k],
        "I": rate[k] * d / math.sqrt(p[k] / density),
        "pk_pc": pk[k] / pc[k],
    }


def load_snapshot(args, name, density, first):
    snap = cg.Snapshot(args.run, name, density, first)
    disp = float(snap.comment[0].split("shear displacement =")[1].split()[0])
    snap.disp = disp
    snap.v_source = "snapshot velocities"
    if args.states:
        mm = 1e3 * disp
        before = "snap_d%05.2fmm" % (mm - args.state_offset)
        after = "snap_d%05.2fmm" % (mm + args.state_offset)
        if os.path.exists(os.path.join(args.states, before + ".bin")) and os.path.exists(os.path.join(args.states, after + ".bin")):
            snap.set_mean_velocity(args.states, before, after, first)
            snap.v_source = "displacements between %s and %s (%.3g s)" % (before, after, snap.v_mean_interval)
    return snap


def main():
    ap = argparse.ArgumentParser(description="Coarse-grained fields of a square-box direct shear snapshot")
    ap.add_argument("--run", required=True)
    ap.add_argument("--snap", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--states", default="")
    ap.add_argument("--state_offset", type=float, default=0.04, help="mm")
    ap.add_argument("--map_widths", default="0.5,2")
    ap.add_argument("--table_widths", default="0.5,1,2")
    ap.add_argument("--table_snaps", default="")
    ap.add_argument("--x0", type=float, default=0.0)
    ap.add_argument("--grid", type=float, default=0.5e-3)
    args = ap.parse_args()

    meta = cg.read_meta(args.run)
    density = meta["density"]
    first = meta["first_particle_index"]
    d = 2 * meta["radius"]
    Lx = meta["box_length"]
    Lz = meta["box_width"]
    gap = meta["gap_between_boxes"]
    map_widths = [float(s) for s in args.map_widths.split(",")]
    table_widths = [float(s) for s in args.table_widths.split(",")]

    snap = load_snapshot(args, args.snap, density, first)
    y_floor = snap.x[:, 1].min() - snap.r.max()
    y_lid = (snap.x[:, 1] + snap.r).max()
    box = Box(Lx, Lz, y_floor, y_lid, snap.disp, gap / 2)
    ys = np.arange(y_floor - 0.002, y_lid + 0.002, args.grid)
    zs = np.arange(-Lz / 2 - 0.002, Lz / 2 + snap.disp + 0.002, args.grid)
    yp = np.arange(y_floor, y_lid, 0.25e-3)

    maps = {wd: compute_map(snap, density, wd * d, ys, zs, args.x0, box) for wd in map_widths}
    profs = {wd: compute_profile(snap, density, wd * d, yp, box) for wd in sorted(set(map_widths + table_widths))}

    dsplot.set_style()

    # Map fields
    Y, Z = np.meshgrid(ys, zs, indexing="ij")
    outside = ~box.inside(Y, Z)

    def map_fields(F):
        pc, _ = cg.invariants(F["sc"])
        Vm = F["Vm"]
        rate = np.abs(np.gradient(Vm[2], ys, axis=0) + np.gradient(Vm[1], zs, axis=1))
        return {
            "nu": np.ma.array(F["nu"], mask=outside),
            "pressure": np.ma.array(1e-3 * pc, mask=outside),
            "rate": np.ma.array(rate, mask=outside),
        }

    mf = {wd: map_fields(maps[wd]) for wd in map_widths}
    labels = {"nu": "solid fraction", "pressure": "contact pressure (kPa)", "rate": "shear rate (1/s)"}
    wmax = max(map_widths)
    limits = {
        "nu": (0.45, 0.70),
        "pressure": (0.0, 2.0 * float(np.percentile(mf[wmax]["pressure"].compressed(), 99))),
        "rate": (0.0, 1.5 * float(np.percentile(mf[wmax]["rate"].compressed(), 99.5))),
    }
    extent = [1e3 * zs[0], 1e3 * zs[-1], 1e3 * ys[0], 1e3 * ys[-1]]

    def draw_walls(ax):
        # Upper box (fixed), lower box (moved by the shear displacement along +z), the lid
        gray = dsplot.EXP_GRAY
        z1, z2, dz = -Lz / 2, Lz / 2, snap.disp
        ax.plot(1e3 * np.array([z1, z1, z2, z2]), 1e3 * np.array([gap, y_lid, y_lid, gap]), color=gray, lw=0.8)
        ax.plot(1e3 * np.array([z1 + dz, z1 + dz, z2 + dz, z2 + dz]), 1e3 * np.array([0.0, y_floor, y_floor, 0.0]), color=gray, lw=0.8)

    def map_panel(ax, wd, key):
        vmin, vmax = limits[key]
        im = ax.imshow(mf[wd][key], origin="lower", extent=extent, cmap=ORANGE_MAP, vmin=vmin, vmax=vmax, interpolation="nearest", aspect="equal")
        draw_walls(ax)
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("y (mm)")
        cb = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
        cb.set_label(labels[key])

    fig, axes = plt.subplots(len(map_widths), 3, figsize=(10.5, 2.15 * len(map_widths) + 0.3), squeeze=False)
    for i, wd in enumerate(map_widths):
        for j, key in enumerate(["nu", "pressure", "rate"]):
            map_panel(axes[i, j], wd, key)
            if j > 0:
                axes[i, j].set_ylabel("")
        axes[i, 0].text(-0.32, 0.5, width_label(wd), transform=axes[i, 0].transAxes, rotation=90, va="center", ha="center", fontsize=11)
    fig.tight_layout()
    dsplot.save(fig, args.out)
    plt.close(fig)

    for wd in map_widths:
        for key in ["nu", "pressure", "rate"]:
            fig, ax = plt.subplots(figsize=(4.47, 2.3))
            map_panel(ax, wd, key)
            fig.tight_layout()
            dsplot.save(fig, "%s_%s_w%s" % (args.out, key, width_tag(wd)))
            plt.close(fig)

    # Profiles over horizontal planes of the box
    shades = [dsplot.DEM_RAMP[1], dsplot.DEM_RAMP[3]]
    colors = {wd: shades[min(i, 1)] for i, wd in enumerate(sorted(map_widths))}
    fig, axes = plt.subplots(1, 5, figsize=(10.5, 3.2), sharey=True)
    for wd in map_widths:
        P = profs[wd]
        c = colors[wd]
        Vz = P["Vm"][2]
        rate = np.abs(np.gradient(Vz, yp))
        p, tau = cg.invariants(P["sc"] + P["sk"])
        with np.errstate(invalid="ignore", divide="ignore"):
            I = rate * d / np.sqrt(np.where(p > 0, p, np.nan) / density)
            mu = tau / np.where(p > 0, p, np.nan)
        yy = 1e3 * yp
        axes[0].plot(P["nu"], yy, color=c, lw=1.2, label=width_label(wd))
        axes[1].plot(1e3 * Vz, yy, color=c, lw=1.2)
        axes[2].plot(rate, yy, color=c, lw=1.2)
        axes[3].plot(mu, yy, color=c, lw=1.2)
        axes[4].plot(I, yy, color=c, lw=1.2)
    axes[0].set_xlabel("solid fraction")
    axes[1].set_xlabel("velocity $V_z$ (mm/s)")
    axes[2].set_xlabel("shear rate (1/s)")
    axes[3].set_xlabel(r"$\tau / p$")
    axes[4].set_xlabel("inertial number $I$")
    axes[4].set_xscale("log")
    axes[0].set_ylabel("y (mm)")
    axes[0].legend(frameon=False, loc="lower left")
    for ax in axes:
        ax.axhspan(0, 1e3 * gap, color=dsplot.LIGHT_GRAY, alpha=0.4, lw=0)
    fig.tight_layout()
    dsplot.save(fig, args.out + "_profiles")
    plt.close(fig)

    # Contact and kinetic parts of the stress, side by side
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.2), sharey=True)
    for wd in map_widths:
        P = profs[wd]
        c = colors[wd]
        pc, tauc = cg.invariants(P["sc"])
        pk, tauk = cg.invariants(P["sk"])
        yy = 1e3 * yp
        axes[0].plot(pc, yy, color=c, lw=1.2, label="contact, " + width_label(wd))
        axes[0].plot(np.abs(pk), yy, color=c, lw=1.2, ls="--", label="kinetic, " + width_label(wd))
        axes[1].plot(tauc, yy, color=c, lw=1.2)
        axes[1].plot(tauk, yy, color=c, lw=1.2, ls="--")
    for ax in axes:
        ax.set_xscale("log")
        ax.axhspan(0, 1e3 * gap, color=dsplot.LIGHT_GRAY, alpha=0.4, lw=0)
    axes[0].set_xlabel("pressure (Pa)")
    axes[1].set_xlabel(r"shear stress $\tau$ (Pa)")
    axes[0].set_ylabel("y (mm)")
    axes[0].legend(frameon=False, fontsize=8, loc="center")
    fig.tight_layout()
    dsplot.save(fig, args.out + "_stress_parts")
    plt.close(fig)

    # Table: values at the center of the shear band
    table_snaps = [s for s in args.table_snaps.split(",") if s] or [args.snap]
    lines = [
        "# Coarse-grained values at the center of the shear band",
        "",
        "Run `%s`. Averages over horizontal planes of the box (%.3f x %.3f m), Gaussian kernel of standard deviation w "
        "(renormalized inside the box). Band center: the height of the largest |dVz/dy|; band width: where the shear rate "
        "is above half of that value. p and tau: invariants of the total (contact plus kinetic) stress; "
        "I = rate d / sqrt(p / rho_s), d = %.4f m." % (args.run, Lx, Lz, d),
        "",
        "| snapshot | velocities | width | band center y (mm) | band width (mm) | solid fraction | shear rate (1/s) | p (Pa) | "
        "tau (Pa) | tau/p | s_yz/s_yy | I | kinetic/contact pressure |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for name in table_snaps:
        s = snap if name == args.snap else load_snapshot(args, name, density, first)
        for wd in table_widths:
            P = profs[wd] if name == args.snap else compute_profile(s, density, wd * d, yp, box)
            b = band_values(P, yp, d, density, y_floor + d, y_lid - d)
            lines.append(
                "| %s | %s | %s | %.1f | %.1f | %.3f | %.3g | %.0f | %.0f | %.3f | %.3f | %.2e | %.1e |"
                % (name, s.v_source, width_label(wd), 1e3 * b["y"], 1e3 * b["width"], b["nu"], b["rate"], b["p"], b["tau"], b["mu"], b["syz_syy"], b["I"], b["pk_pc"])
            )
    with open(args.out + "_table.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
