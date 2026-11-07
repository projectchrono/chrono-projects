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
# Love-Weber stress of triaxial specimens from their particle and contact
# snapshots, checked against the boundary stresses:
# - over the whole specimen: sigma = -(1/V) sum over spheres p, sum over their
#   contacts c, of (x_c - x_p) (x) f_c (compression positive), V the volume
#   enclosed by the membrane. Its axial component must equal sigma1, the mean of
#   the top and bottom platen stresses (the mean includes the specimen's
#   weight). The membrane forces act at the sphere centers, so the average
#   misses their moment arm from the center to the membrane: its lateral
#   component (sigma_xx + sigma_zz)/2 must equal sigma3 (1 - r / R), with r the
#   mean sphere radius and R the membrane radius;
# - in a central cylinder (half the membrane radius, the middle half of the
#   height), from the sphere-sphere contacts whose midpoint lies in it:
#   sigma = (1/V_c) sum l (x) f (l from b1 to b2, f the force on b2), in
#   cylindrical components, against sigma1 and sigma3.
#
# Usage:
#   python tri_love_weber.py --out <file.md> <run folder> [<run folder> ...]
#
# =============================================================================

import argparse
import csv
import json
import math
import os
import re

import numpy as np


def read_meta(run):
    with open(os.path.join(run, "meta.json")) as f:
        return json.load(f)


def read_timeseries(run):
    with open(os.path.join(run, "timeseries.csv")) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    return {k: np.array([float(r[k]) if r[k] not in ("", "nan") else math.nan for r in rows]) for k in rows[0]}


def read_csv(path):
    with open(path) as f:
        header = [line for line in f if line.startswith("#")]
    with open(path) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    data = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    m = re.search(r"time = ([0-9.eE+-]+)", "".join(header))
    return data, float(m.group(1)) if m else math.nan


def check(run, snap, meta, ts):
    P, t = read_csv(os.path.join(run, "particles", snap + ".csv"))
    C, _ = read_csv(os.path.join(run, "contacts", snap + ".csv"))
    first = meta["first_particle_index"]
    row = int(np.argmin(np.abs(np.array(ts["time"]) - t)))
    V = ts["volume"][row]
    H = ts["height"][row]
    R = ts["radius"][row]
    s1 = ts["sigma1"][row]
    s3 = ts["sigma3"][row]
    r_mean = float(np.mean(P["r"]))

    pos = np.zeros((int(P["id"].max()) + 1, 3))
    pos[P["id"].astype(int)] = np.stack([P["x"], P["y"], P["z"]], 1)
    b1 = C["b1"].astype(int)
    b2 = C["b2"].astype(int)
    xc = np.stack([C["px"], C["py"], C["pz"]], 1)
    F = C["fn"][:, None] * np.stack([C["nx"], C["ny"], C["nz"]], 1) + np.stack([C["ftx"], C["fty"], C["ftz"]], 1)

    # Whole specimen: each sphere's contacts, the force on b2 and its opposite on b1
    S = np.zeros((3, 3))
    on2 = b2 >= first
    on1 = b1 >= first
    S += np.einsum("ij,ik->jk", xc[on2] - pos[b2[on2]], F[on2])
    S += np.einsum("ij,ik->jk", xc[on1] - pos[b1[on1]], -F[on1])
    sig = -S / V
    axial = sig[1, 1]
    lateral = 0.5 * (sig[0, 0] + sig[2, 2])

    # Central cylinder: sphere-sphere contacts with their midpoint inside
    pp = on1 & on2
    L = pos[b2[pp]] - pos[b1[pp]]
    mid = 0.5 * (pos[b2[pp]] + pos[b1[pp]])
    Fp = F[pp]
    rho = np.hypot(mid[:, 0], mid[:, 2])
    sel = (rho < 0.5 * R) & (mid[:, 1] > 0.25 * H) & (mid[:, 1] < 0.75 * H)
    Vc = math.pi * (0.5 * R) ** 2 * 0.5 * H
    er = np.stack([mid[sel, 0], np.zeros(sel.sum()), mid[sel, 2]], 1) / np.maximum(rho[sel], 1e-12)[:, None]
    et = np.stack([-er[:, 2], np.zeros(sel.sum()), er[:, 0]], 1)
    core_rr = np.sum(np.einsum("ij,ij->i", Fp[sel], er) * np.einsum("ij,ij->i", L[sel], er)) / Vc
    core_tt = np.sum(np.einsum("ij,ij->i", Fp[sel], et) * np.einsum("ij,ij->i", L[sel], et)) / Vc
    core_yy = np.sum(Fp[sel, 1] * L[sel, 1]) / Vc

    return {
        "eps_a": ts["eps_a"][row], "s1": s1, "s3": s3, "axial": axial, "lateral": lateral,
        "lateral_expected": s3 * (1 - r_mean / R), "core_yy": core_yy, "core_rr": core_rr, "core_tt": core_tt,
        "core_contacts": int(sel.sum()), "ke": ts["kinetic_energy"][row],
    }


def main():
    ap = argparse.ArgumentParser(description="Love-Weber stress of triaxial specimens against the boundary stresses")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("| run | snapshot | ε_a (%) | σ1 (kPa) | axial LW / σ1 | σ3 (kPa) | lateral LW / σ3(1 − r/R) | core axial / σ1 | "
                "core radial / σ3 | core hoop / σ3 | core contacts |\n|---|---|---|---|---|---|---|---|---|---|---|\n")
        worst = 0.0
        for run in args.runs:
            meta = read_meta(run)
            ts = read_timeseries(run)
            snaps = sorted(os.path.splitext(n)[0] for n in os.listdir(os.path.join(run, "contacts")) if n.startswith("snap_"))
            for snap in snaps:
                c = check(run, snap, meta, ts)
                a = c["axial"] / c["s1"]
                b = c["lateral"] / c["lateral_expected"]
                worst = max(worst, abs(a - 1), abs(b - 1))
                f.write("| %s | %s | %.2f | %.1f | %.3f | %.1f | %.3f | %.3f | %.3f | %.3f | %d |\n" % (
                    os.path.basename(os.path.normpath(run)), snap, 100 * c["eps_a"], 1e-3 * c["s1"], a, 1e-3 * c["s3"], b,
                    c["core_yy"] / c["s1"], c["core_rr"] / c["s3"], c["core_tt"] / c["s3"], c["core_contacts"]))
                print("%s %s: axial %.3f, lateral %.3f, core %.3f %.3f %.3f" % (
                    os.path.basename(os.path.normpath(run)), snap, a, b, c["core_yy"] / c["s1"], c["core_rr"] / c["s3"],
                    c["core_tt"] / c["s3"]))
        f.write("\nWhole specimen: σ = −(1/V) Σ_p Σ_c (x_c − x_p) ⊗ f_c over every sphere's contacts (platen contacts "
                "included), V the volume enclosed by the membrane; axial component against σ1 (the mean of the top and "
                "bottom platen stresses on the current area), lateral component (σ_xx + σ_zz)/2 against σ3 (1 − r/R) (the "
                "membrane forces act at the sphere centers, so the average misses their moment arm r). Core: the "
                "sphere-sphere contacts whose midpoint lies within half the membrane radius and the middle half of the "
                "height. Largest deviation of the whole-specimen ratios from 1: %.1f %%.\n" % (100 * worst))


if __name__ == "__main__":
    main()
