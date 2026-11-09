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
# Coarse-graining of particle and contact snapshots (Goldhirsch 2010; Weinhart et
# al. 2012 for the contacts with walls), with a Gaussian kernel phi of standard
# deviation w (the coarse-graining width):
#   density         rho(r) = sum_i m_i phi(r - r_i)
#   solid fraction  nu(r)  = sum_i V_i phi(r - r_i)
#   momentum        j(r)   = sum_i m_i v_i phi(r - r_i),   velocity V = j / rho
#   kinetic stress  sk(r)  = sum_i m_i v_i v_i phi(r - r_i) - rho V V
#   contact stress  sc(r)  = sum_c f_c l_c int_0^1 phi(r - a_c - s l_c) ds
# Stresses are positive in compression. For a contact between two particles, f_c
# is the force on the second particle and l_c the branch vector from the first
# center (a_c) to the second; for a contact with a wall, f_c is the force on the
# particle and l_c runs from the contact point (a_c) to the particle center. The
# line integral is a midpoint rule with nq points along the branch.
#
# Fields are evaluated on a plane (maps), averaged over planes normal to an axis
# (profiles: the kernel integrated over the plane, divided by its area), or
# averaged over the azimuth about the y axis (fields on an (r, y) grid: the kernel
# in the meridian plane divided by 2 pi r, with vectors and tensors in cylindrical
# components).
#
# =============================================================================

import json
import math
import os
import struct

import numpy as np

SQRT_2PI = np.sqrt(2 * np.pi)
_erf = np.vectorize(math.erf)


def gauss(u, w):
    return np.exp(-0.5 * (u / w) ** 2) / (SQRT_2PI * w)


def interval_mass(c, lo, hi, w):
    """Fraction of a Gaussian of standard deviation w centered at c that lies in [lo, hi] (boundary correction)."""
    c = np.asarray(c, dtype=float)
    return 0.5 * (_erf((hi - c) / (math.sqrt(2) * w)) - _erf((lo - c) / (math.sqrt(2) * w)))


def read_state_positions(path, first_particle):
    """Time and particle positions (bodies from first_particle on) of a state file written by dem_utils.h SaveState."""
    with open(path, "rb") as f:
        if f.read(8) != b"CHDEMST1":
            raise ValueError("not a state file: " + path)
        time = struct.unpack("d", f.read(8))[0]
        n = struct.unpack("q", f.read(8))[0]
        data = np.frombuffer(f.read(13 * 8 * n), dtype=np.float64).reshape(n, 13)
    return time, data[first_particle:, 0:3].copy()


def read_csv(path):
    """Read a snapshot CSV into a dict of float arrays; also return its comment lines."""
    with open(path) as f:
        lines = f.readlines()
    comment = [line[1:].strip() for line in lines if line.startswith("#")]
    data = [line for line in lines if not line.startswith("#")]
    names = data[0].strip().split(",")
    arr = np.loadtxt(data[1:], delimiter=",", ndmin=2)
    return {n: arr[:, k] for k, n in enumerate(names)}, comment


def read_meta(run_dir):
    with open(os.path.join(run_dir, "meta.json")) as f:
        return json.load(f)


class Snapshot:
    """Particles and contacts of one snapshot, with the contacts as branch segments."""

    def __init__(self, run_dir, name, density, first_particle):
        P, self.comment = read_csv(os.path.join(run_dir, "particles", name + ".csv"))
        C, _ = read_csv(os.path.join(run_dir, "contacts", name + ".csv"))
        self.name = name
        ids = P["id"].astype(int)
        self.r = P["r"]
        self.x = np.column_stack([P["x"], P["y"], P["z"]])
        self.v = np.column_stack([P["vx"], P["vy"], P["vz"]])
        self.w = np.column_stack([P["wx"], P["wy"], P["wz"]])
        self.vol = 4.0 / 3.0 * np.pi * self.r ** 3
        self.m = density * self.vol
        self.ids = ids
        self.v_mean = None  # optional velocities of the mean motion (e.g. from displacements between two states)

        lookup = -np.ones(ids.max() + 1, dtype=int)
        lookup[ids] = np.arange(len(ids))

        def row(b):
            out = np.full(len(b), -1)
            ok = (b >= first_particle) & (b < len(lookup))
            out[ok] = lookup[b[ok]]
            return out

        i1 = row(C["b1"].astype(int))
        i2 = row(C["b2"].astype(int))
        n = np.column_stack([C["nx"], C["ny"], C["nz"]])
        F = C["fn"][:, None] * n + np.column_stack([C["ftx"], C["fty"], C["ftz"]])  # force on body 2
        pc = np.column_stack([C["px"], C["py"], C["pz"]])

        pair = (i1 >= 0) & (i2 >= 0)
        wall2 = (i1 < 0) & (i2 >= 0)  # wall is body 1, particle body 2
        wall1 = (i1 >= 0) & (i2 < 0)  # particle body 1, wall body 2
        a = np.concatenate([self.x[i1[pair]], pc[wall2], pc[wall1]])
        l = np.concatenate([self.x[i2[pair]] - self.x[i1[pair]], self.x[i2[wall2]] - pc[wall2], self.x[i1[wall1]] - pc[wall1]])
        f = np.concatenate([F[pair], F[wall2], -F[wall1]])
        self.c_start = a
        self.c_branch = l
        self.c_force = f
        self.c_is_wall = np.concatenate([np.zeros(pair.sum(), bool), np.ones(wall2.sum() + wall1.sum(), bool)])
        self.c_fn = np.concatenate([C["fn"][pair], C["fn"][wall2], C["fn"][wall1]])
        self.c_mu = np.concatenate([C["mu"][pair], C["mu"][wall2], C["mu"][wall1]])
        self.c_ratio = np.concatenate([C["ratio"][pair], C["ratio"][wall2], C["ratio"][wall1]])
        self.num_pairs = int(pair.sum())
        self.num_wall_contacts = int(wall2.sum() + wall1.sum())

    def set_mean_velocity(self, states_dir, before, after, first_particle):
        """Velocities of the mean motion: the displacements between two saved states over their time interval."""
        t0, x0 = read_state_positions(os.path.join(states_dir, before + ".bin"), first_particle)
        t1, x1 = read_state_positions(os.path.join(states_dir, after + ".bin"), first_particle)
        rows = self.ids - first_particle
        self.v_mean = (x1[rows] - x0[rows]) / (t1 - t0)
        self.v_mean_interval = t1 - t0

    def particle_sources(self):
        """Positions and weights of the particle sums: mass, volume, momentum (3), second moment of velocity (9), and
        the momentum of the mean motion (3) if set."""
        m, v = self.m, self.v
        wts = [m, self.vol, m * v[:, 0], m * v[:, 1], m * v[:, 2]]
        for a in range(3):
            for b in range(3):
                wts.append(m * v[:, a] * v[:, b])
        if self.v_mean is not None:
            wts += [m * self.v_mean[:, 0], m * self.v_mean[:, 1], m * self.v_mean[:, 2]]
        return self.x, np.column_stack(wts)

    def contact_sources(self, nq=8):
        """Positions and weights of the contact stress sums (9 components f_a l_b), nq points per branch."""
        fl = (self.c_force[:, :, None] * self.c_branch[:, None, :]).reshape(-1, 9) / nq
        pos = np.concatenate([self.c_start + (q + 0.5) / nq * self.c_branch for q in range(nq)])
        return pos, np.tile(fl, (nq, 1))

    def particle_sources_cyl(self):
        """As particle_sources, with the velocities in cylindrical components (r, phi, y) about the y axis."""
        m = self.m
        v = cylindrical(self.x, self.v)
        wts = [m, self.vol, m * v[:, 0], m * v[:, 1], m * v[:, 2]]
        for a in range(3):
            for b in range(3):
                wts.append(m * v[:, a] * v[:, b])
        return self.x, np.column_stack(wts)

    def contact_sources_cyl(self, nq=8):
        """As contact_sources, with force and branch in cylindrical components (r, phi, y) at each point."""
        pos_all, wts_all = [], []
        for q in range(nq):
            pos = self.c_start + (q + 0.5) / nq * self.c_branch
            f = cylindrical(pos, self.c_force)
            l = cylindrical(pos, self.c_branch)
            pos_all.append(pos)
            wts_all.append((f[:, :, None] * l[:, None, :]).reshape(-1, 9) / nq)
        return np.concatenate(pos_all), np.concatenate(wts_all)


def splat_plane(pos, weights, w, grid_u, grid_v, axis_n, x0, axes_uv, chunk=20000):
    """Sum of weights phi(r - pos) at the nodes (u, v) of the plane x[axis_n] = x0. Returns (K, len(u), len(v))."""
    sel = np.abs(pos[:, axis_n] - x0) < 5 * w
    pos = pos[sel]
    weights = weights[sel]
    gn = gauss(pos[:, axis_n] - x0, w)
    out = np.zeros((weights.shape[1], len(grid_u), len(grid_v)))
    for s in range(0, len(pos), chunk):
        p = pos[s:s + chunk]
        Au = gauss(grid_u[None, :] - p[:, axes_uv[0]][:, None], w)
        Av = gauss(grid_v[None, :] - p[:, axes_uv[1]][:, None], w)
        g = gn[s:s + chunk]
        for k in range(weights.shape[1]):
            out[k] += (Au * (weights[s:s + chunk, k] * g)[:, None]).T @ Av
    return out


def splat_profile(pos, weights, w, grid, axis, area, chunk=50000):
    """Plane averages of sum weights phi(r - pos) over planes normal to axis, of the given area. Returns (K, len(grid))."""
    out = np.zeros((weights.shape[1], len(grid)))
    for s in range(0, len(pos), chunk):
        A = gauss(grid[None, :] - pos[s:s + chunk, axis][:, None], w)
        out += weights[s:s + chunk].T @ A
    return out / area


def cylindrical(pos, vec):
    """Components (r, phi, y) of vectors vec at positions pos, about the y axis; phi increases with a positive
    rotation about +y (e_phi = e_y x e_r)."""
    th = np.arctan2(pos[:, 2], pos[:, 0])
    c, s = np.cos(th), np.sin(th)
    vr = vec[:, 0] * c + vec[:, 2] * s
    vp = vec[:, 0] * s - vec[:, 2] * c
    return np.column_stack([vr, vp, vec[:, 1]])


def splat_axisym(pos, weights, w, grid_r, grid_y, chunk=20000):
    """Azimuthal averages of sum weights phi(r - pos) on an (r, y) grid about the y axis: the kernel in the meridian
    plane divided by the circumference 2 pi r (exact as w / r goes to zero). Returns (K, len(grid_r), len(grid_y))."""
    rho = np.hypot(pos[:, 0], pos[:, 2])
    out = np.zeros((weights.shape[1], len(grid_r), len(grid_y)))
    for s in range(0, len(pos), chunk):
        Ar = gauss(grid_r[None, :] - rho[s:s + chunk, None], w)
        Ay = gauss(grid_y[None, :] - pos[s:s + chunk, 1][:, None], w)
        for k in range(weights.shape[1]):
            out[k] += (Ar * weights[s:s + chunk, k][:, None]).T @ Ay
    return out / (2 * np.pi * grid_r[None, :, None])


def fields(part, cont, density, norm=None):
    """Derived fields from the particle sums (14 rows, 17 with the mean motion) and the contact sums (9 rows), on any
    grid shape. norm, if given, is the kernel mass inside the material domain at each node (boundary correction)."""
    if norm is not None:
        part = part / norm[None]
        cont = cont / norm[None]
    rho = part[0]
    nu = part[1]
    safe = np.where(rho > 0, rho, 1)
    V = np.where(rho > 0, part[2:5] / safe, 0.0)
    M2 = part[5:14].reshape((3, 3) + rho.shape)
    sk = M2 - rho[None, None] * V[:, None] * V[None, :]
    sc = cont.reshape((3, 3) + rho.shape)
    out = {"rho": rho, "nu": nu, "V": V, "sk": sk, "sc": sc, "density": density}
    out["Vm"] = np.where(rho > 0, part[14:17] / safe, 0.0) if part.shape[0] >= 17 else V
    return out


def invariants(s):
    """Pressure p = tr(s)/3 and tau = sqrt(dev(s):dev(s)/2) of the symmetric part of a stress field s (3, 3, ...)."""
    sym = 0.5 * (s + np.swapaxes(s, 0, 1))
    p = (sym[0, 0] + sym[1, 1] + sym[2, 2]) / 3
    dev = sym.copy()
    for a in range(3):
        dev[a, a] = dev[a, a] - p
    tau = np.sqrt(0.5 * np.einsum("ab...,ab...->...", dev, dev))
    return p, tau
