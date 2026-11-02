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
# Force networks side by side: images written by test_VDEM_directShear_squareBox
# --render (one per run, at the same shear displacement), cropped to a common
# frame around their content, with a label under each.
#
# Usage:
#   python ds_networks.py --out <file without extension> --images <img> [<img> ...] --labels <text> [<text> ...]
#
# Optional:
#   --size W H   figure size in inches (default 6.19 2.4)
#   --margin F   margin around the common content box, as a fraction of its size (default 0.03)
#
# =============================================================================

import argparse
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dsplot  # noqa: E402
from dsplot import plt  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Force networks side by side")
    ap.add_argument("--out", required=True)
    ap.add_argument("--images", nargs="+", required=True)
    ap.add_argument("--labels", nargs="+", required=True)
    ap.add_argument("--size", nargs=2, type=float, default=[6.19, 2.4])
    ap.add_argument("--margin", type=float, default=0.03)
    args = ap.parse_args()
    if len(args.images) != len(args.labels):
        raise SystemExit("one label per image")

    images = [np.asarray(Image.open(f).convert("RGB")) for f in args.images]

    # Common crop: the box holding every non-white pixel of every image
    mask = np.zeros(images[0].shape[:2], dtype=bool)
    for im in images:
        mask |= (im < 250).any(axis=2)
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    r0, r1, c0, c1 = rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
    dr = int(args.margin * (r1 - r0))
    dc = int(args.margin * (c1 - c0))
    r0, r1 = max(0, r0 - dr), min(mask.shape[0], r1 + dr)
    c0, c1 = max(0, c0 - dc), min(mask.shape[1], c1 + dc)

    dsplot.set_style()
    fig, axes = plt.subplots(1, len(images), figsize=args.size)
    if len(images) == 1:
        axes = [axes]
    for ax, im, label in zip(axes, images, args.labels):
        ax.imshow(im[r0:r1, c0:c1], interpolation="lanczos")
        ax.set_axis_off()
        ax.set_title(label, y=-0.12, fontsize=10)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.1, wspace=0.04)
    dsplot.save(fig, args.out)
    print("%s.png: %d images, crop rows %d-%d, columns %d-%d" % (args.out, len(images), r0, r1, c0, c1))


if __name__ == "__main__":
    main()
