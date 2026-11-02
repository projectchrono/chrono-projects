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
# Clip and stills from the images written by test_VDEM_directShear_squareBox
# --render: the force network layer composited over the grains layer, frame by
# frame (H.264 MP4, loopable, with a poster frame).
#
# Usage:
#   python ds_render.py --render <render folder> --out <file without extension>
#
# Optional:
#   --fps F           frame rate of the clip (default 30)
#   --hold S          seconds the last frame is held before the clip loops (default 1)
#   --fade S          seconds of cross-fade from the last frame back to the first (default 0.5)
#   --poster STEM     the frame used as poster (default: the last frame)
#   --stills STEM...  also write the composite of these frames (<out>_<stem>.png)
#   --lighten F       blend the grains layer toward white by this fraction (default 0.3)
#   --layers          also write a clip of each layer alone (<out>_grains.mp4, <out>_network.mp4)
#
# Pixels that differ between the network and walls layers belong to the
# network; all other pixels come from the grains layer (lightened, so that the
# network stands out; the walls are lightened too).
#
# =============================================================================

import argparse
import os
import shutil
import subprocess
import tempfile

import numpy as np
from PIL import Image


def load(path):
    return np.asarray(Image.open(path).convert("RGB"), dtype=np.float32)


def composite(render_dir, stem, lighten):
    """Network layer over the lightened grains layer, with the network alpha from its difference to the walls layer."""
    grains = load(os.path.join(render_dir, "grains", stem + ".bmp"))
    network = load(os.path.join(render_dir, "network", stem + ".bmp"))
    walls = load(os.path.join(render_dir, "walls", stem + ".bmp"))
    alpha = np.clip(np.abs(network - walls).max(axis=2, keepdims=True) / 48.0, 0.0, 1.0)
    grains = (1.0 - lighten) * grains + lighten * 255.0
    return grains * (1.0 - alpha) + network * alpha


def save_png(image, path):
    Image.fromarray(np.clip(image + 0.5, 0, 255).astype(np.uint8)).save(path)


def encode(frames, fps, hold, fade, out_mp4):
    """Write a loopable clip: the frames, the last frame held, then a cross-fade back to the first frame."""
    n_hold = int(round(hold * fps))
    n_fade = int(round(fade * fps))
    tmp = tempfile.mkdtemp(prefix="ds_render_")
    try:
        k = 0

        def put(image):
            nonlocal k
            save_png(image, os.path.join(tmp, "%05d.png" % k))
            k += 1

        first = None
        last = None
        for frame in frames:
            last = frame()
            if first is None:
                first = last
            put(last)
        for _ in range(n_hold):
            put(last)
        for i in range(1, n_fade + 1):
            s = i / (n_fade + 1)
            put((1 - s) * last + s * first)
        subprocess.run(["ffmpeg", "-v", "error", "-y", "-framerate", str(fps), "-i", os.path.join(tmp, "%05d.png"),
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", out_mp4],
                       check=True)
        return k
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description="Clip and stills from the render layers of the direct shear program")
    ap.add_argument("--render", required=True, help="render folder (with grains/, network/, walls/)")
    ap.add_argument("--out", required=True, help="output file name without extension")
    ap.add_argument("--fps", type=float, default=30)
    ap.add_argument("--hold", type=float, default=1.0)
    ap.add_argument("--fade", type=float, default=0.5)
    ap.add_argument("--poster", default=None)
    ap.add_argument("--stills", nargs="*", default=[])
    ap.add_argument("--lighten", type=float, default=0.3)
    ap.add_argument("--layers", action="store_true")
    args = ap.parse_args()

    stems = sorted(os.path.splitext(f)[0] for f in os.listdir(os.path.join(args.render, "grains")) if f.endswith(".bmp"))
    if not stems:
        raise SystemExit("No images in " + os.path.join(args.render, "grains"))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    frames = [lambda s=s: composite(args.render, s, args.lighten) for s in stems]
    n = encode(frames, args.fps, args.hold, args.fade, args.out + ".mp4")
    print("%s.mp4: %d frames (%.1f s at %g fps) from %d states" % (args.out, n, n / args.fps, args.fps, len(stems)))

    poster = args.poster or stems[-1]
    save_png(composite(args.render, poster, args.lighten), args.out + "_poster.png")
    print("%s_poster.png: %s" % (args.out, poster))
    for stem in args.stills:
        save_png(composite(args.render, stem, args.lighten), "%s_%s.png" % (args.out, stem))
        print("%s_%s.png" % (args.out, stem))

    if args.layers:
        for layer in ("grains", "network"):
            frames = [lambda s=s: load(os.path.join(args.render, layer, s + ".bmp")) for s in stems]
            encode(frames, args.fps, args.hold, args.fade, "%s_%s.mp4" % (args.out, layer))
            print("%s_%s.mp4" % (args.out, layer))


if __name__ == "__main__":
    main()
