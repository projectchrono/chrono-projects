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
# Launch one run of a DEM validation program in its own run folder:
# - copies the executable into the run folder (later rebuilds do not affect the run)
# - records the git state of the Chrono and chrono-projects source trees (commit, and a fingerprint of any
#   uncommitted changes) and passes it to the program, which writes it to meta.json
# - runs the program with the run folder as working directory, at below-normal priority, with its console output
#   in stdout.txt
#
# Usage:
#   python run_case.py --exe <program.exe> --run <run folder> --chrono <Chrono source> --projects <chrono-projects
#                      source> [--path <DLL folder>] -- <program options>
#
# =============================================================================

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time


def git_state(repo):
    def git(*args):
        return subprocess.run(["git", "-C", repo] + list(args), capture_output=True, text=True).stdout.strip()

    head = git("rev-parse", "HEAD")
    status = git("status", "--porcelain", "--untracked-files=no")
    untracked = [f for f in git("ls-files", "--others", "--exclude-standard").splitlines() if not f.endswith(".pyc")]

    # Fingerprint of the uncommitted changes: the diff of the tracked files and the content of the untracked ones
    h = hashlib.sha1()
    h.update(subprocess.run(["git", "-C", repo, "diff", "HEAD"], capture_output=True).stdout)
    for f in sorted(untracked):
        h.update(f.encode())
        with open(os.path.join(repo, f), "rb") as fin:
            h.update(fin.read())
    changed = (len(status.splitlines()) if status else 0) + len(untracked)
    return {
        "commit": head,
        "uncommitted_files": changed,
        "untracked_files": untracked,
        "changes_sha1": h.hexdigest() if changed else "",
    }


def main():
    argv = sys.argv[1:]
    prog_args = []
    if "--" in argv:
        k = argv.index("--")
        prog_args = argv[k + 1:]
        argv = argv[:k]
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--chrono", required=True)
    ap.add_argument("--projects", required=True)
    ap.add_argument("--path", default="")
    args = ap.parse_args(argv)

    os.makedirs(args.run, exist_ok=True)
    exe = os.path.join(args.run, os.path.basename(args.exe))
    shutil.copy2(args.exe, exe)

    prov = {
        "chrono": git_state(args.chrono),
        "chrono_projects": git_state(args.projects),
        "executable": os.path.abspath(args.exe),
        "command": " ".join([os.path.basename(args.exe)] + prog_args),
        "launched": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    cmd = [exe] + prog_args + ["--provenance", json.dumps(prov)]

    env = dict(os.environ)
    if args.path:
        env["PATH"] = args.path + os.pathsep + env.get("PATH", "")
    flags = 0x00004000 if os.name == "nt" else 0  # BELOW_NORMAL_PRIORITY_CLASS
    with open(os.path.join(args.run, "stdout.txt"), "w") as out:
        p = subprocess.Popen(cmd, cwd=args.run, stdout=out, stderr=subprocess.STDOUT, env=env, creationflags=flags)
        rc = p.wait()
    with open(os.path.join(args.run, "exit.txt"), "w") as f:
        f.write("exit code %d at %s\n" % (rc, time.strftime("%Y-%m-%d %H:%M:%S")))
    sys.exit(rc)


if __name__ == "__main__":
    main()
