#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Guard: nothing may be tracked under the HPC's symlinked directories.

On TRUBA, $AMR_HOME/{data,results,runs,models,logs} are symlinks to /arf/scratch.
Git will not write a tracked path through a symlink — it deletes the symlink and
puts a real directory there instead, detaching the repo from every genome and
result on scratch. That is not hypothetical: on 2026-07-15 a `git reset --hard`
severed data/ and runs/ because card_nt/*.n* and runs/.gitkeep were committed.

This test fails the moment someone re-adds a file under those prefixes, which is
cheaper than rediscovering it mid-deploy.
"""

import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Symlinked to scratch on the HPC — see .gitignore
SYMLINKED_PREFIXES = ("data/", "results/", "runs/", "models/", "logs/")


def test_nothing_tracked_under_symlinked_dirs():
    out = subprocess.run(["git", "ls-files"], cwd=PROJECT_ROOT,
                         capture_output=True, text=True, check=True).stdout
    offenders = [f for f in out.splitlines()
                 if f.startswith(SYMLINKED_PREFIXES)]
    assert not offenders, (
        "These files are tracked under a directory that is a symlink to scratch on "
        "the HPC. `git reset --hard` will replace the symlink with a real directory "
        "and detach the repo from the data:\n  " + "\n  ".join(offenders)
    )
