#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run metadata & reproducibility capture (SCALE_MLOPS_PLAN.md §7.1).

Every pipeline run can stamp a reproducible fingerprint — git commit, seeds,
tool versions, data hashes — so a model file, log, result and (future) KB row
all trace back to a single run_id. All helpers are best-effort and must never
crash the pipeline: callers wrap them in try/except and treat failure as
non-fatal.

run_id format (SCALE_MLOPS_PLAN.md §2):
    {organism}__{antibiotic}__{UTC:%Y%m%dT%H%M}__{git7}
"""

import datetime
import hashlib
import json
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _git(*args):
    """Run a git command in the project root; return stripped stdout or None."""
    try:
        out = subprocess.run(
            ["git", *args],
            cwd=str(PROJECT_ROOT),
            capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return None


def git_commit_hash(short=False):
    return _git("rev-parse", "--short", "HEAD") if short else _git("rev-parse", "HEAD")


def git_is_dirty():
    """True if the working tree has uncommitted changes (None if git unavailable)."""
    status = _git("status", "--porcelain")
    if status is None:
        return None
    return bool(status.strip())


def _tool_version(cmd, flag="--version"):
    try:
        out = subprocess.run([cmd, flag], capture_output=True, text=True, check=False)
        text = (out.stdout or out.stderr or "").strip().splitlines()
        line = text[0].strip() if text else ""
        # A wrong flag makes some tools echo an error as their first line (bcalm with
        # --version -> "ERROR: Unknown parameter '--version'"). Storing that as the
        # version is worse than a NULL, so reject error-shaped output and record an
        # honest "unknown" instead.
        low = line.lower()
        if (not line
                or low.startswith(("error", "unknown", "unrecognized", "invalid", "usage"))
                or "unknown parameter" in low):
            return None
        return line
    except Exception:
        return None


def _pkg_version(name):
    try:
        mod = __import__(name)
        return getattr(mod, "__version__", None)
    except Exception:
        return None


def collect_versions(config=None):
    """Best-effort capture of language / library / external-tool versions.

    The tools that DEFINE the science are recorded, not just the incidental ones.
    This used to capture kmc (a QC-only tool for the abandoned k-mer baseline)
    while omitting unitig-caller — which builds the features — and PopPUNK —
    which defines the cross-validation groups. A KB that cannot say which PopPUNK
    produced its lineage labels has a provenance chain broken exactly where it
    matters.

    graph-tool is here for a hard-won reason: on 2026-07-15 a container rebuild
    held poppunk pinned at 2.7.8 while graph-tool — PopPUNK's network backend —
    re-solved 2.98 -> 3.0 underneath it, and the E. coli clustering moved (324 ->
    397 lineages, ARI 0.990 vs the old labels). PopPUNK's own version does NOT
    pin its behaviour; graph-tool's is part of the answer.

    Best-effort by design: a tool living in another container (pyseer is in
    amr-tools.sif, not amr.sif where populate runs) simply returns None here and
    is recorded by the step that actually runs it — see 14_pyseer_lmm, which
    writes its pyseer version into its own summary.
    """
    versions = {
        "python": sys.version.split()[0],
        "xgboost": _pkg_version("xgboost"),
        "scikit_learn": _pkg_version("sklearn"),
        "numpy": _pkg_version("numpy"),
        "scipy": _pkg_version("scipy"),
        "kmc": _tool_version("kmc", "-h") or _tool_version("kmc"),
        "blastn": _tool_version("blastn", "-version"),
        # ── the tools the results actually depend on ──────────────────────────
        "unitig_caller": _tool_version("unitig-caller"),   # builds the features
        # bcalm (unitig-caller's internal cDBG builder) has NO usable version CLI:
        # --version errors, -version prints nothing, and the banner omits it. So this
        # is an honest None; bcalm's version is implied by unitig-caller (captured
        # above) and pinned exactly in environment.lock.yml. Kept as a probe in case
        # a future bcalm build adds one — _tool_version rejects the error line if not.
        "bcalm": _tool_version("bcalm", "-version"),
        "poppunk": _tool_version("poppunk"),               # defines the CV groups
        "graph_tool": _pkg_version("graph_tool"),          # changes PopPUNK's clustering
        "pyseer": _tool_version("pyseer"),                 # None from amr.sif; see docstring
    }
    # The KB schema version comes from CODE, never from config. As a config key it
    # rotted to 0.6.1 while the schema was 0.7.1, and this function copied that
    # stale value into run_metadata.json — i.e. the KB misreported its own schema.
    try:
        from lib.kb_schema import KB_SCHEMA_VERSION
        versions["kb_schema_version"] = KB_SCHEMA_VERSION
    except Exception:
        versions["kb_schema_version"] = None
    if config is not None:
        prov = config.get("provenance", {}) or {}
        versions["card_version"] = prov.get("card_version")
    return versions


def make_run_id(organism, antibiotic):
    """Build a run_id: {org}__{ab}__{UTC}__{git7}."""
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M")
    git7 = git_commit_hash(short=True) or "nogit"
    return f"{organism}__{antibiotic}__{ts}__{git7}"


def hash_files(paths):
    """
    Return a stable fingerprint of a set of files (sha256 of names + sizes +
    per-file sha256). Used to detect when the training data changed between runs.
    """
    h = hashlib.sha256()
    for p in sorted(Path(x) for x in paths):
        try:
            h.update(p.name.encode())
            data = p.read_bytes()
            h.update(str(len(data)).encode())
            h.update(hashlib.sha256(data).digest())
        except Exception:
            h.update(b"<missing>")
    return h.hexdigest()


def build_run_metadata(organism, antibiotic, run_id, seed=None,
                       params=None, data_files=None, config=None, extra=None):
    """Assemble the run_metadata.json payload (SCALE_MLOPS_PLAN §7.1)."""
    meta = {
        "run_id": run_id,
        "organism": organism,
        "antibiotic": antibiotic,
        "git_commit_hash": git_commit_hash(),
        "git_dirty": git_is_dirty(),
        "random_seed": seed,
        "versions": collect_versions(config),
        "params": params or {},
        "started_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    if data_files:
        meta["data_fingerprint"] = {
            "n_files": len(list(data_files)),
            "sha256": hash_files(data_files),
        }
    if extra:
        meta.update(extra)
    return meta


def write_json(path, payload):
    """Write a JSON payload, creating parent dirs. Best-effort (returns bool)."""
    try:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, default=str)
        return True
    except Exception as e:
        print(f"  ⚠ Could not write {path}: {e}")
        return False
