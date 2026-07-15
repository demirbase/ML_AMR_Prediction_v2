#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the shared library (SCALE_MLOPS_PLAN.md §7.5).

Run with:
    pytest tests/
or directly:
    python tests/test_lib.py
"""

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import registry, run_metadata           # noqa: E402
from lib.chunking import get_y_chunk              # noqa: E402
from lib.config import load_config, resolve_path, get_target  # noqa: E402


# ---------------------------------------------------------------------------
# chunking
# ---------------------------------------------------------------------------
def test_get_y_chunk_basic():
    data = list(range(10))
    assert get_y_chunk(data, 0, 3, 10) == [0, 1, 2]
    assert get_y_chunk(data, 1, 3, 10) == [3, 4, 5]


def test_get_y_chunk_last_partial():
    data = list(range(10))
    # chunk 3 with size 3 -> indices 9..10 clamped to total_len
    assert get_y_chunk(data, 3, 3, 10) == [9]


def test_get_y_chunk_out_of_range():
    data = list(range(10))
    assert get_y_chunk(data, 5, 3, 10) == []


# ---------------------------------------------------------------------------
# registry
# ---------------------------------------------------------------------------
def test_registry_classes_structure():
    classes = registry.load_antibiotic_classes()
    assert "Aminoglycosides" in classes
    assert "gentamicin" in classes["Aminoglycosides"]
    # legacy structure: {DisplayName: [members]}
    assert all(isinstance(v, list) for v in classes.values())


def test_registry_reverse_index():
    assert registry.antibiotic_to_class("gentamicin") == "aminoglycosides"
    assert registry.antibiotic_to_class("ciprofloxacin") == "quinolones"
    # carbapenems are their own class (schema 2.0), not lumped with beta-lactams
    assert registry.antibiotic_to_class("meropenem") == "carbapenems"
    assert registry.antibiotic_to_class("not_a_real_drug") is None


def test_registry_metadata_accessors():
    # mechanism_type resolves via class; who_aware is per-antibiotic (registry = source)
    assert registry.antibiotic_mechanism_type("ciprofloxacin") == "target_snp"
    assert registry.antibiotic_mechanism_type("meropenem") == "acquired"
    assert registry.antibiotic_mechanism_type("not_a_real_drug") is None
    assert registry.antibiotic_who_aware("colistin") == "Reserve"
    assert registry.antibiotic_who_aware("ampicillin") == "Access"
    assert registry.antibiotic_who_aware("not_a_real_drug") is None


def test_registry_targets_and_validation():
    targets = registry.list_targets(enabled_only=True)
    # ecoli + kpneumoniae are status: done -> active targets (schema 2.0)
    assert ("ecoli", "gentamicin") in targets
    assert ("kpneumoniae", "meropenem") in targets
    # eskapee_phase filter: pseudomonas is phase 2, excluded from phase-1 list
    phase1 = registry.list_targets(phase=1)
    assert ("ecoli", "gentamicin") in phase1
    assert all(org != "pseudomonas_aeruginosa" for org, _ab in phase1)
    assert registry.validate_target("ecoli", "gentamicin") is True
    assert registry.validate_target("ecoli", "meropenem") is False


# ---------------------------------------------------------------------------
# config / path resolution
# ---------------------------------------------------------------------------
def test_resolve_path_organism_antibiotic():
    # Force the k-mer (base) layout so this templating check is independent of
    # config.yaml's feature_repr default (now 'unitig'); the redirect itself is
    # covered by test_resolve_path_feature_repr_switch below.
    base = load_config()
    cfg = {**base, "preprocessing": {**base.get("preprocessing", {}), "feature_repr": "kmer"}}
    p = resolve_path("matrix_dir", organism="ecoli", antibiotic="gentamicin", config=cfg)
    assert p.as_posix().endswith("data/processed/ecoli/gentamicin/matrix")


def test_resolve_path_feature_repr_switch():
    # The unitig pivot switch (ROADMAP §0 M12): feature_repr redirects ONLY the
    # matrix_dir key, leaving every other path untouched. Use a synthetic config
    # so the test is independent of the repo config.yaml's current value.
    base = load_config()
    cfg_kmer = {**base, "preprocessing": {**base.get("preprocessing", {}), "feature_repr": "kmer"}}
    cfg_unitig = {**base, "preprocessing": {**base.get("preprocessing", {}), "feature_repr": "unitig"},
                  "unitig": {"out_subdir": "matrix_unitig"}}

    p_kmer = resolve_path("matrix_dir", organism="ecoli", antibiotic="ampicillin", config=cfg_kmer)
    p_unitig = resolve_path("matrix_dir", organism="ecoli", antibiotic="ampicillin", config=cfg_unitig)
    assert p_kmer.name == "matrix"
    assert p_unitig.name == "matrix_unitig"
    assert p_unitig.parent == p_kmer.parent          # same {antibiotic} dir, only leaf differs
    # Non-matrix keys must be unaffected by the switch.
    assert resolve_path("models_dir", organism="ecoli", antibiotic="ampicillin",
                        config=cfg_unitig).name == "ampicillin"

    # AMR_FEATURE_REPR env overrides config (HPC convenience).
    import os
    prev = os.environ.get("AMR_FEATURE_REPR")
    try:
        os.environ["AMR_FEATURE_REPR"] = "unitig"
        assert resolve_path("matrix_dir", organism="ecoli", antibiotic="ampicillin",
                            config=cfg_kmer).name == "matrix_unitig"
    finally:
        if prev is None:
            os.environ.pop("AMR_FEATURE_REPR", None)
        else:
            os.environ["AMR_FEATURE_REPR"] = prev


def test_resolve_path_run_id():
    p = resolve_path("run_dir", organism="ecoli", antibiotic="gentamicin", run_id="RID123")
    assert p.name == "RID123"


def test_resolve_path_global_key_no_placeholder():
    # A global key with no placeholder (kmc_bin) resolves directly, with or
    # without organism/antibiotic supplied.
    p = resolve_path("kmc_bin")
    assert p.name == "kmc"


def test_resolve_path_unknown_key():
    raised = False
    try:
        resolve_path("definitely_not_a_real_key")
    except KeyError:
        raised = True
    assert raised, "resolve_path should raise KeyError for an unknown key"


def test_get_target_defaults_from_config():
    org, ab = get_target()
    assert org == "ecoli"
    assert ab  # a non-empty antibiotic from config


# ---------------------------------------------------------------------------
# run metadata
# ---------------------------------------------------------------------------
def test_make_run_id_format():
    rid = run_metadata.make_run_id("ecoli", "gentamicin")
    parts = rid.split("__")
    assert len(parts) == 4
    assert parts[0] == "ecoli" and parts[1] == "gentamicin"


def test_hash_files_stable(tmp_path):
    f = tmp_path / "a.txt"
    f.write_text("hello")
    h1 = run_metadata.hash_files([f])
    h2 = run_metadata.hash_files([f])
    assert h1 == h2 and len(h1) == 64


if __name__ == "__main__":
    # Minimal runner so the file works without pytest installed.
    import tempfile, traceback
    funcs = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = failed = 0
    for fn in funcs:
        try:
            if "tmp_path" in fn.__code__.co_varnames:
                with tempfile.TemporaryDirectory() as d:
                    fn(Path(d))
            else:
                fn()
            print(f"  PASS  {fn.__name__}")
            passed += 1
        except Exception:
            print(f"  FAIL  {fn.__name__}")
            traceback.print_exc()
            failed += 1
    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


# ---- compute-resource env overrides ---------------------------------------
# kmc_mem/threads describe machine size, not science. They used to be a hand-edit
# on the HPC that `git reset --hard` wiped on every deploy — a forgotten re-edit
# meant jobs ran with laptop resources, slowly or OOM, with no clue in the logs.

def test_resource_keys_default_to_laptop_safe_values(monkeypatch):
    from lib.config import load_config
    monkeypatch.delenv("AMR_KMC_MEM", raising=False)
    monkeypatch.delenv("AMR_THREADS", raising=False)
    pre = load_config()["preprocessing"]
    assert pre["kmc_mem"] == 16 and pre["threads"] == 10


def test_env_overrides_resource_keys(monkeypatch):
    from lib.config import load_config
    monkeypatch.setenv("AMR_KMC_MEM", "128")
    monkeypatch.setenv("AMR_THREADS", "20")
    pre = load_config()["preprocessing"]
    assert pre["kmc_mem"] == 128 and pre["threads"] == 20


def test_bad_resource_env_raises_instead_of_silently_defaulting(monkeypatch):
    from lib.config import load_config
    monkeypatch.setenv("AMR_THREADS", "twenty")
    with pytest.raises(ValueError, match="AMR_THREADS"):
        load_config()


def test_env_int_helper(monkeypatch):
    from lib.config import env_int
    monkeypatch.delenv("AMR_TEST_INT", raising=False)
    assert env_int("AMR_TEST_INT", 7) == 7
    monkeypatch.setenv("AMR_TEST_INT", "42")
    assert env_int("AMR_TEST_INT", 7) == 42
    monkeypatch.setenv("AMR_TEST_INT", "  ")      # blank -> default, not a crash
    assert env_int("AMR_TEST_INT", 7) == 7


# ---- organism status vocabulary -------------------------------------------
# is_active() just tests set membership, so an unknown status silently drops the
# organism from the panel with nothing raising. validate_registry pins it shut.

def test_status_vocabulary_is_closed():
    from lib import registry
    assert registry.VALID_STATUS == {
        "done", "in_progress", "planned", "excluded_insufficient_data"}


def test_excluded_organism_is_not_an_active_target():
    from lib import registry
    ent = registry.get_organism("enterobacter_cloacae")
    assert ent["status"] == "excluded_insufficient_data"
    assert not registry.is_active(ent)          # recorded negative finding, not pending work
    active = {o for o, _ in registry.list_targets(enabled_only=True)}
    assert "enterobacter_cloacae" not in active


def test_validate_registry_rejects_an_unknown_status(monkeypatch):
    """A typo'd status must fail loudly, not quietly deactivate the organism."""
    import importlib.util
    from lib import registry
    spec = importlib.util.spec_from_file_location(
        "vr", PROJECT_ROOT / "scripts" / "validate_registry.py")
    vr = importlib.util.module_from_spec(spec); spec.loader.exec_module(vr)

    orgs = {k: dict(v) for k, v in registry.load_organisms().items()}
    orgs["ecoli"]["status"] = "in_progres"      # plausible typo
    monkeypatch.setattr(registry, "load_organisms", lambda: orgs)
    monkeypatch.setattr(vr.registry, "load_organisms", lambda: orgs)

    errors, warnings = [], []
    vr._check_registry(errors, warnings)
    assert any("in_progres" in e for e in errors), errors
