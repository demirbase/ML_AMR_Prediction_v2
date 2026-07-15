#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for 16_external_concordance.py parsers (M13), using the REAL output
formats captured on TRUBA from AMRFinderPlus 4.2.7 (DB 2026-05-15.1) and
ResFinder 4.5.0. Note the software/DB distinction: "2026-05-15.1" is the
DATABASE, which this file used to name as if it were the software version."""

from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# real AMRFinderPlus header (v2026-05-15.1)
AFP_HEADER = ("Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\t"
             "Scope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\t"
             "Reference sequence length\t% Coverage of reference\t% Identity to reference\t"
             "Alignment length\tClosest reference accession\tClosest reference name\t"
             "HMM accession\tHMM description")


def _afp_row(symbol, typ, cls, subcls):
    cols = ["NA", "contig", "1", "2", "+", symbol, "name", "core", typ, "AMR",
            cls, subcls, "EXACTX"] + ["NA"] * 9
    return "\t".join(cols)


def _write_afp(path, rows):
    path.write_text(AFP_HEADER + "\n" + "\n".join(rows) + "\n", encoding="utf-8")


@pytest.fixture
def mod(load_script):
    return load_script("16_external_concordance.py")


def test_amrfinder_tem_only(mod, tmp_path):
    # narrow-spectrum TEM-1 + a marR POINT mutation (MULTIDRUG w/ QUINOLONE)
    f = tmp_path / "afp_g1.tsv"
    _write_afp(f, [
        _afp_row("blaTEM-1", "AMR", "BETA-LACTAM", "BETA-LACTAM"),
        _afp_row("marR_S3N", "AMR", "MULTIDRUG",
                 "AMPICILLIN/CHLORAMPHENICOL/QUINOLONE/RIFAMPIN/TETRACYCLINE"),
        _afp_row("some_vir", "VIRULENCE", "NA", "NA"),   # must be ignored
    ])
    calls = mod.parse_amrfinder(f)
    assert calls["ampicillin"] == 1        # BETA-LACTAM (TEM) + AMPICILLIN (marR)
    assert calls["cefotaxime"] == 0        # no CEPHALOSPORIN -> narrow TEM ≠ cefotaxime
    assert calls["ciprofloxacin"] == 1     # QUINOLONE token (marR)


def test_amrfinder_esbl(mod, tmp_path):
    f = tmp_path / "afp_g2.tsv"
    _write_afp(f, [_afp_row("blaCTX-M-15", "AMR", "BETA-LACTAM", "CEPHALOSPORIN")])
    calls = mod.parse_amrfinder(f)
    assert calls["ampicillin"] == 1        # β-lactamase
    assert calls["cefotaxime"] == 1        # CEPHALOSPORIN ESBL
    assert calls["ciprofloxacin"] == 0


def test_resfinder_pheno_table(mod, tmp_path):
    f = tmp_path / "pheno_table_escherichia_coli.txt"
    f.write_text(
        "# ResFinder phenotype results for escherichia coli.\n"
        "# comment lines ignored\n"
        "# Antimicrobial\tClass\tWGS-predicted phenotype\tMatch\tGenetic background\n"
        "ampicillin\tbeta-lactam\tResistant\t3\tblaTEM-1A (blaTEM-1A_HM749966)\n"
        "cefotaxime\tbeta-lactam\tNo resistance\t0\n"
        "ciprofloxacin\tquinolone\tNo resistance\t0\n"
        "streptomycin\taminoglycoside\tResistant\t3\taadA1\n",
        encoding="utf-8")
    calls = mod.parse_resfinder(f)
    assert calls == {"ampicillin": 1, "cefotaxime": 0, "ciprofloxacin": 0}


def test_resfinder_species_file_preferred(mod, tmp_path):
    d = tmp_path / "rf_g1"
    d.mkdir()
    (d / "pheno_table.txt").write_text(
        "# Antimicrobial\tClass\tWGS-predicted phenotype\tMatch\n"
        "ampicillin\tbeta-lactam\tNo resistance\t0\n", encoding="utf-8")
    (d / "pheno_table_escherichia_coli.txt").write_text(
        "# Antimicrobial\tClass\tWGS-predicted phenotype\tMatch\n"
        "ampicillin\tbeta-lactam\tResistant\t3\tblaCTX\n", encoding="utf-8")
    chosen = mod._resfinder_pheno_file(d)
    assert chosen.name == "pheno_table_escherichia_coli.txt"
    assert mod.parse_resfinder(chosen)["ampicillin"] == 1


def test_head_to_head_shared_genomes(mod):
    # 4 model test genomes; tools + phenotype available for all. Model perfect,
    # AFP over-calls one S->R (a false resistant), RF perfect.
    genomes = ["g1", "g2", "g3", "g4", "gX"]   # gX has no model pred -> excluded
    pheno = {g: {"ampicillin": v} for g, v in
             zip(genomes, [1, 1, 0, 0, 1])}
    afp = {"g1": {"ampicillin": 1}, "g2": {"ampicillin": 1},
           "g3": {"ampicillin": 1}, "g4": {"ampicillin": 0}, "gX": {"ampicillin": 1}}
    rf = {"g1": {"ampicillin": 1}, "g2": {"ampicillin": 1},
          "g3": {"ampicillin": 0}, "g4": {"ampicillin": 0}, "gX": {"ampicillin": 1}}
    model_calls = {"ampicillin": {"g1": 1, "g2": 1, "g3": 0, "g4": 0}}  # no gX
    h = mod.head_to_head(genomes, pheno, afp, rf, model_calls, ["ampicillin"])
    amp = h["ampicillin"]
    assert amp["n_common_test_genomes"] == 4                  # gX dropped
    assert amp["model"]["balanced_accuracy"] == 1.0           # model perfect
    assert amp["resfinder"]["balanced_accuracy"] == 1.0
    assert amp["amrfinderplus"]["major_error_rate"] == pytest.approx(0.5)  # g3 S->R
    # model vs resfinder agree perfectly here
    assert amp["model_vs_resfinder"]["cohen_kappa"] == 1.0


def test_head_to_head_skips_antibiotic_without_model(mod):
    h = mod.head_to_head(["g1"], {"g1": {"cefotaxime": 1}},
                         {"g1": {"cefotaxime": 1}}, {"g1": {"cefotaxime": 1}},
                         {}, ["cefotaxime"])
    assert h == {}


def test_amrfinder_keywords_from_registry(mod):
    # Issue 9: keywords + default antibiotics are registry-driven, not hardcoded.
    from lib.registry import load_amrfinder_keywords
    reg = load_amrfinder_keywords()
    assert reg["cefotaxime"] == {"CEFOTAXIME", "CEPHALOSPORIN"}   # UPPER token set
    assert "BETA-LACTAM" in reg["ampicillin"]
    assert set(mod.AFP_KEYWORDS) == set(reg)                       # module uses the registry
    assert mod.DEFAULT_ANTIBIOTICS == list(mod.AFP_KEYWORDS)


def test_tokens_helper(mod):
    assert mod._tokens("BETA-LACTAM") == {"BETA-LACTAM"}
    assert mod._tokens("AMPICILLIN/QUINOLONE") == {"AMPICILLIN", "QUINOLONE"}
    assert mod._tokens("NA") == set()
    assert mod._tokens("") == set()


def test_write_kb_evidence(mod, tmp_path):
    import sqlite3
    from lib.kb_schema import create_schema
    from lib import concordance as C
    from lib.logging_utils import get_logger
    db = tmp_path / "amrk.db"
    conn = sqlite3.connect(str(db))
    create_schema(conn)
    conn.execute("INSERT INTO pipeline_runs(run_id, organism, antibiotic) VALUES ('R1','ecoli','ampicillin')")
    conn.execute("INSERT INTO antibiotics(antibiotic) VALUES ('ampicillin')")
    conn.execute("INSERT INTO models(model_id, run_id, antibiotic) VALUES (1,'R1','ampicillin')")
    conn.commit(); conn.close()

    yt = [1, 1, 0, 0]
    summary = {"antibiotics": {"ampicillin": {
                   "amrfinderplus": C.score_pair(yt, [1, 0, 0, 0]),
                   "resfinder": C.score_pair(yt, [1, 1, 0, 0])}},
               "head_to_head_model_test_genomes": {"ampicillin": {
                   "n_common_test_genomes": 4, "model": C.score_pair(yt, [1, 1, 0, 0])}}}
    mod.write_kb_evidence(db, summary, get_logger("test"))

    conn = sqlite3.connect(str(db))
    rows = conn.execute("SELECT evidence_type, pipeline_run_id FROM validation_evidence "
                        "ORDER BY evidence_type").fetchall()
    types = [r[0] for r in rows]
    assert "concordance_amrfinderplus" in types
    assert "concordance_resfinder" in types
    assert "head_to_head_model" in types
    assert all(r[1] == "R1" for r in rows)          # linked to the model's run
    # idempotent: second call does not duplicate
    mod.write_kb_evidence(db, summary, get_logger("test"))
    n = conn.execute("SELECT COUNT(*) FROM validation_evidence").fetchone()[0]
    assert n == 3
    conn.close()


# ---- M13 baseline provenance: software AND database versions ---------------
# The headline is "our model beats AMRFinderPlus" (K. pneu cipro: 0.926 vs
# 0.538). That means nothing unless the KB says which AMRFinderPlus — and the
# software version (4.2.7) and the DB version (2026-05-15.1) are different facts.

def _reload(monkeypatch, **env):
    import importlib.util
    for k in ("AMR_AFP_VERSION", "AMR_AFP_DB_VERSION", "AMR_RF_VERSION"):
        monkeypatch.delenv(k, raising=False)
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    spec = importlib.util.spec_from_file_location(
        "conc16", PROJECT_ROOT / "scripts" / "16_external_concordance.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_evidence_source_carries_software_and_db_version(monkeypatch):
    m = _reload(monkeypatch, AMR_AFP_VERSION="AMRFinderPlus 4.2.7",
                AMR_AFP_DB_VERSION="2026-05-15.1", AMR_RF_VERSION="ResFinder 4.5.0")
    assert m.AFP_SOURCE == "AMRFinderPlus 4.2.7 (DB 2026-05-15.1)"
    assert m.RF_SOURCE == "ResFinder 4.5.0"


def test_missing_tool_reports_unknown_not_an_error_message(monkeypatch):
    """`python -m resfinder --version` without resfinder exits 1 and prints
    "No module named resfinder" to stderr. Stamping that into the KB as a version
    is worse than admitting it is unknown."""
    m = _reload(monkeypatch)   # no env, and the tools are not installed locally
    assert "unknown" in m.RF_SOURCE.lower()
    assert "No module named" not in m.RF_SOURCE
    assert "unknown" in m.AFP_SOURCE.lower()
