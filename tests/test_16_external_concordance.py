#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for 16_external_concordance.py parsers (M13), using the REAL
AMRFinderPlus 2026-05-15.1 and ResFinder 4.5.0 output formats captured on TRUBA."""

import pytest

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


def test_tokens_helper(mod):
    assert mod._tokens("BETA-LACTAM") == {"BETA-LACTAM"}
    assert mod._tokens("AMPICILLIN/QUINOLONE") == {"AMPICILLIN", "QUINOLONE"}
    assert mod._tokens("NA") == set()
    assert mod._tokens("") == set()
