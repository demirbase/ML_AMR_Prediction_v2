#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 00 — Metadata preparation (raw download -> binary phenotype matrix).

PURPOSE (to be implemented):
    Convert the raw AMR metadata as downloaded from the source database
    (e.g. BV-BRC: data/external/{organism}/metadata/BVBRC_genome_amr.csv) into
    the binary phenotype matrix the rest of the pipeline expects:

        data/external/{organism}/metadata/amr_phenotypes.csv
        columns: "Genome ID" + one column per antibiotic, values in {0, 1}
                 (1 = Resistant, 0 = Susceptible; NaN/blank = untested)

    Steps 01/03 consume this file directly; the per-genome labels (y_*.csv) are
    materialised later in 03_matrix_construction.py — NOT here and NOT in 01.

PLANNED BEHAVIOUR (content to be written later):
    - read the raw long/wide source file (genome id, antibiotic, phenotype text)
    - normalise antibiotic names (lowercase; "/" -> "_") against the registry
      (config/registry/antibiotics.yaml)
    - map phenotype text (Resistant/Susceptible/Intermediate/…) -> {1, 0, NaN}
      with an explicit, documented mapping (Intermediate handling configurable)
    - pivot to wide: one row per Genome ID, one column per antibiotic
    - keep only antibiotics listed for this organism in organisms.yaml
    - write amr_phenotypes.csv (organism-scoped, via resolve_path)

NOTE: This is a SKELETON. The transformation logic will be added once the exact
raw source schema is confirmed.
"""

import sys


def main() -> None:
    """TODO: implement raw-metadata -> amr_phenotypes.csv conversion."""
    print("00_prepare_metadata.py is a placeholder — implementation pending.")
    print("It will convert the raw AMR download into the binary amr_phenotypes.csv.")
    sys.exit(0)


if __name__ == "__main__":
    main()
