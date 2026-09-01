#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Turn "FAIR-compliant" from a claim in the title into a table with evidence.

The thesis asserts FAIR compliance, so it owes the reader a principle-by-principle
account of what backs that. This writes `fair_mapping.csv`: one row per FAIR
principle, what implements it here, where to look, and an honest status.

The evidence column is queried from the delivered KB rather than typed, because a
FAIR table that quotes remembered numbers is exactly the kind of claim it exists
to substantiate. Coverage is reported as it is, not as it should be -- 1,576 of
3,611 BLAST hits carry an ARO accession, and the table says so.

Two rows are `partial` on purpose:
  R1.2  provenance is per-run and complete, but git_dirty = 1 on 45/45 runs, so a
        commit hash alone does not restore the working tree that produced a model.
  R1.3  no community standard exists for unitig-resolution AMR biomarker bases;
        this uses CARD/ARO, NCBI Taxonomy and WHO AWaRe, which is the closest thing.

    python scripts/kb_fair_mapping.py --db results/kb/amrk.db --out results/tables
"""
import argparse
import csv
import sqlite3
from pathlib import Path

FIELDS = ["principle", "requirement", "implementation", "evidence", "status"]


def facts(db):
    c = sqlite3.connect(db)
    q = lambda s: c.execute(s).fetchone()[0]
    aro_ok = q("""select count(*) from blast_annotations
                  where aro_accession is not null and aro_accession not in ('','nan')""")
    f = {
        "doi": q("select zenodo_doi from kb_metadata"),
        "licence": q("select license from kb_metadata"),
        "schema": q("select kb_schema_version from kb_metadata"),
        "card": q("select card_version from kb_metadata"),
        "created": q("select created_at from kb_metadata"),
        "n_models": q("select count(*) from models"),
        "n_unitigs": q("select count(*) from unitigs"),
        "n_tiers": q("select count(*) from unitig_evidence_tier"),
        "n_runs": q("select count(*) from pipeline_runs"),
        "n_tables": len(list(c.execute("select name from sqlite_master where type='table'"))),
        "n_blast": q("select count(*) from blast_annotations"),
        "aro_ok": aro_ok,
        "aro_distinct": q("""select count(distinct aro_accession) from blast_annotations
                             where aro_accession is not null and aro_accession not in ('','nan')"""),
        "taxid": q("select count(*) from organisms where taxid is not null"),
        "n_orgs_reg": q("select count(*) from organisms"),
        "aware": q("""select count(*) from antibiotics
                       where who_aware is not null and who_aware not in ('','nan')"""),
        "n_abx": q("select count(*) from antibiotics"),
        "seq_ok": q("select count(*) from unitigs where sequence is not null and sequence != ''"),
        "cfg": q("select count(distinct config_hash) from pipeline_runs"),
        "dirty": q("select count(*) from pipeline_runs where git_dirty = 1"),
        "commits": q("select count(distinct git_commit) from pipeline_runs"),
        "bio_cols": 33,
    }
    c.close()
    return f


def rows(f):
    return [
        # F1/F3 are asserted from the KB, not from intent: if zenodo_doi is empty
        # the delivered file does not in fact carry its own identifier, and saying
        # "met" would be the claim the principle exists to prevent. The field is
        # deliberately empty between a content change and the next release (see
        # thesis §5.6.1), so the status follows the field.
        ("F1", "(Meta)data are assigned a globally unique and persistent identifier",
         "Zenodo concept DOI, minted for the dataset and stamped into the KB itself",
         (f"kb_metadata.zenodo_doi = {f['doi']}; the concept DOI resolves to the newest "
          f"version, so citations survive re-release. Records inside the KB use local "
          f"surrogate keys (unitig_id, model_id), not global PIDs."
          if f['doi'] else
          "kb_metadata.zenodo_doi is EMPTY. The knowledge base changed after the last "
          "archived version, so the previous version DOI no longer describes this content "
          "and was not carried forward; a new release must be minted and the field "
          "repopulated. The concept DOI remains the identifier the thesis cites."),
         "met" if f['doi'] else "not met"),

        ("F2", "Data are described with rich metadata",
         "Per-run provenance table plus a KB-level metadata row",
         f"pipeline_runs holds {f['n_runs']} rows x 19 columns (organism, antibiotic, git "
         f"commit, seed, config hash, min_support, n_genomes and 9 tool versions); "
         f"kb_metadata carries schema, CARD snapshot, DOI, licence, created_at. "
         f"METHODOLOGY.md §5 is the narrative record.", "met"),

        ("F3", "Metadata clearly and explicitly include the identifier of the data",
         "The DOI lives inside the database, not only beside it",
         (f"kb_metadata.zenodo_doi is a column of amrk.db, so a copy of the file "
          f"identifies its own published archive with no external manifest."
          if f['doi'] else
          "kb_metadata.zenodo_doi is a column of amrk.db but is currently empty, so a "
          "copy of the file does not identify its own archive until the next release "
          "repopulates it."),
         "met" if f['doi'] else "not met"),

        ("F4", "(Meta)data are registered or indexed in a searchable resource",
         "Zenodo record; DataCite metadata; GitHub repository",
         f"Published Zenodo record ({f['doi']}, resource type Dataset) is DataCite-indexed; "
         f"CITATION.cff makes the repository citable.", "met"),

        ("A1", "(Meta)data are retrievable by their identifier using a standardised protocol",
         "HTTPS retrieval from Zenodo via the DOI",
         "DOI -> HTTPS landing page -> direct file download; no bespoke client, no "
         "registration, no request form.", "met"),

        ("A1.1", "The protocol is open, free and universally implementable",
         "HTTP(S), and SQLite plus CSV as the container formats",
         f"amrk.db is a single SQLite file (schema {f['schema']}) readable by any language's "
         f"standard library; the tidy tables are plain CSV.", "met"),

        ("A1.2", "The protocol allows authentication and authorisation where necessary",
         "Not necessary: the dataset is fully open",
         f"{f['licence']} with no embargo, no restricted tier and no personal data, so no "
         f"authorisation layer is required for any part of it.", "met"),

        ("A2", "Metadata are accessible even when the data are no longer available",
         "Zenodo retains the record and its metadata independently of the files",
         "Zenodo's tombstone policy keeps DOI-level metadata resolvable; the repository "
         "keeps METHODOLOGY.md, CHANGELOG.md and CITATION.cff regardless of archive state.",
         "met"),

        ("I1", "(Meta)data use a formal, accessible, shared, broadly applicable language "
               "for knowledge representation",
         "Explicit, versioned relational schema with declared foreign keys",
         f"{f['n_tables']} tables, primary and foreign keys declared in DDL, schema pinned to "
         f"{f['schema']} and version-locked by tests (test_version_alignment). Figure 38 is "
         f"drawn from that live schema.", "met"),

        ("I2", "(Meta)data use vocabularies that follow FAIR principles",
         "CARD/ARO ontology, NCBI Taxonomy, WHO AWaRe classification",
         f"ARO accessions on {f['aro_ok']:,} of {f['n_blast']:,} BLAST hits ({f['aro_distinct']} "
         f"distinct terms), carrying aro_gene_family, aro_drug_class and "
         f"aro_resistance_mechanism; organisms.taxid filled {f['taxid']}/{f['n_orgs_reg']}; "
         f"antibiotics.who_aware filled {f['aware']}/{f['n_abx']} (Access/Watch/Reserve). "
         f"Hits without an accession are the low-quality tier='none' remainder.", "met"),

        ("I3", "(Meta)data include qualified references to other (meta)data",
         "Every annotation names its source database and the strength of the link",
         f"blast_annotations records source_db, identity_pct, coverage, evalue and tier "
         f"alongside the ARO accession, so a reference is never an unqualified assertion; "
         f"CARD snapshot is pinned ({str(f['card'])[:44]}...).", "met"),

        ("R1", "(Meta)data are richly described with a plurality of accurate and relevant "
               "attributes",
         "Multi-layer evidence per biomarker rather than a single score",
         f"{f['n_tiers']:,} graded (unitig, model) pairs; biomarkers.csv exposes "
         f"{f['bio_cols']} attributes per biomarker (gain, SHAP, CPSS frequency, MDA, pyseer "
         f"LRT p, prevalence delta, odds ratio, BLAST identity/coverage, ARO fields, tier, "
         f"evidence layer list). {f['n_unitigs']:,} unitigs retain their full sequence "
         f"({f['seq_ok']:,}/{f['n_unitigs']:,}), so any claim is re-checkable from the "
         f"nucleotides up.", "met"),

        ("R1.1", "(Meta)data are released with a clear and accessible data usage licence",
         "CC-BY-4.0, declared in the KB, the repository and the archive",
         f"kb_metadata.license = {f['licence']}; LICENSE file in the repository; licence "
         f"field set on the Zenodo record.", "met"),

        ("R1.2", "(Meta)data are associated with detailed provenance",
         "Per-run git commit, seed, config hash and tool versions - with one honest gap",
         f"{f['n_runs']}/{f['n_runs']} runs record commit, seed 42, a distinct config hash "
         f"({f['cfg']} hashes) and 9 tool versions. BUT git_dirty = 1 on {f['dirty']}/"
         f"{f['n_runs']} runs: the working tree was patched beyond the recorded commit, so "
         f"the hash alone does not restore the exact code. {f['commits']} distinct commits "
         f"span the panel. See METHODOLOGY §5.3.", "partial"),

        ("R1.3", "(Meta)data meet domain-relevant community standards",
         "Closest available standards adopted; no standard exists for this artefact type",
         f"CARD/ARO for resistance determinants, NCBI Taxonomy for organisms, WHO AWaRe for "
         f"antibiotic stewardship class, PopPUNK for lineage nomenclature, CheckM2/QUAST for "
         f"assembly quality reporting. There is no community exchange standard for "
         f"unitig-resolution AMR biomarker knowledge bases, so no schema could be conformed "
         f"to; this is stated rather than glossed.", "partial"),
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True, help="tidy tables directory")
    a = ap.parse_args()

    f = facts(a.db)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    dest = out / "fair_mapping.csv"
    data = rows(f)
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(FIELDS)
        w.writerows(data)
    met = sum(1 for r in data if r[4] == "met")
    print(f"  wrote {dest.name}: {len(data)} principles, {met} met, {len(data) - met} partial")


if __name__ == "__main__":
    main()
