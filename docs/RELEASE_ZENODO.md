# Release & Zenodo deposit checklist (M10 / FAIR)

How to cut a versioned, DOI-archived release of the **AMRK-DB** knowledge base.
Local prep (this repo already carries `.zenodo.json`, `CITATION.cff`, the README
*Data availability* section, and a `zenodo_doi`-preserving `populate_database.py`)
is done; the steps below are the **external** actions that need your Zenodo
account. Zenodo's *concept DOI* always resolves to the latest version, and each
version gets its own DOI — so re-depositing after adding an antibiotic is cheap.

## What goes in the deposit

The FAIR artifact is the KB + its evidence, **not** the multi-GB raw matrices
(those are regenerated from BV-BRC assemblies by the pipeline). Bundle:

- `results/ecoli/kb/amrk.db` — the SQLite knowledge base (schema 0.4.0)
- `results/ecoli/kb/15_cross_antibiotic_*` — cross-antibiotic overlap outputs
- `results/ecoli/<antibiotic>/` — per-antibiotic candidate + evidence CSVs/JSON
  (07/09/10/11/12/12b/13/14 outputs) and `runs/.../run_metadata.json`
- `config/config.yaml` + `config/experiments/ecoli/config_*.yaml` (the tuned splits)
- `CITATION.cff`, `README.md`, this file

```bash
# from repo root — build the bundle (adjust globs to the antibiotics present)
tar czf amrk-db_v0.4.0.tar.gz \
  results/ecoli/kb/amrk.db results/ecoli/kb/15_cross_antibiotic_* \
  results/ecoli/*/07_* results/ecoli/*/09_* results/ecoli/*/1[0-4]_* \
  runs/ecoli/*/run_metadata.json \
  config/config.yaml config/experiments/ecoli/config_*.yaml \
  CITATION.cff README.md docs/RELEASE_ZENODO.md
```

## Steps

1. **Decide the version.** `kb_schema_version` = schema shape (0.4.0). The *release*
   version tracks content too: schema-only-unchanged content growth (e.g. adding
   cefotaxime) → bump the **release** to `v0.5.0` even though the schema stays
   0.4.0. Set the same string in `.zenodo.json` (`version`) and `CITATION.cff`.
2. **Create the Zenodo deposit** (https://zenodo.org → New upload). Upload the
   tarball. Zenodo reads `.zenodo.json`-style metadata; verify title, creators
   (add affiliation/ORCID), license = CC-BY-4.0, keywords, related identifier
   (the GitHub repo).
3. **Reserve the DOI** (Zenodo "Reserve DOI" before publishing) so you can write
   it into the artifacts *before* the deposit is frozen.
4. **Stamp the DOI into the KB** — the reserved DOI survives future re-populates
   (populate now preserves it), so set it once:
   ```bash
   # via the env override on the next populate run …
   AMR_ZENODO_DOI="10.5281/zenodo.XXXXXXX" \
     python scripts/populate_database.py --antibiotic <ab>
   # … or directly, without re-populating:
   sqlite3 results/ecoli/kb/amrk.db \
     "UPDATE kb_metadata SET zenodo_doi='10.5281/zenodo.XXXXXXX' WHERE id=1;"
   ```
5. **Write the DOI into the docs:** add it to `CITATION.cff` (`identifiers:` block,
   already stubbed as a comment), `.zenodo.json` is fine as-is, and replace the
   README *"reserved — added on first deposit"* line with the concept DOI badge.
6. **Tag the release** (matches the version):
   ```bash
   git tag -a v0.5.0 -m "AMRK-DB v0.5.0 — <antibiotics>, schema 0.4.0, Zenodo 10.5281/zenodo.XXXXXXX"
   git push origin v0.5.0        # PAT, when ready
   ```
7. **Publish** the Zenodo deposit → the DOI becomes permanent. Put the concept
   DOI in the thesis Methods (data availability) paragraph.

## Methods paragraph (thesis — data availability)

> The AMRK-DB knowledge base is openly available under CC-BY-4.0, archived on
> Zenodo (DOI: 10.5281/zenodo.XXXXXXX; concept DOI resolving to the latest
> version) and versioned by `kb_schema_version` (0.4.0). Each biomarker links to
> a provenance record (git commit, random seed, configuration hash, CARD v4.0.1)
> enabling exact reproduction. Genome assemblies were obtained from BV-BRC and
> the analysis pipeline (MIT-licensed, github.com/demirbase/ML_AMR_Prediction_v2)
> regenerates all features; confidence tiers denote statistical evidence, not
> asserted biological causation.

## Note on the GitHub↔Zenodo webhook

Zenodo can auto-archive GitHub *releases*, but that archives the code repo, not
the KB data bundle, and the repo is currently a fork (see the repo settings). A
**manual dataset upload** of the bundle above is simpler and repo-independent —
prefer it for the KB. Use the GitHub webhook only if you also want the code
snapshot archived.
