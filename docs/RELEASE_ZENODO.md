# Release & Zenodo deposit checklist (M10 / FAIR)

How to cut a versioned, DOI-archived release of the **AMRK-DB** knowledge base.
Local prep (this repo already carries `.zenodo.json`, `CITATION.cff`, the README
*Data availability* section, and a `zenodo_doi`-preserving `populate_database.py`)
is done; the steps below are the **external** actions that need your Zenodo
account. Zenodo's *concept DOI* always resolves to the latest version, and each
version gets its own DOI — so re-depositing after adding an antibiotic is cheap.

## What goes in the deposit

The FAIR artifact is the KB + its evidence + the figures/tables that summarise it —
**not** the multi-GB unitig matrices or the raw assemblies (both are regenerated from
BV-BRC by the pipeline; `unitigs.rtab` alone is ~36 GB for E. coli).

Contents (45 models · 6 ESKAPEE organisms · 14 drug classes · schema 0.7.1):

- `results/kb/amrk.db` — the unified SQLite knowledge base (all organisms in one file;
  `pipeline_runs.organism` separates them)
- `results/tables/` — `models_summary.csv` (45 models, lineage-CV AUC + provenance),
  `kb_overview.csv`, `biomarkers.csv` (3 571 biomarkers × 7 evidence layers),
  `mechanisms.csv`, `cv_comparison.csv` (random-vs-lineage CV),
  `h3_gene_family_overlap.csv` + `h3_summary.json`
- `results/figures/` — the 36 figures (PNG + PDF)
- `results/{organism}/{antibiotic}/` — per-model candidate/evidence outputs
  (07/09/10/11/12/12b/13/14) and `04_evaluation/*.csv`
- `runs/{organism}/{antibiotic}/*/run_metadata.json` — git commit, seed, tool versions
- `config/config.yaml`, `config/registry/*.yaml`, `config/experiments/*/config_*.yaml`
- `environment.lock.yml`, `environment-tools.lock.yml`, `environment-checkm2.lock.yml`
  and `containers/*.def` — the pinned software that produced all of it
- `CITATION.cff`, `README.md`, `METHODOLOGY.md`, this file

```bash
# from the repo root (or $AMR_WORK on TRUBA, where results/ is populated)
tar czf amrk-db_v0.7.1.tar.gz \
  results/kb/amrk.db results/tables results/figures \
  results/*/*/0[4-9]_* results/*/*/1[0-4]_* \
  runs/*/*/*/run_metadata.json \
  config/config.yaml config/registry/*.yaml config/experiments/*/config_*.yaml \
  environment*.lock.yml containers/*.def \
  CITATION.cff README.md METHODOLOGY.md docs/RELEASE_ZENODO.md
```

> **Version:** the release version tracks `kb_schema_version` (**0.7.1**) — `tests/
> test_version_alignment.py` fails the build if `.zenodo.json`, `CITATION.cff`,
> `pyproject.toml` and `config.yaml` disagree with the schema constant. So tag
> `v0.7.1`; do not invent a separate release number.

## Steps

1. **Decide the version.** the release version IS `kb_schema_version` (**0.7.1**) — the alignment test
   requires every declared version to equal the schema constant, so a separate
   release number would break the build. Bump the schema when the shape changes.
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
     python scripts/populate_database.py --organism <org> --antibiotic <ab>
   # … or directly, without re-populating:
   sqlite3 results/kb/amrk.db \
     "UPDATE kb_metadata SET zenodo_doi='10.5281/zenodo.XXXXXXX' WHERE id=1;"
   ```
5. **Write the DOI into the docs:** add it to `CITATION.cff` (`identifiers:` block,
   already stubbed as a comment), `.zenodo.json` is fine as-is, and replace the
   README *"reserved — added on first deposit"* line with the concept DOI badge.
6. **Tag the release** (matches the version):
   ```bash
   git tag -a v0.7.1 -m "AMRK-DB v0.7.1 — 45 models, 6 ESKAPEE organisms, 14 classes, Zenodo 10.5281/zenodo.XXXXXXX"
   git push origin v0.7.1        # PAT, when ready
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
