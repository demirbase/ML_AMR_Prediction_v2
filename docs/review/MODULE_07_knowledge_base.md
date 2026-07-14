# Modül 7 — Bilgi Tabanı (KB)

> Mercek: Python+sertleştir · su geçirmez · reproducibility
> İncelenen: `lib/kb_schema.py`, `populate_database.py`, `migrate_kb_050.py`, `lib/kb_queries.py`. Tarih: 2026-07-13

## 1. Genel Değerlendirme
KB'nin şema + doldurma katmanı. **Sağlam:** 13 tablo, provenance çıpası (`pipeline_runs`: git+seed+CARD+config_hash), FK bütünlüğü + idempotency (unique index + INSERT OR IGNORE), additive migrate deseni (FK re-insert'i atlar), `kb_queries` saf sqlite3 (test-edilebilir). M1'de populate registry-meta okumaya çevrildi; M7'de `cv_method` eklendi → auc_mean_seeds'in dürüstlüğü artık makine-okunur.

## 2. Güçlü Yanlar
- **Provenance-first:** her kayıt `pipeline_runs`'a (git_commit + random_seed + config_hash + card_version) bağlı → tam tekrarlanabilirlik.
- **Idempotent populate** (`ensure_unique_indexes` + INSERT OR IGNORE) — 09+13b çift-insert audit'i kapalı.
- **Additive migrate** (`migrate_kb_050`): model re-insert etmeden backfill (FK sorununu atlar).
- **Registry tek-kaynak (M1):** populate_organisms/antibiotics_meta artık registry'den (gram/phylum/mechanism/AWaRe) okuyor — hardcode yok.

## 3. Problemler
### Critical/High: yok.
### Medium
- **M7-1 [ÇÖZÜLDÜ] — `test_kb_queries` ön-var kırığı:** fixture overlap'i `organism`suz insert ediyordu (0.6.0 NOT NULL) → 6 ERROR. **Fix:** fixture'a `organism='ecoli'` eklendi → **tüm suite 109 passed, 0 kırık.**
- **M7-2 [ÇÖZÜLDÜ] — cv_method KB'de yoktu (M5/M3-3):** `models.cv_method` kolonu eklendi (schema 0.6.0→0.6.1); populate 07b summary'den okuyor; Mac KB ALTER+backfill (21 model = `lineage_group_kfold_5fold`, hepsi PopPUNK'lı). Artık her model "honest lineage-CV mi fallback mı" KB'den sorgulanabilir.
### Low
- **M7-3 — KB organisms tablosu 2 organizma:** populate registry'den okuyor ama mevcut Mac KB re-populate edilene kadar sadece ecoli/kpneu içeriyor (7 ESKAPEE organizması registry'de, KB'ye bir sonraki populate'te girecek). Build-artifact, beklenen.
- **M7-4 — `kb_queries.get_overlap` organism-filtresiz:** overlap 0.6.0'da organism-aware ama get_overlap sorgusu organism süzmüyor → çok-organizmalı KB'de aynı ilaç çifti organizmalar arası karışabilir. **→ M8** (API/query katmanı).
- **M7-5 — Public sürüm drifti sürüyor:** schema 0.6.1 ama CITATION.cff/.zenodo.json hâlâ 0.4.0 (M0 bulgusu) → **M11** yayın-hazırlığında senkronla.
- **M7-6 — TRUBA/Drive KB kopyaları 0.6.0** (cv_method'suz) → re-populate/migrate ile senkron (deploy).

## 4. Düzeltilecek
1. M7-1 ✓ 2. M7-2 ✓ 3. M7-4 → M8 4. M7-5 → M11 5. M7-6 deploy.

## 5. Refactor
- `get_overlap` + kb_queries'e organism parametresi (M8).
- migrate scriptleri tek `migrate_kb.py`'de birleştirilebilir (0.6.x zinciri).

## 6. Bilimsel Eksikler
- Provenance zengin (git+seed+CARD); **unitig-caller/pyseer/bcalm sürümleri** pipeline_runs'ta mı — eklenmezse tam tekrarlanabilirlik için not (kmc/xgboost var). Yayın öncesi (M11) doğrula.
- cv_method sayesinde artık makalede "tüm rapor edilen AUC lineage-CV" iddiası KB'den kanıtlanabilir.

## 7. Literatür: gerekmiyor.

## Uygulama durumu (2026-07-13) — UYGULANDI
- M7-1: test_kb_queries fixture (overlap.organism) → **suite 109 passed, 0 kırık** (tüm 9 ön-var kırık kapandı).
- M7-2: `models.cv_method` (schema 0.6.1); populate okuyor; Mac KB backfill (21 = lineage_group_kfold_5fold).
- Sürüm: KB_SCHEMA_VERSION + config + KB metadata → 0.6.1.

## Sonraki modüllere
- **M8:** kb_queries/API organism-filtresi (M7-4); kb_app/kb_figures/kb_report cv_method'u gösterebilir.
- **M11:** CITATION/zenodo → 0.6.1 (M7-5); provenance tool sürümleri (M7-6).
- **M9:** TRUBA re-populate ile KB senkron.
