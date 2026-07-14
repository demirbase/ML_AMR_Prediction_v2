# Modül 10 — Testler & CI

> Mercek: Python+sertleştir · su geçirmez
> İncelenen: `tests/` (17 dosya), `pytest.ini`, `.github/workflows/ci.yml`. Tarih: 2026-07-13

## 1. Genel Değerlendirme
Kalite kapıları. **Sağlam ve artık TAMAMEN YEŞİL:** 17 test dosyası (unit/smoke/integration), marker'lı hızlı-varsayılan koşu, çok-sürümlü CI (3.10-3.12). M1-M9 boyunca bulunan 9 ön-var kırığın hepsi kapandı → **suite 109 passed, 0 kırık.** M10'da CI'a registry-guard + mypy eklendi.

## 2. Güçlü Yanlar
- **Kapsam:** saf-mantık (lib, bvbrc cleaning, lineage CV, concordance), smoke (her numaralı script import), integration (sentetik uçtan-uca), KB queries, populate idempotency.
- **Marker'lar** (unit/smoke/integration/slow) → `pytest` hızlı varsayılan; ağır testler opt-in (config'i mutasyona uğratmaz).
- **CI** her push/PR, 3 Python sürümü, ruff + pytest.
- **Bu inceleme boyunca 9 ön-var kırık düzeltildi** (test_15 organism, test_kb_queries overlap.organism) + yeni testler (bvbrc nanargmax/policy, registry meta accessor'ları).

## 3. Problemler
### Critical/High: yok.
### Medium
- **M10-1 [ÇÖZÜLDÜ] — CI boşlukları (M0):** CI yalnız ruff+pytest'ti. **Eklendi:** `validate_registry.py` (blocking — registry↔iç tutarlılık bekçisi) + `mypy` (advisory, continue-on-error). Artık registry driftı CI'da yakalanır.
### Low
- **M10-2 — Coverage raporu yok:** opsiyonel; `pytest-cov` + eşik eklenebilir.
- **M10-3 — integration/slow CI'da değil:** KMC/xgboost gerektiriyor → maintainer/HPC'de koşuyor (belgeli, kabul).
- **M10-4 — validate_registry unit-test'i yok:** CI step'i var; istenirse `tests/`e de eklenebilir (çift güvence).

## 4. Düzeltilecek
1. M10-1 ✓ 2. M10-2 coverage (opsiyonel) 3. gerisi kabul.

## 5. Refactor
- `pytest-cov` + `--cov-fail-under` eşiği (yayın-hazırlık).
- `validate_registry` için bir smoke test (import + 0-error assert).

## 6. Bilimsel Eksikler: yok (mühendislik). Suite yeşilliği = "su geçirmez" temel kanıtı.

## 7. Literatür: gerekmiyor.

## Uygulama durumu (2026-07-13) — UYGULANDI
- M10-1: CI'a `validate_registry` (blocking) + `mypy` (advisory) step'leri.
- Doğrulama: ci.yml valid, validate_registry 0 hata, suite **109 passed**.

## Sonraki modüllere
- **M11:** README/CITATION/zenodo 0.6.1 + nextflow-mention temizliği (M9-5) + cv_method/organism surface + METHODOLOGY/ROADMAP 0.6.x + 21-model paneli.
