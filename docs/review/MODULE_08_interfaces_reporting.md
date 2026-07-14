# Modül 8 — Arayüz & Rapor

> Mercek: Python+sertleştir · su geçirmez · yayın-hazırlık
> İncelenen: `kb_api.py` (FastAPI), `kb_app.py` (Streamlit), `kb_tables.py` (tidy CSV), `kb_figures.py` (tez figürleri), `kb_report.py` (Markdown özet), `lib/kb_queries.py`. Tarih: 2026-07-13

## 1. Genel Değerlendirme
KB'nin dış yüzleri (FAIR erişim). **Sağlam:** REST API (FastAPI, OpenAPI /docs, CORS), Streamlit explorer, tidy-CSV export, 9 tez figürü, tek-komut Markdown rapor; hepsi `kb_queries` (saf sqlite3) üzerinden. En kritik iki hata (bayat DB path + organism-karışması) düzeltildi.

## 2. Güçlü Yanlar
- **`kb_queries` saf sqlite3** → API/app/report/tables tek sorgu katmanını paylaşıyor (DRY, test-edilebilir).
- **FastAPI:** `/kmers`, `/kmers/{sequence}` (tam kanıt zinciri), `/overlap`, `/stats`, `/metadata` (FAIR makine-okunur); OpenAPI auto-docs; `AMR_KB_DB` env ile DB seçimi.
- **Streamlit:** 5+ sekme (biyobelirteç/kanıt/provenance/H3/M13 + Ham tablolar), varsayılan doğru birleşik KB.
- **kb_figures/tables/report:** tez çıktıları tek komutla üretilebilir.

## 3. Problemler
### Critical/High: yok.
### Medium
- **M8-1 [ÇÖZÜLDÜ] — `kb_api` bayat DB path:** default `results/ecoli/kb/amrk.db` (eski per-organizma) idi → API birleşik KB'yi bulamaz/yanlış KB açardı. **Fix:** `results/kb/amrk.db` (kb_app ile uyumlu, 2 yer: kod + docstring).
- **M8-2 [ÇÖZÜLDÜ] — `get_overlap` organism-karışması (M7-4):** overlap 0.6.0'da organism-aware ama sorgu süzmüyordu → çok-organizmalı KB'de aynı ilaç çifti (ör. cipro/levo ecoli+kpneu) karışırdı. **Fix:** `get_overlap(...,organism=None)` + `/api/v1/overlap?...&organism=` param (backward-compat: None=tümü).
### Low
- **M8-3 — CORS `allow_origins=["*"]`** (kb_api) — lokal araştırma API'si için sorun değil; production dağıtımında kısıtla (not).
- **M8-4 — cv_method arayüzde gösterilmiyor:** yeni `models.cv_method` (M7) kb_app/kb_report/figürlerde surface edilirse "bu AUC honest lineage-CV" görünür olur (opsiyonel iyileştirme).

## 4. Düzeltilecek
1. M8-1 ✓ 2. M8-2 ✓ 3. M8-3 production notu 4. M8-4 opsiyonel (cv_method surface).

## 5. Refactor
- Streamlit/API'ye organism seçici (overlap artık organism-aware) — çok-organizma UX'i netleşir.
- cv_method'u models tablosu görünümüne ekle (kb_app/kb_report).

## 6. Bilimsel Eksikler
- FAIR erişim tam (API+metadata); cv_method surface edilince "tüm AUC lineage-CV" iddiası UI'dan da görünür olur.

## 7. Literatür: gerekmiyor.

## Uygulama durumu (2026-07-13) — UYGULANDI
- M8-1: kb_api DB path → `results/kb/amrk.db` (kod+docstring).
- M8-2: `get_overlap` + `/overlap` endpoint organism-param; doğrulandı (kpneu-only sorgu çalışıyor, `organism` alanı dönüyor).
- Test: syntax OK, suite **109 passed** (0 kırık).

## Sonraki modüllere
- **M9:** orkestrasyon (slurm/, 08.nf sil, container/lock, reproducibility); TRUBA re-populate.
- **M11:** cv_method/organism surface (kb_app/report), CORS production notu.
