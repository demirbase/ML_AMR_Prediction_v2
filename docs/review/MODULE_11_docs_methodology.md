# Modül 11 — Dokümantasyon & Bilimsel Metodoloji

> Mercek: yayın-hazırlık · doğru kapsam/sürüm · su geçirmez
> İncelenen: `README`, `CITATION.cff`, `.zenodo.json`, `pyproject`, `METHODOLOGY`, `ROADMAP`, `docs/`. Tarih: 2026-07-13

## 1. Genel Değerlendirme
Yayın/atıf yüzeyi. **Altyapı olgun** (cff+zenodo+DOI slotu, METHODOLOGY, ROADMAP, KB_ACIKLAMA/KAVRAMLAR/AKIS_SEMASI). Ana sorun M0'dan beri bilinen **sürüm+kapsam drifti**: public metadata "0.4.0 / E. coli only" derken gerçek "0.6.1 / 2 organizma / 21 model". M11'de public-facing sürüm+kapsam düzeltildi; kalan derin doc-rewrite'lar deposit-finalize'a bırakıldı.

## 2. Güçlü Yanlar
- Yayın altyapısı hazır: `CITATION.cff`, `.zenodo.json` (DOI slotu), MIT+CC-BY, ROADMAP/METHODOLOGY, danışman dokümanları (KB_ACIKLAMA/KAVRAMLAR + akış şeması), inceleme raporları (docs/review/).
- METHODOLOGY reification-safe wording (S10), ROADMAP §0 literatür-pivot kararları.

## 3. Problemler
### Critical/High: yok.
### Medium
- **M11-1 [ÇÖZÜLDÜ] — Public sürüm+kapsam drifti (M0/M7 carry-over):** CITATION/zenodo/pyproject → **version 0.6.1** + başlık/açıklama "Escherichia coli" → **"ESKAPEE pathogens (E. coli, K. pneumoniae)"**; zenodo notes schema 0.6.1 (21 model/2 org). Atıf yapan ilk gördüğü alanlar artık doğru.
### Low
- **M11-2 — Derin prose güncellemesi (deferred):** README abstract ("E. coli across four antibiotics") + zenodo `description` gövdesi ("in E. coli", "schema 0.4.0") + METHODOLOGY/ROADMAP'e 0.6.x şema + 21-model paneli + nextflow-mention temizliği (M9-5) → **Zenodo deposit-finalize (M10 roadmap) + tek doc-pass'te** yapılacak (çok satırlı prose, tez-yazımıyla eşzamanlı).
- **M11-3 — SCALE_MLOPS_PLAN.md** referansları geçerli (M1'de doğrulandı); güncel kalsın.
- **M11-4 — provenance tool sürümleri** (unitig-caller/pyseer) pipeline_runs'a (M7-6) — yayın Availability için.

## 4. Düzeltilecek
1. M11-1 ✓ (version+scope). 2. M11-2 deposit-finalize doc-pass. 3. M11-4 provenance (deploy).

## 5. Refactor
- Tek `VERSION` kaynağı (pyproject) + release-time otomasyon (bump-my-version) → CITATION/zenodo elle-drift biter (M0 §5).
- README'yi 0.6.1 / çok-organizma / Python-orkestrasyon (Nextflow-suz) gerçeğine göre tek seferde yenile.

## 6. Bilimsel Eksikler
- Kapsam beyanı (çok-organizma/ESKAPEE, cross-phylum hedefi) artık başlıklarda doğru → reviewer "kod E. coli diyor" çelişkisi kapandı (public-facing).
- Container+lock (M9-3) olmadan "fully reproducible" iddiası → Availability bölümünde must-fix (deploy).
- cv_method (M7) sayesinde "tüm rapor edilen AUC lineage-CV" iddiası KB'den kanıtlanabilir → Methods'a işle.

## 7. Literatür: gerekmiyor (kapsam kararı ESKAPEE1.md ile verildi).

## Uygulama durumu (2026-07-13) — UYGULANDI
- M11-1: CITATION.cff + .zenodo.json + pyproject → version 0.6.1 + başlık/kapsam ESKAPEE (E. coli, K. pneumoniae) + zenodo notes schema 0.6.1.
- Doğrulama: 3 metadata dosyası parse OK, suite 109 passed.
- Deferred: M11-2 derin prose (deposit-finalize), M11-4 provenance sürümleri (deploy).
