# Modül 6 — Biyolojik Doğrulama (08–16)

> Mercek: Python+sertleştir · ESKAPEE · su geçirmez · bilimsel doğruluk
> İncelenen: `08_blast_annotation` (CARD/NCBI), `09_biological_summary` (tier), `10_kmer_background_frequency` (R/S prevalans), `11_variant_snp_check` (SNP), `12/12b_permutation` (MDA/label), `13/13b_stability` (CPSS), `14_pyseer_lmm` (GWAS), `15_cross_antibiotic` (H3), `16_external_concordance` (M13). Tarih: 2026-07-13

## 1. Genel Değerlendirme
7 katmanlı kanıt zincirinin üretildiği modül — tezin "açıklanabilirlik" iddiasının kalbi. **Sağlam ve tutarlı:** çoğu script env-parametrik (get_target; 12b/13 05'in globallerini miras alıyor, 15 organism-CLI — beklenen), **hiç hardcoded path yok** (08 organizmayı registry'den, 09 NCBI e-mail'i asla hardcode etmiyor = ToS/ban koruması, 14 config-resolved paths), tier eşikleri config'te (citeable), 16 leakage-free head-to-head. Bilimsel katmanlar figür incelemelerinde de doğrulandı (05_significance, 06_evidence_layers, 07_external).

## 2. Güçlü Yanlar
- **Env-parametrik zincir** (08/09/10/11/12/13b/14/16 get_target) → ESKAPEE paralel-koşu.
- **08 BLAST:** blastn subprocess (CARD local + NCBI remote); task auto-seçimi (kısa unitig→blastn, uzun→blastn), organizma registry'den; SIGXCPU workaround belgeli.
- **09 tier:** confirmed/candidate/weak eşikleri config'te; reification güvencesi.
- **14 pyseer:** iki-container (amr.sif prep/post + amr-tools.sif CLI); config-resolved paths (hardcode yok).
- **16 concordance:** amrfinder/resfinder vs model, aynı held-out test, FDA ME/VME + κ; amrfinder_keywords registry-driven.
- **NCBI e-mail loud-required** (placeholder ban riskini engelliyor).

## 3. Problemler
### Critical/High: yok — biyolojik zincir doğru ve figürlerle tutarlı.
### Medium
- **M6-1 [ÇÖZÜLDÜ] — `test_15` ön-var kırığı:** `populate_overlap` imzası 0.6.0'da `organism` aldı, test 4-arg çağırıyordu (3 FAIL). **Fix:** testin 4 çağrısına `organism="ecoli"` eklendi → `test_15` 9/9 yeşil. (9 ön-var kırıktan 3'ü kapandı.)
- **M6-2 — 08_blast_pipeline.nf öksüz:** 08.py BLAST'ı **doğrudan blastn subprocess** ile koşuyor; `.nf` çağrılmıyor ama README/HANDOFF/docs + 08.py docstring referans veriyor. "Python'da kal" kararıyla `.nf` gereksiz. **→ M9** (Nextflow kararı + README + `environment.yml` nextflow-dep + 08.py docstring birlikte temizlenecek).
### Low
- **M6-3 — 08.py docstring bayat** ("runs via Nextflow pipeline") — gerçekte blastn subprocess. M6-2 ile birlikte düzeltilecek.
- **M6-4 — 15 `_drug_family`** artık carbapenems/monobactams'ı doğru topluyor (M1'de güncellendi); test yeşil.

## 4. Düzeltilecek
1. M6-1 ✓ (test_15). 2. M6-2/M6-3 → M9 (nf + docs + dep). 3. gerisi temiz.

## 5. Refactor
- 08.nf sil + nextflow bağımlılığını `environment.yml`'den çıkar (M9) → "Python'da kal" kararıyla tam uyum, container küçülür.

## 6. Bilimsel Eksikler
- Zincir bilimsel olarak tam (7 katman); tek not: mechanism_type/who_aware artık registry'den (M1) → populate M7'de registry'den okumalı (M1'de kod hazır, KB re-populate/migrate M7).
- unitig-caller/pyseer/CARD sürümleri provenance'ta mı → M7'de doğrula.

## 7. Literatür: gerekmiyor.

## Uygulama durumu (2026-07-13) — UYGULANDI
- M6-1: `test_15` 4 çağrıya `organism` eklendi → 9/9 yeşil (ön-var kırıklar 9→6).
- Test: suite **104 passed** (test_15 dahil), kalan 6 kırık = `test_kb_queries` (overlap.organism NOT NULL → M7).

## Sonraki modüllere
- **M7:** KB `models.cv_method` kolonu (M5) + populate registry-meta okuma (M1) + `test_kb_queries` fixture (overlap.organism) + provenance (tool sürümleri).
- **M9:** 08.nf sil + README/env/docstring (M6-2/M6-3).
