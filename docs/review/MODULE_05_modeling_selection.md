# Modül 5 — Modelleme & Seçim (04/05/06/07/07b)

> Mercek: Python+sertleştir · ESKAPEE · su geçirmez · bilimsel doğruluk
> İncelenen: `04_optimization` (HPO), `05_model_training` (full-data), `06_evaluation` (metrik+CI+head-to-head), `07_explainability` (gain top-N), `07b_feature_stability` (lineage-CV + CPSS). Tarih: 2026-07-13

## 1. Genel Değerlendirme
Tahmin motoru + rapor edilen metriğin üretildiği modül. **Sağlam:** 04-07b hepsi env-parametrik (`get_target` ✓), 07b tezin dürüst metriğini (lineage-aware StratifiedGroupKFold → auc_mean_seeds) üretiyor, 06 operating-point + bootstrap CI + leakage-free model-vs-araç head-to-head. **İki-metrik tasarımı bilinçli:** 06 tek-split `roc_auc` (operating point) vs 07b lineage-CV `auc_mean_seeds` (rapor edilen). En kritik risk = soy-CV fallback'in şeffaflığıydı; bu modülde sertleştirildi.

## 2. Güçlü Yanlar
- **Env-parametrik** 04/05/06/07/07b → ESKAPEE paralel-koşu.
- **07b lineage-aware CV** (lib/lineage StratifiedGroupKFold) = rapor edilen dürüst metrik; fallback loud.
- **04 paralel Optuna** (SLURM CPU-verimi), **05 full-data early-stopping** (ExtMem, 384GB node), **06 bootstrap %95 CI + FDA ME/VME + leakage-free head-to-head** (test genomları eğitimde değil).
- **CPSS kararlılık** (07b) — PFER-sınırlı biyobelirteç seçimi.

## 3. Problemler
### Critical/High: yok.
### Medium
- **M5-1 [ÇÖZÜLDÜ] — Soy-CV fallback şeffaflığı (M3-3):** 07b `cv_method`'u yalnız stdout'a yazıyordu; summary CSV'ye ve KB'ye girmiyordu → auc_mean_seeds "lineage-CV mi fallback mı" ayırt edilemiyordu (schema comment bile belirsiz). **Fix:** `cv_method` summary CSV'ye eklendi + fallback'te **prominent ⚠ uyarı** (AUC lineage-inflated olabilir, 02c çalıştır). **Kalan:** KB `models.cv_method` kolonu → **M7**.
- **M5-2 [ÇÖZÜLDÜ] — `no_group_leakage` çağrılmıyordu:** helper tanımlıydı ama leakage-free özelliği yalnız StratifiedGroupKFold'a güveniliyordu. **Fix:** her fold için `no_group_leakage` assert'i eklendi (leakage'de RuntimeError).
### Low
- **M5-3 — 06/HPO chunk-split, lineage-aware değil (audit #22):** tasarım gereği (operating point). **Tek rapor edilen genelleme metriği = auc_mean_seeds (07b)** — Methods'ta net vurgulanmalı (05_significance figürü de 12b'nin chunk-split AUC'sini gösteriyor, bkz. figür incelemesi).
- **M5-4 — `10_repeated_holdout_summary` dosya adı bayat:** artık lineage-CV içeriyor; yeniden adlandırma populate'e sıçrar → ertelendi/not.

## 4. Düzeltilecek
1. M5-1 ✓ (KB kolonu → M7). 2. M5-2 ✓. 3. M5-3 Methods notu. 4. M5-4 ertelendi.

## 5. Refactor
- KB'ye `cv_method` kolonu (M7) + populate summary'den okusun → her modelin metriği "honest lineage-CV mi" makine-okunur olur.
- `10_repeated_holdout_summary` → `10_cv_summary` (M7/populate ile birlikte).

## 6. Bilimsel Eksikler
- **Rapor edilen metrik = auc_mean_seeds (lineage-CV)** her yerde net olmalı; roc_auc/12b-AUC (tek-split) ile karıştırılmamalı (figür incelemelerinde de not edildi).
- Fallback modeller (varsa) makalede "lineage-CV değil" diye işaretlenmeli — artık cv_method ile mümkün.

## 7. Literatür: gerekmiyor.

## Uygulama durumu (2026-07-13) — UYGULANDI
- M5-1: `cv_method` → summary CSV + prominent fallback uyarısı.
- M5-2: `no_group_leakage` assert'i (07b lineage yolu).
- Test: 07b syntax OK, suite 101 passed, 0 yeni kırık.

## Sonraki modüllere
- **M6:** 08-16 biyolojik doğrulama; 08.nf vs 08.py (M4-3); ön-var `test_15` kırığı (populate_overlap logger).
- **M7:** KB `models.cv_method` kolonu + populate; `10_cv_summary` rename; ön-var `test_kb_queries` (overlap.organism).
