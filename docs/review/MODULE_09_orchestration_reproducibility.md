# Modül 9 — Orkestrasyon & Reproducibility

> Mercek: **Python'da kal + sertleştir** · su geçirmez · tam reprodüksiyon
> İncelenen: `slurm/` (Drive rescue), container `*.def`, `environment.yml`, provenance akışı, Nextflow kararı. Tarih: 2026-07-13

## 1. Genel Değerlendirme
Pipeline'ın koşum + tekrarlanabilirlik katmanı. **Karar (kullanıcı): Python'da kal — Nextflow yok.** Bu modülde öksüz Nextflow kaldırıldı; SLURM env-parametrik reçete + container'lar sağlam ama iki reproducibility boşluğu var (slurm git-dışı, lock yok).

## 2. Güçlü Yanlar
- **Env-parametrik SLURM reçetesi:** `sbatch --export=ALL,AMR_ORGANISM=,AMR_ANTIBIOTIC=` → M3'te 02/03 de bu mekanizmaya çekildi = tüm zincir tek desen.
- **Container'lar:** `amr.def` (M0'da geri yüklendi, environment.yml'den türer), `amr-tools.def` (pyseer/quast/amrfinder/resfinder), `amr-checkm2.def`. İki-container deseni (14/02d) farklı env ihtiyaçlarını çözüyor.
- **Provenance:** git_commit+seed+config_hash+card_version her pipeline_runs'ta.

## 3. Problemler
### Critical/High: yok.
### Medium
- **M9-1 [ÇÖZÜLDÜ] — Öksüz Nextflow (M6-2/M4-3):** `08_blast_pipeline.nf` silindi; `environment.yml` nextflow bağımlılığı kaldırıldı (08.py blastn'i subprocess ile koşuyor); 08.py docstring + README badge düzeltildi. Container küçülür, "Python'da kal" ile tam uyum.
- **M9-2 — `slurm/` git-dışı (reproducibility boşluğu):** 33 SLURM scripti sadece TRUBA'da (Drive'a rescue edildi). 4 kanonik `*_env.slurm` (03u/ml/bio/pyseer) repoya girmeli; ~29 tek-seferlik varyant çöp. **Deferred:** Drive'dan çek → repoya `slurm/` (kanonik 4) commit.
- **M9-3 — Lock dosyası yok (M0 carry-over):** `environment.lock.yml`/`requirements.lock.txt` commit'li değil + `amr.def` base `:latest` + pinsiz dep → byte-reproducible değil. **Deferred:** deploy anında conda env resolve edip lock üret + commit.
### Low
- **M9-4 — `amr-gpu.def` ölü** (GPU reddedildi) → `archive/`'e.
- **M9-5 — README/requirements nextflow doc-mention'ları** → **M11** (docs).
- **M9-6 — TRUBA re-populate:** KB 0.6.1 + yeni registry organizmaları (7 ESKAPEE) KB'ye bir sonraki populate'te girecek; TRUBA/Drive KB senkron (deploy).

## 4. Düzeltilecek
1. M9-1 ✓ 2. M9-2 slurm commit (Drive'dan) 3. M9-3 lock (deploy) 4. M9-4 amr-gpu arşiv 5. M9-5 → M11 6. M9-6 deploy.

## 5. Refactor
- `slurm/` kanonik 4 env-script + bir `slurm/README.md` (per-antibiyotik reçete). Tek-seferlik varyantları çöpe.
- `containers/` klasörü + lock dosyaları + base imaj digest-pin → tam reprodüksiyon.

## 6. Bilimsel Eksikler
- Lock + digest-pin olmadan "fully reproducible" iddiası (makale Availability) savunulamaz → M9-3 must-fix (deploy).
- provenance'a unitig-caller/pyseer sürümü (M7-6).

## 7. Literatür: gerekmiyor.

## Uygulama durumu (2026-07-13) — UYGULANDI
- M9-1: `08_blast_pipeline.nf` silindi + environment.yml nextflow kaldırıldı + 08.py docstring + README badge.
- Deferred (Mac-dışı bağımlı): M9-2 slurm commit, M9-3 lock, M9-6 re-populate → deploy.
- Test: 08.py syntax OK, suite yeşil.

## Sonraki modüllere
- **M10:** tests/CI — CI'a mypy + validate_registry + coverage; ön-var kırık yok artık.
- **M11:** README/requirements nextflow-mention temizliği (M9-5), CITATION/zenodo 0.6.1, cv_method surface.
