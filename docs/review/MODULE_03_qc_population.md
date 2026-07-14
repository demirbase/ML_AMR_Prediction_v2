# Modül 3 — QC & Popülasyon Yapısı

> Mercek: **Python'da kal + sertleştir** · **ESKAPEE ölçeklenebilirlik** · **su geçirmez, bilimsel doğruluk**
> İncelenen: `02_kmer_extraction.py`, `02b_global_qc_analysis.py`, `02c_lineage_poppunk.py`, `02d_genome_qc.py`, `02p_kmer_parallel.py`, `03_matrix_construction.py`, `03b_matrix_validation_qc.py`, `scripts/lib/lineage.py`
> Tarih: 2026-07-13

---

## 1. Genel Değerlendirme

Bu modül iki işi yapar: (a) **veri/genom kalite kontrolü** (k-mer spektrum outlier'ı, CheckM2/QUAST) ve (b) tezin en kritik bilimsel savunması olan **popülasyon yapısı → soy-farkındalıklı CV** (PopPUNK + GroupKFold). **Bilimsel çekirdek çok sağlam:** `lib/lineage.py` örnek kalitede (StratifiedGroupKFold, leakage-guard, dtype=str), 02c PopPUNK un-mangling'i loud-fail ediyor (sessiz CV bozulması imkânsız), 02d CheckM2+QUAST must-have'i tam. Adaptive min_support veri boyutuna ölçekleniyor.

**Ana açık = ESKAPEE ölçeklenmesi ve bir bilimsel-dürüstlük riski.** Organizma-seviyesi QC giriş noktaları (02/02b/02p) `--organism` argümanı taşımıyor — organizmayı config.yaml'dan doğrudan okuyorlar (02c/02d ise `--organism` alıyor). Yani yeni bir ESKAPEE organizması için config düzenlemek gerekiyor, paralel per-organizma QC zor. Ayrıca: PopPUNK cluster'ı olmayan bir organizmada 07b'nin **fallback**'e düşüp soy-düzeltmesiz AUC'yi "lineage-CV" diye raporlaması riski var (HANDOFF bunu doğruluyor) — bu **loud + KB'de işaretli** olmalı.

---

## 2. Güçlü Yanlar

- **`lib/lineage.py` (bilimsel pillar):** `StratifiedGroupKFold` (grup bütünlüğü + R/S denge), `no_group_leakage` doğrulama helper'ı, `load_lineage` missing-genome loud-guard + `dtype=str` (audit Issue 5), `collapse_rare_clusters` (belgeli trade-off), `n_splits > n_groups` guard.
- **02c PopPUNK:** organism-level (bir kez, tüm antibiyotikler), `normalize_clusters` PopPUNK'ın `.`→`_` mangling'ini geri çevirip **eşleşmeyende loud-fail** (sessiz CV bozulması yok), `--organism`, config-driven model/refine, `--reuse-db`, `n_clusters < n_splits` uyarısı.
- **02d genom QC (M15 must-have):** CheckM2 (completeness/contamination) + QUAST (N50/contig); prep/post container-chain deseni (pyseer gibi), advisory exclusion (matrisi otomatik bozmaz), pass-rate summary = Methods "data quality" beyanı, `--organism` + eşik CLI'ları.
- **02b spektrum QC:** KMC unique-k-mer üzerinden **IQR outlier** (kontamine/eksik assembly) tespiti, advisory.
- **03 adaptive min_support:** `max(floor, ceil(prevalence × n_genomes))` — küçük/büyük veri setine ölçekleniyor (audit-bilgili, biyolojik olarak güvenli).
- **Tümü `resolve_path` (organism-aware paths)** + `resolve_tool` (poppunk PATH-önce) + `dtype=str` hazard'ı her yerde kapalı.

---

## 3. Problemler (önem sırasıyla)

### Critical / High
- (yok) — bilimsel çekirdek doğru; aşağıdaki riskler süreç/tutarlılık.

### Medium
- **M3-1 — 02/02b/02p `--organism` argümanı YOK:** organizmayı `config['project']['organism']`'dan **doğrudan** (import-time, get_target'sız) okuyorlar. 02c/02d `--organism` alırken bunlar almıyor → yeni ESKAPEE organizması için **config.yaml düzenlemek şart**, paralel per-organizma QC zor. **Fix:** 02/02b/02p'ye `--organism` ekle (02c/02d ile uyum) veya `get_target` kullan.
- **M3-2 — 03/03b `config['project']['target_antibiotic']`'i doğrudan okuyor** (get_target değil) → k-mer baseline yolu env-parametrik değil. Unitig (03u, asıl yol, feature_repr=unitig varsayılan) env-parametrikti; bu yüzden etki **baseline ile sınırlı** ama tutarsız. **Fix:** get_target'e geçir (parallelism refactor 03u'yu aldı, 03'ü atlamış).
- **M3-3 — Soy-CV fallback dürüstlük riski (cross-module → M5'te doğrula):** PopPUNK cluster'ı yoksa 07b non-lineage split'e düşebilir ve AUC'yi yine "lineage-CV" olarak raporlayabilir → yeni organizmada **şişirilmiş, yanlış-etiketli** metrik. HANDOFF "02c'yi çalıştır ki fallback olmasın" diyor = fallback VAR. **Fix (M5):** fallback loud-warn + KB'ye `lineage_cv=false` bayrağı; auc_mean_seeds'i fallback'te "lineage-CV" diye damgalama.

### Low
- **M3-4 — 02c docstring/default drifti:** docstring "refine REQUIRED (bgmm/dbscan alone under-cluster)" + default model=bgmm/refine=on der; ama **validated config `dbscan` + `refine:false`** (config yorumu: bgmm çöktü, refine NaN-fail). Kod config'i doğru okuyor → sadece **docstring + hardcoded default yanıltıcı**. Reconcile et.
- **M3-5 — legacy `paths:` bloğu ölü:** hiçbir script `config['paths']`'i doğrudan okumuyor (hepsi resolve_path); `paths_organism` tüm key'leri taşıdığı için legacy fallback hiç tetiklenmiyor → **silinebilir** (M1 M3-item). Silmeden önce her key'in paths_organism'da olduğunu doğrula.
- **M3-6 — Import-time modül-globalleri:** 02/02b/03/03b organism/target'ı **modül seviyesinde** (import anında) hesaplıyor → test/parametrize zor, config-direct okuma import'ta oluyor. Fonksiyon-kapsamına taşımak (refactor, düşük öncelik).

---

## 4. Düzeltilmesi Gerekenler (madde madde)

1. **02/02b/02p'ye `--organism` ekle** (M3-1) — 02c/02d ile uyumlu; default config'ten.
2. **03/03b'yi `get_target`'e geçir** (M3-2) — k-mer baseline de env-parametrik olsun.
3. **02c docstring + default'ları config gerçeğiyle uyumla** (M3-4) — "dbscan, no-refine (bu veride bgmm+refine başarısız)".
4. **Legacy `paths:` bloğunu kaldır** (M3-5) — paths_organism kapsamını doğrulayıp.
5. (M5'te) **Soy-CV fallback'i loud + KB-işaretli yap** (M3-3).
6. (Refactor, ertelenebilir) modül-global'leri fonksiyona taşı (M3-6).

---

## 5. Refactor Önerileri

- **Organizma çözümünde tek desen:** tüm 02x/03x `get_target(args)` kullansın (CLI `--organism/--antibiotic` > env > config). Şu an karışık (02c/02d CLI, 02/02b/02p config-direct, 03 config-direct, 03u get_target). Tek desen ESKAPEE paralel-koşuyu ve testi netleştirir.
- **`config.yaml` sadeleştirme:** legacy `paths:` kalkınca ~%20 küçülür (M1 ile birleşir).
- **QC orkestrasyonu:** 02/02b/02c/02d organism-level ve bağımsız; bir `make qc ORG=...` hedefi (Makefile) tek komutla per-organizma QC'yi sürebilir (ESKAPEE onboarding'i kolaylaşır).

---

## 6. Bilimsel Eksikler (makale açısından)

- **Soy-CV fallback şeffaflığı (M3-3) — en kritik:** rapor edilen AUC'nin gerçekten soy-düzeltmeli olduğu her model için garanti edilmeli; fallback varsa açıkça belirtilmeli. Aksi halde "lineage-aware" iddiası bir organizmada boşa düşer.
- **PopPUNK model/refine seçiminin organizmaya genellenmesi:** dbscan+no-refine E. coli/K. pneu'da valide; **gram-pozitif S. aureus / A. baumannii'de popülasyon yapısı farklı** → yeni organizma koşulurken PopPUNK ayarı (model, refine, çözünürlük) yeniden doğrulanmalı (n_clusters ≥ n_splits kontrolü zaten var). Run-time metodoloji notu.
- **QC eşiklerinin (CheckM2/QUAST) raporlanması:** pass-rate + dağılımlar Methods'a (02d summary JSON'da var → tabloya taşınmalı).

---

## 7. Literatür Gereksinimi

Bu modül **yerleşik yöntemler (PopPUNK, CheckM2/QUAST, GroupKFold) + mühendislik → derin literatür GEREKMİYOR.** Tek bir **run-time metodoloji notu** (literatür değil, uygulama sırasında karar): yeni ESKAPEE organizmaları için PopPUNK model/refine ayarı organizma-özel doğrulanmalı (Faz 1: A. baumannii, S. aureus). Bu M3 düzeltmelerini bloke etmez.

M3 düzeltmeleri (M3-1,2,4,5) literatür beklemiyor — onayınla uygulanır.

---

## Uygulama durumu (2026-07-13) — UYGULANDI

- **M3-1/M3-2:** `02`, `02b`, `02p`, `03`, `03b` artık `get_target(config=cfg)` kullanıyor → `AMR_ORGANISM`/`AMR_ANTIBIOTIC` env override'ı tanıyor (03u ile aynı; `sbatch --export=...` ile paralel per-organizma/antibiyotik koşu, config düzenlemeden). Doğrulandı: `AMR_ORGANISM=acinetobacter_baumannii` → doğru çözülüyor.
- **M3-4:** 02c `run_poppunk` docstring'i validated karar ile uyumlandı (dbscan alone yeterli; bgmm+refine bu veride başarısız; yeni organizmada doğrula).
- **M3-5:** `config.yaml` legacy `paths:` bloğu **kaldırıldı** (paths_organism superset, resolve_path fallback hiç tetiklenmiyordu). Doğrulandı: config yükleniyor, resolve_path/resolve_tool legacy'siz çalışıyor.
- **Ertelendi:** M3-3 (soy-CV fallback loud+KB-flag) → **M5 (07b)**; M3-6 (import-time global refactor) → düşük öncelik.

**Test durumu:** 6 script syntax OK, `validate_registry` 0 hata, tüm suite **101 passed**, M3 **0 yeni kırık** (9 ön-var kırık M6/M7/M10).

## Sonraki modüllere taşınan notlar
- **M4 (03u unitig):** 03u get_target kullanıyor (env-parametrik ✓) — asıl feature yolu; doğrula.
- **M5 (07b):** **M3-3 soy-CV fallback** burada — loud + KB-flag yapılacak; `no_group_leakage` gerçekten çağrılıyor mu doğrula.
- **M1 carry-over:** legacy `paths:` silme M3-5 ile birleşti.
- **Ön-var 9 test kırığı** hâlâ M6/M7/M10'da.
