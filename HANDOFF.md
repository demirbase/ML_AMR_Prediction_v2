# AMR k-mer Knowledge Base — Project Handoff Document

> **Repo:** `ML_AMR_Prediction_v2` · branch **`main`** · **HEAD `5e3101f`+** (pushed to `github.com/demirbase/ML_AMR_Prediction_v2`). ✅ 45-model KB COMPLETE — **READ §0.-8 FIRST.**
> **Local (Mac) path:** `~/Desktop/IU_master/projects/ML_project_kopyasi`

---

# §0.-8 — LATEST STATE (2026-07-31) — ✅ 45-MODEL KB COMPLETE — READ FIRST, supersedes ALL below

> **Repo HEAD `5e3101f`**, hepsi push'lu. KB `results/kb/amrk.db` (TRUBA) schema **0.7.1**, **populate 45/45 FAIL=0**. §0.-7 ve altı tarihsel arka plan.

### ✅ NİHAİ KB (tezin çekirdek deliverable'ı)
- **45 model · 6 ESKAPEE · 14 sınıf · hepsi `lineage_group_kfold_5fold`.** lineage-CV AUC: min 0.429 · ort **0.842** · max 0.975.
- **Organizma:** Efm 0.933 · Ec 0.917 · Sa 0.836 · Kp 0.834 · Ab 0.759 · Pa 0.755. (Ab/Pa düşük = intrinsik/klonal; lineage-CV dürüstçe ifşa.)
- **⭐ MANŞET — A. baumannii ceftazidime: holdout 0.985 vs lineage-CV 0.429** (Jaccard 0.04) = ders-kitabı klonal-confounding, lineage-CV'nin var-oluş kanıtı (E3 §8 random-vs-lineage tablosunun canlı örneği).
- **evidence_tier:** confirmed 349 · candidate 942 · weak 1920 · none 337 · **strong_novel 23**.
- **provenance:** 8/9 kolon 45/45 dolu (card·xgboost·unitig-caller·poppunk·graph-tool·blast·pyseer). **bcalm 0/45 = dürüst NULL** (bcalm sürüm CLI'si yok; `_tool_version` hata-string'i reddediyor).
- **Biyoloji doğrulandı:** MRSA→mecA · VRE→vanA kümesi · ESBL→CTX-M · KPC/NDM · A.b OXA hepsi geldi. CARD-conf=0 olanlar = mutasyon-tabanlı/intrinsik (gyrA/pbp5/porin) → CARD-nükleotid yakalamaz, tutarlı.

### 🐛 FAZ-B'DE ÖĞRENİLENLER (kritik)
1. **bio, Faz2a'ya BAĞIMLI** — bio adım 10 `07_kb_candidates`'i (09'un çıktısı) okur; paralel-dal DEĞİL. Sıra: 03u→ML→**Faz2a(07-09)**→bio(10-13b)→pyseer(14)→populate.
2. **`5e3101f` — store-subset 03u `unitigs.rtab` yazmıyordu** (fallback yazıyordu) → 45 pyseer "unitig Rtab not found" ile patladı. Fix eklendi. **rtab BÜYÜK** (ecoli ~36GB, yazımı ~50dk python-loop yavaş — opt. adayı) ama disk 3.5PB, sorun yok.
3. **SLURM "launch failed requeued held"** tekrarladı (geçici bozuk node) — `scontrol release` çalışmadı ("interaktif iş" hatası) → **scancel + taze resubmit** (63 idle node varken düzgün node'a düşer) çözdü.
4. Faz2a NCBI seri ~saatler, **detached screen**'de gözetimsiz koştu (`faz2a_all.sh` + `faz2a_retry.sh` 3-geçiş self-healing, `07_kb_candidates` eksik olanı yeniden koşar). 09 markdown-rapor efetch'i kırılgan (kozmetik, kb_candidates rapordan önce yazılır).

### ⏭️ YAYIN-ÖNCESİ KALAN
- **Kanonik slurm'leri COMMIT et** (TRUBA'da untracked): `run_03uml_env` · `run_03ubio_env` · `run_m15_qc_env` · `run_03u_build_env` · `faz2a_all.sh` · `faz2a_retry.sh` · `run_pyseer_env`(AMR_PYSEER_VERSION patch'li).
- KB figürleri/tabloları: `scripts/kb_figures.py` + `kb_tables.py` (tez showcase).
- Analiz kalemleri: MIC duyarlılık · coğrafi hold-out · random-vs-lineage-CV tablosu (**A.b ceftazidime hazır örnek**) · H3 hypergeometric (E2 çerçevesi).
- Zenodo deposit + DOI. Suite/`test_version_alignment` yeşil tut.

---

# §0.-7 — LATEST STATE (2026-07-19) — PİLOT UÇTAN UCA GEÇTİ · 4 GİZLİ BUG DÜZELTİLDİ · TAM KOŞUYA HAZIR — READ FIRST, supersedes ALL below

> **Repo HEAD `06a38f8`** (+ bu handoff commit'i), hepsi push'lu. TRUBA senkron (`git fetch && git reset --hard origin/main`). **KB şeması 0.7.1.** §0.-6 ve altı geçerli arka plan; bu bölüm pilot-sonrası GÜNCEL durum.

### ✅ PİLOT (A. baumannii amikacin) UÇTAN UCA KOŞTU — 3 DOĞRULAMANIN HEPSİ GEÇTİ
Zincir: 03u → ML(04-07b) → 07 → 08 → 09 → bio(10-13b) → pyseer(14) → populate.
- **#1 lineage-CV ✓** — `models.auc_mean_seeds=0.916`, `cv_method='lineage_group_kfold_5fold'` (fallback DEĞİL). roc_auc=0.972 ayrı kolon (holdout, KB metriği DEĞİL). 5 fold 0.735/0.962/0.979/0.925/0.980 — GC2 klonunu tutan fold dürüst düşük.
- **#2 provenance ✓** — pipeline_runs 8/9 gerçek sürüm (poppunk 2.7.8, unitig-caller 1.3.2, graph-tool 3.0, blast 2.17, kmc, xgboost, card, pyseer 1.4.1). bcalm = dürüst None (sürüm CLI'si yok; unitig-caller'a gömülü + lockfile'da pinli).
- **#3 evidence_tier ✓** — weak 50/confirmed 21/candidate 8/**strong_novel 1**/none 1. Biyoloji: APH(3')-VIa/VIb aminoglikozit %100 CARD, recovery 0.64, H2_pass, novel_fraction 0.36, 10 stabil novel.

### 🐛 PİLOTTA BULUNAN + DÜZELTİLEN 4 GİZLİ BUG (hepsi tezi vururdu; testler göremiyordu)
1. **`7ef05f7` — 03u lineage kesişimi.** `--qc-db` (02c) sketch DB'den 38 genom budar → kümesiz kalırlar; 03u model genomlarını AMR etiketinden seçtiği için 34 amikacin genomu kümesiz → **07b SESSİZCE lineage-CV'den holdout'a düşüyordu**. `--no-qc` denendi → A.b **%96.7 mega-kümeye ÇÖKTÜ** (qc ŞART, redundant değil). Fix: 03u seçtiği genomları `poppunk_clusters.csv` ile kesiştirir (kümesiz=QC-fail → modele girmez). **config `lineage.qc: true` KALIR.**
2. **`c3ecf9c` — 08 nextflow → subprocess.** 08 main() `nextflow run <olmayan .nf>` çağırıyordu; nextflow yeni amr.sif'te yok → 08 hiç koşamıyordu. Fix: CARD local + NCBI remote blastn doğrudan subprocess (outfmt-6, 09'un okuduğu kolonlar). nextflow bağımlılığı kalktı.
3. **`06a38f8` — provenance card/bcalm/pyseer.** bcalm probe `--version` (yanlış flag) hata string'i saklıyordu → `-version` + hata-reddi (bcalm yine None, dürüst). card_version NULL → `card_db_version.txt`'ten okunuyor. pyseer_version NULL → 14 artık `AMR_PYSEER_VERSION` env'ini tercih ediyor.
4. **`run_pyseer_env.slurm` PATCHED (TRUBA)** — pyseer sürümünü yakalayıp `AMR_PYSEER_VERSION` ile 14 post'a geçiyor. ⚠️ **UNTRACKED** — kanonik slurm'lerle commit edilmeli.

### ⚠️ AYRICA BULUNANLAR (tam koşu öncesi)
- **`amr_phenotypes.csv`**: `00_prepare_metadata` yeni organizmalarda (Ab/Pa/Efm) koşmamıştı (pilotta A.b üretildi). **Ec/Kp/Sa phenotypes Jul-5 BAYAT** (filtre-düzeltmesi öncesi).
- **02d/M15 CheckM2 QC koşmadı** — 03u "global_qc_outliers.csv yok" dedi; registry'nin QC gate'i, pilot atladı.
- **09 kırılgan** — NCBI Entrez efetch (markdown rapor) uzun + internet; SSH koparsa ölür. `config.ncbi.entrez_email` SET DEĞİL (NCBI ban riski). KB-kritik çıktılar (07_kb_candidates, 08_metrics) rapordan ÖNCE yazılır → populate etkilenmez.

### ⏭️ TAM KOŞU ÖNCESİ CHECKLIST (temiz koşu)
1. **6 organizmada `00_prepare_metadata`** yeniden (3 eski bayat).
2. **02d/M15 CheckM2 QC** → her organizma `global_qc_outliers.csv` (03u eler; 03u-kesişim kümesizleri ayrıca eler → çift güvenli).
3. **Lineage etiketleri:** 6'sı da qc:true ÜRETİLDİ + geçerli (A.b refine 39.5%). Re-cluster GEREKMEZ.
4. **`03u --build-db`** organizma başına BİR kez (hiç kurulmadı). Build slurm YOK → yazılmalı (`-c40 --mem=300G`, --build-db, yalnız AMR_ORGANISM).
5. **09'u nohup/screen'de** + `config.ncbi.entrez_email` set.
6. **Kanonik slurm'leri commit** (run_03u_env, run_ml_env, run_bio_env, run_pyseer_env[patched], build) — TRUBA'da untracked.
7. **45-model per-antibiyotik zincir** env-parametrik, KESİN SIRALI (memory `amr-truba-gotchas`).

### 🔧 OPERASYONEL
- **barbun QOS: -c 20-40 ZORUNLU** (-c16 reddedildi). node 40-çekirdek/384GB/9.5GB-core/3-gün. Dolu cluster'da -c20 --mem=180G küçük işleri hızlı sokar (ML MaxRSS ~21G). Alt kuyruklar orfoz/hamsi/orkinos/sardalya/palamut. (memory güncel.)
- populate, NULL sürüm kolonlarını `run_metadata.json`'dan okur (04 yazar) → kod fix'i ancak 04 yeniden koşunca yansır (tam koşuda otomatik). card/pyseer populate-anında ayrı kaynaktan okunur.

---

# §0.-6 — LATEST STATE (2026-07-16) — DEPLOY-ÖNCESİ HAZIR · CONTAINER'LAR KİLİTLİ · GENOMLAR İNDİ · PANEL KÜRE · LINEAGE AYARLANDI — READ FIRST, supersedes ALL below

> **Repo HEAD `b9c0720`**, hepsi push'lu. **TRUBA senkron** (`$AMR_HOME` = origin/main; `git fetch && git reset --hard origin/main` ile; eski `git checkout origin/main -- <file>` yöntemi ARTIK KULLANILMIYOR).
> Suite **145 passed**, `validate_registry` 0 hata, config'te 0 ölü anahtar. **KB şeması 0.7.1.**

### ⏭️ SIRADAKI İŞ (tek cümle): pilotun ML-yarısı — bir antibiyotik zincirini (A. baumannii amikacin) baştan sona koş, `cv_method='lineage_group_kfold'` + `pipeline_runs`'ta yeni sürüm kolonlarının dolduğunu teyit et. Sonra Faz-2 organizmalarını + tam koşuyu başlat. AYRINTI: en aşağıdaki "PILOT ML-YARISI" bölümü.

### ⚡ KANONİK LINEAGE ETİKETLERİ HENÜZ ÜRETİLMEDİ — İLK İŞ BU OLABİLİR
02c pilot testleri `--out-name poppunk_clusters_*.csv` ile koştu (karşılaştırma için); pipeline'ın okuduğu kanonik `poppunk_clusters.csv` **6 organizma için de yeniden üretilmeli** (yeni container + yeni genomlar + kürasyonlu ayarlar). Mevcut `poppunk_clusters.csv`'ler: A. baumannii + E. coli **eski parametre/refine testlerinden kalma** (kanonik ayarla üretilmedi), K. pneu/Sa/Pa/Efm'de belki hiç yok ya da bayat. **`--out-name` VERMEDEN** (kanonik isme yazsın), `--reuse-db` ile sketch'i tekrar kullanarak, registry'nin çözdüğü ayarlarla (A. baumannii otomatik refine alır) her organizma için 02c koş. SLURM ya da debug node — **login node'da ASLA** (sketching CPU-ulimit'ine takılıp ölür).

### ✅ CONTAINER'LAR: üçü de pinli + kurulu + kilitli + terfi etmiş
Kanonik isimler artık YENİ imajları gösteriyor; eskiler `*.RETIRED_20260715` olarak duruyor (tam koşu bitince silinir; Mac'te `backup/containers_20260715/` + SHA256SUMS var).

| container | ne koşar | kilit | kritik sürümler |
|---|---|---|---|
| `amr.sif` | 00a-13b + **02c PopPUNK** | `environment.lock.yml` (355 pkg) | poppunk 2.7.8 · unitig-caller 1.3.2 · **graph-tool 3.0** · python 3.12.13 · blast 2.17.0 |
| `amr-tools.sif` | **14 pyseer** + **16 M13** + quast | `environment-tools.lock.yml` (206) | pyseer 1.4.1 · **amrfinder 4.2.7** · resfinder 4.5.0 · kma 1.6.11 · **python 3.8.20** · blast 2.16.0 |
| `amr-checkm2.sif` | 02d/M15 QC | `environment-checkm2.lock.yml` (127) | checkm2 1.1.0 · diamond 2.1.11 · python 3.12.13 |

**Ortamlar arası sürüm ayrışması KASITLI ve belgelenmiş** (python 3.8 vs 3.12, numpy 1.22 vs 2.5, pandas 2.0 vs 3.0, blast 2.16 vs 2.17): pyseer/resfinder/AMRFinderPlus python 3.8'e çözülüyor, **tek container bu yüzden imkânsız**. Ortamlar birbirine düz metinle veri veriyor (rtab → pyseer → CSV → populate), güvenli. **Methods'ta "hangi adım hangi ortamda" tablosu yazılmalı.** Def dosyalarının eski gerekçesi ("CheckM2 python<3.9'a çakılı") **tam tersiydi**, düzeltildi.

### 🔬 CONTAINER DOĞRULAMASI YAPILDI — graph-tool kümelemeyi DEĞİŞTİRDİ, bilerek kabul edildi
`amr.sif`'in doğrudan pinleri birebir tuttu ama **transitif bağımlılıklar kaydı**: **graph-tool 2.98 → 3.0** (PopPUNK'ın ağ backend'i), mandrake 1.2.5 → 1.2.4, libprotobuf 6.33.5 → 7.35.1. İş `6098297` E. coli'yi **ESKİ parametrelerle** yeniden kümeledi (tek değişken: container) → `compare_lineage.py`:

```
ARI 0.990368  |  324 → 397 küme  |  singleton 157 → 234
en büyük küme: 817 (%14.9) — İKİSİNDE DE AYNI, top5'in 4'ü aynı
```

Yani fark **saf kuyruk parçalanması**; lineage-CV'nin dayandığı büyük soylar bozulmamış, sızıntı yok, `config.yaml`'ın asıl iddiası ("dbscan bgmm gibi %94'lük mega-kümeye çökmüyor") 3.0'da birebir ayakta. **KARAR (kullanıcı): 3.0 kabul edildi** — zaten her şey yeniden koşuluyor, `environment.lock.yml` 3.0'ı pinliyor. Ölçüm `config.yaml` lineage bloğunda kayıtlı. **Her container rebuild'inden sonra tekrar ölç** (`scripts/compare_lineage.py`).

> **DERS:** build "complete" dedi, 143 test yeşildi, araçlar çalıştı — sürüklenme yalnızca *eski parametrelerle koşup eski etiketlerle karşılaştırdığımız için* göründü. Bunu atlasaydık tüm ESKAPEE panelini farklı bir kümelemeyle koşar, farkı asla bilmezdik.

### ✅ ARAÇ SÜRÜMÜ PROVENANCE'I (şema 0.7.1)
`pipeline_runs` **kmc'yi** (terk edilmiş k-mer baseline'ının QC aracı) kaydediyordu ama **unitig-caller'ı** (özellikleri üretir) ve **PopPUNK'ı** (CV gruplarını tanımlar) kaydetmiyordu — KB kendi soy etiketlerini neyin ürettiğini söyleyemiyordu. Eklendi: `unitig_caller_version` · `bcalm_version` · `poppunk_version` · **`graph_tool_version`** · `blast_version` · `pyseer_version`. graph-tool ayrı kolon çünkü **PopPUNK'ı pinlemek davranışını pinlemiyor** (yukarıdaki ölçüm). pyseer'i **14 kendi summary'sine yazıyor** (amr-tools.sif'te yaşıyor, populate amr.sif'te — sadece aracı çağıran adım dürüstçe raporlayabilir). Migrasyon yolu var: 0.7.1-öncesi KB kolonları kazanır, eski satırlar NULL kalır (dürüst "bilinmiyor").
**M13 provenance düzeltildi:** `AFP_SOURCE` **"AMRFinderPlus 2026-05-15.1"** diyordu — o **veritabanı** sürümü, yazılım **4.2.7**. Artık araca sorulup `AMRFinderPlus 4.2.7 (DB 2026-05-15.1)` yazılıyor; probe exit-code kontrol ediyor (yoksa `"No module named resfinder"` hata mesajını sürüm diye KB'ye damgalıyordu).
**Sürüm drifti kalıcı olarak kapatıldı:** `tests/test_version_alignment.py` — 5 dosyayı `KB_SCHEMA_VERSION`'a çiviler, `.zenodo.json`'ın **prose'unda** gizli şema numarasını yakalar, config'in kopyası geri gelirse kırılır. (config'de **6. bir kopya** vardı: `provenance.kb_schema_version: "0.6.1"` — iki minor bayat ve `collect_versions` onu **KB'ye yazıyordu**; kaldırıldı, kod sabiti tek kaynak.)

### 📥 GENOMLAR İNDİ (yeni filtreyle)
| organizma | indi / hedef | kayıp | not |
|---|---|---|---|
| E. coli | 5681 / 6087 | %6.7 | eski hedef 5470 → **+617** |
| K. pneumoniae | 4791 / 5167 | %7.3 | eski 4615 → **+552** |
| S. aureus | 2639 / 2803 | %5.8 | eski 2494 → **+309** |
| E. faecium | 2078 / 2275 | %8.7 | yeni |
| A. baumannii | 1171 / 1251 | %6.4 | yeni |
| P. aeruginosa | ~1486 | — | ilk koşu **yarım kaldı** (rapor yazılmadan öldü), `--workers 8` ile yeniden koşuldu |

Kayıpların **tamamı** `empty or non-FASTA response` — BV-BRC'de o assembly'ler gerçekten yok (bilinen davranış, %6-9 tutarlı). **Filtre düzeltmesinin gerçek kazancı +%11-12** (ilk yazdığım +%23-108 YANLIŞTI: temizleyici-öncesi sayıyı sonrasıyla kıyaslamışım).

### Bu session ne yaptı — "belge bir şey diyor, kod başka şey yapıyor" sınıfı
12-modül audit'i kodu inceledi ama **kodun kendisi hakkında söylediklerini** incelemedi. Testler bunları göremez (hiçbiri kod hatası değil). ~13 bulgu, üçü doğrudan tezi vuracaktı:
- **`lineage.min_cluster_size: 10`** hiçbir şey tarafından okunmuyordu, `collapse_rare_clusters` hiç çağrılmıyordu (sadece testi vardı; MODULE_03 raporu onu "bilimsel pillar" diye övmüş) → Methods'a **yanlış cümle** yazdıracaktı. **Silindi**, politika belgelendi: PopPUNK kümeleri **olduğu gibi** kullanılır (havuzlama daha kötü: akraba olmayan singleton'ları tek grupta toplayıp tek fold'a atar). E3'ün parent-merge önerisi bize kapalı (iterative-PopPUNK hiyerarşisi gerekir; biz düz dbscan).
- **`card_nt` + 08** → 08 CARD DB yoksa sadece **uyarıp devam ediyordu**; ama tier/gene_symbol/ARO oradan gelir → **içi boş KB, exit 0, kimse fark etmez.** Artık hard-fail (`AMR_ALLOW_MISSING_CARD_DB=true` kaçış).
- **symlink tuzağı** (aşağıda) → deploy ortasında pipeline koptu.

Diğerleri: `shap>=0.44` ölü (TreeSHAP XGBoost built-in, `import shap` yok) · 8 ölü config anahtarı (`encoding`, `optimization_metric`, `booster`, `registry.*_file`, `organism_display`, `blast_db_dir`, `analysis_results_dir` — çoğu **var olmayan yetenek vaat ediyordu**) · `kb_api` "0.4.0" + "for E. coli" · `08` docstring'inde "Why Nextflow?" · `resolve_path` docstring'i silinmiş `paths:` bloğuna atıf · `config.project.version` 0.6.0 ama yorumu "KB şemasıyla hizalı" diyor.

### Yapılanlar (hepsi push'lu)
- **`evidence_tier` (schema 0.7.0)** — `unitig_evidence_tier` tablosu + `classify_evidence_tier` (confirmed / **strong_novel** / candidate / weak / none). CPSS+pyseer = novelty backbone. `/api/v1/novel` ucu. 7 test. *Geçen session yapılmış ama commit edilmemişti.*
- **E2/E3 okundu** (kararlar aşağıda), `docs/literature/` commit'lendi.
- **`environment.yml` TAM PİNLENDİ** + `environment.lock.yml` (355 paket, build-string'li) commit'lendi. `amr.def` tek-container tasarımını belgeliyor (`amr.sif`/`amr-pp.sif` ikiliği **tasarım değil, build-tarihi kazası**: amr-pp = aynı reçete 3 saat sonra, poppunk eklendikten sonra; ortak paketlerin build string'leri birebir aynıydı).
- **02c PopPUNK parametreleri açıldı** — E3 ayarları panel geneline (k 15-35 step 2), S. aureus `sketch_size: 10^5` registry override'ı, `--qc-db` (daha önce **hiç** koşmuyordu). CLI flag'leri (`--min-k/--max-k/--k-step/--sketch-size/--no-qc/--out-name`) + `scripts/compare_lineage.py`.
- **HPC kaynakları env'e** — `AMR_KMC_MEM` / `AMR_THREADS` (`load_config`'te, tek noktada). Artık `git reset` HPC ayarlarını silmiyor; HANDOFF'un "reset sonrası elle tekrar uygula" notu **geçersiz**.
- **Zenodo/CITATION/pyproject/kb_app** → hepsi **0.7.1** + ESKAPEE (artık `test_version_alignment.py` çiviliyor). `.zenodo.json` "E. coli / schema 0.4.0" diyordu — **DOI ile kalıcılaşacaktı.** `notes` artık sayı tekrarlamıyor, "amrk.db'yi oku" diyor.
- **TRUBA temizliği: 766 GB → 216 GB.** `blast_db/core_nt` **318G** (yarıda bırakılmış NCBI indirmesi, kod referansı sıfır) · `data/interim` 155G · `ecoli/ampicillin/matrix` 49G (eski ham k-mer baseline) · staph rtab 3.6G · smoke artıkları. **KORUNDU:** genomlar 58G, card/amrfinder/resfinder/checkm2 DB'leri, containers.
- **Yedekler:** 4 sif + SHA256SUMS → `backup/containers_20260715/` (Mac). `$AMR_WORK/backup/predeploy_20260715/` → KB + figures + tables + runs + **eski `poppunk_clusters.csv`** (= container testinin kontrol grubu).

### ⚠️ SYMLINK TUZAĞI (yaşandı, çözüldü — bir daha kurma)
`$AMR_HOME/{data,results,runs,models,logs}` **symlink** → `/arf/scratch`. Git izlenen bir yolu symlink'in içine **yazmaz**: symlink'i silip yerine gerçek dizin koyar. `data/external/blast_db/card_nt/*` (8 dosya) + `runs/.gitkeep` izlendiği için **`git reset --hard` `data` ve `runs` symlink'lerini kopardı** — 156 GB görünmez oldu (veri kaybı yok; symlink hedefi silinmez). Kurtarma: `rm -rf data runs && ln -s $AMR_WORK/data data && ln -s $AMR_WORK/runs runs`.
**Kalıcı çözüm:** 9 dosya untrack + `.gitignore`'a `data/external/blast_db/` + `!runs/.gitkeep` kaldırıldı + **`tests/test_repo_symlink_safety.py`** (symlink'li önek altına dosya eklenirse kırmızı). `card_nt` git'te tutulmasının tek gerekçesi ("08 kutudan çıkar çıkmaz çalışsın") 08'in hard-fail'iyle karşılandı.

### E2 KARARI (H3 istatistiksel testi) — `docs/literature/E2.md`
**Ham unitig üzerinde Fisher/hypergeometric = hakem reddi** (LD bağımsızlığı yıkar; tek plazmit binlerce korele unitig sokar → sahte mikroskobik p). Doğru çerçeve: **(1)** önce bileşenleştir (cDBG alt-grafı ya da Pearson>0.95 kümeleme) → **(2)** null evren = iki modelin ön-filtre sonrası girdi uzaylarının **KESİŞİMİ** (bizim "union-universe" caveat'ımızın cevabı; pan-genom **değil**, bilinen AMR genleri **değil**) → **(3)** Fisher exact (bileşenlerde) → **(4)** 500-1000 etiket permütasyonu ile temsili çiftlerde doğrula → **(5)** **Benjamini-Yekutieli** (BH değil — çapraz-direnç negatif korelasyon yaratır, PRDS ihlal) → **(6)** Overlap Coefficient + Fold Enrichment (Jaccard boyut-dengesizliğinde yanıltır; setlerimiz 36 vs 83) → **(7)** UpSet + kümelenmiş heatmap. **İyi haber:** `15_cross_antibiotic.py`'nin ARO gen-ailesi katmanı zaten bileşenleştirme → orada Fisher bugün savunulabilir. Minor-revision kalemi.

### E3 KARARI (PopPUNK) — `docs/literature/E3.md`
Literatür **BGMM K=2 + refine** diyor, biz **dbscan + refine yok** kullanıyoruz → **bilinçli sapma, ampirik gerekçeli** (`config.yaml:167-172`: bgmm E. coli'de ~%94'lük mega-kümeye çöktü, refine dejenere NaN sınırda öldü) **ve E3'ün önerisi iterative-PopPUNK varsayıyor, biz düz fit koştuk.** Methods'ta aynen böyle yazılmalı. E3'ün HDBSCAN-yakınsama uyarısı yüksek-rekombinasyonlu türleri hedefliyor → **Enterobacter (Faz 2) patlayabilir**; loud fail var (`n_clusters >= n_splits`, "try --model bgmm"). ✅ StratifiedGroupKFold zaten altın standart (`lib/lineage.py:122`). **Hâlâ eksik ve hakem-beklentisi: random-vs-lineage-aware CV karşılaştırma tablosu** (ucuz, yüksek değer).

### VERİ KAYNAĞI + PANEL — 2026-07-15 ikinci yarı (E4 sonrası, hepsi ÖLÇÜMLE)

**Panel = 6 organizma** (kodda `registry.is_active` ile çözülüyor, iddia değil): E. coli · K. pneumoniae · S. aureus · A. baumannii · P. aeruginosa · E. faecium. **Enterobacter DIŞLANDI** (`status: excluded_insufficient_data`) — tam ECC kompleksi (6 tür, 466 genom) tarandı, hiçbir antibiyotik minority ≥150'yi geçmiyor (en iyi: gentamicin 89). E4 "taksonomi sorunu, kompleksi aç, kurtulur" dedi → **test edildi, çürüdü**: BV-BRC gönderenin etiketini koruyor, 550 büyük kova (11935 satır), hormaechei küçük (1283) — E4'ün varsaydığının tersi. Gerekçe registry bloğunda tam ölçümle yazılı. Diğer 5 organizmanın taxid'leri de kompleks-kontrolünden geçti: **temiz** (K. pneu 85291 vs variicola 1444; A. baumannii 28237 vs pittii 265 — ayrı türler, birleştirilmemeli).

**🔴 FİLTRE HATASI BULUNDU VE DÜZELTİLDİ.** Eski sorgu `evidence="Laboratory Method"` İSTİYORDU. Ama BV-BRC gerçek CLSI/EUCAST ölçümlerinin çoğunda `evidence` alanını **boş** bırakıyor → o satırlar atılıyordu. Doğru kural: **"computational olanı at"**, "lab olanı iste" değil. **Gerçek kazanç +%11-12** (E. coli 5470→6087, K. pneu 4615→5167, S. aureus 2494→2803, P. aeruginosa 1312→1486 + 6 antibiyotik). ⚠️ Bu satır önce "+%108 / +%55 / +%62 / +%23" diyordu — **yanlıştı**: temizleyici-ÖNCESİ satır sayılarını temizleyici-SONRASI genom sayılarıyla kıyaslamıştım; `testing_standard` (EUCAST/CLSI) filtresi devreye girince rakam düşüyor. İKİ katmanda birden düzeltildi (`00a` sorgusu + `lib/bvbrc` temizleyicisi — cleaner'ın kendi `contains("laborator")` kontrolü vardı, düzeltilmeseydi 00a düzeltmesi no-op olurdu). Ayrıca "computational" **İKİ kolonda** saklanıyor: `evidence` + `laboratory_typing_method` (K. pneu'da 9814 satır). İkincisi eski filtreden sızıyordu — etiketlere ULAŞMIYOR (hepsinin `testing_standard`'ı boş, standart filtresi eliyor, ölçüldü: 0 geçiyor) ama artık **açıkça** kapatıldı, şansa bırakılmadı.

**Computational satırlar ASLA etikete girmez** — E. coli'de 6.9M(!) vs 243K lab. Girerse model AMRFinder'ı taklit eder ve M13 flagship'i (0.926 vs 0.538) kendi kendini kıyaslamaya döner.

**MIC YENİDEN-ETİKETLEME: YAPILMIYOR.** E4 §4 "ham MIC'i tek standarda çevir" diyor, ama ölçtük: mg/L doluluk **S. aureus'ta %9** (A.b %76, Kp %66, Pa %58, Efm %44, Ec %28). S. aureus'un %91'ini kaybetmek = gram-pozitif/çapraz-filum iddiasını kaybetmek. E4 kendi §5'inde zaten "BV-BRC hazır etiketleri endüstri standardı, savunulabilir" diyor. **Karar:** BV-BRC etiketleri kullanılır, heterojenlik Methods'ta limitation olarak yazılır, ve **MIC'i olan alt kümede duyarlılık analizi** koşulur (ör. "Kp etiketlerinin %66'sı ham MIC taşıyor; EUCAST 2024'e göre yeniden yorumlayınca AUC 0.94→0.93" = heterojenliğin sonucu etkilemediğinin KANITI). Ham MIC kolonları (`measurement`/`sign`/`value`/`unit`) SELECT_FIELDS'e eklendi — şimdi kullanılmıyor, sonra yeniden indirmemek için. NOT: ölçüm birimleri karışık (E. coli'de %69 birimsiz, %4 'mm'=zon çapı) ve `laboratory_typing_method`'da çöp var (`'2014,2015'`, `Vitek_2-P607_card`) → duyarlılık analizi ciddi veri temizliği gerektirir.

**⛔ E1 TERSİNE DÖNDÜ — TEMPORAL VALIDATION İMKÂNSIZ.** Detay + rakamlar: `FINAL_AUDIT §E1`. Özet: etiketlenebilir genomlarda **≥2023 → E. coli 28, K. pneu 13, S. aureus 0, A. baumannii 11.** BV-BRC'nin AMR verisi 2021'de bitiyor. E1 doluluğa bakıp "fizibıl" demiş, **dağılıma bakmamış**. **Yıl (collection_year + testing_standard_year) tamamen işlemden çıktı.** Yeni dış-doğrulama stratejisi = **coğrafi + soy hold-out**; `isolation_country` %96-99 dolu. Uyarılar: ülke dominansı (E. coli %58 Norveç, A.b %63 USA → hold-out olamaz; dengeli olanlar S. aureus UK/Çin/USA ve K. pneu), ve **ülke-soy confounding** (coğrafinin PopPUNK-CV üstüne ne kattığı ölçülmeli). **Veri yetmezse coğrafiyi de katma** (kullanıcı kararı). M13 concordance yerinde duruyor.

### ✅ 2026-07-16 BİTENLER (bu maddeler ARTIK YAPILDI — yeniden yapma)
- **P. aeruginosa indirmesi:** 1382/1486 (%7 kayıp, normal — ilk koşu login-node/rapor sorunuyla 966'da ölmüştü, `--workers 8` ile toparlandı).
- **Faz B temizlik YAPILDI:** disk **230→102 GB** (128 GB açıldı). `data/processed/*` (matrisler + lineage), `results/{ecoli,kpneumoniae,staph,kb,figures,tables}`, `models/*`, `runs/*` içerikleri silindi (DİZİNLER değil — symlink hedefleri). KORUNDU: `data/raw` genomlar 79G, `data/external` DB'leri 3.3G, `containers` 6.7G, `backup` 6.2G.
- **Panel KÜRE EDİLDİ — 45 model, 14/14 sınıf** (`organisms.yaml` commit `1d11cc3`): Ec 8 · Kp 11 · Sa 8 · Ab 8 · Pa 4 · Efm 6. Ham minority≥150 listesi 73 çiftti; sınıf-showcase için kırpıldı (Ec/Kp'de 6 sefalosporin → 2-3). Tüm çapraz-organizma ilaçlar (cipro/tetra 5x, gent/tmp-smx 4x, mero/ceftaz 3x) + flagshipler korundu (carbapenem mero+imi Kp/Ab, MRSA cefoxitin+oxacillin, VRE vanco+teico). Dışlanan: "extended spectrum beta lactamase" (fenotip etiketi, ilaç değil), linezolid Efm (minority 67). **BUG düzeltildi:** `ampicillin/sulbactam` slash-kanonik kayıtlıydı (yol kırıcı, M9 taraması atlamış) → underscore-kanonik yapıldı (A. baumannii'nin kilit ilacı, minority 405).
- **`--qc-db` DOĞRULANDI:** A. baumannii pilotunda hatasız koştu, genom eliyor (1171→1133/1138). Kör yazılmıştı, artık kanıtlı.
- **⭐ LINEAGE model/refine per-organism AYARLANDI** (commit `b9c0720`) — aşağıdaki ayrı bölüme bak. **A. baumannii `refine: true` registry override aldı**, gerisi global dbscan/no-refine. Registry-override-gölgeleme bug'ı düzeltildi.

### DEPLOY — KALAN ADIMLAR (sıralı — GÜNCEL)
1. **⏭️ KANONİK lineage etiketlerini üret** — 6 organizma için `poppunk_clusters.csv` (yukarıdaki "KANONİK LINEAGE" kutusu). `--out-name` VERME, `--reuse-db` kullan (sketch DB'ler `data/processed/{org}/lineage/_poppunk_work/db`'de duruyor → dakikalar). A. baumannii otomatik refine alır. **Login node'da ASLA.**
2. **⏭️ PILOT ML-YARISI** — bir antibiyotik zinciri (A. baumannii **amikacin**, minority 467) baştan sona: `03u → 04-07b → (07-09) → 10-13b → 14 → populate`. DOĞRULA: (a) `07b` logunda `cv_method='lineage_group_kfold_5fold'` (fallback DEĞİL), (b) `pipeline_runs`'ta yeni sürüm kolonları dolu (`poppunk_version`, `unitig_caller_version`, `graph_tool_version`, `pyseer_version` vb.), (c) `unitig_evidence_tier` tablosu doluyor + `strong_novel` adaylar var mı. Detay: EN ALT "PILOT ML-YARISI".
3. **Tam ESKAPEE koşusu.** ⚡ **`unitig_all` store hiç kurulmamış** → **`03u --build-db` ile organizma başına BİR kez kur**, sonra her antibiyotik subset eder. 45 model, per-antibiyotik zincir env-parametrik SLURM (memory `amr-truba-gotchas`: zincir KESİN SIRALI, pyseer bio'yu okur, PARALEL DEĞİL).
4. **A. baumannii `length_range`** (opsiyonel rafinman) — genomlar indi, gözlenen uzunluk dağılımından türet (E3 §5); global `length_sigma:5` şu an adaptif filtreliyor, blocker değil.
5. `slurm/` **kanonik 4'ü** repoya commit — `run_lineage.slurm` hâlâ emekli `amr-pp.sif`'i çağırıyor, `amr.sif`'e çevir.
6. Zenodo deposit + DOI.

**Analiz kalemleri (modeller ÇIKTIKTAN sonra, koşuyu bloke etmez):** MIC duyarlılık analizi (§ yukarı) · coğrafi hold-out + soy-CV'ye ne kattığının ölçümü · **random-vs-lineage-aware CV karşılaştırma tablosu** (E3 §8: hakem açıkça bekliyor, hâlâ yok) · H3 hypergeometric (E2 çerçevesi).

**Açık, küçük:** `amr-gpu.def` ölü (M9-4) · 4 def'te de base imaj `condaforge/miniforge3:latest` (digest-pin TODO) · `*.RETIRED_20260715` sif'ler (3 tane: amr/amr-pp/amr-tools/amr-checkm2'nin eskileri; tam koşu bitince sil) · biriken `screen` oturumları · `$AMR_WORK/`'te pilot log'ları + `minority.py`/`geo.py`/`leak.py`/`measure.py` (silinebilir).

### ⭐ LINEAGE model/refine — PER-ORGANIZMA, ÖLÇÜMLE (2026-07-16, commit `b9c0720`)
Pilot A. baumannii'de dbscan'in **%76 mega-küme** verdiğini yakaladı (GC2 klonu). Kullanıcı doğru sordu: "değiştirirsek hepsinde bakmak lazım değil mi?" → tüm panel ölçüldü, **en büyük soy oranı** (CV-denge kriteri), dbscan vs dbscan+refine:

| organizma | no-refine | refine | SEÇİLEN |
|---|---|---|---|
| E. coli | 15.4% | 15.4% | no-refine (berabere) |
| K. pneumoniae | **22.3%** | 58.6% | no-refine — **refine BOZUYOR** (ST258 klonlarını birleştiriyor) |
| S. aureus | 20.6% | 20.6% | no-refine (berabere) |
| P. aeruginosa | 12.4% | 19.3% | no-refine |
| E. faecium | 11.4% | 4.2% | no-refine (refine aşırı böler: 872 soy/668 singleton) |
| **A. baumannii** | 76.3% | **52.1%** | **REFINE** — tek override |

**KARAR: ayar değil, KURAL tek tip** — "dbscan varsayılan; refine yalnızca varsayılan tek soyu domine bırakırsa (>%40)". Nesnel kriter, hepsine aynı uygulandı, sadece A. baumannii'de farklı çıktı. Hakem-savunulabilir (gerekçesiz per-organism = cherry-picking). Global default (dbscan/no-refine) 6'nın 5'i için doğru; A. baumannii `organisms.yaml`'da `lineage: refine: true`. **Kp %22 + Ab %52 en dengesizler ama BİYOLOJİ** (yüksek-riskli klonlar CG258/ST258, GC2) — Methods'ta "klonal-domine, lineage-CV en muhafazakâr" diye yaz. sketch 10⁵ A. baumannii'de yardım etmedi (76.3→76.4), bgmm de (76.3, dbscan'le aynı). **Sketch/model/refine bir organizma için değişirse: tümünde yeniden ölç** (`compare_lineage.py` + en-büyük-soy).

### OPERASYONEL (bugün öğrenilenler — memory `amr-truba-gotchas`)
- **Container build VE PopPUNK sketching giriş düğümünde ÖLÜR** — CPU-time ulimit'i uzun/ağır işlemi keser. Build'de `mksquashfs`, 02c'de `poppunk --create-db` çakılıyor ("Command failed", sebep gizli). **`srun -p debug -N1 -c8 --time=02:00:00 --pty bash`** ile compute node'a geç (2026-07-16: 4 organizmanın 02c'si login node'da hep `--create-db`'de öldü; debug'da sorunsuz). İlk satırda `hostname` ile teyit et — `arf-ui1` ÇIKMAMALI.
- **Yeni `srun --pty` oturumu ENV'i SIFIRLAR** — `AMR_HOME`/`AMR_WORK`/`APPTAINER_BINDPATH` gitmiş olur; her yeni debug oturumunda yeniden `export` et yoksa komutlar sessizce düşer (2026-07-16 bunu bir kez yaşadık).
- **SLURM barbun: `--nodes=1` ŞART** — yoksa "node başına 20 çekirdek" QOS kontrolü hesaplayamaz ve reddeder. **`export APPTAINER_BINDPATH=/arf` ŞART** — yoksa container `/arf`'ı görmez. Altın örnek: `$AMR_HOME/slurm/run_lineage.slurm` (PopPUNK'ı daha önce başarıyla koşan script).
- **Container'a `conda list` ile sorma** — prefix vermezsen miniforge **base** env'ini listeler ve "paket yok" yanılgısı verir. Prefix'ler def dosyalarında yazılı: `amr.sif`→`/opt/amr-env`, `amr-tools.sif`→`/opt/amr-tools-env`, `amr-checkm2.sif`→`/opt/amr-checkm2-env`. `which` bu container'larda "Illegal option" veriyor — PATH'i `<tool> --version` ile yokla ya da prefix'i doğrudan ver. Bu yüzden iki kez yanlış teşhis kondu.
- **`00a` indirmesi sessizce ölebilir** — P. aeruginosa 966/1486'da `download_report.csv` yazmadan öldü ve döngü (`set -e` yok) sessizce sonrakine geçti. **Rapor dosyasının varlığı = koşunun bittiğinin kanıtı**; sadece `.fna` saymak yanıltır. Ölürse `--workers` düşürüp yeniden koş (resume-safe, mevcutları atlar).
- TRUBA'ya `rsync`/`ssh` **IP ile**: `172.16.6.14` (`.11`, `.16` de var); `arf-ui1` küme-içi ad, dışarıdan çözülmez. `.sif` çekerken `-z` **kullanma** (zaten sıkıştırılmış).
- Biriken `screen` oturumları (`screen -ls`) — eski koşulardan ölü kabuklar, süpürülebilir.

### 🎯 PILOT ML-YARISI — sıradaki oturumun ilk ML işi (ayrıntı)
02c/`--qc-db`/lineage tarafı bitti; şimdi bir antibiyotik zincirini uçtan uca koşup ML+KB tarafını doğrula. **Hedef: A. baumannii amikacin** (en küçük organizma 1171 genom + güçlü split minority 467 + refine override'ını test eder). Zincir env-parametrik SLURM (memory `amr-truba-gotchas`'taki reçete): `sbatch --chdir=$AMR_WORK --export=ALL,AMR_ORGANISM=acinetobacter_baumannii,AMR_ANTIBIOTIC=amikacin,AMR_FEATURE_REPR=unitig $AMR_HOME/slurm/<script>`. Sıra KESİN: Faz1(03u→04-07b) → Faz2a(07→08→09, screen, 08 NCBI ister) → bio(10-13b) → pyseer(14) → populate → rm rtab. **DOĞRULANACAK 3 ŞEY:**
  1. **`07b` lineage-CV gerçek mi:** logda `cv_method='lineage_group_kfold_5fold'` (fallback `repeated_holdout_5seed` DEĞİL). Kanonik `poppunk_clusters.csv` üretilmiş olmalı (adım 1) yoksa fallback'e düşer.
  2. **Yeni provenance kolonları dolu mu:** populate sonrası `pipeline_runs`'ta `poppunk_version`, `unitig_caller_version`, `bcalm_version`, `graph_tool_version`, `blast_version`, `pyseer_version` NULL DEĞİL. (pyseer'i 14 kendi summary'sine yazıyor → populate oradan okuyor; container sqlite3 CLI yok → `python -c "import sqlite3..."`.)
  3. **`unitig_evidence_tier` (0.7.1) doluyor mu:** `strong_novel` adaylar (CPSS+pyseer geçen, CARD geni olmayan) görünüyor mu — evidence_tier feature'ının asıl amacı.
- populate NOT: **idempotent değil** model child-row'ları varsa (memory'de detay) — zincir pyseer'e KADAR koşulup öyle populate edilmeli, TEK kez. Env inline şart (`AMR_ORGANISM=... AMR_ANTIBIOTIC=... AMR_FEATURE_REPR=unitig apptainer exec ...`) yoksa config default'a (ecoli/gentamicin) düşer.

---

# §0.-5 — (2026-07-14) — 12-MODÜL PRODUCTION-HARDENING İNCELEMESİ + ESKAPEE PİVOTU — **§0.-6 TARAFINDAN DEVRALINDI** (tarihsel)

> ⚠️ Bu bölümün "DEPLOY ÖNCESİ YAPILACAKLAR" listesi **BİTTİ** (evidence_tier ✅, E2/E3 ✅ okundu, temizlik ✅, container ✅ kuruldu) ve şema artık **0.7.0**. Aşağıdaki "yapılacak"ları **uygulama** — güncel durum §0.-6'da. Bölüm, kararların gerekçesi için duruyor.

> **Repo HEAD `3d92bb6`** (pushed to `github.com/demirbase/ML_AMR_Prediction_v2`). Bu session = **sistematik 12-modül audit + sertleştirme** (yeni özellik değil; tez+makale+production+TRUBA hazırlığı). Suite artık **TAMAMEN YEŞİL (109 passed / 0 kırık)** — HANDOFF §0.-4'teki "117 pass" bayattı, suite kırmızıydı; düzelttik.

### Proje yeni yönü (danışman toplantısı 2026-07-14)
- **2 makale** tezden çıkacak. **Araç/web-app en son.** Odak = **organizma değil, antibiyotik SINIFI.** Organizmalar = **ESKAPEE.** Pipeline production-ready + Nur'un pipeline'ına entegre (ayrı makale).
- **ESKAPEE paneli (E1/ESKAPEE1.md literatürüyle):** **Faz 1** (tez+1. makale): K. pneumoniae, E. coli, S. aureus, A. baumannii. **Faz 2**: P. aeruginosa, E. faecium, Enterobacter. (Veri hep BV-BRC'den; uygunluk kapı değil.)

### 12-modül incelemede yapılanlar (hepsi push'lu, her modülün raporu `docs/review/MODULE_00..11_*.md` + `FINAL_AUDIT_pre_deploy.md`)
- **M0 Repo:** `amr.def` TRUBA'dan kurtarılıp geri yüklendi (container reçetesi environment.yml'den türer). `backup/` gitignore. 9 commit push.
- **M1 Config/registry:** **registry schema 2.0** — mekanizma-bazlı 19 sınıf (**carbapenems ayrı**, `others` dağıtıldı: glycopeptides/polymyxins/oxazolidinones/phenicols/rifamycins/fosfomycins/lipopeptides/nitrofurans/glycylcyclines). `organisms.yaml`: gerçekle senkron (ecoli 7, kpneu 14 `status:done`) + 7 ESKAPEE organizması + gram/phylum/eskapee_phase/priority_classes. **mechanism_type/who_aware registry'ye taşındı** (populate artık registry'den okuyor — hardcode yok). `feature_repr: unitig` **varsayılan** (artık env-override ŞART değil). `scripts/validate_registry.py` bekçisi (registry↔KB, CI'da).
- **M2 Veri:** **BUG FIX** `lib/bvbrc._resolve_group` `np.argmax→nanargmax` (kısmi-NaN yılda yanlış R/S etiketi). Bilinmeyen-antibiyotik raporu + `intermediate_policy` config switch (drop varsayılan). API-truncation loud.
- **M3 QC/lineage:** 02/02b/02p/03/03b → `get_target` (env-parametrik, ESKAPEE paralel). legacy `paths:` bloğu silindi. 02c PopPUNK docstring uyumlandı. `lib/lineage.py` mükemmel.
- **M4 Unitig:** **S. AUREUS BLOCKER ÇÖZÜLDÜ** — `03u.select_genomes`'a seed-42 shuffle (fenotip-bloklu sıra→tek-sınıf-chunk→XGBoost NaN). deprecated `np.fromstring`→`frombuffer`.
- **M5 Modelleme:** **`07b` fallback şeffaflığı** — `cv_method` summary'ye yazılıyor + non-lineage-CV loud uyarı + `no_group_leakage` assert.
- **M6 Biyoloji:** `test_15` düzeltildi. **`08_blast_pipeline.nf` öksüzdü → M9'da silindi.**
- **M7 KB:** **`models.cv_method` kolonu** (schema **0.6.1**) — her modelin lineage-CV mi fallback mı olduğu KB'de; Mac KB backfill (21=lineage_group_kfold_5fold). `test_kb_queries` fixture düzeltildi → **tüm suite yeşil.**
- **M8 Arayüz:** **BUG FIX** `kb_api` DB path `results/ecoli/kb`→`results/kb/amrk.db`. `get_overlap`+`/overlap` organism-aware.
- **M9 Orkestrasyon:** **Nextflow tamamen kaldırıldı** (08.nf + environment.yml dep + docstring + README badge) — "Python'da kal" kararı.
- **M10 Test/CI:** CI'a `validate_registry` (blocking) + `mypy` (advisory).
- **M11 Docs:** CITATION/zenodo/pyproject → **version 0.6.1** + "Escherichia coli"→"ESKAPEE (E. coli, K. pneumoniae)".

### KB durumu (Mac `results/kb/amrk.db`, schema 0.6.1)
21 model, 2 organizma, cv_method='lineage_group_kfold_5fold' (hepsi). **TRUBA/Drive KB kopyaları 0.6.0 (bayat) → re-populate ile senkron olacak.** Mac KB drug_class re-sync'li (meropenem/imipenem→carbapenems).

### DEPLOY ÖNCESİ YAPILACAKLAR (sıralı, deployment planı `FINAL_AUDIT_pre_deploy.md §G`)
0. **⭐ Tier sistemini evidence-weighted yap (kullanıcı isteği, deploy-öncesi):** tier şu an SADECE BLAST (identity/coverage/evalue, `09`). 7 katmanı da (prevalans/SNP/MDA/label-perm/CPSS/pyseer) katan bileşik güven notu → özellikle `tier=none` ama CPSS+pyseer geçen **novel adaylar** görünür olur. BLAST-tier'ı yanında ayrı `evidence_tier` olarak ekle (M1 raporu §5 + memory `amr-future-ideas`). Re-populate öncesi yapılırsa yeni KB doğrudan içerir.
1. **TRUBA temizliği** (konservatif — `FINAL_AUDIT §F`): SAFE-DELETE = `data/interim/*` (155G KMC), staph cefoxitin `unitigs.rtab` (3.6G bayat), smoke/log/cache. KEEP = genomlar/KB/figures/models/runs/poppunk_clusters/unitig_all/containers/slurm.
2. **Bayat çıktı temizliği:** kod değişti (shuffle+registry2.0) → re-populate edilecek her (org,ab) için `data/processed/{o}/{a}/matrix_unitig`, `models/{o}/{a}`, `results/{o}/{a}` sil (resume-logic bayat sonuç almasın). `unitig_all` store KORU.
3. **Backup** (rclone, alt-dizin, arf-ui4/screen) → NEEDS-BACKUP seti.
4. **TRUBA repo reset:** `git -C $AMR_HOME fetch origin && git reset --hard origin/main` (untracked amr.def/slurm etkilenmez) + HPC config (kmc_mem=128/threads=20) tekrar uygula + **`slurm/` kanonik 4'ü (`run_03u_env`,`run_ml_env`,`run_bio_env`,`run_pyseer_env`) repoya commit.**
5. **Container lock** (M9-3, publication-must): `apptainer build amr.sif amr.def` + `conda env export --no-builds > environment.lock.yml` commit + base imaj digest-pin.
6. **Test:** pytest (109 passed) + validate_registry (0 hata) + pilot dry-run (1 ab, --max-genomes 50, cv_method='lineage_group_kfold' teyit).
7. **Tam koşu:** ecoli/kpneu re-populate (0.6.1) → S. aureus (blocker çözülü) → A. baumannii (yeni). Env-parametrik zincir.
8. **Zenodo deposit (M10)** + DOI damgası.

### BİLİMSEL GAP'LER (submission-öncesi, literatür kararlı)
- ~~**E1 — Temporal external validation YAPILACAK**~~ → **2026-07-15 ÖLÇÜMLE ÇÜRÜDÜ, YAPILMIYOR.** BV-BRC'de ≥2023 AMR-etiketli genom yok denecek kadar az (S. aureus 0, K. pneu 13, E. coli 28). Yerine **coğrafi + soy hold-out**. §0.-6 ve `FINAL_AUDIT §E1`'e bak — aşağıdaki eski plan geçersiz.
- **E2 + E3: ✅ OKUNDU (2026-07-15), kararları §0.-6'da.** (Bu satır "HENÜZ OKUNMADI" diyordu.)

### Reviewer riski (FINAL_AUDIT §D): reject YOK. Major-revision riski = temporal validation (E1, çözülüyor) + container-lock. Minor = H3 test (E2), cefotaxime(Kp) head-to-head 0.495 açıklaması, provenance tool sürümleri.

> **Repo HEAD `b413f84`** (pushed to `github.com/demirbase/ML_AMR_Prediction_v2`). Advisor meeting prep session — KB brought to a polished, presentable "final" state for a Friday meeting.

### What this session did (all pushed unless noted)
- **21 models COMPLETE** (E. coli 7 + K. pneu 14). Finished the 4 in-flight K. pneu antibiotics: **cefotaxime m18, levofloxacin m19, ceftazidime m20, amikacin m21** (each full chain 03u→…→14→populate). All 7 validation layers present in every model.
- **KB schema 0.4.0 → 0.6.0** (`scripts/lib/kb_schema.py`): added tables **`organisms`** (gram/phylum), **`external_concordance`** (M13), + columns `antibiotics.mechanism_type`/`who_aware` (WHO AWaRe), `models.n_features`; **`unitig_antibiotic_overlap` made organism-aware** (0.6.0, organism in PK — unified KB must not merge same drug across organisms).
- **`scripts/migrate_kb_050.py`** — additive backfill (does NOT re-insert models → avoids the FK issue). Run with `AMR_FEATURE_REPR=unitig AMR_CARD_VERSION=4.0.1`. **This is how the KB is updated now, not full re-populate.**
- **15_cross_antibiotic** refactored organism-aware (per-organism stable_sets/overlap; `--db results/kb/amrk.db`). Ran for ecoli+kpneu → `unitig_antibiotic_overlap` filled (ecoli 18 + kpneu 13).
- **16_external_concordance** run for BOTH organisms → `external_concordance` = 48 rows (model + AMRFinderPlus + ResFinder, leakage-free head-to-head on held-out test genomes). **K. pneu afp/resfinder produced via `$AMR_WORK/run_16_kpneu.slurm`** (40-core parallel, DB paths from `data/external/{amrfinder_db,resfinder_db,pointfinder_db}`; amrfinder `--organism Klebsiella_pneumoniae`, resfinder `--db_path_res/--db_path_point`). Flagship result: **K. pneu ciprofloxacin model bACC 0.926 vs AMRFinderPlus 0.538 / ResFinder 0.540** (unitig model catches gyrA/parC SNP that gene-based tools miss).
- **9 thesis figures** (`scripts/kb_figures.py`): 00 overview, 01 performance, 02 cpss/pfer, 03 cross-org, 04 mechanism-heatmap, 05 significance (real-vs-null), 06 evidence-layers (7-layer backbone), 07 external (model-vs-tools) + `why_unitigs_not_kmers` SVG concept diagram (visualize tool). All in `results/figures/` (Mac + TRUBA).
- **`scripts/kb_app.py`** fixed: default `results/kb/amrk.db` (was old per-org), tab5 reads `external_concordance` table, tab3 organisms panel, **+ new "Ham tablolar" tab (all 13 raw tables)**. ⚠️ **kb_app.py "Ham tablolar" tab NOT yet committed/pushed** (local only).
- **Registry curation** (commit cfd9dc5): oxacillin (penicillins), macrolides+lincosamides classes (7→9), AFP keywords for amikacin/levofloxacin/oxacillin/erythromycin/clindamycin.
- **2 docs NOT yet committed** (local Mac only): `docs/KB_ACIKLAMA.md` (every table/column explained) + `docs/KB_KAVRAMLAR.md` (concepts + real KB examples, for the advisor).
- **Backups on Drive** (`gdrive:TRUBA_25626/scratch_amr/backup/`): `amr_results_20260709_*.tar.gz` (genome_qc + rtab + amrfinder/resfinder raw excluded → ~40MB meaningful backup).

### KB final state (verified, `results/kb/amrk.db`, schema 0.6.0)
21 models · 2 organisms · 7 drug classes · 2363 unitigs · 13 tables all populated · external_concordance 48 · overlap organism-aware · card 4.0.1 · every model 7 evidence layers. **Known caveat:** cefotaxime(Kp) external head-to-head model bACC≈0.495 (small n=200 test slice; lineage-CV is 0.77) — have the explanation ready for the advisor.

### NEXT SESSION — resume order
1. **Commit + push the uncommitted local work** (needs a fresh fine-grained PAT, no Co-Authored-By): `scripts/kb_app.py` (Ham tablolar tab) + `docs/KB_ACIKLAMA.md` + `docs/KB_KAVRAMLAR.md`. Then optionally `git checkout origin/main -- scripts/kb_app.py` on TRUBA.
2. **S. aureus (3rd organism) — STILL PENDING** (paused before the meeting). The blocker + fix are known: 03u writes genomes in phenotype-blocked order (clonal MRSA) → 04's chunk-split gives single-class folds → nan. **Fix ready but NOT applied:** add a deterministic seed-42 shuffle of `valid_genomes` in `scripts/03u_unitig_matrix.py:select_genomes` (before return), then re-run cefoxitin 03u (unitig_all store exists → fast) → Faz1 → full chain. Then the other 7 staph antibiotics (registry already curated: cefoxitin/oxacillin/ciprofloxacin/gentamicin/tetracycline/tmp-smx/erythromycin/clindamycin). S. aureus completes the gram-positive cross-phylum claim. (Also run 02c PopPUNK for staph so 07b lineage-CV is real, not fallback.)
3. **M10 Zenodo deposit** — last must-have; deposit the 21-model 0.6.0 KB, stamp `kb_metadata.zenodo_doi`.
4. **Optional polish:** external_concordance currently head-to-head (test-set); the broader csv-based afp/resfinder-vs-phenotype (all genomes, more antibiotics) is in `16_concordance_{org}.csv` if a wider table is wanted. Docs (METHODOLOGY/ROADMAP) fold-in of the 0.6.0 schema + 21-model panel.

### Operational reminders (also in memory `amr-truba-gotchas`)
- SLURM: submit `sbatch --chdir=$AMR_WORK --export=ALL,... $AMR_HOME/slurm/<s>.slurm` (WorkDir must be scratch; barbun needs ≥20 cores → use `-c40`). Interactive Faz2a needs `export AMR_FEATURE_REPR=unitig`. Per-antibiotic chain is STRICTLY SEQUENTIAL (bio 10-13b fully done BEFORE pyseer 14; verify `14_pyseer_significant_*.csv` exists before populate + before `rm rtab`). `populate`/`migrate` need target via inline env; run from `$AMR_HOME`. `git checkout origin/main -- <file>` needs a `git fetch origin` first. No `sqlite3` CLI in container → use python. amr-tools.sif DBs live in `data/external/{amrfinder_db,resfinder_db,pointfinder_db}` (container has none).

---

# §0.-3 — LATEST STATE (2026-07-07 evening) — 3rd ORGANISM (S. aureus) + K. pneu BREADTH IN FLIGHT — superseded by §0.-4 above

> **Repo HEAD `4946883`** (pushed). Two commits since §0.-2 (HEAD `2528437`): `ac041e4` (populate now fills `antibiotics.drug_class` from registry — was hardcoded `None`) + `4946883` (`scripts/kb_tables.py` tidy-CSV export + `scripts/kb_figures.py` thesis figures). Local untracked: `backup/`, `figures/`; `TRUBA_Proje_Kurulum_Rehberi.md` staged as deleted.

**This session = SCALE-OUT: adding a 3rd organism (Staphylococcus aureus) + widening K. pneumoniae to more antibiotics.** Everything in §0.-2 (17-model KB, env-parametric SLURM workflow, per-antibiotic recipe) still holds — this section only tracks the NEW in-flight work. KB is still the unified `results/kb/amrk.db`; new models APPEND (`populate_database.py` with env `AMR_ORGANISM`/`AMR_ANTIBIOTIC`, no `rm`).

### IN-FLIGHT jobs — UPDATED 2026-07-07 late-evening (all survive laptop shutdown)
- **K. pneu 4× `03u` ALL DONE** (COMPLETED 0:0): ceftazidime (6041988), amikacin (6042214), cefotaxime (6042215), levofloxacin (6042216). Matrices in `data/processed/kpneumoniae/{ab}/matrix_unitig/` (rtabs 8–15 GB each, still on disk → `rm` after each pyseer).
- **K. pneu 4× Faz1 (04→05→06→07b) RUNNING**: cefotaxime (6043613), levofloxacin (6043614), amikacin (6043616), ceftazidime (6043617). Next per-ab step after each leaves `squeue` = **Faz2a 07→08→09** (UI/`screen`, NCBI internet).
- **S. aureus FASTA download DONE**: 2494 usable `.fna` in `data/raw/staphylococcus_aureus/genomes/` (162 failed, the low-ID `1280.9xx` cluster). `00_prepare_metadata.py` built `amr_phenotypes.csv` = 2494 genomes × 41 antibiotics.
- **S. aureus cefoxitin FLAGSHIP `03u` RUNNING** (job **6043622**). First staph antibiotic → builds the organism-level `unitig_all` store once (all 2494 genomes; later staph antibiotics reuse it cheaply), then subsets to cefoxitin's 1692.
- **Registry curation DONE + pushed** (commit **`cfd9dc5`**, deployed to TRUBA via `git checkout origin/main -- config/registry/{antibiotics,organisms}.yaml`): added `oxacillin` (penicillins); split **macrolides** + **lincosamides** classes out of `others` (KB classes **7→9**); AFP keywords for amikacin/levofloxacin/oxacillin `[OXACILLIN,METHICILLIN]`/erythromycin/clindamycin + `METHICILLIN` on cefoxitin; staph `antibiotics:` = curated 8-set `[cefoxitin, oxacillin, ciprofloxacin, gentamicin, tetracycline, trimethoprim_sulfamethoxazole, erythromycin, clindamycin]`, `enabled: true`.
- **⚠️ SLURM submit gotcha (NEW):** TRUBA rejects jobs whose WorkDir is not under `/arf/scratch` (`$AMR_HOME` is `/arf/home/...` → "Lütfen /arf/scratch altında" surfaced as a QOS error). `run_03u_env.slurm` carries an internal scratch chdir, but `run_ml_env.slurm` does not. **Fix: submit with `sbatch --chdir=$AMR_WORK --export=ALL,... $AMR_HOME/slurm/<script>.slurm`** (full path to the script in `$AMR_HOME`, WorkDir forced to scratch; the script `cd`s to `$AMR_HOME` internally at runtime).

### 3rd organism = Staphylococcus aureus (taxid 1280, slug `staphylococcus_aureus`) — GRAM-POSITIVE, cross-phylum
Chosen over P. aeruginosa for the strongest "organism-agnostic" claim: gram-positive vs the two gram-negatives (E. coli, K. pneu) → best cross-phylum generalisation story.
- **`organisms.yaml` was MISSING the staph block on TRUBA** (only had ecoli+kpneumoniae → "Unknown organism" error). Fixed with `git fetch origin && git checkout origin/main -- config/registry/organisms.yaml` (origin/main HAS the staph+pseudo blocks; the `enabled` flag is irrelevant to our env-parametric workflow — 00a only needs the taxid). Verified `grep -c staphylococcus_aureus` = 3.
- **Metadata fetched** (00a `--skip-download`): `data/external/staphylococcus_aureus/metadata/amr_cleaned_long.csv` = **2656 genomes × 41 antibiotics, 24857 R/S pairs** (raw 45876 rows, 114 conflicts dropped).
- **S. aureus R/S (minority = min(R,S), with registry class):**
  | antibiotic | R | S | minority | class | note |
  |---|---|---|---|---|---|
  | ciprofloxacin | 1097 | 1119 | **1097** | quinolone | **cross-org** (E.coli+K.pneu) |
  | erythromycin | 1155 | 860 | 860 | macrolide (others) | **NEW class** (erm/msr) |
  | oxacillin | 588 | 592 | 588 | **UNREG** | **direct MRSA/mecA marker** — add to antibiotics.yaml |
  | clindamycin | 734 | 567 | 567 | lincosamide (others) | new class |
  | cefoxitin | 1258 | 543 | 543 | cephalosporin | **MRSA surrogate → mecA (FLAGSHIP, registered)** |
  | tetracycline | 454 | 1503 | 454 | tetracycline | **cross-org** (K.pneu) |
  | gentamicin | 314 | 1636 | 314 | aminoglycoside | **cross-org** (both) |
  | trimethoprim_sulfamethoxazole | 264 | 943 | 264 | folate | **cross-org** (both) |
  | penicillin | 965 | 80 | 80 | penicillin | near-universal R (imbalanced) |
  (also: fusidic acid 261 UNREG, chloramphenicol 122, daptomycin 78, rifampin 41, tigecycline 12…)
- **Planned S. aureus targets:** flagship **cefoxitin** (mecA/MRSA, registered, ready) + **oxacillin** (add to registry = direct methicillin marker) + cross-organism **ciprofloxacin / gentamicin / tetracycline / trimethoprim_sulfamethoxazole** (3-organism concordance) + new-class **erythromycin** (macrolide) / **clindamycin** (lincosamide).
- **Registry curation NEEDED before the S. aureus pipeline** (do locally + push, then `git checkout` on TRUBA): (1) `organisms.yaml` staph `antibiotics:` list is a placeholder — replace with the curated target set; (2) `antibiotics.yaml` ADD **oxacillin** (class `penicillins`, `amrfinder_keywords: [OXACILLIN, METHICILLIN, mecA]`) + confirm erythromycin/clindamycin mapping (currently "others").

### K. pneu breadth — remaining candidates (minority ≥150, NOT yet done)
Already done (11): gentamicin, tobramycin, meropenem, ciprofloxacin, imipenem, cefoxitin, cefepime, tetracycline, piperacillin_tazobactam, trimethoprim_sulfamethoxazole, **+ ceftazidime (03u running now)**. Plus amikacin/cefotaxime/levofloxacin (03u running). **Next tier after these finish** (slash-free, minority): colistin 263, cefazolin 241, tigecycline 234, ceftriaxone 225, aztreonam 216, amoxicillin_clavulanic_acid 188, trimethoprim 187. **DROP junk** `extended spectrum beta lactamase` (170 — phenotype label, not a drug).

### NEXT SESSION — resume order
1. **Check the 4 K. pneu 03u jobs** (`squeue -u $USER`; 6041988/6042214/6042215/6042216). As each leaves the queue → run its per-antibiotic chain (§0.-2 / §0.0 recipe): Faz1 `04→05→06→07b` → Faz2a `07→08→09` (UI/screen, internet) → bio `10→11` → `12→12b→13→13b` → pyseer `14` → `populate_database.py` (env `AMR_ORGANISM=kpneumoniae AMR_ANTIBIOTIC=<ab>`) → **`rm` that antibiotic's `unitigs.rtab`** to free disk. **Guard: never start an antibiotic's next phase until its SLURM job leaves `squeue`.**
2. **Check S. aureus download** (`screen -r sa_dl`; expect ~2600 usable `.fna`). When done → `00_prepare_metadata.py --organism staphylococcus_aureus` → curate registry (oxacillin + staph antibiotics list; push; checkout on TRUBA) → `03u` for cefoxitin (flagship) via env → full chain → populate.
3. **Disk:** 4 K. pneu `unitigs.rtab` in flight (~50-70 GB each). `lfs quota -u $USER /arf/scratch` before/after; delete each rtab right after its pyseer/populate. Only start the next K. pneu tier once these clear.
4. Fold new models into docs (ROADMAP §0.5 showcase, METHODOLOGY multi-organism) + eventually the Zenodo deposit (M10, still the last must-have).

---

# §0.-2 — LATEST STATE (2026-07-07) — MULTI-ORGANISM KB COMPLETE — supersedes §0.-1 and all below

**The KB is now a unified, 2-organism, 17-model AMR biomarker knowledge base.** Scaled from 3 E. coli antibiotics to **E. coli (7) + Klebsiella pneumoniae (10) = 17 models**, one unified DB at `results/kb/amrk.db` (`models.organism` distinguishes rows). All 17 ran the full pipeline (03u→04→05→06→07b→07→08→09→10→11→12→12b→13→13b→14→populate) with lineage-CV, CPSS+PFER, MDA, label-permutation, pyseer LMM, and CARD/NCBI ARO biology. Every mechanism recovered is biologically correct (below).

### KB — 17 models (unified `results/kb/amrk.db`, schema 0.4.0)
| # | organism | antibiotic | class | confirmed mechanism |
|---|---|---|---|---|
| m1 | ecoli | ampicillin | penicillin | TEM |
| m2 | ecoli | ciprofloxacin | quinolone | gyrA/parC (SNP, step 11) |
| m3 | ecoli | cefotaxime | cephalosporin | CTX-M / CMY (ESBL/AmpC) |
| m4 | ecoli | gentamicin | aminoglycoside | AAC(3)-II |
| m10 | ecoli | trimethoprim_sulfamethoxazole | folate | sul2 + dfrA15 |
| m12 | ecoli | ceftazidime | cephalosporin | CTX-M (ESBL) |
| m17 | ecoli | amoxicillin_clavulanic_acid | penicillin+inh | OXA-1 |
| m5 | kpneumoniae | gentamicin | aminoglycoside | AAC(3)-II |
| m6 | kpneumoniae | tobramycin | aminoglycoside | (AAC/ANT) |
| m7 | kpneumoniae | meropenem | carbapenem | **KPC** |
| m8 | kpneumoniae | ciprofloxacin | quinolone | gyrA/parC |
| m9 | kpneumoniae | imipenem | carbapenem | **KPC** |
| m11 | kpneumoniae | cefoxitin | cephalosporin(cephamycin) | (AmpC) |
| m13 | kpneumoniae | cefepime | cephalosporin | (ESBL/AmpC; hardest, lineage-CV ~0.74) |
| m14 | kpneumoniae | tetracycline | tetracycline | tet(A) |
| m15 | kpneumoniae | piperacillin_tazobactam | penicillin+inh | TEM |
| m16 | kpneumoniae | trimethoprim_sulfamethoxazole | folate | dfrA14 |

**KB stats (verified 2026-07-07):** 17 models, 2008 unitigs, 1902 blast_annotations, 4832 validation_evidence, 17 pipeline_runs; DB 2.2 MB. **7 drug classes** (penicillin, cephalosporin, carbapenem, quinolone, aminoglycoside, folate, tetracycline).
**Cross-organism concordance (same drug, both species → same mechanism):** gentamicin=AAC(3)-II (both), ciprofloxacin=gyrA/parC (both), trimethoprim/sulfa=dfr (both). **K. pneumoniae flagship: KPC carbapenemase** recovered for both meropenem & imipenem (100% id, E≈1e-23…1e-80).
**Signal concentration (PFER) tracks biology:** carbapenem/quinolone/ESBL = concentrated (cip PFER 0.10, mero 2.96, amox-clav 1.03); aminoglycoside = diffuse/co-carried (gent 50.6/35.7). K. pneu genomes = 4615 (2nd organism, this session), E. coli = 5470.

### Infrastructure built this session (all pushed to `main`, HEAD `2528437`)
1. **Parallelism refactor (`a23dc40`):** 15 pipeline scripts (03u,04,05,06,07,07b,08,09,10,11,12,13b,14,16,populate) now resolve (organism,antibiotic) via `lib.config.get_target()` — precedence **CLI-arg > `AMR_ORGANISM`/`AMR_ANTIBIOTIC` env > config.yaml**, backward-compatible. 12b/13 inherit 05's globals. **Removes the config.yaml mutex → many (organism,antibiotic) pipelines run in PARALLEL via per-job env**, no config edits. 117 pytest pass.
2. **Unified multi-organism KB (`2528437`):** `populate_database.py` default DB → `results/kb/amrk.db` (was per-organism `results/{org}/kb`); schema already carries `models.organism`. Pass `--db` to override.
3. **Slash-safe antibiotics (`2528437`):** `antibiotics.yaml` slash canonicals → underscore (`trimethoprim_sulfamethoxazole`, `amoxicillin_clavulanic_acid`, `piperacillin_tazobactam`) with slash spellings kept as aliases; class mapping unchanged. Fixes path/filename breakage so combo drugs run as ML targets. **After deploy, metadata was rebuilt** (re-normalise `amr_cleaned_long.csv` → `00_prepare_metadata.py`) for both organisms so columns are underscore.
4. **Genome backup:** E. coli + K. pneu assemblies (~52 GB, 10085 `.fna`) tar.gz'd on a compute node (gzip needs no CPU-ulimit-limited login node → do it via SLURM, NOT `tar czf|rclone rcat` on a UI/transfer node which gets killed at ~394 MB) → `rclone copy` the 15.5 GB tarball to `gdrive:TRUBA_25626/scratch_amr/backup/`.

### Env-parametric SLURM workflow (the per-antibiotic recipe — REUSE THIS)
Generic env-driven scripts in `$AMR_HOME/slurm/`: `run_03u_env.slurm`, `run_ml_env.slurm` (04-07b), `run_bio_env.slurm` (10-13b), `run_pyseer_env.slurm` (two-container 14). Submit with `sbatch --export=ALL,AMR_ORGANISM=<o>,AMR_ANTIBIOTIC=<a> <script>`. Interactive-only step = **2a (07→08→09)** on the UI node (08 NCBI needs internet; run in `screen` — it survives SSH drops; NCBI is SLOW ~20 min for some, not hung). `populate` + `rm unitigs.rtab` after each antibiotic's pyseer (rtab is pyseer-only, ~40-70 GB, regenerable → delete to free disk). **Guard: never start an antibiotic's next phase until its current SLURM job leaves `squeue`** (same `results/{org}/{ab}/` dir). Full audit anytime: KB `SELECT model_id,antibiotic,run_id FROM models` + filesystem stage scan (features.txt/config_{ab}.yaml/05_final_biological_report.md/13_stability_summary/14_pyseer_summary).

### Remaining / next
- **Docs:** fold the 17-model panel + cross-org story into `docs/ROADMAP.md` §0.5 (showcase) and METHODOLOGY (multi-organism scale-out, unified KB, per-antibiotic PFER/mechanism table). Novelty reframe now supports "**cross-organism, multi-class, PFER-bounded, lineage-validated open AMR biomarker KB**".
- **M13/M10 still open:** external concordance (16) for the new antibiotics needs their `amrfinder_keywords` in `antibiotics.yaml` (only amp/cef/cip present); Zenodo deposit (M10) still the last must-have — deposit the unified 17-model KB.
- **Disk:** all per-antibiotic `unitigs.rtab` deleted after pyseer. K. pneu raw `.fna` present; genome backup on Drive.

---

# §0.-1 — LATEST STATE (2026-07-03) — superseded by §0.-2 above (historical)

**The thesis is essentially research-complete.** 3 antibiotics in the KB, all must-haves done except the Zenodo deposit, full audit + fixes done. What follows below (§0, §0.0…) is earlier/historical detail; where it conflicts, THIS section wins.

### ⏳ IN-FLIGHT (2026-07-03 session — laptop-independent, resume on return)
Two long jobs launched this session; both survive laptop shutdown (verified).
1. **4th antibiotic = gentamicin (E. coli)** — aminoglycoside = **NEW drug class** (KB had 2 β-lactam + 1 FQ; chosen over trimethoprim/sulfa for clean single-drug mechanism). `config.yaml target_antibiotic=gentamicin` set on TRUBA (via `sed`). **`03u` unitig-caller RUNNING** (SLURM job 6022783, barbun55, `-c40 --mem300G`, `AMR_FEATURE_REPR=unitig`, `--threads 40`) on **4194 genomes (646 R / 3548 S, 15.4% R)** — 80 QC-outliers excluded. On completion → `data/processed/ecoli/gentamicin/matrix_unitig/`. Then per-antibiotic playbook (§0.0): Faz 1 `04→05→06→07b` (full node, `AMR_EXTERNAL_MEMORY=false AMR_OPTUNA_PATIENCE=15`) → Faz 2a/2b/3/3c → `populate_database.py --antibiotic gentamicin` (**no `rm`**, appends model_id 4). Add `amrfinder_keywords: gentamicin: [GENTAMICIN, AMINOGLYCOSIDE]` to `antibiotics.yaml` before its M13.
2. **2nd ORGANISM = K. pneumoniae (taxid 573, slug `kpneumoniae`) — SCALE-OUT.** Metadata fetched: **4976 genomes × 70 antibiotics** (`data/external/kpneumoniae/metadata/amr_cleaned_long.csv`). **FASTA download RUNNING** in `screen -r kpneu_dl` on **arf-ui1** (BV-BRC API, `--raw-csv` reuse, `--workers 12`; ~10% "empty/non-FASTA" fails concentrated in the low-ID 573.12xx cluster then thinning → expect **~4480 usable**; **resume-safe** — skips existing `.fna`, re-run same cmd or `--retry-failed`). ⚠️ Compute nodes have NO internet → download MUST stay on a UI/transfer node.
   - **Target = gentamicin (DECIDED — cross-organism story).** Same drug as E. coli's 4th → "does the unitig-KB / stable biomarkers generalise across species?" (S1-style cross-organism overlap). K. pneu gentamicin on the 4615-genome set: **1268 R / 2055 S (3323 tested, 38% R)** vs E. coli's 646/3548 (15% R) — nice balance contrast, both class-weighted. (Rejected alternatives: meropenem 1045/2133 = carbapenem/new-class; trimethoprim/sulfa 2108/1047.) **`03u` submitted** `--organism kpneumoniae --antibiotic gentamicin` (job PENDS behind E. coli gentamicin Faz 1 — 40-core MS limit; auto-starts when it frees) → `data/processed/kpneumoniae/gentamicin/matrix_unitig/`.
   - **Registry curation BEFORE its pipeline:** `organisms.yaml` set `kpneumoniae enabled: true` + a **curated** `antibiotics:` list (current is placeholder `[meropenem,ciprofloxacin,gentamicin,colistin]`); `antibiotics.yaml` add `amrfinder_keywords` for the chosen target. K. pneu UNREGISTERED names seen: **DROP junk** `extended spectrum beta lactamase` / `fluoroquinolones` / `aminogycosides` (phenotype/class labels, not drugs — never make them ML targets); optionally add real minors `cefiderocol`→cephalosporins, `spectinomycin`, alias `cefalexin`→cephalexin (all tiny/zero-R, not targets). Unknown names are harmless (`normalize_antibiotic` passes through, `antibiotic_to_class`→None).
   - Next: download done → `00_prepare_metadata.py` (builds `amr_phenotypes.csv` from present `.fna`) → pick target → curate registry → full playbook.
- **Disk watch:** `/arf/scratch` ~802 GB / 1 TB used, inodes 118.8K / 200K. Gentamicin + K. pneu unitig `unitigs.rtab` (~70-90 GB each) will pressure quota → clean regenerable amp/cip/cef `unitigs.rtab` first when tight.

### 🚀 PARALLELISM: env-driven (organism, antibiotic) — config-mutex REMOVED (2026-07-04)
The pipeline scripts used to read `config['project']['target_antibiotic']`/`organism` **directly**, so only ONE (organism, antibiotic) could run at a time (config.yaml was a global mutex; env override "leaked" because scripts bypassed `get_target`). **FIXED:** 15 scripts now resolve the target via **`lib.config.get_target()`** — precedence **CLI-arg > `AMR_ORGANISM`/`AMR_ANTIBIOTIC` env > config.yaml**, backward-compatible (no env ⇒ config, so existing runs are unchanged). Refactored: `03u,04,05,06,07,07b,08,09,10,11,12,13b,14,16,populate_database` (12b/13 inherit 05's globals). Verified: 117 pytest pass, functional env-override test passes. **Now: set `export AMR_ORGANISM=… AMR_ANTIBIOTIC=…` per SLURM job → run many (organism,antibiotic) pipelines in PARALLEL without touching config.yaml.** (Organism-level QC 01/02/02b/02p/03/03b still read config directly but take `--organism` where used, or are baseline/skippable; refactor them too if a new organism needs parallel QC.) **Deploy to TRUBA:** push, then `git checkout origin/main -- scripts/{04,05,06,07,07b,08,09,10,11,12,13b,14,03u,16,populate}*.py` (none are in the forbidden list — 03 untouched). Safe to checkout anytime (backward-compat); the currently-running gentamicin run keeps working.

### KB (authoritative, TRUBA `results/ecoli/kb/amrk.db`, schema 0.4.0, deduped)
3 models, **acquired-gene vs target-SNP showcase across two β-lactam mechanisms + a fluoroquinolone**:
| antibiotic (model) | genomes R/S | lineage-CV AUC | CPSS stable/PFER | pyseer sig | mechanism (confirmed) |
|---|---|---|---|---|---|
| ampicillin (1) | 4373 (2717/1729) | 0.9511±0.011 | 36 / 2.73 | 25/36 | **acquired** TEM β-lactamase; H2 TRUE 47%; SNP 0 |
| ciprofloxacin (2) | 4150 (1324/2826) | 0.9496±0.007 | 70 / 8.4 | 5/70 | **target SNP** gyrA S83L + parC S80I (step 11) |
| cefotaxime (3) | 3788 (1009/2779) | 0.9546±0.020 | 83 / 13.2 | 17/83 | **acquired** CTX-M-276/278 + CMY ESBL/AmpC; H2 FALSE 25% (ok) |
KB provenance complete (git 5b76f47 + seed 42 + config_hash + CARD 4.0.1; n_genomes/min_support backfilled). `validation_evidence` deduped 980→902.

### Roadmap status: **must-have 15/16, should-have ~all**
- **DONE:** M1-M9, M11, M12, M14, M16, **M13** (concordance + head-to-head), **M15** (CheckM2+QUAST 97.1% pass 5312/5470), M5(minus Zenodo). Should: S1(cross-antibiotic/H3), S4, S5, **S7**(ResFinder), **S8**(FastAPI `scripts/kb_api.py`), **S9**(FAIR /metadata), **S10**(METHODOLOGY §4.4).
- **REMAINING:** **M10** Zenodo deposit (only external step left — `docs/RELEASE_ZENODO.md`; `.zenodo.json`/CITATION ready; stamp DOI via `AMR_ZENODO_DOI` env; deposit ideally v0.5.0 with all 3 antibiotics). **M13 Track A** temporal/geo hold-out (deferred by user, needs BV-BRC year/geo metadata fetch). **Thesis writing** of audit notes 10/21/22/14 (below).

### Key results (thesis headlines)
- **H3 REJECTED (biologically substantive negative finding):** within-β-lactam ampicillin~cefotaxime share **no** stable unitig **and no gene family** (ampicillin=TEM vs cefotaxime=CTX-M/CMY — same class, distinct enzymes); the only β-lactamase overlap is cross-class cefotaxime~ciprofloxacin CTX-M **co-carriage**. `scripts/15_cross_antibiotic.py` (unitig + ARO gene-family level). H3 hypergeometric is `--with-test` (deferred, union-universe caveat).
- **M13 head-to-head (leakage-free, held-out test genomes):** unitig model bACC amp 0.873 / cef 0.925 / cip 0.928 — **matches ResFinder (cef,cip), beats AMRFinderPlus (cip,amp)**. Tool-vs-phenotype κ: ResFinder amp 0.86/cef 0.89/cip 0.76; AMRFinderPlus over-calls amp+cip (naive determinant→class map incl. intrinsic β-lactamase/efflux → ME ~38-40%; cef clean 0.87). All in KB `validation_evidence` (concordance_amrfinderplus/resfinder/head_to_head_model).

### New tooling (this session, all on main)
`scripts/16_external_concordance.py` (M13, prep/post + `--write-kb`), `scripts/lib/concordance.py` (bACC/sens/spec/κ/McNemar/FDA ME-VME), `scripts/kb_api.py` + `scripts/lib/kb_queries.py` (S8/S9 REST API, `uvicorn scripts.kb_api:app`), `scripts/kb_report.py` (one-command thesis results Markdown), `scripts/kb_app.py` (Streamlit + H3/M13 tabs), `scripts/02d_genome_qc.py` (M15). Tests: 117 pass, 1 skip (fastapi). `pip install fastapi uvicorn` to run the API/its test.

### AUDIT (independent technical audit done 2026-07-03) — fixed vs deferred
**Fixed+pushed:** High 1 (METHODOLOGY reconciliation banner — it still documents the k-mer *baseline*; unitig/lineage-CV/CPSS is canonical), 2 (README Results), 3/24 (**populate idempotency** — `kb_schema.ensure_unique_indexes()` + `INSERT OR IGNORE`; blast/evidence were double-inserted by 09+13b); Medium 4 (run_metadata n_genomes/min_support), 5 (Genome-ID `dtype=str` in 06/16/lineage — "562.10"→float hazard), 6 (run_pipeline `--antibiotic auto` was broken: value_counts on 0/1 matrix), 8 (16 tool versions env-overridable), 9 (AFP keywords → `antibiotics.yaml` registry); Low 17/23/25/28. **Deferred (gerekçeli):** Medium 7 (`io_utils.run_command` string→list refactor, risky), Medium 10/21/22 + Low 14 = **thesis-writing doc notes** (10: AMRFinderPlus naive determinant→class map over-calls → frame in Methods that ResFinder is the curated-phenotype reference; 21: "recovery rate" = precision-of-stable-set not recall, define precisely; 22: HPO/06 use chunk-split not lineage-aware — only 07b's reported AUC is lineage-CV, note the defense; 14: git_dirty=1 on all runs — document the TRUBA manual patches to 02p/02b/03/config). Low deferred: 11 (kmer column name = unitig seq, rename is breaking), 12 (`same_class` = registry class, add `same_drug_family`?), 13 (API rate-limit), 15 (ampicillin 12b nulls csv missing), 18/19/26/29 (documented low-impact). **No critical/data-corrupting bug found; code quality high.**

### Operational (unchanged, critical)
- **Push:** ONLY when the user asks + pastes a fresh fine-grained PAT (Contents:write). **No `Co-Authored-By` trailer.** Push via `git push "https://demirbase:<PAT>@github.com/demirbase/ML_AMR_Prediction_v2.git" main`. The auto-mode classifier blocks pushes unless the user explicitly authorized this turn. **The repo is a FORK of `iumobg/ML_AMR_Prediction_v2` → commits DON'T count toward the contribution graph** (make it standalone / detach fork to fix — pending user decision).
- **TRUBA:** no `git pull` — targeted `git checkout origin/main -- <file>` (never `config.yaml`/`02p`/`02b`/`03`). Submit from `$AMR_WORK`, `APPTAINER_BINDPATH=/arf`, **`--exclude=barbun45,barbun46`** (Apptainer `lookup userid` glitch nodes) + **`apptainer exec --no-home`** (the real fix for the glitch). barbun min 20 cores/node. Containers: `amr.sif` (core), `amr-tools.sif` (pyseer/quast/amrfinderplus 4.2.7/resfinder 4.5.0), `amr-checkm2.sif`, `amr-pp.sif` (poppunk). Env overrides: `AMR_FEATURE_REPR=unitig AMR_EXTERNAL_MEMORY=false AMR_OPTUNA_PATIENCE=15 AMR_CARD_VERSION=4.0.1 AMR_ENTREZ_EMAIL AMR_ZENODO_DOI`.
- **Drive backup:** `rclone` to `gdrive:TRUBA_25626/scratch_amr` on transfer host **arf-ui4 (172.16.6.14)** in `screen`. **NEVER `rclone copy $AMR_WORK` whole-tree** (OOM-killed ~1h walking 100GB+). Copy small subdirs individually (KB/results/models/runs); exclude regenerable big data (`unitigs.rtab`, `matrix_unitig/X_*.npz`, `kmc_outputs`, `genomes/*.fna`, `*_db`). Latest science outputs already backed up.
- **CheckM2 DB:** at `$AMR_WORK/data/external/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd` (pass `--database_path`; the download's json-write fails on read-only FS but the .dmnd lands fine). AMRFinderPlus DB `amrfinder_update -d $AMR_WORK/data/external/amrfinder_db`; ResFinder DBs cloned to `$AMR_WORK/data/external/{resfinder,pointfinder}_db`.

---

# 0. TRUBA (ARF) deployment — LIVE STATE (resume here)

> A full real run is in progress on the **TRUBA ARF** cluster (user `edemirbas`). New session: continue from here.

**Where we are:** Data acquired (5470 *E. coli* genomes). Pipeline on TRUBA: `00a✓ 00✓ 01✓ 02 KMC✓ 02b QC✓ 03 matrix✓` — **FEATURES DONE**: `data/processed/ecoli/ampicillin/matrix/` has **22 `.npz` chunks + `features.txt` (1.27 GB)** = **4446 ampicillin genomes × 50.8M k-mers** (~90% sparse; full matrix ≈ **109 GB** decompressed / 21.8B nnz).

**Training regime refactored to full-data boosting (DONE & re-run on TRUBA 2026-06-21).** The old run (2026-06-18, incremental 1-tree/chunk) gave test ROC-AUC 0.903 / MCC 0.693 / acc 0.84. Refactored to **standard full-data gradient boosting** (`scripts/lib/xgb_data.py`); 04/05/07b updated and pushed to `main` (HEAD `a5b9ddc`). Key engineering hurdles solved on TRUBA: 04 HPO now runs **parallel Optuna trials** over a `QuantileDMatrix` subset (`training.optuna_threads_per_trial`, default 2) to use cores without OOM; 05/07b use **`ExtMemQuantileDMatrix`** (external memory, pages spilled to scratch) because an in-core full-train QuantileDMatrix peaked **>400 GB** and OOM-killed the 384 GB barbun node.

**New result (2026-06-21, 4373 genomes, 50 HPO trials, full-data boosting):** test **ROC-AUC 0.930 (CI 0.914–0.945), PR-AUC 0.965, MCC 0.739, acc 0.866, balanced-acc 0.885** (threshold 0.5). Improved over the old regime across the board (MCC 0.693→0.739). NB: 04+05+06 jobs ran to completion; the TRUBA **low-efficiency warning DID fire** during ExtMem training (Eff ~10%, I/O-bound) but **did NOT kill the job** — confirms the warning is a nag, not an auto-killer (every actual death this project was OOM or manual `scancel`).

**UNDERFITTING HYPOTHESIS TESTED → DISPROVEN (2026-06-21/22).** Added full-data early stopping to 05 (`training.max_boost_rounds`, `validation_fraction` for the ES split) and ran it: early stopping found **best tree count = 29** (≈ the HPO subset's 30), val AUC peaks there, more trees overfit. So the model is **well-fit, not underfit**; ROC-AUC ~0.93 / MCC ~0.73 is the **realistic signal ceiling** for ampicillin k-mers (consistent with literature), not a tuning bug. The compressed probability range is just low-lr × few-trees, not capacity-limited. **No ML tweak (early stopping, more trees, min_support) will raise this ceiling — it lives in the feature representation.** Note also: `min_support` removes rare k-mers but barely reduces nnz (dominated by common k-mers), so it does NOT shrink RAM (only fewer GENOMES via `training.max_train_chunks` does); its real value is lineage/noise de-confounding of candidate k-mers, a KB-quality lever, not a prediction lever.

**STRATEGIC PIVOT (agreed): stop optimising the prediction engine; move to the thesis contribution.** Per `docs/ROADMAP.md` the research question is the **queryable, biologically-validated AMR k-mer Knowledge Base**, NOT prediction accuracy — AUC 0.93 is already more than enough to rank k-mers for the KB. ROADMAP M2 accepts **5-seed repeated holdout OR 5-fold CV** (we use 5-seed = 07b); no must-have requires a prediction-accuracy bar. So the training regime is not a thesis lever; we keep **standard full-data boosting** (more conventional/defensible than the old incremental 1-tree/chunk, which a reviewer would question).

**GPU evaluated and REJECTED (2026-06-22).** Built `amr-gpu.sif` (CUDA xgboost, `amr-gpu.def`, USE_CUDA True) and smoke-tested on a Tesla **V100 16 GB** (akya-cuda/barbun-cuda, AllowAccounts=ALL). The 50.8M-feature ultra-wide sparse matrix is GPU-hostile: in-core GPU OOMs (3 chunks alone wanted 11 GB — GPU ELLPACK ≈4 bytes/nnz → full set ≈88 GB ≫ 16 GB) and ExtMem-GPU did not even finish building in 25 min. **CPU is the path.** (`amr-gpu.def` kept in repo for the record; don't pursue GPU on this hardware for this data.)

**ENGINE SETTLED — feature filter made DATA-ADAPTIVE (the real lever for speed + KB quality + generality).** Since the AUC ceiling is feature-bound and full-data training is slow only because of the 50.8M features, 03 now derives min_support adaptively: `min_support = max(min_support_floor=5, ceil(min_prevalence=0.01 * n_genomes))` (config `preprocessing`, with `min_support: null` = auto, or an int to force). This scales across antibiotics/organisms (ampicillin 4373 → 44; 1788 → 18; ≤500 → floor 5), so a small dataset is never over-filtered and a large one gets de-confounding + ~2× faster training. Biologically safe: 1% prevalence is far below step-10's ~10%-prevalence discriminativeness threshold, so no individually-strong marker is dropped (only the rare/lineage/error tail — exactly the confounders ROADMAP §⚠️/S3 flag).

**LITERATURE REVIEW DONE (2026-06-22) → MAJOR METHODOLOGICAL PIVOT (see `docs/ROADMAP.md` §0).** Two rounds of systematic literature research were completed and distilled into binding decisions now in ROADMAP **§0** (authoritative; supersedes older ROADMAP sections). The current raw-k-mer + adaptive-min_support engine WORKS and is a valid baseline, BUT the literature mandates a publication-grade overhaul. Key decisions:
- **Unitigs replace raw k-mers** (`bcalm2` + `unitig-caller`): ~10M→~730k features, ~212→~18 GB, ~7h→~50min, BLAST-mappable, GWAS-standard. Downstream XGBoost unchanged (binary matrix). **This dissolves the min_support/speed/GPU pain we fought all session** — so the raw-k-mer 03 rebuild + GPU work are now mostly moot.
- **Lineage-aware CV** (PopPUNK clusters → `GroupKFold`) replaces random/chunk split. Random CV inflates AUC 20-30%. **Biggest reviewer-blocker.** Final model still trained on all data.
- **Stability: CPSS (B=100, 50%, π≥0.6) + Chi²/MI staged prefilter + SHAP** replaces 5-seed 07b. Importance: **SHAP** not Gain.
- **MUST add:** external validation (temporal/geographic hold-out + AMRFinderPlus/ResFinder concordance: Kappa/McNemar/bACC); pyseer LMM+Bonferroni; CheckM2+QUAST QC; BH-FDR on step 10; ARO/CARD ontology mapping in the KB; reification-safe wording.
- **Confirmed (no change):** k=21, binary+max_bin=2, class-weight (no SMOTE), BV-BRC+EUCAST/CLSI, 4373 genomes adequate, AUC ~0.93 consistent with literature.
- **Novelty reframe:** NOT "first ML AMR DB" (BV-BRC exists) → "first PFER-bounded (CPSS), lineage-validated, k-mer/unitig-resolution, transparent+FAIR open AMR biomarker KB." Target: Database(Oxford)/Briefings in Bioinformatics. Showcase: ciprofloxacin (gyrA/parC SNP) + a β-lactam (acquired gene).

**TRUBA state at handoff:** raw-k-mer 03 matrix rebuild (job 5952577, adaptive min_support=44 → 21.4M features) was running/likely done — **but it is now superseded by the unitig pivot, so do NOT keep building on it.** KB schema (`scripts/lib/kb_schema.py`) + `scripts/populate_database.py` (M8 foundation) are built, tested, committed (HEAD `main`). CARD 4.0.1 full bundle downloaded (step 11 active). TRUBA scratch cleaned (GPU `.sif`, kmc_curve, smoke removed); pending: `rm -rf $AMR_WORK/data/interim/ecoli/kmc_outputs/tmp/*` after 03.

**IMMEDIATE NEXT (resume here):** The full pipeline (04→populate) is DONE end-to-end and **re-run cleanly ("canonical") for TWO antibiotics** with full provenance. KB at `results/ecoli/kb/amrk.db` is now **multi-antibiotic, unitig-named, schema 0.4.0**.
- **Ampicillin (canonical, 2026-06-28):** lineage-CV **0.9511±0.011**, H2 **TRUE** (recovery 47%, 28 confirmed); CPSS 36 stable / **PFER 2.73** (TEM-256/257/258, tet(A)); MDA 0-sig (redundancy); label-perm p=0.0196; pyseer **25/36** lineage-sig (3/3 TEM pass, aminoglycoside co-resistance does NOT); SNP 0 resistant-allele (expected — acquired-gene mechanism). 303 evidence rows.
- **Ciprofloxacin (canonical, 2026-06-29):** lineage-CV **0.9496±0.007**; CPSS 70 stable / PFER 8.4; label-perm p=0.0196; pyseer 5/70 (signal more lineage-entangled — clonal FQ spread). **SNP showcase ✓: gyrA S83L + parC S80I = resistant_allele** (the textbook FQ mutations) — CARD-homolog recovery ~0 (expected: SNP not acquired gene), validated by step 11 instead. The two antibiotics = the **acquired-gene vs target-SNP** showcase pair.
- **Reproducibility:** every run stamps git_commit (5b76f47) + seed + config_hash + CARD 4.0.1 in `pipeline_runs` (run_id no longer `…__unknown`). Local **Streamlit KB explorer** `scripts/kb_app.py` (queryable interface, S8/N1).
- **→ NEXT:** (a) 3rd antibiotic **cefotaxime** (β-lactam, balanced) then gentamicin — same chain (set `target_antibiotic` in config.yaml, run 03u → 04→14 → populate WITHOUT `rm` so it appends to the KB); (b) thesis figures (null histogram, mechanism comparison, evidence-chain); (c) M13 external validation, M15 CheckM2/QUAST QC; S1 cross-antibiotic overlap (now 2 antibiotics in KB), S8 FastAPI.

**Per-antibiotic playbook (cefotaxime/gentamicin):** `sed -i 's/target_antibiotic: .*/target_antibiotic: cefotaxime/' config/config.yaml` → 03u (unitig-caller, ~hours, full node) → Faz1 `04 05 06 07b` (full node, `AMR_EXTERNAL_MEMORY=false AMR_OPTUNA_PATIENCE=15`) → Faz2a `07 08 09` (UI, internet, `AMR_ENTREZ_EMAIL` + `NXF_ANSI_LOG=false`; NCBI remote may need one retry) → Faz2b `10 11` (compute) → Faz3 `12`(`--n-permutations 100`) `12b`(`--n-permutations 50`, full node) `13`(`--base-trees 10`) `13b` → Faz3c `14` (two-container) → `populate_database.py --antibiotic <ab>` (`AMR_CARD_VERSION=4.0.1`, no rm). **Always `sbatch --exclude=barbun45`** (Apptainer `lookup userid` glitch node).

## §0.0 CANONICAL multi-antibiotic state (2026-06-29) — READ FIRST

The whole pipeline was **re-run cleanly from scratch (after the unitig matrix) for two antibiotics** so every artefact is saved with provenance. This supersedes the earlier ampicillin numbers in §0.1/§0.2 (those were the first, partly hand-patched run). The KB (`results/ecoli/kb/amrk.db`) is now **multi-antibiotic, unitig-named, schema 0.4.0**.

### Why re-run
The first ampicillin run had hand-written bits: `config_ampicillin.yaml` was written by hand (04 HPO was cut short, n_estimators=8), so it was not git-reproducible, and `run_metadata.json` was missing → KB run_id was `…__unknown` (no git hash). The canonical re-run starts at **04 (HPO included)** → reproducible config + run_metadata + git provenance.

### Pipeline phases (per antibiotic — the playbook)
Set the target once, then run the phases. **Always `sbatch --exclude=barbun45`** (that node throws the Apptainer `Couldn't determine user account information: lookup userid` glitch — env, not code).
- **Switch target:** `sed -i 's/target_antibiotic: .*/target_antibiotic: <ab>/' config/config.yaml` (universal — scripts read config both via `load_config` and direct `yaml`, so config.yaml is the only reliable switch; an env override would leak).
- **Faz 0 — `03u`** unitig-caller on the antibiotic's genome set → `data/processed/ecoli/<ab>/matrix_unitig/` (no organism `unitig_all` store exists, so each antibiotic runs unitig-caller fresh, ~hours). Full node `-c40 --mem 300G`, `AMR_FEATURE_REPR=unitig`.
- **Faz 1 — `04 05 06 07b`** (full node, in-core): `AMR_FEATURE_REPR=unitig AMR_EXTERNAL_MEMORY=false AMR_OPTUNA_PATIENCE=15`. 04 HPO is the long pole (~hours; high CPU load = good; the brief low-eff dip is the serial 05 DMatrix build — a nag, not fatal).
- **Faz 2a — `07 08 09`** on the **UI node** (internet for NCBI/Entrez): `AMR_FEATURE_REPR=unitig AMR_ENTREZ_EMAIL=… NXF_ANSI_LOG=false`. 08 runs CARD (local, blastn-short/word7) + NCBI (remote, blastn/word11, `txid562[Organism:exp]`). **NCBI remote can transiently fail** (`Connection stream is in bad state`) → just re-run the NCBI blastn directly (idempotent), it's a network blip not a code/param issue.
- **Faz 2b — `10 11`** (compute, offline).
- **Faz 3 — `12`(`--n-permutations 100`, fast MDA) · `12b`(`--n-permutations 50`, full node, slow: N=50 is now standard) · `13`(`--n-candidates 5000 --B 100 --pi 0.6 --base-trees 10`) · `13b`** (compute).
- **Faz 3c — `14`** pyseer LMM (two-container SLURM: prep/post in `amr.sif`, `similarity_pyseer`+`pyseer --lmm` in `amr-tools.sif`; kinship from a subsampled Rtab (~every 100th unitig), LMM on the 5000 Chi² candidates — genome-wide pyseer over the full 60 GB Rtab is impractical (16 h, OOM)).
- **Faz 4 — `populate_database.py --antibiotic <ab>`** (`AMR_CARD_VERSION=4.0.1`). **Do NOT `rm` the db for the 2nd+ antibiotic** — populate appends (new model_id, new unitigs/evidence). Only `rm` for a from-scratch single-antibiotic rebuild.

### Results (both canonical)
| | ampicillin (model_id 1) | ciprofloxacin (model_id 2) |
|---|---|---|
| genomes (R/S) | 4373 (2717/1729) | 4150 (1324/2826) |
| unitigs | 4.94M | 4.62M |
| **lineage-CV ROC-AUC** | **0.9511 ± 0.011** | **0.9496 ± 0.007** |
| 06 single-split AUC | 0.924 | 0.980 |
| trees (early stop) | 146 | 118 |
| CARD recovery / H2 | 47% / **TRUE** | ~0% / FALSE *(expected — SNP not acquired gene)* |
| CPSS stable / PFER | 36 / **2.73** | 70 / 8.4 |
| MDA significant | 0 (redundancy) | 0 (redundancy) |
| label-perm p | 0.0196 | 0.0196 |
| pyseer lineage-sig | **25/36** (3/3 TEM pass; aminogly. co-res NOT) | 5/70 (signal more lineage-entangled — clonal FQ) |
| **SNP (step 11)** | **0 resistant-allele** (expected) | **gyrA S83L + parC S80I = resistant_allele** ✓ |
| confirmed biomarkers | TEM-256/257/258, tet(A), APH(6)-Id, AAC(6')-Ib7 | efflux/co-res (gyrA/parC SNPs aren't CARD homologs) |

**The pair is the thesis showcase:** ampicillin = **acquired gene** (β-lactamase, homolog-BLAST + H2), ciprofloxacin = **target-gene SNP** (gyrA/parC, validated by step 11 not H2). Each mechanism validated by the appropriate tool.

### KB schema renamed k-mer→unitig (0.3.0→0.4.0)
`kb_schema.py`: `kmers`→`unitigs`, `kmer_id`→`unitig_id`, `kmer_model_scores`→`unitig_model_scores`, `kmer_background_frequency`→`unitig_background_frequency`, `kmer_antibiotic_overlap`→`unitig_antibiotic_overlap`, `kb_metadata.n_kmers`→`n_unitigs`. `populate_database.py` matches. NB: candidate-CSV **column** reads stay `kmer` (that's the on-disk column name from 07/10/13), and output **filenames** are unchanged (e.g. `10_kmer_background_frequency_<ab>.csv`).

### Local KB access — `scripts/kb_app.py` (Streamlit, S8/N1)
`pip install streamlit pandas` then `streamlit run scripts/kb_app.py` → point at `amrk.db`. Tabs: biomarkers (filter by method/tier/stability/gene), per-unitig **evidence chain** (BLAST+ARO, discriminativeness, CPSS, permutation, pyseer LMM), model + provenance. Get `amrk.db` to the Mac via `scp` from the transfer host (it's tiny) or the Drive backup.

### Bugs fixed during the canonical re-run (don't reintroduce)
1. `populate_run` used wrong run_metadata keys (`git_commit`/`seed`/`created_at`) and bound `data_fingerprint` (a dict) → crashed once run_metadata.json existed. Fixed to `git_commit_hash`/`random_seed`/`started_at` + store `data_fingerprint.sha256` in `config_hash`.
2. `12b` aligned test labels via `y[test_mask]` but loaded X_test in the config's (non-ascending) `test_files` order → REAL AUC collapsed to ~0.49. Fixed: build the test matrix in **ascending chunk order**.
3. `13` CPSS base learner inherited the model's 66/146 trees → ~317 features/fit → **PFER blew up to ~100**. Fixed: `--base-trees 10` (sparse base selector, decoupled from the final model) → PFER ~2.7. The final SHAP model still uses the full tree count.

### Backup / provenance
Drive backup via `rclone` to `gdrive:TRUBA_25626/scratch_amr` (incremental; transfer host arf-ui4, in `screen`). Repo is on GitHub. Every `pipeline_runs` row stamps git_commit 5b76f47 + seed + config_hash(sha256) + CARD 4.0.1.

## §0.1 Unitig + lineage-CV + biology — STATUS (2026-06-24, ampicillin FIRST run — see §0.0 for the canonical re-run)

### M9 — permutation significance (DONE; one resubmit pending) (2026-06-24)
Two complementary permutation tests written (commits `2306f61`→`7b69af4` on `main`), both reusing existing infra (06's exact held-out split, 07b's frozen-HP loader):
- **`scripts/12_permutation_test.py` — MDA (per-feature permutation importance, ROADMAP §0.2).** Model fixed (no retrain); permute each candidate unitig's column in the held-out test set; measure ROC-AUC drop; BH-FDR. **RAN & SAVED** (`12_permutation_test_ampicillin.csv` + `12_permutation_summary_ampicillin.json`): baseline AUC **0.9534** (reproduces headline), **51/60** candidates model-used, **0 significant at Q<0.05**. **This is expected, not failure:** unitigs are redundant (a β-lactamase gene → many overlapping unitigs), so permuting one is compensated by correlated partners → low individual MDA. Positive: the **top-MDA unitigs are the real genes** (CTX-M-260, CMY-198, OXA-1042). Interpretation = motivation for M4 CPSS (grouped/stability). Impl notes: binary feature → redistribute the 1s keeping count; only changed rows re-predicted; column set via sparse `maximum`/`eliminate_zeros` (no `tolil` — that tripped the **UI-node CPU ulimit**; run MDA on compute or it dies).
- **`scripts/12b_label_permutation_test.py` — label-permutation null (model-level significance, ROADMAP §1.7).** Shuffle ALL labels, retrain frozen-HP 8 trees, build the null ROC-AUC distribution; empirical p. **RESULT ESTABLISHED:** the first (naive) run printed **REAL AUC = 0.9534 vs null ~0.48–0.53** (30 perms, running max 0.5345) → real ≫ null → **model is highly significant (p → 1/(N+1) ≈ 0.0099)**. That naive run was slow (~17 min/perm: it rebuilt the DMatrix every permutation → I/O-bound, **TRUBA low-eff killer** risk) so it was cancelled; rewritten to **build the train matrix ONCE as a streamed in-core `QuantileDMatrix`** (max_bin=2, ~211G; supports `set_label`/`set_weight`, so perms just swap labels + refit — verified locally, ~100× faster). **SAVED (2026-06-24, job 5970202):** `12b_label_permutation_summary_ampicillin.json` → **REAL 0.9534, null mean 0.4994 / max 0.5521 / std 0.021, 0/100 ≥ real, empirical_p = 0.0099, significant=true** + `12b_label_permutation_nulls_ampicillin.csv` (for the thesis null-distribution histogram). Ran on a full node `-c40 --mem 300G` (in-core QuantileDMatrix ~211G), ~minutes. **Apptainer node gotcha** (cost 2 failed submits): some barbun nodes throw `FATAL: Couldn't determine user account information: user: lookup userid …` at container start (env, not code; node-specific) — just resubmit, or `sbatch --exclude=<node>` / `apptainer exec --no-home`.
- **TODO after 12b saves:** add `validation_evidence(evidence_type='permutation', evidence_score=mda_auc_drop / null-p)` rows to `populate_database.py` (M11), re-run populate.

### Block 1 + KB — DONE & REPRODUCIBLE (2026-06-24)
Block 1 biology ran end-to-end on unitigs and the **KB is populated**. Commits `1e09219`→`853143f` on `main`.
- **NCBI remote BLAST fixed (the session's main fight).** The public NCBI server **kills `blastn-short` + `word_size 7` over nt with SIGXCPU** (CPU-usage limit) — short 7-base seeds explode across nt **even when restricted to one species** (tested, still SIGXCPU). So the remote pass is now **DECOUPLED from CARD** (`08_blast_pipeline.nf` + `08_blast_annotation.py`): CARD pass keeps `blastn-short`/word7 (local, fine); **NCBI pass uses `blastn` + `word_size 11` + `-entrez_query txid<taxid>[Organism:exp]` + `-max_target_seqs 50`**. The entrez_query is **taxid-based from the registry** (not the scientific name — a space breaks the Nextflow CLI launcher → "Illegal option --"). Result: CARD 3605 + NCBI 4522 hits (E. coli-restricted). `09` got an `AMR_ENTREZ_EMAIL` env override (no config.yaml edit).
- **Nextflow-under-nohup gotcha:** Nextflow's ANSI console does terminal ioctls; a backgrounded JVM gets **SIGTTOU and STOPS (state `T`)** before submitting processes. Fix baked in: 08 sets **`NXF_ANSI_LOG=false`** (+ strips `NXF_OPTS` to kill the benign "Illegal option --"). HPC background runs now unblocked.
- **M8 KB populated** (`results/ecoli/kb/amrk.db`, **schema 0.2.0**): 65 kmers, 60 model-scores, 60 blast_annotations, 60 background-freq, 11 SNP, 131 validation_evidence. `populate_database.py` prefers `10_kmer_background_frequency` (superset of `07_kb_candidates` — has all candidate cols + prevalence/fisher/discriminative).
- **M16 ARO in KB** (schema bumped 0.1.0→0.2.0): `blast_annotations` gained `aro_accession/aro_gene_family/aro_drug_class/aro_resistance_mechanism`; **13/60 ARO-mapped** (= the 13 confirmed CARD hits: CMY-198 `ARO:3008132`, OXA-1042, CTX-M-260/278, TEM-258 `ARO:3009077`, sul1).
- **M6 CARD version recorded** via `AMR_CARD_VERSION` env override → `kb_metadata.card_version = 4.0.1`.
- **Validation metrics (M7/H2):** 13 confirmed / 47 none; **recovery 32% → H2 FALSE** (<40%); novel fraction 68% (19 stable-novel unitigs → H4). H2 failing is OK per ROADMAP §0.4 reframe (KB value = PFER-bound + lineage-validated + transparent + novel, not recovery %); revisit under M4 CPSS+SHAP.
- **Env overrides added this session:** `AMR_ENTREZ_EMAIL`, `AMR_CARD_VERSION` (both bypass config.yaml). To repro Block 1: `export AMR_FEATURE_REPR=unitig AMR_ENTREZ_EMAIL=… NXF_ANSI_LOG=false` then `08_blast_annotation.py` (UI node, internet) → `09` → `populate_database.py` (`AMR_CARD_VERSION=4.0.1`).
- **Known KB gaps (→ M10):** `run_metadata.json` MISSING (04 HPO was cut short) → `pipeline_runs` has no git_commit/timestamp (run_id `…__unknown`); `blast_annotations.coverage`/`description` NULL (not emitted to candidate CSV); `delta_prevalence` not in CSV.

### RESULTS (the thesis headline)
- **Unitig matrix:** 4,938,938 unitigs × 4373 genomes (median len 34 bp, min 31, p90 54, max 10030), 22 chunks, ~8 GB (vs 21.4M k-mers / 49 GB).
- **06 chunk-split test (unitig):** ROC-AUC **0.9534** (CI 0.936–0.969), MCC **0.8185**, bal-acc 0.913 — **beats k-mer baseline (0.930 / 0.739)**.
- **07b lineage-aware 5-fold GroupKFold (HONEST headline):** ROC-AUC **0.9505 ± 0.0102** (folds 0.934/0.943/0.961/0.957/0.957), **28 stable unitigs** (freq≥0.6), Jaccard 0.236. Lineage-CV ≈ chunk-split → **near-zero lineage leakage → signal is mechanism-driven, generalises across lineages** (answers the Yu-2024 reviewer-blocker).
- **05 early stopping:** best tree count = **8** (lr 0.024, depth 10). `n_estimators: 8` is in the experiment config.
- **Biology (08 fixed):** candidate unitigs map full-length (cov=1, 100% id) to β-lactamases **TEM-258/257, CTX-M-260/278, OXA-1042, CMY-198** (+ catA1, sul1, aadA24) — the real ampicillin mechanisms. rank-1 (32 bp, best 14 bp) is likely novel (H4).

### Containers on TRUBA (`$AMR_WORK/containers/`) — all built
- `amr.sif` — core + `unitig-caller 1.3.2 (Bifrost)` + `bcalm`. (NOTE: rebuild pulled **Nextflow 26.04.4**, strict parser — see 08 fix below.)
- `amr-pp.sif` — core + `poppunk 2.7.8` (env pins `setuptools<81`; PopPUNK needs legacy `pkg_resources`).
- `amr-tools.sif` — `pyseer 1.4.1 + quast 5.3.0 + ncbi-amrfinderplus 4.2.7 + resfinder` (M13/M14). DBs NOT baked: `amrfinder -u`, ResFinder DB clone — download to scratch when needed.
- `amr-checkm2.sif` — `checkm2 1.1.0` (separate: it pins python<3.9). DB: `checkm2 database --download --path $AMR_WORK/data/external/checkm2_db`.
- Build on a **debug node**: `unset APPTAINER_BINDPATH` (build sandbox can't bind `/arf`); login node hits a CPU-time ulimit on `mksquashfs`. `def`s defined in `amr.def`/`amr-pp.def`/`amr-tools.def`/`amr-checkm2.def`.

### Scripts / wiring (all on `main`, HEAD ~`400657f`)
- **`03u_unitig_matrix.py`** — `unitig-caller --call --rtab` → 03's exact chunk contract → `matrix_unitig/`. `--build-db` = organism-level store (`processed/{org}/unitig_all/`); per-antibiotic then SUBSETS it (no re-run). Config `unitig:` (`out_subdir`, `min_support: 10`, `db_min_support: 2`, `threads`).
- **`02c_lineage_poppunk.py`** + **`lib/lineage.py`** — PopPUNK **dbscan** (bgmm degenerate→refine NaN) → `processed/{org}/lineage/poppunk_clusters.csv`. `group_kfold_masks` → StratifiedGroupKFold sample masks. Config `lineage:` (`model: dbscan`, `refine: false`, `n_splits: 5`).
- **`07b`** — `build_cv_splits()`: lineage GroupKFold if `poppunk_clusters.csv` exists, else 5-seed fallback. n_estimators (=8) read from experiment config best_params.
- **`08`/`08_blast_pipeline.nf`** — FIXED: (a) `def OUTFMT`→`params.outfmt` (strict Nextflow 26 rejected top-level `def`); (b) BLAST task by **median query length** (`blastn-short` if median<50, with **word_size 7** — word_size 11 / `blastn` truncated short-unitig hits to ~14 bp noise). `blast.task` overrides.
- **`09`** — ARO mapping (M16): `aro_from_sseqid` + `load_aro_index` (`data/external/card/aro_index.tsv`) → KB cols `aro_accession/aro_gene_family/aro_drug_class/aro_resistance_mechanism`; coverage now = aln/`qlen` (08 emits qlen); reification note (S10).
- **`10`** — BH-FDR (§0.2): `fisher_q` + `discriminative_fdr`.
- **ENV OVERRIDES (no config edit on TRUBA):** `AMR_FEATURE_REPR=unitig` (matrix_dir→matrix_unitig), `AMR_EXTERNAL_MEMORY=false` (in-core 05/07b — faster but ~211 GB RAM for 4.9M feats), `AMR_OPTUNA_PATIENCE=15` (04 early stop), `AMR_POPPUNK_BIN`, `AMR_UNITIG_CALLER_BIN`.

### TRUBA artefacts present (verified)
matrix_unitig (features.txt + 22 X chunks + y + genomes + 70 GB unitigs.rtab) · lineage/poppunk_clusters.csv (324) · models/ecoli/ampicillin/xgboost_ampicillin_final_v2.json (8-tree unitig model) · config/experiments/ecoli/config_ampicillin.yaml (MANUAL — built from Optuna trial 20 after HPO was cut short; `n_estimators: 8`) · full CARD at data/external/card/ (card.json, aro_index.tsv, variant_model fasta) · CARD homolog DB at data/external/blast_db/card_nt/.

### Manual-config caveat (why)
The unitig 04 HPO (50 trials, ~5 h on 4.9M feats, slow) was **cancelled at trial 35**; best = **trial 20** (subset AUC 0.9903). `config_ampicillin.yaml` was written BY HAND (snippet) with trial-20 params + a manual_linspace chunk split (test = parts 0,7,14,21). 04 itself was NOT completed for unitigs. To redo properly later: re-run 04 with `AMR_OPTUNA_PATIENCE=15` (now coded) so it stops ~trial 35 cleanly.

### TRUBA rules that bit us
- **Compute nodes have NO outbound internet** → `blastn -remote` (NCBI) + Entrez FAIL there (errorStrategy 'ignore' keeps CARD alive). Run **08 + 09 on the UI node (arf-ui1)** for NCBI/Entrez; 10/11 on compute (offline OK). Local nt (~200 GB, N6) is the alternative.
- **barbun = whole-node feel:** `-c40 --mem 300G` waits for a fully-free node + fairshare (heavy daily usage lowers priority); `-c20 --mem 100G` (half node) schedules instantly. In-core 05 used **211 GB** (4.9M feats) — request ≥256 G if `AMR_EXTERNAL_MEMORY=false`.
- Pull code with **targeted `git checkout origin/main -- <file>`**; NEVER `config.yaml` (manual HPC tuning) or `02p`/`02b`/`03` (manual parallel patches). The TRUBA `config.yaml` has hand-added `unitig:`/`lineage:` sections, but env overrides make config edits unnecessary.

## §0.2 SLURM templates (copy-paste; submit from `$AMR_WORK`)

Common header rules: submit from `/arf/scratch` (`cd $AMR_WORK && sbatch ...`); `export APPTAINER_BINDPATH=/arf`; env overrides set the behaviour (no config edit); `2>&1` so Optuna/Nextflow stderr lands in the `.out`. `H=/arf/home/edemirbas/ML_AMR_Prediction_v2`, `SIF=$AMR_WORK/containers/amr.sif`. **Half node (`-c20 --mem 100G`) schedules instantly; full node (`-c40 --mem 300G`) waits.** Internet-needing steps (08-NCBI, 09-Entrez) run on the **UI node** (not SLURM).

**A) Unitig ML chain 04→05→06→07b** (`-c40 --mem 300G` in-core, or `-c20 --mem 120G` + `AMR_EXTERNAL_MEMORY=true` ExtMem):
```bash
#SBATCH -J amr-ml -p barbun -N1 -c20 --mem=120G --time=1-00:00:00 -o amr-ml-%j.out -e amr-ml-%j.err
set -euo pipefail; export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig AMR_EXTERNAL_MEMORY=true AMR_OPTUNA_PATIENCE=15
cd $H; for s in 04_optimization 05_model_training 06_evaluation 07b_feature_stability; do apptainer exec $SIF python -u scripts/$s.py 2>&1; done
```
**B) 07b only (in-core, fast)** — `-c40 --mem 300G`, `AMR_EXTERNAL_MEMORY=false`; needs the experiment config from 04.
**C) Biology 10+11 (compute, offline)** — `-c20 --mem 100G`, `AMR_FEATURE_REPR=unitig`; run `10_kmer_background_frequency` + `11_variant_snp_check`.
**D) Biology 07+08+09 on the UI node (interactive, internet for NCBI/Entrez):** `export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig; cd $H; apptainer exec $SIF python -u scripts/08_blast_annotation.py` then `09_biological_summary.py`.
**E) Organism unitig store (once/org, long)** — `-c40 --mem 300G`: `apptainer exec $SIF python -u scripts/03u_unitig_matrix.py --build-db --db-min-support 10`.
**F) PopPUNK lineage (once/org)** — `-c20 --mem 120G`, `SIF=amr-pp.sif`: `python -u scripts/02c_lineage_poppunk.py` (dbscan; `--reuse-db` to skip re-sketch).
**G) Container build (debug node)**: `unset APPTAINER_BINDPATH; export APPTAINER_TMPDIR=/tmp/apptmp APPTAINER_CACHEDIR=/tmp/apcache; cd $H; apptainer build --fakeroot $AMR_WORK/containers/<name>.sif <name>.def`.
**H) Block-2 tool runs** (`SIF=amr-tools.sif`/`amr-checkm2.sif`): pyseer LMM (M14) ✓, **AMRFinderPlus/ResFinder concordance (M13) ✓** (`16_external_concordance.py`; amrfinder DB via `amrfinder_update -d $AMR_WORK/data/external/amrfinder_db`, ResFinder DBs cloned to `$AMR_WORK/data/external/{resfinder,pointfinder}_db`; per-genome loop `xargs -P20` inside one `amr-tools.sif` exec over `16_genomes.txt`, then `--mode post` in `amr.sif`). **CheckM2+QUAST QC (M15) ✓ DONE (2026-07-02)** (`02d_genome_qc.py`): CheckM2 DB at `$AMR_WORK/data/external/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd` (download's json-write is read-only-FS-blocked but the .dmnd lands fine → pass `--database_path`); SLURM chains prep(amr.sif) → `checkm2 predict --database_path $DMND --tmpdir` (amr-checkm2.sif) → `quast.py *.fna --no-plots/html/icarus` (amr-tools.sif) → post(amr.sif). **Result: 5312/5470 = 97.1% pass** (fail 158: N50 131, completeness 29, contigs 26, contamination 13). Outputs in `results/ecoli/global_exploration/genome_qc/`. Advisory (not retrained; fails <3%).

**M13 DONE — concordance (2026-07-02).** `16_external_concordance.py` (+`lib/concordance`): AMRFinderPlus 2026-05-15.1 + ResFinder 4.5.0 on **5468 genomes** vs EUCAST/CLSI + **leakage-free model-vs-tool head-to-head** on the config held-out test split (06 saves `16_model_preds_{ab}.csv`). **Tool vs phenotype:** ResFinder κ amp 0.86 / cef 0.89 / cip 0.76 (solid); AMRFinderPlus over-calls amp+cip (naive Class/Subclass→drug map catches intrinsic β-lactamase + efflux → ME ~38-40%; cef clean κ 0.87). **Head-to-head (held-out):** unitig model bACC amp 0.873 / cef 0.925 / cip 0.928 — matches ResFinder (cef,cip), beats AMRFinderPlus (cip,amp). Outputs: `results/ecoli/external_validation/16_concordance_{ecoli}.{csv,json}`. **Still open:** temporal/geo hold-out (needs BV-BRC year/geo metadata fetch — Track A).
**I) M9 permutation (DONE/pending).** MDA `12_permutation_test.py`: **compute node only** (dies on UI CPU ulimit), offline, `AMR_FEATURE_REPR=unitig`, ~minutes. Label-perm `12b_label_permutation_test.py`: **full node** `-c40 --mem 300G` (in-core QuantileDMatrix ~211G), offline; builds once then 100 fast refits. **Apptainer node gotcha:** some barbun nodes throw `FATAL: Couldn't determine user account information: user: lookup userid …` at container start (env, not code) — resubmit and/or `sbatch --exclude=<node>`, or `apptainer exec --no-home …`.

**Connection / layout (all on TRUBA):**
- Login: OpenVPN → `ssh edemirbas@172.16.6.11` (UI = `arf-ui1`; transfer hosts `arf-ui4/5` = .14/.15).
- `$AMR_HOME=/arf/home/edemirbas/ML_AMR_Prediction_v2` (git clone of `main` = **code**).
- `$AMR_WORK=/arf/scratch/edemirbas/amr` (**data/outputs/container**). `$SIF=$AMR_WORK/containers/amr.sif`.
- All three + `APPTAINER_BINDPATH=/arf` are in `~/.bashrc`; `~/.bash_profile` does `source ~/.bashrc`.
- In the repo, `data/ results/ logs/ runs/ models/` are **symlinks → `$AMR_WORK/…`**. CARD homolog DB copied to `$AMR_WORK/data/external/blast_db/card_nt/`.

**Environment = Apptainer container** (TRUBA forbids conda/pip on the shared FS):
- `apptainer` is at `/usr/bin/apptainer` (v1.3.6, **no `module load` needed**).
- `amr.sif` built from `$AMR_HOME/amr.def` (Bootstrap docker `condaforge/miniforge3` + `environment.yml`) on an **interactive debug node** with `apptainer build --fakeroot` (set `APPTAINER_TMPDIR/CACHEDIR=/tmp/...`). Contains python+xgboost+sklearn+biopython+certifi, KMC 3.2.4, BLAST+ 2.17, Nextflow.
- Run everything as `apptainer exec $SIF python scripts/...`.

**Hard TRUBA rules learned (critical):**
1. **Submit jobs from `/arf/scratch`** (`cd $AMR_WORK` before `sbatch`/`srun`) — else `srun: error: Lutfen islerinizi /arf/scratch/ ...`.
2. **`APPTAINER_BINDPATH=/arf`** required, else the container can't see scratch (symlinks → `mkdir`/FileNotFound errors).
3. **`ftp.bv-brc.org` is FIREWALL-BLOCKED** on TRUBA. Genome FASTAs are fetched from the **BV-BRC Data API** (`www.bv-brc.org/api/genome_sequence`, dna+fasta) — already the repo default in `00a` (commit `5d0c9a3`).
4. **Queues:** `barbun` **min 20 cores/node** (hamsi 28, orfoz 56). **MS-student limit = 40 cores.** → use `barbun -c 20` (or `-c 40` for the parallel ML step). `debug` ≤4h for tests.
5. **Low CPU efficiency → TRUBA warns (`Eff:%…`) and may auto-cancel + cut your core quota.** barbun's 20-core minimum means any single-threaded step looks ~5%. Mitigations applied & on `main`: **`scripts/02p_kmer_parallel.py`** (parallel KMC, 5470 genomes in ~2.5 min), **parallel 02b** spectra extraction (`1cd0119`), **parallel 03** per-genome KMC dump (`3ddc476`). **Caveat:** 03's per-genome *parse* (matching 5M k-mers/genome against the ~8 GB `kmer_to_index` dict) is **GIL-bound**, so thread parallelism only partly helps → 03 still ran ~5% and got warned. **03 is resume-safe** (skips existing `*.npz` + `features.txt`), so a kill is harmless — just resubmit. **True 03 fix (future):** 2-bit integer-encode k-mers + numpy `searchsorted` (drop the Python dict) or multiprocessing. **ML steps (04/05) use cores via `n_jobs=40`** → expected high efficiency.
6. Quotas (banner): `/arf/home` 100 GB / 100K inode; `/arf/scratch` 1 TB / 200K inode; **no backup**; `/arf` is NVMe Lustre (fast — no separate `/tmp` needed for KMC).

**`config.yaml` on TRUBA (NOT committed — TRUBA-specific tuning):** `kmc_mem:128`, `threads:20`, `n_jobs:40`, `chunk_size:200`, `n_trials:30`, `target_antibiotic:ampicillin`.

**SLURM scripts in `$AMR_HOME/slurm/`:** `00_test.slurm` (debug sanity ✓), `run_features.slurm` (02p→02b→03, barbun `-c 20 --mem 120G`, **done**), `run_matrix.slurm` (03-only spare), `run_ml.slurm` (04→05→06, `-c 40 --mem 300G --time 3-00:00:00`, **current**). All use `export APPTAINER_BINDPATH=/arf`, `set -euo pipefail`, mail to `eren0demirbas@gmail.com`, submit from `$AMR_WORK`.

**Data state:** 5470 genomes downloaded (395 had no API sequence, dropped). `amr_phenotypes.csv` = 5470×72. ampicillin 4446 tested (recommended target). KMC dbs for all 5470 in `$AMR_WORK/data/interim/ecoli/kmc_outputs/`.

**Commits made during deployment (all on `main`):** `bvbrc` NaN-safe filters (`17b1301`); `00a` FTP→API download + new `02p` parallel KMC (`5d0c9a3`); `02b` str-path KMC check (`c352c48`); **parallel 02b spectra (`1cd0119`)**; **parallel 03 dump (`3ddc476`)**. **Do NOT `git pull` on TRUBA** — its working copy is manually patched (02p, parallel 02b/03) + TRUBA-specific `config.yaml`; pull would conflict. `main` already contains all the code fixes.

**Immediate next steps (new session):**
1. Check `run_ml.slurm` (`amr-ml`): `squeue -u $USER`; `tail $AMR_WORK/amr-ml-*.out`. 04 HPO (30 Optuna trials × out-of-core over 22 chunks × 50.8M features) is the long phase. On success: `config/experiments/ecoli/config_ampicillin.yaml` + model + `06_evaluation` metrics appear.
2. If `run_ml` is killed for low efficiency: 04 is **not** resume-safe (restarts HPO). Re-submit; if it keeps getting flagged, note that 04/05's per-chunk DMatrix *load* is serial between trees — acceptable for a one-off, or reduce scope.
3. After ML → **biology job** (`barbun -c 20`): `07b → 07 → 09 → 10`. `08` CARD-local BLAST works on compute nodes; `08` NCBI-remote + `09` Entrez need internet → run those on the **UI** (compute nodes may lack outbound internet) or skip.
4. Then download results: `tar` `results/ models/ runs/` from scratch → home or `rsync` to laptop (scratch auto-purges in 30 days).
5. Token: pushes use a fine-grained PAT pasted in chat (expires fast; ask for a fresh one with Contents:write). Full step-by-step + real-run corrections in `docs/TRUBA_Kurulum_ve_Calistirma_Rehberi.md` (Appendix).

---

# 1. Project Overview — what & why

- **Goal.** Predict antimicrobial resistance (AMR) in *E. coli* from whole-genome assemblies using **alignment-free k-mer features + out-of-core XGBoost**, then *reverse-translate* the most important k-mers back into biology (BLAST vs CARD/NCBI). The thesis-level goal is a **queryable, confidence-tiered, cross-antibiotic AMR k-mer Knowledge Base (AMRK-DB)**.
- **The real research gap (why this is novel).** Many papers predict AMR from k-mers; almost none turn the ML feature-importance output into a *reproducible, stability-filtered, biologically validated, discriminativeness-checked* knowledge base. The contribution chain is: **Gain → seed stability → BLAST confidence tier → discriminativeness (R vs S) → SNP-allele check → cross-antibiotic overlap → KB**. See `docs/ROADMAP.md` for the full thesis framing (hypotheses H1–H4, must-haves M1–M11).
- **Current stage.** Pipeline is implemented end-to-end (`00a → 11` + `07b`), multi-organism/antibiotic, audit-clean, **cross-environment (macOS/Linux/HPC)**, and **verified on the real 1788-genome dataset**. The KB persistence layer (SQLite/Postgres + API, M8/M10/M11) and the cross-antibiotic hypergeometric test (S1) are **not yet built** — those are the next big items.

---

# 2. Repository Structure

```
config/
  config.yaml                       # global config: organism, target_antibiotic, params, tiers, tool/data paths
  registry/organisms.yaml           # organism -> taxid, data paths, antibiotic set (single source of truth)
  registry/antibiotics.yaml         # antibiotic classes + ALIASES (name normalisation single source)
  experiments/{organism}/config_{antibiotic}.yaml   # AUTO-generated by step 04 (data split + best HPO params)
scripts/
  00a_download_bvbrc.py             # BV-BRC download (API/CLI/--raw-csv) + clean + parallel {id}.fna fetch
  00_prepare_metadata.py            # cleaned long table -> wide binary amr_phenotypes.csv (∩ present .fna)
  01_data_validation.py / 01b_*     # phenotype validation + EDA plots; ML-target recommendation
  02_kmer_extraction.py             # KMC k-mer counting (k=21) -> per-genome .kmc_pre/.kmc_suf
  02b_global_qc_analysis.py         # global QC: complexity outliers (IQR) + min_support elbow advisory
  03_matrix_construction.py         # global vocab (KMC) -> sparse binary CSR .npz chunks (+ y, genomes, features.txt)
  03u_unitig_matrix.py              # ROADMAP §0 M12: unitig-caller rtab -> SAME chunked matrix (-> matrix_unitig/); replaces raw k-mers downstream
  03b_matrix_validation_qc.py       # matrix QC (sparsity/prevalence)
  04_optimization.py                # Optuna HPO -> experiment config + run_metadata.json
  05_model_training.py              # full-data boosting over streaming QuantileDMatrix -> model + manifest.json + threshold
  lib/xgb_data.py                   # ChunkDMatrixIter + build_quantile_dmatrix (streaming DMatrix; used by 05/07b)
  06_evaluation.py                  # metrics, ROC/PR, calibration, bootstrap 95% CIs, error analysis
  07b_feature_stability.py          # 5-seed repeated holdout: AUC mean±std, selection freq, Jaccard, stable set
  07_explainability.py              # Gain top-N ∪ 07b stable set -> CSV + FASTA (flagged)
  08_blast_annotation.py / .nf      # Nextflow: CARD local + NCBI remote BLAST (blastn-short)
  09_biological_summary.py          # tiered report + KB candidates + recovery/composite/novel metrics
  10_kmer_background_frequency.py   # R-vs-S prevalence + Fisher exact -> discriminativeness
  11_variant_snp_check.py           # CARD variant-model SNP allele check (resistant vs wildtype)
  12_permutation_test.py            # M9: MDA permutation importance (model fixed; per-candidate AUC drop + BH-FDR)
  12b_label_permutation_test.py     # M9: label-permutation null (shuffle labels, retrain frozen-HP -> null AUC dist)
  13_stability_selection.py         # M4: CPSS (Chi² prefilter + B=100 complementary pairs, π≥0.6, PFER) + TreeSHAP
  13b_stable_annotation.py          # M4: BLAST stable set vs CARD, reuse 09 tier+ARO -> KB-ready stable candidates
  14_pyseer_lmm.py                  # M14: pyseer LMM (--mode prep|post; pyseer CLIs run in amr-tools.sif via SLURM)
  kb_app.py                         # Streamlit KB explorer (local queryable UI; pip install streamlit, reads amrk.db)
  populate_database.py              # M8: load all step outputs -> results/{org}/kb/amrk.db (multi-antibiotic, unitig schema 0.4.0)
  lib/                              # shared package (config, registry, chunking, io_utils, run_metadata, bvbrc)
  constants.py, utils.py            # thin backward-compat shims -> lib/ (kept for old imports/tests)
  run_pipeline.py                   # orchestrator: runs the numbered steps in order (subprocess + logging)
  lib/logging_utils.py              # standard timestamped logger factory (orchestrator + new code)
  migrate_to_organism_layout.py     # reversible data-layout migration (already applied)
tests/                              # pytest: smoke / unit / integration + README
docs/
  TECHNICAL_REVIEW.md               # consolidated audit findings + resolution status
  SCALE_MLOPS_PLAN.md               # multi-organism + KB + MLOps plan
  ROADMAP.md                        # thesis roadmap (hypotheses, M1-M11, 6-month plan)
data/  models/  results/  logs/  runs/   # generated — only the CARD homolog BLAST DB is version-controlled
# Research-software-engineering scaffolding:
LICENSE (MIT) · CITATION.cff · pyproject.toml (PEP 621 + ruff/mypy/pytest config)
.github/workflows/ci.yml (ruff + unit/smoke on py3.10-3.12) · .pre-commit-config.yaml
Makefile · CONTRIBUTING.md · CHANGELOG.md
requirements.txt, environment.yml, pytest.ini, README.md, QUICKSTART.md, METHODOLOGY.md
```

---

# 3. Pipeline — step by step (what it does, why, how)

Run order (config-driven; each reads `config.yaml` for organism/antibiotic):

```
00a  download + clean BV-BRC AMR  ─►  00  binary phenotype matrix
01   validate / pick target           02  KMC k-mers     02b  global QC
03   sparse matrix chunks             03b  matrix QC
04   Optuna HPO  ─►  05  train  ─►  06  evaluate
07b  5-seed stability  ─►  07  candidate k-mers (gain ∪ stable)
08   BLAST (CARD + NCBI)  ─►  09  tiered report + KB candidates + metrics
10   discriminativeness (R vs S)       11  variant-model SNP allele check
```

- **00a `download_bvbrc`** — *why:* get a reproducible, cleaned AMR label set + the matching assemblies. *How:* fetches the BV-BRC `genome_amr` table (`--backend api` default, `cli` via `p3-*`, or `--raw-csv` from the website), cleans it via `lib/bvbrc.py` (EUCAST/CLSI only, Lab-Method evidence, R/S→1/0, antibiotic-name normalisation, duplicate-conflict resolution), then downloads each surviving genome as `{genome_id}.fna` in parallel (retry/resume/`--retry-failed`/`--max-genomes`). Writes `amr_cleaned_long.csv`, `download_manifest.json`, logs + reports.
- **00 `prepare_metadata`** — *why:* the numbered pipeline needs a wide binary label matrix. *How:* pivots the cleaned long table to `amr_phenotypes.csv` (`Genome ID` + one 0/1 column per antibiotic, blank = untested), intersected with the `.fna` files actually present. Genomes with AMR labels but no downloadable assembly are dropped (a few % is normal, e.g. `562.1`).
- **01 / 01b `data_validation`** — class balance, missingness, EDA plots, and a scientific ML-target recommendation (minority count/ratio per antibiotic). Uses antibiotic classes from the registry.
- **02 `kmer_extraction`** — KMC counts canonical 21-mers per genome (`min_count=1`); outputs binary KMC DBs. Re-runnable (skips genomes already counted).
- **02b `global_qc_analysis`** — scans all KMC DBs for genome-complexity outliers (IQR on unique-k-mer count) and computes a `min_support` "elbow" advisory (does **not** change config).
- **03 `matrix_construction`** — builds the **global k-mer vocabulary** with one KMC pass over all genomes (`-ci min_support` rare filter, `-cx max_support` drops core-genome k-mers), dumps it to `features.txt`, then writes the genome×k-mer presence/absence matrix as **CSR `.npz` chunks** of `chunk_size=200` genomes + `y_{ab}.csv` + `genomes_{ab}.csv`.
- **04 `optimization`** — Optuna HPO (`n_trials=25`, `eval_metric=auc`). Stratified chunk split into train/test/optuna subsets. `colsample_bytree` searched on a **√p-anchored log range**; `n_estimators = best_iteration+1`; `base_score=0.5` pinned. Writes `config/experiments/{organism}/config_{antibiotic}.yaml` (the data split + best params) and `run_metadata.json`.
- **05 `model_training`** — **standard full-data gradient boosting** over a single, streaming **`ExtMemQuantileDMatrix`** (external memory): `lib/xgb_data.ChunkDMatrixIter` feeds the chunks to XGBoost one at a time and the quantised pages are **spilled to fast scratch** (`cache_prefix`), so the matrix never has to fit in RAM. (An in-core `QuantileDMatrix` of the full train set peaked **>400 GB** and OOM-killed the 384 GB node — external memory keeps RAM bounded to ~one page + histograms.) Trains the Optuna-tuned `n_estimators` trees on the whole training set (every tree sees all training rows), with a single **global `neg/pos` instance weight** for class imbalance, `base_score=0.5`. Saves the model + `manifest.json`, and an **operating threshold fixed at 0.5** (global weighting; no test-set tuning → no leakage). The disk cache (`models/.../_xgb_cache_train`) is removed after training. *Replaces the previous 1-tree-per-chunk incremental regime (weaker fit + very low HPC CPU efficiency).* **07b uses the same external-memory regime per seed.**
- **06 `evaluation`** — single stratified split; ROC/PR curves, calibration, confusion matrix, **bootstrap 95% CIs**, MCC/κ, error analysis. **Does not overwrite** the config threshold (leakage fix).
- **07b `feature_stability`** — runs **before** 07. 5 seeds `[42,123,777,1024,2025]`, stratified 80/20, **fixed HPO across seeds** (no per-seed retune → leakage-safe; "repeated holdout", Mahé 2018 resampling, not true k-fold because of the out-of-core constraint). Each seed trains with the **same full-data boosting regime as 05** — a streaming `QuantileDMatrix` built from the seed's train rows (sample-level `row_mask` over the chunks), so it stays out-of-core (one chunk in RAM at a time) while every tree sees the whole train split. Reports AUC mean±std, per-k-mer **selection_frequency** (`stable` if ≥ `stability_threshold=0.6`), and mean pairwise **Jaccard** of the top-N sets.
- **07 `explainability`** — extracts the single model's Gain top-N k-mers, **then merges in the 07b stable set** (k-mers reproducible across seeds but not in the gain top-N), flagging each row `in_gain_topN` / `stable` / `selection_frequency`. Emits the candidate CSV + FASTA. So BLAST/biology (08–11) covers **both** the gain and stability views.
- **08 `blast_annotation`** — Nextflow runs two parallel BLASTs of the candidate FASTA: **CARD local** and **NCBI nt remote**, both with `-task blastn-short -dust no` (correct for 21-mer queries). Records the CARD DB version (`blastdbcmd -info`).
- **09 `biological_summary`** — grades every hit into **confirmed / candidate / weak / none** using **identity + coverage** (alignment length / k) as the primary, DB-size-independent criteria, with E-value a loose secondary gate (E-value is *not* comparable between CARD and NCBI). Joins 07b stability, resolves NCBI gene names via Entrez, and writes: the Markdown report (with quantitative summary + reification-fallacy + gyrA/SNP caveats), `07_kb_candidates_{ab}.csv` (per-k-mer KB record incl. **composite_score = stability × log10(1/E) × identity**), and `08_validation_metrics_{ab}.json` (**M7 known-mechanism recovery rate**, **H2 pass/fail ≥40%**, **H4 novel-candidate fraction**, tier counts).
- **10 `kmer_background_frequency`** — *why:* BLAST says *which gene*, not *whether the k-mer discriminates*. Streams the matrix once and computes each candidate's prevalence in **resistant vs susceptible** genomes + **Fisher's exact** + a `discriminative` flag (|Δprev| ≥ 0.10 AND p < 0.05). Flags BLAST-confirmed-but-ubiquitous k-mers (likely wildtype/lineage). Output `10_kmer_background_frequency_{ab}.csv`.
- **11 `variant_snp_check`** — *why:* a homolog hit to gyrA only proves "gyrA region present", not "resistance SNP present". BLASTs candidates against CARD's **protein-variant-model** sequences, parses each model's resistance SNPs from `card.json` (protein position, e.g. `S83L`), maps protein pos → CDS codon, reads the k-mer's strand-aware codon, translates, and classifies **resistant_allele / wildtype / other_variant / ambiguous**. Needs the full CARD download; **skips cleanly** with instructions if absent. Output `11_variant_snp_check_{ab}.csv`.

---

# 4. Last real-data run (2026-06-15, ampicillin, 1788 genomes)

Proof the whole chain works on real data, not just synthetic smoke:

| Step | Result |
|---|---|
| 00a / 00 | 1788 genomes; ampicillin 758 R / 828 S / 202 untested |
| 02 KMC | 1788/1788, 0 failures |
| 03 matrix | 1552 genomes × **30,082,953** k-mers, 8 chunks, `features.txt` ≈ 757 MB |
| 04 HPO | best ROC-AUC **0.935** (10-trial test run; config now back to 25) |
| 06 eval | **test ROC-AUC 0.862** (CI 0.81–0.91), PR-AUC 0.90, MCC 0.64, acc 0.82 |
| 07b | AUC **0.750 ± 0.051**, Jaccard 0.23, **5 stable** k-mers |
| 07 | 10 gain ∪ 3 added stable = 13 candidates |
| 08 BLAST | CARD 440 hits (after the blastn-short fix; was 23), NCBI 1.4 MB |
| 09 | tiers: 2 confirmed / 4 weak / 7 none; recovery 20%, novel 20% |
| 10 | **9/13 discriminative**; APH(6)-Id real marker (prev_R 0.53 vs prev_S 0.11, p=4e-76); **OXA-1238 confirmed-but-ubiquitous** (0.63 vs 0.53) |
| 11 | 0 resistant-allele SNPs (expected — ampicillin = β-lactamase **acquisition**, a homolog mechanism, not a point mutation; the SNP machinery's payoff is on ciprofloxacin/gyrA) |

---

# 5. Cross-environment design (macOS / Linux / HPC)

The project must run on the user's Mac **and** a remote HPC pulled from GitHub. Key mechanisms:

- **PATH-aware tool resolution** — `lib.config.resolve_tool(config_key, command, …)` finds KMC/BLAST in this order: **(1)** env override `AMR_<TOOL>_BIN`, **(2)** `shutil.which` (conda/module on PATH), **(3)** the bundled macOS binary under `bin/bin/` **only on Darwin** (it is a Mach-O arm64 build that would mis-fire on Linux). Used by **02, 02b, 03** (`kmc`/`kmc_tools`) and **11** (`blastn`/`makeblastdb`); **08** prepends the conda bin to PATH so Nextflow's `blastn` resolves. *Why:* the old hardcoded `bin/bin/kmc` path broke every Linux/HPC run.
- **HTTPS / SSL** — `00a` builds an SSL context from `certifi` for all HTTPS (BV-BRC API + FTP); conda Python otherwise fails cert verification.
- **No CWD assumptions** — every script anchors paths to `PROJECT_ROOT = Path(__file__).resolve().parent.parent` and imports `lib` because Python puts the script's dir on `sys.path[0]`; running `python scripts/XX.py` works from any directory.
- **`.gitignore` + data hygiene** — all generated data is ignored (`data/raw`, `data/interim`, `data/processed`, `data/external/**/metadata`, the full CARD bundle `data/external/card`, `*.npz`, `features.txt`, models, results, logs, runs, nextflow `work/`). The **only** committed data is the CARD homolog BLAST DB (`data/external/blast_db/card_nt/card.*`, ~8.5 MB) so step 08 works out-of-the-box. A fresh clone reproduces everything by running `00a → …`.
- **Dependencies** — `environment.yml` (conda, installs KMC/BLAST/Nextflow via bioconda) or `requirements.txt` (pip Python deps; tools installed separately). `certifi` is listed. `QUICKSTART.md` has the full fresh-machine + HPC (SLURM) setup.

---

# 6. Key config knobs (`config/config.yaml`)

- `project.organism` = `ecoli`, `project.target_antibiotic` = **`ampicillin`** (was gentamicin; changed because gentamicin is all-susceptible in the current sample).
- `preprocessing`: `k_length=21`, `min_support=5`, `chunk_size=200`, `kmc_mem`, `threads`.
- `training.n_trials=25` (was temporarily 10 for a fast test run; restored).
- `analysis`: `top_n_features=50`, `stability_threshold=0.6`, and **`confidence_tiers`** (identity+coverage+evalue per tier) + `report_max_evalue`.
- `blast`: `card_db_dir/name`, `evalue`, `word_size`, `threads`, and **`card_variant_fasta` / `card_json`** (step 11; full CARD download).
- `ncbi.entrez_email` — empty in config; set it before step 09 to avoid NCBI rate-limit warnings, **or** `export AMR_ENTREZ_EMAIL=…` (env override).
- Tool overrides via env: `AMR_KMC_BIN`, `AMR_KMC_TOOLS_BIN`, `AMR_BLASTN_BIN`, etc.

---

# 7. Important decisions (rationale)

1. **Conflict resolution** for duplicate (genome, antibiotic): majority vote → tie: newest `testing_standard_year` → still tied: drop (NaN) + log.
2. **Antibiotic normalisation** via registry aliases (`co-trimoxazole`→`trimethoprim/sulfamethoxazole`; canonical = registry spelling).
3. **Cleaning filters:** EUCAST/CLSI standards, Lab-Method evidence, R/S phenotypes only.
4. **HPO once, fixed across seeds** in 07b (leakage-safe repeated holdout, not k-fold — out-of-core constraint).
5. **Confidence tiers grade on identity + coverage**, not raw E-value, because E-value depends on DB size and is not comparable across CARD vs NCBI. **Weak hits are kept and flagged**, never silently dropped (so e.g. a gyrA-type partial hit stays visible).
6. **07 carries both** the gain top-N and the 07b stable set forward so biology covers both.
7. **Discriminativeness (step 10) is separate from BLAST** — a k-mer can hit a known ARG yet be non-discriminative (ubiquitous); both facts are recorded.
8. **Organism-scoping** of every path incl. the auto-generated experiment config.
9. **Generated outputs are not version-controlled**; only the CARD homolog DB is.
10. **FASTA filenames must be `{genome_id}.fna`** (pipeline globs `*.fna`, uses the stem as Genome ID).

---

# 8. Audit findings — all resolved (see `docs/TECHNICAL_REVIEW.md`)

- **P-01** data leakage (Youden's J fit on test set) — FIXED (threshold on train/val only; 06 doesn't overwrite config).
- **P-02** shell injection (`shell=True`) — FIXED (`lib.io_utils.run_command`, shlex, never shell=True).
- **P-03** `eval_metric` hardcode — FIXED (from config).
- **P-04** colsample √p mismatch — FIXED (dynamic √p-anchored range).
- R/S double-count (01), hardcoded `top_50` filenames (08/09), O(N²) Gram SVD (03b), Entrez email/api_key, PR-AUC inconsistency, KMC resume, NCBI errorStrategy, CARD-version record, duplicated code → `lib/` — ALL FIXED.
- Integration-test-caught runtime bugs: `base_score must be in (0,1)`, `n_estimators=0`, empty feature-importance crash, empty-stability table — ALL FIXED.

---

# 9. Recent work log (newest first)

- **2026-06-28/29 CANONICAL clean re-run + 2nd antibiotic; KB renamed to unitig (0.4.0).** Re-ran 04→populate from scratch (kept the unitig matrix) so everything is saved with provenance. Fixed on the way: `populate_run` key mismatch + dict-bind crash once run_metadata.json actually existed (git_commit_hash/random_seed/started_at/data_fingerprint.sha256); `12b` test-set misalignment with the new non-ascending config split (REAL AUC collapsed to 0.49 → build the test matrix in ascending chunk order); `13` CPSS PFER blew up to ~100 because the 66/146-tree model over-selected → added `--base-trees 10` (sparse base selector → PFER back to ~2.7). **Renamed the whole KB schema k-mer→unitig** (`unitigs`, `unitig_id`, `unitig_model_scores`, `unitig_background_frequency`, `unitig_antibiotic_overlap`, `n_unitigs`; schema 0.3.0→**0.4.0**); CSV column reads stay `kmer` (on-disk name), output filenames unchanged. Added `scripts/kb_app.py` (Streamlit KB explorer). **Ampicillin** canonical (lineage-CV 0.9511, H2 TRUE 47%, 303 evidence) + **ciprofloxacin** appended to the same KB (model_id=2; **SNP showcase: gyrA S83L / parC S80I resistant_allele** — the FQ positive control; CARD-homolog recovery ~0 as expected). KB now multi-antibiotic. Commits `dc35d99`→`1b7d004`. **Gotchas:** 12b is slow with the bigger tree count (N=50 standard now); barbun45 = bad Apptainer node (`--exclude`); `AMR_CARD_VERSION` typo on cipro fixed in-DB via UPDATE.

- **2026-06-25 M14 pyseer LMM done & in KB.** Genome-wide pyseer was impractical (`similarity_pyseer` over the full 66 GB Rtab ran 16 h at ~2.5% eff, OOM-approaching → cancelled). Switched to a **targeted** design: kinship from a genome-wide **subsample** (~80k unitigs), LMM on the **candidate** unitigs (5000 Chi²). Two-container SLURM (prep/post in `amr.sif`, `similarity_pyseer`+`pyseer --lmm` in `amr-tools.sif`; `amr-tools.sif` has no PyYAML). Bonferroni 1.09e-5: **26/39 CPSS-stable significant; 3/3 TEM β-lactamase pass, both aminoglycoside co-resistance genes do NOT** — LMM cleanly separates the ampicillin mechanism from plasmid-linked co-resistance. Loaded into KB as `pyseer_lmm` evidence (TEM p≈1e-109). Commits `eec562e`→`8343ebb`.
- **2026-06-25 M4 CPSS + SHAP done & in KB; M9 evidence wired in.** `13_stability_selection.py`: Chi² prefilter (top-5000) → CPSS B=100 (200 fits) → **39 stable (π≥0.6), PFER bound 5.4** + built-in TreeSHAP. CPSS picked a *different* set from gain-top50 (0 exact overlap) but the same genomic regions — redundancy made visible. `13b_stable_annotation.py`: BLAST stable set vs CARD, reuse 09's identity+coverage tiering (kills 14 bp noise) + ARO → **5 confirmed full-coverage: TEM-256/257/258, APH(6)-Id, AAC(6')-Ib7**. Loaded into KB (`populate_database.py`, schema **0.3.0**: `kmer_model_scores` += `mean_abs_shap`+`selection_method`; `gain_seed`=60 / `cpss`=39). Also wired step 12/12b permutation results into `validation_evidence` (60 MDA + 1 label-perm p=0.0099). Commits `2223eb1`→`c29f972`.
- **2026-06-24 M9 permutation tests (MDA + label-permutation) — see §0.1 "M9".** Commits `2306f61`→`7b69af4` on `main`. `12_permutation_test.py` (MDA): ran & saved, baseline 0.9534, 0/51 significant at BH-FDR Q<0.05 — expected under unitig redundancy (top-MDA = the real β-lactamases); motivates M4 CPSS. `12b_label_permutation_test.py` (label-perm null): REAL 0.9534 ≫ null ~0.50 → model highly significant (p≈1/(N+1)); rewritten to build the train QuantileDMatrix once + swap labels (set_label/set_weight) after the per-perm rebuild was too slow (low-eff killer) and a plain-CSR build OOM'd 120G. 12b saved (job 5970202): REAL 0.9534, null ~0.50, **p=0.0099 significant**. Gotchas logged: MDA dies on the UI CPU ulimit (run on compute); 12b needs a full node (-c40 --mem 300G, ~211G in-core); some barbun nodes throw an Apptainer `lookup userid` error (resubmit/exclude).
- **2026-06-24 Block 1 biology made reproducible + KB populated (M8/M16/M6) — see §0.1 "Block 1 + KB DONE".** Commits `1e09219`→`853143f` on `main`. Fixed NCBI remote BLAST (SIGXCPU): decoupled the remote pass from CARD → `blastn`/word11 + taxid `-entrez_query` + `-max_target_seqs`; taxid (not scientific name — space breaks Nextflow launcher); `NXF_ANSI_LOG=false` so backgrounded Nextflow doesn't stall (SIGTTOU); `AMR_ENTREZ_EMAIL`/`AMR_CARD_VERSION` env overrides. Re-ran 08 (CARD 3605 + NCBI 4522, E. coli) → 09 → `populate_database.py` → `amrk.db` schema **0.2.0** (M8). Added ARO ontology cols to `blast_annotations` (**M16**, 13/60 mapped) + recorded CARD 4.0.1 (**M6**). H2 still FALSE (recovery 32%) — fine per §0.4. **Next: Block 2, starting M9 MDA permutation.**
- **2026-06-23 FULL unitig pivot executed end-to-end (ROADMAP §0) — see §0.1/§0.2.** Unitig matrix (4.94M unitigs) + PopPUNK dbscan lineages (324) + lineage-aware 07b → **honest headline ROC-AUC 0.9505 ± 0.0102** (5-fold GroupKFold; ≈ chunk-split 0.9534 → near-zero lineage leakage; beats k-mer 0.930). 04 HPO cut short at trial 35 → experiment config written by hand from trial-20 params (`n_estimators: 8`). Biology fixes: Nextflow 26 strict-parser (`def`→`params.outfmt`), BLAST task by median length + **word_size 7** → candidate unitigs hit **TEM/CTX-M/OXA/CMY** β-lactamases at cov=1. Added 09 ARO mapping (M16), 10 BH-FDR, S10 reification, env overrides (`AMR_FEATURE_REPR`/`AMR_EXTERNAL_MEMORY`/`AMR_OPTUNA_PATIENCE`), 4 containers (amr/amr-pp/amr-tools/amr-checkm2). Commits `0909589`→`400657f` on `main`. **IN PROGRESS:** Block 1 biology (08/09 on UI for NCBI, 10/11 on compute).
- **2026-06-22 unitig pipeline — step 1 started (ROADMAP §0 / IMMEDIATE NEXT #1):** added `bcalm`+`unitig-caller` to `environment.yml`; wrote `03u_unitig_matrix.py`; added the config `unitig:` section. (Superseded by the 2026-06-23 entry above.)
- **2026-06-22 literature review (2 rounds) → ROADMAP §0 methodological pivot:** systematic review (Sections A–F + an implementation "how" round) concluded that the raw-k-mer + random-CV approach, while a working baseline, is not publication-grade. Binding decisions written to `docs/ROADMAP.md` **§0** (and must-have table M2/M4 revised + M12–M16 added): switch to **unitigs** (bcalm2/unitig-caller), **lineage-aware CV** (PopPUNK+GroupKFold), **CPSS stability selection + SHAP**, add **external validation + AMRFinderPlus/ResFinder concordance**, **pyseer LMM**, **CheckM2/QUAST QC**, **BH-FDR (step 10)**, **ARO/CARD ontology mapping**, reification-safe wording; novelty reframed away from "first ML AMR DB". Confirmed unchanged: k=21, binary+max_bin=2, class-weight/no-SMOTE, BV-BRC, AUC~0.93. **This pivot makes the day's raw-k-mer/min_support/GPU work largely moot (unitigs dissolve the speed/memory problem).** Next session implements §0.
- **2026-06-22 KB layer started (M8)** (`d49377f`, `54178a1`): `scripts/lib/kb_schema.py` (SQLite DDL, 11 tables per ROADMAP §1.1, stdlib `sqlite3` — no new dep / no container rebuild) + `scripts/populate_database.py` (loads run_metadata, manifest, 06 metrics, 07b holdout, 09/10 candidate+background, 11 SNP → `results/{org}/kb/amrk.db`; idempotent, multi-antibiotic, graceful on missing inputs; writes `validation_evidence` per result (M11) + `kb_metadata` FAIR row incl. CARD 4.0.1). Functionally tested on synthetic inputs. **Next:** run the pipeline to produce real outputs, then `populate_database.py` for real + validate column mapping; then FastAPI endpoints (S8) + cross-antibiotic overlap (S1) + permutation test (M9). **Also done 2026-06-22:** data-adaptive min_support (03), GPU evaluated+rejected (V100 16 GB), TRUBA scratch cleaned (~2.4 GB), full CARD 4.0.1 downloaded (step 11 now active for all antibiotics).
- **2026-06-22 prediction-engine consolidation/audit:** after many incremental patches, audited 04/05/07b/lib for consistency, reproducibility, no-hardcoding. Fixes: seeded the Optuna TPE sampler (`random_seed`) for reproducible HPO; moved the 05 early-stopping val fraction to config (`training.validation_fraction`, was hardcoded 0.15); 07b now respects `training.external_memory` (in-core vs ExtMem) like 05; manifest records the actual `n_trees` + corrected stale `threshold_type` (global neg/pos weight, not per-chunk "Dynamic Instance Weighting"); fixed stale "incremental/epoch" docstrings in 05/07b. No hardcoded abs paths remain in the engine. 57 unit/smoke + 1 integration pass. Config knobs now: `optuna_threads_per_trial`, `max_boost_rounds`, `external_memory`, `max_train_chunks` (all documented in `config.yaml`).
- **2026-06-18→22 full-data boosting refactor** (`a5b9ddc`→`7e5c22b`): replaced incremental 1-tree/chunk with standard boosting over a streaming `(Ext)QuantileDMatrix` (`lib/xgb_data.py`); parallel Optuna trials (04); ExtMem to survive the >400 GB in-core peak; early stopping in 05. Result on 4373 genomes: ROC-AUC 0.930 / MCC 0.739. Early stopping disproved the underfit hypothesis (best ≈ 29–30 trees) → metrics are at the k-mer signal ceiling; pivot to the KB contribution (see §0).
- **2026-06-15 RSE / open-science scaffolding:** added `LICENSE` (MIT), `CITATION.cff`, `pyproject.toml` (packaging metadata + ruff/mypy/pytest config), GitHub Actions CI (`ruff` + unit/smoke on py3.10–3.12), `.pre-commit-config.yaml`, `CONTRIBUTING.md`, `CHANGELOG.md`, a `Makefile`, a `run_pipeline.py` orchestrator, and `lib/logging_utils.py`. Added type hints to `lib/config.py`. Untracked the auto-generated experiment configs (now gitignored). Smoke test now covers steps 10/11 → **56 tests pass**.
- **2026-06-15 cross-env cleanup** (`647c0d9`): `resolve_tool` adopted in 03 + 02b (was hardcoded `bin/bin/kmc`); `.gitignore` hardened + generated artifacts untracked (only CARD homolog DB kept); README pipeline list updated to `00a→11`. Verified: all scripts py_compile (3.10), 54 unit/smoke + 1 integration test pass, no hardcoded/CWD-dependent paths.
- **2026-06-15 step 11** (`306a18f`): `11_variant_snp_check.py` — CARD variant-model SNP allele check (pure helpers unit-tested incl. ± strand; fixed a `card.json` parse bug where `param_value` entries are bare strings).
- **2026-06-15 step 10** (`541f90d`): `10_kmer_background_frequency.py` — R-vs-S prevalence + Fisher + discriminativeness.
- **2026-06-15 BLAST/tier fix** (`e7e2e62`): CARD `blastn-short -dust no` (CARD hits 23→440); tiers → identity+coverage; hits sorted by E-value; gyrA/SNP + cross-DB caveats.
- **2026-06-14 biology pass** (`84dc94d`): M7 recovery rate, composite score, H4 novel fraction in 09; tiers moved to config; 07 merges 07b stable set; order 07b→07.
- **2026-06-13 portability + BV-BRC fixes** (`64b0fab`, `ab68d02`): `resolve_tool`, certifi SSL, API URL-encode (HTTP 400 fix), dry-run AMR sampling, batched CLI fetch, `base_score=0.5` in 04, rewritten QUICKSTART.

---

# 10. Known issues / technical debt

- **`config ncbi.entrez_email` is empty** → step 09 Entrez warns / may rate-limit (does NOT crash — falls back to stitle parsing). Set it via `config.yaml` **or** `export AMR_ENTREZ_EMAIL=…` (env override, added 2026-06-24; no config edit needed on HPC). Same for `AMR_ENTREZ_API_KEY`.
- **Step 11 needs the full CARD download** (`data/external/card/card.json` + `nucleotide_fasta_protein_variant_model.fasta`); not shipped (ignored). Skips cleanly if absent.
- **BV-BRC API deep-pagination cap** on the full 243k-row table — `--backend api` works and is fast, but if truncated fall back to `--backend cli` (batched, slow) or website `--raw-csv`.
- **No pipeline orchestrator** (`run_pipeline.py`) — steps run manually.
- ~~**04 vs 05 training regimes differ** (single early-stopped fit vs 1-tree/chunk incremental)~~ — **RESOLVED 2026-06-18**: 05 (and 07b) now use the same standard full-data boosting as 04's HPO (streaming `QuantileDMatrix`, `lib/xgb_data.py`). 04 still tunes `n_estimators` on a representative chunk subset; 05 trains that many trees on the full set. Minor remaining nuance: the budget is tuned on a subset, not the full train set (acceptable, keeps HPO fast).
- **BV-BRC env footgun:** `source /Applications/BV-BRC.app/user-env.sh` shadows `python` with BV-BRC's Python 2.7 → run scripts with `python3` or the explicit conda path while it's sourced (so `p3-*` stay on PATH).
- **Only ampicillin** has been run on real data so far; cefotaxime/ciprofloxacin/gentamicin matrices not regenerated at full scale.

---

# 11. Next priority tasks

**High priority (thesis novelty)**
- `10_cross_antibiotic_analysis.py` + **hypergeometric / Fisher test** (S1 / H3) — uses 07b stable-k-mer sets across antibiotics (β-lactam-internal vs cross-class overlap).
- **Permutation test** (M9): null Gain distribution from shuffled labels vs the real model.
- Run the full chain on **ciprofloxacin** to exercise step 11 on a genuine SNP mechanism (gyrA/parC) and validate `resistant_allele` calls.

**KB layer (M8/M10/M11)**
- `populate_database.py` → SQLite (then Postgres) with tables incl. `kmer_background_frequency`, `validation_evidence`, `pipeline_runs`; `kb_schema_version` + Zenodo DOI; minimal FastAPI endpoints (`/kmers`, `/overlap`, `/metadata`).

**Infra / optional**
- `run_pipeline.py` orchestrator + Makefile; temporal validation (S2); MLST/phylogenetic bias flagging (S3).

---

# 12. Working conventions & notes for the next session

- **Branch** `fix/amr-audit-remediation` (NOT main); **Conventional Commits**. The user has asked that pushes be attributed to them only — **do not add a `Co-Authored-By: Claude` trailer** to commits that will be pushed.
- **Push** only when asked; uses a GitHub PAT pasted in chat (fine-grained, `Contents: Read and write`). Remind the user to revoke/rotate exposed tokens.
- **Verification discipline:** the assistant sandbox is Python 3.13 but the user runs **3.10** — always `py_compile` with `/opt/anaconda3/envs/bitirme_vol2/bin/python`. The user runs xgboost/KMC/BLAST/BV-BRC locally; **`pytest` (54 tests) + `pytest -m integration` (synthetic 02→07b) are the primary validators** and have caught several real runtime bugs.
- **Single sources of truth:** `lib/` (helpers), `config/registry/` (organisms + antibiotics + aliases), `config/config.yaml` (params/tiers/paths).
- **Naming:** organism slug `ecoli` ↔ taxid 562; antibiotic ids lowercase with `/` preserved; genome FASTA `{genome_id}.fna`; experiment config `config/experiments/{organism}/config_{antibiotic}.yaml`; run_id `{org}__{ab}__{UTC}__{git7}`.
- **Do not rewrite working code.** Keep changes incremental, organism-scoped, cross-environment, and tested. Read `docs/ROADMAP.md`, `docs/TECHNICAL_REVIEW.md`, `config/registry/*`, `lib/bvbrc.py`, and `lib/config.py` (esp. `resolve_path` / `resolve_tool`) before starting.
