# FİNAL AUDIT — TRUBA Deploy Öncesi (M0–M11 + ROADMAP sentezi)

> Amaç: geliştirme DEĞİL — tez + makale + production + TRUBA için hazırlık doğrulaması.
> Tarih: 2026-07-13 · Repo HEAD `7d3822f` (12 modül incelemesi push'lu) · Suite **109 passed / 0 kırık**

---

## A. KOD KALİTESİ

| Soru | Durum |
|---|---|
| Eksik refactor | **Kalanlar (düşük öncelik):** `scripts/lib` → düzgün `amr_lib` paketi (sys.path hack'leri, M0 §5) · 02/02b/03/03b import-time modül-globalleri (M3-6). İşlevsel değil, yayın-blokör değil. |
| İsimlendirme tutarlılığı | ✅ İyi. Tek not: `10_repeated_holdout_summary` dosya adı artık lineage-CV içeriyor (M5-4, populate ile birlikte rename → ertelendi). |
| Gereksiz/dead code | ✅ Temizlendi: legacy `paths:` bloğu (M3-5), `_ORG_META`/`_AWARE`/`_MECH` hardcode (M1), `08_blast_pipeline.nf` + nextflow dep (M9). **Kalan:** `amr-gpu.def` ölü (arşive → M9-4). |
| Kullanılmayan script | `02p_kmer_parallel` (parallel KMC alternatifi) — kullanılıyor mu netleştir; `01b` EDA opsiyonel. Silme, KEEP. |
| Duplicate | ✅ registry tek-kaynak (M1); `kb_queries` tek sorgu katmanı. |
| Yanlış/hardcoded path | ✅ `kb_api` bayat DB path düzeltildi (M8-1); `bin/bin/kmc` sadece Darwin fallback (resolve_tool PATH-önce); başka hardcoded yok. |
| Config doğru mu | ✅ organism-aware templates tek blok; `feature_repr=unitig` varsayılan; `intermediate_policy` switch; sürümler 0.6.1. |

**Sonuç:** Kod kalitesi **production-ready'e yakın.** Kalan refactor'lar kozmetik (amr_lib paketleme, modül-global), yayın/koşu blokörü değil.

## B. PIPELINE TUTARLILIĞI
Zincir: `00a→00→(01/01b)→02/02b/02c/02d→03u→04→05→06→07→07b→08→09→10→11→12/12b→13/13b→14→15→16→populate`.
- ✅ **Kırık bağlantı yok:** 03u→downstream aynı `.npz` sözleşmesi (features/y/genomes/X_part); env-parametrik hedef her adımda `get_target`.
- ✅ **Girdi/çıktı:** metadata `amr_phenotypes.csv` → matris → model → KB; dosya adları tutarlı.
- ⚠️ **k-mer baseline (03/03b)** korunuyor ama asıl yol unitig (03u); baseline artık env-parametrik (M3-2).
- ⚠️ **slurm/** orkestrasyonu git-dışı (M9-2, Drive'da) → deploy'da repoya alınmalı.

## C. BİLİMSEL AUDIT (ROADMAP uyumu)

**Must-Have (M1–M16):** M1(leakage-fix ✅ satır778), M2(lineage-CV ✅), M3(tier ✅), M4(CPSS ✅), M5(repro ✅, Zenodo hariç), M6(CARD ver ✅), M7(recovery ✅), M9(permütasyon ✅), M11(evidence ✅), M12(unitig ✅), M14(pyseer ✅), M15(genom QC ✅), M16(ARO ✅).

**GAP'LER (yayın öncesi karar gerektiren):**
1. **M10 — Zenodo DOI deposit YAPILMADI** (son must-have). FAIR "Findable" için şart. → deploy sonrası.
2. **M13/S2 — Temporal/coğrafi external hold-out YAPILMADI.** ROADMAP bunu **MUST** diyor ("zamansal veya coğrafi hold-out, hakemlerce kabul"). Şu an sadece araç-concordance var. `00a` yıl/coğrafya metadata'sını çekmiyor (SELECT'te sadece `testing_standard_year` var, collection year/geo yok). → **literatür + karar gerekli (§E-1).**
3. **S1 — H3 cross-antibiotic overlap istatistiksel testi ertelenmiş** (`--with-test`, union-universe caveat). ROADMAP: "bu olmadan H3 sadece betimsel." → **literatür + karar (§E-2).**
4. **M8 — PostgreSQL yerine SQLite.** ROADMAP kendisi "önce SQLite, production'da Postgres" diyor → tez için **kabul**, production/public-erişimde tekrar değerlendir.
5. **S3 — MLST klonal-kompleks bias analizi** açık (PopPUNK lineage-CV çekirdeği karşılıyor; ek MLST betimlemesi opsiyonel).

**Eksik validation/benchmark/istatistik:** temporal hold-out (2), H3 hypergeometric (3). Diğer tüm katmanlar (7-layer, concordance, permütasyon, pyseer) tam.

## D. PUBLICATION READINESS (reviewer gözüyle)

**REJECT riski:** yok (yöntem sağlam, leakage yok, lineage-CV var, concordance var).

**MAJOR REVISION riski:**
- **Temporal/geografik external validation eksikliği** — WGS-AMR external-validation için hakemler zamansal hold-out bekleyebilir (özellikle Database/Briefings). *En yüksek major-revision riski.* (§E-1)
- **Container + lock yok** → "fully reproducible" iddiası (Availability) savunulamaz (M9-3). Deploy'da lock + digest-pin şart.

**MINOR REVISION riski:**
- H3 "desteklendi" iddiası istatistiksel testsiz betimsel (S1) — hypergeometric ekle (§E-2).
- cefotaxime(Kp) head-to-head bACC≈0.495 (küçük n) — açıklama Methods'ta net olmalı (auc_mean_seeds 0.97).
- Public metadata derin prose hâlâ "E. coli" izleri (M11-2, deposit-finalize).
- provenance'a unitig-caller/pyseer sürümü (M7-6).
- Zenodo DOI (M10).

**Güçlü savunmalar:** flagship cipro(Kp) 0.926 vs araçlar 0.54 (SNP mekanizması) · cv_method ile "tüm AUC lineage-CV" kanıtı · 7-katman ortogonal kanıt · cross-organism concordance.

## E. LİTERATÜR ARAŞTIRMA İSTEKLERİ (Deep Research — varsayım yapmadım)

**E-1 — Temporal/coğrafi external validation: gerekli mi, nasıl?**
- *Neden:* En yüksek major-revision riski; ROADMAP MUST diyor ama yapılmadı. `00a` yıl/geo metadata çekmiyor.
- *Hangi kararı etkiler:* (a) submission öncesi temporal hold-out koşacak mıyız? (b) `00a`'ya BV-BRC `collection_year`/`geographic_location` alanları eklenip yeniden mi indirilecek? (c) train ≤2019 / test ≥2020 split mi, ülke-bazlı mı?
- *İncelenecek:* WGS-AMR ML'de temporal validation standardı (2023–2026); Database(Oxford)/Briefings in Bioinformatics AMR-ML makalelerinin external-validation beklentisi; BV-BRC'de collection-year kapsamı yeterli mi (ESKAPEE başına).

**E-2 — H3 cross-antibiotic overlap istatistiksel testi (union-universe sorunu).**
- *Neden:* H3 iddiası şu an betimsel; ROADMAP hypergeometric istiyor; "union-universe" caveat'ı var (evren tanımı testi belirler).
- *Hangi kararı etkiler:* hypergeometric mi Fisher mi; evren = tüm unitig mi, iki-antibiyotik-birleşimi mi; çoklu-karşılaştırma düzeltmesi.
- *İncelenecek:* pan-genom/unitig overlap significance testleri (DBGWAS/pyseer ekosistemi); ortak-özellik enrichment için doğru null evreni.

**E-3 — Yeni ESKAPEE organizmaları için PopPUNK ayarı (gram-pozitif + A. baumannii).**
- *Neden:* dbscan+no-refine E. coli/K. pneu'da valide; gram-pozitif S. aureus / A. baumannii popülasyon yapısı farklı (M3 notu).
- *Hangi kararı etkiler:* her yeni organizmada PopPUNK model/refine/çözünürlük (n_clusters≥n_splits) — lineage-CV'nin gerçek olması buna bağlı.
- *İncelenecek:* PopPUNK'ın S. aureus / A. baumannii / E. faecium'da önerilen ayarları ve benchmark'ları (2023–2026).

## F. TRUBA READINESS (KONSERVATİF — hiçbir şey silinmeyecek, sadece kategori)

> **Önce doğrula:** `du -sh $AMR_WORK/* $AMR_WORK/data/* $AMR_WORK/results/* 2>/dev/null | sort -rh | head -40` + `du -sh $AMR_WORK/data/processed/*/*/* 2>/dev/null | sort -rh | head` çıktısını paylaş → kesin liste ona göre.

**SAFE TO DELETE** (yeniden üretilebilir, tüm 21 model + pyseer bittikten SONRA):
- `data/processed/*/*/matrix_unitig/*unitigs.rtab` — pyseer-only, ~40-70GB/adet, regenerable (03u). *Tüm 14/populate bitti → güvenli.*
- `data/interim/*/kmc_outputs/` — KMC QC çıktısı, regenerable.
- `**/__pycache__`, `.pytest_cache`, `*.pyc`, `nohup.out`, `_retry_tmp.csv`, `_poppunk_work/` (ara; **poppunk_clusters.csv HARİÇ**).

**KEEP** (bilimsel değer / pahalı yeniden üretim):
- `results/kb/amrk.db` (KB — ama Mac 0.6.1 daha yeni; deploy'da re-populate) · `results/figures/` · `models/*/*/` (eğitilmiş) · `runs/` (provenance).
- `data/processed/*/lineage/poppunk_clusters.csv` (soy etiketleri) · `data/processed/*/unitig_all/` (organizma-store, pahalı) · `containers/*.sif` · `amr.def` · `slurm/` (kanonik).

**NEEDS BACKUP** (silmeden önce Drive'a):
- `results/kb/amrk.db`, `results/figures/`, `models/`, `runs/`, `poppunk_clusters.csv`, `16_concordance_*.csv`, `download_manifest.json`, `cleaning_report.json`, tüm `*_summary_*.json/csv`.
- (Zaten Drive'da: kod rescue + `amr_results_20260709_*.tar.gz`.)

**UNKNOWN** (elle bak, silme):
- `results/ecoli/kb/` gibi ESKİ per-organizma dizinleri (birleşik `results/kb` ile aşıldı — ama içerik doğrulanmadan silinmez).
- Eski experiment/log dizinleri (pre-0.6.0), eski conda env'leri (`conda env list` ile bak), `data/raw/*/genomes/*.fna` (Drive'da yedekli ama re-download pahalı → KEEP+BACKUP).

## G. DEPLOYMENT PLANI (uygulanmayacak — sadece plan)

1. **TRUBA temizliği** — §F survey çıktısını al → SAFE-TO-DELETE'i onayla → `lfs quota` boşluğu gör → sadece rtab/kmc/cache sil (KEEP'e dokunma).
2. **Backup** — NEEDS-BACKUP setini `rclone` ile Drive'a (arf-ui4, screen; whole-tree DEĞİL, alt-dizin) → `amr_predeploy_YYYYMMDD.tar.gz`.
3. **Repo senkronizasyonu** — `git -C $AMR_HOME fetch origin` → **`git reset --hard origin/main`** (untracked `amr.def`/`slurm/` etkilenmez; hazineler zaten Drive'da) → HPC config'i (`kmc_mem=128/threads=20`) tekrar uygula (ya da env) → `slurm/` kanonik 4'ü repoya commit et.
4. **Environment kurulumu** — `apptainer build --fakeroot amr.sif amr.def` (nextflow-suz environment.yml) → amr-tools/amr-checkm2 sif'leri doğrula → **lock üret** (`conda env export --no-builds > environment.lock.yml`) + commit (M9-3).
5. **Test çalıştırmaları** — `pytest -m "not integration and not slow"` (109 passed bekle) + `validate_registry.py` (0 hata) + bir smoke import.
6. **Küçük pilot run** — 1 antibiyotik, `--max-genomes 50` dry-run 00a → 03u `--min-support 1` → Faz1 → 07b (lineage-CV method logunu doğrula) → KB'ye populate → cv_method='lineage_group_kfold' teyit.
7. **Tam ölçekli antibiyotik çalışmaları** — ESKAPEE Faz-1 sırası: mevcut ecoli/kpneu re-populate (KB 0.6.1) → **S. aureus** (blocker çözülü, seed-42 shuffle) → **A. baumannii** (yeni). Her biri per-antibiyotik zincir (Faz1→Faz2a→bio→pyseer→populate→rm rtab), env-parametrik.
8. **Database oluşturulması** — populate + migrate (registry-meta + cv_method dolu) → `validate_registry --db` (0 hata) → figürler/tablolar yeniden üret.
9. **Son kalite kontrolü** — KB'de 21+ model, her modelde 7 katman + cv_method; `validate_registry` 0 hata; suite yeşil; Zenodo deposit (M10) + DOI damgası; docs prose-finalize (M11-2).

---

## ~~E1 KARARI — TEMPORAL VALIDATION YAPILACAK~~ → **TERSİNE DÖNDÜ (2026-07-15): TEMPORAL İMKÂNSIZ**

> **Bu bölümün kararı ÖLÇÜMLE ÇÜRÜTÜLDÜ. Aşağıdaki eski plan uygulanmayacak; yeni karar bunun altında.**

**Eski karar (2026-07-14):** "yayın öncesi ZORUNLU (TRIPOD+AI/DOME); BV-BRC Collection Year ~%85-89 dolu → fizibıl". Plan: `00a`'ya collection_year ekle, sızıntı-güvenli split (unitig sözlüğü sadece train'de), train ≤2021 / test 2023+.

**ÇÜRÜTME (2026-07-15, ölçüm).** E1 doluluk konusunda haklıydı (%73-85 ölçtük) ama **dağılıma hiç bakmadı.** Etiketlenebilir genomlarda (R/S + CLSI/EUCAST) `collection_year` ≥2023 olanlar:

| organizma | ≤2021 | **≥2023** |
|---|---|---|
| E. coli | 1362 | **28** |
| K. pneumoniae | 3726 | **13** |
| S. aureus | 2326 | **0** |
| A. baumannii | 991 | **11** |

**BV-BRC'nin AMR-etiketli genomları pratikte 2021'de bitiyor.** 0-28 genomluk bir dilim test seti değil, gürültü. Temporal hold-out bu veriyle **yapılamaz** — istemediğimizden değil, veri olmadığından. *Doluluk ≠ fizibilite:* E1'in atladığı ayrım buydu.

**YENİ KARAR — dış doğrulama = COĞRAFİ + SOY hold-out, yıl tamamen işlemden çıkar.**
ROADMAP zaten "zamansal **veya coğrafi** hold-out" diyordu; veri coğrafiyi çok daha iyi destekliyor: `isolation_country` **%96-99 dolu** (collection_year'ın %73-85'ine karşı, üstelik onda 1885/1800/1905 gibi çöp değerler var).

1. **`collection_year` / `testing_standard_year` KULLANILMIYOR.** 00a onları çekebilir (provenance), ama hiçbir split'e girmezler.
2. **Coğrafi hold-out** — `isolation_country` ile. UYARILAR: (a) **ülke dominansı** — E. coli'nin %58'i Norveç, A. baumannii'nin %63'ü USA; onları hold-out yapmak train/test'i ters çevirir. En dengeli adaylar **S. aureus** (UK 773 / Çin 665 / USA 659) ve K. pneumoniae (USA 2045 / Norveç 693 / İtalya 484). (b) **ülke ve soy karışık (confounded)** — tek bir çalışmadan gelen ülke, zaten bir soyu hold-out etmek demek; PopPUNK'la onu zaten yapıyoruz. Coğrafinin soy-CV'nin ÜSTÜNE ne kattığı **ölçülmeli**; katmıyorsa süslü bir tekrardır.
3. **Veri yetersizse coğrafiyi de katma** (kullanıcı kararı) — zorlama bir dış doğrulama, hiç olmamasından kötüdür.
4. **M13 concordance (AMRFinderPlus/ResFinder head-to-head) yerinde duruyor** ve şu an tek gerçek dış doğrulamamız.

**Not:** Bu ölçümdeki E. coli sayıları `cap=150000` ile budanmış olabilir (E. coli'nin lab satırı 243K+); coğrafi dağılım ORANLARI temsili ama mutlak sayılar için yeniden ölçülmeli. Temporal bulgusu bundan etkilenmez — 0/11/13/28 o kadar uç ki dört katı bile test seti etmez.

## ÖZET KARAR NOKTALARI
- **Deploy'a hazır** (kod/pipeline/test yeşil). ~~2 bilimsel gap~~ → **E-1 kapandı** (temporal imkânsız, ölçüldü; yerine coğrafi+soy hold-out — yukarı bak). Kalan: **H3 hypergeometric (E-2)**, submission-öncesi, deploy'u bloke etmez.
- **Deploy-must:** container lock + digest-pin (M9-3), slurm commit (M9-2), Zenodo (M10).
- **TRUBA:** sadece rtab/kmc/cache SAFE-TO-DELETE; gerisi KEEP/BACKUP; survey çıktısıyla kesinleştir.
