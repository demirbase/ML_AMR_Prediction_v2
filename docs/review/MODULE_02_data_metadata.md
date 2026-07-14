# Modül 2 — Veri Toplama & Metadata

> Mercek: **Python'da kal + sertleştir** · **ESKAPEE ölçeklenebilirlik** · **su geçirmez, sıfır-hata, tam reprodüksiyon**
> İncelenen: `scripts/00a_download_bvbrc.py`, `scripts/00_prepare_metadata.py`, `scripts/lib/bvbrc.py`, `scripts/01_data_validation.py`, `scripts/01b_data_validation.py`, `tests/test_bvbrc.py`
> Tarih: 2026-07-13

---

## 1. Genel Değerlendirme

Bu modül ham veriyi (BV-BRC genom + AMR fenotip) çekip **modele girecek temiz binary fenotip matrisine** dönüştürür — pipeline'ın "girdi kalitesi" katmanı. Genel izlenim: **mimari olarak olgun ve ESKAPEE'ye hazır**. Organizma registry `taxid`'inden çözülüyor (yeni organizma = registry kaydı + koşu, kod değil — **M0/M1'in "registry entry, not code" iddiası burada DOĞRULANDI**). İndirme resume-safe + retry'li, provenance manifest'i var, cleaning saf/test-edilebilir (`lib/bvbrc.py`) ve `test_bvbrc.py` iyi kapsıyor. Audit'in "Genome-ID float hazard"ı `dtype=str` ile çözülmüş.

**Ama fenotip-etiket doğruluğunu etkileyen doğrulanmış bir bug ve iki veri-hijyeni açığı var.** Çakışma çözümündeki `np.argmax` NaN'lı yıllarda yanlış etiket üretiyor (verified); bilinmeyen antibiyotik adları matrise sızabiliyor; ve API sayfalama ortada hata verirse veri sessizce kırpılabiliyor. Bunlar "su geçirmez" için düzeltilmeli.

---

## 2. Güçlü Yanlar

- **Registry-güdümlü, ESKAPEE-hazır:** `00a`/`00`/`01`/`01b` organizmayı registry'den (taxid + resolve_path) çözüyor → A. baumannii vb. eklemek = registry kaydı + `00a --organism ...`. Tüm antibiyotikler indiriliyor (ML hedefi sonra seçiliyor — "hepsini indir" yaklaşımıyla uyumlu).
- **Provenance:** `download_manifest.json` (kaynak, filtre, tarih, sayımlar) + `cleaning_report.json` (adım-adım satır/çift sayıları) + ham tablo saklanıyor (`BVBRC_genome_amr.csv`) → re-clean edilebilir.
- **Dayanıklı indirme:** ThreadPool + retry + resume (var olan `.fna` atlanır) + `--retry-failed` + per-genome durum raporu. API/FTP fallback (FTP HPC'de firewall'lu → `www.bv-brc.org` API). certifi SSL (conda cert sorununa karşı).
- **Çok-backend:** BV-BRC CLI (p3-*), HTTP API, web-export CSV (`--raw-csv`) — büyük tablolarda deep-pagination limitine çözüm.
- **Saf, test-edilebilir cleaning** (`lib/bvbrc.py`, no network): evidence→standard(EUCAST/CLSI)→phenotype(R/S)→normalize→conflict-resolution. `test_bvbrc.py` alias/typo/cotrimoxazole/passthrough, filtreler, 3 çakışma senaryosu, CLI-header, pivot'u kapsıyor.
- **Genome-ID `dtype=str`** her yerde (00/retry/pivot) → "562.10"→float audit hazard'ı kapalı.

---

## 3. Problemler (önem sırasıyla)

### Critical
- (yok) — ölçekte veri-bozan bir hata yok; aşağıdaki etiket-bug'ı kenar-durum frekansında.

### High
- (yok)

### Medium
- **M2-1 — Çakışma tie-break bug'ı (DOĞRULANDI):** `lib/bvbrc.py:_resolve_group` satır 89: tie (eşit R/S) durumunda `np.argmax(yr.values)` en yeni yılı seçmesi beklenirken, yıllar **kısmen NaN** ise `np.argmax` NaN indeksini döndürüyor → yanlış satırın etiketi seçiliyor. Doğrulandı: `np.argmax([2015, NaN]) = 1` (NaN), etiket yanlış. **Fix:** `np.nanargmax` + labels'ı yıl-hizalı seç. Kenar-durum (tied conflict + kısmi yıl) ama gerçek bir yanlış eğitim-etiketi.
- **M2-2 — Bilinmeyen antibiyotik adları matrise sızıyor:** `clean_amr_table` (satır 145-146) `normalize_antibiotic` uygular; bilinmeyen adlar **aynen korunur** (sadece boş/None atılır). Yani "extended spectrum beta lactamase", "fluoroquinolones" gibi fenotip-etiketleri/çöp "antibiyotik" sütunu olur (HANDOFF'un işaretlediği). M1-L1'in buradaki tezahürü. **Fix:** cleaning'de `antibiotic_to_class(x) is None` olanları **rapora yaz (warn)** + opsiyonel `strict` ile at (yeni ilaçları kaybetmemek için varsayılan warn).
- **M2-3 — Sessiz API sayfalama kırpılması:** `fetch_amr_table` (satır 127-130) bir sayfa isteği ortada patlarsa `break` edip **kısmi frame** döndürüyor (hard-fail yok). 5 sayfanın 3.'sü düşerse sessizce eksik veri. Reprodüksiyon için tehlikeli. **Fix:** mid-fetch hatasını "tamamlandı"dan ayır; manifest'e `complete: false` yaz veya loud-fail.

### Low
- **M2-4 — Atılan Intermediate/NS/SDD sayımı raporda yok:** binary R/S için `Intermediate`/`Non-susceptible`/`SDD` sessizce düşüyor (tasarım). Methods şeffaflığı için `cleaning_report`'a "phenotype_dropped: {intermediate: n, ...}" eklenmeli.
- **M2-5 — 01/01b docstring/legacy path drifti:** docstring'ler `genome_amr_matrix.csv` (legacy) diyor ama kod `resolve_path('metadata_file')` = `amr_phenotypes.csv` okuyor; 01'de ölü legacy fallback path (`BASE_DIR/metadata/genome_amr_matrix.csv`). Doküman + ölü kol temizlenmeli (M1'deki legacy `paths:` ile birlikte).
- **M2-6 — registry `download_date` manuel:** `download_manifest.json` tarihi tutuyor ama `organisms.yaml.download_date` elle. `00a` çalışınca registry'ye yazılabilir (opsiyonel otomasyon).

---

## 4. Düzeltilmesi Gerekenler (madde madde)

1. **`_resolve_group` `np.nanargmax`'e geç** (M2-1) + `test_bvbrc`'ye kısmi-NaN-yıl tie testi ekle.
2. **Cleaning'de bilinmeyen antibiyotikleri raporla** (M2-2): `clean_amr_table` report'una `unknown_antibiotics: [...]`; opsiyonel `strict=True` ile at. `00a` bunu logla.
3. **API kırpılmasını loud yap** (M2-3): mid-fetch hatasında manifest `complete=false` + uyarı; `--raw-csv`'ye yönlendir.
4. **Cleaning report'a phenotype-drop sayıları** (M2-4).
5. **01/01b docstring + ölü legacy path temizliği** (M2-5).
6. (Ops.) `00a` registry `download_date`'i güncellesin (M2-6).

---

## 5. Refactor Önerileri

- **`clean_amr_table`'a `strict`/`report_unknown` parametresi** — normalize sonrası registry-bilinmeyenlerini görünür kıl (M1'in "strict normalize" önerisiyle birleşir; tek yerde, `lib/bvbrc`).
- **Fetch tamlık sözleşmesi:** `fetch_amr_table` `(df, complete: bool)` döndürsün; `main` `complete=False` ise manifest'e işaretlesin ve exit-nonzero opsiyonu.
- **Intermediate politikası tek switch:** `_PHENOTYPE_MAP`'i config'e taşı (drop vs I→R) → Methods kararı kod-dışı, tek yerden (bkz. §7 karar).
- **01/01b'yi tek "00c_data_report" altında birleştir** (opsiyonel) — validation + EDA aynı girdi; iki script yerine bir rapor adımı. (Şimdilik düşük öncelik.)

---

## 6. Bilimsel Eksikler (makale açısından)

- **Intermediate ('I') politikası belgelenmeli:** binary R/S için I'yı düşürmek yaygın ama bazı çalışmalar CLSI'ye göre I→R (veya I→S) yapar. Seçim + atılan sayılar Methods'ta olmalı (şu an sessiz).
- **Etiket köken/standart şeffaflığı:** EUCAST/CLSI-only filtresi + çakışma-çözüm (majority→newest-year→drop) iyi; ama çakışma/atılan sayıları makalede raporlanmalı (report'ta var, figüre/tabloya taşınmalı).
- **Çöp-antibiyotik sızıntısı** temizlenmezse fenotip matrisine sözde-hedef girer (bilimsel veri hijyeni).

---

## 7. Literatür/Karar Gereksinimi

Derin literatür GEREKMİYOR; tek bir **metodoloji kararı** senin onayını bekliyor (varsayım yapmıyorum):

- **Intermediate ('I') fenotipi ne olsun?** Şu an: **düşürülüyor** (sadece R/S). Alternatifler: I→R (CLSI temkinli yaklaşımı, klinik "işe yaramayabilir") veya I→S. Literatürde (ESKAPEE1.md §3: NCBI Antibiogram + CLSI/EUCAST) her ikisi de görülür; **binary R/S + I-drop savunulabilir varsayılan.** Onaylıyor musun yoksa I→R mi yapalım? (Kod tek switch'le hazırlanacak.)

Bunun dışında M2 düzeltmeleri (M2-1..6) literatür beklemiyor — onayınla uygulanır.

---

## Uygulama durumu (2026-07-13) — UYGULANDI

Kararlar: Intermediate → **config switch** (varsayılan drop); **tüm M2 düzeltmeleri uygulandı**.
- **M2-1 (bug):** `lib/bvbrc._resolve_group` → `np.nanargmax` + labels/years hizalama. Regresyon testi (`test_clean_conflict_tie_partial_nan_year`) eklendi ve **geçti** (argmax NaN'ı seçme hatası kapandı).
- **M2-2:** `clean_amr_table` → `report["unknown_antibiotics"]` + `n_unknown_antibiotic_names`; opsiyonel `strict_antibiotics` ile atma. Test eklendi.
- **M2-3:** `fetch_amr_table` → `(df, complete)`; `00a` mid-fetch kırpılmasında loud-warn + `manifest["fetch_complete"]`.
- **M2-4:** `clean_amr_table` → `report["phenotype_dropped"]` (atılan phenotype sayıları) + `intermediate_policy` rapora/manifeste.
- **Intermediate switch:** `config.yaml` → `metadata.intermediate_policy` (drop|resistant|susceptible); `00a` okuyup `clean_amr_table`'a geçiriyor + manifest'e damgalıyor. `drop`/`resistant` testleri eklendi.
- **M2-5:** 01b docstring + 01 ölü legacy fallback path `amr_phenotypes.csv`/organism-layout'a güncellendi.
- **M2-6** (registry download_date otomasyonu): ertelendi (opsiyonel, düşük değer).

**Test durumu:** `test_bvbrc` 12/12; tüm suite **101 passed**; M2 **0 yeni kırık** (9 ön-var kırık M6/M7/M10'da).

## Sonraki modüllere taşınan notlar
- **M3:** 01/01b'nin legacy fallback path'i + `paths:` legacy bloğu birlikte temizlenmeli; 02/02b/02c/02d/03 organism-aware mı (paralel ESKAPEE).
- **M6/M7/M10:** M1'de bulunan 9 ön-var test kırığı hâlâ açık.
- **Genel:** `_PHENOTYPE_MAP` config'e taşınırsa provenance'a "intermediate_policy" eklenmeli.
