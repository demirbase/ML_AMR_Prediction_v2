# AMRK-DB — Bilgi Tabanı (Knowledge Base) Ayrıntılı Rehber

> **Dosya:** `results/kb/amrk.db` (SQLite) · **Şema sürümü:** `0.6.0` · **Lisans:** CC-BY-4.0
> **Kapsam:** 21 model (E. coli 7 + K. pneumoniae 14), 2 organizma, 7 ilaç sınıfı, 2363 unitig
> Bu belge KB'nin **her tablosunu, her kolonunu, ne işe yaradığını ve ne anlama geldiğini** tek tek açıklar. Danışman toplantısında gelebilecek soruların cevapları en sonda (§7).

---

## 1. KB Nedir, Neden Var?

Bu bilgi tabanı, antibiyotik direncini (AMR) **açıklanabilir biyobelirteçlerle** haritalayan, sorgulanabilir bir veritabanıdır. Amaç yüksek doğrulukla direnç **tahmini** değil (o zaten ~0.93 AUC ile çözülmüş); asıl katkı: **"hangi genomik imza hangi antibiyotiğe direnci sürüklüyor, ve bu ne kadar güvenilir kanıtlanabiliyor?"** sorusunu, her kayıt için **çok katmanlı kanıt zinciri + tam köken (provenance)** ile yanıtlayan **açık, FAIR-uyumlu bir AMR biyobelirteç bilgi tabanı** kurmak.

**Üç tasarım ilkesi:**
1. **Özellik birimi = unitig** (ham k-mer değil). Unitig = de Bruijn grafiğinin dallanmayan yollarının birleştirilmiş hali → daha az, daha uzun, **BLAST'lanabilir (gen-uzunluğu)**, GWAS-standardı özellikler. (Neden unitig? → §7.1)
2. **Her kayıt bir pipeline koşusuna bağlı** (`pipeline_runs`) → git commit + seed + CARD sürümü + config parmak izi ile **tam tekrarlanabilir**.
3. **Her biyobelirteç 7 ortogonal doğrulama katmanından** geçer (`validation_evidence`) → tek bir istatistiğe değil, çok katmanlı kanıta dayanır (reification/aşırı-yorum riskine karşı yapısal güvence).

---

## 2. Genel Mimari — 13 Tablo

Tablolar üç gruba ayrılır:

| Grup | Tablolar | Rol |
|---|---|---|
| **Referans / köken** | `pipeline_runs`, `organisms`, `antibiotics`, `models`, `kb_metadata` | Kim, ne, ne zaman, hangi sürümle üretti |
| **Çekirdek biyobelirteç** | `unitigs`, `unitig_model_scores`, `blast_annotations` | Unitig sözlüğü + model önemi + biyolojik kimlik |
| **Doğrulama / kanıt** | `unitig_background_frequency`, `variant_snp_check`, `validation_evidence`, `unitig_antibiotic_overlap`, `external_concordance` | 7 katmanlı kanıt + dış doğrulama + çapraz-antibiyotik |

**İlişki zinciri:**
`pipeline_runs` → `models` → (`unitig_model_scores`, `blast_annotations`, `unitig_background_frequency`, `variant_snp_check`) → hepsi `unitigs`'e bağlanır. `validation_evidence` her şeyi run_id ile köken'e bağlar. `external_concordance` modeli dış araçlarla kıyaslar.

---

## 3. Referans / Köken Tabloları

### 3.1 `pipeline_runs` — Köken çıpası (provenance anchor)
KB'deki **her satırın nereden geldiğini** garanti eden tablo. Bir antibiyotik pipeline'ının bir koşusu = bir satır. (21 satır)

| Kolon | Anlam |
|---|---|
| `run_id` | Birincil anahtar. Format: `{organizma}__{antibiyotik}__{UTC-zaman}__{git7}` (örn. `kpneumoniae__cefotaxime__20260707T1720__5b76f47`). Kayıtları bu koşuya bağlar. |
| `organism` | Organizma slug'ı (`ecoli`, `kpneumoniae`). `organisms` tablosuna anahtar. |
| `antibiotic` | Hedef antibiyotik (kanonik ad). |
| `git_commit` | O koşuyu üreten kodun 40-karakter git commit hash'i → **kod tekrarlanabilirliği**. |
| `git_dirty` | 0/1 — çalışma ağacı temiz miydi (0) yoksa yerel yamalı mı (1). |
| `card_version` | Kullanılan CARD veritabanı sürümü (**4.0.1**) → BLAST anotasyonlarının kaynağı. |
| `kmc_version`, `xgboost_version` | Araç sürümleri. |
| `random_seed` | Rastgelelik tohumu (42) → **istatistiksel tekrarlanabilirlik**. |
| `config_hash` | Veri + konfigürasyon parmak izi (sha256) → veri sürümü değişirse fark edilir. |
| `min_support` | O koşuda kullanılan **uyarlanabilir özellik filtresi** eşiği (`max(5, %1 × n_genom)`). |
| `n_genomes` | Modele giren genom sayısı. |
| `created_at` | ISO-8601 UTC zaman damgası. |

> **Neden önemli:** Reviewer "bu kaydı nasıl yeniden üretirim?" derse → git_commit + seed + config_hash + card_version dördü birlikte tam reçeteyi verir.

### 3.2 `organisms` — Organizma referansı *(şema 0.5.0'da eklendi)*
Cross-phylum (gram +/−) genelleme iddiasını yapılandırır. (4 satır)

| Kolon | Anlam |
|---|---|
| `organism` | Slug (birincil anahtar). |
| `display_name` | Tam ad (`Escherichia coli`). |
| `taxid` | NCBI taksonomi ID'si (562, 573, 1280). |
| `gram_stain` | `negative` / `positive` → gram-negatif vs gram-pozitif ayrımı. |
| `phylum` | Filum (`Pseudomonadota`, `Bacillota`). |

> **Anlam:** E. coli + K. pneumoniae = gram-negatif (Pseudomonadota); S. aureus (planlanan 3.) = gram-pozitif (Bacillota) → **"organizma-agnostik, cross-phylum"** genelleme hikâyesinin altyapısı.

### 3.3 `antibiotics` — Antibiyotik referansı
Her antibiyotiğin sınıfı ve klinik/mekanik meta bilgisi. (16 satır)

| Kolon | Anlam |
|---|---|
| `antibiotic` | Kanonik ad (birincil anahtar). |
| `drug_class` | İlaç sınıfı (`penicillins`, `cephalosporins`, `quinolones`, `aminoglycosides`, `tetracyclines`, `folate_pathway_inhibitors`, `beta_lactams_carbapenems_others`; ayrıca `macrolides`, `lincosamides` — registry'de 9 sınıf). |
| `mechanism_type` *(0.5.0)* | Baskın direnç mekanizması tipi: **`acquired`** (edinilmiş gen — TEM, CTX-M, AAC…) vs **`target_snp`** (hedef-gen nokta mutasyonu — gyrA/parC, quinolone'lar). Tez'in "edinilmiş-gen vs hedef-SNP showcase"ının anahtarı. |
| `who_aware` *(0.5.0)* | DSÖ (WHO) AWaRe kategorisi: **Access** (birinci basamak) / **Watch** (dikkat, direnç riski yüksek) / **Reserve** (son çare). Klinik önem bağlamı. |

### 3.4 `models` — Eğitilmiş modeller ve performansları
Her antibiyotik-organizma modeli = bir satır, tutulan-dışı (held-out) değerlendirme metrikleriyle. (21 satır)

| Kolon | Anlam |
|---|---|
| `model_id` | Birincil anahtar (otomatik artan). |
| `run_id` | Köken bağı (`pipeline_runs`). |
| `antibiotic` | Hedef. |
| `n_trees` | XGBoost ağaç sayısı. |
| `operating_threshold` | Karar eşiği (0/1 sınıflandırma noktası). |
| `roc_auc` | Tek-split ROC-AUC. |
| `roc_auc_ci_low`, `roc_auc_ci_high` | ROC-AUC %95 güven aralığı. |
| `pr_auc` | Precision-Recall AUC (dengesiz veri için bilgilendirici). |
| `mcc` | Matthews korelasyon katsayısı (dengeli tek-sayı metriği). |
| `balanced_accuracy` | Dengeli doğruluk (sınıf dengesizliğine dayanıklı). |
| `accuracy` | Ham doğruluk. |
| `auc_mean_seeds` | **Lineage-aware CV ROC-AUC ortalaması** (5-seed, PopPUNK kümeleriyle StratifiedGroupKFold) → **rapor edilen asıl performans** (soy/popülasyon-yapısı confounding'e karşı dürüst metrik). |
| `auc_std_seeds` | Aynısının standart sapması. |
| `n_features` *(0.5.0)* | Modelin matrisindeki toplam unitig (özellik) sayısı. |

> **Kritik ayrım:** `roc_auc` (tek split) ≠ `auc_mean_seeds` (lineage-CV). **Tezde savunulan metrik `auc_mean_seeds`'tir** çünkü aynı soydan genomların hem eğitim hem teste sızmasını (lineage leakage) önler. `01_performance` figürü bunu gösterir.

### 3.5 `kb_metadata` — FAIR meta veri (tek satır)
API `/metadata` ucunun döndürdüğü, KB'nin künyesi.

| Kolon | Anlam |
|---|---|
| `kb_schema_version` | Şema sürümü (**0.6.0**) — anlamsal sürümleme. |
| `card_version` | CARD sürümü (4.0.1). |
| `zenodo_doi` | Kalıcı DOI (yayın/Zenodo deposit sonrası dolar — FAIR "Findable"). |
| `license` | CC-BY-4.0. |
| `created_at`, `n_unitigs`, `n_models` | Üretim zamanı + toplam unitig (2363) + model (21). |

---

## 4. Çekirdek Biyobelirteç Tabloları

### 4.1 `unitigs` — Unitig sözlüğü (deduplike)
Tüm modeller arasında paylaşılan, tekilleştirilmiş unitig havuzu. Diğer her şey buraya bağlanır. (2363 satır)

| Kolon | Anlam |
|---|---|
| `unitig_id` | Birincil anahtar. |
| `sequence` | Unitig'in DNA dizisi (TEKİL). Değişken uzunlukta (gen-uzunluğu → BLAST'lanabilir). |
| `k` | Unitig'i kuran de Bruijn k değeri (**21**). |

### 4.2 `unitig_model_scores` — Model önemi + kararlılık skorları
Her (unitig, model) çifti için XGBoost önemi ve seçim kararlılığı. (2250 satır)

| Kolon | Anlam |
|---|---|
| `unitig_id`, `model_id` | Hangi unitig, hangi model. |
| `gain` | XGBoost **Gain** önem skoru (unitig'in modele katkısı). |
| `in_gain_topn` | 0/1 — tek modelin gain top-N listesinde mi (adım **07**). |
| `selection_frequency` | **CPSS seçim sıklığı** (adım **13**): B=100 alt-örnekte unitig kaç kez seçildi (0–1). |
| `stable` | 0/1 — `selection_frequency ≥ 0.6` mı (**kararlı biyobelirteç**). |
| `composite_score` | Bileşik skor: `stability × log10(1/E) × identity` (biyoloji + istatistik birleşimi). |
| `mean_abs_shap` | Ortalama |TreeSHAP| (CPSS satırları için; adım 13) → model-agnostik önem. |
| `selection_method` | `gain_seed` (07/07b, tek-model gain) veya `cpss` (adım 13, kararlılık seçimi). PK'nın parçası. |

> **CPSS (Complementary Pairs Stability Selection):** Meinshausen-Bühlmann / Shah-Samworth. Bir unitig birçok alt-örnekte tekrar seçiliyorsa (π≥0.6) "kararlı" sayılır → **rastgele/gürültü özelliklerini eler, PFER (yanlış-pozitif beklentisi) sınırlı**.

### 4.3 `blast_annotations` — Biyolojik kimlik (CARD + NCBI + ARO)
Her unitig'in gerçek dünyadaki karşılığı: hangi direnç geni? (2250 satır)

| Kolon | Anlam |
|---|---|
| `annotation_id` | Birincil anahtar. |
| `unitig_id`, `model_id` | Hangi unitig, hangi model bağlamında. |
| `source_db` | `card` (yerel CARD) veya `ncbi` (uzak NCBI nt). |
| `gene_symbol` | Gen sembolü (örn. `TEM-1`, `CTX-M-15`, `KPC-2`, `gyrA`). |
| `description` | Hit açıklaması. |
| `identity_pct` | Dizi kimliği yüzdesi. |
| `coverage` | Hizalama kapsamı (hizalama uzunluğu / unitig uzunluğu). |
| `evalue` | BLAST E-değeri (istatistiksel anlamlılık). |
| `tier` | **Güven seviyesi**: `confirmed` > `candidate` > `weak` > `none` (identity + coverage eşiklerine göre; adım 09). Reification güvencesi: "confirmed" ≠ nedensellik iddiası, yüksek-güven eşleşme. |
| `aro_accession` | CARD/ARO ontoloji erişim numarası. |
| `aro_gene_family` | ARO gen ailesi (örn. `CTX-M beta-lactamase`) → biyolojik olarak anlamlı seviye. |
| `aro_drug_class` | ARO ilaç sınıfı. |
| `aro_resistance_mechanism` | ARO direnç mekanizması (örn. antibiyotik inaktivasyonu, hedef değişimi). |

---

## 5. Doğrulama / Kanıt Tabloları

### 5.1 `unitig_background_frequency` — R vs S prevalans (adım 10)
Bir unitig gerçekten dirençli izolatlarda mı zenginleşmiş? (1203 satır)

| Kolon | Anlam |
|---|---|
| `unitig_id`, `model_id` | Hangi unitig/model. |
| `prevalence_resistant` | Dirençli (R) izolatlarda görülme oranı. |
| `prevalence_susceptible` | Duyarlı (S) izolatlarda görülme oranı. |
| `prevalence_overall` | Genel görülme oranı. |
| `delta_prevalence` | R − S prevalans farkı (ne kadar ayırt edici). |
| `odds_ratio` | Odds oranı (R'de olma bahsi). |
| `fisher_p` | Fisher exact test p-değeri. |
| `discriminative` | 0/1 — `|delta| ≥ eşik VE p < alfa` (istatistiksel olarak ayırt edici mi). |

### 5.2 `variant_snp_check` — Hedef-gen SNP kontrolü (adım 11)
Nokta mutasyonu mekanizmaları (quinolone gyrA/parC gibi) için CARD varyant-model allel kontrolü. (193 satır)

| Kolon | Anlam |
|---|---|
| `unitig_id`, `model_id` | Hangi unitig/model. |
| `card_model` | CARD varyant modeli. |
| `snp` | Nokta mutasyonu (örn. `S83L`, `S80I`). |
| `allele_class` | `resistant_allele` (direnç alleli) / `wildtype` / `other` / `ambiguous`. |

> **Anlam:** Edinilmiş gen (BLAST ile bulunur) ile **hedef-SNP mekanizmasını ayırır**. E. coli ciprofloxacin = gyrA S83L + parC S80I → BLAST homolog recovery ~0 (beklenen), ama SNP kontrolü `resistant_allele` yakalar. İki farklı biyoloji, iki farklı araç.

### 5.3 `validation_evidence` — Kanıt defteri (M11) — **7 KATMANIN KALBİ**
Her doğrulama sonucunun tek, birleşik defteri. Her satır = bir kanıt parçası. (~5850 satır)

| Kolon | Anlam |
|---|---|
| `evidence_id` | Birincil anahtar. |
| `unitig_id` | Hangi unitig. |
| `evidence_type` | **7 katmandan biri** (aşağıda). |
| `evidence_source` | Kaynak + sürüm (örn. `CARD 4.0.1`, `CPSS (B=100)`, `pyseer LMM lineage-corrected`). |
| `evidence_score` | İlgili skor (E-değeri / delta-AUC / Fisher p / seçim sıklığı / p-değeri…). |
| `pipeline_run_id` | Köken bağı. |

**7 doğrulama katmanı (`evidence_type`):**

| # | evidence_type | Ne test eder | Kaynak adım |
|---|---|---|---|
| 1 | `blast` | Biyolojik kimlik — bilinen direnç geni mi? | 08 (CARD+NCBI) |
| 2 | `background_frequency` | R'de S'ye göre zenginleşmiş mi? | 10 |
| 3 | `snp` | CARD varyant-modelinde direnç alleli mi? | 11 |
| 4 | `permutation_mda` | Özellik önemi rastgele-permütasyon null'una karşı anlamlı mı? | 12 |
| 5 | `label_permutation` | **Model** AUC'si etiket-karıştırma null'unu geçiyor mu? | 12b |
| 6 | `stability_selection` | CPSS'te kararlı mı (π≥0.6, PFER-sınırlı)? | 13 |
| 7 | `pyseer_lmm` | **Soy-düzeltmeli** GWAS anlamlılığı (popülasyon yapısı çıkarılınca hâlâ anlamlı mı)? | 14 |

> **Neden 7 katman?** Tek bir metrik yanıltıcı olabilir (ör. bir unitig sadece belirli bir soyla birlikte gidiyordur). Yedi ortogonal test aynı unitig'i doğruladığında, o biyobelirteç **gerçek direnç sinyali** olma olasılığı çok yüksektir — "istatistiksel sinyal ≠ biyolojik nedensellik" reviewer itirazına yapısal cevap.

### 5.4 `unitig_antibiotic_overlap` — Çapraz-antibiyotik örtüşme (adım 15, S1/H3) — *(0.6.0'da organizma-farkındalıklı yapıldı)*
Aynı unitig birden çok antibiyotikte kararlı mı? (organizma başına ayrı) (31 satır: ecoli 18 + kpneu 13)

| Kolon | Anlam |
|---|---|
| `unitig_id` | Paylaşılan unitig. |
| `organism` *(0.6.0)* | **Hangi organizma** — E. coli gentamicin ile K. pneu gentamicin'i **birbirine karıştırmamak** için (unified KB'de kritik). |
| `antibiotic_a`, `antibiotic_b` | Örtüşen antibiyotik çifti. |
| `same_class` | 0/1 — aynı registry ilaç sınıfı mı. |

> **H3 hipotezi:** "Sınıf-içi (β-laktam) örtüşme > sınıf-arası örtüşme?" Bulgu: ampicillin~cefotaxime (aynı sınıf) **ortak kararlı unitig yok** — çünkü ampicillin=TEM, cefotaxime=CTX-M/CMY (aynı sınıf, farklı enzimler). **Biyolojik olarak anlamlı negatif bulgu.**

### 5.5 `external_concordance` — Dış doğrulama (M13) — *(0.5.0/0.6.0)*
Bizim modelimiz, referans genotip araçlarına (AMRFinderPlus, ResFinder) karşı, **aynı held-out test genomlarında**, EUCAST/CLSI fenotipine göre. (48 satır)

| Kolon | Anlam |
|---|---|
| `model_id` | Hangi model. |
| `caller` | **`model`** (bizimki) / **`amrfinderplus`** / **`resfinder`**. |
| `reference` | Fenotip standardı (`EUCAST/CLSI (held-out test)`). |
| `n_test` | Test genomu sayısı (üç tahminci de aynı sette — leakage-free). |
| `sensitivity`, `specificity` | Duyarlılık / özgüllük. |
| `balanced_accuracy` | **Dengeli doğruluk** (asıl kıyas metriği). |
| `cohen_kappa` | Cohen κ (fenotiple uyum). |
| `major_error_rate` | **ME** (FDA): yanlış-dirençli oranı. |
| `very_major_error_rate` | **VME** (FDA): yanlış-duyarlı oranı (klinik olarak en tehlikeli hata). |

> **En güçlü bulgu:** K. pneumoniae **ciprofloxacin**'de modelimiz bACC **0.926**, AMRFinderPlus 0.538, ResFinder 0.540. Quinolone direnci nokta-mutasyon (gyrA/parC) kaynaklı; gen-tabanlı araçlar bunu kaçırıyor, **unitig modeli yakalıyor**. (`07_external_concordance` figürü)
> **Dürüst not:** cefotaxime (Kp) model bACC ≈ 0.495 — bu, o modelin küçük held-out test dilimindeki (n=200) head-to-head değeri; lineage-CV genelleme AUC'si (Figür 1) 0.77'dir.

---

## 6. KB Nasıl Sorgulanır?

- **Streamlit arayüzü:** `streamlit run scripts/kb_app.py` → 5 sekme: Biyobelirteçler / Kanıt zinciri / Model & Provenance / Çapraz-antibiyotik (H3) / Dış doğrulama (M13). Varsayılan DB: `results/kb/amrk.db`.
- **REST API:** `uvicorn scripts.kb_api:app` → `/api/v1/kmers`, `/kmers/{sequence}` (tam kanıt zinciri), `/overlap`, `/stats`, `/metadata` (FAIR).
- **Tidy tablolar (analiz/figür için):** `scripts/kb_tables.py` → `models_summary.csv`, `kb_overview.csv`, `biomarkers.csv`, `mechanisms.csv`.
- **Tez figürleri:** `scripts/kb_figures.py` → 8 figür (kapak, performans, PFER, mekanizma-ısıharitası, konkordans, significance, evidence-katmanları, external).

---

## 7. Danışman Sorularına Hazır Cevaplar

**7.1 "Neden ham k-mer değil de unitig?"**
Ham 21-mer'ler: ~50.8M özellik, aşırı redundant (tek SNP → onlarca örtüşen k-mer), ve **21 bp çok kısa → BLAST E-değeri anlamsız → geni tanıyamazsın**. Unitig'ler: de Bruijn grafiğinin dallanmayan yolları birleştirilir → ~0.8–3.5M özellik (≈25× az), **değişken uzunluk = gen-uzunluğu = CARD/NCBI'a eşlenebilir**, redundancy çözülür, DBGWAS/pyseer standardı. Yani unitig hem boyut hem **biyolojik yorumlanabilirlik** kazandırır.

**7.2 "Bir biyobelirtecin gerçek olduğunu nereden biliyorsun?"**
7 ortogonal katman (`validation_evidence`): BLAST kimliği + R/S prevalans + SNP allel + MDA permütasyon + etiket-permütasyon + CPSS kararlılık + soy-düzeltmeli pyseer LMM. Bir unitig yedisini de geçiyorsa, tek bir soya bağlı gürültü değil, gerçek direnç sinyalidir.

**7.3 "Popülasyon yapısı (lineage) sonuçları çarpıtmıyor mu?"**
İki katmanda korunuyoruz: (a) rapor edilen performans `auc_mean_seeds` = **lineage-aware CV** (PopPUNK kümeleriyle GroupKFold — aynı soy hem train hem teste sızmaz); (b) `pyseer_lmm` = **soy-düzeltmeli GWAS** (kinship matrisi ile popülasyon yapısı çıkarılır). Aminoglikozit ko-direnci pyseer'de ELENİR (mekanizmayı ayırır) ama TEM β-laktamaz GEÇER → yöntem çalışıyor.

**7.4 "İstatistiksel sinyali biyolojik nedensellikle karıştırmıyor musun?"**
Hayır — `tier` sistemi (confirmed/candidate/weak) ve "associational-not-causal" politikası bunu açıkça ayırır. Bir unitig "confirmed" olsa bile bu **yüksek-güven eşleşme** demektir, nedensellik iddiası değil. Kanıt-öncelikli tasarım (her iddia `validation_evidence`'ta izlenebilir).

**7.5 "Dış doğrulama var mı?"**
Evet — `external_concordance` (M13): modelimiz + AMRFinderPlus + ResFinder, **aynı held-out test genomlarında**, EUCAST/CLSI fenotipine göre, bACC/κ/FDA ME-VME ile. Model quinolone (gyrA SNP) gibi araçların kaçırdığı mekanizmalarda onları geçiyor.

**7.6 "Bu kayıt tekrarlanabilir mi?"**
Evet — her kayıt `pipeline_runs` üzerinden git_commit + random_seed + config_hash + card_version'a bağlı. Aynı reçete → aynı sonuç.

**7.7 "Cross-organism / organizma-agnostik iddian ne?"**
Aynı ilaç iki organizmada aynı gen ailesini kurtarıyor mu? Evet: gentamicin=AAC(3)/AAC(6′) (her ikisi), ciprofloxacin=gyrA/parC (her ikisi), trimethoprim/sulfa=sul (her ikisi), sefalosporinler=CTX-M/SHV (her ikisi). (`03_cross_organism` figürü). S. aureus (gram-pozitif) eklenince cross-phylum iddiası tamamlanacak.

---

*Bu belge `results/kb/amrk.db` (şema 0.6.0) içeriğine dayanır. Şema değişirse `scripts/lib/kb_schema.py`'deki `KB_SCHEMA_VERSION` ile birlikte güncellenmelidir.*
