# AMR K-mer Knowledge Base — Tez ve Proje Yol Haritası

**Hazırlayan:** Sentetik analiz (tüm materyaller entegre)  
**Tarih:** Haziran 2026  
**Kapsam:** ML_AMR_Prediction_v2 GitHub Repo · AMR_Technical_Report.html · AMR_KmerDB_Roadmap.docx · AMR_Project_Analysis.md · AMR_Literature_Review.md  
**Araştırma Sorusu:** *"ML modellerinin öğrendiği yüksek önem skoruna sahip k-mer'lerden biyolojik olarak anlamlı, sorgulanabilir ve yeniden kullanılabilir bir AMR Knowledge Base oluşturulabilir mi?"*

---

## ⚠️ Önceki Analizlere Katkı, Düzeltme ve Eleştiri

Bu raporu okumadan önce, mevcut analizlerin güçlü ve zayıf yönlerini açıkça kayıt altına almak gerekiyor.

### AMR_Project_Analysis.md — En Güvenilir Kaynak, Ama İki Eksik Var

Bu rapor doğru yapılmış ve en güvenilir belgedir. Tutarsızlık tablosu (10 satır) ve kritik problem listesi (P-01 ile P-15) bilimsel analizin sağlam bir çerçevesidir. Ancak iki önemli eksik var:

**Eksik 1:** "Aşama A / Aşama B" ayrımı doğru ama hangi Aşama A maddelerinin tamamlanmadan tez savunmasına *girilmemesi* gerektiği net belirtilmemiş. P-01 (data leakage) çözülmeden tez teslim edilemez; bu sadece "yüksek öncelikli" değil, submission blocker'dır.

**Eksik 2:** "5-fold CV out-of-core mimaride zor ama zorunlu" deniyor, somut çözüm önerilmiyor. Bu raporda ilerleyen bölümlerde pratik çözüm sunulacak.

### AMR_KmerDB_Roadmap.docx — Vizyon Güçlü, Ama Üç Kritik Hata Var

**Hata 1 — "Doktora Yeterlilik Seviyesi" iddiası:** HTML raporu "Yüksek Lisans Tezi" diyor. DOCX "Doktora Yeterlilik" diyor. Mevcut durumda KB sıfır implementasyon, P-01'den P-07'ye kritik bug'lar açık. Bu durum yüksek lisans tezine dahi ancak kritik düzeltmelerle yeterli olabilir; doktora yeterlilik seviyesinden çok uzak.

**Hata 2 — Yanlış BLAST bulgusu:** DOCX Bölüm 1.2: "Ciprofloxacin Rank 3'te gyrA QRDR bölgesi (E-value 1.5)." Repodaki gerçek BLAST çıktısı: gyrA, **gentamicin Rank 3'te** (E-value 1.5). Ciprofloxacin Rank 3'te CARD'da blaS1 ve AAC6'-Ib (E-val 2.8) var. Bu hata tez metodoloji bölümüne aynen taşınırsa reviewer hemen yakalar.

**Hata 3 — NAR Database Issue hedefi 6 ay için gerçekçi değil:** NAR DB Issue, erişilebilir public web UI, production hosting ve external validation zorunlu kılıyor. Bunlar mevcut durumla 6 ayda yapılamaz. Birincil hedef **Database (Oxford)** veya **Briefings in Bioinformatics** olmalı.

### AMR_Literature_Review.md — Kaliteli, Bir Doğrulama Gerekiyor

Laing et al. 2025 (Çalışma 11) için verilen DOI `10.1016/j.anucene.2025` Nükleer Enerji Dergisi'ne ait bir DOI formatı görünüyor; AMR biyoinformatik çalışması için şüpheli. Bu referansı tez yazımı öncesinde PubMed veya Google Scholar üzerinden doğrulayıp correct DOI ile güncellemek gerekiyor.

Güncellenmiş inceleme Çalışma 14–17 (Meinshausen stability selection teorik temeli, AntiMicrobial-KG, BRIDGE çerçevesi, DeepMDC) eklenerek 17 çalışmaya yükseldi. Diğer 16 referans tutarlı ve iyi seçilmiş.

Yeni incelemenin üç katkısı KB tasarımını doğrudan etkiliyor ve ilgili bölümlere entegre edildi: (1) FAIR prensipleri ve provenance/evidence tracking'in AMR k-mer veritabanlarında henüz hayata geçirilmediğine dair açık literatür boşluğu — Bölüm 1.3'e `kb_schema_version`, `zenodo_doi` alanları ve FAIR metadata gereksinimleri olarak yansıtıldı; (2) Reviewer perspektifinden öngörülen altı somut itiraz — "Feature reproducibility nasıl güvence altına alınıyor?" ve "FAIR uyumluluğu nasıl belgeleniyor?" soruları must-have maddelerine (M10, M11) ve ayrı bir reviewer riskleri alt bölümüne taşındı; (3) Knowledge Base yayınlanabilirlik kriterleri — validation evidence yapısı ve reification fallacy'ye karşı yapısal güvencenin metodoloji notuna dönüştürülmesi gereksinimi Bölüm 2 Should Have'e (S9, S10) eklendi.

### BLAST Sonuçlarına Dair Yeni Bir Gözlem

Mevcut analizlerde yeterince vurgulanmamış: Gentamicin Rank 4 ve 5'teki k-mer'ler ve pek çok "Intergenic Region (Esc11, Esc16...)" NCBI hit'i —  aynı izolat adları tekrar tekrar çıkıyor. Bu, o k-mer'lerin gerçek direnç sinyali taşımayabileceğini; belirli klonal izolatlara özgü lineage sinyali olabileceğini kuvvetle düşündürüyor. Phylogenetic bias kontrolü teorik bir gereklilik değil, mevcut veri için zaten görünen bir problem.

---

## 0. METODOLOJİK REVİZYON — Literatür Taraması Kararları (Haziran 2026)

> İki turluk sistematik literatür taraması (Bölüm A–F + uygulama "nasıl" turu) sonucu kesinleşen, pipeline'ı ve tez çerçevesini bağlayan kararlar. **Bu bölüm aşağıdaki §1–§6'nın önüne geçer; çelişki olursa burası geçerlidir.** Kaynaklar: Yu/Barquist 2024 (PLOS Biol, popülasyon-yapısı confounding), Jaillard 2018 (DBGWAS/unitig), Lees 2018/2019 (pyseer, PopPUNK), Meinshausen-Bühlmann 2010 + Shah-Samworth 2013 (stability selection/CPSS), Lundberg 2020 (TreeSHAP), Chklovski 2023 (CheckM2), Alcock 2023 (CARD/ARO), Bortolaia 2020 (ResFinder), Feldgarden 2021 (AMRFinderPlus), Nguyen 2018/2019 (XGBoost AMR).

### 0.1 MUST-HAVE mimari değişiklikler (yayın öncesi zorunlu)
1. **Soya-duyarlı çapraz doğrulama — rastgele/chunk CV TAMAMEN TERK.** **PopPUNK** küme etiketleri (5470 genom ~1-2h) → **`GroupKFold`** (5–10 fold, stratified-group; <10 izolatlı ST'ler elenir/"diğer"e; ST131 gibi baskın kladlar downsample veya PopPUNK alt-kümelemesi). Doğruladığı: **raporlanan genelleme-AUC + feature stability** (final model yine TÜM veride eğitilir). Neden: random CV, AUC'yi **%20–30 şişirir**; model gerçek mekanizma yerine soyu ezberler (Yu 2024). **En büyük reviewer-blocker.**
2. **⭐ Unitig temsiline geçiş — ham k-mer TERK.** `bcalm2` + `unitig-caller` → ikili unitig varlık/yokluk matrisi; **downstream XGBoost aynı kalır** (satır=izolat, sütun=unitig 0/1). ~10M k-mer → **~730k unitig**, ~212GB→**~18GB** RAM, ~7h→**~50dk**. Faydalar: boyut indirgeme + yorumlanabilirlik (uzun, BLAST-eşlenebilir) + GWAS-standart + min_support baskısını azaltır. KMC yalnız QC/spektrum (02/02b) için kalır; matris girdisi unitig olur. **Bu, önceki min_support/hız/GPU sancılarını kökten çözer.**
3. **Stabilite: 5-seed → CPSS (Complementary Pairs Stability Selection, B=100, %50 alt-örnek, yerine koymadan).** Aşamalı (staged) seçim: (a) sıfır-varyans + frekans filtresi → (b) univariate **Chi²/Mutual Information** ile aday sete in (≈ binlerce) → (c) adaylarda CPSS B=100 → (d) seçim frekansı eşiği **π≥0.6** (M&B/Shah-Samworth PFER sınırı; 0.6 teorik olarak doğrulandı). 5-seed istatistiksel garanti vermez.
4. **Feature importance: Gain → SHAP (TreeSHAP), yalnız son kısıtlı (aday) model üzerinde.** Tam matriste asla; saniyeler sürer.
5. **Dış doğrulama (MUST).** BV-BRC'den **zamansal (son ~2 yıl) veya coğrafi hold-out** (ikisi de hakemlerce kabul). **AMRFinderPlus** (tercih — SNP+edinilmiş) + ResFinder ile head-to-head. Metrikler: **balanced accuracy, sensitivity, specificity, Cohen's Kappa, McNemar**; FDA bandı ME≤%3, VME≤%1.5–7.5.
6. **pyseer LMM + Bonferroni** — bir k-mer/unitig'in mekanizmayla ilişkili olduğu iddiası için popülasyon-yapısı düzeltmeli p-değeri (SHAP üstü k-mer'leri çapraz-kontrol).

### 0.2 Eklenecek doğrulamalar / bileşenler
- **Permütasyon/anlamlılık: MDA (Mean Decrease in Accuracy)** — modeli bir kez eğit, test feature'ını karıştır, doğruluk düşüşünü ölç (yüzlerce yeniden-eğitim YOK). Alternatif: donmuş-HP + early-stop ile etiket permütasyonu. Eşik p<0.05 / **BH-FDR Q<0.05**.
- **Discriminativeness'a BH-FDR:** |Δprev|≥0.10 + Fisher exact **+ Benjamini-Hochberg Q<0.05** (step 10).
- **Genom QC (hibrit):** BV-BRC metadata ön-filtre (EvalG/EvalCon — savunulabilir) **+ yerel CheckM2** (completeness ≥95–99%, contamination ≤5%, sıkı ≤2–3%) **+ QUAST** (N50 ≥50kb, contig üst sınırı). CheckM1 değil **CheckM2** (~saatler / 5470 genom).
- **ARO/CARD ontoloji eşlemesi (KB):** unitig→CARD BLAST → ARO accession → `aro_index.tsv` + `card.json` parse → 5'li şema (**ARO ID, gen adı, gen ailesi, direnç mekanizması, ilaç sınıfı**). `kb_schema`'ya alanlar.
- **Reifikasyon dili (S10):** "neden olur / belirler / tetikler" YOK → "ilişkilidir / yüksek tahmin gücü / istatistiksel sinyal". SHAP ≠ nedensellik (Takefuji 2025).

### 0.3 Doğrulananlar (değişiklik YOK — literatür bizi destekledi)
- **k=21** (hesaplama yönetilebilir + BLAST-eşlenebilir); **ikili presence/absence + max_bin=2** (topluluk standardı, count'tan üstünlük yok); **sınıf-ağırlık (neg/pos), SMOTE YOK**, 0.5 threshold raporla + PR-eğrisi tartış; **BV-BRC + EUCAST/CLSI + Intermediate çıkar** (standart); **4373 genom yeterli** (modern bant 1000–5000+).
- **AUC ~0.93 literatürle UYUMLU** (per-antibiyotik raporla): ciprofloxacin 0.83–0.99, gentamisin 0.89–0.99, ampisilin volatil 0.51–0.98, TMP-SMX 0.90–0.96, tetrasiklin 0.88–0.98, sefotaksim 0.72–0.81.

### 0.4 Özgünlük çerçevesi (B1 — KRİTİK)
**"İlk ML-AMR veritabanı" İDDİA ETME** (BV-BRC zaten XGBoost ile "AMR Regions" sunuyor). Yeni çerçeve: *"istatistiksel hata-sınırlı (PFER/M&B-CPSS), soy-farkında doğrulanmış, k-mer/unitig çözünürlüklü, şeffaf + FAIR ilk AÇIK AMR biyobelirteç bilgi tabanı."* Farklılaştırıcılar: **CPSS stability selection + lineage-aware validation + ARO ontolojisi + açık (Docker/DOI/API) altyapı.** Bu çerçevede özgünlük güveni **Yüksek**; "ilk DB" çerçevesinde **Orta-Düşük**.

### 0.5 Showcase antibiyotik + hedef dergi
- **Showcase:** ciprofloxacin (gyrA/parC **SNP** — step 11 pozitif kontrol) + β-laktam (**edinilmiş gen** pozitif kontrol). Merkeze "**test edilebilir yeni biyolojik keşif**" koy (Database Oxford şartı).
  - ✅ **GERÇEKLEŞTİ (2026-06-29):** **ampicillin** = edinilmiş gen (TEM-256/257/258; CARD recovery %47/H2 TRUE; SNP 0). **ciprofloxacin** = hedef-gen SNP (**gyrA S83L + parC S80I = resistant_allele**, step 11; CARD homolog recovery ~0 = beklenen). İki mekanizma, her biri doğru araçla doğrulanmış — showcase çifti hazır. KB çok-antibiyotikli (schema 0.4.0).
- **Hedef dergi: Database (Oxford)** veya **Briefings in Bioinformatics** (NAR Database Issue 6 ayda gerçekçi değil). FAIR minimum: Zenodo DOI, REST API, Docker/Conda <30dk kurulum, 2-yıl erişim garantisi, ARO eşlemesi.

### 0.6 Pipeline'a etki — değişecek/eklenecek scriptler
- **00a/00:** BV-BRC QC metadata + yıl/coğrafya alanlarını çek (external hold-out + QC için).
- **YENİ unitig adımı (02c/03 yerine):** `bcalm2` + `unitig-caller` → unitig matrisi. KMC 02/02b QC'de kalır.
- **YENİ lineage adımı:** PopPUNK küme etiketleri → GroupKFold; 04/05/07b'nin chunk-level split'i bununla değişir.
- **07b:** CPSS (B=100) + Chi² aday filtresi + SHAP (5-seed yerine).
- **07:** SHAP importance (Gain yerine).
- **10:** BH-FDR ekle.
- **YENİ:** pyseer LMM; MDA permütasyon; external validation + AMRFinderPlus/ResFinder concordance; CheckM2/QUAST QC.
- **09 / KB:** ARO eşleme alanları; `populate_database.py` ARO + provenance + lineage + CPSS skorları.

### 0.7 Major reviewer saldırı noktaları (savunma hazır olmalı)
1. *Soy sızıntısı:* "model ST131 k-mer'lerini ezberledi" → lineage-aware GroupKFold + pyseer LMM ile çözülür.
2. *Mükerrerlik:* "BV-BRC zaten yapıyor" → PFER-sınırlı CPSS + şeffaflık + ARO ile farklılaş.
3. *Feature kararsızlığı:* "5 seed gürültü" → CPSS B=100 + SHAP.
4. *Agresif filtreleme:* "%1 nadir plazmidi siler" → mutlak ≥10 + unitig collapse (oransal eşik değil).

### 0.8 GERÇEKLEŞTİ — Çok-organizmalı KB tamamlandı (2026-07-07)

Showcase, tek-organizma iki-antibiyotikten **iki-organizma 17-modele** ölçeklendi. Tek unified KB (`results/kb/amrk.db`, `models.organism` ayırır): **E. coli 7 + Klebsiella pneumoniae 10 = 17 model, 7 ilaç sınıfı** (penisilin, sefalosporin, karbapenem, kinolon, aminoglikozid, folat, tetrasiklin). Her model tam boru hattından geçti (soy-CV + CPSS/PFER + MDA + etiket-permütasyon + pyseer LMM + CARD/NCBI ARO).

**Özgünlük çerçevesi güncellendi (§0.4'ü genişletir):** artık *"**çok-organizmalı, çok-sınıflı**, PFER-sınırlı (CPSS), soy-farkında doğrulanmış, unitig-çözünürlüklü, şeffaf + FAIR ilk açık AMR biyobelirteç bilgi tabanı."* Cross-organism konkordans yeni bir farklılaştırıcı.

**Cross-organism mekanizma konkordansı (aynı ilaç, iki tür → aynı biyolojik mekanizma — hesaplı-tekrarlanan pozitif kontrol):**
- gentamicin → **AAC(3)-II** (her ikisi); ciprofloxacin → **gyrA/parC** (her ikisi); trimethoprim/sulfa → **dfr** (E. coli sul2+dfrA15, K. pneu dfrA14).
- K. pneumoniae amiral: **KPC karbapenemaz** — hem meropenem hem imipenem için (%100 id, E≈1e-23…1e-80); E. coli 3. kuşak sefalosporin: **CTX-M/CMY** (cefotaxime, ceftazidime).

**PFER biyolojiyi yansıtıyor (yeni tez-gözlemi):** konsantre mekanizmalar düşük PFER (cip 0.10, meropenem 2.96, amox-clav 1.03 = OXA-1); dağınık/ko-taşınan aminoglikozid direnci yüksek PFER (gentamicin 50.6/35.7). Bu, "PFER-sınırlı KB" iddiasını her antibiyotik için sayısal ve biyolojik olarak destekler.

**Altyapı (hepsi `main`'de):** (1) env-parallelism refactor — `get_target()` ile CLI>env>config, config-mutex kalktı, çok-antibiyotik/organizma paralel; (2) unified çok-organizma KB; (3) slash-güvenli kombo antibiyotikler (underscore canonical); (4) genom yedeği (Drive). Kalan: METHODOLOGY'ye panel + Zenodo (M10) deposit; yeni antibiyotikler için M13 concordance (`amrfinder_keywords` eklendi, AFP çıktısına karşı doğrulanmalı).

---

## 1. Knowledge Base Mimarisi

DOCX'teki veri modeli temelden sağlam; 8 tablo yapısı kabul edilebilir. Aşağıda bileşen bazında kısa değerlendirme ve önemli eklemeler sunuluyor.

### 1.1 Veri Modeli

DOCX'te tanımlanan `antibiotics → models → kmer_model_scores → kmers → blast_annotations → kmer_antibiotic_overlap → validation_evidence` zinciri doğru. **İki ekleme gerekiyor:**

**`pipeline_runs` tablosu:** Her scriptın çalışma metadatası (git commit hash, Python versiyonu, çalışma tarihi, kullanılan seed, KMC versiyonu) ayrı bir tabloda tutulmalı. `models.pipeline_version` tek bir string yetersiz; commit hash'e bağlanmadan reproducibility claim yapılamaz.

**`kmer_background_frequency` tablosu:** Her k-mer'in tüm genomlarda (dirençli + duyarlı) görülme sıklığı. Bir k-mer stabil ve yüksek Gain'li olsa bile tüm genomlarda bulunuyorsa discriminative değil. Bu arka plan sıklığı yoksa "novel marker" ile "ubiquitous conserved sequence" ayrımı yapılamaz. Bu tablo, hem H4 hipotezi hem de Lineage-specific k-mer tanımı için kritik.

### 1.2 Database Yapısı

PostgreSQL tercih doğru ve gerekçesi sağlam. Önerilen ek: **İlk 3 ay için SQLite ile başla, Month 4'te PostgreSQL'e migrate et.** Sebebi: PostgreSQL kurulum + server yönetimi tez sürecinin bant genişliğini tüketmemeli; aynı SQLAlchemy ORM kodu her iki veritabanıyla çalışır. Sadece production aşamasında (publication için public erişim gerektiğinde) PostgreSQL zorunlu hale gelir.

**İndeks stratejisi:** `kmers.sequence` üzerinde B-tree (exact match), `blast_annotations.gene_symbol` üzerinde B-tree, `kmer_model_scores.stability_score` üzerinde GiST. pg_trgm extension 21-mer substring araması için faydalı ama zorunlu değil.

### 1.3 Provenance Sistemi

DOCX'teki `models` tablosu iyi tasarlanmış. Kritik eklemeler:
- `git_commit_hash` (char[40]): Pipeline versiyonunu tam olarak sabitleyen tek güvenilir referans
- `card_version` (varchar): Şu an hiçbir yerde kayıtlı değil — tüm BLAST sonuçları referanssız
- `bvbrc_download_date` + `bvbrc_filter_criteria` (text): Hangi genome kalite filtreleri uygulandı?
- `n50_threshold` + `max_contig_count`: Assembly kalite eşikleri kayıt altında olmalı
- `kb_schema_version` (varchar): KB şema versiyonu — semantic versioning (v0.1.0, v0.2.0); her şema değişikliğinde artırılır, API ve veri yayımında bu sürüme atıfta bulunulur
- `zenodo_doi` (varchar): Her versiyonlanmış release için Zenodo kalıcı DOI — FAIR Findable gereksinimi; Methods bölümünde veri erişim referansı olarak kullanılır

**FAIR uyumluluğu notu:** AMR_Literature_Review.md boşluk analizi (madde 6–7), FAIR standartlarında ML-derived AMR k-mer KB'sinin henüz hayata geçirilmediğini doğruluyor. Minimum FAIR karşılama: (F) Zenodo DOI, (A) public REST API, (I) JSON-LD metadata veya Dublin Core header, (R) CC-BY 4.0 lisans + `kb_schema_version`. Makine-okunabilir metadata şeması (`/api/v1/metadata` endpoint üzerinden) yayın öncesi hedef; tez için Zenodo DOI + `/api/v1/metadata` JSON çıktısı yeterli.

### 1.4 Confidence Scoring Sistemi

DOCX'teki üç tier (confirmed/candidate/weak) bilimsel olarak doğru. **Bir kompozit skor eklenebilir:**

`composite_score = stability_score × log10(1/E-value) × (identity_pct/100)`

Bu formül; stabil, yüksek homoloji ve düşük E-value'ya sahip k-mer'leri doğal olarak üste taşır. Tek sayı ile KB sıralaması yapmayı kolaylaştırır. Bu, tez Tables bölümüne koyulacak "Tablo 2: Top-20 KB Entries by Composite Score" için kullanışlı.

### 1.5 Feature Stability Sistemi

> **⚠ REVİZE (bkz. §0.1):** 5-seed repeated holdout istatistiksel olarak yetersiz bulundu. Yerine **CPSS (Complementary Pairs Stability Selection, B=100, %50 alt-örnek, π≥0.6)** + aşamalı Chi²/MI aday filtresi + SHAP kullanılacak (Meinshausen-Bühlmann 2010 / Shah-Samworth 2013 PFER sınırları). Aşağıdaki "güçlendirme" notu tarihsel bağlam için bırakıldı.

5-seed selection frequency yaklaşımı doğru metodoloji (Mahé & Tournoud 2018 referanslı). **Bir güçlendirme:** Bağımsız seed'lerin yanı sıra **5-fold CV ile fold bazlı stability** hesaplamak metodolojik açıdan daha güçlü bir kombinasyon:

- `seed_stability`: 5 farklı random seed ile eğitim → kaç seed'de top-50'ye girdi?
- `cv_stability`: 5-fold CV'de kaç fold'da top-50'ye girdi?
- `combined_stability`: (seed_stability + cv_stability) / 2

Bu çift-boyutlu stabilite, tek bir metrik olarak kullanılabilir ve Mahé (2018) üzerine net bir metodolojik katkı sağlar.

### 1.6 Cross-Antibiotic Analiz

DOCX'te kmer_antibiotic_overlap tablosu var. **Eksik: istatistiksel anlamlılık.** İki antibiyotik arasındaki k-mer overlap'inin şans eseri mi yoksa biyolojik olarak anlamlı mı olduğunu test etmek için **hypergeometric test** kullanılabilir:

- Toplam k-mer havuzu N (tüm stable k-mer'ler)
- Antibiyotik A'nın stable seti: K
- Antibiyotik B'nin stable setinde A'yı da geçen k-mer sayısı: k
- Beklenen overlap (null hipotez): K × |set_B| / N

Fisher's exact test veya hypergeometric p-değeri, H3 hipotezini (beta-laktam overlap > cross-class overlap) istatistiksel olarak destekler. Bu olmadan "H3 desteklendi" iddiası sadece betimsel kalır.

### 1.7 Validation Katmanı

DOCX Section 8 detaylı. Önceliklendirme açısından önerilen sıra:

1. **Known mechanism recovery rate** (1 tablo, hesaplanabilir): Confirmed tier k-mer'lerin kaçı bilinen ARG? Bu proof-of-concept için minimum gereklilik.
2. **Permutation test** (1 script, ~2 gün): Null dağılım → gerçek model Gain vs. shuffled label Gain
3. **Temporal split validation** (1 hafta): 2019 öncesi train, 2020+ test
4. **ResFinder cross-validation** (1 hafta): Confirmed k-mer'ler içeren genomları ResFinder'a yükle, concordance raporla

**Validation evidence tablosu gereksinimi:** Her validation sonucunun `validation_evidence` tablosuna yazılması zorunlu. Minimum alan yapısı:
- `evidence_type`: ENUM('blast', 'resfinder', 'temporal_split', 'permutation') — hangi validation türü?
- `evidence_source`: CARD v3.X.X / ResFinder 4.0 / shuffled_labels — tam versiyon referansı
- `evidence_score`: float — BLAST için E-value, temporal/permutation için AUC delta
- `pipeline_run_id`: `pipeline_runs` tablosuna FK — hangi model çalışması üretti?

Bu tablo olmadan KB'deki her kaydın nereden geldiği izlenemez ve "Validation evidence nasıl belgeleniyor?" reviewer sorusuna yanıt verilemez. M11 olarak must-have listesine taşındı.

### 1.8 API

FastAPI + SQLAlchemy tercih doğru. **Minimum viable API (yayın için yeterli):**

```
GET /api/v1/kmers?antibiotic=gentamicin&min_stability=0.6&tier=confirmed
GET /api/v1/kmers/{sequence}
GET /api/v1/overlap?ab1=ampicillin&ab2=cefotaxime
GET /api/v1/stats
GET /api/v1/metadata
```

Rate limiting (100 req/min per IP), CORS headers, OpenAPI docs otomatik — bunlar FastAPI ile 1-2 günlük iş. Yayın için "API erişilebilir" olması yeterli; production SLA gerekmez.

**FAIR veri erişim stratejisi:** `/api/v1/metadata` endpoint KB şema versiyonu (`kb_schema_version`), Zenodo DOI, toplam kayıt sayısı ve lisans bilgisini (CC-BY 4.0) makine-okunabilir JSON olarak döndürür — FAIR Accessible ve Reusable gereksinimlerini karşılar. Versiyonlanmış endpoint yapısı (`/api/v1/kmers?version=0.1`) KB'nin farklı sürümlerini aynı anda erişilebilir kılar; bu Database (Oxford) submission için gereken "bilimsel veritabanı yönetim" kriteri olarak Methods bölümüne girilmeli. Tez için bu endpoint'in çalışır olması yeterli; external hosting yayın aşamasına kadar ertelenebilir.

### 1.9 Web Arayüzü

**Tez için:** Streamlit ile basit bir arama formu yeterli. 3-4 günlük iş. `pip install streamlit` ile mevcut ortamda çalışır.

**Publication için:** Streamlit yeterli değil; GitHub Pages üzerinde statik HTML + JavaScript (Fetch API) daha güvenilir ve ücretsiz hosting ile çalışır. CARD, ResFinder web arayüzleri ile benzer görünüm akademik standart sağlar.

---

## 2. Yayınlanabilir Hale Gelmek İçin Gerekenler

### Must Have (Bunlar Olmadan Paper Submission Yapılmaz)

| # | Gereklilik | Mevcut Durum | Tahmini Süre |
|---|---|---|---|
| M1 | **P-01 fix: Data leakage** — Youden's J test setten kaldır, sadece train/val üzerinden hesapla | Aktif bug | 2 saat |
| M2 | **Soya-duyarlı çapraz doğrulama (§0.1)** — PopPUNK küme etiketleri + GroupKFold. *(REVİZE: rastgele/5-seed CV artık YETERSİZ; lineage-aware ZORUNLU.)* | ✅ DONE (2026-06-23) — 07b StratifiedGroupKFold, 324 PopPUNK kümesi, AUC 0.9505±0.01 | 1-2 hafta |
| M3 | **E-value confidence tier sistemi** — Mevcut BLAST sonuçlarını confirmed/candidate/weak olarak yeniden sınıflandır | Yok | 3 gün |
| M4 | **Feature stability (07b) — CPSS, B=100, %50 alt-örnek, π≥0.6 + SHAP (§0.1)** *(REVİZE: 5-seed yerine Meinshausen-Bühlmann/Shah-Samworth)* | ✅ DONE (2026-06-25) — `13` CPSS (Chi² ön-filtre→200 fit→π≥0.6) + yerleşik TreeSHAP; **39 kararlı, PFER≤5.4**; `13b` CARD-tier+ARO: 5 confirmed (TEM-256/257/258, APH(6)-Id, AAC(6')-Ib7, tam-boy); KB'de (schema 0.3.0, method='cpss') | 1 hafta |
| M12 | **Unitig temsiline geçiş (§0.1)** — bcalm2 + unitig-caller (ham k-mer yerine) | ✅ DONE (2026-06-23) — 4.94M unitig matrisi, `03u_unitig_matrix.py` | 1 hafta |
| M13 | **Dış doğrulama (§0.1)** — zamansal/coğrafi hold-out + AMRFinderPlus/ResFinder concordance (Kappa, McNemar, bACC) | ✅ **CONCORDANCE DONE (2026-07-02)** — `16_external_concordance.py` + `lib/concordance`: AMRFinderPlus 2026-05-15.1 + ResFinder 4.5.0 on **5468 genomes** vs EUCAST/CLSI (bACC/κ/**FDA ME-VME**) + **leakage-free model-vs-tool head-to-head** on held-out test genomes. Model bACC **amp 0.873 / cef 0.925 / cip 0.928** — matches ResFinder (cef,cip) & beats AMRFinderPlus (cip,amp). Temporal/geo hold-out (needs year/geo metadata) still open. | 1 hafta |
| M14 | **pyseer LMM + Bonferroni (§0.1)** — popülasyon-yapısı düzeltmeli k-mer/unitig anlamlılığı | ✅ DONE (2026-06-25) — `14` pyseer LMM (kinship alt-örnekten, LMM adaylarda); eşik 1.09e-5; **26/39 CPSS-kararlı soy-düzeltmeli anlamlı; 3/3 TEM β-laktamaz geçti, aminoglikozid ko-direnç GEÇMEDI** (LMM mekanizmayı ayırdı); KB'de (pyseer_lmm evidence, TEM p≈1e-109) | 3 gün |
| M15 | **Genom QC (§0.2)** — BV-BRC metadata ön-filtre + yerel CheckM2 + QUAST | ✅ **DONE (2026-07-02)** — `02d_genome_qc.py` (CheckM2 amr-checkm2.sif + QUAST amr-tools.sif) run on all 5470 assemblies: **pass 5312/5470 = 97.1%** (fail 158: N50 131, completeness 29, contigs 26, contamination 13) at completeness≥95 / contamination≤5 / N50≥50kb / contigs≤500. Outputs `02d_genome_qc_ecoli.csv` + summary JSON + advisory exclusion list (158 IDs). Advisory: fails <3% → not retrained; a data-quality statement for Methods. (02b IQR advisory also in place.) | 3 gün |
| M16 | **ARO/CARD ontoloji eşlemesi (§0.2)** — KB'de ARO ID + gen ailesi + mekanizma + ilaç sınıfı | ✅ DONE (2026-06-24) — 09 ARO mapping + `blast_annotations` ARO kolonları (kb_schema 0.2.0), 13/60 eşlendi | 3 gün |
| M5 | **Reproducibility fix** — Ampicillin ve ciprofloxacin matrislerini Zenodo'ya yükle veya yeniden üretim pipeline'ını belgele | ✅ Büyük kısmı DONE (2026-06-29) — kanonik 04→populate yeniden üretim koşusu; her `pipeline_runs` git_commit+seed+config_hash(sha256)+CARD 4.0.1 damgalıyor; Drive yedeği (`TRUBA_25626`). Zenodo DOI kaldı (M10) | 3 gün |
| M6 | **CARD version kaydı** — CARD hangi versiyonu? config.yaml ve Methods bölümüne | ✅ DONE (2026-06-24) — `AMR_CARD_VERSION` env → `kb_metadata.card_version=4.0.1` | 2 saat |
| M7 | **Known mechanism recovery rate tablosu** — Confirmed tier k-mer'lerin kaçı bilinen ARG? | ✅ DONE (2026-06-24) — 09 `08_validation_metrics`: recovery %32 (9/28 stable), **H2 FALSE** (<%40; §0.4'e göre sorun değil) | 2 gün |
| M8 | **PostgreSQL KB + populate script** — En az gentamicin ve cefotaxime için çalışan KB | ✅ SQLite DONE (2026-06-24) — ampicillin `amrk.db` schema 0.2.0 dolu (65 kmer); Postgres migration + 2. antibiyotik kaldı | 2 hafta |
| M9 | **Permutation test** — Null dağılım vs. gerçek model | ✅ DONE & saved (2026-06-24) — `12` MDA (0/51 FDR-sig: unitig redundancy) + `12b` label-perm null (REAL 0.9534 ≫ null mean 0.4994/max 0.5521, **p=0.0099 significant**) | 3 gün |
| M10 | **KB versiyonlama + Zenodo DOI** — `kb_schema_version` alanı DB'ye eklenmeli; git semver tag (v0.1.0) + Zenodo deposit; kalıcı DOI Methods bölümüne ve README'ye yazılmalı | Yok | 2 gün |
| M11 | **Validation evidence tablosu** — `validation_evidence` şeması + pipeline entegrasyonu; her BLAST, ResFinder ve permutation test sonucu `evidence_type`, `evidence_source`, `evidence_score`, `pipeline_run_id` alanlarıyla kayıt altına alınmalı | Sıfır | 3 gün |

### Should Have (Bunlar Paper'ı Güçlendirir, Zorunlu Değil)

| # | Gereklilik | Tahmini Süre |
|---|---|---|
| S1 | Cross-antibiotic overlap analizi (10_cross_antibiotic.py) + hypergeometric test | 1 hafta |
| S2 | Temporal external validation (2019 öncesi vs 2020+) | 1 hafta |
| S3 | Phylogenetic bias kontrolü — MLST tabanlı klonal kompleks analizi | 2 hafta |
| S4 | P-02 fix: shell injection güvenliği | 2 saat |
| S5 | P-03 fix: eval_metric config'den oku | 2 saat |
| S6 | P-04 fix: colsample METHODOLOGY çelişkisi — kodu ya da dokümanı düzelt | 1 saat |
| S7 | ResFinder cross-database validation | 1 hafta |
| S8 | FastAPI minimal REST endpoint ✅ **DONE (2026-07-02)** — `scripts/kb_api.py` (FastAPI + CORS + OpenAPI /docs) over `scripts/lib/kb_queries.py` (pure sqlite3, unit-tested): `/api/v1/kmers` (filter antibiotic/tier/min_stability/stable_only), `/kmers/{sequence}` (full evidence chain), `/overlap`, `/stats`, `/metadata`. `uvicorn scripts.kb_api:app` | 3 gün |
| S9 | FAIR metadata şeması — `/api/v1/metadata` endpoint üzerinden makine-okunabilir KB metadata (schema.org veya Dublin Core JSON-LD); CC-BY 4.0 lisans beyanı; Database (Oxford) submission için "veri erişim stratejisi" kriteri ✅ **DONE (2026-07-02)** — `/api/v1/metadata` returns `kb_schema_version` + `zenodo_doi` + `license` (CC-BY-4.0) + `card_version` + n_unitigs/n_models + antibiotics as JSON | 2 gün |
| S10 | Reification fallacy güvencesi notu — confidence tier sisteminin "istatistiksel sinyal ≠ biyolojik nedensellik" ayrımını açıkça belgeleyen metodoloji paragrafı; Methods ve Discussion bölümlerine eklenmeli; Takefuji 2025 eleştirisine yanıt ✅ **DONE (2026-07-02)** — `METHODOLOGY.md §4.4` (associational-not-causal wording policy + 3 yapısal güvence: layered orthogonal evidence, measured confounding/H3 negative finding, provenance-over-assertion) | 1 gün |

### Nice To Have (Bunlar Gelecek Çalışmaya Kalabilir)

| # | Gereklilik | Not |
|---|---|---|
| N1 | Streamlit veya statik web arayüzü | NAR için zorunlu, Database Oxford için "strongly recommended" |
| N2 | Multi-organism genişleme (K. pneumoniae) | V2 hedefi |
| N3 | SHAP entegrasyonu (Gain'e ek) | Takefuji 2025 eleştirisi bu ile yanıtlanır |
| N4 | Neo4j graph layer | k-mer→gen→yol analizi için uzun vade |
| N5 | Docker container | Reproducibility için ideal ama tez için zorunlu değil |
| N6 | NCBI local nt veritabanı (>200GB) | Remote BLAST yerine — HPC gerektiriyor |

---

## 3. Tez Kapsamı

### 3.1 Tez Başlığı Önerileri

**Önerilen Birincil Başlık (Briefings in Bioinformatics / Database Oxford uyumlu):**  
*"AMRK-DB: A Stability-Filtered, Cross-Antibiotic Knowledge Base of Machine Learning-Derived Genomic Signatures for Antimicrobial Resistance in Escherichia coli"*

**Alternatif — Metodoloji odaklı:**  
*"From Feature Importance to Knowledge Base: Constructing a Queryable Repository of ML-Derived AMR k-mer Signatures with Confidence Tiering and Cross-Antibiotic Analysis"*

**Türkçe tez başlığı (YÖK formatı):**  
*"Makine Öğrenmesi Tabanlı Özellik Önemi Çıktılarından Kararlılık Filtreli ve Antibiyotiklerarası Analize Sahip Sorgulanabilir Bir Antimikrobiyal Direnç K-mer Bilgi Tabanının Oluşturulması: Escherichia coli Örneği"*

### 3.2 Temel Araştırma Soruları

**Ana soru (ana tez sorusuyla uyumlu):**  
Alignment-free k-mer tabanlı ML modellerinden türetilen feature importance çıktıları, E. coli'de AMR mekanizmalarına dair biyolojik olarak anlamlı, yeniden üretilebilir ve sorgulanabilir bir bilgi deposu oluşturmak için yeterli bilgi içeriyor mu?

**Alt sorular:**

1. XGBoost Gain skoru yüksek k-mer'ler, farklı random başlangıç koşullarında istatistiksel olarak tutarlı mı? *(H1 — stability)*
2. Stability filtresi geçen k-mer'ler, bilinen AMR mekanizmaları ile anlamlı biyolojik örtüşme gösteriyor mu ve kataloglenmamış direnç bölgelerini işaret ediyor mu? *(H2 — biological validity)*
3. Aynı antibiyotik sınıfı için eğitilmiş modeller, farklı sınıflar arasındaki modellerden istatistiksel olarak anlamlı daha yüksek k-mer örtüşmesi gösteriyor mu? *(H3 — mechanism specificity)*
4. Bu bilgileri yapılandıran bir knowledge base, araştırmacılar için yeniden kullanılabilir değer üretiyor mu? *(H4 — utility)*

### 3.3 Test Edilebilir Hipotezler

| Hipotez | Test Yöntemi | Kabul Kriteri | Reddedilirse |
|---|---|---|---|
| **H1 (Stability):** Top-50 k-mer setinin ≥%60'ı 5 seed'de tutarlı | Seed-based selection frequency | Median stability ≥ 0.6 | KB geçerliliği sorgulanır; metodoloji zayıflar |
| **H2 (Biological validity):** Stability ≥ 0.6 olan k-mer'lerin ≥%40'ı CARD'da bilinen ARG'ye E≤1e-3, id≥95% ile haritaslıyor | BLAST + confidence tier | ≥40% confirmed rate | Proof-of-concept zayıflar; ama H1 geçerliyse metodoloji katkısı devam eder |
| **H3 (Mechanism specificity):** Beta-laktam içi overlap (AMP–CEF) > cross-class overlap (AMP–GEN) | Hypergeometric test | p < 0.05 | Cross-antibiotic analiz anlamsızlaşır; ama model kalitesine dair önemli negative finding |
| **H4 (Novel candidates):** Stability ≥ 0.6 k-mer'lerin >0 fraksiyonu CARD'da eşleşme bulamıyor | BLAST miss rate | >0% novel candidate | Bu hiç gerçekleşmez — toplam "confirmed" k-mer sayısı az bile olsa novel adaylar çıkacak |

**Methodolojik not H4 için:** H4 "ret edilemez" bir hipotez olarak formüle edilmemelidir. Daha savunulabilir hali: "novel candidate fraksiyonunun istatistiksel dağılımı ve genomik bağlamı karakterize edilebilir mi?" — bu bir keşif sorusu, bir hipotez testi değil.

### 3.4 Beklenen Bilimsel Katkılar

**Katkı 1 — Metodolojik:** Feature importance'ı knowledge base'e dönüştüren sistematik bir pipeline tanımlanması: Gain → seed stability → CV stability → BLAST confidence tier → cross-antibiotic overlap. Bu işlem zinciri literatürde tanımlanmamış.

**Katkı 2 — Analitik:** E. coli'de dört antibiyotik sınıfı için ML-derived stable k-mer'lerin sistematik karşılaştırması ve cross-class MDR sinyal adaylarının tanımlanması. Gentamicin modelindeki gyrA sinyali bu katkının somut örneği.

**Katkı 3 — Kaynak:** AMRK-DB: API ile sorgulanabilir, versiyonlanmış, confidence-tiered bir AMR k-mer kaynağı. Diğer araştırmacılar kendi ML modellerinin çıktılarıyla karşılaştırabileceği bir referans noktası.

**Katkı 4 — Sınırlılık şeffaflığı:** Confidence tier sistemi ve lineage-specific flagging ile "ML'nin bulduğu = gerçek" yanılgısına karşı metodolojik dürüstlüğün somutlaştırılması.

---

## 4. 6 Aylık Yol Haritası

| Ay | Odak | Yapılacaklar | Beklenen Çıktı |
|---|---|---|---|
| **Ay 1** | **Pipeline Sağlamlaştırma** | • P-01 fix (data leakage, 2 saat) • P-02 fix (shell injection) • P-03/P04 fix (eval_metric, colsample) • BioPython requirements.txt ekle, SHAP kaldır • CARD version kaydet • Entrez email config'e taşı | Submission-blocker buglar kapalı; temiz, denetlenebilir pipeline |
| **Ay 1-2** | **Validation Altyapısı** | • 5-seed repeated holdout implementasyonu (out-of-core uyumlu alternatif: 5 farklı seed ile 80/20 split × 5, CV değil; ancak her split stratified) • `07b_feature_stability.py` — selection frequency hesaplama • Permutation test scripti • Ampicillin ve ciprofloxacin matrislerini Zenodo'ya yükle | ROC-AUC güvenilirlik kanıtı; feature stability skorları hazır |
| **Ay 2** | **BLAST Tier Sistemi** | • Mevcut 4 antibiyotik BLAST sonuçlarını confidence tier'larına yeniden sınıflandır • confirmed / candidate / weak dağılımını raporla • Known mechanism recovery rate hesapla (Tablo 2 adayı) • Hypergeometric test ile H1 ve H2 ön analiz | Savunulabilir BLAST annotation; ilk quantitative results |
| **Ay 3** | **Knowledge Base Build** | • `11_populate_database.py` — SQLite ile başla • En az 2 antibiyotik için KB'i doldur • `10_cross_antibiotic_analysis.py` + H3 test • `validation_evidence` tablosu şeması + BLAST sonuçlarını tabloya yaz • `kb_schema_version` alanını DB'ye ekle + git tag v0.1.0 + Zenodo ön kaydı (DOI rezerve et) | Çalışan SQLite KB; cross-antibiotic overlap sonuçları; validation evidence altyapısı hazır; v0.1.0 Zenodo DOI |
| **Ay 4** | **Validation & Tez Yazımı** | • Temporal external validation (2019–2020 split) → sonuçları `validation_evidence` tablosuna yaz • ResFinder cross-database validation (en az 2 antibiyotik) → concordance oranı raporla • PostgreSQL migration (production deployment için) • FastAPI minimal endpoint (`/kmers`, `/overlap`, `/metadata`) • FAIR metadata şeması: `/api/v1/metadata` JSON çıktısı + CC-BY 4.0 lisans beyanı • S10: Reification fallacy metodoloji notu — Methods taslağına ekle | Public-accessible KB; validation evidence tablosu dolu; FAIR minimum karşılama tamamlanmış; tez Methods ve Results taslakları |
| **Ay 5** | **Tez Yazımı + Phylogenetic Bias** | • MLST analizi ile lineage-specific k-mer flagleme (E. coli ST131 kontrolü öncelikli) • Tez Introduction + Discussion bölümleri | Bias-corrected KB; tez taslak tamamlanmış |
| **Ay 6** | **Finalizasyon** | • Danışman revizyonu • Web arayüzü (Streamlit veya statik HTML, en az 1 haftalık iş) • GitHub README + Zenodo DOI (v0.1.0 veya v1.0.0 release) • KB yayınlanabilirlik kriteri kontrolü: public API erişim, versiyonlanmış schema, validation evidence tablosu dolu, FAIR metadata endpoint çalışır, reification fallacy metodoloji notu tez metninde • Tez teslim hazırlığı + paper preprint (bioRxiv) | Tez teslim hazır; KB yayınlanabilirlik kriterleri karşılanmış; paper submission |

### Out-of-Core Mimari ile CV Çözümü

Önceki analizlerde "CV zorlu" denilip geçildi. Pratik çözüm:

```python
# 5-seed repeated holdout (gerçek CV yerine — out-of-core uyumlu)
SEEDS = [42, 123, 777, 1024, 2025]
for seed in SEEDS:
    X_train_idx, X_test_idx = stratified_split(genome_list, seed=seed)
    # Mevcut pipeline zaten bu yapıyı destekliyor
    # Her seed için ayrı model eğit, feature importance çıkar
    # stability_score = sum(selected_in_seed_i) / len(SEEDS)
```

Bu yaklaşım gerçek k-fold CV değil ama Mahé & Tournoud (2018)'in resampling yaklaşımıyla metodolojik tutarlı ve mevcut out-of-core mimariye entegrasyonu 1-2 günlük iş.

---

## 5. Nihai Tavsiye

**Seçim: B — AMR K-mer Knowledge Base yönüne evrilmek.**

Ama bu karar, sık yapılan iki yanlış anlama ile yüklü olmamalı:

---

### Neden A yeterli değil?

"Sadece mevcut AMR prediction pipeline'ını geliştirmek" yolunda gittiğinizi varsayalım. E-value tier düzeltiyorsunuz, CV ekliyorsunuz, tez yazıyorsunuz. Ne elde edilir?

Bilimsel sonuç: "E. coli'de 21-mer + XGBoost ile AMR tahmin edilebilir (ROC-AUC ~0.93–0.99)."

Bu sonuç, Laing et al. 2025 (34 antibiyotik, 4.300 izolat), Jain et al. 2025 (P. aeruginosa, SHAP), Kover (12 organizma, 56 antibiyotik) gibi çalışmalarla zaten kapsanmış bir alana yeni bir tekrar katkısı olur. Yayın yapılabilir ama özgünlük iddiası zayıf. Tek farklılaştırıcı nokta out-of-core mimarinin ölçeklenebilirlik yönü olur — bu ise ana yenilik iddiası için çok teknik ve dar kalır.

---

### Neden B mantıklı?

KB yönüne evrilmenin üç bilimsel gerekçesi var:

**Gerekçe 1 — Gerçek bir literatür boşluğu dolduruluyor.** "ML feature importance → stability filter → cross-antibiotic overlap → confidence-tiered queryable knowledge base" işlem zinciri literatürde sistematik olarak tanımlanmamış. Laing 2025 novel marker keşfediyor ama KB'ye dönüştürmüyor. Kover model görselleştiriyor ama sorgulayabilir veritabanı üretmiyor. ValizadehAslani 2020 biological insight çıkarıyor ama stability analizi ve cross-antibiotic analiz yok. Bu gap özgündür ve gerçektir.

**Gerekçe 2 — Mevcut biyolojik sinyal bu iddiayı destekleyecek kadar güçlü.** AAC(3)-IId (gentamicin, E-val 2.8e-05, %100 identity) ve TEM-1 (ampicillin, E-val 2.8e-05, %100 identity) confidence tier sistemi uygulandıktan sonra "confirmed" kalıyor. Bunlar referanssız keşfedilmiş, literatürle uyumlu ve mekanistik olarak açıklanabilir sinyaller. Proof-of-concept için yeterli temel var.

**Gerekçe 3 — Gentamicin modeli beklenmedik bir gyrA sinyali taşıyor (Rank 3, E-val 1.5, %100 identity).** Bu, aminoglikozid modelinde bir fluorokinolon hedef geni —  bilinen ko-direnç ilişkisi ama ML modeli bunu referanssız öğrenmiş. Bu tam olarak "ML'nin ne öğrendiğini sorgulanabilir yapıya dönüştürmenin" değerini gösteren somut bir örnek. Bu sinyal A yolunda footnote olarak kalır; B yolunda merkezi bir bulgu haline gelir.

---

### B yolunu seçerken kaçınılması gereken tuzaklar

**Tuzak 1:** "KB planlandı = KB var" varsayımı. DOCX'teki 10 bölümlük mimari planı sıfır implementasyona eşit. Ay 1-2'de pipeline kritik buglarını kapatmadan KB'ye geçilmesi durumunda KB'ye koyulan verinin kalitesi güvensiz.

**Tuzak 2:** Confirmed tier'daki k-mer sayısı beklenenden az çıkabilir. Confidence tier sistemi uygulandıktan sonra "confirmed" kategori TEM-1 ve AAC(3)-IId ile sınırlı kalabilir. Bu tezi geçersiz kılmaz — az sayıda ama güvenilir confirmed kayıt + çok sayıda candidate/weak kayıt, metodolojinin değerini gösterir. Önemli olan tablo sayısı değil, şeffaflık.

**Tuzak 3:** KB'yi 6 ay içinde NAR Database Issue'ya submission yapmak için production kalitesinde bir servis olarak tasarlamaya çalışmak. Bu tez kapsamını aşar. Tez için: çalışan API + Zenodo hosted data yeterli. Production servis için zaman dilimi +6 ay daha.

---

### Öngörülen Reviewer Riskleri ve Yanıt Stratejisi

Literatür incelemesi ve proje analizi birleştirildiğinde altı somut reviewer itirazı öne çıkıyor. Bunlar hem tez savunması hem yayın sürecinde hazırlıklı olunması gereken noktalardır.

**Risk 1 — "Threshold test seti üzerinden optimize edilmiş (Youden's J). Bu data leakage'dır."**  
Yanıt stratejisi: P-01 fix (M1) Ay 1'de tamamlanmadan tez savunmasına girilmez. Bu submission blocker'dır. Çözülünce yanıt kendi içinde oluşur.

**Risk 2 — "Cross-validation yok. Tek split ile ROC-AUC=0.99 iddiası kabul edilemez."**  
Yanıt stratejisi: 5-seed repeated holdout (M2) sonuçlarını Methods'a ekle; her seed için AUC, stability skoru ve seçilen k-mer setinin Jaccard benzerliğini raporla. "Gerçek k-fold değil ama Mahé 2018 resampling yaklaşımıyla metodolojik tutarlı" argümanını Methods'ta açıkça yaz.

**Risk 3 — "Feature reproducibility nasıl güvence altına alınıyor?"**  
Yanıt stratejisi: 5-seed selection frequency (M4) + validation evidence tablosu (M11) + git commit hash + Zenodo DOI (M10) kombinasyonu bu soruyu kapatır. Her k-mer kaydına `seed_list`, `model_version`, `pipeline_run_id` eklendiğinde yanıt doğrudan KB'den gösterilebilir.

**Risk 4 — "BLAST annotation yüzeysel — E-value eşiği savunulamaz."**  
Yanıt stratejisi: Confidence tier sistemi (M3) ile E-value ≤ 1e-3 / ≤ 1.0 / ≤ 50 ayrımını uygula. Methods'ta 21-mer için E-value yorumunun sınırlarını açıkla. gyrA (E-val 1.5) "candidate" olarak etiketle, "confirmed" olarak sunma.

**Risk 5 — "FAIR uyumluluğu nasıl belgeleniyor?"**  
Yanıt stratejisi: `/api/v1/metadata` endpoint (S9) + Zenodo DOI (M10) + CC-BY 4.0 lisans + `kb_schema_version` alanı birlikte FAIR minimum karşılama olarak Methods bölümüne girer. "FAIR-F: Zenodo DOI, FAIR-A: public REST API, FAIR-I: JSON-LD metadata, FAIR-R: versiyonlanmış schema + lisans" yapısında sunulabilir.

**Risk 6 — "Yüksek Gain skoru nedensellik kanıtlamaz."**  
Yanıt stratejisi: Bu itiraz doğru ve beklenen. Reification fallacy güvencesi (S10) metodoloji notunda açıkça belgelenmeli: çalışma nedensellik iddia etmiyor, stability filtresiyle desteklenmiş istatistiksel sinyal kataloğu sunuyor. Her kaydın confidence tier sistemi "istatistiksel korelasyon" ile "biyolojik doğrulama" katmanlarını yapısal olarak ayırıyor — bu, AMR k-mer KB'lerinde tasarım ilkesi olarak ilk kez hayata geçirilen bir güvence mekanizması.

---

### Önerilen çalışma çerçevesi özeti

```
[Ay 1-2] Critical bug fix → Validation foundation
    → P-01, stability analizi, E-value tier
    → Bunlar bitmeden KB'ye tek satır kod yazılmaz

[Ay 3-4] KB Build → Cross-antibiotic analysis → FAIR altyapı
    → SQLite → PostgreSQL migration
    → validation_evidence tablosu + kb_schema_version + Zenodo DOI (M10, M11)
    → FastAPI minimal endpoint (/kmers, /overlap, /metadata)

[Ay 5-6] Validation → Tez + Preprint
    → Temporal validation, phylogenetic bias
    → KB yayınlanabilirlik kriteri kontrolü (FAIR minimum, reification fallacy notu)
    → Briefings in Bioinformatics veya Database (Oxford) submission
```

**Tek cümle:** Mevcut pipeline üretim kalitesine ve validation depth'e ulaşmadan yayın yoktur; ama validation tamamlandıktan sonra KB yönüne gitmemek bu projenin özgün katkısını ortadan kaldırır.

---

## Ek: Dergi Hedefi Değerlendirmesi

| Dergi | Impact | KB Gereksinim | 6 Ayda Mümkün? | Öneri |
|---|---|---|---|---|
| **Database (Oxford)** | IF ~4.7 | Çalışan DB + minimal web UI + description paper + FAIR-uyumlu versiyonlanmış veri yayımı (Zenodo DOI, M10) + validation evidence tablosu (M11) + provenance kaydı (M6 + `kb_schema_version`) | **Evet, gerçekçi** | ✅ Birincil hedef |
| **Briefings in Bioinformatics** | IF ~9.5 | Full methodology + validation + novel biological finding + provenance kaydı + confidence tier metodoloji notu (S10: reification fallacy güvencesi) | **Sınırda, zorlu** | ⚠️ Eğer phylogenetic bias + temporal val tamamlanırsa |
| **NAR Database Issue** | IF ~14.9 | Production hosting + public web UI + comprehensive external val + FAIR production infrastructure + external validation dataset | **Hayır, 6 ayda değil** | ❌ 12+ ay sonrası |
| **Bioinformatics (Oxford)** | IF ~5.8 | Strong methodology + adequate validation + reproducibility (M5 + M10) | **Evet, Aşama A+B ile** | ✅ Alternatif |
| **PLoS Computational Biology** | IF ~4.1 | Solid biology + reproducibility + FAIR minimum (Zenodo DOI + public API) | **Evet** | ✅ Güvenli alternatif |

---

*Tüm materyaller entegre analiz: ML_AMR_Prediction_v2 repo (263 dosya) + 5 proje belgesi + 17 literatür çalışması. Tez yazımı sürecinde bu belgenin 3. ve 4. bölümleri Methods, Introduction ve Discussion taslakları için doğrudan girdi olarak kullanılabilir.*

