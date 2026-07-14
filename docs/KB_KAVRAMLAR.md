# AMRK-DB — Kavramlar ve Örneklerle Açıklama

> Bu belge `KB_ACIKLAMA.md`'yi tamamlar: oradaki teknik terimleri **sade dille + KB'den gerçek örneklerle** açıklar. "Bu ne, ne işe yarıyor, ne anlama geliyor, nerede kullanılıyor?" sorularının cevabı. Danışman toplantısı için hazırlık.

İçindekiler:
1. [k-mer vs unitig — unitig varken k-mer neden var?](#1)
2. [GWAS nedir? (pyseer LMM)](#2)
3. [Provenance (köken) nedir?](#3)
4. [FAIR nedir?](#4)
5. [CPSS ve "stable" nasıl belirleniyor?](#5)
6. [SHAP ve Gain nedir? (özellik önemi)](#6)
7. [Bir sürü ROC var — hangisi önemli?](#7)
8. [CARD / NCBI BLAST hit'leri, E-değeri, identity, coverage](#8)
9. [Güven seviyeleri: confirmed / candidate / weak / none](#9)
10. [unitig_background_frequency nedir?](#10)
11. [7 doğrulama katmanı — bir unitig örneği üzerinden](#11)
12. [External concordance — model vs AMRFinderPlus vs ResFinder (ME/VME)](#12)
13. ["0 / None / tek değer" neden görünüyor?](#13)

---

<a name="1"></a>
## 1. k-mer vs unitig — unitig varken k-mer neden oluşturuluyor?

**k-mer:** Bir DNA dizisini sabit uzunlukta (bizde **21 bp**) parçalara bölmek. "ATTGC...CG" dizisini 1'er kaydırarak çıkan tüm 21 harfli pencereler. Bir bakteri genomu milyonlarca k-mer üretir.

**unitig:** k-mer'leri bir **de Bruijn grafiğinde** birleştirip, **dallanmayan yolları tek bir uzun parçaya** çökertmek. Yani art arda gelen ve hep birlikte giden k-mer'ler tek bir "unitig" olur (değişken uzunlukta, çoğu gen-uzunluğunda).

**Neden ikisi de var?** İkisi farklı amaç için:
- **k-mer** aşaması (02/02b, KMC) → **kalite kontrol ve spektrum analizi** için. Genomların k-mer profillerini karşılaştırıp aykırı/bozuk assembly'leri yakalamak. Modelin *özelliği* değil, bir **QC aracı**.
- **unitig** aşaması (03u) → modelin **asıl özellikleri**. KB'de sakladığımız, skorladığımız, BLAST'ladığımız şey unitig.

**Neden özellik olarak k-mer değil de unitig?**
- **Sayı:** ham k-mer ~**50.8 milyon** özellik; unitig ~**0.8–3.5 milyon** (≈25× az). *(KB örneği: ampicillin/E. coli modeli `n_features = 4,938,938` unitig)*
- **Yorumlanabilirlik:** 21 bp çok kısadır → BLAST E-değeri anlamsız çıkar → "bu k-mer hangi gen?" diye soramazsın. Unitig gen-uzunluğunda → **CARD/NCBI'a eşlenir, gerçek gen adı çıkar** (TEM, CTX-M…).
- **Redundancy:** Tek bir mutasyon onlarca örtüşen k-mer üretir (hepsi aynı şeyi söyler). Unitig bunları tek parçada birleştirir.

> **Tek cümle:** k-mer = kalite kontrol için ham malzeme; unitig = modele giren, biyolojik olarak okunabilir özellik.

---

<a name="2"></a>
## 2. GWAS nedir? (ve pyseer LMM)

**GWAS = Genome-Wide Association Study** (Genom Boyu İlişkilendirme Çalışması). "Genomdaki hangi konum/varyant, bir özellikle (burada: antibiyotik direnci) istatistiksel olarak ilişkili?" sorusunu tarar.

**Sorun:** Bakteriler klonaldır — aynı **soydan (lineage)** gelenler binlerce özelliği paylaşır. Bir unitig aslında dirençle değil, sadece "dirençli bir soyla birlikte gitmekle" ilişkili olabilir (yalancı ilişki / popülasyon-yapısı confounding).

**Çözüm — pyseer LMM (Linear Mixed Model):** Genomlar arası akrabalık (kinship) matrisini modele katıp **soy etkisini çıkarır**. Bir unitig soy etkisi çıkarıldıktan *sonra* hâlâ anlamlıysa → gerçek direnç sinyali.

- KB'de: `validation_evidence.evidence_type = 'pyseer_lmm'`, kaynak *"pyseer LMM lineage-corrected (Bonferroni 1.09e-05)"*.
- **Doğrulama:** Aminoglikozit ko-direnci pyseer'de ELENİR (soyla giden, mekanizma değil), ama TEM β-laktamaz GEÇER → yöntem gerçek mekanizmayı ayırıyor.
- *Örnek:* ceftazidime/E. coli modelinde **117** unitig pyseer-anlamlı çıktı.

---

<a name="3"></a>
## 3. Provenance (köken) nedir?

**Provenance = bir verinin nereden, nasıl, hangi koşullarla üretildiğinin kaydı.** "Bu sonucu birebir yeniden üretebilir miyim?" sorusunun cevabı.

KB'de `pipeline_runs` tablosu bunu tutar. Her model için:
- `git_commit` → sonucu üreten **kodun tam sürümü**
- `random_seed` (42) → **rastgelelik sabitlenmiş**
- `card_version` (4.0.1) → hangi referans veritabanı
- `config_hash` → veri + ayar parmak izi
- `created_at` → ne zaman

> **Neden önemli:** Bilimsel tekrarlanabilirlik. Reviewer "bu k-mer'i nasıl buldun, tekrar üretilir mi?" derse → bu dört alan tam reçeteyi verir. FAIR'in "Reusable" (yeniden kullanılabilir) ayağı.

---

<a name="4"></a>
## 4. FAIR nedir?

Bilimsel veriyi yayınlanabilir/kullanılabilir yapan 4 ilke:
- **F — Findable (Bulunabilir):** Kalıcı kimlik (Zenodo DOI — `kb_metadata.zenodo_doi`).
- **A — Accessible (Erişilebilir):** Açık API (`scripts/kb_api.py`, `/api/v1/...`).
- **I — Interoperable (Birlikte çalışabilir):** Standart formatlar + ARO/CARD ontolojisi (`blast_annotations.aro_*`).
- **R — Reusable (Yeniden kullanılabilir):** Açık lisans (CC-BY-4.0) + sürüm (`kb_schema_version`) + provenance.

> KB `kb_metadata` tablosu ve `/api/v1/metadata` ucu tam bu bilgiyi makine-okunur döndürür. Literatürde ML-türevi AMR k-mer veritabanlarında FAIR henüz nadir → **tezin özgün katkılarından biri**.

---

<a name="5"></a>
## 5. CPSS ve "stable" nasıl belirleniyor?

**CPSS = Complementary Pairs Stability Selection** (Tamamlayıcı Çiftler Kararlılık Seçimi; Meinshausen-Bühlmann / Shah-Samworth).

**Fikir:** Bir özellik gerçekten önemliyse, veriyi biraz değiştirsen bile (rastgele alt-örnekler alsan) **tekrar tekrar seçilmeli.** Sadece bir şansa seçilenler gürültüdür.

**Nasıl:** Veriden B=100 kez rastgele yarı-örnek alınır, her seferinde model önemli özellikleri seçer. Bir unitig'in **seçim sıklığı** (`selection_frequency`) = 100 denemenin kaçında seçildiği (0–1 arası).

**"stable" kriteri:** `selection_frequency ≥ 0.6` → `stable = 1`. Yani en az %60 tutarlılıkla seçilenler kararlı biyobelirteç.
- KB'de: `unitig_model_scores.selection_frequency`, `stable`, `selection_method='cpss'`.
- *Örnek:* `selection_frequency = 1.0, stable = 1` → bir unitig 100 denemenin 100'ünde de seçilmiş (çok güçlü).

**Bonus — PFER:** CPSS "beklenen yanlış-pozitif sayısını (PFER)" da sınırlar. Düşük PFER = temiz/konsantre sinyal (ör. carbapenem, quinolone); yüksek PFER = dağınık/ko-taşınan sinyal (ör. aminoglikozit). `02_cpss_pfer` figürü bunu gösterir.

---

<a name="6"></a>
## 6. SHAP ve Gain nedir? (özellik önemi)

İkisi de "bir unitig modele ne kadar katkı yapıyor?" sorusunu ölçer ama farklı yöntemle:

- **Gain** (`unitig_model_scores.gain`): XGBoost'un kendi iç ölçüsü — bu unitig ağaç bölünmelerinde ne kadar "kazanç" (hata azaltma) sağladı. *Örnek: gain = 54.9.*
- **SHAP** (`mean_abs_shap`): **SHapley Additive exPlanations** — oyun teorisinden gelen, model-agnostik bir yöntem. Her unitig'in her tahmine katkısını adil biçimde dağıtır; ortalama mutlak değeri o unitig'in genel önemidir. *Örnek: mean_abs_shap = 0.029.*

**Neden ikisi?** Gain hızlı ama modele özgü; SHAP daha ilkeli/karşılaştırılabilir. İkisinin birbirini desteklemesi güveni artırır.

> **Not:** Gain, `gain_seed` yöntemiyle seçilen unitig'lerde dolu; SHAP, `cpss` yöntemiyle işlenen (adım 13) unitig'lerde dolu. Bu yüzden bir satırda biri, öbür satırda diğeri boş görünebilir (→ §13).

---

<a name="7"></a>
## 7. Bir sürü ROC var — hangisi önemli?

**ROC-AUC** = modelin R/S ayırma gücü (0.5 = rastgele, 1.0 = mükemmel). KB'de birkaç ROC var; **hangisine bakılacağı kritik:**

| Kolon | Ne | Güvenilirlik |
|---|---|---|
| `roc_auc` | **Tek** rastgele train/test split'inde AUC | İyimser olabilir (soy sızıntısı riski) |
| `auc_mean_seeds` | **Lineage-aware CV** (PopPUNK kümeleriyle GroupKFold, 5-seed ortalama) | ✅ **RAPOR EDİLEN, GÜVENİLİR METRİK** |
| `auc_std_seeds` | Yukarıdakinin standart sapması | Kararlılık göstergesi |

**Neden `auc_mean_seeds` önemli?** Bakteriler klonal → aynı soydan genom hem eğitim hem teste düşerse model "ezberler", AUC yapay yüksek çıkar. Lineage-aware CV bunu engeller (aynı soy tek tarafta kalır) → **dürüst genelleme performansı.**

**KB örneği (fark çarpıcı):**
- ampicillin/E. coli: `roc_auc = 0.924` ama `auc_mean_seeds = 0.951` (bu modelde yakın)
- **ciprofloxacin/K. pneu: `roc_auc = 0.967` ama `auc_mean_seeds = 0.868`** → tek-split iyimser; lineage-CV daha düşük ama **dürüst**. Sunumda **0.868'i** savun.

> `01_performance` figürü hep `auc_mean_seeds ± auc_std_seeds` gösterir.

---

<a name="8"></a>
## 8. CARD / NCBI BLAST hit'leri — E-değeri, identity, coverage

Bir unitig'in **hangi gerçek gen olduğunu** bulmak için iki veritabanına BLAST'lanır:
- **CARD** (yerel): antibiyotik direnci genlerinin küratörlü veritabanı (ARO ontolojisiyle).
- **NCBI nt** (uzak): tüm bilinen diziler.

`blast_annotations` tablosundaki her satır bir hit. Kolonlar ne anlama gelir:

| Kolon | Anlam | Örnek |
|---|---|---|
| `gene_symbol` | Eşleşen gen | `AAC(3)-IId`, `CTX-M-276`, `TEM-256` |
| `identity_pct` | Dizi kimliği % (ne kadar birebir) | `100.0` = tam eşleşme |
| `coverage` | Unitig'in ne kadarı hizalandı | `1.0` = tamamı |
| `evalue` | **E-değeri**: bu eşleşmenin **şans eseri** olma beklentisi. Ne kadar **küçük** o kadar anlamlı. | `1.35e-56` = neredeyse imkânsız tesadüf (çok güçlü) |

**"Hata" değil, çoğu bir sonuç:** `evalue = 1e-56` = mükemmel hit. Yüksek E-değeri (ör. 5) = zayıf/anlamsız eşleşme. NCBI/CARD "error" görürsen genelde: (a) o unitig hiçbir bilinen gene uymuyor (yeni/novel olabilir), (b) organizmaya özgü panel yok (ör. AMRFinderPlus K. pneu için nokta-mutasyon paneli sunmaz — bu bir **uyarı**, hata değil).

> **Neden coverage bazen boş?** 2250 hit'in ~1235'inde coverage NULL — bu NCBI hit'lerinde coverage'ın her zaman hesaplanmaması ya da kısa yerel hizalamalardan. Eksik veri normaldir; E-değeri + identity zaten güveni verir (→ §13).

---

<a name="9"></a>
## 9. Güven seviyeleri: confirmed / candidate / weak / none

Her BLAST hit bir **güven seviyesine (`tier`)** atanır — identity + coverage eşiklerine göre (adım 09):

| tier | Anlam | KB'de sayı |
|---|---|---|
| `confirmed` | Yüksek-güven eşleşme (yüksek identity + coverage + düşük E) | 256 |
| `candidate` | Orta-güven, dikkate değer | 11 |
| `weak` | Zayıf eşleşme | 60 |
| `none` | Anlamlı hit yok (muhtemelen novel/bilinmeyen) | 1923 |

> **Önemli (reification / aşırı-yorum güvencesi):** "confirmed" **nedensellik iddiası değildir** — "bu unitig bilinen X geniyle yüksek-güvenle eşleşiyor" demektir. İstatistiksel sinyali biyolojik nedensellikten ayırmak için bu kademeli sistem kullanılır. `none` çokluğu (1923) doğaldır: unitig'lerin çoğu tek tek "bilinen gen" değil, sinyalin dağıldığı parçalardır; asıl önemli olan **kararlı + confirmed** olanlardır.

---

<a name="10"></a>
## 10. unitig_background_frequency nedir?

"Bu unitig gerçekten **dirençli** izolatlarda mı zenginleşmiş, yoksa her yerde mi var?" sorusunu ölçer (adım 10). Bir unitig önemli görünse bile, R ve S izolatlarda **eşit** sıklıktaysa ayırt edici değildir.

| Kolon | Anlam |
|---|---|
| `prevalence_resistant` | Dirençli (R) izolatlarda görülme oranı |
| `prevalence_susceptible` | Duyarlı (S) izolatlarda görülme oranı |
| `delta_prevalence` | R − S farkı (ne kadar ayırt edici) |
| `odds_ratio` | R'de bulunma bahsi (yüksek = R'ye özgü) |
| `fisher_p` | Fisher exact test p-değeri (anlamlılık) |
| `discriminative` | 0/1 — fark yeterince büyük **ve** p yeterince küçük mü |

**KB örneği:** `prevalence_resistant = 0.66`, `prevalence_susceptible = 0.05`, `odds_ratio = 36.9`, `fisher_p ≈ 0`, `discriminative = 1` → bu unitig dirençlilerin %66'sında, duyarlıların sadece %5'inde; 37 kat daha olası → **çok güçlü ayırt edici.**

---

<a name="11"></a>
## 11. 7 doğrulama katmanı — tek bir unitig örneği

`validation_evidence` tablosu, aynı unitig için farklı testlerin sonuçlarını toplar. Bir unitig ne kadar çok katmandan geçerse o kadar güvenilir. **KB'den gerçek bir unitig'in kanıt satırları:**

| evidence_type | evidence_source | evidence_score | Ne diyor |
|---|---|---|---|
| `blast` | CARD 4.0.1 | — | Bilinen bir gene eşleşiyor |
| `background_frequency` | R-vs-S Fisher exact | — | R'de zenginleşmiş |
| `permutation_mda` | MDA test ROC-AUC drop (100 perms, BH-FDR) | 0.0015 | Rastgele-permütasyondan anlamlı |
| `stability_selection` | CPSS (B=100, π≥0.6, PFER-bounded) | 1.0 | 100/100 kararlı |
| `pyseer_lmm` | pyseer LMM lineage-corrected (Bonferroni 1.09e-05) | — | Soy düzeltmesinden sonra da anlamlı |

Buna `snp` (varsa hedef-mutasyon) ve `label_permutation` (model-seviye) eklenir = **7 katman.** Bir unitig bu ortogonal testlerin çoğunu geçtiyse, tek bir soya bağlı gürültü değil, gerçek direnç sinyalidir. (`06_evidence_layers` figürü tüm modelleri bu 7 katmanla gösterir.)

---

<a name="12"></a>
## 12. External concordance — model vs AMRFinderPlus vs ResFinder

`external_concordance` = bizim modelimizi, sahadaki **standart genotip araçlarıyla** (AMRFinderPlus, ResFinder) kıyaslar. Üçü de **aynı held-out test genomlarında**, gerçek laboratuvar fenotipine (EUCAST/CLSI) göre puanlanır (adım 16).

| Kolon | Anlam |
|---|---|
| `balanced_accuracy` | Dengeli doğruluk (asıl kıyas) |
| `cohen_kappa` | Fenotiple uyum (0=şans, 1=tam) |
| `major_error_rate` (ME) | **Yanlış-dirençli** oranı (FDA) |
| `very_major_error_rate` (VME) | **Yanlış-duyarlı** oranı — klinik olarak **en tehlikeli** hata (hastaya işe yaramayan ilaç verdirir) |

**KB örneği — K. pneu ciprofloxacin (n=600 test genomu):**

| caller | bACC | κ | VME |
|---|---|---|---|
| **model (bizim)** | **0.926** | **0.853** | 0.094 |
| amrfinderplus | 0.538 | 0.074 | 0.0 |
| resfinder | 0.540 | 0.077 | 0.0 |

> **Yorum:** Quinolone direnci nokta-mutasyon (gyrA/parC) kaynaklı. Gen-tabanlı araçlar (AMRFinderPlus/ResFinder) edinilmiş gen ararlar → bu SNP'yi kaçırırlar (bACC ~0.54, rastgeleye yakın). **Unitig modelimiz mutasyonu yakalar (0.926).** Bu, sunumun en güçlü tek bulgusu. (`07_external_concordance` figürü)
>
> *(Not: AMRFinderPlus/ResFinder VME=0 çünkü neredeyse hiçbir izolatı "dirençli" demiyorlar — hepsine "duyarlı" derlerse yanlış-duyarlı olur ama onların ölçümünde farklı düşüyor; asıl gösterge düşük bACC/κ = araç bu mekanizmada başarısız.)*

---

<a name="13"></a>
## 13. "0 / None / tek değer" neden görünüyor?

Bunlar hata değil, **anlamlı boşluklar**:

- **`gain` boş ama `mean_abs_shap` dolu (veya tersi):** İki farklı seçim yöntemi var. `selection_method='gain_seed'` satırlarında **gain** dolu, SHAP boş; `selection_method='cpss'` satırlarında **SHAP** dolu, gain boş. Aynı unitig iki farklı yöntemle işlendiği için iki ayrı satırı olabilir. *(KB: gain_seed 1203 satır / 459'u stable; cpss 1047 satır / hepsi stable.)*
- **`coverage = None`:** Bazı BLAST hit'lerinde (özellikle NCBI) coverage hesaplanmaz; güven yine E-değeri + identity'den gelir.
- **`evidence_score = 0`:** Bazı kanıt türlerinde skor alanı kullanılmaz (ör. blast için asıl bilgi `blast_annotations`'ta; validation_evidence'ta sadece "bu katman geçti" işareti) → 0 "başarısız" değil, "skor bu katmanda taşınmıyor" demek.
- **`snp` satırı az (193) / bazı modellerde yok:** Sadece **nokta-mutasyon mekanizması olan** antibiyotiklerde (quinolone gyrA/parC) SNP satırı olur. Aminoglikozit/β-laktam edinilmiş-gen mekanizması olduğu için SNP tablosunda görünmez — **bu doğru biyoloji** (`antibiotics.mechanism_type` = `acquired` vs `target_snp` bunu ayırır).
- **`tier='none'` çokluğu (1923):** Unitig'lerin çoğu tek başına "bilinen gen" değil; asıl değerli olanlar **kararlı + confirmed** olan alt kümedir.

---

## Hızlı hatırlatma tablosu (soru → cevap)

| Soru | Kısa cevap |
|---|---|
| Unitig varken k-mer neden? | k-mer = QC aracı; unitig = model özelliği (BLAST'lanabilir) |
| GWAS? | Genom boyu direnç-ilişki taraması; pyseer LMM soy etkisini çıkarır |
| Provenance? | Kaydın kökeni (git+seed+CARD+config) → tekrarlanabilirlik |
| FAIR? | Bulunabilir/Erişilebilir/Birlikte çalışabilir/Yeniden kullanılabilir |
| CPSS / stable? | 100 alt-örnekte ≥%60 seçilen unitig = kararlı |
| SHAP? | Model-agnostik adil özellik-önemi ölçüsü |
| Hangi ROC? | `auc_mean_seeds` (lineage-CV) — dürüst genelleme |
| E-değeri? | BLAST eşleşmesinin tesadüf olma beklentisi (küçük=iyi) |
| confirmed? | Yüksek-güven eşleşme (nedensellik iddiası DEĞİL) |
| background_frequency? | R'de S'ye göre zenginleşme (ayırt edicilik) |
| external concordance? | Model vs AMRFinderPlus/ResFinder, aynı test genomlarında |
| 0/None neden? | Anlamlı boşluk (yöntem farkı / mekanizma yok / skor taşınmıyor) |

---

*Bu belge `results/kb/amrk.db` (şema 0.6.0) verilerinden gerçek örneklerle hazırlanmıştır. Tablo/kolon referansı: `docs/KB_ACIKLAMA.md`.*
