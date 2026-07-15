# E4 — ARAŞTIRMA İSTEĞİ: ESKAPEE için laboratuvar-doğrulamalı AMR fenotip veri kaynakları

> Tarih: 2026-07-15 · İstek sahibi: proje · Format: E1/E2/E3 ile aynı (deep research)
> **Bu bir varsayım listesi DEĞİL — cevaplanması gereken sorular listesi.** Aşağıdaki
> sayılar bizim ölçtüklerimiz; kaynak önerileri araştırmadan çıkacak.

---

## 0. Bağlam — neden soruyoruz

Tez paneli **ESKAPEE (7 organizma)**. Veri kaynağı bugüne kadar **yalnızca BV-BRC**
(`genome_amr` tablosu, API). 2026-07-15'te panelin tamamı için BV-BRC kapsamını ölçtük:

| organizma | taxid | **lab-doğrulamalı** satır | computational satır | evidence-boş satır | model-uygun genom |
|---|---|---|---|---|---|
| E. coli | 562 | 243.124 | 6.894.231 | 3.067 | 5.470 |
| K. pneumoniae | 573 | 85.291 | 1.728.894 | 26.608 | 4.615 |
| S. aureus | 1280 | 45.876 | 473.467 | 6.298 | 2.494 |
| A. baumannii | 470 | 28.237 | 456.826 | 2.226 | 1.126 |
| P. aeruginosa | 287 | 11.121 | 308.382 | 1.496 | 1.312 |
| E. faecium | 1352 | 25.635 | 42.390 | 20 | 2.270 |
| **Enterobacter cloacae** | 550 | 11.935 | 15.911 | 62 | **320** |

**Kritik ayrım:** `evidence="Computational Method"` kayıtları (her organizmada ezici
çoğunluk) **kullanılamaz** — bunlar bir yazılımın genomdan çıkardığı tahminler. ML
etiketi olarak kullanmak döngüsel olur ve projenin M13 iddiasını (modelimiz
AMRFinderPlus'tan iyi: bACC 0.926 vs 0.538) anlamsızlaştırır. Sadece gerçek AST
(CLSI/EUCAST, MIC/broth dilution/disk diffusion) kullanılabilir.

**Antibiyotik başına azınlık sınıfı (minority = min(R,S)), eşik ≥150:**
- E. faecium: 8 antibiyotik (vancomycin 931 — VRE bayrak hedefi)
- A. baumannii: 7 (imipenem 400, meropenem 185 — carbapenem)
- P. aeruginosa: 5 (meropenem 404, ceftazidime 368)
- **Enterobacter cloacae: 0** (en yüksek minority **44**) → mevcut kaynakla model eğitilemez

**Karar:** Enterobacter tezde **negatif bulgu** olarak raporlanacak. Bu araştırma,
o kararın (ve genel olarak "tek kaynak yeterli mi" sorusunun) sağlamasını yapacak.

---

## 1. ANA SORU

**ESKAPEE patojenleri için, bir genom assembly'sine bağlı laboratuvar-doğrulamalı
AMR fenotipi (AST) içeren, BV-BRC dışındaki hangi açık kaynaklar var; kapsamları
BV-BRC'nin üzerine ne katar; ve mükerrer izolatlar nasıl ayıklanır?**

---

## 2. İNCELENMESİ İSTENEN KAYNAKLAR (liste kapalı değil)

Her biri için: **(a)** ESKAPEE başına gerçek-AST'li izolat sayısı, **(b)** assembly/read
erişimi var mı, **(c)** programatik indirme (API/FTP) var mı, **(d)** lisans ve
yeniden-dağıtım (Zenodo deposit'e konabilir mi), **(e)** BV-BRC ile örtüşme oranı.

- **NCBI Pathogen Detection** — AST tabloları (`Antibiogram`) ne kadar kapsıyor? Isolates
  browser'daki AST verisi toplu indirilebiliyor mu? BioSample'a bağlanıyor mu?
- **NCBI BioSample `antibiogram` / AST paket alanları** — kaç ESKAPEE BioSample'ında dolu?
- **EnteroBase** (Escherichia/Salmonella/Klebsiella) — fenotip taşıyor mu, yoksa sadece
  genom+MLST mi?
- **PubMLST / BIGSdb** — AMR fenotip alanları hangi türlerde dolu?
- **CARD Prevalence / Resistomes & Variants** — fenotip mi, tahmin mi?
- **AllTheBacteria / 661k genome collection** — assembly bol ama fenotip var mı?
- **NARMS (FDA/CDC)**, **EARS-Net/EUCAST surveillance**, **ATLAS (Pfizer)**,
  **SENTRY (JMI)**, **GEARS (Merck)** — genom + AST birlikte açık mı, yoksa sadece
  toplu istatistik mi?
- **Enterobacter'e özel** herhangi bir koleksiyon/çalışma var mı? (320 genom çıkmazını
  aşacak bir kaynak varsa, tezin panelini 7'ye tamamlar.)
- Büyük **kurumsal/ulusal koleksiyonlar** (ör. NCTC 3000, Wellcome/Sanger koleksiyonları,
  CRyPTIC benzeri konsorsiyumlar) — ESKAPEE muadili var mı?

---

## 3. MÜKERRER İZOLAT (DEDUP) — bu kritik

BV-BRC büyük ölçüde **NCBI'ı toplayan** bir kaynak (PATRIC mirası). Yani ikinci bir
kaynaktan çekeceğimizin **çoğu muhtemelen zaten elimizde**. Aynı izolatı iki kez almak
veri sızıntısı yaratır (aynı genom hem train hem test'e düşebilir) ve lineage-CV'yi
bozar — projenin tüm bilimsel omurgası buna dayanıyor.

- Kaynaklar arası **ortak anahtar** ne? (`biosample_accession`, `assembly_accession`
  GCA_/GCF_, SRA run, BioProject) — hangisi ne kadar güvenilir dolu?
- Aynı izolatın **farklı assembly'si** iki kaynakta varsa nasıl tespit edilir?
  İçerik-bazlı dedup (Mash/ANI eşiği?) literatürde nasıl yapılıyor, hangi eşik kabul görüyor?
- **Yayımlanmış çalışmalar** çok-kaynaklı veri setlerini nasıl dedup ediyor? Somut
  protokol arıyoruz, genel tavsiye değil.

---

## 4. ETİKET UYUMLAŞTIRMA (label harmonisation)

Farklı kaynaklar farklı standart/eşik kullanıyor. Bizde şu an `intermediate_policy: drop`
ve BV-BRC'nin `resistant_phenotype` alanı kullanılıyor; ham MIC'e inmiyoruz.

- **CLSI vs EUCAST** breakpoint farkları R/S etiketini ne kadar değiştiriyor? Aynı MIC
  farklı standartta farklı sınıfa düşebiliyor mu — bu ML etiketi için ne kadar gürültü?
- **Breakpoint yılı** (`testing_standard_year`) kayıyor — 2015 CLSI ile 2023 CLSI aynı
  MIC'i farklı etiketliyor olabilir. Literatür bunu nasıl ele alıyor? Ham MIC'ten
  **tek bir standartla yeniden etiketleme** yapılıyor mu, yapılmalı mı?
- **Intermediate / "Susceptible-dose dependent" / "Not defined"** kategorileri: drop mu,
  R'ye mi katılıyor? Güncel en iyi uygulama ne?
- Ham **MIC değerini** (measurement + sign + unit) kullanıp kendi eşiğimizi uygulamak
  daha savunulabilir mi? (BV-BRC bu alanları taşıyor.)

---

## 5. PRESEDANS — belki en önemli soru

**2021-2026 arası WGS-AMR ML makaleleri veri setlerini nereden alıyor?**
- Sadece PATRIC/BV-BRC kullananların oranı ne? (Eğer alan standardı buysa, tek kaynak
  kullanmamız savunulabilir ve bu araştırma "kontrol edildi, gerek yok" diye kapanır.)
- Çok-kaynaklı olanlar hangi kaynakları birleştiriyor, dedup'ı nasıl anlatıyorlar?
- **Organizma başına tipik n** ne? Bizim A. baumannii 1126 / P. aeruginosa 1312
  sayılarımız yayımlanmış çalışmalara göre düşük mü, normal mi? (Hakem "veriniz az"
  derse cevabımız olmalı.)
- Bir organizmayı **veri yetersizliğinden dışlamak** (bizim Enterobacter kararımız)
  literatürde nasıl raporlanıyor? Kabul gören bir dil/format var mı?

---

## 6. HANGİ KARARLARI ETKİLER

1. **Enterobacter panelde kalır mı?** Başka kaynak onu kurtarıyorsa 7 organizma; yoksa
   6 + negatif bulgu.
2. **İkinci kaynak eklenecek mi?** Eklenecekse: hangisi, hangi organizma için, dedup
   protokolü ne, ve `00a`'ya nasıl entegre edilir (yeni bir backend mi?).
3. **Ham MIC'ten yeniden etiketleme yapılacak mı?** Yapılacaksa hangi standart/yıl.
4. **Mevcut n'ler savunulabilir mi?** (A. baumannii 1126 ile carbapenem modeli — hakem
   ne der?)

---

## 7. NE İSTEMİYORUZ

- Computational/in-silico tahmin içeren kaynaklar (AMRFinder/ResFinder çıktısı fenotip
  yerine geçemez — projenin karşılaştırma iddiasını yok eder).
- Genom olmadan sadece AST içeren sürveyans özetleri (ATLAS benzeri toplu istatistikler
  ML için kullanılamaz — izolat düzeyinde genom + fenotip eşleşmesi şart).
- "Kayıt istenirse erişilir" tipi kapalı kaynaklar (FAIR/yeniden-üretilebilirlik iddiası
  ve Zenodo deposit'i ile bağdaşmaz).
