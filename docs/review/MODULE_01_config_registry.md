# Modül 1 — Config & Registry

> Mercek: **Python'da kal + sertleştir** · **ESKAPEE / antibiyotik-sınıfı odağı** · **su geçirmez, sıfır-hata, repodan tam reprodüksiyon**
> İncelenen: `config/config.yaml`, `config/registry/organisms.yaml`, `config/registry/antibiotics.yaml`, `scripts/lib/config.py`, `scripts/lib/registry.py`
> Tarih: 2026-07-13

---

## 1. Genel Değerlendirme

Bu modül pipeline'ın **tek doğruluk kaynağı** olması gereken katman: hangi organizma, hangi antibiyotik, hangi sınıf, hangi path. **Kod tarafı (config.py / registry.py) gerçekten iyi** — `get_target` (CLI > env > config precedence), `resolve_path` (organism-aware + feature_repr switch), `resolve_tool` (PATH önce, bundled binary sadece Darwin fallback), `registry.py` (lru_cache, temiz API, alias normalizasyonu). Bunlar düşünülmüş, test edilebilir kod.

**Ama veri tarafı gerçeğin gerisinde ve iki font-of-truth çatışıyor.** `organisms.yaml` KB'nin gerçek durumuyla **ciddi çelişki içinde** (kpneumoniae `enabled: false` iken elde 14 model var), ve organizma/antibiyotik biyolojik meta verisi (gram/phylum, mechanism_type, who_aware) registry'de DEĞİL, KB şemasında yaşıyor — yani "single source of truth" iddiası kısmen yanlış. Ek olarak, tezin asıl yöntemi (unitig) **her koşuda bir env override gerektiriyor** (`feature_repr` varsayılanı hâlâ `kmer`), bu da yanlışlıkla k-mer koşusu riski doğuran bir ayak-kapanı. ESKAPEE hedefi için registry hem eksik (3 patojen yok) hem de sınıf taksonomisi antibiyotik-sınıfı odaklı bir makale için fazla kaba.

Kod sağlam; **veri + varsayılanlar + tek-kaynak disiplini** düzeltilmeli.

---

## 2. Güçlü Yanlar

- **`config.py` mimarisi örnek:** `get_target` precedence (CLI > `AMR_ORGANISM`/`AMR_ANTIBIOTIC` env > config) geriye-uyumlu paralel koşuyu mümkün kılıyor; `resolve_path` `paths_organism` → legacy `paths` fallback + `{organism}/{antibiotic}/{run_id}` şablonları; `feature_repr` switch'i tek noktada (03b/04/05/06/07/07b hiç değişmeden unitig'e döner).
- **`resolve_tool` taşınabilir:** env override → PATH (conda/module) → sadece Darwin'de bundled binary. Yani `bin/bin/kmc` HPC'de hiç kullanılmıyor (M0'daki "vendored binary" endişesini yumuşatır — kod PATH'i tercih ediyor).
- **`registry.py` temiz:** `lru_cache`, `antibiotic_to_class` ters-index, `normalize_antibiotic` alias katlaması (typo→canonical), `list_targets`/`validate_target` API. Kopyalanan `ANTIBIOTIC_CLASSES` sözlüğünü tek dosyaya indirmiş (audit B09/B10).
- **9 ilaç sınıfı + alias + amrfinder_keywords** registry-güdümlü → yeni antibiyotik "kod değil, YAML" (doğru tasarım hedefi).
- **confidence_tiers config'de** (citeable, Methods'a konur).

---

## 3. Problemler (önem sırasıyla)

### Critical
- (yok) — env-parametrik koşu (`AMR_ORGANISM`/`AMR_ANTIBIOTIC`) stale registry'yi baypas ettiği için production sessizce kırılmadı. Ama aşağıdaki High'lar "su geçirmez" çıtası için bloklayıcı.

### High
- **H1 — `organisms.yaml` KB gerçeğiyle ağır çelişkili (single-source-of-truth yalanı):**
  - `ecoli.antibiotics: [ampicillin, cefotaxime, ciprofloxacin, gentamicin]` = **4**, ama KB'de **7** E. coli modeli var (+ trimethoprim/sulfa, ceftazidime, amoxicillin/clavulanic_acid).
  - `kpneumoniae: enabled: false` + `[meropenem, ciprofloxacin, gentamicin, colistin]` = 4 placeholder — ama KB'de **14** K. pneu modeli var (ve colistin hiç yapılmadı).
  - `staphylococcus_aureus: enabled: true` + 8 antibiyotik — ama **0 model** (henüz koşulmadı; enabled=true yanıltıcı).
  - Sonuç: `list_targets(enabled_only=True)` **yanlış küme** döndürür (14 kpneu'yu atlar, olmayan 8 staph'ı sayar). Registry "tek doğruluk kaynağı" değil.
- **H2 — Yöntem varsayılanı yanlış (reproducibility ayak-kapanı):** `preprocessing.feature_repr: "kmer"`. Tezin ASIL yöntemi unitig; onu koşmak için HER seferinde `AMR_FEATURE_REPR=unitig` override şart (HANDOFF defalarca "ŞART" diye uyarıyor). Override unutulursa sessizce **yanlış (k-mer baseline) koşu** olur. Varsayılan `unitig` olmalı, `kmer` opt-in baseline.

### Medium
- **M1 — Bölünmüş doğruluk kaynağı (biyolojik meta):** Organizma `gram_stain`/`phylum` ve antibiyotik `mechanism_type`(acquired/target_snp)/`who_aware`(AWaRe) **registry'de YOK** — yalnızca KB şemasında (0.5.0/0.6.0) var. Yani cross-phylum/ESKAPEE iddiasının biyolojik verisi registry'de değil KB'de. "Single source of truth" için bunlar registry'de curate edilip KB'ye AKMALI.
- **M2 — config.yaml sürüm/köken drifti:** `project.version: "M4-pro"` (anlamsız iç etiket), `provenance.kb_schema_version: "v0.1.0"` (gerçek 0.6.0), `provenance.card_version: null` (gerçek 4.0.1, env'den geliyor). M0'daki sürüm-drifti temasının devamı.
- **M3 — Legacy `paths:` bloğu köreldi:** `paths:` (yalnız `{antibiotic}`) ile `paths_organism:` (`{organism}/{antibiotic}`) iki ayrı düzen; `kmc_bin`/`data_dir` çift tanımlı. Migration bittiyse legacy blok vestigial → bakım tehlikesi. Tüm scriptler `resolve_path` kullanıyor mu M3'te doğrulanacak; kullanıyorsa legacy blok silinmeli.
- **M4 — Sınıf taksonomisi antibiyotik-sınıfı odağı için fazla kaba (bilimsel):** `beta_lactams_carbapenems_others` **karbapenemleri** (klinik olarak kritik, ayrı bir sınıf) aztreonam + genel "beta-lactam" ile lumpluyor. `others` bir çöp-torbası: polimiksinler (colistin, polymyxin B), glikopeptidler (vancomycin, teicoplanin), oksazolidinon (linezolid), fenikol (chloramphenicol) — hepsi farklı mekanizma. Sınıf-odaklı bir makale için bu granülerlik savunulamaz.
- **M5 — pseudomonas slug büyük-harf uyumsuzluğu:** registry key `Pseudomonas_aeruginosa`, ama `kb_figures._abbr` `paeruginosa`→"Pa" bekliyor (diğer sluglar küçük harf: `ecoli`, `kpneumoniae`, `staphylococcus_aureus`). Etkinleştirilirse slug eşleşmez → path/figure kırılır.

### Low
- **L1 — `normalize_antibiotic` bilinmeyeni sessizce geçiriyor:** typo/alias-dışı isim aynen döner → fenotip-etiketi ("extended spectrum beta lactamase") yanlışlıkla "yeni antibiyotik/target" olabilir. Bir `strict`/warn modu gerekli (su-geçirmez için).
- **L2 — config yorum typo'su:** `encoding: "binary"  # ... frequency-basedß` (başıboş `ß`).
- **L3 — Dangling doküman referansı:** `SCALE_MLOPS_PLAN.md` config.yaml/config.py/registry.py/environment.yml'de defalarca anılıyor — dosya repoda var mı? Yoksa kırık referans.
- **L4 — `lru_cache` registry:** çalışma-zamanı YAML mutasyonu cache'lenir; registry'yi değiştiren testler `cache_clear` gerektirir.
- **L5 — ESKAPEE eksik:** registry'de Enterococcus faecium, Acinetobacter baumannii, Enterobacter yok (tam ESKAPEE paneli için gerekli).

---

## 4. Düzeltilmesi Gerekenler (madde madde)

1. **`organisms.yaml`'ı gerçekle senkronla** (H1): ecoli `antibiotics` → gerçek 7; kpneumoniae `enabled: true` + gerçek 14 antibiyotik; saureus `enabled` durumunu netleştir (data yoksa `false` veya `planned` bayrağı); pseudomonas slug'ını `paeruginosa`'ya düzelt (M5).
2. **`feature_repr` varsayılanını `unitig` yap** (H2); `kmer`'i açık opt-in baseline olarak bırak; HANDOFF'taki "ŞART" uyarısı gereksizleşir.
3. **Biyolojik meta veriyi registry'ye taşı** (M1): organisms.yaml'a `gram_stain`/`phylum`; antibiotics.yaml'a `mechanism_type`/`who_aware`/(class'a) — populate bunları registry'den okusun (KB tek yerden beslensin).
4. **config.yaml sürümlerini düzelt** (M2): `project.version`'ı gerçek sürüme; `provenance.kb_schema_version` → 0.6.0; `card_version`'ı pinle ya da "env-driven" yorumu ekle.
5. **Legacy `paths:` bloğunu değerlendir** (M3): M3 modülünde tüm scriptlerin `resolve_path` kullandığı doğrulanınca sil; çift-tanımlı key'leri tekilleştir.
6. **Sınıf taksonomisini yeniden yapılandır** (M4): carbapenems ayrı sınıf; `others`'ı polymyxins / glycopeptides / oxazolidinones / phenicols olarak böl. *(literatür girdisi gerekiyor — §7)*
7. **`normalize_antibiotic`'e strict/warn modu** (L1); config typo (L2); `SCALE_MLOPS_PLAN.md` referanslarını doğrula/düzelt (L3).
8. **ESKAPEE panelini tamamla** (L5): E. faecium, A. baumannii, Enterobacter blokları. *(literatür girdisi — §7)*

---

## 5. Refactor Önerileri

- **Registry'yi tek gerçek kaynak yap:** organizma bloğu şeması = `{display_name, taxid, gram_stain, phylum, source, filter_criteria, antibiotics, status}`; antibiyotik = `{class, mechanism_type, who_aware, aliases, amrfinder_keywords}`. `populate_database.py` + `migrate_kb_050.py` bu alanları **registry'den** okusun (şu an muhtemelen hardcoded). Böylece KB, registry'nin türevi olur — çift-bakım biter.
- **`status` alanı** (`enabled` yerine): `planned | in_progress | done` → "enabled ama 0 model" belirsizliğini çözer; `list_targets` gerçek durumu yansıtır.
- **Class-first görünüm:** antibiyotik-sınıfı odağı için `list_targets_by_class()` gibi bir yardımcı (sınıf → [(organism, antibiotic)]) — makalenin sınıf-bazlı analizini registry'den sürmek için.
- **Config sadeleştirme:** legacy `paths:` kalkınca config.yaml ~%20 küçülür; `paths_organism` → tek `paths`.
- **Doğrulama scripti:** `scripts/validate_registry.py` (CI'da koşar) — registry ↔ KB tutarlılığı (her enabled target'ın KB'de modeli var mı, her KB modelinin registry'de kaydı var mı), slug lowercase, her member tekil sınıfta. "Su geçirmez" için bu bekçi çok değerli.

---

## 6. Bilimsel Eksikler (makale açısından)

- **ESKAPEE paneli eksik ve tanımsız:** "ESKAPEE odaklı" iddia için hedef patojen kümesi netleşmeli (E. faecium, S. aureus, K. pneumoniae, A. baumannii, P. aeruginosa, Enterobacter, + E. coli). Şu an sadece 2 tamam (ecoli, kpneu), 1 planlı (saureus). Hangi ESKAPEE üyeleri, hangi veri kaynağıyla (BV-BRC kapsamı)?
- **Antibiyotik-sınıfı taksonomisi savunulabilir olmalı:** makalenin ODAĞI sınıf iken, carbapenem'lerin lumplanması ve "others" çöp-torbası bir reviewer'ın ilk hedefi. Standart bir şemaya (WHO AWaRe + farmakolojik sınıf) oturtulmalı.
- **Organizma/ilaç biyolojik meta (gram/phylum, mechanism_type, AWaRe) curate + kaynaklı olmalı** — cross-phylum ve mechanism-showcase iddialarının temeli; registry'de referanslarıyla durmalı.

---

## 7. Literatür Gereksinimleri — ÇÖZÜLDÜ (docs/literature/ESKAPEE1.md, 2026-07-13)

İki bilimsel karar literatür raporuyla netleşti:

**Karar 1 — ESKAPEE paneli (faz'lı):**
- **Faz 1** (tez + 1. SCI makale): **K. pneumoniae, E. coli, S. aureus, A. baumannii** — en yüksek WGS+AST hacmi, en net genotip-fenotip.
- **Faz 2** (2. makale / optimizasyon): **P. aeruginosa** (polijenik/epistatik, zor), **E. faecium** (glikopeptid, orta veri), **Enterobacter** (taksonomik kirlilik, düşük veri).
- Veri-uygunluğu artık KAPI değil (hepsi BV-BRC'den çekiliyor); faz sadece önceliklendirme metadata'sı.

**Karar 2 — Sınıf taksonomisi:**
- Temel = **CARD ARO ontolojisi** (ML AMR altın standardı); fenotip normalizasyonu **NCBI Antibiogram + CLSI/EUCAST**. (WHO AWaRe/ATC/MeSH ML için uygunsuz — mekanizma bazlı değil.)
- **Carbapenemler AYRI sınıf** (beta-laktam altında lumplanmaz → aksi halde blaKPC/blaNDM sinyali genel β-laktamaz gürültüsünde kaybolur, model çöker).
- **"Others" tamamen dağıtılır:** glycopeptides / polymyxins / oxazolidinones / phenicols / rifamycins / fosfomycins / lipopeptides / nitrofurans ayrı mekanizma sınıfları; verisi yetmeyen sınıf **kapsam dışı**, asla lumplanmaz.

Bu kararlara göre `antibiotics.yaml` + `organisms.yaml` yeniden yapılandırma önerisi hazır (aşağıdaki "Önerilen değişiklikler" — sohbette sunuldu, onay bekliyor). Kalan düzeltmeler (H1 senkron, H2 varsayılan, M2/M3/M5, L'ler) literatür beklemiyordu.

### Önerilen değişiklikler (uygulama onay bekliyor)
- **antibiotics.yaml:** 10 sınıf → ~19 mekanizma-bazlı sınıf (carbapenems, monobactams, glycopeptides, polymyxins, oxazolidinones, phenicols, rifamycins, fosfomycins, lipopeptides, nitrofurans, glycylcyclines ayrıştı; `others` silindi; `quinolones`→`fluoroquinolones`).
- **organisms.yaml:** yeni şema alanları `gram_stain`, `phylum`, `eskapee_phase`, `status` (done|planned|not_started), `priority_classes`; `enabled` yerine `status`. ecoli→gerçek 7 abx, kpneumoniae→gerçek 14 + status:done; slug'lar lowercase; +A. baumannii/P. aeruginosa/E. faecium/Enterobacter blokları.
- **Sıralama uyarısı:** taksonomi değişikliği M7'ye (populate/KB `antibiotics` tablosu, `antibiotic_to_class`) sıçrar → registry + populate + `validate_registry.py` **birlikte** güncellenmeli (yoksa geçici tutarsızlık).

---

## Uygulama durumu (2026-07-13) — UYGULANDI

M1 değişiklikleri onayla uygulandı (hepsi şimdi):
- `antibiotics.yaml` → schema 2.0, mekanizma-bazlı 19 sınıf (carbapenems/monobactams/glycopeptides/polymyxins/oxazolidinones/phenicols/rifamycins/fosfomycins/lipopeptides/nitrofurans/glycylcyclines ayrı; `others` silindi; `beta-lactam`+`sulbactam` unclassified).
- `organisms.yaml` → schema 2.0, `status`/`eskapee_phase`/`gram_stain`/`phylum`/`priority_classes`; ecoli→7, kpneumoniae→14 (status: done); +A. baumannii (Faz1) +P. aeruginosa/E. faecium/Enterobacter (Faz2); sluglar lowercase.
- `registry.py` → `list_targets` `status`-tabanlı (+`phase` filtresi, `is_active` helper), geriye-uyumlu.
- `config.yaml` → `feature_repr: unitig` (varsayılan), `kb_schema_version 0.6.0`, `card_version 4.0.1`, `version 0.6.0`, typo.
- `15_cross_antibiotic.py` → `_BETA_LACTAM_CLASSES`: `beta_lactams_carbapenems_others` → `carbapenems`.
- **KB re-sync** (Mac `results/kb/amrk.db`): meropenem/imipenem `drug_class` → carbapenems (2 satır); `same_class` etkilenmedi.
- Testler: `test_lib` (2) + `test_15` `_drug_family` testi güncellendi.

**M1 tam bitirme (single source of truth + watertight):**
- `antibiotics.yaml` → `class_mechanism_type` (acquired/target_snp, per-class) + `who_aware` (WHO AWaRe, per-antibiotic) eklendi → biyolojik meta artık registry'de (KB'de değil).
- `registry.py` → `antibiotic_mechanism_type` / `antibiotic_who_aware` / `clear_cache` accessor'ları.
- `populate_database.py` → hardcoded `_ORG_META`/`_AWARE`/`_MECH` **silindi**; `populate_organisms` + `populate_antibiotics_meta` artık registry'den okuyor (silinmiş `beta_lactams_carbapenems_others` referansı da böylece düzeldi).
- **`scripts/validate_registry.py`** (YENİ) — watertight bekçi: registry iç tutarlılık (üye tekilliği, alias↔member, slug lowercase, priority_class geçerliliği, mechanism/aware değer aralığı) + `--db` ile registry↔KB senkron. **Çıktı: 0 hata / 0 uyarı** (hem iç hem KB).
- Testler: `test_registry_metadata_accessors` + carbapenems assertion eklendi (`test_lib` 15/15 yeşil).

**Test durumu:** M1 değişiklikleri **0 yeni kırık** ekledi (97 passed). Kalan kırıklar M1-ÖNCESİ (stash ile HEAD'de doğrulandı):

### ⚠️ Ön-var test kırıkları (bu session'dan ÖNCE kırıktı — HANDOFF "117 pass" bayat)
- **`test_15_cross_antibiotic.py` (3 FAIL)** — `populate_overlap() missing 'logger'`: script 15'in `populate_overlap` imzası `logger` aldı, test güncellenmedi. → **M6**.
- **`test_kb_queries.py` (6 ERROR)** — `NOT NULL: unitig_antibiotic_overlap.organism`: 0.6.0 organism-aware overlap kolonu eklendi, test fixture'ı güncellenmedi. → **M7 (şema) / M10 (test)**.
- **KB re-sync kalıcılığı:** Mac KB güncellendi; **TRUBA/Drive KB kopyaları drug_class'ta bayat** (regenerate/populate ile senkron olacak — M7).

## Sonraki modüllere taşınan notlar
- **M3 (QC):** legacy `paths:` bloğunu kimse kullanıyor mu? (kullanılmıyorsa sil). Organizma-seviyesi scriptler (01/02/02b/02p/03) `get_target`/`resolve_path` mi yoksa config'i doğrudan mı okuyor (ESKAPEE paralel koşu için).
- **M7 (KB):** `mechanism_type`/`who_aware`/`gram`/`phylum` populate'te nereden geliyor (hardcoded mı registry mi)? → registry'ye taşınınca populate güncellenir.
- **M10 (Test/CI):** `validate_registry.py` bekçisini CI'a ekle.
- **Genel:** `SCALE_MLOPS_PLAN.md` var mı — yoksa tüm referanslar kırık (M0/M11 doküman denetimi).
