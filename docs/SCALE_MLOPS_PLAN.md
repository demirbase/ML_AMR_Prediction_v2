# Çok-Organizmalı Ölçeklenme & MLOps Hazırlık Planı

**Proje:** ML AMR Prediction Framework v2 → AMRK-DB
**Amaç:** Mevcut tek-organizma (E. coli) / antibiyotik-bazlı pipeline'ı, **çok organizma + her organizmanın kendi antibiyotik seti + veritabanı (KB) + MLOps** yapısına *geriye uyumlu ve kademeli* şekilde hazırlamak.
**Yaklaşım:** Mevcut scriptler ve `config.yaml` çalışmaya devam eder; yeni katmanlar üstüne eklenir. Hiçbir adım mevcut çıktıları kırmaz.
**İlke:** Önce yapısal hazırlık (klasör + kimlik + registry + run metadata), sonra DB, sonra MLOps otomasyonu.

> Bu bir **plan dökümanıdır** — kod/klasör henüz oluşturulmadı. Onayladıktan sonra fazları sırayla uygularız.

---

## 0. Mevcut Durumun Kısa Teşhisi (neden refactor gerekli)

Çok-organizmaya geçişi bugün engelleyen yapısal noktalar:

1. **Kimlik tekil:** `config.yaml → project.target_antibiotic` ve `organism` tek string. Her şey tek organizma varsayıyor.
2. **Path şeması yalnız `{antibiotic}`:** `matrix_dir: data/processed/{antibiotic}/matrix`. Organizma boyutu yok → iki organizmada aynı antibiyotik (örn. gentamicin) çakışır.
3. **KMC çıktıları global:** `data/interim/global_kmc_outputs/` tüm antibiyotiklerde paylaşılıyor. Bu *tek organizma için* doğru ama **farklı organizmaların genomları karışır**.
4. **Metadata tek dosya:** `data/external/metadata/genome_amr_matrix.csv` — tüm antibiyotikler kolon. Farklı organizmanın farklı antibiyotik seti ve farklı genom kümesi var.
5. **Antibiyotik sınıfları kopyalanmış:** `ANTIBIOTIC_CLASSES` sözlüğü `01` ve `01b` içinde birebir tekrar (denetimde B09). Tek kaynak (registry) yok.
6. **Run metadata yok:** git commit, seed, KMC/CARD versiyonu, çalışma tarihi hiçbir yerde kalıcı değil → reproducibility + KB provenance imkânsız (Roadmap M6, M10, M11).
7. **Orkestrasyon manuel:** 9 script elle, sırayla çalışıyor; organizma/antibiyotik parametresi yok.

Hedef: bu 7 noktayı **kırmadan** kapatmak.

---

## 1. Hedef Dizin Yapısı (çok-organizma)

Mevcut `{antibiotic}` hiyerarşisinin **bir üstüne `{organism}`** ekliyoruz. Organizma = kısa, kararlı bir slug (örn. `ecoli`, `kpneumoniae`, `paeruginosa`, `saureus`).

```
ML_project_kopyasi/
├── config/
│   ├── config.yaml                     # GLOBAL varsayılanlar (organizmadan bağımsız)
│   ├── registry/
│   │   ├── organisms.yaml              # YENİ — organizma kayıt defteri (tek kaynak)
│   │   └── antibiotics.yaml            # YENİ — antibiyotik→sınıf eşlemesi (tek kaynak; B09/B10 çözümü)
│   └── experiments/
│       └── {organism}/
│           └── config_{antibiotic}.yaml   # 04_optimization.py'nin ürettiği auto-config (taşındı)
│
├── data/
│   ├── raw/
│   │   └── {organism}/genomes/*.fna          # organizma-bazlı ham genomlar
│   ├── external/
│   │   ├── {organism}/metadata/amr_phenotypes.csv   # organizma-bazlı fenotip matrisi
│   │   └── blast_db/card_nt/                 # PAYLAŞILAN (organizmadan bağımsız)
│   ├── interim/
│   │   └── {organism}/kmc_outputs/           # organizma-bazlı KMC DB (global yerine)
│   └── processed/
│       └── {organism}/{antibiotic}/matrix/   # X_*.npz, y_*.csv, genomes_*.csv, features.txt
│
├── models/
│   └── {organism}/{antibiotic}/
│       ├── xgboost_{antibiotic}_final.json
│       ├── optuna_study_{antibiotic}.pkl
│       └── manifest.json                     # YENİ — model kartı (metrik, param, veri hash, git)
│
├── results/
│   └── {organism}/{antibiotic}/01_..05_...   # mevcut yapı korunur
│
├── logs/
│   └── {organism}/{antibiotic}/NN.log
│
├── runs/                                      # YENİ — MLOps run metadata (her pipeline koşusu)
│   └── {organism}/{antibiotic}/{run_id}/
│       ├── run_metadata.json                  # git hash, seed, sürümler, tarih, parametreler
│       └── metrics.json
│
├── db/                                        # YENİ — Knowledge Base katmanı
│   ├── schema/
│   │   ├── 001_core.sql                       # antibiotics, organisms, models, kmers...
│   │   ├── 002_provenance.sql                 # pipeline_runs, validation_evidence
│   │   └── 003_indexes.sql
│   ├── migrations/                            # şema sürüm geçişleri (Alembic veya düz SQL)
│   ├── amrkb.sqlite                           # geliştirme DB (Ay 3) → Postgres'e migrate (Ay 4)
│   └── README.md
│
├── scripts/
│   ├── lib/                                   # YENİ — paylaşılan modüller (kopya kod sonu)
│   │   ├── __init__.py
│   │   ├── config.py                          # registry + config yükleyici, path çözücü
│   │   ├── registry.py                        # organisms.yaml / antibiotics.yaml erişimi
│   │   ├── chunking.py                        # get_y_chunk (tek kopya)
│   │   ├── run_metadata.py                    # git hash, sürüm yakalama, run_id üretimi
│   │   └── io_utils.py                        # run_command (güvenli), npz yükleme yardımcıları
│   ├── db/                                    # YENİ — KB katmanı kodu
│   │   ├── models.py                          # SQLAlchemy ORM (SQLite & Postgres ortak)
│   │   ├── 11_populate_database.py            # k-mer skorları + BLAST → KB
│   │   └── 12_export_fair_metadata.py         # /api/v1/metadata JSON üreteci
│   ├── 01_..09_...                            # mevcut scriptler (küçük düzenlemelerle)
│   ├── 07b_feature_stability.py               # YENİ (Roadmap M4) — 5-seed stability
│   ├── 10_cross_antibiotic_analysis.py        # YENİ (Roadmap S1)
│   └── run_pipeline.py                        # YENİ — orkestratör (--organism --antibiotic)
│
├── tests/                                     # YENİ — birim testler (lib + DB)
├── environment.yml / requirements.txt         # sabitlenmiş sürümler
├── Makefile veya dvc.yaml                      # YENİ — tekrarlanabilir komutlar
└── SCALE_MLOPS_PLAN.md                         # bu dosya
```

**Geriye uyumluluk notu:** Yeni `{organism}` katmanı eklenirken, mevcut E. coli verisi `ecoli` slug'ı altına taşınır (Faz 1, Bölüm 8). Path şablonları `config.yaml`'de değiştiği için scriptler kod değişmeden yeni düzeni okur (bkz. Bölüm 4).

---

## 2. Kimlik & İsimlendirme Standardı

Her şeyin temeli kararlı kimliklerdir. KB'de FK'ler ve dosya yolları bunlara dayanır.

| Varlık | Kimlik (slug) | Örnek | Kural |
|---|---|---|---|
| Organizma | `organism_id` | `ecoli`, `kpneumoniae` | küçük harf, boşluksuz, tür kısaltması; değişmez |
| Antibiyotik | `antibiotic_id` | `gentamicin`, `cefotaxime` | küçük harf; `/` → `_` (örn. `trimethoprim_sulfamethoxazole`) |
| Antibiyotik sınıfı | `class_id` | `aminoglycosides` | registry'den |
| Model koşusu | `run_id` | `ecoli__gentamicin__20260610T1432__a1b2c3d` | `{org}__{ab}__{UTC}__{git7}` |
| KB şema sürümü | `kb_schema_version` | `v0.1.0` | semantic versioning |

`run_id` formatı sayesinde bir model dosyası, log, sonuç ve KB kaydı tek bir koşuya kadar izlenebilir (Roadmap Risk 3 yanıtı).

---

## 3. Registry Sistemi (yeni — tek kaynak)

Çok-organizmanın kalbi: hangi organizmanın hangi antibiyotikleri, hangi veri yolları, hangi metadata kolonu olduğunu **tek yerde** tanımlamak. Bu, kopya `ANTIBIOTIC_CLASSES` sözlüklerini de ortadan kaldırır (B09/B10).

### 3.1 `config/registry/organisms.yaml`

```yaml
# Organizma kayıt defteri — yeni organizma eklemek = buraya bir blok eklemek
schema_version: "1.0"
organisms:
  ecoli:
    display_name: "Escherichia coli"
    taxid: 562
    metadata_file: "data/external/ecoli/metadata/amr_phenotypes.csv"
    genome_id_column: "Genome ID"
    genomes_dir: "data/raw/ecoli/genomes"
    source: "BV-BRC"
    download_date: "2026-06-01"
    filter_criteria: "N50>=20000; max_contigs<=500; species=Escherichia coli"
    # Bu organizma için ML hedefi olan antibiyotikler (alt küme):
    antibiotics: [ampicillin, cefotaxime, ciprofloxacin, gentamicin]
    enabled: true

  kpneumoniae:               # GELECEK — şablon
    display_name: "Klebsiella pneumoniae"
    taxid: 573
    metadata_file: "data/external/kpneumoniae/metadata/amr_phenotypes.csv"
    genome_id_column: "Genome ID"
    genomes_dir: "data/raw/kpneumoniae/genomes"
    source: "BV-BRC"
    download_date: null
    filter_criteria: null
    antibiotics: [meropenem, ciprofloxacin, gentamicin, colistin]
    enabled: false           # veri gelene kadar kapalı
```

Yeni organizma eklemek = bu dosyaya bir blok + veri klasörlerini doldurmak. **Kod değişmez.**

### 3.2 `config/registry/antibiotics.yaml`

`01` ve `01b`'deki kopyalanmış `ANTIBIOTIC_CLASSES` buraya taşınır; her iki script bunu import eder.

```yaml
schema_version: "1.0"
classes:
  penicillins:
    display_name: "Penicillins"
    members: [ampicillin, amoxicillin, "amoxicillin/clavulanic acid", penicillin, ...]
  cephalosporins:
    display_name: "Cephalosporins"
    members: [cefotaxime, ceftazidime, ceftriaxone, cefepime, ...]
  aminoglycosides:
    members: [gentamicin, amikacin, tobramycin, streptomycin, ...]
  quinolones:
    members: [ciprofloxacin, levofloxacin, norfloxacin, "nalidixic acid", ...]
  # ... mevcut 8 sınıf birebir taşınır
# Yardımcı: antibiyotik → sınıf ters indeksi lib/registry.py içinde otomatik üretilir.
```

### 3.3 `scripts/lib/registry.py` (yeni)

Sağlayacağı fonksiyonlar:
- `load_organisms() -> dict`
- `load_antibiotic_classes() -> dict`
- `antibiotic_to_class(ab_id) -> class_id` (ters indeks)
- `list_targets() -> [(organism_id, antibiotic_id), ...]` (enabled olanlar)
- `validate_target(organism_id, antibiotic_id)` (registry'de var mı?)

Böylece `01_data_validation.py` ve `01b` artık sınıf sözlüğünü `registry.load_antibiotic_classes()`'tan alır — **tek kaynak**, kopya yok.

---

## 4. `config.yaml` Düzenlemeleri (geriye uyumlu)

Mevcut `config.yaml`'i **silmeden**, organizma boyutunu ve registry referanslarını ekliyoruz. Path şablonlarına `{organism}` ekleniyor.

### 4.1 Eklenecek/değişecek alanlar

```yaml
project:
  organism: "ecoli"            # DEĞİŞTİ: artık registry slug'ı (eski "Escherichia coli" yerine)
  target_antibiotic: "gentamicin"
  # ... mevcut alanlar korunur

registry:                      # YENİ
  organisms_file: "config/registry/organisms.yaml"
  antibiotics_file: "config/registry/antibiotics.yaml"

paths:
  # GLOBAL/paylaşılan (değişmez):
  blast_db_dir: "data/external/blast_db/card_nt"

  # ORGANİZMA-BAZLI (şablona {organism} eklendi):
  genomes_dir:      "data/raw/{organism}/genomes"                  # eski raw_genomes_dir
  metadata_file:    "data/external/{organism}/metadata/amr_phenotypes.csv"
  kmc_outputs_dir:  "data/interim/{organism}/kmc_outputs"          # eski global_kmc_outputs

  # ORGANİZMA + ANTİBİYOTİK-BAZLI ({organism}/{antibiotic}):
  matrix_dir:           "data/processed/{organism}/{antibiotic}/matrix"
  models_dir:           "models/{organism}/{antibiotic}"
  logs_dir:             "logs/{organism}/{antibiotic}"
  analysis_results_dir: "results/{organism}/{antibiotic}"
  experiment_config:    "config/experiments/{organism}/config_{antibiotic}.yaml"  # auto-config taşındı
  run_dir:              "runs/{organism}/{antibiotic}/{run_id}"    # YENİ — MLOps

provenance:                    # YENİ — KB/reproducibility için (Roadmap M6, M10)
  card_version: null           # örn. "CARD v3.2.9" — BLAST öncesi doldurulmalı
  kb_schema_version: "v0.1.0"
  entrez_email: "eren.demirbas@<kurum>.edu.tr"   # 09'daki user@example.com yerine (B07)
  ncbi_api_key: null
```

### 4.2 Path çözüm yardımcısı — `scripts/lib/config.py` (yeni)

Tek bir fonksiyon tüm scriptlerin path kurma kodunu merkezîleştirir:

```python
def resolve_path(key, organism=None, antibiotic=None, run_id=None) -> Path: ...
def load_config() -> dict: ...        # global config.yaml
def get_target() -> (organism, antibiotic):   # CLI arg > env var > config.yaml
```

**Geriye uyumluluk:** `get_target()` önce `--organism/--antibiotic` CLI argümanına, yoksa `config.yaml`'deki tekil değerlere düşer. Yani mevcut "config'i elle değiştir, scripti çalıştır" akışı aynen çalışır; ek olarak parametreli çağrı mümkün hale gelir.

---

## 5. Mevcut Scriptlerde Yapılacak Küçük Düzenlemeler

Hepsi **geriye uyumlu**; davranış aynı kalır, yalnız organizma boyutu + tek kaynak + run metadata kazanır.

| Script | Düzenleme | Gerekçe |
|---|---|---|
| `lib/*` (yeni) | `get_y_chunk`, `run_command`, config/path yükleme buraya taşınır | Kopya kod sonu (B01 güvenli `run_command` tek yerde) |
| `01_data_validation.py` | `ANTIBIOTIC_CLASSES` → `registry.load_antibiotic_classes()`; metadata yolu `{organism}` ile çözülür; **çift sayım bug'ı düzeltilir** (`get(1.0)+get(1)`) | B09 + denetimdeki kritik çift-sayım |
| `01b_data_validation.py` | Aynı registry import'u; kopya sözlük silinir | B09 |
| `02_kmer_extraction.py` | KMC çıktısı `{organism}/kmc_outputs`; `MIN_COUNT` config'ten; resume kontrolü `.kmc_pre` **ve** `.kmc_suf` | Organizma izolasyonu + bozuk DB resume |
| `02b/03b` (QC görseller) | Path'ler `{organism}/{antibiotic}`; `print("\n="*60)` typo düzeltilir | İzolasyon + kozmetik |
| `03_matrix_construction.py` | Path'ler `{organism}`; `run_command` lib'den (shell=True kaldırılır, B01); KMC `-ci` "prevalence" yorumu netleştirilir | İzolasyon + güvenlik |
| `04_optimization.py` | Auto-config `config/experiments/{organism}/`; `eval_metric` config'ten; metrik etiketi düzeltilir; **run_metadata.json yazılır** | MLOps provenance |
| `05_model_training.py` | Path'ler `{organism}`; model adı sabit-`v2` yerine `run_id`; `manifest.json` yazılır; config'i in-place ezme yerine ayrı `evaluation` bloğu | Model registry |
| `06_evaluation.py` | Youden's J **test'ten** kaldırılır → train/val'den (M1, data leakage); metrikler `runs/.../metrics.json`'a | Submission blocker fix |
| `07_explainability.py` | Path `{organism}`; sabit "ciprofloxacin" metni dinamikleştirilir; `02_top_{TOP_N}` dosya adı 08/09 ile hizalanır | Tutarlılık |
| `08_blast_annotation.py` / `.nf` | FASTA adı `{TOP_N}` ile dinamik (sabit `_50_` kaldır); `card_version` config'ten loglanır; NCBI process'e `errorStrategy 'ignore'` | Denetim bulguları |
| `09_biological_summary.py` | `Entrez.email`/`api_key` config'ten (B07); dosya adı `{TOP_N}`; e-value tier eşiği netleştirilir | Provenance + tutarlılık |

> Not: Bu tablo, daha önce ürettiğim `AUDIT_ISSUES.md` bulgularıyla bilinçli olarak örtüşür — refactor sırasında o düzeltmeler "ücretsiz" yapılır.

---

## 6. Veritabanı (Knowledge Base) Katmanı

Roadmap Bölüm 1 ile birebir hizalı. **SQLite ile başla (Ay 3) → PostgreSQL'e migrate et (Ay 4)**; aynı SQLAlchemy ORM her ikisinde çalışır.

### 6.1 Çekirdek tablolar (`db/schema/001_core.sql`)

```
organisms(organism_id PK, display_name, taxid, source, download_date, filter_criteria)
antibiotics(antibiotic_id PK, class_id, display_name)
models(model_id PK, organism_id FK, antibiotic_id FK, run_id FK,
       roc_auc, pr_auc, mcc, threshold, threshold_type, created_at)
kmers(kmer_id PK, sequence UNIQUE, k_length)
kmer_model_scores(id PK, kmer_id FK, model_id FK, gain_score,
                  seed_stability, cv_stability, combined_stability)
kmer_background_frequency(kmer_id FK, organism_id FK, freq_resistant, freq_susceptible)   # Roadmap 1.1
blast_annotations(id PK, kmer_id FK, source, gene_symbol, pident, evalue, sstart, send, tier)
kmer_antibiotic_overlap(ab1 FK, ab2 FK, overlap_count, expected, hypergeom_p)             # Roadmap 1.6
```

### 6.2 Provenance & validation (`db/schema/002_provenance.sql`) — Roadmap M10/M11

```
pipeline_runs(run_id PK, organism_id, antibiotic_id, git_commit_hash, seed,
              python_version, kmc_version, card_version, xgboost_version,
              kb_schema_version, started_at, finished_at, status)
validation_evidence(id PK, kmer_id FK, evidence_type ENUM(blast,resfinder,temporal_split,permutation),
                    evidence_source, evidence_score, pipeline_run_id FK)
```

### 6.3 İndeksler (`db/schema/003_indexes.sql`)
- `kmers.sequence` B-tree (exact match)
- `blast_annotations.gene_symbol` B-tree
- `kmer_model_scores.combined_stability` (sıralama)
- (Postgres) gerekirse `pg_trgm` 21-mer substring araması

### 6.4 ORM & populate
- `scripts/db/models.py` — SQLAlchemy modelleri (engine URL config'ten: `sqlite:///db/amrkb.sqlite` ↔ `postgresql://...`).
- `scripts/db/11_populate_database.py` — `--organism --antibiotic` alır; `07_explainability` CSV + `07b_stability` + `08` BLAST TSV → tabloları doldurur; her satıra `pipeline_run_id` bağlar.
- `composite_score = stability × log10(1/E-value) × (identity/100)` KB sıralaması için (Roadmap 1.4).

### 6.5 Migration
- Geliştirmede düz SQL dosyaları; ölçek büyüyünce **Alembic** (`db/migrations/`).
- `kb_schema_version` her şema değişikliğinde artar; Zenodo release ile eşlenir.

---

## 7. MLOps Yapısı

### 7.1 Run metadata & reproducibility (`scripts/lib/run_metadata.py`)
Her pipeline koşusu başında otomatik yakalanır ve `runs/{organism}/{antibiotic}/{run_id}/run_metadata.json`'a yazılır:
- `git_commit_hash` (`git rev-parse HEAD`), dirty flag
- `random_seed`, çözülmüş hiperparametreler
- Sürümler: Python, xgboost, scikit-learn, KMC (`kmc -version`), CARD (config'ten)
- Veri parmak izi: kullanılan chunk dosyalarının + `features.txt`'in hash'i / satır sayısı
- Başlangıç/bitiş zaman damgası, durum

Bu, KB `pipeline_runs` tablosuna birebir akar (provenance tek kaynak).

### 7.2 Model registry (`models/{organism}/{antibiotic}/manifest.json`)
Her eğitilen model için "model kartı":
```json
{
  "run_id": "ecoli__gentamicin__20260610T1432__a1b2c3d",
  "metrics": {"roc_auc": 0.97, "pr_auc": 0.95, "mcc": 0.82},
  "params": {...}, "data_split_hash": "...", "threshold": 0.5,
  "git_commit": "a1b2c3d", "created_at": "..."
}
```
"En iyi model" seçimi timestamp tahmini yerine manifest metriğinden yapılır (mevcut `_final_v2.json` / `_v2` sabit-adlandırma karmaşası çözülür).

### 7.3 Veri versiyonlama
- **DVC önerilir** (`dvc.yaml` + `.dvc` pointer'ları): büyük `data/` ve `models/` git dışında, hash ile sürümlenir.
- Alternatif (hafif): her organizma metadata + matris için `manifest.csv` (dosya → sha256, satır sayısı).
- Yayın için: versiyonlanmış matris/KB → **Zenodo DOI** (Roadmap M5/M10).

### 7.4 Orkestrasyon (`scripts/run_pipeline.py` + `Makefile`)
Tek komutla, parametreli, bağımlılık-sıralı çalıştırma:
```bash
python scripts/run_pipeline.py --organism ecoli --antibiotic gentamicin --from 01 --to 09
python scripts/run_pipeline.py --organism ecoli --all-antibiotics     # registry'den döner
make pipeline ORG=ecoli AB=gentamicin
```
- Adım bağımlılık kontrolü (önceki çıktı var mı?), checkpoint/resume, hata durumunda temiz çıkış.
- `--all-antibiotics` / `--all-organisms`: registry'deki `enabled` hedefler üzerinde döner → çok-antibiyotik batch (Roadmap N2 alt yapısı).

### 7.5 Ortam & CI
- `environment.yml` + sabitlenmiş `requirements.txt` (BioPython ekle, kullanılmayanları temizle).
- `tests/`: `lib/chunking.py`, `registry.py`, `config.resolve_path`, `validate_dataset_scientific` için birim testler (denetimdeki verification step ihtiyacı).
- (İleride) GitHub Actions: lint (ruff) + pytest + `dvc status`.
- (İleride, Roadmap N5) Dockerfile — tez için zorunlu değil.

---

## 8. Kademeli Geçiş Planı (faz faz, kırılmadan)

Her faz bağımsız; bir faz bitmeden diğeri başlamaz, ama her faz sonunda **pipeline çalışır durumda** kalır.

### Faz 0 — Hazırlık (yarım gün, risksiz)
- `SCALE_MLOPS_PLAN.md` (bu dosya) onayı.
- `config/registry/` ve `scripts/lib/`, `runs/`, `db/`, `tests/` boş iskeletleri açılır.
- Mevcut `AUDIT_ISSUES.md` düzeltmeleri ile bu refactor'ın kesişimi işaretlenir.

### Faz 1 — Tek kaynak + `lib/` (1–2 gün)
- `antibiotics.yaml` + `organisms.yaml` oluştur; mevcut `ANTIBIOTIC_CLASSES` taşı.
- `lib/registry.py`, `lib/config.py`, `lib/chunking.py`, `lib/io_utils.py`, `lib/run_metadata.py` yaz.
- `01`, `01b` registry'i kullanacak şekilde güncelle (kopya sözlük sil, çift-sayım fix).
- **Çıktı:** Davranış aynı; kod kopyası bitti. Pipeline E. coli'de aynen çalışır.

### Faz 2 — `{organism}` boyutu (1–2 gün)
- `config.yaml` path şablonlarına `{organism}` ekle; `project.organism = ecoli`.
- Mevcut veriyi yeni düzene taşı (tek seferlik, geri alınabilir taşıma scripti):
  - `data/raw/raw_genomes` → `data/raw/ecoli/genomes`
  - `data/external/metadata/genome_amr_matrix.csv` → `data/external/ecoli/metadata/amr_phenotypes.csv`
  - `data/interim/global_kmc_outputs` → `data/interim/ecoli/kmc_outputs`
  - `data/processed/{ab}` → `data/processed/ecoli/{ab}`; `models/`, `results/`, `logs/` aynı şekilde.
- Scriptlerde path kurma → `lib/config.resolve_path()`.
- **Çıktı:** Tek komutla ikinci organizma eklemeye hazır altyapı. E. coli sonuçları birebir reprodüke edilir.

### Faz 3 — Run metadata + model registry (1–2 gün)
- `04/05/06` → `run_id` üretir, `run_metadata.json` + `manifest.json` + `metrics.json` yazar.
- `06` Youden's J leakage fix (M1).
- **Çıktı:** Her koşu izlenebilir; KB provenance verisi hazır.

### Faz 4 — KB (SQLite) (Roadmap Ay 3, ~2 hafta)
- `db/schema/*.sql`, `db/models.py`, `11_populate_database.py`.
- En az `ecoli/gentamicin` + `ecoli/cefotaxime` doldur.
- `07b_feature_stability.py` (M4) + `10_cross_antibiotic_analysis.py` (S1).
- **Çıktı:** Çalışan SQLite KB + provenance + validation_evidence iskeleti.

### Faz 5 — Orkestrasyon + ikinci organizma denemesi (1 hafta)
- `run_pipeline.py` + `Makefile`; `--all-antibiotics`.
- Küçük bir K. pneumoniae alt kümesiyle **uçtan uca çoklu-organizma duman testi** (registry'de `enabled: true`).
- **Çıktı:** Yeni organizma = registry bloğu + veri; kod değişmeden tam pipeline.

### Faz 6 — Postgres + FAIR/API (Roadmap Ay 4)
- SQLite → PostgreSQL migration (aynı ORM).
- `12_export_fair_metadata.py` + FastAPI minimal endpoint (`/kmers`, `/overlap`, `/metadata`).
- Zenodo DOI + `kb_schema_version` (M10).

---

## 9. "Yeni Organizma Ekleme" Akışı (hedef son durum)

Refactor bittiğinde, yeni bir organizma eklemek şu kadar basit olmalı:

1. `config/registry/organisms.yaml`'e bir blok ekle (`enabled: true`).
2. `data/raw/{organism}/genomes/` içine `.fna` dosyalarını koy.
3. `data/external/{organism}/metadata/amr_phenotypes.csv` fenotip matrisini koy.
4. (Gerekirse) yeni antibiyotikleri `antibiotics.yaml` sınıflarına ekle.
5. `python scripts/run_pipeline.py --organism {organism} --all-antibiotics`.

**Hiçbir Python dosyası değiştirilmez.** KB, run metadata, modeller, sonuçlar otomatik olarak izole klasör/tablolarda toplanır.

---

## 10. Özet Kontrol Listesi (uygulama sırası)

- [ ] Faz 0: iskelet klasörler + plan onayı
- [ ] Faz 1: `antibiotics.yaml`, `organisms.yaml`, `lib/` modülleri, 01/01b registry'e geçiş + çift-sayım fix
- [ ] Faz 2: `config.yaml` `{organism}` şablonları + veri taşıma + `resolve_path` entegrasyonu
- [ ] Faz 3: `run_id` + run_metadata + manifest + Youden leakage fix
- [ ] Faz 4: `db/` şema + ORM + `11_populate` + `07b_stability` + `10_cross_antibiotic`
- [ ] Faz 5: `run_pipeline.py` + Makefile + ikinci organizma duman testi
- [ ] Faz 6: Postgres migration + FAIR metadata + FastAPI + Zenodo DOI

---

*Bu plan, mevcut `AMR_Tez_Roadmap_2026.md` (Ay 1–6, M/S/N maddeleri) ile uyumludur: buradaki Faz 1–3 Roadmap Ay 1–2'yi, Faz 4–5 Ay 3'ü, Faz 6 Ay 4'ü besler. Plandaki her yapısal düzeltme `AUDIT_ISSUES.md` bulgularıyla bilinçli kesişir; böylece ölçeklenme ve bug-fix tek refactor turunda tamamlanır.*
