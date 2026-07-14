# Modül 0 — Repository & Temel (Foundation)

> İnceleme merceği: **Python'da kal + sertleştir** (Nextflow rewrite yok) · **ESKAPEE/çok-organizma ölçeklenebilirliği** · **publication-ready**
> İncelenen dosyalar: `pyproject.toml`, `requirements.txt`, `environment.yml`, `environment-tools.yml`, `environment-checkm2.yml`, `.pre-commit-config.yaml`, `.github/workflows/ci.yml`, `Makefile`, `pytest.ini`, `README.md`, `CITATION.cff`, `.zenodo.json`, `CHANGELOG.md`, `LICENSE`, `.gitignore`, `*.def` (containers), dal/git yapısı.
> Tarih: 2026-07-13 · Repo HEAD `b413f84` (main, origin'in 7 commit önünde)

---

## 1. Genel Değerlendirme

Bu modül projenin **paketleme, bağımlılık, CI/kalite-kapıları, lisans, atıf/Zenodo ve git** zeminidir — yani "bu bir yayınlanabilir yazılım mı?" sorusunun altyapısı. Genel izlenim: **temel şaşırtıcı derecede olgun** (pyproject + ruff + mypy + pre-commit + çok-sürümlü CI + CITATION.cff + .zenodo.json + MIT/CC-BY ikili lisans + Makefile). Bir yüksek lisans projesinden beklenenin üstünde.

Ama zemin **bilimin gerisinde kalmış**: metadata hâlâ "*E. coli*, 4 antibiyotik, şema 0.4.0" derken gerçek durum "2 organizma, 21 model, şema 0.6.0". Üç ayrı bağımlılık dosyası birbiriyle çelişiyor ve **çekirdek compute container'ın (`amr.def`) reçetesi repoda yok** → şu an proje "kağıt üstünde reproducible ama pratikte değil". Bunlar mekanik/düşük-riskli düzeltmeler; bilim sağlam, kabuk güncellenmemiş.

---

## 2. Güçlü Yanlar

- **PEP 621 pyproject** + ruff/mypy/pytest tek dosyada konfigüre; `pip install -e ".[dev]"` ile QA zinciri.
- **Çok-sürümlü CI** (3.10/3.11/3.12), her push/PR'da ruff + pytest (unit+smoke); integration/slow bilinçli opt-in → CI hızlı ve config'i mutasyona uğratmıyor.
- **pre-commit** (trailing-ws, eof, yaml, large-files, merge-conflict, ruff --fix) — hijyen otomasyonu kurulu.
- **İkili lisans temiz:** kod MIT (LICENSE + classifier + cff tutarlı), veri CC-BY-4.0 (.zenodo). FAIR "Reusable" ayağı.
- **CITATION.cff + .zenodo.json** hazır (DOI slotu açık) — çoğu tez projesinde olmayan yayın altyapısı.
- **Makefile** insan-dostu hedeflerle (`setup/dev-install/lint/test/pipeline/data/features/train/biology`), `run_pipeline.py` orkestratörü mevcut.
- **Registry-güdümlü, organizma-agnostik tasarım iddiası** (README: "yeni organizma = registry kaydı, kod değil") → ESKAPEE'ye ölçeklenme için doğru zemin (M1/M3'te doğrulanacak).
- **Sızıntı yok:** PAT ne çalışma ağacında ne git geçmişinde; `.gitignore` kapsamlı.

---

## 3. Problemler (önem sırasıyla)

### Critical
- (yok) — foundation'da build-bozan / veri-bozan / sır-sızdıran sorun yok.

### High
- **H1 — [ÇÖZÜLDÜ 2026-07-13] Çekirdek container reçetesi `amr.def` repoda yoktu.** TRUBA'dan kurtarıldı ve repoya geri yüklendi (`amr.def`). İçeriği: `Bootstrap: docker / From: condaforge/miniforge3:latest`, `%post: mamba env create -f environment.yml`, `%runscript: exec python "$@"` (sif'i "çalıştırınca" Python REPL açılmasının nedeni budur). **Kritik gözlem:** container TAMAMEN `environment.yml`'den türüyor → M1'deki `environment.yml` sorunları (shap, pinsizlik, `:latest` base) doğrudan container reproducibility'sini belirliyor. Kalan iş: commit + base imajı digest'e pinle + lock. **Ayrıca `slurm/` (33 dosya) ve TRUBA-only kod farkları da incelendi → hiçbirinde repoda-eksik bir hotfix/bugfix YOK** (00a/02b/config farkları yalnızca origin'in daha yeni/yorumlu olması; TRUBA çalışma ağacı bayat kopya). Yani TRUBA sapması **kayıpsız** çözülebilir.
- **H2 — Public metadata bayat ve tutarsız (publication-blocking).** DÖRT farklı sürüm/kapsam dolaşıyor:
  - `pyproject.toml` → `version 0.1.0`, description "*for E. coli*"
  - `CITATION.cff` → `0.4.0`, başlık/abstract "*for Escherichia coli*"
  - `.zenodo.json` → `0.4.0`, "*for Escherichia coli*", "schema 0.4.0"
  - Gerçek KB → **şema 0.6.0, 2 organizma, 21 model**
  - `README` abstract "multi-organism" diyor ama devamında "*current dataset is E. coli across four antibiotics*"; `CHANGELOG [Unreleased]` "5470 genom / 3. antibiyotik" (E. coli-only dönemi).
- **H3 — `main`, `origin/main`'in 7 commit önünde (push edilmemiş).** 0.5.0/0.6.0 şema, figürler, external_concordance yalnızca yerelde. Laptop kaybı = kod kaybı. (Ortam-güvence adımı; onay+PAT ile push.)

### Medium
- **M1 — Üç bağımlılık kaynağı çelişiyor, lock yok.** `shap` → requirements/pyproject "ölü, kaldırıldı" der ama `environment.yml` hâlâ `shap>=0.44` içeriyor. API/UI (streamlit/fastapi/uvicorn/httpx) requirements'ta var, pyproject optional-deps'te ve environment.yml'de yok. `poppunk/bcalm/unitig-caller` sadece environment.yml'de. Docs `requirements.lock.txt`/`environment.lock.yml`'den bahsediyor ama **lock dosyaları commit'li değil** → tam reproducibility yok.
- **M2 — Vendor'lanmış ikili dosyalar** (`bin/bin/kmc`, `kmc_tools`, `kmc_dump`, `libkmc_core.a` ≈ 15 MB) git'te. Kendi pre-commit'inin `maxkb=2048` kuralını ihlal ediyor (biri 5.2 MB) ve platforma-bağlı binary reproducibility'yi bozuyor. `kmc` zaten `environment.yml`'de bioconda'dan geliyor → git'ten çıkarılmalı.
- **M3 — Kimlik/isim drifti: "k-mer" vs "unitig".** `pyproject name = amr-kmer-kb`, README abstract "alignment-free **k-mer**" merkezli; oysa yöntem artık **unitig** (k-mer sadece QC). CITATION/zenodo "unitig" diyor → dışa mesaj tutarsız. Ek: README `pipeline-Nextflow` rozeti orkestrasyonu abartıyor (gerçekte SLURM + 1 adım .nf).
- **M4 — setuptools paket-adı karışık:** `package-dir = {"amr_lib" = "scripts/lib"}` ama `packages.find ... include=["lib*"]` → kurulabilir paketin adı `amr_lib` mi `lib` mi belirsiz; scriptler `from lib ...` (sys.path) kullanıyor. pip-installable kimlik netleştirilmeli.

### Low
- **L1** — pytest config iki yerde (`pytest.ini` + `pyproject [tool.pytest]`) → tek kaynağa indir.
- **L2** — ruff kural seti çok dar (`E9,F63,F7,F82`); `F401` (kullanılmayan import), `I` (import sıralama) kademeli açılmalı.
- **L3** — mypy konfigüre ama CI'da koşmuyor; coverage raporu yok.
- **L4** — `amr-gpu.def` ölü artefakt (GPU reddedildi) → `archive/`'e taşı veya README'de "kayıt için" notu.
- **L5** — `environment.yml`'de `setuptools<81` iki kez listelenmiş.
- **L6** — `backup/` (6 GB) untracked **ve** `.gitignore`'da yok → `git add .` tehlikesi.
- **L7** — 3 eski dal (`feature/mlops-bio-summary`, `fix/amr-audit-remediation`, `fix/code-review-remediation`) — merge durumu doğrulanıp temizlenmeli.

---

## 4. Düzeltilmesi Gerekenler (madde madde)

1. **`amr.def`'i repoya ekle** (H1) — çekirdek container reçetesi + build komutu (`apptainer build --fakeroot amr.sif amr.def`) README/docs'a.
2. **Tek sürüm kaynağı belirle ve senkronla** (H2) — bir `VERSION`/pyproject `version` gerçeğe (`0.6.0` veya yeni bir yayın etiketi) çekilsin; `CITATION.cff`, `.zenodo.json`, README abstract, CHANGELOG hepsi "2 organizma / 21 model / şema 0.6.0"a güncellensin. "E. coli" → "multi-organism (E. coli, K. pneumoniae; ESKAPEE hedefi)".
3. **Bağımlılıkları tek kaynağa indir** (M1) — pyproject'i tek doğruluk kaynağı yap (core + `[project.optional-dependencies]` `api`, `ui`, `bio`, `dev`); requirements.txt'i pyproject'ten türet ya da kaldır; `shap` kararını netleştir (environment.yml'den çıkar ya da requirements/pyproject'e geri ekle). **Lock dosyalarını commit et** (`environment.lock.yml`, `requirements.lock.txt`).
4. **`bin/bin/*` ikilerini git'ten çıkar** (M2) — `git rm --cached`, bioconda/container'a bırak, `.gitignore`'a `bin/bin/` ekle (veya gerçekten gerekliyse git-lfs).
5. **`.gitignore`'a `backup/` ekle** (L6).
6. **push + commit bekleyen yerel iş** (H3) — onay + taze PAT.
7. **Dal temizliği** (L7) — `git branch --merged main` ile doğrula, merge edilmişleri sil; benzersiz iş varsa origin'e it.
8. **"k-mer" → "unitig" kimliğini netleştir** (M3) — en azından README abstract + pyproject description; paket adı (`amr-kmer-kb`) yayın öncesi kararı.
9. **CI'ya mypy + (opsiyonel) coverage ekle** (L3); ruff kurallarını kademeli genişlet (L2).
10. **pytest config'i tek yerde tut** (L1); `amr-gpu.def`'i arşivle (L4); `environment.yml` çift satırı temizle (L5).

---

## 5. Refactor Önerileri

- **Bağımlılık mimarisi:** `pyproject.toml` = tek kaynak. Ekstralar: `[project.optional-dependencies] api = [fastapi,uvicorn,httpx]`, `ui = [streamlit]`, `dev = [...]`. `environment.yml` yalnızca **conda-only** (kmc/blast/poppunk/bcalm/unitig-caller/nextflow) + `pip: -e .` ile pyproject'i çeksin. Böylece "üç dosya üç gerçek" sorunu biter.
- **Sürüm otomasyonu:** `version`'ı tek yerde tut (pyproject) ve CITATION/zenodo'ya release adımında script/`hatch`/`bump-my-version` ile yay; elle senkron drifti öldürür.
- **Container düzeni:** `containers/` klasörü aç, tüm `.def`'leri (amr, amr-tools, amr-checkm2) oraya taşı + her biri için tek satır build talimatı; `amr-gpu.def`'i `archive/`'e.
- **`bin/` kaldırımı:** KMC bioconda'dan geldiğinden `bin/bin/` tamamen silinebilir; scriptlerdeki `bin/bin/kmc` path referansları PATH/conda'ya çevrilmeli (M3 modülünde path-audit ile birlikte).
- **Paket kimliği:** `scripts/lib` → düzgün bir `amr_lib` paketi (ya da `src/amr_lib/`) haline getirilip scriptler `from amr_lib import ...` kullansın; `sys.path` hack'leri azalır (bu daha büyük bir refactor, M5-M8'e yayılır — burada sadece işaretliyorum).

---

## 6. Bilimsel Eksikler (makale açısından)

- **İsimlendirme = bilimsel doğruluk:** "k-mer tabanlı" ile "unitig tabanlı" ayrımı yayında net olmalı; şu an public kimlik ikisi arasında salınıyor. Makale başlığı/anahtar-kelimeleri unitig'i öne almalı (k-mer = QC).
- **Kapsam beyanı:** metadata "E. coli" derken bilim 2 organizma. Makale/depo "çok-organizmalı, cross-phylum (ESKAPEE hedefli)" olarak yeniden çerçevelenmeli — aksi halde reviewer "kod E. coli diyor, iddia çok-organizma" çelişkisini yakalar.
- **Reproducibility beyanı (software paper için):** çekirdek container reçetesi + lock dosyaları olmadan "fully reproducible" iddiası savunulamaz; bunlar makalenin "Availability" bölümünün ön-koşulu.

---

## 7. Literatür Gereksinimleri

Bu modül **mühendislik ağırlıklı → literatür taraması gerektirmiyor.** Yalnızca bir **karar** girdisi lazım (literatür değil, senin tercihin):

- **Hedef dergi(ler)in yazılım/veri erişilebilirlik politikası.** Öneri: 2 makaleden biri **yazılım/database makalesi** (ör. *Bioinformatics* Application Note, *GigaScience*, *Database (Oxford)*). O derginin "software/data availability + reproducibility checklist"i (container? lock? test coverage? DOI?) M0 düzeltmelerinin kapsamını belirler. Hedef dergiyi netleştirir misin? (Ben M11'de bu checklist'e göre publication-readiness'ı denetlerim.)

Başka literatür girdisi bu modül için gerekmez.

---

## Sonraki modüllere taşınan notlar
- **M1 (Config/registry):** README'nin "yeni organizma = registry kaydı, kod değil" iddiasını doğrula; ESKAPEE için registry şeması yeterli mi?
- **M3 (QC):** organizma-seviyesi scriptlerin (01/02/02b/02p/03) config'i doğrudan okuduğu (env-parametrik değil) HANDOFF'ta yazıyor → paralel ESKAPEE koşuları için darboğaz olabilir.
- **M3/M4 (path-audit):** `bin/bin/kmc` hard-coded path referanslarını scriptlerde ara.
- **M9 (repro):** TRUBA manuel yamaları (`config.yaml/02p/02b/03`) git'e girmemiş — `git diff --stat` çıktısı beklenecek.
