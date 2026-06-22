# TRUBA (ARF) Üzerinde ML_AMR_Prediction_v2 Projesini Kurma ve Çalıştırma Rehberi

> **Kime göre yazıldı:** HPC'ye yeni başlayan, temel Linux bilen bir Moleküler Biyoteknoloji yüksek lisans öğrencisi.
> **Kapsam:** TRUBA'nın **ARF** kümesine SSH ile giriş yaptıktan *sonraki* andan, projenin tam çalıştırılmasına ve sonuçların indirilmesine kadar.
> **Dayanak:** Yalnızca sağladığınız TRUBA dokümantasyonu (`docs.truba.gov.tr`) ve bu repodaki gerçek proje yapısı.
> **Önemli ön kabul:** OpenVPN bağlantınız aktif ve bir kullanıcı arayüzü sunucusuna (`arf-ui1..5`, `172.16.6.11–15`) `ssh kullanici_adiniz@172.16.6.11` ile bağlanabiliyorsunuz. Bu rehber o noktadan başlar.

> **Yayın teşekkürü (zorunlu):** TRUBA kaynaklarını kullanan tüm tez/makale/bildirilerde şu metni ekleyin:
> *"The numerical calculations reported in this paper were fully/partially performed at TUBITAK ULAKBIM, High Performance and Grid Computing Center (TRUBA resources)."*

---

## 0. Projeyi ve TRUBA'yı bir cümlede eşleştirmek

Bu proje (`ML_AMR_Prediction_v2`) *E. coli* genomlarından alignment-free **21-mer** özellikleri çıkarıp **out-of-core XGBoost** ile antibiyotik direnci tahmin eder, sonra önemli k-mer'leri CARD/NCBI BLAST + kararlılık + ayırt edicilik analizleriyle biyolojiye çevirir. Adımlar `scripts/00a … 11` numaralıdır ve `scripts/run_pipeline.py` ile sırayla çalıştırılır.

Kaynak ihtiyaçları (gerçek ölçümler, bu repodaki son çalıştırmadan):

| Aşama | Darboğaz | Tipik ihtiyaç |
|---|---|---|
| 02 KMC k-mer sayımı | Disk I/O (binlerce küçük dosya) | yüksek I/O → **node-yerel NVMe `/tmp`** |
| 03 matris kurulumu | RAM + I/O | 1552 genom × ~30M özellik, `features.txt` ≈ 0.75 GB |
| 04/05 HPO + eğitim | CPU (çok çekirdek) + RAM | çok çekirdek, GB'larca RAM |
| 07b kararlılık | CPU (5 tohum) | çok çekirdek |
| 00a indirme, 08 NCBI BLAST | **İnternet** | yalnızca UI node / yerel makinede |

Bu tablo rehberdeki tüm kararların nedenidir: **kod home'da, büyük veri scratch'te, geçici KMC dosyaları node `/tmp`'sinde, ortam Apptainer konteynerinde, internet gerektiren adımlar UI'da.**

---

# 1. TRUBA'ya giriş sonrası ilk yapılması gerekenler

## 1.1 Nerede olduğunuzu anlayın

Bağlandığınız sunucu bir **kullanıcı arayüzü (UI) sunucusudur** (`arf-ui1..5`). Doküman açıkça belirtir: UI sunucularında **yalnızca** derleme, dosya düzenleme, kısa testler ve iş gönderme yapılır. **Ağır hesap UI'da çalıştırılmaz** — sistem yöneticisi bu tür işleri sonlandırır, ısrar eden hesaplar askıya alınır. Gerçek hesaplar **kuyruğa (SLURM)** gönderilir.

İlk doğrulama komutları:

```bash
whoami                 # TRUBA kullanıcı adınız (SBATCH -A için gerekli)
hostname               # arf-ui1 ... gibi olmalı
pwd                    # /arf/home/<kullanici_adiniz> olmalı
id                     # grup bilginiz
```

## 1.2 Dizin yapısı: home ve scratch

ARF'te iki ana dosya sistemi vardır (doküman: *Dosya Sistemleri ve Depolama*):

| Dosya sistemi | Yol | Amaç | Yaşam süresi | inode limiti |
|---|---|---|---|---|
| **Ev dizini** | `/arf/home/$USER` | betikler, config, küçük girdi/çıktı, kurulumlar | kullanıcı kontrolünde (yedek YOK) | (genel) |
| **Scratch** | `/arf/scratch/$USER` | aktif hesaplama, **büyük veri**, paralel I/O | **maks. 30 gün, otomatik silinir** | **500K dosya** |

> **Kritik uyarılar (dokümandan):**
> - Hiçbir dosya sistemi **yedeklenmez**. Önemli verinin yedeği sizin sorumluluğunuzdadır (yerel makineye indirin).
> - Scratch **kalıcı değildir**; 30 günde temizlenir. Sonuçlarınızı bitince home'a veya yerel makineye taşıyın.
> - `/arf/scratch` için kullanıcı başına **500.000 dosya (inode)** sınırı vardır. Bu proje KMC ile **binlerce küçük dosya** üretebileceğinden inode yönetimi önemlidir (bkz. §5, §11).

Bu projeye uyarlanmış strateji:

```
/arf/home/$USER/ML_AMR_Prediction_v2     <- GİT REPOSU (kod, config, docs) — küçük, kalıcı
/arf/scratch/$USER/amr/                   <- BÜYÜK VERİ + ÇIKTILAR + KONTEYNER (geçici, hızlı)
        ├── containers/amr.sif            <- Apptainer ortamı (tek dosya)
        ├── data/                         <- genomlar, matrisler, metadata
        ├── results/ logs/ runs/ models/  <- pipeline çıktıları
```

## 1.3 Veri depolama stratejisi (neden böyle?)

- **Kod `/arf/home`'da:** Git reposu küçüktür (yalnızca kod/config/docs + ufak CARD homolog BLAST DB), kalıcı tutulmalı.
- **Büyük veri ve çıktılar `/arf/scratch`'te:** 77 GB'a varabilen genom/matris verisi ve yoğun I/O scratch'e aittir (paralel dosya sistemi, yüksek hız).
- **KMC geçici dosyaları node-yerel `/tmp` (NVMe)'de:** Doküman tüm sunucularda `/tmp`'nin NVMe olduğunu ve yüksek I/O için kullanılması gerektiğini söyler. KMC'nin geçici dosyaları için bu idealdir (bkz. §8).
- **Konteyner scratch'te:** `.sif` tek bir büyük dosyadır → inode dostudur; scratch'in hızından yararlanır.

## 1.4 Kota ve dosya (inode) yönetimi

Kota bilginize **UI'ya giriş yaparak** ulaşırsınız (doküman: *"Kullanım kota bilgilerinize arf-ui1 … giriş yaparak ulaşabilirsiniz"*). Genel komutlar:

```bash
# Disk kullanımınız (home ve scratch)
du -sh /arf/home/$USER
du -sh /arf/scratch/$USER

# Bir klasördeki DOSYA SAYISI (inode) — KMC sonrası kontrol edin
find /arf/scratch/$USER/amr -type f | wc -l

# En büyük dosyaları bul (dokümandaki öneri)
find /arf/home/$USER -type f -size +500M -exec ls -lh {} \;

# 30 günden eski scratch dosyalarını bul (temizlik için)
find /arf/scratch/$USER -type f -atime +30
```

İnode/kota baskısına karşı (doküman önerileri): merkezi yazılımları/konteyner kullanın, küçük dosyaları arşivleyin (`tar`), **shared dosya sistemine asla conda/pip ile binlerce küçük dosya kurmayın** (§4).

## 1.5 İlk işte mutlaka faydalı Linux komutları

```bash
ls -lah            # dosyaları izin/boyutla listele
cd / pwd           # gezinme
mkdir -p a/b/c     # iç içe klasör
cp -r / mv / rm -r # kopyala / taşı / sil (rm -r dikkatli!)
nano dosya.sh      # basit metin düzenleyici (vim alternatifi)
cat / less dosya   # dosya içeriğini gör (less: q ile çık)
tail -f dosya.log  # log dosyasını canlı izle
chmod +x betik.sh  # çalıştırılabilir yap
df -h . ; du -sh * # disk durumu / klasör boyutları
tar czf ad.tar.gz klasor/   # arşivle (inode tasarrufu)
tar xzf ad.tar.gz           # arşivi aç
module avail ; module list  # mevcut / yüklü yazılım modülleri
```

> **Güvenlik (dokümandan):** Diğer kullanıcıların verinize erişmesini engellemek için ev dizini izinlerini düzeltin: `chmod 700 $HOME`.

---

# 2. Proje için uygun klasör yapısının oluşturulması

Tek seferlik kurulum:

```bash
# 1) Scratch çalışma alanı (büyük veri + çıktılar + konteyner)
mkdir -p /arf/scratch/$USER/amr/{containers,data,results,logs,runs,models}

# 2) Kolaylık için kısa değişkenler (.bashrc'ye de eklenebilir)
export AMR_HOME=/arf/home/$USER/ML_AMR_Prediction_v2     # KOD
export AMR_WORK=/arf/scratch/$USER/amr                    # VERİ/ÇIKTI
echo "export AMR_HOME=$AMR_HOME" >> ~/.bashrc
echo "export AMR_WORK=$AMR_WORK" >> ~/.bashrc
source ~/.bashrc
```

Bu ayrım, repodaki `.gitignore` felsefesiyle birebir uyumludur: **kod sürümlenir, üretilen veri/çıktı sürümlenmez** — TRUBA'da da kod home'da (git), üretilen her şey scratch'te durur.

---

# 3. GitHub repository'sinin klonlanması

Kodu **home**'a klonlayın:

```bash
cd /arf/home/$USER
git clone https://github.com/demirbase/ML_AMR_Prediction_v2.git
cd ML_AMR_Prediction_v2
git log --oneline -3        # doğru sürümde olduğunuzu doğrulayın
ls scripts/                  # 00a_... 11_... + run_pipeline.py görmelisiniz
```

> Repo özelken `git clone` kullanıcı adı/parola (veya PAT) ister. Public ise doğrudan klonlanır. SSH anahtarı kullanmak isterseniz `ssh-keygen` ile UI'da anahtar üretip GitHub'a ekleyebilirsiniz, ama HTTPS klonlama en basitidir.

UI'da yalnızca **klonlama/düzenleme** yapın; ağır iş yok.

---

# 4. Python ortamının hazırlanması (Apptainer konteyneri ile — TRUBA'nın önerdiği yol)

## 4.1 Neden conda'yı doğrudan kurmuyoruz?

TRUBA dokümanı **kesin** bir kural koyar (*Python Kullanımı* ve *Dosya Sistemleri*):

> *"Depolama sistemine anaconda, miniconda, conda veya herhangi bir Python kütüphanesi kesinlikle yüklenmemelidir, pip ve türevleri kullanılmamalıdır. Küçük boyutlu yüz binlerce dosyadan oluştuğu için … dosya sistemlerinin performanslarını büyük ölçüde düşürür."*

Bu proje conda tabanlıdır (`environment.yml`: KMC, BLAST+, Nextflow + Python paketleri). Bu yüzden `conda env create` ile `/arf/home`'a kurmak **hem yasak hem de inode patlamasına yol açar.** Dokümanın önerdiği iki çözüm:

1. **Merkezi modüller** (`module load`) — hızlı ama bu projenin tüm araçlarını (KMC, BLAST+, Nextflow) içermez.
2. **Apptainer (Singularity) konteyneri** — tüm ortamı **tek `.sif` dosyasına** paketler: inode dostu, taşınabilir, **yeniden üretilebilir**. Bu projenin `environment.yml`'i zaten bunun için biçilmiş kaftandır.

**Bu rehber Apptainer konteynerini birincil yol olarak kullanır.** (Saf-Python ML adımları için merkezi `apps/truba-ai/cpu-2024.0` modülü alternatif olarak §4.6'da verilmiştir.)

## 4.2 Modül sistemi (temel komutlar)

```bash
module avail                       # tüm modüller
module list                        # yüklü modüller
module avail apptainer             # apptainer modülü var mı?
module load apptainer              # (gerekiyorsa) apptainer'ı yükle
apptainer --version                # doğrula
```

Doküman: *"Sistemlerimizde konteyner platformunu kullanmak için apptainer (önceki adıyla singularity) mevcuttur."* Eğer `apptainer` doğrudan PATH'te değilse `module load apptainer` ile yükleyin (`module avail apptainer` ile tam adı görün).

## 4.3 Konteyner tanım dosyası (repo'dan üretilen, yeniden üretilebilir)

> **Kural (dokümandan):** Konteyner **oluşturma** işlemi UI sunucularında YAPILMAZ. Önce `srun` ile interaktif bir hesap sunucusuna geçilir (bkz. §7.3), sonra build yapılır.

`environment.yml`'den bir tanım dosyası oluşturun. Home'da:

```bash
cd $AMR_HOME
cat > amr.def <<'EOF'
Bootstrap: docker
From: condaforge/miniforge3:latest

%files
    environment.yml /opt/environment.yml

%post
    # Tüm bilimsel yığını (KMC, BLAST+, Nextflow + Python paketleri) tek seferde kur
    mamba env create -f /opt/environment.yml -p /opt/amr-env
    mamba clean -a -y

%environment
    export PATH=/opt/amr-env/bin:$PATH
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

%runscript
    exec python "$@"

%labels
    Project ML_AMR_Prediction_v2
EOF
```

Bu tanım, `environment.yml` içindeki **KMC ≥3.2, BLAST+ ≥2.12, Nextflow ≥22.10** ve tüm Python bağımlılıklarını (xgboost, optuna, scikit-learn, biopython, certifi …) tek konteynere koyar.

## 4.4 Konteyneri inşa etmek (interaktif hesap sunucusunda)

```bash
# 1) Önce interaktif bir hesap sunucusu al (debug kuyruğu, kısa süre) — bkz. §7.3
srun -p debug -A $USER -N 1 -n 1 -c 4 --time=02:00:00 --pty /usr/bin/bash

# 2) Artık bir HESAP sunucusundasınız. Apptainer'ı hazırlayın
module load apptainer 2>/dev/null
# Build sırasında çok dosya açılır; geçici alanı NVMe /tmp'ye verin (inode/hız)
export APPTAINER_TMPDIR=/tmp/$USER-apptainer
export APPTAINER_CACHEDIR=/tmp/$USER-apptainer-cache
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

cd $AMR_HOME
# 3) Konteyneri SCRATCH'e inşa et (tek .sif dosyası, inode dostu)
#    İKİ KRİTİK GOTCHA (2026-06 deneyiminden):
#    (a) Build MUTLAKA hesap (compute) sunucusunda yapılmalı. UI'da son adım
#        "creating squashfs" CPU-yoğundur ve login-node CPU-time limitine takılıp
#        "CPU time limit exceeded (core dumped)" ile ölür. srun debug şart.
#    (b) Build'den ÖNCE `unset APPTAINER_BINDPATH`. .bashrc'deki BINDPATH=/arf,
#        build sırasında container'a /arf bind etmeye çalışır ama inşa edilen
#        fs'de /arf yoktur -> "destination /arf doesn't exist in container".
#        (BINDPATH yalnızca ÇALIŞTIRMA zamanında gerekli, build'de zararlı.)
unset APPTAINER_BINDPATH
apptainer build --fakeroot $AMR_WORK/containers/amr.sif amr.def

# GPU SÜRÜMÜ (opsiyonel, denendi ve bu veri/donanım için REDDEDİLDİ):
#   amr-gpu.def, conda CPU xgboost'u CUDA'lı PyPI wheel ile değiştirir
#   (USE_CUDA True). Ama 50.8M-feature ultra-geniş seyrek matris V100 16GB'a
#   sığmaz (GPU ELLPACK ~4 byte/nnz -> ~88 GB) ve ExtMem-GPU pratik değil.
#   Bu veri için CPU kullanın; GPU def repoda kayıt amaçlı durur.

# 4) İnteraktif oturumdan çık
exit
```

> İnternet: build sırasında `condaforge` ve bioconda paketleri indirilir; bu indirme interaktif hesap sunucusunda yapılır (dokümandaki Python konteyner örneği de aynı şekilde `conda install` ile paket çeker). İnşa uzun sürebilir; `--time` bol tutun.

## 4.5 Doğrulama testleri (konteyner çalışıyor mu?)

```bash
SIF=$AMR_WORK/containers/amr.sif

# Python ve kütüphaneler
apptainer exec $SIF python -c "import pandas,numpy,scipy,sklearn,xgboost,optuna,Bio,yaml,certifi; print('python deps OK')"

# Dış araçlar PATH'te mi?
apptainer exec $SIF kmc      | head -1
apptainer exec $SIF blastn  -version | head -1
apptainer exec $SIF nextflow -version | head -2

# Projenin kendi test paketi (kod home'da, konteyner ortamıyla)
cd $AMR_HOME
apptainer exec $SIF python -m pytest -q          # 56 birim/smoke testi
```

`pytest` yeşil geçiyorsa ortam hazırdır. (Ağır `integration` testi xgboost+KMC ister; konteynerde `apptainer exec $SIF python -m pytest -m integration -q` ile çalıştırılabilir, ama bunu küçük bir interaktif oturumda yapın.)

## 4.6 Alternatif: merkezi modül (yalnızca saf-Python adımlar için)

Sadece ML adımları (04/05/06/07b) için, KMC/BLAST gerekmediğinde, dokümandaki hazır yapay zeka modülü kullanılabilir (XGBoost + scikit-learn içerir):

```bash
module load apps/truba-ai/cpu-2024.0
python -c "import xgboost, sklearn; print('truba-ai OK')"
```

Ancak **00a/02/03/08/11** adımları KMC/BLAST/Nextflow/biopython istediğinden, **tutarlılık ve yeniden üretilebilirlik için tüm pipeline'ı tek konteynerle çalıştırmanız önerilir.** Karışık ortam kullanmayın.

---

# 5. Proje verilerinin yerleştirilmesi

Bu projenin beklediği yerleşim (repodaki `config/config.yaml` ve `lib/config.py` ile birebir):

```
data/raw/ecoli/genomes/*.fna                                  # girdi: genom assembly'leri
data/external/ecoli/metadata/amr_phenotypes.csv               # girdi: ikili fenotip matrisi
data/external/blast_db/card_nt/card.*                         # CARD homolog DB (repoda gelir)
data/interim/ecoli/kmc_outputs/                               # ara: KMC çıktıları (çok dosya!)
data/processed/ecoli/<antibiyotik>/matrix/                    # ara: .npz + features.txt + y
results/  logs/  runs/  models/                               # çıktı
```

## 5.1 Strateji: kod home'da, veri scratch'te (sembolik bağ ile)

`data/`, `results/`, `logs/`, `runs/`, `models/` klasörlerini scratch'e yönlendirin ki büyük veri home kotasını/ inode'unu yemesin:

```bash
cd $AMR_HOME
# Scratch'te gerçek klasörler
mkdir -p $AMR_WORK/{data,results,logs,runs,models}
# Repodaki yolları scratch'e sembolik bağla (repo bu adları kullanır)
for d in data results logs runs models; do
    rm -rf "$d" 2>/dev/null
    ln -s "$AMR_WORK/$d" "$d"
done
ls -l data results logs runs models      # -> hepsi scratch'e işaret etmeli
```

> **CARD homolog DB istisnası:** Bu DB repoda gelir (`data/external/blast_db/card_nt/`). `data`'yı scratch'e bağladığınız için, klonladıktan sonra bu DB'yi scratch'teki `data/external/blast_db/`'ye kopyalamanız gerekir:
> ```bash
> mkdir -p $AMR_WORK/data/external/blast_db
> git -C $AMR_HOME show HEAD:.gitignore >/dev/null 2>&1   # repo sağlam mı
> cp -r $AMR_HOME/.git/../data/external/blast_db/card_nt $AMR_WORK/data/external/blast_db/ 2>/dev/null || \
>   echo "Not: card_nt'yi repodaki konumdan scratch/data/external/blast_db/ altına kopyalayın."
> ```
> (Pratikte: klonladıktan hemen sonra, `data`'yı symlink yapmadan ÖNCE `cp -r data/external/blast_db $AMR_WORK/data/external/` yapıp sonra symlink kurmak en temizidir.)

## 5.2 Girdi verisinin scratch'e taşınması (internet/yerel makine üzerinden)

İki seçenek:

**(a) Veriyi yerel makinenizde hazırlayıp TRUBA'ya kopyalamak (önerilen).** İnternet gerektiren `00a`/`00` adımlarını yerel makinenizde (veya UI'da) çalıştırıp hazır `*.fna` + `amr_phenotypes.csv` üretin, sonra `rsync` ile scratch'e atın. Yerel makinenizden:

```bash
# Yerel makinede, proje kökünde — büyük veriyi UI üzerinden scratch'e yolla
rsync -avzP data/raw/ecoli/genomes/ \
  kullanici_adiniz@172.16.6.14:/arf/scratch/kullanici_adiniz/amr/data/raw/ecoli/genomes/

rsync -avzP data/external/ecoli/metadata/amr_phenotypes.csv \
  kullanici_adiniz@172.16.6.14:/arf/scratch/kullanici_adiniz/amr/data/external/ecoli/metadata/
```

> Doküman: dosya transferi için **arf-ui4 (172.16.6.14)** veya **arf-ui5 (172.16.6.15)** tercih edilir; `rsync -avzP` yalnızca değişen dosyaları aktarır, kesilirse kaldığı yerden devam eder. Paralel transferde eşzamanlı bağlantı sayısını **4–8** ile sınırlayın.

**(b) Veriyi UI sunucusunda indirmek.** UI internet erişimine sahiptir ve `00a` "hafif" bir indirme işidir (ağır hesap değil). UI'da konteynerle:

```bash
cd $AMR_HOME
apptainer exec $AMR_WORK/containers/amr.sif python scripts/00a_download_bvbrc.py --organism ecoli --backend api
apptainer exec $AMR_WORK/containers/amr.sif python scripts/00_prepare_metadata.py --organism ecoli
```

> `00a`/`00` ve `08`'in NCBI uzak BLAST kısmı **internet ister**. Hesaplama (batch) sunucularında dışa internet olmayabilir; bu yüzden bu adımları **UI'da** veya **yerel makinede** yapın, ağır hesabı (02–07b, 10) batch'e bırakın. CARD yerel BLAST internet istemez (DB yerel), batch'te çalışır.

## 5.3 Geçici dosyalar ve cache (inode kritik)

- **KMC geçici dosyaları:** node-yerel `/tmp` (NVMe). SLURM betiğinde `export TMPDIR=/tmp/$SLURM_JOB_ID` (bkz. §8). Doküman: yüksek I/O için `/tmp` NVMe kullanın; `export TMPDIR=/arf/scratch/...` da kabul edilebilir ama node-yerel `/tmp` daha hızlıdır.
- **KMC çıktıları (`data/interim/.../kmc_outputs/`):** binlerce küçük `.kmc_pre/.kmc_suf` → scratch'te durur, **iş bitince arşivleyin** (`tar`), inode'u boşaltın.
- **Apptainer cache:** `export APPTAINER_CACHEDIR=/tmp/...` (build sırasında), kalıcı tutmaya gerek yok.

## 5.4 Büyük veri dosyalarının yönetimi

- `features.txt` (~0.75 GB) ve `.npz` matris parçaları **scratch**'te kalmalı (repo `.gitignore`'u bunları zaten dışlar).
- İş bitince sonuç tablolarını (`results/.../05_explainability/*.csv`, `08_validation_metrics_*.json`, modeller, `runs/.../run_metadata.json`) **home'a veya yerel makineye** kopyalayın; ham matrisleri arşivleyip silin (scratch 30 günde temizlenir).

---

# 6. TRUBA (ARF) kaynaklarının verimli kullanılması

## 6.1 Kuyruklar (partition) ve donanım — dokümandan

| Kuyruk | Çekirdek/sunucu (min) | RAM/sunucu | Maks. süre | Not |
|---|---|---|---|---|
| **barbun** | min **20** | 384 GB | 3 gün | Xeon Gold 6248R — **YL için en uygun** |
| **hamsi** | min **28** | 192 GB | 3 gün | |
| **orfoz** | min **56** | 256 GB | 3 gün | Xeon Platinum 8480+ (112 çekirdek/sunucu) |
| barbun-cuda | min 20 + 1 GPU | 384 GB | 3 gün | 2× P100 (GPU işleri için) |
| akya-cuda | min 10 + 1 GPU | 384 GB | 3 gün | 4× V100, `/tmp` altında **1.4 TB NVMe** |
| **debug** | (çeşitli) | çeşitli | **4 saat** | **kısa testler** için |

## 6.2 Yüksek lisans öğrencisi kota sınırı (çok önemli)

Doküman (SSS): **Yüksek lisans öğrencileri için aynı anda en fazla 40 çekirdek** tanımlıdır (lisans 4, doktora/akademik 160). Standart disk kotası 1000 GB'a kadardır.

**Sonuç:**
- **orfoz** minimum 56 çekirdek ister → **40 çekirdek sınırıyla orfoz'a TEK iş bile gönderemezsiniz.** Kullanmayın.
- **barbun (min 20)** ve **hamsi (min 28)** sizin için uygundur. Bu projede **barbun'u, `-c 20` ile** kullanmanızı öneririm (hem alt sınırı karşılar hem 40 çekirdek bütçenizde rahat kalır, bol RAM verir).
- Kısa testler için **debug** (≤4 saat).

> Danışmanınızın bir **proje hesabı** varsa (ör. ARDEB/TBAG/BAP), `-A proje_hesabi` ile daha yüksek çekirdek ve öncelik kullanabilirsiniz. Kendi kullanıcı hesabınız için `-A $USER`.

## 6.3 CPU / RAM / disk / paralelleştirme — bu projeye uyarlama

Bu pipeline **çok-süreçli (multi-thread) tek-node** çalışır (MPI değil). Yani `-N 1` (tek node), `-n 1` (tek görev), `-c <çekirdek>` (görev başına çekirdek) kullanın ve uygulama içi paralelliği bu sayıyla eşitleyin:

`config/config.yaml` içinde TRUBA'ya göre ayarlayın:

```yaml
preprocessing:
  kmc_mem: 64          # GB — barbun 384GB/sunucu; -c 20 ile rahat (varsayılan 16'dan yükseltin)
  threads: 20          # KMC iş parçacığı = SBATCH -c ile aynı
  chunk_size: 200      # out-of-core RAM kontrolü; RAM darsa düşürün
xgboost_params:
  n_jobs: 20           # XGBoost iş parçacığı = SBATCH -c ile aynı
```

> **Kural:** `threads` ve `n_jobs` değerleri SLURM `-c` ile **aynı** olmalı. Daha fazlası kaynak israfı/yavaşlık, daha azı boş çekirdek demektir. Bir sunucudaki çekirdek sayısından fazla `-c` veremezsiniz.

## 6.4 Kaynak israfını önleme

- İşi göndermeden önce **küçük bir test** (az genom, debug kuyruğu) çalıştırın (§8.2).
- `--time`'ı gerçekçi tutun; bittiğinde işi salıverin. Aşırı `--time` kuyrukta daha uzun beklemenize yol açar.
- Gerekenden fazla çekirdek/RAM istemeyin; YL bütçeniz 40 çekirdektir, aynı anda çok iş açarsanız hepsi bu bütçeyi paylaşır.
- KMC ara dosyalarını iş bitince arşivleyip silin (inode).

---

# 7. İş gönderme sistemi (SLURM)

## 7.1 Kuyruk mantığı ve batch kavramı

UI'dan işinizi bir **betik (job script)** ile **SLURM**'a (`sbatch`) gönderirsiniz. SLURM işi uygun bir hesap sunucusunda, sıraya göre çalıştırır. İş bitene kadar UI'da beklemek zorunda değilsiniz; çıktılar dosyalara yazılır.

## 7.2 Temel SLURM komutları (dokümandan)

```bash
sbatch betik.slurm        # işi kuyruğa gönder
squeue -u $USER           # kendi işlerinizin durumu
scancel <JOBID>           # işi iptal et
sinfo                     # kuyrukların boş/dolu durumu
kuyruk                    # (TRUBA komutu) kuyruk doluluk özeti
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS,ReqMem,NCPUS   # bitmiş iş raporu
sstat -j <JOBID> --format=AveCPU,MaxRSS                              # çalışan iş anlık kaynak
scontrol show job <JOBID> # iş detayları
man sbatch                # parametre yardımı
```

## 7.3 İnteraktif iş (test, derleme, konteyner build için)

UI'da ağır iş yasak olduğundan, deneme/derleme için **interaktif** bir hesap sunucusu alın (doküman: `srun` ile debug kuyruğundan):

```bash
# 30 dk'lık interaktif kabuk, debug kuyruğu, 4 çekirdek
srun -p debug -A $USER -N 1 -n 1 -c 4 --time=00:30:00 --pty /usr/bin/bash
# (gerekirse belirli sunucu tipi: -C barbun gibi constraint eklenebilir)
```

Buradan çıkmak için `exit`. İnteraktif oturum konteyner build, hızlı `pytest`, küçük denemeler için idealdir.

## 7.4 SLURM betiğinin anatomisi (dokümandan)

Her betik 3 bölümdür: (1) `#SBATCH` tanımları, (2) modül/ortam yükleme, (3) komut. Temel `#SBATCH` anahtarları:

| Anahtar | Anlamı |
|---|---|
| `-p, --partition` | kuyruk (barbun/hamsi/debug) |
| `-A, --account` | TRUBA hesabı (`$USER` veya proje hesabı) |
| `-J, --job-name` | iş adı |
| `-N, --nodes` | sunucu sayısı (bu proje: 1) |
| `-n, --ntasks` | görev sayısı (bu proje: 1) |
| `-c, --cpus-per-task` | görev başına çekirdek (paralellik) |
| `--time` | `GG-SS:DD:SS` süre limiti |
| `--output` / `--error` | stdout / stderr log dosyaları |
| `--mem` | düğümden istenen RAM (ör. `120G`) |
| `--mail-user` / `--mail-type` | e-posta bildirimi (BEGIN,END,FAIL) |
| `--no-requeue` | sunucu arızasında baştan başlatma (doküman önerisi) |

---

# 8. Bu proje için gerekli job scriptlerinin oluşturulması

Aşağıdaki betikleri `$AMR_HOME/slurm/` altında oluşturun:

```bash
mkdir -p $AMR_HOME/slurm
```

## 8.1 Ortak başlık mantığı

Her betik: konteyneri kullanır (`apptainer exec`), `TMPDIR`'i node NVMe `/tmp`'ye verir, kodu home'dan çalıştırır (veri scratch'e symlink'li).

## 8.2 Küçük test işi (önce bunu çalıştırın — debug, ≤4 saat)

`slurm/00_test.slurm`:

```bash
#!/bin/bash
#SBATCH -p debug
#SBATCH -A kullanici_adiniz
#SBATCH -J amr-test
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=00:30:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail
module load apptainer 2>/dev/null || true
SIF=/arf/scratch/$USER/amr/containers/amr.sif
cd /arf/home/$USER/ML_AMR_Prediction_v2

echo "[$(date)] Node: $(hostname)  JobID: $SLURM_JOB_ID"
# Ortam + testler (gerçek pipeline'a girmeden önce sağlık kontrolü)
apptainer exec "$SIF" python -c "import xgboost,sklearn,Bio,certifi; print('env OK')"
apptainer exec "$SIF" kmc | head -1
apptainer exec "$SIF" python -m pytest -q
echo "[$(date)] TEST DONE"
```

Gönder ve izle:

```bash
cd $AMR_HOME/slurm
sbatch 00_test.slurm
squeue -u $USER
# bitince:
cat amr-test-*.out
```

## 8.3 Tam ölçekli analiz çekirdeği (barbun, çok çekirdek) — ana iş

Bu betik k-mer sayımı → matris → HPO → eğitim → değerlendirme → kararlılık → açıklanabilirlik → ayırt edicilik adımlarını (`01→10`, CARD-only) tek seferde çalıştırır. İnternet gerektiren `00a`/`00` ve NCBI BLAST'ı önce UI'da yaptığınızı varsayar (§5.2).

`slurm/run_core.slurm`:

```bash
#!/bin/bash
#SBATCH -p barbun
#SBATCH -A kullanici_adiniz
#SBATCH -J amr-core
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eposta_adresiniz
#SBATCH --no-requeue

set -euo pipefail
module load apptainer 2>/dev/null || true

SIF=/arf/scratch/$USER/amr/containers/amr.sif
cd /arf/home/$USER/ML_AMR_Prediction_v2

# KMC/geçici dosyalar için node-yerel NVMe /tmp (yüksek I/O + inode home/scratch'i yormaz)
export TMPDIR=/tmp/$SLURM_JOB_ID
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "[$(date)] Node $(hostname)  CPUs=$SLURM_CPUS_PER_TASK  Job=$SLURM_JOB_ID"

# Tüm adımları konteyner ortamında, orchestrator ile sırayla çalıştır.
# (run_pipeline.py varsayılan planı 01->10; ağ gerektiren 00a/00/08 hariçtir.)
apptainer exec --bind "$TMPDIR":"$TMPDIR" "$SIF" \
    python scripts/run_pipeline.py --organism ecoli --antibiotic ampicillin

echo "[$(date)] CORE PIPELINE DONE"
```

> `config.yaml`'da `threads`/`n_jobs`'ı **20** yaptığınızdan emin olun (= `-c 20`). RAM darsa `--mem`'i artırın veya `chunk_size`'ı düşürün.

## 8.4 Adımları ayrı işlere bölmek (daha güvenli, yeniden başlatılabilir)

Uzun pipeline'ı tek işte çalıştırmak yerine, ağır adımları ayrı betiklere bölebilirsiniz (`run_pipeline.py --only` ile). Bu, bir adım hata verirse baştan başlamamanızı sağlar:

```bash
# Sadece özellik çıkarımı (KMC + matris): I/O ağır
apptainer exec "$SIF" python scripts/run_pipeline.py --organism ecoli --only 02 02b 03
# Sadece eğitim hattı
apptainer exec "$SIF" python scripts/run_pipeline.py --organism ecoli --only 04 05 06
# Sadece biyoloji (07b->10; CARD yerel)
apptainer exec "$SIF" python scripts/run_pipeline.py --organism ecoli --only 07b 07 09 10
```

`run_pipeline.py --list` planı gösterir; `--from 02 --to 06` aralık çalıştırır. (Repo'daki `Makefile` hedefleri de aynı mantığı sarmalar: `features`, `train`, `biology`.)

## 8.5 İnternet gerektiren adımlar (UI'da, batch DEĞİL)

```bash
# UI'da (arf-ui), konteynerle — internet var
cd $AMR_HOME
apptainer exec $AMR_WORK/containers/amr.sif python scripts/00a_download_bvbrc.py --organism ecoli --backend api
apptainer exec $AMR_WORK/containers/amr.sif python scripts/00_prepare_metadata.py --organism ecoli
# 08 CARD yerel BLAST batch'te çalışır; NCBI uzak BLAST yalnız internetli ortamda (UI) çalışır.
```

## 8.6 Log dosyaları

`--output=%x-%j.out` ve `--error=%x-%j.err` → `işadı-jobid.out/.err`. Ayrıca `run_pipeline.py` kendi logunu `logs/run_pipeline_<organism>.log`'a yazar (scratch'e symlink'li). Pipeline adımlarının kendi ayrıntılı logları `logs/ecoli/` altında oluşur.

---

# 9. Çalışmaların izlenmesi

```bash
# Kuyruktaki/çalışan işleriniz
squeue -u $USER
# Belirli işin durumu (PD=bekliyor, R=çalışıyor, CG=bitiyor)
squeue -j <JOBID>
# Kuyruk doluluğu (boş node'a yönlenmek için)
sinfo ; kuyruk
# Canlı log izleme
tail -f $AMR_HOME/slurm/amr-core-<JOBID>.out
tail -f $AMR_WORK/logs/ecoli/*.log
# Çalışan işin anlık RAM/CPU kullanımı
sstat -j <JOBID> --format=JobID,AveCPU,MaxRSS,NTasks
# Bitmiş işin özeti (gerçekte ne kadar RAM/çekirdek kullandı?)
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem,NCPUS,ExitCode
```

Hata mesajları için **önce `.err` dosyasına** bakın, sonra pipeline'ın kendi loguna. `State=FAILED` veya `OUT_OF_MEMORY` görürseniz §11'e bakın.

> Doküman notu: `.err` içindeki *"task/cgroup: plugin not compiled with hwloc support, skipping affinity"* bir **uyarıdır**, işi etkilemez, yok sayın.

---

# 10. Sonuçların yönetimi (reproducibility + arşivleme)

## 10.1 Çıktıların düzenlenmesi

Pipeline çıktıları (scratch'te):

```
results/ecoli/ampicillin/04_evaluation/      # metrikler, ROC/PR, kalibrasyon, bootstrap CI
results/ecoli/ampicillin/05_explainability/  # top-k-mer CSV/FASTA, KB adayları,
                                             #   08_validation_metrics_*.json, 10_*_background_*.csv
models/ecoli/ampicillin/                     # eğitilmiş model + manifest.json
runs/ecoli/ampicillin/<run_id>/              # run_metadata.json (git hash, sürümler, seed)
```

## 10.2 Yeniden üretilebilirlik

Bu repo her çalıştırmada `runs/.../run_metadata.json` içine **git commit hash + kütüphane sürümleri + seed** yazar. TRUBA'da ek olarak:

```bash
# Hangi kod sürümüyle çalıştığınızı kaydedin
cd $AMR_HOME && git rev-parse HEAD > $AMR_WORK/runs/CODE_COMMIT.txt
# Konteyner ortamını kaydedin (tam paket listesi)
apptainer exec $AMR_WORK/containers/amr.sif conda list -p /opt/amr-env > $AMR_WORK/runs/ENV_FROZEN.txt
```

Aynı `amr.sif` + aynı git commit + aynı `config.yaml` = aynı sonuç. `.sif` dosyasını saklarsanız ortamı yıllar sonra bile birebir tekrar üretebilirsiniz (FAIR ilkesi).

## 10.3 Versiyonlama

- **Kod:** Git (home'daki repo). Değişiklik yaptıysanız commit edin; sonuç klasör adına commit hash'ini ekleyin.
- **Config:** Çalıştırmada kullandığınız `config.yaml`'ı sonuç klasörüne kopyalayın: `cp config/config.yaml $AMR_WORK/runs/config_used_<JOBID>.yaml`.

## 10.4 Arşivleme ve indirme (scratch 30 günde silinir!)

```bash
# Sonuçları + provenance'ı arşivle (inode dostu, tek dosya)
cd $AMR_WORK
tar czf results_ecoli_ampicillin_$(date +%Y%m%d).tar.gz results/ecoli/ampicillin runs models/ecoli/ampicillin
# Önemli arşivi home'a taşı (kalıcı) veya yerel makineye indir
mv results_ecoli_*.tar.gz /arf/home/$USER/      # ya da:
# (yerel makinede) rsync -avzP kullanici_adiniz@172.16.6.14:/arf/scratch/.../results_*.tar.gz ./
# Ham KMC ara dosyalarını ARŞİVLE ve SİL (inode boşalt)
tar czf kmc_outputs.tar.gz data/interim/ecoli/kmc_outputs && rm -rf data/interim/ecoli/kmc_outputs/*
```

---

# 11. Performans optimizasyonu ve olası hatalar

## 11.1 Disk I/O (en kritik darboğaz)

- KMC binlerce küçük dosya yazar → **node-yerel NVMe `/tmp`** kullanın (`export TMPDIR=/tmp/$SLURM_JOB_ID`). Doküman tüm node'larda `/tmp`'nin NVMe olduğunu belirtir.
- Büyük matrisleri **scratch**'te tutun (home'da değil). Home yüksek I/O için değildir (doküman).
- İş bitince ara dosyaları arşivleyip silin.

## 11.2 RAM yetersizliği (`OUT_OF_MEMORY` / job killed)

- 30M özellikli matriste XGBoost RAM yiyebilir. `--mem`'i artırın (barbun 384 GB/sunucu) veya `config.yaml: chunk_size`'ı düşürün (örn. 200→100; bir seferde RAM'de tek chunk tutulur).
- `sacct -j <JOBID> --format=MaxRSS,ReqMem` ile gerçek kullanımı görüp `--mem`'i ona göre ayarlayın.
- `max_bin=2` zaten ikili veriye göre RAM'i düşürür (repoda ayarlı).

## 11.3 CPU kullanımı

- `threads`/`n_jobs` = `-c` eşitliğini koruyun. CPU %100'e çıkmıyorsa adım I/O-bağlıdır (KMC) → `/tmp` çözümü.
- YL sınırı 40 çekirdek: tek işte `-c 20` mantıklı; iki işi paralel açarsanız `-c 20 + -c 20 = 40` sınırına denk gelir.

## 11.4 Kuyrukta uzun bekleme

- Doküman: ya çekirdek bütçenizi (40) doldurdunuz ya da kuyrukta yer yok. `sinfo`/`kuyruk` ile boş kuyruğa yönlenin (barbun yoğunsa hamsi deneyin, ama hamsi min 28 çekirdek ister).
- `--time`'ı küçültmek işin daha erken başlamasına yardımcı olur.

## 11.5 Yaygın hatalar (dokümandan + bu projeye özel)

| Belirti | Neden / Çözüm |
|---|---|
| `Ev dizinime dosya kopyalayamıyorum` | Disk kotası dolmuş (standart 1000 GB). Büyük veriyi scratch'e alın, home'u temizleyin. |
| Çok dosya hatası / yavaşlama | inode (500K) sınırına yaklaştınız → KMC çıktılarını `tar`'layıp silin; conda'yı shared FS'e kurmayın. |
| `AssocGrpCPUMinutesLimit` / `AssociationJobLimit` | Çekirdek/zaman bütçeniz dolu; iş bekler. Daha az çekirdek/iş veya proje hesabı (`-A`) kullanın. |
| İş `R`→ baştan başlıyor | Sunucu arızası → SLURM requeue. `#SBATCH --no-requeue` ekleyin. |
| `squeue: Socket timed out` | Sistem yoğun; birkaç dk sonra tekrar deneyin (geçici). |
| UI'da işim sonlanıyor | UI'da ağır iş yasak → `srun` interaktif veya `sbatch` kullanın. |
| `hwloc support, skipping affinity` (.err) | Zararsız uyarı, yok sayın. |
| `CERTIFICATE_VERIFY_FAILED` (00a) | Konteynerde `certifi` var; bu adımı internetli UI'da çalıştırın. |
| `KMC executable not found` | Konteyner kullanmıyorsunuz; tüm komutları `apptainer exec $SIF …` ile çağırın. |
| NCBI uzak BLAST takılıyor (08) | Batch node'unda internet yok → 08'i UI'da çalıştırın veya CARD-yerel ile yetinin. |
| Konteyner build UI'da reddediliyor | Build'i `srun` interaktif hesap sunucusunda yapın (doküman kuralı). |

---

# 12. Bu projeye özel öneriler (RSE / FAIR)

- **Repo yapısına saygı:** Kod `/arf/home`'da git ile; üretilen her şey (`data/processed`, `results`, `logs`, `runs`, `models`) `.gitignore`'da olduğu gibi scratch'te. TRUBA'da symlink ile bu ayrımı koruyun (§5.1).
- **Pipeline mimarisi:** `run_pipeline.py` orchestrator'ı SLURM içinde kullanın; ağ-gerektiren adımları (00a/00/08-NCBI) UI'da, hesap-yoğun adımları (02–07b,10) batch'te ayırın. Bu, TRUBA'nın "UI hafif / node ağır / node internetsiz" modeline tam oturur.
- **Veri boyutu:** Şu an 1788 genom × ~30M özellik. Tek node + out-of-core + `chunk_size` ile yönetiliyor. Veri büyürse: (a) `chunk_size`'ı sabit tutup `--mem`/`-c`'yi artırın; (b) farklı antibiyotikleri **ayrı SLURM işleri** olarak paralel gönderin (her biri `--antibiotic X`), 40 çekirdek bütçesini bölüştürün.
- **Gelecekteki büyük veri:** Daha çok genomda KMC I/O darboğazı büyür → mutlaka node-yerel `/tmp` NVMe; gerekirse yüksek-I/O `akya-cuda` (`/tmp`'de 1.4 TB NVMe) düğümü değerlendirilebilir. GPU gerekirse `apps/truba-ai/gpu-2024.0` + GPU kuyruğu (ama bu proje şu an CPU-XGBoost).
- **Yeniden üretilebilirlik (FAIR):** `amr.sif` + git commit + `config.yaml` üçlüsünü her sonuç arşiviyle birlikte saklayın. Konteyner = ortamın dondurulmuş, taşınabilir kopyası. Yayında TRUBA teşekkür metnini ekleyin.
- **Sürdürülebilirlik:** Yeni organizma/antibiyotik eklemek için kod değiştirmeyin — `config/registry/*.yaml`'a blok ekleyip veriyi `data/raw/{organism}/`'a koyun (repo CONTRIBUTING.md). TRUBA'da yalnızca yeni veriyi scratch'e atıp aynı SLURM betiğini `--organism/--antibiotic` ile yeniden gönderin.

---

# 13. Sorun giderme — hızlı başvuru

```bash
# "İşim neden çalışmıyor / bekliyor?"
squeue -u $USER ; sinfo ; kuyruk
scontrol show job <JOBID> | grep -i reason

# "İşim çöktü, neden?"
cat slurm/amr-core-<JOBID>.err          # önce hata dosyası
tail -50 logs/ecoli/*.log                # pipeline logu
sacct -j <JOBID> --format=State,MaxRSS,ReqMem,Elapsed,ExitCode

# "Kotam/inode doldu mu?"
du -sh /arf/home/$USER /arf/scratch/$USER
find /arf/scratch/$USER/amr -type f | wc -l     # 500K'ya yaklaşıyorsa arşivle

# "Ortam bozuk mu?"
apptainer exec $AMR_WORK/containers/amr.sif python -m pytest -q

# "Diğerleri verime erişiyor mu?"
chmod 700 $HOME

# "Şifre / hesap"
passwd                    # UI'da şifre değiştir
# hesap askıya alındıysa: trubadestek@tubitak.gov.tr (kullanıcı adı + proje belirtin)
```

**Destek:** `trubadestek@tubitak.gov.tr` (e-postada TRUBA kullanıcı adınızı ve varsa proje hesabınızı belirtin).

---

## Özet akış (tek bakışta)

```
[UI] OpenVPN+SSH  →  git clone (home)  →  amr.def hazırla
[srun debug]      →  apptainer build amr.sif (scratch'e)
[UI/yerel]        →  00a + 00 (internet) → veriyi scratch'e rsync
[home]            →  data/results/... scratch'e symlink; config.yaml (-c=threads=n_jobs=20, kmc_mem=64)
[sbatch debug]    →  00_test.slurm (sağlık kontrolü)
[sbatch barbun]   →  run_core.slurm  (01→10, TMPDIR=/tmp)   ← ANA İŞ
[UI/internet]     →  08 NCBI BLAST (gerekirse)
[home/yerel]      →  sonuçları tar + indir; KMC arası arşivle+sil
```

Bu rehberdeki tüm yollar, kuyruk adları, sınırlar ve kurallar sağladığınız TRUBA dokümantasyonundan; tüm komut uyarlamaları bu repodaki gerçek `scripts/`, `config/config.yaml`, `environment.yml` ve `run_pipeline.py` yapısından alınmıştır.

---

# Ek: Gerçek Dağıtımdan Düzeltmeler (ARF, 2026-06 — `edemirbas`)

> Bu bölüm, rehberi ARF üzerinde **gerçekten çalıştırırken** öğrenilen ve yukarıdaki genel anlatımı **geçersiz kılan/güncelleyen** noktaları içerir. Çelişki olursa **bu bölüm geçerlidir.**

1. **Apptainer modül değil:** `/usr/bin/apptainer` (v1.3.6) doğrudan PATH'te. `module load apptainer` **çalışmaz/gereksiz**. Direkt `apptainer` kullan.
2. **Konteyner inşası:** interaktif debug node'da, `apptainer build --fakeroot $AMR_WORK/containers/amr.sif amr.def`. Öncesinde `export APPTAINER_TMPDIR=/tmp/$USER-ap APPTAINER_CACHEDIR=/tmp/$USER-cache` (inode/hız). İnternet build sırasında çalışır.
3. **`APPTAINER_BINDPATH=/arf` ZORUNLU:** Apptainer `/arf/scratch`'i otomatik bağlamaz; `data/` symlink'leri scratch'e gittiği için bind olmadan konteyner içinde `FileNotFoundError`/`mkdir` hatası alırsın. `~/.bashrc`'ye ekle **ve her SLURM betiğinde `export APPTAINER_BINDPATH=/arf`** yaz.
4. **İşler `/arf/scratch`'ten gönderilmeli:** `sbatch`/`srun` öncesi `cd $AMR_WORK`. Aksi halde `srun: error: Lutfen islerinizi /arf/scratch/ dizini altinda calistiriniz`. (Betik içinde `cd $AMR_HOME` yapmak serbest; kural **gönderim dizini** içindir.)
5. **`ftp.bv-brc.org` ENGELLİ:** TRUBA'dan FTP'ye bağlantı timeout (HTTP=000). Genom FASTA'ları **BV-BRC Data API**'sinden (`www.bv-brc.org/api/genome_sequence`, `Accept: application/dna+fasta`) indirilir — bu artık `00a`'nın **repo varsayılanı**. API host'u erişilebilir. (Bazı genomların API'de dizisi yok → "empty/non-FASTA", elenir; 5865 adaydan **5470 indi**.)
6. **`barbun` min 20 çekirdek/node** (hamsi 28, orfoz 56). **YL limiti 40 çekirdek.** `-c 8` reddedilir. KMC/QC/matris için `-c 20`, paralel ML (04/05) için `-c 40`.
7. **Düşük CPU verimi cezalıdır:** TRUBA, çekirdekleri boşa kullanan işi uyarır ve **otomatik iptal edip çekirdek hakkını düşürebilir** (`Eff:%2` uyarısı görüldü). KMC sıralı olduğu için boşa çekirdek yakar → **`scripts/02p_kmer_parallel.py`** ile genom-paralel KMC kullan (eşzamanlılık = `preprocessing.threads` = `-c`). 5470 genom **~2.5 dk**'da bitti.
8. **Kotalar (banner):** `/arf/home` 100 GB / 100K inode, `/arf/scratch` 1 TB / 200K inode, **yedek yok**, `/arf` NVMe Lustre (hızlı; KMC için ayrı `/tmp` şart değil).
9. **conda/pip shared FS'e KURULMAZ** → her şey Apptainer konteynerinde. Merkezi `apps/truba-ai` modülü KMC/BLAST/Nextflow içermez; konteyner kullan.
10. **Gerçek SLURM betikleri** `$AMR_HOME/slurm/`'de: `00_test.slurm` (debug sağlık), `run_features.slurm` (`02p`→`02b`→`03`, `-c 20 --mem 120G`), `run_matrix.slurm` (sadece `03`, `-c 20 --mem 240G`). ML için sıradaki: `04`→`05`→`06` `-c 40` ile; sonra `07b`→`07`→`09`→`10`.
11. **config.yaml TRUBA ayarı (repoya commit ETME):** `kmc_mem:128 threads:20 n_jobs:40 chunk_size:200 n_trials:30 target:ampicillin`.
12. **`git pull` YAPMA (TRUBA'da):** çalışma kopyası elle yamalı + `02p` + TRUBA'ya özel config içerir; pull çakışır. Kod düzeltmeleri zaten `main`'de.
13. **Paralel 02b ve 03 (verim):** 02b'nin 5470 KMC veritabanından spektrum çıkarması ve 03'ün genom-başına `kmc_tools dump`'ı serialdi → düşük verim. İkisi de thread-pool ile paralelleştirildi (`workers = preprocessing.threads`; `main` commit'leri `1cd0119`, `3ddc476`). **Uyarı:** 03'ün genom-başına *parse* kısmı (5M k-mer'i ~8 GB'lık Python sözlüğüyle eşleştirme) **GIL'e bağlı** — thread'ler bu kısmı tam hızlandıramaz, 03 yine ~%5 verimde görünebilir ve TRUBA uyarısı gelebilir. **03 resume-safe** (mevcut `*.npz` + `features.txt` atlanır), o yüzden iş öldürülse bile `sbatch run_features.slurm` ile kaldığı chunk'tan devam eder — panik yok. **Gerçek 03 çözümü (gelecek):** k-mer'leri 2-bit integer kodlayıp numpy `searchsorted` (Python dict yerine) veya çoklu-süreç.
14. **Gerçek ölçek (ampicillin, 5470 genom):** matris = **22 chunk × ~50.8M k-mer**, `features.txt` ≈ **1.27 GB**, ~%90 seyrek. 02 KMC paralel ~2.5 dk; 02b+03 birkaç saat (03 tek-thread parse yüzünden). Tüm 5470 KMC db'si `data/interim/ecoli/kmc_outputs/`'ta.
15. **ML işi (`run_ml.slurm`):** 04 HPO → 05 eğitim → 06, **`-c 40 --mem 300G --time 3-00:00:00`**, `n_jobs=40`. XGBoost ağaç kurarken çekirdekleri kullanır → 03'ten yüksek verim. 04, 30 Optuna trial'ını 22 chunk üzerinde out-of-core çalıştırır (en uzun faz). 04 resume-safe **değil** (yeniden başlar). Çıktı: `config/experiments/ecoli/config_ampicillin.yaml` + model + 06 metrikleri.
16. **Biyoloji işi (ML'den sonra):** `07b → 07 → 09 → 10`, `-c 20`. `08` CARD-yerel BLAST hesap node'unda çalışır; `08` NCBI-remote + `09` Entrez **internet** ister → bunları **UI'da** çalıştır (hesap node'larında dış internet olmayabilir) ya da atla.
17. **Sonuçları indir (scratch 30 günde silinir):** bitince `tar czf` ile `results/ models/ runs/` → `/arf/home` veya yerel makineye `rsync`.
