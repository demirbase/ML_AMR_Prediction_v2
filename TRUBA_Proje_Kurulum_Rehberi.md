# TRUBA (HPC) Proje Kurulum ve Çalıştırma Rehberi: AMR Prediction Pipeline

Bu rehber, Moleküler Biyoteknoloji yüksek lisans öğrencileri için sıfırdan TRUBA (Yüksek Başarım Merkezi) üzerinde çalışmaya başlamak ve mevcut **Antimikrobiyal Direnç (AMR) Tahmin** projesini (Makine Öğrenmesi pipeline'ı) güvenilir ve tekrarlanabilir bir şekilde yürütmek amacıyla hazırlanmıştır.

## İçindekiler
1. TRUBA'ya Giriş ve Temel Kavramlar
2. Dosya Sistemini Anlamak (Home vs Scratch)
3. Projeyi TRUBA'ya Aktarmak
4. TRUBA Modül Sistemi ve Miniconda
5. Proje Bağımlılıklarının (Conda) Kurulumu
6. Slurm Kuyruk Sistemine Giriş
7. Pipeline İçin Slurm İş Betiği (Job Script) Hazırlama
8. İşleri Kuyruğa Gönderme (Job Submission)
9. İş Durumunu Takip Etme ve Yönetme
10. Çıktıların (Logların) İncelenmesi
11. GPU Kullanımı ve Kısıtlamalar
12. Performans Optimizasyonu ve İyi Uygulamalar
13. Sık Karşılaşılan Sorunlar ve Çözümleri

---

## 1. TRUBA'ya Giriş ve Temel Kavramlar

TRUBA sistemine erişim için SSH (Secure Shell) kullanılır. Doğrudan hesaplama düğümlerine (node) bağlanılamaz; sadece **Kullanıcı Arayüzü (User Interface - UI/Login)** sunucularına giriş yapılabilir.

Terminali (veya MobaXterm, PuTTY gibi bir aracı) açın ve aşağıdaki komutla bağlanın:
```bash
ssh kullanici_adiniz@172.16.7.1 # VPN ile bağlanıyorsanız veya TRUBA'nın güncel giriş adresi (genellikle levrek1, barbun vs. UI sunucuları)
```
> **Önemli:** Giriş sunucularında (login nodes) **kesinlikle** ağır hesaplamalar, kod derleme veya model eğitimi yapmayın. Bu sunucular sadece dosya transferi, iş betiği (script) hazırlama ve işleri kuyruğa göndermek (submit) içindir.

---

## 2. Dosya Sistemini Anlamak (Home vs Scratch)

TRUBA'da kullanıcılar için iki ana disk alanı bulunur. Bu projede büyük veri setleri ve modeller kullanılacağı için bu ayrım hayatidir:

- **`/truba/home/kullanici_adiniz` (Ev Dizini):** 
  - Kapasitesi sınırlıdır ve yedeklenir.
  - Sadece kaynak kodları, metin dosyaları, küçük konfigürasyon dosyalarını (örn. proje repo'nuzu) tutmak içindir.
  - **DİKKAT:** Büyük paket kurulumları (özellikle Conda/Pip environment'ları) veya yoğun I/O gerektiren işlemler burada **yapılmamalıdır**. Sistem performansını düşürür.

- **`/truba/scratch/kullanici_adiniz` (Çalışma Dizini):**
  - Yüksek kapasitelidir ve hızlı disk okuma/yazma (I/O) sağlar.
  - Makine öğrenmesi eğitim süreçleri (training), geçici dosyalar, büyük FASTA dosyaları, KMC analiz sonuçları burada olmalıdır.
  - **DİKKAT:** Yedeklenmez. İşiniz bittiğinde önemli sonuçları kendi bilgisayarınıza veya `home` dizinine almalısınız.

---

## 3. Projeyi TRUBA'ya Aktarmak

Projeyi klonlamak için `home` dizininde bir klasör oluşturun:
```bash
cd /truba/home/$USER
git clone <projenizin-github-linki> ML_AMR_Prediction
cd ML_AMR_Prediction
```
> *Eğer GitHub'da proje private (gizli) ise Personal Access Token kullanarak veya doğrudan bilgisayarınızdan `scp` komutuyla dosyaları TRUBA'ya kopyalayarak aktarabilirsiniz.*

---

## 4. TRUBA Modül Sistemi ve Miniconda

TRUBA'da yazılımlar "Modül" sistemiyle yönetilir. Python, R, GCC vb. programlar önceden kuruludur ancak sizin bunları aktif etmeniz gerekir.

Sistemde yüklü modülleri görmek için:
```bash
module avail
```

Projemiz `Conda` kullandığı için TRUBA üzerindeki merkezi miniconda modülünü yükleyeceğiz:
```bash
module load miniconda3
```
> *Not: Bu komutu her TRUBA'ya girdiğinizde ve her Slurm betiğinizin içinde çalıştırmanız gerekecektir.*

---

## 5. Proje Bağımlılıklarının (Conda) Kurulumu

Projenin bağımlılıkları `environment.yml` dosyasında (Python paketleri, KMC, BLAST, Nextflow) tanımlanmıştır. 
**Kritik Kural:** Conda ortamını `home` dizininde oluşturursanız, çok sayıda küçük dosya (metadata) yüzünden TRUBA kota/performans limitlerine takılabilirsiniz. Bu yüzden Conda ortamının sadece projenize gereken paketlerini kurmaya dikkat edin. Gerekirse conda prefix kullanarak scratch dizinini hedef gösterin. 

TRUBA önerilerine uygun olarak, ortam kurulumunu aşağıdaki gibi login node'da (veya interaktif oturumda) yapabilirsiniz:

```bash
# Modülü yükle
module load miniconda3

# Proje klasörüne gir
cd /truba/home/$USER/ML_AMR_Prediction

# Conda ortamını ortam dosyasına dayanarak kurulum yapın
# (Bu işlem internet yoğunluğuna bağlı olarak vakit alabilir)
conda env create -f environment.yml
```
*(Conda ortamınız `amr-prediction` adıyla kurulacaktır).*

---

## 6. Slurm Kuyruk Sistemine Giriş

Hesaplama işlerimizi doğrudan çalıştırmak yerine, işlerimizi bir dosya (betik) halinde yazıp Slurm'a (iş yönetim sistemi) iletmeliyiz. Slurm, boş bir sunucu (node) bulduğunda kodunuzu çalıştırır.

Temel komutlar:
- `sbatch <betik.sh>`: İşi kuyruğa gönderir.
- `squeue -u $USER`: Kuyruktaki kendi işlerinizi gösterir.
- `scancel <jobID>`: Çalışan veya bekleyen işi iptal eder.

**Partition (Kuyruk) Tipleri:**
- CPU işleri için genelde `sardalya`, `barbun`, `hamsi` gibi kuyruklar kullanılır.
- GPU (XGBoost GPU ivmelendirmesi gerekiyorsa) için `akya-cuda` veya `kolyoz-cuda` partition'ları kullanılmalıdır.

---

## 7. Pipeline İçin Slurm İş Betiği (Job Script) Hazırlama

Projenizde süreci çalıştırmak için bir `Makefile` bulunuyor. Makine öğrenmesi eğitim aşaması uzun süreceği için bunu Slurm üzerinden çalıştıracağız.

`ML_AMR_Prediction` dizini içine `run_pipeline.slurm` adında bir dosya oluşturun:

```bash
#!/bin/bash
#SBATCH --job-name=amr_pipeline        # İşin adı
#SBATCH --partition=barbun              # ARF CPU kuyruğu (barbun, hamsi veya orfoz)
#SBATCH --nodes=1                      # Kaç sunucu kullanılacak
#SBATCH --ntasks=1                     # Görev sayısı
#SBATCH --cpus-per-task=20             # İşlemci çekirdek sayısı (Barbun için min 20 zorunludur)
#SBATCH --time=48:00:00                # Maksimum çalışma süresi (Sa-Dk-Sn)
#SBATCH --output=logs/slurm-%j.out     # Standart çıktı dosyası (%j JobID ile değişir)
#SBATCH --error=logs/slurm-%j.err      # Hata çıktı dosyası
#SBATCH --mail-type=ALL                # İşe başlama ve bitişte mail at
#SBATCH --mail-user=mail@adresiniz.edu # Mail adresiniz

echo "=========================================================="
echo "İşlem Başlıyor..."
date
echo "Çalışan Node: $SLURM_NODELIST"

# Modülleri yükle
module purge
module load miniconda3

# Conda ortamını aktif et
source activate amr-prediction

# Proje dizinine git (Kendi kullanıcı adınızı yazın)
cd /arf/home/$USER/ML_AMR_Prediction

# Pipeline'ı otomatik antibiyotik seçimiyle çalıştır
make pipeline ORG=ecoli AB=auto

echo "=========================================================="
echo "İşlem Bitti."
date
```

---

## 8. İşleri Kuyruğa Gönderme (Job Submission)

Betik dosyasını hazırladıktan sonra işi sisteme verin:
```bash
sbatch run_pipeline.slurm
```
Ekranda size bir iş numarası (Job ID) dönecektir (Örn: `Submitted batch job 1234567`). Bu ID ile işinizi takip edeceksiniz.

---

## 9. İş Durumunu Takip Etme ve Yönetme

İşin durumunu kontrol etmek için:
```bash
squeue -u $USER
```
Eğer `PD` (Pending) görüyorsanız iş sıradadır. `R` (Running) görüyorsanız hesaplama düğümünde kodunuz çalışıyor demektir.

> **İpucu:** Ağ yoğunluğu kaynaklı `squeue` hataları alırsanız paniğe kapılmayın, kısa bir süre bekleyip tekrar deneyin. TRUBA dokümantasyonuna göre bu olağan bir durumdur.

İşinizi iptal etmeniz gerekirse:
```bash
scancel 1234567
```

---

## 10. Çıktıların (Logların) İncelenmesi

Job script'imizde çıktı dosyalarını `logs/` klasörüne kaydetmesini söylemiştik. İş çalışırken eşzamanlı olarak çıktıları okuyabilirsiniz:
```bash
tail -f logs/slurm-1234567.out
```
Eğer işlemlerde hata çıkarsa `logs/slurm-1234567.err` dosyasını inceleyin. 

> **Not:** TRUBA'da `err` dosyalarında sıklıkla `task/cgroup` uyarıları görebilirsiniz (kaynak yönetimi ile ilgilidir). Kodunuzla ilgili gerçek bir hata satırı (Python Traceback, Out of Memory vb.) görmüyorsanız bu cgroup uyarılarını genelde göz ardı edebilirsiniz.

---

## 11. GPU Kullanımı ve Kısıtlamalar

XGBoost modelinizin eğitiminde GPU kullanmak isterseniz (veri büyükse bu süreyi çok kısaltır), Slurm betiğinizi şu şekilde değiştirmelisiniz:

```bash
#SBATCH --partition=akya-cuda          # veya kolyoz-cuda
#SBATCH --gres=gpu:1                   # 1 adet GPU talep et
```
Ayrıca projenizin XGBoost konfigürasyonunda `tree_method='hist', device='cuda'` ayarlarının kullanıldığından emin olun.
TRUBA'da CUDA modülünü de betiğe eklemelisiniz:
```bash
module load cuda/11.8 # (Veya sistemdeki uyumlu sürüm)
```

---

## 12. Performans Optimizasyonu ve İyi Uygulamalar

- **I/O Yoğunluğu (KMC, Nextflow):** K-mer sayma (KMC) ve Nextflow işlemleri disk üzerinde çok sayıda geçici dosya üretir. Çok yüksek I/O operasyonları gerektiren bir adımınız olursa, bu işlem adımını kod içerisinden geçici olarak `/truba/scratch/$USER/` dizininde gerçekleşecek şekilde yönlendirmeniz tavsiye edilir.
- **MPI / UCX Optimizasyonu:** Eğer birden fazla node kullanacak şekilde dağıtık (distributed) eğitim yapmıyorsanız (ki standart XGBoost ve Scikit-learn genelde tek node kullanır), MPI ayarlarına şimdilik girmeyin. `nodes=1` her zaman daha hızlı çalışmaya başlar.
- **Konteyner Kullanımı (Apptainer):** TRUBA belgelerinde önerildiği üzere Conda ortamı sorun çıkartırsa (veya farklı paket bağımlılığı oluşursa), Apptainer (Singularity) ile projenizi Docker konteyneri içinde çalıştırabilirsiniz. Şu anki `environment.yml` yapısı Conda için uygundur.

---

## 13. Sık Karşılaşılan Sorunlar ve Çözümleri

**S: İşimi gönderdim ancak hep PD (Pending) durumunda kalıyor.**
- **Ç:** İstediğiniz kaynaklar (özellikle GPU veya çok sayıda node) şu an dolu olabilir. `--time`, `--cpus-per-task` gibi değerleri sadece ihtiyacınız kadar talep edin. Ne kadar az kaynak isterseniz o kadar çabuk çalışmaya başlar.

**S: Disk kotası (Quota Exceeded) hatası alıyorum.**
- **Ç:** `home` dizininde büyük dosyalar veya çok sayıda küçük dosya (.conda cache) birikmiş olabilir. `du -sh *` komutu ile boyutları kontrol edin. Geçici dosyaları temizleyin veya `scratch` diskine taşıyın.

**S: Out Of Memory (OOM) hatası veya işin KILLED olması.**
- **Ç:** Betikte talep ettiğiniz RAM yetmemiştir. TRUBA'da bellek yönetimi katıdır. `#SBATCH --mem=64G` (veya daha fazla) parametresi ekleyerek daha fazla RAM isteyin.

**S: Conda ortamı çok yavaş kuruluyor.**
- **Ç:** TRUBA internet çıkışları bazen kısıtlıdır. Kurulum uzun sürse de bir kere yapılacaktır, bekleyin.

---

*Bu rehber, projenizin yapısı (`Makefile`, `environment.yml`) ve güncel TRUBA HPC standartları baz alınarak yeni başlayan bir araştırmacı için adım adım hazırlanmıştır.*
