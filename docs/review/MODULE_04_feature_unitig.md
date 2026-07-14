# Modül 4 — Özellik Çıkarımı (Unitig)

> Mercek: **Python'da kal + sertleştir** · **ESKAPEE ölçeklenebilirlik** · **su geçirmez, bilimsel doğruluk**
> İncelenen: `scripts/03u_unitig_matrix.py` (bcalm2/unitig-caller → matris), matris sözleşmesi (`features.txt` / `y_{ab}.csv` / `genomes_{ab}.csv` / `X_{ab}_part_*.npz`), `08_blast_pipeline.nf` (sadece varlık; BLAST → M6)
> Tarih: 2026-07-13

---

## 1. Genel Değerlendirme

Bu modül tezin **asıl özelliğini** üretir: unitig × genom ikili matrisi (03u, unitig-caller/Bifrost). **Tasarım güçlü ve ESKAPEE-verimli:** env-parametrik (get_target ✓), organizma-seviyesi **unitig store** (unitig-caller BİR kez tüm organizma için → her antibiyotik store'u ucuza subset'ler, saatlerce yeniden çağırma yok), absolute min_support (ROADMAP risk-4), ve iki yol (subset / rtab) **aynı matris sözleşmesini** üretiyor. rtab ayrıştırma loud-fail (sample mismatch, malformed satır).

**Ama bir HIGH blocker var: gram-pozitif/klonal organizmalarda Faz1'i çökerten satır-sıralaması sorunu** — ESKAPEE genişlemesinin (projenin yeni yönü) önündeki tek somut engel. HANDOFF'ta doğrulanmış (S. aureus'ta yaşandı), fix hazır ama **uygulanmamış.**

---

## 2. Güçlü Yanlar

- **Env-parametrik:** `get_target` (CLI > AMR_ORGANISM/AMR_ANTIBIOTIC env > config) → paralel per-(organizma,antibiyotik) koşu.
- **Organizma-seviyesi unitig store (`--build-db`):** unitig-caller ALL genomes üzerinde bir kez → `unitig_all/` store; her antibiyotik `subset_store_to_antibiotic` ile store'u subset'liyor (unitig-caller re-run YOK). ESKAPEE'de çok-antibiyotik-per-organizma maliyetini çökertir.
- **İki yol, tek sözleşme:** `rtab_to_chunks` (fallback) ve `subset_store_to_antibiotic` (fast) ikisi de aynı çıktı sözleşmesini (`features.txt`/`y_{ab}.csv`/`genomes_{ab}.csv`/`X_{ab}_part_*.npz`) satırları **valid_genomes sırasında** yazıyor → tutarlı, downstream 04/05/06/07b hiç değişmeden okuyor.
- **Absolute min_support** (proportional değil — risk-4): `min_support ≤ support ≤ n-1` (zero-variance core + singleton dropped). Nadir-ama-gerçek plazmid unitig'i korur.
- **Sağlam rtab parsing:** rtab kolon sırasından bağımsız (header'dan `row_of_rtabcol`), sample-set mismatch + malformed satırda **loud sys.exit** (sessiz matris bozulması yok). `dtype=str` Genome ID.
- **resolve_tool** unitig-caller için (PATH-önce + env override).

---

## 3. Problemler (önem sırasıyla)

### Critical
- (yok)

### High
- **M4-1 — Fenotip-bloklu satır sırası → chunk-split tek-sınıflı fold → XGBoost NaN (ESKAPEE BLOCKER):** `select_genomes` (satır 124) genomları **metadata sırasında** döndürüyor. Klonal organizmalarda (ör. MRSA — genom ID'leri fenotipe göre kümeleniyor) bu sıra **fenotip-bloklu** (önce hep R, sonra hep S). 04/05/06'nın **contiguous chunk-split**'i o zaman tek-sınıflı chunk üretiyor → XGBoost NaN → Faz1 çöküyor. **Doğrulandı:** HANDOFF S. aureus'ta bunu yaşamış; **gram-pozitif ESKAPEE genişlemesini (projenin yönü) bloke ediyor.**
  - **Fix (hazır):** `select_genomes` return'ünden ÖNCE `valid_genomes` + `valid_labels`'a deterministik (seed 42) `np.random` permütasyon uygula (ikisini birlikte). **Güvenli:** her iki matris yolu satırları valid_genomes sırasında yazdığı için shuffle her iki yolda da propagate olur; lineage-CV **etkilenmez** (o `genomes_{ab}.csv` sırası + PopPUNK gruplarıyla yeniden hizalanır, satır-sırası-değişmez). Deterministik → tekrarlanabilir. Mevcut 21 modeli **etkilemez** (yeniden kurulmuyor).

### Low
- **M4-2 — `np.fromstring(..., sep="\t")` (satır 212) DEPRECATED:** modern numpy'de uyarı veriyor, ileride kaldırılabilir. Su-geçirmez/gelecek-uyum için `np.frombuffer` / `np.array(line.split("\t"), dtype=np.int8)` gibi bir alternatife çevir.
- **M4-3 — 08_blast_pipeline.nf tek Nextflow parçası:** BLAST orkestrasyonu burada (asıl inceleme M6). Python'da-kal kararıyla, bu tek .nf'in Python-orkestrasyona mı taşınacağı yoksa kalıp mı kalacağı M6/M9'da netleşecek.

---

## 4. Düzeltilmesi Gerekenler (madde madde)

1. **M4-1 shuffle fix'i uygula** — `select_genomes`'a seed-42 permütasyon (ESKAPEE gram-pozitif'i açar). Bir smoke/unit test ekle: fenotip-bloklu girdi → shuffle sonrası hiçbir contiguous chunk tek-sınıflı değil.
2. **M4-2** `np.fromstring` → non-deprecated parse.
3. (M6/M9) 08_blast_pipeline.nf kararı.

---

## 5. Refactor Önerileri

- **Shuffle'ı stratified yapma seçeneği:** düz seed-42 shuffle yeterli (sınıfları chunk'lara dağıtır); istenirse stratified-interleave daha da garanti eder ama gerekmez. Basit tut.
- **Matris sözleşmesini tek yerde belgele** (docstring'de var; bir `docs/DATA_CONTRACT.md`'ye taşımak downstream modülleri netleştirir — opsiyonel).

---

## 6. Bilimsel Eksikler (makale açısından)

- **Shuffle yalnızca yardımcı adımları (04/05/06 chunk-split) düzeltir; rapor edilen metrik (auc_mean_seeds = lineage-CV) satır-sırasından bağımsızdır** — bu Methods'ta netleştirilmeli ki "shuffle sonucu değiştirdi mi?" sorusu önden kapansın. (Cevap: hayır, sadece dejenere tek-sınıflı fold'u önler.)
- unitig-caller/Bifrost sürümü provenance'a (pipeline_runs) — CARD gibi araç sürümü kayıt altında mı, M7'de doğrula.

---

## 7. Literatür Gereksinimi

**Gerekmiyor** — mühendislik + yerleşik yöntem (unitig-caller, sparse chunking). M4-1 fix'i literatür beklemiyor; onayınla uygulanır.

---

## Uygulama durumu (2026-07-13) — UYGULANDI
- **M4-1 [ÇÖZÜLDÜ]:** `select_genomes` return öncesi seed-42 deterministik shuffle (`np.random.RandomState(42).permutation`; genome+label hizalı) uygulandı → klonal/fenotip-bloklu sırada tek-sınıflı chunk önlendi. ESKAPEE gram-pozitif (S. aureus) + A. baumannii önü açıldı. Lineage-CV etkilenmez (satır-sırası-değişmez); mevcut 21 model yeniden kurulmadığı için değişmez.
- **M4-2 [ÇÖZÜLDÜ]:** deprecated `np.fromstring(...sep=)` → `np.frombuffer(strip(tabs).encode) - ord('0')` (fast, non-deprecated). Aynı 0/1 int8 semantiği (sum/nonzero korunur).
- **M4-3 → M6** (08.nf konsolidasyonu), **M4-4** run-time (küçük set min_support).

**Test durumu:** 03u syntax OK, shuffle determinizmi + frombuffer parse doğrulandı, tüm suite yeşil (0 yeni kırık).

## Sonraki modüllere taşınan notlar
- **M5 (04/05/06/07b):** M4-1 shuffle 04/05/06'nın chunk-split'ini düzeltir; 07b lineage-CV'yi doğrula (+ M3-3 soy-CV fallback loud+KB-flag burada).
- **M7 (KB/provenance):** unitig-caller sürümü pipeline_runs'ta kayıtlı mı.
- **M6:** 08_blast_pipeline.nf incelemesi.
