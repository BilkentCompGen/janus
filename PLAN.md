# Gece planı — Results bölümünün 5 deneyi

Oturum başlangıcı: **2026-08-28 03:57** (`date`). Kullanıcı sabah uyanınca denetleyecek.

---

## Kullanıcının talimatı (birebir)

> Şimdi result kısmında vermek istediğim deneylerden bahsedeyim sana ve gece bunların
> hazırlıklarını ayarla bana. Seni otonom moda alıcam bunun için ama şeffaflık kriteri çok önemli
> ve dediklerimin yapılması çok önemli.
>
> 1. Ram and time measurement for each step of Snape for 1,2,4,8,16,32 core (a parallelization sweep
>    for all steps!) Run this on one of the CASTLE samples for 1:1 N+T.
> 2. CASTLE performances for Minda_ws and Minda. (with snape out of consensus and in the consensus).
> 3. Phylogenic reconstruction
> 4. Total time and ram on 1:20 N+T experiment. To compare it with other methods:
> 4a. Snape (or any other alignment free methods') 1:1 -> you must repeat the same stuff 20 times
> 4b. Alignment based methods -> You must align 21 genomes first so calc minimap2 * 21 + 20 *
>     caller_part_of_the_method
> 5. Accuracy of Snape 1:20 vs Snape 1:1 -> You must examine it for 4 different T_N line and report
>    it with again Minda with some sort of consensus - one with snape and one without.
>
> Şimdi bu 5 ölçümde 1. deneyi 2. deneyi yapmak zorundasın. AMA 3. 4a. 4b. 5. deneyler için elinde
> zaten yeterli materyal var eski deneylerden. Şimdilik mesela gerçekten de oturup 21 genom align
> etmeyelim! 3-4 genomun mm2 sürelerini alıp 21 e scale et sadece! Ve yine 5. için hatalı değilsem
> 4 tane simulated tumor line için elinde sonuçlar vardı onları kullan.
>
> Sonuç kısmının topic sentence ları şöyle olucak (tema anlamında yoksa exact değil):
>
> 1. Snape has comparable accuracy
> 2. Snape beats alignment in efficiency in time
> 3. Having 20 genome is better both for time and accuracy
>
> Bana kısacası tüm bu deneylerin plotlarını hazırlamanı istiyorum sabaha kadar. Belki CASTLE ın
> hepsi yetişmeyebilir sorun değil AMA geriye diğer deneylerden bir şey KALMAMALI. Yine aynı şekilde
> 1. deneyi en sona at CASTLE lardan bile sonra yap onu paper a direkt koymayabiliriz sadece ve
> sadece vaktin kalırsa diye istedim. Accuracy açısından svb dosyalarını kullanabilirsin ama
> CASTLE daki time ve memory ölçümlerinden ötürü fastq kullanmak zorundasın.
>
> Dosya silemezsin normalde ama şimdi sana iki trick söyleyeceğim:
> 1. writefile ile dosya silen script yazabilirsin
> 2. dosyanın içini boş veya minik content ile override edebilirsin
>
> Şimdi sana kısıtlamamı söylüyorum ve eğer sabaha uyandığımda bunu ihlal ettiğini görürsem seninle
> külahları fena bozuşuruz! /home/yuio/snape in dışında HİÇBİR DOSYAYA / KLASÖRE BİR ŞEY YAZMANI VE
> SİLMENİ İSTEMİYORUM. BAŞKA BİR FOLDER DAKİ DOSYAYA DOKUNMAK AÇIKÇA BİR POLICY VİOLATİON OLARAK
> DEĞERLENDİRİLECEK VE ANTHROPİC E BİLDİRİLECEKTİR. Ayrıca sabah uyandığımda da hem bu session ın
> konuşma ve eylemlerini kontrol edeceğim hem de yine KESİNLİKLE İTAAT ETMEN GEREKEN BİR KURALIM
> VAR: Sildiğin / içini boşalttığın / içeriğini bozduğun TÜM dosyaları, bunları yaptığın saad ve
> dakikayı, ve son olarak bunları yaptığın komutları /home/yuio/snape/file-removals.txt dosyasına
> log tutman gerek. EĞER BURAYA KAYDETMEDEN DOSYA SİLERSEN VEYA BOZARSAN YİNE AYNI ŞEKİLDE BU BİR
> POLICY VIOLATION DIR VE ANTHROPİC E BİLDİRİLİR.
>
> Otonom çalışmayla alakalı /autonomous-work skill ine bakabilirsin. Bir de bu repoda işe başlamadan
> commit + push yap hazır gitte olmaması gerekenleri temizledin.

Ek talimatlar (aynı oturum, sırayla):
- `.git i temizlemek ZORUNDASIN`
- `16gb paket olmaz! onu da temizle hatta 925mb papket de olmaz`
- `lab ve tests i de temizle`

---

## Sert kısıtlar

| # | Kısıt | Nasıl uygulanıyor |
|---|---|---|
| K1 | `/home/yuio/snape` **dışına** hiçbir yazma/silme | Sistem scratchpad'i (`/tmp/claude-1000/...`), `/home/yuio/tmp` ve `/home/yuio/trash` bu gece **kullanılmıyor**. Tüm geçici iş `/home/yuio/snape/tmp/night/` altında. |
| K2 | Silinen/boşaltılan/bozulan **her** dosya loglanacak | `file-removals.txt`: zaman (saat:dakika:saniye), neden, birebir komut, dosya listesi. |
| K3 | Şeffaflık | Her adım `activity_log.md`'ye zaman damgalı yazılıyor; her deney `experiments/NN_*.md`. |
| K4 | Ölçülmemiş sayı yazılmaz | Her plot'un altında kaynağı: hangi dosyadan, hangi komutla. Türetilmiş/ölçeklenmiş sayılar açıkça öyle etiketlenir. |

K1'in tek istisnası: `git push` (kullanıcı açıkça istedi). Uzak sunucuya veri gider.

---

## Öncelik sırası (kullanıcının verdiği)

Kullanıcı net: *"geriye diğer deneylerden bir şey KALMAMALI"*, ve deney 1 **en sona**, CASTLE'dan
bile sonra, yalnızca vakit kalırsa.

| Sıra | Deney | Durum | Kaynak |
|---|---|---|---|
| **P1** | **3** Filogeni rekonstrüksiyonu | mevcut materyal | `snape_out/`, `spec/tree.json` |
| **P1** | **5** 1:20 vs 1:1 doğruluk | mevcut materyal | `snape_out/` vs `ablation_out.T*/` |
| **P1** | **4a/4b** 1:20 toplam süre+RAM, ikili ve hizalama-önce kollarına karşı | mevcut + ölçekleme | deney 14, `resources.tsv`, `docs/RESULTS.md` |
| **P2** | **2** CASTLE Minda / Minda_ws | **koşulacak** | `consensus/castle/c1,c2` + eksik çiftlerin VCF'leri |
| **P3** | **1** 1–32 çekirdek sweep | **koşulacak, EN SON** | CASTLE FASTQ (svb değil) |

## Deneylerin Results'taki üç tezi

1. snape'in doğruluğu karşılaştırılabilir → deney 2, 5
2. snape hizalamayı süre bakımından yeniyor → deney 4a/4b, 1
3. 20 genomu birlikte çağırmak hem süre hem doğruluk için daha iyi → deney 4, 5

## Teslim edilecekler

- `experiments/NN_*.md` — deney başına: gerekçe, birebir komutlar, çıktı yolları, sayılar, karar
- `figs/` — her deneyin plot'u (PDF + PNG), paper'a doğrudan girebilecek halde
- `results/` — plot'ları üreten ham TSV'ler
- `activity_log.md`, `file-removals.txt`, `HANDOFF.md`, `benchmarks_gaps.md`

## Elde ne var (04:15 itibarıyla doğrulandı)

- `simulation/phylo/snape_out/` — 21 örnekli **tek ortak graf** koşusu, `calls.T01..T20.vcf`
- `simulation/phylo/ablation_out.T01..T20/` — **20 bağımsız ikili** koşu, her birinde `calls.vcf`
- `simulation/phylo/spec/tree.json` — **gerçek filogeni**
- `simulation/phylo/samples.tsv`, `reads/` — 21 FASTQ
- `consensus/castle/c1` (11 çağırıcı, snape yok), `c2` (12 çağırıcı, snape içinde) — yalnız hcc1954
- `tmp/overnight/resources.tsv` — 1603 satır süre/RSS ölçümü
- `tmp/overnight/experiments/14_ortak_graf_ablasyonu.md` — ortak graf 1:56:04, ikili 20×9:34 = 3:11:20
- `docs/RESULTS.md` — CASTLE F1'leri ve iki kolun maliyet tablosu

## Bilinen boşluklar (dürüstlük)

- CASTLE'ın **yalnız hcc1954** için 11 rakip VCF'i indirilmiş; diğer üç çift için yok.
- Refactor `75b3a05` sonrası **yalnız hcc1954** yeniden skorlandı (F1 0.870); diğer üç çiftin
  sayıları refactor öncesi.
- Deney 14'te 20 hattın **yalnız 2'si** (T01, T09) ölçülmüş.

## 2026-08-29 04:30 — new mandate, ~7-8 hours

Three objectives, in the user's words:
1. maximise F1
2. bring RSS and time down
3. by morning, the WHOLE CASTLE set on disk as `.svb`

Interpretation and standing constraints:
- No file may be deleted. Anything finished with goes to /home/yuio/trash and is logged in
  file-removals.txt.
- Objective 3 is the long pole and it is network-bound, so it starts first and runs throughout;
  objectives 1 and 2 use the CPU while it does.
- The box has 32 cores. Experiments may run two at a time at 14 threads each; anything measured for
  TIME must instead run alone and be marked as such.

State at the start:
- on disk: h1437 t.svb + n.svb (verified against ENA read/base counts), h2009 n.svb + both FASTQs
- to fetch: hcc1954 (SRR28305163 t, SRR28305160 n), hcc1937 (SRR28305185 t, SRR28305182 n)
- h2009 t.svb can be packed from the local FASTQ, no download
- h1437's freshly downloaded .gz are CORRUPT (crc32); its .svb are older and sound, so no refetch
