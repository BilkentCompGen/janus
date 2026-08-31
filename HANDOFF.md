# HANDOFF — 2026-08-31 gecesi

Bu dosya gecenin tamamini ozetler. Her sayi olculmustur; olculmemis olanlar acikca oyle yaziyor.
Ayrinti: `experiments/06-09`, ham log'lar `e2e/` ve `tmp/` altinda, zaman cizelgesi
`activity_log.md`.

## Bir bakista

| | durum |
|---|---|
| Bes CASTLE cifti, ucdan uca | **BITTI**, besi de. 3 sa 34 dk. Tablo asagida |
| Ters indeks | **kuruldu ve olculdu**. Tuketicilere (junctions/pick) **baglanmadi**, sebebi asagida |
| graph iki gecisli | **bitti**: tepe RSS 32.13 -> 10.24 GB, duvar 12:13 -> 16:59 |
| graph determinizmi | **duzeltildi ve dogrulandi**: iki tam kosu bayt bayt ayni |
| minimap2 | `.mmi` kuruldu ve dogrulandi; thread + kirpma deneyleri kosuyor |
| Kod | `clean` dalinda, push edildi |

## 1. graph iki gecise ayrildi

Kullanicinin fikri. Sayim ile kenar/dizi uretimi ayrildi; ornek basina kenar kumesi (olculen
8 GB) tamamen ortadan kalkti, paketleme emisyon basina degil dugum basina calisiyor
(7.9 milyar yerine 70.4 milyon).

| hcc1954, 16 thread | tek gecis | iki gecis |
|---|---|---|
| tepe RSS | 32.13 GB | **10.24 GB** (-68%) |
| duvar | 12:13 | 16:59 |
| gecis A / B | — | 8.1 dk / 8.9 dk |

Bedelin tamami tek bir sey: ikinci bir `extract2`. Parse disi is %45 azaldi (-2,007 CPU-s) ama
parse +6,467 CPU-s ekledi. Yeni kosuda `extract2` toplam CPU'nun %84'u.

400k altornekte cikti tek gecisliyle **bayt bayt ayni**.

Detay: `experiments/06_iki_gecisli_graph.md`.

## 2. graph.bin neden deterministik degildi

`extract.cpp`'nin `emit()`'i son cekirdek icin `map[end+1]` okuyordu; `map`'in boyutu `hn+1`,
yani **sinir disi**. Okunan sey ayni thread'in onceki okumasindan kalan sayi, ve o sayi `en[i]`
olarak margin testine giriyor. Olculdu: 20,000 okumada 158 cekirdek (%0.8 okuma, tek iplikte).

Tek satirlik duzeltme (indeksi son gecerli girdiye kirp). Dogrulama: ayni girdi, iki tam kosu,
**bayt bayt ayni** (17:10.78 / 10.12 GB ve 17:07.97 / 10.21 GB).

Detay: `experiments/07_determinizm.md`.

## 3. Ters indeks

`snape svbidx` -> `.svb.idx` (okuma basina offset; on ornek, 880 MB, ornek basina 2.2 s).
`snape graph --index N` -> `index.bin` (CSR postings, kayit 16 bayt, dugum basina en kucuk N
okuma -- thread sirasindan bagimsiz).

Olculdu (hcc1954, tam olcek, N=8): **554,221,232 kayit, dugumlerin %97.3'u dolu, 8.52 GB**,
alanina sigmayan kayit 0.

**junctions icin yetiyor**: dugumlerin %97.95'inde capa testini (iki yanda >=300 bp) gecen >=5
aday var, `per_chunk` kapagi 5.

**pick icin acik**: `pick` bugun lokusta ~44x goruyor; N=8 ile ~8x'e duser. Olculmus zincir:
3x -> 290 yanlis cagri, 11x -> 136, 44x -> 116. Yani hiz icin precision odenir.

**Bu yuzden tuketicilere baglamadim.** Uctan uca dogrulama yapmadan baglamak cagrilari sessizce
bozardi. Yapilmasi gereken siradaki is bu: junctions'i indeksle kosup contig kumesini bugunkuyle
karsilastirmak, ve pick icin ya daha buyuk N ya da taramayi korumak.

Detay: `experiments/08_ters_indeks.md`.

## 4. minimap2

Boru hattinin **%52.3'u**. Olculen baslangic: contig'ler 2.04 Gbp / 4:00 = 8.50 Mbp/s
(ortalama sorgu 9,173 bp); normal okumalar 3.79 Gbp / 16:25 = 3.85 Mbp/s (17,814 bp). Baz
basina 2.21 kat fark, uzunluk orani 1.94 -- maliyet sorgu uzunlugunda **superlineer**.

- **Paralellik**: `align()` `-t (threads-2)` kullaniyor, `--threads 16` ile 14. Kutuda 32
  cekirdek var, 16'si hic kullanilmiyor. Olculmemis, ama olcekleme kayitli (%1349 CPU / 16 thread).
- **`.mmi`**: her cagri hg38'i bastan indeksliyor (~25-33 s). Bir kez kurulup paylasildi;
  **cikti birebir ayni oldugu dogrulandi** (2,097 hizalama, bayt bayt).
- **Kirpma**: `pick` okumalari butun yaziyor (4.551 Gbp, ortalama 17,800 bp) ama germline testi
  cagrinin +-500 bp'sine bakiyor. `pick`'e `trim_bp` eklendi; sweep (0/3000/1500/500) **kosuyor**,
  her kol icin hem sure hem F1 olculuyor.
- **Yapmadigim, ve neden**: referansin yalniz aday lokuslarindan kucuk bir indeks kurmak.
  MAPQ indeksteki alternatiflere gore hesaplaniyor; kucuk indekste sisirilir ve germline testinin
  hem `qual < 10` suzgeci hem "normal burayi yerlestiremiyor" maddesi bozulur. Bu maddenin
  olculmus degeri var. Once kirpma ve thread, bu en sona ve F1 olcumuyle.

Detay: `experiments/09_minimap2.md`.

## 5. Bes cift, ucdan uca

Skorlama notu: **Minda bu kutuda kurulu degil.** Ayni iki kural (min_size 50, tolerance 500)
kendi eslestiricimizle uygulaniyor (`e2e/score.py`). Kalibrasyon: Agustos'un hcc1954 cagrilari
bu betikle 0.891/0.862/0.876, Minda ile kayitli 0.895/0.847/0.870.

Truth setleri diskte yoktu; Zenodo 14541057 arsivi (15 GB) indirildi, `minda_4_2` ensemble'lari
`truth/` altinda (1106 / 925 / 703 / 221 kayit, `docs/RESULTS.md` ile birebir).

Tablo `e2e/collect.sh` ile uretiliyor; sonuclar bittikce buraya yazilacak.

### Bes cift, hepsi bitti (04:02:59 -> 07:36:56, 3 sa 34 dk)

Asama basina duvar suresi:

| cift | graph | bubbles | junctions | mm2 ctg | call aday | pick | mm2 normal | call | TOPLAM |
|---|---|---|---|---|---|---|---|---|---|
| hcc1954 | 19:43.91 * | 0:06.58 | 2:21.99 | 7:49.33 * | 0:09.64 | 2:17.66 | 19:58.14 * | 0:21.05 | **52:48** |
| h2009 | 23:38.00 * | 0:09.36 | 2:37.95 | 17:45.03 | 0:23.06 | 2:06.68 | 12:28.70 | 0:37.04 | **59:45** |
| hcc1937 | 15:57.15 | 0:05.95 | 1:48.95 | 4:44.72 | 0:08.90 | 1:43.57 | 8:39.04 | 0:13.43 | **33:21** |
| h1437 | 18:38.81 | 0:05.90 | 1:56.98 | 9:08.74 | 0:11.93 | 2:11.76 | 9:41.69 | 0:21.45 | **42:17** |
| hs578 | 11:35.76 | 0:05.86 | 1:43.41 | 8:34.53 | 0:11.01 | 0:49.67 | 2:25.97 | 0:13.65 | **25:39** |

`*` es zamanli bir isle olculdu (tam olcekli indeks kurulumu / kirpma deneyi). Temiz karsiliklar:
hcc1954 graph 17:10.78, hcc1937 graph 15:57.15, hs578 graph 11:35.76.

Asama basina tepe RSS (GB):

| cift | graph | bubbles | junctions | mm2 ctg | pick | mm2 normal | ZIRVE |
|---|---|---|---|---|---|---|---|
| hcc1954 | 10.30 | 6.10 | 1.84 | 18.30 | 0.11 | 19.71 | **19.71** |
| h2009 | 10.41 | 6.11 | 1.87 | 18.96 | 0.10 | 19.17 | **19.17** |
| hcc1937 | 10.64 | 6.06 | 1.84 | 18.04 | 0.10 | 19.26 | **19.26** |
| h1437 | 9.77 | 6.09 | 1.85 | 18.36 | 0.11 | 19.43 | **19.43** |
| hs578 | 11.32 | 6.03 | 1.83 | 18.27 | 0.09 | 18.35 | **18.35** |

**Zirve artik minimap2'nin.** Iki gecisli graph 32.8 -> ~10.3 GB indirdi; boru hattinin tavani
19.7 GB'da ve o minimap2. Daha asagi inmek isteniyorsa kaldirac orada.

Sure dagilimi: iki minimap2 cagrisi toplam surenin **%52.4'u** (kayitli %52.3 ile ayni).

Cikti buyuklukleri:

| cift | graph.bin | ctg.t.fa | ctg.t.bam | normal.bam | cagri |
|---|---|---|---|---|---|
| hcc1954 | 4.35 GB | 3.10 GB | 0.61 GB | 0.97 GB | 1066 |
| h2009 | 4.36 GB | 8.09 GB | 1.48 GB | 1.42 GB | 960 |
| hcc1937 | 4.33 GB | 3.12 GB | 0.62 GB | 0.59 GB | 605 |
| h1437 | 4.35 GB | 4.08 GB | 0.76 GB | 1.09 GB | 265 |
| hs578 | 4.30 GB | 3.22 GB | 0.67 GB | 0.25 GB | 132 |

### Dogruluk

| cift | cagri | truth | P | R | **F1** | tarihsel F1 | fark |
|---|---|---|---|---|---|---|---|
| hcc1954 | 1066 | 1106 | 0.894 | 0.862 | **0.878** | 0.876 | +0.002 |
| h2009 | 960 | 925 | 0.848 | 0.880 | **0.864** | 0.865 | -0.001 |
| hcc1937 | 605 | 703 | 0.830 | 0.714 | **0.768** | 0.765 | +0.003 |
| h1437 | 265 | 221 | 0.668 | 0.801 | **0.728** | 0.754 | **-0.026** |
| hs578 | 132 | — | — | — | — | — | CASTLE'da truth yok |

Uc ciftte tarihselle ayni. h1437'deki -0.026, gecenin basinda konusulan bilinen sey: ters iplik
dilim duzeltmesinin bu cifte maliyeti. Yani iki gecisli graph, `emit()` duzeltmesi, yeni bubble
kurali ve yorumlarin silinmesi **dogrulugu bozmadi**.

Gorsel rapor: `results/gece_20260831.html` -> http://192.168.1.174:8799/gece_20260831.html

### hcc1954 asama detayi (ilk bitirilen)


| asama | duvar | tepe RSS |
|---|---|---|
| graph | 19:43.91 * | 10.30 GB |
| bubbles | 0:06.58 | 6.10 GB |
| junctions | 2:21.99 | 1.84 GB |
| mm2 contig | 7:49.33 * | 18.30 GB |
| call (aday) | 0:09.64 | 0.05 GB |
| pick | 2:17.66 | 0.11 GB |
| mm2 normal | 19:58.14 * | 19.71 GB |
| call | 0:21.05 | 0.05 GB |

`*` = tam olcekli indeks kurulumu es zamanli kosarken olculdu, cekismeli. graph'in temiz
karsiligi 17:10.78 / 10.12 GB (determinizm kosusu).

**P 0.894 / R 0.862 / F1 0.878** (Agustos, eski binary, ayni skorlayici: 0.891/0.862/0.876).

## Yarin ilk yapilacaklar

1. `bash e2e/collect.sh` -- bes cift icin asama basina sure/RSS ve skor tablosu.
2. Kirpma sweep'inin F1 sonucu: kazanc precision'a mal oluyor mu.
3. Thread taramasi (bos kutu gerek): `-t 14` -> `-t 30`, beklenen ~%26 boru hatti.
4. Indeksi `junctions`'a baglayip contig kumesini bugunkuyle karsilastirmak.
