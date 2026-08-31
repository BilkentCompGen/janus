# Gece plani — 2026-08-31 03:13 baslangic

## Kullanicinin talimati (birebir)

> Ben yatiyorum! /autonomous-work skill i uyarinca gece boyu bahsettiklerimi yapmani ve denemeni
> ve 5 CASTLE sample inin da tum step lerini tamamlamani istiyorum. Yarin uyandigimda eksik kalan
> unutulmus yarida birakilmis is gormeyecegim. Ayrica bana soru da sormayacaksin sabah
> uyandigimda reverse index + 5 castle sample i e2e bitmis ve sureler ram ler etc. olculmus
> gorucem. Hatta sonrasinda vaktin kalirsa minimap2 yi hizlandirmanin yollarina bak
> (1. dogru flagler ve parallelism 2. daha az sayida veya daha kucuk okumalar 3. minimap2 lerde
> tekrar eden isler varsa eger gerekirse alignment i kendin yazip o ortak isleri ortak halletmek)

Onceki oturumda kararlastirilanlar da gecerli:
- Ters indeks: her dugum icin N=8 kayit, kayit = (sample, read_no, ord, st, corelen, dleft, dright),
  16 bayt, CSR (`off[K+1]` + postings). Ayri dosya (`index.bin`), `.svb` yaninda `.svb.idx`
  okuma-basi offset tablosu. Yalniz ILERI iplik kayitlari (junctions/pick tek iplik ayristiriyor).
- Determinizm: dugum basina EN KUCUK N okuma numarasi tutulur, ilk gorulen N degil.

## Yorum ve oncelik

Teslim edilmesi sart olan, azalan onemde:

1. **Bes CASTLE cifti, ucdan uca, olculmus.** Iki kez adiyla soylendi. Guvenilirligi her seyin
   onunde: dogrulanmis binary ile, sirayla, bos kutuda, her asama `/usr/bin/time -v` altinda.
2. **Ters indeks kurulmus ve olculmus.**
3. Vakit kalirsa minimap2.

Cakisma: temiz sure/RSS olcumu bos kutu ister; gelistirme CPU yer. Cozum, e2e once ve
kesintisiz kossun, indeks kodu onun uzerine yazilsin ve yalniz 400k'lik altornekte
(~40 s) denensin. E2e bitince indeks tam olcekte kosulur.

## Zaman butcesi (03:13 baslangic, ~08:00 hedef)

| | |
|---|---|
| 03:13-03:35 | bootstrap, truth arsivi indirmeye alindi, svbidx tamamlandi, e2e betigi |
| 03:35-07:45 | **bes cift e2e, sirayla, bos kutu** (~50 dk/cift) |
| paralel, CPU'suz | indeks kodu, docs, experiments |
| 07:45- | indeks tam olcekte, skorlama, HANDOFF |

## Durum notlari

- `graph` iki gecisli: tepe RSS 32.13 -> 10.24 GB, duvar 12:13 -> 16:59 (hcc1954, olculdu).
- `extract.cpp`'de `emit()`'in sinir disi okumasi duzeltildi; determinizm testi kosuyor.
- Truth setleri diskte yok (yalniz hcc1954). Zenodo 14541057 arsivi (15 GB) indiriliyor.
- Minda kurulu degil. Skorlama icin ya kurulur ya da ayni kurala (min_size 50, tolerance 500)
  gore kendi eslestiricim kullanilir ve boyle oldugu ACIKCA yazilir.
- hs578'in CASTLE'da truth'u yok; e2e kosulur, skorlanmaz.
