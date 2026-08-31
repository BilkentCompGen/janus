# Deney 6 — graph'i iki gecise ayirmak

2026-08-31. Kullanicinin fikri: sayimi, kenar/dizi uretiminden ayirip iki gecis yapmak.

## Gerekce

Tek gecisli kurulumda, tepe anindaki bayt dokumu (hcc1954, `[bytes]` satiri, olculdu):

```
sample-map    10.00 GB     bu ornegin tam sayim tablosu
sample-edges   8.00 GB     bu ornegin kenar kumesi
sample-seq     3.71 GB
graph-map      3.00 GB
graph-seq      3.97 GB
graph-seqref   2.00 GB
graph-edges    3.00 GB
graph-cov      0.25 GB
------------
              33.93 GB     (RSS 30.24, zirve 32.13)
```

`sample-edges`'in buyuk kismi hic yasamayacak cekirdeklerin kenarlari: hayatta kalma o asamada
belli degil, o yuzden hepsi ornek bitene kadar tutuluyor. Ayni sekilde `seq` her ilk-goruste
uretiliyor, yani emisyonlarin %97.2'sinde.

## Ne yapildi

**Gecis A** yalniz sayiyor. Ornek basina tablonun degeri `cov_t`, `SeqRef` degil (kayit 32 -> 16
bayt); ornek basina kenar kumesi hic yok.

**Gecis B** hayatta kalanlar belliyken kenarlari ve cekirdek baytlarini uretiyor. Kenar dogrudan
global kumeye giriyor. Bir kenarin o ornekte yasayip yasamadigi `cov[col] != 0` ile soruluyor --
tek gecisli surumun ornek basina `survives()` testiyle ayni sey, cunku `cov[s][i]` yalniz o ornek
esigi gectiyse yaziliyor. Kenar anahtari `(u32,u32)` ciftinin u64 paketi, `SEdge`'in 16 bayti
yerine 8. Paketleme emisyon basina degil DUGUM basina calisiyor.

## Olculen (hcc1954, 16 thread, bos kutu)

| | tek gecis | iki gecis |
|---|---|---|
| tepe RSS | 32.13 GB | **10.24 GB** (-68%) |
| duvar | 12:13 | 16:59 (+39%) |
| gecis A | — | 8.1 dk, zirve 8.24 GB |
| gecis B | — | 8.9 dk, zirve 10.24 GB |
| graph.bin | 4.35 GB | ayni |

CPU dagilimi, farkin nereden geldigi:

```
                     tek gecis   iki gecis
toplam CPU-s             10863       15324
  extract2                6399       12866    +6467
  geri kalan              4464        2457    -2007
net                                   +4461 CPU-s = +4.8 dk
```

Parse disi is **%45 azaldi**. Ama ikinci gecis bir tam `extract2` daha odetiyor, kazanilanin uc
kati. Yeni kosuda `extract2` toplam CPU'nun **%84'u** (eskiden %59).

## Dogrulama

- 400k okumalik altornekte cikti tek gecisli surumle **bayt bayt ayni**.
- Tam olcekte 70.4 M dugumun 5'i farkliydi, hepsi `min_mult=3` sinirindaki sayimlarda. Sebep
  `extract.cpp`'de bulundu ve duzeltildi (Deney 7). Duzeltmeden sonra iki tam kosu **bayt bayt
  ayni** (17:10.78 ve 17:07.97, 10.12 ve 10.21 GB).

## Ne kazandirmadi

Boru hattinin tepe RSS'i graph'in degil artik: iki minimap2 cagrisi olculmus 19.1 ve 20.9 GB.
Yani zirve 32.8 -> 20.9 GB (-%36), graph'in -%68'i degil. 10 GB'in altina inmek icin asil
kaldirac minimap2'de.
