# Deney 7 — graph.bin neden deterministik degildi

2026-08-31. Iki gecisli graph tam olcekte tek gecisliyle **birebir ayni cikmadi**: 70.4 M dugumun
5'i farkliydi ve `K` 1 fazlaydi. Once bunun benim degisikligim oldugunu varsaydim; degildi.

## Belirti

hcc1954, iki gecisli cikti ile o geceki tek gecisli referans:

```
K            70,395,067  vs  70,395,068
ov_s, ov_c   birebir ayni
okuma sayilari, ornek basina ayrik imza sayilari  birebir ayni

yalniz yenide 2 dugum, yalniz eskide 1
ortak dugumlerde cov farki: duzlem 0'da 4 dugum, duzlem 1'de 1
farkli degerler: (6,5) (3,4) (3,0) (3,4) ve (3,0)
```

Hepsi 0..6 araliginda, yani `min_mult = 3` sinirinda. Doygunluk degil.

## Sebep

`extract.cpp`'nin `emit()`'i:

```c
const uint64_t a = cores[i].start, b = cores[i].end + 1;
s.st[i] = offset + map[a];
s.en[i] = offset + map[b];     // <-- b == hn+1 olabilir
```

`map`'in boyutu `hn + 1` (indeksler 0..hn, `map[hn] = len` sentineli). Bir strandin **son
cekirdeginin** `end`'i `hn` olunca `b = hn + 1` ve `map[b]` **sinir disi**. Vektor `reserve(len+1)`
ile ayrildigi icin cokme yok; okunan sey ayrilmis ama yazilmamis alanda duran, **ayni thread'in
onceki okumasindan kalan** sayi. O sayi `en[i]` olarak margin testine giriyor
(`en[i] + margin > len`), yani sinirdaki cekirdegin sayilip sayilmadigini thread'e is dagitimi
belirliyor.

Dosyanin kendi basligi bunu zaten yaziyordu ("one entry past the run map"), ama sonucunun
determinizm oldugu fark edilmemis.

## Olcum

`tmp/prof/t_lastcore.cpp`, 20,000 gercek h1437 okumasi, seviye 5, tek iplik:

```
3,351,577 cekirdek
map[end+1] sinir disi: 158 cekirdek (%0.0047), 158 okumada (%0.8)
```

`extract2` iki iplik yaptigi icin gercekte ~%1.6 okuma. hcc1954'un 22.9 M okumasinda kosu basina
~366 bin sinir disi okuma; cogu cop deger margin testini zaten gecemeyecek kadar buyuk, birkaci
kucuk cikip cekirdegi tutuyor.

## Duzeltme

```c
s.en[i] = offset + map[b < map.size() ? b : map.size() - 1];
```

Kirpilinca `en = len` oluyor ve `en + margin > len` her zaman dogru: okumanin sonuna dayanan
cekirdek her zaman eleniyor. Zaten margini olmayan bir cekirdek.

## Dogrulama

Ayni girdi, ayni binary, iki tam kosu (hcc1954, 16 thread):

```
kosu 1  17:10.78  10.12 GB
kosu 2  17:07.97  10.21 GB
graph.bin BAYT BAYT AYNI
```

Referanstan (duzeltmeden onceki kosu) farkli, beklendigi gibi sinirdaki birkac dugum kadar.

## Not

`docs/NOTES.md`'deki "graph tekrar uretilebilir" testi 600k'lik bir **altornekte** yapilmisti;
tam olcekte hic test edilmemisti. Kullanici da tam olcekte deterministik olmadigini zaten
biliyordu. Bu deney sebebini buldu ve kapatti.
