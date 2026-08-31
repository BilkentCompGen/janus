# Deney 8 — ters indeks

2026-08-31. Kullanicinin fikri: `junctions` ve `pick` okumalara donuyor; hangi okumanin hangi
cekirdegi tasidigini saklarsak ikisi de taramaz.

## Ne yazildi

**`snape svbidx <in.svb>` -> `<in.svb>.idx`**: okuma basina dosya offseti. Kayitlar degisken
uzunlukta ve `.svb`'nin indeksi yoktu, yani `(ornek, okuma_no)` adreslenebilir degildi. Sadece
`u32 tag`'leri okuyup atlayarak kuruluyor. Olculdu: on ornegin hepsi, **880 MB toplam**, ornek
basina 2.2 s. En uzun okuma **63,407 bp** (hs578/n).

**`snape graph --index N --pad-cores P` -> `<dir>/index.bin`**: CSR postings.

```
char[4] "SNPI" | u32 ver | u32 S | u32 pad
u64 K | u32 N | u32 reserved
u32 off[K+1]
u64 post[2*off[K]]      // 16 bayt/kayit
```

Kayit paketi:

| bit | alan | neden yeter |
|---|---|---|
| 0-27 | read_no | 268 M; en buyuk ornek 13.4 M |
| 28-33 | sample | S <= 64, CLI'da zorlaniyor |
| 34-47 | ord | 16,383; 63 kb okumada ~650 cekirdek |
| 48-63 | corelen | 65,535; bir cekirdek okumadan uzun olamaz |
| w1 0-31 | st | okuma ici ham offset |
| w1 32-47 | dleft = st - ws | okuma uzunlugundan kucuk |
| w1 48-63 | dright = we - en | ayni |

`ws`/`we`, cekirdekten P=50 pozisyon oncesi/sonrasindaki cekirdegin sinirlari, yani
`junctions`'in kesim penceresi dogrudan kayittan cikiyor: `ws(ilk)`, `we(son)`,
`lf = dleft(ilk)`, `rf = dright(son)`.

**Yalniz ILERI iplik.** `junctions` ve `pick` okumalari tek iplik ayristiriyor; ters iplik
cekirdeklerinin koordinatlari rc cercevesinde ve o okumada kesim icin kullanilamaz. Ters aile
dugumleri, karsi yonde dizilenmis okumalarin ileri ayristirmasindan bulunuyor -- bugunku
davranisin aynisi.

**Determinizm.** Dugum basina EN KUCUK N okuma numarasi tutuluyor, ilk gorulen N degil. Shard
kilidi altinda sirali ekleme; dugum dolduktan sonra gelenler tek karsilastirmada eleniyor.
Sonuc thread sirasindan bagimsiz.

## Bir hata ve dersi

Ilk surumde `corelen`'e 12 bit (max 4,095) verilmisti; `docs/NOTES.md`'de olculen en uzun ham
cekirdek 2,474 bp idi. Tam olcekte asiliyor. Daha kotusu, tasma `abort()` ediyordu: **e2e'nin
bes graph asamasi da rc=134 ile dustu**, 40 dakika kayip.

Iki duzeltme: alan genislikleri yeniden dagitildi (`ord` 16->14, `corelen` 12->16 -- ikincisi
artik okuma uzunlugu ile sinirli oldugu icin **garantili**), ve tasma artik oldurmuyor, kayit
atlanip sayiliyor. Bir optimizasyon asla ana isi dusurmemeli.

## Olculen (400k okumalik altornek, K=32.8 M, N=8)

```
92,459,297 kayit, dugum basina ortalama 2.82, atlanan 0
dolu (cnt==8) dugum: 450,429 (%1.4)
capa testi (iki yanda >=300 bp) gecen kayit: 88,599,292 / 92,459,297 = %95.8
dugum basina >=5 aday: %11.73
somatik dugum: 12,689,415 (%38.7), bunlarin %7.84'unde >=5 aday, %5.38'inde 0 aday
```

Bu, tahminimi duzeltti: **capa testi kayitlarin %95.8'inde geciyor**, daha once soyledigim
%34.3 degil. O sayiyi yanlis turetmistim -- 355,026 kisa-capa reddi, capasi hic olmayan
chunk'lardan geliyor, kayit basina degil.

Bagliyaci kisit capa degil **derinlik**: bu altornekte kapsam ~2x, o yuzden dugumlerin ancak
%1.4'u N=8'i dolduruyor.

## Olculen (TAM OLCEK, hcc1954, K=70.4 M, N=8)

`tmp/idxfull`, 48:45 (nice 19, -t 8, e2e ile cekisirken), tepe RSS 18.52 GB, index.bin **8.52 GB**.

```
554,221,232 kayit, atlanan 0
dolu (cnt==8) dugum: 68,500,438 (%97.3)
capa testi gecen kayit: 531,049,256 / 554,221,232 = %95.8

dugum basina capa testini gecen aday:
   8 aday: 48,746,230 (%69.25)
   7 aday: 16,990,930 (%24.14)
   6 aday:  2,783,807 ( %3.95)
  >=5 aday: 68,954,665 (%97.95)

somatik dugum: 674,079 (%1.0)
  bunlarin %26.83'unde >=5 aday, %5.72'sinde 0 aday
```

**Cevap: N=8 `junctions` icin yetiyor.** Dugumlerin %97.95'inde capa testini gecen >=5 aday var,
`per_chunk` kapagi ise 5. Somatik dugumlerde oran daha dusuk (%26.83) ama `junctions` chunk
basina soruyor ve bir chunk ~10 cekirdek, yani birlesimde aday sayisi cok daha yuksek.

Bir de dikkat cekici: tam derinlikte somatik dugum yalnizca **%1.0**. `docs/NOTES.md`'deki
%28.17, 600k'lik bir altornekten geliyordu -- orada normal kapsami dusuk oldugu icin cok sayida
cekirdek "normal 3 kez gormedi" diye somatik sayiliyordu. Tam derinlikte gercek somatik %1.

**Yine de baglamadim.** `pick` tarafi acik: bugun lokusta ~44x goruyor, N=8 ile ~8x'e duser ve
olculmus zincir (11x -> 136 yanlis cagri, 44x -> 116) precision kaybi soyluyor. Uctan uca
dogrulama yapilmadan baglanmamali.

## Acik kalan

`junctions` ve `pick` henuz indeksi KULLANMIYOR. Kod yazilmadan once cevaplanmasi gereken:
tam olcekte N=8, chunk basina 5 contig icin yetiyor mu, ve `pick`'in ihtiyaci olan lokus
derinligini veriyor mu. `pick` bugun lokusta ~44x goruyor; olculmus: 11x -> 136 yanlis cagri,
44x -> 116. N=8 ile ~8x'e duserse precision'a mal olur. Bu, olculmeden baglanmamali.

## `junctions`'i indekse baglamanin tasarimi (yazilmadi, cikarilan gereksinimler)

Kodu yazmadan once bugunku davranisi tam olarak okudum; indeksle ayni ciktiyi vermek icin
gerekenler:

1. **Kayittan pencere dogrudan cikiyor.** Bir dizi `[a..b]` icin
   `ws = st(a) - dleft(a)`, `we = st(b) + corelen(b) + dright(b)`,
   `lf = dleft(a)`, `rf = dright(b)`. Dordu de kayitli alanlar, `pad` header'da.

2. **Ama yon ters cevrilmeli.** Bugun `junctions` okuma -> chunk yonunde gidiyor:
   bir okumanin butun somatik cekirdeklerini bulup maksimal diziyi cikariyor, ve
   `nxt - b > 2*PAD` testiyle **farkli chunk'lara ait** dizileri de birlestiriyor
   (`cc[b+1] >= 0` herhangi bir chunk demek, `ch = cc[a]` ilkini aliyor). Mobil eleman
   insersiyonunun iki birlesme noktasini ayni pencerede tutan sey bu.
   Indeks chunk -> okuma yonunde. Ayni davranisi vermek icin kayitlar once **okumaya gore**
   gruplanmali, sonra `ord`'a gore siralanip dizi cikarilmali. Yani indeks okunduktan sonra
   bir kez ters cevirme gerekiyor: (okuma -> [ord, chunk, alanlar]).

3. **Kapak ve isimlendirme.** `taken[chunk*L + line] < per_chunk` ve `<okuma>_<n>` numaralandirmasi
   okuma icinde artiyor; ikisi de yeni yonde kurulabilir ama sira bugunkunden farkli olacagi
   icin **hangi 5 contig'in secildigi degisebilir**. Contig kumesi bugunkuyle karsilastirilmali,
   birebir ayni olmasi beklenmemeli.

4. **Okuma baytlari.** `.svb.idx` ile dogrudan seek; okumalarin %70.1'i hic acilmaz
   (olculmus: 600k'lik altornekte contig ureten okuma orani %29.9).

Beklenen kazanc, kayitli butceye gore: `junctions` 2:06 -> ~1:00. Bu gecenin olcumunde
`junctions` 2:21.99 (hcc1954) ve 2:37.95 (h2009), yani kazanc ~1:15/cift mertebesinde.
Boru hattinin ~%3'u. `pick`'in kazanci daha buyuk olurdu (2:17.66 -> ~0:20) ama derinlik
sorunu yuzunden baglanamiyor.

**Sonuc: indeks kuruldu ve olculdu; baglanmasi bir sonraki isin konusu ve tek basina
boru hattinin ~%3'unu kazandiriyor. minimap2 %52.3, oncelik orada.**
