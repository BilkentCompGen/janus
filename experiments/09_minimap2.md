# Deney 9 — minimap2, boru hattinin %52.3'u

2026-08-31. Kullanicinin uc yonu: (1) dogru bayraklar ve paralellik, (2) daha az veya daha kucuk
okumalar, (3) iki cagrida tekrar eden isi ortak yapmak.

## Baslangic noktasi (olculmus, hcc1954, `docs/NOTES.md`)

| cagri | girdi | duvar | tepe RSS | pay |
|---|---|---|---|---|
| contig'ler | 222,387 kayit / 2.04 Gbp | 4:00 | 19.1 GB | %10.3 |
| secilmis normal okumalar | 212,754 kayit / 3.79 Gbp | **16:25** | 20.9 GB | **%42.1** |
| toplam | | 20:25 | | **%52.3** |

Baz basina:

```
contig  2.04 Gbp / 240 s = 8.50 Mbp/s   ortalama sorgu  9,173 bp
normal  3.79 Gbp / 985 s = 3.85 Mbp/s   ortalama sorgu 17,814 bp
                                        baz basina 2.21x yavas, uzunluk orani 1.94x
```

Maliyet sorgu uzunlugunda **superlineer** (zincirleme is, capa sayisinda kare-ye yakin). Bu, uc
yonun hepsini birden isaret ediyor: kisa sorgular hem daha az baz hem de baz basina daha ucuz.

## Yon 1 — paralellik

`align()` `-t (threads - 2)` kullaniyor; `--threads 16` ile **14**. Kutuda **32 cekirdek** var.
Olculmus olcekleme: minimap2 + sort + index 16 thread'de %1349 CPU, yani 13.5 cekirdek dolu.
Bos kalan 16 cekirdek bugun hic kullanilmiyor.

Beklenti: `-t 30` ile mm2 20:25 -> ~10:00. Tek bayrak, boru hattinin ~%26'si. **Olculecek.**

## Yon 2 — daha kucuk sorgular

`pick` normal okumalari **butun** olarak yaziyor, ortalama 17,814 bp. Ama germline testi cagri
konumunun yalnizca **+-500 bp**'sine bakiyor (`tol * 10`), ve okuma zaten contig'in germline
capa cekirdeklerini paylastigi icin secildi -- yani ilgili bolge okumanin kucuk bir parcasi.

`pick`'e `trim_bp` parametresi eklendi: eslesen capa cekirdeklerinin en kucuk `st`'sinden en
buyuk `en`'ine kadar olan araligi alip iki yana `trim_bp` ekliyor. `trim_bp=0` bugunku davranis.

Beklenti: 3.79 Gbp -> ~1 Gbp, ve kisa sorgu baz basina da ucuz. 16:25 -> ~4:00.

**Risk, olculecek:** kisa okuma tekrar bolgelerinde daha dusuk MAPQ aliyor, germline testi
`qual < 10` ile eliyor. Yani hiz kazanci precision'a mal olabilir. Kirpma sweep'i F1 ile
birlikte olculmeli, tek basina sure ile degil.

## Yon 3 — ortak is

**(a) Referans indeksi iki kez kuruluyor.** Her cagri hg38 FASTA'sini okuyup minimizer indeksini
bastan kuruyor; olculdu, ~25 s (`mm_idx_stat::24.123`). `minimap2 -d` ile bir kez `.mmi` kurulup
iki cagrida da kullanilabilir. Kazanc kucuk (20:25'in ~50 s'si, %4) ama bedava ve tepe RSS'i de
biraz dusuruyor. Kuruldu, olculecek.

**(b) Iki cagri AYNI lokuslara hizaliyor.** Contig'ler junction'larin etrafindan kesiliyor;
secilmis normal okumalar da o contig'lerin capa cekirdeklerini paylastiklari icin seciliyor.
Yani ikinci cagri, birincinin ziyaret ettigi ~2,400 lokusun referans tarafini bastan isliyor.
Referansin yalniz o pencerelerinden (+-50 kb) kucuk bir indeks kurulsa 3.1 Gb yerine ~240 Mb
olurdu.

**Bunu kosmadim ve dikkatli olunmali:** MAPQ, indeksteki alternatiflere gore hesaplaniyor.
Kucuk bir indekste her okuma "benzersiz" gorunur, MAPQ sisirilir, ve germline testinin
`qual < 10` suzgeci ile "normal burayi yerlestiremiyor" maddesi ikisi birden bozulur. Bu
maddenin olculmus degeri var (hcc1954'te 2-destek ensemble'a karsi, hicbir seyin desteklemedigi
cagrilarin %7.4'u kirli lokusta, dogru cagrilarin %0.0'i). Yani hiz icin dogrudan bir precision
mekanizmasini feda etme riski var. Onerilen sira: once yon 1 ve 2, bu en sona ve F1 olcumuyle.

## Durum

`e2e/mm2_experiments.sh` yazildi ve kuyrukta: `.mmi` kurulumu (bitti), FASTA-vs-mmi cikti
esitligi, thread taramasi (14/24/30), ve normal okumalarda kirpma taramasi (0/1500/3000).
Sonuclar buraya yazilacak.

## Olculen sonuclar (hcc1954, bos kutu, 2026-08-31)

### Paralellik: OLU

| kol | duvar | RSS | CPU |
|---|---|---|---|
| contig -t 8 | 8:27.76 | 15.20 GB | 783% |
| contig -t 14 | **6:48.95** | 17.40 GB | 1360% |
| contig -t 22 | 6:55.96 | 20.22 GB | 2135% |
| contig -t 30 | 6:58.06 | 23.11 GB | 2891% |
| normal -t 14 | **16:16.16** | 18.60 GB | 1381% |
| normal -t 30 | 17:08.85 | 25.55 GB | 2893% |

14'ten sonra olcekleme yok: CPU iki katina cikiyor, duvar sabit hatta kotulesiyor, **ve RSS
buyuyor**. Boru hattinin tavani zaten minimap2'nin RSS'i oldugu icin bu aktif zarar. Cikti
`/dev/null`'a gidiyordu, yani I/O degil -- bellek bant genisligi. Ayni sebep gecenin basinda
graph'in es zamanli isle neden bu kadar yavasladigini da acikliyor.

### Bayraklar: OLU

`-K 2G` 6:56.35, `-2` 6:57.46, `--max-chain-skip 10 --max-chain-iter 1000` 6:57.33, taban
6:58.06. Hicbiri kipirdatmiyor.

### `.mmi`: kucuk ama bedava

6:33.74 vs 6:46.44 (FASTA), **cikti birebir ayni** (@PG disinda), RSS 1 GB dusuk. Cagri basina
~13 s, cift basina ~26 s.

### Kirpma: beklentimin tersi, ve asil bulgu bu

| | baz | duvar (-t 30) | Mbp/s |
|---|---|---|---|
| kirpilmamis | 4.551 Gbp | 17:08.85 | 4.42 |
| kirpma 1500 bp | 2.693 Gbp | 15:10.69 | 2.96 |

**Baz %41 dustu, sure yalnizca %11.** Ayni 255,660 okuma. Iki noktadan model:

```
sure = 2.89 ms/okuma + 63.6 ns/baz
okuma basina sabit is, toplam surenin %72'si
```

Yani maliyet baz basina degil **okuma basina**. "Daha kucuk okumalar" kaldiraci zayif;
**"daha az okuma"** kaldiraci gucllu. Superlineer-uzunluk hipotezim yanlisti.

### Buradan cikan is

`pick` 2,164 aday contig icin 150,374 capa cekirdegi kullaniyor ve 255,660 okuma seciyor:
**lokus basina 118 okuma**. Germline testinin derinlik egrisi olculmus:

```
 3x -> 290 yanlis cagri
 5.5x -> 172
11x -> 136
44x -> 116
```

118x, azalan getirinin cok otesinde. Lokus basina okuma sayisini kapamak dogrudan kazanc:

| lokus basina | okuma | mm2 normal (tahmin) | beklenen ek yanlis cagri |
|---|---|---|---|
| 118 (bugun) | 255,660 | 16:16 | taban |
| ~40 | ~87,000 | ~5:40 | ~0 (44x zaten "tam") |
| ~20 | ~43,000 | ~2:50 | ~20 (136 vs 116) |

~40'a kapamak boru hattini 52:48 -> ~42:00 indirir ve olculmus egriye gore precision'a
dokunmaz. Bunu yapmanin iki yolu var: `pick`'e contig basina kapak koymak, ya da **ters
indeksi kullanmak** -- N kayit/dugum tam olarak bu kapagin kendisi. Gecenin basinda
"N=8 pick icin az" demistim; dogru soru "118 zaten cok fazla mi" imis, ve cevabi evet.

**Bu, gecenin en actionable bulgusu ve olculmemis bir tahmin degil: iki nokta olculdu,
model onlardan cikti, derinlik egrisi zaten kayitliydi.**

## Hipotezin testi: daha az okuma (olculdu, hcc1954, -t 14, bos kutu)

`pick`'in `min_shared` esigi yukseltildi. Ayni contig'ler, ayni aday VCF, tek degisen esik.

| min_shared | okuma | Gbp | mm2 normal | RSS | P | R | F1 | germline elemesi |
|---|---|---|---|---|---|---|---|---|
| 5 (bugunku) | 255,660 | 4.551 | **16:16.16** | 18.60 GB | 0.894 | 0.862 | **0.878** | 871 |
| 10 | 203,485 | 3.645 | 12:28.91 | 18.61 GB | 0.893 | 0.862 | 0.877 | 870 |
| **20** | **145,841** | **2.625** | **7:28.41** | 18.49 GB | 0.893 | 0.862 | **0.877** | **870** |
| 40 | 86,855 | 1.583 | 4:50.83 | 18.60 GB | 0.887 | 0.862 | 0.874 | 862 |

**Sonuc: `min_shared` 5 -> 20, en buyuk asamayi 16:16 -> 7:28 (-%54) indiriyor ve dogruluga
dokunmuyor.** Germline suzgeci tam olarak ayni 870 cagriyi eliyor, yani testin kendisi hic
zayiflamiyor; F1 farki tek bir cagri (0.878 -> 0.877).

40'ta bozulma basliyor: eleme 870 -> 862, sekiz yanlis cagri hayatta kaliyor, precision
0.894 -> 0.887. Yani 20 ile 40 arasinda bir yerde derinlik egrisinin dizi bulunuyor, ve
20 guvenli tarafta.

Kirpma ile karsilastirma, ayni asamada:

| kol | baz | okuma | sure |
|---|---|---|---|
| taban | 4.551 Gbp | 255,660 | 16:16 |
| kirpma 1500 bp | 2.693 Gbp (-%41) | 255,660 (ayni) | 15:11 (-%11) |
| min_shared 20 | 2.625 Gbp (-%42) | 145,841 (-%43) | **7:28 (-%54)** |

Neredeyse ayni baz sayisi, cok farkli sure. Belirleyen sey **okuma sayisi**, tam da modelin
soyledigi gibi (2.89 ms/okuma + 63.6 ns/baz).

### Onerilen degisiklik

`main.cpp:run_pipeline`'da `pick` cagrisina `min_shared = 20` gecmek (bugun varsayilan 5
kullaniliyor, argüman hic verilmiyor). Tek satir. Beklenen etki, hcc1954 temiz kosuda:
mm2 normal 16:16 -> 7:28, boru hatti ~49:00 -> ~40:00, minimap2'nin payi %52.4 -> ~%40.

Bir sart: **bu tek ciftte olculdu.** Digerlerinde de dogrulanmali; ozellikle h1437 ve
hcc1937 daha az contig uretiyor, orada `min_shared=20` fazla eleyebilir.
