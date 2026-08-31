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
