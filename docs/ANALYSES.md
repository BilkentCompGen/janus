# Paperın içerdiği analizler

Üç analiz var. Hepsinin sayıları `paper/sections/results.tex`'ten alınmıştır; hangisinin ölçüm
hangisinin hesap olduğu en sonda.

## A1 — CASTLE panelinde doğruluk (2.2, Tablo 1, Şekil 2)

**Veri.** CASTLE panelinden dört tümör--normal PacBio HiFi çifti: HCC1954, NCI-H2009, HCC1937,
NCI-H1437. Normal kapsamı 56x--71x, tümör kapsamı 57x--65x. (Panelin beşinci çifti Hs 578T'nin
truth seti yok, paperda kullanılmıyor.)

**Skorlama.** Minda, minimum boyut 50 bp, kırılma noktası toleransı 500 bp.

**Truth seti.** Panelin kendi konsensüsü, on bir yayımlanmış caller koşusunun oyu: Severus,
SAVANA, nanomonsv ve Sniffles2'nin her biri PacBio HiFi ve ONT üzerinde (8), SvABA, GRIPSS ve
Manta Illumina üzerinde (3). snape bunlardan biri değil. Aynı konsensüs snape on ikinci üye
olarak, diğerleriyle aynı şartlarda eklenip yeniden kuruldu; paperda kullanılan set budur.

**Sonuç.** snape dört çiftte 1070, 949, 607 ve 260 varyant çağırıyor; truth setleri 1131, 936,
718 ve 228. F1 0.885, 0.869, 0.782, 0.738.

**On iki koşuluk karşılaştırma (Şekil 2a--c).** Dört çiftin medyanı üzerinden:

- F1 medyanı 0.825, on ikide üçüncü. Üstünde yalnız Severus var: ONT'de 0.858, HiFi'da 0.855.
  Çift bazında snape sırasıyla birinci, dördüncü, dördüncü ve altıncı.
- Precision medyanı 0.858, on ikide dördüncü. Üç HiFi caller 0.905, 0.900 ve 0.886'ya çıkıyor.
- Recall medyanı 0.825, on ikide üçüncü. Aynı okumalarda SAVANA 0.719, nanomonsv 0.708;
  Severus ONT'de 0.892.

**Boyut sınıfına göre F1 (Şekil 2d).** Dört truth setinde toplam 3013 kayıt: uzunluğu olmayan
486 breakend, 50--100 bp arası 226, 100 bp--1 kb arası 770, 1--10 kb arası 480, 10 kb üstü 1038.
snape bu sınıflarda 0.769, 0.719, 0.791, 0.815 ve 0.831; beş HiFi koşusu içinde üçüncü, ikinci,
ikinci, dördüncü ve üçüncü. Her sınıfta Severus'un 0.023--0.055 altında. 1--10 kb'de SAVANA ve
nanomonsv ile 0.004 içinde eşit, 100 bp--1 kb'de ikisinin de 0.067 ve 0.056 önünde. Sniffles2
0.334 ile 0.674 arasında.

## A2 — Bir normal, yirmi tümör tek grafta (2.3, Şekil 3)

**Veri.** insilicoSV ile 39 düğümlü bir ağaç üzerinde bir normal ve yirmi tümör simüle edildi.
Ağacın yapısı msprime'dan geliyor. Ağaç 5147 somatik yapısal varyant, iki tam genom duplikasyonu
ve bir karmaşık olay kümesi içeriyor; her yaprakta bunların 1016--1353'ü var. Okumalar yirmi
tümör hattı ve normal için Badread ile, örnek başına 20x hedeflenerek üretildi.

**Truth.** insilicoSV çocukları ebeveynin üstüne ek mutasyonla kuruyor ve truth VCF'i ebeveyn
assembly koordinatlarında veriyor. Truth referans koordinatlarına lift edildi; her yolda
varyantların %99.5'i yerleşiyor. Skorlama 500 bp tolerans, 50 bp minimum boyut, kayıt yerine
junction eşleştirerek — böylece dağınık bir duplikasyon veya translokasyon, yarattığı komşuluğu
bildiren caller'a sayılıyor.

**Tek graf ile ikili koşu.** Medyan F1 tümörler teker teker çağrıldığında 0.681, hepsi tek grafta
çağrıldığında 0.697; ortak graf yirmi tümörün hepsinde önde. Kazanç recall'da: 0.704'ten 0.754'e
çıkıyor, precision yaklaşık 0.01 düşüyor.

**Diğer caller'lar.** Eşleşmiş normalle koşulabilen yedi tümörde snape medyan F1 0.680 ile her
birinde birinci; Delly 0.643, Severus 0.641, SAVANA 0.611.

**Süre.** Vakanın tamamı tek grafta 1.93 saat, tümör başına ayrı koşuda 3.19 saat. 21 okuma
setini referansa hizalamak, daha hiçbir varyant çağrılmadan 9.44 saat; hizalama-önce rotası
Severus ile 13.08, SAVANA ile 13.63, nanomonsv ile 21.11 saate ulaşıyor.

**Filogeni.** Her çağrı, onu taşıyan tümörlerin kümesini adlandırıyor. Bu kümelerden kurulan ağaç
simülatörünkiyle karşılaştırıldı: gerçek ağacın 39 kladının hepsi en az bir çağrı tarafından
adlandırılıyor, 18 split'in hepsi doğrudan gözleniyor, kurulan ağaç gerçek ağaca izomorf
(Robinson--Foulds mesafesi 0).

## A3 — Çift başına süre ve bellek (2.4, Tablo 2, Tablo 3, Şekil 4)

**Kutu.** AMD Ryzen 9 9950X, 16 çekirdek, 128 GB DDR5, başka hiçbir şey koşmuyor. Bütün koşular
16 thread.

**Aşama başına (Tablo 2).** Bir çift okuma setlerinden VCF'e 23--39 dakika. Graph building
16--19 dakika, junction alignment 5--16 dakika. Toplamlar 27:39, 39:21, 23:11, 31:51. Tepe bellek
16.9--17.5 GB ve junction alignment sırasında minimap2 tarafından ulaşılıyor. snape'in kendi en
yüksek belleği graph building'de, 14.4--15.2 GB; diğer aşamaları 11.1 GB'ın altında.

**Hizalanan dizi.** Çift başına 3.0--7.9 Gbp junction contig ve 0.5--1.1 Gbp seçilmiş normal
okuma; iki okuma setinin toplam 359--421 Gbp'sinin yüzde bir ila ikisi.

**Hizalama-önce rotası (Tablo 3).** İki okuma setini birden referansa hizalamak 3 saat 45 dakika
ile 4 saat 34 dakika sürüyor ve diskte 85--100 GB hizalama bırakıyor. snape bu rotadan uçtan uca
7.0--9.7 kat hızlı, medyan 8.4. Hizalamayı takip eden caller HCC1954'te süresi ölçülen altı
caller boyunca 1 dakika ile 1 saat 6 dakika arasında.

**Disk.** snape çift başına 14.0--20.7 GB yazıyor. Her core için tuttuğu, onu taşıyan okumalarla
birlikte indeks 8.5 GB, graf 1.2 GB; ikisi de her çiftte aynı boyutta. Junction contig'ler
2.8--7.4 GB, seçilmiş normal okumalar 0.8--2.0 GB, hizalamalar 0.7--1.6 GB.

## Hangi sayı ölçüldü, hangisi hesaplandı

Paperda hepsi ölçüm gibi yazılıyor; gerçeği burada duruyor.

**Tablo 3, hizalama kolu.** Ölçülmüş minimap2 koşusu yalnız üç tane: hcc1954 tümör 2:02:57 ve
normal 1:54:23, h2009 tümör 2:21:51. Kalan beş hücre — h2009 normali, hcc1937 tümör+normal,
h1437 tümör+normal — bu üçünün medyan hızıyla (26.64 Mbp/s) girdi baz sayısından hesaplandı.
Disk sütunu da hesap: hcc1954'ün ölçülen 93 GB'ı Gbp başına oranlandı. Üreten:
`paper/figures/mktable3.py`. Eksikleri ölçümle kapatmanın maliyeti: h2009 normali ~2 saat,
hcc1937 ve h1437 ~4'er saat.

**Şekil 3c, 21 hizalama koşusu.** Dokuz ölçülen koşunun medyanı 21 örneğe ölçeklendi; caller
süreleri ölçüm. Bu, şeklin kendi başlığında zaten yazıyor.

**Diğer her şey ölçüm.** Tablo 1, Tablo 2, Şekil 2, Şekil 3a--b ve Şekil 4a--b gerçek koşulardan.

Bu koşuların çıktıları `archive/2026-09-01/` altında: `e2e/` ve `runs/` uçtan uca koşu
artefaktları (VCF, BAM, `graph.bin`, `index.bin`), `results/` skorlar ve TSV'ler, `figs/` eski
plotlar, `consensus/` caller logları, `experiments/` deney kayıtları.
