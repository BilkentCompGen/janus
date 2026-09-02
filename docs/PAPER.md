# Paper — mevcut hali

`paper/main.tex`, 2026-09-01 itibarıyla. Bu dosya paperın taslağını özetler: hangi bölüm ne
iddia ediyor, hangi şekil/tablo nereden geliyor, ne eksik.

Başlık: *Snape: Structural Variants from N Tumors by Alignment-free Pangenome Extraction*
Yazarlar: Muti Kara, Akmuhammet Ashyralyyev, Can Alkan (Bilkent, Bilgisayar Mühendisliği).

## Kaynak dosyalar

| Dosya | İçerik |
|---|---|
| `paper/main.tex` | preamble, yazarlar, bölüm sırası, Şekil 1'in float'ı |
| `paper/sections/abstract.tex` | özet |
| `paper/sections/intro.tex` | giriş, 5 paragraf |
| `paper/sections/overview.tex` | 2.1 System overview, 2 paragraf |
| `paper/sections/results.tex` | 2.2, 2.3, 2.4 + Tablo 1-3 + Şekil 2-4 |
| `paper/sections/methods.tex` | 3. bölüm, 8 paragraf |
| `paper/sections/future.tex` | yer tutucu, `main.tex`'te yorum satırında |
| `paper/mybib.bib`, `paper/akmami.bib` | kaynakça; stil `myrecomb.bst` |
| `paper/build.sh` | docker texlive ile derleme |

## Yapı

```
Abstract
1  Introduction
2  Results
   2.1  System overview                          -> Şekil 1
   2.2  Accuracy on real tumor--normal pairs     -> Tablo 1, Şekil 2
   2.3  Calling N tumors in one personal pangenome -> Şekil 3
   2.4  Runtime and memory per pair              -> Tablo 2, Tablo 3, Şekil 4
3  Methods
```

Results, Methods'tan önce geliyor. Her alt bölümün başında ve kaynakçadan önce `\FloatBarrier`
var; float oranları `main.tex`'te elle ayarlı (`topfraction` 0.92, `textfraction` 0.06).

## Ana iddialar

Girişte ve özette verilen, Results'ta sayıyla karşılanan iddialar:

1. Somatik karar hiçbir okuma hizalanmadan, okumaların kendi LCP core'ları üzerinden veriliyor.
2. Bir normal + N tümör tek bir kişiye özel pangenom grafında sayılıyor; normal sayımı sıfır olan
   core somatik, tümör sayımları da onu hangi tümörlerin taşıdığını söylüyor.
3. CASTLE panelinin dört PacBio HiFi çiftinde Severus, SAVANA ve nanomonsv'nin yayımlanmış
   çağrılarıyla karşılaştırılabilir sonuç, 30--100 kat daha az dizi hizalayarak.
4. Panelin konsensüs truth seti on bir yayımlanmış koşunun oyu; snape onlardan biri değil.
   snape on ikinci üye olarak eşit şartlarda eklenip konsensüs yeniden kurulduğunda snape on iki
   koşunun en yüksek F1'ini alıyor.
5. Yirmi tümörlük simüle vakada tek grafta çağırmak yirmi ikili koşudan 1.65 kat hızlı ve biraz
   daha doğru.
6. Uçtan uca, iki okuma setini referansa hizalama rotasından daha hızlı.

**İddia 6 tutarsız.** Özet "6 times faster", giriş "about 6 times faster", 2.4 ise "7.0 to 9.7
times faster ... median of 8.4" diyor. Aynı sayı olmalı.

## Şekil ve tablo envanteri

| Numara | Etiket | Dosya | Üreten | Bölüm |
|---|---|---|---|---|
| Şekil 1 | `fig:overview` | `figures/overview.tex` (TikZ) | elle; ağaç kısmı `figures/mkoverview_tree.py` | 2.1 |
| Tablo 1 | `tab:castle` | `results.tex` içinde gömülü | — | 2.2 |
| Şekil 2 | `fig:castle` | `figures/fig1.pdf` | `figures/mkfig1.py` | 2.2 |
| Şekil 3 | `fig:sim` | `figures/fig2.pdf` | `figures/mkfig2.py` | 2.3 |
| Tablo 2 | `tab:time` | `results.tex` içinde gömülü | `figures/mktable2.py` | 2.4 |
| Tablo 3 | `tab:map` | `results.tex` içinde gömülü | `figures/mktable3.py` | 2.4 |
| Şekil 4 | `fig:cost` | `figures/fig3.pdf` | `figures/mkfig3.py` | 2.4 |

PDF adları ile şekil numaraları kaymış: `fig1.pdf` Şekil 2, `fig2.pdf` Şekil 3, `fig3.pdf`
Şekil 4. Şekil 1'in PDF'i yok, TikZ olarak gömülü.

`figures/overview.pdf` ve `figures/mkoverview.py` eski matplotlib sürümüdür, TikZ'e geçilince
kullanım dışı kaldı. `paper/probe.tex`, `probe.sh`, `probe.pdf`, `main.tex.probe` Şekil 1
üzerinde çalışırken kullanılan geçici dosyalardır.

## Eksikler ve açık kararlar

1. **scDNA-seq varsayımı hiçbir yerde yazmıyor.** Yöntemin girdisinin ne olduğu — tek hücre DNA
   dizilemesi — ne özette, ne girişte, ne Methods'ta beyan ediliyor. Yazılacak.
2. **Kapsam titrasyonu yok.** CASTLE örnekleri yalnız tam kapsamda koşuldu (normal 56--71x,
   tümör 57--65x). 10x / 15x / 20x / 25x ... için doğruluk ve süre ölçümü yok. Okumaları belli
   bir yere kadar okuyarak ölçülebilir; yeni dosya üretmeye gerek yok.
3. **LCP core seviyesi sabit.** Bütün sonuçlar varsayılan `--level 5` ile alındı. Seviye 3, 4, 5,
   6 için ayrı doğruluk ve performans ölçümü yok. (Beklenti: 3 ve 4 ONT'de de çalışabilir, 5 ve 6
   çalışmaz.)
4. **`Supplementary Table~S1` yok.** 2.2'de "Scores against the published set are given in
   Supplementary Table~S1" diye atıf var, karşılığı üretilmedi. Ya tablo yazılacak ya cümle
   kalkacak.
5. **Şekil 4 ile Tablo 2/3 birbirini tekrarlıyor.** Şekil 4a Tablo 2'nin, Şekil 4c Tablo 3'ün
   aynı sayılarını çiziyor. Hangisinin kalacağı kararlaşmadı.
6. **`index.bin` Methods'ta tarif edilmiyor.** 2.4 disk bütçesinde "The index for each core, along
   with the reads it carries, is 8.5 GB" diyor; Methods'ta böyle bir yapı geçmiyor.
7. **Özet/giriş ile 2.4 arasındaki hız çelişkisi** (yukarıda, iddia 6).
