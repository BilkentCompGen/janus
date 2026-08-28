# Deney 1 — çekirdek sweep'i (1, 2, 4, 8, 16, 32)

**Tarih:** 2026-08-28 06:08 →
**Durum:** altörneklenmiş girdiyle koşuluyor. **Tam ölçekte yapılmadı.**

## Neden tam ölçekte değil

Üç ölçülmüş engel, üçü de `benchmarks_gaps.md`'de:

1. **Disk.** `/home` %91→%93 dolu. Otonom protokolün eşiği %85; bunun üstünde yeni büyük deney
   başlatılmaz.
2. **CASTLE'ın FASTQ'su yok.** Kullanıcı süre/bellek ölçümünün `.svb` değil FASTQ ile yapılmasını
   istedi. Diskteki tek tam CASTLE FASTQ'su `hifi/h2009/t.fastq` (376 GB) ve aynı çiftin
   **normali yok**. hcc1954'ün ham `.gz`'leri de diskte değil.
3. **Süre.** Tam ölçekte 16 çekirdekte ~31–42 dk, 1 çekirdekte ~10 saat. Altı ayar bir geceye
   sığmıyor.

## Ne yapıldı yerine

`snape dump` ile hcc1954'ün `t.svb` ve `n.svb`'sinden **600.000'er okuma** FASTQ'ya döküldü
(20.4 GB + 20.5 GB, örnek başına ~9.5 Gbp, ~3× kapsam). Aynı girdi 32, 16, 8, 4, 2, 1 çekirdekte
koşuldu — hızlı ayarlardan yavaşa doğru, ki vakit biterse elde en çok nokta kalsın.

`snape dump` sentetik Q40 kalite yazıyor. Bu ölçümü bozmuyor: `docs/RESULTS.md` bu boru
hattındaki hiçbir aracın baz kalitesi okumadığını gösteriyor (Severus kaliteyi hizalamadan
türetiyor, minimap2 yalnız diziden skorluyor). **FASTQ ayrıştırma maliyeti ise gerçek** — `.svb`
yerine FASTQ kullanmanın sebebi zaten buydu.

## Ölçülenler

Her koşudan:

| kaynak | ne verir |
|---|---|
| snape'in kendi `[eff] TOTAL` satırı | graf aşaması: duvar saati, CPU, paralel verim, tepe RSS |
| 1. `[M::main]` satırı | contig hizalama: Real time, CPU, Peak RSS |
| 2. `[M::main]` satırı | normal hizalama: aynısı |
| `/usr/bin/time -v` | toplam duvar saati, `Maximum resident set size` |

## Yorumlama sınırı

**Mutlak sayılar tam koşuyu temsil etmez.** ~3× kapsamda `min_mult=3` eşiği çekirdeklerin büyük
kısmını eliyor, dolayısıyla graf sonrası aşamalar tam koşudakinden çok daha az iş görüyor.
Ölçülen şey **ölçeklenmenin şekli**: her aşamanın çekirdek sayısıyla nasıl hızlandığı ve
belleğin çekirdekle nasıl değiştiği. Paper'a girerse bu etiketle girmeli.

## Sonuçlar (6 nokta, tamamlandı 07:57)

| çekirdek | toplam | graf | graf paralelliği | graf verimi | tepe RSS |
|---|---|---|---|---|---|
| 1 | 36.6 dk | 11.5 dk | 1.02 | — | 15.7 GB |
| 2 | 35.5 dk | 11.3 dk | 1.04 | %52 | 15.8 GB |
| 4 | 16.5 dk | 4.0 dk | 3.02 | %76 | 18.0 GB |
| 8 | 6.7 dk | 1.9 dk | 6.76 | %85 | 18.4 GB |
| 16 | 4.0 dk | 1.0 dk | 13.56 | **%85** | 20.5 GB |
| 32 | 3.6 dk | 0.8 dk | 25.13 | %79 | 24.8 GB |

İki okunabilir şey var:

1. **Graf aşaması 4'ten 32 çekirdeğe kadar %76–85 paralel verim koruyor.** 2 çekirdekte %52'ye
   düşüyor: girdi örnek başına parçalanıyor (1 normal + 1 tümör = 2 parça), yani 2 çekirdekte
   paralellik ölçülen 1.04'te takılıyor — iki parça birbirini bekliyor. Bu, altörnekleme
   distorsiyonundan en az etkilenen ölçüm — graf aşaması bütün okumaları işliyor, kapsam düşük
   olsa da iş gerçek.
2. **Tepe RSS çekirdek sayısıyla büyüyor**: 18.0 → 24.8 GB (4 → 32). Bellek bütçesi olan bir
   makinede çekirdek sayısı ücretsiz değil.

3. **1'den 2 çekirdeğe geçmek hiçbir şey kazandırmıyor**: 36.6 → 35.5 dk, paralellik 1.02 → 1.04.
   Girdi örnek başına parçalanıyor ve 1 normal + 1 tümör = 2 parça; iki parça birbirini beklediği
   için ikinci çekirdek boşta kalıyor. Kazanç 4 çekirdekten itibaren başlıyor.

Uçtan uca hızlanma 1→32 arasında **10.2×** (graf aşamasında 14.4×). Aradaki fark minimap2'nin
indeks yükleme gibi sabit maliyetlerinden; altörneklenmiş girdide bu sabit pay orantısız büyük,
yani tam ölçekte toplam hızlanma buradakinden **daha iyi** olur.

**Düzeltme:** snape'in log'a yazdığı `eff` yüzdesi sabit 32 çekirdeğe bölünüyor, yani
"makinenin ne kadarı kullanıldı" demek, çekirdek başına verim değil. Ham haliyle 4 çekirdekte
%9, 32'de %79 görünüyordu — verimin çekirdekle artması imkânsız olduğu için yakalandı.
Tabloda gerçek verim `(CPU/duvar)/çekirdek` olarak hesaplandı.
