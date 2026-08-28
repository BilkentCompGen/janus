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

## Sonuçlar

`results/exp1_sweep.tsv` ve `figs/exp1_core_sweep.{pdf,png}` — koşu bittikçe dolduruluyor.
