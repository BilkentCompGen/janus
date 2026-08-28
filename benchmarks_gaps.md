# Ölçülmeyenler ve neden

Bu gecenin işi beş deneyin plot'larını hazırlamaktı. Aşağıdakiler **kasıtlı olarak** dışarıda
kaldı; sessizce atlanan bir şey yok.

## Deney 1 — çekirdek sweep'i (1, 2, 4, 8, 16, 32)

**Yapılmadı.** Kullanıcı zaten en sona koymuştu ("sadece ve sadece vaktin kalırsa"). Üç ölçülmüş
engel:

1. **Disk.** `/home` %91 dolu, 327 GB boş. Otonom protokolün eşiği %85; bu eşiğin üstünde yeni
   büyük deney başlatılmaz.
2. **FASTQ yok.** Kullanıcı süre/bellek ölçümünün `.svb` değil FASTQ ile yapılmasını istedi.
   Diskteki tek CASTLE FASTQ'su `hifi/h2009/t.fastq` (376 GB); aynı çiftin **normali yok**
   (`n.svb` var, FASTQ'su yok, açmak ~200 GB daha ister). hcc1954'ün ham `.gz`'leri de diskte değil.
3. **Süre.** Altı çekirdek ayarında tam koşu: 16 çekirdekte ~40 dk, 1 çekirdekte ~10 saat.
   Toplam bir geceye sığmıyor.

Yapılabilir olan ve kalan zamanda denenen: **altörneklenmiş** bir sweep — paralelleşmenin
*şeklini* (hızlanma ↔ çekirdek) verir, mutlak RAM'i tam ölçekte vermez. Öyle etiketlenecek.

## Deney 2 — CASTLE'ın dört çifti

**Yalnız hcc1954 tamamlandı**, hcc1937 bu gece koşuldu. İki çift dışarıda:

- **h2009**: tümör `.svb`'si diskte **yok** (`t.svb.tmp` 45 GB, yarım kalmış paketleme).
  `t.fastq` 376 GB duruyor, paketleme yeniden koşulabilir ama bu gece disk elvermiyor.
- **h1437**: `t.svb` ve `n.svb` var, koşulabilirdi; zaman yetmedi.

## Deney 4 — snape'in 21 örnekli koşusunun tepe RSS'i

**Ölçülmedi.** Elde 1:1 için 32.8 GB var (`docs/RESULTS.md`, hcc1954). 21 örnek için plan
belgesinde ~75 GB tahmini vardı ama bu bir tahmin; tabloda "ölçülmedi" yazıyor, tahmin
yazılmadı.

## Deney 4 — delly

Dışarıda bırakıldı. Koşuları 30 dakikalık boşluk eşiğiyle ayrılamıyor (4 "koşu", medyan 380 dk),
yani ölçüm güvenilir değil. Diğer beş çağırıcı raporlandı.

## Deney 5 — koordinat tabanlı TP/FP

Simülasyonda **hazır truth VCF yok**: `mutations.tsv` mutasyonları ebeveyn genom koordinatlarında
tutuyor, referansta değil. `eval_phylo.py`'nin klade üzerinden çalışmasının sebebi de bu. Bu
yüzden 5a'da TP/FP değil **kayıt farkı** raporlandı; TP/FP tarafı 5b'de Minda konsensüsüyle
verildi.

## Colored de Bruijn karşılaştırması

**Hiç koşulmadı.** Paper'ın intro'su Cortex'in üç sınırını sayıyor ama hiçbir dBG çağırıcısı
çalıştırılmadı. Bu, introda ampirik bir iddia gibi okunursa hakem oradan girer — ya P5'in tonu
tasarım tespitine çekilmeli ya da gerçekten bir dBG koşulmalı. Bu gecenin kapsamında değildi.

## Çoklu tohum / güven aralığı

Bütün sayılar tek koşudan. `docs/RESULTS.md` snape'in koşular arası oynaklığını ölçmüş
(çağrı kümesi dört panelde de birebir aynı), ama bu gecenin konsensüs koşuları için ayrı bir
tekrar yapılmadı.
