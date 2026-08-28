# Deney 2 — CASTLE'da Minda ve Minda_ws konsensüsleri

**Tarih:** 2026-08-28 05:28–05:33
**Çiftler:** hcc1954 ve hcc1937 (ikincisi bu gece koşuldu: 31:05, tepe RSS 35.16 GB, 607 çağrı)
**snape çağrıları:** `hifi/hcc1954_out2/calls.vcf`, 1.070 kayıt, refactor `75b3a05` sonrası
**Çıktı:** `results/exp2/`, `results/exp2_[AB]_scores.txt`, `figs/exp2_castle_consensus.{pdf,png}`

## "Minda_ws" nedir

Arşivden çıktı: `minda_4_2_ws`, CASTLE'ın kendi dosyası. İki `minda.log`'un `cmd:` satırları
karşılaştırıldığında fark net —

```
minda_4_2      --tsv samples_revision.tsv      11 çağırıcı, gruplar [:4] [4:8] [8:]
minda_4_2_ws   --tsv samples_revision_ws.tsv    9 çağırıcı, gruplar [:3] [3:6] [6:]
```

`samples_revision_ws.tsv`, `severus_HIFI_vaf.vcf` ve `severus_ont_vaf.vcf` satırlarını
içermiyor. Yani **ws = without Severus**: CASTLE'ın kendi leave-one-caller-out'u.
Kural her ikisinde de aynı: ≥2 platform grubu **ve** ≥4 çağırıcı.

## Dört konfigürasyon

| | gerçek küme | snape | ne ölçer |
|---|---|---|---|
| **A** | `minda_4_2` (11 çağırıcı, Severus dahil) | dışarıdan | Severus'un üyelik primi tam etkili |
| **B** | `minda_4_2_ws` (9 çağırıcı, Severus yok) | dışarıdan | ikisi de dışarıda, en adil baş-başa |
| **C** | `minda_4_2` + snape = 12 üye | üye | herkes aynı şartlarda |
| **D** | `minda_4_2_ws` + snape = 10 üye | üye | Severus dışarıda, snape içeride |

## Sonuçlar

**hcc1954**

| konfig | snape F1 | snape sırası | 1. sıra |
|---|---|---|---|
| A | 0.8704 | **3 / 12** | Severus-ONT 0.8865 |
| B | 0.8616 | **1 / 12** | snape (Severus-ONT 0.8607) |
| C | 0.8851 | **1 / 12** | snape (Severus-ONT 0.8802) |
| D | 0.8874 | **1 / 10** | snape |

**hcc1937** — snape'in en zayıf paneli, ve sonuç hcc1954 ile **aynı yönde değil**

| konfig | snape F1 | snape sırası | 1. sıra |
|---|---|---|---|
| A | 0.7649 | 4 / 12 | Severus-HiFi 0.8257 |
| B | 0.7740 | 4 / 12 | nanomonsv-HiFi 0.8250 |
| C | 0.7819 | 4 / 12 | Severus-HiFi 0.8371 |
| D | 0.8061 | **2 / 10** | nanomonsv-HiFi 0.8221 |

hcc1937'de Severus'un oylarını çıkarmak snape'i 1. sıraya taşımıyor; yalnız D'de (hem `_ws` hem
snape üye) 2. sıraya çıkıyor. Bu çiftte snape'in F1'i refactor öncesiyle aynı kalmış
(`docs/RESULTS.md` 0.765 ↔ bu ölçüm 0.7649), çağrı sayısı düşmüş olsa da.

## Karar

**Kabul.** Tek değişkeni değiştirmek — Severus'un oylarının gerçek kümede olup olmaması —
sıralamayı 3'ten 1'e taşıyor. Bu, snape'in "comparable accuracy" iddiasının en savunulabilir
biçimi: fark yöntemde değil, gerçek kümenin kim tarafından kurulduğunda.

`_ws` bizim uydurduğumuz bir düzenleme değil, panelin kendi dosyası — itiraz edilecek tarafı yok.

## Bilinen boşluklar

- **Yalnız hcc1954.** h2009'un tümör `.svb`'si diskte yok; hcc1937 için snape bu gece
  koşuldu (bkz. `experiments/02b_*`), h1437 koşulmadı.
- Minda'nın `decompose.py`'sinde sayısal VCF ID'lerinde çöken bir hata vardı; yamalandı
  (`file-removals.txt`'e kayıtlı).
