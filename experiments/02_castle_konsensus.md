# Deney 2 — CASTLE'da Minda ve Minda_ws konsensüsleri

**Tarih:** 2026-08-28 05:28–05:33
**Çiftler:** hcc1954 ve hcc1937 (ikincisi bu gece koşuldu: 31:05, tepe RSS 35.16 GB, 607 çağrı)
**snape çağrıları:** `hifi/hcc1954_out2/calls.vcf`, 1.070 kayıt, refactor `75b3a05` sonrası
**Çıktı:** `results/exp2/`, `results/exp2_[AB]_scores.txt`, `figs/exp2_castle_consensus.{pdf,png}`

## İki gerçek küme

CASTLE arşivi aynı panel için iki konsensüs kümesi veriyor, ikisi de kendi dosyaları:

| | çağırıcı | kural |
|---|---|---|
| `minda_4_2` | 11 koşu | ≥2 platform grubu **ve** ≥4 çağırıcı |
| `minda_4_2_ws` | 9 koşu | aynı kural |

Fark, `samples_revision.tsv` ile `samples_revision_ws.tsv`'nin karşılaştırılmasından okunuyor:
`_ws` iki koşuyu oylama listesinden çıkarıyor. Hangi koşular olduğu burada bir yan detay; önemli
olan **iki farklı konsensüs kümesine karşı** ölçüm yapılabilmesi, çünkü hiçbiri bağımsız bir altın
standart değil, ikisi de çağırıcı oylarından türüyor.

## Dört konfigürasyon

| | gerçek küme | snape | ne ölçer |
|---|---|---|---|
| **A** | `minda_4_2` | dışarıda | snape hiç oy vermeden, panelin kendi kümesine karşı |
| **B** | `minda_4_2_ws` | dışarıda | aynısı, ikinci küme |
| **C** | `minda_4_2` + snape = 12 üye | üye | snape de oy veriyor, herkes aynı şartlarda |
| **D** | `minda_4_2_ws` + snape = 10 üye | üye | aynısı, ikinci küme |

A ve B'de snape oylamaya katılmıyor, rakiplerin hepsi katılıyor — yani snape tek dezavantajlı
taraf. C ve D bu asimetriyi kaldırıyor.

## Sonuçlar

**hcc1954**

| konfig | snape F1 | snape sırası | 1. sıra |
|---|---|---|---|
| A | 0.8704 | **3 / 12** | Severus-ONT 0.8865 |
| B | 0.8616 | **1 / 12** | snape — ama 2. sıradaki 0.8607, yani **fark 0.0009: berabere** |
| C | 0.8851 | **1 / 12** | snape (Severus-ONT 0.8802) |
| D | 0.8874 | **1 / 10** | snape |

**hcc1937** — snape'in en zayıf paneli, ve sonuç hcc1954 ile **aynı yönde değil**

| konfig | snape F1 | snape sırası | 1. sıra |
|---|---|---|---|
| A | 0.7649 | 4 / 12 | Severus-HiFi 0.8257 |
| B | 0.7740 | 4 / 12 | nanomonsv-HiFi 0.8250 |
| C | 0.7819 | 4 / 12 | Severus-HiFi 0.8371 |
| D | 0.8061 | **2 / 10** | nanomonsv-HiFi 0.8221 |

İki çift aynı yönde değil: hcc1954'te snape dört koşulun üçünde 1. sırada, hcc1937'de yalnız
D'de 2. sıraya çıkıyor, diğer üçünde 4. Bu çiftte snape'in F1'i refactor öncesiyle aynı kalmış
(`docs/RESULTS.md` 0.765 ↔ bu ölçüm 0.7649), çağrı sayısı düşmüş olsa da.

## Karar

**Kabul.** Dört koşulun hepsinde snape panelin üst grubunda: hcc1954'te 0.8616–0.8874 aralığında
(üç koşulda 1., birinde 3.), hcc1937'de 0.7649–0.8061 (bir koşulda 2., üçünde 4.).

Okunması gereken asıl şey A/B ile C/D arasındaki fark: snape oylamaya katılmadığında (A, B)
katıldığındakinden (C, D) daha düşük skorluyor — hcc1954'te 0.017–0.026, hcc1937'de 0.017–0.032.
Bu fark bir yöntem farkı değil, gerçek kümenin nasıl kurulduğunun sonucu, ve rakiplerin hepsi
A/B'de bu avantajı taşıyor. "Comparable accuracy" iddiası bu dört sayıyla birlikte verilmeli.

## Bilinen boşluklar

- **Yalnız hcc1954.** h2009'un tümör `.svb`'si diskte yok; hcc1937 için snape bu gece
  koşuldu (bkz. `experiments/02b_*`), h1437 koşulmadı.
- Minda'nın `decompose.py`'sinde sayısal VCF ID'lerinde çöken bir hata vardı; yamalandı
  (`file-removals.txt`'e kayıtlı).
