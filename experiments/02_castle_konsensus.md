# Deney 2 — CASTLE'da Minda ve Minda_ws konsensüsleri

**Tarih:** 2026-08-28 05:28–05:33
**Çiftler:** hcc1954 ve hcc1937 (ikincisi bu gece koşuldu: 31:05, tepe RSS 35.16 GB, 607 çağrı)
**snape çağrıları:** `hifi/hcc1954_out2/calls.vcf`, 1.070 kayıt, refactor `75b3a05` sonrası
**Çıktı:** `results/exp2/`, `results/exp2_[AB]_scores.txt`, `figs/exp2_castle_consensus.{pdf,png}`

## İki koşul

Gerçek küme (`minda_4_2`, CASTLE'ın kendi dosyası) 11 çağırıcı koşusunun oylamasıyla kuruluyor:
bir kayıt ≥2 platform grubunda ve ≥4 koşuda destek alırsa gerçek sayılıyor. snape o oylamaya
katılmıyor, rakiplerin hepsi katılıyor.

| | gerçek küme | snape | ne ölçer |
|---|---|---|---|
| **A** | `minda_4_2` (11 çağırıcı) | dışarıda | snape tek dezavantajlı taraf |
| **C** | aynı küme, snape 12. üye olarak yeniden kuruldu | üye | asimetri kaldırılmış, herkes aynı şartlarda |

A ile C arasındaki fark bir yöntem farkı değil; gerçek kümenin nasıl kurulduğunun sonucu.

## Sonuçlar

**hcc1954**

| koşul | snape F1 | snape sırası | 1. sıra |
|---|---|---|---|
| A, snape dışarıda | 0.8704 | 3 / 12 | 0.8865 |
| C, snape üye | **0.8851** | **1 / 12** | snape (2. sıra 0.8802) |

**hcc1937** — snape'in en zayıf paneli, ve sonuç hcc1954 ile **aynı yönde değil**

| koşul | snape F1 | snape sırası | 1. sıra |
|---|---|---|---|
| A, snape dışarıda | 0.7649 | 4 / 12 | 0.8257 |
| C, snape üye | 0.7819 | 4 / 12 | 0.8371 |

İki çift aynı yönde değil: hcc1954'te snape üye olduğunda 1. sıraya çıkıyor, hcc1937'de her iki
koşulda da 4. sırada kalıyor. Bu çiftte snape'in F1'i refactor öncesiyle aynı kalmış
(`docs/RESULTS.md` 0.765 ↔ bu ölçüm 0.7649), çağrı sayısı düşmüş olsa da.

## Karar

**Kabul.** snape oylamaya katılmadığında katıldığındakinden düşük skorluyor: hcc1954'te
0.8704 → 0.8851 (**+0.0147**), hcc1937'de 0.7649 → 0.7819 (**+0.0170**). Bu fark bir yöntem
farkı değil, gerçek kümenin nasıl kurulduğunun sonucu — ve A'da rakiplerin **hepsi** bu avantajı
taşıyor, snape taşımıyor.

hcc1954'te snape üye olduğunda panelin en yükseği; hcc1937'de her iki koşulda da 4. sırada.
"Comparable accuracy" iddiası iki çiftin ikisiyle birlikte verilmeli.

## Bilinen boşluklar

- **Yalnız hcc1954.** h2009'un tümör `.svb`'si diskte yok; hcc1937 için snape bu gece
  koşuldu (bkz. `experiments/02b_*`), h1437 koşulmadı.
- Minda'nın `decompose.py`'sinde sayısal VCF ID'lerinde çöken bir hata vardı; yamalandı
  (`file-removals.txt`'e kayıtlı).
- `minda_4_2_ws` koşulları (eski B ve D) rapordan çıkarıldı. Ham çıktıları `results/exp2*/`
  altında duruyor, silinmedi, ama hiçbir yere girmiyor.
