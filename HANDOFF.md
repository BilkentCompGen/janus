# HANDOFF — 2026-08-28 gecesi

**Oturum:** 03:57 → (devam ediyor). Amaç: Results bölümünün beş deneyi için plot hazırlamak.

## Tek cümlede

Deney 2, 3, 4 ve 5 **bitti ve plot'landı**; deney 1 (çekirdek sweep'i) tam ölçekte yapılamadı,
altörneklenmiş hali koşuluyor. Ayrıca depo `.git`'i 141 GB'dan 7.6 MB'a indirildi.

## Üretilen figürler (`figs/`)

| dosya | deney | ne gösteriyor |
|---|---|---|
| `exp2_castle_consensus` | 2 | CASTLE hcc1954, dört konsensüs konfigürasyonu, snape'in sırası |
| `exp2_castle_consensus_hcc1937` | 2 | aynısı hcc1937 için |
| `exp3_phylogeny` | 3 | ağaç kurma yöntemleri (RF), bootstrap desteği, klad oranı |
| `exp4_total_time` | 4 | 1+20 vakanın toplam süresi, üç kol |
| `exp5_joint_vs_pair` | 5a | 20 hatta ortak graf ↔ ikili kayıt farkı |
| `exp5b_consensus` | 5b | 7 yaprakta snape konsensüs içinde/dışında |

Her figürün altında kaynak notu var: hangi dosyadan, hangi yöntemle, ölçüldü mü ölçeklendi mi.

## Başlıca sayılar

**Deney 2 — CASTLE.** `minda_4_2_ws` CASTLE'ın kendi dosyası ve **"without Severus"** demek
(iki `minda.log`'un `cmd:` satırından doğrulandı).

| konfig | hcc1954 | hcc1937 |
|---|---|---|
| A `minda_4_2`, snape dışarıda | 3/12, F1 0.8704 | 4/12, F1 0.7649 |
| B `minda_4_2_ws`, snape dışarıda | **1/12**, 0.8616 | 4/12, 0.7740 |
| C `minda_4_2` + snape üye | **1/12**, 0.8851 | 4/12, 0.7819 |
| D `_ws` + snape üye | **1/10**, 0.8874 | 2/10, 0.8061 |

hcc1954'te Severus'un oyları gerçek kümeden çıkınca snape 3.'lükten 1.'liğe geçiyor.
**hcc1937'de bu olmuyor — snape orada 4. sırada.** İki çift aynı yönde değil.

**Deney 3 — filogeni.** `greedy-split` gerçek topolojiyi **birebir** kuruyor (RF = 0, kapsam ve
kesinlik 1.000). Klad oranı %88.0 (null %70.7); bilgilendirici örüntülerde %59.6 (null **%1.0**).
39/39 klad tanıklandı. 17 kladın 13'ü %100 bootstrap.

**Deney 4 — süre.** (C kolu **ölçeklendi**, 21 genom gerçekten hizalanmadı.) Ortak graf 21 örnek **1.93 sa**; ikili ×20 3.19 sa (1.65×);
21×minimap2 + 20×severus 13.08 sa (**6.76×**); +sniffles2 9.83 sa (5.08×); +nanomonsv 21.11 sa.

**Deney 5.** Ortak graf 20 hatta ikiliden %10.7 daha çok çağrı (+2.641), kaybettiği 367.
Konsensüste snape dışarıda 0.773–0.820 (3/8), içeride 0.900–0.928 (6/7 yaprakta 1/8).

## Sabah yapılması gerekenler

1. **`git push --force origin changes`** — yerel `.git` temizlendi (141 GB → 7.6 MB) ama uzak hâlâ
   eski, 16 GB BAM'lı tarihçeyle. Force-push classifier tarafından engellendi (doğru bir engel:
   yayımlanmış tarihçeyi sen yokken ezmek geri alınamaz).
2. Deney 5b'nin "içeride" sayıları **tek başına kullanılmamalı**: üyelik primi orada ~0.12 F1,
   CASTLE'da ~0.015. Sebep konsensüsün küçük olması (7–8 çağırıcı, ≥3 oy).
3. hcc1937'nin sonucu hcc1954'ünkiyle aynı yönde değil. Paper'da "comparable" iddiası bu iki
   çiftin ikisini de göstermeli.

## Ne denendi ve bırakıldı

- **Deney 1 tam ölçekte**: disk %91 dolu (protokol eşiği %85), CASTLE'ın normal FASTQ'su yok,
  6 çekirdek ayarında tam koşu ~10+ saat. Altörneklenmiş sürüm koşuluyor (600k okuma/örnek,
  ~3× kapsam) ve **öyle etiketlenecek** — mutlak RAM'i tam ölçekte vermez.
- **h2009**: tümör `.svb`'si diskte yok (`t.svb.tmp` yarım kalmış). **h1437**: koşulabilirdi,
  zaman yetmedi.
- **delly** deney 4'ten çıkarıldı: koşuları zaman boşluğuyla ayrılamıyor, ölçüm güvenilmez.

Kalanların tam listesi: `benchmarks_gaps.md`.

## Düzeltilen hatalar

| ne | nerede | yedek |
|---|---|---|
| `eval_phylo.py` değişken gölgeleme (`got` int'e dönüyor, JSON dökümü patlıyor) | `tmp/overnight/eval_phylo.py:296` | `tmp/night/eval_phylo.py.orig` |
| Minda sayısal VCF ID'lerinde çöküyor | `tmp/.venv/.../minda/decompose.py:357` | `tmp/night/minda_decompose.py.orig` |
| Arşivdeki HCC1937 `minda_4_2` VCF'inde `#` önekli olmayan `cmd:` satırı | temiz kopya | `tmp/night/clean/` |
| İlk `true_burden` yazımım yanlıştı (T01 için 391, doğrusu 1277) | `tmp/night/exp5_joint_vs_pair.py` | — |

Silinen/bozulan her şey `file-removals.txt`'te; zaman, gerekçe, komut.

## Nerede ne var

```
PLAN.md               talimat + yorumu
activity_log.md       zaman damgalı olay günlüğü
file-removals.txt     silme/bozma kaydı
benchmarks_gaps.md    ölçülmeyenler ve neden
experiments/          deney başına kayıt (02, 03, 04, 05)
figs/                 PDF + PNG
results/              plot'ları üreten ham TSV/JSON
tmp/night/            bu gecenin betikleri, yedekleri, ara çıktıları
```
