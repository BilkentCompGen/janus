# Deney 4 — 1 normal + 20 tümörlü bir vakanın toplam maliyeti

**Tarih:** 2026-08-28 04:20–04:36
**Çıktı:** `results/exp4_runs.tsv`, `results/exp4_summary.{tsv,json}`, `figs/exp4_total_time.{pdf,png}`

## Üç kol

| kol | ne | kaynak |
|---|---|---|
| **A** | snape, tek ortak graf, 21 okuma kümesi birlikte | **ölçüldü** — deney 14 |
| **B** | snape, ikili, aynı iş 20 kez (4a) | **ölçüldü** — deney 14 |
| **C** | hizalama-önce: 21 × minimap2 + 20 × çağırıcı (4b) | minimap2 **ölçeklendi**, çağırıcı **ölçüldü** |

Kullanıcının talimatı gereği 21 genom **gerçekten hizalanmadı**. `tmp/overnight/resources.tsv`'de
`map-hifi` ile ölçülmüş 9 simülasyon örneği var (germline, T02–T08, T04 iki kez); medyanı 21'e
çarpıldı.

## Ölçüm hattı

`resources.tsv` her **prosesi** ayrı satır yazıyor — savana tek koşuda onlarca alt-proses açıyor,
o yüzden ham satır sayısı koşu sayısı değil. `tmp/night/agg_resources.py` satırları (araç, örnek)
başına koşulara indiriyor: duvar saati = max(bitiş) − min(başlangıç), tepe RSS = max(VmHWM),
CPU = Σ. İki koşuyu ayıran eşik 1800 s boşluk.

`severus` ve `cuteSV` satırları `multiprocessing` forkserver çocukları ve komut satırları örnek
adını **taşımıyor**. Her yaprak bir kez koşulduğu için koşu sürelerinin **dağılımı** yine de yaprak
başına çağırıcı maliyetidir; hangi yaprak olduğunu bilmek gerekmiyor.

## Sonuçlar

minimap2, örnek başına medyan **26.97 dk** (9 ölçüm), tepe RSS 26.1 GB.

| kol | toplam | ortak grafa göre |
|---|---|---|
| **A snape 1:20, ortak graf** | **1.93 sa** (116.1 dk) | — |
| B snape 1:1 × 20 | 3.19 sa (191.3 dk) | 1.65× yavaş |
| C 21×mm2 + 20×sniffles2 | 9.83 sa | 5.08× |
| C 21×mm2 + 20×cuteSV | 10.03 sa | 5.19× |
| C 21×mm2 + 20×severus | 13.08 sa | 6.76× |
| C 21×mm2 + 20×savana | 13.63 sa | 7.05× |
| C 21×mm2 + 20×nanomonsv | 21.11 sa | 10.91× |

Çağırıcı aşaması medyanları (koşu dağılımından): sniffles2 1.18 dk, cuteSV 1.82, severus 10.92,
savana 12.58, nanomonsv 35.02.

## Karar

**Kabul.** İki tez de destekleniyor: ortak graf ikili düzenden **1.65×**, hizalama-önce kolundan
**5.1–10.9×** hızlı.

## Bilinen boşluklar

- **snape'in 21 örnekli koşusunun tepe RSS'i ölçülmedi.** Elimizde 1:1 için 32.8 GB var
  (`docs/RESULTS.md`, hcc1954). 21 örnek için tahmin ~75 GB idi ama **ölçülmedi**, o yüzden
  tabloda "ölçülmedi" yazıyor.
- `delly-lr` dışarıda bırakıldı: 30 dk boşluk eşiği onun koşularını birleştirdiği için ölçümü
  güvenilir değil (4 "koşu", medyan 380 dk).
- Ölçümlerin bir kısmı makine yüklüyken alındı; duvar saatleri bu yüzden yer yer şişkin
  (delly T02 177 dk ↔ T04 6.5 dk). Bu, C kolunu **iyimser değil kötümser** tarafa itiyor değil —
  her iki kol da aynı makinede ölçüldü, ama koşular eşzamanlı değildi.
