# Deney 3 — Filogeni rekonstrüksiyonu

**Tarih:** 2026-08-28 04:16–04:36
**Girdi:** `simulation/phylo/snape_out/calls.vcf` — 21 örnekli tek ortak graf koşusu, 6.425 çağrı, 20 hat sütunu
**Gerçek:** `simulation/phylo/genomes/tumor/plan.json` (39 düğümlü ağaç), `mutations.tsv` (66.070 mutasyon, `birth_node` sütunuyla)

## Gerekçe

Ortak grafın verdiği şey koordinat değil, **destek örüntüsü**: bir çekirdeğin hangi tümör
hatlarında görüldüğü. Simülatör her mutasyonu doğduğu düğümle kaydettiği için, X düğümünde doğan
bir mutasyonu tam olarak X'in altındaki yapraklar taşır — yani gerçek destek örüntüsü her zaman
ağacın bir **kladı**dır, asla rastgele bir alt küme değil. Bu, çok-örnekli iddiayı koordinat
gerektirmeyen bir teste çeviriyor.

## Komutlar

```
tmp/.venv/bin/python tmp/overnight/eval_phylo.py \
    simulation/phylo/snape_out/calls.vcf --json results/exp3_eval_phylo.json
tmp/.venv/bin/python tmp/overnight/build_tree.py \
    simulation/phylo/snape_out/calls.vcf --boot 500 --newick results/exp3_tree.nwk
tmp/.venv/bin/python tmp/night/exp3_plot.py
```

## Sonuçlar

| ölçüm | değer | null (etiket permütasyonu, 200 tekrar) |
|---|---|---|
| klad oranı, tüm çağrılar | **5.655/6.425 = %88.0** | %70.7 [70.4–73.7] |
| klad oranı, bilgilendirici örüntüler (1 < \|destek\| < 20) | **1.135/1.905 = %59.6** | **%1.0** [0.0–11.3] |
| klad kapsamı | **39/39** klad en az bir çağrıyla tanıklandı | — |
| split geri kazanımı | **18/18** gerçek split doğrudan gözlendi | 316 sahte |

Ağaç kurma, beş yöntem, gerçek topolojiye Robinson–Foulds:

| yöntem | klad | ortak | RF | kapsam | kesinlik |
|---|---|---|---|---|---|
| **greedy-split** | 17 | 17 | **0** | **1.000** | **1.000** |
| nj-jaccard | 18 | 17 | 2 | 0.944 | 0.944 |
| greedy-clade | 18 | 16 | 4 | 0.889 | 0.889 |
| upgma | 18 | 14 | 8 | 0.778 | 0.778 |
| nj-hamming | 18 | 13 | 10 | 0.722 | 0.722 |

Bootstrap (500 tekrar, greedy-split): 17 gerçek kladın **13'ü %100**, 16'sı %70 üzeri destekli.
En zayıf klad `T15,T16,T17,T18,T19,T20` (%60.8).

## Karar

**Kabul.** Destek örüntüleri tek başına gerçek topolojiyi birebir kuruyor (RF = 0), hiçbir
koordinat ya da hizalama kullanılmadan. Bilgilendirici örüntülerdeki %59.6'ya karşı null'un %1.0
olması, bunun tesadüf olmadığını gösteren asıl sayı.

**Dürüstlük notu:** `eval_phylo.py`'de bir değişken gölgeleme hatası vardı (296. satırdaki döngü
`got` adını yeniden bağlıyor, sonra JSON dökümü `len(got)` diyip patlıyordu). Konsol çıktısı
etkilenmiyordu. Döngü değişkeni `ngot` yapıldı; orijinali `tmp/night/eval_phylo.py.orig`.

**Çıktılar:** `results/exp3_eval_phylo.{txt,json}`, `results/exp3_build_tree.txt`,
`results/exp3_tree.nwk`, `figs/exp3_phylogeny.{pdf,png}`
