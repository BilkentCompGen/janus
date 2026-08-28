# Hakem turu 2 — editör: dil, tutarlılık, metin↔tablo↔figür uyumu

Rol: sonucu değil sunumu denetleyen editör. Sayılar üç yerde geçiyor (deney dosyası, HANDOFF,
HTML rapor) — üçünün birbirini tutması gerekiyor.

---

**E1. Aynı sayı üç yerde, hepsi tutuyor mu?**

Kontrol edildi: deney 2'nin dört F1'i (`0.8704 / 0.8616 / 0.8851 / 0.8874`),
deney 3'ün `RF = 0` ve `%88.0`'i, deney 4'ün `1.93 / 3.19 / 13.08 sa`'i,
deney 5'in `+%10.7` ve `367 / 3.008`'i — `experiments/*.md`, `HANDOFF.md` ve
`results/gece_raporu.html`'de aynı. **Uyumlu.**

**E2. "Ölçüldü" ile "ölçeklendi" her yerde ayrılmış mı?**

`results/exp4_summary.tsv`'nin `kaynak` sütunu satır satır `OLCULDU` / `OLCEKLENDI` diyor;
figürün altındaki not 21 genomun gerçekten hizalanmadığını yazıyor; deney dosyası aynı şeyi
tekrar ediyor. Deney 1'in figür notu `ALTORNEKLENMIS` ile başlıyor. **Uyumlu.**

Tek eksik: `HANDOFF.md`'nin "Başlıca sayılar" tablosunda deney 4 satırında `ölçeklendi` kelimesi
geçmiyor, yalnız süreler var. Okuyan biri bunun ölçüm sandığı olabilir.
**Düzeltilmeli.**

**E3. Dil tutarlılığı.**

Belgeler Türkçe, kod yorumları Türkçe (mevcut repo kuralına uygun — `tmp/overnight` altındaki
betikler de öyle), figür etiketleri Türkçe. Paper İngilizce; **figürler paper'a girerken
etiketlerin çevrilmesi gerekecek.** Bu bir kusur değil, sabaha not.

**E4. Sınırlar her iddianın yanında mı, yoksa yalnız sonunda mı?**

`benchmarks_gaps.md` toplu listeyi veriyor ama asıl mesele iddianın yanında durması.
Kontrol: deney 5'in "Karar" bölümü üyelik primini iddiayla aynı paragrafta veriyor ✓;
deney 4 RSS'in ölçülmediğini "Bilinen boşluklar" başlığıyla veriyor ✓; deney 1'in yorumlama
sınırı kendi başlığında ✓; deney 2 hcc1937'nin ters yönde olduğunu tablonun hemen altında
söylüyor ✓. **Uyumlu.**

**E5. Figürler tek başına okunabiliyor mu?**

Yedi figürün hepsinde altta kaynak notu var: hangi dosya, hangi yöntem, ölçüldü mü ölçeklendi mi.
Başlıklar ne gösterdiğini söylüyor. **Uyumlu.**

Eksik: `exp5b_consensus` figüründe üyelik priminin büyüklüğü (~0.12) yazılı değil, yalnız
çubuklardan okunuyor. Notta var ama sayı yok. **Küçük, düzeltilmedi.**

**E6. HANDOFF sabah ilk işi net söylüyor mu?**

`git push --force origin changes` en üstte, gerekçesiyle. ✓

---

## Bu turdan çıkan ve YAPILAN düzeltme

- E2: `HANDOFF.md`'nin deney 4 satırına ölçekleme uyarısı eklendi.

## Yapılmayan

- E5: figüre sayı eklemek. Küçük, ve sweep hâlâ koşarken figür üretimi tekrarlanacak.
- E3: figür etiketlerinin İngilizceye çevrilmesi — paper'a girerken yapılacak iş, bu gecenin değil.
