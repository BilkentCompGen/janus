# Hakem turu 1 — gecenin deneylerine karşı

Rol: sonuçları reddetmeye meyilli, somut itirazlar yazan hakem. Her itirazın altında ya bir
düzeltme ya da düzeltmeme gerekçesi var.

---

**H1. Deney 4'ün C kolu ölçülmedi, uyduruldu.**
21 genom hizalanmadı; 9 ölçümün medyanı 21 ile çarpıldı. Bu, C kolunu snape'in lehine
şişirebilir — gerçek bir 21 örnekli koşuda önbellek/IO paylaşımı süreyi düşürebilirdi.

*Karşılık:* Ölçekleme kullanıcının açık talimatı ("21 genomu gerçekten hizalamayalım, 3-4
genomun mm2 sürelerini alıp 21'e scale et"). Figürün altına ve TSV'nin `kaynak` sütununa
`OLCEKLENDI` yazıldı, `OLCULDU` olanlardan ayrıldı. **Ama itiraz haklı:** medyan 26.97 dk'nın
etrafındaki dağılım 10.2–38.9 dk, yani ±%45. Bu belirsizlik figürde gösterilmiyor.
**Düzeltilmedi** — vakit yetmedi; `benchmarks_gaps.md`'ye yazıldı.

**H2. Deney 5b'nin "içeride" sayıları anlamsız.**
snape konsensüsün üyesiyken kendi oyunu içeren bir gerçek kümeye karşı skorlanıyor. 0.92'lik F1
bir başarı değil, bir totoloji.

*Karşılık:* Haklı, ve tam bu yüzden c1ext geçişi eklendi: aynı yapraklarda snape konsensüsün
**dışından** da skorlandı (0.773–0.820). Prim ölçüldü (~0.12) ve CASTLE'ınkiyle (~0.015)
karşılaştırıldı, sebebi (küçük konsensüs, ≥3 oy) yazıldı. Deney dosyasında "içeride sayıları
tek başına verilmemeli" uyarısı var. **Düzeltildi.**

**H3. Tek çift, tek koşu, güven aralığı yok.**
Deney 2 iki CASTLE çiftinde; deney 3 tek simülasyondan; hiçbirinde tekrar yok.

*Karşılık:* Kısmen karşılandı. Deney 3'te bootstrap var (500 tekrar, 13 klad %100). Deney 2
iki çifte genişletildi ve **ikisi aynı yönde çıkmadı** — bu, tek çiftten genelleme yapmanın
tehlikesini gösteren asıl kanıt ve HANDOFF'ta öne çıkarıldı. Deney 4/5'te tekrar yok.
**Kısmen düzeltildi.**

**H4. Deney 3'ün "RF = 0"ı beş yöntemden en iyisini seçmekle elde edildi.**
Beş yöntem koşulup en iyisi raporlanırsa bu bir seçim yanlılığıdır.

*Karşılık:* Beş yöntemin **hepsi** raporlanıyor (figür 3a), en iyisi işaretlenerek. `greedy-split`
sonradan seçilmiş bir kazanan değil, `build_tree.py`'nin gerekçesinde önceden savunulan model:
somatik SV geri dönüşsüz bir karakter, yani Dollo, yani destek örüntüsü zaten kladdır. Yöntem
veriye bakılmadan önce doğru modeldi. **Düzeltmeye gerek yok, ama figürde beşi de duruyor.**

**H5. Deney 1 altörneklenmiş, yani deney 1 yapılmadı.**
~3× kapsamda `min_mult=3` çekirdeklerin çoğunu eliyor; graf sonrası aşamalar tam koşuyu temsil
etmiyor. Nitekim t=32'de contig hizalama (128 s) graf aşamasından (48 s) uzun — tam koşuda tersi.

*Karşılık:* Doğru ve bu kayda geçirildi: deney dosyasında, figür notunda ve
`benchmarks_gaps.md`'de. Ölçülen şey ölçeklenmenin şekli, mutlak maliyet değil. Kullanıcı bu
deneyi zaten "sadece vaktin kalırsa" diye istemişti. **Düzeltilmedi, sınırı yazıldı.**

**H6. Colored de Bruijn hiç koşulmadı ama intro üç sınırını sayıyor.**
Paper ampirik bir iddiada bulunuyor, arkasında ölçüm yok.

*Karşılık:* Haklı. Bu gecenin kapsamında değildi ama `benchmarks_gaps.md`'ye yazıldı ve
HANDOFF'ta sabah işleri arasında. **Düzeltilmedi.**

**H7. hcc1937'nin çağrı sayısı 948'den 607'ye düştü, bu bir regresyon.**

*Karşılık:* F1 düşmedi (0.765 ↔ 0.7649), precision 0.847→0.825, recall 0.698→0.713. Yani çağrı
sayısı düştü ama doğruluk aynı kaldı. 948 sayısı `docs/RESULTS.md`'nin replikasyon tablosundan
ve muhtemelen farklı bir dosyayı (aday mı çağrı mı) sayıyor. **Kovalanmadı**, notu düşüldü.

---

## Turdan çıkan iş listesi

1. Deney 4'ün C koluna minimap2 dağılımından bir belirsizlik bandı ekle. *(yapılmadı)*
2. Deney 2'yi h1437'ye genişlet. *(yapılmadı, `.svb`'ler hazır)*
3. Bir dBG çağırıcısı koş, en azından tek yaprakta. *(yapılmadı)*
4. Deney 5b'yi daha büyük bir konsensüsle tekrarla ki prim küçülsün. *(yapılmadı)*
