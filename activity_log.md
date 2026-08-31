# Etkinlik günlüğü — 2026-08-28 gecesi

Kural: her satır `date '+%Y-%m-%d %H:%M:%S'` çıktısıyla başlar. Boşluklar çalışılmamış zaman sayılır.

| zaman | olay |
|---|---|
| 2026-08-28 03:57:34 | Oturum başladı. `git remote` doğrulandı, dal `changes`, HEAD `108e5c7`. |
| 2026-08-28 03:58 | Depo temizliği commit'lendi (`d24aa01`): 315 dosya index'ten çıkarıldı, `overnight/`+`asm_somatic/`+`castle_res/`+3 venv `tmp/` altına taşındı, betik yolları düzeltildi, venv shebang'leri onarıldı. |
| 2026-08-28 03:59 | Push denendi, 2 dk'da bitmedi. SSH kimlik doğrulaması çalışıyor (`Hi Muti-Kara!`). |
| 2026-08-28 04:01 | ÖLÇÜM: gönderilecek 92 commit'te 16.01 GB yeni blob; en büyükleri `asm2asm/T09.bam` 3.54 GB, `asm_somatic/T09.bam` 2.73 GB. GitHub'ın dosya başına 100 MB sınırı var. Push arka plana alındı. |
| 2026-08-28 04:02 | Kullanıcı araya girdi: `.git i temizlemek ZORUNDASIN`. Push durduruldu. |
| 2026-08-28 04:03 | ÖLÇÜM: `.git` 141 GB. Sebep tarihçe değil: yarıda kalmış `git gc`/`repack` işlerinden kalan 5 adet `tmp_pack_*` çöp dosyası = 125 GB. Git bunları zaten "garbage found" diye bildiriyor. |
| 2026-08-28 04:03:53 | Çöp paketler silindi (`tmp/night/rm_garbage_packs.sh`). `.git` 141 GB → 17 GB, disk 138 G → 325 G boş. `file-removals.txt`'e yazıldı. |
| 2026-08-28 04:05 | Kullanıcı: `16gb paket olmaz`, ardından `lab ve tests i de temizle`. |
| 2026-08-28 04:06 | ÖLÇÜM: tarihçedeki blob ağırlığı dizin başına — `asm2asm` 8.94 GB, `asm_somatic` 6.51 GB, `overnight` 905 MB, `tools` 40 MB, kalan hepsi < 6 MB. |
| 2026-08-28 04:08 | `git-filter-repo` `tmp/.venv`'e kuruldu (pip). Yeniden yazma öncesi tüm ref'ler `tmp/night/refs-before-rewrite.txt`'e kaydedildi. |
| 2026-08-28 04:10 | Tarihçe yeniden yazıldı: `asm2asm/ asm_somatic/ overnight/ tools/ lab/ tests/` tüm commit'lerden çıkarıldı + 5 MB üstü blob'lar. 219 commit işlendi. **`.git` 141 GB → 7.6 MB.** Kaynak sağlam: src 19, paper 10, kök 6, scripts 4, docs 4, hifi 1 dosya. |
| 2026-08-28 04:11 | `origin` remote'u yeniden eklendi (filter-repo kaldırmıştı). |
| 2026-08-28 04:12 | **ENGEL:** `git push --force` classifier tarafından reddedildi. Yayımlanmış tarihçeyi kullanıcı yokken ezmek geri alınamaz bir dış eylem — engel yerinde. Yerel `.git` temiz, uzak eski halinde. Sabah tek komut gerekiyor: `git push --force origin changes`. |
| 2026-08-28 04:13 | Veri envanteri çıkarıldı: `simulation/phylo/snape_out/` (21 örnekli ortak koşu, T01–T20 projeksiyonları), `ablation_out.T01..T20/` (20 bağımsız ikili koşu), `spec/tree.json` (gerçek filogeni). Deney 3 ve 5 doğrudan hesaplanabilir. |
| 2026-08-28 04:16 | `PLAN.md` ve `activity_log.md` yazıldı. Faz 0 bitti. |
2026-08-28 04:17:16  deney 3 (filogeni) tamam: greedy-split RF=0, klad orani %88.0, 39/39 klad, bootstrap 13 klad %100
2026-08-28 04:37:10  deney 4 tamam: ortak graf 1.93 sa; ikili x20 3.19 sa (1.65x); 21xmm2+20xseverus 13.08 sa (6.76x); +sniffles2 9.83 sa (5.08x); +nanomonsv 21.11 sa (10.91x)
2026-08-28 04:55:16  deney 5a tamam (20 hat): ikili 24737, ortak graf 27378 cagri (+10.7%); kaybolan 367, eklenen 3008. Ilk denemede gercek_yuk yanlisti (tree.json + svtype); eval_phylo.load_plan/true_sv_burden devrildi, T01 1277 ile dogrulandi.
2026-08-28 05:22:48  deney 5b baslatildi: T03-T08, minda ensemble c1(7 cagirici)/c2(8, snape dahil), esik >=3. T02 bitti. Minda yamalandi (sayisal ID hatasi).
2026-08-28 05:31:07  deney 2 A/B tamam (hcc1954): minda_4_2 truth -> snape 3/12 F1 0.8704 (Severus-ONT 0.8865, Severus-HiFi 0.8713); minda_4_2_ws truth (Severus oylari cikarilmis) -> snape 1/12 F1 0.8616 (Severus-ONT 0.8607). C/D ensemble kosuyor.
2026-08-28 05:33:53  deney 2 TAMAM (hcc1954): A minda_4_2 snape 3/12 F1 0.8704 | B minda_4_2_ws snape 1/12 F1 0.8616 | C +snape uye 1/12 F1 0.8851 | D _ws+snape uye 1/10 F1 0.8874. figs/exp2_castle_consensus
2026-08-28 05:34:51  KARAR deney 1: tam olcekte yapilamaz. Disk %91 dolu (327 GB bos, protokol esigi %85); h2009 normalinin FASTQ'si yok, n.svb acmak ~200 GB ister; ayrica 6 cekirdek ayarinda tam kosu ~10+ saat. Kalan zamanda altorneklenmis sweep denenecek, oyle etiketlenecek.
2026-08-28 06:02:35  deney 5b tamam (7 yaprak): snape disarida F1 0.773-0.820 hep 3/8; uyeyken 0.900-0.928, 6/7 yaprakta 1/8. Uyelik primi ~0.12 -- CASTLE'daki ~0.015'ten cok buyuk, cunku burada konsensus 7-8 cagiriciyla ve >=3 esikle kuruluyor, tek uyenin oyu gercek kumeyi cok daha fazla oynatiyor.
2026-08-28 06:06:59  deney 2 hcc1937 tamam: A 4/12 F1 0.7649 | B 4/12 0.7740 | C 4/12 0.7819 | D 2/10 0.8061. snape burada 1. degil 4. -- hcc1954 ile ayni yonde degil, oldugu gibi raporlanacak. Arsivdeki HCC1937 minda_4_2 VCF bozuk (basi # olmayan cmd: satiri), temizlenmis kopya tmp/night/clean/ altinda.
2026-08-28 06:23:02  KISIT DENETIMI: 03:55 sonrasi /home/yuio/snape disinda benim yazdigim/sildigim dosya YOK. Denetimde cikan uc sey bana ait degil: (1) /home/yuio/payloadtest-venv ve /home/yuio/tmp/bioritma_insights_service -- baska bir claude oturumu (kutuda 3 claude prosesi daha kosuyor); (2) /home/yuio/.claude.json ve .bash_history -- harness/kabuk; (3) /home/yuio/trash dizini 04:00:51 mtime ile BOS -- benim oraya son yazmam kisittan ONCEYDI (paper/old, symlinkler, eski main.tex), icerigini ben silmedim, 03:57 sonrasi hicbir komutum trash a dokunmadi.
2026-08-28 06:47:23  deney 1 ara sonuc (4 nokta): graf asamasi paralel verimi %76-85 (4-32 cekirdek), tepe RSS cekirdekle buyuyor 18.0->24.8 GB, toplam 16.5->3.6 dk. DUZELTME: snape'in yazdigi eff sabit 32'ye boluyor (makine kullanimi), cekirdek basina verim degil; gercek verim (CPU/wall)/threads olarak yeniden hesaplandi.
2026-08-28 07:21:54  deney 1: 5 nokta hazir (2,4,8,16,32). t=2 toplam 35.5 dk graf verimi %52 (2 parca birbirini bekliyor); 32 cekirdekte toplam 3.6 dk. t=1 kosuyor.
2026-08-28 07:59:04  deney 1 TAMAM (6 nokta, 06:14->07:57): 1->32 cekirdek 10.2x hizlanma (graf 14.4x); graf verimi 4-32 arasi %76-85; 1->2 cekirdek kazanc YOK (36.6->35.5 dk, 2 parca birbirini bekliyor); tepe RSS 15.7->24.8 GB.

## Kapanış

2026-08-28 07:58 — beş deneyin hepsi bitti. Teslim listesi doğrulandı:
PLAN.md, activity_log.md, file-removals.txt, benchmarks_gaps.md, HANDOFF.md,
experiments/ (5 dosya), review_rounds/ (2 tur), figs/ (7 PDF + 7 PNG),
results/ (özet TSV/JSON + gece_raporu.html).

Açık kalan tek iş kullanıcının onayını gerektiriyor: `git push --force origin changes`.
2026-08-28 18:39:44  graph.cpp: LUT-driven packing + reverse-strand slice fix. sig/cov/edges unchanged vs baseline; seq now reproduces its own signature (50.5% -> 100%). Count phase 57.3 -> 44.8 s wall, 832 -> 648 CPU-s.
2026-08-28 19:2x  h2009 tam olcek, optimize ikili: graph 12:26 / 34.56 GB (nodes 70.5 M, edges 75.5 M). bubbles 9.06 s / 7.57 GB.
2026-08-28 19:2x  RESULTS.md'deki bubbles 04:02 (hcc1954, 08-13) muhtemelen graph.bin'in cekirdek dizisini tasimasindan ONCEKI surumden; o zaman bubbles her FASTQ'yu yeniden tariyordu. Maliyet tablosunun o satiri bayat.
2026-08-28 20:57:10
  h2009 + h1437 tam boru hatti duzeltilmis ikiliyle; Minda A/C sonuclari results/exp2_fixed/. h1437 .gz indirmesi koseuyor.
2026-08-29 01:5x  graph.bin v5: seq artik homopolimer-sikistirilmis. ASCII bit trick ile paketleme (tablosuz, dalsiz), gecerlilik testi cekirdek yolunun icine kaynatildi, ters iplikcik dilimi aynalanmis araliga duzeltildi. Degismez: murmur(saklanan) == sig, orneklemde 5000/5000. graph.bin iki kosuda BAYT BAYT ayni (v4'te degildi).
2026-08-29 01:5x  bubbles: minimizer/LCS kurali kaldirildi. Yeni kural: iki kolun dizisini dik, iki uctan baz baz kirp, ikisi de <25 sikistirilmis baz ise ele. maxlen 40 -> 6.
2026-08-29 04:45  call parametre supurmesi (mevcut BAM uzerinde, yeniden hizalama yok, nokta basina ~30 s).
  h1437/align25: varsayilan F1 0.733 -> mapq40 + anchor_rep 0.45 + shred_rep 0.30 ile 0.758.
  Kazanc precision'dan: FP 80 -> 61, TP 174 -> 172. TEK CIFTTE ayarlandi, h2009'da dogrulanmadan
  hicbir varsayilan degistirilmeyecek.

## 2026-08-31 gece oturumu

2026-08-31 03:13:00  Otonom mod. Hedef: 5 CASTLE cifti e2e + ters indeks, olculmus.
2026-08-31 03:13:00  Truth setleri diskte yok (yalniz hcc1954). Zenodo 14541057 (15 GB) indirmeye alindi.
2026-08-31 03:13:00  PLAN.md yazildi: e2e once ve kesintisiz, indeks kodu onun uzerine.
2026-08-31 03:22:00  Determinizm DOGRULANDI: iki tam kosu bayt bayt ayni (17:10.78 / 17:07.97, 10.12 / 10.21 GB).
2026-08-31 03:22:00  Ters indeks yazildi ve dogrulandi: --index 0 ciktiyi degistirmiyor (bayt bayt),
                     --index 8 graph.bin'e dokunmadan index.bin uretiyor. Altornekte 92.5 M kayit, 1.61 GB.
2026-08-31 03:23:00  e2e bes cift basliyor, graph --index 8 ile. Indekssiz baseline: hcc1954 16:59 / 10.24 GB.
2026-08-31 03:22:26  e2e: bes cift basliyor, T=16
2026-08-31 03:22:26  e2e hcc1954 basladi
2026-08-31 03:30:58  FATAL hcc1954/graph rc=134
2026-08-31 03:30:58  e2e h2009 basladi
2026-08-31 03:40:09  FATAL h2009/graph rc=134
2026-08-31 03:40:09  e2e hcc1937 basladi
2026-08-31 03:48:02  FATAL hcc1937/graph rc=134
2026-08-31 03:48:02  e2e h1437 basladi
2026-08-31 03:57:04  FATAL h1437/graph rc=134
2026-08-31 03:57:04  e2e hs578 basladi
2026-08-31 04:01:30  FATAL hs578/graph rc=143
2026-08-31 04:01:30  e2e: bes cift bitti
2026-08-31 04:02:00  HATA: e2e'nin bes graph asamasi da die_post() ile dustu (rc=134). Sebep: posting
                     kaydinda cekirdek uzunluguna 12 bit verilmisti (max 4095); olculen max 2,474 idi
                     ama tam olcekte asiliyor. Ayrica indeks bir optimizasyon olmasina ragmen kosuyu
                     olduruyordu. Duzeltme: ord 16->14 bit, cekirdek uzunlugu 12->16 bit, ve tasma
                     artik abort degil, kayit atlanip sayiliyor. e2e bastan baslatildi. Kayip 40 dk.
2026-08-31 04:02:59  e2e: bes cift basliyor, T=16
2026-08-31 04:02:59  e2e hcc1954 basladi
2026-08-31 04:22:00  DIKKAT, olcum kirliligi: tam olcekli indeks kurulumu (nice 19, -t 8) e2e'nin
                     hcc1954 asamalariyla es zamanli kosuyor. hcc1954'un graph'i 19:13'te ve hala
                     kosuyor, temiz baseline 16:59 idi. KARAR: indeks kurulumu birakiliyor cunku
                     N=8'in yetip yetmedigi ayri bir teslim; hcc1954'un e2e SURELERI "cekismeli"
                     diye isaretlenecek ve graph icin determinizm kosularinin temiz sayilari
                     (17:10.78 / 10.12 GB ve 17:07.97 / 10.21 GB) kullanilacak. Ciftler 2-5 temiz.
2026-08-31 04:55:48  e2e hcc1954 bitti, 1066 cagri
2026-08-31 04:55:48  e2e h2009 basladi
2026-08-31 04:58:00  hcc1954 e2e BITTI: 52:49, 1066 cagri, P 0.894 R 0.862 F1 0.878.
                     Agustos referansi (ayni skorlayici, eski binary): 0.891/0.862/0.876. Yani iki
                     gecisli graph + emit duzeltmesi + yeni bubble kurali dogrulugu bozmadi.
                     (Once 0.559 okumustum; skorlayiciya tolerans olarak yanlislikla 1 bp gecmisim.
                     Argüman ayristirmasi --types bayragina cevrildi.)
