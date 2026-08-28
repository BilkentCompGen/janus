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
