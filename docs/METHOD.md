# Paperın tarif ettiği yöntem

`paper/sections/methods.tex`'in özeti. Parametre varsayılanları `src/main.cpp:170-171` ve
`src/main.cpp:26-37`'deki kullanımdan doğrulanmıştır; paper metniyle koddaki değerin ayrıldığı
yerler işaretlidir.

## Boru hattı

```
Extraction -> Graph building -> Bubble extraction -> Junction extraction -> Junction alignment
           -> Candidate calling -> Normal extraction -> Normal alignment -> Variant calling
```

Tablo 2'nin aşama adları bunlardır; Extraction ayrı bir aşama olarak zamanlanmıyor, graph
building'in içinde.

## Girdi

Üç sütunlu, boşlukla ayrılmış bir örnek tablosu: ad, rol, yol. Her okuma seti bir satır. Tam bir
satırın rolü `normal`, en az bir satırınki `tumor` olmalı. En fazla 63 tümör. Ad VCF'te
raporlamada kullanılıyor, satır başına benzersiz olmalı. Yollar FASTQ dosyasını göstermeli.

## Extraction

Bütün okumalar tek geçişte core dizisine çevriliyor: homopolimer sıkıştırma, LCP, MurmurHash64A.
Bu hash grafın düğüm kimliği. FASTQ strand bilgisi taşıyamadığı için aynı adımlar okumanın ters
tümleyenine de uygulanıyor; yani her okuma grafta iki kez temsil ediliyor.

LCP seviyesi varsayılan `--level 5`, margin varsayılan 10. Paperın bütün sonuçları bu
varsayılanlarla alınmıştır ve metin seviyeyi bir parametre olarak tartışmıyor.

## Graph building

Graf örnek örnek kuruluyor. Her örnek için core'lar yukarıdaki adımlarla çıkarılıyor, sonra o
örnekte 3 kereden az görülen core'lar okuma hatası sayılıp atılıyor (`--min-mult 3`). Kalan
core'lar ve aralarındaki kenarlar örnek başına sayımlarıyla grafa ekleniyor. Sayımlar bir bayt ve
255'te doyuyor (`snape.h`: `using cov_t = uint8_t; constexpr cov_t COV_MAX = 255;` — paper metni
kodla uyumlu). Her düğüm ayrıca core'unun baytlarını homopolimer-sıkıştırılmış ve iki-bit paketli
olarak tutuyor; baytları, o core'u ilk terfi ettiren örnek bağışlıyor.

## Bubble extraction

Çıkış derecesi en az iki olan her core için bütün çıkan kenarlar geziliyor. Her kol, core'ların
çıkış ve giriş dereceleri 1 olduğu sürece derinleştiriliyor. Bir yapının kabarcık sayılması için:
kollar aynı sink'e dönmeli, tek renkli olmalı, 7 core'dan uzun olmamalı ve en az iki kolun rengi
farklı olmalı.

Her kolun dizisi, ardışık core'lar en uzun sonek-önek örtüşmelerinden birleştirilerek yeniden
kuruluyor. Örtüşme yoksa core'lar olduğu gibi bırakılıyor. Varsa iki kol birbirine karşı iki
uçtan, baz baz, her iki taraftaki ilk uyumsuzluğa kadar kırpılıyor ve kalan diziler birbirine
hizalanıyor. Diziler arasındaki fark 25 bp'den azsa SV değil diye atılıyor.

## Junction extraction

Bir somatik olay tek bir core değil, bir somatik core zinciri üretiyor; bu yüzden snape tek tek
core'lara değil bağlı öbeklere bakıyor. Drop listesinden sağ çıkan, birbirine bağlı somatik
core'lar bir junction adayı oluşturuyor. Birden çok tümör hattı olduğunda junction ortak olduğu
için yalnız bir kez tanımlanıyor.

## Junction alignment

snape girdi FASTQ'larını bir kez daha işliyor ve bir junction'ı taşıyan okumaları, junction artı
her iki yanda 50 germline core olacak şekilde kırpıyor (`pad_cores=50`). Aynı okumadaki iki
junction 100 core'dan yakınsa tek contig'de birleştiriliyor. Bir contig, iki tarafta da 300 bp
flank tutturuyorsa yazılıyor (`min_flank=300`) ve junction başına hat başına en fazla 5 tane
(`per_chunk=5`). Referansa hizalama minimap2 ile.

## Normal extraction ve alignment

Germline sorusu yalnız aday lokuslarda soruluyor, bu yüzden snape caller'ı önce normalsiz
koşturup lokusları buluyor (Candidate calling). Her aday kendi contig'ini adlandırıyor; o
contig'in germline flankları normalin de taşıdığı dizi olduğundan, o lokustaki normal okumalar
onlarla core paylaşan okumalar. Yalnız onlar seçiliyor, normalin geri kalanı hiç hizalanmıyor.
Seçilen normal okumalar minimap2 ile referansa hizalanıyor.

## Variant calling

snape bütün contig hizalamalarını okuyup MAPQ ve sübstitüsyon oranlarına bakıyor: MAPQ en az 10,
sübstitüsyon oranı %1'in altında. Bu aşama standart bir uzun okuma kırılma noktası caller'ı.
Contig'ler üzerinde çalıştığı için minimal: bir contig kümesi pileup değil, modellenecek kapsam
yok. Graf adayları zaten somatik diziye daralttığı için bu aşama onları aramak yerine yerleştirip
doğruluyor.

Bir SV imzası, her biri bir pozisyon ve bir yön taşıyan bir kırılma noktası çifti olarak
tanımlanıyor; yön, hizalamanın referanstan ayrılma biçimi. P konumunda L uzunluğunda bir delesyon
(P,+):(P+L,-) çifti; çift tipi de veriyor. Bu, Severus'un da kullandığı gösterim. Bir
junction'ı kapsayan contig tek parça hizalanamayacağı için snape segmentlerini contig sırasına
dizip her ardışık çiftten bir imza alıyor. Contig'ler CIGAR'larında en az 50 bp'lik delesyon veya
insersiyon taşıdıklarında da imza veriyor. Contig'in kendi flankıyla eşleşen bir insersiyon tandem
duplikasyon, küçük indel patlaması ise — contig orada tekrarlı değilse — inversiyon.

Aynı olayı taşıyan contig'lerin her biri bir imza ürettiği için, iki kırılma noktasında da 50 bp
içinde uyuşan imzalar tek bir çağrıya atanıp kümenin medyanına yerleştiriliyor (`cluster=50`).
Destek, çağrıdaki farklı contig sayısı, hat başına sayılıyor; tek bir hat 5'e ulaşırsa çağrı
raporlanıyor (`--supp-discover 5`). Bir çağrı, flanklayan dizisi tandem tekrarsa düşürülüyor —
aligner'ın oraya koyduğu indel bir olay değil. Ayrıca yakında benzer uzunlukta bir indeli normal
de taşıyorsa, ya da oradaki normal okumaların dörtte biri %2'nin üstünde uyumsuzlukla
hizalanıyorsa düşürülüyor. Sağ kalan, tümör hattı başına bir sütunla tek bir VCF'e yazılıyor.

## Metinde tarif edilmeyen parametreler

Kodda varsayılanı olan ama Methods'ta geçmeyenler:

| Parametre | Varsayılan | Nerede |
|---|---|---|
| `--margin` | 10 | `main.cpp:170` |
| `--supp-assign` | 2 | `main.cpp:170`, hat başına genotip eşiği |
| `--index` | 8 | `main.cpp:171` |
| `min_vaf` | 0 | `junctions` kullanımı, süzgeç fiilen kapalı |
| `min_aln` | 0 | `call` kullanımı |

`index.bin` — core başına, o core'u taşıyan okumaların listesi — Methods'ta hiç tarif
edilmiyor; buna karşılık 2.4'te disk bütçesinin 8.5 GB'lık kalemi olarak geçiyor.
