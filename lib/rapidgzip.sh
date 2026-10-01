#!/usr/bin/env bash
# rapidgzip.sh: fetches rapidgzip 0.16.0 (parallel gzip decompression, used by snape pack) and builds its bundled ISA-L
# as lib/rapidgzip/libisal.a. Needs nasm. Run once from the repository root before make.
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p rapidgzip; cd rapidgzip
[ -d src ] || { curl -sSL https://files.pythonhosted.org/packages/source/r/rapidgzip/rapidgzip-0.16.0.tar.gz | tar xz; mv rapidgzip-0.16.0/src .; }
I=src/external/isa-l; O=obj; mkdir -p $O
for f in igzip/igzip_inflate.c igzip/igzip.c igzip/hufftables_c.c igzip/encode_df.c igzip/igzip_icf_base.c igzip/igzip_icf_body.c \
         igzip/igzip_base.c igzip/adler32_base.c igzip/flatten_ll.c igzip/huff_codes.c crc/crc_base.c; do
    gcc -O3 -fPIC -c $I/$f -I$I/include -I$I/igzip -I$I -o $O/$(basename $f .c).o
done
for f in include/reg_sizes.asm include/multibinary.asm igzip/igzip_inflate_multibinary.asm igzip/igzip_decode_block_stateless_01.asm \
         igzip/igzip_decode_block_stateless_04.asm igzip/rfc1951_lookup.asm igzip/stdmac.asm igzip/igzip_deflate_hash.asm igzip/igzip_body.asm \
         igzip/igzip_multibinary.asm igzip/igzip_update_histogram_01.asm igzip/igzip_update_histogram_04.asm igzip/encode_df_04.asm \
         igzip/encode_df_06.asm igzip/igzip_finish.asm igzip/igzip_gen_icf_map_lh1_04.asm igzip/igzip_gen_icf_map_lh1_06.asm \
         igzip/igzip_icf_body_h1_gr_bt.asm igzip/igzip_icf_finish.asm igzip/igzip_set_long_icf_fg_04.asm igzip/igzip_set_long_icf_fg_06.asm \
         igzip/proc_heap.asm igzip/adler32_avx2_4.asm igzip/adler32_sse.asm crc/crc_multibinary.asm crc/crc32_gzip_refl_by16_10.asm \
         crc/crc32_gzip_refl_by8_02.asm crc/crc32_gzip_refl_by8.asm; do
    nasm -f elf64 -I$I/include/ -I$I/igzip/ -I$I/ $I/$f -o $O/$(basename $f .asm).asm.o
done
ar rcs libisal.a $O/*.o
echo "lib/rapidgzip/libisal.a built"
