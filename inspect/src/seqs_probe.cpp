// seqs_probe — Stage 0 of the inspect/ pipeline.
//
// Probe `outputs/svp/wgs/seqs.bin` (and also `posmap.bin`, which has conflicting
// schemas across the codebase) to confirm the actual byte layout BEFORE any
// downstream stage depends on it.
//
// Usage:  ./seqs_probe <dir>
//   <dir>  outputs/svp/wgs (contains counts.bin graph.bin posmap.bin seqs.bin)

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <sys/stat.h>

namespace {

[[noreturn]] void die(const char* msg) {
    std::fprintf(stderr, "[FATAL] %s\n", msg);
    std::abort();
}

uint64_t fsize(const std::string& p) {
    struct stat st{};
    if (::stat(p.c_str(), &st) != 0) {
        std::fprintf(stderr, "[FATAL] stat %s\n", p.c_str());
        std::abort();
    }
    return (uint64_t)st.st_size;
}

void hex_dump(const std::string& path, uint64_t offset, size_t nbytes, const char* label) {
    std::printf("\n== %s @ offset %lu (%zu bytes) ==\n", label, (unsigned long)offset, nbytes);
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) die("fopen");
    std::fseek(f, (long)offset, SEEK_SET);
    std::vector<uint8_t> buf(nbytes);
    size_t got = std::fread(buf.data(), 1, nbytes, f);
    std::fclose(f);
    for (size_t i = 0; i < got; i += 16) {
        std::printf("%08lx  ", (unsigned long)(offset + i));
        for (size_t j = 0; j < 16; j++) {
            if (i + j < got) std::printf("%02x ", buf[i + j]);
            else             std::printf("   ");
            if (j == 7) std::printf(" ");
        }
        std::printf(" |");
        for (size_t j = 0; j < 16 && i + j < got; j++) {
            uint8_t c = buf[i + j];
            std::putchar((c >= 32 && c < 127) ? (char)c : '.');
        }
        std::printf("|\n");
    }
}

uint64_t rd_u64(const std::string& path, uint64_t off) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) die("fopen");
    std::fseek(f, (long)off, SEEK_SET);
    uint64_t v;
    if (std::fread(&v, 8, 1, f) != 1) die("read u64");
    std::fclose(f);
    return v;
}

uint32_t rd_u32(const std::string& path, uint64_t off) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) die("fopen");
    std::fseek(f, (long)off, SEEK_SET);
    uint32_t v;
    if (std::fread(&v, 4, 1, f) != 1) die("read u32");
    std::fclose(f);
    return v;
}

void probe_graph_K(const std::string& dir, uint64_t& K) {
    const std::string p = dir + "/graph.bin";
    K = rd_u64(p, 0);
    std::printf("[graph.bin] K = %lu nodes\n", (unsigned long)K);
}

// ---------- posmap.bin ----------
//
// Two schemas exist in the codebase:
//   A) "viz_sv.py" schema   (chrom-names FIRST, then per-sig variable):
//        u64 nc; nc * { u32 name_len; char name[L]; u64 chrom_len; }
//        u64 ns; ns * { u64 sig; u32 npos; npos * { u16 ci; u16 _; u32 ps; u32 pe; u8 strand; u8 _; u16 _; } }
//   B) src/posmap.h schema  (16 B header + flat records, names at trailer):
//        u64 magic = 0x53565050'4D415000
//        u64 n_records
//        n_records * { u64 sig; u32 chrom_id; u32 pos_start; u32 pos_end; u8 strand; u8 pad[7]; }   (32 B)
//        u32 n_chroms; n_chroms * { u32 name_len; char name[L]; u64 chrom_len; }
//
// We probe both. The first u64 distinguishes:
//   - A: a small count (~455 chroms) → small number
//   - B: 0x53565050'4D415000 ≈ 6.0e18 → huge specific magic constant
void probe_posmap(const std::string& dir) {
    const std::string p = dir + "/posmap.bin";
    uint64_t sz = fsize(p);
    std::printf("\n========== posmap.bin (%lu bytes) ==========\n", (unsigned long)sz);
    hex_dump(p, 0, 64, "posmap.bin head");

    const uint64_t u0 = rd_u64(p, 0);
    std::printf("[posmap.bin] u64[0] = 0x%016lx = %lu\n", (unsigned long)u0, (unsigned long)u0);

    constexpr uint64_t MAGIC_B = 0x53565050'4D415000ULL;
    if (u0 == MAGIC_B) {
        std::printf("[posmap.bin] MATCHES schema B (src/posmap.h magic)\n");
        uint64_t n_rec = rd_u64(p, 8);
        std::printf("[posmap.bin]   n_records = %lu\n", (unsigned long)n_rec);
        std::printf("[posmap.bin]   expected size of record region = %lu + 16 = %lu vs actual %lu\n",
                    (unsigned long)(n_rec * 32), (unsigned long)(n_rec * 32 + 16),
                    (unsigned long)sz);
        hex_dump(p, 16, 96, "first 3 records (32 B each)");
    } else if (u0 < 10000) {
        std::printf("[posmap.bin] LIKELY schema A (viz_sv.py format; nc=%lu)\n", (unsigned long)u0);
        uint64_t off = 8;
        for (uint64_t i = 0; i < u0 && i < 5; i++) {
            uint32_t L = rd_u32(p, off); off += 4;
            std::string name(L, '\0');
            FILE* f = std::fopen(p.c_str(), "rb");
            std::fseek(f, (long)off, SEEK_SET);
            if (std::fread(name.data(), 1, L, f) != L) die("name read");
            std::fclose(f);
            off += L;
            uint64_t clen = rd_u64(p, off); off += 8;
            std::printf("[posmap.bin]   chrom[%lu]: name='%s' len=%lu\n",
                        (unsigned long)i, name.c_str(), (unsigned long)clen);
        }
    } else {
        std::printf("[posmap.bin] UNKNOWN — does not match either schema\n");
    }
}

// ---------- seqs.bin ----------
//
// 14 GB / 85M nodes ≈ 165 B/node, consistent with L5/flank=10 median window (170 B).
// Candidate schemas:
//   X) u64 N + per-node { u8 seqlen; u8 seq[seqlen] }   (variable, sequential)
//   Y) u64 N + u64 offset[N+1] + flat u8 seq            (offset-indexed)
//   Z) other (custom magic, mmap-friendly header, etc.)
//
// We probe by:
//   1) dumping first 256 B as hex+ASCII
//   2) reading u64[0] and checking against graph.bin's K
//   3) trying schema X: u8 seqlen + seq[seqlen]; check that first seqlen records' bytes
//      look like ACGTN (uppercase) and the seqlen is plausibly 100..255
//   4) trying schema Y: u64 offsets; check monotonicity
void probe_seqs(const std::string& dir, uint64_t K) {
    const std::string p = dir + "/seqs.bin";
    uint64_t sz = fsize(p);
    std::printf("\n========== seqs.bin (%lu bytes ; %.2f B/node) ==========\n",
                (unsigned long)sz, (double)sz / (double)K);
    hex_dump(p, 0, 256, "seqs.bin head");

    const uint64_t u0 = rd_u64(p, 0);
    std::printf("[seqs.bin] u64[0] = 0x%016lx = %lu  (graph.bin K = %lu)\n",
                (unsigned long)u0, (unsigned long)u0, (unsigned long)K);

    // schema X probe: assume u64 N at offset 0, then variable records.
    if (u0 == K) {
        std::printf("[seqs.bin] u64[0] == K — consistent with header = node count\n");
    }

    // Schema X test: read u8 seqlen at offset 8, then seq[seqlen]
    {
        std::printf("\n-- trying schema X (u8 seqlen + seq[seqlen]) --\n");
        FILE* f = std::fopen(p.c_str(), "rb");
        std::fseek(f, 8, SEEK_SET);
        uint64_t off = 8;
        int acgtn_total = 0, byte_total = 0;
        for (int i = 0; i < 5 && off < sz; i++) {
            uint8_t L;
            if (std::fread(&L, 1, 1, f) != 1) break;
            off += 1;
            std::vector<uint8_t> seq(L);
            if (L) {
                if (std::fread(seq.data(), 1, L, f) != L) break;
                off += L;
            }
            int acgtn = 0;
            for (uint8_t b : seq) if (b == 'A' || b == 'C' || b == 'G' || b == 'T' || b == 'N') acgtn++;
            acgtn_total += acgtn;
            byte_total += L;
            std::printf("[X] rec %d: seqlen=%u  acgtn=%d/%u  bytes='%.*s'\n",
                        i, L, acgtn, L, (int)std::min<size_t>(L, 60u), (const char*)seq.data());
        }
        std::fclose(f);
        std::printf("[X] ACGTN ratio over first 5 records: %d / %d\n", acgtn_total, byte_total);
    }

    // Schema Y test: u64 N + u64 offset[N+1] + payload. Check offset[0]=0 and monotonicity.
    {
        std::printf("\n-- trying schema Y (u64 N + u64 offset[N+1] + payload) --\n");
        FILE* f = std::fopen(p.c_str(), "rb");
        std::fseek(f, 8, SEEK_SET);
        uint64_t off0; if (std::fread(&off0, 8, 1, f) != 1) { std::fclose(f); return; }
        std::printf("[Y] offset[0] = %lu (should be 0)\n", (unsigned long)off0);
        uint64_t off1; if (std::fread(&off1, 8, 1, f) != 1) { std::fclose(f); return; }
        std::printf("[Y] offset[1] = %lu (should be small positive, ~100-200)\n", (unsigned long)off1);
        uint64_t off2; if (std::fread(&off2, 8, 1, f) != 1) { std::fclose(f); return; }
        std::printf("[Y] offset[2] = %lu (delta from [1] = %lu)\n",
                    (unsigned long)off2, (unsigned long)(off2 - off1));
        std::fclose(f);
    }

    // Schema W test:  u64 K + u64 total_payload_size + u16 seqlen[K] + u8 payload[]
    // Hex dump suggested this: u64[1] looks like a total-payload size (matches file_size - 16 - K*2).
    {
        std::printf("\n-- trying schema W (u64 K + u64 total_payload + u16 seqlen[K] + u8 payload[]) --\n");
        uint64_t total = rd_u64(p, 8);
        uint64_t expected_header = 16 + 2 * K;
        uint64_t expected_total = sz - expected_header;
        std::printf("[W] declared total_payload = %lu\n", (unsigned long)total);
        std::printf("[W] expected total (file_size - 16 - 2*K) = %lu\n", (unsigned long)expected_total);
        std::printf("[W] match: %s\n", total == expected_total ? "YES — schema confirmed" : "NO");

        if (total == expected_total) {
            // Read first 8 seqlens
            FILE* f = std::fopen(p.c_str(), "rb");
            std::fseek(f, 16, SEEK_SET);
            uint16_t seqlens[8];
            if (std::fread(seqlens, 2, 8, f) != 8) die("seqlens read");
            std::fclose(f);
            std::printf("[W] first 8 seqlens: ");
            for (int i = 0; i < 8; i++) std::printf("%u ", seqlens[i]);
            std::printf("\n");

            // Read first record's payload (it starts at offset 16 + 2*K)
            uint64_t payload_off = 16 + 2 * K;
            f = std::fopen(p.c_str(), "rb");
            std::fseek(f, (long)payload_off, SEEK_SET);
            std::vector<uint8_t> buf;
            int dumped_records = 0;
            for (int i = 0; i < 8; i++) {
                buf.resize(seqlens[i]);
                if (seqlens[i] && std::fread(buf.data(), 1, seqlens[i], f) != seqlens[i]) {
                    std::fprintf(stderr, "[W] read failed at record %d\n", i);
                    break;
                }
                // Check ACGTN
                int acgtn = 0;
                for (uint8_t b : buf) if (b == 'A' || b == 'C' || b == 'G' || b == 'T' || b == 'N') acgtn++;
                std::printf("[W] rec %d: len=%u  acgtn=%d/%u  bytes='%.*s'\n",
                            i, seqlens[i], acgtn, seqlens[i],
                            (int)std::min<size_t>(seqlens[i], 80u), (const char*)buf.data());
                dumped_records++;
            }
            std::fclose(f);
            std::printf("[W] dumped %d records, all look like ACGTN strings? confirm visually above\n",
                        dumped_records);
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <svp_output_dir>\n", argv[0]);
        return 1;
    }
    const std::string dir = argv[1];

    uint64_t K;
    probe_graph_K(dir, K);
    probe_posmap(dir);
    probe_seqs(dir, K);
    return 0;
}
