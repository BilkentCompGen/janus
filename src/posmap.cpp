// posmap: reference-derived sig -> (chrom_id, pos) oracle. Extracts LCP cores from BOTH
// strands of the reference at the SAME level/flank as the LCP pass, hashes them with the
// same Murmur64 + seed=0 scheme as compute_sigs, and writes a sig-sorted binary for O(log N)
// lookup. RC pass is required because the caller's read pipeline does dual fw+RC extraction
// per read, so super-nodes whose bytes come from RC-strand reads carry sigs that only the
// RC-of-reference extraction reproduces.
//
// Schemas:
//   sig_ref_map.bin = [u64 N][N x PosRec { u64 sig, u32 chrom_id, u32 pos }]   (sorted by sig)
//   chroms.tsv      = lines "chrom_id\tname\tlength"
//
// Forward records have pos = core.start (in chrom coordinates).
// RC records have pos = chrom_len - core.end (in chrom forward coordinates).
#include "graph.h"
#include "svp.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

extern "C" {
#include "encoding.h"
}

namespace {

struct PosRec { uint64_t sig; uint32_t chrom_id; uint32_t pos_start; uint32_t pos_end; uint32_t _pad; };  // 28 B

char COMP[256];
bool COMP_INIT = false;
void init_comp() {
    if (COMP_INIT) return;
    for (int i = 0; i < 256; i++) COMP[i] = 'N';
    COMP[(uint8_t)'A']='T'; COMP[(uint8_t)'T']='A';
    COMP[(uint8_t)'C']='G'; COMP[(uint8_t)'G']='C'; COMP[(uint8_t)'N']='N';
    COMP_INIT = true;
}

uint32_t find_chrom_id(uint64_t global_pos, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end) {
    auto it = std::upper_bound(rec_start.begin(), rec_start.end(), global_pos);
    if (it == rec_start.begin()) { std::fprintf(stderr, "[FATAL] posmap: pos %lu before any record\n", (unsigned long)global_pos); std::abort(); }
    const uint32_t k = (uint32_t)(std::distance(rec_start.begin(), it) - 1);
    if (global_pos >= rec_end[k]) { std::fprintf(stderr, "[FATAL] posmap: pos %lu in barrier after chrom %u\n", (unsigned long)global_pos, k); std::abort(); }
    return k;
}

// pos_start = left edge of the sig flank window (clamped to chrom start)
// pos_end   = right edge of the sig flank window (clamped to chrom end)
// These define the byte range that this sig "owns" — anchor's content extent.
void build_records_fwd(const std::vector<Core>& cores, const std::vector<uint64_t>& sigs,
                       const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end,
                       int flank, std::vector<PosRec>& out) {
    const size_t base = out.size();
    out.resize(base + cores.size());
    const uint64_t mu = (uint64_t)flank;
    for (size_t i = 0; i < cores.size(); i++) {
        const uint32_t c = find_chrom_id(cores[i].start, rec_start, rec_end);
        const uint64_t cs = rec_start[c], ce = rec_end[c];
        const uint64_t cw_start = cores[i].start > mu ? cores[i].start - mu : cs;
        const uint64_t cw_end   = cores[i].end + mu < ce ? cores[i].end + mu : ce;
        const uint64_t ps = cw_start - cs;
        const uint64_t pe = cw_end - cs;
        if (pe > UINT32_MAX) { std::fprintf(stderr, "[FATAL] posmap: pos %lu exceeds u32\n", (unsigned long)pe); std::abort(); }
        out[base + i] = {sigs[i], c, (uint32_t)ps, (uint32_t)pe, 0};
    }
}

void rc_chrom(const std::string& seq, uint64_t s, uint64_t e, std::string& out) {
    out.assign(e - s, 'N');
    for (uint64_t i = 0; i < e - s; i++) out[e - s - 1 - i] = COMP[(uint8_t)seq[s + i]];
}

void run_chrom_rc(uint32_t chrom_id, const std::string& seq, uint64_t s, uint64_t e, int level, int flank, std::vector<PosRec>& out) {
    std::string rc; rc_chrom(seq, s, e, rc);
    std::vector<uint64_t> rs = {0}, re = {rc.size()};
    std::vector<Core> cores; extract_cores(rc, rs, re, level, cores);
    std::vector<uint64_t> sigs; compute_sigs(rc, cores, flank, 0, rs, re, sigs);
    const uint64_t L = rc.size();
    const uint64_t mu = (uint64_t)flank;
    const size_t base = out.size();
    out.resize(base + cores.size());
    for (size_t i = 0; i < cores.size(); i++) {
        // Sig window in RC coords: [start - flank, end + flank). FWD coord of this window:
        //   fwd_start = L - (rc_end + flank) clamped at 0
        //   fwd_end   = L - (rc_start - flank) clamped at L
        const uint64_t rc_w_start = cores[i].start > mu ? cores[i].start - mu : 0;
        const uint64_t rc_w_end   = cores[i].end + mu < L ? cores[i].end + mu : L;
        const uint64_t fwd_start  = L - rc_w_end;
        const uint64_t fwd_end    = L - rc_w_start;
        if (fwd_end > UINT32_MAX) { std::fprintf(stderr, "[FATAL] posmap: rc fwd_end %lu exceeds u32\n", (unsigned long)fwd_end); std::abort(); }
        out[base + i] = {sigs[i], chrom_id, (uint32_t)fwd_start, (uint32_t)fwd_end, 0};
    }
}

void write_sig_ref_map(const std::string& dir, const std::vector<PosRec>& records) {
    const std::string path = dir + "/sig_ref_map.bin";
    const std::string tmp = path + ".tmp";
    FILE* f = xopen(tmp, "wb");
    const uint64_t N = records.size();
    wr(f, &N, 8);
    wr(f, records.data(), N * sizeof(PosRec));
    atomic_close(f, tmp, path);
}

void write_chroms_tsv(const std::string& dir, const std::vector<std::string>& names, const std::vector<uint64_t>& lengths) {
    const std::string path = dir + "/chroms.tsv";
    const std::string tmp = path + ".tmp";
    FILE* f = xopen(tmp, "wb");
    for (size_t i = 0; i < names.size(); i++) {
        std::fprintf(f, "%zu\t%s\t%lu\n", i, names[i].c_str(), (unsigned long)lengths[i]);
    }
    atomic_close(f, tmp, path);
}

uint64_t count_unique_sigs(const std::vector<PosRec>& records) {
    if (records.empty()) return 0;
    uint64_t u = 1;
    for (size_t i = 1; i < records.size(); i++) if (records[i].sig != records[i - 1].sig) u++;
    return u;
}

}  // namespace

void run_posmap(const std::string& dir, const std::string& ref_path, int level, int flank) {
    LCP_INIT();
    init_comp();
    std::string seq;
    std::vector<uint64_t> rec_start, rec_end;
    std::vector<std::string> names;
    std::fprintf(stderr, "[posmap] reading %s ...\n", ref_path.c_str());
    read_fasta_flat_named(ref_path, flank, seq, rec_start, rec_end, names);
    std::fprintf(stderr, "[posmap] read %lu chroms, %.2f Gb total\n", (unsigned long)names.size(), seq.size() / 1e9);
    std::vector<uint64_t> lengths(names.size());
    for (size_t i = 0; i < names.size(); i++) lengths[i] = rec_end[i] - rec_start[i];

    std::vector<PosRec> records;
    {
        std::vector<Core> cores;
        std::fprintf(stderr, "[posmap] forward: extracting cores L%d ...\n", level);
        extract_cores(seq, rec_start, rec_end, level, cores);
        std::fprintf(stderr, "[posmap] forward: %lu cores\n", (unsigned long)cores.size());
        std::vector<uint64_t> sigs;
        compute_sigs(seq, cores, flank, 0, rec_start, rec_end, sigs);
        build_records_fwd(cores, sigs, rec_start, rec_end, flank, records);
    }

    std::fprintf(stderr, "[posmap] RC pass per chromosome ...\n");
    for (size_t ci = 0; ci < names.size(); ci++) {
        run_chrom_rc((uint32_t)ci, seq, rec_start[ci], rec_end[ci], level, flank, records);
        if (ci % 50 == 0 || ci + 1 == names.size())
            std::fprintf(stderr, "[posmap] RC chrom %zu/%lu records=%lu\n", ci + 1, (unsigned long)names.size(), (unsigned long)records.size());
    }
    std::string().swap(seq);

    std::fprintf(stderr, "[posmap] sorting %lu records by sig ...\n", (unsigned long)records.size());
    std::sort(records.begin(), records.end(), [](const PosRec& a, const PosRec& b) { return a.sig < b.sig; });

    const uint64_t uniq = count_unique_sigs(records);
    std::fprintf(stderr, "[posmap] %lu records / %lu unique sigs (%.2f%% collisions)\n",
                 (unsigned long)records.size(), (unsigned long)uniq,
                 100.0 * (records.size() - uniq) / std::max<uint64_t>(1, records.size()));

    std::fprintf(stderr, "[posmap] writing sig_ref_map.bin + chroms.tsv ...\n");
    write_sig_ref_map(dir, records);
    write_chroms_tsv(dir, names, lengths);
    std::fprintf(stderr, "[posmap] done: %lu records / %lu chroms -> %s/sig_ref_map.bin\n",
                 (unsigned long)records.size(), (unsigned long)names.size(), dir.c_str());
}
