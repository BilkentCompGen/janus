// anchor: per-bubble (chrom, s1_left, s1_right, s2_left, s2_right) resolution. Each Mode B
// bubble has two color=3 anchor super-nodes (s1, s2). Their constituent original-graph node
// IDs map to counts.bin sigs, and posmap binds each sig to a (chrom, pos_start, pos_end) in
// the reference. Majority-vote across the anchor's sigs gives the consensus chrom; we then
// record the leftmost and rightmost flank-window extents of BOTH anchors on the consensus
// chrom, so vcf.cpp can compute strand-aware variant POS via:
//     strand=0:  POS = s1_left + pre + 1
//     strand=1:  POS = s2_right - suf
// (verified against truth SNPs: 91.2% exact POS match vs 1.4% with the old pos_start+1 rule.)
//
// bubble_pos.bin schema (52 B per record):
//   u64 B
//   B x BubblePos { u32 bubble_id, u32 chrom_id,
//                   u32 s1_left, u32 s1_right, u32 s2_left, u32 s2_right,
//                   u32 n_anchor_s1, u32 n_resolved_s1,
//                   u32 n_anchor_s2, u32 n_resolved_s2,
//                   u8  strand, u8 _pad[3] }
// chrom_id == UINT32_MAX means unanchored (consensus failed or no hits).
// strand = 1 when s2_left < s1_right on the consensus chrom (s2 sits before s1 in forward,
//   meaning the arm reads as RC of the forward sequence).
#include "graph.h"
#include "svp.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>

namespace {

struct SuperNode { uint32_t chain_len, block_len, runs_len; uint8_t color; };
struct CovRun    { uint32_t end_pos, n_cov, t_cov; };
struct PosRec    { uint64_t sig; uint32_t chrom_id; uint32_t pos_start; uint32_t pos_end; uint32_t _pad; };

struct BubbleTag {
    uint32_t s1, s2;
    uint8_t  tag, bucket; uint16_t _pad;
    uint32_t pre, suf;
    uint32_t n_resid_off, n_resid_len, t_resid_off, t_resid_len;
    uint32_t n_runs_off, n_runs_len, t_runs_off, t_runs_len;
};

struct BubblePos {
    uint32_t bubble_id, chrom_id;
    uint32_t s1_left, s1_right, s2_left, s2_right;
    uint32_t n_anchor_s1, n_resolved_s1, n_anchor_s2, n_resolved_s2;
    uint8_t  strand; uint8_t _pad[3];
};

void load_counts_sigs(const std::string& dir, std::vector<uint64_t>& node_sigs) {
    FILE* f = xopen(dir + "/counts.bin", "rb");
    uint64_t K; rd(f, &K, 8);
    std::vector<SigCount> rows(K); rd(f, rows.data(), K * sizeof(SigCount));
    node_sigs.resize(K);
    for (uint64_t i = 0; i < K; i++) node_sigs[i] = rows[i].sig;
    std::fclose(f);
}

void load_modea_chains(const std::string& dir, std::vector<SuperNode>& snodes,
                       std::vector<uint32_t>& chain_ids, std::vector<uint64_t>& chain_off) {
    FILE* f = xopen(dir + "/modea.bin", "rb");
    uint64_t S; rd(f, &S, 8); snodes.resize(S); rd(f, snodes.data(), S * sizeof(SuperNode));
    uint64_t tc; rd(f, &tc, 8); chain_ids.resize(tc); rd(f, chain_ids.data(), tc * 4);
    std::fclose(f);
    chain_off.assign(S + 1, 0);
    for (uint64_t i = 0; i < S; i++) chain_off[i + 1] = chain_off[i] + snodes[i].chain_len;
}

void load_modeb_tags(const std::string& dir, std::vector<BubbleTag>& tags) {
    FILE* f = xopen(dir + "/modeb.bin", "rb");
    uint64_t B; rd(f, &B, 8); tags.resize(B); rd(f, tags.data(), B * sizeof(BubbleTag));
    std::fclose(f);
}

void load_posmap(const std::string& dir, std::vector<PosRec>& records) {
    FILE* f = xopen(dir + "/sig_ref_map.bin", "rb");
    uint64_t N; rd(f, &N, 8); records.resize(N); rd(f, records.data(), N * sizeof(PosRec));
    std::fclose(f);
}

struct Hit { uint32_t chrom_id, pos_start, pos_end; };

void gather_anchor_hits(uint32_t super_id, const std::vector<uint64_t>& chain_off,
                        const std::vector<uint32_t>& chain_ids, const std::vector<uint64_t>& node_sigs,
                        const std::vector<PosRec>& posmap, std::vector<Hit>& hits) {
    hits.clear();
    const uint64_t lo = chain_off[super_id], hi = chain_off[super_id + 1];
    for (uint64_t i = lo; i < hi; i++) {
        const uint64_t sig = node_sigs[chain_ids[i]];
        auto it = std::lower_bound(posmap.begin(), posmap.end(), sig,
                                    [](const PosRec& r, uint64_t s) { return r.sig < s; });
        for (; it != posmap.end() && it->sig == sig; ++it) hits.push_back({it->chrom_id, it->pos_start, it->pos_end});
    }
}

uint32_t majority_chrom(const std::vector<Hit>& hits) {
    if (hits.empty()) return UINT32_MAX;
    std::unordered_map<uint32_t, uint32_t> counts;
    for (auto& h : hits) counts[h.chrom_id]++;
    uint32_t best = UINT32_MAX, best_n = 0;
    for (auto& kv : counts) if (kv.second > best_n) { best = kv.first; best_n = kv.second; }
    if (best_n * 2 < hits.size()) return UINT32_MAX;
    return best;
}

void chrom_range(uint32_t chrom, const std::vector<Hit>& hits,
                 uint32_t& left_edge, uint32_t& right_edge) {
    left_edge = UINT32_MAX; right_edge = 0;
    for (auto& h : hits) if (h.chrom_id == chrom) {
        if (h.pos_start < left_edge) left_edge = h.pos_start;
        if (h.pos_end > right_edge) right_edge = h.pos_end;
    }
}

uint32_t count_resolved_on_chrom(uint32_t chrom, const std::vector<Hit>& hits) {
    if (chrom == UINT32_MAX) return 0;
    uint32_t n = 0;
    for (auto& h : hits) if (h.chrom_id == chrom) n++;
    return n;
}

BubblePos process_bubble(uint32_t bubble_id, const BubbleTag& t,
                         const std::vector<uint64_t>& chain_off, const std::vector<uint32_t>& chain_ids,
                         const std::vector<uint64_t>& node_sigs, const std::vector<PosRec>& posmap) {
    std::vector<Hit> hits_s1, hits_s2;
    gather_anchor_hits(t.s1, chain_off, chain_ids, node_sigs, posmap, hits_s1);
    gather_anchor_hits(t.s2, chain_off, chain_ids, node_sigs, posmap, hits_s2);
    std::vector<Hit> combined;
    combined.insert(combined.end(), hits_s1.begin(), hits_s1.end());
    combined.insert(combined.end(), hits_s2.begin(), hits_s2.end());
    const uint32_t chrom = majority_chrom(combined);
    BubblePos bp{};
    bp.bubble_id = bubble_id;
    bp.chrom_id = chrom;
    bp.n_anchor_s1 = (uint32_t)(chain_off[t.s1 + 1] - chain_off[t.s1]);
    bp.n_anchor_s2 = (uint32_t)(chain_off[t.s2 + 1] - chain_off[t.s2]);
    bp.n_resolved_s1 = count_resolved_on_chrom(chrom, hits_s1);
    bp.n_resolved_s2 = count_resolved_on_chrom(chrom, hits_s2);
    if (chrom == UINT32_MAX) return bp;
    chrom_range(chrom, hits_s1, bp.s1_left, bp.s1_right);
    chrom_range(chrom, hits_s2, bp.s2_left, bp.s2_right);
    // Strand: if s2 starts before s1 ends in forward, this bubble's arm reads as RC of the
    // forward genome (LCP path was extracted from RC-strand reads).
    bp.strand = (bp.n_resolved_s1 > 0 && bp.n_resolved_s2 > 0 && bp.s2_left < bp.s1_right) ? 1 : 0;
    return bp;
}

void write_bubble_pos(const std::string& dir, const std::vector<BubblePos>& records) {
    const std::string path = dir + "/bubble_pos.bin";
    const std::string tmp = path + ".tmp";
    FILE* f = xopen(tmp, "wb");
    const uint64_t B = records.size(); wr(f, &B, 8);
    wr(f, records.data(), B * sizeof(BubblePos));
    atomic_close(f, tmp, path);
}

}  // namespace

void run_anchor(const std::string& dir) {
    std::vector<SuperNode> snodes;
    std::vector<uint32_t> chain_ids;
    std::vector<uint64_t> chain_off;
    std::vector<uint64_t> node_sigs;
    std::vector<BubbleTag> tags;
    std::vector<PosRec> posmap;
    std::fprintf(stderr, "[anchor] loading inputs ...\n");
    load_modea_chains(dir, snodes, chain_ids, chain_off);
    load_counts_sigs(dir, node_sigs);
    load_modeb_tags(dir, tags);
    load_posmap(dir, posmap);
    std::fprintf(stderr, "[anchor] S=%lu chain_ids=%lu nodes=%lu bubbles=%lu posmap=%lu\n",
                 (unsigned long)snodes.size(), (unsigned long)chain_ids.size(),
                 (unsigned long)node_sigs.size(), (unsigned long)tags.size(), (unsigned long)posmap.size());

    std::vector<BubblePos> out(tags.size());
    uint64_t resolved = 0, strict_both = 0, strand1 = 0;
    for (size_t i = 0; i < tags.size(); i++) {
        out[i] = process_bubble((uint32_t)i, tags[i], chain_off, chain_ids, node_sigs, posmap);
        if (out[i].chrom_id != UINT32_MAX) {
            resolved++;
            const uint32_t need_s1 = out[i].n_anchor_s1 / 2 + 1, need_s2 = out[i].n_anchor_s2 / 2 + 1;
            if (out[i].n_resolved_s1 >= need_s1 && out[i].n_resolved_s2 >= need_s2) strict_both++;
            if (out[i].strand == 1) strand1++;
        }
    }
    write_bubble_pos(dir, out);
    std::fprintf(stderr, "[anchor] resolved %lu/%lu (%.1f%%); strict-both %lu (%.1f%%); strand1 %lu -> %s/bubble_pos.bin\n",
                 (unsigned long)resolved, (unsigned long)tags.size(),
                 100.0 * resolved / std::max<uint64_t>(1, tags.size()),
                 (unsigned long)strict_both, 100.0 * strict_both / std::max<uint64_t>(1, tags.size()),
                 (unsigned long)strand1, dir.c_str());
}
