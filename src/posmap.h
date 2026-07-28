#pragma once
#include <cstdint>
#include <string>

// The position map ("reference linked list"): sig -> every (contig, pos, strand) the core occupies in
// the reference, KEEPING repeats (multi-locus) and BOTH strands. Built merge-free by tools/refpos so a
// posaware graph build (which runs ref_merge OFF) sees byte-identical cores. Records are 17 B, packed,
// sorted by (sig, contig, pos); file = u64 N then Rec[N]. Only consumed under SVP_POSAWARE.
namespace posmap {

#pragma pack(push, 1)
struct Rec { uint64_t sig; uint32_t contig; uint32_t pos; uint8_t strand; };   // 17 bytes
#pragma pack(pop)

constexpr int      BUCKET_BITS = 24;
constexpr uint64_t NBUCKET = 1ULL << BUCKET_BITS;
constexpr int      SHIFT = 64 - BUCKET_BITS;                 // 40

// Packed node locus: strand(1) | contig(23) | pos(40). NOVEL = a core with no reference occurrence
// (tumour-specific sequence / SV-junction-spanning core / read error) -> kept, keyed by sig alone.
// UNRESOLVED = a repeat core whose copy could NOT be determined from the read's anchors. Splitting it
// on a guess is what scatters a repeat's edges across copies and TANGLES the graph (measured: splitting
// on the smallest-locus guess tripled deg>=8 vs merge-free). Collapsing every unresolved instance of a
// sig onto this one sentinel keeps them a single node -- exactly the merge-free behaviour -- so only
// CONFIDENT placements split. Both sentinels are outside the packed range (strand bit + all-ones pos).
constexpr uint64_t NOVEL = UINT64_MAX;
constexpr uint64_t UNRESOLVED = UINT64_MAX - 1;
inline uint64_t pack_locus(uint32_t contig, uint32_t pos, uint8_t strand) {
    return ((uint64_t)(strand & 1) << 63) | ((uint64_t)(contig & 0x7FFFFF) << 40) | ((uint64_t)pos & 0xFFFFFFFFFFULL);
}
inline uint8_t  locus_strand(uint64_t L) { return (uint8_t)(L >> 63); }
inline uint32_t locus_contig(uint64_t L) { return (uint32_t)((L >> 40) & 0x7FFFFF); }
inline uint64_t locus_pos(uint64_t L)    { return L & 0xFFFFFFFFFFULL; }

// Read-only, mmapped straight off disk (like refidx::Index). A per-sig lookup returns the contiguous
// [first, first+count) slice of Recs; a bucket table over sig>>SHIFT bounds the search to ~one cache line.
struct Map {
    int fd = -1;
    const uint8_t* base = nullptr;
    size_t len = 0;
    uint64_t N = 0;
    const Rec* recs = nullptr;
    uint32_t* bstart = nullptr;                             // NBUCKET+1 offsets into recs, by sig>>SHIFT

    bool open(const std::string& path);
    void close();

    // Return the first Rec with .sig==sig and set count = how many share it (0 if none). Recs are
    // sorted, so the slice is contiguous. Misaligned 8-byte reads (packed 17-B recs) are fine on x86.
    inline const Rec* find(uint64_t sig, uint32_t& count) const {
        const uint64_t b = sig >> SHIFT;
        uint32_t i = bstart[b]; const uint32_t hi = bstart[b + 1];
        while (i < hi && recs[i].sig < sig) i++;
        const uint32_t lo = i;
        while (i < hi && recs[i].sig == sig) i++;
        count = i - lo;
        return recs + lo;
    }
};

}  // namespace posmap
