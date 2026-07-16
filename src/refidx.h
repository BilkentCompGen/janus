// refidx -- ONE reference index answering BOTH questions svp asks about a core.
//
// For every core of every read the caller needs to know:
//
//   1. Which strand family is it in?  FW-only / RC-only / BOTH / NEITHER.
//      FW-only and RC-only vote on the read's orientation. BOTH means the core sits in an inverted
//      repeat -- the same sequence occurs in the reference forwards AND backwards -- so it carries no
//      strand information and must ABSTAIN. (Scoring BOTH as FW is precisely the bug that made
//      satellites and palindromes look like they misoriented 13% of their reads; with abstention the
//      rate is 1 in 320,000 on normal reads.) NEITHER means the core is not reference sequence at all:
//      a read error, or an SV junction.
//
//   2. Is it a REPEAT core?  The blacklist: a core occurring more than once in the reference, which
//      graph.cpp merges into its flanks instead of making a node.
//
// These used to be two structures and two lookups. They are one, because the blacklist is by
// construction a SUBSET of the reference cores: tools/build_blacklist parses the reference with the
// same level, the same margin, both strands, owning each locus once, and hashes with the same
// core_sig. Every blacklisted signature is therefore already a key here. It costs a third flag bit.
//
// LAYOUT -- shaped by where the time actually goes.
//
// The first version used a 2^24-bucket offset array (67 MB) plus 6-byte entries (522 MB). Measured:
// 130 ns/lookup, against 92 ns for the two hash sets it replaced. SLOWER, despite touching a third of
// the memory. The cause was not TLB pressure (huge pages were tried and made it 0.96x -- worse) but
// two DEPENDENT cache misses: off[bucket] must be read before ent[i] can be addressed, so the
// latencies serialise.
//
// Fix: make the offset array fit in L3 so the first access stops being a miss.
//
//     bucket    = sig >> 42                2^22 buckets -> offsets are 16.8 MB, resident in L3 (32 MB)
//     entry     = fp29 << 3 | flags        4 bytes, and 20.6 entries per bucket
//     flags     = bit0 FW, bit1 RC, bit2 REPEAT
//
//     86.4M entries x 4 B  = 346 MB
//     4.2M offsets  x 4 B  =  17 MB
//                            ------
//                            363 MB, one real cache miss, one structure
//
// The 40-bit key of v1 was overkill. All the key has to do is keep false positives rare, and a 29-bit
// fingerprint over a 20.6-entry bucket gives P(false match) = 20.6 / 2^29 = 3.8e-8. A whole run makes
// ~460M lookups on error-broken cores (17% of 2.7e9), so that is ~18 spurious votes GENOME-WIDE.
// Collisions between two REAL signatures (same bucket, same fingerprint -- 51 discriminating bits) run
// to about 3 pairs across 86.4M cores; they merge flags and end up abstaining, which is the safe way
// to be wrong.
//
// A Bloom filter would have been smaller still and was rejected: 17% of a read's cores are error-broken
// and land in exactly the bucket a false positive would corrupt, turning noise into a vote.
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace refidx {

constexpr uint64_t MAGIC = 0x3249444649455255ULL;   // "UREIFDI2"
constexpr int      BUCKET_BITS = 22;
constexpr uint64_t NBUCKET = 1ULL << BUCKET_BITS;
constexpr int      SHIFT = 64 - BUCKET_BITS;        // 42
constexpr int      FP_BITS = 29;
constexpr uint32_t FP_MASK = (1u << FP_BITS) - 1;
constexpr int      FLAG_BITS = 3;

enum Flag : uint8_t { FW = 1, RC = 2, REPEAT = 4 };

// what the strand vote sees (the REPEAT bit is asked for separately)
enum Class : uint8_t { NEITHER = 0, FW_ONLY = 1, RC_ONLY = 2, BOTH = 3 };

struct Header {
  uint64_t magic;
  uint64_t level;      // LCP level the index was built at -- a mismatch is a silent disaster
  uint64_t margin;
  uint64_t nentry;
  uint64_t nfw, nrc, nboth, nrepeat;
};

inline uint64_t bucket_of(uint64_t sig) { return sig >> SHIFT; }
// The fingerprint takes the TOP 29 of the 42 bits below the bucket. Signatures are MurmurHash64A
// outputs, so every bit is already uniform: no further mixing is needed or useful.
inline uint32_t fp_of(uint64_t sig) { return (uint32_t)((sig >> (SHIFT - FP_BITS)) & FP_MASK); }
inline uint32_t pack(uint32_t fp, uint8_t fl) { return (fp << FLAG_BITS) | (fl & 7); }

// Read-only view of the index, mmapped straight off disk. Opening it is a syscall, not a parse: the
// hash-set cache this replaces took 4.5 seconds to load on every single run.
struct Index {
  int fd = -1;
  const uint8_t* base = nullptr;
  size_t len = 0;
  Header hdr{};
  const uint32_t* off = nullptr;    // NBUCKET+1 offsets into the entry array
  const uint32_t* ent = nullptr;    // nentry entries, sorted by fingerprint within each bucket

  bool open(const std::string& path, int level) {
    fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) { std::fprintf(stderr, "[refidx] cannot open %s\n", path.c_str()); return false; }
    struct stat st;
    if (fstat(fd, &st) != 0) return false;
    len = (size_t)st.st_size;
    void* m = mmap(nullptr, len, PROT_READ, MAP_PRIVATE, fd, 0);
    if (m == MAP_FAILED) { std::fprintf(stderr, "[refidx] mmap failed\n"); return false; }
    base = (const uint8_t*)m;
    std::memcpy(&hdr, base, sizeof(Header));
    if (hdr.magic != MAGIC) { std::fprintf(stderr, "[refidx] bad magic in %s\n", path.c_str()); return false; }
    if ((int)hdr.level != level) {
      std::fprintf(stderr, "[refidx] level mismatch: index=%d, asked=%d\n", (int)hdr.level, level);
      return false;
    }
    off = (const uint32_t*)(base + sizeof(Header));
    ent = (const uint32_t*)(base + sizeof(Header) + (NBUCKET + 1) * sizeof(uint32_t));
    // Warm the offset array into L3 -- it is 17 MB and every lookup starts there. Leaving it to
    // demand-faulting means the first million lookups each eat a page fault AND a miss.
    { volatile uint32_t s = 0;
      for (uint64_t i = 0; i <= NBUCKET; i += 16) s += off[i];
      (void)s; }
    // the entries are touched at random; suppress read-ahead
    madvise((void*)ent, len - (size_t)((const uint8_t*)ent - base), MADV_RANDOM);
    return true;
  }
  void close() {
    if (base) { munmap((void*)base, len); base = nullptr; }
    if (fd >= 0) { ::close(fd); fd = -1; }
  }

  // THE hot path: one L3 hit on off[], one miss on ent[]. Entries are sorted by fingerprint, so the
  // scan bails out early on a miss instead of walking the whole bucket.
  inline uint8_t flags(uint64_t sig) const {
    const uint64_t b = sig >> SHIFT;
    const uint32_t fp = fp_of(sig);
    const uint32_t lo = off[b], hi = off[b + 1];
    for (uint32_t i = lo; i < hi; i++) {
      const uint32_t e = ent[i];
      const uint32_t f = e >> FLAG_BITS;
      if (f == fp) return (uint8_t)(e & 7);
      if (f > fp) break;
    }
    return 0;
  }
  inline Class classify(uint64_t sig) const { return (Class)(flags(sig) & 3); }
  inline bool  is_repeat(uint64_t sig) const { return (flags(sig) & REPEAT) != 0; }
};

}  // namespace refidx

// ============================================================================================
// What svp actually calls (src/refidx.cpp). ONE index open, ONE lookup per core, TWO answers:
//
//   the STRAND VOTE   -- is this read the reference's forward copy, or its reverse complement?
//   the REPEAT MERGE  -- is this core a reference repeat, to be absorbed into its flanks?
//
// They used to be unrelated: the merge came from a separate blacklist file and the graph simply
// inserted BOTH strands of every read because it could not tell which one was real. That doubled the
// node count, and it hid every inversion (a read's own reverse complement was always in the graph
// beside it). The vote replaces it: parse the read forward, ask the index which family its cores fall
// in, and insert the ONE strand that wins. Below VOTE_PCT the read is undecided and both strands go in
// (the safe fallback -- measured on CHM13 at 30x: 99.965% of normal reads decide, 1 in 320,000 decides
// WRONGLY, 1.0004 core-sets per read instead of 2.0).
//
// Cores in BOTH families ABSTAIN. They sit in an inverted repeat, so the same bytes occur in the
// reference forwards and backwards and they carry no strand information at all. Scoring them as FW is
// what made palindromes and satellites look like they misoriented 13% of their reads; with abstention
// that rate is 0.000%. Cores in NEITHER family abstain too -- they are read errors or SV junctions.
//
// KNOWN AND ACCEPTED: a read lying wholly INSIDE a large inversion votes "reverse-complement" (it is,
// with respect to the reference) and is inserted flipped -- landing exactly on the normal reads' nodes.
// The inversion's INTERIOR therefore reads as shared, not somatic. That is not a regression: the old
// dual-strand build lost it too (normal's RC mirror was already in the graph). What is callable about a
// large inversion is its two BREAKPOINT junctions, and those survive -- the reads crossing them split
// their vote, come out undecided, and go in both ways.
// ============================================================================================

// One effective core after the repeat merge: the byte span, and the signature that names its node.
// sig is authoritative -- for an unmerged core it is extract()'s own sig (no re-hash), for a merged
// (A+X)/(X+B) span it is core_sig over the merged bytes.
struct Span { uint32_t st, en; uint64_t sig; };

void     ref_load(const char* path, int level);   // fail-fast on a missing file or a level mismatch
bool     ref_on();                                // false = no index: legacy dual-strand, no merge
uint64_t ref_dropped();                           // merged spans dropped for exceeding len_cap
void     ref_read_stats(uint64_t& fwd, uint64_t& rev, uint64_t& undecided);
uint64_t ref_flipped();                           // reads whose vote SVP_VOTE_ERR deliberately flipped

// Flags for the margin-kept cores of one parse. fl[i] belongs to keep[i]. This is THE lookup: every
// later question about these cores is answered from fl, never from a second pass over the index.
void ref_flags(const uint64_t* sig, const std::vector<int>& keep, std::vector<uint8_t>& fl);

// The strand vote over one parse's flags: +1 the bytes as given are reference-forward, -1 they are the
// reverse complement, 0 undecided (caller inserts both strands). A pure function of the read's bytes and
// the index, so the count pass and the edge pass reach the same verdict -- only the count pass tallies.
int ref_vote(const std::vector<uint8_t>& fl, bool tally);

// (A+X)->(X+B): a maximal run of REPEAT cores is absorbed into its non-repeat flanks as two OVERLAPPING
// merged cores, so A->X->B becomes A->B and the bare repeat core never becomes a node. With no index
// loaded this is the identity (every kept core passes through). A merged span longer than len_cap is
// dropped and the chain breaks there.
void ref_merge(const char* data, const uint64_t* st, const uint64_t* en, const uint64_t* sig,
               const std::vector<int>& keep, const std::vector<uint8_t>& fl, uint64_t len_cap,
               std::vector<Span>& out);
