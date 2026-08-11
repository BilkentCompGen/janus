// snape call <ctg.bam> <normal.bam|-> <out.vcf> [min_mapq=10] [min_aln=0] [min_sv=50] [cluster=50] [min_supp=5]
//
// Breakpoints from the aligned junction contigs -- the stage after minimap2, in snape instead of
// Severus. It is the part of Severus that was doing work on this input and nothing else: a MAPQ and
// aligned-length gate, split alignments chained in QUERY order, CIGAR indels at or above the SV
// floor, clustering, and a support floor. What Severus carries besides that is idle here by
// construction -- phasing (one haplotype-free contig), coverage histograms and VAF (a contig is not
// a read pileup), multi-sample subtraction (the somatic decision was already made in the graph).
//
// A junction is two reference sides plus the direction each side leaves the reference in: d=+1 the
// alignment ENDS there going right, d=-1 it STARTS there going left. That pair is the entire type
// vocabulary -- (+,-) deletion, (-,+) duplication, equal signs inversion, two contigs a breakend --
// so the caller never needs to know what an SV is in order to type one.
//
// Support counts DISTINCT contig names, so a contig that reaches the same junction twice (once
// through its CIGAR, once through a split) still supports it once, and `snape junctions` naming its
// records <read>_<n> keeps that equal to one supporting READ.
//
// Positions are 1-based from the moment a junction is built (a 0-based exclusive reference end is
// already the 1-based coordinate of the last aligned base; a 0-based start needs the +1), so the VCF
// writer converts nothing.
#include "call.h"

#include "unordered_dense.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <htslib/sam.h>

namespace {

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[call] FATAL: %s\n", m); std::exit(1); }

struct Seg {                                   // one surviving alignment of one contig
    uint32_t q;                                // contig id
    int32_t tid;
    hts_pos_t rs, re;                          // reference span, 0-based half-open
    int32_t qs, qe;                            // query span in the CONTIG's own orientation
    uint8_t rev;
    double rep_s, rep_e;                       // how repetitive the anchor is at each end
};

struct Jn {
    int32_t c1, c2;
    hts_pos_t p1, p2;                          // 1-based
    int8_t d1, d2;
    uint32_t q;
    double rep;                                // the more repetitive of the two anchors
};

struct Ins { int32_t c; hts_pos_t p; int32_t len; uint32_t q; double rep; };

struct Stats {
    uint64_t rec = 0, unmapped = 0, secondary = 0, lowq = 0, shortaln = 0, highmm = 0, germ = 0;
    uint64_t seg = 0, split = 0, cigar_del = 0, cigar_ins = 0;
};

double qident(const uint8_t* s, int32_t a, int32_t b, int32_t len) {
    int m = 0;
    for (int32_t i = 0; i < len; i++) m += bam_seqi(s, a + i) == bam_seqi(s, b + i);
    return (double)m / len;
}

// The largest share of any base agreeing with the base `p` before it, over periods 1..6: how much of
// this stretch is a tandem repeat. Random sequence sits near 0.25 (one base in four agrees by
// chance) and a VNTR runs far above it.
double repeat_score(const uint8_t* sq, int32_t a, int32_t b) {
    const int32_t n = b - a;
    if (n < 12) return 1.0;                       // too short to judge: treat as repeat, i.e. refuse
    double best = 0;
    for (int p = 1; p <= 6 && p < n; p++) {
        int m = 0;
        for (int32_t k = a + p; k < b; k++) m += bam_seqi(sq, k) == bam_seqi(sq, k - p);
        best = std::max(best, (double)m / (n - p));
    }
    return best;
}

// The same measurement at a junction, taken on the ANCHOR side only. Those bases are aligned to the
// reference at under 1% mismatch, so they ARE the reference at that locus and no FASTA is needed.
// Measured on hcc1954 against the 2-support ensemble, over the calls nothing else corroborates:
// 27.2% of them sit above 0.4 against 2.9% of the true calls, and above 0.6 it is 11.0% against
// 0.2%. A breakpoint whose flank is a VNTR is where an aligner puts an indel that is not an event.
constexpr int32_t REP_WIN = 300;

double rep_before(const bam1_t* b, int32_t q) {
    const int32_t a = std::max(0, q - REP_WIN);
    return q > a ? repeat_score(bam_get_seq(b), a, q) : 1.0;
}

double rep_after(const bam1_t* b, int32_t q) {
    const int32_t e = std::min((int32_t)b->core.l_qseq, q + REP_WIN);
    return e > q ? repeat_score(bam_get_seq(b), q, e) : 1.0;
}

// A short tandem duplication does not reach the aligner as a split: minimap2 writes the second copy
// as an insertion and the event silently becomes an INS. The copy it duplicates is right there in the
// contig, on one side or the other of the inserted bases, and those flanking bases are aligned -- so
// they ARE the reference sequence at that locus. Comparing the two needs no reference of our own.
// Which side matches says which copy the aligner kept, and that fixes the duplicated interval.
bool tandem(const bam1_t* b, int32_t qp, int len, hts_pos_t rp, int32_t tid, uint32_t q,
            std::vector<Jn>& jn) {
    const uint8_t* sq = bam_get_seq(b);
    const double aft = qp + 2 * len <= b->core.l_qseq ? qident(sq, qp, qp + len, len) : 0;
    const double bef = qp - len >= 0 ? qident(sq, qp, qp - len, len) : 0;
    if (std::max(aft, bef) < 0.85) return false;
    const double rep = rep_before(b, qp);
    if (aft >= bef) jn.push_back({tid, tid, rp + 1, rp + 1 + len, -1, 1, q, rep});
    else            jn.push_back({tid, tid, rp - len + 1, rp + 1, -1, 1, q, rep});
    return true;
}

// An inversion shorter than the contig does not split-align either: the aligner walks straight
// through it and pays for the inverted bases with a dense run of tiny indels. Measured on `hard`, a
// 194 bp inversion produced the same ~15-operation shred in every contig carrying it, and not one
// operation of it reached the 50 bp floor. A reference window that shredded is a breakpoint pair
// whatever the aligner called the pieces; the type is read from the directions like any other.
//
// On real HiFi the same rule fires on tandem repeats, where a dense run of small indels is what a
// VNTR looks like and not what an inversion looks like: measured on hcc1954 it produced 41 false
// calls and not one true one, while on `hard` it is worth 4 inversions in the 50-1000 bp band. The
// difference is visible in the contig's own bases without a reference, so the window is asked whether
// it is repetitive before it is allowed to be a breakpoint.
void shred(const std::vector<std::pair<hts_pos_t, int32_t>>& ops, const bam1_t* b, int32_t tid,
           uint32_t q, int min_sv, int nmin, hts_pos_t win, double max_rep, std::vector<Jn>& jn) {
    const uint8_t* sq = bam_get_seq(b);
    for (size_t i = 0; i < ops.size();) {
        size_t j = i;
        while (j + 1 < ops.size() && ops[j + 1].first - ops[i].first <= win) j++;
        if ((int)(j - i + 1) >= nmin && ops[j].first - ops[i].first >= min_sv &&
            repeat_score(sq, ops[i].second, ops[j].second) <= max_rep) {
            jn.push_back({tid, tid, ops[i].first, ops[j].first, 1, 1, q, rep_before(b, ops[i].second)});
            i = j + 1;
        } else {
            i++;
        }
    }
}

// One CIGAR walk: the segment's query offsets and reference end, plus every indel big enough to be
// an SV on its own. Leading/trailing clips are query positions the aligner did not place, and they
// are what makes two segments of the same contig orderable.
void walk(const bam1_t* b, uint32_t q, int min_sv, Seg& s, std::vector<Jn>& jn, std::vector<Ins>& ins,
          int nmin, hts_pos_t win, double max_rep) {
    const uint32_t* cg = bam_get_cigar(b);
    hts_pos_t rp = b->core.pos;
    int32_t lead = 0, qaln = 0, qtot = 0, qp = 0;
    std::vector<std::pair<hts_pos_t, int32_t>> small;      // reference position, contig position
    for (uint32_t i = 0; i < b->core.n_cigar; i++) {
        const int op = bam_cigar_op(cg[i]), len = (int)bam_cigar_oplen(cg[i]);
        const int t = bam_cigar_type(op);
        const bool clip = (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP);
        const bool indel = (op == BAM_CDEL || op == BAM_CINS);
        if (op == BAM_CDEL && len >= min_sv)
            jn.push_back({s.tid, s.tid, rp, rp + len, 1, -1, q, std::max(rep_before(b, qp), rep_after(b, qp))});
        if (op == BAM_CINS && len >= min_sv && !tandem(b, qp, len, rp, s.tid, q, jn))
            ins.push_back({s.tid, rp, len, q, std::max(rep_before(b, qp), rep_after(b, qp + len))});
        if (indel && len < min_sv && len >= 2) small.push_back({rp, qp});
        if (clip && qaln == 0) lead += len;
        qtot += (clip || (t & 1)) ? len : 0;
        qaln += (!clip && (t & 1)) ? len : 0;
        qp += (t & 1) ? len : 0;
        rp += (t & 2) ? len : 0;
    }
    shred(small, b, s.tid, q, min_sv, nmin, win, max_rep, jn);
    s.re = rp;
    s.qs = b->core.flag & BAM_FREVERSE ? qtot - lead - qaln : lead;
    s.qe = s.qs + qaln;
    // The bases the aligner placed, at each end of what it placed: the anchor a split leaves from.
    s.rep_s = rep_after(b, s.qs);
    s.rep_e = rep_before(b, s.qe);
}

// Two consecutive pieces of one contig: where each leaves the reference, ordered so that the same
// junction seen from either side lands in the same cluster.
Jn chain(const Seg& a, const Seg& b) {
    Jn j{a.tid, b.tid, a.rev ? a.rs + 1 : a.re, b.rev ? b.re : b.rs + 1,
         a.rev ? (int8_t)-1 : (int8_t)1, b.rev ? (int8_t)1 : (int8_t)-1, a.q,
         std::max(a.rep_e, b.rep_s)};
    if (j.c1 > j.c2 || (j.c1 == j.c2 && j.p1 > j.p2)) {
        std::swap(j.c1, j.c2); std::swap(j.p1, j.p2); std::swap(j.d1, j.d2);
    }
    return j;
}

// An alignment carrying far more mismatch than the reads can contain is not this locus. Measured on
// hcc1954 real HiFi (read error 0.018%): the alignments that carry a TRUE call sit at 0.08% mismatch
// across their germline anchors, the ones that carry a FALSE call at 8.35%, and 89% of the false ones
// are above 1%. The cause is that 47% of those contigs have no match anywhere in hg38 better than 2%
// -- the copy they came from is missing from the reference or is one of a dozen -- and MAPQ cannot
// warn about it, because with a single candidate the aligner is confident and simply wrong.
// SUBSTITUTIONS only. NM counts inserted and deleted bases as well, so using it raw would score a
// contig that carries a real 400 bp deletion at 8% and throw away exactly the events we are looking
// for -- measured: it took the simulation from 542/600 junctions to 310. Subtracting the indel bases
// leaves the quantity that actually says "this is not the locus I came from".
double mm_rate(const bam1_t* b) {
    const uint8_t* nm = bam_aux_get(b, "NM");
    if (!nm) return 0;
    const uint32_t* cg = bam_get_cigar(b);
    int64_t aln = 0, indel = 0;
    for (uint32_t i = 0; i < b->core.n_cigar; i++) {
        const int op = bam_cigar_op(cg[i]), len = (int)bam_cigar_oplen(cg[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) aln += len;
        else if (op == BAM_CINS || op == BAM_CDEL) indel += len;
    }
    const int64_t sub = (int64_t)bam_aux2i(nm) - indel;
    return aln > 0 && sub > 0 ? (double)sub / (double)aln : 0;
}

// min_aln defaults to 0 because a length floor on a SEGMENT is the wrong shape. A deletion needs one
// surviving alignment of a contig and a breakend needs two, so a breakend pays the floor twice -- and
// worse, dropping a contig's MIDDLE segment chains the two that remain, which fabricates a junction
// between pieces that were never adjacent. Measured on hcc1954 with everything else held: 1000 bp
// gives 0.749/0.718, 600 gives 0.761/0.762, 400 gives 0.771/0.797, and 200 gives 0.798/0.863 -- both
// precision AND recall improve as the floor comes off, and false positives fall from 266 to 241.
// Below 200 the numbers are identical to three digits at 150, 100, 50 and 0, because minimap2 does not
// emit anything shorter here; on `hard` 0 costs nothing and recovers 2 duplications.
bool keep(const bam1_t* b, int min_mapq, hts_pos_t min_aln, double max_mm, Stats& st) {
    st.rec++;
    if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP)) { st.unmapped++; return false; }
    if (b->core.flag & BAM_FSECONDARY) { st.secondary++; return false; }
    if (b->core.qual < min_mapq) { st.lowq++; return false; }
    if (bam_endpos(b) - b->core.pos < min_aln) { st.shortaln++; return false; }
    if (max_mm > 0 && mm_rate(b) > max_mm) { st.highmm++; return false; }
    return true;
}

// ---- the normal, read at the reference locus ----
//
// The graph decides somatic BEFORE alignment, by exact core sequence, and in a repeat that is not the
// same question as somatic. Measured on hcc1954: of the calls this caller gets wrong in the 50-1000 bp
// band, 93.1% have an indel of the same size sitting in the NORMAL reads at the same locus (against
// 2.9% of the calls it gets right, size-matched), at a normal allele fraction of 0.75. They are
// germline. The tumour and normal reads simply do not spell that locus with the same 217 bp core --
// different stutter, a different paralogous copy sampled -- so ncov == 0 is true of the STRING while
// the VARIANT is shared. Nothing before alignment can see that; this can.
struct Normal {
    samFile* fp = nullptr;
    hts_idx_t* idx = nullptr;
    sam_hdr_t* hdr = nullptr;
    bam1_t* b = nullptr;
    bool on() const { return fp != nullptr; }
};

Normal open_normal(const char* path) {
    Normal n;
    if (!path || !std::strcmp(path, "-")) return n;
    n.fp = sam_open(path, "r");
    if (!n.fp) die("open normal bam");
    n.hdr = sam_hdr_read(n.fp);
    n.idx = sam_index_load(n.fp, path);
    if (!n.hdr || !n.idx) die("normal bam needs an index");
    n.b = bam_init1();
    return n;
}

// Same event in the normal, or a locus the normal cannot place either: both mean not somatic.
//
// The second clause is why this also runs on breakends, which it did not before. Measured on
// hcc1954 against the 2-support ensemble: of the calls nothing corroborates, 7.4% sit at a locus
// where a quarter or more of the normal reads are above 2% mismatch, against 0.0% of the true calls
// -- and 27 of those 136 were breakends, which the same-chromosome guard was letting past
// unexamined.
bool germline(Normal& n, const char* chrom, hts_pos_t pos, hts_pos_t len, hts_pos_t win, double dirty_f) {
    if (!n.on()) return false;
    const int tid = sam_hdr_name2tid(n.hdr, chrom);
    if (tid < 0) return false;
    hts_itr_t* it = sam_itr_queryi(n.idx, tid, (int)std::max<hts_pos_t>(0, pos - win), (int)(pos + win));
    if (!it) return false;
    int seen = 0, dirty = 0; bool same = false;
    while (sam_itr_next(n.fp, it, n.b) >= 0) {
        if (n.b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY)) continue;
        if (n.b->core.qual < 10) continue;
        seen++;
        dirty += mm_rate(n.b) > 0.02;
        const uint32_t* cg = bam_get_cigar(n.b);
        hts_pos_t rp = n.b->core.pos;
        for (uint32_t i = 0; i < n.b->core.n_cigar && !same; i++) {
            const int op = bam_cigar_op(cg[i]); const hts_pos_t L = bam_cigar_oplen(cg[i]);
            if ((op == BAM_CINS || op == BAM_CDEL) && len > 0 &&
                L * 2 >= len && L <= len * 2 && std::llabs(rp - pos) <= win) same = true;
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
                op == BAM_CEQUAL || op == BAM_CDIFF) rp += L;
        }
    }
    hts_itr_destroy(it);
    return same || (seen >= 4 && (double)dirty >= dirty_f * seen);
}

// Read the BAM once. Coordinate order scatters a contig's alignments, so segments are bucketed per
// contig name and chained only after the file is done.
void collect(const char* bam, int min_mapq, hts_pos_t min_aln, int min_sv, std::vector<Jn>& jn,
             std::vector<Ins>& ins, sam_hdr_t** hdr, Stats& st, int nmin, hts_pos_t win, double max_mm,
             double max_rep, std::vector<std::string>& qname) {
    samFile* fp = sam_open(bam, "r");
    if (!fp) die("open bam");
    *hdr = sam_hdr_read(fp);
    if (!*hdr) die("read bam header");
    bam1_t* b = bam_init1();
    ankerl::unordered_dense::map<std::string, uint32_t> qid;
    std::vector<std::vector<Seg>> byq;
    while (sam_read1(fp, *hdr, b) >= 0) {
        if (!keep(b, min_mapq, min_aln, max_mm, st)) continue;
        const char* nm = bam_get_qname(b);
        auto it = qid.find(nm);
        if (it == qid.end()) {
            it = qid.emplace(nm, (uint32_t)byq.size()).first; byq.emplace_back(); qname.emplace_back(nm);
        }
        Seg s{it->second, b->core.tid, b->core.pos, 0, 0, 0,
              (uint8_t)((b->core.flag & BAM_FREVERSE) != 0), 1.0, 1.0};
        const size_t before = jn.size();
        walk(b, it->second, min_sv, s, jn, ins, nmin, win, max_rep);
        st.cigar_del += jn.size() - before;
        byq[it->second].push_back(s);
        st.seg++;
    }
    bam_destroy1(b);
    sam_close(fp);
    st.cigar_ins = ins.size();
    for (auto& v : byq) {
        std::sort(v.begin(), v.end(), [](const Seg& x, const Seg& y) { return x.qs < y.qs; });
        for (size_t i = 1; i < v.size(); i++) jn.push_back(chain(v[i - 1], v[i]));
        st.split += v.size() ? v.size() - 1 : 0;
    }
}

struct Cl { Jn j; hts_pos_t p1, p2; uint32_t supp; double rep; };

uint32_t distinct(std::vector<uint32_t>& q) {
    std::sort(q.begin(), q.end());
    return (uint32_t)(std::unique(q.begin(), q.end()) - q.begin());
}

// Junctions within `tol` on BOTH sides, with the same contig pair and the same direction pair, are
// one event. The cluster's position is the median, which cannot be dragged by a single ragged end.
std::vector<Cl> cluster(std::vector<Jn>& jn, hts_pos_t tol, uint32_t min_supp) {
    std::sort(jn.begin(), jn.end(), [](const Jn& a, const Jn& b) {
        return std::tie(a.c1, a.c2, a.d1, a.d2, a.p1, a.p2) < std::tie(b.c1, b.c2, b.d1, b.d2, b.p1, b.p2);
    });
    std::vector<Cl> out;
    for (size_t i = 0; i < jn.size();) {
        size_t k = i + 1;
        while (k < jn.size() && jn[k].c1 == jn[i].c1 && jn[k].c2 == jn[i].c2 && jn[k].d1 == jn[i].d1 &&
               jn[k].d2 == jn[i].d2 && jn[k].p1 - jn[i].p1 <= tol && std::llabs(jn[k].p2 - jn[i].p2) <= tol) k++;
        std::vector<uint32_t> q; std::vector<hts_pos_t> a, c; std::vector<double> rp;
        for (size_t m = i; m < k; m++) {
            q.push_back(jn[m].q); a.push_back(jn[m].p1); c.push_back(jn[m].p2); rp.push_back(jn[m].rep);
        }
        const uint32_t s = distinct(q);
        std::sort(a.begin(), a.end()); std::sort(c.begin(), c.end()); std::sort(rp.begin(), rp.end());
        if (s >= min_supp) out.push_back({jn[i], a[a.size() / 2], c[c.size() / 2], s, rp[rp.size() / 2]});
        i = k;
    }
    return out;
}

struct ClIns { Ins v; int32_t len; uint32_t supp; double rep; };

std::vector<ClIns> cluster_ins(std::vector<Ins>& ins, hts_pos_t tol, uint32_t min_supp) {
    std::sort(ins.begin(), ins.end(), [](const Ins& a, const Ins& b) {
        return std::tie(a.c, a.p) < std::tie(b.c, b.p);
    });
    std::vector<ClIns> out;
    for (size_t i = 0; i < ins.size();) {
        size_t k = i + 1;
        while (k < ins.size() && ins[k].c == ins[i].c && ins[k].p - ins[i].p <= tol) k++;
        std::vector<uint32_t> q; std::vector<int32_t> L; std::vector<double> rp;
        for (size_t m = i; m < k; m++) { q.push_back(ins[m].q); L.push_back(ins[m].len); rp.push_back(ins[m].rep); }
        const uint32_t s = distinct(q);
        std::sort(L.begin(), L.end()); std::sort(rp.begin(), rp.end());
        if (s >= min_supp) out.push_back({ins[i], L[L.size() / 2], s, rp[rp.size() / 2]});
        i = k;
    }
    return out;
}

const char* svtype(const Cl& c) {
    if (c.j.c1 != c.j.c2) return "BND";
    if (c.j.d1 == c.j.d2) return "INV";
    return c.j.d1 > 0 ? "DEL" : "DUP";
}

// A breakend names its partner in the ALT field; the bracket says which side of the partner the
// sequence continues on, and whether N leads or trails says the same for this side.
std::string bnd_alt(const Cl& c, sam_hdr_t* h) {
    const std::string m = std::string(sam_hdr_tid2name(h, c.j.c2)) + ":" + std::to_string(c.p2);
    const std::string br = c.j.d2 < 0 ? "[" + m + "[" : "]" + m + "]";
    return c.j.d1 > 0 ? "N" + br : br + "N";
}

void write_vcf(const char* path, sam_hdr_t* h, const std::vector<Cl>& cl, const std::vector<ClIns>& ci,
               const std::vector<std::string>& qname) {
    FILE* f = std::fopen(path, "w");
    if (!f) die("open out.vcf");
    std::fprintf(f, "##fileformat=VCFv4.2\n##source=snape call\n");
    for (int i = 0; i < sam_hdr_nref(h); i++)
        std::fprintf(f, "##contig=<ID=%s,length=%" PRIhts_pos ">\n", sam_hdr_tid2name(h, i), sam_hdr_tid2len(h, i));
    std::fprintf(f, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n"
                    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n"
                    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Second breakpoint\">\n"
                    "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Direction each side leaves the reference\">\n"
                    "##INFO=<ID=SUPP,Number=1,Type=Integer,Description=\"Supporting contigs\">\n"
                    "##INFO=<ID=REP,Number=1,Type=Float,Description=\"Tandem-repeat fraction of the anchor\">\n"
                    "##INFO=<ID=CTG,Number=1,Type=String,Description=\"One junction contig supporting this call; `snape pick` reads the normal at ITS anchors\">\n"
                    "##ALT=<ID=DEL,Description=\"Deletion\">\n##ALT=<ID=DUP,Description=\"Duplication\">\n"
                    "##ALT=<ID=INV,Description=\"Inversion\">\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tt\n");
    int n = 0;
    for (const Cl& c : cl) {
        const char* t = svtype(c);
        const std::string alt = c.j.c1 == c.j.c2 ? "<" + std::string(t) + ">" : bnd_alt(c, h);
        std::fprintf(f, "%s\t%" PRIhts_pos "\tsvp_%s%d\tN\t%s\t.\tPASS\tSVTYPE=%s;", sam_hdr_tid2name(h, c.j.c1),
                     c.p1, t, n++, alt.c_str(), t);
        if (c.j.c1 == c.j.c2)
            std::fprintf(f, "SVLEN=%" PRIhts_pos ";END=%" PRIhts_pos ";", c.p2 - c.p1, c.p2);
        std::fprintf(f, "STRANDS=%c%c;SUPP=%u;REP=%.3f;CTG=%s\tGT:DV\t0/1:%u\n",
                     c.j.d1 > 0 ? '+' : '-', c.j.d2 > 0 ? '+' : '-', c.supp, c.rep,
                     qname[c.j.q].c_str(), c.supp);
    }
    for (const ClIns& c : ci)
        std::fprintf(f, "%s\t%" PRIhts_pos "\tsvp_INS%d\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=%d;SUPP=%u;CTG=%s\tGT:DV\t0/1:%u\n",
                     sam_hdr_tid2name(h, c.v.c), c.v.p, n++, c.len, c.supp, qname[c.v.q].c_str(), c.supp);
    std::fclose(f);
}

}  // namespace

int cmd_call(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: snape call <ctg.bam> <normal.bam|-> <out.vcf> [min_mapq=10] [min_aln=0] "
                             "[min_sv=50] [cluster=50] [min_supp=5] [shred_n=6] [shred_win=300] [max_mm=0.01] "
                             "[shred_rep=0.40] [anchor_rep=0.50] [normal_dirty=0.25]\n");
        return 2;
    }
    const int min_mapq   = argc >  4 ? std::atoi(argv[4])  : 10;
    const hts_pos_t min_aln = argc > 5 ? std::atoll(argv[5]) : 0;
    const int min_sv     = argc >  6 ? std::atoi(argv[6])  : 50;
    const hts_pos_t tol  = argc >  7 ? std::atoll(argv[7]) : 50;
    // 5, swept on all four pairs through the whole shipped pipeline (F1): 3 is worse than 4
    // everywhere, and between 4 and 5 the split is two pairs each -- but 4 wins hcc1954 by 0.002 and
    // hcc1937 by 0.003 while 5 wins h2009 by 0.014 and h1437 by 0.028, so 5 has the better mean
    // (0.808 against 0.799). 6 gives back more recall than it buys.
    const uint32_t min_supp = argc > 8 ? (uint32_t)std::atoi(argv[8]) : 5;
    const int shred_n    = argc >  9 ? std::atoi(argv[9])  : 6;
    const hts_pos_t shred_win = argc > 10 ? std::atoll(argv[10]) : 300;
    const double max_mm  = argc > 11 ? std::atof(argv[11]) : 0.01;
    const double shred_rep = argc > 12 ? std::atof(argv[12]) : 0.40;
    // 0.50 and not 0.45: the knee is between 0.45 and 0.40, where 33 true calls go at once, and the
    // simulation pays 12 of 600 junctions for 0.50 already -- `hard` has no VNTR-driven false calls
    // to win back, so every step down is pure loss there and pure gain on real reads.
    const double max_rep   = argc > 13 ? std::atof(argv[13]) : 0.50;
    const double dirty_f   = argc > 14 ? std::atof(argv[14]) : 0.25;

    std::vector<Jn> jn; std::vector<Ins> ins; sam_hdr_t* hdr = nullptr; Stats st;
    std::vector<std::string> qname;
    collect(argv[1], min_mapq, min_aln, min_sv, jn, ins, &hdr, st, shred_n, shred_win, max_mm, shred_rep, qname);
    // A same-contig junction shorter than the SV floor is a local alignment artefact, not an event.
    jn.erase(std::remove_if(jn.begin(), jn.end(), [&](const Jn& j) {
        return j.c1 == j.c2 && j.p2 - j.p1 < min_sv; }), jn.end());
    std::vector<Cl> cl = cluster(jn, tol, min_supp);
    std::vector<ClIns> ci = cluster_ins(ins, tol, min_supp);

    // The normal, at the reference locus: the one question the graph could not ask.
    // A breakpoint whose anchor is a tandem repeat is where an aligner puts an indel that is not an
    // event. Applied to the CALL, not just inside shred(), because the same population reaches the
    // caller through CIGAR indels and splits as well.
    if (max_rep < 1.0) {
        const size_t bj = cl.size(), bi = ci.size();
        cl.erase(std::remove_if(cl.begin(), cl.end(), [&](const Cl& c) { return c.rep > max_rep; }), cl.end());
        ci.erase(std::remove_if(ci.begin(), ci.end(), [&](const ClIns& c) { return c.rep > max_rep; }), ci.end());
        std::fprintf(stderr, "[call] cipasi tekrarli oldugu icin elenen: %zu (%zu + %zu kaldi)\n",
                     (bj - cl.size()) + (bi - ci.size()), cl.size(), ci.size());
    }

    // The normal, at the reference locus: the one question the graph could not ask. A breakend is
    // checked at BOTH of its ends, and it has no length to match, so only the placement clause bites.
    Normal nrm = open_normal(argv[2]);
    if (nrm.on()) {
        const size_t before_j = cl.size(), before_i = ci.size();
        cl.erase(std::remove_if(cl.begin(), cl.end(), [&](const Cl& c) {
            const char* n1 = sam_hdr_tid2name(hdr, c.j.c1);
            if (c.j.c1 == c.j.c2) return germline(nrm, n1, c.p1, c.p2 - c.p1, tol * 10, dirty_f);
            return germline(nrm, n1, c.p1, 0, tol * 10, dirty_f) ||
                   germline(nrm, sam_hdr_tid2name(hdr, c.j.c2), c.p2, 0, tol * 10, dirty_f); }), cl.end());
        ci.erase(std::remove_if(ci.begin(), ci.end(), [&](const ClIns& c) {
            return germline(nrm, sam_hdr_tid2name(hdr, c.v.c), c.v.p, c.len, tol * 10, dirty_f); }), ci.end());
        st.germ = (before_j - cl.size()) + (before_i - ci.size());
        std::fprintf(stderr, "[call] normalde de bulundugu icin elenen: %llu (%zu junction + %zu insertion kaldi)\n",
                     (unsigned long long)st.germ, cl.size(), ci.size());
    }
    write_vcf(argv[3], hdr, cl, ci, qname);

    std::fprintf(stderr, "[call] kayit %llu | atlanan: unmapped %llu, secondary %llu, MAPQ<%d %llu, <%lldbp %llu, mismatch>%.0f%% %llu\n",
                 (unsigned long long)st.rec, (unsigned long long)st.unmapped, (unsigned long long)st.secondary,
                 min_mapq, (unsigned long long)st.lowq, (long long)min_aln, (unsigned long long)st.shortaln,
                 100*max_mm, (unsigned long long)st.highmm);
    std::fprintf(stderr, "[call] segment %llu -> split %llu + CIGAR DEL %llu | INS %llu\n",
                 (unsigned long long)st.seg, (unsigned long long)st.split,
                 (unsigned long long)st.cigar_del, (unsigned long long)st.cigar_ins);
    std::fprintf(stderr, "[call] aday %zu junction + %zu insertion -> supp>=%u: %zu + %zu -> %s\n",
                 jn.size(), ins.size(), min_supp, cl.size(), ci.size(), argv[3]);
    sam_hdr_destroy(hdr);
    return 0;
}
