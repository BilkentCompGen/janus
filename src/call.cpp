#include "cmd.h"
#include "log.h"
#include "unordered_dense.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <htslib/sam.h>

static std::vector<Ctg> parse_ctg(const char* arg) {
    std::vector<Ctg> v;
    const std::string in(arg);
    for (size_t i = 0; i <= in.size();) {
        const size_t c = std::min(in.find(',', i), in.size());
        const std::string f = in.substr(i, c - i);
        if (f.empty()) die("call", "empty contig-set entry");
        const size_t e = f.find('=');
        if (e == std::string::npos) v.push_back({"t", f});
        else v.push_back({f.substr(0, e), f.substr(e + 1)});
        if (v.back().name.empty() || v.back().bam.empty()) die("call", "contig set needs <name>=<bam>");
        i = c + 1;
    }
    if (v.size() > 64) die("call", "at most 64 contig sets");
    for (size_t i = 1; i < v.size(); i++)
        for (size_t j = 0; j < i; j++)
            if (v[i].name == v[j].name) die("call", "duplicate contig-set name");
    return v;
}

static double qident(const uint8_t* s, int32_t a, int32_t b, int32_t len) {
    int m = 0;
    for (int32_t i = 0; i < len; i++) m += bam_seqi(s, a + i) == bam_seqi(s, b + i);
    return (double)m / len;
}

static double repeat_score(const uint8_t* sq, int32_t a, int32_t b) {
    const int32_t n = b - a;
    if (n < 12) return 1.0;
    double best = 0;
    for (int p = 1; p <= 6 && p < n; p++) {
        int m = 0;
        for (int32_t k = a + p; k < b; k++) m += bam_seqi(sq, k) == bam_seqi(sq, k - p);
        best = std::max(best, (double)m / (n - p));
    }
    return best;
}

constexpr int32_t REP_WIN = 300;

static double rep_before(const bam1_t* b, int32_t q) {
    const int32_t a = std::max(0, q - REP_WIN);
    return q > a ? repeat_score(bam_get_seq(b), a, q) : 1.0;
}

static double rep_after(const bam1_t* b, int32_t q) {
    const int32_t e = std::min((int32_t)b->core.l_qseq, q + REP_WIN);
    return e > q ? repeat_score(bam_get_seq(b), q, e) : 1.0;
}

static bool tandem(const bam1_t* b, int32_t qp, int len, hts_pos_t rp, int32_t tid, uint32_t q,
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

static void shred(const std::vector<std::pair<hts_pos_t, int32_t>>& ops, const bam1_t* b, int32_t tid,
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

static void walk(const bam1_t* b, uint32_t q, int min_sv, Seg& s, std::vector<Jn>& jn, std::vector<Ins>& ins,
                 int nmin, hts_pos_t win, double max_rep) {
    const uint32_t* cg = bam_get_cigar(b);
    hts_pos_t rp = b->core.pos;
    int32_t lead = 0, qaln = 0, qtot = 0, qp = 0;
    std::vector<std::pair<hts_pos_t, int32_t>> small;
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

    s.rep_s = rep_after(b, s.qs);
    s.rep_e = rep_before(b, s.qe);
}

static Jn chain(const Seg& a, const Seg& b) {
    Jn j{a.tid, b.tid, a.rev ? a.rs + 1 : a.re, b.rev ? b.re : b.rs + 1,
         a.rev ? (int8_t)-1 : (int8_t)1, b.rev ? (int8_t)1 : (int8_t)-1, a.q,
         std::max(a.rep_e, b.rep_s)};
    if (j.c1 > j.c2 || (j.c1 == j.c2 && j.p1 > j.p2)) {
        std::swap(j.c1, j.c2);
        std::swap(j.p1, j.p2);
        std::swap(j.d1, j.d2);
    }
    return j;
}

static double mm_rate(const bam1_t* b) {
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

static bool keep(const bam1_t* b, int min_mapq, hts_pos_t min_aln, double max_mm, Stats& st) {
    st.rec++;
    if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP)) { st.unmapped++; return false; }
    if (b->core.flag & BAM_FSECONDARY) { st.secondary++; return false; }
    if (b->core.qual < min_mapq) { st.lowq++; return false; }
    if (bam_endpos(b) - b->core.pos < min_aln) { st.shortaln++; return false; }
    if (max_mm > 0 && mm_rate(b) > max_mm) { st.highmm++; return false; }
    return true;
}

static Normal open_normal(const char* path) {
    Normal n;
    if (!path || !std::strcmp(path, "-")) return n;
    n.fp = sam_open(path, "r");
    if (!n.fp) die("call", "open normal bam");
    n.hdr = sam_hdr_read(n.fp);
    n.idx = sam_index_load(n.fp, path);
    if (!n.hdr || !n.idx) die("call", "normal bam needs an index");
    n.b = bam_init1();
    return n;
}

static bool germline(Normal& n, const char* chrom, hts_pos_t pos, hts_pos_t len, hts_pos_t win, double dirty_f) {
    if (!n.fp) return false;
    const int tid = sam_hdr_name2tid(n.hdr, chrom);
    if (tid < 0) return false;
    hts_itr_t* it = sam_itr_queryi(n.idx, tid, (int)std::max<hts_pos_t>(0, pos - win), (int)(pos + win));
    if (!it) return false;
    int seen = 0, dirty = 0;
    bool same = false;
    while (sam_itr_next(n.fp, it, n.b) >= 0) {
        if (n.b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY)) continue;
        if (n.b->core.qual < 10) continue;
        seen++;
        dirty += mm_rate(n.b) > 0.02;
        const uint32_t* cg = bam_get_cigar(n.b);
        hts_pos_t rp = n.b->core.pos;
        for (uint32_t i = 0; i < n.b->core.n_cigar && !same; i++) {
            const int op = bam_cigar_op(cg[i]);
            const hts_pos_t L = bam_cigar_oplen(cg[i]);
            if ((op == BAM_CINS || op == BAM_CDEL) && len > 0 &&
                L * 2 >= len && L <= len * 2 && std::llabs(rp - pos) <= win) same = true;
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
                op == BAM_CEQUAL || op == BAM_CDIFF) rp += L;
        }
    }
    hts_itr_destroy(it);
    return same || (seen >= 4 && (double)dirty >= dirty_f * seen);
}

static void same_refs(sam_hdr_t* a, sam_hdr_t* b) {
    bool ok = sam_hdr_nref(a) == sam_hdr_nref(b);
    for (int i = 0; ok && i < sam_hdr_nref(a); i++)
        ok = std::strcmp(sam_hdr_tid2name(a, i), sam_hdr_tid2name(b, i)) == 0 &&
             sam_hdr_tid2len(a, i) == sam_hdr_tid2len(b, i);
    if (!ok) die("call", "contig BAMs were aligned to different references");
}

static bool seg_cmp(const Seg& x, const Seg& y) { return x.qs < y.qs; }

static bool resumes(const Seg& a, const Seg& c, hts_pos_t min_sv, hts_pos_t& pos) {
    if (a.tid != c.tid || a.rev != c.rev) return false;
    pos = a.rev ? a.rs + 1 : a.re;
    return std::llabs((long long)(a.rev ? a.rs - c.re : c.rs - a.re)) < min_sv;
}

static void collect(const char* bam, uint8_t line, int min_mapq, hts_pos_t min_aln, int min_sv,
                    std::vector<Jn>& jn, std::vector<Ins>& ins, sam_hdr_t** hdr, Stats& st, int nmin,
                    hts_pos_t win, double max_mm, double max_rep, std::vector<std::string>& qname,
                    std::vector<uint8_t>& qline) {
    samFile* fp = sam_open(bam, "r");
    if (!fp) die("call", "open bam");
    sam_hdr_t* h = sam_hdr_read(fp);
    if (!h) die("call", "read bam header");
    if (!*hdr) *hdr = h;
    else same_refs(*hdr, h);
    bam1_t* b = bam_init1();
    ankerl::unordered_dense::map<std::string, uint32_t> qid;
    std::vector<std::vector<Seg>> byq;
    const uint32_t base = (uint32_t)qname.size();
    while (sam_read1(fp, h, b) >= 0) {
        if (!keep(b, min_mapq, min_aln, max_mm, st)) continue;
        const char* nm = bam_get_qname(b);
        auto it = qid.find(nm);
        if (it == qid.end()) {
            it = qid.emplace(nm, (uint32_t)byq.size()).first;
            byq.emplace_back();
            qname.emplace_back(nm);
            qline.push_back(line);
        }
        Seg s{base + it->second, b->core.tid, b->core.pos, 0, 0, 0,
              (uint8_t)((b->core.flag & BAM_FREVERSE) != 0), 1.0, 1.0};
        const size_t before = jn.size();
        walk(b, base + it->second, min_sv, s, jn, ins, nmin, win, max_rep);
        st.cigar_del += jn.size() - before;
        byq[it->second].push_back(s);
        st.seg++;
    }
    bam_destroy1(b);
    sam_close(fp);
    if (h != *hdr) sam_hdr_destroy(h);
    st.cigar_ins = ins.size();

    for (auto& v : byq) {
        std::sort(v.begin(), v.end(), seg_cmp);
        st.split += v.size() ? v.size() - 1 : 0;
        size_t i = 1;
        while (i < v.size()) {
            hts_pos_t pos;
            if (i + 1 < v.size() && resumes(v[i - 1], v[i + 1], min_sv, pos) &&
                v[i + 1].qs - v[i - 1].qe >= min_sv) {
                ins.push_back({v[i - 1].tid, pos, v[i + 1].qs - v[i - 1].qe, v[i - 1].q,
                               std::max(v[i - 1].rep_e, v[i + 1].rep_s)});
                st.ins_split++;
                i += 2;
                continue;
            }
            if (resumes(v[i - 1], v[i], min_sv, pos) && v[i].qs - v[i - 1].qe >= min_sv) {
                ins.push_back({v[i - 1].tid, pos, v[i].qs - v[i - 1].qe, v[i - 1].q,
                               std::max(v[i - 1].rep_e, v[i].rep_s)});
                st.ins_split++;
                i += 1;
                continue;
            }
            jn.push_back(chain(v[i - 1], v[i]));
            i += 1;
        }
    }
}

static uint32_t distinct(std::vector<uint32_t>& q, const std::vector<uint8_t>& qline, uint32_t* cnt, int L) {
    std::sort(q.begin(), q.end());
    const size_t n = (size_t)(std::unique(q.begin(), q.end()) - q.begin());
    for (int l = 0; l < L; l++) cnt[l] = 0;
    for (size_t i = 0; i < n; i++) cnt[qline[q[i]]]++;
    return (uint32_t)n;
}

static uint32_t best_of(const uint32_t* cnt, int L) {
    uint32_t b = 0;
    for (int l = 0; l < L; l++) if (cnt[l] > b) b = cnt[l];
    return b;
}

static bool jn_cmp(const Jn& a, const Jn& b) {
    if (a.c1 != b.c1) return a.c1 < b.c1;
    if (a.c2 != b.c2) return a.c2 < b.c2;
    if (a.d1 != b.d1) return a.d1 < b.d1;
    if (a.d2 != b.d2) return a.d2 < b.d2;
    if (a.p1 != b.p1) return a.p1 < b.p1;
    return a.p2 < b.p2;
}

static std::vector<Cl> cluster(std::vector<Jn>& jn, const std::vector<uint8_t>& qline, hts_pos_t tol,
                               uint32_t min_supp, int L, std::vector<uint32_t>& sup) {
    std::sort(jn.begin(), jn.end(), jn_cmp);
    std::vector<Cl> out;
    std::vector<uint32_t> cnt((size_t)L);
    for (size_t i = 0; i < jn.size();) {
        size_t k = i + 1;
        while (k < jn.size() && jn[k].c1 == jn[i].c1 && jn[k].c2 == jn[i].c2 && jn[k].d1 == jn[i].d1 &&
               jn[k].d2 == jn[i].d2 && jn[k].p1 - jn[i].p1 <= tol && std::llabs(jn[k].p2 - jn[i].p2) <= tol) k++;
        std::vector<uint32_t> q;
        std::vector<hts_pos_t> a, c;
        std::vector<double> rp;
        for (size_t m = i; m < k; m++) {
            q.push_back(jn[m].q);
            a.push_back(jn[m].p1);
            c.push_back(jn[m].p2);
            rp.push_back(jn[m].rep);
        }
        const uint32_t s = distinct(q, qline, cnt.data(), L);
        std::sort(a.begin(), a.end());
        std::sort(c.begin(), c.end());
        std::sort(rp.begin(), rp.end());

        if (best_of(cnt.data(), L) >= min_supp) {
            out.push_back({jn[i], a[a.size() / 2], c[c.size() / 2], s,
                           (uint32_t)(sup.size() / (size_t)L), rp[rp.size() / 2]});
            sup.insert(sup.end(), cnt.begin(), cnt.end());
        }
        i = k;
    }
    return out;
}

static bool ins_cmp(const Ins& a, const Ins& b) {
    if (a.c != b.c) return a.c < b.c;
    return a.p < b.p;
}

static std::vector<ClIns> cluster_ins(std::vector<Ins>& ins, const std::vector<uint8_t>& qline, hts_pos_t tol,
                                      uint32_t min_supp, int L, std::vector<uint32_t>& sup) {
    std::sort(ins.begin(), ins.end(), ins_cmp);
    std::vector<ClIns> out;
    std::vector<uint32_t> cnt((size_t)L);
    for (size_t i = 0; i < ins.size();) {
        size_t k = i + 1;
        while (k < ins.size() && ins[k].c == ins[i].c && ins[k].p - ins[i].p <= tol) k++;
        std::vector<uint32_t> q;
        std::vector<int32_t> Ln;
        std::vector<double> rp;
        for (size_t m = i; m < k; m++) { q.push_back(ins[m].q); Ln.push_back(ins[m].len); rp.push_back(ins[m].rep); }
        const uint32_t s = distinct(q, qline, cnt.data(), L);
        std::sort(Ln.begin(), Ln.end());
        std::sort(rp.begin(), rp.end());
        if (best_of(cnt.data(), L) >= min_supp) {
            out.push_back({ins[i], Ln[Ln.size() / 2], s, (uint32_t)(sup.size() / (size_t)L),
                           rp[rp.size() / 2]});
            sup.insert(sup.end(), cnt.begin(), cnt.end());
        }
        i = k;
    }
    return out;
}

static const char* svtype(const Cl& c) {
    if (c.j.c1 != c.j.c2) return "BND";
    if (c.j.d1 == c.j.d2) return "INV";
    return c.j.d1 > 0 ? "DEL" : "DUP";
}

static std::string bnd_alt(const Cl& c, sam_hdr_t* h) {
    const std::string m = std::string(sam_hdr_tid2name(h, c.j.c2)) + ":" + std::to_string(c.p2);
    const std::string br = c.j.d2 < 0 ? "[" + m + "[" : "]" + m + "]";
    return c.j.d1 > 0 ? "N" + br : br + "N";
}

static void fmt_cols(FILE* f, const std::vector<uint32_t>& sup, uint32_t slot, int L, uint32_t assign, int only) {
    for (int l = 0; l < L; l++) {
        if (only >= 0 && l != only) continue;
        const uint32_t d = sup[(size_t)slot * (size_t)L + (size_t)l];
        std::fprintf(f, "\t%s:%u", d >= assign ? "0/1" : "0/0", d);
    }
    std::fputc('\n', f);
}

// idea 1: report only the 24 primary chromosomes; drop GL*/KI* alt-contig calls
static bool is_primary(const char* n) {
    if (!n || std::strncmp(n, "chr", 3) != 0) return false;
    const char* d = n + 3;
    if ((*d == 'X' || *d == 'Y' || *d == 'M') && d[1] == '\0') return true;
    if (*d >= '1' && *d <= '9') {
        int v = 0;
        while (*d >= '0' && *d <= '9') v = v * 10 + (*d++ - '0');
        return *d == '\0' && v <= 22;
    }
    return false;
}

// idea 3: mask (flag, never drop) calls inside the hg38 pseudoautosomal regions
static const struct {
    const char* n;
    long s, e;
} kPar[] = {
    {"chrX", 10002, 2800000}, {"chrX", 155700000, 156040895},
    {"chrY", 10002, 2781474},  {"chrY", 56887143, 57217401},
};
static bool in_par(const char* n, hts_pos_t p) {
    for (const auto& r : kPar)
        if (std::strcmp(r.n, n) == 0 && p >= r.s && p <= r.e) return true;
    return false;
}

static void write_vcf(const char* path, sam_hdr_t* h, const std::vector<Cl>& cl, const std::vector<ClIns>& ci,
                      const std::vector<std::string>& qname, const std::vector<Ctg>& sets,
                      const std::vector<uint32_t>& sup, int L, uint32_t assign, int only) {
    const std::string tmp = std::string(path) + ".tmp";
    FILE* f = xopen(tmp.c_str(), "w");
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
                    "##ALT=<ID=INV,Description=\"Inversion\">\n##ALT=<ID=INS,Description=\"Insertion\">\n"
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Supporting contigs from this line\">\n");
    std::fprintf(f, "##INFO=<ID=PAR,Number=0,Type=Flag,Description=\"Call lies within a pseudoautosomal region (masked, not dropped)\">\n");
    std::fprintf(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (int l = 0; l < L; l++) if (only < 0 || l == only) std::fprintf(f, "\t%s", sets[(size_t)l].name.c_str());
    std::fputc('\n', f);
    int n = 0;
    for (const Cl& c : cl) {
        if (only >= 0 && sup[(size_t)c.slot * (size_t)L + (size_t)only] < assign) { n++; continue; }
        const char* nm1 = sam_hdr_tid2name(h, c.j.c1);
        if (!is_primary(nm1)) { n++; continue; }
        bool par = in_par(nm1, c.p1);
        if (c.j.c1 != c.j.c2) {
            const char* nm2 = sam_hdr_tid2name(h, c.j.c2);
            if (!is_primary(nm2)) { n++; continue; }
            par = par || in_par(nm2, c.p2);
        }
        const char* t = svtype(c);
        const std::string alt = c.j.c1 == c.j.c2 ? "<" + std::string(t) + ">" : bnd_alt(c, h);
        std::fprintf(f, "%s\t%" PRIhts_pos "\tsvp_%s%d\tN\t%s\t.\tPASS\tSVTYPE=%s;", nm1,
                     c.p1, t, n++, alt.c_str(), t);
        if (c.j.c1 == c.j.c2)
            std::fprintf(f, "SVLEN=%" PRIhts_pos ";END=%" PRIhts_pos ";", c.p2 - c.p1, c.p2);
        std::fprintf(f, "STRANDS=%c%c;SUPP=%u;REP=%.3f;CTG=%s",
                     c.j.d1 > 0 ? '+' : '-', c.j.d2 > 0 ? '+' : '-', c.supp, c.rep,
                     qname[c.j.q].c_str());
        if (par) std::fprintf(f, ";PAR");
        std::fprintf(f, "\tGT:DV");
        fmt_cols(f, sup, c.slot, L, assign, only);
    }
    for (const ClIns& c : ci) {
        if (only >= 0 && sup[(size_t)c.slot * (size_t)L + (size_t)only] < assign) { n++; continue; }
        const char* nmi = sam_hdr_tid2name(h, c.v.c);
        if (!is_primary(nmi)) { n++; continue; }
        std::fprintf(f, "%s\t%" PRIhts_pos "\tsvp_INS%d\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=%d;SUPP=%u;CTG=%s",
                     nmi, c.v.p, n++, c.len, c.supp, qname[c.v.q].c_str());
        if (in_par(nmi, c.v.p)) std::fprintf(f, ";PAR");
        std::fprintf(f, "\tGT:DV");
        fmt_cols(f, sup, c.slot, L, assign, only);
    }
    atomic_close(f, tmp.c_str(), path);
}

int run_call(const std::vector<Ctg>& sets, const char* normal, const std::string& out, const CallArgs& a) {
    const int L = (int)sets.size();
    std::vector<Jn> jn;
    std::vector<Ins> ins;
    sam_hdr_t* hdr = nullptr;
    Stats st;
    std::vector<std::string> qname;
    std::vector<uint8_t> qline;
    for (int l = 0; l < L; l++)
        collect(sets[(size_t)l].bam.c_str(), (uint8_t)l, a.min_mapq, a.min_aln, a.min_sv, jn, ins, &hdr, st,
                a.shred_n, a.shred_win, a.max_mm, a.shred_rep, qname, qline);

    {
        size_t k = 0;
        for (size_t i = 0; i < jn.size(); i++) {
            const Jn& j = jn[i];
            if (!(j.c1 == j.c2 && j.p2 - j.p1 < a.min_sv)) jn[k++] = j;
        }
        jn.resize(k);
    }
    std::vector<uint32_t> sup;
    std::vector<Cl> cl = cluster(jn, qline, a.tol, a.min_supp, L, sup);
    std::vector<ClIns> ci = cluster_ins(ins, qline, a.tol, a.min_supp, L, sup);

    if (a.max_rep < 1.0) {
        const size_t bj = cl.size(), bi = ci.size();
        {
            size_t k = 0;
            for (size_t i = 0; i < cl.size(); i++) if (cl[i].rep <= a.max_rep) cl[k++] = cl[i];
            cl.resize(k);
        }
        {
            size_t k = 0;
            for (size_t i = 0; i < ci.size(); i++) if (ci[i].rep <= a.max_rep) ci[k++] = ci[i];
            ci.resize(k);
        }
        say("call", "dropped as repetitive: %zu (%zu + %zu left)",
            (bj - cl.size()) + (bi - ci.size()), cl.size(), ci.size());
    }

    Normal nrm = open_normal(normal);
    if (nrm.fp) {
        const size_t before_j = cl.size(), before_i = ci.size();
        {
            size_t k = 0;
            for (size_t i = 0; i < cl.size(); i++) {
                const Cl& c = cl[i];
                const char* n1 = sam_hdr_tid2name(hdr, c.j.c1);
                bool drop;
                if (c.j.c1 == c.j.c2) drop = germline(nrm, n1, c.p1, c.p2 - c.p1, a.tol * 10, a.dirty_f);
                else drop = germline(nrm, n1, c.p1, 0, a.tol * 10, a.dirty_f) ||
                            germline(nrm, sam_hdr_tid2name(hdr, c.j.c2), c.p2, 0, a.tol * 10, a.dirty_f);
                if (!drop) cl[k++] = c;
            }
            cl.resize(k);
        }
        {
            size_t k = 0;
            for (size_t i = 0; i < ci.size(); i++) {
                const ClIns& c = ci[i];
                if (!germline(nrm, sam_hdr_tid2name(hdr, c.v.c), c.v.p, c.len, a.tol * 10, a.dirty_f)) ci[k++] = c;
            }
            ci.resize(k);
        }
        st.germ = (before_j - cl.size()) + (before_i - ci.size());
        say("call", "dropped as present in normal: %llu (%zu junction + %zu insertion left)",
            (unsigned long long)st.germ, cl.size(), ci.size());
    }
    write_vcf(out.c_str(), hdr, cl, ci, qname, sets, sup, L, a.assign, -1);

    std::string base(out);
    if (base.size() > 4 && base.compare(base.size() - 4, 4, ".vcf") == 0) base.resize(base.size() - 4);
    for (int l = 0; l < L; l++)
        write_vcf((base + "." + sets[(size_t)l].name + ".vcf").c_str(), hdr, cl, ci, qname, sets,
                 sup, L, a.assign, l);

    say("call", "records %llu | skipped: unmapped %llu, secondary %llu, MAPQ<%d %llu, <%lldbp %llu, mismatch>%.0f%% %llu",
        (unsigned long long)st.rec, (unsigned long long)st.unmapped, (unsigned long long)st.secondary,
        a.min_mapq, (unsigned long long)st.lowq, (long long)a.min_aln, (unsigned long long)st.shortaln,
        100 * a.max_mm, (unsigned long long)st.highmm);
    say("call", "segment %llu -> split %llu + CIGAR DEL %llu | INS %llu "
               "(%llu of them between the two segments spanning the reference position)",
        (unsigned long long)st.seg, (unsigned long long)st.split,
        (unsigned long long)st.cigar_del, (unsigned long long)st.cigar_ins,
        (unsigned long long)st.ins_split);
    say("call", "aday %zu junction + %zu insertion -> supp>=%u: %zu + %zu -> %s",
        jn.size(), ins.size(), a.min_supp, cl.size(), ci.size(), out.c_str());
    sam_hdr_destroy(hdr);
    return 0;
}

int cmd_call(int argc, const char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: snape call <ctg.bam|name=bam,name=bam,...> <normal.bam|-> <out.vcf> "
                             "[min_mapq=10] [min_aln=0] [min_sv=50] [cluster=50] [min_supp=5] [shred_n=6] "
                             "[shred_win=300] [max_mm=0.01] [shred_rep=0.40] [anchor_rep=0.50] "
                             "[normal_dirty=0.25] [supp_assign=2]\n");
        return 2;
    }
    CallArgs a;
    a.min_mapq   = argi(argc, argv, 4, 10);
    a.min_aln    = (hts_pos_t)argll(argc, argv, 5, 0);
    a.min_sv     = argi(argc, argv, 6, 50);
    a.tol        = (hts_pos_t)argll(argc, argv, 7, 50);
    a.min_supp   = (uint32_t)argi(argc, argv, 8, 5);
    a.shred_n    = argi(argc, argv, 9, 6);
    a.shred_win  = (hts_pos_t)argll(argc, argv, 10, 300);
    a.max_mm     = argd(argc, argv, 11, 0.01);
    a.shred_rep  = argd(argc, argv, 12, 0.40);
    a.max_rep    = argd(argc, argv, 13, 0.50);
    a.dirty_f    = argd(argc, argv, 14, 0.25);
    a.assign     = (uint32_t)argi(argc, argv, 15, 2);
    return run_call(parse_ctg(argv[1]), argv[2], argv[3], a);
}
