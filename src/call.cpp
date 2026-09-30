#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "cmd.h"
#include "extract.h"
#include "files.h"
#include "io.h"
#include "reads.h"

// Every site from select is somatic by construction (its combination has no normal read and at least 3 tumor reads). call
// validates that on the reads themselves: it fetches each site's tumor read, measures the distance between the site's
// innermost anchors on it and keeps the sequence between them, then streams every read of every sample once and, for each
// read that carries a left and a right anchor of a site in the same orientation, measures the same distance. A read whose
// distance differs from the tumor read's by less than DSAME bp carries the same-length allele; by MS bp or more a different
// one. Reads within MS bp are also compared by sequence (bounded edit distance ESAME plus EFRAC of the length). Per site and
// sample this gives spanning, same-length, near, different and same-sequence read counts, symmetric for normal and tumor.
// The records come from the pieces' alignments (gaps of MS bp or more inside an alignment; jumps between consecutive
// alignments give DEL/INS/DUP/INV/BND), attached to a site whose interval is within ROW bp in the read, merged within TOL bp;
// insertions of NEAR bp or more copying the adjacent reference (KM-mer share between ROT and UNQ) are tandem duplications.
// Rules shared by the three sets (thresholds of the published long-read somatic callers): the normal has more than NCOV
// reads at the site and at most CN reads at the tumor allele's length, at most CVAF of its spanning reads (Severus,
// nanomonsv); the tumor has at least MM reads at that length and at least VAF of its depth (Severus, nanomonsv); the
// nearest normal allele (a length cluster of at least 2 reads, CVAF of the normal's reads) is at least MS bp away, so the
// somatic difference is SV-sized; one record per site (the largest event; a junction outranks the gaps next to it); a site
// with no anchor pair to measure needs Sniffles2's coverage-scaled support SB + SM * depth of its graph count. Set 1 also
// drops events in or within VB bp of a tandem repeat (the TRF bed, as Severus' --vntr-bed); set 2 keeps such an event when
// the normal is clean there: at most NC length clusters, NS of its reads in the two largest, the site at most SL bp, VAF2
// and MQ2; set 3 asks VAF2 and MQ2 of every event. Duplicates are judged as the benchmark does: within DUPT bp at both
// ends, insertions by position; the record with more tumor reads stays.
static const uint32_t MS = 50, MQ = 20, MQ2 = 40, MM = 3, NCOV = 3, CN = 1, TOL = 50, ISL = 100, ROW = 200, DUPT = 500, NEAR = 500, KM = 15, DSAME = 20, ESAME = 20, NA = 3, VB = 50, NC = 2, SL = 2000;
static const float VAF = 0.05f, VAF2 = 0.10f, CVAF = 0.03f, EFRAC = 0.02f, ROT = 0.8f, UNQ = 1.25f, SB = 1.5f, SM = 0.1f, NS = 0.9f;
static const int32_t UNSET = INT32_MIN;

struct seg { uint32_t qs, qe, ts, te, mq; char st; std::string t, cg; };

struct ev { std::string ty, c1, c2; uint32_t p1, p2, q, mq, o; char st; };

// a site: interval a..b in its piece, tumor sample s, combination counts n, tumor count t, anchor depths nd/td, the read it
// came from, anchor nodes l/r nearest first, the tumor read's distances d[li][ri], its sequence sq between the innermost
// anchors, and the per-sample measurements m[sample * NM + kind]
enum { SPAN, SAME, NEAR_, DIFF, SEQ, NM };
struct spot { uint32_t a, b, s, t, nd, td; uint64_t read; std::vector<uint32_t> n, l, r, m; std::vector<std::vector<int32_t>> dl; int32_t d[NA][NA]; std::string sq; uint8_t lk; };

struct sv { ev e; std::string nm, why; uint32_t pc, t, mq, st, rank; const spot* r; };

struct bed { std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> m; };

// the read pass: anchors sorted by node (av) with their site/side/index (ai), a flag per node, sites by tumor read (rs, ri)
struct mctx { std::vector<spot>* w; std::vector<uint32_t> av, ai; std::vector<uint8_t> isa; std::vector<uint64_t> rs; std::vector<uint32_t> ri; std::vector<std::vector<uint64_t>> hb; std::vector<int> fs, fx; std::vector<scr> sc; rpass* p; int t; };

static std::vector<uint32_t> nums(const char* s) {
    std::vector<uint32_t> v;
    for (char* p = (char*)s; *p && *p != '\t' && *p != '\n'; p += *p == ',') v.push_back(std::strtoul(p, &p, 10));
    return v;
}

// anchors "node:pos,node:pos,..." nearest first
static void anchors(spot* r, const char* ls, const char* rs) {
    const char* p[2] = {ls, rs};
    for (int k = 0; k < 2; k++)
        for (const char* q = p[k]; *q && *q != '-' && *q != '\t' && *q != '\n';) {
            (k ? r->r : r->l).push_back(std::strtoull(q, (char**)&q, 10));
            std::strtoul(q + 1, (char**)&q, 10);
            if (*q == ',') q++;
        }
}

// depth of a sample at the site: its count at the node the combination enters from and at the node it leaves to (both are
// on every allele's path); the larger of the two
static void depths(spot* r, uint64_t key, const gph* g) {
    uint64_t u = g->e[(key >> 32) - 1] >> 32, w = g->e[key & 0xffffffff] & 0xffffffff;
    r->nd = std::max(g->c[u], g->c[w]);
    r->td = std::max(g->c[r->s * g->K + u], g->c[r->s * g->K + w]);
}

static bool cmp_seg(const seg& x, const seg& y) { return x.qs < y.qs; }

static bool cmp_sv(const sv& x, const sv& y) {
    return x.e.ty != y.e.ty ? x.e.ty < y.e.ty : x.e.c1 != y.e.c1 ? x.e.c1 < y.e.c1 : x.e.c2 != y.e.c2 ? x.e.c2 < y.e.c2 : x.e.o != y.e.o ? x.e.o < y.e.o : x.e.p1 != y.e.p1 ? x.e.p1 < y.e.p1 : x.e.p2 != y.e.p2 ? x.e.p2 < y.e.p2 : x.e.q != y.e.q ? x.e.q < y.e.q : x.nm < y.nm;
}

static bool cmp_pos(const sv& x, const sv& y) {
    return x.e.c1 != y.e.c1 ? x.e.c1 < y.e.c1 : x.e.c2 != y.e.c2 ? x.e.c2 < y.e.c2 : x.e.p1 != y.e.p1 ? x.e.p1 < y.e.p1 : x.e.p2 < y.e.p2;
}

static bool better(const sv& x, const sv& y) { return x.st != y.st ? x.st > y.st : x.mq != y.mq ? x.mq > y.mq : x.pc > y.pc; }

// the benchmark's match: both ends within w bp (either orientation), insertions by position
static bool match(const ev& x, const ev& y, uint32_t w) {
    if (x.c1 != y.c1 || x.p1 + w < y.p1 || y.p1 + w < x.p1) return false;
    return x.ty == "INS" || y.ty == "INS" || (x.c2 == y.c2 && x.p2 + w >= y.p2 && y.p2 + w >= x.p2);
}

// the length clusters of a sample's spanning reads (single linkage within DSAME bp), sizes largest first
static std::vector<uint32_t> clusters(std::vector<int32_t> d) {
    std::vector<uint32_t> o;
    std::sort(d.begin(), d.end());
    for (size_t i = 0, j; i < d.size(); i = j) {
        for (j = i; j < d.size() && (j == i || d[j] - d[j - 1] <= (int32_t)DSAME); j++);
        o.push_back(j - i);
    }
    std::sort(o.begin(), o.end());
    std::reverse(o.begin(), o.end());
    return o;
}

// the nearest allele of a sample to the tumor allele: a length cluster of at least two reads and CVAF of the sample's
// spanning reads; NONE when the sample has no such cluster
static uint32_t nearest(std::vector<int32_t> d) {
    std::sort(d.begin(), d.end());
    uint32_t best = NONE;
    for (size_t i = 0, j; i < d.size(); i = j) {
        double sum = 0;
        for (j = i; j < d.size() && (j == i || d[j] - d[j - 1] <= (int32_t)DSAME); j++) sum += d[j];
        if (j - i >= std::max<size_t>(2, (size_t)std::ceil(CVAF * d.size()))) best = std::min<uint32_t>(best, std::fabs(sum / (j - i)));
    }
    return best;
}

static void load_bed(const char* path, bed* b) {
    char c[256];
    uint32_t s, e;
    FILE* f = std::fopen(path, "r");
    while (std::fscanf(f, "%255s %u %u%*[^\n]", c, &s, &e) == 3) b->m[c].push_back({s, e});
    std::fclose(f);
    for (auto& kv : b->m) std::sort(kv.second.begin(), kv.second.end());
}

// whether position p of chromosome c lies in or within VB bp of a tandem repeat
static bool inbed(const bed* b, const std::string& c, uint32_t p) {
    auto it = b->m.find(c);
    if (it == b->m.end()) return false;
    const std::vector<std::pair<uint32_t, uint32_t>>& v = it->second;
    size_t i = std::upper_bound(v.begin(), v.end(), std::make_pair(p + VB, ~0u)) - v.begin();
    for (size_t j = i > 3 ? i - 3 : 0; j < i; j++) if (v[j].second + VB >= p) return true;
    return false;
}

// bounded edit distance between two sequences (common ends trimmed first); returns X when it is X or more
static uint32_t edits(const std::string& a, const std::string& b, uint32_t X) {
    uint32_t k = std::min(a.size(), b.size()), s = 0, e = 0;
    while (s < k && a[s] == b[s]) s++;
    while (e < k - s && a[a.size() - 1 - e] == b[b.size() - 1 - e]) e++;
    uint32_t n = a.size() - s - e, m = b.size() - s - e;
    if (n > m + X || m > n + X) return X;
    const char *x = a.data() + s, *y = b.data() + s;
    std::vector<uint32_t> p(m + 2, X), c(m + 2, X);
    for (uint32_t j = 0; j <= m && j < X; j++) p[j] = j;
    for (uint32_t i = 1; i <= n; i++) {
        uint32_t lo = i > X ? i - X : 1, hi = std::min(m, i + X), best = X;
        c[lo - 1] = lo == 1 ? std::min(i, X) : X;
        for (uint32_t j = lo; j <= hi; j++) best = std::min(best, c[j] = std::min(std::min(p[j] + 1, c[j - 1] + 1), std::min(p[j - 1] + (x[i - 1] != y[j - 1]), X)));
        if (hi < m) c[hi + 1] = X;
        if (best >= X) return X;
        p.swap(c);
    }
    return std::min(p[m], X);
}

// the anchor hits of one read strand, one 64-bit key each: site, side, anchor index, then the position in the strand's list
static void hits(mctx* x, const rth* c, uint32_t m, std::vector<uint64_t>& h) {
    h.clear();
    for (uint32_t i = 0; i < m; i++) {
        uint32_t v = c->id[i];
        if (!x->isa[v]) continue;
        for (size_t j = std::lower_bound(x->av.begin(), x->av.end(), v) - x->av.begin(); j < x->av.size() && x->av[j] == v; j++) h.push_back((uint64_t)x->ai[j] << 32 | i);
    }
    std::sort(h.begin(), h.end());
}

// the tumor read of its own sites: every left/right pair present gives d[li][ri] (the anchors' order in the strand that
// carries them is the same in every read, since a node is one orientation of a core); the innermost pair gives the sequence
static void own(rpass* p, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) {
    mctx* x = (mctx*)p->ctx;
    std::vector<uint64_t>& h = x->hb[c - p->t.data()];
    hits(x, c, m, h);
    for (size_t i = 0, j; i < h.size(); i = j) {
        uint32_t s = h[i] >> 35;
        for (j = i; j < h.size() && h[j] >> 35 == s; j++);
        spot* w = &(*x->w)[s];
        if (w->read != ((uint64_t)p->f << 40 | rid)) continue;
        for (size_t a = i; a < j; a++)
            for (size_t b = i; b < j; b++) {
                uint32_t la = h[a] >> 32 & 7, rb = h[b] >> 32 & 7, pa = (uint32_t)h[a], pb = (uint32_t)h[b];
                if (la >= 4 || rb < 4 || pa == pb) continue;
                uint32_t li = la, ri = rb - 4, al = c->r.a[c->ix[pa]], ar = c->r.a[c->ix[pb]];
                if (w->d[li][ri] != UNSET) continue;
                w->d[li][ri] = al < ar ? ar - al : al - ar;
                if (!li && !ri) w->sq = bases(&c->r, rec, std::min(al, ar), std::max(al, ar), rc);
            }
    }
}

// every read of every sample: the innermost pair present in strand order gives the distance; count it against the tumor's
static void measure(rpass* p, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t) {
    mctx* x = (mctx*)p->ctx;
    std::vector<uint64_t>& h = x->hb[c - p->t.data()];
    hits(x, c, m, h);
    for (size_t i = 0, j; i < h.size(); i = j) {
        uint32_t s = h[i] >> 35, li = NONE, ri = NONE, pl = 0, pr = 0;
        for (j = i; j < h.size() && h[j] >> 35 == s; j++) {
            uint32_t k = h[j] >> 32 & 7, pos = (uint32_t)h[j];
            if (k < 4 && li == NONE) { li = k; pl = pos; }
            if (k >= 4 && ri == NONE) { ri = k - 4; pr = pos; }
        }
        spot* w = &(*x->w)[s];
        if (li == NONE || ri == NONE || pl == pr || w->d[li][ri] == UNSET) continue;
        uint32_t al = c->r.a[c->ix[pl]], ar = c->r.a[c->ix[pr]], *cnt = w->m.data() + p->f * NM;
        if (al > ar) std::swap(al, ar);
        int32_t dl = (int32_t)(ar - al) - w->d[li][ri], ad = dl < 0 ? -dl : dl;
        __atomic_fetch_add(&cnt[SPAN], 1, __ATOMIC_RELAXED);
        __atomic_fetch_add(&cnt[ad < (int32_t)DSAME ? SAME : ad < (int32_t)MS ? NEAR_ : DIFF], 1, __ATOMIC_RELAXED);
        while (__atomic_exchange_n(&w->lk, 1, __ATOMIC_ACQUIRE));
        w->dl[p->f].push_back(dl);
        __atomic_store_n(&w->lk, 0, __ATOMIC_RELEASE);
        if (ad < (int32_t)MS && !li && !ri && !w->sq.empty()) {
            std::string q = bases(&c->r, rec, al, ar, rc);
            uint32_t X = std::max(ESAME, (uint32_t)(EFRAC * w->sq.size()));
            if (edits(q, w->sq, X) < X) __atomic_fetch_add(&cnt[SEQ], 1, __ATOMIC_RELAXED);
        }
    }
}

static void fetch(const mctx* x, uint64_t rd, std::vector<uint8_t>& rec) {
    uint64_t off;
    uint32_t tag, s = rd >> 40;
    pread(x->fx[s], &off, 8, sizeof(bin_header) + 8 * (rd & 0xffffffffff));
    pread(x->fs[s], &tag, 4, off);
    rec.assign(4 + (tag >> 31 ? tag & 0x7fffffff : ((tag & 0x7fffffff) + 3) / 4) + 8, 0);
    pread(x->fs[s], rec.data(), rec.size() - 8, off);
}

// thread i fetches every t-th tumor read that has sites, recomputes its cores and runs own() on both strands
static void grab(void* a, int i) {
    mctx* x = (mctx*)a;
    rth* c = &x->p->t[i];
    std::vector<uint8_t> rec;
    for (size_t k = i; k + 1 < x->ri.size(); k += x->t) {
        fetch(x, x->rs[k], rec);
        cores(&x->sc[i], &c->r, rec.data(), LEV);
        lookup(x->p, c);
        x->p->f = x->rs[k] >> 40;
        chains(x->p, c, rec.data(), x->rs[k] & 0xffffffffff);
    }
}

static void add(std::vector<ev>& o, uint32_t q, const char* ty, std::string c1, uint32_t p1, std::string c2, uint32_t p2, const seg& x, const seg& y, uint32_t d = 0) {
    if (c1 == c2 && (p1 > p2 ? p1 - p2 : p2 - p1) < MS && std::strcmp(ty, "INS")) return;
    if (c1 > c2 || (c1 == c2 && p1 > p2)) { std::swap(c1, c2); std::swap(p1, p2); d = (d & 1) * 2 + d / 2; }
    o.push_back({ty, c1, c2, p1, p2, q, std::min(x.mq, y.mq), d, x.st});
}

// events of one piece: gaps of MS bp or more inside an alignment, and the jumps between consecutive alignments along the read
static std::vector<ev> events(std::vector<seg>& sg) {
    std::vector<ev> o;
    std::sort(sg.begin(), sg.end(), cmp_seg);
    for (const seg& s : sg) {
        if (s.mq < MQ) continue;
        uint32_t t = s.ts, q = s.st == '+' ? s.qs : s.qe, gt = 0, gq = 0, gl = 0, gap = ISL;
        char gop = 0;
        for (const char* p = s.cg.c_str(); *p;) {
            uint32_t l = std::strtoul(p, (char**)&p, 10);
            char op = *p++;
            if ((op == 'I' || op == 'D') && l >= MS) {
                if (op == gop && gap < ISL) gl += l + (op == 'D' ? gap : 0);
                else {
                    if (gop) add(o, gq, gop == 'I' ? "INS" : "DEL", s.t, gt, s.t, gt + gl, s, s);
                    gop = op; gt = t; gq = op == 'I' && s.st != '+' ? q - l : q; gl = l;
                }
                gap = 0;
            } else if (op != 'I' && op != 'D') gap += l;
            if (op != 'I') t += l;
            if (op != 'D') q = s.st == '+' ? q + l : q - l;
        }
        if (gop) add(o, gq, gop == 'I' ? "INS" : "DEL", s.t, gt, s.t, gt + gl, s, s);
    }
    for (size_t i = 0; i + 1 < sg.size(); i++) {
        const seg &x = sg[i], &y = sg[i + 1];
        if (x.mq < MQ || y.mq < MQ) continue;
        uint32_t p1 = x.st == '+' ? x.te : x.ts, p2 = y.st == '+' ? y.ts : y.te, q = x.qe, b = (x.st == '+') * 2 + (y.st == '-');
        int64_t g = (int64_t)y.qs - x.qe, d = x.st == '+' ? (int64_t)y.ts - x.te : (int64_t)x.ts - y.te;
        if (x.t != y.t) add(o, q, "BND", x.t, p1, y.t, p2, x, y, b);
        else if (x.st != y.st) add(o, q, "INV", x.t, p1, y.t, p2, x, y, b);
        else if (d - g >= MS) add(o, q, "DEL", x.t, p1, y.t, p2, x, y);
        else if (g - d >= MS && d > -(int64_t)MS) add(o, q, "INS", x.t, p1, x.t, p1 + g - d, x, y);
        else if (d <= -(int64_t)MS) add(o, q, "DUP", x.t, p1, y.t, p2, x, y);
    }
    return o;
}

static bool near(const ev& x, const ev& y, uint32_t w) {
    return x.c1 == y.c1 && x.c2 == y.c2 && x.p1 + w >= y.p1 && y.p1 + w >= x.p1 && x.p2 + w >= y.p2 && y.p2 + w >= x.p2;
}

static bool primary(const std::string& c) {
    return c.size() <= 5 && c.find('_') == std::string::npos && c != "chrM" && c != "M";
}

static std::string ref(FILE* f, const std::vector<uint64_t>& v, uint32_t a, uint32_t b) {
    uint64_t s = v[1] + a / v[2] * v[3] + a % v[2], e = v[1] + b / v[2] * v[3] + b % v[2];
    std::string r(e - s, 0), o;
    std::fseek(f, s, SEEK_SET);
    r.resize(std::fread(&r[0], 1, e - s, f));
    for (char ch : r) if (ch != '\n') o += ch & 0xdf;
    return o;
}

static float shared(const std::string& a, const std::string& b) {
    std::unordered_set<std::string> k;
    uint32_t m = 0;
    for (size_t i = 0; i + KM <= b.size(); i++) k.insert(b.substr(i, KM));
    for (size_t i = 0; i + KM <= a.size(); i++) m += k.count(a.substr(i, KM));
    return (float)m / (a.size() - KM + 1);
}

static void tandem(ev& e, const std::string& pc, FILE* fa, const std::vector<uint64_t>& v) {
    uint32_t L = e.p2 - e.p1;
    std::string s = pc.substr(e.q, L), o;
    if (e.st != '+') { for (size_t i = s.size(); i--;) o += s[i] == 'A' ? 'T' : s[i] == 'C' ? 'G' : s[i] == 'G' ? 'C' : s[i] == 'T' ? 'A' : 'N'; s = o; }
    float u = shared(s, ref(fa, v, e.p1 - L, e.p1)), d = shared(s, ref(fa, v, e.p1, std::min<uint64_t>(e.p1 + L, v[0])));
    if (u + d < ROT || u + d > UNQ) return;
    e.ty = "DUP";
    e.p1 -= (uint32_t)(u / (u + d) * L);
    e.p2 = e.p1 + L;
}

static std::string list(const std::vector<uint32_t>& m, uint32_t S, uint32_t kind) {
    std::string s;
    for (uint32_t i = 0; i < S; i++) s += (i ? "," : "") + std::to_string(m[i * NM + kind]);
    return s;
}

// the sorted length differences of a sample's spanning reads, "-" when there are none
static std::string deltas(std::vector<int32_t> d) {
    std::string s;
    std::sort(d.begin(), d.end());
    for (int32_t v : d) s += (s.empty() ? "" : ",") + std::to_string(v);
    return s.empty() ? "-" : s;
}

int cmd_call(const char* in, const char* rd, const char* fasta, const char* gp, const char* out, const char* vb, int t) {
    LCP_INIT();
    std::string pre(in);
    std::vector<spot> w;
    std::unordered_map<std::string, std::vector<uint32_t>> bypc;
    std::unordered_map<std::string, std::vector<seg>> al;
    char line[1 << 16], nm[256], cn[4096], ls[4096], rs[4096], st, tn[256];
    uint32_t qs, qe, ts, te, mq, z, a, b, s;
    unsigned long key;
    bed vn;
    if (vb) load_bed(vb, &vn);
    gph g;
    read_gph(gp, &g);
    uint64_t K = g.K, S = g.S;
    FILE* fi = std::fopen((pre + ".tsv").c_str(), "r");
    while (std::fgets(line, sizeof line, fi)) {
        std::sscanf(line, "%255s %u %u %*u %4095s %lu %u %4095s %4095s", nm, &a, &b, cn, &key, &s, ls, rs);
        spot r = {a, b, s, 0, 0, 0, 0, nums(cn), {}, {}, std::vector<uint32_t>(S * NM, 0), std::vector<std::vector<int32_t>>(S), {}, "", 0};
        r.t = r.n[s];
        for (uint32_t i = 0; i < NA; i++) for (uint32_t j = 0; j < NA; j++) r.d[i][j] = UNSET;
        unsigned long rs_, ss_;
        std::sscanf(nm, "s%lu_%lu_", &ss_, &rs_);
        r.read = ss_ << 40 | rs_;
        anchors(&r, ls, rs);
        depths(&r, key, &g);
        bypc[nm].push_back(w.size());
        w.push_back(r);
    }
    std::fclose(fi);
    g.c = {};
    g.e = {};
    // the read pass: which nodes are anchors of which sites, the tumor reads of the sites, then every read of every sample
    mctx x = {&w, {}, {}, std::vector<uint8_t>(K, 0), {}, {}, std::vector<std::vector<uint64_t>>(t), {}, {}, std::vector<scr>(t), nullptr, t};
    std::vector<std::pair<uint32_t, uint32_t>> an;
    for (uint32_t i = 0; i < w.size(); i++) {
        for (uint32_t k = 0; k < w[i].l.size() && k < NA; k++) an.push_back({w[i].l[k], i << 3 | k});
        for (uint32_t k = 0; k < w[i].r.size() && k < NA; k++) an.push_back({w[i].r[k], i << 3 | 4 | k});
        x.rs.push_back(w[i].read);
    }
    std::sort(an.begin(), an.end());
    for (auto& p : an) { x.av.push_back(p.first); x.ai.push_back(p.second); x.isa[p.first] = 1; }
    std::sort(x.rs.begin(), x.rs.end());
    x.rs.erase(std::unique(x.rs.begin(), x.rs.end()), x.rs.end());
    x.ri.resize(x.rs.size() + 1);
    for (std::string& f : split(rd)) { x.fs.push_back(fileno(std::fopen(f.c_str(), "rb"))); x.fx.push_back(fileno(std::fopen((f + ".idx").c_str(), "rb"))); }
    rpass p = {std::vector<rth>(t), &g.g, {}, 0, own, &x, 0, 0};
    make_map(&p);
    x.p = &p;
    par(t, grab, &x);
    x.sc = {};
    std::fprintf(stderr, "[call] %zu sites on %zu tumor reads measured\n", w.size(), x.rs.size());
    each_read(&p, g.g, rd, t, measure, &x);
    p = {};
    g = {};
    // records from the pieces' alignments
    fi = std::fopen((pre + ".paf").c_str(), "r");
    while (std::fgets(line, sizeof line, fi)) {
        std::sscanf(line, "%255s %u %u %u %c %255s %u %u %u %u %u %u", nm, &z, &qs, &qe, &st, tn, &z, &ts, &te, &z, &z, &mq);
        const char* cg = std::strstr(line, "cg:Z:");
        al[nm].push_back({qs, qe, ts, te, mq, st, tn, cg ? std::string(cg + 5, std::strcspn(cg + 5, "\t\n")) : ""});
    }
    std::fclose(fi);
    std::unordered_map<std::string, std::vector<uint64_t>> fai;
    FILE* fx = std::fopen((std::string(fasta) + ".fai").c_str(), "r");
    unsigned long f4[4];
    while (std::fscanf(fx, "%255s %lu %lu %lu %lu", nm, &f4[0], &f4[1], &f4[2], &f4[3]) == 5) fai[nm] = {f4[0], f4[1], f4[2], f4[3]};
    std::fclose(fx);
    std::vector<sv> som;
    std::unordered_map<std::string, std::string> sq;
    for (auto& pc : bypc) {
        if (!al.count(pc.first)) continue;
        for (const ev& e : events(al[pc.first])) {
            const spot* r = nullptr;
            for (uint32_t i : pc.second) if (e.q + ROW >= w[i].a && e.q <= w[i].b + ROW && (!r || w[i].t > r->t)) r = &w[i];
            if (!r || !primary(e.c1) || !primary(e.c2) || r->t < MM) continue;
            som.push_back({e, pc.first, "", 1, r->t, e.mq, 0, 0, r});
            if (e.ty == "INS" && e.p2 - e.p1 >= NEAR) sq[pc.first];
        }
    }
    FILE *fa = std::fopen(fasta, "r"), *fp = std::fopen((pre + ".fa").c_str(), "r");
    char* ln = nullptr;
    size_t cap = 0;
    std::string* cur = nullptr;
    for (ssize_t n; (n = getline(&ln, &cap, fp)) > 0;)
        if (ln[0] == '>') { auto it = sq.find(std::string(ln + 1, n - 2)); cur = it == sq.end() ? nullptr : &it->second; }
        else if (cur) cur->assign(ln, n - 1);
    std::free(ln);
    std::fclose(fp);
    for (sv& v : som) if (v.e.ty == "INS" && v.e.p2 - v.e.p1 >= NEAR && v.e.p1 > v.e.p2 - v.e.p1) tandem(v.e, sq[v.nm], fa, fai[v.e.c1]);
    std::fclose(fa);
    // one record per event: same type and orientation, both ends within TOL; the best piece (most tumor reads) represents it
    std::sort(som.begin(), som.end(), cmp_sv);
    std::vector<sv> res;
    for (const sv& v : som) {
        if (res.empty() || res.back().e.ty != v.e.ty || res.back().e.o != v.e.o || !near(res.back().e, v.e, TOL)) { res.push_back(v); continue; }
        sv& q = res.back();
        uint32_t pc = q.pc + 1, m = std::max(q.mq, v.mq);
        if (better(v, q)) q = v;
        q.pc = pc;
        q.mq = m;
    }
    // the measurements behind the rules: tumor reads at the tumor allele's length (the graph count when the site has no
    // anchor pair to measure), the normal's; per site the number of records and which is the largest event
    std::unordered_map<const spot*, std::vector<uint32_t>> bysite;
    for (uint32_t i = 0; i < res.size(); i++) {
        sv& v = res[i];
        const uint32_t* mn = v.r->m.data(), *mt = v.r->m.data() + v.r->s * NM;
        v.st = mn[SPAN] || mt[SPAN] ? mt[SAME] : v.r->t;
        v.rank = v.e.ty == "BND" || v.e.ty == "INV" ? 3 : v.e.p2 - v.e.p1 >= NEAR ? 2 : 1;
        bysite[v.r].push_back(i);
    }
    for (int mode = 1; mode <= 3; mode++) {
        for (sv& v : res) {
            const uint32_t* mn = v.r->m.data(), *mt = v.r->m.data() + v.r->s * NM;
            const std::vector<uint32_t>& g_ = bysite[v.r];
            const std::vector<int32_t>& dn = v.r->dl[0];
            std::vector<uint32_t> cn_ = clusters(dn);
            uint32_t nn = nearest(dn), L = v.e.p2 - v.e.p1, nclu = 0;
            for (uint32_t c_ : cn_) nclu += c_ >= 2;
            bool top = true, measured = mn[SPAN] || mt[SPAN], tr = inbed(&vn, v.e.c1, v.e.p1);
            for (uint32_t i : g_) top &= res[i].rank < v.rank || (res[i].rank == v.rank && res[i].e.p2 - res[i].e.p1 <= L);
            float vaf = (float)v.st / std::max(v.r->td, v.st), ns = dn.empty() ? 0 : (float)(cn_[0] + (cn_.size() > 1 ? cn_[1] : 0)) / dn.size();
            bool clean = mt[SPAN] && nclu <= NC && ns >= NS && v.r->b - v.r->a <= SL && vaf >= VAF2 && v.mq >= MQ2;
            v.why = v.r->nd <= NCOV ? "LOWNORMAL" : mn[SAME] > CN || mn[SAME] > CVAF * mn[SPAN] ? "GERMLINE" : nn != NONE && nn < MS ? "SMALL" : v.st < MM || vaf < VAF ? "LOWVAF"
                  : g_.size() > 1 && !top ? "COMPLEX" : !measured && v.st < std::lround(SB + SM * v.r->td) ? "LOWSUPPORT"
                  : tr && (mode == 1 || (mode == 2 && !clean) || mode == 3) ? "VNTR" : mode == 3 && (vaf < VAF2 || v.mq < MQ2) ? "LOWVAF" : "PASS";
        }
        // records that describe one event within the benchmark tolerance: the better one stays
        std::sort(res.begin(), res.end(), cmp_pos);
        for (size_t i = 0; i < res.size(); i++)
            for (size_t j = i + 1; j < res.size() && res[j].e.c1 == res[i].e.c1 && res[j].e.p1 <= res[i].e.p1 + DUPT; j++)
                if (res[i].why == "PASS" && res[j].why == "PASS" && match(res[i].e, res[j].e, DUPT)) res[better(res[i], res[j]) ? j : i].why = "DUPL";
        FILE* o = std::fopen((std::string(out) + "." + std::to_string(mode) + ".vcf").c_str(), "w");
        std::fprintf(o, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (uint64_t i = 0; i < S; i++) std::fprintf(o, "\tS%lu", (unsigned long)i);
        std::fprintf(o, "\n");
        uint64_t np = 0;
        for (size_t i = 0; i < res.size(); i++) {
            const sv& v = res[i];
            const spot* r = v.r;
            std::string dls;
            for (uint64_t k = 0; k < S; k++) dls += (k ? ";D" : "D") + std::to_string(k) + "=" + deltas(r->dl[k]);
            np += v.why == "PASS";
            std::string j = (v.e.o & 1 ? "]" : "[") + v.e.c2 + ":" + std::to_string(v.e.p2) + (v.e.o & 1 ? "]" : "["), alt = v.e.ty != "BND" ? "<" + v.e.ty + ">" : v.e.o & 2 ? "N" + j : j + "N";
            std::fprintf(o, "%s\t%u\tsv%zu\tN\t%s\t.\t%s\tSVTYPE=%s;END=%u;SVLEN=%u;O=%u;PIECES=%u;PIECE=%s;QPOS=%u;SITE=%u-%u;MQ=%u;TDP=%u;NDP=%u;SPAN=%s;SAME=%s;NEAR=%s;DIFF=%s;SEQ=%s;ST=%u;%s\tGT:DV", v.e.c1.c_str(), v.e.p1, i, alt.c_str(), v.why.c_str(), v.e.ty.c_str(),
                         v.e.ty == "BND" || v.e.ty == "INS" ? v.e.p1 : v.e.p2, v.e.c1 == v.e.c2 ? v.e.p2 - v.e.p1 : 0, v.e.o, v.pc, v.nm.c_str(), v.e.q, r->a, r->b, v.mq, r->td, r->nd,
                         list(r->m, S, SPAN).c_str(), list(r->m, S, SAME).c_str(), list(r->m, S, NEAR_).c_str(), list(r->m, S, DIFF).c_str(), list(r->m, S, SEQ).c_str(), v.st, dls.c_str());
            for (uint32_t n : r->n) std::fprintf(o, "\t%s:%u", n >= MM ? "0/1" : "0/0", n);
            std::fprintf(o, "\n");
        }
        std::fclose(o);
        std::fprintf(stderr, "[call] rule set %d: %lu PASS\n", mode, (unsigned long)np);
    }
    std::fprintf(stderr, "[call] %zu events at sites, %zu records\n", som.size(), res.size());
    return 0;
}
