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

// call reads the sites (ev.tsv), the pieces (ev.fa) and their alignments (ev.paf), the graph and the packed reads, makes
// its own pass over every read of both samples, and writes the records with their evidence features for the verdict.
// The pass: each site's tumor read is fetched, the distance between the site's innermost anchors is measured on it and
// the read's nodes between them are counted (tumor-only, branch); then every read is streamed once, and each read that
// carries a left and a right anchor of a site gives its distance, so every site gets per-sample spanning reads, their
// length differences from the tumor allele, reads within DSAME of it, and the same over consistent reads (their other
// anchors at the tumor read's spacing). Records come from the alignments (gaps of MS bp or more inside an alignment;
// jumps between consecutive alignments give DEL/INS/DUP/INV/BND), attached to a site within ROW bp in the piece, merged
// within TOL bp; insertions of NEAR bp or more copying the adjacent reference (KM-mer share between ROT and UNQ) become
// tandem duplications. Every record gets NF features from the measurements, the reference repeat content in WIN bp,
// the tandem-repeat bed, and the alignment; its error probability is a logistic with weights W0 and PASS is the PCUT
// cut; duplicates within DUPT bp keep the likelier record. The verdict (bench/verdict.py) reads the F= features.
static const uint32_t MS = 50, MQ = 10, TOL = 50, NEAR = 500, ISL = 100, ROW = 200, WIN = 150, NF = 36, DUPT = 500, KM = 15, DSAME = 20, SAME = 20, NA = 3, CL = 3;
static const float PCUT = 0.5f, ROT = 0.8f, UNQ = 1.25f, W0[NF] = {-6, 1.5f, 1, 0.5f, -0.3f, 1, 0.5f, -4, -3, -1.5f, -1.5f, -1, -0.05f, -0.05f, 3, -0.5f, -1.5f, -1, 0, 0, 0, 0, 0, -0.2f, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static const int32_t UNSET = INT32_MIN;

struct seg { uint32_t qs, qe, ts, te, mq; char st; std::string t, cg; float de; };

struct ev { std::string ty, c1, c2; uint32_t p1, p2, q, mq, ln, l1, l2; float de; uint32_t o; char st; };

// the evidence row of a site: interval, site index, tumor sample, tumor length clusters and the piece's cluster size,
// spanning tumor/normal reads, normal reads at the tumor length, shadow (none here), tumor-only and branch nodes and
// bases between the anchors, the normal's median/min/max length difference, the same counts over consistent reads,
// the combination counts n and the depths d per sample
struct arow { uint32_t a, b, u, s, ncl, cl, tspan, nspan, nsame, shadow, tonly, nbr, span; int32_t lest, dmin, dmax; uint32_t nspc, n20c, tspc; std::vector<uint32_t> n, d; };

struct sv { ev e; std::string nm; uint32_t pc; arow r; uint32_t nf, ms, mn; };

struct bed { std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> m; };

// a site as measured: interval a..b in its piece, tumor sample s, tumor count t, depths, tumor-only and branch nodes
// between the innermost anchors, the read and combination it came from, counts n, anchor nodes l/r nearest first,
// per-sample measurements m[sample * NM + kind] and length differences dl, the tumor read's distances d[li][ri]
enum { SPAN, SAME_, OKS, OK20, NM };
struct spot { uint32_t a, b, s, t, nd, td, tonly, nbr; uint64_t read, key; std::vector<uint32_t> n, l, r, m, dep; std::vector<std::vector<int32_t>> dl; int32_t d[NA][NA]; std::string nm; uint8_t lk; };

// the read pass: anchors sorted by node (av) with their site/side/index (ai), flags per node (anchor, no normal read,
// branch), the sites' tumor reads (rs), per-thread hit buffers, the read files and indexes, core scratch, the map
struct mctx { std::vector<spot>* w; std::vector<uint32_t> av, ai; std::vector<uint8_t> isa, nz, br; std::vector<uint64_t> rs; std::vector<std::vector<uint64_t>> hb; std::vector<int> fs, fx; std::vector<scr> sc; rpass* p; int t; };

static std::vector<uint32_t> nums(const char* s) {
    std::vector<uint32_t> v;
    for (char* p = (char*)s; *p && *p != '\t' && *p != '\n'; p += *p == ',') v.push_back(std::strtoul(p, &p, 10));
    return v;
}

static uint32_t median(std::vector<uint32_t> v) {
    if (v.empty()) return 0;
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}

// anchors "node:pos,node:pos,..." nearest first; the depth of a sample at the site is the median count of its anchors (unique
// cores, so their counts are the local read depth) when it has any, else its count at the node the combination enters from
// or leaves to (on every allele's path, but a high-copy repeat there counts every copy). The normal's depth test takes the
// larger of the two.
static void anchors(spot* r, const char* ls, const char* rs, uint64_t key, const gph* g) {
    const char* p[2] = {ls, rs};
    std::vector<uint32_t> cn, ct;
    for (int k = 0; k < 2; k++)
        for (const char* q = p[k]; *q && *q != '-' && *q != '\t' && *q != '\n';) {
            uint64_t v = std::strtoull(q, (char**)&q, 10);
            (k ? r->r : r->l).push_back(v);
            cn.push_back(g->c[v]);
            ct.push_back(g->c[r->s * g->K + v]);
            std::strtoul(q + 1, (char**)&q, 10);
            if (*q == ',') q++;
        }
    uint64_t u = g->e[(key >> 32) - 1] >> 32, w = g->e[key & 0xffffffff] & 0xffffffff;
    r->nd = std::max<uint32_t>(std::max(g->c[u], g->c[w]), median(cn));
    r->td = ct.empty() ? std::max(g->c[r->s * g->K + u], g->c[r->s * g->K + w]) : median(ct);
    for (uint64_t k = 0; k < g->S; k++) r->dep.push_back(r->l.empty() || r->r.empty() ? std::max(g->c[k * g->K + u], g->c[k * g->K + w]) : std::min(g->c[k * g->K + r->l[0]], g->c[k * g->K + r->r[0]]));
}

static uint32_t fit(const arow& r, const ev& e) {
    uint32_t o = e.q < r.a || e.q > r.b, L = e.p2 - e.p1;
    return e.ty == "BND" || e.ty == "INV" ? o : o ? 2 * L : r.nspan ? std::abs(std::abs(r.lest) - (int32_t)L) : L;
}

static bool cmp_seg(const seg& x, const seg& y) { return x.qs < y.qs; }

static bool cmp_sv(const sv& x, const sv& y) {
    return x.e.ty != y.e.ty ? x.e.ty < y.e.ty : x.e.c1 != y.e.c1 ? x.e.c1 < y.e.c1 : x.e.c2 != y.e.c2 ? x.e.c2 < y.e.c2 : x.e.o != y.e.o ? x.e.o < y.e.o : x.e.p1 != y.e.p1 ? x.e.p1 < y.e.p1 : x.e.p2 != y.e.p2 ? x.e.p2 < y.e.p2 : x.e.q != y.e.q ? x.e.q < y.e.q : x.nm < y.nm;
}

static uint32_t top(const std::vector<uint32_t>& n) {
    return *std::max_element(n.begin() + 1, n.end());
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

static void load_bed(const char* path, bed* b) {
    char c[256];
    uint32_t s, e;
    FILE* f = std::fopen(path, "r");
    while (std::fscanf(f, "%255s %u %u%*[^\n]", c, &s, &e) == 3) b->m[c].push_back({s, e});
    std::fclose(f);
    for (auto& kv : b->m) {
        std::vector<std::pair<uint32_t, uint32_t>>& v = kv.second;
        std::sort(v.begin(), v.end());
        size_t j = 0;
        for (size_t i = 1; i < v.size(); i++)
            if (v[i].first <= v[j].second) v[j].second = std::max(v[j].second, v[i].second);
            else v[++j] = v[i];
        v.resize(j + 1);
    }
}

static bool inbed(const bed* b, const std::string& c, uint32_t p) {
    auto it = b->m.find(c);
    if (it == b->m.end()) return false;
    auto k = std::upper_bound(it->second.begin(), it->second.end(), std::make_pair(p + 50, ~0u));
    return k != it->second.begin() && (k - 1)->second + 50 >= p;
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
static void own(rpass* p, rth* c, const uint8_t*, uint32_t m, uint32_t, uint64_t rid) {
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
                if (!li && !ri) {
                    w->tonly = w->nbr = 0;
                    for (uint32_t k = std::min(pa, pb); k <= std::max(pa, pb); k++) { w->tonly += x->nz[c->id[k]]; w->nbr += x->br[c->id[k]]; }
                }
            }
    }
}

// every read of every sample: the innermost pair present in strand order gives the distance; count it against the tumor's
static void measure(rpass* p, rth* c, const uint8_t*, uint32_t m, uint32_t, uint64_t) {
    mctx* x = (mctx*)p->ctx;
    std::vector<uint64_t>& h = x->hb[c - p->t.data()];
    hits(x, c, m, h);
    for (size_t i = 0, j; i < h.size(); i = j) {
        uint32_t s = h[i] >> 35, li = NONE, ri = NONE, pl = 0, pr = 0, ps[2][NA];
        for (uint32_t k = 0; k < NA; k++) ps[0][k] = ps[1][k] = NONE;
        for (j = i; j < h.size() && h[j] >> 35 == s; j++) {
            uint32_t k = h[j] >> 32 & 7, pos = (uint32_t)h[j];
            if (k < 4 && li == NONE) { li = k; pl = pos; }
            if (k >= 4 && ri == NONE) { ri = k - 4; pr = pos; }
            if (ps[k >= 4][k & 3] == NONE) ps[k >= 4][k & 3] = pos;
        }
        spot* w = &(*x->w)[s];
        if (li == NONE || ri == NONE || pl == pr || w->d[li][ri] == UNSET) continue;
        uint32_t al = c->r.a[c->ix[pl]], ar = c->r.a[c->ix[pr]], *cnt = w->m.data() + p->f * NM;
        if (al > ar) std::swap(al, ar);
        int32_t dl = (int32_t)(ar - al) - w->d[li][ri], ad = dl < 0 ? -dl : dl;
        // consistent read: its other anchors on each side sit at the tumor read's spacing from the innermost one
        uint32_t ok[2] = {0, 0};
        for (uint32_t sd = 0; sd < 2; sd++)
            for (uint32_t k = 1; k < NA; k++) {
                if (ps[sd][0] == NONE || ps[sd][k] == NONE) continue;
                int32_t t0 = sd ? w->d[0][0] : w->d[0][0], tk = sd ? w->d[0][k] : w->d[k][0];
                if (t0 == UNSET || tk == UNSET) continue;
                int32_t sp = std::labs(t0 - tk), got = std::labs((int32_t)c->r.a[c->ix[ps[sd][k]]] - (int32_t)c->r.a[c->ix[ps[sd][0]]]);
                ok[sd] += std::labs(got - sp) <= (int32_t)DSAME;
            }
        if (ok[0] >= NA - 1 && ok[1] >= NA - 1) { __atomic_fetch_add(&cnt[OKS], 1, __ATOMIC_RELAXED); if (ad < (int32_t)DSAME) __atomic_fetch_add(&cnt[OK20], 1, __ATOMIC_RELAXED); }
        __atomic_fetch_add(&cnt[SPAN], 1, __ATOMIC_RELAXED);
        if (ad < (int32_t)DSAME) __atomic_fetch_add(&cnt[SAME_], 1, __ATOMIC_RELAXED);
        while (__atomic_exchange_n(&w->lk, 1, __ATOMIC_ACQUIRE));
        w->dl[p->f].push_back(dl);
        __atomic_store_n(&w->lk, 0, __ATOMIC_RELEASE);
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
    for (size_t k = i; k < x->rs.size(); k += x->t) {
        fetch(x, x->rs[k], rec);
        cores(&x->sc[i], &c->r, rec.data(), LEV);
        lookup(x->p, c);
        x->p->f = x->rs[k] >> 40;
        chains(x->p, c, rec.data(), x->rs[k] & 0xffffffffff);
    }
}

static void add(std::vector<ev>& o, uint32_t q, const char* ty, std::string c1, uint32_t p1, std::string c2, uint32_t p2, const seg& x, const seg& y, uint32_t d = 0) {
    if (c1 == c2 && (p1 > p2 ? p1 - p2 : p2 - p1) < MS && std::strcmp(ty, "INS")) return;
    uint32_t l1 = x.qe - x.qs, l2 = y.qe - y.qs;
    if (c1 > c2 || (c1 == c2 && p1 > p2)) { std::swap(c1, c2); std::swap(p1, p2); std::swap(l1, l2); d = (d & 1) * 2 + d / 2; }
    o.push_back({ty, c1, c2, p1, p2, q, std::min(x.mq, y.mq), std::min(l1, l2), l1, l2, std::min(x.de, y.de), d, x.st});
}

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

static bool near(const sv& x, const sv& y) {
    return x.e.ty == y.e.ty && x.e.o == y.e.o && x.e.c1 == y.e.c1 && x.e.c2 == y.e.c2 && x.e.p1 + TOL >= y.e.p1 && y.e.p1 + TOL >= x.e.p1 && x.e.p2 + TOL >= y.e.p2 && y.e.p2 + TOL >= x.e.p2;
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

static float rep(const std::string& s) {
    float b = 0;
    for (uint32_t p = 1; p <= 6 && s.size() > p; p++) {
        uint32_t k = 0;
        for (uint32_t i = p; i < s.size(); i++) k += s[i] == s[i - p];
        b = std::max(b, (float)k / (s.size() - p));
    }
    return b;
}

static float prob(const sv& s, float rp, bool vn, float med, const float* w, float* f) {
    const arow& r = s.r;
    uint32_t t = top(r.n), dt = r.d[r.s], L = s.e.p2 - s.e.p1;
    int32_t al = std::abs(r.lest);
    bool sp = r.nspan > 0, ind = s.e.ty != "BND" && s.e.ty != "INV";
    f[0] = 1;
    f[1] = std::log2(t + 1.f);
    f[2] = std::min(1.f, (float)t / std::max(dt, 1u));
    f[3] = std::log2((float)s.pc);
    f[4] = r.ncl;
    f[5] = std::min(1.f, (float)t / std::max(r.tspan, 1u));
    f[6] = sp;
    f[7] = r.nsame > 0;
    f[8] = sp && ind && al < (int32_t)MS;
    f[9] = sp && al >= (int32_t)MS && ind && (2 * al < (int32_t)L || al > 2 * (int32_t)L);
    f[10] = r.nspan >= 3 && r.dmin < -20 && r.dmax > 20;
    f[11] = r.shadow > 0;
    f[12] = 1000.f * r.tonly / std::max(r.span, 1u);
    f[13] = 1000.f * r.nbr / std::max(r.span, 1u);
    f[14] = s.e.mq / 60.f;
    f[15] = std::min(s.e.de, 0.05f) * 100;
    f[16] = rp;
    f[17] = vn;
    f[18] = s.e.ty == "INS";
    f[19] = s.e.ty == "DUP";
    f[20] = s.e.ty == "INV";
    f[21] = s.e.ty == "BND";
    f[22] = s.e.ty == "BND" ? 0 : std::log2(L + 1.f);
    f[23] = std::log2((dt + 1.f) / (med + 1));
    f[24] = r.d[0] == 0;
    f[25] = (float)r.nsame / std::max(r.nspan, 1u);
    f[26] = std::log2(r.nsame + 1.f);
    f[27] = std::log2(r.nspan + 1.f);
    f[28] = f[5] - f[25];
    f[29] = s.ms > 0;
    f[30] = (float)s.ms / std::max(s.mn, 1u);
    f[31] = ind ? std::min(1.f, (float)s.nf / std::max(L, 1u)) : s.nf;
    f[32] = r.n20c > 0;
    f[33] = std::log2(r.n20c + 1.f);
    f[34] = std::log2(r.nspc + 1.f);
    f[35] = std::min(1.f, (float)r.tspc / std::max(r.tspan, 1u));
    float z = 0;
    for (uint32_t i = 0; i < NF; i++) z += w[i] * f[i];
    return 1 / (1 + std::exp(-z));
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

// records in reference order, for the duplicate sweep
struct by_place {
    const std::vector<sv>& r;
    explicit by_place(const std::vector<sv>& v) : r(v) {}
    bool operator()(size_t a, size_t b) const {
        const ev &x = r[a].e, &y = r[b].e;
        return x.c1 != y.c1 ? x.c1 < y.c1 : x.c2 != y.c2 ? x.c2 < y.c2 : x.p1 != y.p1 ? x.p1 < y.p1 : x.p2 != y.p2 ? x.p2 < y.p2 : a < b;
    }
};

// the evidence row of a site, from the pass
static arow row(const spot& e, uint32_t u) {
    const uint32_t *mn = e.m.data(), *mt = e.m.data() + e.s * NM;
    std::vector<int32_t> dn = e.dl[0];
    std::sort(dn.begin(), dn.end());
    int32_t lest = dn.empty() ? 0 : dn[dn.size() / 2], dmin = dn.empty() ? 0 : dn.front(), dmax = dn.empty() ? 0 : dn.back();
    uint32_t nsame = 0, ncl = 0;
    for (int32_t d : dn) nsame += std::labs(d) <= std::max<int32_t>(SAME, std::labs(lest) / 2);
    for (uint32_t c : clusters(e.dl[e.s])) ncl += c > 0;
    return {e.a, e.b, u, e.s, std::min(ncl, CL), mt[SAME_], mt[SPAN], mn[SPAN], nsame, 0, e.tonly, e.nbr, e.b - e.a, lest, dmin, dmax, mn[OKS], mn[OK20], mt[OKS], e.n, e.dep};
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
    float wt[NF];
    bed vn;
    std::copy(W0, W0 + NF, wt);
    if (vb) load_bed(vb, &vn);
    gph g;
    read_gph(gp, &g);
    uint64_t K = g.K, S = g.S;
    FILE* fi = std::fopen((pre + ".tsv").c_str(), "r");
    while (std::fgets(line, sizeof line, fi)) {
        std::sscanf(line, "%255s %u %u %*u %4095s %lu %u %4095s %4095s", nm, &a, &b, cn, &key, &s, ls, rs);
        spot r = {a, b, s, 0, 0, 0, 0, 0, 0, key, nums(cn), {}, {}, std::vector<uint32_t>(S * NM, 0), {}, std::vector<std::vector<int32_t>>(S), {}, nm, 0};
        r.t = r.n[s];
        for (uint32_t i = 0; i < NA; i++) for (uint32_t j = 0; j < NA; j++) r.d[i][j] = UNSET;
        unsigned long rs_, ss_;
        std::sscanf(nm, "s%lu_%lu_", &ss_, &rs_);
        r.read = ss_ << 40 | rs_;
        anchors(&r, ls, rs, key, &g);
        bypc[nm].push_back(w.size());
        w.push_back(r);
    }
    std::fclose(fi);
    // the pass: node flags, which nodes are anchors of which sites, the sites' tumor reads, then every read of every sample
    mctx x = {&w, {}, {}, std::vector<uint8_t>(K, 0), std::vector<uint8_t>(K, 0), std::vector<uint8_t>(K, 0), {}, std::vector<std::vector<uint64_t>>(t), {}, {}, std::vector<scr>(t), nullptr, t};
    for (uint64_t v = 0; v < K; v++) x.nz[v] = g.c[v] == 0;
    {
        std::vector<uint8_t> in(K, 0), on(K, 0);
        for (uint64_t e : g.e) { uint64_t u = e >> 32, v = e & 0xffffffff; if (on[u] < 2) on[u]++; if (in[v] < 2) in[v]++; }
        for (uint64_t v = 0; v < K; v++) x.br[v] = in[v] > 1 || on[v] > 1;
    }
    g.c = {};
    g.e = {};
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
    std::vector<uint32_t> dep;
    for (const spot& e : w) dep.push_back(e.dep[e.s]);
    std::sort(dep.begin(), dep.end());
    float med = dep.empty() ? 0 : dep[dep.size() / 2];
    // records from the pieces' alignments, each attached to the site with most tumor reads within ROW bp of it
    fi = std::fopen((pre + ".paf").c_str(), "r");
    while (std::fgets(line, sizeof line, fi)) {
        std::sscanf(line, "%255s %u %u %u %c %255s %u %u %u %u %u %u", nm, &z, &qs, &qe, &st, tn, &z, &ts, &te, &z, &z, &mq);
        const char *cg = std::strstr(line, "cg:Z:"), *de = std::strstr(line, "de:f:");
        al[nm].push_back({qs, qe, ts, te, mq, st, tn, cg ? std::string(cg + 5, std::strcspn(cg + 5, "\t\n")) : "", de ? std::strtof(de + 5, nullptr) : 0});
    }
    std::fclose(fi);
    std::unordered_map<std::string, std::vector<uint64_t>> fai;
    FILE* fx = std::fopen((std::string(fasta) + ".fai").c_str(), "r");
    unsigned long f4[4];
    while (std::fscanf(fx, "%255s %lu %lu %lu %lu", nm, &f4[0], &f4[1], &f4[2], &f4[3]) == 5) fai[nm] = {f4[0], f4[1], f4[2], f4[3]};
    std::fclose(fx);
    std::unordered_map<std::string, std::vector<ev>> evs;
    std::unordered_map<std::string, std::string> sq;
    for (auto& pc : bypc) if (al.count(pc.first)) for (const ev& e : evs[pc.first] = events(al[pc.first])) if (e.ty == "INS" && e.p2 - e.p1 >= NEAR) sq[pc.first];
    FILE *fa = std::fopen(fasta, "r"), *fp = std::fopen((pre + ".fa").c_str(), "r");
    char* ln = nullptr;
    size_t cap = 0;
    std::string* cur = nullptr;
    for (ssize_t n; (n = getline(&ln, &cap, fp)) > 0;)
        if (ln[0] == '>') { auto it = sq.find(std::string(ln + 1, n - 2)); cur = it == sq.end() ? nullptr : &it->second; }
        else if (cur) cur->assign(ln, n - 1);
    std::free(ln);
    std::fclose(fp);
    std::vector<sv> som;
    for (auto& pc : bypc) {
        if (!al.count(pc.first)) continue;
        for (const ev& e : evs[pc.first]) {
            const spot* r = nullptr;
            uint32_t u = 0;
            for (uint32_t i : pc.second) if (e.q + ROW >= w[i].a && e.q <= w[i].b + ROW && (!r || w[i].t > r->t)) { r = &w[i]; u = i; }
            if (!r || !primary(e.c1) || !primary(e.c2) || r->t < MM) continue;
            arow a_ = row(*r, u);
            som.push_back({e, pc.first, 1, a_, fit(a_, e), a_.nsame, a_.nspan});
        }
    }
    // one record per event: same type and orientation, both ends within TOL; the piece with most tumor reads represents it
    std::sort(som.begin(), som.end(), cmp_sv);
    std::vector<sv> res;
    for (const sv& v : som) {
        if (res.empty() || !near(res.back(), v)) { res.push_back(v); continue; }
        sv& q = res.back();
        arow m = top(v.r.n) > top(q.r.n) ? v.r : q.r;
        if (top(v.r.n) > top(q.r.n)) { q.e = v.e; q.nm = v.nm; }
        if (v.nf < q.nf || (v.nf == q.nf && v.ms > q.ms)) { q.nf = v.nf; q.ms = v.ms; q.mn = v.mn; }
        m.nsame = std::max(v.r.nsame, q.r.nsame);
        m.nspan = std::max(v.r.nspan, q.r.nspan);
        m.nspc = std::max(v.r.nspc, q.r.nspc);
        m.n20c = std::max(v.r.n20c, q.r.n20c);
        m.ncl = std::max(v.r.ncl, q.r.ncl);
        q.r = m;
        q.pc++;
    }
    FILE* o = std::fopen(out, "w");
    std::fprintf(o, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (uint64_t i = 0; i < S; i++) std::fprintf(o, "\tS%lu", (unsigned long)i);
    std::fprintf(o, "\n");
    std::vector<float> pv(res.size()), rv(res.size());
    std::vector<std::vector<float>> fv(res.size(), std::vector<float>(NF));
    std::vector<uint8_t> dup(res.size(), 0);
    std::vector<size_t> ix(res.size());
    for (size_t i = 0; i < res.size(); i++) {
        const sv& v = res[i];
        rv[i] = rep(ref(fa, fai[v.e.c1], v.e.p1 > WIN ? v.e.p1 - WIN : 0, std::min<uint64_t>(v.e.p1 + WIN, fai[v.e.c1][0])));
        pv[i] = prob(v, rv[i], inbed(&vn, v.e.c1, v.e.p1) || inbed(&vn, v.e.c2, v.e.p2), med, wt, fv[i].data());
        if (v.e.ty == "INS" && v.e.p2 - v.e.p1 >= NEAR && v.e.p1 > v.e.p2 - v.e.p1) tandem(res[i].e, sq[v.nm], fa, fai[v.e.c1]);
        ix[i] = i;
    }
    std::sort(ix.begin(), ix.end(), by_place(res));
    for (size_t a_ = 0; a_ < ix.size(); a_++)
        for (size_t b_ = a_ + 1; b_ < ix.size() && res[ix[b_]].e.c1 == res[ix[a_]].e.c1 && res[ix[b_]].e.c2 == res[ix[a_]].e.c2 && res[ix[b_]].e.p1 <= res[ix[a_]].e.p1 + DUPT; b_++)
            if (res[ix[b_]].e.p2 + DUPT >= res[ix[a_]].e.p2 && res[ix[a_]].e.p2 + DUPT >= res[ix[b_]].e.p2 && (res[ix[a_]].e.ty != res[ix[b_]].e.ty || res[ix[a_]].e.o == res[ix[b_]].e.o)) dup[pv[ix[a_]] < pv[ix[b_]] ? ix[a_] : ix[b_]] = 1;
    uint64_t np = 0;
    for (size_t i = 0; i < res.size(); i++) {
        const sv& v = res[i];
        const spot& e = w[v.r.u];
        float pr = pv[i];
        np += pr >= PCUT && !dup[i];
        std::string j = (v.e.o & 1 ? "]" : "[") + v.e.c2 + ":" + std::to_string(v.e.p2) + (v.e.o & 1 ? "]" : "["), alt = v.e.ty != "BND" ? "<" + v.e.ty + ">" : v.e.o & 2 ? "N" + j : j + "N", fs, dls;
        for (uint32_t k = 0; k < NF; k++) { std::snprintf(nm, sizeof nm, "%s%.4g", k ? "," : "", fv[i][k]); fs += nm; }
        for (uint64_t k = 0; k < S; k++) dls += ";D" + std::to_string(k) + "=" + deltas(e.dl[k]);
        std::fprintf(o, "%s\t%u\tsv%zu\tN\t%s\t.\t%s\tSVTYPE=%s;END=%u;SVLEN=%u;O=%u;PIECES=%u;PIECE=%s;ERR=%.4f;F=%s;SPAN=%s;SAME=%s%s\tGT:DV", v.e.c1.c_str(), v.e.p1, i, alt.c_str(), dup[i] ? "DUPL" : pr >= PCUT ? "PASS" : "LOW", v.e.ty.c_str(),
                     v.e.ty == "BND" || v.e.ty == "INS" ? v.e.p1 : v.e.p2, v.e.c1 == v.e.c2 ? v.e.p2 - v.e.p1 : 0, v.e.o, v.pc, v.nm.c_str(), 1 - pr, fs.c_str(), list(e.m, S, SPAN).c_str(), list(e.m, S, SAME_).c_str(), dls.c_str());
        for (uint32_t n : v.r.n) std::fprintf(o, "\t%s:%u", n >= MM ? "0/1" : "0/0", n);
        std::fprintf(o, "\n");
    }
    std::fclose(o);
    std::fclose(fa);
    std::fprintf(stderr, "[call] %zu tumor events at sites, %zu records, %lu PASS\n", som.size(), res.size(), (unsigned long)np);
    return 0;
}
