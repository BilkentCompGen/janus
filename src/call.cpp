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
#include "io.h"

// Every site from select is somatic by construction (its combination has no normal read and at least 3 tumor reads), so call
// only turns the pieces' alignments into SV records and applies the thresholds of the published long-read somatic callers:
// MS 50 bp minimum size (Severus, nanomonsv), MQ 20 mapping quality (nanomonsv's median-MAPQ threshold), VAF 0.05 minimum
// tumor VAF (Severus, nanomonsv), Sniffles2's coverage-scaled minimum support round(SB + SM * coverage) with coverage
// SR local + (1 - SR) genome-wide (its "auto" rule; local and genome-wide tumor depth come from the anchors), NCOV 3: the normal's zero counts only when the normal has more than 3 reads at the site
// (Severus' control coverage threshold), TOL 50 bp breakpoint clustering (Severus), and Severus' VNTR rule: events inside
// one VNTR (VB 25 bp buffer) are one event, and more than VN 3 distinct events in one VNTR fail as complex.
// Geometry: same-type gaps less than ISL aligned bases apart are one gap, an event belongs to a site whose interval is within
// ROW bp of it in the read, records within DUPT bp at both ends are duplicates (the benchmark's tolerance), and insertions of
// NEAR bp or more that copy the adjacent reference (KM-mer share between ROT and UNQ) are tandem duplications.
static const uint32_t MS = 50, MQ = 20, MM = 3, NCOV = 3, TOL = 50, ISL = 100, ROW = 200, DUPT = 500, NEAR = 500, KM = 15, VN = 3, VB = 25;
static const float VAF = 0.05f, SB = 1.5f, SM = 0.1f, SR = 0.75f, ROT = 0.8f, UNQ = 1.25f;

struct seg { uint32_t qs, qe, ts, te, mq; char st; std::string t, cg; };

struct ev { std::string ty, c1, c2; uint32_t p1, p2, q, mq, o; char st; };

struct spot { uint32_t a, b, s, t, nd, td; std::vector<uint32_t> n; };

struct sv { ev e; std::string nm, why; uint32_t pc, t, mq; spot r; };

struct bed { std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> m; };

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

// Depth of sample x at the site: the median graph count of its anchors on each side; nd is the normal's higher side and td the
// tumor's mean of the sides it has.
static void depths(spot* r, const char* ls, const char* rs, const std::vector<uint16_t>& c, uint64_t K) {
    std::vector<uint32_t> side[2][2];
    const char* p[2] = {ls, rs};
    for (int k = 0; k < 2; k++)
        for (const char* q = p[k]; *q && *q != '-' && *q != '\t' && *q != '\n';) {
            uint64_t v = std::strtoull(q, (char**)&q, 10);
            std::strtoul(q + 1, (char**)&q, 10);
            side[k][0].push_back(c[v]);
            side[k][1].push_back(c[r->s * K + v]);
            if (*q == ',') q++;
        }
    uint32_t n0 = median(side[0][0]), n1 = median(side[1][0]), t0 = median(side[0][1]), t1 = median(side[1][1]);
    bool h0 = !side[0][0].empty(), h1 = !side[1][0].empty();
    r->nd = std::max(n0, n1);
    r->td = h0 && h1 ? (t0 + t1) / 2 : h0 ? t0 : t1;
}

static bool cmp_seg(const seg& x, const seg& y) { return x.qs < y.qs; }

static bool cmp_sv(const sv& x, const sv& y) {
    return x.e.ty != y.e.ty ? x.e.ty < y.e.ty : x.e.c1 != y.e.c1 ? x.e.c1 < y.e.c1 : x.e.c2 != y.e.c2 ? x.e.c2 < y.e.c2 : x.e.o != y.e.o ? x.e.o < y.e.o : x.e.p1 != y.e.p1 ? x.e.p1 < y.e.p1 : x.e.p2 != y.e.p2 ? x.e.p2 < y.e.p2 : x.e.q != y.e.q ? x.e.q < y.e.q : x.nm < y.nm;
}

static bool cmp_pos(const sv& x, const sv& y) {
    return x.e.c1 != y.e.c1 ? x.e.c1 < y.e.c1 : x.e.c2 != y.e.c2 ? x.e.c2 < y.e.c2 : x.e.p1 != y.e.p1 ? x.e.p1 < y.e.p1 : x.e.p2 < y.e.p2;
}

static bool better(const sv& x, const sv& y) { return x.t != y.t ? x.t > y.t : x.mq != y.mq ? x.mq > y.mq : x.pc > y.pc; }

static void add(std::vector<ev>& o, uint32_t q, const char* ty, std::string c1, uint32_t p1, std::string c2, uint32_t p2, const seg& x, const seg& y, uint32_t d = 0) {
    if (c1 == c2 && (p1 > p2 ? p1 - p2 : p2 - p1) < MS && std::strcmp(ty, "INS")) return;
    if (c1 > c2 || (c1 == c2 && p1 > p2)) { std::swap(c1, c2); std::swap(p1, p2); d = (d & 1) * 2 + d / 2; }
    o.push_back({ty, c1, c2, p1, p2, q, std::min(x.mq, y.mq), d, x.st});
}

// Events of one piece: gaps of MS bp or more inside an alignment, and the jumps between consecutive alignments along the read.
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

// The VNTR holding both ends of an event on one chromosome (with VB bp of buffer), as its start; 0 when there is none.
static uint32_t vntr(const bed* b, const ev& e) {
    auto it = b->m.find(e.c1);
    if (e.c1 != e.c2 || e.ty == "INV" || it == b->m.end()) return 0;
    auto k = std::upper_bound(it->second.begin(), it->second.end(), std::make_pair(e.p1 + VB, ~0u));
    if (k == it->second.begin()) return 0;
    --k;
    return k->second + VB >= std::max(e.p1, e.p2) ? k->first + 1 : 0;
}

int cmd_call(const char* paf, const char* evi, const char* fasta, const char* gp, const char* out, const char* vb) {
    std::unordered_map<std::string, std::vector<spot>> sites;
    std::unordered_map<std::string, std::vector<seg>> al;
    char line[1 << 16], nm[256], cn[4096], ls[4096], rs[4096], st, tn[256];
    uint32_t qs, qe, ts, te, mq, z, a, b, s;
    bed vn;
    if (vb) load_bed(vb, &vn);
    bin_header h;
    FILE* fg = open_in(gp, &h);
    uint64_t K = h.nreads, S = h.resv[0], U;
    std::vector<uint16_t> c(S * K);
    std::fread(&U, 8, 1, fg);
    std::fseek(fg, 8 * K, SEEK_CUR);
    std::fread(c.data(), 2, S * K, fg);
    std::fclose(fg);
    FILE* fi = std::fopen(evi, "r");
    while (std::fgets(line, sizeof line, fi)) {
        std::sscanf(line, "%255s %u %u %*u %4095s %*s %u %4095s %4095s", nm, &a, &b, cn, &s, ls, rs);
        spot r = {a, b, s, 0, 0, 0, nums(cn)};
        r.t = r.n[s];
        depths(&r, ls, rs, c, K);
        sites[nm].push_back(r);
    }
    std::fclose(fi);
    c = {};
    fi = std::fopen(paf, "r");
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
    for (auto& pc : sites) {
        if (!al.count(pc.first)) continue;
        for (const ev& e : events(al[pc.first])) {
            const spot* r = nullptr;
            for (const spot& x : pc.second) if (e.q + ROW >= x.a && e.q <= x.b + ROW && (!r || x.t > r->t)) r = &x;
            if (!r || !primary(e.c1) || !primary(e.c2) || r->t < MM) continue;
            som.push_back({e, pc.first, "", 1, r->t, e.mq, *r});
            if (e.ty == "INS" && e.p2 - e.p1 >= NEAR) sq[pc.first];
        }
    }
    FILE *fa = std::fopen(fasta, "r"), *fp = std::fopen((std::string(evi, std::strlen(evi) - 3) + "fa").c_str(), "r");
    char* ln = nullptr;
    size_t cap = 0;
    std::string* cur = nullptr;
    for (ssize_t n; (n = getline(&ln, &cap, fp)) > 0;)
        if (ln[0] == '>') { auto it = sq.find(std::string(ln + 1, n - 2)); cur = it == sq.end() ? nullptr : &it->second; }
        else if (cur) cur->assign(ln, n - 1);
    std::free(ln);
    std::fclose(fp);
    for (sv& x : som) if (x.e.ty == "INS" && x.e.p2 - x.e.p1 >= NEAR && x.e.p1 > x.e.p2 - x.e.p1) tandem(x.e, sq[x.nm], fa, fai[x.e.c1]);
    std::fclose(fa);
    // one record per event: same type and orientation, both ends within TOL; the best piece (most tumor reads) represents it
    std::sort(som.begin(), som.end(), cmp_sv);
    std::vector<sv> res;
    for (const sv& x : som) {
        if (res.empty() || res.back().e.ty != x.e.ty || res.back().e.o != x.e.o || !near(res.back().e, x.e, TOL)) { res.push_back(x); continue; }
        sv& q = res.back();
        uint32_t pc = q.pc + 1, m = std::max(q.mq, x.mq);
        if (better(x, q)) q = x;
        q.pc = pc;
        q.mq = m;
    }
    std::vector<uint32_t> td;
    for (const sv& x : res) if (x.r.td) td.push_back(x.r.td);
    float gd = median(td);
    for (sv& x : res)
        x.why = x.r.nd <= NCOV ? "LOWNORMAL" : !x.r.td ? "NODEPTH" : (float)x.t / std::max(x.r.td, x.t) < VAF ? "LOWVAF" : x.t < std::lround(SB + SM * (SR * x.r.td + (1 - SR) * gd)) ? "LOWSUPPORT" : "PASS";
    // Severus' VNTR rule: same-type events inside one VNTR are one event; more than VN distinct events there fail as complex
    std::unordered_map<std::string, std::vector<uint32_t>> inv;
    for (uint32_t i = 0; i < res.size(); i++) {
        uint32_t v = vn.m.empty() ? 0 : vntr(&vn, res[i].e);
        if (v && res[i].why == "PASS") inv[res[i].e.c1 + ":" + std::to_string(v)].push_back(i);
    }
    for (auto& kv : inv) {
        std::vector<uint32_t>& g = kv.second;
        if (g.size() > VN) { for (uint32_t i : g) res[i].why = "VNTR"; continue; }
        for (uint32_t i = 0; i < g.size(); i++)
            for (uint32_t j = i + 1; j < g.size(); j++)
                if (res[g[i]].e.ty == res[g[j]].e.ty && res[g[i]].why == "PASS" && res[g[j]].why == "PASS") res[better(res[g[i]], res[g[j]]) ? g[j] : g[i]].why = "DUPL";
    }
    // records within DUPT bp at both ends describe one event; the better one stays
    std::sort(res.begin(), res.end(), cmp_pos);
    for (size_t i = 0; i < res.size(); i++)
        for (size_t j = i + 1; j < res.size() && res[j].e.c1 == res[i].e.c1 && res[j].e.p1 <= res[i].e.p1 + DUPT; j++)
            if (res[i].why == "PASS" && res[j].why == "PASS" && near(res[i].e, res[j].e, DUPT) && (res[i].e.ty != res[j].e.ty || res[i].e.o == res[j].e.o))
                res[better(res[i], res[j]) ? j : i].why = "DUPL";
    FILE* o = std::fopen(out, "w");
    std::fprintf(o, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (size_t i = 0; i < (res.empty() ? 0 : res[0].r.n.size()); i++) std::fprintf(o, "\tS%zu", i);
    std::fprintf(o, "\n");
    uint64_t np = 0;
    for (size_t i = 0; i < res.size(); i++) {
        const sv& x = res[i];
        np += x.why == "PASS";
        std::string j = (x.e.o & 1 ? "]" : "[") + x.e.c2 + ":" + std::to_string(x.e.p2) + (x.e.o & 1 ? "]" : "["), alt = x.e.ty != "BND" ? "<" + x.e.ty + ">" : x.e.o & 2 ? "N" + j : j + "N";
        std::fprintf(o, "%s\t%u\tsv%zu\tN\t%s\t.\t%s\tSVTYPE=%s;END=%u;SVLEN=%u;O=%u;PIECES=%u;PIECE=%s;QPOS=%u;SITE=%u-%u;MQ=%u;TDP=%u;NDP=%u;VAF=%.3f\tGT:DV", x.e.c1.c_str(), x.e.p1, i, alt.c_str(), x.why.c_str(), x.e.ty.c_str(),
                     x.e.ty == "BND" || x.e.ty == "INS" ? x.e.p1 : x.e.p2, x.e.c1 == x.e.c2 ? x.e.p2 - x.e.p1 : 0, x.e.o, x.pc, x.nm.c_str(), x.e.q, x.r.a, x.r.b, x.mq, x.r.td, x.r.nd, x.r.td ? (float)x.t / std::max(x.r.td, x.t) : 0);
        for (uint32_t n : x.r.n) std::fprintf(o, "\t%s:%u", n >= MM ? "0/1" : "0/0", n);
        std::fprintf(o, "\n");
    }
    std::fclose(o);
    std::fprintf(stderr, "[call] %zu events at sites, %zu records, %lu PASS\n", som.size(), res.size(), np);
    return 0;
}
