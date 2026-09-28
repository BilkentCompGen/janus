#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include "call.h"

static const uint32_t MM = 3, MQ = 10, MS = 50, MSI = 45, TOL = 50, NEAR = 500, JN = 100, ISL = 250, ROW = 200, WIN = 150, NF = 32, DUPT = 500, KM = 15;
static const float PCUT = 0.5f, W0[NF] = {-6, 1.5f, 1, 0.5f, -0.3f, 1, 0.5f, -4, -3, -1.5f, -1.5f, -1, -0.05f, -0.05f, 3, -0.5f, -1.5f, -1, 0, 0, 0, 0, 0, -0.2f, -1, 0, 0, 0, 0, 0, 0, 0};

struct seg { uint32_t qs, qe, ts, te, mq; char st; std::string t, cg; float de; };

struct ev { std::string ty, c1, c2; uint32_t p1, p2, q, mq, ln, l1, l2; float de; char st; uint32_t syn; };

struct arow { uint32_t a, b, u, s, ncl, cl, tspan, nspan, nsame, shadow, tonly, nbr, span, tsame, orphan; int32_t lest, dmin, dmax; std::vector<uint32_t> n, d, lc, rc; };

struct sv { ev e; std::string nm; uint32_t pc; arow r; };

struct bed { std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> m; };

static std::vector<uint32_t> nums(const char* s) {
    std::vector<uint32_t> v;
    for (char* p = (char*)s; *p && *p != '\t' && *p != '\n'; p += *p == ',') v.push_back(std::strtoul(p, &p, 10));
    return v;
}

static bool cmp_seg(const seg& x, const seg& y) { return x.qs < y.qs; }

static bool cmp_sv(const sv& x, const sv& y) {
    return x.e.ty != y.e.ty ? x.e.ty < y.e.ty : x.e.c1 != y.e.c1 ? x.e.c1 < y.e.c1 : x.e.c2 != y.e.c2 ? x.e.c2 < y.e.c2 : x.e.p1 < y.e.p1;
}

static uint32_t top(const std::vector<uint32_t>& n) {
    return *std::max_element(n.begin() + 1, n.end());
}

static void add(std::vector<ev>& o, uint32_t q, const char* ty, std::string c1, uint32_t p1, std::string c2, uint32_t p2, const seg& x, const seg& y, uint32_t ms) {
    if (c1 == c2 && (p1 > p2 ? p1 - p2 : p2 - p1) < ms && std::strcmp(ty, "INS")) return;
    uint32_t l1 = x.qe - x.qs, l2 = y.qe - y.qs;
    if (c1 > c2 || (c1 == c2 && p1 > p2)) { std::swap(c1, c2); std::swap(p1, p2); std::swap(l1, l2); }
    o.push_back({ty, c1, c2, p1, p2, q, std::min(x.mq, y.mq), std::min(l1, l2), l1, l2, std::min(x.de, y.de), x.st, 0});
}

static uint32_t refpos(const seg& s, uint32_t q) {
    uint32_t t = s.ts, x = s.st == '+' ? s.qs : s.qe;
    for (const char* p = s.cg.c_str(); *p;) {
        uint32_t l = std::strtoul(p, (char**)&p, 10);
        char op = *p++;
        if (op != 'D' && (s.st == '+' ? x + l > q : x - l <= q)) return op == 'I' ? t : t + (s.st == '+' ? q - x : x - q);
        if (op != 'I') t += l;
        if (op != 'D') x = s.st == '+' ? x + l : x - l;
    }
    return t;
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
            if ((op == 'I' || op == 'D') && l >= MSI) {
                if (op == gop && gap < ISL) gl += l + (op == 'D' ? gap : 0);
                else {
                    if (gop) add(o, gq, gop == 'I' ? "INS" : "DEL", s.t, gt, s.t, gt + gl, s, s, MSI);
                    gop = op; gt = t; gq = op == 'I' && s.st != '+' ? q - l : q; gl = l;
                }
                gap = 0;
            } else if (op != 'I' && op != 'D') gap += l;
            if (op != 'I') t += l;
            if (op != 'D') q = s.st == '+' ? q + l : q - l;
        }
        if (gop) add(o, gq, gop == 'I' ? "INS" : "DEL", s.t, gt, s.t, gt + gl, s, s, MSI);
    }
    for (size_t i = 0; i + 2 < sg.size(); i++) {
        const seg &x = sg[i], &y = sg[i + 2];
        int64_t d = x.st == '+' ? (int64_t)y.ts - x.te : (int64_t)x.ts - y.te;
        if (x.mq >= MQ && y.mq >= MQ && x.t == y.t && x.st == y.st && d > -(int64_t)MS && d < (int64_t)MS && y.qs > x.qe + MS) add(o, x.qe, "INS", x.t, x.st == '+' ? x.te : x.ts, x.t, (x.st == '+' ? x.te : x.ts) + y.qs - x.qe, x, y, MS);
    }
    for (size_t i = 0; i + 1 < sg.size(); i++) {
        const seg &x = sg[i], &y = sg[i + 1];
        if (x.mq < MQ || y.mq < MQ) continue;
        uint32_t p1 = x.st == '+' ? x.te : x.ts, p2 = y.st == '+' ? y.ts : y.te, q = x.qe;
        int64_t g = (int64_t)y.qs - x.qe, d = x.st == '+' ? (int64_t)y.ts - x.te : (int64_t)x.ts - y.te;
        if (x.t != y.t) add(o, q, "BND", x.t, p1, y.t, p2, x, y, MS);
        else if (x.st != y.st) add(o, q, "INV", x.t, p1, y.t, p2, x, y, MS);
        else if (d - g >= MS) add(o, q, "DEL", x.t, p1, y.t, p2, x, y, MS);
        else if (g - d >= MS && d > -(int64_t)MS) add(o, q, "INS", x.t, p1, x.t, p1 + g - d, x, y, MS);
        else if (d <= -(int64_t)MS) add(o, q, "DUP", x.t, p1, y.t, p2, x, y, MS);
    }
    return o;
}

static bool near(const sv& x, const sv& y) {
    return x.e.ty == y.e.ty && x.e.c1 == y.e.c1 && x.e.c2 == y.e.c2 && x.e.p1 + TOL >= y.e.p1 && y.e.p1 + TOL >= x.e.p1 && x.e.p2 + TOL >= y.e.p2 && y.e.p2 + TOL >= x.e.p2;
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

struct pfa { FILE* f; std::unordered_map<std::string, long> at; };

static void index_fa(const char* path, pfa* x) {
    x->f = std::fopen(path, "r");
    char* line = nullptr;
    size_t cap = 0;
    std::string nm;
    for (ssize_t n; (n = getline(&line, &cap, x->f)) > 0;)
        if (line[0] == '>') { nm.assign(line + 1, n - 2); x->at[nm] = std::ftell(x->f); }
    std::free(line);
}

static std::string piece(pfa* x, const std::string& nm, uint32_t lo, uint32_t n, bool rc) {
    std::string s(n, 0), o;
    std::fseek(x->f, x->at[nm] + lo, SEEK_SET);
    s.resize(std::fread(&s[0], 1, n, x->f));
    if (!rc) return s;
    for (size_t i = s.size(); i--;) o += s[i] == 'A' ? 'T' : s[i] == 'C' ? 'G' : s[i] == 'G' ? 'C' : s[i] == 'T' ? 'A' : 'N';
    return o;
}

static float shared(const std::string& a, const std::string& b) {
    std::unordered_map<uint64_t, uint8_t> k;
    uint32_t n = 0, m = 0;
    for (size_t i = 0; i + KM <= b.size(); i++) k[std::hash<std::string>()(b.substr(i, KM))] = 1;
    for (size_t i = 0; i + KM <= a.size(); i++) { n++; m += k.count(std::hash<std::string>()(a.substr(i, KM))); }
    return n ? (float)m / n : 0;
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

static uint32_t rank(const arow& x, const ev& e, bool ind) {
    uint32_t L = e.p2 - e.p1, al = std::abs(x.lest);
    return !x.nspan ? 1 : ind && 2 * al >= L && al <= 2 * L ? 2 : 0;
}

static float prob(const sv& s, float rp, bool vn, float med, const float* w, float* f) {
    const arow& r = s.r;
    uint32_t t = top(r.n), dt = std::min(r.lc[r.s], r.rc[r.s]), L = s.e.p2 - s.e.p1;
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
    f[24] = std::min(r.lc[0], r.rc[0]) == 0;
    f[25] = (float)r.nsame / std::max(std::min(r.nspan, 32u), 1u);
    f[26] = std::log2(r.nsame + 1.f);
    f[27] = std::log2(r.nspan + 1.f);
    f[28] = f[5] - f[25];
    f[29] = (float)t / (t + r.tsame);
    f[30] = r.orphan;
    f[31] = s.e.syn;
    float z = 0;
    for (uint32_t i = 0; i < NF; i++) z += w[i] * f[i];
    return 1 / (1 + std::exp(-z));
}

int cmd_call(const char* paf, const char* evi, const char* fasta, const char* out, const char* wf, const char* vb) {
    std::unordered_map<std::string, std::vector<arow>> alt;
    std::unordered_map<std::string, std::vector<seg>> al;
    char line[1 << 16], nm[256], cn[4096], dd[4096], lc[4096], rc[4096], st, tn[256];
    uint32_t qs, qe, ts, te, mq, z;
    float w[NF];
    bed vn;
    std::vector<uint32_t> dep;
    std::copy(W0, W0 + NF, w);
    if (wf) { FILE* fw = std::fopen(wf, "r"); for (uint32_t i = 0; i < NF; i++) std::fscanf(fw, "%f", &w[i]); std::fclose(fw); }
    if (vb) load_bed(vb, &vn);
    FILE* fi = std::fopen(evi, "r");
    while (std::fgets(line, sizeof line, fi)) {
        arow r;
        std::sscanf(line, "%255s %u %u %u %4095s %4095s %*s %u %u %u %u %u %u %d %d %d %u %u %u %u %4095s %4095s %u", nm, &r.a, &r.b, &r.u, cn, dd, &r.s, &r.ncl, &r.cl, &r.tspan, &r.nspan, &r.nsame, &r.lest, &r.dmin, &r.dmax, &r.shadow, &r.tonly, &r.nbr, &r.span, lc, rc, &r.tsame);
        r.orphan = 0;
        r.n = nums(cn);
        r.d = nums(dd);
        r.lc = nums(lc);
        r.rc = nums(rc);
        dep.push_back(std::min(r.lc[r.s], r.rc[r.s]));
        alt[nm].push_back(r);
    }
    std::fclose(fi);
    std::sort(dep.begin(), dep.end());
    float med = dep.empty() ? 0 : dep[dep.size() / 2];
    fi = std::fopen(paf, "r");
    while (std::fgets(line, sizeof line, fi)) {
        std::sscanf(line, "%255s %u %u %u %c %255s %u %u %u %u %u %u", nm, &z, &qs, &qe, &st, tn, &z, &ts, &te, &z, &z, &mq);
        const char *cg = std::strstr(line, "cg:Z:"), *de = std::strstr(line, "de:f:");
        al[nm].push_back({qs, qe, ts, te, mq, st, tn, cg ? std::string(cg + 5, std::strcspn(cg + 5, "\t\n")) : "", de ? std::strtof(de + 5, nullptr) : 0});
    }
    std::fclose(fi);
    std::unordered_map<std::string, std::vector<uint64_t>> fai;
    FILE* fx = std::fopen((std::string(fasta) + ".fai").c_str(), "r");
    unsigned long fl_[4];
    while (std::fscanf(fx, "%255s %lu %lu %lu %lu", nm, &fl_[0], &fl_[1], &fl_[2], &fl_[3]) == 5) fai[nm] = {fl_[0], fl_[1], fl_[2], fl_[3]};
    std::fclose(fx);
    FILE* fa = std::fopen(fasta, "r");
    pfa pf;
    index_fa((std::string(evi).substr(0, std::string(evi).size() - 4) + ".fa").c_str(), &pf);
    std::vector<sv> som;
    for (auto& pc : alt) {
        if (!al.count(pc.first)) continue;
        std::vector<ev> es = events(al[pc.first]);
        std::vector<uint8_t> used(pc.second.size(), 0);
        for (ev& e : es) {
            const arow* r = nullptr;
            bool ind = e.ty == "INS" || e.ty == "DEL" || e.ty == "DUP";
            uint32_t L = e.p2 - e.p1, wn = ind ? std::max(ROW, 2 * L) : ROW;
            for (const arow& x : pc.second) if (e.q + wn >= x.a && e.q <= x.b + wn && (!r || rank(x, e, ind) > rank(*r, e, ind) || (rank(x, e, ind) == rank(*r, e, ind) && top(x.n) > top(r->n)))) r = &x;
            arow o;
            if (!r) {
                const arow* q = nullptr;
                for (const arow& x : pc.second) if (!q || std::min(x.a > e.q ? x.a - e.q : 0u, e.q > x.b ? e.q - x.b : 0u) < std::min(q->a > e.q ? q->a - e.q : 0u, e.q > q->b ? e.q - q->b : 0u)) q = &x;
                o = *q; o.orphan = 1; o.nspan = o.nsame = o.ncl = o.cl = 0; o.lest = o.dmin = o.dmax = 0; r = &o;
            } else used[r - pc.second.data()] = 1;
            uint32_t al_ = std::abs(r->lest);
            if (((e.ty == "DEL" && r->lest > 0) || (e.ty == "INS" && r->lest < 0)) && !r->nsame && al_ >= MS && (al_ > L ? al_ - L : L - al_) <= 20) { e.p2 = e.p1 + al_; L = al_; }
            if (ind && L < MS) continue;
            if (e.ty == "INS" && L >= 1000 && e.p1 > L) {
                std::string in = piece(&pf, pc.first, e.q, L, e.st != '+');
                if (shared(in, ref(fa, fai[e.c1], e.p1 - L, e.p1)) >= 0.5f) { e.ty = "DUP"; e.p2 = e.p1; e.p1 -= L; }
                else if (shared(in, ref(fa, fai[e.c1], e.p1, std::min<uint64_t>(e.p1 + L, fai[e.c1][0]))) >= 0.5f) e.ty = "DUP";
            }
            if (!primary(e.c1) || !primary(e.c2) || top(r->n) < MM) continue;
            som.push_back({e, pc.first, 1, *r});
        }
        for (size_t i = 0; i < pc.second.size(); i++) {
            const arow& x = pc.second[i];
            uint32_t al_ = std::abs(x.lest);
            if (used[i] || x.nsame || al_ < MS || x.cl < 3 || x.nspan < 3 || top(x.n) < MM) continue;
            const seg* g = nullptr;
            for (const seg& q : al[pc.first]) if (q.mq >= MQ && q.qs <= x.a && x.b <= q.qe) g = &q;
            if (!g || !primary(g->t)) continue;
            uint32_t p1 = refpos(*g, x.lest > 0 ? x.a : x.b);
            som.push_back({{x.lest > 0 ? "DEL" : "INS", g->t, g->t, p1, p1 + al_, x.a, g->mq, g->qe - g->qs, g->qe - g->qs, g->qe - g->qs, g->de, g->st, 1}, pc.first, 1, x});
        }
    }
    std::sort(som.begin(), som.end(), cmp_sv);
    std::vector<sv> res;
    for (const sv& s : som) {
        if (res.empty() || !near(res.back(), s)) { res.push_back(s); continue; }
        sv& q = res.back();
        if (top(s.r.n) > top(q.r.n)) { q.e = s.e; q.nm = s.nm; q.r = s.r; }
        q.pc++;
    }
    FILE* o = std::fopen(out, "w");
    std::fprintf(o, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (size_t i = 0; i < (res.empty() ? 0 : res[0].r.n.size()); i++) std::fprintf(o, "\tS%zu", i);
    std::fprintf(o, "\n");
    std::vector<float> pv(res.size()), rv(res.size());
    std::vector<std::vector<float>> fv(res.size(), std::vector<float>(NF));
    std::vector<uint8_t> dup(res.size(), 0);
    std::vector<size_t> ix(res.size());
    for (size_t i = 0; i < res.size(); i++) {
        const sv& s = res[i];
        rv[i] = rep(ref(fa, fai[s.e.c1], s.e.p1 > WIN ? s.e.p1 - WIN : 0, std::min<uint64_t>(s.e.p1 + WIN, fai[s.e.c1][0])));
        pv[i] = prob(s, rv[i], inbed(&vn, s.e.c1, s.e.p1) || inbed(&vn, s.e.c2, s.e.p2), med, w, fv[i].data());
        ix[i] = i;
    }
    std::sort(ix.begin(), ix.end(), [&](size_t a, size_t b) { const ev &x = res[a].e, &y = res[b].e; return x.c1 != y.c1 ? x.c1 < y.c1 : x.c2 != y.c2 ? x.c2 < y.c2 : x.p1 < y.p1; });
    for (size_t a = 0; a < ix.size(); a++)
        for (size_t b = a + 1; b < ix.size() && res[ix[b]].e.c1 == res[ix[a]].e.c1 && res[ix[b]].e.c2 == res[ix[a]].e.c2 && res[ix[b]].e.p1 <= res[ix[a]].e.p1 + DUPT; b++)
            if (res[ix[b]].e.ty == res[ix[a]].e.ty && res[ix[b]].e.p2 + DUPT >= res[ix[a]].e.p2 && res[ix[a]].e.p2 + DUPT >= res[ix[b]].e.p2) dup[pv[ix[a]] < pv[ix[b]] ? ix[a] : ix[b]] = 1;
    uint64_t np = 0;
    for (size_t i = 0; i < res.size(); i++) {
        const sv& s = res[i];
        float p = pv[i];
        np += p >= PCUT && !dup[i];
        std::string alt_ = s.e.ty == "BND" ? "N[" + s.e.c2 + ":" + std::to_string(s.e.p2) + "[" : "<" + s.e.ty + ">", fs;
        for (uint32_t k = 0; k < NF; k++) { std::snprintf(nm, sizeof nm, "%s%.4g", k ? "," : "", fv[i][k]); fs += nm; }
        std::fprintf(o, "%s\t%u\tsv%zu\tN\t%s\t.\t%s\tSVTYPE=%s;END=%u;SVLEN=%u;PIECES=%u;PIECE=%s;ERR=%.4f;F=%s\tGT:DV", s.e.c1.c_str(), s.e.p1, i, alt_.c_str(), dup[i] ? "DUPL" : p >= PCUT ? "PASS" : "LOW", s.e.ty.c_str(),
                     s.e.ty == "BND" || s.e.ty == "INS" ? s.e.p1 : s.e.p2, s.e.c1 == s.e.c2 ? s.e.p2 - s.e.p1 : 0, s.pc, s.nm.c_str(), 1 - p, fs.c_str());
        for (uint32_t n : s.r.n) std::fprintf(o, "\t%s:%u", n >= MM ? "0/1" : "0/0", n);
        std::fprintf(o, "\n");
    }
    std::fclose(o);
    std::fclose(fa);
    std::fprintf(stderr, "[call] %zu tumor events at sites, %zu records, %lu PASS\n", som.size(), res.size(), np);
    return 0;
}
