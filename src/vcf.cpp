#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_set>
#include "call.h"

static const uint32_t MS = 50;
static const uint32_t MQ = 10;
static const uint32_t ISL = 100;

static void add(std::vector<Event>& o, uint32_t q, const char* ty, std::string c1, uint32_t p1, std::string c2, uint32_t p2, const Align& x, const Align& y, uint32_t d = 0) {
    if (c1 == c2 && (p1 > p2 ? p1 - p2 : p2 - p1) < MS && std::strcmp(ty, "INS")) return;
    if (c1 > c2 || (c1 == c2 && p1 > p2)) { std::swap(c1, c2); std::swap(p1, p2); d = (d & 1) * 2 + d / 2; }
    o.push_back({ty, c1, c2, p1, p2, q, std::min(x.mq, y.mq), std::min(x.de, y.de), d, x.st, 0});
}

static bool by_start(const Align& x, const Align& y) {
    return x.qs < y.qs;
}

static void runs(std::vector<Event>& o, const Align& s) {
    uint32_t t = s.ts, q = s.st == '+' ? s.qs : s.qe, gt = 0, gq = 0, gl = 0, gap = ISL;
    char gop = 0;
    for (const char* p = s.cg.c_str(); *p;) {
        uint32_t l = std::strtoul(p, (char**)&p, 10);
        char op = *p++;
        if ((op == 'I' || op == 'D') && l >= MS) {
            if (op == gop && gap < ISL) gl += l + (op == 'D' ? gap : 0);
            else {
                if (gop) add(o, gq, gop == 'I' ? "INS" : "DEL", s.t, gt, s.t, gt + gl, s, s);
                gop = op;
                gt = t;
                gq = op == 'I' && s.st != '+' ? q - l : q;
                gl = l;
            }
            gap = 0;
        } else if (op != 'I' && op != 'D') gap += l;
        if (op != 'I') t += l;
        if (op != 'D') q = s.st == '+' ? q + l : q - l;
    }
    if (gop) add(o, gq, gop == 'I' ? "INS" : "DEL", s.t, gt, s.t, gt + gl, s, s);
}

std::vector<Event> events(std::vector<Align>& sg) {
    std::vector<Event> o;
    std::sort(sg.begin(), sg.end(), by_start);
    for (const Align& s : sg)
        if (s.mq >= MQ) runs(o, s);
    for (size_t i = 0; i + 1 < sg.size(); i++) {
        const Align &x = sg[i], &y = sg[i + 1];
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

static const int32_t SAMEMIN = 20;
static const int32_t NEARD = 10;
static const int32_t NEARP = 3;

static uint32_t clusters(const std::vector<int32_t>& d) {
    uint32_t n = 0;
    for (size_t i = 0, j; i < d.size(); i = j, n++)
        for (j = i; j < d.size() && (j == i || d[j] - d[j - 1] <= (int32_t)DSAME); j++);
    return n;
}

Row::Row(const Spot& e, uint32_t u) : a(e.a), b(e.b), u(u), s(e.s), n(e.n), d(e.dep) {
    const uint32_t *mn = e.m.data(), *mt = e.m.data() + e.s * KINDS;
    const std::vector<int32_t>& dn = e.dl[0];
    lest = dn.empty() ? 0 : dn[dn.size() / 2];
    dmin = dn.empty() ? 0 : dn.front();
    dmax = dn.empty() ? 0 : dn.back();
    nsame = 0;
    for (int32_t x : dn) nsame += std::labs(x) <= std::max<int32_t>(SAMEMIN, std::labs(lest) / 2);
    ncl = std::min(clusters(e.dl[e.s]), CL);
    tspan = mt[SPAN];
    nspan = mn[SPAN];
    tonly = e.tonly;
    nbr = e.nbr;
    span = e.b - e.a;
    nspc = mn[CONSISTENT];
    n20c = mn[CONSISTENT_SAME];
    tspc = mt[CONSISTENT];
}

uint32_t misfit(const Row& r, const Event& e) {
    uint32_t o = e.q < r.a || e.q > r.b, L = e.p2 - e.p1;
    return e.ty == "BND" || e.ty == "INV" ? o : o ? 2 * L : r.nspan ? std::abs(std::abs(r.lest) - (int32_t)L) : L;
}

float repeats(const std::string& s) {
    float b = 0;
    for (uint32_t p = 1; p <= 6 && s.size() > p; p++) {
        uint32_t k = 0;
        for (uint32_t i = p; i < s.size(); i++) k += s[i] == s[i - p];
        b = std::max(b, (float)k / (s.size() - p));
    }
    return b;
}

static float lg(float x) {
    return std::log2(x + 1);
}

static uint32_t pick(const Probe* b, uint32_t s) {
    const uint32_t P = NA * NA;
    size_t top[2] = {0, 0}, got = P, best = 0;
    for (uint32_t x = 0; x < P; x++) {
        top[0] = std::max(top[0], b->dl[x].size());
        top[1] = std::max(top[1], b->dl[s * P + x].size());
        if (b->dl[x].size() + b->dl[s * P + x].size() > b->dl[best].size() + b->dl[s * P + best].size()) best = x;
    }
    for (uint32_t t = 0; t < 2 * NA - 1 && got == P; t++)
        for (uint32_t x = 0; x <= t && x < NA && got == P; x++)
            if (t - x < NA && 5 * b->dl[x * NA + t - x].size() >= 4 * top[0] && 5 * b->dl[s * P + x * NA + t - x].size() >= 4 * top[1]) got = x * NA + t - x;
    return got == P ? best : got;
}

static void measured(const Probe* b, uint32_t s, float L, float* f) {
    for (uint32_t i = 0; i < 20; i++) f[i] = 0;
    if (!b || !b->set) return;
    const uint32_t P = NA * NA, got = pick(b, s);
    const std::vector<std::pair<int32_t, int32_t>> &dn = b->dl[got], &dt = b->dl[s * P + got];
    float n = dn.size(), t = dt.size(), c[2][3] = {{0, 0, 0}, {0, 0, 0}};
    int32_t near = INT32_MAX, path = 30;
    for (uint32_t k = 0; k < 2; k++)
        for (const std::pair<int32_t, int32_t>& z : k ? dt : dn) {
            bool d = std::abs(z.first) <= NEARD, p = z.second <= NEARP;
            c[k][0] += d;
            c[k][1] += p;
            c[k][2] += d && p;
            if (!k) { near = std::min(near, std::abs(z.first)); path = std::min(path, z.second); }
        }
    L = std::max(L, 1.f);
    f[0] = 1;
    f[1] = b->set == 2;
    f[2] = lg(n);
    f[3] = lg(t);
    f[4] = c[1][0] / std::max(t, 1.f);
    f[5] = lg(c[0][0]);
    f[6] = c[0][0] / std::max(n, 1.f);
    if (n > 0) {
        f[7] = std::min(near / L, 2.f);
        f[8] = lg(near);
        f[9] = lg(dn.back().first - dn.front().first);
        f[10] = std::min(std::abs(near - L) / L, 2.f);
        f[14] = path;
    }
    f[11] = lg(c[0][1]);
    f[12] = c[0][1] / std::max(n, 1.f);
    f[13] = c[1][1] / std::max(t, 1.f);
    f[15] = lg(c[0][2]);
    f[16] = c[1][2] / std::max(t, 1.f) - c[0][2] / std::max(n, 1.f);
    f[17] = lg(b->gap);
    f[18] = lg(b->absent);
    f[19] = lg(b->mid.size());
}

void features(const Record& s, const Probe* b, float rp, bool vn, float med, float* f) {
    const Row& r = s.r;
    uint32_t t = r.n[r.s], dt = r.d[r.s], L = s.e.p2 - s.e.p1;
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
    f[11] = 0;
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
    measured(b, r.s, ind ? L : 0, f + 36);
}

Bed::Bed(const char* path) {
    if (!path) return;
    char c[256];
    uint32_t s, e;
    FILE* f = std::fopen(path, "r");
    while (std::fscanf(f, "%255s %u %u%*[^\n]", c, &s, &e) == 3) m[c].push_back({s, e});
    std::fclose(f);
    for (std::pair<const std::string, std::vector<std::pair<uint32_t, uint32_t>>>& kv : m) {
        std::vector<std::pair<uint32_t, uint32_t>>& v = kv.second;
        std::sort(v.begin(), v.end());
        size_t j = 0;
        for (size_t i = 1; i < v.size(); i++)
            if (v[i].first <= v[j].second) v[j].second = std::max(v[j].second, v[i].second);
            else v[++j] = v[i];
        v.resize(j + 1);
    }
}

bool Bed::near(const std::string& c, uint32_t p) const {
    uint32_t s, e;
    return span(c, p, &s, &e);
}

bool Bed::span(const std::string& c, uint32_t p, uint32_t* s, uint32_t* e) const {
    std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>>::const_iterator it = m.find(c);
    if (it == m.end()) return false;
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator k = std::upper_bound(it->second.begin(), it->second.end(), std::make_pair(p + 50, ~0u));
    if (k == it->second.begin() || (k - 1)->second + 50 < p) return false;
    *s = (k - 1)->first;
    *e = (k - 1)->second;
    return true;
}

static const uint32_t WIN = 150;
static const uint32_t DUPT = 500;
static const uint32_t KM = 15;
static const float ROT = 0.8f;
static const float UNQ = 1.25f;

static float shared(const std::string& a, const std::string& b) {
    std::unordered_set<std::string> k;
    uint32_t m = 0;
    for (size_t i = 0; i + KM <= b.size(); i++) k.insert(b.substr(i, KM));
    for (size_t i = 0; i + KM <= a.size(); i++) m += k.count(a.substr(i, KM));
    return (float)m / (a.size() - KM + 1);
}

static void tandem(Event& e, const std::string& pc, Fai& fai, FILE* fa) {
    uint32_t L = e.p2 - e.p1;
    std::string s = pc.substr(e.q, L), o;
    if (e.st != '+') {
        for (size_t i = s.size(); i--;) o += s[i] == 'A' ? 'T' : s[i] == 'C' ? 'G' : s[i] == 'G' ? 'C' : s[i] == 'T' ? 'A' : 'N';
        s = o;
    }
    float u = shared(s, fai.range(fa, e.c1, e.p1 - L, e.p1)), d = shared(s, fai.range(fa, e.c1, e.p1, std::min<uint64_t>(e.p1 + L, fai.m[e.c1].len)));
    if (u + d < ROT || u + d > UNQ) return;
    e.ty = "DUP";
    e.p1 -= (uint32_t)(u / (u + d) * L);
    e.p2 = e.p1 + L;
}

static std::string list(const std::vector<uint32_t>& m, uint32_t S, uint32_t kind) {
    std::string s;
    for (uint32_t i = 0; i < S; i++) s += (i ? "," : "") + std::to_string(m[i * KINDS + kind]);
    return s;
}

static std::string deltas(const std::vector<int32_t>& d) {
    std::string s;
    for (int32_t v : d) s += (s.empty() ? "" : ",") + std::to_string(v);
    return s.empty() ? "-" : s;
}

struct Placed { const Event* e; size_t i; };

static bool by_place(const Placed& a, const Placed& b) {
    const Event &x = *a.e, &y = *b.e;
    return x.c1 != y.c1 ? x.c1 < y.c1 : x.c2 != y.c2 ? x.c2 < y.c2 : x.p1 != y.p1 ? x.p1 < y.p1 : x.p2 != y.p2 ? x.p2 < y.p2 : a.i < b.i;
}


void Call::write(std::vector<Record>& res, size_t found) {
    Fai fai(a.need(a.ref, "--ref"));
    FILE* fa = std::fopen(a.ref, "r");
    uint64_t S = w[0].m.size() / KINDS;
    FILE* o = std::fopen(a.need(a.out, "--out"), "w");
    std::fprintf(o, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (uint64_t i = 0; i < S; i++) std::fprintf(o, "\tS%lu", (unsigned long)i);
    std::fprintf(o, "\n");
    std::vector<std::vector<float>> fv(res.size(), std::vector<float>(NF));
    std::vector<uint8_t> dup(res.size(), 0);
    std::vector<Placed> ix(res.size());
    for (size_t i = 0; i < res.size(); i++) {
        const Record& v = res[i];
        float rv = repeats(fai.range(fa, v.e.c1, v.e.p1 > WIN ? v.e.p1 - WIN : 0, std::min<uint64_t>(v.e.p1 + WIN, fai.m[v.e.c1].len)));
        features(v,v.e.pb ? &pbs[v.e.pb - 1] : nullptr, rv, bed.near(v.e.c1, v.e.p1) || bed.near(v.e.c2, v.e.p2), med[v.r.s], fv[i].data());
        if (v.e.ty == "INS" && v.e.p2 - v.e.p1 >= NEAR && v.e.p1 > v.e.p2 - v.e.p1 && sq.count(v.nm)) tandem(res[i].e, sq.at(v.nm), fai, fa);
        ix[i] = {&res[i].e, i};
    }
    std::sort(ix.begin(), ix.end(), by_place);
    for (size_t x = 0; x < ix.size(); x++) {
        const Event& p = *ix[x].e;
        for (size_t y = x + 1; y < ix.size() && ix[y].e->c1 == p.c1 && ix[y].e->c2 == p.c2 && ix[y].e->p1 <= p.p1 + DUPT; y++) {
            const Event& q = *ix[y].e;
            if (res[ix[x].i].r.s == res[ix[y].i].r.s && q.p2 + DUPT >= p.p2 && p.p2 + DUPT >= q.p2 && (p.ty != q.ty || p.o == q.o)) dup[res[ix[x].i].r.n[res[ix[x].i].r.s] < res[ix[y].i].r.n[res[ix[y].i].r.s] ? ix[x].i : ix[y].i] = 1;
        }
    }
    uint64_t nd = 0;
    char num[64];
    for (size_t i = 0; i < res.size(); i++) {
        const Record& v = res[i];
        const Spot& e = w[v.r.u];
        nd += dup[i];
        std::string j = (v.e.o & 1 ? "]" : "[") + v.e.c2 + ":" + std::to_string(v.e.p2) + (v.e.o & 1 ? "]" : "["), alt = v.e.ty != "BND" ? "<" + v.e.ty + ">" : v.e.o & 2 ? "N" + j : j + "N", fs, dls;
        for (uint32_t k = 0; k < NF; k++) { std::snprintf(num, sizeof num, "%s%.4g", k ? "," : "", fv[i][k]); fs += num; }
        for (uint64_t k = 0; k < S; k++) dls += ";D" + std::to_string(k) + "=" + deltas(e.dl[k]);
        std::fprintf(o, "%s\t%u\tsv%zu\tN\t%s\t.\t%s\tSVTYPE=%s;END=%u;SVLEN=%u;O=%u;PIECES=%u;PIECE=%s;F=%s;SPAN=%s;SAME=%s%s\tGT:DV", v.e.c1.c_str(), v.e.p1, i, alt.c_str(), dup[i] ? "DUPL" : ".", v.e.ty.c_str(),
                     v.e.ty == "BND" || v.e.ty == "INS" ? v.e.p1 : v.e.p2, v.e.c1 == v.e.c2 ? v.e.p2 - v.e.p1 : 0, v.e.o, v.pc, v.nm.c_str(), fs.c_str(), list(e.m, S, SPAN).c_str(), list(e.m, S, SAME).c_str(), dls.c_str());
        for (uint64_t k = 0; k < v.r.n.size(); k++) std::fprintf(o, "\t%s:%u", (k ? k == v.r.s : v.r.n[k] >= a.support) ? "0/1" : "0/0", v.r.n[k]);
        std::fprintf(o, "\n");
    }
    std::fclose(o);
    std::fclose(fa);
    std::fprintf(stderr, "[call] %zu tumor events at sites, %zu records, %lu duplicates\n", found, res.size(), (unsigned long)nd);
}
