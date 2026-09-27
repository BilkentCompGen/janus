#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include "call.h"

static const uint32_t MM = 3, MQ = 10, MS = 50, TOL = 50, NEAR = 500, JN = 100, COV = 200, SHORT = 1000, ISL = 100, SUP = 5, SUPI = 8, GMQ = 30, WIN = 150, END = 500000, TINY = 80, SUPT = 15, ROW = 200;
static const float MDE = 0.015f, MREP = 0.5f;

struct seg { uint32_t qs, qe, ts, te, mq; char st; std::string t, cg; float de; };

struct ev { std::string ty, c1, c2; uint32_t p1, p2, q, mq, ln, l1, l2; float de; };

struct arow { uint32_t a, b, u; std::vector<uint32_t> n, d; std::string al; uint64_t k; };

struct crow { std::string pc; uint32_t s; };

struct sv { ev e; std::string nm; uint32_t f, pc, el; std::vector<uint32_t> n, d; std::vector<uint8_t> cov, has; };

static std::vector<uint32_t> nums(const char* s) {
    std::vector<uint32_t> v;
    for (char* p = (char*)s; *p && *p != '\t' && *p != '\n'; p += *p == ',') v.push_back(std::strtoul(p, &p, 10));
    return v;
}

static bool cmp_seg(const seg& x, const seg& y) { return x.qs < y.qs; }

static bool cmp_sv(const sv& x, const sv& y) {
    return x.e.ty != y.e.ty ? x.e.ty < y.e.ty : x.e.c1 != y.e.c1 ? x.e.c1 < y.e.c1 : x.e.c2 != y.e.c2 ? x.e.c2 < y.e.c2 : x.e.p1 < y.e.p1;
}

static void add(std::vector<ev>& o, uint32_t q, const char* ty, const std::string& c1, uint32_t p1, const std::string& c2, uint32_t p2, const seg& x, const seg& y) {
    if (c1 == c2 && (p1 > p2 ? p1 - p2 : p2 - p1) < MS && std::strcmp(ty, "INS")) return;
    uint32_t l1 = x.qe - x.qs, l2 = y.qe - y.qs;
    if (c1 == c2 && p1 > p2) { std::swap(p1, p2); std::swap(l1, l2); }
    o.push_back({ty, c1, c2, p1, p2, q, std::min(x.mq, y.mq), std::min(l1, l2), l1, l2, std::min(x.de, y.de)});
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
        uint32_t p1 = x.st == '+' ? x.te : x.ts, p2 = y.st == '+' ? y.ts : y.te, q = x.qe;
        int64_t g = (int64_t)y.qs - x.qe, d = x.st == '+' ? (int64_t)y.ts - x.te : (int64_t)x.ts - y.te;
        if (x.t != y.t) add(o, q, "BND", x.t, p1, y.t, p2, x, y);
        else if (x.st != y.st) add(o, q, "INV", x.t, p1, y.t, p2, x, y);
        else if (d - g >= MS) add(o, q, "DEL", x.t, p1, y.t, p2, x, y);
        else if (g - d >= MS && d > -(int64_t)MS) add(o, q, "INS", x.t, p1, x.t, p1 + g - d, x, y);
        else if (d <= -(int64_t)MS) add(o, q, "DUP", x.t, p1, y.t, p2, x, y);
    }
    return o;
}

static bool like(const ev& x, const ev& y) {
    uint32_t nr = x.ty == "DEL" || x.ty == "INS" ? std::min(NEAR, std::max(100u, 3 * (x.p2 - x.p1))) : NEAR;
    if (x.ty != y.ty || x.c1 != y.c1 || x.c2 != y.c2 || x.p1 + nr < y.p1 || y.p1 + nr < x.p1) return false;
    if (x.ty == "DEL" || x.ty == "INS") return 2 * (y.p2 - y.p1) >= x.p2 - x.p1 && y.p2 - y.p1 <= 2 * (x.p2 - x.p1);
    return x.p1 + JN >= y.p1 && y.p1 + JN >= x.p1 && x.p2 + JN >= y.p2 && y.p2 + JN >= x.p2;
}

static bool covers(const std::vector<seg>& sg, const std::string& c, uint32_t p) {
    for (const seg& s : sg) if (s.mq >= MQ && s.t == c && s.ts + COV <= p && s.te >= p + COV) return true;
    return false;
}

static bool near(const sv& x, const sv& y) {
    return x.e.ty == y.e.ty && x.e.c1 == y.e.c1 && x.e.c2 == y.e.c2 && x.e.p1 + TOL >= y.e.p1 && y.e.p1 + TOL >= x.e.p1 && x.e.p2 + TOL >= y.e.p2 && y.e.p2 + TOL >= x.e.p2;
}

static bool primary(const std::string& c) {
    return c.size() <= 5 && c.find('_') == std::string::npos && c != "chrM" && c != "M";
}

static uint32_t top(const std::vector<uint32_t>& n);

static std::string ref(FILE* f, const std::vector<uint64_t>& v, uint32_t a, uint32_t b) {
    uint64_t s = v[1] + a / v[2] * v[3] + a % v[2], e = v[1] + b / v[2] * v[3] + b % v[2];
    std::string r(e - s, 0), o;
    std::fseek(f, s, SEEK_SET);
    r.resize(std::fread(&r[0], 1, e - s, f));
    for (char ch : r) if (ch != '\n') o += ch & 0xdf;
    return o;
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

static const char* verdict(const sv& s, FILE* fa, std::unordered_map<std::string, std::vector<uint64_t>>& fai) {
    uint32_t t = top(s.n), L = s.e.p2 - s.e.p1;
    if (!s.cov[0]) return "NOREF";
    if (s.has[0]) return "NORMAL";
    if (t < SUP || (s.e.ty == "INS" && t < SUPI) || ((s.e.ty == "DEL" || s.e.ty == "INS") && L < TINY && t < SUPT) || (s.n[0] && (s.n[0] > 2 || t < 10 * s.n[0]))) return "LOWSUP";
    if (s.e.mq < GMQ) return "MAPQ";
    if (s.e.de > MDE) return "DIVERGE";
    if (t < 20 && s.e.ty != "BND" && s.e.ty != "INV" && L < SHORT && rep(ref(fa, fai[s.e.c1], s.e.p1 > WIN ? s.e.p1 - WIN : 0, std::min<uint64_t>(s.e.p1 + WIN, fai[s.e.c1][0]))) > MREP) return "REPEAT";
    const std::vector<uint64_t>& c = fai[s.e.c1];
    if ((s.e.ty == "DEL" || s.e.ty == "INS") && L < SHORT && (s.e.p1 < END || s.e.p1 + END > c[0])) return "CHREND";
    return "PASS";
}

static uint32_t top(const std::vector<uint32_t>& n) {
    return *std::max_element(n.begin() + 1, n.end());
}

int cmd_call(const char* paf, const char* evi, const char* fasta, const char* out) {
    std::unordered_map<std::string, std::vector<arow>> alt;
    std::unordered_map<uint32_t, std::vector<crow>> cp;
    std::unordered_map<uint64_t, std::vector<crow>> byk;
    std::unordered_map<std::string, std::vector<seg>> al;
    std::unordered_map<std::string, std::vector<ev>> me;
    char line[1 << 16], nm[256], role[8], cn[4096], x[4096], dd[4096], st, tn[256];
    uint32_t a, b, u, qs, qe, ts, te, mq, z;
    FILE* f = std::fopen(evi, "r");
    while (std::fgets(line, sizeof line, f)) {
        unsigned long ky = ~0ul;
        dd[0] = 0;
        std::sscanf(line, "%255s %u %u %7s %u %4095s %4095s %4095s %lu", nm, &a, &b, role, &u, cn, x, dd, &ky);
        if (role[0] != 'a') cp[u].push_back({nm, (uint32_t)std::atoi(cn)});
        else { alt[nm].push_back({a, b, u, nums(cn), nums(dd), x, ky}); byk[ky].push_back({nm, (uint32_t)std::atoi(nm + 1)}); }
    }
    std::fclose(f);
    f = std::fopen(paf, "r");
    while (std::fgets(line, sizeof line, f)) {
        std::sscanf(line, "%255s %u %u %u %c %255s %u %u %u %u %u %u", nm, &z, &qs, &qe, &st, tn, &z, &ts, &te, &z, &z, &mq);
        const char *cg = std::strstr(line, "cg:Z:"), *de = std::strstr(line, "de:f:");
        al[nm].push_back({qs, qe, ts, te, mq, st, tn, cg ? std::string(cg + 5, std::strcspn(cg + 5, "\t\n")) : "", de ? std::strtof(de + 5, nullptr) : 0});
    }
    std::fclose(f);
    std::vector<sv> som;
    for (auto& pc : alt) {
        uint32_t fs = std::atoi(pc.first.c_str() + 1);
        if (!fs || !al.count(pc.first)) continue;
        for (const ev& e : events(al[pc.first])) {
            const arow* r = nullptr;
            std::vector<const arow*> rs;
            for (const arow& w : pc.second)
                if (e.q + ROW >= w.a && e.q <= w.b + ROW) { rs.push_back(&w); if (!r || top(w.n) > top(r->n)) r = &w; }
            if (!r || !primary(e.c1) || !primary(e.c2) || top(r->n) < MM) continue;
            sv s = {e, pc.first, fs, 1, r->b - r->a, r->n, r->d, std::vector<uint8_t>(r->n.size(), 0), std::vector<uint8_t>(r->n.size(), 0)};
            std::vector<uint8_t> c2(r->n.size(), e.ty == "INS" || e.l2 < SHORT);
            for (uint32_t k = 0; k < r->n.size(); k++) s.cov[k] = e.l1 < SHORT;
            s.cov[fs] = c2[fs] = 1;
            for (const arow* w : rs) for (uint32_t k = 0; k < w->n.size(); k++) s.has[k] |= w->al[k] == 'a' || (w->al[k] == '1' && (e.ty == "DEL" || e.ty == "INS") && e.p2 - e.p1 < SHORT);
            std::vector<crow> ms;
            for (const arow& w : pc.second) ms.insert(ms.end(), cp[w.u].begin(), cp[w.u].end());
            for (const arow* w : rs) for (const crow& c : byk[w->k]) if (c.s != fs && w->k != ~0ul) ms.push_back(c);
            for (const crow& c : ms) {
                if (!al.count(c.pc)) continue;
                std::vector<seg>& sg = al[c.pc];
                s.cov[c.s] |= covers(sg, e.c1, e.p1);
                c2[c.s] |= covers(sg, e.c2, e.p2);
                if (!me.count(c.pc)) me[c.pc] = events(sg);
                for (const ev& y : me[c.pc]) s.has[c.s] |= like(e, y);
            }
            for (uint32_t k = 0; k < r->n.size(); k++) s.cov[k] &= c2[k];
            som.push_back(s);
        }
    }
    std::sort(som.begin(), som.end(), cmp_sv);
    std::vector<sv> res;
    for (const sv& s : som) {
        if (!res.empty() && near(res.back(), s)) {
            sv& r = res.back();
            for (size_t i = 0; i < s.n.size(); i++) { r.n[i] = std::max(r.n[i], s.n[i]); r.cov[i] |= s.cov[i]; r.has[i] |= s.has[i]; }
            if (r.d.size() == s.d.size()) for (size_t i = 0; i < s.d.size(); i++) r.d[i] = std::max(r.d[i], s.d[i]);
            r.pc++;
            continue;
        }
        res.push_back(s);
    }
    std::unordered_map<std::string, std::vector<uint64_t>> fai;
    FILE* fx = std::fopen((std::string(fasta) + ".fai").c_str(), "r");
    unsigned long fl_[4];
    while (std::fscanf(fx, "%255s %lu %lu %lu %lu", nm, &fl_[0], &fl_[1], &fl_[2], &fl_[3]) == 5) fai[nm] = {fl_[0], fl_[1], fl_[2], fl_[3]};
    std::fclose(fx);
    FILE *fa = std::fopen(fasta, "r"), *o = std::fopen(out, "w");
    std::fprintf(o, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (size_t i = 0; i < (res.empty() ? 0 : res[0].n.size()); i++) std::fprintf(o, "\tS%zu", i);
    std::fprintf(o, "\n");
    uint64_t np = 0, nr = 0;
    for (size_t i = 0; i < res.size(); i++) {
        const sv& s = res[i];
        const char* fl = verdict(s, fa, fai);
        np += !std::strcmp(fl, "PASS");
        nr += !std::strcmp(fl, "NOREF");
        std::string alt_ = s.e.ty == "BND" ? "N[" + s.e.c2 + ":" + std::to_string(s.e.p2) + "[" : "<" + s.e.ty + ">", cv, hs;
        std::string dp;
        for (size_t k = 0; k < s.n.size(); k++) { cv += (k ? "," : "") + std::to_string(s.cov[k]); hs += (k ? "," : "") + std::to_string(s.has[k]); }
        for (size_t k = 0; k < s.d.size(); k++) dp += (k ? "," : "") + std::to_string(s.d[k]);
        std::fprintf(o, "%s\t%u\tsv%zu\tN\t%s\t.\t%s\tSVTYPE=%s;END=%u;SVLEN=%u;PIECES=%u;PIECE=%s;EVLEN=%u;COV=%s;HAS=%s;MQ=%u;ALN=%u;DE=%.4f;DP=%s\tGT:DV", s.e.c1.c_str(), s.e.p1, i, alt_.c_str(), fl, s.e.ty.c_str(),
                     s.e.ty == "BND" || s.e.ty == "INS" ? s.e.p1 : s.e.p2, s.e.c1 == s.e.c2 ? s.e.p2 - s.e.p1 : 0, s.pc, s.nm.c_str(), s.el, cv.c_str(), hs.c_str(), s.e.mq, s.e.ln, s.e.de, dp.c_str());
        for (uint32_t n : s.n) std::fprintf(o, "\t%s:%u", n >= MM ? "0/1" : "0/0", n);
        std::fprintf(o, "\n");
    }
    std::fclose(o);
    std::fclose(fa);
    std::fprintf(stderr, "[call] %zu tumor events between anchors, %zu records: %lu PASS, %lu without a normal read across the breakpoints\n", som.size(), res.size(), np, nr);
    return 0;
}
