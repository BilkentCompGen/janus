#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include "call.h"

static const uint32_t MM = 3, MQ = 10, MS = 50, TOL = 50, NEAR = 500, COV = 200;

struct seg { uint32_t qs, qe, ts, te, mq; char st; std::string t, cg; };

struct ev { std::string ty, c1, c2; uint32_t p1, p2, q; };

struct arow { uint32_t a, b, u; std::vector<uint32_t> n; std::string al; };

struct crow { std::string pc; uint32_t s; };

struct sv { ev e; std::string nm; uint32_t f, pc, el; std::vector<uint32_t> n; std::vector<uint8_t> cov, has; };

static std::vector<uint32_t> nums(const char* s) {
    std::vector<uint32_t> v;
    for (char* p = (char*)s; *p && *p != '\t' && *p != '\n'; p += *p == ',') v.push_back(std::strtoul(p, &p, 10));
    return v;
}

static bool cmp_seg(const seg& x, const seg& y) { return x.qs < y.qs; }

static bool cmp_sv(const sv& x, const sv& y) {
    return x.e.ty != y.e.ty ? x.e.ty < y.e.ty : x.e.c1 != y.e.c1 ? x.e.c1 < y.e.c1 : x.e.c2 != y.e.c2 ? x.e.c2 < y.e.c2 : x.e.p1 < y.e.p1;
}

static void add(std::vector<ev>& o, uint32_t q, const char* ty, const std::string& c1, uint32_t p1, const std::string& c2, uint32_t p2) {
    if (c1 == c2 && p1 > p2) std::swap(p1, p2);
    o.push_back({ty, c1, c2, p1, p2, q});
}

static std::vector<ev> events(std::vector<seg>& sg) {
    std::vector<ev> o;
    std::sort(sg.begin(), sg.end(), cmp_seg);
    for (const seg& s : sg) {
        if (s.mq < MQ) continue;
        uint32_t t = s.ts, q = s.st == '+' ? s.qs : s.qe;
        for (const char* p = s.cg.c_str(); *p;) {
            uint32_t l = std::strtoul(p, (char**)&p, 10);
            char op = *p++;
            if (op == 'I' && l >= MS) add(o, s.st == '+' ? q : q - l, "INS", s.t, t, s.t, t + l);
            if (op == 'D' && l >= MS) add(o, q, "DEL", s.t, t, s.t, t + l);
            if (op != 'I') t += l;
            if (op != 'D') q = s.st == '+' ? q + l : q - l;
        }
    }
    for (size_t i = 0; i + 1 < sg.size(); i++) {
        const seg &x = sg[i], &y = sg[i + 1];
        if (x.mq < MQ || y.mq < MQ) continue;
        uint32_t p1 = x.st == '+' ? x.te : x.ts, p2 = y.st == '+' ? y.ts : y.te, q = x.qe;
        int64_t g = (int64_t)y.qs - x.qe, d = x.st == '+' ? (int64_t)y.ts - x.te : (int64_t)x.ts - y.te;
        if (x.t != y.t) add(o, q, "BND", x.t, p1, y.t, p2);
        else if (x.st != y.st) add(o, q, "INV", x.t, p1, y.t, p2);
        else if (d - g >= MS) add(o, q, "DEL", x.t, p1, y.t, p2);
        else if (g - d >= MS && d > -(int64_t)MS) add(o, q, "INS", x.t, p1, x.t, p1 + g - d);
        else if (d <= -(int64_t)MS) add(o, q, "DUP", x.t, p1, y.t, p2);
    }
    return o;
}

static bool like(const ev& x, const ev& y) {
    if (x.ty != y.ty || x.c1 != y.c1 || x.c2 != y.c2 || x.p1 + NEAR < y.p1 || y.p1 + NEAR < x.p1) return false;
    if (x.ty == "DEL" || x.ty == "INS") return 2 * (y.p2 - y.p1) >= x.p2 - x.p1 && y.p2 - y.p1 <= 2 * (x.p2 - x.p1);
    return x.p2 + NEAR >= y.p2 && y.p2 + NEAR >= x.p2;
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

static uint32_t top(const std::vector<uint32_t>& n) {
    return *std::max_element(n.begin() + 1, n.end());
}

int cmd_call(const char* paf, const char* evi, const char* out) {
    std::unordered_map<std::string, std::vector<arow>> alt;
    std::unordered_map<uint32_t, std::vector<crow>> cp;
    std::unordered_map<std::string, std::vector<seg>> al;
    std::unordered_map<std::string, std::vector<ev>> me;
    char line[1 << 16], nm[256], role[8], cn[4096], x[4096], st, tn[256];
    uint32_t a, b, u, qs, qe, ts, te, mq, z;
    FILE* f = std::fopen(evi, "r");
    while (std::fgets(line, sizeof line, f)) {
        std::sscanf(line, "%255s %u %u %7s %u %4095s %4095s", nm, &a, &b, role, &u, cn, x);
        if (role[0] == 'a') alt[nm].push_back({a, b, u, nums(cn), x});
        else cp[u].push_back({nm, (uint32_t)std::atoi(cn)});
    }
    std::fclose(f);
    f = std::fopen(paf, "r");
    while (std::fgets(line, sizeof line, f)) {
        std::sscanf(line, "%255s %u %u %u %c %255s %u %u %u %u %u %u", nm, &z, &qs, &qe, &st, tn, &z, &ts, &te, &z, &z, &mq);
        const char* cg = std::strstr(line, "cg:Z:");
        al[nm].push_back({qs, qe, ts, te, mq, st, tn, cg ? std::string(cg + 5, std::strcspn(cg + 5, "\t\n")) : ""});
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
                if (e.q + TOL >= w.a && e.q <= w.b + TOL) { rs.push_back(&w); if (!r || top(w.n) > top(r->n)) r = &w; }
            if (!r || !primary(e.c1) || !primary(e.c2) || top(r->n) < MM) continue;
            sv s = {e, pc.first, fs, 1, r->b - r->a, r->n, std::vector<uint8_t>(r->n.size(), 0), std::vector<uint8_t>(r->n.size(), 0)};
            std::vector<uint8_t> c2(r->n.size(), e.ty == "INS");
            s.cov[fs] = c2[fs] = 1;
            for (const arow* w : rs) {
                for (uint32_t k = 0; k < w->n.size(); k++) s.has[k] |= w->al[k] == '1' || w->al[k] == 'a';
                for (const crow& c : cp[w->u]) {
                    if (!al.count(c.pc)) continue;
                    std::vector<seg>& sg = al[c.pc];
                    s.cov[c.s] |= covers(sg, e.c1, e.p1);
                    c2[c.s] |= covers(sg, e.c2, e.p2);
                    if (!me.count(c.pc)) me[c.pc] = events(sg);
                    for (const ev& y : me[c.pc]) s.has[c.s] |= like(e, y);
                }
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
            r.pc++;
            continue;
        }
        res.push_back(s);
    }
    FILE* o = std::fopen(out, "w");
    std::fprintf(o, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (size_t i = 0; i < (res.empty() ? 0 : res[0].n.size()); i++) std::fprintf(o, "\tS%zu", i);
    std::fprintf(o, "\n");
    uint64_t np = 0, nr = 0;
    for (size_t i = 0; i < res.size(); i++) {
        const sv& s = res[i];
        const char* fl = !s.cov[0] ? "NOREF" : s.has[0] ? "NORMAL" : "PASS";
        np += !s.has[0] && s.cov[0];
        nr += !s.cov[0];
        std::string alt_ = s.e.ty == "BND" ? "N[" + s.e.c2 + ":" + std::to_string(s.e.p2) + "[" : "<" + s.e.ty + ">", cv, hs;
        for (size_t k = 0; k < s.n.size(); k++) { cv += (k ? "," : "") + std::to_string(s.cov[k]); hs += (k ? "," : "") + std::to_string(s.has[k]); }
        std::fprintf(o, "%s\t%u\tsv%zu\tN\t%s\t.\t%s\tSVTYPE=%s;END=%u;SVLEN=%u;PIECES=%u;PIECE=%s;EVLEN=%u;COV=%s;HAS=%s\tGT:DV", s.e.c1.c_str(), s.e.p1, i, alt_.c_str(), fl, s.e.ty.c_str(),
                     s.e.ty == "BND" || s.e.ty == "INS" ? s.e.p1 : s.e.p2, s.e.c1 == s.e.c2 ? s.e.p2 - s.e.p1 : 0, s.pc, s.nm.c_str(), s.el, cv.c_str(), hs.c_str());
        for (uint32_t n : s.n) std::fprintf(o, "\t%s:%u", n >= MM ? "0/1" : "0/0", n);
        std::fprintf(o, "\n");
    }
    std::fclose(o);
    std::fprintf(stderr, "[call] %zu tumor events between anchors, %zu records: %lu PASS, %lu without a normal read across the breakpoints, %lu the normal has\n", som.size(), res.size(), np, nr, res.size() - np - nr);
    return 0;
}
