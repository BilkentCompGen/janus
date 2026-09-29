#include <algorithm>
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include "extract.h"
#include "place.h"

static const uint32_t WIN = 40000, MAR = 2000, TOL = 500, MINC = 3, EXPL = 100;

struct ent { uint64_t h; uint32_t pos, cs; };

struct seqs { std::vector<std::string> name, seq; };

struct hitr { uint32_t cs, st; int64_t d; uint32_t q, pos; };

struct win { uint32_t ci, lo, hi; };

struct rctx { const seqs* f; std::vector<win> job; std::vector<std::vector<ent>> out; std::vector<scr> s; std::vector<ent> idx; std::vector<uint8_t> on; std::vector<std::vector<win>> w; int t; };

static bool by_ent(const ent& a, const ent& b) { return a.h != b.h ? a.h < b.h : a.cs != b.cs ? a.cs < b.cs : a.pos < b.pos; }

static bool by_hit(const hitr& a, const hitr& b) { return a.cs != b.cs ? a.cs < b.cs : a.st != b.st ? a.st < b.st : a.d != b.d ? a.d < b.d : a.q < b.q; }

static bool by_win(const win& a, const win& b) { return a.ci != b.ci ? a.ci < b.ci : a.lo < b.lo; }

static void fasta(const char* path, seqs* f) {
    FILE* in = std::fopen(path, "r");
    char* ln = nullptr;
    size_t cap = 0;
    for (ssize_t n; (n = getline(&ln, &cap, in)) > 0;) {
        while (n && (ln[n - 1] == '\n' || ln[n - 1] == '\r')) n--;
        if (ln[0] == '>') { f->name.push_back(std::string(ln + 1, std::strcspn(ln + 1, " \t\n"))); f->seq.push_back(""); }
        else for (ssize_t i = 0; i < n; i++) f->seq.back() += ln[i] & 0xdf;
    }
    std::free(ln);
    std::fclose(in);
}

static void hashes(scr* s, const std::string& q, uint32_t lo, uint32_t hi, std::vector<uint64_t>& h) {
    std::vector<uint8_t> rec(4 + hi - lo + 8, 0);
    uint32_t tag = (hi - lo) | 1u << 31;
    std::memcpy(rec.data(), &tag, 4);
    std::memcpy(rec.data() + 4, q.data() + lo, hi - lo);
    cores(s, rec.data(), 5);
    std::vector<uint64_t> w(s->r.k / 32 + 2);
    h.resize(s->r.n);
    for (uint32_t i = 0; i < s->r.n; i++) h[i] = core_hash(w.data(), core_cut(&s->r, i, w.data()));
}

static uint32_t fwd(const rdcore* r, uint32_t i) { return i < r->nf ? r->a[i] : r->len - r->b[i]; }

static void ref_job(void* a, int i) {
    rctx* c = (rctx*)a;
    std::vector<uint64_t> h;
    for (size_t j = i; j < c->job.size(); j += c->t) {
        const win& w = c->job[j];
        const std::string& q = c->f->seq[w.ci];
        uint32_t lo = w.lo, e = w.hi;
        while (lo > 0 && w.lo - lo < MAR && q[lo - 1] != 'N') lo--;
        while (e < q.size() && e - w.hi < MAR && q[e] != 'N') e++;
        hashes(&c->s[i], q, lo, e, h);
        for (uint32_t k = 0; k < h.size(); k++) {
            uint32_t p = lo + fwd(&c->s[i].r, k);
            if (p >= w.lo && p < w.hi) c->out[i].push_back({h[k], p, w.ci << 1 | (k >= c->s[i].r.nf)});
        }
    }
}

int cmd_refidx(const char* in, const char* out, int t) {
    LCP_INIT();
    seqs f;
    fasta(in, &f);
    rctx c = {&f, {}, std::vector<std::vector<ent>>(t), std::vector<scr>(t), {}, {}, {}, t};
    for (uint32_t ci = 0; ci < f.seq.size(); ci++)
        for (uint32_t a = 0, b; a < f.seq[ci].size(); a = b) {
            while (a < f.seq[ci].size() && f.seq[ci][a] == 'N') a++;
            for (b = a; b < f.seq[ci].size() && f.seq[ci][b] != 'N'; b++);
            for (uint32_t w = a; w < b; w += WIN) c.job.push_back({ci, w, std::min(b, w + WIN)});
        }
    par(t, ref_job, &c);
    std::vector<ent> all;
    for (auto& v : c.out) all.insert(all.end(), v.begin(), v.end());
    std::sort(all.begin(), all.end(), by_ent);
    FILE* fo = open_tmp(out);
    std::fseek(fo, sizeof(bin_header), SEEK_SET);
    std::fwrite(all.data(), sizeof(ent), all.size(), fo);
    bin_header o = {{'S', 'V', 'B', 'R'}, 1, all.size(), f.seq.size(), {}};
    finish(fo, out, &o);
    std::fprintf(stderr, "[refidx] %zu cores on %zu sequences\n", all.size(), f.seq.size());
    return 0;
}

static void place_job(void* a, int i) {
    rctx* c = (rctx*)a;
    std::vector<uint64_t> h;
    std::vector<hitr> hs;
    for (size_t j = i; j < c->f->seq.size(); j += c->t) {
        const std::string& q = c->f->seq[j];
        hashes(&c->s[i], q, 0, q.size(), h);
        hs.clear();
        for (uint32_t k = 0; k < h.size(); k++) {
            auto r = std::equal_range(c->idx.begin(), c->idx.end(), ent{h[k], 0, 0}, [](const ent& x, const ent& y) { return x.h < y.h; });
            if (r.second - r.first != 1) continue;
            uint32_t st = (r.first->cs & 1) ^ (k >= c->s[i].r.nf), p = fwd(&c->s[i].r, k);
            hs.push_back({r.first->cs >> 1, st, st ? (int64_t)r.first->pos + p : (int64_t)r.first->pos - p, p, r.first->pos});
        }
        std::sort(hs.begin(), hs.end(), by_hit);
        for (size_t x = 0, y; x < hs.size(); x = y) {
            uint32_t lo = hs[x].pos, hi = hs[x].pos;
            for (y = x + 1; y < hs.size() && hs[y].cs == hs[x].cs && hs[y].st == hs[x].st && hs[y].d - hs[y - 1].d <= TOL; y++) { lo = std::min(lo, hs[y].pos); hi = std::max(hi, hs[y].pos); }
            if (y - x < MINC) continue;
            c->on[j] = 1;
            c->w[i].push_back({hs[x].cs, lo > q.size() + MAR ? lo - (uint32_t)q.size() - MAR : 0, hi + (uint32_t)q.size() + MAR});
        }
    }
}

static void write_fa(const char* path, const seqs& f, const std::vector<uint8_t>& on, uint8_t v) {
    FILE* fo = std::fopen(path, "w");
    for (size_t j = 0; j < f.seq.size(); j++) if (on[j] == v) std::fprintf(fo, ">%s\n%s\n", f.name[j].c_str(), f.seq[j].c_str());
    std::fclose(fo);
}

int cmd_place(const char* in, const char* idx, const char* ref, const char* out, int t) {
    LCP_INIT();
    seqs f, g;
    fasta(in, &f);
    bin_header h;
    FILE* fi = open_in(idx, &h);
    rctx c = {&f, {}, {}, std::vector<scr>(t), std::vector<ent>(h.nreads), std::vector<uint8_t>(f.seq.size(), 0), std::vector<std::vector<win>>(t), t};
    std::fread(c.idx.data(), sizeof(ent), h.nreads, fi);
    std::fclose(fi);
    par(t, place_job, &c);
    std::vector<win> w;
    for (auto& v : c.w) w.insert(w.end(), v.begin(), v.end());
    std::sort(w.begin(), w.end(), by_win);
    fasta(ref, &g);
    std::string o(out);
    FILE* fo = std::fopen((o + ".win.fa").c_str(), "w");
    for (size_t x = 0, y; x < w.size(); x = y) {
        uint32_t hi = w[x].hi;
        for (y = x + 1; y < w.size() && w[y].ci == w[x].ci && w[y].lo <= hi; y++) hi = std::max(hi, w[y].hi);
        hi = std::min<uint32_t>(hi, g.seq[w[x].ci].size());
        std::fprintf(fo, ">%s:%u\n%s\n", g.name[w[x].ci].c_str(), w[x].lo, g.seq[w[x].ci].substr(w[x].lo, hi - w[x].lo).c_str());
    }
    std::fclose(fo);
    write_fa((o + ".on.fa").c_str(), f, c.on, 1);
    write_fa((o + ".off.fa").c_str(), f, c.on, 0);
    uint64_t n = 0;
    for (uint8_t v : c.on) n += v;
    std::fprintf(stderr, "[place] %zu pieces, %lu placed by cores\n", f.seq.size(), n);
    return 0;
}

int cmd_lift(const char* paf, const char* on, const char* ref, const char* out) {
    std::map<std::string, uint64_t> len;
    char nm[256];
    unsigned long l, z;
    FILE* fi = std::fopen((std::string(ref) + ".fai").c_str(), "r");
    while (std::fscanf(fi, "%255s %lu %lu %lu %lu", nm, &l, &z, &z, &z) == 5) len[nm] = l;
    std::fclose(fi);
    std::map<std::string, std::vector<std::string>> rows;
    std::map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> cov;
    std::map<std::string, uint32_t> ql;
    fi = std::fopen(paf, "r");
    char* ln = nullptr;
    size_t cap = 0;
    for (ssize_t n; (n = getline(&ln, &cap, fi)) > 0;) {
        char q[256], s[4], tn[512];
        unsigned qn, qs, qe, tl, ts, te;
        int k;
        std::sscanf(ln, "%255s %u %u %u %3s %511s %u %u %u%n", q, &qn, &qs, &qe, s, tn, &tl, &ts, &te, &k);
        char* c = std::strrchr(tn, ':');
        *c = 0;
        uint32_t lo = std::strtoul(c + 1, nullptr, 10);
        rows[q].push_back(std::string(q) + "\t" + std::to_string(qn) + "\t" + std::to_string(qs) + "\t" + std::to_string(qe) + "\t" + s + "\t" + tn + "\t" + std::to_string(len[tn]) + "\t" + std::to_string(ts + lo) + "\t" + std::to_string(te + lo) + (ln + k));
        cov[q].push_back({qs, qe});
        ql[q] = qn;
    }
    std::free(ln);
    std::fclose(fi);
    seqs f;
    fasta(on, &f);
    std::string o(out);
    FILE *fp = std::fopen((o + ".paf").c_str(), "w"), *fm = std::fopen((o + ".more.fa").c_str(), "w");
    for (size_t j = 0; j < f.name.size(); j++) {
        std::vector<std::pair<uint32_t, uint32_t>>& v = cov[f.name[j]];
        std::sort(v.begin(), v.end());
        uint32_t got = 0, end = 0;
        for (auto& x : v) { if (x.second > end) got += x.second - std::max(x.first, end); end = std::max(end, x.second); }
        if (f.seq[j].size() - got >= EXPL) std::fprintf(fm, ">%s\n%s\n", f.name[j].c_str(), f.seq[j].c_str());
        else for (std::string& r : rows[f.name[j]]) std::fputs(r.c_str(), fp);
    }
    std::fclose(fp);
    std::fclose(fm);
    return 0;
}
