#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include "cmd.h"
#include "extract.h"
#include "files.h"
#include "io.h"

static const uint32_t WIN = 40000; // reference bases per refidx job
static const uint32_t MAR = 2000; // bases parsed past a refidx window on each side, and added around a placed piece
static const uint32_t TOL = 500; // largest diagonal step between neighbouring core hits of one placement
static const uint32_t MINC = 3; // unique core hits a placement needs

struct ent { uint64_t h; uint32_t pos, cs; };

struct hitr { uint32_t cs, st; int64_t d; uint32_t q, pos; };

struct win { uint32_t ci, lo, hi; };

struct rctx { const seqs* f; std::vector<win> job; std::vector<std::vector<ent>> out; std::vector<scr> s; int t; };

struct wctx { const seqs* f; std::vector<scr> s; std::vector<ent> idx; std::vector<uint8_t> on; std::vector<std::vector<win>> w; int t; };

static bool by_hash(const ent& a, const ent& b) { return a.h < b.h; }

static bool by_ent(const ent& a, const ent& b) { return a.h != b.h ? a.h < b.h : a.cs != b.cs ? a.cs < b.cs : a.pos < b.pos; }

static bool by_hit(const hitr& a, const hitr& b) { return a.cs != b.cs ? a.cs < b.cs : a.st != b.st ? a.st < b.st : a.d != b.d ? a.d < b.d : a.q < b.q; }

static bool by_win(const win& a, const win& b) { return a.ci != b.ci ? a.ci < b.ci : a.lo < b.lo; }

static void hashes(scr* s, const std::string& q, uint32_t lo, uint32_t hi, std::vector<uint64_t>& h) {
    std::vector<uint8_t> rec(4 + hi - lo + 8, 0);
    uint32_t tag = (hi - lo) | 1u << 31;
    std::memcpy(rec.data(), &tag, 4);
    std::memcpy(rec.data() + 4, q.data() + lo, hi - lo);
    cores(s, &s->r, rec.data(), LEV);
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
    read_fasta(in, &f);
    rctx c = {&f, {}, std::vector<std::vector<ent>>(t), std::vector<scr>(t), t};
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
    bin_header o = header('R', 1, all.size(), f.seq.size());
    save(out, &o, all.data(), all.size() * sizeof(ent));
    std::fprintf(stderr, "[refidx] %zu cores on %zu sequences\n", all.size(), f.seq.size());
    return 0;
}

static void place_job(void* a, int i) {
    wctx* c = (wctx*)a;
    std::vector<uint64_t> h;
    std::vector<hitr> hs;
    for (size_t j = i; j < c->f->seq.size(); j += c->t) {
        const std::string& q = c->f->seq[j];
        hashes(&c->s[i], q, 0, q.size(), h);
        hs.clear();
        for (uint32_t k = 0; k < h.size(); k++) {
            auto r = std::equal_range(c->idx.begin(), c->idx.end(), ent{h[k], 0, 0}, by_hash);
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
    read_fasta(in, &f);
    bin_header h;
    FILE* fi = open_in(idx, &h);
    wctx c = {&f, std::vector<scr>(t), std::vector<ent>(h.nreads), std::vector<uint8_t>(f.seq.size(), 0), std::vector<std::vector<win>>(t), t};
    std::fread(c.idx.data(), sizeof(ent), h.nreads, fi);
    std::fclose(fi);
    par(t, place_job, &c);
    std::vector<win> w;
    for (auto& v : c.w) w.insert(w.end(), v.begin(), v.end());
    std::sort(w.begin(), w.end(), by_win);
    read_fasta(ref, &g);
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
