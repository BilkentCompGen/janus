#include <algorithm>
#include <cstring>
#include "cmd.h"
#include "core.h"
#include "files.h"
#include "io.h"

static const uint32_t MAR = 2000;

struct Entry { uint64_t h; uint32_t pos, cs; };
struct Window { uint32_t ci, lo, hi; };

class Cores {
public:
    std::vector<uint64_t> h;
    void of(const std::string& q, uint32_t lo, uint32_t hi);
    uint32_t forward(uint32_t i) const { return r.forward(i); }
    bool reverse(uint32_t i) const { return i >= r.nf; }

private:
    Read r;
    Parser parser;
};

class RefIndex : public Job {
public:
    explicit RefIndex(const Args& a);
    int run();

private:
    const Args& a;
    Fasta f;
    std::vector<Window> job;
    std::vector<std::vector<Entry>> out;
    std::vector<Cores> cores;
    void task(int i) override;
};

class Place : public Job {
public:
    explicit Place(const Args& a);
    int run();

private:
    const Args& a;
    Fasta f;
    std::vector<Cores> cores;
    std::vector<Entry> idx;
    std::vector<uint8_t> on;
    std::vector<std::vector<Window>> w;
    void task(int i) override;
    void write(const char* path, uint8_t v) const;
};

static const uint32_t WIN = 40000;

void Cores::of(const std::string& q, uint32_t lo, uint32_t hi) {
    std::vector<uint8_t> rec(4 + hi - lo + 8, 0);
    uint32_t tag = (hi - lo) | 1u << 31;
    std::memcpy(rec.data(), &tag, 4);
    std::memcpy(rec.data() + 4, q.data() + lo, hi - lo);
    parser.cores(r, rec.data());
    h.resize(r.n);
    for (uint32_t i = 0; i < r.n; i++) h[i] = r.hash(i);
}

static bool by_entry(const Entry& a, const Entry& b) { return a.h != b.h ? a.h < b.h : a.cs != b.cs ? a.cs < b.cs : a.pos < b.pos; }

RefIndex::RefIndex(const Args& a) : a(a), f(a.in), out(a.threads), cores(a.threads) {
    for (uint32_t ci = 0; ci < f.seq.size(); ci++)
        for (uint32_t x = 0, y; x < f.seq[ci].size(); x = y) {
            while (x < f.seq[ci].size() && f.seq[ci][x] == 'N') x++;
            for (y = x; y < f.seq[ci].size() && f.seq[ci][y] != 'N'; y++);
            for (uint32_t w = x; w < y; w += WIN) job.push_back({ci, w, std::min(y, w + WIN)});
        }
}

void RefIndex::task(int i) {
    Cores& c = cores[i];
    for (size_t j = i; j < job.size(); j += a.threads) {
        const Window& w = job[j];
        const std::string& q = f.seq[w.ci];
        uint32_t lo = w.lo, e = w.hi;
        while (lo > 0 && w.lo - lo < MAR && q[lo - 1] != 'N') lo--;
        while (e < q.size() && e - w.hi < MAR && q[e] != 'N') e++;
        c.of(q, lo, e);
        for (uint32_t k = 0; k < c.h.size(); k++) {
            uint32_t p = lo + c.forward(k);
            if (p >= w.lo && p < w.hi) out[i].push_back({c.h[k], p, w.ci << 1 | c.reverse(k)});
        }
    }
}

int RefIndex::run() {
    par(a.threads);
    std::vector<Entry> all;
    for (std::vector<Entry>& v : out) all.insert(all.end(), v.begin(), v.end());
    std::sort(all.begin(), all.end(), by_entry);
    Output::save(a.out, Header('R', 1, all.size(), f.seq.size()), all.data(), all.size() * sizeof(Entry));
    std::fprintf(stderr, "[refidx] %zu cores on %zu sequences\n", all.size(), f.seq.size());
    return 0;
}

static const uint32_t TOL = 500;
static const uint32_t MINC = 3;

struct Match { uint32_t cs, st; int64_t d; uint32_t pos; };

static bool by_hash(const Entry& a, const Entry& b) { return a.h < b.h; }

static bool by_match(const Match& a, const Match& b) { return a.cs != b.cs ? a.cs < b.cs : a.st != b.st ? a.st < b.st : a.d < b.d; }

static bool by_window(const Window& a, const Window& b) { return a.ci != b.ci ? a.ci < b.ci : a.lo < b.lo; }

Place::Place(const Args& a) : a(a), f(a.in), cores(a.threads), on(f.seq.size(), 0), w(a.threads) {
    Input in(a.reads);
    idx.resize(in.h.n);
    std::fread(idx.data(), sizeof(Entry), in.h.n, in.f);
    std::fclose(in.f);
}

void Place::task(int i) {
    Cores& c = cores[i];
    std::vector<Match> hs;
    for (size_t j = i; j < f.seq.size(); j += a.threads) {
        const std::string& q = f.seq[j];
        c.of(q, 0, q.size());
        hs.clear();
        for (uint32_t k = 0; k < c.h.size(); k++) {
            std::pair<std::vector<Entry>::const_iterator, std::vector<Entry>::const_iterator> r = std::equal_range(idx.begin(), idx.end(), Entry{c.h[k], 0, 0}, by_hash);
            if (r.second - r.first != 1) continue;
            uint32_t st = (r.first->cs & 1) ^ c.reverse(k), p = c.forward(k);
            hs.push_back({r.first->cs >> 1, st, st ? (int64_t)r.first->pos + p : (int64_t)r.first->pos - p, r.first->pos});
        }
        std::sort(hs.begin(), hs.end(), by_match);
        for (size_t x = 0, y; x < hs.size(); x = y) {
            uint32_t lo = hs[x].pos, hi = hs[x].pos;
            for (y = x + 1; y < hs.size() && hs[y].cs == hs[x].cs && hs[y].st == hs[x].st && hs[y].d - hs[y - 1].d <= TOL; y++) { lo = std::min(lo, hs[y].pos); hi = std::max(hi, hs[y].pos); }
            if (y - x < MINC) continue;
            on[j] = 1;
            w[i].push_back({hs[x].cs, lo > q.size() + MAR ? lo - (uint32_t)q.size() - MAR : 0, hi + (uint32_t)q.size() + MAR});
        }
    }
}

void Place::write(const char* path, uint8_t v) const {
    FILE* fo = std::fopen(path, "w");
    for (size_t j = 0; j < f.seq.size(); j++)
        if (on[j] == v) Fasta::put(fo, f.name[j], f.seq[j]);
    std::fclose(fo);
}

int Place::run() {
    par(a.threads);
    std::vector<Window> all;
    for (std::vector<Window>& v : w) all.insert(all.end(), v.begin(), v.end());
    std::sort(all.begin(), all.end(), by_window);
    Fasta g(a.data);
    std::string o(a.out);
    FILE* fo = std::fopen((o + ".win.fa").c_str(), "w");
    for (size_t x = 0, y; x < all.size(); x = y) {
        uint32_t hi = all[x].hi;
        for (y = x + 1; y < all.size() && all[y].ci == all[x].ci && all[y].lo <= hi; y++) hi = std::max(hi, all[y].hi);
        hi = std::min<uint32_t>(hi, g.seq[all[x].ci].size());
        std::fprintf(fo, ">%s:%u\n%s\n", g.name[all[x].ci].c_str(), all[x].lo, g.seq[all[x].ci].substr(all[x].lo, hi - all[x].lo).c_str());
    }
    std::fclose(fo);
    write((o + ".on.fa").c_str(), 1);
    uint64_t n = 0;
    for (uint8_t v : on) n += v;
    std::fprintf(stderr, "[place] %zu pieces, %lu placed by cores\n", f.seq.size(), (unsigned long)n);
    return 0;
}

int cmd_refidx(const Args& a) {
    return RefIndex(a).run();
}

int cmd_place(const Args& a) {
    return Place(a).run();
}
