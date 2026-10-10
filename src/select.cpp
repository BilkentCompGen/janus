#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>
#include "cmd.h"
#include "core.h"
#include "files.h"
#include "io.h"
#include "keys.h"
#include "reads.h"

struct Anchor { double sc; uint32_t v, a, b; };

struct Site {
    uint64_t read, key;
    uint32_t s, lo, hi, a, b, pc, lv, rv;
    std::string n;
    std::vector<Anchor> l, r;
};

class Anchors {
public:
    explicit Anchors(const Graph& g);
    void place(const Chain& c, const std::vector<uint32_t>& ids, uint64_t s, const Hit& h, uint32_t m, uint32_t rc, Site& e) const;

private:
    const Graph& g;
    std::vector<uint32_t> med;
    double fit(const std::vector<uint32_t>& ids, uint32_t v, uint64_t s) const;
    void side(const Chain& c, const std::vector<uint32_t>& ids, uint64_t s, int64_t from, int64_t step, int64_t end, std::vector<Anchor>& o) const;
};

static const uint32_t PAD = 1000;
static const uint32_t FL = 3;
static const uint32_t MULT = 2;
static const double DL = 15000;

Anchors::Anchors(const Graph& g) : g(g) {
    std::vector<uint32_t> smp;
    for (uint64_t s = 0; s < g.S; s++) {
        smp.clear();
        for (uint64_t v = 0; v < g.K; v++) {
            bool with = s && g.reads(s, v) > 0;
            for (uint64_t x = 1; x < g.S && !s && !with; x++) with = g.reads(x, v) > 0;
            if (with && g.reads(0, v)) smp.push_back(g.reads(s, v));
        }
        std::nth_element(smp.begin(), smp.begin() + smp.size() / 2, smp.end());
        med.push_back(smp[smp.size() / 2]);
    }
}

double Anchors::fit(const std::vector<uint32_t>& ids, uint32_t v, uint64_t s) const {
    std::pair<std::vector<uint32_t>::const_iterator, std::vector<uint32_t>::const_iterator> r = std::equal_range(ids.begin(), ids.end(), v);
    uint32_t cn = g.reads(0, v), ct = g.reads(s, v);
    if (!cn || !ct || r.second - r.first != 1 || cn > MULT * med[0] || ct > MULT * med[s]) return 0;
    double n = (double)cn / med[0], t = (double)ct / med[s];
    return std::sqrt(std::min(1.0, n) * std::min(1.0, t));
}

static bool by_cost(const Anchor& p, const Anchor& q) {
    return p.sc < q.sc;
}

void Anchors::side(const Chain& c, const std::vector<uint32_t>& ids, uint64_t s, int64_t from, int64_t step, int64_t end, std::vector<Anchor>& o) const {
    uint32_t at = c.start(from);
    for (int64_t j = from; j != end; j += step) {
        uint32_t p = c.start(j);
        double w = std::max(0.0, 1 - (p > at ? p - at : at - p) / DL), f;
        if (w <= 0 || (o.size() == NA && w <= -o.back().sc)) break;
        if ((f = fit(ids, c.id[j], s)) == 0) continue;
        o.push_back({-(f * w), c.id[j], p, c.end(j)});
        std::sort(o.begin(), o.end(), by_cost);
        if (o.size() > NA) o.pop_back();
    }
    for (Anchor& x : o) x.sc = x.a > at ? x.a - at : at - x.a;
    std::sort(o.begin(), o.end(), by_cost);
}

static void flip(uint32_t len, uint32_t rc, uint32_t* a, uint32_t* b) {
    if (!rc) return;
    uint32_t z = len - *b;
    *b = len - *a;
    *a = z;
}

void Anchors::place(const Chain& c, const std::vector<uint32_t>& ids, uint64_t s, const Hit& q, uint32_t m, uint32_t rc, Site& e) const {
    uint32_t len = c.r.len;
    e.l.clear();
    e.r.clear();
    side(c, ids, s, q.a, -1, -1, e.l);
    side(c, ids, s, q.b, 1, m, e.r);
    e.lv = e.l.empty() ? NONE : e.l[0].v;
    e.rv = e.r.empty() ? NONE : e.r[0].v;
    e.a = e.l.empty() ? c.start(q.a) : e.l[0].a;
    e.b = e.r.empty() ? c.end(q.b) : e.r[0].b;
    uint32_t pd = std::max(PAD, FL * (e.b - e.a)), lo = e.l.empty() ? e.a : e.l.back().a, hi = e.r.empty() ? e.b : e.r.back().b;
    e.lo = lo > pd ? lo - pd : 0;
    e.hi = std::min(len, hi + pd);
    flip(len, rc, &e.lo, &e.hi);
    flip(len, rc, &e.a, &e.b);
    for (Anchor& x : e.l) flip(len, rc, &x.a, &x.b);
    for (Anchor& x : e.r) flip(len, rc, &x.a, &x.b);
    if (rc) e.l.swap(e.r);
}

struct Seen { uint32_t site, d; uint64_t sc; };

class Groups {
public:
    Groups(const std::vector<Site>& w, const Combos& combos, int t);
    void measure(Chain& c, uint32_t m, uint64_t rid, uint32_t sample);
    std::vector<Site> rows(uint32_t support);

private:
    const std::vector<Site>& w;
    const Combos& combos;
    std::unordered_map<uint64_t, std::vector<uint32_t>> bykey;
    std::vector<std::vector<Seen>> seen;
};

Groups::Groups(const std::vector<Site>& w, const Combos& combos, int t) : w(w), combos(combos), seen(t) {
    for (uint32_t i = 0; i < w.size(); i++) bykey[w[i].key].push_back(i);
}

void Groups::measure(Chain& c, uint32_t m, uint64_t rid, uint32_t sample) {
    combos.find(c, m);
    for (const Hit& q : c.hits) {
        std::unordered_map<uint64_t, std::vector<uint32_t>>::const_iterator it = bykey.find(q.key);
        if (it == bykey.end()) continue;
        uint64_t sc = Combos::centrality(c, q, rid);
        for (uint32_t i : it->second) {
            const Site& e = w[i];
            if (e.s != sample) continue;
            int64_t l = q.a, r = q.b;
            if (e.lv != NONE) for (l = (int64_t)q.a - 1; l >= 0 && c.id[l] != e.lv; l--);
            if (e.rv != NONE) for (r = q.b + 1; r < m && c.id[r] != e.rv; r++);
            if (l < 0 || r >= m) continue;
            seen[c.slot].push_back({i, c.start(r) - c.start(l), sc});
        }
    }
}

static bool by_seen(const Seen& p, const Seen& q) {
    return p.site != q.site ? p.site < q.site : p.d < q.d;
}

static bool by_size(const std::pair<uint32_t, uint64_t>& p, const std::pair<uint32_t, uint64_t>& q) {
    return p.first != q.first ? p.first > q.first : p.second > q.second;
}

std::vector<Site> Groups::rows(uint32_t support) {
    std::vector<Seen> o;
    for (std::vector<Seen>& v : seen) { o.insert(o.end(), v.begin(), v.end()); v = {}; }
    std::sort(o.begin(), o.end(), by_seen);
    if (const char* dp = std::getenv("SNAPE_DUMP_D")) {
        FILE* f = std::fopen(dp, "w");
        for (const Seen& x : o) std::fprintf(f, "%u\t%lu\t%lu\t%u\n", x.site, (unsigned long)w[x.site].read, (unsigned long)x.sc, x.d);
        std::fclose(f);
    }
    std::vector<Site> add;
    for (size_t i = 0, j; i < o.size(); i = j) {
        for (j = i + 1; j < o.size() && o[j].site == o[i].site; j++);
        std::vector<std::pair<uint32_t, uint64_t>> gr;
        for (size_t k = i, l; k < j; k = l) {
            uint64_t b = o[k].sc;
            for (l = k + 1; l < j && o[l].d - o[l - 1].d <= DSAME; l++) b = std::max(b, o[l].sc);
            gr.push_back({(uint32_t)(l - k), b});
        }
        std::sort(gr.begin(), gr.end(), by_size);
        const Site& e = w[o[i].site];
        std::vector<uint64_t> got = {e.read};
        for (size_t k = 0; k < gr.size() && k < CL && gr[k].first >= support; k++) {
            uint64_t r = Combos::central_read(e.s, gr[k].second);
            if (std::find(got.begin(), got.end(), r) != got.end()) continue;
            got.push_back(r);
            Site x = {};
            x.read = r;
            x.key = e.key;
            x.s = e.s;
            x.n = e.n;
            add.push_back(x);
        }
    }
    return add;
}

struct Piece { uint64_t read; uint32_t lo, hi, k; std::string s; };

class Pieces : public Job {
public:
    std::vector<Piece> list;
    uint64_t nb;
    Pieces(const ReadFiles& files, int t) : nb(0), files(files), t(t) {}
    void build(std::vector<Site>& w, const std::vector<uint64_t>& rd, const std::vector<uint32_t>& ro);
    void cut(FILE* f);
    static std::string name(const Piece& q);

private:
    const ReadFiles& files;
    int t;
    size_t lo, hi;
    std::vector<uint32_t> gs;
    void task(int i) override;
};

static const size_t BATCH = 1 << 14;

void Pieces::build(std::vector<Site>& w, const std::vector<uint64_t>& rd, const std::vector<uint32_t>& ro) {
    for (uint32_t k = 0; k + 1 < ro.size(); k++) {
        std::vector<std::pair<uint32_t, uint32_t>> iv;
        uint32_t g = list.size();
        for (uint32_t i = ro[k]; i < ro[k + 1]; i++)
            if (w[i].hi > w[i].lo) iv.push_back({w[i].lo, w[i].hi});
        std::sort(iv.begin(), iv.end());
        for (uint32_t j = 0; j < iv.size(); j++) {
            if (j && iv[j].first <= list.back().hi) { list.back().hi = std::max(list.back().hi, iv[j].second); continue; }
            if (!j) gs.push_back(g);
            list.push_back({rd[k], iv[j].first, iv[j].second, j ? list.back().k + 1 : 0, ""});
        }
        for (uint32_t i = ro[k]; i < ro[k + 1]; i++)
            for (w[i].pc = g; w[i].pc + 1 < list.size() && list[w[i].pc + 1].lo <= w[i].lo; w[i].pc++);
    }
    gs.push_back(list.size());
}

void Pieces::task(int i) {
    Read r;
    std::vector<uint8_t> rec;
    for (size_t k = lo + i; k < hi; k += t) {
        files.get(list[gs[k]].read, rec);
        r.squeeze(rec.data());
        for (uint32_t j = gs[k]; j < gs[k + 1]; j++) list[j].s = r.bases(rec.data(), list[j].lo, list[j].hi, 0);
    }
}

void Pieces::cut(FILE* f) {
    for (lo = 0; lo + 1 < gs.size(); lo = hi) {
        hi = std::min(lo + BATCH, gs.size() - 1);
        par(t);
        for (uint32_t j = gs[lo]; j < gs[hi]; j++) {
            Fasta::put(f, name(list[j]), list[j].s);
            nb += list[j].s.size();
            std::string().swap(list[j].s);
        }
    }
}

std::string Pieces::name(const Piece& q) {
    return "s" + std::to_string(read_sample(q.read)) + "_" + std::to_string(read_index(q.read)) + "_" + std::to_string(q.k);
}

class Select : public Walker {
public:
    explicit Select(const Args& a);
    int run();

private:
    const Args& a;
    Graph g;
    Adjacency adj;
    Combos combos;
    Anchors anchors;
    ReadFiles files;
    std::vector<Site> w;
    std::vector<uint64_t> rd;
    std::vector<uint32_t> ro;
    std::vector<std::vector<uint32_t>> ids;
    Groups* groups;
    void load();
    void index_reads();
    void visit(Chain& c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) override;
    void anchor(Chain& c, uint32_t m, uint32_t rc, uint64_t rid);
    void write(Pieces& p, size_t n0) const;
};

Select::Select(const Args& a) : Walker(a.threads), a(a), g(a.in), adj(g), combos(g, adj), anchors(g), files(a.reads), ids(a.threads), groups(nullptr) {
    combos.keep(std::vector<uint8_t>(g.K, 1));
}

void Select::load() {
    Input in(a.data);
    std::vector<uint64_t> r(2 + in.h.total);
    for (uint64_t i = 0; i < in.h.n; i++) {
        std::fread(r.data(), 8, r.size(), in.f);
        Site e = {};
        e.read = r[0];
        e.key = r[1];
        e.s = read_sample(r[0]);
        for (uint64_t s = 0; s < in.h.total; s++) e.n += (s ? "," : "") + std::to_string(r[2 + s]);
        w.push_back(e);
    }
    std::fclose(in.f);
    index_reads();
}

void Select::index_reads() {
    rd.clear();
    ro.clear();
    for (size_t i = 0; i < w.size(); i++)
        if (!i || w[i].read != w[i - 1].read) { rd.push_back(w[i].read); ro.push_back(i); }
    ro.push_back(w.size());
}

void Select::visit(Chain& c, const uint8_t*, uint32_t m, uint32_t rc, uint64_t rid) {
    if (groups) groups->measure(c, m, rid, sample);
    else anchor(c, m, rc, rid);
}

void Select::anchor(Chain& c, uint32_t m, uint32_t rc, uint64_t rid) {
    uint64_t r = read_id(sample, rid), k = std::lower_bound(rd.begin(), rd.end(), r) - rd.begin();
    if (k == rd.size() || rd[k] != r) return;
    combos.find(c, m);
    std::vector<uint32_t>& sorted = ids[c.slot];
    sorted.assign(c.id.begin(), c.id.begin() + m);
    std::sort(sorted.begin(), sorted.end());
    for (uint32_t i = ro[k]; i < ro[k + 1]; i++)
        for (const Hit& q : c.hits)
            if (q.key == w[i].key) anchors.place(c, sorted, sample, q, m, rc, w[i]);
}

static bool by_row(const Site& p, const Site& q) {
    return p.read != q.read ? p.read < q.read : p.key < q.key;
}

static std::string list(const std::vector<Anchor>& v, uint32_t lo) {
    std::string s;
    for (const Anchor& x : v) s += (s.empty() ? "" : ",") + std::to_string(x.v) + ":" + std::to_string(x.a - lo);
    return s.empty() ? "-" : s;
}

void Select::write(Pieces& p, size_t n0) const {
    std::string o(a.out);
    FILE *f = std::fopen((o + ".fa").c_str(), "w"), *v = std::fopen((o + ".tsv").c_str(), "w");
    p.cut(f);
    for (uint32_t i = 0; i < w.size(); i++) {
        const Site& e = w[i];
        if (e.hi <= e.lo) continue;
        const Piece& q = p.list[e.pc];
        std::fprintf(v, "%s\t%u\t%u\t%u\t%s\t%lu\t%u\t%s\t%s\t%u\n", Pieces::name(q).c_str(), e.a - q.lo, e.b - q.lo, i, e.n.c_str(), (unsigned long)e.key, e.s, list(e.l, q.lo).c_str(), list(e.r, q.lo).c_str(), q.lo);
    }
    std::fclose(f);
    std::fclose(v);
    std::fprintf(stderr, "[select] %zu sites, %zu more rows from distance groups, on %zu reads; %zu pieces, %lu bp\n", n0, w.size() - n0, rd.size(), p.list.size(), (unsigned long)p.nb);
}

int Select::run() {
    load();
    index(g.hashes);
    fetch(files, rd);
    size_t n0 = w.size();
    std::vector<Site> add;
    {
        Groups gr(w, combos, a.threads);
        groups = &gr;
        stream(std::vector<std::string>(files.names.begin() + 1, files.names.end()), 1);
        groups = nullptr;
        add = gr.rows(a.support);
    }
    std::vector<Site> old;
    old.swap(w);
    w = add;
    std::sort(w.begin(), w.end(), by_row);
    index_reads();
    fetch(files, rd);
    w.insert(w.end(), old.begin(), old.end());
    std::stable_sort(w.begin(), w.end(), by_row);
    index_reads();
    Pieces p(files, a.threads);
    p.build(w, rd, ro);
    write(p, n0);
    return 0;
}

int cmd_select(const Args& a) {
    return Select(a).run();
}
