#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include "cmd.h"
#include "core.h"
#include "io.h"

class CoreTable : public Job {
public:
    uint64_t z, K;
    uint8_t full;
    CoreTable(int t, uint32_t support);
    ~CoreTable();
    static uint64_t slots(uint64_t cores);
    bool make(uint64_t z);
    void release();
    uint64_t home(uint64_t h) const { return (h >> 48) * z + (((h >> 32) & 0xffff) * z >> 16); }
    void prefetch(uint64_t h) const { __builtin_prefetch(&tab[home(h)], 1); }
    uint64_t add(uint64_t h) {
        uint64_t s = (h >> 48) * z, tag = h << 16 >> 16, at = home(h);
        for (uint64_t i = at;;) {
            uint64_t v = __atomic_load_n(&tab[i], __ATOMIC_RELAXED);
            while (!v || v >> 16 == tag) {
                if ((v & 0xffff) == 0xffff) return i;
                if (__atomic_compare_exchange_n(&tab[i], &v, v ? v + 1 : tag << 16 | 1, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return i;
            }
            i = i + 1 == s + z ? s : i + 1;
            if (i == at) { __atomic_store_n(&full, 1, __ATOMIC_RELAXED); return at; }
        }
    }
    bool filled() const { return __atomic_load_n(&full, __ATOMIC_RELAXED); }
    void rank(int fd);
    uint64_t first(int i) const { return shard(i) * z; }
    uint64_t id(uint64_t slot) const { return tab[slot]; }

private:
    uint64_t* tab;
    std::vector<uint32_t> sc;
    int t, fd;
    uint32_t support;
    bool counted;
    uint64_t shard(int i) const;
    void task(int i) override;
    void count(int i);
    void emit(int i);
};

static const uint32_t NS = 65536;
static const uint64_t ZDIV = 24;
static const uint64_t ZADD = 120000000;
static const uint64_t ZPAD = 64;

CoreTable::CoreTable(int t, uint32_t support) : z(0), K(0), full(0), tab(nullptr), sc(NS + 1), t(t), fd(-1), support(support), counted(false) {}

CoreTable::~CoreTable() {
    release();
}

uint64_t CoreTable::slots(uint64_t cores) {
    return (cores / ZDIV + ZADD) / NS + ZPAD;
}

bool CoreTable::make(uint64_t size) {
    release();
    if (NS * size > (1ull << 31)) return false;
    z = size;
    full = 0;
    tab = (uint64_t*)std::calloc(NS * z, 8);
    return true;
}

void CoreTable::release() {
    std::free(tab);
    tab = nullptr;
}

uint64_t CoreTable::shard(int i) const {
    return share(NS, i, t);
}

void CoreTable::count(int i) {
    for (uint64_t s = shard(i); s < shard(i + 1); s++)
        for (uint64_t k = s * z; k < (s + 1) * z; k++) sc[s + 1] += (tab[k] & 0xffff) >= support;
}

void CoreTable::emit(int i) {
    std::vector<uint64_t> v(z);
    std::vector<uint16_t> n(z);
    for (uint64_t s = shard(i); s < shard(i + 1); s++) {
        uint64_t m = 0, *slot = tab + s * z;
        for (uint64_t k = 0; k < z; k++) {
            if ((slot[k] & 0xffff) >= support) v[m++] = slot[k] >> 16 << 16 | k;
            else slot[k] = 0;
        }
        std::sort(v.begin(), v.begin() + m);
        for (uint64_t j = 0; j < m; j++) {
            uint64_t k = v[j] & 0xffff;
            n[j] = (uint16_t)slot[k];
            v[j] = s << 48 | v[j] >> 16;
            slot[k] = sc[s] + j + 1;
        }
        pwrite(fd, v.data(), m * 8, sizeof(Header) + (uint64_t)sc[s] * 8);
        pwrite(fd, n.data(), m * 2, sizeof(Header) + K * 8 + (uint64_t)sc[s] * 2);
    }
}

void CoreTable::task(int i) {
    if (counted) emit(i);
    else count(i);
}

void CoreTable::rank(int out) {
    fd = out;
    par(t);
    for (uint32_t s = 0; s < NS; s++) sc[s + 1] += sc[s];
    K = sc[NS];
    counted = true;
    par(t);
}

struct Edge { uint64_t ends; uint16_t reads; };

class PairTable : public Job {
public:
    std::vector<std::vector<Edge>> edges;
    std::vector<uint32_t> degree;
    uint8_t full;
    PairTable(CoreTable& cores, int t);
    ~PairTable();
    void make();
    void release();
    static uint64_t key(uint64_t a, uint64_t b) { return 1ull << 63 | a << 31 | b; }
    uint64_t home(uint64_t k) const { return (core_hash(&k, 32) >> 32) * np >> 32; }
    void prefetch(uint64_t k) const { __builtin_prefetch(&pr[home(k)], 1); }
    void add(uint64_t k) {
        uint64_t i = home(k);
        for (uint32_t p = 0; p < 1u << 16; p++, i = i + 1 == np ? 0 : i + 1) {
            uint64_t v = 0;
            if (__atomic_compare_exchange_n(&pr[i], &v, k, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED) || v == k) {
                uint16_t c = __atomic_load_n(&pc[i], __ATOMIC_RELAXED);
                while (c < 0xffff && !__atomic_compare_exchange_n(&pc[i], &c, c + 1, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
                return;
            }
        }
        __atomic_store_n(&full, 1, __ATOMIC_RELAXED);
    }
    bool filled() const { return __atomic_load_n(&full, __ATOMIC_RELAXED); }
    void pick();

private:
    CoreTable& cores;
    uint64_t np, *pr;
    uint16_t* pc;
    int t;
    void task(int i) override;
};

static const uint64_t M31 = (1u << 31) - 1;

PairTable::PairTable(CoreTable& cores, int t) : edges(t), full(0), cores(cores), np(0), pr(nullptr), pc(nullptr), t(t) {}

PairTable::~PairTable() {
    release();
}

void PairTable::make() {
    release();
    np = 65536 * cores.z * 5 / 4;
    full = 0;
    pr = (uint64_t*)std::calloc(np, 8);
    pc = (uint16_t*)std::calloc(np, 2);
}

void PairTable::release() {
    std::free(pr);
    std::free(pc);
    pr = nullptr;
    pc = nullptr;
}

static bool by_ends(const Edge& a, const Edge& b) {
    return a.ends < b.ends;
}

void PairTable::task(int i) {
    std::vector<Edge>& out = edges[i];
    uint64_t lo = cores.first(i), hi = cores.first(i + 1);
    for (uint64_t j = 0; j < np; j++) {
        uint64_t k = pr[j], s1 = k >> 31 & M31;
        if (!k || s1 < lo || s1 >= hi) continue;
        uint64_t r1 = cores.id(s1), r2 = cores.id(k & M31);
        if (r1 && r2) {
            out.push_back({(r1 - 1) << 32 | (r2 - 1), pc[j]});
            degree[r1 - 1]++;
        }
    }
    std::sort(out.begin(), out.end(), by_ends);
}

void PairTable::pick() {
    degree.resize(cores.K);
    par(t);
}

class Count : public Stream {
public:
    explicit Count(const Args& a);
    int run();

private:
    const Args& a;
    CoreTable cores;
    PairTable pairs;
    int w;
    std::vector<Read> read;
    std::vector<std::vector<uint64_t>> slots;
    void block(Block& b) override;
    void chain(Read& r, uint32_t lo, uint32_t hi, uint64_t* h);
    void write(Output& out);
};

static const uint32_t MARGIN = 10;

Count::Count(const Args& a) : a(a), cores(a.threads, a.support), pairs(cores, a.threads), w(0), read(a.threads), slots(a.threads) {}

void Count::chain(Read& r, uint32_t lo, uint32_t hi, uint64_t* h) {
    uint64_t* q = h + r.n;
    uint32_t m = 0, e = 0;
    for (uint32_t i = lo; i < hi; i++)
        if (r.a[i] >= MARGIN && r.b[i] + MARGIN <= r.len) h[m++] = r.hash(i);
    for (uint32_t i = 0; i < m; i++) {
        if (i + PF < m) cores.prefetch(h[i + PF]);
        h[i] = cores.add(h[i]);
    }
    for (uint32_t i = 1; i < m; i++)
        if (h[i - 1] != h[i]) q[e++] = PairTable::key(h[i - 1], h[i]);
    for (uint32_t i = 0; i < e; i++) {
        if (i + PF < e) pairs.prefetch(q[i + PF]);
        pairs.add(q[i]);
    }
}

void Count::block(Block& b) {
    Read& r = read[b.id];
    std::vector<uint64_t>& h = slots[b.id];
    for (size_t o = 0, e = 0; o < b.inlen[0] && !cores.filled() && !pairs.filled();) {
        r.next(b, &o, &e, w);
        if (h.size() < 2 * r.n) h.resize(4 * r.n);
        chain(r, 0, r.nf, h.data());
        chain(r, r.nf, r.n, h.data());
    }
}

static int nbits(uint64_t v) {
    return v ? 64 - __builtin_clzll(v) : 1;
}

void Count::write(Output& out) {
    uint64_t E = 0, most = 0;
    for (uint32_t d : pairs.degree) { E += d; most = std::max<uint64_t>(most, d); }
    std::vector<char> stream;
    std::vector<uint16_t> reads;
    BitWriter bw(&stream);
    int rb = nbits(cores.K - 1), db = nbits(most);
    for (uint32_t d : pairs.degree) bw.put(d, db);
    for (std::vector<Edge>& part : pairs.edges)
        for (Edge& e : part) { bw.put((uint32_t)e.ends, rb); reads.push_back(e.reads); }
    bw.flush();
    std::fseek(out.f, sizeof(Header) + cores.K * 10, SEEK_SET);
    std::fwrite(reads.data(), 2, E, out.f);
    std::fwrite(stream.data(), 1, stream.size(), out.f);
    Header h('C', 1, cores.K, E);
    h.resv[0] = rb;
    h.resv[1] = db;
    h.resv[2] = a.support;
    h.resv[3] = MARGIN;
    out.finish(h);
}

int Count::run() {
    std::string ex = std::string(a.in) + ".exc";
    Input svb(a.in), exc(ex);
    Output out(a.out);
    w = exc.h.resv[1];
    cores.make(CoreTable::slots(exc.h.total));
    for (;;) {
        pairs.make();
        pass(svb.f, exc.f, nullptr, a.threads);
        if (!cores.full && !pairs.full) break;
        uint64_t z = cores.z * 3 / 2;
        std::fprintf(stderr, "[count] the tables filled up at %lu slots per shard; counting again with %lu\n", (unsigned long)cores.z, (unsigned long)z);
        pairs.release();
        if (!cores.make(z)) { std::fprintf(stderr, "[count] too many distinct cores\n"); return 1; }
        svb = Input(a.in);
        exc = Input(ex);
    }
    cores.rank(fileno(out.f));
    pairs.pick();
    cores.release();
    pairs.release();
    write(out);
    return 0;
}

int cmd_count(const Args& a) {
    return Count(a).run();
}
