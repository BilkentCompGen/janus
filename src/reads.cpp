#include <fcntl.h>
#include <unistd.h>
#include "io.h"
#include "reads.h"

ReadFiles::ReadFiles(const char* list) : names(split(list)) {
    for (std::string& f : names) {
        svb.push_back(open(f.c_str(), O_RDONLY));
        idx.push_back(open((f + ".idx").c_str(), O_RDONLY));
    }
}

void ReadFiles::get(uint64_t id, std::vector<uint8_t>& rec) const {
    uint64_t off[2];
    uint32_t s = read_sample(id);
    pread(idx[s], off, 16, sizeof(Header) + 8 * read_index(id));
    rec.assign(off[1] - off[0] + 8, 0);
    pread(svb[s], rec.data(), off[1] - off[0], off[0]);
}

class FetchJob : public Job {
public:
    explicit FetchJob(Walker& w) : w(w) {}

private:
    Walker& w;
    void task(int i) override { w.fetched(i); }
};

Walker::Walker(int t) : t(t), sample(0), chains(t), mask(0), w(0), files(nullptr), reads(nullptr), lo(0), hi(0) {
    for (int i = 0; i < t; i++) chains[i].slot = i;
}

void Walker::index(const std::vector<uint64_t>& nodes) {
    uint64_t z = 1, n = nodes.size();
    while (z < n + n / 2) z <<= 1;
    map.assign(z, {0, NONE});
    mask = z - 1;
    for (uint32_t v = 0; v < n; v++) {
        if (v + PF < n) __builtin_prefetch(&map[nodes[v + PF] & mask]);
        uint64_t i = nodes[v] & mask;
        while (map[i].h) i = (i + 1) & mask;
        map[i] = {nodes[v], v};
    }
}

void Walker::walk(Chain& c, const uint8_t* rec, uint64_t rid) {
    Read& r = c.r;
    if (c.hash.size() < r.n) c.hash.resize(2 * r.n);
    if (c.node.size() < r.n) c.node.resize(2 * r.n);
    for (uint32_t i = 0; i < r.n; i++) {
        c.hash[i] = r.hash(i);
        __builtin_prefetch(&map[c.hash[i] & mask]);
    }
    for (uint32_t i = 0; i < r.n; i++) {
        uint64_t j = c.hash[i] & mask;
        while (map[j].h && map[j].h != c.hash[i]) j = (j + 1) & mask;
        c.node[i] = map[j].v;
    }
    if (c.id.size() < r.n) { c.id.resize(2 * r.n); c.ix.resize(2 * r.n); }
    for (uint32_t s = 0; s < 2; s++) {
        uint32_t m = 0;
        for (uint32_t i = s ? r.nf : 0; i < (s ? r.n : r.nf); i++)
            if (c.node[i] != NONE) { c.id[m] = c.node[i]; c.ix[m++] = i; }
        visit(c, rec, m, s ? 3 : 0, rid);
    }
}

void Walker::block(Block& b) {
    Chain& c = chains[b.id];
    uint64_t rid = b.rec;
    for (size_t o = 0, e = 0; o < b.inlen[0]; rid++) {
        const uint8_t* rec = c.r.next(b, &o, &e, w);
        walk(c, rec, rid);
    }
}

void Walker::stream(const std::vector<std::string>& names, int first) {
    sample = first;
    for (const std::string& f : names) {
        Input svb(f), exc(f + ".exc");
        w = exc.h.resv[1];
        pass(svb.f, exc.f, nullptr, t);
        sample++;
    }
}

void Walker::fetched(int i) {
    Chain& c = chains[i];
    std::vector<uint8_t> rec;
    for (size_t k = lo + i; k < hi; k += t) {
        files->get((*reads)[k], rec);
        c.parser.cores(c.r, rec.data());
        walk(c, rec.data(), read_index((*reads)[k]));
    }
}

void Walker::fetch(const ReadFiles& f, const std::vector<uint64_t>& list) {
    FetchJob job(*this);
    files = &f;
    reads = &list;
    sample = 0;
    for (size_t k = 0; k < list.size(); sample++) {
        lo = k;
        while (k < list.size() && read_sample(list[k]) == (uint64_t)sample) k++;
        hi = k;
        if (hi > lo) job.par(t);
    }
}
