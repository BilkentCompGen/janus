#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "core.h"
#include "io.h"

inline uint64_t read_id(uint64_t sample, uint64_t index) { return sample << 40 | index; }
inline uint32_t read_sample(uint64_t id) { return id >> 40; }
inline uint64_t read_index(uint64_t id) { return id & 0xffffffffff; }

class ReadFiles {
public:
    std::vector<std::string> names;
    explicit ReadFiles(const char* list);
    void get(uint64_t id, std::vector<uint8_t>& rec) const;

private:
    std::vector<int> svb, idx;
};

struct Hit { uint64_t key; uint32_t a, b; };

class Chain {
public:
    Read r;
    Parser parser;
    int slot;
    std::vector<uint64_t> hash;
    std::vector<uint32_t> node, id, ix, ed;
    std::vector<Hit> hits;
    uint32_t start(uint32_t i) const { return r.a[ix[i]]; }
    uint32_t end(uint32_t i) const { return r.b[ix[i]]; }
};

class Walker : public Stream {
public:
    explicit Walker(int t);
    void index(const std::vector<uint64_t>& nodes);
    void stream(const std::vector<std::string>& files, int first);
    void fetch(const ReadFiles& files, const std::vector<uint64_t>& reads);
    void fetched(int i);

protected:
    int t, sample;
    std::vector<Chain> chains;
    virtual void visit(Chain& c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) = 0;

private:
    struct Slot { uint64_t h; uint32_t v; };
    std::vector<Slot> map;
    uint64_t mask;
    int w;
    const ReadFiles* files;
    const std::vector<uint64_t>* reads;
    size_t lo, hi;
    void block(Block& b) override;
    void walk(Chain& c, const uint8_t* rec, uint64_t rid);
};
