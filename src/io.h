#pragma once
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

static const int MAXT = 64;

inline void lock(uint8_t& l) {
    while (__atomic_exchange_n(&l, 1, __ATOMIC_ACQUIRE));
}

inline void unlock(uint8_t& l) {
    __atomic_store_n(&l, 0, __ATOMIC_RELEASE);
}

class Job {
public:
    virtual ~Job() {}
    void par(int n);
    static uint64_t share(uint64_t total, int i, int n) { return total * i / n; }

protected:
    virtual void task(int i) = 0;

private:
    static void* enter(void* seat);
};

struct Header {
    char magic[4];
    uint32_t ver;
    uint64_t n, total;
    uint8_t resv[8];
    Header() {}
    Header(char kind, uint32_t ver, uint64_t n, uint64_t total);
};

class Input {
public:
    FILE* f;
    Header h;
    explicit Input(const std::string& path);
};

class Output {
public:
    FILE* f;
    explicit Output(const char* path);
    void finish(const Header& h);
    static void save(const char* path, const Header& h, const void* p, size_t n);

private:
    std::string path;
};

inline uint64_t be64(const uint8_t* p) {
    uint64_t x;
    std::memcpy(&x, p, 8);
    return __builtin_bswap64(x);
}

class BitReader {
public:
    BitReader(const uint8_t* p, uint64_t pos) : p(p), pos(pos) {}
    uint32_t get(int w) {
        uint32_t v = (uint32_t)((be64(p + (pos >> 3)) << (pos & 7)) >> (64 - w));
        pos += w;
        return v;
    }

private:
    const uint8_t* p;
    uint64_t pos;
};

class BitWriter {
public:
    explicit BitWriter(std::vector<char>* out) : out(out), acc(0), pos(0) {}
    void put(uint64_t v, int w) {
        acc = acc << w | v;
        for (pos += w; pos >= 8; pos -= 8) out->push_back((char)(acc >> (pos - 8)));
    }
    void flush() {
        if (pos) put(0, 8 - pos);
    }

private:
    std::vector<char>* out;
    uint64_t acc;
    int pos;
};

static const uint32_t NONE = ~0u;
static const uint32_t PF = 16;

inline size_t svb_bytes(uint32_t tag) {
    uint32_t len = tag & 0x7fffffff;
    return 4 + (tag >> 31 ? len : (len + 3) / 4);
}

inline size_t svb_size(const char* p, size_t n) {
    if (n < 4) return 0;
    uint32_t tag;
    std::memcpy(&tag, p, 4);
    size_t sz = svb_bytes(tag);
    return sz <= n ? sz : 0;
}

inline size_t exc_size(const char* p, size_t n) {
    if (n < 6) return 0;
    size_t sz = (uint8_t)p[4] << 8 | (uint8_t)p[5];
    return sz <= n ? sz : 0;
}

struct Block {
    int id;
    const char* in[2];
    size_t inlen[2];
    std::vector<char> out;
    uint64_t rec;
};

class Stream : public Job {
public:
    uint64_t pass(FILE* fi, FILE* fe, FILE* fo, int t);

protected:
    virtual void block(Block& b) = 0;
    void task(int i) override;

private:
    struct Side { FILE* f; char* b[2]; size_t len, off, start, next; };
    Side s[2];
    Block b[MAXT];
    int n, nb, cur;
    size_t cap;
    uint64_t rec, nrec;
    bool whole();
    void cut();
    bool grow();
};

std::vector<std::string> split(const char* s);

class Lines {
public:
    char* s;
    ssize_t n;
    explicit Lines(const char* path);
    ~Lines();
    bool next();

private:
    FILE* f;
    size_t cap;
};
