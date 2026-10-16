#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fcntl.h>
#include <filereader/Standard.hpp>
#include <map>
#include <memory>
#include <mutex>
#include <rapidgzip/ParallelGzipReader.hpp>
#include <string>
#include <unistd.h>
#include <vector>
#include "cmd.h"
#include "core.h"
#include "io.h"

struct Gzip;

class Source {
public:
    Source(const char* path, int t);
    ~Source();
    size_t read(char* o, size_t n);

private:
    int fd;
    Gzip* gz;
};

static const uint64_t GZCHUNK = 16u << 20;

struct Gzip {
    rapidgzip::ParallelGzipReader<rapidgzip::ChunkData> reader;
    Gzip(const char* path, int t) : reader(std::make_unique<rapidgzip::StandardFileReader>(path), t, GZCHUNK) {
        reader.setCRC32Enabled(true);
    }
};

Source::Source(const char* path, int t) : fd(open(path, O_RDONLY)), gz(nullptr) {
    unsigned char m[2] = {0, 0};
    if (fd < 0 || pread(fd, m, 2, 0) < 0) { std::fprintf(stderr, "[pack] cannot open %s\n", path); std::exit(1); }
    if (m[0] == 0x1f && m[1] == 0x8b) {
        close(fd);
        fd = -1;
        gz = new Gzip(path, t);
    } else {
        int d = open(path, O_RDONLY | O_DIRECT);
        if (d >= 0) { close(fd); fd = d; }
        else posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
    }
}

Source::~Source() {
    if (fd >= 0) close(fd);
    delete gz;
}

static void corrupt(const char* what) {
    std::fprintf(stderr, "[pack] corrupt input: %s\n", what);
    std::exit(1);
}

size_t Source::read(char* o, size_t n) {
    size_t got = 0;
    ssize_t r = 0;
    try {
        while (got < n && (r = gz ? (ssize_t)gz->reader.read(o + got, n - got) : ::read(fd, o + got, n - got)) > 0) got += r;
    } catch (const std::exception& e) {
        corrupt(e.what());
    }
    if (!gz && r < 0) corrupt("read failed");
    return got;
}

class Chunk {
public:
    uint64_t id;
    size_t len;
    std::vector<uint8_t> out;
    std::vector<uint32_t> size;
    uint64_t bases, longest;
    bool bam;
    Chunk();
    ~Chunk();
    bool fill(Source& src, std::vector<char>& carry);
    void encode();

private:
    char *buf, *txt;
    size_t head;
    void put(const char* s, uint32_t n);
    void fastq();
    void records();
};

static const size_t CHUNK = 64u << 20;
static const size_t HEAD = 4u << 20, AL = 4096;

static char* aligned(size_t head) {
    return (char*)std::aligned_alloc(AL, head + CHUNK + AL);
}

Chunk::Chunk() : id(0), len(0), bases(0), longest(0), bam(false), buf(aligned(HEAD)), txt(nullptr), head(HEAD) {}

Chunk::~Chunk() {
    std::free(buf);
}

static size_t width(const char* p, const size_t* st, size_t i) {
    size_t a = st[i], b = st[i - 1] - 1;
    return (b > a && p[b - 1] == '\r') ? b - a - 1 : b - a;
}

static size_t last_record_end(const char* p, size_t n) {
    size_t st[12], m = 0, e = n;
    while (m < 12) {
        const char* q = e ? (const char*)memrchr(p, '\n', e) : nullptr;
        st[m++] = q ? q - p + 1 : 0;
        if (!q) break;
        e = q - p;
    }
    for (size_t i = 4; i < m; i++)
        if (p[st[i]] == '@' && p[st[i - 2]] == '+' && width(p, st, i - 1) == width(p, st, i - 3)) return st[i - 4];
    return 0;
}

static size_t last_alignment_end(const char* p, size_t n) {
    size_t at = 0;
    uint32_t bs;
    while (at + 4 <= n) {
        std::memcpy(&bs, p + at, 4);
        if (at + 4 + bs > n) break;
        at += 4 + bs;
    }
    return at;
}

static void skip_bam_header(Source& src) {
    char b[8];
    uint32_t n, l;
    if (src.read(b, 8) != 8 || std::memcmp(b, "BAM\1", 4)) corrupt("not a BAM file");
    std::memcpy(&n, b + 4, 4);
    std::vector<char> t(n);
    src.read(t.data(), n);
    src.read((char*)&n, 4);
    for (uint32_t i = 0; i < n; i++) {
        src.read((char*)&l, 4);
        t.resize(l + 4);
        src.read(t.data(), l + 4);
    }
}

bool Chunk::fill(Source& src, std::vector<char>& carry) {
    if (carry.size() > head) {
        std::free(buf);
        head = (carry.size() + AL) / AL * AL * 2;
        buf = aligned(head);
    }
    txt = buf + head - carry.size();
    std::memcpy(txt, carry.data(), carry.size());
    size_t r = src.read(buf + head, CHUNK), n = carry.size() + r, cut;
    bool eof = r < CHUNK;
    if (bam) cut = last_alignment_end(txt, n);
    else if (eof) {
        if (n && txt[n - 1] != '\n') txt[n++] = '\n';
        cut = n;
    } else cut = last_record_end(txt, n);
    carry.assign(txt + cut, txt + n);
    len = cut;
    return eof;
}

static bool line(const char*& p, const char* e, const char*& s, size_t& l) {
    if (p >= e) return false;
    const char* q = (const char*)std::memchr(p, '\n', e - p);
    if (!q) q = e;
    s = p;
    l = q - p;
    p = q + 1;
    while (l && s[l - 1] == '\r') l--;
    return true;
}

void Chunk::put(const char* s, uint32_t n) {
    uint32_t raw = !BASES.plain((const uint8_t*)s, n), tag = n | raw << 31;
    size_t at = out.size(), bytes = svb_bytes(tag);
    out.resize(at + bytes);
    std::memcpy(out.data() + at, &tag, 4);
    if (raw) std::memcpy(out.data() + at + 4, s, n);
    else BASES.pack((const uint8_t*)s, n, out.data() + at + 4);
    size.push_back(bytes);
    bases += n;
    longest = std::max<uint64_t>(longest, n);
}

void Chunk::fastq() {
    const char *p = txt, *e = p + len, *h, *s, *x, *y;
    size_t hl, sl, xl, yl;
    while (line(p, e, h, hl) && line(p, e, s, sl) && line(p, e, x, xl) && line(p, e, y, yl)) put(s, sl);
}

void Chunk::records() {
    static const char CODE[] = "=ACMGRSVTWYHKDBN", COMP[] = "=TGKCYSBAWRDMHVN";
    const uint8_t *p = (const uint8_t*)txt, *e = p + len;
    std::string s;
    while (p < e) {
        uint32_t bs, n;
        uint16_t cigars, flag;
        std::memcpy(&bs, p, 4);
        std::memcpy(&cigars, p + 16, 2);
        std::memcpy(&flag, p + 18, 2);
        std::memcpy(&n, p + 20, 4);
        const uint8_t* q = p + 36 + p[12] + 4 * cigars;
        p += 4 + bs;
        if (flag & 0x900) continue;
        s.resize(n);
        for (uint32_t i = 0; i < n; i++) {
            uint8_t c = (q[i >> 1] >> (i & 1 ? 0 : 4)) & 15;
            if (flag & 16) s[n - 1 - i] = COMP[c];
            else s[i] = CODE[c];
        }
        put(s.data(), n);
    }
}

void Chunk::encode() {
    out.clear();
    size.clear();
    bases = longest = 0;
    out.reserve(len / 7);
    if (bam) records();
    else fastq();
}

class Pack : public Job {
public:
    explicit Pack(const Args& a);
    int run();

private:
    const Args& a;
    int t, encoders;
    Source src;
    std::vector<Chunk> chunks;
    std::mutex m;
    std::condition_variable cv;
    std::deque<Chunk*> todo, pool;
    std::map<uint64_t, Chunk*> done;
    bool fed;
    uint64_t made, reads, total, longest, pos;
    FILE *svb, *idx;
    bool bam;
    void task(int i) override;
    void read();
    void encode();
    void write();
    void shuffle();
};

Pack::Pack(const Args& a) : a(a), t(std::max(1, a.threads)), encoders(std::max(2, t / 5)), src(a.in, t), chunks(encoders + 4), fed(false), made(0), reads(0), total(0), longest(0), pos(sizeof(Header)) {
    std::string in(a.in);
    bam = in.size() > 4 && in.compare(in.size() - 4, 4, ".bam") == 0;
    for (Chunk& c : chunks) { c.bam = bam; pool.push_back(&c); }
}

void Pack::task(int i) {
    if (i == 0) read();
    else if (i == 1) write();
    else encode();
}

void Pack::read() {
    std::vector<char> carry;
    if (bam) skip_bam_header(src);
    for (bool eof = false; !eof;) {
        Chunk* b;
        {
            std::unique_lock<std::mutex> l(m);
            while (pool.empty()) cv.wait(l);
            b = pool.front();
            pool.pop_front();
        }
        eof = b->fill(src, carry);
        std::lock_guard<std::mutex> l(m);
        if (b->len) { b->id = made++; todo.push_back(b); }
        else pool.push_back(b);
        cv.notify_all();
    }
    std::lock_guard<std::mutex> l(m);
    fed = true;
    cv.notify_all();
}

void Pack::encode() {
    for (;;) {
        Chunk* b;
        {
            std::unique_lock<std::mutex> l(m);
            while (todo.empty() && !fed) cv.wait(l);
            if (todo.empty()) return;
            b = todo.front();
            todo.pop_front();
        }
        b->encode();
        std::lock_guard<std::mutex> l(m);
        done[b->id] = b;
        cv.notify_all();
    }
}

void Pack::write() {
    std::vector<uint64_t> off;
    for (uint64_t id = 0;; id++) {
        Chunk* b;
        {
            std::unique_lock<std::mutex> l(m);
            while (!done.count(id) && !(fed && todo.empty() && id >= made)) cv.wait(l);
            if (!done.count(id)) break;
            b = done[id];
            done.erase(id);
        }
        off.clear();
        for (size_t i = 0; i < b->size.size(); i++) { off.push_back(pos); pos += b->size[i]; }
        std::fwrite(off.data(), 8, off.size(), idx);
        std::fwrite(b->out.data(), 1, b->out.size(), svb);
        reads += b->size.size();
        total += b->bases;
        longest = std::max(longest, b->longest);
        std::lock_guard<std::mutex> l(m);
        pool.push_back(b);
        cv.notify_all();
    }
}

void Pack::shuffle() {
    std::string x = std::string(a.out) + ".idx";
    std::vector<uint64_t> off(reads + 1), order(reads);
    {
        Input in(x);
        std::fread(off.data(), 8, off.size(), in.f);
        std::fclose(in.f);
    }
    uint64_t seed = 0x9e3779b97f4a7c15ull;
    for (uint64_t i = 0; i < reads; i++) order[i] = i;
    for (uint64_t i = reads; i > 1; i--) {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        std::swap(order[i - 1], order[seed % i]);
    }
    int fd = open(a.out, O_RDONLY);
    Output so(a.out), xo(x.c_str());
    std::vector<char> rec;
    uint64_t at = sizeof(Header), n = 0, kept = 0, top = 0;
    for (uint64_t i = 0; i < reads; i++) {
        uint64_t o = off[order[i]], len = off[order[i] + 1] - o;
        uint32_t tag;
        rec.resize(len);
        if (pread(fd, rec.data(), len, o) != (ssize_t)len) corrupt("read failed");
        std::memcpy(&tag, rec.data(), 4);
        std::fwrite(&at, 8, 1, xo.f);
        std::fwrite(rec.data(), 1, len, so.f);
        at += len;
        n++;
        kept += tag & 0x7fffffffu;
        top = std::max<uint64_t>(top, tag & 0x7fffffffu);
    }
    close(fd);
    std::fwrite(&at, 8, 1, xo.f);
    reads = n;
    total = kept;
    longest = top;
    so.finish(Header('1', 1, reads, total));
    xo.finish(Header('I', 1, reads, longest));
}

int Pack::run() {
    std::string x = std::string(a.out) + ".idx";
    {
        Output so(a.out), xo(x.c_str());
        svb = so.f;
        idx = xo.f;
        par(encoders + 2);
        std::fwrite(&pos, 8, 1, idx);
        so.finish(Header('1', 1, reads, total));
        xo.finish(Header('I', 1, reads, longest));
    }
    if (bam && a.shuffle) shuffle();
    std::fprintf(stderr, "[pack] %lu reads, %lu bases\n", reads, total);
    return 0;
}

int cmd_pack(const Args& a) {
    return Pack(a).run();
}
