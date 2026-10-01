#include <algorithm>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <filereader/Standard.hpp>
#include <rapidgzip/ParallelGzipReader.hpp>
#include "cmd.h"
#include "io.h"

// pack: a FASTQ (plain or gzip, any number of gzip members) into an svb and its index. Reads are stored in input order:
// a 4-byte tag (length, top bit set for a raw read) and the bases at 2 bits each, or the raw bytes when the read has a
// character other than ACGT (either case). The index holds each read's offset in the svb and the end offset.
// One reader fills CHUNK-sized blocks straight from the source (gzip is inflated in parallel by rapidgzip with ISA-L,
// on t threads, CRC checked) and cuts each block after its last whole record, found by looking back from the block's
// end; max(2, t/5) workers encode the blocks; the caller writes them in input order, so the output does not depend on t.

static const size_t CHUNK = 64u << 20;   // text per block
static const uint64_t GZCHUNK = 16u << 20;   // rapidgzip's compressed chunk per decoding task: 16 MiB kept the peak at ~6 GB

struct pblock { uint64_t id; char* buf; size_t cap, len; std::vector<uint8_t> out; std::vector<uint32_t> size; uint64_t bases, mx; };

// the text source: rapidgzip for gzip (one or more members), plain reads otherwise
struct source {
    FILE* f = nullptr; std::unique_ptr<rapidgzip::ParallelGzipReader<rapidgzip::ChunkData>> gz;
    source(const char* path, int t) {
        f = std::fopen(path, "rb");
        if (!f) { std::fprintf(stderr, "[pack] cannot open %s\n", path); std::exit(1); }
        int a = std::fgetc(f), b = std::fgetc(f); std::rewind(f);
        if (a == 0x1f && b == 0x8b) {
            std::fclose(f); f = nullptr;
            gz = std::make_unique<rapidgzip::ParallelGzipReader<rapidgzip::ChunkData>>(
                std::make_unique<rapidgzip::StandardFileReader>(path), t, GZCHUNK);
            gz->setCRC32Enabled(true);
        } else std::setvbuf(f, nullptr, _IONBF, 0);
    }
    ~source() { if (f) std::fclose(f); }
    size_t read(char* o, size_t n) {   // fills o with up to n bytes; less only at the end
        size_t got = 0, r;
        try {
            while (got < n && (r = gz ? gz->read(o + got, n - got) : std::fread(o + got, 1, n - got, f)) > 0) got += r;
        } catch (const std::exception& e) { std::fprintf(stderr, "[pack] corrupt input: %s\n", e.what()); std::exit(1); }
        return got;
    }
};

// the end of the last whole 4-line record in [p, p + n), or 0: a record is "@..." / sequence / "+..." / quality of the
// sequence's length; a sequence line never starts with '@' or '+', so the pattern cannot start inside a record
static size_t last_record_end(const char* p, size_t n) {
    std::vector<size_t> st;   // line starts, last first; st[0] is the start of the unterminated tail (or n)
    size_t e = n;
    while (st.size() < 12) {
        const char* q = e ? (const char*)memrchr(p, '\n', e) : nullptr;
        st.push_back(q ? q - p + 1 : 0);
        if (!q) break;
        e = q - p;
    }
    auto len = [&](size_t i) { size_t a = st[i], b = st[i - 1] - 1; return (b > a && p[b - 1] == '\r') ? b - a - 1 : b - a; };
    for (size_t i = 4; i < st.size(); i++)   // line st[i] .. st[i-3] are the four lines, st[i-4] the cut
        if (p[st[i]] == '@' && p[st[i - 2]] == '+' && len(i - 1) == len(i - 3)) return st[i - 4];
    return 0;
}

static unsigned char ENC[256], BAD[256];

static void encode(pblock* b) {
    const char *p = b->buf, *e = p + b->len;
    auto next = [&](const char*& s, size_t& l) {   // one line without its trailing \r / \n
        if (p >= e) return false;
        const char* q = (const char*)std::memchr(p, '\n', e - p); if (!q) q = e;
        s = p; l = q - p; p = q + 1;
        while (l && s[l - 1] == '\r') l--;
        return true;
    };
    const char *h, *s, *x, *y; size_t hl, sl, xl, yl;
    b->out.clear(); b->size.clear(); b->bases = b->mx = 0;
    b->out.reserve(b->len / 7);
    while (next(h, hl) && next(s, sl) && next(x, xl) && next(y, yl)) {
        uint32_t L = sl, raw = 0;
        const unsigned char* u = (const unsigned char*)s;
        for (size_t i = 0; i < sl; i++) raw |= BAD[u[i]];
        uint32_t tag = L | raw << 31; size_t n = svb_bytes(tag), at = b->out.size();
        b->out.resize(at + n);
        std::memcpy(&b->out[at], &tag, 4);
        uint8_t* d = &b->out[at + 4];
        if (raw) std::memcpy(d, s, L);
        else {
            size_t i = 0, k = 0;
            for (; i + 4 <= L; i += 4) d[k++] = ENC[u[i]] << 6 | ENC[u[i + 1]] << 4 | ENC[u[i + 2]] << 2 | ENC[u[i + 3]];
            if (i < L) { uint8_t c = 0; for (int j = 0; i < L; i++, j++) c |= ENC[u[i]] << (6 - 2 * j); d[k] = c; }
        }
        b->size.push_back(n); b->bases += L; b->mx = std::max<uint64_t>(b->mx, L);
    }
}

int cmd_pack(const char* in, const char* out, int t) {
    t = std::max(1, t);
    for (int c = 0; c < 256; c++) { ENC[c] = code(c); BAD[c] = code(c) == 3 && (c & 0xdf) != 'T'; }
    // decoding takes ~90 % of the CPU (one encoder does ~4 GB of text per second), so gzip gets all t threads and
    // about a fifth as many encoders run beside it
    int et = std::max(2, t / 5);
    source src(in, t);
    std::string x = std::string(out) + ".idx";
    FILE *fo = open_tmp(out), *fx = open_tmp(x.c_str());
    std::mutex m; std::condition_variable cv;
    std::deque<pblock*> todo, pool; std::map<uint64_t, pblock*> done; bool fed = false; uint64_t made = 0;
    const size_t MAXQ = (size_t)et + 4;   // blocks in flight: bounds memory at about MAXQ x CHUNK x 1.2
    for (size_t i = 0; i < MAXQ; i++) pool.push_back(new pblock{0, new char[CHUNK + 1], CHUNK, 0, {}, {}, 0, 0});

    std::thread reader([&] {
        std::vector<char> carry;
        for (bool eof = false; !eof;) {
            pblock* b;
            { std::unique_lock<std::mutex> l(m); cv.wait(l, [&] { return !pool.empty(); }); b = pool.front(); pool.pop_front(); }
            if (carry.size() * 2 > b->cap) { delete[] b->buf; b->cap = carry.size() * 4; b->buf = new char[b->cap + 1]; }
            std::memcpy(b->buf, carry.data(), carry.size());
            size_t n = carry.size() + src.read(b->buf + carry.size(), b->cap - carry.size());
            eof = n < b->cap;
            size_t cut;
            if (eof) {
                if (n && b->buf[n - 1] != '\n') b->buf[n++] = '\n';   // a last line without its newline still counts
                cut = n;
            } else cut = last_record_end(b->buf, n);
            carry.assign(b->buf + cut, b->buf + n);
            b->len = cut;
            std::lock_guard<std::mutex> l(m);
            if (cut) { b->id = made++; todo.push_back(b); } else pool.push_back(b);
            cv.notify_all();
        }
        std::lock_guard<std::mutex> l(m); fed = true; cv.notify_all();
    });
    std::vector<std::thread> ws;
    for (int i = 0; i < et; i++) ws.emplace_back([&] {
        for (;;) {
            pblock* b;
            {
                std::unique_lock<std::mutex> l(m);
                cv.wait(l, [&] { return !todo.empty() || fed; });
                if (todo.empty()) return;
                b = todo.front(); todo.pop_front();
            }
            encode(b);
            std::lock_guard<std::mutex> l(m); done[b->id] = b; cv.notify_all();
        }
    });
    uint64_t nr = 0, tot = 0, mx = 0, pos = sizeof(bin_header);
    std::vector<uint64_t> off;
    for (uint64_t id = 0;; id++) {   // write the blocks in input order
        pblock* b;
        {
            std::unique_lock<std::mutex> l(m);
            cv.wait(l, [&] { return done.count(id) || (fed && todo.empty() && id >= made); });
            if (!done.count(id)) break;
            b = done[id]; done.erase(id);
        }
        off.clear();
        for (uint32_t s : b->size) { off.push_back(pos); pos += s; }
        std::fwrite(off.data(), 8, off.size(), fx);
        std::fwrite(b->out.data(), 1, b->out.size(), fo);
        nr += b->size.size(); tot += b->bases; mx = std::max(mx, b->mx);
        std::lock_guard<std::mutex> l(m); pool.push_back(b); cv.notify_all();
    }
    reader.join(); for (auto& w : ws) w.join();
    for (pblock* b : pool) { delete[] b->buf; delete b; }
    std::fwrite(&pos, 8, 1, fx);
    bin_header ho = header('1', 1, nr, tot), hx = header('I', 1, nr, mx);
    finish(fo, out, &ho);
    finish(fx, x.c_str(), &hx);
    std::fprintf(stderr, "[pack] %lu reads, %lu bases\n", nr, tot);
    return 0;
}
