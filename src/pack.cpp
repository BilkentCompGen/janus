#include <zlib-ng.h>
#include <algorithm>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include "cmd.h"
#include "io.h"

// pack: a FASTQ (plain or gzip, any number of gzip members) into an svb and its index. Reads are stored in input order:
// a 4-byte tag (length, top bit set for a raw read) and the bases at 2 bits each, or the raw bytes when the read has a
// character other than ACGT (either case). The index holds each read's offset in the svb and the end offset.
// One thread inflates (zlib-ng) and cuts the text into CHUNK-sized blocks that end on a record, T workers encode the
// blocks, and the caller writes them in order, so the output does not depend on the thread count.

static const size_t CHUNK = 64 << 20, IN = 16 << 20;

struct pblock { uint64_t id; std::string text; std::vector<uint8_t> out; std::vector<uint32_t> size; uint64_t bases = 0, mx = 0; };

// the text source: inflates gzip members one after the other, or passes a plain file through
struct source {
    FILE* f; bool gz; zng_stream z{}; std::vector<uint8_t> in; bool eof = false;
    explicit source(const char* path) : f(std::fopen(path, "rb")), in(IN) {
        if (!f) { std::fprintf(stderr, "[pack] cannot open %s\n", path); std::exit(1); }
        int a = std::fgetc(f), b = std::fgetc(f); gz = a == 0x1f && b == 0x8b; std::rewind(f);
        if (gz && zng_inflateInit2(&z, 15 + 32) != Z_OK) { std::fprintf(stderr, "[pack] inflate init failed\n"); std::exit(1); }
    }
    ~source() { if (gz) zng_inflateEnd(&z); std::fclose(f); }
    size_t read(char* o, size_t n) {   // up to n bytes of text; 0 at the end
        if (!gz) return std::fread(o, 1, n, f);
        z.next_out = (uint8_t*)o; z.avail_out = n;
        while (z.avail_out) {
            if (!z.avail_in) {
                if (eof) break;
                z.avail_in = std::fread(in.data(), 1, in.size(), f); z.next_in = in.data();
                if (!z.avail_in) { eof = true; break; }
            }
            int r = zng_inflate(&z, Z_NO_FLUSH);
            if (r == Z_STREAM_END) zng_inflateReset(&z);   // another gzip member may follow
            else if (r != Z_OK && r != Z_BUF_ERROR) { std::fprintf(stderr, "[pack] corrupt gzip input (%d)\n", r); std::exit(1); }
        }
        return n - z.avail_out;
    }
};

// encodes every complete 4-line record of the block's text (a trailing incomplete record is dropped, as before)
static void encode(pblock* b) {
    const char *p = b->text.data(), *e = p + b->text.size();
    auto next = [&](const char*& s, size_t& l) {   // one line without its trailing \r / \n
        if (p >= e) return false;
        const char* q = (const char*)std::memchr(p, '\n', e - p); if (!q) q = e;
        s = p; l = q - p; p = q + 1;
        while (l && s[l - 1] == '\r') l--;
        return true;
    };
    const char *h, *s, *x, *y; size_t hl, sl, xl, yl;
    b->out.reserve(b->text.size() / 7);
    while (next(h, hl) && next(s, sl) && next(x, xl) && next(y, yl)) {
        uint32_t L = sl, raw = 0;
        for (size_t i = 0; i < sl; i++) raw |= code(s[i]) == 3 && (s[i] & 0xdf) != 'T';
        uint32_t tag = L | raw << 31; size_t n = svb_bytes(tag), at = b->out.size();
        b->out.resize(at + n, 0);
        std::memcpy(&b->out[at], &tag, 4);
        uint8_t* d = &b->out[at + 4];
        if (raw) std::memcpy(d, s, L);
        else for (uint32_t i = 0; i < L; i++) d[i >> 2] |= code(s[i]) << (6 - 2 * (i & 3));
        b->size.push_back(n); b->bases += L; b->mx = std::max<uint64_t>(b->mx, L);
    }
    std::string().swap(b->text);
}

int cmd_pack(const char* in, const char* out, int t) {
    t = std::max(1, t);
    source src(in);
    std::string x = std::string(out) + ".idx";
    FILE *fo = open_tmp(out), *fx = open_tmp(x.c_str());
    std::mutex m; std::condition_variable cv;
    std::deque<pblock*> todo; std::map<uint64_t, pblock*> done; bool fed = false; uint64_t made = 0;
    const size_t MAXQ = 2 * (size_t)t + 2;   // blocks in flight: bounds memory at about MAXQ x CHUNK x 2

    std::thread reader([&] {   // cut the text after the last newline that ends a record
        std::string carry; uint64_t lines = 0;   // newlines seen before the current buffer, mod 4
        std::vector<char> buf(CHUNK);
        for (;;) {
            size_t n = src.read(buf.data(), buf.size());
            std::string text = std::move(carry); carry.clear();
            text.append(buf.data(), n);
            size_t cut = 0;
            if (n) {
                uint64_t k = lines; const char* p = text.data(); const char* e = p + text.size();
                for (const char* q; (q = (const char*)std::memchr(p, '\n', e - p)); p = q + 1) if (++k % 4 == 0) cut = q + 1 - text.data();
                lines = 0; carry.assign(text, cut, std::string::npos); text.resize(cut);
            } else {
                if (!text.empty() && text.back() != '\n') text += '\n';   // a last line without its newline still counts
            }
            if (!text.empty()) {
                pblock* b = new pblock; b->text = std::move(text);
                std::unique_lock<std::mutex> l(m);
                cv.wait(l, [&] { return todo.size() + done.size() < MAXQ; });
                b->id = made++; todo.push_back(b); cv.notify_all();
            }
            if (!n) break;
        }
        std::lock_guard<std::mutex> l(m); fed = true; cv.notify_all();
    });
    std::vector<std::thread> ws;
    for (int i = 0; i < t; i++) ws.emplace_back([&] {
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
            b = done[id]; done.erase(id); cv.notify_all();
        }
        off.clear();
        for (uint32_t s : b->size) { off.push_back(pos); pos += s; }
        std::fwrite(off.data(), 8, off.size(), fx);
        std::fwrite(b->out.data(), 1, b->out.size(), fo);
        nr += b->size.size(); tot += b->bases; mx = std::max(mx, b->mx);
        delete b;
    }
    reader.join(); for (auto& w : ws) w.join();
    std::fwrite(&pos, 8, 1, fx);
    bin_header ho = header('1', 1, nr, tot), hx = header('I', 1, nr, mx);
    finish(fo, out, &ho);
    finish(fx, x.c_str(), &hx);
    std::fprintf(stderr, "[pack] %lu reads, %lu bases\n", nr, tot);
    return 0;
}
