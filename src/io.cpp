#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <pthread.h>
#include <unistd.h>
#include "io.h"

struct Seat { Job* job; int i; };

void* Job::enter(void* seat) {
    Seat* s = (Seat*)seat;
    s->job->task(s->i);
    return nullptr;
}

void Job::par(int n) {
    pthread_t th[MAXT + 1];
    Seat s[MAXT + 1];
    for (int i = 0; i < n; i++) { s[i] = {this, i}; pthread_create(&th[i], nullptr, enter, &s[i]); }
    for (int i = 0; i < n; i++) pthread_join(th[i], nullptr);
}

Header::Header(char kind, uint32_t ver, uint64_t n, uint64_t total) : magic{'S', 'V', 'B', kind}, ver(ver), n(n), total(total), resv{} {}

Input::Input(const std::string& path) : f(std::fopen(path.c_str(), "rb")) {
    if (!f || std::fread(&h, 1, sizeof(Header), f) != sizeof(Header)) { std::fprintf(stderr, "janus: cannot read %s\n", path.c_str()); std::exit(1); }
}

Output::Output(const char* path) : f(std::fopen((std::string(path) + ".tmp").c_str(), "wb+")), path(path) {
    if (!f) { std::fprintf(stderr, "janus: cannot write %s\n", path); std::exit(1); }
    std::setvbuf(f, nullptr, _IOFBF, 1 << 22);
    std::fseek(f, sizeof(Header), SEEK_SET);
}

void Output::finish(const Header& h) {
    std::fseek(f, 0, SEEK_SET);
    std::fwrite(&h, 1, sizeof(Header), f);
    std::fclose(f);
    std::rename((path + ".tmp").c_str(), path.c_str());
}

void Output::save(const char* path, const Header& h, const void* p, size_t n) {
    Output o(path);
    std::fwrite(p, 1, n, o.f);
    o.finish(h);
}

static const int NR = 8;

class Fill : public Job {
public:
    static size_t read(char* p, size_t n, FILE* f);

private:
    int fd;
    char* p;
    size_t n, got[NR];
    long o;
    Fill(int fd, char* p, size_t n, long o) : fd(fd), p(p), n(n), got(), o(o) {}
    void task(int i) override;
};

size_t Fill::read(char* p, size_t n, FILE* f) {
    Fill r(fileno(f), p, n, std::ftell(f));
    r.par(NR);
    size_t got = 0;
    for (size_t g : r.got) got += g;
    std::fseek(f, r.o + (long)got, SEEK_SET);
    return got;
}

void Fill::task(int i) {
    size_t at = n * i / NR, end = n * (i + 1) / NR;
    for (ssize_t g = 1; g > 0 && at < end; at += g) g = pread(fd, p + at, end - at, o + at);
    got[i] = at - n * i / NR;
}

static const size_t IBUF = 1 << 26;
void Stream::task(int i) {
    if (i < nb) { block(b[i]); return; }
    for (int j = 0; j < n; j++) {
        size_t k = s[j].len - s[j].off;
        std::memcpy(s[j].b[cur ^ 1], s[j].b[cur] + s[j].off, k);
        s[j].next = k + Fill::read(s[j].b[cur ^ 1] + k, cap - k, s[j].f);
    }
}

bool Stream::whole() {
    size_t z[2];
    for (int i = 0; i < n; i++) {
        const char* p = s[i].b[cur] + s[i].off;
        size_t left = s[i].len - s[i].off;
        if (!(z[i] = i ? exc_size(p, left) : svb_size(p, left))) return false;
    }
    for (int i = 0; i < n; i++) s[i].off += z[i];
    nrec++;
    return true;
}

void Stream::cut() {
    Block& x = b[nb++];
    x.rec = rec;
    rec = nrec;
    for (int i = 0; i < n; i++) {
        x.in[i] = s[i].b[cur] + s[i].start;
        x.inlen[i] = s[i].off - s[i].start;
        s[i].start = s[i].off;
    }
}

bool Stream::grow() {
    int i;
    for (i = 0; i < n && s[i].len < cap; i++);
    if (i == n) return false;
    cap *= 2;
    for (i = 0; i < n; i++) {
        s[i].b[0] = (char*)std::realloc(s[i].b[0], cap + 8);
        s[i].b[1] = (char*)std::realloc(s[i].b[1], cap + 8);
        s[i].len += Fill::read(s[i].b[cur] + s[i].len, cap - s[i].len, s[i].f);
    }
    return true;
}

uint64_t Stream::pass(FILE* fi, FILE* fe, FILE* fo, int t) {
    s[0] = {fi, {}, 0, 0, 0, 0};
    s[1] = {fe, {}, 0, 0, 0, 0};
    n = fe ? 2 : 1;
    cur = 0;
    cap = IBUF;
    rec = nrec = 0;
    t = std::min(t, MAXT);
    for (int i = 0; i < t; i++) {
        b[i].id = i;
        if (fo) b[i].out.reserve(IBUF);
    }
    for (int i = 0; i < n; i++) {
        s[i].b[0] = (char*)std::malloc(IBUF + 8);
        s[i].b[1] = (char*)std::malloc(IBUF + 8);
        s[i].len = Fill::read(s[i].b[0], IBUF, s[i].f);
    }
    while (s[0].len) {
        nb = 0;
        for (int i = 0; i < n; i++) s[i].off = s[i].start = 0;
        while (whole())
            if (s[0].off - s[0].start > s[0].len / t && nb < t - 1) cut();
        if (s[0].off > s[0].start) cut();
        if (!nb) {
            if (!grow()) break;
            continue;
        }
        par(nb + 1);
        for (int i = 0; i < nb; i++) {
            if (fo) std::fwrite(b[i].out.data(), 1, b[i].out.size(), fo);
            b[i].out.clear();
        }
        cur ^= 1;
        for (int i = 0; i < n; i++) s[i].len = s[i].next;
    }
    for (int i = 0; i < n; i++) {
        std::fclose(s[i].f);
        std::free(s[i].b[0]);
        std::free(s[i].b[1]);
    }
    return nrec;
}

std::vector<std::string> split(const char* s) {
    std::vector<std::string> v;
    std::string p = s;
    for (size_t a = 0, b; a < p.size(); a = b + 1) {
        b = std::min(p.find(',', a), p.size());
        v.push_back(p.substr(a, b - a));
    }
    return v;
}

Lines::Lines(const char* path) : s(nullptr), n(0), f(std::fopen(path, "r")), cap(0) {
    if (!f) { std::fprintf(stderr, "janus: cannot read %s\n", path); std::exit(1); }
}

Lines::~Lines() {
    std::free(s);
    std::fclose(f);
}

bool Lines::next() {
    return (n = getline(&s, &cap, f)) > 0;
}
