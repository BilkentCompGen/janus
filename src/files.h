#pragma once
#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

class Graph {
public:
    uint64_t K, E, S;
    std::vector<uint64_t> hashes, edges;
    std::vector<uint16_t> counts;
    std::vector<uint8_t> novel;
    std::vector<uint32_t> starts;
    Graph() : K(0), E(0), S(0) {}
    explicit Graph(const char* path) { read(path); }
    void read(const char* path);
    void write(const char* path) const;
    uint16_t reads(uint64_t s, uint64_t v) const { return counts[s * K + v]; }
    static uint32_t from(uint64_t e) { return e >> 32; }
    static uint32_t to(uint64_t e) { return e & 0xffffffff; }
};

class Adjacency {
public:
    std::vector<uint32_t> first_out, target, first_in, source;
    explicit Adjacency(const Graph& g);
    uint32_t outs(uint32_t v) const { return first_out[v + 1] - first_out[v]; }
    uint32_t ins(uint32_t v) const { return first_in[v + 1] - first_in[v]; }
};

class Fasta {
public:
    std::vector<std::string> name, seq;
    explicit Fasta(const char* path);
    static void put(FILE* f, const std::string& name, const std::string& seq) { std::fprintf(f, ">%s\n%s\n", name.c_str(), seq.c_str()); }
};

struct FaiEntry { uint64_t len, off, lb, lw; };

class Fai {
public:
    std::unordered_map<std::string, FaiEntry> m;
    explicit Fai(const char* fasta);
    std::string range(FILE* fa, const std::string& c, uint32_t a, uint32_t b);
};

class Paf {
public:
    char q[512], st[4], t[512];
    unsigned ql, qs, qe, tl, ts, te, nm, bl, mq;
    int k4, k9;
    int parse(const char* line) { return std::sscanf(line, "%511s %u %u %u%n %3s %511s %u %u %u%n %u %u %u", q, &ql, &qs, &qe, &k4, st, t, &tl, &ts, &te, &k9, &nm, &bl, &mq); }
};
