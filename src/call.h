#pragma once
#include <climits>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include "cmd.h"
#include "files.h"

static const int32_t UNSET = INT32_MIN;
enum { SPAN, SAME, CONSISTENT, CONSISTENT_SAME, KINDS };

struct Spot {
    uint32_t a, b, s, t, tonly, nbr, lo;
    uint64_t read;
    std::vector<uint32_t> n, l, r, m, dep, pb;
    std::vector<std::vector<int32_t>> dl;
    int32_t d[NA][NA];
    std::string nm;
    uint8_t lk;
    void parse(const char* line, const Graph& g);
};

struct Align { uint32_t qs, qe, ts, te, mq; char st; std::string t, cg; float de; };

struct Event {
    std::string ty, c1, c2;
    uint32_t p1, p2, q, mq;
    float de;
    uint32_t o;
    char st;
    uint32_t pb;
};

struct Probe {
    uint32_t qb, tb, v[2][NA], at[2][NA], gap, absent;
    int32_t qa, ta, d[NA][NA];
    std::vector<std::vector<std::pair<int32_t, int32_t>>> dl;
    std::vector<std::pair<uint32_t, uint32_t>> mid;
    uint8_t lk, set;
};

std::vector<Event> events(std::vector<Align>& sg);

static const uint32_t NF = 56;
static const uint32_t NEAR = 500;

struct Row {
    uint32_t a, b, u, s, ncl, tspan, nspan, nsame, tonly, nbr, span;
    int32_t lest, dmin, dmax;
    uint32_t nspc, n20c, tspc;
    std::vector<uint32_t> n, d;
    Row(const Spot& e, uint32_t u);
};

struct Record {
    Event e;
    std::string nm;
    uint32_t pc;
    Row r;
    uint32_t nf, ms, mn;
};

uint32_t top(const std::vector<uint32_t>& n);
uint32_t misfit(const Row& r, const Event& e);
float repeats(const std::string& s);
float prob(const Record& s, const Probe* b, float rp, bool vn, float med, float* f);

class Bed {
public:
    explicit Bed(const char* path);
    bool near(const std::string& c, uint32_t p) const;
    bool span(const std::string& c, uint32_t p, uint32_t* s, uint32_t* e) const;

private:
    std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> m;
};

class Call {
public:
    explicit Call(const Args& a) : a(a), pre(a.in), bed(a.vntr), med(0) {}
    int run();

private:
    const Args& a;
    std::string pre;
    Bed bed;
    std::vector<Spot> w;
    std::vector<Probe> pbs;
    std::unordered_map<std::string, std::vector<uint32_t>> bypc;
    std::unordered_map<std::string, std::vector<Align>> al;
    std::unordered_map<std::string, std::vector<Event>> evs;
    std::unordered_map<std::string, std::string> sq;
    float med;
    void sites();
    void probes(uint64_t S);
    void alignments();
    void sequences();
    std::vector<Record> attach(size_t* found);
    void write(std::vector<Record>& res, size_t found);
};
