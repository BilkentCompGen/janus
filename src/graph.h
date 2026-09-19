#pragma once
#include "samples.h"
#include "snape.h"
#include "unordered_dense.h"
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <queue>
#include <string>
#include <vector>

struct Shard { std::string path; int col; };
struct Batch { std::vector<std::string> seqs; int col = 0; uint64_t first = 0; };

using SigNodeMap = ankerl::unordered_dense::map<uint64_t, uint32_t>;
using SigCntMap = ankerl::unordered_dense::map<uint64_t, cov_t>;
using EdgeKeySet = ankerl::unordered_dense::set<uint64_t>;

struct NodeTab { SigNodeMap m; std::vector<cov_t> a; };
struct ShTab { SigCntMap m; };

struct SigIndex {
    static constexpr int BITS = 24;
    std::vector<uint32_t> off;
    const SigCount* v = nullptr;
    void build(const std::vector<SigCount>& nodes);
    uint32_t find(uint64_t sig) const;
    size_t bytes() const;
};

struct PostIdx {
    uint32_t N = 0, pad = 0;
    std::vector<uint64_t> w;
    std::vector<uint8_t> cnt;
    std::vector<uint64_t> worst;
    std::mutex* mtx = nullptr;
    size_t bytes() const;
};

struct Feed {
    std::mutex qmtx;
    std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    std::atomic<int> readers_left{1};
    std::atomic<size_t> next{0};
    std::vector<uint64_t> reads_per;
    uint64_t max_bases = 0;
};

struct FillCtx {
    const SigIndex* idx = nullptr;
    const std::vector<cov_t>* cov = nullptr;
    PostIdx* px = nullptr;
    int S = 0;
};

struct LocalB { EdgeKeySet e; };

void run_graph_reads(const std::vector<Sample>& samples, const std::string& dir,
                     int level, int margin, int min_mult, int threads, uint64_t max_bases = 0,
                     int index_n = 8, int pad_cores = 50);
