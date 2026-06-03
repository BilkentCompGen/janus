// bubble_census — Stage 1 of inspect/.
//
// Walks every bubble in the svp graph, tracks node IDs (the canonical
// decompose() only records aggregate counts), looks up posmap coords for the
// color=3 anchors adjacent to each bubble, classifies against the truth VCF,
// flags HLA / centromere / segdup / simple-repeat regions from BED files.
//
// Outputs:
//   <dir>/census.tsv         one row per bubble (greppable, sortable)
//   <dir>/census.bin         packed binary for inspect/inspect.cpp
//   <dir>/bubble_summary.md  human-readable WG summary
//
// Usage:  bubble_census <dir> [<truth_vcf>] [<bed_dir>]

#include "common.h"

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace insp;

namespace {

constexpr int MIN_MULT = 3;          // matches make run default + project memory
constexpr int TRUTH_PAD = 500;       // truvari default --refdist
constexpr int LARGE_BUBBLE_THRESHOLD = 100;
constexpr int SINGLETON_THRESHOLD = 2;

// ---------- truth VCF ----------
struct TruthRec {
    std::string chrom;
    uint32_t pos;     // 0-based
    uint32_t end;     // 0-based, exclusive (for INS, == pos)
    int32_t  svlen;
    std::string svtype;
    std::string sv_id;
    int vset;
};

void load_truth_vcf(const std::string& path, std::vector<TruthRec>& out) {
    std::ifstream f(path);
    if (!f) die("cannot open truth VCF");
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> t;
        size_t b = 0;
        for (size_t i = 0; i <= line.size(); i++) {
            if (i == line.size() || line[i] == '\t') {
                t.push_back(line.substr(b, i - b));
                b = i + 1;
            }
        }
        if (t.size() < 8) continue;
        TruthRec r{};
        r.chrom = t[0];
        r.pos = (uint32_t)(std::stoul(t[1]) - 1);  // VCF is 1-based
        r.sv_id = t[2];
        std::string info = t[7];
        std::map<std::string, std::string> kv;
        size_t kvb = 0;
        for (size_t i = 0; i <= info.size(); i++) {
            if (i == info.size() || info[i] == ';') {
                std::string item = info.substr(kvb, i - kvb);
                size_t eq = item.find('=');
                if (eq != std::string::npos) kv[item.substr(0, eq)] = item.substr(eq + 1);
                kvb = i + 1;
            }
        }
        r.svtype = kv.count("SVTYPE") ? kv["SVTYPE"] : "";
        if (r.svtype.empty()) continue;
        r.svlen = kv.count("SVLEN") ? std::stoi(kv["SVLEN"]) : 0;
        r.end = kv.count("END") ? (uint32_t)(std::stoul(kv["END"]) - 1) : r.pos;
        r.vset = kv.count("VSET") ? std::stoi(kv["VSET"]) : 0;
        out.push_back(std::move(r));
    }
}

// ---------- BED region tables (with max-end-prefix for O(log n + k) overlap) ----------
struct BedIv { uint32_t start, end; };
struct BedTable {
    std::unordered_map<std::string, std::vector<BedIv>> by_chrom;
    std::unordered_map<std::string, std::vector<uint32_t>> max_end_prefix;
    bool loaded = false;
    void load(const std::string& path) {
        std::ifstream f(path);
        if (!f) return;
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::vector<std::string> t;
            size_t b = 0;
            for (size_t i = 0; i <= line.size(); i++) {
                if (i == line.size() || line[i] == '\t') {
                    t.push_back(line.substr(b, i - b));
                    b = i + 1;
                }
            }
            if (t.size() < 3) continue;
            BedIv iv{(uint32_t)std::stoul(t[1]), (uint32_t)std::stoul(t[2])};
            by_chrom[t[0]].push_back(iv);
        }
        for (auto& kv : by_chrom) {
            auto& v = kv.second;
            std::sort(v.begin(), v.end(),
                      [](const BedIv& a, const BedIv& b){ return a.start < b.start; });
            auto& m = max_end_prefix[kv.first];
            m.resize(v.size());
            uint32_t mx = 0;
            for (size_t i = 0; i < v.size(); i++) { mx = std::max(mx, v[i].end); m[i] = mx; }
        }
        loaded = true;
    }
    bool overlaps(const std::string& chrom, uint32_t lo, uint32_t hi) const {
        auto it = by_chrom.find(chrom);
        if (it == by_chrom.end()) return false;
        const auto& v = it->second;
        auto mit = max_end_prefix.find(chrom);
        const auto& m = mit->second;
        auto u = std::lower_bound(v.begin(), v.end(), hi,
                                  [](const BedIv& iv, uint32_t k){ return iv.start < k; });
        if (u == v.begin()) return false;
        size_t idx = (size_t)(u - v.begin()) - 1;
        if (m[idx] <= lo) return false;
        // walk back while the prefix-max can still produce an overlap
        while (true) {
            if (v[idx].end > lo && v[idx].start < hi) return true;
            if (idx == 0) return false;
            if (m[idx - 1] <= lo) return false;
            idx--;
        }
    }
};

// ---------- census row (binary on-disk) ----------
struct CensusRow {
    uint32_t bubble_id;
    uint8_t  type;          // 'D' / 'I' / 'C'
    uint8_t  flags;         // bit0=in_hla, 1=in_cent, 2=in_segdup, 3=in_simple, 4=in_gap
    uint8_t  anchorable;    // bool
    uint8_t  mixed_strand;
    uint8_t  mixed_chrom;
    uint8_t  truth_match;   // bool
    uint8_t  singleton;     // bool (n_nodes <= 2)
    uint8_t  large;         // bool (n_nodes > 100)
    uint32_t n_nodes;
    uint32_t n_normal;
    uint32_t n_tumor;
    uint32_t n_anchors;
    uint64_t n_cov_sum;
    uint64_t t_cov_sum;
    uint32_t n_cov_max;
    uint32_t t_cov_max;
    int32_t  chrom_idx;     // -1 if not anchored
    uint8_t  L_strand;
    uint8_t  R_strand;
    uint16_t _pad;
    uint32_t L_pos;
    uint32_t R_pos;
    int32_t  truth_svlen;
    uint32_t truth_pos;
    uint32_t truth_end;
    uint8_t  truth_svtype;  // 'D'/'I'/'P'(dup)/'V'(inv)/'S'(snp) or 0
    uint8_t  _pad2[3];
};
static_assert(sizeof(CensusRow) == 88, "CensusRow size");

// ---------- node_to_sig from counts.bin (streaming, applies MIN_MULT filter) ----------
void build_node_to_sig(const std::string& path, uint64_t K, std::vector<uint64_t>& node_to_sig) {
    FILE* f = xopen(path, "rb");
    uint64_t N; xread(f, &N, 8);
    node_to_sig.resize(K);
    uint64_t kept = 0;
    SigCount rec;
    for (uint64_t i = 0; i < N; i++) {
        xread(f, &rec, sizeof(rec));
        if (std::max(rec.n, rec.t) < (uint32_t)MIN_MULT) continue;
        if (kept >= K) die("counts.bin: more kept sigs than graph.bin K — inconsistent files");
        node_to_sig[kept++] = rec.sig;
    }
    std::fclose(f);
    if (kept != K) {
        std::fprintf(stderr, "[FATAL] counts.bin kept=%lu but graph.bin K=%lu\n",
                     (unsigned long)kept, (unsigned long)K);
        std::abort();
    }
}

// ---------- replay decompose() but track node IDs and anchor IDs ----------
struct Bubble {
    char type;                       // 'D'/'I'/'C'
    std::vector<uint32_t> node_ids;
    std::vector<uint32_t> anchor_ids;
    uint32_t n_normal, n_tumor;
    uint64_t n_cov_sum, t_cov_sum;
    uint32_t n_cov_max, t_cov_max;
};

void decompose_with_ids(const std::vector<std::vector<int32_t>>& adj,
                        const std::vector<uint8_t>& color,
                        const std::vector<uint32_t>& ncov,
                        const std::vector<uint32_t>& tcov,
                        std::vector<Bubble>& out) {
    const uint32_t K = (uint32_t)color.size();
    std::vector<char> vis(K, 0);
    std::queue<uint32_t> q;
    std::vector<uint32_t> anchor_buf;       // reused across bubbles
    anchor_buf.reserve(64);
    for (uint32_t s = 0; s < K; s++) {
        if (color[s] == 3 || vis[s]) continue;
        q.push(s); vis[s] = 1;
        Bubble b{};
        uint8_t seen = 0;
        anchor_buf.clear();
        while (!q.empty()) {
            const uint32_t u = q.front(); q.pop();
            b.node_ids.push_back(u);
            seen |= color[u];
            if (color[u] == 1) b.n_normal++;
            else if (color[u] == 2) b.n_tumor++;
            b.n_cov_sum += ncov[u]; b.t_cov_sum += tcov[u];
            if (ncov[u] > b.n_cov_max) b.n_cov_max = ncov[u];
            if (tcov[u] > b.t_cov_max) b.t_cov_max = tcov[u];
            for (int32_t v : adj[u]) {
                if (color[v] == 3) {
                    anchor_buf.push_back((uint32_t)v);
                } else if (!vis[v]) {
                    vis[v] = 1; q.push((uint32_t)v);
                }
            }
        }
        b.type = (seen == 1 ? 'D' : seen == 2 ? 'I' : 'C');
        std::sort(anchor_buf.begin(), anchor_buf.end());
        anchor_buf.erase(std::unique(anchor_buf.begin(), anchor_buf.end()), anchor_buf.end());
        b.anchor_ids = anchor_buf;
        out.push_back(std::move(b));
    }
}

// ---------- anchor analysis: pick (L, R) pair from posmap ----------
struct AnchorAnalysis {
    bool anchorable = false;
    bool mixed_strand = false;
    bool mixed_chrom = false;
    int32_t chrom_idx = -1;
    uint32_t L_pos = 0, R_pos = 0;
    uint8_t L_strand = 0, R_strand = 0;
    std::string chrom_name;
};

void analyze_anchors(const Bubble& bub,
                     const std::vector<uint64_t>& node_to_sig,
                     const Posmap& pm,
                     AnchorAnalysis& aa) {
    // Collect every posmap entry that any anchor's sig points to
    struct Hit { uint16_t ci; uint32_t ps, pe; uint8_t strand; };
    std::vector<Hit> hits;
    for (uint32_t aid : bub.anchor_ids) {
        uint64_t sig = node_to_sig[aid];
        auto it = pm.sig_pos.find(sig);
        if (it == pm.sig_pos.end()) continue;
        for (const Pos& p : it->second) {
            hits.push_back({p.chrom_idx, p.pos_start, p.pos_end, p.strand});
        }
    }
    if (hits.empty()) return;

    // Group by chrom: pick chrom with most hits
    std::map<uint16_t, std::vector<Hit>> by_chrom;
    for (const Hit& h : hits) by_chrom[h.ci].push_back(h);
    if (by_chrom.size() > 1) aa.mixed_chrom = true;
    auto best = by_chrom.begin();
    for (auto it = by_chrom.begin(); it != by_chrom.end(); ++it)
        if (it->second.size() > best->second.size()) best = it;
    auto& bh = best->second;

    // Strand consensus within picked chrom
    std::map<uint8_t, int> by_strand;
    for (const Hit& h : bh) by_strand[h.strand]++;
    if (by_strand.size() > 1) aa.mixed_strand = true;
    uint8_t best_strand = 0; int best_strand_count = -1;
    for (auto& kv : by_strand) if (kv.second > best_strand_count) { best_strand_count = kv.second; best_strand = kv.first; }

    // Pick L = min pos, R = max pos on best_strand subset
    uint32_t lo = UINT32_MAX, hi = 0;
    for (const Hit& h : bh) {
        if (h.strand != best_strand) continue;
        if (h.ps < lo) lo = h.ps;
        if (h.pe > hi) hi = h.pe;
    }
    if (lo == UINT32_MAX) return;
    if ((uint64_t)hi - (uint64_t)lo > 1'000'000ULL) return;  // span > 1 Mbp: drop

    aa.anchorable = true;
    aa.chrom_idx = (int32_t)best->first;
    aa.L_pos = lo;
    aa.R_pos = hi;
    aa.L_strand = best_strand;
    aa.R_strand = best_strand;
    if ((size_t)aa.chrom_idx < pm.chrom_names.size())
        aa.chrom_name = pm.chrom_names[aa.chrom_idx];
}

// ---------- truth match: best record within ±TRUTH_PAD on anchor chrom ----------
struct TruthHit {
    bool matched = false;
    int32_t  svlen = 0;
    uint32_t pos = 0, end = 0;
    char svtype = 0;
    std::string sv_id;
};

void match_truth(const std::string& chrom, uint32_t lo, uint32_t hi,
                 const std::vector<TruthRec>& truth,
                 std::map<std::string, std::vector<size_t>>& truth_by_chrom,
                 TruthHit& th) {
    auto it = truth_by_chrom.find(chrom);
    if (it == truth_by_chrom.end()) return;
    const auto& idxs = it->second;  // already sorted by pos
    // binary search for first truth.pos >= lo - TRUTH_PAD
    uint32_t lo_pad = lo > TRUTH_PAD ? lo - TRUTH_PAD : 0;
    uint32_t hi_pad = hi + TRUTH_PAD;
    auto first = std::lower_bound(idxs.begin(), idxs.end(), lo_pad,
                                  [&](size_t i, uint32_t k){ return truth[i].pos < k; });
    double best_score = -1;
    for (auto p = first; p != idxs.end(); ++p) {
        const TruthRec& r = truth[*p];
        if (r.pos > hi_pad) break;
        if (r.vset == 0) continue;  // skip germline SNPs (VSET=0)
        // simple score: 1.0 - normalized distance
        double dist = std::min(std::abs((int64_t)r.pos - (int64_t)lo),
                               std::abs((int64_t)r.end - (int64_t)hi));
        double score = 1.0 / (1.0 + dist / 100.0);
        if (score > best_score) {
            best_score = score;
            th.matched = true;
            th.svlen = r.svlen;
            th.pos = r.pos;
            th.end = r.end;
            th.svtype = r.svtype == "DEL" ? 'D' :
                        r.svtype == "INS" ? 'I' :
                        r.svtype == "DUP" ? 'P' :
                        r.svtype == "INV" ? 'V' :
                        r.svtype == "SNP" ? 'S' : '?';
            th.sv_id = r.sv_id;
        }
    }
}

// ---------- summary writer ----------
void write_summary(const std::vector<CensusRow>& rows,
                   const std::vector<std::string>& truth_ids,
                   const Posmap& pm,
                   const std::string& path) {
    FILE* f = xopen(path, "w");
    std::fprintf(f, "# Bubble census summary (svp graph-inspector, Stage 1)\n\n");
    std::fprintf(f, "Total bubbles: %zu\n\n", rows.size());

    int n_D = 0, n_I = 0, n_C = 0;
    int n_anchored = 0, n_mixed_strand = 0, n_mixed_chrom = 0;
    int n_truth = 0, n_singleton = 0, n_large = 0;
    int n_hla = 0, n_cent = 0, n_segdup = 0, n_simple = 0, n_gap = 0;
    std::map<int, int> size_hist;          // 1, 2, 3-5, 6-10, 11-50, 51-500, 501+
    std::map<std::string, int> by_chrom;
    std::map<std::string, int> by_chrom_truth;
    std::map<std::string, int> by_chrom_FP;
    int n_truth_by_type[5] = {0, 0, 0, 0, 0};  // D/I/P/V/S
    const char* type_chars = "DIPVS";

    auto bucket = [](uint32_t n) -> int {
        if (n == 1) return 0;
        if (n == 2) return 1;
        if (n <= 5) return 2;
        if (n <= 10) return 3;
        if (n <= 50) return 4;
        if (n <= 500) return 5;
        return 6;
    };
    const char* bucket_lbl[] = {"1", "2", "3-5", "6-10", "11-50", "51-500", "501+"};

    for (size_t i = 0; i < rows.size(); i++) {
        const CensusRow& r = rows[i];
        if (r.type == 'D') n_D++;
        else if (r.type == 'I') n_I++;
        else if (r.type == 'C') n_C++;
        if (r.anchorable) n_anchored++;
        if (r.mixed_strand) n_mixed_strand++;
        if (r.mixed_chrom) n_mixed_chrom++;
        if (r.truth_match) {
            n_truth++;
            for (int j = 0; j < 5; j++) if (r.truth_svtype == type_chars[j]) n_truth_by_type[j]++;
        }
        if (r.singleton) n_singleton++;
        if (r.large) n_large++;
        if (r.flags & 1) n_hla++;
        if (r.flags & 2) n_cent++;
        if (r.flags & 4) n_segdup++;
        if (r.flags & 8) n_simple++;
        if (r.flags & 16) n_gap++;
        size_hist[bucket(r.n_nodes)]++;
        if (r.chrom_idx >= 0 && (size_t)r.chrom_idx < pm.chrom_names.size()) {
            const std::string& ch = pm.chrom_names[r.chrom_idx];
            by_chrom[ch]++;
            if (r.truth_match) by_chrom_truth[ch]++;
            else by_chrom_FP[ch]++;
        }
    }

    std::fprintf(f, "## Types\n\n| type | count | %% |\n|---|---|---|\n");
    std::fprintf(f, "| D (normal-only) | %d | %.2f%% |\n", n_D, 100.0 * n_D / rows.size());
    std::fprintf(f, "| I (tumor-only / somatic candidate) | %d | %.2f%% |\n", n_I, 100.0 * n_I / rows.size());
    std::fprintf(f, "| C (mixed) | %d | %.2f%% |\n", n_C, 100.0 * n_C / rows.size());

    std::fprintf(f, "\n## Truth matching (vs. somatic VSET>0 within ±%d bp)\n\n", TRUTH_PAD);
    std::fprintf(f, "TP bubbles: %d / %zu = %.2f%%\n", n_truth, rows.size(), 100.0 * n_truth / rows.size());
    std::fprintf(f, "FP (no truth match, anchored): %d\n", n_anchored - n_truth);
    std::fprintf(f, "UNMAPPED (not anchorable): %d\n\n", (int)rows.size() - n_anchored);
    std::fprintf(f, "Truth match by svtype:\n");
    const char* type_names[] = {"DEL", "INS", "DUP", "INV", "SNP"};
    for (int j = 0; j < 5; j++) std::fprintf(f, "  %s: %d\n", type_names[j], n_truth_by_type[j]);

    std::fprintf(f, "\n## Anchorability\n\n");
    std::fprintf(f, "anchorable (L,R same chrom+strand, span <= 1 Mb): %d / %zu = %.2f%%\n",
                 n_anchored, rows.size(), 100.0 * n_anchored / rows.size());
    std::fprintf(f, "mixed_strand (anchors disagree on strand — INV candidate): %d\n", n_mixed_strand);
    std::fprintf(f, "mixed_chrom (anchors on multiple chroms — BND candidate): %d\n", n_mixed_chrom);

    std::fprintf(f, "\n## Size distribution\n\n| n_nodes | count |\n|---|---|\n");
    for (auto& kv : size_hist)
        std::fprintf(f, "| %s | %d |\n", bucket_lbl[kv.first], kv.second);

    std::fprintf(f, "\n## Special sizes\n\n");
    std::fprintf(f, "singleton (n_nodes <= %d): %d (%.2f%%) — suspected RC-pair / residual error\n",
                 SINGLETON_THRESHOLD, n_singleton, 100.0 * n_singleton / rows.size());
    std::fprintf(f, "large (n_nodes > %d): %d (%.2f%%) — suspected tangle\n",
                 LARGE_BUBBLE_THRESHOLD, n_large, 100.0 * n_large / rows.size());

    std::fprintf(f, "\n## Region overlap (from BED assets)\n\n");
    std::fprintf(f, "| region | count | FP rate (within region) |\n|---|---|---|\n");
    auto fpr = [&](uint8_t mask) {
        int total = 0, fp = 0;
        for (const CensusRow& r : rows) if (r.flags & mask && r.anchorable) {
            total++;
            if (!r.truth_match) fp++;
        }
        if (total == 0) return std::string("n/a");
        std::ostringstream os;
        os.precision(2);
        os << std::fixed << (100.0 * fp / total) << "%";
        return os.str();
    };
    std::fprintf(f, "| HLA | %d | %s |\n", n_hla, fpr(1).c_str());
    std::fprintf(f, "| centromere | %d | %s |\n", n_cent, fpr(2).c_str());
    std::fprintf(f, "| segdup | %d | %s |\n", n_segdup, fpr(4).c_str());
    std::fprintf(f, "| simple repeat | %d | %s |\n", n_simple, fpr(8).c_str());
    std::fprintf(f, "| gap | %d | %s |\n", n_gap, fpr(16).c_str());

    // Clean (no region flags) FP rate as baseline
    {
        int total = 0, fp = 0;
        for (const CensusRow& r : rows) if (r.flags == 0 && r.anchorable) { total++; if (!r.truth_match) fp++; }
        if (total)
            std::fprintf(f, "| clean (no flags) | %d | %.2f%% |\n", total, 100.0 * fp / total);
    }

    std::fprintf(f, "\n## Top-20 chroms by FP density (anchored bubbles)\n\n");
    std::vector<std::pair<std::string, double>> fp_rates;
    for (auto& kv : by_chrom) {
        int total = kv.second;
        int fp = by_chrom_FP[kv.first];
        if (total >= 10) fp_rates.emplace_back(kv.first, 100.0 * fp / total);
    }
    std::sort(fp_rates.begin(), fp_rates.end(),
              [](const auto& a, const auto& b){ return a.second > b.second; });
    std::fprintf(f, "| chrom | bubbles | TP | FP | FP%% |\n|---|---|---|---|---|\n");
    for (size_t i = 0; i < fp_rates.size() && i < 20; i++) {
        const std::string& ch = fp_rates[i].first;
        int total = by_chrom[ch];
        int tp = by_chrom_truth[ch];
        int fp = by_chrom_FP[ch];
        std::fprintf(f, "| %s | %d | %d | %d | %.2f%% |\n", ch.c_str(), total, tp, fp, 100.0 * fp / total);
    }

    (void)truth_ids;
    std::fclose(f);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <svp_output_dir> [<truth_vcf>] [<bed_dir>]\n", argv[0]);
        return 1;
    }
    const std::string dir = argv[1];
    const std::string truth_vcf = argc > 2 ? argv[2] : "";
    const std::string bed_dir = argc > 3 ? argv[3] : "";

    std::fprintf(stderr, "[census] loading graph.bin ... ");
    Graph g; load_graph(dir + "/graph.bin", g);
    std::fprintf(stderr, "K=%lu E=%lu\n", (unsigned long)g.K, (unsigned long)g.edges.size());

    std::fprintf(stderr, "[census] building adjacency ... ");
    std::vector<std::vector<int32_t>> adj(g.K);
    for (const auto& e : g.edges) {
        adj[e.first].push_back((int32_t)e.second);
        adj[e.second].push_back((int32_t)e.first);
    }
    std::fprintf(stderr, "done\n");

    std::fprintf(stderr, "[census] building node_to_sig from counts.bin ... ");
    std::vector<uint64_t> node_to_sig;
    build_node_to_sig(dir + "/counts.bin", g.K, node_to_sig);
    std::fprintf(stderr, "done\n");

    std::fprintf(stderr, "[census] loading posmap.bin ... ");
    Posmap pm; load_posmap(dir + "/posmap.bin", pm);
    std::fprintf(stderr, "chroms=%zu sigs_with_pos=%zu\n", pm.chrom_names.size(), pm.sig_pos.size());

    std::vector<TruthRec> truth;
    std::map<std::string, std::vector<size_t>> truth_by_chrom;
    if (!truth_vcf.empty()) {
        std::fprintf(stderr, "[census] loading truth VCF ... ");
        load_truth_vcf(truth_vcf, truth);
        for (size_t i = 0; i < truth.size(); i++) truth_by_chrom[truth[i].chrom].push_back(i);
        for (auto& kv : truth_by_chrom) {
            std::sort(kv.second.begin(), kv.second.end(),
                      [&](size_t a, size_t b){ return truth[a].pos < truth[b].pos; });
        }
        std::fprintf(stderr, "%zu records, %zu chroms\n", truth.size(), truth_by_chrom.size());
    }

    BedTable hla, cent, segdup, simple, gap;
    if (!bed_dir.empty()) {
        std::fprintf(stderr, "[census] loading BED assets ... ");
        hla.load(bed_dir + "/hg38_hla.bed");
        cent.load(bed_dir + "/hg38_centromere.bed");
        segdup.load(bed_dir + "/hg38_segdup.bed");
        simple.load(bed_dir + "/hg38_simpleRepeat.bed");
        gap.load(bed_dir + "/hg38_gap.bed");
        std::fprintf(stderr, "loaded\n");
    }

    std::fprintf(stderr, "[census] running decompose with node-id tracking ... ");
    std::vector<Bubble> bubs;
    decompose_with_ids(adj, g.color, g.ncov, g.tcov, bubs);
    std::fprintf(stderr, "%zu bubbles\n", bubs.size());

    // ---------- per-bubble analysis + write ----------
    const std::string out_bin = dir + "/census.bin";
    const std::string out_tsv = dir + "/census.tsv";
    FILE* fbin = xopen(out_bin + ".tmp", "wb");
    FILE* ftsv = xopen(out_tsv, "w");
    std::fputs("bubble_id\ttype\tn_nodes\tn_normal\tn_tumor\tn_anchors\tn_cov_sum\tt_cov_sum\t"
               "n_cov_max\tt_cov_max\tanchorable\tmixed_strand\tmixed_chrom\tchrom\tL_pos\tR_pos\tL_strand\tR_strand\t"
               "truth_match\ttruth_id\ttruth_svtype\ttruth_svlen\ttruth_pos\ttruth_end\t"
               "in_hla\tin_centromere\tin_segdup\tin_simple_repeat\tin_gap\tsingleton\tlarge\n", ftsv);

    uint64_t nrows = bubs.size();
    xread(fbin, nullptr, 0);  // no-op
    std::fwrite(&nrows, 8, 1, fbin);

    std::vector<CensusRow> rows;
    rows.reserve(bubs.size());
    std::vector<std::string> truth_ids;
    truth_ids.reserve(bubs.size());

    uint32_t bubble_id = 0;
    for (const Bubble& b : bubs) {
        AnchorAnalysis aa;
        analyze_anchors(b, node_to_sig, pm, aa);

        CensusRow r{};
        r.bubble_id = bubble_id;
        r.type = (uint8_t)b.type;
        r.n_nodes = (uint32_t)b.node_ids.size();
        r.n_normal = b.n_normal;
        r.n_tumor = b.n_tumor;
        r.n_anchors = (uint32_t)b.anchor_ids.size();
        r.n_cov_sum = b.n_cov_sum;
        r.t_cov_sum = b.t_cov_sum;
        r.n_cov_max = b.n_cov_max;
        r.t_cov_max = b.t_cov_max;
        r.anchorable = aa.anchorable ? 1 : 0;
        r.mixed_strand = aa.mixed_strand ? 1 : 0;
        r.mixed_chrom = aa.mixed_chrom ? 1 : 0;
        r.chrom_idx = aa.chrom_idx;
        r.L_pos = aa.L_pos;
        r.R_pos = aa.R_pos;
        r.L_strand = aa.L_strand;
        r.R_strand = aa.R_strand;
        r.singleton = (r.n_nodes <= SINGLETON_THRESHOLD) ? 1 : 0;
        r.large = (r.n_nodes > LARGE_BUBBLE_THRESHOLD) ? 1 : 0;

        TruthHit th;
        if (aa.anchorable && !aa.chrom_name.empty())
            match_truth(aa.chrom_name, aa.L_pos, aa.R_pos, truth, truth_by_chrom, th);
        r.truth_match = th.matched ? 1 : 0;
        r.truth_svtype = (uint8_t)th.svtype;
        r.truth_svlen = th.svlen;
        r.truth_pos = th.pos;
        r.truth_end = th.end;
        truth_ids.push_back(th.sv_id);

        if (aa.anchorable && !aa.chrom_name.empty()) {
            if (hla.overlaps(aa.chrom_name, aa.L_pos, aa.R_pos)) r.flags |= 1;
            if (cent.overlaps(aa.chrom_name, aa.L_pos, aa.R_pos)) r.flags |= 2;
            if (segdup.overlaps(aa.chrom_name, aa.L_pos, aa.R_pos)) r.flags |= 4;
            if (simple.overlaps(aa.chrom_name, aa.L_pos, aa.R_pos)) r.flags |= 8;
            if (gap.overlaps(aa.chrom_name, aa.L_pos, aa.R_pos)) r.flags |= 16;
        }

        std::fwrite(&r, sizeof(r), 1, fbin);
        rows.push_back(r);

        // TSV row
        std::fprintf(ftsv,
            "%u\t%c\t%u\t%u\t%u\t%u\t%lu\t%lu\t%u\t%u\t%d\t%d\t%d\t",
            r.bubble_id, (char)r.type, r.n_nodes, r.n_normal, r.n_tumor, r.n_anchors,
            (unsigned long)r.n_cov_sum, (unsigned long)r.t_cov_sum, r.n_cov_max, r.t_cov_max,
            r.anchorable, r.mixed_strand, r.mixed_chrom);
        if (r.chrom_idx >= 0 && (size_t)r.chrom_idx < pm.chrom_names.size())
            std::fprintf(ftsv, "%s\t", pm.chrom_names[r.chrom_idx].c_str());
        else
            std::fprintf(ftsv, ".\t");
        std::fprintf(ftsv, "%u\t%u\t%u\t%u\t%d\t%s\t%c\t%d\t%u\t%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            r.L_pos, r.R_pos, r.L_strand, r.R_strand,
            r.truth_match,
            th.sv_id.empty() ? "." : th.sv_id.c_str(),
            r.truth_svtype ? (char)r.truth_svtype : '.',
            r.truth_svlen, r.truth_pos, r.truth_end,
            (r.flags & 1) ? 1 : 0, (r.flags & 2) ? 1 : 0, (r.flags & 4) ? 1 : 0,
            (r.flags & 8) ? 1 : 0, (r.flags & 16) ? 1 : 0,
            r.singleton, r.large);

        bubble_id++;
        if (bubble_id % 10000 == 0) std::fprintf(stderr, "[census]   %u / %zu\n", bubble_id, bubs.size());
    }
    std::fclose(fbin);
    std::fclose(ftsv);

    if (std::rename((out_bin + ".tmp").c_str(), out_bin.c_str()) != 0)
        die("rename census.bin");

    std::fprintf(stderr, "[census] writing bubble_summary.md ... ");
    write_summary(rows, truth_ids, pm, dir + "/bubble_summary.md");
    std::fprintf(stderr, "done\n");

    std::fprintf(stderr, "[census] outputs in %s/{census.bin,census.tsv,bubble_summary.md}\n", dir.c_str());
    return 0;
}
