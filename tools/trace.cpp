#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1   // POSIX getline for the VCF reader
#endif
// trace: map svp graph nodes <-> simulation truth (FASTA/VCF). Out-of-band analysis tool
// (NOT built by `make`; compile via tools/build_trace.sh). Reuses the production signature
// path verbatim: a node's stored window bytes ARE the bytes that were murmur-hashed, so
//   sig(node_i) = emit_sig(payload_i, seqlen_i, 0)
// and genome/VCF windows hashed the same way match node ids directly.
//
// graph.bin v2 layout (graph.cpp run_graph_reads / compact.cpp build_super_graph):
//   u64 K, u8 color[K], u32 ncov[K], u32 tcov[K], u64 E, IEdge edges[E],
//   u64 total_payload, u16 seqlen[K], u8 payload[total_payload]
//
// Commands:
//   probe <graph.bin> <fa> <chrom> <start0> <len>     -- G2 sig-bridge sanity
//   run   <graph.bin> <ref.fa> <out.fa> <out.vcf> <outdir>
//         -- whole-genome truth parse -> per-node {truth-source, genome pos}, then
//            marking-1 (variant->node via ref->out coord lift), marking-2 (read errors),
//            marking-3 (clean bubbles/cycles via unitig compaction). Resumable: the truth
//            parse checkpoints to <outdir>/annot.bin; later phases reuse it.
//   neigh <graph.bin> <outdir> <seeds.tsv> <radius>   -- per-seed subgraph JSON+GFA
#include "svp.h"            // extract decl, IEdge, SigCount
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

int LEVEL = 5;                          // runtime (SVP_TRACE_LEVEL); MUST match the graph's build --level
constexpr int FLANK = 0, MARGIN = 10;
static inline uint64_t emit_sig(const void* p, int n, int) { return core_sig((const char*)p, n); }   // XXH3 over the bare core (no flank); matches graph.bin
constexpr uint64_t CHUNK = 4u << 20, PAD = 8192;
constexpr uint8_t REF_FWD=1, REF_RC=2, OUT_FWD=4, OUT_RC=8, REF_ANY=REF_FWD|REF_RC, OUT_ANY=OUT_FWD|OUT_RC;
constexpr uint64_t LOC_NONE = ~0ull;
inline uint64_t pack_loc(uint32_t cid, uint64_t pos) { return ((uint64_t)cid << 40) | (pos & ((1ull<<40)-1)); }
inline uint32_t  loc_chrom(uint64_t v) { return (uint32_t)(v >> 40); }
inline uint64_t  loc_pos(uint64_t v)   { return v & ((1ull<<40)-1); }

uint8_t COMP[256];
void init_comp() { for (int i = 0; i < 256; i++) COMP[i] = 'N';
    COMP[(uint8_t)'A']='T'; COMP[(uint8_t)'T']='A'; COMP[(uint8_t)'C']='G'; COMP[(uint8_t)'G']='C'; }

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[FATAL] %s\n", m); std::abort(); }

int env_int(const char* k, int dflt) { const char* v = std::getenv(k); return v ? std::atoi(v) : dflt; }

// ---- mmap'd graph.bin ----
struct Graph {
    int fd = -1; const uint8_t* base = nullptr; size_t fsize = 0;
    uint64_t K = 0, E = 0, total_payload = 0;
    const uint8_t* color = nullptr; const uint32_t* ncov = nullptr; const uint32_t* tcov = nullptr;
    const IEdge* edges = nullptr; const uint16_t* seqlen = nullptr; const uint8_t* payload = nullptr;
    std::vector<uint64_t> seq_off;
};

const uint8_t* rdptr(const uint8_t* p, const uint8_t* end, size_t n, const char* what) {
    if (p + n > end) { std::fprintf(stderr, "[FATAL] graph.bin truncated reading %s\n", what); std::abort(); }
    return p;
}

void load_graph(const std::string& path, Graph& g) {
    g.fd = open(path.c_str(), O_RDONLY);
    if (g.fd < 0) die("open graph.bin");
    struct stat st; if (fstat(g.fd, &st) != 0) die("fstat graph.bin");
    g.fsize = st.st_size;
    void* m = mmap(nullptr, g.fsize, PROT_READ, MAP_PRIVATE, g.fd, 0);
    if (m == MAP_FAILED) die("mmap graph.bin");
    madvise(m, g.fsize, MADV_RANDOM);
    g.base = (const uint8_t*)m; const uint8_t* p = g.base, *end = g.base + g.fsize;
    std::memcpy(&g.K, rdptr(p, end, 8, "K"), 8); p += 8;
    g.color  = rdptr(p, end, g.K, "color");                    p += g.K;
    g.ncov   = (const uint32_t*)rdptr(p, end, g.K*4, "ncov");   p += g.K*4;
    g.tcov   = (const uint32_t*)rdptr(p, end, g.K*4, "tcov");   p += g.K*4;
    std::memcpy(&g.E, rdptr(p, end, 8, "E"), 8); p += 8;
    g.edges  = (const IEdge*)rdptr(p, end, g.E*sizeof(IEdge), "edges"); p += g.E*sizeof(IEdge);
    std::memcpy(&g.total_payload, rdptr(p, end, 8, "total_payload"), 8); p += 8;
    g.seqlen = (const uint16_t*)rdptr(p, end, g.K*2, "seqlen"); p += g.K*2;
    g.payload= rdptr(p, end, g.total_payload, "payload");
    g.seq_off.assign(g.K + 1, 0);
    for (uint64_t i = 0; i < g.K; i++) g.seq_off[i+1] = g.seq_off[i] + g.seqlen[i];
    std::fprintf(stderr, "[trace] graph.bin: K=%lu E=%lu payload=%.2fGB (%.2fGB mmap)\n",
                 (unsigned long)g.K, (unsigned long)g.E, g.total_payload/1e9, g.fsize/1e9);
}

inline uint64_t node_sig(const Graph& g, uint64_t i) {
    return emit_sig(g.payload + g.seq_off[i], (int)g.seqlen[i], 0);
}

void build_index(const Graph& g, std::unordered_map<uint64_t, uint32_t>& idx) {
    idx.reserve(g.K * 2);
    for (uint64_t i = 0; i < g.K; i++) idx.emplace(node_sig(g, i), (uint32_t)i);
    std::fprintf(stderr, "[trace] sig index: %zu entries\n", idx.size());
}

// ---- faidx random access ----
struct FaEntry { uint64_t len, offset; uint32_t linebases, linewidth; };
struct Faidx { int fd = -1; std::unordered_map<std::string, FaEntry> ix; std::vector<std::string> order; };

void load_faidx(const std::string& fa, Faidx& f) {
    f.fd = open(fa.c_str(), O_RDONLY);
    if (f.fd < 0) die("open fasta");
    FILE* fi = std::fopen((fa + ".fai").c_str(), "r");
    if (!fi) die("open .fai");
    char name[512]; uint64_t len, off; uint32_t lb, lw;
    while (std::fscanf(fi, "%511s\t%lu\t%lu\t%u\t%u\n", name, &len, &off, &lb, &lw) == 5) {
        f.ix[name] = {len, off, lb, lw}; f.order.push_back(name);
    }
    std::fclose(fi);
}

std::string fa_fetch(const Faidx& f, const std::string& chrom, uint64_t start0, uint64_t want) {
    auto it = f.ix.find(chrom);
    if (it == f.ix.end()) return "";
    const FaEntry& e = it->second;
    if (start0 >= e.len) return "";
    if (start0 + want > e.len) want = e.len - start0;
    const uint64_t fstart = e.offset + (start0 / e.linebases) * e.linewidth + (start0 % e.linebases);
    const uint64_t nl = e.linewidth - e.linebases;
    const uint64_t fbytes = want + (want / e.linebases + 2) * nl;
    std::string raw(fbytes, '\0');
    const ssize_t got = pread(f.fd, &raw[0], fbytes, (off_t)fstart);
    if (got <= 0) return "";
    raw.resize((size_t)got);
    std::string out; out.reserve(want);
    for (char c : raw) { if (c=='\n'||c=='\r') continue; out.push_back(c>='a'&&c<='z'?(char)(c-32):c); if (out.size()==want) break; }
    return out;
}

std::string revcomp(const std::string& s) {
    std::string r(s.size(), 'N');
    for (size_t j = 0; j < s.size(); j++) r[j] = (char)COMP[(uint8_t)s[s.size()-1-j]];
    return r;
}

bool has_edge(const Graph& g, uint32_t a, uint32_t b) {
    const IEdge key{a, b}; const IEdge *lo = g.edges, *hi = g.edges + g.E;
    while (lo < hi) { const IEdge* mid = lo + (hi - lo)/2; if (*mid < key) lo = mid+1; else hi = mid; }
    return lo != g.edges + g.E && lo->first == a && lo->second == b;
}

// =====================================================================================
// G2 probe
// =====================================================================================
int cmd_probe(int argc, char** argv) {
    if (argc < 6) die("probe <graph.bin> <fa> <chrom> <start0> <len>");
    Graph g; load_graph(argv[1], g);
    std::unordered_map<uint64_t, uint32_t> idx; build_index(g, idx);
    Faidx fa; load_faidx(argv[2], fa);
    std::string seq = fa_fetch(fa, argv[3], std::strtoull(argv[4],nullptr,10), std::strtoull(argv[5],nullptr,10));
    if (seq.empty()) die("fa_fetch empty");
    for (int pass = 0; pass < 2; pass++) {
        const std::string s = pass==0 ? seq : revcomp(seq);
        const int len = (int)s.size();
        int n; const uint64_t *st,*en; const uint64_t *lb; extract(s.data(), len, 0, LEVEL, &n,&st,&en,&lb); (void)lb;
        uint64_t cores=0,hit=0,adj_ok=0,adj_tot=0; bool have=false; uint32_t prev=0;
        for (int i=0;i<n;i++){ if (st[i]<(uint64_t)MARGIN||en[i]+MARGIN>(uint64_t)len) continue;
            const uint64_t sig=emit_sig(s.data()+st[i]-FLANK,(int)(en[i]-st[i]+2*FLANK),0); cores++;
            auto it=idx.find(sig); const bool f=it!=idx.end(); if(f)hit++;
            if(f&&have){adj_tot++; if(has_edge(g,prev,it->second))adj_ok++;} have=f; if(f)prev=it->second; }
        std::fprintf(stderr,"[probe] %s: cores=%lu found=%lu (%.1f%%) adj=%lu/%lu (%.1f%%)\n",
            pass==0?"fwd":"rc",(unsigned long)cores,(unsigned long)hit,cores?100.0*hit/cores:0,
            (unsigned long)adj_ok,(unsigned long)adj_tot,adj_tot?100.0*adj_ok/adj_tot:0);
    }
    return 0;
}

// =====================================================================================
// Per-node truth annotation (the resumable checkpoint)
// =====================================================================================
struct Annot {
    uint64_t K = 0;
    std::vector<uint8_t> truthsrc;     // REF_FWD|REF_RC|OUT_FWD|OUT_RC bits
    std::vector<uint64_t> ref_loc, out_loc;
    std::vector<std::string> chrom_names;   // ref chrom id -> name
};

bool load_annot(const std::string& path, Annot& a) {
    FILE* f = std::fopen(path.c_str(), "rb"); if (!f) return false;
    char magic[4]; if (std::fread(magic,1,4,f)!=4 || std::memcmp(magic,"TRAN",4)) { std::fclose(f); return false; }
    if (std::fread(&a.K,8,1,f)!=1) die("annot K");
    a.truthsrc.resize(a.K); a.ref_loc.resize(a.K); a.out_loc.resize(a.K);
    if (std::fread(a.truthsrc.data(),1,a.K,f)!=a.K) die("annot truthsrc");
    if (std::fread(a.ref_loc.data(),8,a.K,f)!=a.K) die("annot ref_loc");
    if (std::fread(a.out_loc.data(),8,a.K,f)!=a.K) die("annot out_loc");
    uint32_t nc; if (std::fread(&nc,4,1,f)!=1) die("annot nc");
    for (uint32_t i=0;i<nc;i++){ uint16_t L; if(std::fread(&L,2,1,f)!=1)die("annot name"); std::string s(L,0); if(std::fread(&s[0],1,L,f)!=L)die("annot name"); a.chrom_names.push_back(s);}
    std::fclose(f);
    std::fprintf(stderr, "[trace] loaded annot.bin (K=%lu, %zu chroms)\n", (unsigned long)a.K, a.chrom_names.size());
    return true;
}

void write_annot(const std::string& path, const Annot& a) {
    FILE* f = std::fopen((path+".tmp").c_str(), "wb"); if (!f) die("write annot");
    std::fwrite("TRAN",1,4,f); std::fwrite(&a.K,8,1,f);
    std::fwrite(a.truthsrc.data(),1,a.K,f); std::fwrite(a.ref_loc.data(),8,a.K,f); std::fwrite(a.out_loc.data(),8,a.K,f);
    uint32_t nc=(uint32_t)a.chrom_names.size(); std::fwrite(&nc,4,1,f);
    for (auto& s : a.chrom_names){ uint16_t L=(uint16_t)s.size(); std::fwrite(&L,2,1,f); std::fwrite(s.data(),1,L,f);}
    std::fclose(f);
    if (std::rename((path+".tmp").c_str(), path.c_str())) die("rename annot");
}

struct ParseJob { int g; std::string name; uint32_t cid; uint64_t clen, chunk_start; };

// process one chunk: fwd + rc passes; update truth annotation + recall counters (owned cores only)
void parse_chunk(const ParseJob& j, const Faidx& fa, const std::unordered_map<uint64_t,uint32_t>& idx,
                 std::vector<std::atomic<uint8_t>>& truthsrc, std::vector<std::atomic<uint64_t>>& loc,
                 std::atomic<uint64_t>* tot, std::atomic<uint64_t>* fnd, uint8_t bit_fwd, uint8_t bit_rc) {
    const uint64_t own_lo = j.chunk_start, own_hi = std::min(j.clen, j.chunk_start + CHUNK);
    const uint64_t cs_lo = j.chunk_start > PAD ? j.chunk_start - PAD : 0;
    const uint64_t cs_hi = std::min(j.clen, own_hi + PAD);
    std::string s = fa_fetch(fa, j.name, cs_lo, cs_hi - cs_lo);
    if ((int)s.size() < 2*MARGIN+4) return;
    uint64_t ltot=0, lfnd=0;
    for (int pass=0; pass<2; pass++) {
        const std::string d = pass==0 ? s : revcomp(s);
        const int len=(int)d.size(); const uint8_t bit = pass==0?bit_fwd:bit_rc;
        int n; const uint64_t *st,*en; const uint64_t *lb; extract(d.data(),len,0,LEVEL,&n,&st,&en,&lb); (void)lb;
        for (int i=0;i<n;i++){ if (st[i]<(uint64_t)MARGIN||en[i]+MARGIN>(uint64_t)len) continue;
            const uint64_t gstart = pass==0 ? cs_lo+st[i] : cs_lo+(len-en[i]);
            if (gstart<own_lo || gstart>=own_hi) continue;            // owned-once
            ltot++;
            const uint64_t sig=emit_sig(d.data()+st[i]-FLANK,(int)(en[i]-st[i]+2*FLANK),0);
            auto it=idx.find(sig); if (it==idx.end()) continue; lfnd++;
            truthsrc[it->second].fetch_or(bit, std::memory_order_relaxed);
            loc[it->second].store(pack_loc(j.cid, gstart), std::memory_order_relaxed);
        }
    }
    tot[j.cid].fetch_add(ltot, std::memory_order_relaxed);
    fnd[j.cid].fetch_add(lfnd, std::memory_order_relaxed);
}

// strip trailing _hapA/_hapB so out.fa contigs map onto the ref chrom id
std::string base_chrom(const std::string& n) {
    if (n.size()>5){ const std::string s=n.substr(n.size()-5); if (s=="_hapA"||s=="_hapB") return n.substr(0,n.size()-5);} return n;
}

void run_truth_parse(const Graph& g, const std::string& ref_fa, const std::string& out_fa, Annot& a) {
    std::unordered_map<uint64_t,uint32_t> idx; build_index(g, idx);
    Faidx fref, fout; load_faidx(ref_fa, fref); load_faidx(out_fa, fout);
    // chrom id table = ref faidx order; out contigs map via base_chrom
    std::unordered_map<std::string,uint32_t> cid;
    for (auto& nm : fref.order){ cid[nm]=(uint32_t)a.chrom_names.size(); a.chrom_names.push_back(nm);}
    const uint32_t NC=(uint32_t)a.chrom_names.size();
    std::vector<std::atomic<uint8_t>> truthsrc(g.K);  for (auto&x:truthsrc) x.store(0);
    std::vector<std::atomic<uint64_t>> refl(g.K), outl(g.K);
    for (uint64_t i=0;i<g.K;i++){ refl[i].store(LOC_NONE); outl[i].store(LOC_NONE);}
    std::vector<std::atomic<uint64_t>> tot(NC*2), fnd(NC*2); for (uint32_t i=0;i<NC*2;i++){tot[i].store(0);fnd[i].store(0);}
    // build jobs
    std::vector<ParseJob> jobs;
    auto add_jobs=[&](int gi, Faidx& fa){ for (auto& nm:fa.order){ const uint64_t L=fa.ix[nm].len;
        const std::string bc=base_chrom(nm); auto c=cid.find(bc); if(c==cid.end()) continue;
        for (uint64_t cs=0; cs<L; cs+=CHUNK) jobs.push_back({gi,nm,c->second,L,cs}); } };
    add_jobs(0,fref); add_jobs(1,fout);
    std::fprintf(stderr,"[trace] truth parse: %zu chunks over %u chroms\n", jobs.size(), NC);
    std::atomic<size_t> next{0}; const int NT=env_int("TRACE_THREADS",16);
    auto worker=[&](){ Faidx* fas[2]={&fref,&fout}; size_t k;
        while((k=next.fetch_add(1))<jobs.size()){ const ParseJob& j=jobs[k];
            parse_chunk(j,*fas[j.g],idx,truthsrc, j.g==0?refl:outl, j.g==0?tot.data():tot.data()+NC,
                        j.g==0?fnd.data():fnd.data()+NC, j.g==0?REF_FWD:OUT_FWD, j.g==0?REF_RC:OUT_RC);
            if ((k&4095)==0) std::fprintf(stderr,"\r[trace]   chunk %zu/%zu",k,jobs.size()); } };
    std::vector<std::thread> th; for (int i=0;i<NT;i++) th.emplace_back(worker); for (auto&t:th) t.join();
    std::fprintf(stderr,"\r[trace]   chunks done                          \n");
    a.K=g.K; a.truthsrc.resize(g.K); a.ref_loc.resize(g.K); a.out_loc.resize(g.K);
    for (uint64_t i=0;i<g.K;i++){ a.truthsrc[i]=truthsrc[i].load(); a.ref_loc[i]=refl[i].load(); a.out_loc[i]=outl[i].load(); }
    // recall table to stderr (also returned via summary later)
    uint64_t rt=0,rf=0,ot=0,of=0; for (uint32_t c=0;c<NC;c++){rt+=tot[c];rf+=fnd[c];ot+=tot[NC+c];of+=fnd[NC+c];}
    std::fprintf(stderr,"[trace] recall  ref: %lu/%lu (%.2f%%)  out: %lu/%lu (%.2f%%)\n",
        (unsigned long)rf,(unsigned long)rt, rt?100.0*rf/rt:0,(unsigned long)of,(unsigned long)ot, ot?100.0*of/ot:0);
}

}  // namespace

// later phases (marking-1/2/3, neigh, diag) appended in part 2
int trace_run(int argc, char** argv);
int trace_neigh(int argc, char** argv);
int trace_diag(int argc, char** argv);
int trace_probeseq(int argc, char** argv);
int trace_shapes(int argc, char** argv);

int main(int argc, char** argv) {
    init_comp();
    LEVEL = env_int("SVP_TRACE_LEVEL", 5);
    if (LEVEL != 5) std::fprintf(stderr, "[trace] LCP level = %d (SVP_TRACE_LEVEL)\n", LEVEL);
    if (argc < 2) { std::fprintf(stderr, "usage: trace <probe|run|neigh> ...\n"); return 1; }
    const std::string cmd = argv[1];
    if (cmd == "probe") return cmd_probe(argc-1, argv+1);
    if (cmd == "run")   return trace_run(argc-1, argv+1);
    if (cmd == "neigh") return trace_neigh(argc-1, argv+1);
    if (cmd == "diag")  return trace_diag(argc-1, argv+1);
    if (cmd == "probeseq") return trace_probeseq(argc-1, argv+1);
    if (cmd == "shapes")   return trace_shapes(argc-1, argv+1);
    std::fprintf(stderr, "[FATAL] unknown command %s\n", cmd.c_str());
    return 1;
}

// =====================================================================================
// PART 2: marking-1 (variant lift), marking-2 (errors), marking-3 (bubbles), neigh
// =====================================================================================
namespace {

void split_tab(const std::string& s, std::vector<std::string>& out) {
    out.clear(); size_t p = 0;
    while (true) { size_t q = s.find('\t', p); out.push_back(s.substr(p, q==std::string::npos?std::string::npos:q-p)); if (q==std::string::npos) break; p=q+1; }
}
std::string info_get(const std::string& info, const char* key) {
    const std::string k = std::string(key) + "="; size_t p = info.find(k);
    if (p==std::string::npos) return ""; p += k.size(); size_t e = info.find(';', p);
    return info.substr(p, e==std::string::npos?std::string::npos:e-p);
}

struct Variant { uint32_t cid; uint64_t pos0; int64_t dl; uint64_t svlen; uint64_t out_pos;
                 std::string type, svid;
                 uint32_t n_nodes=0,n_c1=0,n_c2=0,n_c3=0,best_node=UINT32_MAX,best_ncov=0,best_tcov=0; uint8_t best_color=0; int64_t pos_err=-1; };

void classify_variant(const std::string& info, Variant& v) {
    const std::string svt = info_get(info, "SVTYPE");
    v.svlen = std::strtoull(info_get(info, "SVLEN").c_str(), nullptr, 10);
    if      (svt=="DEL") { v.type="DEL"; v.dl=-(int64_t)v.svlen; }
    else if (svt=="INS") { v.type="INS"; v.dl=(int64_t)v.svlen; }
    else if (svt=="DUP") { v.type="DUP"; v.dl=(int64_t)v.svlen; }
    else if (svt=="INV") { v.type="INV"; v.dl=0; }
    else                 { v.type="SNP"; v.svlen=0; v.dl=0; }
}

// parse VCF -> per-chrom variants; compute the ref->out coordinate lift (out_pos).
bool parse_vcf(const std::string& path, const std::unordered_map<std::string,uint32_t>& cid, std::vector<Variant>& vs) {
    FILE* f = std::fopen(path.c_str(), "r"); if (!f) return false;
    char* line=nullptr; size_t cap=0; ssize_t n; std::vector<std::string> fld;
    while ((n=getline(&line,&cap,f))>0) {
        if (line[0]=='#') continue;
        std::string ln(line, line[n-1]=='\n'?n-1:n); split_tab(ln, fld);
        if (fld.size()<8) continue;
        auto c = cid.find(fld[0]); if (c==cid.end()) continue;
        Variant v; v.cid=c->second; v.pos0=std::strtoull(fld[1].c_str(),nullptr,10)-1; v.svid=fld[2];
        classify_variant(fld[7], v); vs.push_back(std::move(v));
    }
    std::free(line); std::fclose(f);
    // lift per chrom: stable order by pos, cumulative dl of preceding variants
    std::stable_sort(vs.begin(), vs.end(), [](const Variant&a,const Variant&b){ return a.cid!=b.cid?a.cid<b.cid:a.pos0<b.pos0; });
    uint32_t cur=UINT32_MAX; int64_t off=0;
    for (auto& v : vs){ if (v.cid!=cur){cur=v.cid;off=0;} v.out_pos=(uint64_t)((int64_t)v.pos0+off); off+=v.dl; }
    std::fprintf(stderr,"[trace] VCF: %zu variants\n", vs.size());
    return true;
}

// marking-1: attribute OUT-only (somatic) nodes to variants in REF coordinates. A somatic
// junction node neighbours ref-backed anchors, so it inherits a ref locus from a graph
// neighbour (drift-free) -- far more reliable than a 126k-variant cumulative coordinate lift.
void marking1(const Graph& g, const Annot& a, std::vector<Variant>& vs, const std::string& outdir) {
    auto is_som = [&](uint64_t i){ const uint8_t t=a.truthsrc[i]; return (t&OUT_ANY) && !(t&REF_ANY); };
    std::vector<uint64_t> som_ref(g.K, LOC_NONE);                   // somatic node -> ref locus (from a neighbour)
    for (uint64_t e=0;e<g.E;e++){ const uint32_t u=g.edges[e].first, w=g.edges[e].second;
        if (is_som(u) && a.ref_loc[w]!=LOC_NONE && som_ref[u]==LOC_NONE) som_ref[u]=a.ref_loc[w];
        if (is_som(w) && a.ref_loc[u]!=LOC_NONE && som_ref[w]==LOC_NONE) som_ref[w]=a.ref_loc[u]; }
    std::vector<std::pair<uint64_t,uint32_t>> som;                  // key=(cid<<40)|ref_pos, node
    uint64_t tot_som=0, located=0;
    for (uint64_t i=0;i<g.K;i++){ if(!is_som(i)) continue; tot_som++;
        if (som_ref[i]!=LOC_NONE){ located++; som.emplace_back(som_ref[i], (uint32_t)i); } }
    std::sort(som.begin(), som.end());
    std::fprintf(stderr,"[trace] marking-1: %lu somatic nodes, %lu located via ref-neighbour (%.1f%%)\n",
        (unsigned long)tot_som,(unsigned long)located, tot_som?100.0*located/tot_som:0);
    const uint64_t W=64;
    for (auto& v : vs){
        const bool spans = (v.type=="DEL"||v.type=="INV"||v.type=="DUP");
        const uint64_t lo = v.pos0>W? v.pos0-W : 0;
        const uint64_t hi = v.pos0 + (spans?v.svlen:0) + W;          // [pos,end] for span types; pos±W for SNP/INS
        const uint64_t klo=pack_loc(v.cid,lo), khi=pack_loc(v.cid,hi);
        auto it=std::lower_bound(som.begin(),som.end(),std::make_pair(klo,(uint32_t)0));
        for (; it!=som.end() && it->first<=khi; ++it){ const uint32_t nd=it->second; v.n_nodes++;
            const uint8_t col=g.color[nd]; if(col==1)v.n_c1++; else if(col==2)v.n_c2++; else v.n_c3++;
            const bool better = (v.best_node==UINT32_MAX) || (col==2 && v.best_color!=2) || (col==v.best_color && g.tcov[nd]>v.best_tcov);
            if (better){ v.best_node=nd; v.best_color=col; v.best_ncov=g.ncov[nd]; v.best_tcov=g.tcov[nd]; } }
    }
    // localization precision: bp distance from the true position to the NEAREST attributed somatic node
    for (auto& v : vs){ const uint64_t key=pack_loc(v.cid,v.pos0);
        auto it=std::lower_bound(som.begin(),som.end(),std::make_pair(key,(uint32_t)0)); int64_t best=-1;
        if(it!=som.end()&&loc_chrom(it->first)==v.cid){ int64_t d=(int64_t)loc_pos(it->first)-(int64_t)v.pos0; best=d<0?-d:d; }
        if(it!=som.begin()){ auto p=std::prev(it); if(loc_chrom(p->first)==v.cid){ int64_t d=(int64_t)v.pos0-(int64_t)loc_pos(p->first); if(d<0)d=-d; if(best<0||d<best)best=d; } }
        v.pos_err=best; }
    FILE* f=std::fopen((outdir+"/variants.tsv").c_str(),"w");
    std::fprintf(f,"svid\tchrom\tpos1\ttype\tsvlen\tout_pos\tn_nodes\tn_c1\tn_c2\tn_c3\tbest_node\tbest_color\tbest_ncov\tbest_tcov\tfound\tpos_err\n");
    uint64_t found=0; for (auto& v: vs){ if(v.n_nodes>0)found++;
        std::fprintf(f,"%s\t%s\t%lu\t%s\t%lu\t%lu\t%u\t%u\t%u\t%u\t%d\t%u\t%u\t%u\t%d\t%ld\n",
            v.svid.c_str(), a.chrom_names[v.cid].c_str(), (unsigned long)v.pos0+1, v.type.c_str(),
            (unsigned long)v.svlen,(unsigned long)v.out_pos, v.n_nodes,v.n_c1,v.n_c2,v.n_c3,
            v.best_node==UINT32_MAX?-1:(int)v.best_node, v.best_color, v.best_ncov, v.best_tcov, v.n_nodes>0?1:0, (long)v.pos_err); }
    std::fclose(f);
    std::fprintf(stderr,"[trace] marking-1: %lu/%zu variants have >=1 somatic node in graph (%.1f%%)\n",
        (unsigned long)found, vs.size(), vs.size()?100.0*found/vs.size():0);
}

// ---- union-find ----
uint32_t uf_find(std::vector<uint32_t>& p, uint32_t x){ while(p[x]!=x){p[x]=p[p[x]];x=p[x];} return x; }
void uf_union(std::vector<uint32_t>& p, uint32_t a, uint32_t b){ a=uf_find(p,a); b=uf_find(p,b); if(a!=b)p[a<b?b:a]=a<b?a:b; }

// marking-2: read-error nodes (truthsrc==0); cluster via error-error edges; locate via truth neighbor.
void marking2(const Graph& g, const Annot& a, const std::string& outdir) {
    uint64_t nerr=0; for (uint64_t i=0;i<g.K;i++) if(a.truthsrc[i]==0) nerr++;
    std::vector<uint32_t> par(g.K); for (uint64_t i=0;i<g.K;i++) par[i]=(uint32_t)i;
    for (uint64_t e=0;e<g.E;e++){ const uint32_t u=g.edges[e].first,w=g.edges[e].second;
        if (a.truthsrc[u]==0 && a.truthsrc[w]==0) uf_union(par,u,w); }
    // per-root aggregate
    struct Cl { uint32_t size=0,c1=0,c2=0,c3=0,seed=UINT32_MAX,sncov=0,stcov=0; uint64_t scov=0,nearloc=LOC_NONE; };
    std::unordered_map<uint32_t,Cl> cls; cls.reserve(nerr/2+1);
    for (uint64_t i=0;i<g.K;i++){ if(a.truthsrc[i]!=0) continue; const uint32_t r=uf_find(par,(uint32_t)i);
        Cl& c=cls[r]; c.size++; const uint8_t col=g.color[i]; if(col==1)c.c1++;else if(col==2)c.c2++;else c.c3++;
        const uint64_t cov=(uint64_t)g.ncov[i]+g.tcov[i]; if(cov>c.scov||c.seed==UINT32_MAX){c.scov=cov;c.seed=(uint32_t)i;c.sncov=g.ncov[i];c.stcov=g.tcov[i];} }
    // locate: an edge error<->truth gives the cluster a genome anchor
    for (uint64_t e=0;e<g.E;e++){ const uint32_t u=g.edges[e].first,w=g.edges[e].second;
        const bool eu=a.truthsrc[u]==0, ew=a.truthsrc[w]==0;
        if (eu==ew) continue; const uint32_t er=eu?u:w, tr=eu?w:u;
        if (a.ref_loc[tr]==LOC_NONE) continue; Cl& c=cls[uf_find(par,er)]; if(c.nearloc==LOC_NONE)c.nearloc=a.ref_loc[tr]; }
    std::vector<std::pair<uint32_t,Cl>> v(cls.begin(),cls.end());
    std::sort(v.begin(),v.end(),[](const auto&x,const auto&y){return x.second.size>y.second.size;});
    FILE* f=std::fopen((outdir+"/error_clusters.tsv").c_str(),"w");
    std::fprintf(f,"cluster\tsize\tn_c1\tn_c2\tn_c3\tseed_node\tseed_ncov\tseed_tcov\tnear_chrom\tnear_pos\n");
    const size_t CAP=300000; size_t emitted=0;
    for (auto& kv: v){ if(emitted++>=CAP) break; const Cl& c=kv.second;
        const char* nch = c.nearloc==LOC_NONE?"NA":a.chrom_names[loc_chrom(c.nearloc)].c_str();
        std::fprintf(f,"%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%ld\n", kv.first,c.size,c.c1,c.c2,c.c3,
            c.seed,c.sncov,c.stcov,nch, c.nearloc==LOC_NONE?-1L:(long)loc_pos(c.nearloc)); }
    std::fclose(f);
    std::fprintf(stderr,"[trace] marking-2: %lu error nodes in %zu clusters (%.2f%% of K); emitted top %zu\n",
        (unsigned long)nerr, cls.size(), g.K?100.0*nerr/g.K:0, emitted-1<CAP?emitted-1:CAP);
}

// ---- marking-3: unitig compaction (topology only) + bubble/cycle enumeration ----
struct SuperMeta { uint8_t color, alltruth; uint32_t members, bp, rep, ncov, tcov; uint64_t refloc, boff; };

// KMP longest k with A[la-k:]==B[:k]; spliced super-node bp = sum seqlen minus these overlaps.
static uint32_t sfx_pfx(const uint8_t* A, uint32_t la, const uint8_t* B, uint32_t lb){
    if(!la||!lb)return 0; std::vector<uint32_t> f(lb,0);
    for(uint32_t i=1,j=0;i<lb;i++){ while(j&&B[i]!=B[j])j=f[j-1]; if(B[i]==B[j])j++; f[i]=j; }
    uint32_t j=0; for(uint32_t i=0;i<la;i++){ while(j&&A[i]!=B[j])j=f[j-1]; if(A[i]==B[j])j++; if(j==lb){ if(i+1==la)return lb; j=f[j-1]; } }
    return j;
}

void compact(const Graph& g, const Annot& a, std::vector<uint32_t>& nts,
             std::vector<SuperMeta>& sm, std::vector<IEdge>& sedges, std::vector<uint8_t>* sblocks=nullptr) {
    const uint64_t K=g.K;
    std::vector<uint32_t> ind(K,0),outd(K,0),succ(K,UINT32_MAX),pred(K,UINT32_MAX);
    for (uint64_t e=0;e<g.E;e++){ outd[g.edges[e].first]++; ind[g.edges[e].second]++; succ[g.edges[e].first]=g.edges[e].second; pred[g.edges[e].second]=g.edges[e].first; }
    nts.assign(K,UINT32_MAX); std::vector<char> vis(K,0);
    auto emit=[&](const std::vector<uint32_t>& chain){ const uint32_t sid=(uint32_t)sm.size();
        SuperMeta m; m.color=g.color[chain[0]]; m.alltruth=1; m.members=(uint32_t)chain.size(); m.rep=chain[0];
        m.ncov=0; m.tcov=0; m.refloc=LOC_NONE;
        for (uint32_t nd: chain){ nts[nd]=sid; if(a.truthsrc[nd]==0)m.alltruth=0;
            if(g.ncov[nd]>m.ncov)m.ncov=g.ncov[nd]; if(g.tcov[nd]>m.tcov)m.tcov=g.tcov[nd];
            if(m.refloc==LOC_NONE&&a.ref_loc[nd]!=LOC_NONE)m.refloc=a.ref_loc[nd]; }
        if(sblocks){ m.boff=sblocks->size();                     // build the spliced super-node sequence
            sblocks->insert(sblocks->end(), g.payload+g.seq_off[chain[0]], g.payload+g.seq_off[chain[0]]+g.seqlen[chain[0]]);
            for(size_t i=1;i<chain.size();i++){ const uint32_t cl=chain[i], lt=g.seqlen[chain[i-1]];
                const uint32_t k=sfx_pfx(sblocks->data()+(sblocks->size()-lt),lt,g.payload+g.seq_off[cl],g.seqlen[cl]);
                sblocks->insert(sblocks->end(), g.payload+g.seq_off[cl]+k, g.payload+g.seq_off[cl]+g.seqlen[cl]); }
            m.bp=(uint32_t)(sblocks->size()-m.boff);
        } else { m.boff=0; m.bp=g.seqlen[chain[0]];               // length only (marking-3 doesn't need the bytes)
            for(size_t i=1;i<chain.size();i++){ const uint32_t pl=chain[i-1],cl=chain[i];
                const uint32_t k=sfx_pfx(g.payload+g.seq_off[pl],g.seqlen[pl],g.payload+g.seq_off[cl],g.seqlen[cl]); m.bp+=g.seqlen[cl]-k; } }
        sm.push_back(m); };
    std::vector<uint32_t> chain;
    for (uint32_t u=0;u<(uint32_t)K;u++){ if(vis[u]) continue;
        const bool left_end=(ind[u]!=1)||(outd[pred[u]]!=1)||(g.color[pred[u]]!=g.color[u]);
        if(!left_end) continue;
        chain.assign(1,u); vis[u]=1; uint32_t cur=u;
        while(outd[cur]==1){ uint32_t nx=succ[cur]; if(ind[nx]!=1||g.color[nx]!=g.color[cur]||vis[nx])break; chain.push_back(nx); vis[nx]=1; cur=nx; }
        emit(chain); }
    uint64_t cyc=0; for (uint32_t u=0;u<(uint32_t)K;u++){ if(vis[u])continue; chain.assign(1,u); vis[u]=1; cyc++; emit(chain); }
    for (uint64_t e=0;e<g.E;e++){ const uint32_t su=nts[g.edges[e].first],sv=nts[g.edges[e].second]; if(su!=sv)sedges.push_back({su,sv}); }
    std::sort(sedges.begin(),sedges.end()); sedges.erase(std::unique(sedges.begin(),sedges.end()),sedges.end());
    std::fprintf(stderr,"[trace] compact: K=%lu -> S=%zu (cycle-singletons=%lu) super_edges=%zu\n",
        (unsigned long)K, sm.size(), (unsigned long)cyc, sedges.size());
}

struct CSR { std::vector<uint64_t> off; std::vector<uint32_t> nbr; };
void build_csr(const std::vector<IEdge>& se, uint64_t S, CSR& c){
    c.off.assign(S+1,0); for (auto&e:se) c.off[e.first+1]++;
    for (uint64_t i=0;i<S;i++) c.off[i+1]+=c.off[i];
    c.nbr.resize(se.size()); std::vector<uint64_t> w(c.off.begin(),c.off.begin()+S);
    for (auto&e:se) c.nbr[w[e.first]++]=e.second;
}

// bfs over supers from s1 to s2, interior colour-masked; returns interior super ids in `arm`.
bool bfs_arm(uint32_t s1, uint32_t s2, const CSR& og, const std::vector<uint8_t>& col,
             uint8_t start_mask, uint8_t interior_mask, uint32_t maxhops,
             std::vector<uint32_t>& seen, std::vector<uint32_t>& par, std::vector<uint32_t>& dist, uint32_t gen,
             std::vector<uint32_t>& arm){
    arm.clear(); seen[s1]=gen; par[s1]=s1; dist[s1]=0;
    std::vector<uint32_t> q; q.push_back(s1); size_t h=0;
    while(h<q.size()){ const uint32_t u=q[h++]; if(dist[u]>=maxhops)continue;
        const uint8_t mask=(u==s1)?start_mask:interior_mask;
        for (uint64_t i=og.off[u];i<og.off[u+1];i++){ const uint32_t v=og.nbr[i];
            if(v==s2&&col[v]==3){ if(seen[v]!=gen||u<par[v]){seen[v]=gen;par[v]=u;dist[v]=dist[u]+1;} continue; }
            if(!(mask&(1u<<(col[v]-1))))continue; if(seen[v]==gen)continue; seen[v]=gen;par[v]=u;dist[v]=dist[u]+1; q.push_back(v); } }
    if(seen[s2]!=gen) return false;
    for(uint32_t cur=par[s2]; cur!=s1; cur=par[cur]){ if(par[cur]==cur)return false; arm.push_back(cur);} return true;
}

// dfs from branch start x (colour!=3 interior); record c3 closures + cycle hit.
void dfs_branch(uint32_t s1, uint32_t x, const std::vector<uint8_t>& col, const CSR& og,
                std::vector<uint32_t>& vis, uint32_t gen, std::vector<uint32_t>& closed, bool& cyc){
    closed.clear(); cyc=false; std::vector<uint32_t> stk; stk.push_back(x); vis[x]=gen; uint32_t nodes=0;
    while(!stk.empty()&&nodes<500){ const uint32_t n=stk.back(); stk.pop_back(); nodes++;
        for (uint64_t i=og.off[n];i<og.off[n+1];i++){ const uint32_t nx=og.nbr[i];
            if(nx==s1){cyc=true;continue;} if(col[nx]==3){closed.push_back(nx);continue;}
            if(vis[nx]==gen){cyc=true;continue;} vis[nx]=gen; stk.push_back(nx); } }
}

struct Bubble { uint32_t s1,s2; uint8_t bucket; bool cycle; };

void enumerate(const std::vector<uint8_t>& col, const CSR& og, std::vector<Bubble>& out,
               uint64_t& direct, uint64_t& multi){
    out.clear(); direct=multi=0; const uint64_t S=col.size();
    std::vector<uint32_t> vis(S,0),closed; uint32_t gen=0;
    struct PI{uint32_t c1=0,c2=0; bool direct=false,cyc=false;}; std::unordered_map<uint64_t,PI> pi;
    for (uint32_t s1=0;s1<(uint32_t)S;s1++){ if(col[s1]!=3)continue; pi.clear();
        for (uint64_t i=og.off[s1];i<og.off[s1+1];i++){ const uint32_t x=og.nbr[i]; if(col[x]==3)pi[((uint64_t)s1<<32)|x].direct=true; }
        for (uint64_t i=og.off[s1];i<og.off[s1+1];i++){ const uint32_t x=og.nbr[i]; if(col[x]==3)continue;
            const uint8_t bc=col[x]; bool cy=false; dfs_branch(s1,x,col,og,vis,++gen,closed,cy);
            for (uint32_t s2: closed){ PI& p=pi[((uint64_t)s1<<32)|s2]; if(bc==1)p.c1++; else p.c2++; if(cy)p.cyc=true; } }
        for (auto& kv: pi){ const uint32_t s2=(uint32_t)kv.first; const PI& v=kv.second;
            const bool h1=v.c1>0,h2=v.c2>0; if(!h1&&!h2){direct++;continue;}
            if(v.cyc){ out.push_back({s1,s2,0,true}); continue; }
            if(v.c1>1||v.c2>1){multi++;continue;}
            uint8_t bucket = (h1&&h2)?0:(h2?1:2); out.push_back({s1,s2,bucket,false}); } }
}

void marking3(const Graph& g, const Annot& a, std::vector<Variant>& vs, const std::string& outdir){
    std::vector<uint32_t> nts; std::vector<SuperMeta> sm; std::vector<IEdge> sedges;
    compact(g,a,nts,sm,sedges);
    const uint64_t S=sm.size(); std::vector<uint8_t> col(S); for (uint64_t i=0;i<S;i++) col[i]=sm[i].color;
    CSR og; build_csr(sedges,S,og);
    std::vector<Bubble> bubs; uint64_t direct=0,multi=0; enumerate(col,og,bubs,direct,multi);
    uint64_t ncyc=0,nbub=0; for (auto&b:bubs){ if(b.cycle)ncyc++; else nbub++; }
    std::fprintf(stderr,"[trace] marking-3: bubbles=%lu cycles=%lu (direct=%lu multi=%lu) [expect ~494923 / 120]\n",
        (unsigned long)nbub,(unsigned long)ncyc,(unsigned long)direct,(unsigned long)multi);
    // per-chrom variant index for matching
    std::vector<std::pair<uint64_t,uint32_t>> vix;   // (cid<<40|pos0, variant idx)
    for (uint32_t i=0;i<vs.size();i++) vix.emplace_back(pack_loc(vs[i].cid,vs[i].pos0), i);
    std::sort(vix.begin(),vix.end());
    auto match=[&](uint64_t refloc)->const char*{ if(refloc==LOC_NONE)return "NA";
        const uint64_t cid=loc_chrom(refloc),pos=loc_pos(refloc);
        const uint64_t lo=pack_loc((uint32_t)cid,pos>200?pos-200:0), hi=pack_loc((uint32_t)cid,pos+200);
        auto it=std::lower_bound(vix.begin(),vix.end(),std::make_pair(lo,(uint32_t)0));
        return (it!=vix.end()&&it->first<=hi)? vs[it->second].svid.c_str() : "NA"; };
    std::vector<uint32_t> seen(S,0),par(S,0),dist(S,0); uint32_t gen=0; std::vector<uint32_t> arm;
    FILE* f=std::fopen((outdir+"/bubbles.tsv").c_str(),"w");
    std::fprintf(f,"id\tkind\tbucket\ts1_node\ts2_node\tchrom\tpos\tarm_members\tclean\tncov\ttcov\tmatched_svid\n");
    const char* bnames[3]={"std_both","T_only","N_only"};
    for (size_t i=0;i<bubs.size();i++){ const Bubble& b=bubs[i];
        bool clean = sm[b.s1].alltruth && sm[b.s2].alltruth; uint32_t armmem=0;
        for (int pass=0;pass<2;pass++){ const uint8_t smk=pass==0?0x1:0x2;
            if (bfs_arm(b.s1,b.s2,og,col,smk,0x3,500,seen,par,dist,++gen,arm))
                for (uint32_t s: arm){ armmem+=sm[s].members; if(!sm[s].alltruth)clean=false; } }
        const uint64_t rl=sm[b.s1].refloc;
        std::fprintf(f,"%zu\t%s\t%s\t%u\t%u\t%s\t%ld\t%u\t%d\t%u\t%u\t%s\n", i, b.cycle?"cycle":"bubble",
            bnames[b.bucket], sm[b.s1].rep, sm[b.s2].rep, rl==LOC_NONE?"NA":a.chrom_names[loc_chrom(rl)].c_str(),
            rl==LOC_NONE?-1L:(long)loc_pos(rl), armmem, clean?1:0, sm[b.s1].ncov, sm[b.s1].tcov, match(rl)); }
    std::fclose(f);
}

}  // namespace

int trace_run(int argc, char** argv) {
    if (argc < 6) die("run <graph.bin> <ref.fa> <out.fa> <out.vcf> <outdir>");
    const std::string gp=argv[1], ref=argv[2], out=argv[3], vcf=argv[4], dir=argv[5];
    mkdir(dir.c_str(),0755);
    Graph g; load_graph(gp, g);
    Annot a;
    if (!load_annot(dir+"/annot.bin", a)) {
        run_truth_parse(g, ref, out, a);
        write_annot(dir+"/annot.bin", a);
        std::fprintf(stderr,"[trace] wrote %s/annot.bin\n", dir.c_str());
    }
    // class histogram
    uint64_t err=0,refonly=0,outonly=0,both=0; for (uint64_t i=0;i<g.K;i++){ const uint8_t t=a.truthsrc[i];
        if(t==0)err++; else if((t&REF_ANY)&&(t&OUT_ANY))both++; else if(t&REF_ANY)refonly++; else outonly++; }
    std::fprintf(stderr,"[trace] nodes: total=%lu both=%lu ref-only=%lu out-only(somatic)=%lu read-error=%lu\n",
        (unsigned long)g.K,(unsigned long)both,(unsigned long)refonly,(unsigned long)outonly,(unsigned long)err);
    std::unordered_map<std::string,uint32_t> cid; for (uint32_t i=0;i<a.chrom_names.size();i++) cid[a.chrom_names[i]]=i;
    std::vector<Variant> vs; parse_vcf(vcf, cid, vs);
    marking1(g,a,vs,dir);
    marking2(g,a,dir);
    marking3(g,a,vs,dir);
    // summary.json
    FILE* s=std::fopen((dir+"/summary.json").c_str(),"w");
    std::fprintf(s,"{\"K\":%lu,\"E\":%lu,\"both\":%lu,\"ref_only\":%lu,\"somatic\":%lu,\"read_error\":%lu,\"variants\":%zu}\n",
        (unsigned long)g.K,(unsigned long)g.E,(unsigned long)both,(unsigned long)refonly,(unsigned long)outonly,(unsigned long)err, vs.size());
    std::fclose(s);
    std::fprintf(stderr,"[trace] run done -> %s\n", dir.c_str());
    return 0;
}

// =====================================================================================
// neigh: per-seed subgraph (BFS radius) -> JSON (viewer) + GFA (Bandage)
// =====================================================================================
namespace {
const char* node_class(uint8_t t){ if(t==0)return "error"; const bool r=t&REF_ANY,o=t&OUT_ANY;
    return (r&&o)?"shared":(r?"ref":"somatic"); }
// truth-class palette for Bandage --color CSV + the in-browser viewer
const char* class_hex(const char* cls){
    if(!std::strcmp(cls,"error"))return "#111111"; if(!std::strcmp(cls,"somatic"))return "#e6194b";
    if(!std::strcmp(cls,"shared"))return "#4363d8"; return "#9aa0a6"; /*ref*/ }

void out_edges_range(const Graph& g, uint32_t u, uint64_t& lo, uint64_t& hi){
    const IEdge key{u,0}; const IEdge *L=g.edges,*H=g.edges+g.E;
    while(L<H){const IEdge* m=L+(H-L)/2; if(*m<key)L=m+1; else H=m;} lo=(uint64_t)(L-g.edges);
    const IEdge key2{u,UINT32_MAX}; L=g.edges+lo; H=g.edges+g.E;
    while(L<H){const IEdge* m=L+(H-L)/2; if(*m<key2)L=m+1; else H=m;} hi=(uint64_t)(L-g.edges);
}

void emit_seed(const Graph& g, const Annot& a, const CSR& und, uint32_t seed, const std::string& label,
               uint32_t radius, uint32_t maxn, const std::string& dir){
    std::vector<uint32_t> q; q.push_back(seed); size_t h=0;
    std::unordered_map<uint32_t,uint32_t> depth; depth[seed]=0;
    while(h<q.size() && q.size()<maxn){ const uint32_t u=q[h++]; const uint32_t d=depth[u]; if(d>=radius)continue;
        for (uint64_t i=und.off[u];i<und.off[u+1];i++){ const uint32_t v=und.nbr[i];
            if(depth.count(v))continue; depth[v]=d+1; q.push_back(v); if(q.size()>=maxn)break; } }
    std::unordered_map<uint32_t,uint32_t> inset; for (uint32_t i=0;i<q.size();i++) inset[q[i]]=i;
    FILE* j=std::fopen((dir+"/"+label+".json").c_str(),"w");
    std::fprintf(j,"{\"seed\":%u,\"nodes\":[",seed);
    for (size_t i=0;i<q.size();i++){ const uint32_t n=q[i]; const uint64_t rl=a.ref_loc[n]!=LOC_NONE?a.ref_loc[n]:a.out_loc[n];
        std::string seq((const char*)g.payload+g.seq_off[n], g.seqlen[n]);
        std::fprintf(j,"%s{\"id\":%u,\"color\":%u,\"cls\":\"%s\",\"ncov\":%u,\"tcov\":%u,\"len\":%u,\"chrom\":\"%s\",\"pos\":%ld,\"seq\":\"%s\"}",
            i?",":"", n, g.color[n], node_class(a.truthsrc[n]), g.ncov[n], g.tcov[n], g.seqlen[n],
            rl==LOC_NONE?"NA":a.chrom_names[loc_chrom(rl)].c_str(), rl==LOC_NONE?-1L:(long)loc_pos(rl), seq.c_str()); }
    std::fprintf(j,"],\"edges\":[");
    FILE* gf=std::fopen((dir+"/"+label+".gfa").c_str(),"w"); std::fprintf(gf,"H\tVN:Z:1.0\n");
    FILE* cf=std::fopen((dir+"/"+label+".colors.csv").c_str(),"w"); std::fprintf(cf,"Node,Colour\n");
    for (size_t i=0;i<q.size();i++){ const uint32_t n=q[i]; std::string seq((const char*)g.payload+g.seq_off[n],g.seqlen[n]);
        std::fprintf(gf,"S\t%u\t%s\tKC:i:%u\n", n, seq.c_str(), g.ncov[n]+g.tcov[n]);
        std::fprintf(cf,"%u,%s\n", n, class_hex(node_class(a.truthsrc[n]))); }
    std::fclose(cf);
    bool first=true;
    for (size_t i=0;i<q.size();i++){ const uint32_t u=q[i]; uint64_t lo,hi; out_edges_range(g,u,lo,hi);
        for (uint64_t e=lo;e<hi;e++){ const uint32_t v=g.edges[e].second; if(!inset.count(v))continue;
            std::fprintf(j,"%s[%u,%u]", first?"":",", u, v); std::fprintf(gf,"L\t%u\t+\t%u\t+\t0M\n",u,v); first=false; } }
    std::fprintf(j,"]}\n"); std::fclose(j); std::fclose(gf);
}
}  // namespace

int trace_neigh(int argc, char** argv){
    if (argc < 5) die("neigh <graph.bin> <outdir> <seeds.tsv> <radius>");
    const std::string gp=argv[1], dir=argv[2], seeds=argv[3]; const uint32_t radius=(uint32_t)std::atoi(argv[4]);
    Graph g; load_graph(gp,g); Annot a; if(!load_annot(dir+"/annot.bin",a)) die("need annot.bin (run first)");
    // undirected CSR
    CSR und; und.off.assign(g.K+1,0);
    for (uint64_t e=0;e<g.E;e++){ und.off[g.edges[e].first+1]++; und.off[g.edges[e].second+1]++; }
    for (uint64_t i=0;i<g.K;i++) und.off[i+1]+=und.off[i];
    und.nbr.resize(2*g.E); std::vector<uint64_t> w(und.off.begin(),und.off.begin()+g.K);
    for (uint64_t e=0;e<g.E;e++){ const uint32_t u=g.edges[e].first,v=g.edges[e].second; und.nbr[w[u]++]=v; und.nbr[w[v]++]=u; }
    mkdir((dir+"/neigh").c_str(),0755);
    FILE* f=std::fopen(seeds.c_str(),"r"); if(!f) die("open seeds");
    char* line=nullptr; size_t cap=0; ssize_t n; std::vector<std::string> fld; uint64_t cnt=0;
    while((n=getline(&line,&cap,f))>0){ if(line[0]=='#'||n<2)continue;
        std::string ln(line,line[n-1]=='\n'?n-1:n); split_tab(ln,fld); if(fld.size()<2)continue;
        const std::string label=fld[0]; const uint32_t seed=(uint32_t)std::strtoul(fld[1].c_str(),nullptr,10);
        if(seed>=g.K)continue; emit_seed(g,a,und,seed,label,radius,600,dir+"/neigh"); cnt++; }
    std::free(line); std::fclose(f);
    std::fprintf(stderr,"[trace] neigh: emitted %lu subgraphs -> %s/neigh\n",(unsigned long)cnt,dir.c_str());
    return 0;
}

// =====================================================================================
// diag: WHY is a variant missed? Construct its alt allele from ref+VCF, compute the novel
// window signatures it should create, and probe the graph to bucket the reason. This is
// independent of marking-1's ref-neighbour attribution, so it separates a true caller miss
// from an attribution gap. Robust to breakpoint off-by-ones (tries offsets, unions sigs).
// =====================================================================================
namespace {
constexpr int CTX = 600;   // flank each side of a junction: long enough for stable level-5 LCP
                           // cores (local-consistency neighbourhood) so they match the read parse

struct DVar { uint32_t cid; uint64_t pos0, end0, svlen; std::string type, svid, alt, inseq; bool found; };

// collect interior-core window sigs of s (same window logic as the build)
void collect_sigs(const std::string& s, std::vector<uint64_t>& out) {
    const int len=(int)s.size(); if (len<2*MARGIN+4) return;
    int n; const uint64_t *st,*en; const uint64_t *lb; extract(s.data(),len,0,LEVEL,&n,&st,&en,&lb); (void)lb;
    for (int i=0;i<n;i++){ if (st[i]<(uint64_t)MARGIN||en[i]+MARGIN>(uint64_t)len) continue;
        out.push_back(emit_sig(s.data()+st[i]-FLANK,(int)(en[i]-st[i]+2*FLANK),0)); }
}
std::string sub(const Faidx& fa,const std::string& c,long a,long b){ if(a<0)a=0; if(b<=a)return ""; return fa_fetch(fa,c,(uint64_t)a,(uint64_t)(b-a)); }

// build the set of alt-allele local sequences (>=1) for a variant, plus the ref-local sequence.
void build_alleles(const Faidx& fa, const std::string& chrom, const DVar& v, std::vector<std::string>& alts, std::string& refloc) {
    const long p=(long)v.pos0, e=(long)v.end0; const std::string L=sub(fa,chrom,p-CTX,p);
    refloc = sub(fa,chrom,p-CTX,(v.type=="SNP"||v.type=="INS")?p+CTX:e+CTX);
    if (v.type=="SNP"){ alts.push_back(L + v.alt + sub(fa,chrom,p+1,p+CTX)); }
    else if (v.type=="INS"){ alts.push_back(L + v.inseq + sub(fa,chrom,p,p+CTX)); }
    else if (v.type=="DEL"){ for (int o=-2;o<=2;o++) alts.push_back(L + sub(fa,chrom,e+o,e+o+CTX)); }
    else if (v.type=="DUP"){ for (int o=-2;o<=2;o++) alts.push_back(sub(fa,chrom,e+o-CTX,e+o) + sub(fa,chrom,p,p+CTX)); }
    else if (v.type=="INV"){ // both inversion junctions
        alts.push_back(L + revcomp(sub(fa,chrom,e-CTX,e))); alts.push_back(revcomp(sub(fa,chrom,p,p+CTX)) + sub(fa,chrom,e,e+CTX)); }
}

struct DiagOut { const char* bucket; uint32_t n_novel,n_present_som,n_present_shared,n_absent,flank_tot,flank_present; };

DiagOut diag_one(const Graph& g, const std::unordered_map<uint64_t,uint32_t>& idx, const Annot& a,
                 const Faidx& fa, const std::string& chrom, const DVar& v) {
    std::vector<std::string> alts; std::string refloc; build_alleles(fa,chrom,v,alts,refloc);
    std::vector<uint64_t> rs; collect_sigs(refloc, rs);
    std::unordered_set<uint64_t> refset(rs.begin(),rs.end());
    // flank coverage: how many ref-local windows are present in the graph (region covered at all?)
    uint32_t ftot=0,fpre=0; for (uint64_t s: rs){ ftot++; if (idx.count(s)) fpre++; }
    std::unordered_set<uint64_t> novel;
    for (auto& al : alts){ std::vector<uint64_t> as; collect_sigs(al,as); for (uint64_t s: as) if(!refset.count(s)) novel.insert(s); }
    uint32_t nsom=0,nshared=0,nabs=0;
    for (uint64_t s: novel){ auto it=idx.find(s); if(it==idx.end()){nabs++;continue;}
        const uint8_t t=a.truthsrc[it->second]; if((t&OUT_ANY)&&!(t&REF_ANY))nsom++; else nshared++; }
    DiagOut d{"other",(uint32_t)novel.size(),nsom,nshared,nabs,ftot,fpre};
    if (refloc.empty()||ftot==0)            d.bucket="no_context";
    else if (novel.empty())                 d.bucket="no_novel_window";
    else if (nsom>0)                         d.bucket="sig_present_attribution_gap";
    else if (nshared>0 && nabs==0)           d.bucket="nonnovel_global_repeat";
    else if (fpre*2 < ftot)                  d.bucket="low_coverage_region";
    else                                     d.bucket="novel_absent_undercov";
    return d;
}

// parse VCF with the construction fields; join `found` from variants.tsv (svid -> found).
bool parse_vcf_diag(const std::string& path, const std::unordered_map<std::string,uint32_t>& cid,
                    const std::unordered_map<std::string,uint8_t>& found, std::vector<DVar>& vs) {
    FILE* f=std::fopen(path.c_str(),"r"); if(!f) return false;
    char* line=nullptr; size_t cap=0; ssize_t n; std::vector<std::string> fl;
    while((n=getline(&line,&cap,f))>0){ if(line[0]=='#')continue; std::string ln(line,line[n-1]=='\n'?n-1:n); split_tab(ln,fl);
        if(fl.size()<8)continue; auto c=cid.find(fl[0]); if(c==cid.end())continue;
        DVar v; v.cid=c->second; v.pos0=std::strtoull(fl[1].c_str(),nullptr,10)-1; v.svid=fl[2];
        const std::string svt=info_get(fl[7],"SVTYPE"); v.svlen=std::strtoull(info_get(fl[7],"SVLEN").c_str(),nullptr,10);
        const std::string endv=info_get(fl[7],"END"); v.end0=endv.empty()?v.pos0:std::strtoull(endv.c_str(),nullptr,10); // 0-based first-kept after span
        if(svt=="DEL")v.type="DEL"; else if(svt=="INS"){v.type="INS";v.inseq=info_get(fl[7],"INSSEQ");}
        else if(svt=="DUP")v.type="DUP"; else if(svt=="INV")v.type="INV"; else {v.type="SNP";v.alt=fl[4];}
        auto fd=found.find(v.svid); v.found = fd!=found.end() && fd->second;
        vs.push_back(std::move(v)); }
    std::free(line); std::fclose(f); return true;
}
}  // namespace

int trace_diag(int argc, char** argv){
    if (argc < 6) die("diag <graph.bin> <ref.fa> <out.vcf> <variants.tsv> <outdir>");
    const std::string gp=argv[1], ref=argv[2], vcf=argv[3], vtsv=argv[4], dir=argv[5];
    Graph g; load_graph(gp,g); std::unordered_map<uint64_t,uint32_t> idx; build_index(g,idx);
    Annot a; if(!load_annot(dir+"/annot.bin",a)) die("need annot.bin (run first)");
    Faidx fa; load_faidx(ref,fa);
    std::unordered_map<std::string,uint32_t> cid; for(uint32_t i=0;i<a.chrom_names.size();i++)cid[a.chrom_names[i]]=i;
    // found status from variants.tsv (col1 svid, col15 found)
    std::unordered_map<std::string,uint8_t> found; { FILE* f=std::fopen(vtsv.c_str(),"r"); char* L=nullptr; size_t c=0; ssize_t n; std::vector<std::string> fl; bool hdr=true;
        while((n=getline(&L,&c,f))>0){ if(hdr){hdr=false;continue;} std::string ln(L,L[n-1]=='\n'?n-1:n); split_tab(ln,fl); if(fl.size()>=15)found[fl[0]]=(uint8_t)(fl[14]=="1"); }
        std::free(L); std::fclose(f); }
    std::vector<DVar> vs; parse_vcf_diag(vcf,cid,found,vs);
    // bucket every MISSED variant; also sample 4000 FOUND as a control
    std::unordered_map<std::string,uint64_t> bucket; std::unordered_map<std::string,uint64_t> bucket_ctrl;
    FILE* out=std::fopen((dir+"/diag.tsv").c_str(),"w");
    std::fprintf(out,"svid\tchrom\tpos1\ttype\tsvlen\tfound\tbucket\tn_novel\tn_present_som\tn_present_shared\tn_absent\tflank_present\tflank_tot\n");
    uint64_t done=0, nmiss=0;
    for (auto& v : vs){ const bool ctrl = v.found && (done%29==0); if (v.found && !ctrl){done++;continue;}  // all missed + ~1/29 found as control
        DiagOut d=diag_one(g,idx,a,fa,a.chrom_names[v.cid],v);
        (v.found?bucket_ctrl:bucket)[d.bucket]++; if(!v.found)nmiss++;
        std::fprintf(out,"%s\t%s\t%lu\t%s\t%lu\t%d\t%s\t%u\t%u\t%u\t%u\t%u\t%u\n", v.svid.c_str(), a.chrom_names[v.cid].c_str(),
            (unsigned long)v.pos0+1, v.type.c_str(), (unsigned long)v.svlen, v.found?1:0, d.bucket,
            d.n_novel,d.n_present_som,d.n_present_shared,d.n_absent,d.flank_present,d.flank_tot);
        if(((++done)&8191)==0) std::fprintf(stderr,"\r[diag] %lu processed",(unsigned long)done); }
    std::fclose(out);
    std::fprintf(stderr,"\r[diag] MISSED variants (%lu) bucket breakdown:\n",(unsigned long)nmiss);
    for (auto& kv: bucket) std::fprintf(stderr,"   %-30s %8lu  (%.1f%%)\n",kv.first.c_str(),(unsigned long)kv.second,100.0*kv.second/nmiss);
    uint64_t nc=0; for(auto&kv:bucket_ctrl)nc+=kv.second;
    std::fprintf(stderr,"[diag] FOUND control (%lu sampled):\n",(unsigned long)nc);
    for (auto& kv: bucket_ctrl) std::fprintf(stderr,"   %-30s %8lu  (%.1f%%)\n",kv.first.c_str(),(unsigned long)kv.second,100.0*kv.second/nc);
    return 0;
}

// =====================================================================================
// probeseq: independent verifier. Given a file of `label<TAB>sequence` lines, parse each
// sequence with the production window logic and report how many of its core windows are
// present in the graph, colour-2 (tumor-only), and out-only (somatic truth-class). Lets an
// auditor construct any allele by hand and check whether its signatures live in the graph,
// independently of diag's own construction.
// =====================================================================================
int trace_probeseq(int argc, char** argv){
    if (argc < 4) die("probeseq <graph.bin> <annot_dir> <seqs.tsv>");
    Graph g; load_graph(argv[1], g); std::unordered_map<uint64_t,uint32_t> idx; build_index(g, idx);
    Annot a; if(!load_annot(std::string(argv[2])+"/annot.bin",a)) die("need annot.bin");
    FILE* f=std::fopen(argv[3],"r"); if(!f) die("open seqs.tsv");
    char* L=nullptr; size_t c=0; ssize_t n;
    std::printf("label\tn_cores\tn_present\tn_color1\tn_color2\tn_color3\tn_outonly\tn_error\n");
    while((n=getline(&L,&c,f))>0){ if(L[0]=='#'||n<3)continue; std::string ln(L,L[n-1]=='\n'?n-1:n);
        size_t tab=ln.find('\t'); if(tab==std::string::npos)continue;
        const std::string lab=ln.substr(0,tab); std::string s=ln.substr(tab+1);
        for(char& ch:s) if(ch>='a'&&ch<='z')ch-=32;
        const int len=(int)s.size(); uint32_t nc=0,np=0,c1=0,c2=0,c3=0,oo=0,er=0;
        if(len>=2*MARGIN+4){ int m; const uint64_t *st,*en; const uint64_t *lb; extract(s.data(),len,0,LEVEL,&m,&st,&en,&lb);(void)lb;
            for(int i=0;i<m;i++){ if(st[i]<(uint64_t)MARGIN||en[i]+MARGIN>(uint64_t)len)continue; nc++;
                const uint64_t sig=emit_sig(s.data()+st[i]-FLANK,(int)(en[i]-st[i]+2*FLANK),0);
                auto it=idx.find(sig); if(it==idx.end())continue; np++; const uint32_t nd=it->second;
                const uint8_t col=g.color[nd]; if(col==1)c1++;else if(col==2)c2++;else c3++;
                const uint8_t t=a.truthsrc[nd]; if(t==0)er++; else if((t&OUT_ANY)&&!(t&REF_ANY))oo++; } }
        std::printf("%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",lab.c_str(),nc,np,c1,c2,c3,oo,er);
    }
    std::free(L); std::fclose(f); return 0;
}

// =====================================================================================
// shapes: per-variant EXPECTED-vs-OBSERVED graph shape in the SIMPLIFIED (compacted) graph,
// a structure inventory (shape<->variant both ways), and an INV breakpoint/interior deep-dive.
// Coordinate-free: a variant's somatic nodes come from constructing its alt allele(s) (as in
// diag) and looking up the novel out-only sigs, then mapping nodes -> super-nodes whose local
// topology is classified. Shapes: simple/asym/seam bubbles, superbubble, cycle, tip, tangle.
//   shapes <graph.bin> <ref.fa> <out.vcf> <outdir>     (reuses <outdir>/annot.bin)
// =====================================================================================
namespace {
const char* const SHAPE_NAME[9]={"linear","simple","asym","seam","super","cycle","tip","tangle","orphan"};
struct STag { uint8_t shape=0, longcol=0; uint32_t len1=0, len2=0; uint32_t sid=UINT32_MAX; };
struct Struct { uint32_t s1,s2; uint8_t shape,longcol; uint32_t len1,len2,tcov; uint64_t refloc; uint32_t atc,anc; uint8_t cc; };  // cc: bit0=has normal-only(c1) arm, bit1=has tumor-only(c2) arm

// Classify every anchored region around colour-3 supers and tag arm-member supers; then a
// second pass labels leftover coloured supers as tip / tangle / orphan from local degree.
void tag_shapes(const std::vector<uint8_t>& col, const CSR& og, const CSR& ig,
                const std::vector<SuperMeta>& sm, const std::vector<uint8_t>& sblocks, std::vector<STag>& tag,
                std::vector<Struct>& structs, uint64_t* sc /*[9]*/){
    const uint64_t S=col.size(); tag.assign(S,STag{});
    std::vector<uint32_t> dvis(S,0),closed; uint32_t dgen=0;
    std::vector<uint32_t> seen(S,0),par(S,0),dist(S,0),arm; uint32_t agen=0;
    std::vector<uint8_t> a1,a2,a2rc;                              // reused arm-sequence buffers for the divergent middle
    auto measure=[&](uint32_t s1,uint32_t s2,uint8_t color)->uint32_t{ uint32_t m=0; const uint8_t smk=(uint8_t)(1u<<(color-1));
        if(bfs_arm(s1,s2,og,col,smk,0x3,500,seen,par,dist,++agen,arm)) for(uint32_t s:arm)m+=sm[s].bp; return m; };
    // peak tumor/normal coverage over a colour's divergent arm (read-error arms stay thin; a real
    // GT1|1 allele piles to ~single-copy depth; a paralog copy carries depth on BOTH colours)
    auto armcov=[&](uint32_t s1,uint32_t s2,uint8_t color,uint32_t& tc,uint32_t& nc){ tc=0;nc=0; const uint8_t smk=(uint8_t)(1u<<(color-1));
        if(bfs_arm(s1,s2,og,col,smk,0x3,500,seen,par,dist,++agen,arm)) for(uint32_t s:arm){ if(sm[s].tcov>tc)tc=sm[s].tcov; if(sm[s].ncov>nc)nc=sm[s].ncov; } };
    // splice an arm (s1->s2 order) into `out` from the super-node sequences, KMP-overlapping neighbours.
    auto splice_arm=[&](uint32_t s1,uint32_t s2,uint8_t color,std::vector<uint8_t>& out)->void{ out.clear();
        const uint8_t smk=(uint8_t)(1u<<(color-1));
        if(!bfs_arm(s1,s2,og,col,smk,0x3,500,seen,par,dist,++agen,arm)) return;
        uint32_t prevlen=0;
        for(auto it=arm.rbegin();it!=arm.rend();++it){ const uint32_t s=*it; const uint8_t* sb=sblocks.data()+sm[s].boff; const uint32_t sl=sm[s].bp;
            if(out.empty()) out.insert(out.end(),sb,sb+sl);
            else { const uint32_t lt=prevlen<(uint32_t)out.size()?prevlen:(uint32_t)out.size();
                const uint32_t k=sfx_pfx(out.data()+(out.size()-lt),lt,sb,sl); out.insert(out.end(),sb+k,sb+sl); }
            prevlen=sl; } };
    auto tagarm=[&](uint32_t s1,uint32_t s2,uint8_t color,uint8_t shp,uint8_t lc,uint32_t l1,uint32_t l2){ const uint8_t smk=(uint8_t)(1u<<(color-1));
        if(bfs_arm(s1,s2,og,col,smk,0x3,500,seen,par,dist,++agen,arm)) for(uint32_t s:arm)
            if(tag[s].shape==0||tag[s].shape>=6){tag[s].shape=shp;tag[s].longcol=lc;tag[s].len1=l1;tag[s].len2=l2;tag[s].sid=(uint32_t)structs.size();} };
    std::unordered_map<uint32_t,uint64_t> pic; std::unordered_set<uint32_t> cyc2;
    for(uint32_t s1=0;s1<(uint32_t)S;s1++){ if(col[s1]!=3)continue; pic.clear(); cyc2.clear();
        for(uint64_t i=og.off[s1];i<og.off[s1+1];i++){ const uint32_t x=og.nbr[i]; if(col[x]==3)continue;
            bool cy=false; dfs_branch(s1,x,col,og,dvis,++dgen,closed,cy); const uint8_t bc=col[x];
            for(uint32_t s2:closed){ uint64_t& w=pic[s2]; if(bc==1)w+=((uint64_t)1<<20); else w+=1; if(cy)cyc2.insert(s2);} }
        for(auto& kv:pic){ const uint32_t s2=kv.first; const uint32_t c1=(uint32_t)(kv.second>>20),c2=(uint32_t)(kv.second&0xfffff);
            const bool h1=c1>0,h2=c2>0; if(!h1&&!h2)continue; uint8_t shp=0,lc=0; uint32_t l1=0,l2=0;
            if(cyc2.count(s2)){ shp=5; if(h1)tagarm(s1,s2,1,5,0,0,0); if(h2)tagarm(s1,s2,2,5,0,0,0); }
            else if(c1>1||c2>1){ shp=4; if(h1)tagarm(s1,s2,1,4,0,0,0); if(h2)tagarm(s1,s2,2,4,0,0,0); }
            else if(h1&&h2){ splice_arm(s1,s2,1,a1); splice_arm(s1,s2,2,a2);     // divergent middle = arm minus shared prefix/suffix (the zip)
                a2rc.resize(a2.size()); for(size_t z=0;z<a2.size();z++)a2rc[z]=COMP[a2[a2.size()-1-z]];   // arms may be stored on opposite strands
                uint32_t bm=UINT32_MAX;                                          // pick the orientation that shares the most
                for(int o=0;o<2;o++){ const std::vector<uint8_t>& y=o?a2rc:a2;
                    uint32_t lcp=0,n=(uint32_t)std::min(a1.size(),y.size()); while(lcp<n&&a1[lcp]==y[lcp])lcp++;
                    uint32_t lcs=0,cap=n-lcp; while(lcs<cap&&a1[a1.size()-1-lcs]==y[y.size()-1-lcs])lcs++;
                    const uint32_t q1=(uint32_t)a1.size()-lcp-lcs, q2=(uint32_t)y.size()-lcp-lcs;
                    if(std::max(q1,q2)<bm){bm=std::max(q1,q2); l1=q1; l2=q2;} }
                uint32_t mx=std::max(l1,l2),mn=std::min(l1,l2); if(!mn)mn=1;     // indel = one-sided (ratio huge); SNP = two-sided equal (ratio ~1)
                shp=(mx<=mn*2)?1:2; lc=(l1>=l2)?1:2; tagarm(s1,s2,1,shp,lc,l1,l2); tagarm(s1,s2,2,shp,lc,l1,l2); }
            else { const uint8_t color=h1?1:2; const uint32_t l=measure(s1,s2,color); shp=3; lc=color; l1=h1?l:0; l2=h2?l:0; tagarm(s1,s2,color,3,color,l1,l2); }
            uint32_t atc=0,anc=0,tc=0,nc=0;                                  // tumor-allele / normal-allele arm depth
            if(shp==1||shp==2){ armcov(s1,s2,2,tc,nc); atc=tc; armcov(s1,s2,1,tc,nc); anc=nc; }
            else if(shp==3){ armcov(s1,s2,lc,tc,nc); if(lc==2)atc=tc; else anc=nc; }
            sc[shp]++; structs.push_back({s1,s2,shp,lc,l1,l2,std::max(sm[s1].tcov,sm[s2].tcov),sm[s1].refloc,atc,anc,(uint8_t)((h1?1:0)|(h2?2:0))}); } }
    for(uint32_t u=0;u<(uint32_t)S;u++){ if(col[u]==3||tag[u].shape!=0)continue;
        const uint32_t outd=(uint32_t)(og.off[u+1]-og.off[u]), ind=(uint32_t)(ig.off[u+1]-ig.off[u]);
        bool ina=false,outa=false;
        for(uint64_t i=ig.off[u];i<ig.off[u+1];i++) if(col[ig.nbr[i]]==3){ina=true;break;}
        for(uint64_t i=og.off[u];i<og.off[u+1];i++) if(col[og.nbr[i]]==3){outa=true;break;}
        const bool branchy=ind>1||outd>1; uint8_t shp;
        if(ina&&outa)shp=7; else if(ina||outa)shp=branchy?7:6; else shp=8;
        tag[u].shape=shp; sc[shp]++; }
}
}  // namespace

// Apply the production X-cut (resolve.cpp logic) to trace's super-graph so the shape analysis runs
// on the RESOLVED graph: split a repeat super into one copy per (in,out) pair when the reads (graph
// dir/threads.bin) pin its in-flanks to its out-flanks 1:1 AND cover its topology. Auto-applied if
// threads.bin sits next to graph.bin; otherwise the graph is left un-resolved.
static void apply_xcut(std::vector<SuperMeta>& sm, std::vector<IEdge>& sedges,
                       const std::vector<uint32_t>& nts, const std::string& threads_path) {
    FILE* f=std::fopen(threads_path.c_str(),"rb"); if(!f){ std::fprintf(stderr,"[shapes] no threads.bin -> un-resolved graph\n"); return; }
    uint64_t n; if(std::fread(&n,8,1,f)!=1){std::fclose(f);return;}
    std::vector<RepTriple> raw(n); if(n&&std::fread(raw.data(),sizeof(RepTriple),n,f)!=n){std::fclose(f);return;} std::fclose(f);
    const bool force=std::getenv("SVP_TRACE_RESOLVE_ALL");   // resolve EVERY branching repeat per observed (in,out) pair, even cross-product (measures the blow-up)
    std::unordered_map<uint32_t,std::unordered_map<uint64_t,uint32_t>> bymid;
    for(const RepTriple& t:raw){ const uint32_t smid=nts[t.mid];
        const uint32_t sin=t.in==UINT32_MAX?UINT32_MAX:nts[t.in], sout=t.out==UINT32_MAX?UINT32_MAX:nts[t.out];
        if(sin==smid||sout==smid)continue; bymid[smid][((uint64_t)sin<<32)|sout]++; }
    std::unordered_map<uint32_t,std::vector<uint32_t>> inset,outset;
    for(const IEdge& e:sedges){ outset[e.first].push_back(e.second); inset[e.second].push_back(e.first); }
    std::unordered_set<uint64_t> rm; std::vector<IEdge> add; uint64_t nsplit=0,ncopy=0,nforce=0,maxcopy=0;
    for(const auto& kv:bymid){ const uint32_t C=kv.first;
        auto ii=inset.find(C), oo=outset.find(C); if(ii==inset.end()||oo==outset.end())continue;
        const std::vector<uint32_t>& is=ii->second; const std::vector<uint32_t>& os=oo->second;
        if(is.size()<2||os.size()<2)continue;
        std::unordered_map<uint32_t,uint32_t> ind,outd; std::vector<std::pair<uint32_t,uint32_t>> M;
        for(const auto& p:kv.second){ const uint32_t a=(uint32_t)(p.first>>32),b=(uint32_t)p.first; if(a==UINT32_MAX||b==UINT32_MAX)continue; ind[a]++; outd[b]++; M.push_back({a,b}); }
        bool clean = (M.size()==is.size() && is.size()==os.size());        // perfect 1:1 matching covering the topology
        if(clean) for(uint32_t a:is) if(ind[a]!=1){ clean=false; break; }
        if(clean) for(uint32_t b:os) if(outd[b]!=1){ clean=false; break; }
        if(!clean && !force) continue;                                     // cross-product: refuse unless RESOLVE-ALL
        if(!clean) nforce++;                                               // an otherwise-refused (cross-product) repeat, forced
        nsplit++; if(M.size()>maxcopy) maxcopy=M.size();
        for(const auto& p:M){ const uint32_t ci=(uint32_t)sm.size(); sm.push_back(sm[C]); add.push_back({p.first,ci}); add.push_back({ci,p.second}); ncopy++; }
        for(uint32_t a:is)rm.insert(((uint64_t)a<<32)|C); for(uint32_t b:os)rm.insert(((uint64_t)C<<32)|b);
    }
    std::vector<IEdge> ne; ne.reserve(sedges.size()+add.size());
    for(const IEdge& e:sedges) if(!rm.count(((uint64_t)e.first<<32)|e.second)) ne.push_back(e);
    ne.insert(ne.end(),add.begin(),add.end()); std::sort(ne.begin(),ne.end()); ne.erase(std::unique(ne.begin(),ne.end()),ne.end()); sedges.swap(ne);
    std::string extra = force ? (" | "+std::to_string(nforce)+" cross-product forced, worst single repeat -> "+std::to_string(maxcopy)+" copies") : std::string();
    std::fprintf(stderr,"[shapes] X-cut applied (%s): split %lu repeat supers -> %lu copies (S now %zu)%s\n",
        force?"RESOLVE-ALL":"RESOLVED",(unsigned long)nsplit,(unsigned long)ncopy,sm.size(),extra.c_str());
}

int trace_shapes(int argc, char** argv){
    if(argc<5) die("shapes <graph.bin> <ref.fa> <out.vcf> <outdir>");
    const std::string gp=argv[1],ref=argv[2],vcf=argv[3],dir=argv[4];
    Graph g; load_graph(gp,g); std::unordered_map<uint64_t,uint32_t> idx; build_index(g,idx);
    Annot a; if(!load_annot(dir+"/annot.bin",a)) die("need annot.bin (run `trace run` first)");
    Faidx fa; load_faidx(ref,fa);
    std::vector<uint32_t> nts; std::vector<SuperMeta> sm; std::vector<IEdge> sedges; std::vector<uint8_t> sblocks;
    compact(g,a,nts,sm,sedges,&sblocks);
    apply_xcut(sm, sedges, nts, gp.substr(0, gp.find_last_of('/')+1) + "threads.bin");
    const uint64_t S=sm.size(); std::vector<uint8_t> col(S); for(uint64_t i=0;i<S;i++)col[i]=sm[i].color;
    CSR og; build_csr(sedges,S,og);
    std::vector<IEdge> rse; rse.reserve(sedges.size()); for(auto&e:sedges)rse.push_back({e.second,e.first});
    CSR ig; build_csr(rse,S,ig);
    std::vector<STag> tag; std::vector<Struct> structs; uint64_t sc[9]={0};
    tag_shapes(col,og,ig,sm,sblocks,tag,structs,sc);
    std::fprintf(stderr,"[shapes] structures: simple=%lu asym=%lu seam=%lu super=%lu cycle=%lu | supers tip=%lu tangle=%lu orphan=%lu\n",
        (unsigned long)sc[1],(unsigned long)sc[2],(unsigned long)sc[3],(unsigned long)sc[4],(unsigned long)sc[5],
        (unsigned long)sc[6],(unsigned long)sc[7],(unsigned long)sc[8]);

    std::unordered_map<std::string,uint32_t> cid; for(uint32_t i=0;i<a.chrom_names.size();i++)cid[a.chrom_names[i]]=i;
    std::unordered_map<std::string,uint8_t> nofound; std::vector<DVar> vs; parse_vcf_diag(vcf,cid,nofound,vs);
    std::vector<uint8_t> var_super(S,0);

    // gather the somatic super-nodes implied by a list of constructed alt-allele sequences
    auto somatic_supers=[&](const std::vector<std::string>& alts, const std::unordered_set<uint64_t>& refset,
                            std::unordered_set<uint32_t>& sup, uint32_t& nnodes){
        std::unordered_set<uint32_t> nodes;
        for(const auto& al:alts){ if(al.empty())continue; std::vector<uint64_t> as; collect_sigs(al,as);
            for(uint64_t s:as){ if(refset.count(s))continue; auto it=idx.find(s); if(it==idx.end())continue;
                const uint32_t nd=it->second; const uint8_t t=a.truthsrc[nd]; if(!((t&OUT_ANY)&&!(t&REF_ANY)))continue;
                if(nodes.insert(nd).second) sup.insert(nts[nd]); } }
        nnodes=(uint32_t)nodes.size(); };
    // pick the strongest-shape super for a variant (strong structures beat tip/tangle; tie -> tcov)
    auto best_obs=[&](const std::unordered_set<uint32_t>& sup, uint8_t& obs,uint8_t& lc,uint32_t& l1,uint32_t& l2,
                      uint32_t& bt,uint32_t& bm, bool& hcyc,bool& htip,bool& hsup,bool& hseam){
        obs=0;lc=0;l1=l2=0;bt=0;bm=0;hcyc=htip=hsup=hseam=false; long bestkey=-1;
        for(uint32_t s:sup){ const STag& t=tag[s];
            if(t.shape==5)hcyc=true; if(t.shape==6)htip=true; if(t.shape==4)hsup=true; if(t.shape==3)hseam=true;
            const int rank=(t.shape>=1&&t.shape<=5)?2:(t.shape?1:0);
            const long key=(long)rank*1000000L+(long)std::min<uint32_t>(sm[s].tcov,999999u);
            if(key>bestkey){bestkey=key;obs=t.shape;lc=t.longcol;l1=t.len1;l2=t.len2;bt=sm[s].tcov;bm=sm[s].members;} } };

    FILE* fv=std::fopen((dir+"/variant_shapes.tsv").c_str(),"w");
    std::fprintf(fv,"svid\ttype\tsvlen\tn_som_nodes\tn_som_supers\tobs_shape\tobs_longcol\tlen1\tlen2\thas_cycle\thas_tip\thas_super\thas_seam\tbest_tcov\tbest_members\n");
    FILE* fi=std::fopen((dir+"/inv_breakpoints.tsv").c_str(),"w");
    std::fprintf(fi,"svid\tpos1\tend1\tsvlen\tleft_found\tleft_shape\tright_found\tright_shape\tn_bp_found\tint_cores\tint_present\tint_c3\tint_c2\tint_c1\n");

    for(auto& v:vs){ const std::string& chrom=a.chrom_names[v.cid];
        std::vector<std::string> alts; std::string reflocseq; build_alleles(fa,chrom,v,alts,reflocseq);
        std::vector<uint64_t> rs; collect_sigs(reflocseq,rs); std::unordered_set<uint64_t> refset(rs.begin(),rs.end());
        if(v.type=="INV"){
            std::unordered_set<uint32_t> lsup,rsup; uint32_t ln=0,rn=0;
            if(alts.size()>0){ std::vector<std::string> la{alts[0]}; somatic_supers(la,refset,lsup,ln); }
            if(alts.size()>1){ std::vector<std::string> ra{alts[1]}; somatic_supers(ra,refset,rsup,rn); }
            uint8_t lo,ro,jlc; uint32_t j1,j2,jt,jm; bool jc,jt2,js,jsm;
            best_obs(lsup,lo,jlc,j1,j2,jt,jm,jc,jt2,js,jsm); best_obs(rsup,ro,jlc,j1,j2,jt,jm,jc,jt2,js,jsm);
            const bool lf=!lsup.empty(), rf=!rsup.empty(); const int nbp=(lf?1:0)+(rf?1:0);
            std::string interior=sub(fa,chrom,(long)v.pos0,(long)v.end0);
            std::vector<uint64_t> is; collect_sigs(interior,is); uint32_t ic=0,ip=0,i3=0,i2=0,i1=0;
            for(uint64_t s:is){ ic++; auto it=idx.find(s); if(it==idx.end())continue; ip++; const uint8_t cc=g.color[it->second];
                if(cc==3)i3++; else if(cc==2)i2++; else i1++; }
            std::fprintf(fi,"%s\t%lu\t%lu\t%lu\t%d\t%s\t%d\t%s\t%d\t%u\t%u\t%u\t%u\t%u\n", v.svid.c_str(),
                (unsigned long)v.pos0+1,(unsigned long)v.end0,(unsigned long)v.svlen, lf?1:0,SHAPE_NAME[lo], rf?1:0,SHAPE_NAME[ro], nbp, ic,ip,i3,i2,i1);
            for(uint32_t s:lsup)var_super[s]=1; for(uint32_t s:rsup)var_super[s]=1;
            std::unordered_set<uint32_t> allsup=lsup; for(uint32_t s:rsup)allsup.insert(s);
            uint8_t obs,olc; uint32_t L1,L2,bt,bm; bool hc,ht,hs,hsm; best_obs(allsup,obs,olc,L1,L2,bt,bm,hc,ht,hs,hsm);
            const char* os = allsup.empty()? "absent" : SHAPE_NAME[obs];
            std::fprintf(fv,"%s\tINV\t%lu\t%u\t%zu\t%s\t%u\t%u\t%u\t%d\t%d\t%d\t%d\t%u\t%u\n", v.svid.c_str(),(unsigned long)v.svlen,
                ln+rn,(size_t)allsup.size(), os,olc,L1,L2,hc?1:0,ht?1:0,hs?1:0,hsm?1:0,bt,bm);
            continue;
        }
        std::unordered_set<uint32_t> sup; uint32_t nn=0; somatic_supers(alts,refset,sup,nn);
        for(uint32_t s:sup)var_super[s]=1;
        uint8_t obs,olc; uint32_t L1,L2,bt,bm; bool hc,ht,hs,hsm; best_obs(sup,obs,olc,L1,L2,bt,bm,hc,ht,hs,hsm);
        const char* os = sup.empty()? "absent" : SHAPE_NAME[obs];
        std::fprintf(fv,"%s\t%s\t%lu\t%u\t%zu\t%s\t%u\t%u\t%u\t%d\t%d\t%d\t%d\t%u\t%u\n", v.svid.c_str(),v.type.c_str(),
            (unsigned long)v.svlen, nn,(size_t)sup.size(), os,olc,L1,L2, hc?1:0,ht?1:0,hs?1:0,hsm?1:0, bt,bm);
    }
    std::fclose(fv); std::fclose(fi);

    // structure inventory (shape -> variant?): a struct traces to a variant if either anchor or
    // any of its tagged arm members were touched by a variant's somatic supers.
    FILE* fs=std::fopen((dir+"/structures.tsv").c_str(),"w");
    std::fprintf(fs,"shape\tlongcol\tlen1\tlen2\ttcov\tmatched\tchrom\tpos\ta_tcov\ta_ncov\n");
    // a structure backs a variant iff one of its somatic ARM supers was attributed to a VCF variant
    // (var_super). The color-3 anchors s1/s2 are shared backbone and never somatic, so checking the
    // arm (via tag[s].sid) is the real signal; matched==0 => an "empty" shape that traces to no variant.
    std::vector<uint8_t> struct_matched(structs.size(),0);
    for(uint32_t s=0;s<S;s++) if(var_super[s] && tag[s].sid!=UINT32_MAX) struct_matched[tag[s].sid]=1;
    for(uint32_t si=0; si<structs.size(); si++){ const Struct& st=structs[si]; const uint64_t rl=st.refloc;
        const int m=(struct_matched[si]||var_super[st.s1]||var_super[st.s2])?1:0;
        std::fprintf(fs,"%s\t%u\t%u\t%u\t%u\t%d\t%s\t%ld\t%u\t%u\n", SHAPE_NAME[st.shape],st.longcol,st.len1,st.len2,st.tcov,m,
            rl==LOC_NONE?"NA":a.chrom_names[loc_chrom(rl)].c_str(), rl==LOC_NONE?-1L:(long)loc_pos(rl), st.atc,st.anc); }
    std::fclose(fs);
    // per-shape matched vs empty, AND colour content. A shape is shared (all color-3), normal-only(c1),
    // or tumor-divergent(c2 = somatic candidate). Repeat copies are all-shared and produce NO shape, so
    // 'shared' should be ~0 -- every counted shape is colour-divergent. Empties are bucketed by colour so
    // the somatic-relevant FP surface (c2) is separated from normal/germline (c1).
    uint64_t sc_m[9]={0}, n_c2[9]={0}, n_c1[9]={0}, n_sh[9]={0}, e_c2[9]={0}, e_c1[9]={0}, e_sh[9]={0};
    for(uint32_t si=0;si<structs.size();si++){ const Struct& st=structs[si]; const uint8_t sh=st.shape;
        const bool m=struct_matched[si]||var_super[st.s1]||var_super[st.s2]; if(m)sc_m[sh]++;
        if(st.cc&2){ n_c2[sh]++; if(!m)e_c2[sh]++; } else if(st.cc==1){ n_c1[sh]++; if(!m)e_c1[sh]++; } else { n_sh[sh]++; if(!m)e_sh[sh]++; } }
    for(uint32_t s=0;s<S;s++){ const uint8_t sh=tag[s].shape; if(sh<6)continue; const bool m=var_super[s]; if(m)sc_m[sh]++;
        if(sm[s].color==2){ n_c2[sh]++; if(!m)e_c2[sh]++; } else { n_c1[sh]++; if(!m)e_c1[sh]++; } }
    FILE* fc=std::fopen((dir+"/shape_census.tsv").c_str(),"w");
    std::fprintf(fc,"shape\tgranularity\tcount\tmatched\tempty\tc2_somatic\tc1_normal\tshared\tempty_c2\tempty_c1\n");
    const char* gran[9]={"-","structure","structure","structure","structure","structure","super","super","super"};
    uint64_t tot_all=0,tot_m=0,T_c2=0,T_c1=0,T_sh=0,Te_c2=0,Te_c1=0,Te_sh=0;
    for(int i=1;i<=8;i++){ std::fprintf(fc,"%s\t%s\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",SHAPE_NAME[i],gran[i],
            (unsigned long)sc[i],(unsigned long)sc_m[i],(unsigned long)(sc[i]-sc_m[i]),
            (unsigned long)n_c2[i],(unsigned long)n_c1[i],(unsigned long)n_sh[i],(unsigned long)e_c2[i],(unsigned long)e_c1[i]);
        tot_all+=sc[i]; tot_m+=sc_m[i]; T_c2+=n_c2[i]; T_c1+=n_c1[i]; T_sh+=n_sh[i]; Te_c2+=e_c2[i]; Te_c1+=e_c1[i]; Te_sh+=e_sh[i]; }
    std::fclose(fc);
    std::fprintf(stderr,"[shapes] colour filter (your ask -> drop all-shared): all-shared=%lu (dropped), colour-specific=%lu | tumor-divergent c2=%lu, normal-only c1=%lu\n",
        (unsigned long)T_sh,(unsigned long)(tot_all-T_sh),(unsigned long)T_c2,(unsigned long)T_c1);
    std::fprintf(stderr,"[shapes] EMPTIES by colour: tumor-divergent c2=%lu (the real somatic FP surface), normal-only c1=%lu, all-shared=%lu\n",
        (unsigned long)Te_c2,(unsigned long)Te_c1,(unsigned long)Te_sh);
    std::fprintf(stderr,"[shapes] empty (no variant): %lu of %lu shapes (%.1f%%) trace to no VCF variant\n",
        (unsigned long)(tot_all-tot_m),(unsigned long)tot_all, tot_all? 100.0*(double)(tot_all-tot_m)/(double)tot_all : 0.0);
    std::fprintf(stderr,"[shapes] wrote variant_shapes.tsv inv_breakpoints.tsv structures.tsv shape_census.tsv\n");
    return 0;
}
