#pragma once
#include "postings.h"
#include "unordered_dense.h"
#include <htslib/sam.h>
#include <cstdint>
#include <string>
#include <vector>

int cmd_bubbles(int argc, const char** argv);
int cmd_junctions(int argc, const char** argv);
int cmd_call(int argc, const char** argv);
int cmd_pick(int argc, const char** argv);
int cmd_pack(int argc, const char** argv);
int cmd_dump(int argc, const char** argv);
int cmd_svbidx(int argc, const char** argv);

struct Graph {
    uint64_t K = 0;
    std::vector<uint64_t> sig;
    std::vector<uint8_t> som;
    std::vector<uint64_t> oo, io;
    std::vector<uint32_t> on, in;
    int S = 0;
};

struct Bub { std::vector<uint32_t> som_arm, germ_arm; };

struct Chunks {
    std::vector<uint64_t> carry;
    int lines = 0;
    ankerl::unordered_dense::map<uint64_t, uint32_t> of;
    size_t n = 0;
};

struct JHit { uint32_t read, ord, chunk, st, clen, dl, dr; };

struct Ctg { std::string name, bam; };

struct Seg {
    uint32_t q;
    int32_t tid;
    hts_pos_t rs, re;
    int32_t qs, qe;
    uint8_t rev;
    double rep_s, rep_e;
};

struct Jn {
    int32_t c1, c2;
    hts_pos_t p1, p2;
    int8_t d1, d2;
    uint32_t q;
    double rep;
};

struct Ins { int32_t c; hts_pos_t p; int32_t len; uint32_t q; double rep; };

struct Stats {
    uint64_t rec = 0, unmapped = 0, secondary = 0, lowq = 0, shortaln = 0, highmm = 0, germ = 0;
    uint64_t seg = 0, split = 0, cigar_del = 0, cigar_ins = 0, ins_split = 0;
};

struct Normal {
    samFile* fp = nullptr;
    hts_idx_t* idx = nullptr;
    sam_hdr_t* hdr = nullptr;
    bam1_t* b = nullptr;
};

struct Cl { Jn j; hts_pos_t p1, p2; uint32_t supp, slot; double rep; };
struct ClIns { Ins v; int32_t len; uint32_t supp, slot; double rep; };
