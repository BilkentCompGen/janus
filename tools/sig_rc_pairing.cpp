// sig_rc_pairing.cpp — For each sig X in counts.bin, RC its stored byte window and
// check whether ANY other sig in counts.bin has matching bytes. Pairing is via a
// per-window murmur64 with seed=0 (label-independent — see [[project_label_seed_redundant]]).
// Reports pairing rate stratified by per-sig genomic provenance class (BOTH / N_only
// / T_only / TRUE_PHANTOM), to test whether RC-pairing distinguishes genuine sigs
// from sequencing-error phantoms.
//
// usage: sig_rc_pairing <counts.bin> <provenance.bin> [out_pairing.bin]

#include "svp.h"
#include "MurmurHash3.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

static char comp(char c) {
    switch (c & 0xDF) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
    }
    return 'N';
}

struct Prov { uint32_t n_fwd, n_rc, t_fwd, t_rc; };

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: sig_rc_pairing <counts.bin> <provenance.bin> [out_pairing.bin]\n");
        return 2;
    }
    const char* cp = argv[1]; const char* pp = argv[2];
    const char* op = (argc >= 4) ? argv[3] : nullptr;

    // counts.bin: u64 K | K×SigCount | u64 total | u16 seqlen[K] | u8 payload[total]
    FILE* fc = std::fopen(cp, "rb");
    if (!fc) { std::fprintf(stderr, "[FATAL] open %s\n", cp); std::abort(); }
    uint64_t K; std::fread(&K, 8, 1, fc);
    std::vector<SigCount> sigs(K);
    std::fread(sigs.data(), sizeof(SigCount), K, fc);
    uint64_t total; std::fread(&total, 8, 1, fc);
    std::vector<uint16_t> seqlen(K); std::fread(seqlen.data(), 2, K, fc);
    std::vector<uint8_t> payload(total); std::fread(payload.data(), 1, total, fc);
    std::fclose(fc);
    std::fprintf(stderr, "[rcp] loaded counts.bin: K=%lu, total=%lu\n", (unsigned long)K, (unsigned long)total);

    // provenance.bin
    FILE* fp = std::fopen(pp, "rb");
    if (!fp) { std::fprintf(stderr, "[FATAL] open %s\n", pp); std::abort(); }
    uint64_t Kp; std::fread(&Kp, 8, 1, fp);
    if (Kp != K) { std::fprintf(stderr, "[FATAL] K mismatch %lu vs %lu\n", (unsigned long)Kp, (unsigned long)K); std::abort(); }
    std::vector<Prov> prov(K);
    std::fread(prov.data(), sizeof(Prov), K, fp);
    std::fclose(fp);

    std::vector<uint64_t> seq_off(K + 1, 0);
    for (uint64_t i = 0; i < K; i++) seq_off[i + 1] = seq_off[i] + seqlen[i];

    // Pass 1: hash each sig's bytes with seed=0 (label-independent) and build a map.
    std::fprintf(stderr, "[rcp] hashing %lu sigs with seed=0...\n", (unsigned long)K);
    std::unordered_map<uint64_t, uint32_t> hash_to_id;
    hash_to_id.reserve(K * 2);
    std::vector<uint64_t> seed0_hash(K);
    for (uint32_t i = 0; i < K; i++) {
        const uint8_t* b = payload.data() + seq_off[i];
        seed0_hash[i] = murmur64((const char*)b, seqlen[i], 0);
        hash_to_id.emplace(seed0_hash[i], i);
    }
    std::fprintf(stderr, "[rcp] hash_to_id size = %zu (collisions = %lu)\n",
                 hash_to_id.size(), (unsigned long)K - hash_to_id.size());

    // Pass 2: for each sig, RC its bytes, hash, look up.
    std::fprintf(stderr, "[rcp] RC-pairing pass...\n");
    std::vector<uint8_t> rcbuf(2048);
    std::vector<uint8_t> paired(K, 0);  // 0 = no RC partner, 1 = partner found, 2 = partner = self (palindrome)
    for (uint32_t i = 0; i < K; i++) {
        const uint16_t l = seqlen[i];
        if ((size_t)l > rcbuf.size()) rcbuf.resize(l);
        const uint8_t* b = payload.data() + seq_off[i];
        for (uint16_t j = 0; j < l; j++) rcbuf[j] = (uint8_t)comp((char)b[l - 1 - j]);
        const uint64_t h = murmur64((const char*)rcbuf.data(), l, 0);
        auto it = hash_to_id.find(h);
        if (it == hash_to_id.end()) paired[i] = 0;
        else if (it->second == i)   paired[i] = 2;
        else                        paired[i] = 1;
    }

    // Optional: write pairing flags.
    if (op) {
        FILE* fo = std::fopen(op, "wb");
        if (!fo) { std::fprintf(stderr, "[FATAL] open %s\n", op); std::abort(); }
        std::fwrite(&K, 8, 1, fo);
        std::fwrite(paired.data(), 1, K, fo);
        std::fclose(fo);
        std::fprintf(stderr, "[rcp] wrote %s\n", op);
    }

    // Tabulate by provenance class.
    auto cls = [&](uint64_t i) -> int {
        const Prov& p = prov[i];
        const bool n = (p.n_fwd + p.n_rc) > 0;
        const bool t = (p.t_fwd + p.t_rc) > 0;
        if (n && t) return 0;   // BOTH
        if (n)      return 1;   // N_only
        if (t)      return 2;   // T_only
        return 3;               // TRUE_PHANTOM
    };
    const char* CN[] = {"BOTH", "N_only", "T_only", "TRUE_PHANTOM"};
    uint64_t total_by[4] = {0,0,0,0};
    uint64_t paired_by[4] = {0,0,0,0};
    uint64_t palin_by[4]  = {0,0,0,0};
    uint64_t unpaired_by[4]={0,0,0,0};
    for (uint64_t i = 0; i < K; i++) {
        int c = cls(i);
        total_by[c]++;
        if      (paired[i] == 1) paired_by[c]++;
        else if (paired[i] == 2) palin_by[c]++;
        else                     unpaired_by[c]++;
    }
    std::printf("== RC-pairing by provenance class ==\n");
    std::printf("  %-14s  %-12s  %-18s  %-18s  %s\n",
                "class", "n_sigs", "RC partner found", "palindromic-self", "no RC partner");
    for (int c = 0; c < 4; c++) {
        uint64_t t = total_by[c]; if (!t) continue;
        std::printf("  %-14s  %-12lu  %-18s  %-18s  %s\n",
                    CN[c], (unsigned long)t,
                    (std::to_string(paired_by[c])   + "  (" + std::to_string(100.0 * paired_by[c]   / t).substr(0,5) + "%)").c_str(),
                    (std::to_string(palin_by[c])    + "  (" + std::to_string(100.0 * palin_by[c]    / t).substr(0,5) + "%)").c_str(),
                    (std::to_string(unpaired_by[c]) + "  (" + std::to_string(100.0 * unpaired_by[c] / t).substr(0,5) + "%)").c_str());
    }

    // Cross-tabulate: phantom WITH RC partner — is its partner also a phantom? Could indicate
    // systematic errors on both strands.
    if (total_by[3] > 0) {
        uint64_t ph_partner_in_class[4] = {0,0,0,0};
        for (uint64_t i = 0; i < K; i++) {
            if (cls(i) != 3 || paired[i] != 1) continue;
            const uint16_t l = seqlen[i];
            if ((size_t)l > rcbuf.size()) rcbuf.resize(l);
            const uint8_t* b = payload.data() + seq_off[i];
            for (uint16_t j = 0; j < l; j++) rcbuf[j] = (uint8_t)comp((char)b[l - 1 - j]);
            const uint64_t h = murmur64((const char*)rcbuf.data(), l, 0);
            auto it = hash_to_id.find(h);
            if (it == hash_to_id.end()) continue;
            ph_partner_in_class[cls(it->second)]++;
        }
        std::printf("\n== TRUE_PHANTOM sigs that DO have an RC partner — class of the partner ==\n");
        for (int c = 0; c < 4; c++) {
            std::printf("  partner = %-14s = %lu\n", CN[c], (unsigned long)ph_partner_in_class[c]);
        }
    }
    return 0;
}
