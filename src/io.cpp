#include "svp.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <zlib.h>

void read_fasta_flat(const std::string& path, int flank, std::string& out, std::vector<uint64_t>& rec_start, std::vector<uint64_t>& rec_end) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[FATAL] gzopen %s\n", path.c_str()); std::abort(); }
    const std::string barrier(2 * flank + 1, 'N');
    char line[1 << 16];
    bool in_record = false;
    while (gzgets(f, line, sizeof(line)) != nullptr) {
        if (line[0] == '>') {
            if (in_record) rec_end.push_back(out.size());
            if (!rec_start.empty()) out += barrier;
            rec_start.push_back(out.size());
            in_record = true;
            continue;
        }
        size_t n = std::strlen(line);
        while (n && (line[n - 1] == '\n' || line[n - 1] == '\r')) n--;
        for (size_t k = 0; k < n; k++) if (line[k] >= 'a' && line[k] <= 'z') line[k] -= 32;  // uppercase: signatures must be case-insensitive (soft-masking transparent, like lcptools cores)
        out.append(line, n);
    }
    if (in_record) rec_end.push_back(out.size());
    gzclose(f);
}

// Same as read_fasta_flat but also captures the chromosome name from each '>' header
// (first whitespace-delimited token after '>'). Used by posmap to emit chroms.tsv.
void read_fasta_flat_named(const std::string& path, int flank, std::string& out, std::vector<uint64_t>& rec_start, std::vector<uint64_t>& rec_end, std::vector<std::string>& names) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[FATAL] gzopen %s\n", path.c_str()); std::abort(); }
    const std::string barrier(2 * flank + 1, 'N');
    char line[1 << 16];
    bool in_record = false;
    while (gzgets(f, line, sizeof(line)) != nullptr) {
        if (line[0] == '>') {
            if (in_record) rec_end.push_back(out.size());
            if (!rec_start.empty()) out += barrier;
            rec_start.push_back(out.size());
            size_t k = 1;
            while (line[k] && line[k] != ' ' && line[k] != '\t' && line[k] != '\n' && line[k] != '\r') k++;
            names.emplace_back(line + 1, k - 1);
            in_record = true;
            continue;
        }
        size_t n = std::strlen(line);
        while (n && (line[n - 1] == '\n' || line[n - 1] == '\r')) n--;
        for (size_t k = 0; k < n; k++) if (line[k] >= 'a' && line[k] <= 'z') line[k] -= 32;
        out.append(line, n);
    }
    if (in_record) rec_end.push_back(out.size());
    gzclose(f);
}
