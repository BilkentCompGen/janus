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
        out.append(line, n);
    }
    if (in_record) rec_end.push_back(out.size());
    gzclose(f);
}

static void emit_seq(FILE* f, std::string_view s) {
    if (s.empty()) std::fputc('.', f);
    else std::fwrite(s.data(), 1, s.size(), f);
}

void write_bubbles(const std::string& path, const std::vector<Bubble>& bubs) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) { std::fprintf(stderr, "[FATAL] fopen %s\n", path.c_str()); std::abort(); }
    std::fputs("#type\tn_start\tn_end\tt_start\tt_end\tn_seq\tt_seq\n", f);
    for (const Bubble& b : bubs) {
        std::fprintf(f, "%s\t%lu\t%lu\t%lu\t%lu\t", b.type, (unsigned long)b.n_start, (unsigned long)b.n_end, (unsigned long)b.t_start, (unsigned long)b.t_end);
        emit_seq(f, b.n_seq); std::fputc('\t', f);
        emit_seq(f, b.t_seq); std::fputc('\n', f);
    }
    std::fclose(f);
}
