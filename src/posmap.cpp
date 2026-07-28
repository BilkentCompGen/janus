#include "posmap.h"
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace posmap {

bool Map::open(const std::string& path) {
    fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) { std::fprintf(stderr, "[posmap] cannot open %s\n", path.c_str()); return false; }
    struct stat st;
    if (fstat(fd, &st) != 0) return false;
    len = (size_t)st.st_size;
    void* m = mmap(nullptr, len, PROT_READ, MAP_PRIVATE, fd, 0);
    if (m == MAP_FAILED) { std::fprintf(stderr, "[posmap] mmap failed\n"); return false; }
    base = (const uint8_t*)m;
    std::memcpy(&N, base, 8);
    if (8 + N * sizeof(Rec) != len) {
        std::fprintf(stderr, "[posmap] size mismatch: N=%lu implies %lu B, file is %zu B (rec=%zu)\n",
                     (unsigned long)N, (unsigned long)(8 + N * sizeof(Rec)), len, sizeof(Rec));
        return false;
    }
    recs = (const Rec*)(base + 8);
    // random access over the record array; the bucket table below is the hot index
    madvise((void*)recs, len - 8, MADV_RANDOM);
    // bucket table over sig>>SHIFT: records are sorted by sig, so one linear pass fills it.
    bstart = new uint32_t[NBUCKET + 1];
    uint64_t cur = 0;
    for (uint64_t i = 0; i < N; i++) {
        const uint64_t b = recs[i].sig >> SHIFT;
        while (cur <= b) bstart[cur++] = (uint32_t)i;
    }
    while (cur <= NBUCKET) bstart[cur++] = (uint32_t)N;
    std::fprintf(stderr, "[posmap] %lu records (%zu B), bucket table %llu -> mmapped %s\n",
                 (unsigned long)N, len, (unsigned long long)NBUCKET, path.c_str());
    return true;
}

void Map::close() {
    if (bstart) { delete[] bstart; bstart = nullptr; }
    if (base) { munmap((void*)base, len); base = nullptr; }
    if (fd >= 0) { ::close(fd); fd = -1; }
}

}  // namespace posmap
