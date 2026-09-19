#include "log.h"

#include <stdlib.h>
#include <cstring>

void say(const char* tag, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    std::fprintf(stderr, "[%s] ", tag);
    std::vfprintf(stderr, fmt, ap);
    std::fputc('\n', stderr);
    va_end(ap);
}

[[noreturn]] void die(const char* tag, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    std::fprintf(stderr, "[%s] FATAL: ", tag);
    std::vfprintf(stderr, fmt, ap);
    std::fputc('\n', stderr);
    va_end(ap);
    std::exit(1);
}

FILE* xopen(const char* path, const char* mode) {
    FILE* f = std::fopen(path, mode);
    if (!f) die("io", "open %s", path);
    return f;
}

void wr(FILE* f, const void* p, size_t n) {
    if (n && std::fwrite(p, 1, n, f) != n) die("io", "fwrite");
}

void atomic_close(FILE* f, const char* tmp, const char* dst) {
    std::fclose(f);
    if (std::rename(tmp, dst)) die("io", "rename %s", dst);
}

bool is_svb(const char* p) {
    const size_t n = std::strlen(p);
    return n >= 4 && std::strcmp(p + n - 4, ".svb") == 0;
}

bool is_gz(const char* p) {
    const size_t n = std::strlen(p);
    return n >= 3 && std::strcmp(p + n - 3, ".gz") == 0;
}

std::string dir_of(const std::string& path) {
    const size_t c = path.find_last_of('/');
    return c == std::string::npos ? std::string(".") : path.substr(0, c);
}

int level_env(void) {
    const char* e = std::getenv("SNAPE_LEVEL");
    return e ? std::atoi(e) : 5;
}

size_t flush_env(void) {
    const char* e = std::getenv("SNAPE_FLUSH");
    return e ? (size_t)std::atoll(e) : (size_t)(1u << 20);
}
