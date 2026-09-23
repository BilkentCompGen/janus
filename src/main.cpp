#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "count.h"
#include "extract.h"

static int usage() {
    std::printf("usage: snape <extract|count> -i in [-o out] [-t N] [-l L]\n");
    return 1;
}

int main(int argc, char** argv) {
    const char *in = nullptr, *out = nullptr;
    int t = 16, l = 5;
    for (int i = 2; i + 1 < argc; i += 2)
        if (!std::strcmp(argv[i], "-i")) in = argv[i + 1];
        else if (!std::strcmp(argv[i], "-o")) out = argv[i + 1];
        else if (!std::strcmp(argv[i], "-t")) t = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-l")) l = std::atoi(argv[i + 1]);
        else return usage();
    if (argc < 2 || !in) return usage();
    bool ex = !std::strcmp(argv[1], "extract");
    std::string o = out ? out : std::string(in) + (ex ? ".exc" : ".cnt");
    if (ex) return cmd_extract(in, o.c_str(), t, l);
    if (!std::strcmp(argv[1], "count")) return cmd_count(in, o.c_str(), t);
    return usage();
}
