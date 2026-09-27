#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "count.h"
#include "extract.h"
#include "graph.h"
#include "prune.h"
#include "select.h"
#include "call.h"

static int usage() {
    std::printf("usage: snape <extract|count|graph|prune|select|call> -i in [-d drop] [-r reads] [-o out] [-t N] [-l L] [-m mode]\n");
    return 1;
}

int main(int argc, char** argv) {
    const char *in = nullptr, *out = nullptr, *rd = nullptr, *dr = nullptr;
    int t = 16, l = 5, md = 0;
    for (int i = 2; i + 1 < argc; i += 2)
        if (!std::strcmp(argv[i], "-i")) in = argv[i + 1];
        else if (!std::strcmp(argv[i], "-o")) out = argv[i + 1];
        else if (!std::strcmp(argv[i], "-r")) rd = argv[i + 1];
        else if (!std::strcmp(argv[i], "-t")) t = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-l")) l = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-d")) dr = argv[i + 1];
        else if (!std::strcmp(argv[i], "-m")) md = std::atoi(argv[i + 1]);
        else return usage();
    if (argc < 2 || !in) return usage();
    bool ex = !std::strcmp(argv[1], "extract");
    std::string o = out ? out : std::string(in) + (ex ? ".exc" : ".cnt");
    if (ex) return cmd_extract(in, o.c_str(), t, l);
    if (!std::strcmp(argv[1], "count")) return cmd_count(in, o.c_str(), t);
    if (!std::strcmp(argv[1], "graph")) return cmd_graph(in, out);
    if (!std::strcmp(argv[1], "prune")) return cmd_prune(in, rd, out, t);
    if (!std::strcmp(argv[1], "select")) return cmd_select(in, dr, rd, out, t, md);
    if (!std::strcmp(argv[1], "call")) return cmd_call(in, rd, dr, out);
    return usage();
}
