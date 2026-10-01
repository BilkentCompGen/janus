#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "cmd.h"

uint32_t MM = 3;   // see cmd.h

int main(int argc, char** argv) {
    const char *in = nullptr, *out = nullptr, *rd = nullptr, *dr = nullptr, *gp = nullptr, *vb = nullptr, *c = argc > 1 ? argv[1] : "";
    int t = 16;
    unsigned long long maxb = 0;   // pack -b: stop after this many bases (0: all)
    for (int i = 2; i + 1 < argc; i += 2)
        if (!std::strcmp(argv[i], "-i")) in = argv[i + 1];
        else if (!std::strcmp(argv[i], "-o")) out = argv[i + 1];
        else if (!std::strcmp(argv[i], "-r")) rd = argv[i + 1];
        else if (!std::strcmp(argv[i], "-d")) dr = argv[i + 1];
        else if (!std::strcmp(argv[i], "-t")) t = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-g")) gp = argv[i + 1];
        else if (!std::strcmp(argv[i], "-v")) vb = argv[i + 1];
        else if (!std::strcmp(argv[i], "-b")) maxb = std::strtoull(argv[i + 1], nullptr, 10);
        else if (!std::strcmp(argv[i], "-m")) MM = std::atoi(argv[i + 1]);
    if (!std::strcmp(c, "pack")) return cmd_pack(in, out, t, maxb);
    if (!std::strcmp(c, "extract")) return cmd_extract(in, out, t);
    if (!std::strcmp(c, "count")) return cmd_count(in, out, t);
    if (!std::strcmp(c, "graph")) return cmd_graph(in, out);
    if (!std::strcmp(c, "prune")) return cmd_prune(in, rd, out, t);
    if (!std::strcmp(c, "select")) return cmd_select(in, dr, rd, out, t);
    if (!std::strcmp(c, "refidx")) return cmd_refidx(in, out, t);
    if (!std::strcmp(c, "place")) return cmd_place(in, rd, dr, out, t);
    if (!std::strcmp(c, "lift")) return cmd_lift(in, dr, out);
    if (!std::strcmp(c, "clip")) return cmd_clip(in, rd, out);
    if (!std::strcmp(c, "patch")) return cmd_patch(in, rd, dr, out);
    if (!std::strcmp(c, "call")) return cmd_call(in, rd, dr, gp, out, vb, t);
    std::printf("usage: snape <pack|extract|count|graph|prune|select|refidx|place|lift|clip|patch|call> -i in|ev -o out [-r reads|index|pieces] [-d sites|ref.fa] [-t threads] [-g graph] [-v vntr.bed] [-b max bases (pack)] [-m support reads, default 3]\n");
    return 1;
}
