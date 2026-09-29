#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "count.h"
#include "extract.h"
#include "files.h"
#include "graph.h"
#include "pack.h"
#include "place.h"
#include "prune.h"
#include "reads.h"
#include "select.h"
#include "call.h"

static int usage() {
    std::printf("usage: snape <pack|extract|count|graph|prune|ids|select|refidx|place|lift|call> -i in [-d drop|ref.fa] [-r reads|index|pieces] [-o out] [-t N] [-l L] [-m mode] [-w weights] [-v vntr.bed]\n");
    return 1;
}

static int ids(const char* in, const char* rd, int t) {
    gph g;
    read_gph(in, &g);
    write_ids(g.g, rd, in, t);
    return 0;
}

int main(int argc, char** argv) {
    const char *in = nullptr, *out = nullptr, *rd = nullptr, *dr = nullptr, *wf = nullptr, *vb = nullptr;
    int t = 16, l = 5, md = 0;
    for (int i = 2; i + 1 < argc; i += 2)
        if (!std::strcmp(argv[i], "-i")) in = argv[i + 1];
        else if (!std::strcmp(argv[i], "-o")) out = argv[i + 1];
        else if (!std::strcmp(argv[i], "-r")) rd = argv[i + 1];
        else if (!std::strcmp(argv[i], "-t")) t = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-l")) l = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-d")) dr = argv[i + 1];
        else if (!std::strcmp(argv[i], "-m")) md = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-w")) wf = argv[i + 1];
        else if (!std::strcmp(argv[i], "-v")) vb = argv[i + 1];
        else return usage();
    if (argc < 2 || !in) return usage();
    bool ex = !std::strcmp(argv[1], "extract");
    std::string o = out ? out : std::string(in) + (ex ? ".exc" : ".cnt");
    if (ex) return cmd_extract(in, o.c_str(), t, l);
    if (!std::strcmp(argv[1], "count")) return cmd_count(in, o.c_str(), t);
    if (!std::strcmp(argv[1], "graph")) return cmd_graph(in, out);
    if (!std::strcmp(argv[1], "prune")) return cmd_prune(in, rd, out, t);
    if (!std::strcmp(argv[1], "pack")) return cmd_pack(in, out);
    if (!std::strcmp(argv[1], "ids")) return ids(in, rd, t);
    if (!std::strcmp(argv[1], "select")) return cmd_select(in, dr, rd, out, t, md);
    if (!std::strcmp(argv[1], "refidx")) return cmd_refidx(in, out, t);
    if (!std::strcmp(argv[1], "place")) return cmd_place(in, rd, dr, out, t);
    if (!std::strcmp(argv[1], "lift")) return cmd_lift(in, rd, dr, out);
    if (!std::strcmp(argv[1], "call")) return cmd_call(in, rd, dr, out, wf, vb);
    return usage();
}
