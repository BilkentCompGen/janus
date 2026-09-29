#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "cmd.h"

int main(int argc, char** argv) {
    const char *in = nullptr, *out = nullptr, *rd = nullptr, *dr = nullptr, *wf = nullptr, *vb = nullptr, *c = argc > 1 ? argv[1] : "";
    int t = 16;
    for (int i = 2; i + 1 < argc; i += 2)
        if (!std::strcmp(argv[i], "-i")) in = argv[i + 1];
        else if (!std::strcmp(argv[i], "-o")) out = argv[i + 1];
        else if (!std::strcmp(argv[i], "-r")) rd = argv[i + 1];
        else if (!std::strcmp(argv[i], "-d")) dr = argv[i + 1];
        else if (!std::strcmp(argv[i], "-t")) t = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-w")) wf = argv[i + 1];
        else if (!std::strcmp(argv[i], "-v")) vb = argv[i + 1];
    if (!std::strcmp(c, "pack")) return cmd_pack(in, out);
    if (!std::strcmp(c, "extract")) return cmd_extract(in, out, t);
    if (!std::strcmp(c, "count")) return cmd_count(in, out, t);
    if (!std::strcmp(c, "graph")) return cmd_graph(in, out);
    if (!std::strcmp(c, "prune")) return cmd_prune(in, rd, out, t);
    if (!std::strcmp(c, "select")) return cmd_select(in, dr, rd, out, t);
    if (!std::strcmp(c, "refidx")) return cmd_refidx(in, out, t);
    if (!std::strcmp(c, "place")) return cmd_place(in, rd, dr, out, t);
    if (!std::strcmp(c, "lift")) return cmd_lift(in, rd, dr, out);
    if (!std::strcmp(c, "call")) return cmd_call(in, rd, dr, out, wf, vb);
    std::printf("usage: snape <pack|extract|count|graph|prune|select|refidx|place|lift|call> -i in -o out [-r reads|index|pieces] [-d drop|ref.fa] [-t threads] [-w weights] [-v vntr.bed]\n");
    return 1;
}
