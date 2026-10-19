#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "cmd.h"
#include "encoding.h"

Args::Args(int argc, char** argv) : command(argc > 1 ? argv[1] : ""), in(nullptr), out(nullptr), reads(nullptr), data(nullptr), graph(nullptr), vntr(nullptr), threads(16), support(3), shallow(10), clonal(40), shuffle(false), germline(false) {
    for (int i = 2; i + 1 < argc; i += 2)
        if (!std::strcmp(argv[i], "-i")) in = argv[i + 1];
        else if (!std::strcmp(argv[i], "-o")) out = argv[i + 1];
        else if (!std::strcmp(argv[i], "-r")) reads = argv[i + 1];
        else if (!std::strcmp(argv[i], "-d")) data = argv[i + 1];
        else if (!std::strcmp(argv[i], "-t")) threads = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-g")) graph = argv[i + 1];
        else if (!std::strcmp(argv[i], "-v")) vntr = argv[i + 1];
        else if (!std::strcmp(argv[i], "-m")) support = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-N")) shallow = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-T")) clonal = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-s")) shuffle = std::atoi(argv[i + 1]) != 0;
        else if (!std::strcmp(argv[i], "-G")) germline = std::atoi(argv[i + 1]) != 0;
}

int main(int argc, char** argv) {
    Args a(argc, argv);
    const char* c = a.command;
    LCP_INIT();
    if (!std::strcmp(c, "pack")) return cmd_pack(a);
    if (!std::strcmp(c, "extract")) return cmd_extract(a);
    if (!std::strcmp(c, "count")) return cmd_count(a);
    if (!std::strcmp(c, "graph")) return cmd_graph(a);
    if (!std::strcmp(c, "prune")) return cmd_prune(a);
    if (!std::strcmp(c, "select")) return cmd_select(a);
    if (!std::strcmp(c, "refidx")) return cmd_refidx(a);
    if (!std::strcmp(c, "place")) return cmd_place(a);
    if (!std::strcmp(c, "lift")) return cmd_lift(a);
    if (!std::strcmp(c, "clip")) return cmd_clip(a);
    if (!std::strcmp(c, "patch")) return cmd_patch(a);
    if (!std::strcmp(c, "call")) return cmd_call(a);
    std::printf("usage: janus <pack|extract|count|graph|prune|select|refidx|place|lift|clip|patch|call> -i in|ev -o out [-r reads|index|pieces] [-d sites|ref.fa] [-t threads] [-g graph] [-v vntr.bed] [-m support reads, 0 = automatic (count), default 3] [-N normal depth threshold (count), default 10] [-T tumor depth threshold (count), default 40] [-s 1 shuffle (pack)] [-G 1 germline]\n");
    return 1;
}
