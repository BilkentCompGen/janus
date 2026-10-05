#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "cmd.h"
#include "encoding.h"

Args::Args(int argc, char** argv) : command(argc > 1 ? argv[1] : ""), in(nullptr), out(nullptr), reads(nullptr), data(nullptr), graph(nullptr), vntr(nullptr), threads(16), bases(0), support(3) {
    for (int i = 2; i + 1 < argc; i += 2)
        if (!std::strcmp(argv[i], "-i")) in = argv[i + 1];
        else if (!std::strcmp(argv[i], "-o")) out = argv[i + 1];
        else if (!std::strcmp(argv[i], "-r")) reads = argv[i + 1];
        else if (!std::strcmp(argv[i], "-d")) data = argv[i + 1];
        else if (!std::strcmp(argv[i], "-t")) threads = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "-g")) graph = argv[i + 1];
        else if (!std::strcmp(argv[i], "-v")) vntr = argv[i + 1];
        else if (!std::strcmp(argv[i], "-b")) bases = std::strtoull(argv[i + 1], nullptr, 10);
        else if (!std::strcmp(argv[i], "-m")) support = std::atoi(argv[i + 1]);
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
    std::printf("usage: snape <pack|extract|count|graph|prune|select|refidx|place|lift|clip|patch|call> -i in|ev -o out [-r reads|index|pieces] [-d sites|ref.fa] [-t threads] [-g graph] [-v vntr.bed] [-b max bases (pack)] [-m support reads, default 3]\n");
    return 1;
}
