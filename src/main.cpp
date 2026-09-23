#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "count.h"
#include "extract.h"

struct opts { const char* in; const char* out; int t, l; };

static void usage() { printf("usage: snape <pack|extract|count> [-i in] [-o out] [-t N] [-l L]\n"); }

static bool parse_opts(int argc, char** argv, opts* o) {
    for (int i = 0; i < argc; i++)
        if (!strcmp(argv[i], "-t")) o->t = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-l")) o->l = atoi(argv[++i]);
        else if ((!strcmp(argv[i], "-i") || !strcmp(argv[i], "--in")) && !o->in) o->in = argv[++i];
        else if ((!strcmp(argv[i], "-o") || !strcmp(argv[i], "--out")) && !o->out) o->out = argv[++i];
        else return false;
    return true;
}

void pack(const opts* o) {
    printf("snape pack: %d threads\n", o->t);
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    opts o = {nullptr, nullptr, 16, 5};
    if (!parse_opts(argc - 2, argv + 2, &o)) { usage(); return 1; }
    if (!strcmp(argv[1], "pack")) pack(&o);
    else if (!strcmp(argv[1], "extract")) {
        if (!o.in) { usage(); return 1; }
        return cmd_extract(o.in, o.out, o.t, o.l);
    }
    else if (!strcmp(argv[1], "count")) {
        if (!o.in) { usage(); return 1; }
        return cmd_count(o.in, o.out, o.t);
    }
    else { usage(); return 1; }
    return 0;
}
