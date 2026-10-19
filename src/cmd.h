#pragma once
#include <cstdint>

static const uint32_t NA = 3;
static const uint32_t DSAME = 20;
static const uint32_t CL = 3;

struct Args {
    const char *command, *in, *out, *reads, *data, *graph, *vntr;
    int threads;
    uint32_t support, shallow, clonal;
    bool shuffle, germline;
    Args(int argc, char** argv);
};

int cmd_pack(const Args& a);
int cmd_extract(const Args& a);
int cmd_count(const Args& a);
int cmd_graph(const Args& a);
int cmd_prune(const Args& a);
int cmd_select(const Args& a);
int cmd_refidx(const Args& a);
int cmd_place(const Args& a);
int cmd_lift(const Args& a);
int cmd_clip(const Args& a);
int cmd_patch(const Args& a);
int cmd_call(const Args& a);
