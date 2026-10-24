#pragma once
#include <cstdint>

static const uint32_t NA = 3;
static const uint32_t DSAME = 20;
static const uint32_t CL = 3;

struct Args {
    const char* command;
    const char *reads, *out, *normal, *counts, *graph, *sites, *index, *ref, *pieces, *paf, *aln, *vntr, *model, *hist, *minimum, *calls, *samples;
    int threads;
    uint32_t support, shallow, clonal;
    double raise_, floor_, mu;
    bool shuffle, noclade;
    Args(int argc, char** argv);
    const char* need(const char* value, const char* option) const;
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
int cmd_classify(const Args& a);
int cmd_joint(const Args& a);
int cmd_run(const Args& a);
