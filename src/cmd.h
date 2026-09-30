#pragma once

int cmd_pack(const char* fastq, const char* out, int t);
int cmd_extract(const char* reads, const char* out, int t);
int cmd_count(const char* reads, const char* out, int t);
int cmd_graph(const char* counts, const char* out);
int cmd_prune(const char* graph, const char* reads, const char* out, int t);
int cmd_select(const char* graph, const char* sites, const char* reads, const char* out, int t);
int cmd_refidx(const char* ref, const char* out, int t);
int cmd_place(const char* pieces, const char* idx, const char* ref, const char* out, int t);
int cmd_lift(const char* paf, const char* ref, const char* out);
int cmd_clip(const char* in, const char* paf, const char* out);
int cmd_patch(const char* in, const char* paf, const char* own, const char* out);
int cmd_call(const char* in, const char* rd, const char* fasta, const char* graph, const char* out, const char* vntr, int t);
