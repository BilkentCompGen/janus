#pragma once

int cmd_pack(const char* in, const char* out);
int cmd_extract(const char* in, const char* out, int t);
int cmd_count(const char* in, const char* out, int t);
int cmd_graph(const char* in, const char* out);
int cmd_prune(const char* in, const char* rd, const char* out, int t);
int cmd_select(const char* in, const char* dr, const char* rd, const char* out, int t);
int cmd_refidx(const char* in, const char* out, int t);
int cmd_place(const char* in, const char* idx, const char* ref, const char* out, int t);
int cmd_lift(const char* paf, const char* on, const char* ref, const char* out);
int cmd_call(const char* paf, const char* evi, const char* fasta, const char* out, const char* w, const char* v);
