#pragma once
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string>

void say(const char* tag, const char* fmt, ...);
[[noreturn]] void die(const char* tag, const char* fmt, ...);

FILE* xopen(const char* path, const char* mode);
void wr(FILE* f, const void* p, size_t n);
void atomic_close(FILE* f, const char* tmp, const char* dst);

bool is_svb(const char* p);
bool is_gz(const char* p);
std::string dir_of(const std::string& path);

int level_env(void);
size_t flush_env(void);
