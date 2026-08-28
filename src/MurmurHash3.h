#pragma once
#include <cstdint>

uint64_t murmur64(const void* key, int len, uint32_t seed);
