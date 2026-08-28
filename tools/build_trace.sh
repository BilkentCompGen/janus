#!/usr/bin/env bash
# Build the out-of-band analysis tool tools/trace (NOT part of `make`). Links the SAME
# prebuilt lcptools objects + flags as svp so the LCP parse + signatures are byte-identical.
# -DLARGE_GENOME is load-bearing (must match lib/lcptools build).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p build
g++ -std=c++17 -O3 -march=native -DLARGE_GENOME -Wall -Wextra \
    -Ilib/lcptools -Isrc \
    tools/trace.cpp build/extract.o build/MurmurHash3.o \
    -Llib/lcptools -llcptools -lm -lpthread \
    -o tools/trace
echo "built tools/trace"
