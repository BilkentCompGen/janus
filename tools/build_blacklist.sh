#!/usr/bin/env bash
# Build tools/build_blacklist (NOT part of `make`). Links the SAME extract.o + lcptools as svp so
# the LCP parse + XXH3 core signatures are byte-identical. -DLARGE_GENOME must match lib/lcptools.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p build
g++ -std=c++17 -O3 -march=native -DLARGE_GENOME -Wall -Wextra \
    -Ilib/lcptools -Isrc \
    tools/build_blacklist.cpp build/extract.o \
    -Llib/lcptools -llcptools -lm -lpthread \
    -o tools/build_blacklist
echo "built tools/build_blacklist"
