CXX      := g++
# -DLARGE_GENOME: biolcp struct core uses 64-bit start/end (no offset truncation on a
# whole-genome reference). MUST match the flag the static lib is built with (below).
CXXFLAGS := -std=c++17 -O3 -march=native -flto=auto -DLARGE_GENOME -Wall -Wextra -Wpedantic -MMD
INC      := -Ilib/lcptools -Isrc
# -lhts: `snape call` reads the aligned BAMs (system htslib, headers in /usr/include/htslib).
LDFLAGS  := -Llib/lcptools -llcptools -lm -lhts
# mold linker if installed (faster links); default linker otherwise. -flto link uses CXXFLAGS.
LDFLAGS  += $(shell command -v mold >/dev/null 2>&1 && echo -fuse-ld=mold)

BUILD  := build

# src/ holds the pipeline and nothing else, and all of it links into ONE binary: `snape`.
# scripts/ is measurement, never part of calling.
SRC := $(wildcard src/*.cpp)
OBJ := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SRC))

LEVEL   ?= 5
THREADS ?= 16
MINMULT ?= 3
# box rule: never rm -- move to a timestamped trash dir (resolved once per make run)
TRASH   := /home/yuio/trash/svp-clean-$(shell date +%Y%m%d-%H%M%S)

snape: $(OBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LDFLAGS) -lpthread

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(BUILD):
	mkdir -p $@

# ---------- the pipeline, end to end (reads only -- an assembly is never a caller input) ----------
# make run NORMAL=hifi/hcc1954/n.svb TUMOUR=hifi/hcc1954/t.svb DIR=outputs/svp/hcc1954
NORMAL ?=
TUMOUR ?=
REF    ?= /home/datasets/references/genomes/hg38/hg38.fa
DIR    ?= outputs/svp/run
run: snape
	@test -n "$(NORMAL)" -a -n "$(TUMOUR)" || { echo "usage: make run NORMAL=<n.svb|fq> TUMOUR=<t.svb|fq> [REF=... DIR=...]"; exit 2; }
	@mkdir -p $(DIR)
	taskset -c 0-15 ./snape run $(NORMAL) $(TUMOUR) $(REF) $(DIR) --level $(LEVEL) --min-mult $(MINMULT) --threads $(THREADS)

install-deps: install-deps-lcptools install-deps-badread

# biolcp branch: DNA-optimised lcptools. Its own Makefile installs to /usr/local (needs
# root) so we build the static archive in-place here. -DLARGE_GENOME must match snape's
# CXXFLAGS (64-bit core start/end). -march=native: same box builds + runs.
install-deps-lcptools:
	cd lib/lcptools && gcc -O3 -march=native -DLARGE_GENOME -c encoding.c core.c lps.c && ar rcs liblcptools.a encoding.o core.o lps.o

# Badread HiFi read simulator (rrwick) -- only needed to regenerate the simulation benchmark.
# Single-threaded, so scripts/sim_reads_badread.sh shards contigs across threads. Installed into a
# self-contained venv; the system python is left untouched.
install-deps-badread: lib/Badread/venv/bin/badread
lib/Badread/venv/bin/badread:
	@if [ ! -d lib/Badread ]; then \
	    echo ">> cloning Badread"; \
	    git clone --depth 1 https://github.com/rrwick/Badread.git lib/Badread; \
	fi
	cd lib/Badread && python3 -m venv venv && ./venv/bin/pip install --quiet --upgrade pip && ./venv/bin/pip install --quiet .
	@./lib/Badread/venv/bin/badread --version

# clean: build artifacts only (binary + objects). Run outputs are expensive -- a graph.bin is tens of
# minutes -- so they are NOT swept by a plain `make clean`.
clean:
	@mkdir -p $(TRASH) && mv -f $(BUILD) snape $(TRASH)/ 2>/dev/null || true
	@echo "clean: build/ + snape -> $(TRASH)"

clean-outputs:
	@mkdir -p $(TRASH) && mv -f outputs/svp/* $(TRASH)/ 2>/dev/null || true
	@echo "clean-outputs: outputs/svp/* -> $(TRASH)"

clean-all: clean clean-outputs

.PHONY: install-deps install-deps-lcptools install-deps-badread clean clean-outputs clean-all run

-include $(OBJ:.o=.d)
