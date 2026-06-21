CXX      := g++
# -DLARGE_GENOME: biolcp struct core uses 64-bit start/end (no offset truncation on a
# whole-genome reference). MUST match the flag the static lib is built with (below).
CXXFLAGS := -std=c++17 -O3 -march=native -flto=auto -DLARGE_GENOME -Wall -Wextra -Wpedantic -MMD
INC      := -Ilib/lcptools -Isrc -Ilib/btllib/include
LDFLAGS  := -Llib/lcptools -llcptools -lm
# mold linker if installed (faster links); default linker otherwise. -flto link uses CXXFLAGS.
LDFLAGS  += $(shell command -v mold >/dev/null 2>&1 && echo -fuse-ld=mold)

BUILD  := build
SRC    := $(wildcard src/*.cpp)
OBJ    := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SRC))

# btllib bloom filters (lib/btllib, untracked like lcptools). Compiled with relaxed warnings — it is not written to
# svp's -Wall -Wextra -Wpedantic. Only the plain BloomFilter path is used (status +
# bloom_filter); headers are pulled in by INC above.
BT_CXXFLAGS := -std=c++17 -O3 -march=native -flto=auto -Ilib/btllib/include
BT_OBJ      := $(BUILD)/bt_status.o $(BUILD)/bt_bloom_filter.o

# PGO (profile-guided optimization) — biggest free codegen win for the branch-heavy LCP/count
# pass. Two-step; force a full object rebuild between steps (make does NOT track flag changes):
#   make clean && make PGO=generate && make run        # instrumented run writes build/pgo/*.gcda
#   <remove build/*.o + svp, KEEP build/pgo> && make PGO=use
PGO ?=
ifeq ($(PGO),generate)
  CXXFLAGS    += -fprofile-generate=build/pgo
  BT_CXXFLAGS += -fprofile-generate=build/pgo
endif
ifeq ($(PGO),use)
  CXXFLAGS    += -fprofile-use=build/pgo -fprofile-correction -Wno-missing-profile
  BT_CXXFLAGS += -fprofile-use=build/pgo -fprofile-correction -Wno-missing-profile
endif

LEVEL   ?= 5
OUT     := outputs/svp
# box rule: never rm — move to a timestamped trash dir (resolved once per make run)
TRASH   := /home/yuio/trash/svp-clean-$(shell date +%Y%m%d-%H%M%S)

svp: $(OBJ) $(BT_OBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LDFLAGS)


$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(BUILD)/bt_status.o: lib/btllib/src/status.cpp | $(BUILD)
	$(CXX) $(BT_CXXFLAGS) -c $< -o $@
$(BUILD)/bt_bloom_filter.o: lib/btllib/src/bloom_filter.cpp | $(BUILD)
	$(CXX) $(BT_CXXFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $@

# ---------- run target (whole-genome, reads only — assemblies are never a caller input) ----------
READS_N  ?= simulation/reads/normal/hifi.fastq.gz
READS_T  ?= simulation/reads/wgs/hifi.fastq.gz
MINMULT  ?= 10
THREADS  ?= 16
run: svp
	@mkdir -p $(OUT)/wgs
	/usr/bin/time -v ./svp run $(READS_N) $(READS_T) $(OUT)/wgs --level $(LEVEL) --min-mult $(MINMULT) --threads $(THREADS)

run-list:
	@echo "read-based pipeline (reads only — no assembly in the caller):"
	@echo "  make run    normal=$(READS_N)"
	@echo "              tumor =$(READS_T)"
	@echo "  overrides: LEVEL=$(LEVEL) MINMULT=$(MINMULT) READS_N=... READS_T=..."

install-deps: install-deps-lcptools install-deps-badread install-deps-rapidgzip

# biolcp branch: DNA-optimised lcptools. Its own Makefile installs to /usr/local (needs
# root) so we build the static archive in-place here. -DLARGE_GENOME must match svp's
# CXXFLAGS (64-bit core start/end). -march=native: same box builds + runs.
install-deps-lcptools:
	cd lib/lcptools && gcc -O3 -march=native -DLARGE_GENOME -c encoding.c core.c lps.c && ar rcs liblcptools.a encoding.o core.o lps.o

# Badread HiFi read simulator (rrwick) — single-threaded, so sim_reads_badread.sh
# shards contigs across threads. Installed into a self-contained venv; the system
# python is left untouched.
install-deps-badread: lib/Badread/venv/bin/badread
lib/Badread/venv/bin/badread:
	@if [ ! -d lib/Badread ]; then \
	    echo ">> cloning Badread"; \
	    git clone --depth 1 https://github.com/rrwick/Badread.git lib/Badread; \
	fi
	cd lib/Badread && python3 -m venv venv && ./venv/bin/pip install --quiet --upgrade pip && ./venv/bin/pip install --quiet .
	@./lib/Badread/venv/bin/badread --version

# rapidgzip: parallel single-stream gzip decoder. The `graph` stage popen's it at runtime
# (open_decode in graph.cpp; default lib/rapidgzip/venv/bin/rapidgzip, override SVP_RAPIDGZIP).
# Self-contained venv so the system python is untouched; NOT linked into svp. Lives in the
# repo (lib/), never /tmp — so it survives reboots and `svp run` from the repo root works.
install-deps-rapidgzip: lib/rapidgzip/venv/bin/rapidgzip
lib/rapidgzip/venv/bin/rapidgzip:
	python3 -m venv lib/rapidgzip/venv && ./lib/rapidgzip/venv/bin/pip install --quiet --upgrade pip && ./lib/rapidgzip/venv/bin/pip install --quiet rapidgzip
	@./lib/rapidgzip/venv/bin/rapidgzip --version

# clean: build artifacts only (binary + objects). Keeps outputs/ — run results are expensive
# (a graph.bin is ~18 min) so they are NOT swept by a plain `make clean`.
clean:
	@mkdir -p $(TRASH) && mv -f $(BUILD) svp $(TRASH)/ 2>/dev/null || true
	@echo "clean: build/ + svp -> $(TRASH)"

# clean-outputs: run results under outputs/svp (graph.bin / simplify.bin / logs / experiment dirs).
clean-outputs:
	@mkdir -p $(TRASH) && mv -f $(OUT)/* $(TRASH)/ 2>/dev/null || true
	@echo "clean-outputs: $(OUT)/* -> $(TRASH)"

# clean-all: both of the above (build artifacts + run outputs).
clean-all: clean clean-outputs

.PHONY: install-deps install-deps-lcptools install-deps-badread install-deps-rapidgzip \
        clean clean-outputs clean-all run-list run

-include $(OBJ:.o=.d)
