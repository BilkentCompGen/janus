CXX      := g++
# -DLARGE_GENOME: biolcp struct core uses 64-bit start/end (no offset truncation on a
# whole-genome reference). MUST match the flag the static lib is built with (below).
CXXFLAGS := -std=c++17 -O3 -march=native -DLARGE_GENOME -Wall -Wextra -Wpedantic -MMD
INC      := -Ilib/lcptools -Isrc
LDFLAGS  := -Llib/lcptools -llcptools -lz -lm

BUILD  := build
SRC    := $(wildcard src/*.cpp)
OBJ    := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SRC))

LEVEL   ?= 5
FLANK   ?= 10
OUT     := outputs/svp

svp: $(OBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LDFLAGS)

# Per-bubble TP/FP + sig recall vs a paired assembly oracle (rule-1 compliant).
# Reuses io / cores / MurmurHash3 objects from build/.
ANALYZE_DEPS := $(BUILD)/io.o $(BUILD)/cores.o $(BUILD)/MurmurHash3.o $(BUILD)/extract_fast.o
tools/bubble_recall: tools/bubble_recall.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(BUILD):
	mkdir -p $@

# ---------- run target (whole-genome, reads only — assemblies are never a caller input) ----------
READS_N  ?= simulation/reads/normal/hifi.fastq.gz
READS_T  ?= simulation/reads/wgs/hifi.fastq.gz
MINMULT  ?= 3
THREADS  ?= 16
run: svp
	@mkdir -p $(OUT)/wgs
	/usr/bin/time -v ./svp run $(READS_N) $(READS_T) $(OUT)/wgs --level $(LEVEL) --flank $(FLANK) --min-mult $(MINMULT) --threads $(THREADS)

run-list:
	@echo "read-based pipeline (reads only — no assembly in the caller):"
	@echo "  make run    normal=$(READS_N)"
	@echo "              tumor =$(READS_T)"
	@echo "  overrides: LEVEL=$(LEVEL) FLANK=$(FLANK) MINMULT=$(MINMULT) READS_N=... READS_T=..."

install-deps: install-deps-lcptools install-deps-pbsim3 install-deps-ccs

# biolcp branch: DNA-optimised lcptools. Its own Makefile installs to /usr/local (needs
# root) so we build the static archive in-place here. -DLARGE_GENOME must match svp's
# CXXFLAGS (64-bit core start/end). -march=native: same box builds + runs.
install-deps-lcptools:
	cd lib/lcptools && gcc -O3 -march=native -DLARGE_GENOME -c encoding.c core.c lps.c && ar rcs liblcptools.a encoding.o core.o lps.o

# pbsim3 ships generated autotools but the bundled config.guess/config.sub are
# absent — autoreconf -fi installs them, then standard configure+make.
install-deps-pbsim3: lib/pbsim3/src/pbsim
lib/pbsim3/src/pbsim:
	@if [ ! -d lib/pbsim3 ]; then \
	    echo ">> cloning pbsim3"; \
	    git clone --depth 1 https://github.com/yukiteruono/pbsim3.git lib/pbsim3; \
	fi
	cd lib/pbsim3 && autoreconf -fi && ./configure && $(MAKE)

install-deps-ccs: lib/ccs/ccs
lib/ccs/ccs:
	@mkdir -p lib/ccs
	@echo ">> downloading PacBio ccs v6.4.0 static binary"
	cd lib/ccs && wget -q https://github.com/PacificBiosciences/ccs/releases/download/v6.4.0/ccs.tar.gz && tar xzf ccs.tar.gz && rm -f ccs.tar.gz
	@./lib/ccs/ccs --version

clean:
	mkdir -p /home/yuio/trash/svp-make-clean && \
	mv -f $(BUILD) svp $(OUT) /home/yuio/trash/svp-make-clean/ 2>/dev/null || true

.PHONY: install-deps install-deps-lcptools install-deps-pbsim3 install-deps-ccs \
        clean run-list run

-include $(OBJ:.o=.d)
