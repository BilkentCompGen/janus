CXX      := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -Wpedantic -MMD
INC      := -Ilib/lcptools -Isrc
LDFLAGS  := -Llib/lcptools -llcptools -lz

BUILD  := build
SRC    := $(wildcard src/*.cpp)
OBJ    := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SRC))

LEVEL   ?= 5
FLANK   ?= 10
OUT     := outputs/svp

svp: $(OBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LDFLAGS)

# Analysis-only binary: reuses io / cores / MurmurHash3 objects from build/
ANALYZE_DEPS := $(BUILD)/io.o $(BUILD)/cores.o $(BUILD)/MurmurHash3.o
tools/analyze_cores: tools/analyze_cores.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

tools/analyze_uniqueness: tools/analyze_uniqueness.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

tools/analyze_kmer: tools/analyze_kmer.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

tools/dump_gfa: tools/dump_gfa.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

tools/gfa_collapse: tools/gfa_collapse.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

tools/analyze_tangle: tools/analyze_tangle.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

tools/call: tools/call.cpp $(BUILD)/MurmurHash3.o
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(BUILD)/MurmurHash3.o $(LDFLAGS)

# Per-bubble TP/FP + sig recall vs a paired assembly oracle (rule-1 compliant).
tools/bubble_recall: tools/bubble_recall.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

# Per-block (Mode A super-node) TP/FP + sig recall vs the same assembly oracle.
tools/block_recall: tools/block_recall.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

# LCP-topology experiment binaries (read-only probes; docs/LCP_TOPOLOGY.md)
LCP_EXP := lcp_bags lcp_triplet_bags lcp_collision_anatomy lcp_perturb lcp_rc_invariance lcp_crosslevel lcp_repeat_probe lcp_consistency
tools/lcp_%: tools/lcp_%.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)
lcp-exp: $(addprefix tools/,$(LCP_EXP))

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

install-deps-lcptools:
	$(MAKE) -C lib/lcptools

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
