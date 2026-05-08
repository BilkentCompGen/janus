CXX      := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -Wpedantic -MMD
INC      := -Ilib/lcptools -Isrc
LDFLAGS  := -Llib/lcptools -llcptools -lz

BUILD  := build
SRC    := $(wildcard src/*.cpp)
OBJ    := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SRC))

LEVEL  ?= 5
FLANK  ?= 10
OUT    := outputs/svp
NORMAL ?= simulation/ref.fasta

svp: $(OBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LDFLAGS)

# Analysis-only binary: reuses io / cores / MurmurHash3 objects from build/
ANALYZE_DEPS := $(BUILD)/io.o $(BUILD)/cores.o $(BUILD)/MurmurHash3.o
tools/analyze_cores: tools/analyze_cores.cpp $(ANALYZE_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(ANALYZE_DEPS) $(LDFLAGS)

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(BUILD):
	mkdir -p $@

# ---------- simulation auto-build (delegates to simulation/Makefile) ----------
simulation/truth/%/out.fa: ; $(MAKE) -C simulation $*

# ---------- run targets ----------
# id == configs/<id>.yaml. Output: $(OUT)/<id>/bubbles.tsv.
# N defaults to ref.fasta. Override e.g. NORMAL=simulation/truth/germline/out.fa
SIMS := del dup ins inv mixed snp

define run_rule
run-$(1): svp simulation/truth/$(1)/out.fa $$(NORMAL)
	@mkdir -p $$(OUT)/$(1)
	/usr/bin/time -v ./svp $$(NORMAL) simulation/truth/$(1)/out.fa $$(OUT)/$(1)/bubbles.tsv --level $$(LEVEL) --flank $$(FLANK)
endef

$(foreach s,$(SIMS),$(eval $(call run_rule,$(s))))

run-list:
	@echo "available run-<id> targets (N=$(NORMAL)):"
	@echo "  ids: $(SIMS)"
	@echo "  scope is decided by the yaml's overlap_regions: chr19.bed → chr19-only, genome.bed → whole genome"
	@echo "overrides: LEVEL=$(LEVEL) FLANK=$(FLANK) NORMAL=$(NORMAL)"
	@echo "  e.g. make run-mixed LEVEL=4 FLANK=20"
	@echo "  e.g. make run-mixed NORMAL=simulation/truth/germline/out.fa"

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
        clean run-list $(addprefix run-,$(SIMS))

-include $(OBJ:.o=.d)
