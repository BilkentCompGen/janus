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

install-deps:
	$(MAKE) -C lib/lcptools

clean:
	mkdir -p /home/yuio/trash/svp-make-clean && \
	mv -f $(BUILD) svp $(OUT) /home/yuio/trash/svp-make-clean/ 2>/dev/null || true

.PHONY: install-deps clean run-list $(addprefix run-,$(SIMS))

-include $(OBJ:.o=.d)
