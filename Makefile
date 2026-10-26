# janus: somatic structural variants from LCP cores of HiFi reads.
#   make                 build/janus (the sources of lcptools and rapidgzip are fetched into lib/ once; needs g++, gcc, git, nasm, zlib)
#   make ref             the references in ref/ (ref/build.sh: hg38 with the chrY PARs masked, CHM13v2, their indexes)
#   make run             janus on a sample sheet: make run SAMPLES=datasets/hcc1954.tsv OUT=out/hcc1954 [REF=hg38_parmask] [THREADS=16]
#                        (the simulated cohort: SAMPLES=datasets/simulation.tsv REF=chm13v2; one tumor gives a pairwise run, several a joint run)
#   make datasets        the reads and truth sets of the paper (datasets/fetch.sh)
#   make simulation      the simulated cohort (datasets/simulate.sh)
#   make bench-setup     the benchmark environment (bench/setup.sh)
#   make bench-<step>    the experiments of the paper in order: full, competitors, consensus, coverage, simulation, runtime, supp
#   make figures         the figures, tables, source data and numbers of the paper (bench/figures.py, paper/figures/system)
#   make paper           the figures and the PDF of the paper (paper/main.pdf)
CXX := g++
CC := gcc
CXXFLAGS := -std=c++17 -DWITH_ISAL -O3 -march=native -flto=auto -Wall -Wextra -Wpedantic -MMD -MP
CFLAGS := -O3 -march=native -Wall -Wextra -Wpedantic -MMD -MP
RG := lib/rapidgzip/src
ISAL := $(RG)/external/isa-l
LCP := lib/lcptools
INC := -I$(LCP) -Isrc -isystem $(RG) -isystem $(RG)/core -isystem $(ISAL)/include -isystem $(ISAL)/igzip -isystem $(ISAL)
LDFLAGS := $(ISAL)/bin/isa-l.a -lz -lm -lpthread
LDFLAGS += $(shell command -v mold >/dev/null 2>&1 && echo -fuse-ld=mold)
BUILD := build
SRC := $(wildcard src/*.cpp)
OBJ := $(SRC:src/%.cpp=$(BUILD)/%.o) $(BUILD)/lcp_batch.o
LCPOBJ := $(BUILD)/lcptools/core.o $(BUILD)/lcptools/lps.o $(BUILD)/lcptools/encoding.o
PY := .venv-bench/bin/python
# $(call TEX,dir) runs a command in the latex-tikz image with dir mounted as the working directory
TEX = docker run --rm --log-driver none -u "$$(id -u):$$(id -g)" -e HOME=/tmp -v $(1):/w -w /w latex-tikz

all: $(BUILD)/janus

$(BUILD)/janus: $(OBJ) $(LCPOBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LCPOBJ) $(LDFLAGS)

$(BUILD)/%.o: src/%.cpp $(LCP)/config.h $(ISAL)/bin/isa-l.a | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(BUILD)/lcp_batch.o: src/lcp_batch.c $(LCP)/config.h | $(BUILD)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

# lcptools (BilkentCompGen/lcptools, biolcp branch, pinned): its sources compiled here with 64-bit positions, no labels,
# variable-length cores and the DNA alphabet; src/lcp_batch.c adds janus's batch parsing on top of it
LCPCOMMIT := 7197f492152359dd970348f1dd61ea29fb0f7993
$(LCP)/config.h.in:
	[ -d $(LCP)/.git ] || git clone -q https://github.com/BilkentCompGen/lcptools.git $(LCP)
	git -C $(LCP) checkout -q $(LCPCOMMIT)

$(LCP)/config.h: $(LCP)/config.h.in
	sed -e 's/@LABEL@/0/' -e 's/@POS@/64/' -e 's/@DCT@/1/' -e 's/@FIXED_CORE@/0/' -e 's/@FACTOR@/1.5/' -e 's/@ALPHABET_PROTEIN@/0/' -e 's/@SYMBOL_BITS@/2/' $< > $@

$(BUILD)/lcptools/%.o: $(LCP)/config.h
	mkdir -p $(BUILD)/lcptools
	$(CC) $(CFLAGS) -I$(LCP) -c $(LCP)/$*.c -o $@

# rapidgzip 0.16.0 (parallel gzip and BAM input of janus pack) and its isa-l, fetched and built once
$(ISAL)/bin/isa-l.a:
	[ -d lib/rapidgzip/.git ] || git clone -q --depth 1 --branch rapidgzip-v0.16.0 --recurse-submodules --shallow-submodules https://github.com/mxmlnkn/rapidgzip.git lib/rapidgzip
	$(MAKE) -C $(ISAL) -f Makefile.unx lib

-include $(OBJ:.o=.d)

SAMPLES ?= datasets/hcc1954.tsv
OUT ?= out/$(basename $(notdir $(SAMPLES)))
REF ?= hg38_parmask
THREADS ?= 16
run: $(BUILD)/janus
	$(BUILD)/janus run --samples $(SAMPLES) --ref ref/$(REF) --model models/janus.model --out $(OUT) --threads $(THREADS)

ref: $(BUILD)/janus
	bash ref/build.sh

datasets:
	bash datasets/fetch.sh

simulation: $(BUILD)/janus
	bash datasets/simulate.sh all

bench-setup:
	bash bench/setup.sh

# janus on the seven pairs (Fig. 2a, 2c, Fig. 4)
bench-full: $(BUILD)/janus
	bash bench/full.sh

# the other callers on the pairs (Fig. 2a, 2c, Fig. 4a)
bench-competitors:
	bash bench/competitors.sh align && bash bench/competitors.sh phase && bash bench/competitors.sh call

# the consensus call sets, the classifier with each pair held out, every call set scored (Fig. 2a, 2c, 2d)
bench-consensus:
	bash bench/consensus.sh labels && bash bench/consensus.sh without && bash bench/classifier.sh verdict && bash bench/consensus.sh with && bash bench/consensus.sh loo && bash bench/consensus.sh external

# depth and purity (Fig. 2b)
bench-coverage:
	bash bench/coverage.sh depth && bash bench/coverage.sh purity && bash bench/coverage.sh callers && bash bench/classifier.sh verdict && bash bench/coverage.sh score

# the simulated cohort, pairwise and joint (Fig. 3)
bench-simulation:
	bash bench/simulation.sh truth && bash bench/simulation.sh reads && bash bench/simulation.sh pairs && bash bench/simulation.sh joint && bash bench/simulation.sh pooled && bash bench/simulation.sh score

# the clean timing runs (Fig. 4)
bench-runtime:
	bash bench/runtime.sh

# the supplementary experiments: the 11-call-set labels, the classifier experiments and grid, the count-stage threshold sweep
bench-supp:
	bash bench/classifier.sh labels && bash bench/classifier.sh supp && bash bench/classifier.sh grid && bash bench/sweep.sh runs && bash bench/sweep.sh score

figures:
	docker image inspect latex-tikz > /dev/null 2>&1 || docker build -q -t latex-tikz bench/callers/images/latex-tikz
	docker image inspect janus-bench/phangorn > /dev/null 2>&1 || docker build -q -t janus-bench/phangorn bench/callers/images/phangorn
	$(call TEX,$(CURDIR)/paper/figures/system) sh -c "pdflatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp 20_system.tex > /dev/null && cp /tmp/20_system.pdf /w/" && cp paper/figures/system/20_system.pdf paper/figures/fig_system.pdf
	$(PY) bench/figures.py > bench/figures.log

paper: figures
	$(call TEX,$(CURDIR)/paper) sh -c "pdflatex -interaction=nonstopmode main.tex > /dev/null; bibtex main > /dev/null; pdflatex -interaction=nonstopmode main.tex > /dev/null; pdflatex -interaction=nonstopmode main.tex | { ! grep -E '^!|undefined'; }"

.PHONY: all run ref datasets simulation bench-setup bench-full bench-competitors bench-consensus bench-coverage bench-simulation bench-runtime bench-supp figures paper
