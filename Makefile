CXX      := g++
CXXFLAGS := -std=c++17 -O3 -march=native -flto=auto -Wall -Wextra -Wpedantic -MMD
INC      := -Ilib/lcptools -Isrc
LDFLAGS  := -Llib/lcptools -llcptools -lm -lhts
LDFLAGS  += $(shell command -v mold >/dev/null 2>&1 && echo -fuse-ld=mold)

BUILD  := build

SRC := $(wildcard src/*.cpp)
OBJ := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SRC))

LEVEL   ?= 5
THREADS ?= 16
MINMULT ?= 3
TRASH   := /home/yuio/trash/svp-clean-$(shell date +%Y%m%d-%H%M%S)

snape: $(OBJ)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LDFLAGS) -lpthread

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(BUILD):
	mkdir -p $@

NORMAL ?=
TUMOUR ?=
REF    ?= /home/datasets/references/genomes/hg38/hg38.fa
DIR    ?= outputs/svp/run
run: snape
	@test -n "$(NORMAL)" -a -n "$(TUMOUR)" || { echo "usage: make run NORMAL=<n.svb|fq> TUMOUR=<t.svb|fq> [REF=... DIR=...]"; exit 2; }
	@mkdir -p $(DIR)
	printf 'name\trole\tpath\nN01\tnormal\t$(NORMAL)\nT01\ttumor\t$(TUMOUR)\n' > $(DIR)/sheet.tsv
	taskset -c 0-15 ./snape run --samples $(DIR)/sheet.tsv $(REF) $(DIR) --level $(LEVEL) --min-mult $(MINMULT) --threads $(THREADS)

install-deps: install-deps-lcptools install-deps-badread

install-deps-lcptools:
	cd lib/lcptools && $(MAKE) config.h LABEL=0 POS=64 DCT=1 CORE=var ALPHABET=dna && gcc -O3 -march=native -c encoding.c core.c lps.c && ar rcs liblcptools.a encoding.o core.o lps.o

install-deps-badread: lib/Badread/venv/bin/badread
lib/Badread/venv/bin/badread:
	@if [ ! -d lib/Badread ]; then \
	    echo ">> cloning Badread"; \
	    git clone --depth 1 https://github.com/rrwick/Badread.git lib/Badread; \
	fi
	cd lib/Badread && python3 -m venv venv && ./venv/bin/pip install --quiet --upgrade pip && ./venv/bin/pip install --quiet .
	@./lib/Badread/venv/bin/badread --version

clean:
	@mkdir -p $(TRASH) && mv -f $(BUILD) snape $(TRASH)/ 2>/dev/null || true
	@echo "clean: build/ + snape -> $(TRASH)"

clean-outputs:
	@mkdir -p $(TRASH) && mv -f outputs/svp/* $(TRASH)/ 2>/dev/null || true
	@echo "clean-outputs: outputs/svp/* -> $(TRASH)"

clean-all: clean clean-outputs

.PHONY: install-deps install-deps-lcptools install-deps-badread clean clean-outputs clean-all run

-include $(OBJ:.o=.d)
