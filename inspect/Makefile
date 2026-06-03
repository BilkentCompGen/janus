CXX      := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -MMD
INC      := -I../lib/lcptools -I../src -Isrc
LDFLAGS  := -L../lib/lcptools -llcptools -lz

PARENT_DEPS := ../build/io.o ../build/cores.o ../build/MurmurHash3.o

BINS := seqs_probe bubble_census inspect

all: $(BINS)

$(PARENT_DEPS):
	$(MAKE) -C .. build/io.o build/cores.o build/MurmurHash3.o

seqs_probe: src/seqs_probe.cpp $(PARENT_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(PARENT_DEPS) $(LDFLAGS)

bubble_census: src/bubble_census.cpp $(PARENT_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(PARENT_DEPS) $(LDFLAGS)

inspect: src/inspect.cpp $(PARENT_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(PARENT_DEPS) $(LDFLAGS)

clean:
	rm -f $(BINS) src/*.d

.PHONY: all clean
