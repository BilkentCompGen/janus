CXX = g++
CXXFLAGS = -std=c++20 -O3 -g -Wall -Wextra -pedantic
LDFLAGS = -lz -llcptools -lhts -L$(LIBDIR)/parasail -lparasail -Wl,-rpath,$(LIBDIR)/parasail

TARGET = svp
SRCDIR = src
OBJDIR = obj
LIBDIR = lib

SOURCES = main.cpp $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(SOURCES)))

TESTDIR = tests
SRC_OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(wildcard $(SRCDIR)/*.cpp)))

.PHONY: all clean compare test help

# When 'compare' or 'test' is invoked, treat the extra args as no-op targets
ifeq (compare,$(firstword $(MAKECMDGOALS)))
  COMPARE_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(COMPARE_ARGS):;@:)
endif
ifeq (test,$(firstword $(MAKECMDGOALS)))
  TEST_NAME := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(TEST_NAME):;@:)
endif

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

VPATH = .:$(SRCDIR)

# Pattern rule for source files
$(OBJDIR)/%.o: %.cpp $(wildcard $(SRCDIR)/*.h) $(wildcard $(SRCDIR)/*.hpp)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I. -I$(LIBDIR)/lcptools -I$(LIBDIR)/parasail -c -o $@ $<

test: $(SRC_OBJS)
	$(CXX) $(CXXFLAGS) -I. -I$(LIBDIR)/lcptools -I$(LIBDIR)/parasail -o $(OBJDIR)/$(TEST_NAME) $(TESTDIR)/$(TEST_NAME).cpp $^ $(LDFLAGS) && ./$(OBJDIR)/$(TEST_NAME)

compare:
	@python3 scripts/compare_vcf.py $(COMPARE_ARGS)

help:
	@echo "Usage:"
	@echo "  make              					Build the svp binary"
	@echo "  make clean        					Remove build artifacts and outputs"
	@echo "  make compare <simulated.vcf> <detected.vcf>		Compare VCFs (reports SNP/short-indel, DEL, INS metrics)"
	@echo "  make test <name>  					Build and run a test from tests/"
	@echo "  make help         					Show this help message"

clean:
	rm -rf $(OBJDIR) $(TARGET) *.tsv outputs
