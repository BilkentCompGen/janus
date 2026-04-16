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
DIAGDIR = tests/diagnostics
SRC_OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(wildcard $(SRCDIR)/*.cpp)))

.PHONY: all clean compare test diag help

# When 'compare' or 'test' is invoked, treat the extra args as no-op targets
ifeq (compare,$(firstword $(MAKECMDGOALS)))
  COMPARE_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(COMPARE_ARGS):;@:)
endif
ifeq (test,$(firstword $(MAKECMDGOALS)))
  TEST_NAME := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(TEST_NAME):;@:)
endif
ifeq (diag,$(firstword $(MAKECMDGOALS)))
  DIAG_NAME := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(DIAG_NAME):;@:)
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

diag: $(SRC_OBJS)
	$(CXX) $(CXXFLAGS) -I. -I$(LIBDIR)/lcptools -I$(LIBDIR)/parasail -o $(OBJDIR)/$(DIAG_NAME) $(DIAGDIR)/$(DIAG_NAME).cpp $^ $(LDFLAGS) && ./$(OBJDIR)/$(DIAG_NAME)

compare:
	@./scripts/truvari.sh $(COMPARE_ARGS)

help:
	@echo "Usage:"
	@echo "  make              					Build the svp binary"
	@echo "  make clean        					Remove build artifacts and outputs"
	@echo "  make compare <base.vcf> <comp.vcf> <ref.fa>		Compare VCFs with truvari bench"
	@echo "  make test <name>  					Build and run a test from tests/"
	@echo "  make diag <name>  					Build and run a diagnostic from tests/diagnostics/"
	@echo "  make help         					Show this help message"

clean:
	rm -rf $(OBJDIR) $(TARGET) *.tsv outputs
