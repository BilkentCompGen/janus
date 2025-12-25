CXX = g++
CXXFLAGS = -std=c++20 -O3 -g -Wall -Wextra -pedantic
LDFLAGS = -lz -llcptools -lhts -L$(LIBDIR)/parasail -lparasail -Wl,-rpath,$(LIBDIR)/parasail

TARGET = svp
SRCDIR = src
OBJDIR = obj
LIBDIR = lib

SOURCES = main.cpp $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(SOURCES)))

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

VPATH = .:$(SRCDIR)

# Pattern rule for source files
$(OBJDIR)/%.o: %.cpp $(wildcard $(SRCDIR)/*.h) $(wildcard $(SRCDIR)/*.hpp)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(LIBDIR)/lcptools -I$(LIBDIR)/parasail -c -o $@ $<

images:
	./scripts/generate_images.sh

clean:
	rm -rf $(OBJDIR) $(TARGET) *.tsv outputs
