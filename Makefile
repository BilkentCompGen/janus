CXX = g++
CXXFLAGS = -std=c++20 -O3 -g -Wall -Wextra -pedantic
LDFLAGS = -lz -llcptools

TARGET = svp
SRCDIR = src
OBJDIR = obj
CSVFILES = *.csv
GFAFILES = *.gfa
TSVFILES = *.tsv

SOURCES = main.cpp $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(SOURCES)))

.PHONY: all clean profile

all: $(TARGET)

profile: CXXFLAGS += -g
profile: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

VPATH = .:$(SRCDIR)

$(OBJDIR)/%.o: %.cpp $(wildcard $(SRCDIR)/*.h) $(wildcard $(SRCDIR)/*.hpp)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR) $(TARGET) $(CSVFILES) $(GFAFILES) $(TSVFILES)
