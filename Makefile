CXX = g++
CXXFLAGS = -std=c++20 -O3 -Wall -Wextra -pedantic -Isrc
LDFLAGS = -lz

TARGET = svp
SRCDIR = src
OBJDIR = obj

SOURCES = main.cpp $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(SOURCES)))

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

VPATH = .:$(SRCDIR)

$(OBJDIR)/%.o: %.cpp $(wildcard $(SRCDIR)/*.h) $(wildcard $(SRCDIR)/*.hpp)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR) $(TARGET)
