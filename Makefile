CXX      := g++
CXXFLAGS := -std=c++17 -O3 -march=native -flto=auto -Wall -Wextra -Wpedantic -MMD -MP
INC      := -Ilib/lcptools -Isrc
LDFLAGS  := -Llib/lcptools -llcptools -lm
LDFLAGS  += $(shell command -v mold >/dev/null 2>&1 && echo -fuse-ld=mold)

BUILD := build
SRC := $(wildcard src/*.cpp)
OBJ := $(SRC:src/%.cpp=$(BUILD)/%.o)

$(BUILD)/snape: $(OBJ) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $(OBJ) $(LDFLAGS) -lpthread

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

-include $(OBJ:.o=.d)

clean:
	rm -rf $(BUILD)

.PHONY: clean
