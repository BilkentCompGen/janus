CXX      := g++
CXXFLAGS := -std=c++17 -DWITH_ISAL -O3 -march=native -flto=auto -Wall -Wextra -Wpedantic -MMD -MP
RG       := lib/rapidgzip/src
INC      := -Ilib/lcptools -Isrc -isystem $(RG) -isystem $(RG)/core -isystem $(RG)/external/isa-l/include -isystem $(RG)/external/isa-l/igzip -isystem $(RG)/external/isa-l
LDFLAGS  := -Llib/lcptools -llcptools -Llib/rapidgzip -lisal -lm -lz
LDFLAGS  += $(shell command -v mold >/dev/null 2>&1 && echo -fuse-ld=mold)

BUILD := build
SRC := $(wildcard src/*.cpp)
OBJ := $(SRC:src/%.cpp=$(BUILD)/%.o)

$(BUILD)/janus: $(OBJ) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $(OBJ) $(LDFLAGS) -lpthread

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

-include $(OBJ:.o=.d)

clean:
	rm -rf $(BUILD)

.PHONY: clean
