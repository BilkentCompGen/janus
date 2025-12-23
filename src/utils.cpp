#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdint>

#include "utils.h"
#include "zstr.hpp"
#include "core.h"

using namespace std;

//////////////////////////////////////////////////////////// Logger ////////////////////////////////////////////////////////////
Logger::Logger(const string& name) : name(name) {}

void Logger::log_impl(const string& msg, const char* file, int line) const {
    auto tp_ms = chrono::time_point_cast<chrono::milliseconds>(chrono::system_clock::now());
    cout << format("{:%H:%M:%S}", tp_ms) << " | " << name << " | " << file << ":" << line << " | " << msg << endl;
}

/////////////////////////////////////////////////////////// Gzip I/O ///////////////////////////////////////////////////////////
unique_ptr<istream> open_input_stream(const string& file_path, bool is_gzipped) {
    if (is_gzipped) return make_unique<gzifstream>(file_path);
    else return make_unique<ifstream>(file_path);
}

////////////////////////////////////////////////////////// FASTA Reader ////////////////////////////////////////////////////////
static inline void ltrim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !isspace(ch);
    }));
}

static inline void rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !isspace(ch);
    }).base(), s.end());
}

static inline void trim(string &s) {
    ltrim(s);
    rtrim(s);
}

FastaStreamer::FastaStreamer(const string& file_path, bool is_gzipped) {
    input_stream = open_input_stream(file_path, is_gzipped);
    if (!*input_stream) {
        throw runtime_error("Could not open file: " + file_path);
    }
}

FastaIterator FastaStreamer::begin() {
    return FastaIterator(input_stream.get());
}

FastaIterator FastaStreamer::end() {
    return FastaIterator(nullptr);
}

FastaIterator::FastaIterator(istream* in) : input(in), at_end(in == nullptr) {
    if (!at_end) read_next_sequence();
}

FastaIterator& FastaIterator::operator++() {
    read_next_sequence();
    return *this;
}

const FastaPair& FastaIterator::operator*() const {
    return current_pair;
}

bool FastaIterator::operator!=(const FastaIterator& other) const {
    return at_end != other.at_end;
}

void FastaIterator::read_next_sequence() {
    string line;
    current_pair.first.clear();
    current_pair.second.clear();

    if (input->peek() == EOF) {
        at_end = true;
        return;
    }

    while (getline(*input, line)) {
        trim(line);
        if (line.empty()) continue;
        if (line[0] == '>') {
            current_pair.first = line.substr(1);
            break;
        }
    }

    if (current_pair.first.empty()) {
        at_end = true;
        return;
    }

    while (input->peek() != '>' && input->peek() != EOF) {
        getline(*input, line);
        trim(line);
        if(!line.empty()) current_pair.second += line;
    }
}


vector<FastaPair> read_fasta_file(const string& file_path, bool is_gzipped) {
    vector<FastaPair> fasta_data;
    FastaStreamer streamer(file_path, is_gzipped);
    for (const FastaPair& pair : streamer) fasta_data.push_back(pair);
    return fasta_data;
}

//////////////////////////////////////////////////////// Argument Parser /////////////////////////////////////////////////////
ArgParser::ArgParser(const string& description) : program_description(description) {}

void ArgParser::add(const string& name, const string& flag, const string& short_flag, const string& description, const string& type, bool required) {
    definitions[name] = {flag, short_flag, description, type, required};
    flag_to_name[flag] = name;
    flag_to_name[short_flag] = name;
}

void ArgParser::print_help() {
    cerr << "Usage: " << program_name << " [options]" << endl;
    if (!program_description.empty()) cerr << endl << program_description << endl;
    cerr << endl << "Options:" << endl;
    for (const pair<const string, Argument>& pair : definitions) {
        const Argument& arg = pair.second;
        cerr << "  " << arg.flag << ", " << arg.short_flag << "\t" << arg.description << " (type: " << arg.type << ")" << (arg.required ? " [required]" : "") << endl;
    }
    cerr << "  --help, -h\tShow this help message and exit" << endl;
    exit(0);
}

map<string, any> ArgParser::parse(int argc, char** argv) {
    program_name = argv[0];
    if (argc == 1) print_help();
    
    map<string, any> results;
    
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") print_help();
        
        if (!flag_to_name.count(arg)) continue;
        
        string name = flag_to_name[arg];
        const Argument& def = definitions[name];
        
        if (def.type == "bool") {
            results[name] = true;
        } else {
            if (i + 1 >= argc) throw runtime_error("Missing value for argument: " + arg);
            
            string value_str = argv[++i];
            
            if (def.type == "string") {
                results[name] = value_str;
            } else if (def.type == "int") {
                try {
                    results[name] = stoi(value_str);
                } catch (const invalid_argument&) {
                    throw runtime_error("Invalid integer value for argument: " + arg);
                }
            } else if (def.type == "float") {
                try {
                    results[name] = stof(value_str);
                } catch (const invalid_argument&) {
                    throw runtime_error("Invalid float value for argument: " + arg);
                }
            }
        }
    }
    
    for (const auto& [name, def] : definitions) {
        if (def.required && !results.count(name)) {
            throw runtime_error("Missing required argument: " + def.flag);
        }
    }
    
    return results;
}

//////////////////////////////////////////////////////// LCP Core Hashing ///////////////////////////////////////////////////////////

#define UBLOCK_BIT_SIZE 32
#define ROTATION_SHIFT (UBLOCK_BIT_SIZE - 1)
#define SPARSE_HASHING_CONSTANT 0x9e3779b9u

// Portable 64-bit MurmurHash3 implementation
// Returns h1 (low 64 bits) of MurmurHash3_x64_128
// Fixed: unaligned access, strict aliasing, endianness, overflow issues

static inline uint64_t rotl64(uint64_t x, int r) {
    return (x << r) | (x >> (64 - r));
}

static inline uint64_t fmix64(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

// Safe little-endian load with optional byteswap for big-endian
static inline uint64_t load_u64_le(const void *p) {
    uint64_t v;
    std::memcpy(&v, p, sizeof(v));
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    // Byteswap on big-endian systems
    v = __builtin_bswap64(v);
#elif defined(_WIN32) || defined(_WIN64)
    // Windows: assume little-endian (x86/x64)
    // Could add runtime check if needed
#endif
    return v;
}

// Returns h1 (low 64 bits) of MurmurHash3_x64_128
uint64_t MurmurHash3_64(const void *key, size_t len, uint64_t seed) {
    const uint8_t *data = (const uint8_t *)key;
    size_t nblocks = len / 16;  // Fixed: use size_t instead of int
    
    uint64_t h1 = seed;
    uint64_t h2 = seed;
    
    const uint64_t c1 = 0x87c37b91114253d5ULL;
    const uint64_t c2 = 0x4cf5ad432745937fULL;
    
    // Body: Process blocks of 16 bytes at a time
    // Fixed: use safe memcpy-based loads instead of direct pointer casting
    for (size_t i = 0; i < nblocks; i++) {
        uint64_t k1 = load_u64_le(data + i * 16 + 0);
        uint64_t k2 = load_u64_le(data + i * 16 + 8);
        
        k1 *= c1;
        k1 = rotl64(k1, 31);
        k1 *= c2;
        h1 ^= k1;
        
        h1 = rotl64(h1, 27);
        h1 += h2;
        h1 = h1 * 5 + 0x52dce729;
        
        k2 *= c2;
        k2 = rotl64(k2, 33);
        k2 *= c1;
        h2 ^= k2;
        
        h2 = rotl64(h2, 31);
        h2 += h1;
        h2 = h2 * 5 + 0x38495ab5;
    }
    
    // Tail: Process remaining bytes
    const uint8_t *tail = data + nblocks * 16;
    uint64_t k1 = 0;
    uint64_t k2 = 0;
    
    switch (len & 15) {
    case 15: k2 ^= ((uint64_t)tail[14]) << 48; [[fallthrough]];
    case 14: k2 ^= ((uint64_t)tail[13]) << 40; [[fallthrough]];
    case 13: k2 ^= ((uint64_t)tail[12]) << 32; [[fallthrough]];
    case 12: k2 ^= ((uint64_t)tail[11]) << 24; [[fallthrough]];
    case 11: k2 ^= ((uint64_t)tail[10]) << 16; [[fallthrough]];
    case 10: k2 ^= ((uint64_t)tail[9]) << 8; [[fallthrough]];
    case  9: k2 ^= ((uint64_t)tail[8]) << 0;
             k2 *= c2; k2 = rotl64(k2, 33); k2 *= c1; h2 ^= k2; [[fallthrough]];
    case  8: k1 ^= ((uint64_t)tail[7]) << 56; [[fallthrough]];
    case  7: k1 ^= ((uint64_t)tail[6]) << 48; [[fallthrough]];
    case  6: k1 ^= ((uint64_t)tail[5]) << 40; [[fallthrough]];
    case  5: k1 ^= ((uint64_t)tail[4]) << 32; [[fallthrough]];
    case  4: k1 ^= ((uint64_t)tail[3]) << 24; [[fallthrough]];
    case  3: k1 ^= ((uint64_t)tail[2]) << 16; [[fallthrough]];
    case  2: k1 ^= ((uint64_t)tail[1]) << 8; [[fallthrough]];
    case  1: k1 ^= ((uint64_t)tail[0]) << 0;
             k1 *= c1; k1 = rotl64(k1, 31); k1 *= c2; h1 ^= k1; [[fallthrough]];
    case  0: break;
    }
    
    // Finalization
    h1 ^= (uint64_t)len;
    h2 ^= (uint64_t)len;
    
    h1 += h2;
    h2 += h1;
    
    h1 = fmix64(h1);
    h2 = fmix64(h2);
    
    h1 += h2;
    
    return h1;  // Return only h1 (low 64 bits of x64_128)
}

static const uint32_t pow_table[33] = {
    0x00000001u,  // KBASE^0
    0x01000193u,  // KBASE^1
    0x26027A69u,  // KBASE^2
    0x3EE6B34Bu,  // KBASE^3
    0x502C3F11u,  // KBASE^4
    0x46A747C3u,  // KBASE^5
    0xFC55F7F9u,  // KBASE^6
    0x34555CFBu,  // KBASE^7
    0x5D615F21u,  // KBASE^8
    0x2148C0F3u,  // KBASE^9
    0x5887BE89u,  // KBASE^10
    0xE6B0F1ABu,  // KBASE^11
    0xD38C7031u,  // KBASE^12
    0x37149D23u,  // KBASE^13
    0xD8735E19u,  // KBASE^14
    0xD69D215Bu,  // KBASE^15
    0x345B8241u,  // KBASE^16
    0xAD0E0C53u,  // KBASE^17
    0xC01D66A9u,  // KBASE^18
    0x17489C0Bu,  // KBASE^19
    0xB24DA551u,  // KBASE^20
    0x013B3E83u,  // KBASE^21
    0x73436839u,  // KBASE^22
    0xAC1D11BBu,  // KBASE^23
    0xACC2E961u,  // KBASE^24
    0x57D563B3u,  // KBASE^25
    0xF7EBF2C9u,  // KBASE^26
    0x116F326Bu,  // KBASE^27
    0xDD0C5E71u,  // KBASE^28
    0x6B78ABE3u,  // KBASE^29
    0x11F69659u,  // KBASE^30
    0xA02EAE1Bu,  // KBASE^31
    0x447C1481u,  // KBASE^32
};

uint64_t hash_core(const struct core* c, const char* sequence) {
    if (c == NULL) return 0;

    // If we have the raw DNA sequence, hash it directly
    if (sequence != nullptr && c->end > c->start) {
        uint64_t length = c->end - c->start;
        const char* seq_start = sequence + c->start;
        
        // Hash the raw DNA sequence substring with MurmurHash3-64
        // Use label as seed for additional uniqueness
        uint64_t hash = MurmurHash3_64(seq_start, length, (uint64_t)c->label);
        return hash;
    }
    
    // Fallback: hash bit_rep if sequence not available
    ubit_size num_blocks = (c->bit_size + UBLOCK_BIT_SIZE - 1) / UBLOCK_BIT_SIZE;
    size_t byte_len = num_blocks * sizeof(ublock);
    uint64_t hash = MurmurHash3_64(c->bit_rep, byte_len, (uint64_t)c->label);
    return hash;
}